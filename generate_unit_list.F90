program generate_unit_list

!   This software is in the public domain in the U.S., per U.S. law, and is available under the
! CC0 1.0 license elsewhere in the world.  See LICENSE.md for the license that applies to this file.

!   This program provides a list of unit descriptions with a count of duplicated entries.  This list
! can be used as written to guide the selection of scaling factors, or it could be edited by a human
! to consolidate the list to avoid having multiple variables that scale the same way.

!   Although this program was developed to work with the unit description syntax used in the MOM6
! ocean model and the SIS2 sea-ice model, and some of this code was originally derived from MOM6
! code, it is a stand-alone program.

  implicit none

  !> A type used to create a list of the unit descriptions and the number of times they are used.
  type :: unit_entry
    character(len=128) :: unit_string  !< A unit description, such as [L T-1 ~> m s-1] for a velocity
    integer :: count = 0  !< The number of times this description has been used.
    integer, dimension(:), allocatable :: dim_powers  !< The powers of each of the rescaled units in this descriptor.
    type(unit_entry), pointer :: prev => NULL()  !< The previous entry in the double linked list.
    type(unit_entry), pointer :: next => NULL()  !< The next entry in the double linked list.
  end type unit_entry

  type(unit_entry), pointer :: list_start ! The starting point of the double-linked list
  type(unit_entry), pointer :: item  ! The entry to work on in the double-linked list
  type(unit_entry), pointer :: item_to_discard  ! The entry to remove from the double-linked list

  character(len=1024) :: arg, in_file, line, key_argument
  character(len=80) :: sample_key, unit_string, most_used_str
  character(len=8), dimension(:), allocatable :: key
  integer, dimension(:), allocatable :: dp_here
  integer :: file_arg, key_arg
  integer :: n, narg, arglen, nkey
  integer :: most_used_count, total_count
  integer :: inext, next_pos
  integer :: ios, iounit
  logical :: write_help, simplify

  !-----------------------------------------------------------------------

  ! Parse the command-line arguments
  narg = command_argument_count()
  file_arg = 1 ; key_arg = 0

  write_help = (narg < 1)
  simplify = .false.

  do n=1,narg
    call get_command_argument(n, arg, arglen)
    if (arglen > 1080) call write_errmsg("Excessively long input argument.")
    if (trim(arg) == "--help") write_help = .true.
    if (trim(arg) == "--simplify") simplify = .true.
    if (trim(arg) == "-f") then
      file_arg = n+1
      if (n >= narg) call write_errmsg("Missing file argument")
    endif
    if (trim(arg) == "-k") then
      key_arg = n+1
      if (n >= narg) call write_errmsg("Missing key argument")
    endif
  enddo

  if (narg < file_arg) call write_errmsg("Missing file argument")
  call get_command_argument(file_arg, in_file)

  if (write_help) then
    write(sample_key, '(A)') '["Z", "H", "L", ...]'
    write(*,*) "generate_unit_list [-k '"//trim(sample_key)//"'] [--simplify] [-f] input_file"
    stop
  endif

  if (key_arg < 1) then
    nkey = 8
    allocate(key(nkey), source="        ")
    key(1:nkey) = ["Z", "H", "L", "T", "R", "Q", "C", "S"]
  else
    call get_command_argument(key_arg, key_argument)
    call parse_key_argument(key_argument, key, nkey)
  endif

  open(newunit=iounit, file=trim(in_file), access='SEQUENTIAL', &
         form='FORMATTED', action='READ', position='REWIND', iostat=ios)

  allocate(list_start) ; list_start%count = 0 ; list_start%next => NULL() ; list_start%prev => NULL()

  ! Read the individual lines from the file and add any unit descriptions in them to the list,
  ! either by incrementing the count of similar descriptions or adding a new member of the list.
  do while (.true.)
    read(iounit, '(a)', end=8) line
    inext = 1

    do while (inext < len(line))
      call extract_unit_string(line(inext:), unit_string, next_pos)
      if (len_trim(unit_string) == 0) exit
      inext = inext + next_pos

      call increment_list(list_start, key, unit_string)

      ! For debugging: write(*,'(A)') trim(unit_string)
    enddo

  enddo ! while (.true.)
8 continue ! Jump to this point when read() reaches the end of the file.

  close(iounit)

  ! Uncomment this to write the entire list without consolidation.
  !  call write_list(list_start, key)

  ! Write the list of common terms grouped by common scaling factors, deallocating terms and
  ! contracting the list once they have been written
  allocate(dp_here(nkey), source=0)
  do while (associated(list_start))
    ! Identify the scaling factors that have the same rescaling as the top of the list.
    dp_here(:) = list_start%dim_powers(:)
    total_count = list_start%count
    most_used_str = list_start%unit_string
    most_used_count = list_start%count

    item => list_start%next
    do while (associated(item))
      ! Look for other terms with the same scaling, keeping track of the most frequently used
      ! factor or the one with the shortest description.
      if (same_powers(item%dim_powers, dp_here)) then
        total_count = total_count + item%count
        if ((most_used_count < item%count) .or. &
            ((most_used_count == item%count) .and. &
             (scan(item%unit_string, "~") < scan(most_used_str, "~")))) then
          most_used_str = item%unit_string
          most_used_count = item%count
        endif
      endif
      item => item%next
    enddo

    if (total_count > list_start%count) then
      if (simplify) then
        write(*, '(I4.4," ",A)') total_count, trim(most_used_str)
      else
        write(*, '("!==== Equivalent scaling to : ",I4.4," ",A)') total_count, trim(most_used_str)
      endif
      item => list_start
      list_start => NULL()
      do while (associated(item))
        if (same_powers(item%dim_powers, dp_here)) then
          if (.not.simplify) write(*, '(I4.4," ", A)') item%count, trim(item%unit_string)
          ! Discard this item.
          if (associated(item%next)) item%next%prev => item%prev
          if (associated(item%prev)) item%prev%next => item%next
          item_to_discard => item
          item => item%next
          deallocate(item_to_discard%dim_powers) ; deallocate(item_to_discard)
        else
          if (.not.associated(list_start)) list_start => item
          item => item%next
        endif
      enddo
    else ! There is only the one entry with this scaling.
      if (.not.simplify) write(*, '("!------------------------")')
      write(*, '(I4.4," ", A)') list_start%count, trim(list_start%unit_string)

      ! Discard this item.
      if (associated(list_start%next)) list_start%next%prev => list_start%prev
      if (associated(list_start%prev)) list_start%prev%next => list_start%next
      item_to_discard => list_start
      list_start => list_start%next
      deallocate(item_to_discard%dim_powers) ; deallocate(item_to_discard)
    endif

  enddo

  ! Clean up allocated memory.
  deallocate(dp_here, key)

contains

!> Return true if all of the powers in the dimension lists are the same.
logical function same_powers(dims1, dims2)
  integer, intent(in) :: dims1(:), dims2(:) !< The dimension powers to compare
  ! Local variables
  integer :: n

  same_powers = .true.
  do n=1,max(size(dims1),size(dims2))
    if (item%dim_powers(n) /= dp_here(n)) same_powers = .false.
  enddo
end function same_powers

!> Determine the powers to rescale from an argument provided as a run-time argument and store the array
!! of them in key, and also return the number of such dimensions.
subroutine parse_key_argument(key_argument, key, nkey)
  character(len=*), intent(in)  :: key_argument    !< The argument string to parse for the key
  character(len=*), dimension(:), allocatable, intent(inout)  :: key(:)  !< The key for the unit scaling
  integer, intent(out)  :: nkey     !< The number of rescaled units to categorize

  ! Local variables
  integer :: nquote, isk, iek, ioff, istart, n

  ! Count the quotes to determine the number of keys
  nquote = 0
  istart = 1
  do while (istart < len_trim(key_argument))
    ioff = scan(key_argument(istart:), '"')
    if (ioff > 0) then
      nquote = nquote + 1
      istart = istart + ioff
    else
      istart = istart + len_trim(key_argument)
    endif
  enddo

  nkey = nquote/2

  if (2*nkey /= nquote) call write_errmsg("Apparently mismatched quotes in key "//trim(key_argument))
  if (nkey == 0) call write_errmsg("No quotes found in key "//trim(key_argument))

  ! Populate the key.
  allocate(key(nkey), source="        ")
  istart = 1
  do n=1,nkey
    isk = scan(key_argument(istart:), '"')
    iek = scan(key_argument(istart+isk:), '"') + isk - 2
    if (iek < isk) call write_errmsg("Apparently blank key in "//trim(key_argument))
    if (iek >= isk+8) call write_errmsg("Apparently overly long key (limit 8 char) in "//trim(key_argument))
    key(n) = key_argument(istart+isk:istart+iek)
    istart = istart + iek + 2

    ! For debugging: write(*,'(A," from ",A)') trim(key(n)), trim(key_argument)
  enddo

end subroutine parse_key_argument


!> Add unit_string to the list or increment its count.
subroutine increment_list(list, key, unit_string)
  type(unit_entry), pointer :: list    !< A pointer to the head of a double-linked list.
  character(len=*), dimension(:), intent(in)  :: key(:)    !< The key for the unit scaling
  character(len=*), intent(in)  :: unit_string     !< The unit string to add to the list

  ! Local variables
  type(unit_entry), pointer :: item  ! The entry to work on in the double-linked list

  if (len_trim(unit_string) == 0) return

  item => list
  do while (.true.)
    if (trim(item%unit_string) == trim(unit_string)) then
      item%count = item%count + 1
      return
    elseif (item%count == 0) then  ! This only happens on the first calls
      item%unit_string = trim(unit_string) ; item%count = 1
      allocate(item%dim_powers(size(key)), source=0)
      call encode_dim_powers(unit_string, key, item%dim_powers)
      return
    elseif (associated(item%next)) then
      item => item%next
    else  ! Add a new entry to the double linked list.
      allocate(item%next)
      item%next%prev => item
      item => item%next
      allocate(item%dim_powers(size(key)), source=0)
      item%next => NULL()
      item%unit_string = trim(unit_string) ; item%count = 1
      call encode_dim_powers(unit_string, key, item%dim_powers)
      return
    endif
  enddo

  write (*,*) "Did nothing with entry "//trim(unit_string)
end subroutine increment_list

!> Add unit_string to the list or increment its count.
subroutine write_list(list, key)
  type(unit_entry), pointer :: list    !< A pointer to the head of a double-linked list.
  character(len=*), dimension(:), intent(in)  :: key(:)    !< The key for the unit scaling

  ! Local variables
  type(unit_entry), pointer :: item  ! The entry to work on in the double-linked list
  character(len=128) :: mesg, tmp_mesg

  item => list
  do while (associated(item))
    mesg = " "
    do n=1,size(key)
      write(tmp_mesg, '(" ",A,I2)') trim(key(n)), item%dim_powers(n)
      mesg = trim(mesg)//trim(tmp_mesg)
    enddo
    write(*, '(I4.4," ", A, " ", A)') item%count, trim(item%unit_string), trim(mesg)
    item => item%next
  enddo

end subroutine write_list


!> Return the first properly formatted unit string in a line fragment, along with the offset to the
!| next position to start looking for the next such string.
subroutine extract_unit_string(frag, unit_string, next_pos)
  character(len=*), intent(in)  :: frag     !< The fragment of the input line to be parsed.
  character(len=*), intent(out) :: unit_string !< A unit string extracted from the line or a blank
  integer,          intent(out) :: next_pos !< The offset to the next position to start looking
                                            !! for another unit string

  ! Local variables
  integer :: iconv, istart, iend

  iconv = index(frag, "~>")

  istart = 0 ; iend = 0
  if (iconv > 2) then
    istart = index(frag(:iconv-1), "[", back=.true.)
    iend = index(frag(iconv+2:), "]")
  endif

  if ((iconv==0) .or. (istart==0) .or. (iend==0)) then
    unit_string = ""
    next_pos = len_trim(frag)
  else
    unit_string = frag(istart:iconv+iend+1)
    next_pos = iconv+iend + 1
  endif

end subroutine extract_unit_string

!> Convert a unit scaling descriptor into an array of the dimensions of powers given in the key
subroutine encode_dim_powers(scaling, key, dim_powers)

  character(len=*),               intent(in)  :: scaling   !< The unit description that will be converted
  character(len=*), dimension(:), intent(in)  :: key(:)    !< The key for the unit scaling
  integer, dimension(size(key)),  intent(out) :: dim_powers !< The dimensions in scaling of each
                                                           !! element of they key.

  ! Local variables
  character(len=:), allocatable :: actstr ! The full active remaining string to be parsed.
  character(len=:), allocatable :: fragment ! The space-delimited fragment being parsed.
  character(len=:), allocatable :: dimnm  ! The probable dimension name
  character(len=11) :: numbers ! The list of characters that could make up the exponent.
  integer :: istart, iend, ieq, ief, ipow  ! Positions in strings.
  integer :: dp   ! The power for this dimension.
  integer :: ndim ! The number of dimensional scaling factors to consider.
  integer :: n

  dim_powers(:) = 0

  iend = index(scaling, "~>") - 1
  if (iend < 1) return

  ! Parse the key.
  ndim = size(key)
  numbers = "-0123456789"

  ! Strip away any leading square brace.
  istart = index(scaling(:iend), "[") + 1
  ! If there is an "=" in the string, start after this.
  ieq = index(scaling(istart:iend), "=", back=.true.)
  if (ieq > 0) istart = istart + ieq

  ! Set up the active string to work on.
  actstr = trim(adjustl(scaling(istart:iend)))
  do  ! Loop over each of the elements in the unit scaling descriptor.
    if (len_trim(actstr) == 0) exit
    ief = index(actstr, " ") - 1
    if (ief <= 0) ief = len_trim(actstr)
    fragment = actstr(:ief)
    ipow = scan(fragment, "-")
    if (ipow == 0) ipow = scan(fragment, numbers)

    if (ipow == 0) then ! There is no exponent
      dimnm = fragment
      dp = 1
    else
      if (verify(fragment(ipow:), numbers) == 0) then
        read(fragment(ipow:),*) dp
        dimnm = fragment(:ipow-1)
      else
        dimnm = fragment
        dp = 1
      endif
    endif

    do n=1,ndim
      if (trim(dimnm) == trim(key(n))) then
        dim_powers(n) = dim_powers(n) + dp
        exit
      endif
    enddo

    ! Remove the leading fragment that has been parsed from actstr
    actstr = trim(adjustl(actstr(ief+1:)))
  enddo

end subroutine encode_dim_powers

!> This subroutine writes an error message and stops the run.
subroutine write_errmsg(mesg)
  character(len=*), intent(in)  :: mesg !< The error message to write

  write(*,*) "Error: "//trim(mesg)
  stop

end subroutine write_errmsg

!> \namespace generate_unit_list
!!
!!   The stand-alone tool generate_unit_list in this file extacts the rescaled unit
!! descriptions in the syntax like [L T-1 ~> m s-1] from the file given by a
!! command line argument, and writes to stdout a list of the various units that have
!! distinct dimensional scaling and the counts of their occurance, accoording to a
!! key of rescaled variables that can also be provided as a run-time parameter.
!!
!!   The syntax for generate_unit_list is:
!!
!!    generate_unit_list [-k <key>] [--simplify] [-f] <input_file>
!!
!!  <input_file> is the name of the file to parse for the unit descriptions.
!!
!!  <key> is the key of rescaled units in the double-quote delimited syntax, but
!!    if omitted the key defaults to the full list currently used by MOM6:
!!    '["Z", "H", "L", "T", "R", "Q", "C", "S"]'
!!
!!  --simplify specifies that all units with common rescaling are collapsed to a
!!    single entry giving the most frequently used unit description.  If it is
!!    omitted all unique unit descriptions are listed, but grouped by common
!!    rescaling and with added lines separating groups with different rescaling to
!!    aid in the manual consolidation of descriptions with common rescaling.
!!

end program generate_unit_list
