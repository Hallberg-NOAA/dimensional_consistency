# dimensional_consistency
This repository contains standalone tools that help to manage the dimensional consistency testing as is carried out in the
MOM6 ocean model (github.com/mom-ocean/MOM6) and the SIS2 sea ice model (https://github.com/NOAA-GFDL/SIS2).

`generate_unit_list()` extacts the rescaled unit descriptions in the syntax like [L T-1 ~> m s-1] from the file given by a
command line argument, and writes to stdout a list of the various units that have distinct dimensional scaling and the
counts of their occurance, accoording to a key of rescaled variables that can also be provided as a run-time parameter.
This output can in turn be used to identify rescaling factors that avoid aliasing different units together, and hence
allow for a single test to evaluate the dimensional consitency of the MOM6 or SIS2 code for multiple dimensions at once.

   The syntax for generate_unit_list is:

    generate_unit_list [-k <key>] [--simplify] [-f] <input_file>

  <input_file> is the name of the file to parse for the unit descriptions.

  <key> is the key of rescaled units in the double-quote delimited syntax, but
      if omitted the key defaults to the full list currently used by MOM6:
      '["Z", "H", "L", "T", "R", "Q", "C", "S"]'

  --simplify specifies that all units with common rescaling are collapsed to a
      single entry giving the most frequently used unit description.  If it is
      omitted all unique unit descriptions are listed, but grouped by common
      rescaling and with added lines separating groups with different rescaling to
      aid in the manual consolidation of descriptions with common rescaling.
