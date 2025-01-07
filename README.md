## Version of LM4.1 with GLASS snow model

Built based on ppa_20220812_litt_slm

Recent version of LM4.1 with CENTURY soil carbon

With GLASS v1.0 snow model

Search for tag EZSNOW-2022SC for changes within the snow module necessary for back compatibility
The entire snowlayers_io.F90 and snowlayers_io.inc have been changed


# Note on restart

The restart for GLASS snow was developed for an older version of the model
Here it was updated

# Description of GLASS snow namelists

# snowpack_nml
# snicar_nml
# snow_data_nml
# snow_evolution_nml
# cm_snow_nml


There is a change in namelist: use_mcm_masking, depth_crit, csw, clw, cpw are now in snowpack_nml

# version Jan 6 2025:
# branching out from 
# add an albedo modification for snow over glaciers, similar to the albedo change by Chris Milly
# decrease snow NIR albedo of the same amount as in Milly snow scheme.

