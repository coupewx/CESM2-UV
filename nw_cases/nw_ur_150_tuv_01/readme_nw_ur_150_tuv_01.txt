
# This document contains info about how the original case was created on CU Research Computing's blanca.
# Note that you will need to complete porting steps to run CESM on your own machine.


# Copy all files in this folder to the CASEDIR created from ./create_newcase to be able to run.


cd /home/joco6825/cesm_tags/cesm2.1.5-UVphyto/cime/scripts/


./create_newcase --compset 2000_CAM60%WCCM_CLM50%BGC-CROP_CICE_POP2%ECO_MOSART_CISM2%NOEVOLVE_WW3_BGC%BDRD --res f19_g17 --machine blanca-curc --run-unsupported --case /projects/joco6825/nw_cases/nw_ur_150_tuv_01 --output-root /scratch/alpine/joco6825/cesm/nw_cases/ --project blanca-toon


cd /projects/joco6825/nw_cases/nw_ur_150_tuv_01


./case.setup



# Update CAM_CONFIG_OPTS in env_build.xml to:

CAM_CONFIG_OPTS: -phys cam4 -rad rrtmg -age_of_air_trcs -chem waccm_ma_sulfur -nadv 172 -usr_mech_infile /projects/joco6825/nw_cases/nw_ur_150_tuv_01/SourceMods/src.cam/chem_mech.in -carma meteor_impact_all -max_n_rad_cnst=172


# Copy restart files to RUNDIR (restart files are available on data archive)

cp -r /scratch/alpine/joco6825/cesm/nw_cases/archive/nw_cntrl_tuv_01/rest/0005-01-01-00000/* /scratch/alpine/joco6825/cesm/nw_cases/nw_ur_150_tuv_01/run/


# Copy aerosol property files to RUNDIR (files are available on data archive)

cp -r /projects/joco6825/nw_cases/nw_ip_46p8_tuv_ens01/aerosol_files/combined/* /scratch/alpine/joco6825/cesm/nw_cases/nw_ur_150_tuv_01/run/


./xmlchange RUN_REFCASE=nw_cntrl_tuv_01
./xmlchange RUN_REFDATE=0005-01-01
./xmlchange RUN_STARTDATE=0005-01-01
./xmlchange CONTINUE_RUN=FALSE
./xmlchange RUN_TYPE=hybrid
./xmlchange STOP_N=1
./xmlchange STOP_OPTION=nmonths

./xmlchange GET_REFCASE=FALSE



# Other notes:

user_nl_cam must have the following for emission of 150 Tg black carbon on May 15th

carma_emis_soot                = 1.5e11
carma_emis_startdate           = 05135
carma_emis_stopdate            = 05142
carma_emis_starttime           = 0
carma_emis_stoptime            = 0
carma_fractal_soot             = .true.




./case.build --skip-provenance-check

./case.submit
