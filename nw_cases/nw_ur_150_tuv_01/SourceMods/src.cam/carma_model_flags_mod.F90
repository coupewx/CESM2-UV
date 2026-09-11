!! This module handles reading the namelist and provides access to some other flags
!! that control a specific CARMA model's behavior.
!!
!! By default the specific CARMA model does not have any unique namelist values. If
!! a CARMA model wishes to have its own namelist, then this file needs to be copied
!! from physics/cam to physics/model/<model_name> and the code needed to read in the
!! namelist values added there. This file will take the place of the one in
!! physics/cam. 
!!
!! It needs to be in its own file to resolve some circular dependencies.
!!
!! @author  Chuck Bardeen
!! @version Mar-2011
module carma_model_flags_mod

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use spmd_utils,     only: masterproc

  ! Flags for integration with CAM Microphysics
  public carma_model_readnl                   ! read the carma model namelist
  

  ! Namelist flags
  !
  ! Create a public definition of any new namelist variables that you wish to have,
  ! and default them to an inital value.
  real(r8), public               :: carma_emis_dust           = 0._r8     !! Total dust emission for the event (kg)
  real(r8), public               :: carma_emis_soot           = 0._r8     !! Total soot emission for the event (kg)
  real(r8), public               :: carma_emis_co2_splash     = 0._r8     !! Splashed co2 emission for the event (kg)
  real(r8), public               :: carma_emis_co2_fire       = 0._r8     !! Fire generated co2 emission for the event (kg)
  real(r8), public               :: carma_emis_co2_impact     = 0._r8     !! Ballistic co2 emission for the event (kg)
  real(r8), public               :: carma_emis_h2o_splash     = 0._r8     !! Splashed h2o emission for the event (kg)
  real(r8), public               :: carma_emis_h2o_fire       = 0._r8     !! Fire generated h2o emission for the event (kg)
  real(r8), public               :: carma_emis_h2o_impact     = 0._r8     !! Ballistic h2o emission for the event (kg)
  real(r8), public               :: carma_emis_so2_splash     = 0._r8     !! Splashed so2 emission for the event (kg)
  real(r8), public               :: carma_emis_so2_fire       = 0._r8     !! Fire generated so2 emission for the event (kg)
  real(r8), public               :: carma_emis_so2_impact     = 0._r8     !! Ballistic so2 emission for the event (kg)
  real(r8), public               :: carma_emis_hbr_splash     = 0._r8     !! Splashed hbr emission for the event (kg)
  real(r8), public               :: carma_emis_hbr_fire       = 0._r8     !! Fire generated hbr emission for the event (kg)
  real(r8), public               :: carma_emis_hbr_impact     = 0._r8     !! Ballistic hbr emission for the event (kg)
  real(r8), public               :: carma_emis_hcl_splash     = 0._r8     !! Splashed hcl emission for the event (kg)
  real(r8), public               :: carma_emis_hcl_fire       = 0._r8     !! Fire generated hcl emission for the event (kg)
  real(r8), public               :: carma_emis_hcl_impact     = 0._r8     !! Ballistic hcl emission for the event (kg)
  real(r8), public               :: carma_emis_heat_fire      = 0._r8     !! Fire generated heat emission for the event (J)
  real(r8), public               :: carma_emis_no_splash      = 0._r8     !! Splashed no emission for the event (kg)
  real(r8), public               :: carma_emis_no_fire        = 0._r8     !! Fire generated no emission for the event (kg)
  real(r8), public               :: carma_emis_no_impact      = 0._r8     !! Ballistic no emission for the event (kg)  
  real(r8), public               :: carma_emis_ctrlat         = 0._r8     !! impact latitude
  real(r8), public               :: carma_emis_ctrlon         = 0._r8     !! impact longitude
  real(r8), public               :: carma_emis_dust_radius    = 20000._r8 !! distance from impact affected by dust (km)
  real(r8), public               :: carma_emis_soot_radius    = 20000._r8 !! distance from impact affected by soot (km)
  real(r8), public               :: carma_emis_splash_radius  = 20000._r8 !! distance from impact affected by soot (km)
  real(r8), public               :: carma_emis_fine_trop_frac = 0.5_r8    !! fraction of fine soot placed at tropopasue (fraction)
  real(r8), public               :: carma_emis_coarse_trop_frac = 0.5_r8  !! fraction of coarse soot placed at tropopasue (fraction)
  real(r8), public               :: carma_emis_fire_fine_frac = 0.5_r8    !! fraction of soot placed at tropopasue (fraction)
  real(r8), public               :: carma_emis_fire_trop_frac = 0.5_r8    !! fraction of soot placed at tropopasue (fraction)
  integer, public                :: carma_emis_startdate      = 1     !! start year and day of year (yyyyddd)
  integer, public                :: carma_emis_stopdate       = 1     !! stop year and day of year (yyyyddd)
  integer, public                :: carma_emis_starttime      = 0     !! start time of day (s)
  integer, public                :: carma_emis_stoptime       = 0     !! stop time of day (s)
  integer, public                :: carma_emis_dust_dtime     = 0     !! duration of dust and gas emissions (s)
  integer, public                :: carma_emis_z_ballistic    = 50.   !! altitude of vaporized dust and gas emissions (km)
  integer, public                :: carma_emis_heat_fire_levs = 1     !! number of surface levels for fire heating.
  logical, public                :: carma_fractal_soot        = .false. !! fractal Soot
  logical, public                :: carma_emis_fire_land_only = .false.   !! emit from fires only over land
  logical, public                :: carma_emis_fire_veg_scale = .false.   !! scale fires by vegetation
  character(len=256), public     :: carma_mdust_file          = "index_Mg01Fe09O.nc" !! dust refractive index file

contains


  !! Read the CARMA model runtime options from the namelist
  !!
  !! @author  Chuck Bardeen
  !! @version Mar-2011
  subroutine carma_model_readnl(nlfile)
  
    ! Read carma namelist group.
  
    use cam_abortutils,  only: endrun
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
  
    ! args
  
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
  
    ! local vars
  
    integer :: unitn, ierr
  
    ! read namelist for CARMA
    namelist /carma_model_nl/ &
      carma_emis_dust, carma_emis_soot, carma_emis_h2o_splash, carma_emis_h2o_fire, carma_emis_h2o_impact, &
      carma_emis_co2_splash, carma_emis_co2_fire, carma_emis_co2_impact, &
      carma_emis_so2_splash, carma_emis_so2_fire, carma_emis_so2_impact, &
      carma_emis_co2_splash, carma_emis_co2_fire, carma_emis_co2_impact, &
      carma_emis_no_splash, carma_emis_no_fire, carma_emis_no_impact, &
      carma_emis_hbr_splash, carma_emis_hbr_fire, carma_emis_hbr_impact, &
      carma_emis_hcl_splash, carma_emis_hcl_fire, carma_emis_hcl_impact, &
      carma_emis_no_splash, carma_emis_no_fire, carma_emis_no_impact, &
      carma_emis_heat_fire, carma_emis_heat_fire_levs, &
      carma_emis_startdate, carma_emis_stopdate, carma_emis_dust_dtime, &
      carma_emis_starttime, carma_emis_stoptime, carma_emis_ctrlat, carma_emis_ctrlon, &
      carma_emis_dust_radius, carma_emis_soot_radius, carma_emis_splash_radius, carma_fractal_soot, &
      carma_emis_fire_trop_frac,carma_emis_z_ballistic, carma_mdust_file, carma_emis_fine_trop_frac, carma_emis_coarse_trop_frac, &
      carma_emis_fire_fine_frac, carma_emis_fire_land_only, carma_emis_fire_veg_scale
      
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'carma_model_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, carma_model_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('carma_model_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if
  
#ifdef SPMD
    call mpibcast(carma_emis_dust,           1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_soot,           1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_co2_splash,     1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_co2_fire,       1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_co2_impact,     1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_h2o_splash,     1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_h2o_fire,       1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_h2o_impact,     1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_hbr_splash,     1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_hbr_fire,       1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_hbr_impact,     1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_hcl_splash,     1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_hcl_fire,       1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_hcl_impact,     1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_so2_splash,     1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_so2_fire,       1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_so2_impact,     1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_heat_fire,      1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_no_splash,      1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_no_fire,        1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_no_impact,      1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_startdate,      1,  mpiint,  0, mpicom)
    call mpibcast(carma_emis_stopdate,       1,  mpiint,  0, mpicom)
    call mpibcast(carma_emis_starttime,      1,  mpiint,  0, mpicom)
    call mpibcast(carma_emis_stoptime,       1,  mpiint,  0, mpicom)
    call mpibcast(carma_emis_dust_dtime,     1,  mpiint,  0, mpicom)
    call mpibcast(carma_emis_heat_fire_levs, 1,  mpiint,  0, mpicom)
    call mpibcast(carma_emis_ctrlat,         1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_ctrlon,         1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_dust_radius,    1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_soot_radius,    1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_splash_radius,  1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_z_ballistic,    1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_fine_trop_frac, 1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_coarse_trop_frac, 1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_fire_fine_frac, 1,  mpir8,   0, mpicom)
    call mpibcast(carma_emis_fire_trop_frac, 1,  mpir8,   0, mpicom)
    call mpibcast(carma_fractal_soot,        1,  mpilog,  0, mpicom)
    call mpibcast(carma_emis_fire_land_only, 1,  mpilog,  0, mpicom)
    call mpibcast(carma_emis_fire_veg_scale, 1,  mpilog,  0, mpicom)
    call mpibcast(carma_mdust_file,          len(carma_mdust_file),  mpichar,   0, mpicom)
#endif
  
  end subroutine carma_model_readnl

end module carma_model_flags_mod
