!! This module is used to define a particular CARMA microphysical model. For 
!! simple cases, this may be the only code that needs to be modified. This module
!! defines several constants and has three methods:
!!
!!   - CARMA_DefineModel()
!!   - CARMA_EmitParticle()
!!   - CARMA_InitializeParticle()
!!
!! These methods define the microphysical model, the particle emissions and
!! the initial conditions of the particles. Each realization of CARMA
!! microphysics has its own version of this file.
!!
!! This file is used to model a meteor impact upon the land. This model is
!! preliminary. Please talk to Chuck Bardeen (bardeenc@ucar.edu) if you are
!! interested in this model.
!!
!! @version Oct-2012 
!! @author  Chuck Bardeen 
module carma_model_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmaelement_mod
  use carmagas_mod
  use carmagroup_mod
  use carmasolute_mod
  use carmastate_mod
  use carma_mod
  use carma_flags_mod
  use carma_model_flags_mod
  
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun
  use physics_types,  only: physics_state, physics_ptend
  use pmgrid,         only: plat, plon
  use ppgrid,         only: pcols, pver
  use physics_buffer, only: physics_buffer_desc

  implicit none

  private

  ! Declare the public methods.
  public CARMA_DefineModel
  public CARMA_Detrain
  public CARMA_DiagnoseBins
  public CARMA_DiagnoseBulk
  public CARMA_EmitGas
  public CARMA_EmitHeat
  public CARMA_EmitParticle
  public CARMA_InitializeModel
  public CARMA_InitializeParticle
  public CARMA_WetDeposition
  
  ! Declare public constants
  integer, public, parameter      :: NGROUP   = 4               !! Number of particle groups
  integer, public, parameter      :: NELEM    = 4               !! Number of particle elements
  integer, public, parameter      :: NBIN     = 21              !! Number of particle bins
  integer, public, parameter      :: NSOLUTE  = 0               !! Number of particle solutes
  integer, public, parameter      :: NGAS     = 7               !! Number of gases


  !! Relative humidities for mie and radiation calculations. The RRTMG radiation code will interpolate
  !! based upon the current relative humidity from a table built using the specified relative
  !! humidities.
  integer, public, parameter      :: NMIE_RH  = 8               !! Number of relative humidities for mie calculations
  real(kind=f), public            :: mie_rh(NMIE_RH) = (/ 0._f, 0.5_f, 0.7_f, 0.8_f, 0.9_f, 0.95_f, 0.98_f, 0.99_f /)
  
  ! Defines whether the groups should undergo deep convection in phase 1 or phase 2.
  ! Water vapor and cloud particles are convected in phase 1, while all other constituents
  ! are done in phase 2.
  logical, public                 :: is_convtran1(NGROUP) = .false.  !! Should the group be transported in the first phase?

  ! Define any particle compositions that are used. Each composition type
  ! should have a unique number.
  integer, public, parameter      :: I_DUST         = 1         !! dust composition
  integer, public, parameter      :: I_SOOT         = 2         !! soot composition
  integer, public, parameter      :: I_H2SO4        = 3         !! sulfate aerosol composition
  integer, public, parameter      :: I_WATER        = 4         !! water

  ! Define group, element, solute and gas indexes.
  integer, public, parameter      :: I_GRP_DUST     = 1         !! dust aerosol group
  integer, public, parameter      :: I_GRP_CSOOT    = 2         !! coarse soot aerosol group
  integer, public, parameter      :: I_GRP_SULFATE  = 3         !! sulfate aerosol
  integer, public, parameter      :: I_GRP_SOOT     = 4         !! fine soot aerosol group

  integer, public, parameter      :: I_ELEM_DUST    = 1         !! dust aerosol element
  integer, public, parameter      :: I_ELEM_CSOOT   = 2         !! coarse aerosol element
  integer, public, parameter      :: I_ELEM_SULFATE = 3         !! sulfate aerosol
  integer, public, parameter      :: I_ELEM_SOOT    = 4         !! fine soot aerosol

  integer, public, parameter      :: I_GAS_H2O       = 1        !! water vapor
  integer, public, parameter      :: I_GAS_H2SO4     = 2        !! sulphuric acid
  integer, public, parameter      :: I_GAS_SO2       = 3        !! sulfur dioxide
  integer, public, parameter      :: I_GAS_CO2       = 4        !! carbon dioxide
  integer, public, parameter      :: I_GAS_HBR       = 5        !! hydrobromic acid
  integer, public, parameter      :: I_GAS_HCL       = 6        !! hydrochloric acid
  integer, public, parameter      :: I_GAS_NO        = 7        !! nitric oxide

  real(kind=f), public, parameter :: WTMOL_H2SO4    = 98.078479_f    !! molecular weight of sulphuric acid  
  real(kind=f), public, parameter :: WTMOL_HBR      = 80.91_f        !! molecular weight of hydrobromic acid  
  real(kind=f), public, parameter :: WTMOL_HCL      = 36.46094_f     !! molecular weight of hydrochloric acid  
  real(kind=f), public, parameter :: WTMOL_NO       = 30.01_f        !! molecular weight of nitric oxide  

  ! Physics buffer index for sulfate surface area density
  integer                         :: ipbuf4sad, ipbuf4reff, ipbuf4so4mmr
  
  integer                         :: carma_dustmap(NBIN)        !! mapping of the CARMA dust bins to the surface dust bins.
  real(kind=f)                    :: carma_sootbinfactor(NBIN)  !! bin weighting factor for soot emissions
  real(kind=f)                    :: carma_csootbinfactor(NBIN) !! bin weighting factor for coarse soot emissions
  real(kind=f)                    :: carma_emis_dtime           !! duration of the fire event (s)
  real(kind=f)                    :: carma_emis_soot_area       !! surface area where soot emissions are happening (m2)
  real(kind=f)                    :: carma_emis_dust_area       !! surface area where dust emissions are happening (m2)
  real(kind=f)                    :: carma_emis_splash_area     !! surface area where splash emissions are happening (m2)
  complex(kind=f)                 :: carma_refidx_dust(NWAVE)   !! dust refractive indicies
  real(kind=f)                    :: carma_emis_tot_veg         !! total vegetation fraction where soot emissions are happening, max of 1 if all veg burning
  real(kind=f)                    :: carma_vegfrac_map(plon, plat) !! fraction of total vegetation per area in a gridbox (1/m2)


  ! At some point, these should be moved to CARMA flags and added to the namelist and the
  ! number of regions should be increased.
  integer, parameter              :: maxRegions = 2             !! Maximum number of regions that can be defined
  real(r8)     :: carma_emis_rgnFraction(maxRegions) = (/ 0.45_r8, 0.55_r8 /)
  real(r8)     :: carma_emis_rgnMinLat(maxRegions)   = (/ 28._r8, 48._r8 /)
  real(r8)     :: carma_emis_rgnMaxLat(maxRegions)   = (/ 56._r8, 76._r8 /)
  real(r8)     :: carma_emis_rgnMinLon(maxRegions)   = (/ -122.5_r8, 27.5_r8 /)
  real(r8)     :: carma_emis_rgnMaxLon(maxRegions)   = (/ -72.5_r8, 157.5_r8 /)
  real(r8)     :: carma_emis_rgnMinP(maxRegions)     = (/ 150._r8, 150._r8 /)
  real(r8)     :: carma_emis_rgnMaxP(maxRegions)     = (/ 300._r8, 300._r8 /)

  ! This part stays here.
  integer                         :: nRegions   = 0             !! Number of regions that are defined.
  real(r8)     :: carma_emis_rgnArea(maxRegions)                !! Emission surface area within the region. 


contains


  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DefineModel(carma, rc)
    use ioFileMod,                     only: getfil
    use wrap_nf
    use mpishorthand
    use physics_buffer, only: pbuf_add_field, dtype_r8

    type(carma_type), intent(inout)    :: carma     !! the carma object
    integer, intent(out)               :: rc        !! return code, negative indicates failure
    
    ! Local variables
    real(kind=f), parameter            :: RHO_DUST = 2.7_f     ! density of dust particles (g/cm)
    real(kind=f)                       :: RHO_SOOT             ! density of soot particles (g/cm)
    real(kind=f)                       :: RHO_CSOOT = 1.0_f    ! density of coarse soot particles (g/cm)
    real(kind=f), parameter            :: RHO_SULFATE = 1.923_f    ! dry density of sulfate particles (g/cm3)
    real(kind=f), parameter            :: dust_rmin = 20.e-7_f ! dust minimum radius (cm)
    real(kind=f), parameter            :: dust_vmrat = 2.49_f  ! dust volume ratio
    real(kind=f), parameter            :: soot_rmin = 20.e-7_f ! soot minimum radius (cm)
    real(kind=f), parameter            :: soot_vmrat = 2.49_f  ! soot volume ratio
    real(kind=f), parameter            :: csoot_rmin = 1.e-4_f ! coarse soot minimum radius (cm)
    real(kind=f), parameter            :: csoot_vmrat = 2._f   ! coarse soot volume ratio
    !  Set radius of smallest bin such that mass is that of 2 molecules of H2SO4.   
    real(kind=f), parameter            :: so4_rmin = 3.43230298e-8_f  ! sulfate minimum radius (cm)
    real(kind=f), parameter            :: so4_vmrat = 5.5_f    ! sulfate volume ratio
    complex(kind=f)                    :: refidx(NWAVE)        ! refractice indices

    integer                            :: LUNOPRT               ! logical unit number for output
    logical                            :: do_print              ! do print output?
    real(kind=f)                       :: soot_rmon = 30.e-7_f  ! soot monomer radius (cm)
    real(kind=f)                       :: soot_df(NBIN) = 2.2_f ! soot fractal dimension
    real(kind=f)                       :: soot_falpha = 1._f    ! soot fractal packing coefficient

    integer                            :: i
    integer                            :: j
    integer                            :: irgn
    real(kind=f)                       :: wave(NWAVE)           ! CAM band wavelength centers (cm)
    integer                            :: fid
    integer                            :: wave_did
    integer                            :: wave_vid
    integer                            :: real_vid
    integer                            :: imag_vid
    character(len=256)                 :: efile                 ! vaporized dust refractive index file name
    real(kind=f)                       :: interp
    integer                             :: dust_nwave          ! number of dust wavelengths in file
    real(r8), allocatable, dimension(:) :: dust_wave           ! dust wavelengths
    real(r8), allocatable, dimension(:) :: dust_real           ! dust, real part of m
    real(r8), allocatable, dimension(:) :: dust_imag           ! dust, imag part of m

    ! Default return code.
    rc = RC_OK



    ! Determine how many emission regions are defined.
    nRegions = 0
    do irgn = 1, maxRegions
      if (carma_emis_rgnMinLat(irgn) .ne. carma_emis_rgnMaxLat(irgn)) then
        nRegions = nRegions + 1
      end if
    end do


    ! Report model specific namelist configuration parameters.
    if (masterproc) then
      call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT, wave=wave)
      if (rc < 0) call endrun("CARMA_InitializeModel: CARMA_Get failed.")
    
      if (do_print) write(LUNOPRT,*) ''
      if (do_print) write(LUNOPRT,*) 'CARMA ', trim(carma_model), ' specific settings :'
      if (do_print) write(LUNOPRT,*) '  carma_emis_startdate      = ', carma_emis_startdate
      if (do_print) write(LUNOPRT,*) '  carma_emis_starttime      = ', carma_emis_starttime
      if (do_print) write(LUNOPRT,*) '  carma_emis_stopdate       = ', carma_emis_stopdate
      if (do_print) write(LUNOPRT,*) '  carma_emis_stoptime       = ', carma_emis_stoptime
      if (do_print) write(LUNOPRT,*) '  carma_emis_ctrlat         = ', carma_emis_ctrlat
      if (do_print) write(LUNOPRT,*) '  carma_emis_ctrlon         = ', carma_emis_ctrlon
      if (do_print) write(LUNOPRT,*) '  carma_emis_dust           = ', carma_emis_dust, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_soot           = ', carma_emis_soot, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_co2_fire       = ', carma_emis_co2_fire, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_co2_impact     = ', carma_emis_co2_impact, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_co2_splash     = ', carma_emis_co2_splash, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_h2o_fire       = ', carma_emis_h2o_fire, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_h2o_impact     = ', carma_emis_h2o_impact, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_h2o_splash     = ', carma_emis_h2o_splash, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_hbr_fire       = ', carma_emis_hbr_fire, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_hbr_impact     = ', carma_emis_hbr_impact, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_hbr_splash     = ', carma_emis_hbr_splash, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_hcl_fire       = ', carma_emis_hcl_fire, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_hcl_impact     = ', carma_emis_hcl_impact, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_hcl_splash     = ', carma_emis_hcl_splash, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_no_fire        = ', carma_emis_no_fire, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_no_impact      = ', carma_emis_no_impact, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_no_splash      = ', carma_emis_no_splash, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_so2_fire       = ', carma_emis_so2_fire, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_so2_impact     = ', carma_emis_so2_impact, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_so2_splash     = ', carma_emis_so2_splash, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_heat_fire      = ', carma_emis_heat_fire, ' (kg)'
      if (do_print) write(LUNOPRT,*) '  carma_emis_soot_radius    = ', carma_emis_soot_radius
      if (do_print) write(LUNOPRT,*) '  carma_emis_dust_radius    = ', carma_emis_dust_radius
      if (do_print) write(LUNOPRT,*) '  carma_emis_splash_radius  = ', carma_emis_splash_radius
      if (do_print) write(LUNOPRT,*) '  carma_emis_fire_fine_frac   = ', carma_emis_fire_fine_frac
      if (do_print) write(LUNOPRT,*) '  carma_emis_fine_trop_frac   = ', carma_emis_fine_trop_frac
      if (do_print) write(LUNOPRT,*) '  carma_emis_coarse_trop_frac = ', carma_emis_coarse_trop_frac
      if (do_print) write(LUNOPRT,*) '  carma_emis_fire_trop_frac = ', carma_emis_fire_trop_frac
      if (do_print) write(LUNOPRT,*) '  carma_emis_fire_land_only = ', carma_emis_fire_land_only
      if (do_print) write(LUNOPRT,*) '  carma_emis_fire_veg_scale = ', carma_emis_fire_veg_scale
      if (do_print) write(LUNOPRT,*) '  carma_fractal_soot        = ', carma_fractal_soot
      if (do_print) write(LUNOPRT,*) '  carma_mdust_file          = ', carma_mdust_file
      if (do_print) write(LUNOPRT,*) '  carma_hetchem_feedback    = ', carma_hetchem_feedback
    
      if (do_print) then
        do irgn = 1, nRegions
          write(LUNOPRT,*) ''
          write(LUNOPRT,*) '  region                    = ', irgn
          write(LUNOPRT,*) '    carma_emis_rgnFraction    = ', carma_emis_rgnFraction(irgn)
          write(LUNOPRT,*) '    carma_emis_rgnMinLat      = ', carma_emis_rgnMinLat(irgn)
          write(LUNOPRT,*) '    carma_emis_rgnMaxLat      = ', carma_emis_rgnMaxLat(irgn)
          write(LUNOPRT,*) '    carma_emis_rgnMinLon      = ', carma_emis_rgnMinLon(irgn)
          write(LUNOPRT,*) '    carma_emis_rgnMaxLon      = ', carma_emis_rgnMaxLon(irgn)
          write(LUNOPRT,*) '    carma_emis_rgnMinP        = ', carma_emis_rgnMinP(irgn)
          write(LUNOPRT,*) '    carma_emis_rgnMaxP        = ', carma_emis_rgnMaxP(irgn)
          write(LUNOPRT,*) ''
        end do
      end if


    end if

    ! Define the Groups
    !
    ! NOTE: If NWAVE > 0 then the group should have refractive indices defined.
    !
    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.

    ! Use the same refractive index at all wavelengths. This value is typical of soot and
    ! is recommended by Toon et al. 2012. TBD Wagner et al. 2011 shows variability in the
    ! real part (0.003 (IR) to 0.05 (UV)).
;    refidx(:) = (1.53_f, 0.008_f)
        
    ! Use Mg0.1Fe0.9 for dust refractive indicies from Hervig et al. [2009] and personal
    ! communication. Vaporized dust is mostly from the impactor. 
    !
    ! NOTE: These values probably should be a band average, but for now just do band centers.
    if (masterproc) then 
    
      ! Open the netcdf file (read only)
      call getfil(carma_mdust_file, efile, fid)
      if (do_print) write(LUNOPRT,*) 'carma_init(): Reading dust refractive indexes from ', efile

      call wrap_open(efile, 0, fid)

      ! Alocate the table arrays
      call wrap_inq_dimid(fid, "wvlen", wave_did)
      call wrap_inq_dimlen(fid, wave_did, dust_nwave)
  
      allocate(dust_wave(dust_nwave))
      allocate(dust_real(dust_nwave))
      allocate(dust_imag(dust_nwave))
  
      ! Read in the tables.
      call wrap_inq_varid(fid, 'wvlen', wave_vid)
      call wrap_get_var_realx(fid, wave_vid, dust_wave)
      dust_wave = dust_wave * 1e-4_f          ! um -> cm

      call wrap_inq_varid(fid, 'm_r', real_vid)
      call wrap_get_var_realx(fid, real_vid, dust_real)

      call wrap_inq_varid(fid, 'm_i', imag_vid)
      call wrap_get_var_realx(fid, imag_vid, dust_imag)

      ! Close the file.
      call wrap_close(fid)
        
      ! Interpolate the values.
      if (do_print) write(LUNOPRT,*) 'wvlen(um)    m_dust'

      do i = 1, NWAVE
        do j = 1, dust_nwave
          if (wave(i) <= dust_wave(j)) then
            if ((j > 1) .and. (wave(i) /= dust_wave(j))) then
              interp = (wave(i) - dust_wave(j-1)) / (dust_wave(j) - dust_wave(j-1))
              carma_refidx_dust(i) = cmplx(dust_real(j-1) + interp*(dust_real(j) - dust_real(j-1)), dust_imag(j-1) + interp*(dust_imag(j) - dust_imag(j-1)))
            else
              carma_refidx_dust(i) = cmplx(dust_real(j), dust_imag(j))
            endif
          
            exit
          else if (j == dust_nwave) then
            carma_refidx_dust(i) = cmplx(dust_real(j), dust_imag(j))
          end if
        end do

        if (do_print) write(LUNOPRT,*) wave(i)*1e4_f, carma_refidx_dust(i)
      end do
    end if
    
#if ( defined SPMD )
    call mpibcast(carma_refidx_dust,  NWAVE, mpic16, 0, mpicom)
#endif

    call CARMAGROUP_Create(carma, I_GRP_DUST, "Dust", dust_rmin, dust_vmrat, I_SPHERE, 1._f, .false., &
                           rc, do_wetdep=.true., do_drydep=.true., solfac=0.3_f, &
                           scavcoef=0.1_f, shortname="CRDUST", refidx=carma_refidx_dust, &
                           imiertn=I_MIERTN_BOHREN1983, do_mie=.true.)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')

    ! Use the same refractive index at all wavelengths. This value is typical of soot and
    ! is recommended by Toon et al. 2012.
    refidx(:) = (1.8_f, 0.67_f)
    
    if (carma_fractal_soot) then
      RHO_SOOT = 1.8_f

      ! This matches the df profile used by Wolf and Toon [2010].
      soot_df(:) = (/ 3.0000_f, 3.0000_f, 1.5033_f, 1.5082_f, 1.5494_f, 1.6168_f, 1.7589_f, &
                      1.9957_f, 2.2519_f, 2.3840_f, 2.4000_f, 2.4000_f, 2.4000_f, 2.4000_f, &
                      2.4000_f, 2.4000_f, 2.4000_f, 2.4000_f, 2.4000_f, 2.4000_f, 2.4000_f /)

      call CARMAGROUP_Create(carma, I_GRP_SOOT, "Soot", soot_rmin, soot_vmrat, I_SPHERE, 1._f, .false., &
                             rc, do_wetdep=.true., do_drydep=.true., solfac=0.1_f, &
                             scavcoef=0.1_f, shortname="CRSOOT", refidx=refidx, do_mie=.true., &
                             is_fractal=.true., rmon=soot_rmon, df=soot_df, falpha=soot_falpha, &
                             imiertn=I_MIERTN_BOTET1997)
    else
      RHO_SOOT = 1.0_f
      call CARMAGROUP_Create(carma, I_GRP_SOOT, "Soot", soot_rmin, soot_vmrat, I_SPHERE, 1._f, .false., &
                             rc, do_wetdep=.true., do_drydep=.true., solfac=0.1_f, &
                             scavcoef=0.1_f, shortname="CRSOOT", refidx=refidx, do_mie=.true.)
    end if
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')
   
    call CARMAGROUP_Create(carma, I_GRP_SULFATE, "sulfate", so4_rmin, so4_vmrat, I_SPHERE, 1._f, .false., &
                           rc, irhswell=I_WTPCT_H2SO4, do_wetdep=.true., do_drydep=.true., solfac=1.0_f, &
                           scavcoef=0.1_f, is_sulfate=.true., shortname="PURSUL")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')

    call CARMAGROUP_Create(carma, I_GRP_CSOOT, "Coarse Soot", csoot_rmin, csoot_vmrat, I_SPHERE, 1._f, .false., &
                           rc, do_wetdep=.true., do_drydep=.true., solfac=0.1_f, &
                           scavcoef=0.1_f, shortname="CRCSOT", refidx=refidx, do_mie=.true.)

    
    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.
    call CARMAELEMENT_Create(carma, I_ELEM_DUST, I_GRP_DUST, "Dust", RHO_DUST, I_INVOLATILE, I_DUST, rc, shortname="CRDUST")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_SOOT, I_GRP_SOOT, "Soot", RHO_SOOT, I_INVOLATILE, I_SOOT, rc, shortname="CRSOOT")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_SULFATE, I_GRP_SULFATE, "Sulfate", RHO_SULFATE, &
         I_VOLATILE, I_H2SO4, rc, shortname="PURSUL")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_CSOOT, I_GRP_CSOOT, "Soot", RHO_CSOOT, I_INVOLATILE, I_SOOT, rc, shortname="CRCSOT")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')


    
    ! Define the Solutes

    
    ! Define the Gases
    !
    ! NOTE: Gases are defined just for emissions. They don't currently interact with the
    ! particles.
    call CARMAGAS_Create(carma, I_GAS_H2O, "Water Vapor", WTMOL_H2O, &
         I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, rc, shortname="Q")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')

    call CARMAGAS_Create(carma, I_GAS_H2SO4, "Sulfuric Acid", WTMOL_H2SO4, I_VAPRTN_H2SO4_AYERS1980, &
                         I_GCOMP_H2SO4, rc, shortname = "H2SO4", ds_threshold=-0.2_f)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')

    ! NOTE: We are only using this for emission, so it doesn't really matter that CARMA
    ! doesn't understand its properties.
    call CARMAGAS_Create(carma, I_GAS_SO2, "Sulfur Dioxide", WTMOL_SO2, &
         I_VAPRTN_NONE, I_GCOMP_SO2, rc, shortname="SO2")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')

    ! NOTE: We are only using this for emission, so it doesn't really matter that CARMA
    ! doesn't understand its properties.
    call CARMAGAS_Create(carma, I_GAS_CO2, "Carbon Dioxide", WTMOL_CO2, &
         I_VAPRTN_NONE, I_GCOMP_CO2, rc, shortname="CO2")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')

    ! NOTE: We are only using this for emission, so it doesn't really matter that CARMA
    ! doesn't understand its properties.
    call CARMAGAS_Create(carma, I_GAS_HBR, "Hydrobromic Acid", WTMOL_HBR, &
         I_VAPRTN_NONE, I_GCOMP_HBR, rc, shortname="HBR")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')

    ! NOTE: We are only using this for emission, so it doesn't really matter that CARMA
    ! doesn't understand its properties.
    call CARMAGAS_Create(carma, I_GAS_HCL, "Hydrochloric Acid", WTMOL_HCL, &
         I_VAPRTN_NONE, I_GCOMP_HCL, rc, shortname="HCL")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')

    ! NOTE: We are only using this for emission, so it doesn't really matter that CARMA
    ! doesn't understand its properties.
    call CARMAGAS_Create(carma, I_GAS_NO, "Nitric Oxide", WTMOL_NO, &
         I_VAPRTN_NONE, I_GCOMP_NO, rc, shortname="NO")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')

    
    ! Define the Processes
    call CARMA_AddCoagulation(carma, I_GRP_DUST, I_GRP_DUST, I_GRP_DUST, I_COLLEC_DATA, rc)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_SOOT, I_GRP_SOOT, I_GRP_SOOT, I_COLLEC_DATA, rc)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_CSOOT, I_GRP_CSOOT, I_GRP_CSOOT, I_COLLEC_DATA, rc)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    ! For sulfates, set H2SO4 to be the condensing gas, water vapor is assumed to be in equilibrium
    ! and will be used to define the wet particle radius.
    call CARMA_AddGrowth(carma, I_ELEM_SULFATE, I_GAS_H2SO4, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

    call CARMA_AddNucleation(carma, I_ELEM_SULFATE, I_ELEM_SULFATE, I_HOMNUC, 0._f, rc, igas=I_GAS_H2SO4)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_SULFATE, I_GRP_SULFATE, I_GRP_SULFATE, I_COLLEC_FUCHS, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    call pbuf_add_field('SADSULF', 'global', dtype_r8, (/pcols, pver/), ipbuf4sad)
    
    if (carma_rad_feedback) then
       call pbuf_add_field('VOLC_RAD_GEOM', 'global', dtype_r8, (/pcols, pver/), ipbuf4reff)
       call pbuf_add_field('VOLC_MMR', 'global', dtype_r8, (/pcols, pver/), ipbuf4so4mmr)
    endif

    return
  end subroutine CARMA_DefineModel


  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen 
  !!
  !!  @see CARMASTATE_SetDetrain
  subroutine CARMA_Detrain(carma, cstate, cam_in, dlf, state, icol, dt, rc, rliq, prec_str, snow_str, &
      tnd_qsnow, tnd_nsnow)
    use camsrfexch,         only: cam_in_t
    use physconst,          only: latice, latvap, cpair

    implicit none
    
    type(carma_type), intent(in)         :: carma            !! the carma object
    type(carmastate_type), intent(inout) :: cstate           !! the carma state object
    type(cam_in_t),  intent(in)          :: cam_in           !! surface input
    real(r8), intent(in)                 :: dlf(pcols, pver) !! Detraining cld H20 from convection (kg/kg/s)
    type(physics_state), intent(in)      :: state            !! physics state variables
    integer, intent(in)                  :: icol             !! column index
    real(r8), intent(in)                 :: dt               !! time step (s)
    integer, intent(out)                 :: rc               !! return code, negative indicates failure
    real(r8), intent(inout), optional    :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional    :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional    :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    real(r8), intent(out), optional      :: tnd_qsnow(pcols,pver) !! snow mass tendency (kg/kg/s)
    real(r8), intent(out), optional      :: tnd_nsnow(pcols,pver) !! snow number tendency (#/kg/s)
    
    ! Default return code.
    rc = RC_OK
        
    return
  end subroutine CARMA_Detrain


  !! For diagnostic groups, sets up up the CARMA bins based upon the CAM state.
  !!
  !!  @version July-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DiagnoseBins(carma, cstate, state, pbuf, icol, dt, rc, rliq, prec_str, snow_str)
    use time_manager,     only: is_first_step

    implicit none
    
    type(carma_type), intent(in)          :: carma        !! the carma object
    type(carmastate_type), intent(inout)  :: cstate       !! the carma state object
    type(physics_state), intent(in)       :: state        !! physics state variables
    type(physics_buffer_desc), pointer    :: pbuf(:)      !! physics buffer
    integer, intent(in)                   :: icol         !! column index
    real(r8), intent(in)                  :: dt           !! time step
    integer, intent(out)                  :: rc           !! return code, negative indicates failure
    real(r8), intent(in), optional        :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional     :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional     :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    
    real(r8)                             :: mmr(pver) !! elements mass mixing ratio
    integer                              :: ibin      !! bin index
    
    ! Default return code.
    rc = RC_OK
    
    ! By default, do nothing. If diagnosed groups exist, this needs to be replaced by
    ! code to determine the mass in each bin from the CAM state.
    
    return
  end subroutine CARMA_DiagnoseBins
  
  
  !! For diagnostic groups, determines the tendencies on the CAM state from the CARMA bins.
  !!
  !!  @version July-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DiagnoseBulk(carma, cstate, cam_out, state, pbuf, ptend, icol, dt, rc, rliq, prec_str, snow_str, &
    prec_sed, snow_sed, tnd_qsnow, tnd_nsnow, re_ice)
    use camsrfexch,       only: cam_out_t
    use physics_buffer,   only: pbuf_get_field

    implicit none
    
    type(carma_type), intent(in)         :: carma     !! the carma object
    type(carmastate_type), intent(inout) :: cstate    !! the carma state object
    type(cam_out_t),      intent(inout)  :: cam_out   !! cam output to surface models
    type(physics_state), intent(in)      :: state     !! physics state variables
    type(physics_buffer_desc), pointer   :: pbuf(:)   !! physics buffer
    type(physics_ptend), intent(inout)   :: ptend     !! constituent tendencies
    integer, intent(in)                  :: icol      !! column index
    real(r8), intent(in)                 :: dt        !! time step
    integer, intent(out)                 :: rc        !! return code, negative indicates failure
    real(r8), intent(inout), optional    :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional    :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional    :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    real(r8), intent(inout), optional    :: prec_sed(pcols)       !! total precip from cloud sedimentation (m/s)
    real(r8), intent(inout), optional    :: snow_sed(pcols)       !! snow from cloud ice sedimentation (m/s)
    real(r8), intent(inout), optional    :: tnd_qsnow(pcols,pver) !! snow mass tendency (kg/kg/s)
    real(r8), intent(inout), optional    :: tnd_nsnow(pcols,pver) !! snow number tendency (#/kg/s)
    real(r8), intent(out), optional      :: re_ice(pcols,pver)    !! ice effective radius (m)
    
    integer                              :: ielem     ! element index
    integer                              :: ibin      ! bin index
    real(r8)                             :: mmr(pver) ! mass mixing ration (kg/kg)
    real(r8)                             :: sflx      ! surface flux (kg/m2/s)
    real(r8)                             :: numberDensity(pver)
    real(r8)                             :: ad(pver)       ! stratospheric aerosol wet surface area density (cm2/cm3)
    real(r8)                             :: reff(pver)     ! stratospheric wet effective radius (m)
    real(r8)                             :: md(pver)       ! bin integrated stratospheric mass mixing ratio (kg/kg)
    real(r8)                             :: r_wet(pver)    ! Sulfate aerosol bin wet radius (cm)
    real(r8), pointer, dimension(:,:)    :: sadsulf_ptr           ! Sulfate surface area density pointer
    real(r8), pointer, dimension(:,:)    :: reffsulf_ptr          ! Sulfate effective radius pointer
    real(r8), pointer, dimension(:,:)    :: mmrsulf_ptr           ! Sulfate mass mixing ratio pointer
    integer                              :: igroup

    ! Default return code.
    rc = RC_OK

    ! Add the sedimentation and dry deposition fluxes to the hydrophilic black carbon.
    !
    ! NOTE: Don't give the surface model negative values for the surface fluxes.
    ielem = I_ELEM_SOOT
    do ibin = 1, NBIN
    
      call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, sedimentationFlux=sflx)
      if (rc < 0) call endrun('CARMA_DiagnoseBulk::CARMA_GetBin failed.')
      
      cam_out%bcphidry(icol) = cam_out%bcphidry(icol) + max(sflx, 0._r8)
    end do

    ielem = I_ELEM_DUST
    do ibin = 1, NBIN
    
      call CARMASTATE_GetBin(cstate, ielem, ibin, mmr, rc, sedimentationFlux=sflx)
      if (rc < 0) call endrun('CARMA_DiagnoseBulk::CARMA_GetBin failed.')
      
      if (carma_dustmap(ibin) == 1) then
        cam_out%dstdry1(icol) = cam_out%dstdry1(icol) + max(sflx, 0._r8)
      else if (carma_dustmap(ibin) == 2) then
        cam_out%dstdry2(icol) = cam_out%dstdry2(icol) + max(sflx, 0._r8)
      else if (carma_dustmap(ibin) == 3) then
        cam_out%dstdry3(icol) = cam_out%dstdry3(icol) + max(sflx, 0._r8)
      else if (carma_dustmap(ibin) == 4) then
        cam_out%dstdry4(icol) = cam_out%dstdry4(icol) + max(sflx, 0._r8)
      end if
    end do
    
    ! Calcualte the radiative properties for the sulfates.
    call CARMAELEMENT_Get(carma, I_ELEM_SULFATE, rc, igroup=igroup)
    if (rc < 0) call endrun('CARMA_DiagnoseBulk::CARMAELEMENT_Get failed.')

    ad(:)  = 0.0_r8     ! stratospheric wet aerosol surface area density (cm2/cm3)
    md(:)  = 0.0_r8     ! bin integrated stratospheric mass mixing ratio (kg/kg)
    reff(:)  = 0.0_r8   ! stratospheric effective radius (m)

    do ibin = 1, NBIN
      call CARMASTATE_GetBin(cstate, I_ELEM_SULFATE, ibin, mmr(:), rc, &
                             numberDensity=numberDensity, r_wet=r_wet)
      if (rc < 0) call endrun('CARMA_DiagnoseBulk::CARMASTATE_GetBin failed.')

      ! Calculate the total densities.
      !
      ! NOTE: Calculate AD in cm2/cm3.
      if (numberDensity(1) /= CAM_FILL) then
        ad(:)  = ad(:)  + numberDensity(:) * (r_wet(:)**2)
        reff(:) = reff(:) + numberDensity(:) * (r_wet(:)**3)
        md(:)  = md(:)  + mmr(:)  ! bin integrated stratospheric mass mixing ratio (kg/kg)
      end if
    end do
    
    reff(:) = reff(:) / ad(:) ! wet effective radius in cm
    reff(:) = reff(:) / 100.0_r8 ! cm -> m
    ad(:)  = ad(:) * 4.0_r8 * PI ! surface area density in cm2/cm3
    
    ! Add the computed fields to the physics buffer.
    call pbuf_get_field(pbuf, ipbuf4sad, sadsulf_ptr)
    sadsulf_ptr(icol, :cstate%f_NZ) = ad(:cstate%f_NZ)    ! stratospheric aerosol wet surface area density (cm2/cm3)

    if (carma_rad_feedback) then
      call pbuf_get_field(pbuf, ipbuf4reff, reffsulf_ptr)
      reffsulf_ptr(icol, :cstate%f_NZ) = reff(:cstate%f_NZ) ! stratospheric wet effective radius (m)

      call pbuf_get_field(pbuf, ipbuf4so4mmr, mmrsulf_ptr)
      mmrsulf_ptr(icol, :cstate%f_NZ) = md(:cstate%f_NZ)    ! bin integrated stratospheric mass mixing ratio (kg/kg)
    end if    
    
    return
  end subroutine CARMA_DiagnoseBulk


  !! Calculates the emissions for gases associated with the CARMA model. By default,
  !! there is no emission, but this routine can be overridden for models that wish
  !! to have a gas emission related to the aerosol emission.
  !!
  !! @author  Chuck Bardeen
  !! @version August-2016
  subroutine CARMA_EmitGas(carma, igas, icnst, dt, state, cam_in, tendency, surfaceFlux, rc)
    use shr_kind_mod,  only: r8 => shr_kind_r8
    use ppgrid,        only: pcols, pver
    use physics_types, only: physics_state
    use phys_grid,     only: get_rlon_all_p, get_rlat_all_p
    use time_manager,  only: get_curr_date, get_perp_date, get_curr_calday, &
                             is_perpetual, is_first_step
    use camsrfexch,    only: cam_in_t
    use tropopause,    only: tropopause_find
    use physconst,     only: gravit
    
    implicit none
    
    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: igas                  !! gas index
    integer, intent(in)                :: icnst                 !! consituent index
    real(r8), intent(in)               :: dt                    !! time step (s)
    type(physics_state), intent(in)    :: state                 !! physics state
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    real(r8), intent(out)              :: tendency(pcols, pver) !! constituent tendency (kg/kg/s)
    real(r8), intent(out)              :: surfaceFlux(pcols)    !! constituent surface flux (kg/m^2/s)
    integer, intent(out)               :: rc                    !! return code, negative indicates failure
    
    real(r8), parameter                :: mu_dust = 6.6_r8  ! width parameter, dust, tropopause (km)
    real(r8), parameter                :: mu_soot_gnd  = 1._r8  ! width parameter, soot, ground (km)
    real(r8), parameter                :: mu_soot_trop = 3._r8  ! width parameter, soot, tropopause (km)

    integer       :: tropLev(pcols)           ! tropopause level index   
    real(r8)      :: tropP(pcols)             ! tropopause pressure (Pa)  
    real(r8)      :: tropT(pcols)             ! tropopause temperature (K) 
    real(r8)      :: tropZ(pcols)             ! tropopause height (m) 

    real(r8)     :: lon(pcols)              ! longitude
    real(r8)     :: lat(pcols)              ! latitude
    integer      :: igroup                  ! group index
    integer      :: lchnk                   ! chunk identifier
    integer      :: ncol                    ! number of columns in chunk
    integer      :: icol                    ! column index
    integer      :: k                       ! vertical index
    real(r8)     :: calday                  ! current calendar day
    integer      :: currentDate             ! current date (yyyydoy)
    integer      :: yr                      ! year
    integer      :: mon                     ! month
    integer      :: day                     ! day of month
    integer      :: ncsec                   ! time of day (seconds)
    integer      :: doy                     ! day of year
    real(r8)     :: startyear               ! start year
    real(r8)     :: stopyear                ! stop year
    real(r8)     :: startdoy                ! start year
    real(r8)     :: stopdoy                 ! stop year
    integer      :: emis_time               ! length of time for emission
    real(r8)     :: vfunc(pver)             ! scaling factor to preserve total emission
    character(len=32) :: shortname          ! the shortname of the group
    real(r8)     :: zmid                    ! layer midpoint altitude (km)
    real(r8)     :: ztrop                   ! tropopause altitude (km)
    real(r8)     :: rate                    ! emission rate (kg/s/m)
    real(r8)     :: massflux                ! mass flux (kg/m3/s)
    real(r8)     :: thickness               ! layer thickness (m)
    real(r8)     :: dist                    ! great circle distance
    real(r8)     :: tendency_splash(pcols, pver) !! constituent tendency - splash (kg/kg/s)
    real(r8)     :: tendency_fire(pcols, pver)   !! constituent tendency - fire (kg/kg/s)
    real(r8)     :: tendency_impact(pcols, pver) !! constituent tendency - ballistic (kg/kg/s)
    integer      :: elapsed_dtime           ! time since start of event (sec)
    real(r8)     :: dlat
    real(r8)     :: dlon
    integer      :: ilat
    integer      :: ilon

    ! Default return code.
    rc = RC_OK

    lchnk = state%lchnk
    ncol = state%ncol
    
    ! Add any surface flux here.
    ncol = state%ncol
    surfaceFlux(:ncol) = 0.0_f
    
    ! For emissions into the atmosphere, put the emission here.
    !
    ! Use Toon et al. [2016] as the source function for gases from a
    ! meteor impact.
    ! 
    ! For water vapor, there are three locations for emission:
    !
    !   - from the impact (with the dust)
    !   - from the fires  (with the soot)
    !   - splashed sea water (uniform mixing ratio above the tropopause).
    !
    ! See EmitParticle for how the soot and dust regions are defined.
    tendency(:ncol, :pver) = 0.0_r8
    
    ! We only do enissions for certain gases, and the gases must be defined as
    ! CARMA gases.
    call CARMAGAS_GET(carma, igas, rc, shortname=shortname)
    if (RC < RC_ERROR) return
    
    if ((shortname /= "CO2") .and. (shortname /= "Q")   .and. (shortname /= "SO2") .and. &
        (shortname /= "HBR") .and. (shortname /= "HCL") .and. (shortname /= "NO")) then
      return
    end if    

    ! Determine the latitude and longitude of each column.
    lchnk = state%lchnk
    ncol = state%ncol
    
    lat(:ncol)  = state%lat(:ncol)
    lon(:ncol)  = state%lon(:ncol)    

    ! Assume a regular grid.
    dlat = PI / (plat-1)
    dlon = 2._f*PI / plon
    
    ! Determine the day of year.
    calday = get_curr_calday()
    if ( is_perpetual() ) then
      call get_perp_date(yr, mon, day, ncsec)
    else
      call get_curr_date(yr, mon, day, ncsec)
    end if
    doy = floor(calday)
    
    ! Determine the start and stop year and day of year from the namelist
    ! variables.
    currentDate = yr * 1000 + doy
    startyear = carma_emis_startdate / 1000
    stopyear  = carma_emis_stopdate  / 1000
    
    startdoy  = mod(carma_emis_startdate, 1000)
    stopdoy   = mod(carma_emis_stopdate, 1000)

    ! Make sure to emit for at least one timestep and in multiples of the time
    ! step length.
    ! TBD - This has a leap year problem, but works otherwise ...
    carma_emis_dtime = INT((((stopyear - startyear) * 365._f + (stopdoy - startdoy)) * 24._f * 3600._f + &
                 (carma_emis_stoptime - carma_emis_starttime)) / dt) * dt
    elapsed_dtime = INT((((yr - startyear) * 365._f + (doy - startdoy)) * 24._f * 3600._f + &
                 (ncsec - carma_emis_starttime)) / dt) * dt

    ! Find the tropopause using the default algorithm backed by the climatology.
    call tropopause_find(state, tropLev, tropZ=tropZ)

    
    ! Do the water vapor from the impact site and from splashed sea water. This
    ! is uniformly mixed above the tropopause and is all emitted at the start of
    ! the event.
    tendency_splash(:ncol, :pver) = 0._r8
!    if ((currentDate == carma_emis_startdate) .and. &
!        ((ncsec >= carma_emis_starttime) .and. (ncsec < (carma_emis_starttime + dt)))) then
    if (((currentDate > carma_emis_startdate) .or. &
         ((currentDate == carma_emis_startdate) .and. (ncsec >= carma_emis_starttime))) .and. &
        ((currentDate < carma_emis_stopdate) .or. &
         ((currentDate == carma_emis_stopdate) .and. (ncsec < carma_emis_stoptime))) .and. &
        (elapsed_dtime < carma_emis_dust_dtime)) then

      ! Make sure to emit for at least one timestep and in multiples of the time
      ! step length.
      ! TBD - This has a leap year problem, but works otherwise ...
      carma_emis_dtime = INT((((stopyear - startyear) * 365._f + (stopdoy - startdoy)) * 24._f * 3600._f + &
                   (carma_emis_stoptime - carma_emis_starttime)) / dt) * dt

      ! Loop over all of the columns.
      do icol = 1, ncol

        ! Calculate the great circle distance (in km).
        dist = REARTH  / 1e5_f * &
               abs(acos(sin(lat(icol))*sin(carma_emis_ctrlat*DEG2RAD) + &
               cos(lat(icol))*cos(carma_emis_ctrlat*DEG2RAD)*cos(lon(icol)-(carma_emis_ctrlon*DEG2RAD))))
                
        ! Is the column one of the ones over which there should be emissions>
        if (dist <= carma_emis_splash_radius) then

          ! Set tendencies for any sources or sinks in the atmosphere.
          do k = 1, pver
        
            ! Get the cell midpoint and height
            zmid  = state%zm(icol, k) / 1000._f

            ! Get the tropopause height.
            ztrop = tropZ(icol) / 1000._f

            ! These are uniformly mixed above the tropopause, with a constant
            ! mass mixing ratio.
            if (zmid > ztrop) then
              vfunc(k) = state%pdel(icol, k) / gravit
            else
              vfunc(k) = 0._r8
            end if
            
            if (shortname == "CO2") then
              rate = carma_emis_co2_splash
            else if (shortname == "Q") then
              rate = carma_emis_h2o_splash
            else if (shortname == "SO2") then
              rate = carma_emis_so2_splash
            else if (shortname == "HBR") then
              rate = carma_emis_hbr_splash
            else if (shortname == "HCL") then
              rate = carma_emis_hcl_splash
            else if (shortname == "NO") then
              rate = carma_emis_no_splash
            end if
          
            ! Calculate a rate by dividing by total emission time.
            rate = rate  * vfunc(k) / carma_emis_dust_dtime
            
            ! Scale for the fraction of the total surface area that is emitting and
            ! convert to kg/m2/s
            massflux = rate / carma_emis_splash_area
        
            ! Convert the mass flux to a tendency on the mass mixing ratio.
            tendency_splash(icol, k) = massflux / (state%pdel(icol, k) / gravit)
          end do
      
          ! Now normalize in the vertical to preserve the total mass.
          tendency_splash(icol, :) = tendency_splash(icol, :) / sum(vfunc(:))
          end if
      end do
    end if
    
      
    ! Do the ballistically distributed water vapor. It is assumed to all be injected
    ! in the first time step after the impact.
    tendency_impact(:ncol, :pver) = 0._r8
    if (((currentDate > carma_emis_startdate) .or. &
         ((currentDate == carma_emis_startdate) .and. (ncsec >= carma_emis_starttime))) .and. &
        ((currentDate < carma_emis_stopdate) .or. &
         ((currentDate == carma_emis_stopdate) .and. (ncsec < carma_emis_stoptime))) .and. &
        (elapsed_dtime < carma_emis_dust_dtime)) then
!    if ((currentDate == carma_emis_startdate) .and. &
!        ((ncsec >= carma_emis_starttime) .and. (ncsec < (carma_emis_starttime + dt)))) then
  
      ! Make sure to emit for at least one timestep and in multiples of the time
      ! step length.
      ! TBD - This has a leap year problem, but works otherwise ...
      carma_emis_dtime = INT((((stopyear - startyear) * 365._f + (stopdoy - startdoy)) * 24._f * 3600._f + &
                   (carma_emis_stoptime - carma_emis_starttime)) / dt) * dt

      ! Loop over all of the columns.
      do icol = 1, ncol

        ! Calculate the great circle distance (in km).
        dist = REARTH  / 1e5_f * &
               abs(acos(sin(lat(icol))*sin(carma_emis_ctrlat*DEG2RAD) + &
               cos(lat(icol))*cos(carma_emis_ctrlat*DEG2RAD)*cos(lon(icol)-(carma_emis_ctrlon*DEG2RAD))))
                  
        ! Is the column one of the ones over which there should be emissions>
        if (dist <= carma_emis_dust_radius) then

          ! Set tendencies for any sources or sinks in the atmosphere.
          do k = 1, pver
          
            ! Get the cell midpoint and height
            zmid  = state%zm(icol, k) / 1000._f
  
            ! Determine the total emission rate for this grid box using equation 2
            ! from Toon et al. [2012] and also adjust for the fraction of the
            ! mass that goes into the specified bin based on the assumed size
            ! distribution also from Toon et al. [2012]. This is a Gaussian
            ! centered at 70 km.
            vfunc(k) = 1._f / (mu_dust * sqrt(2._f * PI)) * &
                       exp(-0.5_f * (((zmid - carma_emis_z_ballistic) / mu_dust)**2)) * &
                       (state%zi(icol, k) - state%zi(icol, k+1))
                     
            if (shortname == "CO2") then
              rate = carma_emis_co2_impact
            else if (shortname == "Q") then
              rate = carma_emis_h2o_impact
            else if (shortname == "SO2") then
              rate = carma_emis_so2_impact
            else if (shortname == "HBR") then
              rate = carma_emis_hbr_impact
            else if (shortname == "HCL") then
              rate = carma_emis_hcl_impact
            else if (shortname == "NO") then
              rate = carma_emis_no_impact
            end if
            
            ! Calculate a rate by dividing by total emission time.
            rate = rate  * vfunc(k) / carma_emis_dust_dtime
            
            ! Scale for the fraction of the total surface area that is emitting and
            ! convert to kg/m2/s
            massflux = rate / carma_emis_dust_area
          
            ! Convert the mass flux to a tendency on the mass mixing ratio.
            tendency_impact(icol, k) = massflux / (state%pdel(icol, k) / gravit)
          end do
        
          ! Now normalize in the vertical to preserve the total mass.
          tendency_impact(icol, :) = tendency_impact(icol, :) / sum(vfunc(:))
        end if
      end do
    end if
    
    
    ! Do the fire generated water vapor. It is assumed to be injected during
    ! the entire time period.
    tendency_fire(:ncol, :pver) = 0._r8
    if (((currentDate > carma_emis_startdate) .or. &
         ((currentDate == carma_emis_startdate) .and. (ncsec >= carma_emis_starttime))) .and. &
        ((currentDate < carma_emis_stopdate) .or. &
         ((currentDate == carma_emis_stopdate) .and. (ncsec < carma_emis_stoptime)))) then

      ! Make sure to emit for at least one timestep and in multiples of the time
      ! step length.
      ! TBD - This has a leap year problem, but works otherwise ...
      carma_emis_dtime = INT((((stopyear - startyear) * 365._f + (stopdoy - startdoy)) * 24._f * 3600._f + &
                   (carma_emis_stoptime - carma_emis_starttime)) / dt) * dt

      ! Loop over all of the columns.
      do icol = 1, ncol

        ! Calculate the great circle distance (in km)
        dist = REARTH / 1e5_f * &
               abs(acos(sin(lat(icol))*sin(carma_emis_ctrlat*DEG2RAD) + &
               cos(lat(icol))*cos(carma_emis_ctrlat*DEG2RAD)*cos(lon(icol)-(carma_emis_ctrlon*DEG2RAD))))
                
        ! Is the column one of the ones over which there should be emissions>
        if (dist <= carma_emis_soot_radius) then

          ! NOTE: Assume a regular grids for now.
          ilat = nint((lat(icol) + (PI / 2._r8)) / dlat) + 1
          ilon = nint(lon(icol) / dlon) + 1
    
          ! Set tendencies for any sources or sinks in the atmosphere.
          do k = 1, pver
          
            ! Get the cell midpoint and height
            zmid  = state%zm(icol, k) / 1000._f
  
            ! Get the tropopause height.
            ztrop = tropZ(icol) / 1000._f

            ! Determine the total emission rate for this grid box using equation 2
            ! from Toon et al. [2012] and also adjust for the fraction of the
            ! mass that goes into the specified bin based on the assumed size
            ! distribution also from Toon et al. [2012].
            vfunc(k) = &
                  (1._f / sqrt(2._f * PI)) * ( &
                    carma_emis_fine_trop_frac / mu_soot_trop * exp(-0.5_f * (((zmid - ztrop) / mu_soot_trop)**2)) + &
                    (1._f - carma_emis_fine_trop_frac) / mu_soot_gnd  * exp(-0.5_f * ((zmid / mu_soot_gnd)**2))) * &
                  (state%zi(icol, k) - state%zi(icol, k+1))
            
            if (shortname == "CO2") then
              rate = carma_emis_co2_fire
            else if (shortname == "Q") then
              rate = carma_emis_h2o_fire
            else if (shortname == "SO2") then
              rate = carma_emis_so2_fire
            else if (shortname == "HBR") then
              rate = carma_emis_hbr_fire
            else if (shortname == "HCL") then
              rate = carma_emis_hcl_fire
            else if (shortname == "NO") then
              rate = carma_emis_no_fire
            end if

            ! Apply a linearly decreasing trend on the rate.
            rate = rate * 2._r8 * max(0._r8, (1._r8 - elapsed_dtime / (carma_emis_dtime-dt)))
              
            ! Calculate a rate by dividing by total emission time.
            rate = rate  * vfunc(k) / carma_emis_dtime
              
            ! Scale for the fraction of the total surface area that is emitting and
            ! convert to kg/m2/s
            massflux = rate / carma_emis_soot_area
            
            ! If only emitting proportional to vegetation, then scale by veg fraction.
            if (carma_emis_fire_veg_scale) then
              massflux = massflux * carma_emis_soot_area / carma_emis_tot_veg * carma_vegfrac_map(ilon, ilat)
            ! If only emitting over land, then scale by land fraction.
            else if (carma_emis_fire_land_only) then
              massflux = massflux * cam_in.landfrac(icol)
            end if

            ! Convert the mass flux to a tendency on the mass mixing ratio.
            tendency_fire(icol, k) = massflux / (state%pdel(icol, k) / gravit)
          end do
          
          ! Now normalize in the vertical to preserve the total mass.
          tendency_fire(icol, :) = tendency_fire(icol, :) / sum(vfunc(:))
        end if
      end do
    end if
    
    ! Combine the tendencies from different sources.
    tendency(:ncol, :) = tendency_splash(:ncol, :) + tendency_fire(:ncol, :) + tendency_impact(:ncol, :)
    
    return
  end subroutine CARMA_EmitGas


  !! Calculates the heat emissions associated with the CARMA model. By default,
  !! there is no emission, but this routine can be overridden for models that wish
  !! to have a gas emission related to the aerosol emission.
  !!
  !! @author  Chuck Bardeen
  !! @version June-2017
  subroutine CARMA_EmitHeat(carma, dt, state, cam_in, tendency, rc)
    use shr_kind_mod,  only: r8 => shr_kind_r8
    use ppgrid,        only: pcols, pver
    use physics_types, only: physics_state
    use phys_grid,     only: get_rlon_all_p, get_rlat_all_p
    use time_manager,  only: get_curr_date, get_perp_date, get_curr_calday, &
                             is_perpetual, is_first_step
    use camsrfexch,    only: cam_in_t
    use physconst,     only: gravit
    
    implicit none
    
    type(carma_type), intent(in)       :: carma                 !! the carma object
    real(r8), intent(in)               :: dt                    !! time step (s)
    type(physics_state), intent(in)    :: state                 !! physics state
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    real(r8), intent(out)              :: tendency(pcols, pver) !! constituent tendency (kg/kg/s)
    integer, intent(out)               :: rc                    !! return code, negative indicates failure
    
    integer      :: heat_fire_ktop          ! top layer for fire heating
    
    real(r8)     :: lon(pcols)              ! longitude
    real(r8)     :: lat(pcols)              ! latitude
    integer      :: igroup                  ! group index
    integer      :: lchnk                   ! chunk identifier
    integer      :: ncol                    ! number of columns in chunk
    integer      :: icol                    ! column index
    integer      :: k                       ! vertical index
    real(r8)     :: calday                  ! current calendar day
    integer      :: currentDate             ! current date (yyyydoy)
    integer      :: yr                      ! year
    integer      :: mon                     ! month
    integer      :: day                     ! day of month
    integer      :: ncsec                   ! time of day (seconds)
    integer      :: doy                     ! day of year
    real(r8)     :: startyear               ! start year
    real(r8)     :: stopyear                ! stop year
    real(r8)     :: startdoy                ! start year
    real(r8)     :: stopdoy                 ! stop year
    integer      :: emis_time               ! length of time for emission
    real(r8)     :: vfunc(pver)             ! scaling factor to preserve total emission
    real(r8)     :: mass(pver)              ! air mass
    real(r8)     :: tmass                   ! total column air mass for emission
    real(r8)     :: rate                    ! emission rate (J/s)
    real(r8)     :: heatflux                ! mass flux (kg/m3/s)
    real(r8)     :: dist                    ! great circle distance
    real(r8)     :: tendency_fire(pcols, pver)   !! heat tendency - fire (J/s)
    integer      :: elapsed_dtime           ! time since start of event (sec)
    real(r8)     :: dlat
    real(r8)     :: dlon
    integer      :: ilat
    integer      :: ilon

    ! Default return code.
    rc = RC_OK

    lchnk = state%lchnk
    ncol = state%ncol
    
    ! For emissions into the atmosphere, put the emission here.
    !
    ! Use Toon et al. [2016] as the source function for heat from a
    ! meteor impact.
    ! 
    ! For heat, there is currently one location for emission:
    !
    !   - from the fires  (at the surface)
 
    ! See EmitParticle for how the soot and dust regions are defined.
    tendency(:ncol, :pver) = 0.0_f
    
    
    ! Determine the day of year.
    calday = get_curr_calday()
    if ( is_perpetual() ) then
      call get_perp_date(yr, mon, day, ncsec)
    else
      call get_curr_date(yr, mon, day, ncsec)
    end if
    doy = floor(calday)

    ! Determine the latitude and longitude of each column.
    lat(:ncol)  = state%lat(:ncol)
    lon(:ncol)  = state%lon(:ncol)    

    ! Assume a regular grid.
    dlat = PI / (plat-1)
    dlon = 2*PI / plon
    
    ! Determine the start and stop year and day of year from the namelist
    ! variables.
    currentDate = yr * 1000 + doy
    startyear = carma_emis_startdate / 1000
    stopyear  = carma_emis_stopdate  / 1000
    
    startdoy  = mod(carma_emis_startdate, 1000)
    stopdoy   = mod(carma_emis_stopdate, 1000)

    ! Make sure to emit for at least one timestep and in multiples of the time
    ! step length.
    ! TBD - This has a leap year problem, but works otherwise ...
    carma_emis_dtime = INT((((stopyear - startyear) * 365._f + (stopdoy - startdoy)) * 24._f * 3600._f + &
                 (carma_emis_stoptime - carma_emis_starttime)) / dt) * dt
    elapsed_dtime = INT((((yr - startyear) * 365._f + (doy - startdoy)) * 24._f * 3600._f + &
                 (ncsec - carma_emis_starttime)) / dt) * dt
    
    
    ! Do the fire generated heat. It is assumed to be injected during
    ! the entire time period.
    tendency_fire(:ncol, :pver) = 0._r8
    heat_fire_ktop = pver - carma_emis_heat_fire_levs + 1

    if (((currentDate > carma_emis_startdate) .or. &
         ((currentDate == carma_emis_startdate) .and. (ncsec >= carma_emis_starttime))) .and. &
        ((currentDate < carma_emis_stopdate) .or. &
         ((currentDate == carma_emis_stopdate) .and. (ncsec < carma_emis_stoptime)))) then

      ! Loop over all of the columns.
      do icol = 1, ncol

        ! Calculate the great circle distance (in km)
        dist = REARTH / 1e5_f * &
               abs(acos(sin(lat(icol))*sin(carma_emis_ctrlat*DEG2RAD) + &
               cos(lat(icol))*cos(carma_emis_ctrlat*DEG2RAD)*cos(lon(icol)-(carma_emis_ctrlon*DEG2RAD))))
                
        ! Is the column one of the ones over which there should be emissions>
        if (dist <= carma_emis_soot_radius) then
        
          ! NOTE: Assume a regular grids for now.
          ilat = nint((lat(icol) + (PI / 2._r8)) / dlat) + 1
          ilon = nint(lon(icol) / dlon) + 1
    
          ! Deposit energy in the same fraction as the mass fraction in the column.
          mass(:) = state%pdel(icol, :) / gravit

          ! Emit in the surface layers.
          tmass = sum(mass(heat_fire_ktop:pver))
          
          ! Set tendencies for any sources or sinks in the atmosphere.
          do k = heat_fire_ktop, pver
          
            ! Heat emission is proportional to the mass fraction of this
            ! gridbox compared to the total column over which heat emissions
            ! are happening.
            vfunc(k) = mass(k) / tmass

            rate = carma_emis_heat_fire
            
            ! Apply a linearly decreasing trend on the rate.
            rate = rate * 2._r8 * max(0._r8, (1._r8 - elapsed_dtime / carma_emis_dtime))

            ! Calculate a rate by dividing by total emission time.
            rate = rate  * vfunc(k) / carma_emis_dtime
              
            ! Scale for the fraction of the total surface area that is emitting and
            ! convert to J/m2/s
            heatflux = rate / carma_emis_soot_area
            
            ! If only emitting proportional to vegetation, then scale by veg fraction.
            if (carma_emis_fire_veg_scale) then
              heatflux = heatflux * carma_emis_soot_area / carma_emis_tot_veg * carma_vegfrac_map(ilon, ilat)
            ! If only emitting over land, then scale by land fraction.
            else if (carma_emis_fire_land_only) then
              heatflux = heatflux * cam_in.landfrac(icol)
            end if

            ! Convert the mass flux to a tendency on the dry static energy (J/kg/s)
            tendency_fire(icol, k) = heatflux / (state%pdel(icol, k) / gravit)
          end do
        end if
      end do
    end if
    
    ! Combine the tendencies from different sources.
    tendency(:ncol, :) = tendency_fire(:ncol, :)
    
    return
  end subroutine CARMA_EmitHeat  


  !! Calculates the emissions for CARMA aerosol particles. By default, there is no
  !! emission, but this routine can be overridden for models that wish to have
  !! an aerosol emission.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_EmitParticle(carma, ielem, ibin, icnst, dt, state, cam_in, tendency, surfaceFlux, rc)
    use shr_kind_mod,  only: r8 => shr_kind_r8
    use ppgrid,        only: pcols, pver
    use physics_types, only: physics_state
    use phys_grid,     only: get_rlon_all_p, get_rlat_all_p
    use time_manager,  only: get_curr_date, get_perp_date, get_curr_calday, &
                             is_perpetual, is_first_step
    use camsrfexch,    only: cam_in_t
    use tropopause,    only: tropopause_find
    use physconst,     only: gravit
    
    implicit none
    
    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: ielem                 !! element index
    integer, intent(in)                :: ibin                  !! bin index
    integer, intent(in)                :: icnst                 !! consituent index
    real(r8), intent(in)               :: dt                    !! time step (s)
    type(physics_state), intent(in)    :: state                 !! physics state
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    real(r8), intent(out)              :: tendency(pcols, pver) !! constituent tendency (kg/kg/s)
    real(r8), intent(out)              :: surfaceFlux(pcols)    !! constituent surface flux (kg/m^2/s)
    integer, intent(out)               :: rc                    !! return code, negative indicates failure
    
    real(r8), parameter                :: mu_dust = 6.6_r8  ! width parameter, dust, tropopause (km)
    real(r8), parameter                :: mu_soot_gnd  = 1._r8  ! width parameter, soot, ground (km)
    real(r8), parameter                :: mu_soot_trop = 3._r8  ! width parameter, soot, tropopause (km)

    integer       :: tropLev(pcols)           ! tropopause level index   
    real(r8)      :: tropP(pcols)             ! tropopause pressure (Pa)  
    real(r8)      :: tropT(pcols)             ! tropopause temperature (K) 
    real(r8)      :: tropZ(pcols)             ! tropopause height (m) 

    real(r8)     :: lon(pcols)              ! longitude
    real(r8)     :: lat(pcols)              ! latitude
    integer      :: igroup                  ! group index
    integer      :: lchnk                   ! chunk identifier
    integer      :: ncol                    ! number of columns in chunk
    integer      :: icol                    ! column index
    integer      :: irgn                    ! emission region index
    integer      :: k                       ! vertical index
    real(r8)     :: calday                  ! current calendar day
    integer      :: currentDate             ! current date (yyyydoy)
    integer      :: yr                      ! year
    integer      :: mon                     ! month
    integer      :: day                     ! day of month
    integer      :: ncsec                   ! time of day (seconds)
    integer      :: doy                     ! day of year
    real(r8)     :: startyear               ! start year
    real(r8)     :: stopyear                ! stop year
    real(r8)     :: startdoy                ! start year
    real(r8)     :: stopdoy                 ! stop year
    integer      :: emis_time               ! length of time for emission
    real(r8)     :: vfunc(pver)             ! scaling factor to preserve total emission
    character(len=32) :: shortname          ! the shortname of the group
    real(r8)     :: zmid                    ! layer midpoint altitude (km)
    real(r8)     :: ztrop                   ! tropopause altitude (km)
    real(r8)     :: rate                    ! emission rate (kg/s/m)
    real(r8)     :: massflux                ! mass flux (kg/m3/s)
    real(r8)     :: thickness               ! layer thickness (m)
    real(r8)     :: dist                    ! great circle distance
    integer      :: elapsed_dtime           ! time since start of event (sec)
    real(r8)     :: dlat
    real(r8)     :: dlon
    integer      :: ilat
    integer      :: ilon
    real(r8)     :: mnlon
    real(r8)     :: mxlon
    ! Default return code.
    rc = RC_OK

    ! Determine the day of year.
    calday = get_curr_calday()
    if ( is_perpetual() ) then
      call get_perp_date(yr, mon, day, ncsec)
    else
      call get_curr_date(yr, mon, day, ncsec)
    end if
    doy = floor(calday)

    ! Determine the latitude and longitude of each column.
    lchnk = state%lchnk
    ncol = state%ncol
    
    lat(:ncol)  = state%lat(:ncol)
    lon(:ncol)  = state%lon(:ncol)

    ! Assume a regular grid.
    dlat = PI / (plat-1)
    dlon = 2._f*PI / plon

    ! Add any surface flux here.
    surfaceFlux(:ncol) = 0.0_f
    
    ! For emissions into the atmosphere, put the emission here.
    !
    ! Use Toon et al. [2012] as the source function for soot and dust
    ! from a 1 km meteor impact.
    !
    ! For soot, it is assumed that the soot is emitted in one column
    ! containing the impact and that there are two gaussian
    ! distributions: one centered at the surface and one centered at
    ! the tropopause. The emission rate of soot is given as g/s/km and
    ! we assume that the total mass is delivered in one time step.
    !
    ! NOTE: Perhaps some of these fields should end up in the CARMA
    ! model namelist, so different experiments can be run more easily.
    tendency(:ncol, :pver) = 0.0_r8
    
    ! Determine the start and stop year and day of year from the namelist
    ! variables.
    currentDate = yr * 1000 + doy
    startyear = carma_emis_startdate / 1000
    stopyear  = carma_emis_stopdate  / 1000
    
    startdoy  = mod(carma_emis_startdate, 1000)
    stopdoy   = mod(carma_emis_stopdate, 1000)

    ! Make sure to emit for at least one timestep and in multiples of the time
    ! step length.
    ! TBD - This has a leap year problem, but works otherwise ...
    carma_emis_dtime = INT((((stopyear - startyear) * 365._f + (stopdoy - startdoy)) * 24._f * 3600._f + &
                 (carma_emis_stoptime - carma_emis_starttime)) / dt) * dt

    elapsed_dtime = INT((((yr - startyear) * 365._f + (doy - startdoy)) * 24._f * 3600._f + &
               (ncsec - carma_emis_starttime)) / dt) * dt
  
    call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup)
    if (RC < RC_ERROR) return
    
    call CARMAGROUP_GET(carma, igroup, rc, shortname=shortname)
    if (RC < RC_ERROR) return
      









!
!    ! Do dust and soot separately.
!    if ((shortname == "CRDUST") .and. (carma_emis_dust > 0._f)) then
!
!      ! For dust, it is assumed to all be injected in the first time step after the
!      ! impact.
!      if (((currentDate > carma_emis_startdate) .or. &
!           ((currentDate == carma_emis_startdate) .and. (ncsec >= carma_emis_starttime))) .and. &
!          ((currentDate < carma_emis_stopdate) .or. &
!           ((currentDate == carma_emis_stopdate) .and. (ncsec < carma_emis_stoptime))) .and. &
!          (elapsed_dtime < carma_emis_dust_dtime)) then
!  
!        ! For vaporized dust, all of the mass goes into the first size bin (20 nm).
!        if (ibin == 1) then
!          ! Loop over all of the columns.
!          do icol = 1, ncol
!  
!            ! Calculate the great circle distance (in km).
!            dist = REARTH  / 1e5_f * &
!                   abs(acos(sin(lat(icol))*sin(carma_emis_ctrlat*DEG2RAD) + &
!                   cos(lat(icol))*cos(carma_emis_ctrlat*DEG2RAD)*cos(lon(icol)-(carma_emis_ctrlon*DEG2RAD))))
!                    
!            ! Is the column one of the ones over which there should be emissions>
!            if (dist <= carma_emis_dust_radius) then
!  
!              ! Set tendencies for any sources or sinks in the atmosphere.
!              do k = 1, pver
!            
!                ! Get the cell midpoint and height
!                zmid  = state%zm(icol, k) / 1000._f
!    
!                ! Determine the total emission rate for this grid box using equation 2
!                ! from Toon et al. [2012] and also adjust for the fraction of the
!                ! mass that goes into the specified bin based on the assumed size
!                ! distribution also from Toon et al. [2012]. This is a Gaussian
!                ! centered at 70 km.
!                vfunc(k) = 1._f / (mu_dust * sqrt(2._f * PI)) * &
!                           exp(-0.5_f * (((zmid - carma_emis_z_ballistic) / mu_dust)**2)) * &
!                       (state%zi(icol, k) - state%zi(icol, k+1))
!                       
!                rate = carma_emis_dust
!              
!                ! Calculate a rate by dividing by total emission time.
!                rate = rate  * vfunc(k) / carma_emis_dust_dtime
!              
!!                ! Scale for the fraction of the total surface area that is emitting and
!                ! convert to kg/m2/s
!                massflux = rate / carma_emis_dust_area
!            
!                ! Convert the mass flux to a tendency on the mass mixing ratio.
!                tendency(icol, k) = massflux / (state%pdel(icol, k) / gravit)
!!              end do
!          
!              ! Now normalize in the vertical to preserve the tota mass.
!              tendency(icol, :) = tendency(icol, :) / sum(vfunc(:))
!            end if
!          end do
!        end if 
!      end if


   if (shortname == "CRSOOT") then!  .and. (carma_emis_soot > 0._f)) then

      ! For soot it is assumed to be injected during the entire time period.
      if (((currentDate > carma_emis_startdate) .or. &
           ((currentDate == carma_emis_startdate) .and. (ncsec >= carma_emis_starttime))) .and. &
          ((currentDate < carma_emis_stopdate) .or. &
           ((currentDate == carma_emis_stopdate) .and. (ncsec < carma_emis_stoptime)))) then

        ! Loop over all of the columns.
        !write(*,*) "WE GOT THIS FAR",currentDate
        do icol = 1, ncol

          ! Check each region. Regions should not overlap.
          do irgn = 1, nRegions

!write(*,*) "tend", tendency(icol,1)
!write(*,*) "vfunc", sum(vfunc(:))
            ! grid lons are 0-360, but specification can be -180 to +180, so convert lons
            mnlon =  carma_emis_rgnMinLon(irgn)
            if (mnlon .lt. 0._r8) mnlon = 360._r8 + mnlon
            mnlon = mnlon*DEG2RAD

            mxlon =  carma_emis_rgnMaxLon(irgn)
            if (mxlon .lt. 0._r8) mxlon = 360._r8 + mxlon
            mxlon = mxlon*DEG2RAD

            ! Determine if this column is in the range of the regions.
            if ((lat(icol) .ge. carma_emis_rgnMinLat(irgn)*DEG2RAD) .and. (lat(icol) .le. carma_emis_rgnMaxLat(irgn)*DEG2RAD) .and. &
                (lon(icol) .ge. mnlon) .and. (lon(icol) .le. mxlon)) then

              ! In the vertical, split the mass evenly between levels, so scale mixing
              ! ratio by 1/density
              vfunc(:) = 0._r8

              do k = 1, pver
                ! state is in Pa, and namelist is in hPa
                if ((state%pmid(icol, k) .ge. carma_emis_rgnMinP(irgn) * 100._r8) .and. &
                    (state%pmid(icol, k) .le. carma_emis_rgnMaxP(irgn) * 100._r8)) then
                  vfunc(k) = 1._r8 / state%pdel(icol, k)
                end if

                ! The rate is a combination of the overall soot, the fraction of soot in
                ! this region, the temporal scaling of the soot, and the bin factor.
!                rate = carma_emis_soot * &
!                       carma_emis_rgnFraction(irgn) * &
!                       2._r8 * max(0._r8, (1._r8 - elapsed_dtime / carma_emis_dtime)) * &
!                       carma_sootbinfactor(ibin)
                rate = carma_emis_soot * &
                       carma_emis_rgnFraction(irgn) * &
                       2._r8 * max(0._r8, (1._r8 - (2._r8 * elapsed_dtime + dt) / (2._r8 * carma_emis_dtime))) * &
                       carma_sootbinfactor(ibin)

                ! Calculate a rate by dividing by total emission time.
                rate = rate  * vfunc(k) / carma_emis_dtime

                ! Scale for the fraction of the total surface area that is emitting and
                ! convert to kg/m2/s
                massflux = rate / carma_emis_rgnArea(irgn)

                ! If only emitting over land, then scale by land fraction.
                if (carma_emis_fire_land_only) then
                  massflux = massflux * cam_in.landfrac(icol)
                end if

                ! Convert the mass flux to a tendency on the mass mixing ratio.
                tendency(icol, k) = massflux / (state%pdel(icol, k) / gravit)
              end do

              ! Now normalize in the vertical to preserve the total mass.
              tendency(icol, :) = tendency(icol, :) / sum(vfunc(:))
            end if
          end do
        end do
      end if
    end if






! BELOW WAS ASTEROID IMPACT CODE

!    else if ((shortname == "CRSOOT")  .and. (carma_emis_soot > 0._f)) then
!            
!      ! For soot it is assumed to be injected during the entire time period.
!      if (((currentDate > carma_emis_startdate) .or. &
!           ((currentDate == carma_emis_startdate) .and. (ncsec >= carma_emis_starttime))) .and. &
!          ((currentDate < carma_emis_stopdate) .or. &
!           ((currentDate == carma_emis_stopdate) .and. (ncsec < carma_emis_stoptime)))) then
!
!        ! Find the tropopause using the default algorithm backed by the climatology.
!        call tropopause_find(state, tropLev, tropZ=tropZ)
!  
!        ! Loop over all of the columns.
!        do icol = 1, ncol
!
!          ! Calculate the great circle distance (in km)
!          dist = REARTH / 1e5_f * &
!                 abs(acos(sin(lat(icol))*sin(carma_emis_ctrlat*DEG2RAD) + &
!                 cos(lat(icol))*cos(carma_emis_ctrlat*DEG2RAD)*cos(lon(icol)-(carma_emis_ctrlon*DEG2RAD))))
!                 
!         ! Is the column one of the ones over which there should be emissions>
!          if (dist <= carma_emis_soot_radius) then
!
!            ! NOTE: Assume a regular grids for now.
!            ilat = nint((lat(icol) + (PI / 2._f)) / dlat) + 1
!            ilon = nint(lon(icol) / dlon) + 1
!
!            ! Set tendencies for any sources or sinks in the atmosphere.
!            do k = 1, pver
!            
!              ! Get the cell midpoint and height
!              zmid  = state%zm(icol, k) / 1000._f
!    
!              ! Get the tropopause height.
!              ztrop = tropZ(icol) / 1000._f
!
!              ! Determine the total emission rate for this grid box using equation 2
!              ! from Toon et al. [2012] and also adjust for the fraction of the
!              ! mass that goes into the specified bin based on the assumed size
!              ! distribution also from Toon et al. [2012].
!! Ballistic soot
!!              vfunc(k) = 1._f / (mu_dust * sqrt(2._f * PI)) * &
!!                         exp(-0.5_f * (((zmid - carma_emis_z_ballistic) / mu_dust)**2)) * &
!!                     (state%zi(icol, k) - state%zi(icol, k+1))
!              vfunc(k) = 1._f / sqrt(2._f * PI) * &
!                    ((1._f - carma_emis_fine_trop_frac) / mu_soot_gnd  * exp(-0.5_f * ((zmid / mu_soot_gnd)**2)) + &
!                     carma_emis_fine_trop_frac / mu_soot_trop * exp(-0.5_f * (((zmid - ztrop) / mu_soot_trop)**2))) * &
!                     (state%zi(icol, k) - state%zi(icol, k+1))
!              
!              rate = carma_emis_soot * carma_emis_fire_fine_frac * carma_sootbinfactor(ibin)
!                
!              ! Apply a linearly decreasing trend on the rate.
!              rate = rate * 2._r8 * max(0._r8, (1._r8 - elapsed_dtime / (carma_emis_dtime-dt)))
!
!              ! Calculate a rate by dividing by total emission time.
!              rate = rate  * vfunc(k) / carma_emis_dtime
!                
!              ! Scale for the fraction of the total surface area that is emitting and
!              ! convert to kg/m2/s
!              massflux = rate / carma_emis_soot_area
!              
!              ! If only emitting proportional to vegetation, then scale by veg fraction.
!              if (carma_emis_fire_veg_scale) then
!                massflux = massflux * carma_emis_soot_area / carma_emis_tot_veg * carma_vegfrac_map(ilon, ilat)
!             ! If only emitting over land, then scale by land fraction.
!             else if (carma_emis_fire_land_only) then
!                massflux = massflux * cam_in.landfrac(icol)
!              end if
!              
!              ! Convert the mass flux to a tendency on the mass mixing ratio.
!              tendency(icol, k) = massflux / (state%pdel(icol, k) / gravit)
!            end do
!            
!            ! Now normalize in the vertical to preserve the total mass.
!            tendency(icol, :) = tendency(icol, :) / sum(vfunc(:))
!         end if
!        end do
!      end if


!    else if (shortname == "CRCSOT") then
!            
!      ! For soot it is assumed to be injected during the entire time period.
!      if (((currentDate > carma_emis_startdate) .or. &
!           ((currentDate == carma_emis_startdate) .and. (ncsec >= carma_emis_starttime))) .and. &
!          ((currentDate < carma_emis_stopdate) .or. &
!           ((currentDate == carma_emis_stopdate) .and. (ncsec < carma_emis_stoptime)))) then
!
!        ! Find the tropopause using the default algorithm backed by the climatology.
!        call tropopause_find(state, tropLev, tropZ=tropZ)
!  
!        ! Loop over all of the columns.
!        do icol = 1, ncol
!
!          ! Calculate the great circle distance (in km)
!          dist = REARTH / 1e5_f * &
!                 abs(acos(sin(lat(icol))*sin(carma_emis_ctrlat*DEG2RAD) + &
!                 cos(lat(icol))*cos(carma_emis_ctrlat*DEG2RAD)*cos(lon(icol)-(carma_emis_ctrlon*DEG2RAD))))
!                  
!          ! Is the column one of the ones over which there should be emissions>
!          if (dist <= carma_emis_soot_radius) then
!
!            ! NOTE: Assume a regular grids for now.
!            ilat = nint((lat(icol) + (PI / 2._f)) / dlat) + 1
!            ilon = nint(lon(icol) / dlon) + 1
!
!            ! Set tendencies for any sources or sinks in the atmosphere.
!            do k = 1, pver
!            
!              ! Get the cell midpoint and height
!              zmid  = state%zm(icol, k) / 1000._f
!    
!              ! Get the tropopause height.
!              ztrop = tropZ(icol) / 1000._f
!
!              ! Determine the total emission rate for this grid box using equation 2
!              ! from Toon et al. [2012] and also adjust for the fraction of the
!              ! mass that goes into the specified bin based on the assumed size
!              ! distribution also from Toon et al. [2012].
!! Ballistic soot
!!              vfunc(k) = 1._f / (mu_dust * sqrt(2._f * PI)) * &
!!                         exp(-0.5_f * (((zmid - carma_emis_z_ballistic) / mu_dust)**2)) * &
!!                     (state%zi(icol, k) - state%zi(icol, k+1))
!              vfunc(k) = 1._f / sqrt(2._f * PI) * &
!                    ((1._f - carma_emis_coarse_trop_frac) / mu_soot_gnd  * exp(-0.5_f * ((zmid / mu_soot_gnd)**2)) + &
!                     carma_emis_coarse_trop_frac / mu_soot_trop * exp(-0.5_f * (((zmid - ztrop) / mu_soot_trop)**2))) * &
!                     (state%zi(icol, k) - state%zi(icol, k+1))
!              
!              rate = carma_emis_soot * (1._f - carma_emis_fire_fine_frac) * carma_csootbinfactor(ibin)
!                
!              ! Apply a linearly decreasing trend on the rate.
!              rate = rate * 2._r8 * max(0._r8, (1._r8 - elapsed_dtime / (carma_emis_dtime-dt)))
!
!              ! Calculate a rate by dividing by total emission time.
!              rate = rate  * vfunc(k) / carma_emis_dtime
!                
!              ! Scale for the fraction of the total surface area that is emitting and
!              ! convert to kg/m2/s
!              massflux = rate / carma_emis_soot_area
!              
!              ! If only emitting proportional to vegetation, then scale by veg fraction.
!              if (carma_emis_fire_veg_scale) then
!                massflux = massflux * carma_emis_soot_area / carma_emis_tot_veg * carma_vegfrac_map(ilon, ilat)
!
!              ! If only emitting over land, then scale by land fraction.
!              else if (carma_emis_fire_land_only) then
!                massflux = massflux * cam_in.landfrac(icol)
!              end if
!              
!              ! Convert the mass flux to a tendency on the mass mixing ratio.
!              tendency(icol, k) = massflux / (state%pdel(icol, k) / gravit)
!            end do
!            
!            ! Now normalize in the vertical to preserve the total mass.
!            tendency(icol, :) = tendency(icol, :) / sum(vfunc(:))
!          end if
!        end do
!      end if
!    end if
    
    return
  end subroutine CARMA_EmitParticle


  !! Allows the model to perform its own initialization in addition to what is done
  !! by default in CARMA_init.
  !!
  !! NOTE: If CARMA constituents appear in the initial condition file, then those
  !! values will override anything set here.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeModel(carma, lq_carma, rc)
    use constituents, only: pcnst
    use dyn_grid, only: get_horiz_grid_dim_d, get_horiz_grid_d
    use time_manager, only: is_first_step

    implicit none
    
    type(carma_type), intent(in)       :: carma                 !! the carma object
    logical, intent(inout)             :: lq_carma(pcnst)       !! flags to indicate whether the constituent
                                                                !! could have a CARMA tendency
    integer, intent(out)               :: rc                    !! return code, negative indicates failure
    
    ! NOTE: The dust distribution has not been specified yet, but it should be different
    ! from the soot.
    real(kind=f), parameter            :: rm_dust    = 0.11_f     ! dust mean radius (um)
    real(kind=f), parameter            :: sigma_dust = 1.6_f      ! dust variance
    real(kind=f), parameter            :: rm_soot    = 0.11_f     ! soot mean radius (um)
    real(kind=f), parameter            :: sigma_soot = 1.6_f      ! soot variance
    !real(kind=f), parameter            :: rm_csoot   = 10.0_f     ! soot mean radius (um)
    !real(kind=f), parameter            :: sigma_csoot= 1.46_f     ! soot variance

    real(kind=f), parameter            :: rm_csoot   = 1.85_f     ! soot mean radius (um)
    real(kind=f), parameter            :: sigma_csoot= 1.8_f     ! soot variance

    integer                            :: i
    integer                            :: irgn
    integer                            :: hdim1_d
    integer                            :: hdim2_d
    integer                            :: ngcols
    real(kind=f)                       :: dist                  ! great circle distance
    real(kind=f)                       :: r(NBIN)
    real(kind=f)                       :: dr(NBIN)
    real(kind=f)                       :: rmass(NBIN)
    real(kind=f)                       :: dM(NBIN)
    real(kind=f), allocatable          :: lat(:)
    real(kind=f), allocatable          :: lon(:)
    real(kind=f), allocatable          :: colarea(:)
    real(kind=f), allocatable          :: landfrac(:)
    real(kind=f), allocatable          :: vegfrac(:)
    character(len=32)                  :: shortname             ! the shortname of the group
    real(r8)                           :: mnlon
    real(r8)                           :: mxlon

    integer                            :: LUNOPRT               ! logical unit number for output
    logical                            :: do_print              ! do print output?

  1 format(i3,5x,i3,4x,e10.3,4x,e10.3) 

    ! Default return code.
    rc = RC_OK

    ! Create a mapping of the CARMA dust bins to the dust sizes assumed at the
    ! surface. The sizes of the dust bins at the surface are from Mahowald et al.
    ! [2006].
    !
    !   1 :  0.1 - 1.0 um
    !   2 :  1.0 - 2.5 um
    !   3 :  2.5 - 5.0 um
    !   4 :  5.0 - 10.0 um
    call CARMAGROUP_GET(carma, I_GRP_DUST, rc, r=r)
    if (RC < RC_ERROR) return
    
    do i = 1, NBIN
      if (r(i) .le. 1e-4_f) then
        carma_dustmap(i)  = 1
      else if (r(i) .le. 2.5e-4_f) then
        carma_dustmap(i) = 2
      else if (r(i) .le. 5e-4_f) then
        carma_dustmap(i) = 3
      else
        carma_dustmap(i) = 4
      end if
    end do
    
    ! Determine the weight of mass in each bin based upon the size distribution specified
    ! in Toon et al. [2012], for soot.
    
    call CARMAGROUP_GET(carma, I_GRP_SOOT, rc, shortname=shortname, r=r, dr=dr, rmass=rmass)
    if (RC < RC_ERROR) return
    
    dM(:) = rmass(:) * &
         exp(-(log(r(:) * 1e4_f / rm_soot) ** 2) / (2._f * (log(sigma_soot) ** 2))) / &
         log(sigma_soot) * (dr(:) / r(:))
    carma_sootbinfactor(:)  = dM / sum(dM)

    ! Same thing for the coarse soot.
    call CARMAGROUP_GET(carma, I_GRP_CSOOT, rc, shortname=shortname, r=r, dr=dr, rmass=rmass)
    if (RC < RC_ERROR) return
    
    dM(:) = rmass(:) * &
         exp(-(log(r(:) * 1e4_f / rm_csoot) ** 2) / (2._f * (log(sigma_csoot) ** 2))) / &
         log(sigma_soot) * (dr(:) / r(:))
    carma_csootbinfactor(:)  = dM / sum(dM)

    ! Determine the total area in which debris will be emitted. This is used to scale
    ! the emission per column, based upon the fraction of surface area. This assumes a
    ! regular physics grid.
    call get_horiz_grid_dim_d(hdim1_d, hdim2_d)
  
    ngcols = hdim1_d*hdim2_d
  
    allocate(lat(ngcols))
    allocate(lon(ngcols))
    allocate(colarea(ngcols))
  
    call get_horiz_grid_d(ngcols, clat_d_out=lat, clon_d_out=lon, area_d_out=colarea)

    ! If scaling fires by vegetation, then get vegfrac to the global column list.
    !
    ! NOTE: VEGFRAC must be in the IC file, and the IC is only available on an
    ! initial run. Need to put these results into the restart file for the
    ! general case. This will work for now if emissions are completed in the
    ! initial run.
    if (carma_emis_fire_veg_scale) then
      allocate(vegfrac(ngcols))
    
      if (is_first_step()) then      
        call carma_getVEGFRAC(ngcols, lat, lon, vegfrac, carma_vegfrac_map)
      else
        vegfrac(:) = 0._f
      end if      

    ! If only doing fires on land, then get landfrac to the global column list.
    !
    ! NOTE: LANDFRAC must be in the IC file, and the IC is only available on an
    ! initial run. Need to put these results into the restart file for the
    ! general case. This will work for now if emissions are completed in the
    ! initial run.
    else if (carma_emis_fire_land_only) then
      allocate(landfrac(ngcols))
    
      if (is_first_step()) then      
        call carma_getLANDFRAC(ngcols, lat, lon, landfrac)
      else
        landfrac(:) = 0._f
      end if
    end if
  
    ! rad2 -> m2
    colarea = colarea * REARTH * REARTH / 1e4


    ! Include grid boxes that are within a certain radius of the impact center.
    carma_emis_soot_area = 0._f
    carma_emis_rgnArea(:) = 0._f


    ! Include grid boxes that are within a certain radius of the impact center.
    carma_emis_soot_area = 0._f
    carma_emis_dust_area = 0._f
    carma_emis_splash_area = 0._f
    carma_emis_tot_veg = 0._f
  
    do i = 1, ngcols

      ! Calculate the great circle distance (in km).
      dist = REARTH / 1e5_f * &
             abs(acos(sin(lat(i))*sin(carma_emis_ctrlat*DEG2RAD) + &
                  cos(lat(i))*cos(carma_emis_ctrlat*DEG2RAD)*cos(lon(i)-(carma_emis_ctrlon*DEG2RAD))))
                  
      if (dist <= carma_emis_soot_radius) then
        if (carma_emis_fire_veg_scale) then
          if (vegfrac(i) .gt. 0._f) then
            carma_emis_soot_area = carma_emis_soot_area + colarea(i)
            carma_emis_tot_veg   = carma_emis_tot_veg   + colarea(i) * vegfrac(i)
          end if
        else if (carma_emis_fire_land_only) then
          carma_emis_soot_area = carma_emis_soot_area + colarea(i) * landfrac(i)
        else 
          carma_emis_soot_area = carma_emis_soot_area + colarea(i)
        end if
      end if
    
      if (dist <= carma_emis_dust_radius) then
        carma_emis_dust_area = carma_emis_dust_area + colarea(i)
      end if

      if (dist <= carma_emis_splash_radius) then
        carma_emis_splash_area = carma_emis_splash_area + colarea(i)
      end if

      ! Calculate the soot area in the regions.
      do irgn = 1, nRegions

        ! grid lons are 0-360, but specification can be -180 to +180, so convert lons
        mnlon =  carma_emis_rgnMinLon(irgn)
        if (mnlon .lt. 0._r8) mnlon = 360._r8 + mnlon
        mnlon = mnlon*DEG2RAD

        mxlon =  carma_emis_rgnMaxLon(irgn)
        if (mxlon .lt. 0._r8) mxlon = 360._r8 + mxlon
        mxlon = mxlon*DEG2RAD

        if ((lat(i) .ge. carma_emis_rgnMinLat(irgn)*DEG2RAD) .and. (lat(i) .le. carma_emis_rgnMaxLat(irgn)*DEG2RAD) .and. &
            (lon(i) .ge. mnlon) .and. (lon(i) .le. mxlon)) then

          if (carma_emis_fire_land_only) then
            carma_emis_rgnArea(irgn) = carma_emis_rgnArea(irgn) + colarea(i) * landfrac(i)
          else
            carma_emis_rgnArea(irgn) = carma_emis_rgnArea(irgn) + colarea(i)
          end if
        end if
      end do  
  
      end do
  
    deallocate(lat)
    deallocate(lon)
    deallocate(colarea)
  
    if (carma_emis_fire_veg_scale) then
      deallocate(vegfrac)
    else if (carma_emis_fire_land_only) then
      deallocate(landfrac)
    end if
 
 
    ! Report model specific namelist configuration parameters.
    if (masterproc) then
      call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
      if (rc < 0) call endrun("CARMA_InitializeModel: CARMA_Get failed.")
    
      
      if (do_print) then
        write(LUNOPRT,*) ''
        write(LUNOPRT,*) 'CARMA Initialization ...'
        
        write(LUNOPRT,*) ''
        write(LUNOPRT,*) 'ibin  dustmap  sootfactor csootfactor'

        do i = 1, NBIN
          write(LUNOPRT,1) i, carma_dustmap(i), carma_sootbinfactor(i), carma_csootbinfactor(i)
        end do

        write(LUNOPRT,*) ''
        write(LUNOPRT,*) '  Dust area         :  ', carma_emis_dust_area / 1e6_f, ' (km^2)'
        write(LUNOPRT,*) '  Soot area         :  ', carma_emis_soot_area / 1e6_f, ' (km^2)'
        write(LUNOPRT,*) '  Splash area       :  ', carma_emis_splash_area / 1e6_f, ' (km^2)'
        write(LUNOPRT,*) ''
        
        if (carma_emis_fire_veg_scale) then
          write(LUNOPRT,*) '  Total Veg       :  ', carma_emis_tot_veg
        end if

        write(LUNOPRT,*) ''
        do irgn = 1, nRegions
          write(LUNOPRT,*) '  Region area         :  ', irgn, carma_emis_rgnArea(irgn) / 1e6_f, ' (km^2)'
        end do
        write(LUNOPRT,*) ''




      end if
    end if
    
    return
  end subroutine CARMA_InitializeModel


  !! Sets the initial condition for CARMA aerosol particles. By default, there are no
  !! particles, but this routine can be overridden for models that wish to have an
  !! initial value.
  !!
  !! NOTE: If CARMA constituents appear in the initial condition file, then those
  !! values will override anything set here.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeParticle(carma, ielem, ibin, latvals, lonvals, mask, q, rc)
    use shr_kind_mod,   only: r8 => shr_kind_r8
    use pmgrid,         only: plat, plev, plon

    implicit none
    
    type(carma_type), intent(in)  :: carma      !! the carma object
    integer,          intent(in)  :: ielem      !! element index
    integer,          intent(in)  :: ibin       !! bin index
    real(r8),         intent(in)  :: latvals(:) !! lat in degrees (ncol)
    real(r8),         intent(in)  :: lonvals(:) !! lon in degrees (ncol)
    logical,          intent(in)  :: mask(:)    !! Only initialize where .true.
    real(r8),         intent(out) :: q(:,:)     !! mass mixing ratio (gcol, lev)
    integer,          intent(out) :: rc         !! return code, negative indicates failure

    ! Default return code.
    rc = RC_OK

    ! Add initial condition here.
    
    return
  end subroutine CARMA_InitializeParticle

    
  !!  Called after wet deposition has been performed. Allows the specific model to add
  !!  wet deposition of CARMA aerosols to the aerosols being communicated to the surface.
  !!
  !!  @version July-2011 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_WetDeposition(carma, ielem, ibin, sflx, cam_out, state, rc)
    use camsrfexch,       only: cam_out_t

    implicit none
    
    type(carma_type), intent(in)         :: carma       !! the carma object
    integer, intent(in)                  :: ielem       !! element index
    integer, intent(in)                  :: ibin        !! bin index
    real(r8), intent(in)                 :: sflx(pcols) !! surface flux (kg/m2/s)
    type(cam_out_t), intent(inout)       :: cam_out     !! cam output to surface models
    type(physics_state), intent(in)      :: state       !! physics state variables
    integer, intent(out)                 :: rc          !! return code, negative indicates failure
    
    integer    :: icol
 
    ! Default return code.
    rc = RC_OK
    
    ! Add the wet deposition fluxes to the hydrophilic black carbon.
    !
    ! NOTE: Don't give the surface model negative values for the surface fluxes.
    if (ielem == I_ELEM_SOOT) then
      do icol = 1, state%ncol
        cam_out%bcphiwet(icol) = cam_out%bcphiwet(icol) + max(sflx(icol), 0._r8)
      end do
    end if

    if (ielem == I_ELEM_DUST) then
      do icol = 1, state%ncol
        if (carma_dustmap(ibin) == 1) then
          cam_out%dstwet1(icol) = cam_out%dstwet1(icol) + max(sflx(icol), 0._r8)
        else if (carma_dustmap(ibin) == 2) then
          cam_out%dstwet2(icol) = cam_out%dstwet2(icol) + max(sflx(icol), 0._r8)
        else if (carma_dustmap(ibin) == 3) then
          cam_out%dstwet3(icol) = cam_out%dstwet3(icol) + max(sflx(icol), 0._r8)
        else if (carma_dustmap(ibin) == 4) then
          cam_out%dstwet4(icol) = cam_out%dstwet4(icol) + max(sflx(icol), 0._r8)
        end if
      end do
    end if
    
    return
  end subroutine CARMA_WetDeposition 
  
end module
