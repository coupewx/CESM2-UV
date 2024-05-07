module radxfr_cam
  use shr_kind_mod,   only : r8 => shr_kind_r8
  use ppgrid,         only : pcols, pver
  use ppgrid,         only : begchunk, endchunk
  use cam_logfile,    only : iulog
  use cam_abortutils, only : endrun
  use spmd_utils,     only : masterproc
  use infnan, only : nan, assignment(=)

  use tuv_radiation_transfer, only: tuv_radiation_transfer_init
  use tuv_radiation_transfer, only: tuv_radiation_transfer_run
  use molec_ox_xsect, only: molec_ox_xsect_init
  use molec_ox_xsect, only: molec_ox_xsect_run
  use wavelength_grid, only: nwave
  use interpolate_data, only: lininterp_init, lininterp, interp_type
  use radconstants, only: nswbands
  use tuv_subs, only: fery, futr, addpnt, inter2

  implicit none

  logical :: do_radxfr = .false.
  logical :: has_aer_ra_feedback = .true.
  logical :: do_tuv_photo = .false.
  logical :: do_tuv_cloud = .false.
  logical :: do_sfc_export = .false.
  logical :: do_actinic_co2 = .false.
  logical :: do_actinic_h2o = .false.
  logical :: do_actinic_hbr = .false.
  logical :: do_actinic_hcl = .false.
  logical :: do_actinic_no  = .false.

  real(r8), protected, allocatable :: actinic_fluxes(:,:,:,:) ! (nwave, pver, pcols, begchunk:endchunk )
  character(len=256) :: weights_data_root
  
  integer :: swaertau_idx   = -1
  integer :: swaertauw_idx  = -1
  integer :: swaertauwg_idx = -1
  
  real(r8) :: rrtmg_wavelength(nswbands-1)
  
  type (interp_type) :: interp_wgts
  
  integer, parameter  :: nintegrals = 23
  real(r8), protected, allocatable :: swgt(:,:) ! (nintegrals, nwave)
  character(len=16) :: sname(nintegrals)

  integer, parameter  :: nfluxes = 9
  real(r8), parameter :: fwc(nfluxes) = (/ 0., 121., 180., 210., 240., 308., 400., 600., 999. /)
  integer             :: ituv600 = -1
  integer             :: fwv(nfluxes) = -1
  character(len=16)   :: fname(nfluxes)
  
  
contains
  subroutine radxfr_cam_readnl(nlfile)
    use spmd_utils, only : mpicom, masterprocid, mpi_character, mpi_logical
    use units,      only : getunit, freeunit
    use namelist_utils, only : find_group_name  
    use wavelength_grid, only: wavelength_grid_init

    use params_mod, only: input_data_root

    character(len=*), intent(in) :: nlfile

    character(len=64)  :: radxfr_wavelength_grid_file = 'NONE'
    character(len=265) :: radxfr_input_data_root = 'NONE'
    character(len=512) :: wavelen_grid_filepath
    character(len=512) :: errmsg
    logical            :: radxfr_has_aer_ra_feedback = .false.
    logical            :: radxfr_do_tuv_photo = .false.
    logical            :: radxfr_do_tuv_cloud = .false.
    logical            :: radxfr_do_sfc_export = .false.
    logical            :: radxfr_do_actinic_co2 = .false.
    logical            :: radxfr_do_actinic_h2o = .false.
    logical            :: radxfr_do_actinic_hbr = .false.
    logical            :: radxfr_do_actinic_hcl = .false.
    logical            :: radxfr_do_actinic_no = .false.
    integer :: errflg
    integer :: unitn, ierr

    namelist /photo_radxfr_opts/ radxfr_wavelength_grid_file, radxfr_input_data_root, \
      radxfr_has_aer_ra_feedback, radxfr_do_tuv_photo, radxfr_do_tuv_cloud, \
      radxfr_do_sfc_export, \
      radxfr_do_actinic_co2, radxfr_do_actinic_h2o, radxfr_do_actinic_hbr, \
      radxfr_do_actinic_hcl, radxfr_do_actinic_no

    errflg=0
    errmsg=' '

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       write(*,*) 'read : '//trim(nlfile)
       call find_group_name(unitn, 'photo_radxfr_opts', status=ierr)
       if (ierr == 0) then
          read(unitn, photo_radxfr_opts, iostat=ierr)
          if (ierr /= 0) then
             call endrun('radxfr_cam_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

    call mpi_bcast(radxfr_wavelength_grid_file, len(radxfr_wavelength_grid_file), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(radxfr_input_data_root, len(radxfr_input_data_root), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(radxfr_has_aer_ra_feedback, 1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(radxfr_do_tuv_photo, 1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(radxfr_do_tuv_cloud, 1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(radxfr_do_sfc_export, 1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(radxfr_do_actinic_co2, 1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(radxfr_do_actinic_h2o, 1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(radxfr_do_actinic_hbr, 1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(radxfr_do_actinic_hcl, 1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(radxfr_do_actinic_no, 1, mpi_logical, masterprocid, mpicom, ierr)
 
    do_radxfr = radxfr_wavelength_grid_file /= 'NONE'

    if (.not.do_radxfr) return

    if (masterproc) then    
       write(iulog,*) 'radxfr_cam_readnl: radxfr_input_data_root = '//trim(radxfr_input_data_root)
       write(iulog,*) 'radxfr_cam_readnl: radxfr_wavelength_grid_file = '//trim(radxfr_wavelength_grid_file)
       write(iulog,*) 'radxfr_cam_readnl: radxfr_has_aer_ra_feedback = ', radxfr_has_aer_ra_feedback
       write(iulog,*) 'radxfr_cam_readnl: radxfr_do_tuv_photo = ', radxfr_do_tuv_photo
       write(iulog,*) 'radxfr_cam_readnl: radxfr_do_tuv_cloud = ', radxfr_do_tuv_cloud
       write(iulog,*) 'radxfr_cam_readnl: radxfr_do_sfc_export = ', radxfr_do_sfc_export
       write(iulog,*) 'radxfr_cam_readnl: radxfr_do_actinic_co2 = ', radxfr_do_actinic_co2
       write(iulog,*) 'radxfr_cam_readnl: radxfr_do_actinic_h2o = ', radxfr_do_actinic_h2o
       write(iulog,*) 'radxfr_cam_readnl: radxfr_do_actinic_hbr = ', radxfr_do_actinic_hbr
       write(iulog,*) 'radxfr_cam_readnl: radxfr_do_actinic_hcl = ', radxfr_do_actinic_hcl
       write(iulog,*) 'radxfr_cam_readnl: radxfr_do_actinic_no = ', radxfr_do_actinic_no
    end if

    input_data_root = trim(radxfr_input_data_root)
    wavelen_grid_filepath = trim(input_data_root)//'/'//trim(radxfr_wavelength_grid_file)

    ! call this here since nwave needs to be known earlier than the init phase
    call wavelength_grid_init(wavelen_grid_filepath, errmsg, errflg)
    if(errflg/=0) then
       call endrun('radxfr_cam_readnl: '//trim(errmsg))
    end if

    has_aer_ra_feedback  = radxfr_has_aer_ra_feedback
    
    do_tuv_photo = radxfr_do_tuv_photo
    do_tuv_cloud = radxfr_do_tuv_cloud
    do_sfc_export = radxfr_do_sfc_export
    do_actinic_co2 = radxfr_do_actinic_co2
    do_actinic_h2o = radxfr_do_actinic_h2o
    do_actinic_hbr = radxfr_do_actinic_hbr
    do_actinic_hcl = radxfr_do_actinic_hcl
    do_actinic_no = radxfr_do_actinic_no

    if (has_aer_ra_feedback .and. .not. do_tuv_photo) then
       call endrun('radxfr_cam_readnl: Including aerosols in photolysis requires using in-line TUV.')
    end if

  end subroutine radxfr_cam_readnl

  subroutine radxfr_cam_register
    use physics_buffer, only: pbuf_add_field, dtype_r8
    
    integer :: uva_idx, uvb_idx, par_idx, ephyto1_idx, ephyto2_idx, ephyto4_idx, ephyto5_idx
    
    if (do_sfc_export) then
      call pbuf_add_field('TUV_UVA', 'physpkg', dtype_r8, (/pcols/), uva_idx)
      call pbuf_add_field('TUV_UVB', 'physpkg', dtype_r8, (/pcols/), uvb_idx)
      call pbuf_add_field('TUV_PARE', 'physpkg', dtype_r8, (/pcols/), par_idx)
      call pbuf_add_field('TUV_PHYTOPHAEO92', 'physpkg', dtype_r8, (/pcols/), ephyto1_idx)
      call pbuf_add_field('TUV_PHYTOPRORO92', 'physpkg', dtype_r8, (/pcols/), ephyto2_idx)
      !call pbuf_add_field('TUV_PHYTO94', 'physpkg', dtype_r8, (/pcols/), ephyto3_idx)
      call pbuf_add_field('TUV_PHYTO17', 'physpkg',dtype_r8, (/pcols/),ephyto4_idx)
      call pbuf_add_field('TUV_PHYTO19', 'physpkg',dtype_r8, (/pcols/),ephyto5_idx)
    end if
  	
  end subroutine radxfr_cam_register
  
  subroutine radxfr_cam_init
    use spmd_utils,      only: mpicom, mpir8
    use physics_buffer,  only: pbuf_get_index
    use radconstants,    only: get_sw_spectral_boundaries
    use cam_history,     only: addfld, horiz_only
    use wavelength_grid, only: nwave, wc, wl
    use ioFileMod,       only: getfil
    use params_mod,      only: input_data_root
    use units,           only : getunit, freeunit
    use module_xsections, only: no2xs_a, no2xs_b, o2_xs, so2_xs, co2_xs, hbr_xs, hcl_xs, h2o_xs, no_xs

    integer, parameter  :: kdata = 1000
    real(r8), parameter :: deltax = 1.E-5_r8
    real(r8) :: x1(kdata)
    real(r8) :: y1(kdata)

    integer :: errflg, ierr, kin, fid
    integer :: i, n, iw, j
    character(len=1200) :: errmsg
    real(r8) :: wvl_low(nswbands)
    real(r8) :: wvl_high(nswbands)
    real(r8) :: wdiff
    character(len=288) :: efile
    integer  :: idum
    real(r8) :: dum1, dum2
    
    if (.not.do_radxfr) return

    ! physic buffer fields for aerosol optical properties.
    swaertau_idx   = pbuf_get_index('SWAERTAU')
    swaertauw_idx  = pbuf_get_index('SWAERTAUW')
    swaertauwg_idx = pbuf_get_index('SWAERTAUWG')

    allocate( actinic_fluxes(nwave, pver, pcols, begchunk:endchunk ) )
    actinic_fluxes = nan

    allocate( swgt(nintegrals, nwave ) )
    
    errflg=0
    errmsg=' '

    call molec_ox_xsect_init( errmsg, errflg )
    if (errflg/=0) then
       call endrun('radxfr_cam_init: '//trim(errmsg))
    end if

    call tuv_radiation_transfer_init( r8, errmsg, errflg )
    if (errflg/=0) then
       call endrun('radxfr_cam_init: '//trim(errmsg))
    end if
    
    ! Output the cross sections used for actinic flux calculations.
    if (masterproc) then
        write(iulog,*) "Cross-sections of TUV Actinic flux"
        write(iulog,*) "WV, NO2a, NO2b, O2, SO2, CO2, HBR, HCL, H2O, NO"
        do i = 1, nwave
          write(iulog,'(f6.2,2x,9(e11.4,2x))') wc(i), no2xs_a(i), no2xs_b(i), o2_xs(i), so2_xs(i), co2_xs(i), hbr_xs(i), hcl_xs(i), h2o_xs(i), no_xs(i)
        end do
    end if
    
    
    ! Get the RRTMG wavenumber edges and convert to a wavelength center.
    !
    ! NOTE: Last band is a broadband that overlaps the other bands, so skip it.
    call get_sw_spectral_boundaries(wvl_low, wvl_high, "nm")
    rrtmg_wavelength = (wvl_low(1:nswbands-1) + wvl_high(1:nswbands-1)) / 2._r8

    ! Calculate weights needed to interpolate from the RRTMG wavelengths to the
    ! radxfr wavelengths.
    call lininterp_init(rrtmg_wavelength, nswbands-1, wc, nwave, 1, interp_wgts)

    write(iulog, *) 'radxfr: rrtmg wavelengths: ', rrtmg_wavelength
    write(iulog, *) 'radxfr: tuv wavelengths: ', wc 
    
    ! Find the point closest to 600 nm on the TUV wavelength grid
    fname(1) = 'TUV_AFLXLOW'
    fwv(1) = 1
    
    fname(nfluxes) = 'TUV_AFLXHIGH'
    fwv(nfluxes) = nwave
    
    fname(2) = 'TUV_AFLX121'
    fname(3) = 'TUV_AFLX180'
    fname(4) = 'TUV_AFLX210'
    fname(5) = 'TUV_AFLX240'
    fname(6) = 'TUV_AFLX308'
    fname(7) = 'TUV_AFLX400'
    fname(8) = 'TUV_AFLX600'
    
    do j = 2, nfluxes-1
      wdiff = 1000._r8
      fwv(j) = -1
      do i = 1, nwave
        if (abs(wc(i) - fwc(j)) .lt. wdiff) then
          fwv(j) = i
          wdiff = abs(wc(i) - fwc(j))
        end if
      end do

      write(iulog, *) 'radxfr: closest band to ', fwc(j), ' nm ', fwv(j), wc(fwv(j))
    end do
    
    ituv600 = fwv(8)

    ! Add output fields for the tuv aerosol optics in each of the bands: 300, 400, 600, 999 nm.
    call addfld('TUV_TAULOW', (/ 'lev' /), 'A', '1', 'TUV aerosol optical depth, shortest', flag_xyfill=.true.)
    call addfld('TUV_WLOW',   (/ 'lev' /), 'A', '1', 'TUV single scatter albedo, shortest',flag_xyfill=.true.)
    call addfld('TUV_GLOW',   (/ 'lev' /), 'A', '1', 'TUV asymmetry factor, shortest', flag_xyfill=.true.)

    call addfld('TUV_TAU600', (/ 'lev' /), 'A', '1', 'TUV aerosol optical depth, 600 nm', flag_xyfill=.true.)
    call addfld('TUV_W600',   (/ 'lev' /), 'A', '1', 'TUV single scatter albedo, 600 nm', flag_xyfill=.true.)
    call addfld('TUV_G600',   (/ 'lev' /), 'A', '1', 'TUV asymmetry factor, 600 nm', flag_xyfill=.true.)

    call addfld('TUV_TAUHIGH',(/ 'lev' /), 'A', '1', 'TUV aerosol optical depth, longest', flag_xyfill=.true.)
    call addfld('TUV_WHIGH',  (/ 'lev' /), 'A', '1', 'TUV single scatter albedo, longest', flag_xyfill=.true.)
    call addfld('TUV_GHIGH',  (/ 'lev' /), 'A', '1', 'TUV asymmetry factor, longest', flag_xyfill=.true.)

    call addfld(trim(fname(1)), (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, shortest', flag_xyfill=.true.)
    call addfld(trim(fname(2)), (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, 121 nm', flag_xyfill=.true.)
    call addfld(trim(fname(3)), (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, 180 nm', flag_xyfill=.true.)
    call addfld(trim(fname(4)), (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, 210 nm', flag_xyfill=.true.)
    call addfld(trim(fname(5)), (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, 240 nm', flag_xyfill=.true.)
    call addfld(trim(fname(6)), (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, 308 nm', flag_xyfill=.true.)
    call addfld(trim(fname(7)), (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, 400 nm', flag_xyfill=.true.)
    call addfld(trim(fname(8)),  (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, 600 nm', flag_xyfill=.true.)
    call addfld(trim(fname(nfluxes)), (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, longest', flag_xyfill=.true.)

    
    ! Create the weights for the spectral integrals and add the fields

    ! UV-A (315-400 nm)
    j = 1
    sname(j) = 'TUV_UVA'
    call addfld(trim(sname(j)), horiz_only, 'A', 'W m-2', 'UV-A, 315-400 nm', flag_xyfill=.true.)
    do iw = 1, nwave
       if (wc(iw) .gt. 315._r8 .and. wc(iw) .lt. 400._r8) then
          swgt(j,iw) = 1._r8
       else
          swgt(j,iw) = 0._r8
       end if
    end do

    ! UV-B (280-330 nm)
    j = j + 1
    sname(j) = 'TUV_UVB'
    call addfld(trim(sname(j)), horiz_only, 'A', 'W m-2', 'UV-B (280-315 nm)', flag_xyfill=.true.)
    do iw = 1, nwave
       if (wc(iw) .gt. 280._r8 .and. wc(iw) .lt. 315._r8) then
          swgt(j,iw) = 1._r8
       else
          swgt(j,iw) = 0._r8
       end if
    end do

    ! UV-B* (280-315 nm)
    j = j + 1
    sname(j) = 'TUV_UVBSTAR'
    call addfld(trim(sname(j)), horiz_only, 'A', 'W m-2', 'UV-B* (280-320 nm)', flag_xyfill=.true.)
    do iw = 1, nwave
       if (wc(iw) .gt. 280._r8 .and. wc(iw) .lt. 320._r8) then
          swgt(j,iw) = 1._r8
       else
          swgt(j,iw) = 0._r8
       end if
    end do
    
    ! UV-C (100-280 nm)
    ! NOTE: Our current minimum wavelength in CESM is 121 nm, so this won't cover  the
    ! full UV-C range.
    j = j + 1
    sname(j) = 'TUV_UVC'
    call addfld(trim(sname(j)), horiz_only, 'A', 'W m-2', 'UV-C (100-280 nm)', flag_xyfill=.true.)
    do iw = 1, nwave
       if (wc(iw) .gt. 100._r8 .and. wc(iw) .lt. 280._r8) then
          swgt(j,iw) = 1._r8
       else
          swgt(j,iw) = 0._r8
       end if
    end do

    ! Visible (>400 nm)
    j = j + 1
    sname(j) = 'TUV_VIS'
    call addfld(trim(sname(j)), horiz_only, 'A', 'W m-2', 'Visible (>400 nm)', flag_xyfill=.true.)
    do iw = 1, nwave
       if (wc(iw) .gt. 400._r8) then
          swgt(j,iw) = 1._r8
       else
          swgt(j,iw) = 0._r8
       end if
    end do

   ! Photosynthetic Active Radiation (400 < PAR < 700 nm)
   ! conversion to micro moles m-2 s-1:
   !   s = s * (1e6/6.022142E23)(w/1e9)/(6.626068E-34*2.99792458E8)
    j = j + 1
    sname(j) = 'TUV_PAR'
    call addfld(trim(sname(j)), horiz_only, 'A', 'umol m-2 s-1', 'Photosynthetic Active Radiation (400 < PAR < 700 nm)', flag_xyfill=.true.)
    do iw = 1, nwave
       if (wc(iw) .gt. 400._r8 .and. wc(iw) .lt. 700._r8) then
          swgt(j,iw) = 8.36e-3_r8 * wc(iw)
       else
          swgt(j,iw) = 0._r8
       end if
    end do

    ! Photosynthetic Active Radiation (400 < PAR < 700 nm), but in energy units
    j = j + 1
    sname(j) = 'TUV_PARE'
    call addfld(trim(sname(j)), horiz_only, 'A', 'W m-2', 'Photosynthetic Active Radiation (400 < PAR < 700 nm)', flag_xyfill=.true.)
    do iw = 1, nwave
       if (wc(iw) .gt. 400._r8 .and. wc(iw) .lt. 700._r8) then
          swgt(j,iw) = 1._r8
       else
          swgt(j,iw) = 0._r8
       end if
    end do

    ! UV index (Canadian - WMO/WHO)
    ! Report of the WMO Meeting of experts on UV-B measurements, data quality 
    ! and standardization of UV indices, World Meteorological Organization
    ! (WMO), report No. 95, Geneva, 1994.
    ! based on the CIE erythema weighting, multiplied by 40.
    j = j + 1
    sname(j) = 'TUV_UVINDEX'
    call addfld(trim(sname(j)), horiz_only, 'A', 'W m-2', 'UV index (Canadian - WMO/WHO)', flag_xyfill=.true.)
    do iw = 1, nwave
       swgt(j,iw) = 40._r8 * fery(wc(iw))
    end do

    ! UV index (Canadian - WMO/WHO)
    ! Report of the WMO Meeting of experts on UV-B measurements, data quality 
    ! and standardization of UV indices, World Meteorological Organization
    ! (WMO), report No. 95, Geneva, 1994.
    ! based on the CIE erythema weighting, multiplied by 40.
    j = j + 1
    sname(j) = 'TUV_UVINDEXMX'
    call addfld(trim(sname(j)), horiz_only, 'X', '1', 'UV index (Canadian - WMO/WHO), Maximum', flag_xyfill=.true.)
    do iw = 1, nwave
       swgt(j,iw) = 40._r8 * fery(wc(iw))
    end do
    
    ! skin cancer in mice,  Utrecht/Phildelphia study
    ! from de Gruijl, F. R., H. J. C. M. Sterenborg, P. D. Forbes, 
    ! R. E. Davies, C. Cole, G. Kelfkens, H. van Weelden, H. Slaper,
    ! and J. C. van der Leun, Wavelength dependence of skin cancer 
    ! induction by ultraviolet irradiation of albino hairless mice, 
    ! Cancer Res., 53, 53-60, 1993.
    ! Calculate with function futr(w), normalize at 300 nm.
    j = j + 1
    sname(j) = 'TUV_SCUPMICE93'
    call addfld(trim(sname(j)), horiz_only, 'A', '1', 'SCUP-mice (de Gruijl et al., 1993)', flag_xyfill=.true.)
    do iw = 1, nwave
       swgt(j,iw) = futr(wc(iw)) / futr(300._r8)
    end do

    ! CIE standard human erythema action spectrum
    ! from:
    ! McKinlay, A. F., and B. L. Diffey, A reference action spectrum for 
    ! ultraviolet induced erythema in human skin, in Human Exposure to 
    ! Ultraviolet Radiation: Risks and Regulations, W. R. Passchler 
    ! and B. F. M. Bosnajokovic, (eds.), Elsevier, Amsterdam, 1987.
    j = j + 1
    sname(j) = 'TUV_ERYTCIE87'
    call addfld(trim(sname(j)), horiz_only, 'A', '1', 'CIE human erythema (McKinlay and Diffey, 1987)', flag_xyfill=.true.)
    do iw = 1, nwave
       swgt(j,iw) = fery(wc(iw))
    end do

    ! phytoplankton, Boucher et al. (1994) 
    ! from Boucher, N., Prezelin, B.B., Evens, T., Jovine, R., Kroon, B., Moline, M.A.,
    ! and Schofield, O., Icecolors '93: Biological weighting function for the ultraviolet
    ! inhibition  of carbon fixation in a natural antarctic phytoplankton community, 
    ! Antarctic Journal, Review 1994, pp. 272-275, 1994.
    ! In original paper, value of b and m (em below are given as positive.  Correct values
    ! are negative. Also, limit to positive values.
    !j = j + 1
    !sname(j) = 'TUV_PHYTO94'
    !call addfld(trim(sname(j)), horiz_only, 'A', '1', 'Phytoplankton (Boucher et al., 1994)', flag_xyfill=.true.)
    !do iw = 1, nwave
    !   if (wc(iw) .gt. 290._r8 .and. wc(iw) .lt. 400._r8) then
    !      swgt(j,iw) = -3.17e-6_r8 + EXP(112.5_r8 + -6.223e-01_r8 * wc(iw) + 7.670E-04_r8 * wc(iw)*wc(iw))
    !   else
    !      swgt(j,iw) = 0._r8
    !   endif
    !   swgt(j,iw) = max(swgt(j,iw), 0._r8)
    !end do


    ! prochlorococcus, Neale and Thomas (2017)
    ! from paper examining high and mid latitude picophytoplankton
    j = j+1
    sname(j)='TUV_PHYTO17'
    call addfld(trim(sname(j)), horiz_only, 'A', '1', 'Phytoplankton (Neale and Thomas, 2017)',flag_xyfill=.true.)
    do iw = 1, nwave
        if (wc(iw) .gt. 280._r8 .and. wc(iw) .lt. 400._r8) then
           swgt(j,iw) =  -3.10572006e-15_r8*wc(iw)**6 + 6.27284524e-12_r8*wc(iw)**5 - 5.22236845e-9_r8*wc(iw)**4 + 2.28848111e-6_r8*wc(iw)**3-5.54843743e-4_r8*wc(iw)**2+7.02209062e-2_r8*wc(iw) - 3.59639531
        else
          swgt(j,iw) = 0._r8
       endif
       swgt(j,iw) = max(swgt(j,iw), 0._r8)
    end do
    swgt(j,:) = swgt(j,:) * 1000._r8 ! equation is in units (mW m^-2)^-1 and spectral irradiance is in units W m^-2 ... wait, does this make sense?

    !coccolithophores, Lorenzo et al. (2019)
    ! from paper examining Emiliana Huxleyi
    j = j+1
    sname(j)='TUV_PHYTO19'
    call addfld(trim(sname(j)), horiz_only, 'A', '1', 'Coccolithophores (Lorenzo et al., 2019)',flag_xyfill=.true.)
    do iw = 1, nwave
        swgt(j,iw) = 0._r8 ! by default set this to 0 
        if (wc(iw) .gt. 280._r8 .and. wc(iw) .lt. 320._r8) then !UV-B
           swgt(j,iw) = -1.40530706e-12_r8*wc(iw)**6 + 2.39353416e-09_r8*wc(iw)**5 - 1.68928198e-06_r8*wc(iw)**4 + 6.31865068e-04_r8*wc(iw)**3 - 1.31975456e-01_r8*wc(iw)**2 + 1.45753274e+01_r8*wc(iw) - 6.63797953e+02_r8
       else if (wc(iw) .gt. 320._r8 .and. wc(iw) .lt. 400_r8) then  !UV-A
           swgt(j,iw) = 8.89581784e-15_r8*wc(iw)**6 - 1.93953328e-11_r8*wc(iw)**5 + 1.76199028e-08_r8*wc(iw)**4 - 8.53798419e-06_r8*wc(iw)**3 + 2.32766669e-03_r8*wc(iw)**2 - 3.38550975e-01_r8*wc(iw) + 2.05261692e+01_r8
       endif
       swgt(j,iw) = max(swgt(j,iw), 0._r8)
    end do
    swgt(j,:) = swgt(j,:) * 1000._r8

    ! Plant damage - Caldwell 1971
    ! Caldwell, M. M., Solar ultraviolet radiation and the growth and 
    ! development of higher plants, Photophysiology 6:131-177, 1971.
    ! Alternative fit to Caldwell (1971) by 
    ! Micheletti, M. I. and R. D. Piacentini, Photochem. Photobiol.,
    ! 76, pp.?, 2002.
    j = j + 1
    sname(j) = 'TUV_PLANT71'
    call addfld(trim(sname(j)), horiz_only, 'A', '1', 'Plant damage (Caldwell, 1971)', flag_xyfill=.true.)
    do iw = 1, nwave
       swgt(j,iw) = 570.25_r8 + -4.70144_r8*wc(iw) + 0.01274_r8 *wc(iw)**2  + -1.13118E-5_r8*wc(iw)**3
       if (swgt(j,iw) .lt. 0._r8 .or. wc(iw) .gt. 313._r8) then
          swgt(j,iw) = 0._r8
       end if
    end do

    ! Plant damage - Flint & Caldwell 2003
    ! Flint, S. D. and M. M. Caldwell, A biological spectral weigthing
    ! function for ozone depletion research with higher plants, Physiologia
    ! Plantorum, in press, 2003.
    ! Data available to 366 nm
    j = j + 1
    sname(j) = 'TUV_PLANT03'
    call addfld(trim(sname(j)), horiz_only, 'A', '1', 'Plant damage (Flint & Caldwell, 2003)', flag_xyfill=.true.)
    do iw = 1, nwave
       swgt(j,iw) = EXP( 4.688272_r8 * EXP( -EXP(0.1703411_r8 * (wc(iw)-307.867_r8) / 1.15_r8))+ \
          ((390._r8-wc(iw))/121.7557_r8 - 4.183832_r8))

       ! put on per joule (rather than per quantum) basis:
       swgt(j,iw) = swgt(j,iw) * wc(iw)/300._r8

       if (swgt(j,iw) .lt. 0._r8 .or. wc(iw) .gt. 366._r8) then
          swgt(j,iw) = 0._r8
       end if
    end do
    
    
    ! All of the following require reading in weights from a file. TUV defines them as text files,
    ! but perhaps for the future they should be converted to NETCDF files. Also, perhaps this should
    ! be done in the masterproc and the broadcast out to the other nodes
    
    ! DNA damage action spectrum
    ! from: Setlow, R. B., The wavelengths in sunlight effective in 
    ! producing skin cancer: a theoretical analysis, Proceedings 
    ! of the National Academy of Science, 71, 3363 -3366, 1974.
    !  normalize to unity at 300 nm
    ! Data read from original hand-drawn plot by Setlow
    ! received from R. Setlow in May 1995
    ! data is per quantum (confirmed with R. Setlow in May 1995).  
    ! Therefore must put on energy basis if irradiance is is energy
    ! (rather than quanta) units.
    j = j + 1
    sname(j) = 'TUV_DNA74'
    call addfld(trim(sname(j)), horiz_only, 'A', '1', 'DNA damage, in vitro (Setlow, 1974)', flag_xyfill=.true.)
    
    ! Read from file.
    if (masterproc) then
       call getfil(trim(input_data_root)//'/'//trim('DATAS1/dna.setlow.new'), efile, fid)
       kin = getunit()
       open(unit=kin,file=efile,status='old')
       do i = 1, 11
          read(kin,*)
       enddo
       n = 55
       do i = 1, n
          read(kin,*) x1(i), y1(i)
          y1(i) = y1(i) / 2.4e-02_r8  *  x1(i)/300._r8
       enddo
       close(kin)
       call freeunit(kin)
       ! Interpolate onto wavelength grid.
       call addpnt(x1,y1,kdata,n,x1(1)*(1._r8-deltax),y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,0._r8,y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,x1(n)*(1._r8+deltax),0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,1.e+38_r8,0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call inter2(nwave,wl,swgt(j,:nwave),n,x1,y1,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to interpolate for TUV_DNA74.')
       endif
       do i=1, nwave
         write(*,*) "TUV_DNA74: ", i, wc(i), swgt(j,i)
       end do
    end if
    
    call mpibcast(swgt(j,:), nwave, mpir8, 0, mpicom)


    ! Utrecht/Philadelphia mice spectrum corrected for humans skin.
    ! From de Gruijl, F.R. and J. C. van der Leun, Estimate of the wavelength 
    ! dependency of ultraviolet carcinogenesis and its relevance to the
    ! risk assessment of a stratospheric ozone depletion, Health Phys., 4,
    ! 317-323, 1994.
    j = j + 1
    sname(j) = 'TUV_SCUPMAN94'
    call addfld(trim(sname(j)), horiz_only, 'A', '1', 'SCUP-human (de Gruijl and van der Leun, 1994)', flag_xyfill=.true.)
    
    ! Read from file.
    if (masterproc) then
       call getfil(trim(input_data_root)//'/'//trim('DATAS1/SCUP-h'), efile, fid)
       kin = getunit()
       open(unit=kin,file=efile,status='old')
       n = 28
       do i = 1, n
          read(kin,*) x1(i), y1(i)
       enddo
       close(kin)
       call freeunit(kin)
       ! Interpolate onto wavelength grid.
       call addpnt(x1,y1,kdata,n,x1(1)*(1._r8-deltax),y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,0._r8,y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,x1(n)*(1._r8+deltax),0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,1.e+38_r8,0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call inter2(nwave,wl,swgt(j,:nwave),n,x1,y1,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to interpolate for TUV_DNA74.')
       endif
    end if
    
    call mpibcast(swgt(j,:), nwave, mpir8, 0, mpicom)
    

    ! Human erythema - Anders et al.
    ! from:
    ! Anders, A., H.-J. Altheide, M. Knalmann, and H. Tronnier,
    ! Action spectrum for erythema in humands investigated with dye lasers, 
    ! Photochem. and Photobiol., 61, 200-203, 1995.
    ! for skin types II and III, Units are J m-2.
    j = j + 1
    sname(j) = 'TUV_ERYH95'
    call addfld(trim(sname(j)), horiz_only, 'A', '1', 'Erythema, humans (Anders et al., 1995)', flag_xyfill=.true.)
    
    ! Read from file.
    if (masterproc) then
       call getfil(trim(input_data_root)//'/'//trim('DATAS1/ery.anders'), efile, fid)
       kin = getunit()
       open(unit=kin,file=efile,status='old')
       do i = 1, 5
         read(kin,*)
       end do
       n = 28
       do i = 1, n
          read(kin,*) x1(i), y1(i)
          y1(i) = 1._r8 / y1(i)
       enddo
       close(kin)
       call freeunit(kin)
       ! Interpolate onto wavelength grid.
       call addpnt(x1,y1,kdata,n,x1(1)*(1._r8-deltax),y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,0._r8,y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,x1(n)*(1._r8+deltax),0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,1.e+38_r8,0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call inter2(nwave,wl,swgt(j,:nwave),n,x1,y1,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to interpolate for TUV_DNA74.')
       endif
    end if
    
    call mpibcast(swgt(j,:), nwave, mpir8, 0, mpicom)


    ! 1991-92 ACGIH threshold limit values
    ! from
    ! ACGIH, 1991-1992 Threshold Limit Values, American Conference 
    ! of Governmental and Industrial Hygienists, 1992.
    j = j + 1
    sname(j) = 'TUV_OTLV92'
    call addfld(trim(sname(j)), horiz_only, 'A', '1', 'Occupational TLV (ACGIH, 1992)', flag_xyfill=.true.)
    
    ! Read from file.
    if (masterproc) then
       call getfil(trim(input_data_root)//'/'//trim('DATAS1/acgih.1992'), efile, fid)
       kin = getunit()
       open(unit=kin,file=efile,status='old')
       n = 56
       do i = 1, n
          read(kin,*) x1(i), y1(i)
       enddo
       close(kin)
       call freeunit(kin)
       ! Interpolate onto wavelength grid.
       call addpnt(x1,y1,kdata,n,x1(1)*(1._r8-deltax),y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,0._r8,y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,x1(n)*(1._r8+deltax),0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,1.e+38_r8,0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call inter2(nwave,wl,swgt(j,:nwave),n,x1,y1,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to interpolate for TUV_DNA74.')
       endif
    end if
    
    call mpibcast(swgt(j,:), nwave, mpir8, 0, mpicom)


    ! phytoplankton, Cullen et al.
    ! Cullen, J.J., Neale, P.J., and Lesser, M.P., Biological weighting function for the  
    ! inhibition of phytoplankton photosynthesis by ultraviolet radiation, Science, 25,
    ! 646-649, 1992.
    ! phaeo
    j = j + 1
    sname(j) = 'TUV_PHYTOPHAEO92'
    call addfld(trim(sname(j)), horiz_only, 'A', '1', 'Phytoplankton, phaeo (Cullen et al., 1992', flag_xyfill=.true.)
    
    ! Read from file.
    if (masterproc) then
       !call getfil(trim(input_data_root)//'/'//trim('DATAS1/phaeo.bio'), efile, fid)
       call getfil(trim('/glade/work/jcoupe/misc_files/phaeo_dataextended.bio'), efile, fid) !extended file to 280-400 nm
       kin = getunit()
       open(unit=kin,file=efile,status='old')
       n = 106
       do i = 1, n
          read(kin,*) idum, dum1, dum2, y1(i)
          x1(i) = (dum1+dum2)/2._r8
       enddo
       close(kin)
       call freeunit(kin)
       ! Interpolate onto wavelength grid.
       call addpnt(x1,y1,kdata,n,x1(1)*(1._r8-deltax),y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,0._r8,y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,x1(n)*(1._r8+deltax),0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,1.e+38_r8,0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call inter2(nwave,wl,swgt(j,:nwave),n,x1,y1,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to interpolate for TUV_DNA74.')
       endif
    end if
    
    ! The units in the files for the spectral weights are (mW/m2)-1, but the spectral
    ! irradiance is in W/m2/nm. Thus we need to scale the weights by 10^3 so that the
    ! resulting integral will be unitlessthe 
    swgt(j,:) = swgt(j,:) * 1000._r8
    
    call mpibcast(swgt(j,:), nwave, mpir8, 0, mpicom)


    ! phytoplankton, Cullen et al.
    ! Cullen, J.J., Neale, P.J., and Lesser, M.P., Biological weighting function for the  
    ! inhibition of phytoplankton photosynthesis by ultraviolet radiation, Science, 25,
    ! 646-649, 1992.
    ! proro
    j = j + 1
    sname(j) = 'TUV_PHYTOPRORO92'
    call addfld(trim(sname(j)), horiz_only, 'A', '1', 'Phytoplankton, proro (Cullen et al., 1992)', flag_xyfill=.true.)
    
    ! Read from file.
    if (masterproc) then
       call getfil(trim(input_data_root)//'/'//trim('DATAS1/proro.bio'), efile, fid)
       kin = getunit()
       open(unit=kin,file=efile,status='old')
       n = 100
       do i = 1, n
          read(kin,*) idum, dum1, dum2, y1(i)
          x1(i) = (dum1+dum2)/2._r8
       enddo
       close(kin)
       call freeunit(kin)
       ! Interpolate onto wavelength grid.
       call addpnt(x1,y1,kdata,n,x1(1)*(1._r8-deltax),y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,0._r8,y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,x1(n)*(1._r8+deltax),0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,1.e+38_r8,0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call inter2(nwave,wl,swgt(j,:nwave),n,x1,y1,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to interpolate for TUV_DNA74.')
       endif
    end if
        
    ! The units in the files for the spectral weights are (mW/m2)-1, but the spectral
    ! irradiance is in W/m2/nm. Thus we need to scale the weights by 10^3 so that the
    ! resulting integral will be unitlessthe 
    swgt(j,:) = swgt(j,:) * 1000._r8
    
    call mpibcast(swgt(j,:), nwave, mpir8, 0, mpicom)


    ! Damage to lens of pig eyes, from 
    ! Oriowo, M. et al. (2001). Action spectrum for in vitro
    ! UV-induced cataract using whole lenses. Invest. Ophthalmol. & Vis. Sci. 42,
    ! 2596-2602.  For pig eyes. Last two columns computed by L.O.Bjorn.
    j = j + 1
    sname(j) = 'TUV_LENS01'
    call addfld(trim(sname(j)), horiz_only, 'A', '1', 'Cataract, pig (Oriowo et al., 2001)', flag_xyfill=.true.)
    
    ! Read from file.
    if (masterproc) then
       call getfil(trim(input_data_root)//'/'//trim('DATAS1/cataract_oriowo'), efile, fid)
       kin = getunit()
       open(unit=kin,file=efile,status='old')
       do i = 1, 7
         read(kin,*)
       end do
       n = 18
       do i = 1, n
          read(kin,*) x1(i), dum1, dum1, y1(i)
       enddo
       close(kin)
       call freeunit(kin)
       ! extrapolation to 400 nm (has very little effect on raf):
!      do i = 1, 30
!         n = n + 1
!         x1(n) = x1(n-1) + 1.
!         y1(n) = 10**(5.7666 - 0.0254*x1(n))
!      enddo
       ! Interpolate onto wavelength grid.
       call addpnt(x1,y1,kdata,n,x1(1)*(1._r8-deltax),y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,0._r8,y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,x1(n)*(1._r8+deltax),0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,1.e+38_r8,0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call inter2(nwave,wl,swgt(j,:nwave),n,x1,y1,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to interpolate for TUV_DNA74.')
       endif
    end if
    
    call mpibcast(swgt(j,:), nwave, mpir8, 0, mpicom)


    ! Vitamin D - CIE 2006
    ! Action spectrum for the production fo previtamin-D3 in human skin, 
    ! CIE Techincal Report TC 6-54, Commission Internatinale del'Eclairage, 2006.
    ! Wavelength range of data is 252-330 nm, but Values below 260 nm and beyond 
    ! 315 nm were interpolated by CIE using a spline fit.
    ! TUV also assigns the 252nm value to shorter wavelengths, and zero 
    ! beyond 330nm.
    j = j + 1
    sname(j) = 'TUV_VITDCIE06'
    call addfld(trim(sname(j)), horiz_only, 'A', '1', 'Previtamin-D3 (CIE 2006)', flag_xyfill=.true.)
    
    ! Read from file.
    if (masterproc) then
       call getfil(trim(input_data_root)//'/'//trim('DATAS1/vitamin_D.txt'), efile, fid)
       kin = getunit()
       open(unit=kin,file=efile,status='old')
       do i = 1, 7
         read(kin,*)
       end do
       n = 79
       do i = 1, n
          read(kin,*) x1(i), y1(i)
       enddo
       close(kin)
       call freeunit(kin)
       ! extrapolation to 400 nm (has very little effect on raf):
       ! Interpolate onto wavelength grid.
       call addpnt(x1,y1,kdata,n,x1(1)*(1._r8-deltax),y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,0._r8,y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,x1(n)*(1._r8+deltax),0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,1.e+38_r8,0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call inter2(nwave,wl,swgt(j,:nwave),n,x1,y1,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to interpolate for TUV_DNA74.')
       endif
    end if
    
    call mpibcast(swgt(j,:), nwave, mpir8, 0, mpicom)


    ! Non-melanoma skin cancer, CIE 2006.
    ! Action spectrum for the induction of non-melanoma skin cancer. From:
    ! Photocarcinogenesis Action Spectrum (Non-Melanoma Skin Cancers), 
    ! CIE S 019/E:2006, Commission Internationale de l'Eclairage, 2006.
    ! 1 nm spacing from 250 to 400 nm. Normalized at maximum, 299 nm.
    ! Set constanta at 3.94E-04 between 340 and 400 nm.
    ! Assume zero beyond 400 nm.
    ! Assume constant below 250 nm.
    j = j + 1
    sname(j) = 'TUV_NMSCCIE06'
    call addfld(trim(sname(j)), horiz_only, 'A', '1', 'NMSC (CIE 2006)', flag_xyfill=.true.)
    
    ! Read from file.
    if (masterproc) then
       call getfil(trim(input_data_root)//'/'//trim('DATAS1/nmsc_cie.txt'), efile, fid)
       kin = getunit()
       open(unit=kin,file=efile,status='old')
       do i = 1, 7
         read(kin,*)
       end do
       n = 151
       do i = 1, n
          read(kin,*) x1(i), y1(i)
       enddo
       close(kin)
       call freeunit(kin)
       ! extrapolation to 400 nm (has very little effect on raf):
       ! Interpolate onto wavelength grid.
       call addpnt(x1,y1,kdata,n,x1(1)*(1._r8-deltax),y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,0._r8,y1(1),ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,x1(n)*(1._r8+deltax),0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call addpnt(x1,y1,kdata,n,1.e+38_r8,0._r8,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to addpnt for TUV_DNA74.')
       endif
       call inter2(nwave,wl,swgt(j,:nwave),n,x1,y1,ierr)
       if (ierr .ne. 0) then
          call endrun('radxfr_cam_init - ERROR: unable to interpolate for TUV_DNA74.')
       endif
    end if
    
    call mpibcast(swgt(j,:), nwave, mpir8, 0, mpicom)

   end subroutine radxfr_cam_init

  subroutine radxfr_cam_update(ncol, lchnk, esfact, zenith, albedo, press_mid, alt, temp, o2vmr, o3vmr, so2vmr, no2vmr, novmr, co2vmr, h2ovmr, hbrvmr, hclvmr, cldfrac, cldw, pbuf)
    use physconst,       only : rairv
    use ref_pres,        only : press_top=>ptop_ref
    use physics_buffer,  only : pbuf_get_field, pbuf_get_index, physics_buffer_desc
    use cam_history,     only : outfld
    use wavelength_grid, only : nwave, wl
    use cam_abortutils,  only : endrun
    
    integer,  intent(in) :: ncol, lchnk
    real(r8), intent(in) :: esfact
    real(r8), intent(in) :: zenith(:)
    real(r8), intent(in) :: albedo(:)
    real(r8), intent(in) :: press_mid(:,:)
    real(r8), intent(in) :: alt(:,:) !  kilometers
    real(r8), intent(in) :: temp(:,:)
    real(r8), intent(in) :: o2vmr(:,:)
    real(r8), intent(in) :: o3vmr(:,:)
    real(r8), intent(in) :: so2vmr(:,:)
    real(r8), intent(in) :: no2vmr(:,:)
    real(r8), intent(in) :: novmr(:,:)
    real(r8), intent(in) :: co2vmr(:,:)
    real(r8), intent(in) :: h2ovmr(:,:)
    real(r8), intent(in) :: hbrvmr(:,:)
    real(r8), intent(in) :: hclvmr(:,:)
    real(r8), intent(in) :: cldfrac(:,:)
    real(r8), intent(in) :: cldw(:,:) ! kg/kg
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer :: i, j, k, iwv
    integer :: errflg
    character(len=512) :: errmsg

    real(r8) :: alt_meters(pver)
    real(r8) :: watdens(pver) ! g/m3 <-- cldw (kg/kg)
    real(r8) :: cloudfr(pver) ! 

    real(r8) :: dto2(pver,nwave)
    real(r8) :: srb_o2_xs(nwave,pver)
    
    real(r8), pointer, dimension(:,:,:) :: swaertau   ! shortwave aerosol tau
    real(r8), pointer, dimension(:,:,:) :: swaertauw  ! shortwave aerosol tau * w
    real(r8), pointer, dimension(:,:,:) :: swaertauwg ! shortwave aerosol tau * w * g

    real(r8) :: swaerw(pcols, pver, nswbands)
    real(r8) :: swaerg(pcols, pver, nswbands)

    real(r8) :: tauaer(pcols, pver, nwave)
    real(r8) :: waer(pcols, pver, nwave)
    real(r8) :: gaer(pcols, pver, nwave)
    
    real(r8) :: spectral_irradiance(pver, nwave, pcols)
    real(r8) :: sintegral(nintegrals, pcols)
    real(r8) :: tmp(pcols, pver)
    
    integer  :: ierr, tmp_idx
    real(r8), pointer, dimension(:) :: tmp_field
    
    if (.not.do_radxfr) return

    errflg=0
    errmsg=' '
    
    ! Get the aerosol optical properties.
    if (has_aer_ra_feedback) then
      call pbuf_get_field(pbuf, swaertau_idx,   swaertau)
      call pbuf_get_field(pbuf, swaertauw_idx,  swaertauw)
      call pbuf_get_field(pbuf, swaertauwg_idx, swaertauwg)
    
      ! Need to convert tau*w to w and tau*w*g to g for the radiation code.
      where(swaertau .ne. 0._r8)
        swaerw = swaertauw / swaertau
      elsewhere
        swaerw = 1._r8
      end where

      where(swaertauw .ne. 0._r8)
        swaerg = swaertauwg / swaertauw
      elsewhere
        swaerg = 0._r8
      end where
    
      ! The CESM wavelengths to the wavelength grid used by TUV.
      do i = 1, ncol
        do k = 1, pver
          call lininterp(swaertau(i,k,1:nswbands-1), nswbands-1, tauaer(i,k,:), nwave, interp_wgts)
          call lininterp(swaerw(i,k,1:nswbands-1), nswbands-1, waer(i,k,:), nwave, interp_wgts)
          call lininterp(swaerg(i,k,1:nswbands-1), nswbands-1, gaer(i,k,:), nwave, interp_wgts)
        end do
      end do
    else
      tauaer(:,:,:) = 0._r8
      waer(:,:,:) = 1._r8
      gaer(:,:,:)   = 0._r8
    end if
    
    ! For DEBUG, output the aerosol optical properties for the tuv bands.
    call outfld('TUV_TAULOW', tauaer(:ncol,:,1), ncol, lchnk)
    call outfld('TUV_WLOW',   waer(:ncol,:,1), ncol, lchnk)
    call outfld('TUV_GLOW',   gaer(:ncol,:,1), ncol, lchnk)
    call outfld('TUV_TAU600', tauaer(:ncol,:,ituv600), ncol, lchnk)
    call outfld('TUV_W600',   waer(:ncol,:,ituv600), ncol, lchnk)
    call outfld('TUV_G600',   gaer(:ncol,:,ituv600), ncol, lchnk)
    call outfld('TUV_TAUHIGH', tauaer(:ncol,:,nwave), ncol, lchnk)
    call outfld('TUV_WHIGH',   waer(:ncol,:,nwave), ncol, lchnk)
    call outfld('TUV_GHIGH',   gaer(:ncol,:,nwave), ncol, lchnk)    
    
    do i = 1,ncol    

       alt_meters(:) = alt(i,:)*1.e3_r8 ! km --> m

       call molec_ox_xsect_run( pver, zenith(i), alt_meters, temp(i,:), press_mid(i,:), press_top, o2vmr(i,:), dto2, srb_o2_xs, errmsg, errflg )
       if (errflg/=0) then
          call endrun('radxfr_cam_update: '//trim(errmsg))
       end if

       !            1.e3_r8 * kg/kg * Pa / ((J/K/kg) * K) --> g/m3
       watdens(:) = 1.e3_r8 * cldw(i,:)*press_mid(i,:)/(rairv(i,:,lchnk)*temp(i,:)) 
       cloudfr(:) = cldfrac(i,:)

       call tuv_radiation_transfer_run( pver, nwave, &
            zenith(i), albedo(i), press_mid(i,:), press_top, alt_meters, temp(i,:), &
            o3vmr(i,:), so2vmr(i,:), no2vmr(i,:), novmr(i,:), co2vmr(i,:), h2ovmr(i,:), hbrvmr(i,:), hclvmr(i,:), cloudfr, watdens, dto2, &
            has_aer_ra_feedback, tauaer(i,:,:), waer(i,:,:), gaer(i,:,:), &
            actinic_fluxes(:,:,i,lchnk), spectral_irradiance(:,:,i), errmsg, errflg )
       if (errflg/=0) then
          call endrun('radxfr_cam_update: '//trim(errmsg))
       end if
       
               
      ! Check for small negative values. Actinic flux should be positive.
      do k = 1, pver
        do iwv = 1, nwave
          if (actinic_fluxes(iwv, k, i, lchnk) .lt. 0._r8) then
!            write(iulog, *) "radxfr: negative actinic flux ... reseting to 0 for ", iwv, k, i, lchnk
            actinic_fluxes(iwv, k, i, lchnk) = 0._r8
          end if

          if (spectral_irradiance(k, iwv, i) .lt. 0._r8) then
!            write(iulog, *) "radxfr: negative spectral irradiance ... reseting to 0 for ", iwv, k, i, lchnk
            spectral_irradiance(k, iwv, i) = 0._r8
          end if          
        end do
      end do
    end do
    
    ! Only need this at the surface for now.
    !
    ! NOTE: The actinic flux and the spectral irradiance are defined with k=1 as the bottom,
    ! not the top as is normal in CAM.
    sintegral(:,:) = 0._r8
    do j = 1, nintegrals
      do iwv = 1, nwave-1
        sintegral(j, :ncol) = sintegral(j, :ncol) + &
          esfact * swgt(j, iwv) * spectral_irradiance(1, iwv, :ncol) * abs(wl(iwv+1) - wl(iwv))
      end do
      
      ! Output the diagnostic.
      call outfld(trim(sname(j)), sintegral(j, :ncol), ncol, lchnk)
      
      ! If required, output the data to a pbuf field
      tmp_idx = -1
      tmp_idx = pbuf_get_index(sname(j), ierr)
      if (tmp_idx > 0) then
        call pbuf_get_field(pbuf, tmp_idx, tmp_field)
        tmp_field(:ncol) = sintegral(j, :ncol)
      end if 
    end do
    
    ! output actinic flux at certain wavelengths import to stratospheric chemistry
    !
    ! NOTE: Need to flip the fluxes for output and also reverse in the vertical.
    do j = 1, nfluxes
      do i = 1, ncol
        do k = 1, pver
          tmp(i,k) = esfact * actinic_fluxes(fwv(j),pver-k+1,i,lchnk)
        end do
      end do
      call outfld(trim(fname(j)), tmp(:ncol, :), ncol, lchnk)
    end do   

  end subroutine radxfr_cam_update
  
  end module radxfr_cam
