module tuv_radiation_transfer
  use phot_kind_mod, only: rk => kind_phot
  use phot_kind_mod, only: kind_phys => kind_phot

  implicit none

  private
  public :: tuv_radiation_transfer_init
  public :: tuv_radiation_transfer_run
  
contains
  
!> \section arg_table_tuv_radiation_transfer_init Argument Table
!! \htmlinclude tuv_radiation_transfer_init.html
!!
subroutine tuv_radiation_transfer_init( realkind, errmsg, errflg )
    use wavelength_grid,  only: nwave, wl
    use module_xsections, only: rdxs_init

    integer,          intent(in)  :: realkind
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
    
    errmsg = ''
    errflg = 0

    if ( realkind/=kind_phys ) then
       errmsg = 'tuv_radiation_transfer_init: realkind does not match kind_phot'
       errflg = 1
       return
    end if
    
    call rdxs_init( nwave, wl, errmsg, errflg )
    
    if (errflg.ne.0) return
   
  end subroutine tuv_radiation_transfer_init
  
!> \section arg_table_tuv_radiation_transfer_run Argument Table
!! \htmlinclude tuv_radiation_transfer_run.html
!!
  subroutine tuv_radiation_transfer_run( nlev, tuv_n_wavelen, zenith, albedo, press_mid, press_top, alt, temp, &
       o3vmr, so2vmr, no2vmr, novmr, co2vmr, h2ovmr, hbrvmr, hclvmr, cldfrc, cldwat, dto2, has_aer_ra_feedback, tauaerin, waerin, gaerin, radfld, efld, errmsg, errflg )

    use tuv_subs,         only: tuv_radfld
    use wavelength_grid,  only: nwave, wl, wc
    use module_xsections, only: o2_xs, so2_xs, o3xs, no2xs_jpl06a, no_xs, co2_xs, h2o_xs, hbr_xs, hcl_xs
    use params_mod
 
    integer,          intent(in)  :: nlev

    !! NOTE THIS VARIABLE WILL GO AWAY - FOR NOW IS REQUIRED WORKAROUND FOR CPF
    integer,          intent(in)    :: tuv_n_wavelen

    real(kind_phys),  intent(in)  :: zenith
    real(kind_phys),  intent(in)  :: albedo
    real(kind_phys),  intent(in)  :: press_mid(:)
    real(kind_phys),  intent(in)  :: press_top
    real(kind_phys),  intent(in)  :: alt(:)  ! m
    real(kind_phys),  intent(in)  :: temp(:) ! K
    real(kind_phys),  intent(in)  :: o3vmr(:)
    real(kind_phys),  intent(in)  :: so2vmr(:)
    real(kind_phys),  intent(in)  :: no2vmr(:)
    real(kind_phys),  intent(in)  :: novmr(:)
    real(kind_phys),  intent(in)  :: co2vmr(:)
    real(kind_phys),  intent(in)  :: h2ovmr(:)
    real(kind_phys),  intent(in)  :: hbrvmr(:)
    real(kind_phys),  intent(in)  :: hclvmr(:)
    real(kind_phys),  intent(in)  :: cldfrc(:)
    real(kind_phys),  intent(in)  :: cldwat(:) ! cld water content (g/m3)
    real(kind_phys),  intent(in)  :: dto2(:,:)
    logical,          intent(in)  :: has_aer_ra_feedback
    real(kind_phys),  intent(in)  :: tauaerin(:,:)
    real(kind_phys),  intent(in)  :: waerin(:,:)
    real(kind_phys),  intent(in)  :: gaerin(:,:)
    real(kind_phys),  intent(out) :: radfld(:,:) ! /sec
    real(kind_phys),  intent(out) :: efld(:,:)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer, parameter :: nlambda_start=1
    integer, parameter :: cld_od_opt=1
    real(rk), parameter :: dobsi = 0._rk
    
    real(rk) :: zen
    real(rk) :: alb(nwave)
    real(rk) :: zlev(nlev+1) ! km 
    real(rk) :: tlev(nlev)
    real(rk) :: cldfrclev(nlev)
    real(rk) :: cldwatlev(nlev)
    real(rk) :: aircol(nlev)  ! # molecules / cm2 in each layer
    real(rk) :: o3col(nlev) 
    real(rk) :: so2col(nlev)
    real(rk) :: no2col(nlev)
    real(rk) :: nocol(nlev)
    real(rk) :: co2col(nlev)
    real(rk) :: h2ocol(nlev)
    real(rk) :: hbrcol(nlev)
    real(rk) :: hclcol(nlev)
    real(rk) :: dpress(nlev)

    real(rk) :: dtaer(nlev,nwave), omaer(nlev,nwave), gaer(nlev,nwave)
    real(rk) :: dtcld(nlev,nwave), omcld(nlev,nwave), gcld(nlev,nwave)
    real(rk) :: dt_cld(nlev)
    
    real(rk) :: e_dir(nlev,nwave)
    real(rk) :: e_dn(nlev,nwave)
    real(rk) :: e_up(nlev,nwave)
    real(rk) :: dir_fld(nlev,nwave)
    real(rk) :: dwn_fld(nlev,nwave)
    real(rk) :: up_fld(nlev,nwave)

    integer :: k, kk, nlyr

    real(rk) :: o3_xs(nwave,nlev)
    real(rk) :: no2_xs(nwave,nlev)
    real(rk) :: o3_xs_tpose(nlev,nwave)
    real(rk) :: no2_xs_tpose(nlev,nwave)
    real(rk) :: delz_km, delz_cm
    
    errmsg = ''
    errflg = 0

    nlyr = nlev -1
    
    dpress(nlyr:1:-1) = press_mid(2:nlyr+1) - press_mid(1:nlyr)
    do k=1,nlyr
       kk=nlyr-k+1
       aircol(k) = 10._rk*dpress(k)*R/(kboltz*g)
       o3col(k)  = 0.5_rk*(o3vmr(kk)+o3vmr(kk+1))*aircol(k)
       so2col(k) = 0.5_rk*(so2vmr(kk)+so2vmr(kk+1))*aircol(k)
       no2col(k) = 0.5_rk*(no2vmr(kk)+no2vmr(kk+1))*aircol(k)
       nocol(k)  = 0.5_rk*(novmr(kk)+novmr(kk+1))*aircol(k)
       co2col(k) = 0.5_rk*(co2vmr(kk)+co2vmr(kk+1))*aircol(k)
       h2ocol(k) = 0.5_rk*(h2ovmr(kk)+h2ovmr(kk+1))*aircol(k)
       hbrcol(k) = 0.5_rk*(hbrvmr(kk)+hbrvmr(kk+1))*aircol(k)
       hclcol(k) = 0.5_rk*(hclvmr(kk)+hclvmr(kk+1))*aircol(k)
    end do

    ! inputs need to be bottom up vert coord
    tlev(nlev:1:-1) = temp(1:nlev)
    cldfrclev(nlev:1:-1) = cldfrc(1:nlev)
    cldwatlev(nlev:1:-1) = cldwat(1:nlev)
    zlev(nlev:1:-1) = alt(1:nlev)*1.e-3_rk ! m -> km

    delz_km = zlev(nlev) - zlev(nlev-1)
    delz_cm = delz_km*1.e5_rk ! cm

    zlev(nlev+1) = zlev(nlev) + delz_km ! km
    aircol(nlev) = delz_cm * 10._rk * press_top / ( kboltz * tlev(nlev) ) ! molecules / cm2
    o3col(nlev)  = o3vmr(1)  * aircol(nlev)
    so2col(nlev) = so2vmr(1) * aircol(nlev)
    no2col(nlev) = no2vmr(1) * aircol(nlev)
    nocol(nlev)  = novmr(1)  * aircol(nlev)
    co2col(nlev) = co2vmr(1) * aircol(nlev)
    h2ocol(nlev) = h2ovmr(1) * aircol(nlev)
    hbrcol(nlev) = hbrvmr(1) * aircol(nlev)
    hclcol(nlev) = hclvmr(1) * aircol(nlev)

    if (has_aer_ra_feedback) then
      dtaer(nlev:1:-1,:) = tauaerin(1:nlev,:)
      omaer(nlev:1:-1,:) = waerin(1:nlev,:)
      gaer(nlev:1:-1,:)  = gaerin(1:nlev,:)
    else
      dtaer(:,:) = 0.0_rk
      omaer(:,:) = 1.0_rk
      gaer(:,:)  = 0.0_rk
    end if

    zen = zenith
    alb(:) = albedo

    o3_xs_tpose = 0.0_rk
    no2_xs_tpose= 0.0_rk

    call o3xs( nlev,tlev,nwave,wl,o3_xs_tpose )
    call no2xs_jpl06a( nlev,tlev,nwave,wl,no2_xs_tpose )

    o3_xs  = transpose( o3_xs_tpose )
    no2_xs = transpose( no2_xs_tpose )

    call tuv_radfld( nlambda_start, cld_od_opt, cldfrclev, nlev, nwave, &
         zen, zlev, alb, &
         aircol, o3col, so2col, no2col, nocol, co2col, h2ocol, hbrcol, hclcol, &
         dtaer, omaer, gaer, dtcld, omcld, gcld, &
         cldwatlev, dobsi, o3_xs, no2_xs, no_xs, o2_xs, &
         so2_xs, co2_xs, h2o_xs, hbr_xs, hcl_xs, wl(1), wc, tlev, dto2, radfld, efld, &
         e_dir, e_dn, e_up, &
         dir_fld, dwn_fld, up_fld, dt_cld, errmsg, errflg )

  end subroutine tuv_radiation_transfer_run

end module tuv_radiation_transfer
