module ocn_import_export

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_FieldMod
   use POP_GridHorzMod
   use POP_HaloMod
   use POP_IOUnitsMod
   use POP_MCT_vars_mod
   use POP_CplIndices

   use seq_flds_mod
   use seq_timemgr_mod
   use shr_file_mod 
   use shr_cal_mod,       only : shr_cal_date2ymd
   use shr_sys_mod

   use perf_mod
   use kinds_mod,         only: int_kind, r8
   use communicate,       only: my_task, master_task
   use constants
   use blocks
   use domain,            only: distrb_clinic, POP_haloClinic
   use exit_mod
   use forcing_shf,       only: SHF_QSW
   use forcing_sfwf,      only: lsend_precip_fact, precip_fact
   use forcing_fields,    only: EVAP_F, PREC_F, SNOW_F, MELT_F, ROFF_F, IOFF_F
   use forcing_fields,    only: SALT_F
   use forcing_fields,    only: SENH_F, LWUP_F, LWDN_F, MELTH_F
   use forcing_fields,    only: ATM_CO2_PROG_nf_ind, ATM_CO2_DIAG_nf_ind
   use forcing_fields,    only: ATM_NHx_nf_ind, ATM_NOy_nf_ind
   use forcing_fields,    only: IFRAC, U10_SQR, ATM_PRESS
   use forcing_fields,    only: LAMULT, USTOKES, VSTOKES
   use forcing_fields,    only: ATM_FINE_DUST_FLUX, ATM_COARSE_DUST_FLUX, SEAICE_DUST_FLUX
   use forcing_fields,    only: ATM_BLACK_CARBON_FLUX, SEAICE_BLACK_CARBON_FLUX
   use forcing_fields,    only: ATM_PAR,ATM_EPHYTO1,ATM_EPHYTO1_NET, ATM_EPHYTO4,ATM_EPHYTO4_NET,ATM_EPHYTO5, ATM_EPHYTO5_NET
   use mcog,              only: lmcog, mcog_ncols, import_mcog
   use forcing_coupled,   only: ncouple_per_day,  &
                                update_ghost_cells_coupler_fluxes, &
                                rotate_wind_stress, pop_set_coupled_forcing, &
                                pop_init_coupled,  &
                                orb_eccen, orb_obliqr, orb_lambm0, orb_mvelpp
   use ice,               only: tfreez, tmelt, liceform,QFLUX, QICE, AQICE, &
                                tlast_ice
   use global_reductions, only: global_sum_prod
   use io_tools,          only: document
   use named_field_mod,   only: named_field_register, named_field_get_index, &
                                named_field_set, named_field_get
   use prognostic
   use time_management
   use registry
   ! QL, 150526, ocn<->wav
   use vmix_kpp,          only: KPP_HBLT      ! ocn -> wav, bounadry layer depth

   implicit none
   public
   save

   ! accumulated sum of send buffer quantities for averaging before being sent
   real (r8), dimension(:,:,:,:), allocatable ::  SBUFF_SUM 

   real (r8) :: tlast_coupled 

contains

!***********************************************************************
!BOP
! !IROUTINE: ocn_import
! !INTERFACE:

  subroutine ocn_import(x2o, ldiag_cpl, errorCode)

! !DESCRIPTION:
!-----------------------------------------------------------------------
!  This routine receives message from cpl7 driver
!
!    The following fields are always received from the coupler:
! 
!    o  taux   -- zonal wind stress (taux)                 (W/m2   )
!    o  tauy   -- meridonal wind stress (tauy)             (W/m2   )
!    o  snow   -- water flux due to snow                   (kg/m2/s)
!    o  rain   -- water flux due to rain                   (kg/m2/s)
!    o  evap   -- evaporation flux                         (kg/m2/s)
!    o  meltw  -- snow melt flux                           (kg/m2/s)
!    o  salt   -- salt                                     (kg(salt)/m2/s)
!    o  swnet  -- net short-wave heat flux                 (W/m2   )
!    o  sen    -- sensible heat flux                       (W/m2   )
!    o  lwup   -- longwave radiation (up)                  (W/m2   )
!    o  lwdn   -- longwave radiation (down)                (W/m2   )
!    o  melth  -- heat flux from snow&ice melt             (W/m2   )
!    o  ifrac  -- ice fraction
!    o  rofl   -- river runoff flux                        (kg/m2/s)
!    o  rofi   -- ice runoff flux                          (kg/m2/s)
! 
!    The following fields are sometimes received from the coupler,
!      depending on model options:
! 
!    o  pslv   -- sea-level pressure                       (Pa)
!    o  duu10n -- 10m wind speed squared                   (m^2/s^2)
!    o  co2prog-- bottom atm level prognostic co2
!    o  co2diag-- bottom atm level diagnostic co2
! 
!-----------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

    real(r8)           , intent(inout) :: x2o(:,:)
    logical (log_kind) , intent(in)    :: ldiag_cpl
    integer (POP_i4)   , intent(out)   :: errorCode  ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (char_len) ::   &
      label,                 &
      message
 
   integer (int_kind) ::  &
      i,j,k,n,ncol,iblock

   real (r8), dimension(nx_block,ny_block) ::  &
      WORKB

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::   &
      WORK1, WORK2        ! local work space

   real (r8), dimension(mcog_ncols) ::  &
      frac_col_1pt, fracr_col_1pt, qsw_fracr_col_1pt

   real (r8) ::  &
      m2percm2,  &
      gsum

   type (block) :: this_block ! local block info

   integer (int_kind) :: nrecv

!-----------------------------------------------------------------------
!
!  zero out padded cells 
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   WORK1 = c0
   WORK2 = c0

!-----------------------------------------------------------------------
!
!  unpack and distribute wind stress, then convert to correct units
!  and rotate components to local coordinates
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)

      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         WORK1(i,j,iblock) = x2o(index_x2o_Foxx_taux,n)
         WORK2(i,j,iblock) = x2o(index_x2o_Foxx_tauy,n)
      enddo
      enddo
   enddo ! iblock

   !***
   !*** do NOT perform halo updates now, because vector updates must
   !***   be done after the rotation is completed.
   !***

!-----------------------------------------------------------------------
!
!  rotate true zonal/meridional wind stress into local coordinates,
!  convert to dyne/cm**2, and shift SMFT to U grid
!
!  halo updates are performed in subroutine rotate_wind_stress, 
!  following the rotation
!
!-----------------------------------------------------------------------

      call rotate_wind_stress(WORK1, WORK2)

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)

!-----------------------------------------------------------------------
!
!  unpack and distribute fresh water flux and salt flux
!
!  NOTE: if there are code changes associated with changing the names or
!        the number of fluxes received from the coupler, then subroutine
!        update_ghost_cells_coupler_fluxes will need to be modified also
!
!-----------------------------------------------------------------------


      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         SNOW_F(i,j,iblock) = x2o(index_x2o_Faxa_snow,n)
         WORKB (i,j       ) = x2o(index_x2o_Faxa_rain,n)
         EVAP_F(i,j,iblock) = x2o(index_x2o_Foxx_evap,n)
         MELT_F(i,j,iblock) = x2o(index_x2o_Fioi_meltw,n)
         ROFF_F(i,j,iblock) = x2o(index_x2o_Foxx_rofl,n)
         IOFF_F(i,j,iblock) = x2o(index_x2o_Foxx_rofi,n)
         SALT_F(i,j,iblock) = x2o(index_x2o_Fioi_salt,n)

         PREC_F(i,j,iblock) = WORKB(i,j) + SNOW_F(i,j,iblock)    ! rain + snow

         WORKB(i,j        ) = x2o(index_x2o_Foxx_swnet,n)
         SHF_QSW(i,j,iblock) = WORKB(i,j)*  &
            RCALCT(i,j,iblock)*hflux_factor  !  convert from W/m**2
         SENH_F(i,j,iblock)  = x2o(index_x2o_Foxx_sen,n)
         LWUP_F(i,j,iblock)  = x2o(index_x2o_Foxx_lwup,n)
         LWDN_F(i,j,iblock)  = x2o(index_x2o_Faxa_lwdn,n)
         MELTH_F(i,j,iblock) = x2o(index_x2o_Fioi_melth,n)

         WORKB(i,j       ) = x2o(index_x2o_Si_ifrac,n)
         IFRAC(i,j,iblock) = WORKB(i,j) * RCALCT(i,j,iblock)

         !***  converting from Pa to dynes/cm**2
         WORKB(i,j       ) = x2o(index_x2o_Sa_pslv,n)
         ATM_PRESS(i,j,iblock) = c10 * WORKB(i,j) * RCALCT(i,j,iblock)

         !***  converting from m**2/s**2 to cm**2/s**2
         WORKB(i,j       ) = x2o(index_x2o_So_duu10n,n)
         U10_SQR(i,j,iblock) = cmperm * cmperm * WORKB(i,j) * RCALCT(i,j,iblock)
         
         ! QL, 150526, langmuir mixing variables
         ! QL, 150908, only apply enhancement factor to ice-free region
         if (IFRAC(i,j,iblock) <= 0.05_r8) then
            ! import enhancement factor (unitless)
            WORKB(i,j       ) = x2o(index_x2o_Sw_lamult,n)
            LAMULT(i,j,iblock) = WORKB(i,j)*RCALCT(i,j,iblock)
         else
            LAMULT(i,j,iblock) = c1
         endif
         ! QL, 150706, DEBUG
         !print*, "LAMULT = ", LAMULT(i,j,iblock)  
         !print*, "RCALCT = ", RCALCT(i,j,iblock)  
         !print*, "WORKB = ", WORKB(i,j)  
         !! DEBUG
         ! import surface Stokes drift (m/s)
         WORKB(i,j       ) = x2o(index_x2o_Sw_ustokes,n)
         USTOKES(i,j,iblock) = WORKB(i,j)*RCALCT(i,j,iblock)
         WORKB(i,j       ) = x2o(index_x2o_Sw_vstokes,n)
         VSTOKES(i,j,iblock) = WORKB(i,j)*RCALCT(i,j,iblock)

         ! convert dust flux from MKS (kg/m^2/s) to CGS (g/cm^2/s)
         ATM_FINE_DUST_FLUX(i,j,iblock) = 0.1_r8 * RCALCT(i,j,iblock) * ( &
            x2o(index_x2o_Faxa_dstwet1,n) + x2o(index_x2o_Faxa_dstdry1,n))

         ! convert dust flux from MKS (kg/m^2/s) to CGS (g/cm^2/s)
         ATM_COARSE_DUST_FLUX(i,j,iblock) = 0.1_r8 * RCALCT(i,j,iblock) * ( &
            x2o(index_x2o_Faxa_dstwet2,n) + x2o(index_x2o_Faxa_dstwet3,n) + x2o(index_x2o_Faxa_dstwet4,n) + &
            x2o(index_x2o_Faxa_dstdry2,n) + x2o(index_x2o_Faxa_dstdry3,n) + x2o(index_x2o_Faxa_dstdry4,n))

         ! convert dust flux from MKS (kg/m^2/s) to CGS (g/cm^2/s)
         SEAICE_DUST_FLUX(i,j,iblock) = 0.1_r8 * RCALCT(i,j,iblock) * ( &
            x2o(index_x2o_Fioi_flxdst,n))

         ! convert black carbon flux from MKS (kg/m^2/s) to CGS (g/cm^2/s)
         ATM_BLACK_CARBON_FLUX(i,j,iblock) = 0.1_r8 * RCALCT(i,j,iblock) * ( &
            x2o(index_x2o_Faxa_bcphidry,n) + x2o(index_x2o_Faxa_bcphodry,n) + &
            x2o(index_x2o_Faxa_bcphiwet,n))

         ! convert black carbon flux from MKS (kg/m^2/s) to CGS (g/cm^2/s)
         SEAICE_BLACK_CARBON_FLUX(i,j,iblock) = 0.1_r8 * RCALCT(i,j,iblock) * ( &
            x2o(index_x2o_Fioi_bcpho,n) + x2o(index_x2o_Fioi_bcphi,n))
         
         !  (Optional) Surface Radiation fields from ATM via TUV   
         if (index_x2o_Faxa_par > 0) then
            ATM_PAR(i,j,iblock) = x2o(index_x2o_Faxa_par,n)
         end if
         !if (index_x2o_Faxa_par_fake > 0) then
         !   ATM_PAR_FAKE(i,j,iblock) = x2o(index_x2o_Faxa_par_fake,n) ! want to convert this to swnet with ATM_PAR_FAKE
         !end if
         if (index_x2o_Sa_ephyto1 > 0) then
            ATM_EPHYTO1(i,j,iblock) = x2o(index_x2o_Sa_ephyto1,n) ! original - does not consider albedo
         end if
         if (index_x2o_Sa_ephyto1_net > 0) then
            ATM_EPHYTO1_NET(i,j,iblock) = x2o(index_x2o_Sa_ephyto1_net,n) ! replaced ephyto1 with ephyto1_net
         end if
         !if (index_x2o_Sa_ephyto2 > 0) then
         !   ATM_EPHYTO2(i,j,iblock) = x2o(index_x2o_Sa_ephyto2,n)
         !end if
        ! if (index_x2o_Sa_ephyto3 > 0) then
        !    ATM_EPHYTO3(i,j,iblock) = x2o(index_x2o_Sa_ephyto3,n)
        ! end if
         if (index_x2o_Sa_ephyto4 > 0) then
            ATM_EPHYTO4(i,j,iblock) = x2o(index_x2o_Sa_ephyto4,n)
         end if
         if (index_x2o_Sa_ephyto4_net > 0) then
            ATM_EPHYTO4_NET(i,j,iblock) = x2o(index_x2o_Sa_ephyto4_net,n)
         end if
         if (index_x2o_Sa_ephyto5 > 0) then
            ATM_EPHYTO5(i,j,iblock) = x2o(index_x2o_Sa_ephyto5,n)
         end if
         if (index_x2o_Sa_ephyto5_net > 0) then
            ATM_EPHYTO5_NET(i,j,iblock) = x2o(index_x2o_Sa_ephyto5_net,n)
         end if
      enddo
      enddo

   enddo

!-----------------------------------------------------------------------
!
!  optional fields per mcog column
!
!-----------------------------------------------------------------------

   if (lmcog) then

      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1

            ! extract fields for each column and pass to import_mcog

            do ncol = 1, mcog_ncols
               frac_col_1pt(ncol)  = max(c0, min(c1, x2o(index_x2o_frac_col(ncol), n)))
               fracr_col_1pt(ncol) = max(c0, min(c1, x2o(index_x2o_fracr_col(ncol), n)))
               qsw_fracr_col_1pt(ncol) = x2o(index_x2o_qsw_fracr_col(ncol), n)
            enddo

            call import_mcog(frac_col_1pt, fracr_col_1pt, qsw_fracr_col_1pt, &
                             x2o(index_x2o_Foxx_swnet,n), iblock, i, j)

         enddo ! do j
         enddo ! do i
      enddo ! do iblock = 1, nblocks_clinic

   else ! if (lmcog) then

      ! if mcog is off, fill its arrays with cpl aggregated full cell means

      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1

            ncol = 1
            frac_col_1pt(ncol)  = c1
            fracr_col_1pt(ncol) = c1
            qsw_fracr_col_1pt(ncol) = x2o(index_x2o_Foxx_swnet,n)

            call import_mcog(frac_col_1pt, fracr_col_1pt, qsw_fracr_col_1pt, &
                             x2o(index_x2o_Foxx_swnet,n), iblock, i, j)

         enddo ! do j
         enddo ! do i
      enddo ! do iblock = 1, nblocks_clinic

   endif ! if (lmcog) then

!-----------------------------------------------------------------------
!
!  update ghost cells for fluxes received from the coupler
!
!-----------------------------------------------------------------------

   call update_ghost_cells_coupler_fluxes(errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'ocn_import: error in update_ghost_cells_coupler_fluxes')
      return
   endif

!-----------------------------------------------------------------------
!
!  unpack atmospheric CO2
!
!-----------------------------------------------------------------------

   if (index_x2o_Sa_co2prog > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Sa_co2prog,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating PROG CO2 halo')
         return
      endif

      call named_field_set(ATM_CO2_PROG_nf_ind, WORK1)
   endif

   if (index_x2o_Sa_co2diag > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Sa_co2diag,n)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating DIAG CO2 halo')
         return
      endif

      call named_field_set(ATM_CO2_DIAG_nf_ind, WORK1)
   endif
 
   if (index_x2o_Faxa_nhx > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         ! Note - the input units are kgN/m2/s to nmolN/cm2/s
         ! TODO: Keith has pointed out might want to use 14.007_r8 instead of 14.0_r8 for more
         ! consistency when bringing in N isotopes into the code
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Faxa_nhx,n) * (1.0e-1_r8 * (c1/14.0_r8) * 1.0e9_r8) 
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating DIAG NHx halo')
         return
      endif

      call named_field_set(ATM_NHx_nf_ind, WORK1)
   endif

   if (index_x2o_Faxa_noy > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         ! Note - the input units are kgN/m2/s to nmolN/cm2/s
         ! TODO: Keith has pointed out might want to use 14.007_r8 instead of 14.0_r8 for more
         ! consistency when bringing in N isotopes into the code

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            WORK1(i,j,iblock) = x2o(index_x2o_Faxa_noy,n) * (1.0e-1_r8 * (c1/14.0_r8) * 1.0e9_r8)
         enddo
         enddo
      enddo

      call POP_HaloUpdate(WORK1,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'ocn_import_mct: error updating DIAG NOy halo')
         return
      endif

      call named_field_set(ATM_NOy_nf_ind, WORK1)
   endif

!-----------------------------------------------------------------------
!
!  diagnostics
!
!-----------------------------------------------------------------------

   if (ldiag_cpl) then

     write(message,'(6a,1x,5a)')  &
         ' Global averages of fluxes received from cpl at ',  &
           cyear,'/',cmonth ,'/',cday,  chour,':',cminute,':',csecond
     call document ('pop_recv_from_coupler', trim(message))
 
     m2percm2  = mpercm*mpercm
     nrecv = size(x2o, dim=1)
     do k = 1,nrecv

         n = 0
         do iblock = 1, nblocks_clinic
            this_block = get_block(blocks_clinic(iblock),iblock)

            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
               n = n + 1
               WORK1(i,j,iblock) = x2o(k,n)  ! mult. by TAREA in global_sum_prod
            enddo
            enddo
         enddo

         gsum = global_sum_prod(WORK1 , TAREA, distrb_clinic, &
                                 field_loc_center, RCALCT)*m2percm2
         if (my_task == master_task) then
            call seq_flds_getField(label,k,seq_flds_x2o_fields)
            write(stdout,1100)'ocn','recv', label ,gsum
            call shr_sys_flush(stdout)
         endif
      enddo
   endif


1100  format ('comm_diag ', a3, 1x, a4, 1x, a8, 1x, es26.19:, 1x, a6)

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_import

!***********************************************************************
!BOP
! !IROUTINE: ocn_export_mct
! !INTERFACE:

 subroutine ocn_export(o2x, ldiag_cpl, errorCode)   

! !DESCRIPTION:
!  This routine calls the routines necessary to send pop fields to
!  the CCSM cpl7 driver
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

   real(r8)           , intent(inout) :: o2x(:,:)
   logical (log_kind) , intent(in)    :: ldiag_cpl
   integer (POP_i4)   , intent(out)   :: errorCode  ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: n, iblock
           
   character (char_len)    :: label
 
   integer (int_kind) ::  &
        i,j,k

   real (r8), dimension(nx_block,ny_block) ::   &
        WORK1, WORK2,      &! local work space
        WORK3, WORK4

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::   &
        WORKA               ! local work space with full block dimension

   real (r8) ::   &
      m2percm2,   &
      gsum

   type (block) :: this_block ! local block info

   integer (int_kind) :: nsend

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  initialize control buffer
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!
!     interpolate onto T-grid points and rotate on T grid
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
     this_block = get_block(blocks_clinic(iblock),iblock)

     call ugrid_to_tgrid(WORK3,SBUFF_SUM(:,:,iblock,index_o2x_So_u),iblock)
     call ugrid_to_tgrid(WORK4,SBUFF_SUM(:,:,iblock,index_o2x_So_v),iblock)

     WORK1 = (WORK3*cos(ANGLET(:,:,iblock))+WORK4*sin(-ANGLET(:,:,iblock)))  &
            * mpercm/tlast_coupled
     WORK2 = (WORK4*cos(ANGLET(:,:,iblock))-WORK3*sin(-ANGLET(:,:,iblock)))  &
            * mpercm/tlast_coupled

     do j=this_block%jb,this_block%je
     do i=this_block%ib,this_block%ie
        n = n + 1
        o2x(index_o2x_So_u,n) = WORK1(i,j)
        o2x(index_o2x_So_v,n) = WORK2(i,j)
     enddo
     enddo
  enddo

!-----------------------------------------------------------------------
!
!     convert and pack surface temperature
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         o2x(index_o2x_So_t,n) =   &
             SBUFF_SUM(i,j,iblock,index_o2x_So_t)/tlast_coupled + T0_Kelvin
      enddo
      enddo
   enddo

!-----------------------------------------------------------------------
!
!     convert and pack salinity
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         o2x(index_o2x_So_s,n) =   &
             SBUFF_SUM(i,j,iblock,index_o2x_So_s)*salt_to_ppt/tlast_coupled
      enddo
      enddo
   enddo

!-----------------------------------------------------------------------
!
!    QL, 150526, convert and pack boundary layer depth
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         o2x(index_o2x_So_bldepth,n) =   &
             SBUFF_SUM(i,j,iblock,index_o2x_So_bldepth)/100./tlast_coupled
      enddo
      enddo
   enddo

!-----------------------------------------------------------------------
!
!     interpolate onto T-grid points, then rotate on T grid
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      call ugrid_to_tgrid(WORK3,SBUFF_SUM(:,:,iblock,index_o2x_So_dhdx),iblock)
      call ugrid_to_tgrid(WORK4,SBUFF_SUM(:,:,iblock,index_o2x_So_dhdy),iblock)
 
      WORK1 = (WORK3*cos(ANGLET(:,:,iblock)) + WORK4*sin(-ANGLET(:,:,iblock)))  &
              /grav/tlast_coupled
      WORK2 = (WORK4*cos(ANGLET(:,:,iblock)) - WORK3*sin(-ANGLET(:,:,iblock)))  &
              /grav/tlast_coupled

      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         o2x(index_o2x_So_dhdx,n) = WORK1(i,j)
         o2x(index_o2x_So_dhdy,n) = WORK2(i,j)
      enddo
      enddo
   enddo

!-----------------------------------------------------------------------
!
!     pack heat flux due to freezing/melting (W/m^2)
!     QFLUX computation and units conversion occurs in ice.F
!
!-----------------------------------------------------------------------

   n = 0
   do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         n = n + 1
         o2x(index_o2x_Fioo_q,n) = QFLUX(i,j,iblock)
      enddo
      enddo
   enddo

   tlast_ice = c0
   AQICE     = c0
   QICE      = c0

!-----------------------------------------------------------------------
!
!     pack co2 flux, if requested (kg CO2/m^2/s)
!     units conversion occurs where co2 flux is computed
!
!-----------------------------------------------------------------------

   if (index_o2x_Faoo_fco2_ocn > 0) then
      n = 0
      do iblock = 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            n = n + 1
            o2x(index_o2x_Faoo_fco2_ocn,n) = &
               SBUFF_SUM(i,j,iblock,index_o2x_Faoo_fco2_ocn)/tlast_coupled
         enddo
         enddo
      enddo
   endif

!-----------------------------------------------------------------------
!
!     diagnostics
!
!-----------------------------------------------------------------------

   if (ldiag_cpl) then
      call ccsm_char_date_and_time
      !DEBUG      write(message,'(6a,1x,5a)')' Global averages of fluxes sent to cpl at ', &
      !DEBUG           cyear,'/',cmonth, '/',cday,  chour,':',cminute,':',csecond
      !DEBUG      call document ('pop_send_to_coupler', message)
      write(stdout,*)'pop_send_to_coupler'

      m2percm2  = mpercm*mpercm
      nsend = size(o2x,dim=1)
      do k = 1,nsend
        n = 0
        do iblock = 1, nblocks_clinic
           this_block = get_block(blocks_clinic(iblock),iblock)
           do j=this_block%jb,this_block%je
           do i=this_block%ib,this_block%ie
              n = n + 1
              WORKA(i,j,iblock) = o2x(k,n)
           enddo
           enddo
        enddo

        call POP_HaloUpdate(WORKA,POP_haloClinic,          &
                       POP_gridHorzLocCenter,          &
                       POP_fieldKindScalar, errorCode, &
                       fillValue = 0.0_POP_r8)
       
         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'ocn_export_mct: error updating halo for state')
            return
         endif

        gsum = global_sum_prod(WORKA , TAREA, distrb_clinic, &
                                   field_loc_center, RCALCT)*m2percm2
        if (my_task == master_task) then
           call seq_flds_getField(label,k,seq_flds_o2x_fields)
           write(stdout,1100)'ocn','send', label ,gsum
        endif
      enddo ! k
      if (my_task == master_task) call shr_sys_flush(stdout)
   endif

1100 format ('comm_diag ', a3, 1x, a4, 1x, a8, 1x, es26.19:, 1x, a6)
    
    tlast_coupled = c0

!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_export

!***********************************************************************

!BOP
! !IROUTINE: POP_sum_buffer
! !INTERFACE:

 subroutine POP_sum_buffer

! !DESCRIPTION:
!  This routine accumulates sums for averaging fields to
!  be sent to the coupler
!
! !REVISION HISTORY:
!  same as module
! 
!EOP
!BOC

#ifdef CCSMCOUPLED
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::  &
      WORK                ! local work arrays

   real (r8) ::   &
      delt,             & ! time interval since last step
      delt_last           ! time interval for previous step

   integer (int_kind) :: &
      iblock,           & ! block index
      sflux_co2_nf_ind = 0! named field index of fco2

   logical (log_kind) :: &
      first = .true.      ! only true for first call

   save first

!-----------------------------------------------------------------------
!
!  zero buffer if this is the first time after a coupling interval
!
!-----------------------------------------------------------------------

   if (tlast_coupled == c0) SBUFF_SUM = c0
   WORK = c0

!-----------------------------------------------------------------------
!
!  update time since last coupling
!
!-----------------------------------------------------------------------

   if (avg_ts .or. back_to_back) then
      delt = p5*dtt
   else
      delt =    dtt
   endif
   tlast_coupled = tlast_coupled + delt

!-----------------------------------------------------------------------
!
!  allow for fco2 field to not be registered on first call
!     because init_forcing is called before init_passive_tracers
!  use weight from previous timestep because flux used here is that
!     computed during the previous timestep
!
!-----------------------------------------------------------------------

   if (index_o2x_Faoo_fco2_ocn > 0) then
      if (sflux_co2_nf_ind == 0) then
         call named_field_get_index('SFLUX_CO2', sflux_co2_nf_ind, &
                                    exit_on_err=.not. first)
      endif

      if (avg_ts .or. back_to_back) then
         delt_last = p5*dtt
      else
         delt_last =    dtt
      endif
   endif

!-----------------------------------------------------------------------
!
!  accumulate sums of U,V,T,S and GRADP
!  accumulate sum of co2 flux, if requested
!     implicitly use zero flux if fco2 field not registered yet
!  ice formation flux is handled separately in ice routine
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1, nblocks_clinic
   SBUFF_SUM(:,:,iblock,index_o2x_So_u) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_u) + delt*  &
                                   UVEL(:,:,1,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2x_So_v) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_v) + delt*  &
                                   VVEL(:,:,1,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2x_So_t ) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_t ) + delt*  &
                                   TRACER(:,:,1,1,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2x_So_s ) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_s ) + delt*  &
                                   TRACER(:,:,1,2,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2x_So_dhdx) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_dhdx) + delt*  &
                                   GRADPX(:,:,curtime,iblock)

   SBUFF_SUM(:,:,iblock,index_o2x_So_dhdy) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_dhdy) + delt*  &
                                   GRADPY(:,:,curtime,iblock)

   ! QL, 150526, boundary layer depth
   SBUFF_SUM(:,:,iblock,index_o2x_So_bldepth) =   &
      SBUFF_SUM(:,:,iblock,index_o2x_So_bldepth) + delt*  &
                                   KPP_HBLT(:,:,iblock)

   if (index_o2x_Faoo_fco2_ocn > 0 .and. sflux_co2_nf_ind > 0) then
      call named_field_get(sflux_co2_nf_ind, iblock, WORK(:,:,iblock))
      SBUFF_SUM(:,:,iblock,index_o2x_Faoo_fco2_ocn) = &
         SBUFF_SUM(:,:,iblock,index_o2x_Faoo_fco2_ocn) + delt_last*WORK(:,:,iblock)
   endif

   enddo
   !$OMP END PARALLEL DO

   first = .false.

#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_sum_buffer
 
!***********************************************************************

end module ocn_import_export
