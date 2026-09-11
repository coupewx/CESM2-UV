! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine calculates new gas concentrations.
!!
!! @author Andy Ackerman, Bill McKie, Chuck Bardeen
!! @version Dec-1995, Sep-1997, Nov-2009
subroutine gsolve(carma, cstate, iz, previous_ice, previous_liquid, scale_threshold, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
  use shr_infnan_mod, only: shr_infnan_isnan

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(in)                  :: iz      !! z index
  real(kind=f), intent(in)             :: previous_ice(NGAS)      !! total ice at the start of substep
  real(kind=f), intent(in)             :: previous_liquid(NGAS)   !! total liquid at the start of substep
  real(kind=f)                         :: scale_threshold !! Scaling factor for convergence thresholds
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local Variables
  integer                              :: igas    !! gas index
  real(kind=f)                         :: gc_cgs
  real(kind=f)                         :: rvap
  real(kind=f)                         :: total_ice(NGAS)      ! total ice
  real(kind=f)                         :: total_liquid(NGAS)   ! total liquid
  real(kind=f)                         :: threshold            ! convergence threshold
  
  
  1 format(/,'gsolve::ERROR - negative gas concentration for ',a,' : iz=',i4,',lat=', &
              f7.2,',lon=',f7.2,',gc=',e10.3,',gasprod=',e10.3,',supsati=',e10.3, &
              ',supsatl=',e10.3,',t=',f6.2)
  2 format('gsolve::ERROR - conditions at beginning of the step : gc=',e10.3,',supsati=',e17.10, &
              ',supsatl=',e17.10,',t=',f6.2,',d_gc=',e10.3,',d_t=',f6.2)
  3 format(/,'microfast::WARNING - gas concentration change exceeds threshold: ',a,' : iz=',i4,',lat=', &
              f7.2,',lon=',f7.2, ', (gc-gcl)/gcl=', e10.3)
  4 format(/,'Oh No! Gas Destroyed:',a,' : iz=',i4,',lat=', &
              f7.2,',lon=',f7.2,',gc=',e10.3,',gasprod=',e10.3,',supsati=',e10.3, &
              ',supsatl=',e10.3,',t=',f6.2)
  

  ! Determine the total amount of condensate for each gas.
  call totalcondensate(carma, cstate, iz, total_ice, total_liquid, rc)
  
  do igas = 1,NGAS
  
    ! We do not seem to be conserving mass and energy, so rather than relying upon gasprod
    ! and rlheat, recalculate the total change in condensate to determine the change
    ! in gas and energy.
    !
    ! This is because in the old scheme, the particles were solved for implicitly, but the
    ! gas and latent heat were solved for explicitly using the same rates.

    ! NOTE: This is a hack for now, but if the gas is used for oxidation, the
    ! result can't be claculated as a difference and the claculated gasprod needs
    ! to be used. It is currently not allowed for a particle to be both growing
    ! and oxidizing if you want to exactly conserve the gas.
!    if (.not. is_oxid_gas(igas)) then

!    if ((igas /= 7) .or. (igas /= 8) .or. (igas /= 9)) then
    if (igas /= 2) then
      gasprod(igas) = ((previous_ice(igas) - total_ice(igas)) + (previous_liquid(igas) - total_liquid(igas))) / dtime
      rlprod        = rlprod - ((previous_ice(igas) - total_ice(igas)) * (rlhe(iz,igas) + rlhm(iz,igas)) + &
                       (previous_liquid(igas) - total_liquid(igas)) * (rlhe(iz,igas))) / (CP * rhoa(iz) * dtime)
      gc(iz,igas) = gc(iz,igas) + dtime * gasprod(igas)
!    else if (((igas == 7) .or. (igas == 8) .or. (igas == 9)) .and. (gc(iz,igas) > 1.0e-24_f)) then
    else if ((igas == 2) .and. (gc(iz,igas) > 1.0e-23_f)) then
!    if (gc(iz,igas) > 1.0e-24_f) then
      gasprod(igas) = ((previous_ice(igas) - total_ice(igas)) + (previous_liquid(igas) - total_liquid(igas))) / dtime
      rlprod        = rlprod - ((previous_ice(igas) - total_ice(igas)) * (rlhe(iz,igas) + rlhm(iz,igas)) + &
                       (previous_liquid(igas) - total_liquid(igas)) * (rlhe(iz,igas))) / (CP * rhoa(iz) * dtime)
      gc(iz,igas) = gc(iz,igas) + dtime * gasprod(igas)
    else
      gc(iz,igas)   = 0.0_f
      gasprod(igas) = 0.0_f
      rlprod        = 0.0_f
    end if

    ! Don't let the gas concentration go negative.

    if (gc(iz,igas) < 0.0_f) then
      if (do_substep) then
        if (nretries == maxretries) then
!         gc(iz,igas)   = gcl(iz,igas)
          gc(iz,igas)   = 0.0_f
          gasprod(igas) = 0.0_f
          rlprod        = 0.0_f
        else
          rc = RC_WARNING_RETRY
          if (do_grow_limit .and. (nretries >= glim_retries)) then
            grow_scale(iz, igas) = grow_scale(iz, igas) * glim_scale
          end if
        end if
      end if
    end if

    if (shr_infnan_isnan(gc(iz,igas))) then
      if (do_substep) then
        if (nretries == maxretries) then
!         gc(iz,igas)   = gcl(iz,igas)
          gc(iz,igas)   = 0.0_f
          gasprod(igas) = 0.0_f
          rlprod        = 0.0_f
        else
          rc = RC_WARNING_RETRY
          if (do_grow_limit .and. (nretries >= glim_retries)) then
            grow_scale(iz, igas) = grow_scale(iz, igas) * glim_scale
          end if
        end if
      end if
    end if

    ! If gas changes by too much, then retry the calculation.
    !    threshold = dgc_threshold(igas) / 0.5_f !CRT do not think it is actuall  used
    threshold = dgc_threshold(igas) / scale_threshold
    
    if (threshold /= 0._f) then
      if ((dtime * gasprod(igas) / gc(iz,igas)) > threshold) then
        if (do_substep) then
          if (nretries == maxretries) then 
            if (do_print) write(LUNOPRT,3) trim(gasname(igas)), iz, lat, lon, dtime * gasprod(igas) / gc(iz,igas)
            if (do_print) write(LUNOPRT,2) gcl(iz,igas), supsatiold(iz,igas), supsatlold(iz,igas), told(iz), d_gc(iz,igas), d_t(iz)
          end if
        else
          if (do_print) write(LUNOPRT,3) trim(gasname(igas)), iz, lat, lon, dtime * gasprod(igas) / gc(iz,igas)
        end if
  
        rc = RC_WARNING_RETRY
      end if
    end if
  end do

  ! Return to caller with new gas concentrations.
  return
end
