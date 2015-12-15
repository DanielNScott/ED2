module isotopes
use ed_max_dims , only : n_pft

!==========================================================================================!
! This module contains the nelder-mead simplicial optimization algorith, associated error  !
! functions, and the quicksort function.                                                   !
!==========================================================================================!

! Logical switches whether or not print detailed non-fatal messages.
logical :: close_nonfatalmessage = .false.

! Logical whether or not flush memory.
logical :: flushing_memory = .true.

!---------------------------------------------------------------------------------------!
! ISO_ALLOC_SCHEME -- How should the isotopes be distributed throughout the model?      !
!                      -2.  Today's resp. d13C = d13C of yesterday's assimilate.        !
!                      -1.  Overlay and propagate, first quasi-mechanistic scheme.      !
!                       0.  Use an optimization alg. to 'enforce'                       !
!                             dl - dr ~ P1                                              !
!                             dl - ds ~ P2                                              !
!                       1.  Enforce d13C bleaf = d13C assimilates, solve                !
!                             ds - dr = P1 - P2                                         !
!                                                                                       !
! ISO_P1 - ISO_P2  -- Parameters for the schemes that fix things.                       !
!                      P1:  (default: -1.50) (initially: -3.0)                          !
!                      P2:  (default: -3.20) (initially: -1.5)                          !
!                                                                                       !
! ISO_R_FLG        -- Determines what (if any) respiration fractionation occurs.        ! 
!                       0:  No respiration fractionation.                               !
!                       1:  Apparent frac. at leaves and roots from Bowling 2008        !
!---------------------------------------------------------------------------------------!
integer  :: c13af
real     :: initial_d13C
real     :: iso_P1 
real     :: iso_P2
integer  :: iso_r_flg
!---------------------------------------------------------------------------------------!

!----- Pee Dee belemnite ratio 13C:12C -------------------------------------------------!
real, parameter     :: R_std = 0.0112372 
!---------------------------------------------------------------------------------------!

!type(simtime)  :: trench_date   ! Date for trenching exp. to begin.
!real(kind=8)   :: trench_time   ! Model trench time ala grid_coms: 'time'


!---------------------------------------------------------------------------------------!
! C_ALLOC_FLG -- The C Allocation scheme.                                               !
!                   0.  Original / Normal ED 2.1 C allocation.                          !
!                   1.  New C allocation, respiration as loss from 'logical' places.    !
!---------------------------------------------------------------------------------------!
integer :: c_alloc_flg
!---------------------------------------------------------------------------------------!

!=======================================================================================!
!============================== !!!DSC!!! ==============================================!
! Tissue (c13):(C) ratios for new recruits ("C" = Total carbon!) 						         	       !
!---------------------------------------------------------------------------------------!
real, dimension(n_pft) :: cri_bdead
real, dimension(n_pft) :: cri_bleaf
real, dimension(n_pft) :: cri_broot
real, dimension(n_pft) :: cri_bsapwooda
real, dimension(n_pft) :: cri_bsapwoodb
real, dimension(n_pft) :: cri_bstorage
!=======================================================================================!
!=======================================================================================!


!-------- DS Additional parameters -----------------------------------------------------!
! RTRMULT       -- Factor multiplying the root turnover rate for PFTs 6 through 11      !
!                  (1.0 = default).                                                     !
!---------------------------------------------------------------------------------------!
real(kind=4)               :: rtrfact  ! root turnover rate multiplier
real(kind=4)               :: larprop  ! leaf assimilate respiration proportion

real(kind=4)               :: iso_lrf  ! leaf respiration fractionation
real(kind=4)               :: iso_rrf  ! root respiration fractionation
real(kind=4)               :: iso_grf  ! growth respiration fractionation
real(kind=4)               :: iso_strf ! storage respiration fractionation
real(kind=4)               :: iso_vlrf ! virtual leaf respiration fractionation
real(kind=4)               :: iso_hrf  ! heterotrophic respiration fractionation

end module isotopes
