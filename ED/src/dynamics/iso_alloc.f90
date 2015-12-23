module iso_alloc
implicit none
contains

!==========================================================================================!
subroutine pheninit_iso(bleaf, broot, bsapwooda, bsapwoodb, balive, bstorage     &
                       ,bleaf_c13, broot_c13, bsapwooda_c13, bsapwoodb_c13       &
                       ,balive_c13, bstorage_c13)
   use isotopes    	, only : cri_bleaf		     & ! intent(in)
                            , cri_broot	        & ! intent(in)
                            , cri_bsapwooda       & ! intent(in)
                            , cri_bsapwoodb       ! ! intent(in)
   implicit none
   !----- Arguments --------------------------------------------------------------------!
   real                     , intent(in)  :: bleaf             ! Leaf biomass
   real                     , intent(in)  :: broot             ! Root biomass
   real                     , intent(in)  :: bsapwooda         ! AG Sapwood biomass 
   real                     , intent(in)  :: bsapwoodb         ! BG Sapwood biomass 
   real                     , intent(in)  :: balive            ! Living tissue biomass
   real                     , intent(in)  :: bstorage          ! Storage biomass
   
   real                     , intent(out) :: bleaf_c13         ! Carbon-13 Analogues
   real                     , intent(out) :: broot_c13         !
   real                     , intent(out) :: bsapwooda_c13     !
   real                     , intent(out) :: bsapwoodb_c13     !
   real                     , intent(out) :: balive_c13        !
   real                     , intent(out) :: bstorage_c13      !
   !------------------------------------------------------------------------------------!

   !------------------------------------------------------------------------------------!
   !  Apply -27 permil ratio to everything. Recruits are small so even if this is not   !
   ! an ideal initialization it will change quickly.                                    !
   !------------------------------------------------------------------------------------!
      bleaf_c13      = bleaf     * 0.010931/(1.0 + 0.010931)
      broot_c13      = broot     * 0.010931/(1.0 + 0.010931)
      bsapwooda_c13  = bsapwooda * 0.010931/(1.0 + 0.010931)
      bsapwoodb_c13  = bsapwoodb * 0.010931/(1.0 + 0.010931)
      balive_c13     = balive    * 0.010931/(1.0 + 0.010931)
      bstorage_c13   = bstorage  * 0.010931/(1.0 + 0.010931)
   !------------------------------------------------------------------------------------!

end subroutine pheninit_iso
!==========================================================================================!


!==========================================================================================!
!==========================================================================================!
! This subroutine computes the leaf and root respiration 13C composition.                  !
!------------------------------------------------------------------------------------------!
subroutine leaf_root_resp_c13(csite,ipa)
   use ed_state_vars  , only : sitetype                  & ! structure
                             , patchtype                 ! ! structure
   use ed_misc_coms   , only : dtlsm                     & ! intent(in)
                             , current_time              & ! intent(in)
                             , frqsum                    ! ! intent(in)
   use isotopes       , only : c13af                     ! ! intent(in)
   use iso_utils      , only : hotc                      & ! intent(in)
                             , check_patch_c13           & ! intent(in)
                             , htIsoDelta                & ! intent(in)
                             , check_site_c13            ! ! intent(in)
   use consts_coms    , only : tiny_num                  & ! intent(in)
                             , umols_2_kgCyr             & ! intent(in)
                             , umol_2_kgC                ! ! intent(in)
   use consts_coms    , only : yr_day                 & ! intent(in)
                             , day_sec                ! ! intent(in)
   use growth_balive  , only : get_maintenance        & ! intent(in)
                             , apply_maintenance      & ! intent(in)
                             , apply_maintenance_c13  & ! intent(in)
                             , update_growth_resp_co  & ! intent(in)
                             , get_storage_resp       ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype)            , intent(inout) :: csite    ! Current site
   integer                   , intent(in)    :: ipa      ! Current patch number
   !----- Local variables. ----------------------------------------------------------------!
   type(patchtype), pointer :: cpatch                     ! Current patch
   integer                  :: ico                        ! Current cohort number
   real                     :: resp_loss                  ! Respiration excl. storage_resp
   real                     :: carbon_balance             ! Respiration excl. storage_resp
   real                     :: carbon_debt                ! Respiration excl. storage_resp
   real                     :: bloss_max                  ! Respiration excl. storage_resp
   real                     :: lf_bloss                   ! Respiration excl. storage_resp
   real                     :: rf_bloss                   ! Respiration excl. storage_resp
   real                     :: ff_gpp                     ! Fraction of above from gpp.
   real                     :: ratio_gpp                  ! heavy:total carbon in GPP
   real                     :: ratio_stor                 ! heavy:total carbon in storage
   real                     :: ratio_resp                 ! heavy:total carbon in resp
   real                     :: ratio_leaf                 ! heavy:total carbon in resp
   real                     :: ratio_root                 ! heavy:total carbon in resp
   real                     :: gpp_loss                   ! heavy:total carbon in resp
   real                     :: stor_loss                  ! heavy:total carbon in resp
   real                     :: leaf_loss                  ! heavy:total carbon in resp
   real                     :: root_loss                  ! heavy:total carbon in resp
   real                     :: f_leaf                     ! heavy:total carbon in resp
   real                     :: f_root                     ! heavy:total carbon in resp
   real                     :: f_lgr                      ! heavy:total carbon in resp
   real                     :: f_rgr                      ! heavy:total carbon in resp
   real                     :: f_sagr                     ! heavy:total carbon in resp
   real                     :: f_sbgr                     ! heavy:total carbon in resp
   real                     :: leaf_r_sum                     ! heavy:total carbon in resp
   real                     :: root_r_sum                     ! heavy:total carbon in resp
   real                     :: lrf                     ! heavy:total carbon in resp
   real                     :: lgf                     ! heavy:total carbon in resp
   real                     :: rrf                     ! heavy:total carbon in resp
   real                     :: rgf                       ! heavy:total carbon in resp
   real                     :: gen_loss                     ! heavy:total carbon in resp
   real                     :: leaf_resp                  ! heavy:total carbon in resp
   real                     :: root_resp                  ! heavy:total carbon in resp
   real                     :: dtlsm_o_frqsum             ! heavy:total carbon in resp
   real                     :: dtlsm_o_daysec             ! heavy:total carbon in resp
   real                     :: tfact             ! heavy:total carbon in resp
   real                     :: cb_decrement             ! heavy:total carbon in resp
   real                     :: lloss_aim             ! heavy:total carbon in resp
   real                     :: rloss_aim    ! heavy:total carbon in resp
   real                     :: lloss_overrun         ! heavy:total carbon in resp
   real                     :: rloss_overrun         ! heavy:total carbon in resp
   real                     :: leaf_int_loss         ! heavy:total carbon in resp
   real                     :: root_int_loss    ! heavy:total carbon in resp
   real                     :: leaf_respiration_c13         ! heavy:total carbon in resp
   real                     :: root_respiration_c13         ! heavy:total carbon in resp
   real                     :: carbon13_balance         ! heavy:total carbon in resp
   real                     :: flux_fact
   real                     :: flux_fact_inv
   real                     :: dtlsm_c_gain
   real                     :: dtlsm_c13_gain

   real         , dimension(30) :: patch_check_vals         ! heavy:total carbon in resp
   character(18), dimension(30) :: patch_check_labs       ! heavy:total carbon in resp
   !---------------------------------------------------------------------------------------!
   dtlsm_o_frqsum = dtlsm/frqsum  
   dtlsm_o_daysec = dtlsm/day_sec
   tfact = 1.0/yr_day
   if (c13af == 0) then
      return
   end if
      
   !------------------------------------------------------------------------------------!
   ! This subroutine has two major components:                                          !
   ! First, we determine how much respiration there is and how much each what fraction  !
   ! of the total each component comprises.                                             !
   ! Next, we determine what source (i.e. storage, leaf biomass, or root biomass) is    !
   ! being metabolised for each respiration term.                                       !
   !------------------------------------------------------------------------------------!
   
   !---------------------------------------------------------------------------------------!
   !    Loop over all cohorts.                                                             !
   !---------------------------------------------------------------------------------------!
   cpatch => csite%patch(ipa)
   cohortloop: do ico = 1,cpatch%ncohorts     
      flux_fact     = umol_2_kgC *dtlsm /cpatch%nplant(ico)
      flux_fact_inv = cpatch%nplant(ico) /umol_2_kgC /dtlsm
      
      if (current_time%time < dtlsm) then
         !write(*,*) 'New day!'
         call update_growth_resp_co(cpatch,ico)
      end if
      
      call get_maintenance(cpatch,ico,tfact*dtlsm_o_daysec,csite%can_temp(ipa),cb_decrement)             
      call apply_maintenance(cpatch,ico)
      call apply_maintenance_c13(cpatch,ico)
      call get_storage_resp(cpatch,ico,tfact)
      
      cpatch%bstorage(ico) = cpatch%bstorage(ico)                                          & 
                           - cpatch%leaf_storage_resp(ico)                                 &
                           - cpatch%root_storage_resp(ico)                                 &
                           - cpatch%sapa_storage_resp(ico)                                 &
                           - cpatch%sapb_storage_resp(ico)
      cpatch%bstorage_c13(ico) = cpatch%bstorage_c13(ico)                      &
                               - cpatch%leaf_storage_resp_c13(ico)             &
                               - cpatch%root_storage_resp_c13(ico)             &
                               - cpatch%sapa_storage_resp_c13(ico)             &
                               - cpatch%sapb_storage_resp_c13(ico)

      !------------------------------------------------------------------------------------!
      ! Standardize respiration units and get total respiration                            !
      !------------------------------------------------------------------------------------!
      leaf_resp = cpatch%leaf_respiration(ico) *flux_fact
      root_resp = cpatch%root_respiration(ico) *flux_fact

      resp_loss = cpatch%leaf_growth_resp(ico) + cpatch%root_growth_resp(ico)              &
                + cpatch%sapa_growth_resp(ico) + cpatch%sapb_growth_resp(ico)              &
                + leaf_resp                    + root_resp

      !------------------------------------------------------------------------------------!
      ! Find the proportion of total respiration each term comprises:                      !
      !------------------------------------------------------------------------------------!
      f_leaf = max(min(hotc(leaf_resp                   ,resp_loss),1.0),0.0)
      f_root = max(min(hotc(root_resp                   ,resp_loss),1.0),0.0)
      f_lgr  = max(min(hotc(cpatch%leaf_growth_resp(ico),resp_loss),1.0),0.0)
      f_rgr  = max(min(hotc(cpatch%root_growth_resp(ico),resp_loss),1.0),0.0)
      f_sagr = max(min(hotc(cpatch%sapa_growth_resp(ico),resp_loss),1.0),0.0)
      f_sbgr = max(min(hotc(cpatch%sapb_growth_resp(ico),resp_loss),1.0),0.0)
      
      !------------------------------------------------------------------------------------!
      ! Find the proportions to partition leaf-based and root-based resp.                  !
      !------------------------------------------------------------------------------------!
      leaf_r_sum = leaf_resp + cpatch%leaf_growth_resp(ico)
      root_r_sum = root_resp + cpatch%root_growth_resp(ico)
                      
      lrf = hotc(leaf_resp                   ,leaf_r_sum)
      lgf = hotc(cpatch%leaf_growth_resp(ico),leaf_r_sum)
      rrf = hotc(root_resp                   ,root_r_sum)
      rgf = hotc(cpatch%root_growth_resp(ico),root_r_sum)
      
      dtlsm_c_gain   =  cpatch%gpp(ico)*flux_fact - leaf_resp - root_resp
      carbon_balance = (cpatch%gpp(ico)*flux_fact) - resp_loss
      carbon_debt    = 0.0

      !------------------------------------------------------------------------------------!
      ! When we call get_c_xfers() in growth_balive.f90, we will have 1 of three cases:    !
      ! 1) (carbon_balance > 0)                         => Add new carbon to plant pools   !
      ! 2) (carbon_balance < 0 & plants should grow)    => Lose storage first              !
      ! 3) (carbon_balance < 0 & plants shouldn't grow) => Lose tissues first              !
      !                                                                                    !
      ! We determine which will occur and set flux 13C content appropriately. Note the     !
      ! following variables which determine branching in get_c_xfers() later:              !
      !                                                                                    !
      ! available_carbon = cpatch%bstorage(ico) + carbon_balance                           !
      ! time_to_flush    = carbon_balance > 0.0 .or.                                     & !
      !                 ( available_carbon > 0.0 .and. cpatch%phenology_status(ico) == 1 ) !
      !------------------------------------------------------------------------------------!
      ratio_gpp  = hotc(cpatch%gpp_c13     (ico),cpatch%gpp     (ico))
      ratio_leaf = hotc(cpatch%bleaf_c13   (ico),cpatch%bleaf   (ico))
      ratio_root = hotc(cpatch%broot_c13   (ico),cpatch%broot   (ico))
      ratio_stor = hotc(cpatch%bstorage_c13(ico),cpatch%bstorage(ico))
      
      gpp_loss  = 0.0
      stor_loss = 0.0
      leaf_loss = 0.0
      root_loss = 0.0

      bloss_max   = cpatch%bleaf(ico) + cpatch%broot(ico)
      lf_bloss    = cpatch%bleaf(ico) /bloss_max
      rf_bloss    = cpatch%broot(ico) /bloss_max
      
      lloss_aim     = 0.0
      rloss_aim     = 0.0
      leaf_int_loss = 0.0
      root_int_loss = 0.0
      lloss_overrun = 0.0
      rloss_overrun = 0.0

      if (carbon_balance > 0.0) then
         !---------------------------------------------------------------------------------!
         ! In get_c_xfers we will simply put carbon_balance into some pool(s), and since   !
         ! carbon_balance = gpp - non_stor_resp_loss > 0.0                                 !
         !                                                                                 !
         ! We should set the signature of the respiration to the signature of gpp.         !
         ! Note that this also implies carbon13_balance > 0.0.                             !
         !---------------------------------------------------------------------------------!
         gpp_loss  = ratio_gpp *resp_loss
         
      !elseif (cpatch%phenology_status(ico) == 1 .and. (carbon_balance <= 0.0)) then
         !---------------------------------------------------------------------------------!
         ! "time_to_flush" will be .true. in get_c_xfers and we will move storage into     !
         ! plant pools despite having lost all gpp and some storage to respiration.        !
         !---------------------------------------------------------------------------------!
      !   gpp_loss  = cpatch%gpp_c13(ico) *flux_fact
      !   stor_loss = -1.0*carbon_balance *ratio_stor
         
      !elseif (cpatch%phenology_status(ico) /= 1 .and. (carbon_balance <= 0.0)) then
      elseif (carbon_balance <= 0.0) then
         !---------------------------------------------------------------------------------!
         ! In this case we aren't growing any tissues, we are only taking carbon out.      !
         !                                                                                 !
         ! Recall the phenology_status codes...           ... which later imply:           !
         !  0 - plant has the maximum LAI, given its size  => Lose storage first           !
         !  1 - plant is growing leaves                    => Lose storage first           !
         ! -1 - plant is actively dropping leaves          => Lose biomass first           !
         ! -2 - plant has no leaves                        => Lose biomass first           !
         !---------------------------------------------------------------------------------!
         carbon_debt = -1.0 *carbon_balance
         gpp_loss    = cpatch%gpp_c13(ico) *flux_fact

         select case (cpatch%phenology_status(ico))
         case (0,1)
            if (cpatch%bstorage(ico) > carbon_debt) then
               !--------------------------------------------------------------------------!
               ! Only some storage loss necessary, no tissue loss necessary.              !
               !--------------------------------------------------------------------------!
               stor_loss = carbon_debt *ratio_stor
               
            else
               !--------------------------------------------------------------------------!
               ! All storage lost, then some leaves + roots.                              !
               ! 13C_removed = var_13C:C * proportion_from_var * total_C_being_removed    !
               !--------------------------------------------------------------------------!
               carbon_debt = carbon_debt - cpatch%bstorage(ico)
               stor_loss   = cpatch%bstorage_c13(ico)
               
               if (bloss_max > carbon_debt) then
                  !--------------------------------------------------------------------------!
                  ! All respiration has leaf or root substrate.                              !
                  ! 13C_removed = var_13C:C * proportion_from_var * total_C_being_removed    !
                  !--------------------------------------------------------------------------!               
                  lloss_aim   = f_leaf*carbon_debt
                  rloss_aim   = f_root*carbon_debt

                  ! Only one of these will be > 0, otherwise carbon_debt > bloss_max
                  lloss_overrun = max(lloss_aim - cpatch%bleaf(ico),0.0)
                  rloss_overrun = max(rloss_aim - cpatch%broot(ico),0.0)

                  leaf_int_loss = min(cpatch%bleaf(ico),lloss_aim) 
                  root_int_loss = min(cpatch%broot(ico),rloss_aim)
    
                  leaf_loss = leaf_int_loss *ratio_leaf + rloss_overrun *ratio_root
                  root_loss = root_int_loss *ratio_root + lloss_overrun *ratio_leaf
               else
                  leaf_loss = cpatch%bleaf_c13(ico)
                  root_loss = cpatch%broot_c13(ico)
               end if
            end if
         case (-2,-1)
            if (bloss_max > carbon_debt) then
               !--------------------------------------------------------------------------!
               ! All respiration has leaf or root substrate.                              !
               ! 13C_removed = var_13C:C * proportion_from_var * total_C_being_removed    !
               !--------------------------------------------------------------------------!               
               lloss_aim   = f_leaf*carbon_debt
               rloss_aim   = f_root*carbon_debt

               ! Only one of these will be > 0, otherwise carbon_debt > bloss_max
               lloss_overrun = max(lloss_aim - cpatch%bleaf(ico),0.0)
               rloss_overrun = max(rloss_aim - cpatch%broot(ico),0.0)

               leaf_int_loss = min(cpatch%bleaf(ico),lloss_aim) 
               root_int_loss = min(cpatch%broot(ico),rloss_aim)
 
               leaf_loss = leaf_int_loss *ratio_leaf + rloss_overrun *ratio_root
               root_loss = root_int_loss *ratio_root + lloss_overrun *ratio_leaf
               
               !leaf_loss = ratio_leaf * lf_bloss * min(bloss_max,carbon_debt)
               !root_loss = ratio_root * rf_bloss * min(bloss_max,carbon_debt)
            else
               !--------------------------------------------------------------------------!
               ! All leaves and roots will go, then some storage.                         !
               ! 13C_removed = var_13C:C * proportion_from_var * total_C_being_removed    !
               !--------------------------------------------------------------------------!
               stor_loss = ratio_stor * min(cpatch%bstorage(ico), carbon_debt - bloss_max)
               leaf_loss = cpatch%bleaf_c13(ico)
               root_loss = cpatch%broot_c13(ico)
               !--------------------------------------------------------------------------!
            end if
         end select
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! Now we partition the total 13C loss sensibly among the respiration terms by...     !      
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! 3) Dividing up storage and gpp based resp evenly and dividing leaf loss and root   !
      !    loss evenly into leaf-based partitions and root-based partitions respectively.  ! 
      !------------------------------------------------------------------------------------!
      gen_loss = stor_loss + gpp_loss

      leaf_respiration_c13 = f_leaf*(gen_loss) + leaf_loss*lrf
      root_respiration_c13 = f_root*(gen_loss) + root_loss*rrf
      
      cpatch%leaf_respiration_c13(ico) = leaf_respiration_c13 * flux_fact_inv
      cpatch%root_respiration_c13(ico) = root_respiration_c13 * flux_fact_inv

      cpatch%leaf_growth_resp_c13(ico) = f_lgr *(gen_loss) + leaf_loss*lgf
      cpatch%root_growth_resp_c13(ico) = f_rgr *(gen_loss) + root_loss*rgf
      cpatch%sapa_growth_resp_c13(ico) = f_sagr*(gen_loss)
      cpatch%sapb_growth_resp_c13(ico) = f_sbgr*(gen_loss)
      
      dtlsm_c13_gain   = cpatch%gpp_c13(ico) * flux_fact                                   &
                       - leaf_respiration_c13                                              &
                       - root_respiration_c13
      carbon13_balance = cpatch%gpp_c13(ico) * flux_fact                                   &
                       - leaf_respiration_c13                                              &
                       - root_respiration_c13                                              &
                       - cpatch%leaf_growth_resp_c13(ico)                                  &
                       - cpatch%root_growth_resp_c13(ico)                                  &
                       - cpatch%sapa_growth_resp_c13(ico)                                  &
                       - cpatch%sapb_growth_resp_c13(ico)
      !------------------------------------------------------------------------------------!

      
      !------------------------------------------------------------------------------------!
      ! Save the loss partitioning so it can be accounted for in get_c13_xfers later.      !
      !------------------------------------------------------------------------------------!
      cpatch%bleaf_c13_loss(ico) = max((leaf_int_loss + lloss_overrun)*ratio_leaf,0.0)
      cpatch%broot_c13_loss(ico) = max((root_int_loss + rloss_overrun)*ratio_root,0.0)
      cpatch%bstor_c13_loss(ico) = max(stor_loss,0.0)
      !------------------------------------------------------------------------------------!

      
      !----- The output variables must be in [kgC/plant/yr]. ---------------------!
      cpatch%fmean_leaf_resp_c13(ico) = cpatch%fmean_leaf_resp_c13 (ico)          &
                                      + cpatch%leaf_respiration_c13(ico)          &
                                      * dtlsm_o_frqsum * umols_2_kgCyr            &
                                      / cpatch%nplant          (ico)
      cpatch%fmean_root_resp_c13(ico) = cpatch%fmean_root_resp_c13 (ico)          &
                                      + cpatch%root_respiration_c13(ico)          &
                                      * dtlsm_o_frqsum * umols_2_kgCyr            &
                                      / cpatch%nplant          (ico)
      !---------------------------------------------------------------------------!
      !write(*,*) htIsoDelta(cpatch%leaf_respiration_c13(ico),cpatch%leaf_respiration(ico)) &
      !          ,htIsoDelta(cpatch%root_respiration_c13(ico),cpatch%root_respiration(ico))
      patch_check_vals = &
     (/carbon_balance, carbon13_balance,cpatch%gpp(ico)*flux_fact,cpatch%gpp_c13(ico)*flux_fact &
      ,     stor_loss,         gpp_loss,                leaf_loss,                    root_loss &
      ,        f_leaf,           f_root,                    f_lgr,                        f_rgr &
      ,           lrf,              rrf,                      lgf,                          rgf &
      ,        f_sagr,           f_sbgr,                 lf_bloss,                     rf_bloss &
      , leaf_int_loss,    root_int_loss,            lloss_overrun,                rloss_overrun &
      ,     bloss_max, real(cpatch%phenology_status(ico))                                       &
      , htIsoDelta(abs(dtlsm_c13_gain)  ,abs(dtlsm_c_gain))     &
      , htIsoDelta(abs(carbon13_balance),abs(carbon_balance)), dtlsm_c_gain, dtlsm_c13_gain      /)

      patch_check_labs = &
     (/'    carbon_balance','  carbon13_balance','               gpp','           gpp_c13'&
      ,'         stor_loss','          gpp_loss','         leaf_loss','         root_loss'&
      ,'   gen_frac_leaf_r','   gen_frac_root_r','  gen_frac_leaf_gr','  gen_frac_root_gr'&
      ,'  leaf_frac_leaf_r','  root_frac_root_r',' leaf_frac_leaf_gr',' root_frac_root_gr'&
      ,'  gen_frac_sapa_gr','  gen_frac_sapb_gr','          lf_bloss','          rf_bloss'&
      ,'     leaf_int_loss','     root_int_loss','     lloss_overrun','     rloss_overrun'&
      ,'         bloss_max','  phenology_status','          dcg_d13C','   carbon_bal_d13C'&
      ,'      dtlsm_c_gain','    dtlsm_c13_gain'/)

 
      call check_patch_c13(cpatch,ico,'leaf_root_resp_c13','iso_alloc.f90',patch_check_vals&
                          ,patch_check_labs &
                          ,(/1,1,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,3/))
   end do cohortloop
   call check_site_c13(csite,ipa,'leaf_root_resp_c13','iso_alloc.f90')

end subroutine leaf_root_resp_c13
!==========================================================================================!
!==========================================================================================!



!==========================================================================================!
!==========================================================================================!
! This subroutine computes the leaf and root respiration 13C composition.                  !
!------------------------------------------------------------------------------------------!
subroutine leaf_root_resp_alt_c13(csite,ipa)
   use ed_state_vars  , only : sitetype                  & ! structure
                             , patchtype                 ! ! structure
   use ed_misc_coms   , only : dtlsm                     & ! intent(in)
                             , frqsum                    ! ! intent(in)
   use isotopes       , only : c13af                     ! ! intent(in)
   use consts_coms    , only : tiny_num                  & ! intent(in)
                             , umols_2_kgCyr             ! ! intent(in)
   use iso_utils      , only : hotc                      & ! intent(in)
                             , check_patch_c13           & ! intent(in)
                             , check_site_c13            ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype)            , intent(inout) :: csite    ! Current site
   integer                   , intent(in)    :: ipa      ! Current patch number
   !----- Local variables. ----------------------------------------------------------------!
   type(patchtype), pointer :: cpatch                     ! Current patch
   integer                  :: ico                        ! Current cohort number
   real                     :: non_stor_resp              ! Respiration excl. storage_resp
   real                     :: ff_gpp                     ! Fraction of above from gpp.
   real                     :: ratio_gpp                  ! heavy:total carbon in GPP
   real                     :: ratio_stor                 ! heavy:total carbon in storage
   real                     :: ratio_resp                 ! heavy:total carbon in resp
   real                     :: dtlsm_o_frqsum             ! heavy:total carbon in resp
   !---------------------------------------------------------------------------------------!
   dtlsm_o_frqsum = dtlsm/frqsum  

 
   !---------------------------------------------------------------------------------------!
   !    Loop over all cohorts.                                                             !
   !---------------------------------------------------------------------------------------!
   cpatch => csite%patch(ipa)
   cohortloop: do ico = 1,cpatch%ncohorts
      !------------------------------------------------------------------------------------!
      !     Determine if the sum of respiration, excluding [var]_storage_resp terms, is    !
      ! greater than gpp. If so, some fraction of this respiration must come from plant    !
      ! storage, in which case we implement a 2-pool mixing model.                         ! 
      !------------------------------------------------------------------------------------!
      if (c13af > 0) then
         non_stor_resp = cpatch%leaf_respiration(ico) + cpatch%root_respiration(ico)      &
                       + cpatch%leaf_growth_resp(ico) + cpatch%root_growth_resp(ico)      &
                       + cpatch%sapa_growth_resp(ico) + cpatch%sapb_growth_resp(ico)
         
         if (non_stor_resp > tiny_num) then
            ff_gpp     = min(cpatch%gpp(ico) /non_stor_resp, 1.0)
            ratio_gpp  = hotc(     cpatch%gpp_c13(ico),     cpatch%gpp(ico))
            ratio_stor = hotc(cpatch%bstorage_c13(ico),cpatch%bstorage(ico))
            ratio_resp = ff_gpp*ratio_gpp + (1.0 - ff_gpp)*ratio_stor
         else
            ratio_resp = 0.0
         end if
         
         cpatch%leaf_respiration_c13(ico) = cpatch%leaf_respiration(ico) *ratio_resp
         cpatch%root_respiration_c13(ico) = cpatch%root_respiration(ico) *ratio_resp
         cpatch%leaf_growth_resp_c13(ico) = cpatch%leaf_growth_resp(ico) *ratio_resp
         cpatch%root_growth_resp_c13(ico) = cpatch%root_growth_resp(ico) *ratio_resp
         cpatch%sapa_growth_resp_c13(ico) = cpatch%sapa_growth_resp(ico) *ratio_resp
         cpatch%sapb_growth_resp_c13(ico) = cpatch%sapb_growth_resp(ico) *ratio_resp

         !----- The output variables must be in [kgC/plant/yr]. ---------------------!
         cpatch%fmean_leaf_resp_c13(ico) = cpatch%fmean_leaf_resp_c13 (ico)          &
                                         + cpatch%leaf_respiration_c13(ico)          &
                                         * dtlsm_o_frqsum * umols_2_kgCyr            &
                                         / cpatch%nplant          (ico)
         cpatch%fmean_root_resp_c13(ico) = cpatch%fmean_root_resp_c13 (ico)          &
                                         + cpatch%root_respiration_c13(ico)          &
                                         * dtlsm_o_frqsum * umols_2_kgCyr            &
                                         / cpatch%nplant          (ico)
         !---------------------------------------------------------------------------!
         
         call check_patch_c13(cpatch,ico,'leaf_root_resp_c13','iso_alloc.f90')
         !---------------------------------------------------------------------------!
      else
         return
      end if
      !------------------------------------------------------------------------------!
   end do cohortloop
   call check_site_c13(csite,ipa,'leaf_root_resp_c13','iso_alloc.f90')

end subroutine leaf_root_resp_alt_c13
!==========================================================================================!
!==========================================================================================!

end module iso_alloc
