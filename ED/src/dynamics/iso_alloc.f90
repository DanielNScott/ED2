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
                             , simtime                   & ! intent(in)
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
   type(simtime)            :: dummytime                  ! Current time after next update
   real                     :: resp_loss                  ! Respiration excl. storage_resp
   real                     :: carbon_balance             ! Respiration excl. storage_resp
   real                     :: carbon_debt                ! Respiration excl. storage_resp
   real                     :: bloss_max                  ! Respiration excl. storage_resp
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

   real                     :: f_lr                       ! heavy:total carbon in resp
   real                     :: f_rr                       ! heavy:total carbon in resp
   real                     :: f_lgr                      ! heavy:total carbon in resp
   real                     :: f_rgr                      ! heavy:total carbon in resp
   real                     :: f_sagr                     ! heavy:total carbon in resp
   real                     :: f_sbgr                     ! heavy:total carbon in resp

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

   real                     :: gpp_2_lrr
   real                     :: gpp_2_gr

   real                     :: stor_2_lrr
   real                     :: stor_2_gr

   real                     :: leaf_loss_2_lr
   real                     :: root_loss_2_rr
   real                     :: leaf_loss_2_lgr
   real                     :: root_loss_2_rgr

   real                     :: l2lrsum
   real                     :: r2lrsum
   real                     :: lgr2gr
   real                     :: rgr2gr
   real                     :: sagr2gr
   real                     :: sbgr2gr

   real                     :: f_debt_2_lrr
   real                     :: growth_resp
   real                     :: lrresp_residual

   real                     :: lgr2lgr_plus_sagr
   real                     :: sagr2lgr_plus_sagr
   real                     :: rgr2rgr_plus_sbgr
   real                     :: sbgr2rgr_plus_sbgr

   real         , dimension(46) :: patch_check_vals         ! heavy:total carbon in resp
   character(18), dimension(46) :: patch_check_labs       ! heavy:total carbon in resp
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
      
      dummytime = current_time
      call update_model_time_dm(dummytime,dtlsm)
      if (dummytime%time < dtlsm) then
         call update_growth_resp_co(cpatch,ico)
      end if
      
      call get_maintenance(cpatch,ico,tfact*dtlsm_o_daysec,csite%can_temp(ipa),cb_decrement)             
      call apply_maintenance(cpatch,ico)
      call apply_maintenance_c13(cpatch,ico)
      call get_storage_resp(cpatch,ico,tfact*dtlsm_o_daysec)
      
      cpatch%bstorage(ico)     = cpatch%bstorage(ico)                                          & 
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
      
      growth_resp= cpatch%leaf_growth_resp(ico) + cpatch%root_growth_resp(ico)              &
                 + cpatch%sapa_growth_resp(ico) + cpatch%sapb_growth_resp(ico)

      resp_loss = growth_resp + leaf_resp + root_resp

      !------------------------------------------------------------------------------------!
      ! Find the proportion of total respiration each term comprises:                      !
      !------------------------------------------------------------------------------------!
      l2lrsum = max(min(hotc(leaf_resp ,leaf_resp+root_resp),1.0),0.0)
      r2lrsum = max(min(hotc(root_resp ,leaf_resp+root_resp),1.0),0.0)
      
      lgr2gr  = max(min(hotc(cpatch%leaf_growth_resp(ico),growth_resp),1.0),0.0)
      rgr2gr  = max(min(hotc(cpatch%root_growth_resp(ico),growth_resp),1.0),0.0)
      sagr2gr = max(min(hotc(cpatch%sapa_growth_resp(ico),growth_resp),1.0),0.0)
      sbgr2gr = max(min(hotc(cpatch%sapb_growth_resp(ico),growth_resp),1.0),0.0)
      
      
      !------------------------------------------------------------------------------------!
      ! Find the proportions to partition leaf-based and root-based resp.                  !
      !------------------------------------------------------------------------------------!
      
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
      
      gpp_loss   = 0.0
      gpp_2_lrr  = 0.0
      gpp_2_gr   = 0.0
      stor_loss  = 0.0
      stor_2_lrr = 0.0
      stor_2_gr  = 0.0
      leaf_loss  = 0.0
      root_loss  = 0.0
            
      leaf_loss_2_lr  = 0.0
      root_loss_2_rr  = 0.0
      leaf_loss_2_lgr = 0.0
      root_loss_2_rgr = 0.0

      bloss_max   = cpatch%bleaf(ico) + cpatch%broot(ico)
      
      lloss_aim     = 0.0
      rloss_aim     = 0.0
      leaf_int_loss = 0.0
      root_int_loss = 0.0
      lloss_overrun = 0.0
      rloss_overrun = 0.0
      
      ! if dtlsm_c_gain < 0 then leaf and root resp are using more than just gpp.                           
      f_debt_2_lrr = max(min(dtlsm_c_gain/carbon_balance,1.0),0.0)
      
      if (dtlsm_c_gain > 0.0) then
         ! None of leaf or root biomass is going to leaf_respiration or root_respiration.
         f_lr = 0.0
         f_rr = 0.0

         if (carbon_balance > 0.0) then
            ! There will be no leaf loss or root loss, so these fractions don't mean
            ! anything and we can set them to 0.
            f_lgr  = 0.0 
            f_rgr  = 0.0
            f_sagr = 0.0 
            f_sbgr = 0.0 
         else
            ! Leaf or root loss will be split between growth resp terms evenly
            f_lgr  = lgr2gr
            f_rgr  = rgr2gr
            f_sagr = sagr2gr
            f_sbgr = sbgr2gr
         end if
      else
         ! Leaf/Root loss will be split evenly between the fraction of leaf_respiration not
         ! already accounted for and the other above/below-ground growth resps.
         lrresp_residual = abs(dtlsm_c_gain)
         f_lr = l2lrsum *hotc(abs(dtlsm_c_gain),growth_resp + abs(dtlsm_c_gain))
         f_rr = r2lrsum *hotc(abs(dtlsm_c_gain),growth_resp + abs(dtlsm_c_gain))
        
         f_lgr  = hotc(cpatch%leaf_growth_resp(ico),growth_resp + abs(dtlsm_c_gain)) 
         f_rgr  = hotc(cpatch%root_growth_resp(ico),growth_resp + abs(dtlsm_c_gain)) 
         f_sagr = hotc(cpatch%sapa_growth_resp(ico),growth_resp + abs(dtlsm_c_gain)) 
         f_sbgr = hotc(cpatch%sapb_growth_resp(ico),growth_resp + abs(dtlsm_c_gain)) 

      end if

      if (carbon_balance > 0.0) then
         !---------------------------------------------------------------------------------!
         ! In get_c_xfers we will simply put carbon_balance into some pool(s), and since   !
         ! carbon_balance = gpp - non_stor_resp_loss > 0.0                                 !
         !                                                                                 !
         ! We should set the signature of the respiration to the signature of gpp.         !
         ! Note that this also implies carbon13_balance > 0.0.                             !
         !---------------------------------------------------------------------------------!
         gpp_loss  = ratio_gpp *resp_loss

         gpp_2_lrr = ratio_gpp *(leaf_resp + root_resp)
         gpp_2_gr  = ratio_gpp *growth_resp
         
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
                  lloss_aim   = (f_lr + f_lgr + f_sagr)*carbon_debt
                  rloss_aim   = (f_rr + f_rgr + f_sbgr)*carbon_debt

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
               lloss_aim   = (f_lr + f_lgr + f_sagr)*carbon_debt
               rloss_aim   = (f_rr + f_rgr + f_sbgr)*carbon_debt

               ! Only one of these will be > 0, otherwise carbon_debt > bloss_max
               lloss_overrun = max(lloss_aim - cpatch%bleaf(ico),0.0)
               rloss_overrun = max(rloss_aim - cpatch%broot(ico),0.0)

               leaf_int_loss = min(cpatch%bleaf(ico),lloss_aim) 
               root_int_loss = min(cpatch%broot(ico),rloss_aim)
 
               leaf_loss = leaf_int_loss *ratio_leaf + rloss_overrun *ratio_root
               root_loss = root_int_loss *ratio_root + lloss_overrun *ratio_leaf
               
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
         
         ! Either dtlsm_c_gain > 0 and therefore it covers leaf and root resp, or 
         if (dtlsm_c_gain > 0.0) then
            

            gpp_2_lrr  = gpp_loss *max(min(     (leaf_resp+root_resp)/(cpatch%gpp(ico)*flux_fact) ,1.0),0.0)
            gpp_2_gr   = gpp_loss *max(min((1.0-(leaf_resp+root_resp)/(cpatch%gpp(ico)*flux_fact)),1.0),0.0)
            stor_2_lrr = 0.0
            stor_2_gr  = stor_loss

            leaf_loss_2_lr  = 0.0
            root_loss_2_rr  = 0.0
            leaf_loss_2_lgr = leaf_loss
            root_loss_2_rgr = root_loss

         else
            gpp_2_lrr  = gpp_loss
            gpp_2_gr   = 0.0
            stor_2_lrr = stor_loss *       f_debt_2_lrr
            stor_2_gr  = stor_loss *(1.0 - f_debt_2_lrr)
           
            leaf_loss_2_lr  = leaf_loss * hotc(           f_lr , (f_lr + f_sagr + f_lgr))
            root_loss_2_rr  = root_loss * hotc(           f_rr , (f_rr + f_sbgr + f_rgr))
            leaf_loss_2_lgr = leaf_loss * hotc((f_lgr + f_sagr), (f_lr + f_sagr + f_lgr))
            root_loss_2_rgr = root_loss * hotc((f_rgr + f_sbgr), (f_rr + f_sbgr + f_rgr))
         end if 
      end if
      !------------------------------------------------------------------------------------!
      lgr2lgr_plus_sagr  = max(min(cpatch%leaf_growth_resp(ico)/(cpatch%leaf_growth_resp(ico)+cpatch%sapa_growth_resp(ico)),1.0),0.0)
      sagr2lgr_plus_sagr = max(min(cpatch%sapa_growth_resp(ico)/(cpatch%leaf_growth_resp(ico)+cpatch%sapa_growth_resp(ico)),1.0),0.0)
          
      rgr2rgr_plus_sbgr  = max(min(cpatch%root_growth_resp(ico)/(cpatch%root_growth_resp(ico)+cpatch%sapb_growth_resp(ico)),1.0),0.0)
      sbgr2rgr_plus_sbgr = max(min(cpatch%sapb_growth_resp(ico)/(cpatch%root_growth_resp(ico)+cpatch%sapb_growth_resp(ico)),1.0),0.0)


      !------------------------------------------------------------------------------------!
      ! Now we partition the total 13C loss sensibly among the respiration terms by...     !      
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! 3) Dividing up storage and gpp based resp evenly and dividing leaf loss and root   !
      !    loss evenly into leaf-based partitions and root-based partitions respectively.  ! 
      !------------------------------------------------------------------------------------!
      ! if dtlsm_c_gain > 0 then gpp covers leaf & root resp.
      ! if dtlsm_c_gain < 0 then gpp does not cover leaf and root resp
      
      leaf_respiration_c13 = (gpp_2_lrr + stor_2_lrr)*l2lrsum + leaf_loss_2_lr
      root_respiration_c13 = (gpp_2_lrr + stor_2_lrr)*r2lrsum + root_loss_2_rr

      cpatch%leaf_respiration_c13(ico) = leaf_respiration_c13 * flux_fact_inv
      cpatch%root_respiration_c13(ico) = root_respiration_c13 * flux_fact_inv

      cpatch%leaf_growth_resp_c13(ico) = (gpp_2_gr + stor_2_gr)*lgr2gr  + leaf_loss_2_lgr *lgr2lgr_plus_sagr
      cpatch%root_growth_resp_c13(ico) = (gpp_2_gr + stor_2_gr)*rgr2gr  + root_loss_2_rgr *rgr2rgr_plus_sbgr
      cpatch%sapa_growth_resp_c13(ico) = (gpp_2_gr + stor_2_gr)*sagr2gr + leaf_loss_2_lgr *sagr2lgr_plus_sagr
      cpatch%sapb_growth_resp_c13(ico) = (gpp_2_gr + stor_2_gr)*sbgr2gr + root_loss_2_rgr *sbgr2rgr_plus_sbgr
     
      if (f_debt_2_lrr < tiny(1.0) .and. &
             abs(carbon_debt) > tiny(1.0) .and. &
          htIsoDelta(cpatch%leaf_growth_resp_c13(ico),cpatch%leaf_growth_resp(ico)) /= ratio_leaf) then
          dtlsm_c13_gain = 10000
      end if
 
      !------------------------------------------------------------------------------------!
      ! 3) Dividing up storage and gpp based resp evenly and dividing leaf loss and root   !
      !    loss evenly into leaf-based partitions and root-based partitions respectively.  ! 
      !------------------------------------------------------------------------------------!
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
      ,        f_lr,           f_rr,                    f_lgr,                        f_rgr &
      ,        f_sagr,           f_sbgr                 &
      , leaf_int_loss,    root_int_loss,            lloss_overrun,                rloss_overrun &
      ,     bloss_max, real(cpatch%phenology_status(ico))                                       &
      , htIsoDelta(abs(dtlsm_c13_gain)  ,abs(dtlsm_c_gain))     &
      , htIsoDelta(abs(carbon13_balance),abs(carbon_balance)), dtlsm_c_gain, dtlsm_c13_gain     &
      , gpp_2_lrr    ,         gpp_2_gr,               stor_2_lrr,                    stor_2_gr &
      ,       l2lrsum,          r2lrsum,                   lgr2gr,                       rgr2gr &
      ,       sagr2gr,          sbgr2gr,           leaf_loss_2_lr,               root_loss_2_rr &
      ,leaf_loss_2_lgr, root_loss_2_rgr,              growth_resp,          f_debt_2_lrr &
      ,             dtlsm_c_gain,               dtlsm_c13_gain &
      ,lgr2lgr_plus_sagr, sagr2lgr_plus_sagr,rgr2rgr_plus_sbgr,sbgr2rgr_plus_sbgr/)

      patch_check_labs = &
     (/'    carbon_balance','  carbon13_balance','               gpp','           gpp_c13'&
      ,'         stor_loss','          gpp_loss','         leaf_loss','         root_loss'&
      ,'   gen_frac_leaf_r','   gen_frac_root_r','  gen_frac_leaf_gr','  gen_frac_root_gr'&
      ,'  gen_frac_sapa_gr','  gen_frac_sapb_gr'&
      ,'     leaf_int_loss','     root_int_loss','     lloss_overrun','     rloss_overrun'&
      ,'         bloss_max','  phenology_status','          dcg_d13C','   carbon_bal_d13C'&
      ,'      dtlsm_c_gain','    dtlsm_c13_gain'                                          &
      ,'         gpp_2_lrr','          gpp_2_gr','        stor_2_lrr','         stor_2_gr'&
      ,'           l2lrsum','           r2lrsum','            lgr2gr','            rgr2gr'&
      ,'           sagr2gr','           sbgr2gr','    leaf_loss_2_lr','    root_loss_2_rr'&
      ,'   leaf_loss_2_lgr','   root_loss_2_rgr','       growth_resp','    r_c_debt_2_lrr'&
      ,'      daily_c_gain','    daily_c13_gain'&
      ,' lgr2lgr_plus_sagr','sagr2lgr_plus_sagr',' rgr2rgr_plus_sbgr','sbgr2rgr_plus_sbgr'/)

 
      call check_patch_c13(cpatch,ico,'leaf_root_resp_c13','iso_alloc.f90',patch_check_vals&
                          ,patch_check_labs &
                          ,(/1,1,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,3   &
                            ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,3/))
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
