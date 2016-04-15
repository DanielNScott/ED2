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

   real                     :: resp_sum                   ! Respiration excl. storage_resp

   real                     :: carbon_balance             ! Respiration excl. storage_resp
   real                     :: carbon_debt                ! Respiration excl. storage_resp

   real                     :: bloss_max                  ! Respiration excl. storage_resp

   real                     :: ratio_gpp                  ! GPP
   real                     :: ratio_stor                 ! storage
   real                     :: ratio_leaf                 ! resp
   real                     :: ratio_root                 ! resp

   real                     :: gpp_loss                   ! resp
   real                     :: stor_loss                  ! resp
   real                     :: leaf_loss                  ! resp
   real                     :: root_loss                  ! resp

   real                     :: f_lr                       ! resp
   real                     :: f_rr                       ! resp
   real                     :: f_lgr                      ! resp
   real                     :: f_rgr                      ! resp
   real                     :: f_sagr                     ! resp
   real                     :: f_sbgr                     ! resp

   real                     :: leaf_resp                  ! In kgC/pl not umolC/m2/s
   real                     :: root_resp                  ! In kgC/pl not umolC/m2/s
   real                     :: leaf_resp_c13              ! In kgC/pl not umolC/m2/s
   real                     :: root_resp_c13              ! In kgC/pl not umolC/m2/s
   
   real                     :: dtlsm_o_frqsum             ! resp
   real                     :: dtlsm_o_daysec             ! resp
   
   real                     :: tfact                      ! resp
   real                     :: cb_decrement             ! resp

   real                     :: lloss_aim             ! resp
   real                     :: rloss_aim             ! resp

   real                     :: lloss_overrun         ! resp
   real                     :: rloss_overrun         ! resp
   real                     :: leaf_int_loss         ! resp
   real                     :: root_int_loss    ! resp

   real                     :: carbon13_balance         ! resp

   real                     :: flux_fact
   real                     :: flux_fact_inv
   real                     :: dtlsm_c_gain
   real                     :: dtlsm_c13_gain

   real                     :: gloss_via_lrr
   real                     :: gloss_via_gr

   real                     :: sloss_via_lrr
   real                     :: sloss_via_gr

   real                     :: lloss_via_lr
   real                     :: rloss_via_rr
   real                     :: lloss_via_lgr
   real                     :: rloss_via_rgr

   real                     :: lr_o_lrr
   real                     :: rr_o_lrr
   real                     :: lgr_o_gr
   real                     :: rgr_o_gr
   real                     :: sagr_o_gr
   real                     :: sbgr_o_gr
   real                     :: lgr_o_lgr_sagr
   real                     :: rgr_o_rgr_sbgr
   real                     :: sagr_o_lgr_sagr
   real                     :: sbgr_o_rgr_sbgr

   real                     :: f_debt_2_lrr
   real                     :: growth_resp
   real                     :: lrresp_residual


   real         , dimension(46) :: patch_check_vals         ! resp
   character(18), dimension(46) :: patch_check_labs       ! resp
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
      call get_storage_resp(cpatch,ico,tfact)
      
      cpatch%bstorage(ico)     = cpatch%bstorage(ico)                                      & 
                               - ( cpatch%leaf_storage_resp(ico)                           &
                                 + cpatch%root_storage_resp(ico)                           &
                                 + cpatch%sapa_storage_resp(ico)                           &
                                 + cpatch%sapb_storage_resp(ico)     ) *dtlsm_o_daysec

      cpatch%bstorage_c13(ico) = cpatch%bstorage_c13(ico)                                  &
                               - ( cpatch%leaf_storage_resp_c13(ico)                       &
                                 + cpatch%root_storage_resp_c13(ico)                       &
                                 + cpatch%sapa_storage_resp_c13(ico)                       &
                                 + cpatch%sapb_storage_resp_c13(ico) ) *dtlsm_o_daysec

      !------------------------------------------------------------------------------------!
      ! Standardize respiration units and get total respiration                            !
      !------------------------------------------------------------------------------------!
      leaf_resp = cpatch%leaf_respiration(ico) *flux_fact
      root_resp = cpatch%root_respiration(ico) *flux_fact
      
      growth_resp= cpatch%leaf_growth_resp(ico) + cpatch%root_growth_resp(ico)             &
                 + cpatch%sapa_growth_resp(ico) + cpatch%sapb_growth_resp(ico)

      resp_sum = growth_resp + leaf_resp + root_resp

      !------------------------------------------------------------------------------------!
      ! Find the proportion of total respiration each term comprises:                      !
      !------------------------------------------------------------------------------------!
      lr_o_lrr = max(min(hotc(leaf_resp ,leaf_resp+root_resp),1.0),0.0)
      rr_o_lrr = max(min(hotc(root_resp ,leaf_resp+root_resp),1.0),0.0)
      
      lgr_o_gr  = max(min(hotc(cpatch%leaf_growth_resp(ico),growth_resp),1.0),0.0)
      rgr_o_gr  = max(min(hotc(cpatch%root_growth_resp(ico),growth_resp),1.0),0.0)
      sagr_o_gr = max(min(hotc(cpatch%sapa_growth_resp(ico),growth_resp),1.0),0.0)
      sbgr_o_gr = max(min(hotc(cpatch%sapb_growth_resp(ico),growth_resp),1.0),0.0)
      
      
      !------------------------------------------------------------------------------------!
      ! Find the proportions to partition leaf-based and root-based resp.                  !
      !------------------------------------------------------------------------------------!
      
      dtlsm_c_gain   =  cpatch%gpp(ico)*flux_fact - leaf_resp - root_resp
      carbon_balance = (cpatch%gpp(ico)*flux_fact) - resp_sum
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
      gloss_via_lrr  = 0.0
      gloss_via_gr   = 0.0
      stor_loss  = 0.0
      sloss_via_lrr = 0.0
      sloss_via_gr  = 0.0
      leaf_loss  = 0.0
      root_loss  = 0.0
            
      lloss_via_lr  = 0.0
      rloss_via_rr  = 0.0
      lloss_via_lgr = 0.0
      rloss_via_rgr = 0.0

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
            f_lgr  = lgr_o_gr
            f_rgr  = rgr_o_gr
            f_sagr = sagr_o_gr
            f_sbgr = sbgr_o_gr
         end if
      else
         ! Leaf/Root loss will be split evenly between the fraction of leaf_respiration not
         ! already accounted for and the other above/below-ground growth resps.
         ! (Remember: abs(dtlsm_c_gain) is the total flux loss minus gpp)
         f_lr = lr_o_lrr *hotc(abs(dtlsm_c_gain),growth_resp + abs(dtlsm_c_gain))
         f_rr = rr_o_lrr *hotc(abs(dtlsm_c_gain),growth_resp + abs(dtlsm_c_gain))
        
         f_lgr  = hotc(cpatch%leaf_growth_resp(ico),growth_resp + abs(dtlsm_c_gain)) 
         f_rgr  = hotc(cpatch%root_growth_resp(ico),growth_resp + abs(dtlsm_c_gain)) 
         f_sagr = hotc(cpatch%sapa_growth_resp(ico),growth_resp + abs(dtlsm_c_gain)) 
         f_sbgr = hotc(cpatch%sapb_growth_resp(ico),growth_resp + abs(dtlsm_c_gain)) 

      end if

      if (carbon_balance > 0.0) then
         !---------------------------------------------------------------------------------!
         ! In get_c_xfers we will simply put carbon_balance into some pool(s), and since   !
         ! carbon_balance = gpp - resp_sum > 0.0                                           !
         !                                                                                 !
         ! We should set the signature of the respiration to the signature of gpp.         !
         ! Note that this also implies carbon13_balance > 0.0.                             !
         !---------------------------------------------------------------------------------!
         gpp_loss  = ratio_gpp *resp_sum

         gloss_via_lrr = ratio_gpp *(leaf_resp + root_resp)
         gloss_via_gr  = ratio_gpp *growth_resp
         
      !elseif (cpatch%phenology_status(ico) == 1 .and. (carbon_balance <= 0.0)) then
         !---------------------------------------------------------------------------------!
         ! "time_to_flush" will be .true. in get_c_xfers and we will move storage into     !
         ! plant pools despite having lost all gpp and some storage to respiration.        !
         !---------------------------------------------------------------------------------!
         !gpp_loss  = cpatch%gpp_c13(ico) *flux_fact
         !stor_loss = -1.0*carbon_balance *ratio_stor
         
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
            

            gloss_via_lrr  = gpp_loss *max(min(     (leaf_resp+root_resp)/(cpatch%gpp(ico)*flux_fact) ,1.0),0.0)
            gloss_via_gr   = gpp_loss *max(min((1.0-(leaf_resp+root_resp)/(cpatch%gpp(ico)*flux_fact)),1.0),0.0)
            sloss_via_lrr = 0.0
            sloss_via_gr  = stor_loss

            lloss_via_lr  = 0.0
            rloss_via_rr  = 0.0
            lloss_via_lgr = leaf_loss
            rloss_via_rgr = root_loss

         else
            gloss_via_lrr  = gpp_loss
            gloss_via_gr   = 0.0
            sloss_via_lrr = stor_loss *       f_debt_2_lrr
            sloss_via_gr  = stor_loss *(1.0 - f_debt_2_lrr)
           
            lloss_via_lr  = leaf_loss * hotc(           f_lr , (f_lr + f_sagr + f_lgr))
            rloss_via_rr  = root_loss * hotc(           f_rr , (f_rr + f_sbgr + f_rgr))
            lloss_via_lgr = leaf_loss * hotc((f_lgr + f_sagr), (f_lr + f_sagr + f_lgr))
            rloss_via_rgr = root_loss * hotc((f_rgr + f_sbgr), (f_rr + f_sbgr + f_rgr))
         end if 
      end if
      !------------------------------------------------------------------------------------!
      lgr_o_lgr_sagr  = max(min(hotc(cpatch%leaf_growth_resp(ico),(cpatch%leaf_growth_resp(ico)+cpatch%sapa_growth_resp(ico))),1.0),0.0)
      sagr_o_lgr_sagr = max(min(hotc(cpatch%sapa_growth_resp(ico),(cpatch%leaf_growth_resp(ico)+cpatch%sapa_growth_resp(ico))),1.0),0.0)
          
      rgr_o_rgr_sbgr  = max(min(hotc(cpatch%root_growth_resp(ico),(cpatch%root_growth_resp(ico)+cpatch%sapb_growth_resp(ico))),1.0),0.0)
      sbgr_o_rgr_sbgr = max(min(hotc(cpatch%sapb_growth_resp(ico),(cpatch%root_growth_resp(ico)+cpatch%sapb_growth_resp(ico))),1.0),0.0)


      !------------------------------------------------------------------------------------!
      ! Now we partition the total 13C loss sensibly among the respiration terms by...     !      
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! 3) Dividing up storage and gpp based resp evenly and dividing leaf loss and root   !
      !    loss evenly into leaf-based partitions and root-based partitions respectively.  ! 
      !------------------------------------------------------------------------------------!
      ! if dtlsm_c_gain > 0 then gpp covers leaf & root resp.
      ! if dtlsm_c_gain < 0 then gpp does not cover leaf and root resp
      
      leaf_resp_c13 = (gloss_via_lrr + sloss_via_lrr)*lr_o_lrr + lloss_via_lr
      root_resp_c13 = (gloss_via_lrr + sloss_via_lrr)*rr_o_lrr + rloss_via_rr

      cpatch%leaf_respiration_c13(ico) = leaf_resp_c13 * flux_fact_inv
      cpatch%root_respiration_c13(ico) = root_resp_c13 * flux_fact_inv

      cpatch%leaf_growth_resp_c13(ico) = (gloss_via_gr + sloss_via_gr)*lgr_o_gr  + lloss_via_lgr *lgr_o_lgr_sagr
      cpatch%root_growth_resp_c13(ico) = (gloss_via_gr + sloss_via_gr)*rgr_o_gr  + rloss_via_rgr *rgr_o_rgr_sbgr
      cpatch%sapa_growth_resp_c13(ico) = (gloss_via_gr + sloss_via_gr)*sagr_o_gr + lloss_via_lgr *sagr_o_lgr_sagr
      cpatch%sapb_growth_resp_c13(ico) = (gloss_via_gr + sloss_via_gr)*sbgr_o_gr + rloss_via_rgr *sbgr_o_rgr_sbgr
     
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
                       - leaf_resp_c13                                                     &
                       - root_resp_c13
      carbon13_balance = cpatch%gpp_c13(ico) * flux_fact                                   &
                       - leaf_resp_c13                                                     &
                       - root_resp_c13                                                     &
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
      patch_check_vals = &
     (/carbon_balance, carbon13_balance,cpatch%gpp(ico)*flux_fact,cpatch%gpp_c13(ico)*flux_fact &
      ,     stor_loss,         gpp_loss,                leaf_loss,                    root_loss &
      ,        f_lr,           f_rr,                    f_lgr,                        f_rgr &
      ,        f_sagr,           f_sbgr                 &
      , leaf_int_loss,    root_int_loss,            lloss_overrun,                rloss_overrun &
      ,     bloss_max, real(cpatch%phenology_status(ico))                                       &
      , htIsoDelta(abs(dtlsm_c13_gain)  ,abs(dtlsm_c_gain))     &
      , htIsoDelta(abs(carbon13_balance),abs(carbon_balance)), dtlsm_c_gain, dtlsm_c13_gain     &
      , gloss_via_lrr    ,         gloss_via_gr,               sloss_via_lrr,                    sloss_via_gr &
      ,       lr_o_lrr,          rr_o_lrr,                   lgr_o_gr,                       rgr_o_gr &
      ,       sagr_o_gr,          sbgr_o_gr,           lloss_via_lr,               rloss_via_rr &
      ,lloss_via_lgr, rloss_via_rgr,              growth_resp,          f_debt_2_lrr &
      ,             dtlsm_c_gain,               dtlsm_c13_gain &
      ,lgr_o_lgr_sagr, sagr_o_lgr_sagr,rgr_o_rgr_sbgr,sbgr_o_rgr_sbgr/)

      patch_check_labs = &
     (/'    carbon_balance','  carbon13_balance','               gpp','           gpp_c13'&
      ,'         stor_loss','          gpp_loss','         leaf_loss','         root_loss'&
      ,'   gen_frac_leaf_r','   gen_frac_root_r','  gen_frac_leaf_gr','  gen_frac_root_gr'&
      ,'  gen_frac_sapa_gr','  gen_frac_sapb_gr'&
      ,'     leaf_int_loss','     root_int_loss','     lloss_overrun','     rloss_overrun'&
      ,'         bloss_max','  phenology_status','          dcg_d13C','   carbon_bal_d13C'&
      ,'      dtlsm_c_gain','    dtlsm_c13_gain'                                          &
      ,'         gloss_via_lrr','          gloss_via_gr','        sloss_via_lrr','         sloss_via_gr'&
      ,'           lr_o_lrr','           rr_o_lrr','            lgr_o_gr','            rgr_o_gr'&
      ,'           sagr_o_gr','           sbgr_o_gr','    lloss_via_lr','    rloss_via_rr'&
      ,'   lloss_via_lgr','   rloss_via_rgr','       growth_resp','    r_c_debt_2_lrr'&
      ,'      daily_c_gain','    daily_c13_gain'&
      ,' lgr_o_lgr_sagr','sagr_o_lgr_sagr',' rgr_o_rgr_sbgr','sbgr_o_rgr_sbgr'/)

 
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
   real                     :: ratio_gpp                  ! GPP
   real                     :: ratio_stor                 ! storage
   real                     :: ratio_resp                 ! resp
   real                     :: dtlsm_o_frqsum             ! resp
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
