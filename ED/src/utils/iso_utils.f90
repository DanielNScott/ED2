module iso_utils
implicit none
contains

!==========================================================================================!
!     This fn returns the ratio 13C/12C of carbon fixed under specified conditions, in     !
!     accordance with the farquhar 1983 photosynthesis fractionation model                 !
!------------------------------------------------------------------------------------------!
! Model: Farquhar et al. 1982,1991 ; Implementation sim. to: Ehleringer, Law et al. 2006 } !
!  Note: Proportionality of partial pressures with concentrations == Henry's Law           !
!  Note: Leaf respiration fractionation (sixth term in F1982) is NOT included here since   !
!        leaf resp. is subtracted from photosynth. output in the model broadly             !
!------------------------------------------------------------------------------------------!
! Details:
!------------------------------------------------------------------------------------------!
real function photo_h2tc(d13C_atm,can_co2,fs_open,lsfc_co2_open,lsfc_co2_closed            &
                        ,lint_co2_open,lint_co2_closed)
   use isotopes, only  : R_std  ! ! intent(in)
   implicit none

   !------ Arguments. ------------------------------------------------------------------------!
   real(kind=4), intent(in)   :: d13C_atm          ! delta 13C for atm           [ per mil ]
   real(kind=4), intent(in)   :: fs_open           ! Open proportion of stom.    [ dim.less]
   real(kind=4), intent(in)   :: can_co2           ! Canopy air CO2              [ µmol/mol]

   real(kind=4), intent(in)   :: lsfc_co2_open     ! Leaf sfc. CO2    (op.)      [ µmol/mol]
   real(kind=4), intent(in)   :: lsfc_co2_closed   ! Leaf sfc. CO2    (cl.)      [ µmol/mol]
   real(kind=4), intent(in)   :: lint_co2_open     ! Intercell. CO2   (op.)      [ µmol/mol]
   real(kind=4), intent(in)   :: lint_co2_closed   ! Intercell. CO2   (cl.)      [ µmol/mol]

   !------ Local vars. -----------------------------------------------------------------------!
   ! Note: All concentrations are divided by other concentrations, hence units cancel in eq   !
   !       for f. Local variables are derived from stomatal props.                            !
   !------------------------------------------------------------------------------------------!
   real(kind=4)            :: d  ! d = (d13C_atm - d13C_assim)/(1+d13C_atm/1000) [ per mil ]
   real(kind=4)            :: ca ! CO2 conc in atm                               [ µmol/mol]
   real(kind=4)            :: cs ! CO2 conc on leaf surface                      [ µmol/mol]
   real(kind=4)            :: ci ! CO2 conc intercellular                        [ µmol/mol]
   real(kind=4)            :: cc ! CO2 conc chloroplast                          [ µmol/mol]
   
   real(kind=4)            :: d13C_assim  ! d13C_atm - (1 + d13C_atm/1000) *d    [ per mil ]
   real                    :: R_assim     ! Heavy to light carbon                [         ]

   !----- Fractionation constants, Ehlerigner, Law et al. 2006 -------------------------------!
   real(kind=4), parameter :: ab = 2.9 ! diff across bdry layer                  [ per mil ]
   real(kind=4), parameter :: a  = 4.4 ! diff across stomatal cavity             [ per mil ]
   real(kind=4), parameter :: es = 1.1 ! diff into solution                      [ per mil ]
   real(kind=4), parameter :: as = 0.7 ! diff w/in solution                      [ per mil ]
   real(kind=4), parameter :: b  = 28.2! by rubisco, assuming 5% PEP fix.        [ per mil ]
   
   ca = can_co2
   cs = (fs_open * lsfc_co2_open) + ( (1.0-fs_open) * lsfc_co2_closed )
   ci = (fs_open * lint_co2_open) + ( (1.0-fs_open) * lint_co2_closed )

   !-------------------------------------------------------------------!
   !- RIGHT NOW THIS IS TOTALLY ARBITRARY, NEED A DEFENSIBLE ESTIMATE -!
   !-------------------------------------------------------------------!
   cc = 0.995 * ci  
   !-------------------------------------------------------------------!

   !f = ( ab*(ca - cs)  +  a*(cs - ci)  +  (es + as)*(ci-cc)  +  b*cc )/ca
   d  = a + (27.0 - a)*ci/ca
   d13C_assim = d13C_atm - (1.0 + d13C_atm/1000.0) *d
   
   R_assim = (d13C_assim/1000.0 + 1.0)*R_std
   photo_h2tc = R_assim/(R_assim + 1.0)
   return
end function photo_h2tc
!==========================================================================================!


!==========================================================================================!
real function hotc(heavy,total)
   use consts_coms, only : tiny_num    ! intent(in)
   implicit none
   !------ Arguments ----------------------------------------------------------------------!
   real              :: heavy       ! C-13
   real              :: total       ! C-12 + C-13
   !---------------------------------------------------------------------------------------!

   ! This function computes the safe ratio of heavy to total carbon.
   if (total > tiny_num) then
      hotc = heavy/total
   else
      hotc = 0.0
   end if
   
   return
end function hotc
!==========================================================================================!

!==========================================================================================!
real (kind=8) function sdiv8(heavy,total)
   use consts_coms, only : tiny_num    ! intent(in)
   implicit none
   !------ Arguments ----------------------------------------------------------------------!
   real              :: heavy       ! C-13
   real              :: total       ! C-12 + C-13
   !---------------------------------------------------------------------------------------!

   ! This function computes the safe ratio of heavy to total carbon.
   if (total > tiny_num) then
      sdiv8 = dble(heavy)/dble(total)
   else
      sdiv8 = dble(0.0)
   end if
   
   return
end function sdiv8
!==========================================================================================!

!==========================================================================================!
real(kind=8) function resp_h2tc(rtype,heavy,total)
   use isotopes,        only : R_std      & ! intent(in)
                             , iso_lrf    & ! intent(in)
                             , iso_rrf    & ! intent(in)
                             , iso_grf    & ! intent(in)
                             , iso_strf   & ! intent(in)
                             , iso_vlrf   & ! intent(in)
                             , iso_hrf    ! ! intent(in)
   implicit none

   !------ Arguments. ------------------------------------------------------------------------!
   character(len=*)  :: rtype       ! respiration type
   real              :: heavy       ! 13C
   real              :: total       ! 12C + 13C
   
   !------ Local Var. ------------------------------------------------------------------------!
   real(kind=8)               :: frac     ! permil "apparent fractionation" by respiring item
   real(kind=8)               :: heavy8   ! Double precision heavy carbon
   real(kind=8)               :: total8   ! Double precision total carbon
   real(kind=8)               :: h2lc     ! heavy to light carbon
   real(kind=8)               :: rh2lc    ! respired heavy to light carbon
   !------------------------------------------------------------------------------------------!

   !------------------------------------------------------------------------------------------!
   ! Select fractionation associated with each respiration type.                              !
   ! 1: Leaf   2: Root  3: Growth  4: Storage   5: Virtual Leaf   6: Heterotrophic            !
   !------------------------------------------------------------------------------------------!
   select case(rtype)
   case('leaf')
      frac = dble(iso_lrf)
   case('root')
      frac = dble(iso_rrf)
   case('growth')
      frac = dble(iso_grf)
   case('stor')
      frac = dble(iso_strf)
   case('het')
      frac = dble(iso_hrf)
   case default
      write(*,*) 'Subroutine resp_h2tc has been misused.'
      stop
   end select

   heavy8 = dble(heavy)
   total8 = dble(total)
   
   if (total > tiny(1.0)) then
      !--- Maximize accuracy when fractionation is set to 0 ----------------------------------!
      if (abs(frac) < tiny(1.0)) then
         resp_h2tc = heavy8/total8
         return
      end if
      rh2lc        = (dble(frac + htIsoDelta(heavy,total))/1000.d0 + 1.d0)*R_std
      resp_h2tc    = rh2lc/(rh2lc + 1.d0)
   else
      resp_h2tc = 0.0
   end if
   
   return
end function resp_h2tc
!==========================================================================================!

!==========================================================================================!
real function htIsoDelta(heavy,total)
   use isotopes, only  : R_std  ! ! intent(in)
   implicit none
   
   !---------------------------------------------------------------------------------------!
   ! Returns VPD referenced delta 13C value
   !------ Arguments. ---------------------------------------------------------------------!
   real     :: heavy       ! Heavy Carbon
   real     :: total       ! Total Carbon
   !------ Local (Double Precision) -------------------------------------------------------!
   real(kind=8)     :: heavy8       ! Heavy Carbon
   real(kind=8)     :: total8       ! Total Carbon
   
   heavy8 = dble(heavy)
   total8 = dble(total)
      
   if (abs(total) > tiny(1.0)) then
      if (sign(1.0,heavy) /= sign(1.0,total)) then
         htIsoDelta = huge(1.0)
      else
         htIsoDelta = sngl((abs(heavy8)/(abs(total8)-abs(heavy8)))/R_std - 1.d0) * 1000.0
      end if
   else
      if (abs(heavy) < tiny(1.0)) then
         htIsoDelta = 0.0
      else
         htIsoDelta = huge(1.0)
      end if      
   end if
   
   return
end function htIsoDelta
!==========================================================================================!   

!==========================================================================================!
real function htIsoDelta8(heavy,total)
   use isotopes, only  : R_std  ! ! intent(in)
   implicit none
   
   !---------------------------------------------------------------------------------------!
   ! Returns VPD referenced delta 13C value
   !------ Arguments. ---------------------------------------------------------------------!
   real(kind=8)     :: heavy       ! Heavy Carbon
   real(kind=8)     :: total       ! Total Carbon
   
   if (total > tiny(1.0)) then
      htIsoDelta8 = ((heavy/(total-heavy))/R_std - 1.0) * 1000.0
   else
      if (heavy < tiny(1.0)) then
         htIsoDelta8 = 0.0
      else
         htIsoDelta8 = huge(1.0)
      end if      
   end if
   
   return
end function htIsoDelta8
!==========================================================================================!   

!==========================================================================================!
real function d13C2Ratio(delta)
   use isotopes, only  : R_std  ! ! intent(in)
   implicit none
   
   !---------------------------------------------------------------------------------------!
   ! Returns ratio of heavy to total carbon for a given delta value.
   !------ Arguments. ---------------------------------------------------------------------!
   real     :: delta       ! Delta value
   real     :: R_hl        ! Ratio of heavy to light carbon
   
   R_hl        = (delta/1000.0 + 1.0) * R_std
   d13C2Ratio  = R_hl/(R_hl + 1.0)
   
   return
end function d13C2Ratio
!==========================================================================================!   

!==========================================================================================!
real function d13C2Ratio8(delta)
   use isotopes, only  : R_std  ! ! intent(in)
   implicit none
   
   !---------------------------------------------------------------------------------------!
   ! Returns ratio of heavy to total carbon for a given delta value.
   !------ Arguments. ---------------------------------------------------------------------!
   real(kind=8)     :: delta       ! Delta value
   real(kind=8)     :: R_hl        ! Ratio of heavy to light carbon
   
   R_hl        = (delta/1000.0 + 1.0) * R_std
   d13C2Ratio8 = R_hl/(R_hl + 1.0)
   
   return
end function d13C2Ratio8
!==========================================================================================!



!==========================================================================================!
!==========================================================================================!
!     This sub-routine does a pretty comprehensive sanity check on C-13 variables.         !
!------------------------------------------------------------------------------------------!
subroutine check_patch_c13(cpatch,ico,call_loc,fname,aux_vals,aux_labs,aux_pair)
   use ed_max_dims    , only : str_len            ! ! intent(in)
   use ed_state_vars  , only : patchtype          ! ! structure
   use isotopes       , only : c13af              & ! intent(in)
                             , c_alloc_flg        ! ! intent(in)
   use isotopes       , only : larprop            ! ! intent(in)
   use ed_misc_coms   , only : current_time       ! ! intent(in)
   use consts_coms    , only : umol_2_kgC         ! ! intent(in)
   use ed_misc_coms   , only : dtlsm              ! ! intent(in)
  implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(patchtype)           , intent(in) :: cpatch      ! Current patch
   integer                   , intent(in) :: ico         ! Current cohort number
   character(*)              , intent(in) :: call_loc    ! What routine called this check?
   character(*)              , intent(in) :: fname       ! What file is that routine in?
   real         , dimension(:), optional, intent(in) :: aux_vals
   character(18), dimension(:), optional, intent(in) :: aux_labs
   integer      , dimension(:), optional, intent(in) :: aux_pair
   !----- Local variables. ----------------------------------------------------------------!
   logical        :: error_found = .false.               ! Is there a problem?
   logical        :: check_delta = .true.               ! Is there a problem?
   character(40)  :: reason                              ! Error or warning diagnosis reason
   character(7)   :: Cfmt                                ! Character format, for strings
   character(9)   :: Rfmt                                ! Real format, for reals.
   real           :: gpp_delta                           ! Delta C-13 vals for variables.
   real           :: lr_delta                            ! ...
   real           :: rr_delta                            ! ...
   real           :: lg_delta                            ! ...
   real           :: rg_delta                            ! ...
   real           :: sag_delta                           ! ...
   real           :: sbg_delta                           ! ...
   real           :: ls_delta                            ! ...
   real           :: rs_delta                            ! ...
   real           :: sas_delta                           ! ...
   real           :: sbs_delta                           ! ...
   real           :: leaf_delta                          ! ...
   real           :: root_delta                          ! ...
   real           :: bsa_delta                           ! ...
   real           :: bsb_delta                           ! ...
   real           :: bst_delta                           ! ...
   real           :: non_stor_resp                       ! Sum of non storage resps
   real           :: nsr_c13                             ! C-13 content of above
   logical        :: valid                               ! Valid C-13 to C ratio logical
   logical        :: gpp_override                               ! Valid C-13 to C ratio logical
   integer        :: aux_size                            ! Length of auxiliary input
   real           :: leaf_resp                             ! C-13 content of above
   real           :: root_resp                            ! C-13 content of above
   real           :: leaf_resp_c13                             ! C-13 content of above
   real           :: root_resp_c13                             ! C-13 content of above
   integer        :: loop_ind                            ! Aux var printout loop index.
   integer        :: loop_ind2                            ! Aux var printout loop index.
   real           :: flux_fact                       ! Aux var printout loop index.
   real           :: loop_delta                       ! Aux var printout loop index.
   !---------------------------------------------------------------------------------------!
   Cfmt = '(11A18)'
   Rfmt = '(11E18.2)'

   flux_fact = umol_2_kgC*dtlsm/cpatch%nplant(ico)
   !------------------------------------------------------------------------------!
   ! Perform a great big sanity check. Today_XXX vars are included to potentially !
   ! detect problems generated elsewhere in the code. Clearly it is the case in   !
   ! this routine that if a today var is messed up it is because it's sub-daily   !
   ! analogue is to blame.                                                        !
   !------------------------------------------------------------------------------!
   ! Assimilated C-13 should not exceed assimilated C.
   reason = ''

   ! If gpp has a crazy delta, we should allow the model to proceed with
   ! matching fluxes and even matching pool values if small enough. If we have a
   ! real problem, this will manifest at a time when gpp has a normal delta too.  !
   gpp_override = .false.

   call check_c13(cpatch%gpp_c13(ico) &
                 ,cpatch%gpp    (ico) &
                 ,check_delta,gpp_delta,valid,'gpp',reason)
   if (not(valid)) then
      gpp_override = .true.
   end if
   
   call check_c13(cpatch%leaf_respiration_c13(ico) &
                 ,cpatch%leaf_respiration    (ico) &
                 ,check_delta,lr_delta,valid,'leaf_resp',reason)
   
   call check_c13(cpatch%root_respiration_c13(ico) &
                 ,cpatch%root_respiration    (ico) &
                 ,check_delta,rr_delta,valid,'root_resp',reason)
   
   call check_c13(cpatch%leaf_growth_resp_c13(ico) &
                 ,cpatch%leaf_growth_resp    (ico) &
                 ,check_delta,lg_delta,valid,'leaf_grow_resp',reason)
                 
   call check_c13(cpatch%root_growth_resp_c13(ico) &
                 ,cpatch%root_growth_resp    (ico) &
                 ,check_delta,rg_delta,valid,'root_grow_resp',reason)
   
   call check_c13(cpatch%sapa_growth_resp_c13(ico) &
                 ,cpatch%sapa_growth_resp    (ico) &
                 ,check_delta,sag_delta,valid,'sapa_grow_resp',reason)
   
   call check_c13(cpatch%sapb_growth_resp_c13(ico) &
                 ,cpatch%sapb_growth_resp    (ico) &
                 ,check_delta,sbg_delta,valid,'sapb_grow_resp',reason)
   
   call check_c13(cpatch%leaf_storage_resp_c13(ico) &
                 ,cpatch%leaf_storage_resp    (ico) &
                 ,check_delta,ls_delta,valid,'leaf_stor_resp',reason)
   
   call check_c13(cpatch%root_storage_resp_c13(ico) &
                 ,cpatch%root_storage_resp    (ico) &
                 ,check_delta,rs_delta,valid,'root_stor_resp',reason)
   
   call check_c13(cpatch%sapa_storage_resp_c13(ico) &
                 ,cpatch%sapa_storage_resp    (ico) &
                 ,check_delta,sas_delta,valid,'sapa_stor_resp',reason)
   
   call check_c13(cpatch%sapb_storage_resp_c13(ico) &
                 ,cpatch%sapb_storage_resp    (ico) &
                 ,check_delta,sbs_delta,valid,'sapb_stor_resp',reason)
   
   call check_c13(cpatch%bleaf_c13(ico) &
                 ,cpatch%bleaf    (ico) &
                 ,check_delta,leaf_delta,valid,'bleaf',reason)
   
   call check_c13(cpatch%broot_c13(ico) &
                 ,cpatch%broot    (ico) &
                 ,check_delta,root_delta,valid,'broot',reason)
   
   call check_c13(cpatch%bsapwooda_c13(ico) &
                 ,cpatch%bsapwooda    (ico) &
                 ,check_delta,bsa_delta,valid,'sapwooda',reason)
                 
   call check_c13(cpatch%bsapwoodb_c13(ico) &
                 ,cpatch%bsapwoodb    (ico) &
                 ,check_delta,bsb_delta,valid,'sapwoodb',reason)
   
   call check_c13(cpatch%bstorage_c13(ico) &
                 ,cpatch%bstorage    (ico) &
                 ,check_delta,bst_delta,valid,'storage',reason)

   leaf_resp = cpatch%leaf_respiration(ico) *flux_fact
   root_resp = cpatch%root_respiration(ico) *flux_fact

   leaf_resp_c13 = cpatch%leaf_respiration_c13(ico) *flux_fact
   root_resp_c13 = cpatch%root_respiration_c13(ico) *flux_fact

   non_stor_resp = cpatch%leaf_growth_resp(ico) + cpatch%root_growth_resp(ico)              &
                 + cpatch%sapa_growth_resp(ico) + cpatch%sapb_growth_resp(ico)              &
                 + leaf_resp                    + root_resp   

   nsr_c13 = cpatch%leaf_growth_resp_c13(ico) + cpatch%root_growth_resp_c13(ico)      &
           + cpatch%leaf_growth_resp_c13(ico) + cpatch%root_growth_resp_c13(ico)      &
           + leaf_resp_c13                    + root_resp_c13
   
   if (present(aux_pair)) then
      loop_ind = 1
      aux_size = size(aux_vals)
      outerloop: do loop_ind = 1,aux_size
         if (aux_pair(loop_ind) /= 0 .and. loop_ind < aux_size) then
            do loop_ind2 = loop_ind+1,aux_size
               if (aux_pair(loop_ind2) == aux_pair(loop_ind)) then
                  !write(*,*) 'Checking ', aux_labs(loop_ind2), ' ', aux_labs(loop_ind)
                  !write(*,*) aux_vals(loop_ind2), aux_vals(loop_ind)
                  call check_c13(aux_vals(loop_ind2) &
                                ,aux_vals(loop_ind)  &
                                ,check_delta,loop_delta,valid,aux_labs(loop_ind),reason)
                  !valid = .false.
                  exit outerloop
               end if
            end do
         end if
      end do outerloop
   end if
   
   if (not(valid) .and. not(gpp_override)) then
      write(*,*) '======================================================================='
      write(*,*) ' ', reason
      write(*,*) ' C-13 sanity check error in ', call_loc, '!'
      write(*,*) '======================================================================='      
      write(*,'(A13,4I4)') ' Model Time: ', current_time%date, current_time%hour         &
                                          , current_time%min , current_time%sec 
      write(*,*) ' cohort, pft, larprop : ', ico, cpatch%pft(ico), larprop
      write(*,*) ' '
      write(*,*) '-----------------------------------------------------------------------'
      write(*,*) 'Flux Diagnostics. Row 1: Totals. Row 2: Heavy C. Row 3: d13C.'
      write(*,*) '-----------------------------------------------------------------------'
      write(*,*) 'Units: kgC/pl/yr'
      write(*,'(4A18)')   '               GPP'          &
                         ,'  LEAF_RESPIRATION'        &
                         ,'  ROOT_RESPIRATION'        &
                         ,'        **RESP_SUM'
      write(*,'(4ES18.8)')  cpatch%gpp(ico)*flux_fact                &
                          ,leaf_resp   &
                          ,root_resp   &
                          ,non_stor_resp
      write(*,'(4ES18.8)')  cpatch%gpp_c13(ico)*flux_fact             &
                          ,leaf_resp_c13  &
                          ,root_resp_c13  &
                          ,nsr_c13
      write(*,'(4ES18.8)') gpp_delta,lr_delta,rr_delta, htIsoDelta(nsr_c13,non_stor_resp)
      write(*,*) ' '
      write(*,*) 'Units: kgC/pl/yr'
      write(*,'(5A18)') '  LEAF_GROWTH_RESP','  ROOT_GROWTH_RESP','  SAPA_GROWTH_RESP' &
                       ,'  SAPB_GROWTH_RESP'
      write(*,'(5ES18.8)') cpatch%leaf_growth_resp(ico) ,cpatch%root_growth_resp(ico)       &
                       ,cpatch%sapa_growth_resp(ico) ,cpatch%sapb_growth_resp(ico)
      write(*,'(5ES18.8)') cpatch%leaf_growth_resp_c13(ico)    &
                       ,cpatch%root_growth_resp_c13(ico)    &
                       ,cpatch%sapa_growth_resp_c13(ico)    &
                       ,cpatch%sapb_growth_resp_c13(ico)
      write(*,'(5ES18.8)') lg_delta,rg_delta,sag_delta,sbg_delta
      write(*,*) ' '
      write(*,*) 'Units: kgC/pl/yr'
      write(*,'(5A18)') ' LEAF_STORAGE_RESP',' ROOT_STORAGE_RESP',' SAPA_STORAGE_RESP' &
                         ,' SAPB_STORAGE_RESP'
      write(*,'(5ES18.8)') cpatch%leaf_storage_resp(ico)       &
                       ,cpatch%root_storage_resp(ico)       &
                       ,cpatch%sapa_storage_resp(ico)       &
                       ,cpatch%sapb_storage_resp(ico)
      write(*,'(5ES18.8)') cpatch%leaf_storage_resp_c13(ico)   &
                       ,cpatch%root_storage_resp_c13(ico)   &
                       ,cpatch%sapa_storage_resp_c13(ico)   &
                       ,cpatch%sapb_storage_resp_c13(ico)
      write(*,'(5ES18.8)') ls_delta,rs_delta,sas_delta,sbs_delta

      write(*,*) ' '
      write(*,*) ' **RESP_SUM is the sum of non-storage respiration terms.'
      write(*,*) '-----------------------------------------------------------------------'
      write(*,*) 'Pool Diagnostics. Row 1: Totals. Row 2: Heavy C. Row 3: d13C.'
      write(*,*) '-----------------------------------------------------------------------'
      write(*,*) 'Units: kgC/pl'
      write(*,'(5A14)') '         BLEAF','         BROOT','     BSAPWOODA'               &
                        ,'    BSAPWOODB','      BSTORAGE'
      write(*,'(5ES14.6)') cpatch%bleaf(ico),cpatch%broot(ico),cpatch%bsapwooda(ico)      &
                         ,cpatch%bsapwoodb(ico), cpatch%bstorage(ico)

      write(*,'(5ES14.6)') cpatch%bleaf_c13(ico),cpatch%broot_c13(ico)                    &
                         ,cpatch%bsapwooda_c13(ico),cpatch%bsapwoodb_c13(ico)            &
                         ,cpatch%bstorage_c13(ico)

      write(*,'(5ES14.6)') leaf_delta,root_delta,bsa_delta,bsb_delta,bst_delta
      write(*,*) ' '
      write(*,*) '-----------------------------------------------------------------------'
      if (present(aux_vals)) then
         aux_size = size(aux_vals)
         write(*,*) 'Auxiliary Diagnostics (see labels). Row format as above.'
         write(*,*) '-----------------------------------------------------------------------'
         write(*,'(4A18)')    aux_labs(1:min(4,aux_size))
         write(*,'(4ES18.8)') aux_vals(1:min(4,aux_size))
         loop_ind = 5;
         do while (aux_size > loop_ind .and. loop_ind < 100)
            write(*,*) ' '
            write(*,'(4A18)')    aux_labs(loop_ind:min(loop_ind+3,aux_size))
            write(*,'(4ES18.8)') aux_vals(loop_ind:min(loop_ind+3,aux_size))
            loop_ind = loop_ind + 4
         end do
         write(*,*) '-----------------------------------------------------------------------'
      end if
      !if (present(assim_h2tc)) then
      !   write(*,*) '-----------------------------------------------------------------------'
      !   write(*,*) ' Pieces of leaf C-13 respiration computation'
      !   write(*,Cfmt) 'R Leaf Tissue', 'R Leaf Assim',                        &
      !                 'Leaf C-13/C'  , '"leaf_h2tc" '
      !   write(*,Rfmt) cpatch%leaf_respiration(ico) - cpatch%lassim_resp(ico), &
      !                 cpatch%lassim_resp     (ico)                          , &
      !                 cpatch%bleaf_c13       (ico) / cpatch%bleaf      (ico), &
      !                 leaf_h2tc
      !   write(*,*)
      !   write(*,*) ' Leaf assimilate C-13/C ratio and recalculation from GPP and GPP_C13'
      !   write(*,*) ' If this msg is not from canopy_photosynthesis, ignore them.        '
      !   write(*,*) ' assim_h2tc : ', assim_h2tc, cpatch%gpp_c13(ico) / cpatch%gpp(ico)
      !   write(*,*) '-----------------------------------------------------------------------'
      !end if
      call fatal_error(reason,call_loc,fname)
   end if
   
   
   return
end subroutine check_patch_c13
!==========================================================================================!
!==========================================================================================!



!==========================================================================================!
!==========================================================================================!
!     This sub-routine does a pretty comprehensive sanity check on C-13 variables.         !
!------------------------------------------------------------------------------------------!
subroutine check_site_c13(csite,ipa,call_loc,fname,aux_vals,aux_labs,aux_pair)
   use ed_max_dims    , only : str_len            ! ! intent(in)
   use ed_state_vars  , only : sitetype           ! ! structure
   use isotopes       , only : c13af              & ! intent(in)
                             , larprop            ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype)            , intent(in) :: csite       ! Current site
   integer                   , intent(in) :: ipa         ! Current patch number
   character(*)              , intent(in) :: call_loc    ! What routine called this check?
   character(*)              , intent(in) :: fname       ! What file is that routine in?
   real         , dimension(:), optional, intent(in) :: aux_vals
   character(18), dimension(:), optional, intent(in) :: aux_labs
   integer      , dimension(:), optional, intent(in) :: aux_pair
   !----- Local variables. ----------------------------------------------------------------!
   logical        :: error_found = .false.               ! Is there a problem?
   logical        :: check_delta = .true.               ! Is there a problem?
   character(40)  :: reason                              ! Error or warning diagnosis reason
   character(6)   :: Cfmt                                ! Character format, for strings
   character(8)   :: Rfmt                                ! Real format, for reals.
   real           :: fsc_delta                           ! Delta C-13 vals for variables.
   real           :: ssc_delta                           ! ...
   real           :: stsc_delta                          ! ...
   real           :: stsl_delta                          ! ...
   real           :: can_co2_delta                       ! ...
   real           :: rh_delta                            ! ...
   real           :: cwd_delta                           ! ...
   logical        :: valid                               ! Valid C-13 to C ratio logical
   integer        :: aux_size                            ! Aux var printout loop index.
   integer        :: loop_ind                            ! Aux var printout loop index.
   integer        :: loop_ind2                            ! Aux var printout loop index.
   real           :: flux_fact                       ! Aux var printout loop index.
   real           :: loop_delta                       ! Aux var printout loop index.
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   Cfmt = '(4A18)'
   Rfmt = '(4E18.2)'

   reason = ''
   call check_c13(csite%fast_soil_c13(ipa) &
                 ,csite%fast_soil_c  (ipa) &
                 ,check_delta,fsc_delta,valid,'fast_soil_c',reason)
   
   call check_c13(csite%slow_soil_c13(ipa) &
                 ,csite%slow_soil_c  (ipa) &
                 ,check_delta,ssc_delta,valid,'slow_soil_c',reason)
   
   call check_c13(csite%structural_soil_c13(ipa) &
                 ,csite%structural_soil_c  (ipa) &
                 ,check_delta,stsc_delta,valid,'struc_soil_c',reason)
   
   call check_c13(csite%structural_soil_L_c13(ipa) &
                 ,csite%structural_soil_L    (ipa) &
                 ,check_delta,stsl_delta,valid,'struc_soil_L',reason)
   
   call check_c13(csite%can_co2_c13(ipa) &
                 ,csite%can_co2    (ipa) &
                 ,check_delta,can_co2_delta,valid,'can_co2',reason)
   
   call check_c13(csite%rh_c13(ipa) &
                 ,csite%rh    (ipa) &
                 ,check_delta,rh_delta,valid,'rh',reason)
   
   call check_c13(csite%cwd_rh_c13(ipa) &
                 ,csite%cwd_rh    (ipa) &
                 ,check_delta,cwd_delta,valid,'cwd_rh',reason)
   
   if (present(aux_pair)) then
      loop_ind = 1
      aux_size = size(aux_vals)
      outerloop: do loop_ind = 1,aux_size
         if (aux_pair(loop_ind) /= 0 .and. loop_ind < aux_size) then
            do loop_ind2 = loop_ind+1,aux_size
               if (aux_pair(loop_ind2) == aux_pair(loop_ind)) then
                  !write(*,*) 'Checking ', aux_labs(loop_ind2), ' ', aux_labs(loop_ind)
                  !write(*,*) aux_vals(loop_ind2), aux_vals(loop_ind)
                  call check_c13(aux_vals(loop_ind2) &
                                ,aux_vals(loop_ind)  &
                                ,check_delta,loop_delta,valid,aux_labs(loop_ind),reason)
                  !valid = .false.
                  exit outerloop
               end if
            end do
         end if
      end do outerloop
   end if
   
   if (not(valid)) then
      write(*,*) '======================================================================='
      write(*,*) ' C-13 sanity check error in ', call_loc, '!'
      write(*,*) reason
      write(*,*) '======================================================================='
      write(*,*) ' patch : ', ipa
      write(*,*) ' '
      write(*,*) '-----------------------------------------------------------------------'
      write(*,*) 'Pool Diagnostics. Row 1: Totals. Row 2: Heavy C. Row 3: d13C.'
      write(*,*) '-----------------------------------------------------------------------'
      write(*,Cfmt) '       FAST_SOIL_C','       SLOW_SOIL_C',' STRUCTURAL_SOIL_C'       &
                   ,' STRUCTURAL_SOIL_L'
      write(*,Rfmt) csite%fast_soil_C(ipa),csite%slow_soil_C(ipa)                        &
                   ,csite%structural_soil_C(ipa),csite%structural_soil_L(ipa)
                   
      write(*,Rfmt) csite%fast_soil_c13(ipa),csite%slow_soil_c13(ipa)                    &
                   ,csite%structural_soil_c13(ipa),csite%structural_soil_L_c13(ipa)
                   
      write(*,Rfmt) fsc_delta, ssc_delta, stsc_delta, stsl_delta
      write(*,*)    ' '
      write(*,Cfmt) '           CAN_CO2','               RH','             CWD_RH'
      write(*,Rfmt) csite%can_co2(ipa),csite%rh(ipa),csite%cwd_rh(ipa)
      write(*,Rfmt) csite%can_co2_c13(ipa),csite%rh_c13(ipa),csite%cwd_rh_c13(ipa)
      write(*,Rfmt) can_co2_delta,rh_delta,cwd_delta
      write(*,*) '-----------------------------------------------------------------------'
      if (present(aux_vals)) then
         aux_size = size(aux_vals)
         write(*,*) 'Auxiliary Diagnostics (see labels). Row format as above.'
         write(*,*) '-----------------------------------------------------------------------'
         write(*,'(4A18)')    aux_labs(1:min(4,aux_size))
         write(*,'(4ES18.8)') aux_vals(1:min(4,aux_size))
         loop_ind = 5;
         do while (aux_size > loop_ind .and. loop_ind < 100)
            write(*,*) ' '
            write(*,'(4A18)')    aux_labs(loop_ind:min(loop_ind+3,aux_size))
            write(*,'(4ES18.8)') aux_vals(loop_ind:min(loop_ind+3,aux_size))
            loop_ind = loop_ind + 4
         end do
         write(*,*) '-----------------------------------------------------------------------'
      end if
      

      call fatal_error(reason,call_loc,fname)
   end if
   
   return
end subroutine check_site_c13
!==========================================================================================!
!==========================================================================================!


!==========================================================================================!
!==========================================================================================!
! Check if a C-13 var too high, ignoring tiny-value comparison, and get it's delta value.  !
!------------------------------------------------------------------------------------------!
subroutine check_c13(heavy,total,check_delta,delta,valid,vname,reason)
   use consts_coms, only : tiny_num    ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real         , intent(in)    :: heavy
   real         , intent(in)    :: total
   logical      , intent(in)    :: check_delta
   real         , intent(out)   :: delta
   logical      , intent(out)   :: valid
   character(*) , intent(in)    :: vname
   character(*) , intent(out)   :: reason
   !----- Local Vars ----------------------------------------------------------------------!
   logical                    :: all_bets_are_off
   !---------------------------------------------------------------------------------------!
   delta = htIsoDelta(heavy,total)
   
   all_bets_are_off = abs(total) < 0.00000001 ! 10E-8
   if (all_bets_are_off) then
      valid = .true.
      return
   end if
   
   if (heavy > tiny_num) then
      if (total > tiny_num) then
         if (heavy > total) then
            valid = .false.
         else
            valid = .true.
         end if
      else
         valid = .false.
      end if
   else
      valid = .true.
   end if
   
   if (not(valid)) then
      reason = 'Bad ' // trim(vname) // ' total C-13'
      return
   end if
   
   ! Last condition here, delta /= delta is a NaN check.
   if ( ( abs(total) > tiny_num) .and. check_delta .and. &
        ( delta <  -200.0   .or. tiny_num < delta .or. delta /= delta)) then
      valid = .false.
      reason = 'Bad ' // trim(adjustl(vname)) // ' delta C-13'
   end if

end subroutine check_c13
!==========================================================================================!
!==========================================================================================!

end module iso_utils
