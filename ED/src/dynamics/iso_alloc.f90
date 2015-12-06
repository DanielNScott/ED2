module iso_alloc
implicit none
contains

!==========================================================================================!
! alloc_c13 is basically a wrapper function for calling the different allocation schemes,  !
! but it includes a lot of error checking.                                                 !
!==========================================================================================!
subroutine alloc_c13 (cpatch,ico,tr_bleaf,tr_broot,tr_bsapwooda,tr_bsapwoodb,daily_C_gain, &
                      daily_c13_gain,carbon_balance,carbon13_balance,st_h2tc_in)
   use ed_state_vars,   only : patchtype        ! ! intent(in)
   use isotopes    ,    only : c13af            ! ! intent(in)
   implicit none
   
   !------ Arguments. ---------------------------------------------------------------------!
   type(patchtype) , target         :: cpatch
   integer         , intent(in)     :: ico
   real            , intent(in)     :: tr_bleaf          ! Xfer of total C from calling fn
   real            , intent(in)     :: tr_broot          ! Xfer of total C from calling fn
   real            , intent(in)     :: tr_bsapwooda      ! Xfer of total C from calling fn
   real            , intent(in)     :: tr_bsapwoodb      ! Xfer of total C from calling fn
   real            , intent(in)     :: daily_C_gain
   real            , intent(in)     :: daily_c13_gain
   real            , intent(in)     :: carbon_balance
   real            , intent(in)     :: carbon13_balance
   real            , intent(in)     :: st_h2tc_in
   
   !------ Local Vars. --------------------------------------------------------------------!
   real(kind=8)   :: sum_c13_in              ! Pre-allocation C-13
   real(kind=8)   :: sum_c13_out             ! Post-allocation C-13
   real(kind=8)   :: sum_c13_dif             ! The difference of the two
   real           :: sum_tr                  ! Sum of xfers from storage to tissues
   real           :: sum_C                   ! Sum of plant C excluding storage.
   real           :: bl_c13_in, bl_c13_out   ! Leaf C-13 content before & after allocation
   real           :: br_c13_in, br_c13_out   ! Root C-13 content before & after allocation
   real           :: sa_c13_in, sa_c13_out   ! SapA C-13 content before & after allocation
   real           :: sb_c13_in, sb_c13_out   ! SapB C-13 content before & after allocation
   real           :: st_c13_in, st_c13_out   ! Stor C-13 content before & after allocation
   integer        :: print_details           !
   logical        :: stop_run                ! Flag to stop the run if an error is serious.
   character(70)  :: alloc_msg               ! Tells us where a problem was if encountered.
   character(20)  :: reason   
   !---------------------------------------------------------------------------------------!

   stop_run       = .false.
   print_details  = .false.
   
   !---------------------------------------------------------------------------------------!
   ! Please see the note in plant_cbal_sanity_check, this is not a meaningful check.       !
   !---------------------------------------------------------------------------------------!
   !if (abs(carbon13_balance) > abs(carbon_balance) .and. check_erythang) then
   !   write (*,*) '------------------------------------------------------------------------'
   !   write (*,*) 'Input check in alloc_c13: abs(carbon13_balance) > abs(carbon_balance)   !'
   !   write (*,*) '------------------------------------------------------------------------'
   !   write (*,'(A13, I4, I3)') 'Cohort, PFT: ', ico, cpatch%pft(ico)
   !   write (*,*) 'carbon13_balance, carbon_balance : ', carbon13_balance, carbon_balance
   !   write (*,*) ''
   !end if
   

   !---------------------------------------------------------------------------------------!
   ! Save the input C-13 fields.                                                           !
   !---------------------------------------------------------------------------------------!
   bl_c13_in = cpatch%bleaf_c13(ico)
   br_c13_in = cpatch%broot_c13(ico)
   sa_c13_in = cpatch%bsapwooda_c13(ico)
   sb_c13_in = cpatch%bsapwoodb_c13(ico)
   st_c13_in = cpatch%bstorage_c13(ico)

   sum_c13_in = dble(cpatch%bleaf_c13    (ico))                                            &
              + dble(cpatch%broot_c13    (ico))                                            &
              + dble(cpatch%bsapwooda_c13(ico))                                            &
              + dble(cpatch%bsapwoodb_c13(ico))                                            &
              + dble(cpatch%bstorage_c13 (ico))                                            &
              + dble(carbon13_balance)
   !---------------------------------------------------------------------------------------!

   

   !---------------------------------------------------------------------------------------!
   ! Allocate the C-13.                                                                    !
   !---------------------------------------------------------------------------------------!
   !select case(c13af)
   !case(1)
      call mech_alloc_2 (cpatch,ico,carbon_balance,carbon13_balance,tr_bleaf,tr_broot,     &
                         tr_bsapwooda,tr_bsapwoodb,st_h2tc_in,alloc_msg)
   !case(2)
   !   call nm_alloc_1(cpatch,ico,carbon13_balance)
   !case(3)
   !   call quad_alloc(cpatch,ico,carbon13_balance,lh2tc_in,sth2tc_in)
   !case(4)
   !   call l_r_diff(cpatch,ico,carbon13_balance,lh2tc_in,rh2tc_in,sth2tc_in,bleaf_in)
   !case default
   !   fatal_error('c13af Not in (1,2,4) but c13_alloc is called!','alloc_c13','iso_alloc.f90')
   !end select
   
   ! Update the living biomass.
   cpatch%balive_c13(ico) = cpatch%bleaf_c13    (ico) + cpatch%broot_c13    (ico)          &
                          + cpatch%bsapwooda_c13(ico) + cpatch%bsapwoodb_c13(ico)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! Alias the output C-13 fields.                                                         !
   !---------------------------------------------------------------------------------------!
   bl_c13_out = cpatch%bleaf_c13(ico)
   br_c13_out = cpatch%broot_c13(ico)
   sa_c13_out = cpatch%bsapwooda_c13(ico)
   sb_c13_out = cpatch%bsapwoodb_c13(ico)
   st_c13_out = cpatch%bstorage_c13(ico)
 
   sum_c13_out = dble(cpatch%bleaf_c13    (ico))                                           &
               + dble(cpatch%broot_c13    (ico))                                           &
               + dble(cpatch%bsapwooda_c13(ico))                                           &
               + dble(cpatch%bsapwoodb_c13(ico))                                           &
               + dble(cpatch%bstorage_c13 (ico))
   !---------------------------------------------------------------------------------------!
   

   
   !---------------------------------------------------------------------------------------!
   ! Check no C-13 was lost.                                                               !
   !---------------------------------------------------------------------------------------!
   if (( abs(sum_c13_dif / sum_c13_in ) > 0.000001 .or.                                    & 
         abs(sum_c13_dif / sum_c13_out) > 0.000001) ) then
      print_details = .true.
      stop_run      = .true.
      reason        = 'Greater than 1 millionth of total C-13 lost.'
   end if
   !---------------------------------------------------------------------------------------!

   
   !---------------------------------------------------------------------------------------!
   ! Check no field was set < 0.                                                           !
   !---------------------------------------------------------------------------------------!   
   if ((cpatch%bleaf_c13(ico)   < 0.0         .or. cpatch%broot_c13(ico)     < 0.0 .or.   &
      cpatch%bsapwooda_c13(ico) < 0.0         .or. cpatch%bsapwoodb_c13(ico) < 0.0 .or.   &
      cpatch%bstorage_c13(ico)  < -0.000000001 )) then
      print_details = .true.
      stop_run      = .true.
      reason        = 'A C-13 field is negative!'
   end if
   !---------------------------------------------------------------------------------------!

      
      
   !---------------------------------------------------------------------------------------!
   ! Check C-13 <= total C                                                                 !
   !---------------------------------------------------------------------------------------!
   if ((cpatch%bleaf_c13(ico)   > cpatch%bleaf(ico)     + tr_bleaf     .or.   & 
      cpatch%broot_c13(ico)     > cpatch%broot(ico)     + tr_broot     .or.   &
      cpatch%bsapwooda_c13(ico) > cpatch%bsapwooda(ico) + tr_bsapwooda .or.   &
      cpatch%bsapwoodb_c13(ico) > cpatch%bsapwoodb(ico) + tr_bsapwoodb        )) then
      !cpatch%bstorage_c13(ico)  > cpatch%bstorage(ico)         ) .and. check_erythang) then
      print_details = .true.
      stop_run      = .true.
      reason        = 'Total C-13 is greater than total C!'
   end if
   !---------------------------------------------------------------------------------------!
   

   
   !---------------------------------------------------------------------------------------!
   ! Print any issue we may have found.                                                    !
   !---------------------------------------------------------------------------------------!
   if (print_details) then
      sum_tr = tr_bleaf + tr_broot + tr_bsapwooda + tr_bsapwoodb
      sum_C  = cpatch%bleaf(ico) + cpatch%broot(ico) + cpatch%bsapwooda(ico)               &
             + cpatch%bsapwoodb(ico)
      
      write(*,*) '-----------------------------------------------------------------'
      write(*,*) ' C-13 was not conserved in alloc_c13!'
      write(*,*) ' ', reason
      write(*,*) '-----------------------------------------------------------------'
      write(*,'(A18, I4, I3)') ' Cohort, PFT     :', ico, cpatch%pft(ico)
      write(*,*) ' '
      write(*,*) '----------- Inputs ----------------------------------------------'
      write(*,*) ' carbon_balance  :', carbon_balance, carbon13_balance
      write(*,*) ' bstorage        :', cpatch%bstorage(ico)
      write(*,*) ' bstorage d13C   :', htIsoDelta(st_c13_in,cpatch%bstorage(ico))
      write(*,*) ' bstorage ratio  :', st_h2tc_in
      write(*,*) ' '
      write(*,*) '- Context Vars:  ----- Total C ----- 13C Value ---------------------'
      write(*,*) ' daily_C_gain    :', daily_C_gain                  , daily_c13_gain
!      write(*,*) ' today_gpp       :', cpatch%today_gpp(ico)         , cpatch%today_gpp_c13(ico)
!      write(*,*) ' today_leaf_resp :', cpatch%today_leaf_resp(ico)   , cpatch%today_leaf_resp_c13(ico)
!      write(*,*) ' today_root_resp :', cpatch%today_root_resp(ico)   , cpatch%today_root_resp_c13(ico)
      write(*,*) ' '
      write(*,*) '- Variable Name: ---- Input Val ---- Var Will Be --- Xfer to Var ---'
      write(*,*) ' bleaf           :', cpatch%bleaf(ico)    , cpatch%bleaf(ico)     + tr_bleaf    , tr_bleaf    
      write(*,*) ' broot           :', cpatch%broot(ico)    , cpatch%broot(ico)     + tr_broot    , tr_broot    
      write(*,*) ' bsapwooda       :', cpatch%bsapwooda(ico), cpatch%bsapwooda(ico) + tr_bsapwooda, tr_bsapwooda
      write(*,*) ' bsapwoodab      :', cpatch%bsapwoodb(ico), cpatch%bsapwoodb(ico) + tr_bsapwoodb, tr_bsapwoodb
      write(*,*) ' TOTAL           :', sum_C                , sum_C + sum_tr                      , sum_tr
      write(*,*) '' 
      write(*,*) '---------- Subroutine Results -----------------------------------'
      write(*,*) '- Variable Name: ---- Input Val ----- Out Val ---- Difference ---'
      write(*,*) ' bleaf_c13       :', bl_c13_in  , bl_c13_out  , bl_c13_out - bl_c13_in
      write(*,*) ' broot_c13       :', br_c13_in  , br_c13_out  , br_c13_out - br_c13_in
      write(*,*) ' bsapwooda_c13   :', sa_c13_in  , sa_c13_out  , sa_c13_out - sa_c13_in
      write(*,*) ' bsapwoodb_c13   :', sb_c13_in  , sb_c13_out  , sb_c13_out - sb_c13_in
      write(*,*) ' bstorage_c13    :', st_c13_in  , st_c13_out  , st_c13_out - st_c13_in
      write(*,*) ' TOTAL           :', sum_c13_in ,sum_c13_out  ,sum_c13_out -sum_c13_in
      write(*,*) ''
      write(*,*) alloc_msg
   end if
   
   if (stop_run) then
      call fatal_error(reason,'alloc_c13','iso_alloc.f90')
   end if
   !---------------------------------------------------------------------------------------!

end subroutine alloc_c13
!==========================================================================================!

   
!==========================================================================================!
subroutine mech_alloc_2(cpatch   , ico       , carbon_balance, carbon13_balance,         &
                        tr_bleaf , tr_broot  , tr_bsapwooda  , tr_bsapwoodb    ,         &
                        sth2tc_in, alloc_msg)
   use ed_state_vars, only : patchtype
   implicit none
   
   !------ Arguments. ---------------------------------------------------------------------!
   type(patchtype) , target      :: cpatch
   integer         , intent(in)  :: ico
   real            , intent(in)  :: carbon_balance
   real            , intent(in)  :: carbon13_balance
   real            , intent(in)  :: tr_bleaf
   real            , intent(in)  :: tr_broot
   real            , intent(in)  :: tr_bsapwooda
   real            , intent(in)  :: tr_bsapwoodb
   real            , intent(in)  :: sth2tc_in
   character(70)   , intent(out) :: alloc_msg
   
   !------ Local Var. ---------------------------------------------------------------------!
   real                          :: sum_tr
   real                          :: sum_tr_c13
   real                          :: tr_bleaf_c13
   real                          :: tr_broot_c13
   real                          :: tr_bsapwooda_c13
   real                          :: tr_bsapwoodb_c13
   real                          :: tr_bstorage_c13
   !---------------------------------------------------------------------------------------!
   
   sum_tr = tr_bleaf + tr_broot + tr_bsapwooda + tr_bsapwoodb
    
   if (sum_tr > 0.0) then
   !---------------------------------------------------------------------------------------!
   !     We're about to transfer tr_xx vars into pools. (They're positive.)                !
   !                                                                                       !
   ! The fraction of total carbon being moved into living pools that comes from carbon     !
   ! balance will be 1 if carbon_balance > sum_tr and carbon_balance / sum_tr otherwise.   !
   !                                                                                       !
   ! Note that carbon_balance will always be positive because growth_resp is directly prop !
   ! to daily_C_gain by growth_resp_factor; GRF is in [0,1) by necessity else there would  !
   ! be no carbon gain and nothing would ever grow.                                        !
   !                                                                                       !
   ! Similarly, carbon13_balance must be >= 0 by the same proportionality. (If no frac.)   !
   ! Here we transfer 13C in proportion to the xfer coming from it's source, and if there  !
   ! wasn't enough we then rectify the matter.                                             !
   !---------------------------------------------------------------------------------------!
      if (carbon_balance >= sum_tr) then
         tr_bleaf_c13     = carbon13_balance / carbon_balance * tr_bleaf
         tr_broot_c13     = carbon13_balance / carbon_balance * tr_broot
         tr_bsapwooda_c13 = carbon13_balance / carbon_balance * tr_bsapwooda
         tr_bsapwoodb_c13 = carbon13_balance / carbon_balance * tr_bsapwoodb
         
         sum_tr_c13 = sum_tr * carbon13_balance / carbon_balance
         
         !Truncating round-off error... a bad thing here, or a VERY bad thing...?
         !It's typically less than 10E-9, and we really want to keep things >= 0...
         tr_bstorage_c13  = max(0.0,carbon13_balance - sum_tr_c13)
         
         alloc_msg = 'Allocation Conditional (see src): carbon_balance >= sum_tr'
         
      elseif (carbon_balance < sum_tr) then
      !------------------------------------------------------------------------------------!
      ! Formula ( 13C_bal + sth2tc * (sum_tr - c_bal) ) / sum_tr * bleaf simplified from:  !
      ! ((13C_bal / C_bal) * (C_bal/sum_tr) + sth2tc * (sum_tr - c_bal)/ sum_tr )* bleaf   !
      !------------------------------------------------------------------------------------!
         tr_bleaf_c13     = ( carbon13_balance + sth2tc_in *(sum_tr - carbon_balance) )    &
                            /sum_tr *tr_bleaf
         tr_broot_c13     = ( carbon13_balance + sth2tc_in *(sum_tr - carbon_balance) )    &
                            /sum_tr *tr_broot
         tr_bsapwooda_c13 = ( carbon13_balance + sth2tc_in *(sum_tr - carbon_balance) )    &
                            /sum_tr *tr_bsapwooda
         tr_bsapwoodb_c13 = ( carbon13_balance + sth2tc_in *(sum_tr - carbon_balance) )    &
                            /sum_tr *tr_bsapwoodb
         !Truncating round-off error... a bad thing here, or a VERY bad thing...?
         !It's typically less than 10E-9, and we really want to keep things >= 0...
         tr_bstorage_c13 = - min( sth2tc_in* (sum_tr - carbon_balance)                     &
                                , cpatch%bstorage_c13(ico))

      alloc_msg = 'Allocation Conditional (see src): carbon_balance < sum_tr'
      
      else
         tr_bleaf_c13     = 0.0
         tr_broot_c13     = 0.0
         tr_bsapwooda_c13 = 0.0
         tr_bsapwoodb_c13 = 0.0
         tr_bstorage_c13  = 0.0
      end if
   
      cpatch%bleaf_c13(ico)     = cpatch%bleaf_c13(ico)     + tr_bleaf_c13
      cpatch%broot_c13(ico)     = cpatch%broot_c13(ico)     + tr_broot_c13
      cpatch%bsapwooda_c13(ico) = cpatch%bsapwooda_c13(ico) + tr_bsapwooda_c13
      cpatch%bsapwoodb_c13(ico) = cpatch%bsapwoodb_c13(ico) + tr_bsapwoodb_c13
      cpatch%bstorage_c13(ico)  = cpatch%bstorage_c13(ico)  + tr_bstorage_c13
      
      !----- If there wasn't enough in storage, 'undo' some proportion of the xfers -------!
      ! if (cpatch%bstorage_c13(ico) < 0.0) then
      !    call dist_c13_prop(cpatch%bstorage_c13(ico),cpatch%bleaf_c13(ico)                 &
      !                      ,cpatch%broot_c13(ico),cpatch%bsapwooda_c13(ico)                &
      !                      ,cpatch%bsapwoodb_c13(ico))
      !    cpatch%bstorage_c13(ico) = 0.0  
      ! end if
   else
      !---------------------------------------------------------------------------------------!
      ! We just put carbon_balance into storage in main routine, only reached if flphen == 1. !
      !---------------------------------------------------------------------------------------!
      cpatch%bstorage_c13(ico) = cpatch%bstorage_c13(ico) + carbon13_balance
      alloc_msg = 'Allocation Conditional (see src): sum_tr is not greater than 0.0'
      !---------------------------------------------------------------------------------------!
   end if
end subroutine mech_alloc_2
!=======================================================================================!      


!==========================================================================================!
!  Here we 'guarantee' fixed permil relationships between pools and let the whole    !
! structure float as necessary to maintain that. Also, we don't care about loop flg. !
! The reason this is 'zero' and the next is 'first' is because here we pay NO heed   !
! to where our carbon 13 just was. We just redistribute willy nilly. It's a party.   !
!==========================================================================================!   
subroutine nm_alloc_1(cpatch,ico,carbon13_balance)
   !use isotope_utils, only : nelder_mead           & ! subroutine
   !                        , iso_err_fn            & ! function
   !                        , htIsoDelta            ! ! function
!   use iso_checks          , check_patch_c13      ! ! subroutine
   use ed_state_vars, only : patchtype             ! ! intent(inout)
   implicit none
   
   !------ Arguments. ---------------------------------------------------------------------!
   type(patchtype) , target        :: cpatch
   integer         , intent(in)    :: ico
   real            , intent(in)    :: carbon13_balance
   
   !------ Local Var. ---------------------------------------------------------------------!   
   real, dimension(5)            :: pools
   real, dimension(5)            :: labels
   integer                       :: fn_flg

   real                          :: tc13    ! Total 13C (w/o bdead)
   real, dimension(3)            :: tC      ! Total C, by pool (w/o bdead)
   real, dimension(2)            :: sol     ! Best solution found
   real, dimension(4,2)          :: x       ! Used for passing in centroid etc.
   real, dimension(4)            :: p       ! Params, see below.
   real, dimension(3,2)          :: pts     ! Simplex
   real, dimension(2,2)          :: c       ! Constraints on domains of coordinates
   integer                       :: step    ! Branch of the algorithm being executed
   integer                       :: cnt     ! Current recursion depth

   labels = (/1.,2.,3.,4.,5./)
  
   !---------------------------------------------------------------------------------------!
   ! Assign 'total' c13 (Doesn't include bdead, which we always just grow from storage)    !
   ! and 'total' cohort carbon (Also ignoring bdead)                                       !
   !---------------------------------------------------------------------------------------!
   tc13  =  max(0.0,   cpatch%bleaf_c13    (ico) + cpatch%broot_c13    (ico)               &
                     + cpatch%bsapwooda_c13(ico) + cpatch%bsapwoodb_c13(ico)               &
                     + cpatch%bstorage_c13 (ico) + carbon13_balance )
   
   tC(1) =  cpatch%bleaf(ico)
   tC(2) =  cpatch%broot(ico)
   tC(3) =  cpatch%bsapwooda(ico) + cpatch%bsapwoodb(ico) + cpatch%bstorage(ico)

   !----- If cohort carbon is 0, we should stop right here... -----------------------------!
   if (sum(tC) < tiny(1.0)) then
      cpatch%bleaf_c13(ico)     = 0.0
      cpatch%broot_c13(ico)     = 0.0
      cpatch%bsapwooda_c13(ico)     = 0.0
      cpatch%bsapwoodb_c13(ico)     = 0.0
      cpatch%bstorage_c13 (ico)     = 0.0
      return
   end if
   
   !----- Make sure tc13 isn't larger than tC ---------------------------------------------!
   if (tc13 >= sum(tC)) then
      write (*,*) '-----------------------------------------------------------'
      write (*,*) 'ERROR? tc13 >= tC'
      write (*,*) 'tC   = ', sum(tC)
      write (*,*) 'tc13 = ', tc13
      write (*,*) 'carbon13_balance = ', carbon13_balance
      write (*,*) 'If this is "tiny" i.e. fortran ~ 0, it may not be an issue.'
      cpatch%bleaf_c13(ico) = cpatch%bleaf(ico)
      cpatch%broot_c13(ico) = cpatch%broot(ico)
      cpatch%bsapwooda_c13(ico) = cpatch%bsapwooda(ico)
      cpatch%bsapwoodb_c13(ico) = cpatch%bsapwoodb(ico)
      cpatch%bstorage_c13 (ico) = cpatch%bstorage(ico)
      return
   end if
   
   !---------------------------------------------------------------------------------------!
   ! If we have no leaves or roots we should use a different error function, and if both   !
   ! are zero we should give everything in 's' the same ratio. Also, the initial simplex   !
   ! depends on this because we want to know if we're reducing the dimensionality, i.e.    !
   ! putting in zeros for leaf or root coords.                                             !
   !---------------------------------------------------------------------------------------!
   if (cpatch%bleaf(ico) < tiny(1.0) .and. cpatch%broot(ico) < tiny(1.0)) then
      cpatch%bleaf_c13(ico)     = 0.0
      cpatch%broot_c13(ico)     = 0.0
      cpatch%bsapwooda_c13(ico) = tc13 * cpatch%bsapwooda(ico) / tC(3)
      cpatch%bsapwoodb_c13(ico) = tc13 * cpatch%bsapwoodb(ico) / tC(3)
      cpatch%bstorage_c13 (ico) = tc13 * cpatch%bstorage (ico) / tC(3)
      return
      
   else if (cpatch%bleaf(ico) < tiny(1.0) .and. cpatch%broot(ico) > tiny(1.0)) then
      cpatch%bleaf_c13(ico)   = 0.0
      fn_flg = 1
      pts = reshape  ( (/               0.0 ,                             0.0              &
                         ,              0.0 , tc13*tC(2)/sum(tC)                           &
                         ,tc13*tC(2)/sum(tC), tc13*(tC(2)-0.001*tC(2))/sum(tC) /)          &
                     ,(/3,2/) )      
      
   else if (cpatch%bleaf(ico) > tiny(1.0) .and. cpatch%broot(ico) < tiny(1.0)) then
      cpatch%broot_c13(ico)   = 0.0
      fn_flg = 2
      pts = reshape  ( (/ tc13*tC(1)/sum(tC), tc13*(tC(1)-0.001*tC(1))/sum(tC)             &
                         ,tc13*tC(1)/sum(tC),                             0.0              &
                         ,              0.0 ,                             0.0  /)          &
                     ,(/3,2/) )
   else
      fn_flg = 0
      pts = reshape  ( (/ tc13*tC(1)/sum(tC), tc13*(tC(1)-0.001*tC(1))/sum(tC)             &
                         ,tc13*tC(1)/sum(tC), tc13*tC(2)/sum(tC)                           &
                         ,tc13*tC(2)/sum(tC), tc13*(tC(2)-0.001*tC(2))/sum(tC) /)          &
                     ,(/3,2/) )
   end if
   !---------------------------------------------------------------------------------------!
   
   !----- Initialize the ones driving the nm_alg (except pts)------------------------------!
   sol(:)   = 0.                       ! Best solution found (init to 0)
   x(:,:)   = 0.                       ! Used for passing in centroid etc. (init to 0)
   p        = (/1.0,2.0,-0.5,0.5/)     ! Params, see subroutine nelder_mead
   c(:,1)   = 0.                       ! Constraints on domains of coordinates
   c(:,2)   = tc13                     ! Constraints on domains of coordinates
   step     = 1.                       ! Branch of the algorithm being executed
   cnt      = 0.                       ! Current recursion depth
   
   !------- Find and assign the carbon 13 distribution ------------------------------------!
   call nelder_mead(tc13,tC,sol,x,p,pts,c,step,cnt,fn_flg)

   cpatch%bleaf_c13(ico)     = sol(1)
   cpatch%broot_c13(ico)     = sol(2)
   cpatch%bsapwooda_c13(ico) = (tc13 - sol(1) - sol(2)) * cpatch%bsapwooda(ico) / tC(3)
   cpatch%bsapwoodb_c13(ico) = (tc13 - sol(1) - sol(2)) * cpatch%bsapwoodb(ico) / tC(3)
   cpatch%bstorage_c13 (ico) = (tc13 - sol(1) - sol(2)) * cpatch%bstorage (ico) / tC(3)
   
   cpatch%balive_c13(ico) = cpatch%bleaf_c13    (ico) + cpatch%broot_c13    (ico)          &
                          + cpatch%bsapwooda_c13(ico) + cpatch%bsapwoodb_c13(ico)
   
!   call check_patch_c13(cpatch,ico,"Leaving alloc_c13")

end subroutine nm_alloc_1
!=======================================================================================!



!=======================================================================================!
subroutine l_r_diff(cpatch,ico,carbon13_balance,lh2tc_in,rh2tc_in,sth2tc_in,bleaf_in)
   !use isotope_utils, only : iso_err_fn            & ! function
   !                        , htIsoDelta            & ! function
   !                        , htIsoDelta8           & ! function
   !                        , d13C2Ratio      & ! function
   !                        , d13C2Ratio8     & ! function
   !                        , get_lhc_target        ! ! subroutine
!   use iso_checks          , check_patch_c13      ! ! subroutine
   use isotopes,  only : iso_P1                & ! intent(in)
                           , iso_P2                & ! intent(in)
                           , R_std                 ! ! intent(in)
   use ed_state_vars, only : patchtype             ! ! intent(inout)
   implicit none
   
   !------ Arguments. ---------------------------------------------------------------------!
   type(patchtype) , target      :: cpatch
   integer         , intent(in)  :: ico
   real            , intent(in)  :: carbon13_balance
   real            , intent(in)  :: lh2tc_in           ! Pre pre alloc_c_bal. leaf ratio
   real            , intent(in)  :: rh2tc_in           ! Pre pre alloc_c_bal. root ratio
   real            , intent(in)  :: sth2tc_in          ! Pre pre alloc_c_bal. stor. ratio
   real            , intent(in)  :: bleaf_in           ! Pre alloc_c_bal leaf C
   
   !------ Local Var. ---------------------------------------------------------------------!
   real(kind=8)                  :: lh2tc                ! Leaf 13C:12C that we assign
   real(kind=8)                  :: rh2tc                ! root 13C:12C that we assign
   real(kind=8)                  :: lhc_target           ! Aim for bleaf 13C
   real                          :: lhc_target4          ! Aim for bleaf 13C
   real(kind=8)                  :: rhc_target           ! Aim for broot 13C
   real(kind=8)                  :: lhc_diff             ! bleaf 13C - lhc_target
   real(kind=8)                  :: rhc_diff             ! broot 13C - lhc_target

   real(kind=8)                  :: lc , sac , sbc , stc , rc  ! Aliases for C pools
   real(kind=8)                  :: lhc, sahc, sbhc, sthc, rhc ! Heavy C pools
   real(kind=8)                  :: lhc_in, sahc_in, sbhc_in, sthc_in, rhc_in ! Heavy C pools
   
   real(kind=8)                  :: tc13                 ! Total 13C (leaf, root, s-)
   real(kind=8)                  :: dist_c13             ! 13C to distribute

   real(kind=8)                  :: d13Cleaf             ! d13C of leaves 
   integer                       :: prtflg               ! Something's weird? Print it!
   
   prtflg = 0
   
   !---------------------------------------------------------------------------------------!
   ! Assign 'total' c13 (Doesn't include bdead, which we always just grow from storage)    !
   ! and 'total' cohort carbon (Also ignoring bdead)                                       !
   !---------------------------------------------------------------------------------------!
   tc13     = max(0.0, cpatch%bleaf_c13    (ico) + cpatch%broot_c13    (ico)               &
                     + cpatch%bsapwooda_c13(ico) + cpatch%bsapwoodb_c13(ico)               &
                     + cpatch%bstorage_c13 (ico) + carbon13_balance )
                     
   dist_c13 = tc13
   if (dist_c13 < tiny(1.0)) then
      write(*,*) 'dist_c13 < 0, init'
   end if
   
   !----------------- Double precision aliases --------------------------------------------!
   lc     =  dble(cpatch%bleaf(ico))
   rc     =  dble(cpatch%broot(ico))
   sac    =  dble(cpatch%bsapwooda(ico))
   sbc    =  dble(cpatch%bsapwoodb(ico))
   stc    =  dble(cpatch%bstorage (ico))

   lhc    =  dble(cpatch%bleaf_c13(ico))
   rhc    =  dble(cpatch%broot_c13(ico))
   sahc   =  dble(cpatch%bsapwooda_c13(ico))
   sbhc   =  dble(cpatch%bsapwoodb_c13(ico))
   sthc   =  dble(cpatch%bstorage_c13 (ico))
   
   lhc_in    =  dble(cpatch%bleaf_c13(ico))
   rhc_in    =  dble(cpatch%broot_c13(ico))
   sahc_in   =  dble(cpatch%bsapwooda_c13(ico))
   sbhc_in   =  dble(cpatch%bsapwoodb_c13(ico))
   sthc_in   =  dble(cpatch%bstorage_c13 (ico))
      
   !----- If cohort carbon is 0, we should stop right here... ----------------------------!
   control: do
   if (lc + rc + sac + sbc + stc < tiny(1.0)) then
      write(*,*) '!---------------------- tCl + tCr + tCs < tiny -------------------------! '
      lhc    = 0.d0
      rhc    = 0.d0
      sahc   = 0.d0
      sbhc   = 0.d0
      sthc   = 0.d0
      exit control
   end if
      
   !---------------------------------------------------------------------------------------!
   ! First check if this cohort has leaves. We have the following cases to consider:       !
   !  - We have leaves and gpp, so set leaf ratio = gpp ratio                              !
   !  - We have leaves but no gpp (LAI unresolvable, etc.)                                 !
   !     - Just flushing so there is no previous ratio - make the ratio the storage ratio  !
   !     - There is some previous ratio, which we want to keep                             !
   !  - No leaves, so leaf heavy carbon = 0.0                                              !
   !  Why don't we have to worry about storage = 0.0?                                      !
   !---------------------------------------------------------------------------------------!
   lhc_target4 = 0.0
   !call get_lhc_target(cpatch%pft(ico),cpatch%bleaf(ico),cpatch%bleaf_c13(ico),cpatch%today_gpp(ico)  &
   !                   ,cpatch%today_gpp_c13(ico),lh2tc_in,sth2tc_in,bleaf_in,lhc_target4)
   lhc_target = dble(lhc_target4)

   ! if (cpatch%bleaf(ico) > tiny(1.0)) then
      ! if (cpatch%today_gpp(ico) > tiny(1.0)) then
         ! lh2tc = cpatch%today_gpp_c13(ico)/cpatch%today_gpp(ico)
      ! else
         ! if (lh2tc_in > tiny(1.0)) then
            ! lh2tc = lh2tc_in
         ! else
            ! lh2tc = sth2tc_in
         ! end if
      ! end if
   ! else
      ! lh2tc  = 0.0
   ! end if
   
   ! lhc_target = lh2tc * lc
   lhc_diff   = lhc - lhc_target
   
   if (lhc_diff > tiny(1.0) .or. -lhc_diff < tc13) then
      lhc       = lhc_target
      dist_c13  = tc13 - lhc
      if (dist_c13 < 0.0) then
         write(*,*) 'dist_c13 < 0, lhc assignment'
      end if
   else
      !--------------------- THIS SHOULDN'T HAPPEN -------------------------------------!
      write(*,*) '!------------------------ lhc_diff < -tc13 ---------------------------! '
      lhc   = tc13
      rhc   = 0.0
      sahc  = 0.0
      sbhc  = 0.0
      sthc  = 0.0
      exit control
   end if
   
   !---------------------------------------------------------------------------------------!
   ! Backtrack for a second: dist_c13 at this point shouldn't be allowed to be > tCr + tCs !
   ! if we want to distribute it. So if it is, we have no choice about leaf ratio.         !
   !---------------------------------------------------------------------------------------!
   if (dist_c13 > rc + sac + sbc + stc) then
      write(*,*) '!---------------- dist_c13 > rc + sac + sbc + stc -----------------------! '
      lhc  = lhc + rc + sac + sbc + stc - (rhc + sahc + sbhc + sthc)
      rhc  = rc
      sahc = sac
      sbhc = sbc
      sthc = stc
      exit control
   end if
   
   !---------------------------------------------------------------------------------------!
   ! At this point dist_c13 is everything we want to distribute.                           !
   ! If we have no roots everything goes in s-. Otherwise we determine the ratio we want   !
   ! the roots to have and try to assign it.                                               !
   !---------------------------------------------------------------------------------------!
   if (cpatch%broot(ico) < tiny(1.0)) then
      rhc_target = 0.0
      rhc        = rhc_target
      
      sahc   = dist_c13 * sac / (sac + sbc + stc)
      sbhc   = dist_c13 * sbc / (sac + sbc + stc)
      sthc   = dist_c13 * stc / (sac + sbc + stc)
      exit control
   end if

   !--------------  Attempt to fix iso_P1 = d13Cleaf - d13Croot ---------------------------!
   if (cpatch%bleaf(ico) > tiny(1.0)) then
      d13Cleaf   = htIsoDelta8(lhc,lc)
      rhc_target = d13C2Ratio8(d13Cleaf - iso_P1) * rc
      rhc_diff   = rhc - rhc_target
   else
      rhc_target = dist_c13 * rc  / (sac + sbc + stc + rc)
      rhc_diff   = rhc - rhc_target
      rhc        = rhc_target
      
      sahc  = dist_c13 * sac / (sac + sbc + stc + rc)
      sbhc  = dist_c13 * sbc / (sac + sbc + stc + rc)
      sthc  = dist_c13 * stc / (sac + sbc + stc + rc)

      !rhc_target = rh2tc_in * rc
      !rhc_diff   = rhc - rhc_target
      exit control
   end if
   
   !-------- These should only be reached if we never exited control loop -----------------!
   rhc      =            min(rhc_target,dist_c13)
   dist_c13 = dist_c13 - min(rhc_target,dist_c13)
   if (dist_c13 < 0.0) then
      write(*,*) '!-------------- dist_c13 < 0, rhc assignment ----------------------------!'
   end if

   
   if (dist_c13 < sac + sbc + stc) then
      sahc     = dist_c13 * sac / (sac + sbc + stc)
      sbhc     = dist_c13 * sbc / (sac + sbc + stc)
      sthc     = dist_c13 * stc / (sac + sbc + stc)
   else
      dist_c13 = dist_c13 - sahc - sbhc  - sthc
      if (dist_c13 < 0.0) then
      write(*,*) '!---------------- dist_c13 < 0, dist_c13 > s- assignment ----------------!'
      end if

      rhc      = rhc + dist_c13
      sahc     = sac
      sbhc     = sbc
      sthc     = stc
   end if

   exit control
   end do control
   
   ! Convert Back.... ---------------------------------------------------------------------!
   cpatch%bleaf_c13(ico)     = lhc
   cpatch%broot_c13(ico)     = rhc
   cpatch%bsapwooda_c13(ico) = sahc
   cpatch%bsapwoodb_c13(ico) = sbhc
   cpatch%bstorage_c13 (ico) = sthc
   
   ! if (htIsoDelta8(lhc,lc) > 0.0 .or. htIsoDelta8(lhc,lc) < -65.0) then
      ! write(*,*) '!--------------------------------------------------------------------------!'
      ! write(*,*) ' d13C_leaf is outside [-65.0 , 0.0]... :', htIsoDelta8(lhc,lc)
      ! prtflg = 1
   ! end if
   ! if (htIsoDelta8(sahc,sac) > 0.0 .or. htIsoDelta8(sahc,sac) < -65.0) then
      ! write(*,*) '!--------------------------------------------------------------------------!'
      ! write(*,*) ' d13C_sapwooda is outside [-65.0 , 0.0]'
      ! prtflg = 1
   ! end if
   ! if (abs(htIsoDelta8(sahc,sac) - htIsoDelta8(sthc,stc)) > 0.1 .and. stc > tiny(1.0)) then
      ! write(*,*) '!--------------------------------------------------------------------------!'
      ! write(*,*) '|d13C_sa - d13C_st| > 1'
      ! prtflg = 1
   ! end if
   ! if (abs(htIsoDelta8(rhc_target,rc) - d13Cleaf) > 3.0001 .and. lc > tiny(1.0)) then
      ! write(*,*) '!--------------------------------------------------------------------------!'
      ! write(*,*) '|d13C rhc_target - d13C_leaf| > 3.0001'
      ! prtflg = 1
   ! end if   
   ! if (abs(htIsoDelta8(rhc_target,rc)) > 65.0 ) then
      ! write(*,*) '!--------------------------------------------------------------------------!'
      ! write(*,*) '|d13C rhc_target| > 65'
      ! prtflg = 1
   ! end if
   ! if ((lhc   < 0.0 .and.   lc > tiny(1.0)) .or.                                       &
       ! (rhc   < 0.0 .and.  rhc > tiny(1.0)) .or.                                       &
       ! (sahc  < 0.0 .and.  sac > tiny(1.0)) .or.                                       &
       ! (sbhc  < 0.0 .and.  sbc > tiny(1.0)) .or.                                       &
       ! (sthc  < 0.0 .and.  stc > tiny(1.0))      ) then
      ! write(*,*) '!--------------------------------------------------------------------------!'
      ! write(*,*) ' some hc field is less than zero... '
      ! prtflg = 1
   ! end if   
   ! if (prtflg == 1) then
      ! write(*,*) '!--------------------------------------------------------------------------!'
      ! write(*,*) ' Cohort, PFT    : ', ico , cpatch%pft(ico)
      ! write(*,*) ' C in           : ' 
      ! write(*,*) '   - lc         : ', lc
      ! write(*,*) '   - rc         : ', rc
      ! write(*,*) '   - sac        : ', sac
      ! write(*,*) '   - sbc        : ', sbc
      ! write(*,*) '   - stc        : ', stc
      ! write(*,*) ''
      ! write(*,*) ' c13 in, out    : ' 
      ! write(*,*) '   - lhc        : ', lhc_in  , lhc
      ! write(*,*) '   - rhc        : ', rhc_in  , rhc
      ! write(*,*) '   - sahc       : ', sahc_in , sahc
      ! write(*,*) '   - sbhc       : ', sbhc_in , sbhc
      ! write(*,*) '   - sthc       : ', sthc_in , sthc
      ! write(*,*) ''
      ! write(*,*) 'tc13            : ', tc13
      ! write(*,*) 'dist_c13        : ', dist_c13
      ! write(*,*) ''
      ! write(*,*) 'lhc_target      : ', lhc_target
      ! write(*,*) 'lhc_diff        : ', lhc_diff
      ! write(*,*) 'd13Cleaf        : ', d13Cleaf
      ! write(*,*) ''
      ! write(*,*) 'rhc_target      : ', rhc_target
      ! write(*,*) 'd13C rhc_target : ', htIsoDelta8(rhc_target,rc)
      ! write(*,*) 'rhc_diff        : ', rhc_diff
      ! write(*,*) ''
      ! write(*,*) 'lh2tc_in        : ', lh2tc_in
      ! write(*,*) 'sth2tc_in       : ', sth2tc_in
      ! write(*,*) ''
      ! write(*,*) 'dl , ds , dr    : ', htIsoDelta8(lhc,lc), htIsoDelta8(sahc + sbhc + sthc, sac + sbc + stc ), htIsoDelta8(rhc,rc)
      ! write(*,*) 'dsa, dsb, dst   : ', htIsoDelta8(sahc,sac), htIsoDelta8(sbhc,sbc), htIsoDelta8(sthc,stc)
      ! write(*,*) 'dsa - dst       : ', htIsoDelta8(sahc,sac) - htIsoDelta8(sthc,stc)
      ! write(*,*) ''
      ! write(*,*) 'today_gpp       : ', cpatch%today_gpp(ico)
      ! write(*,*) 'today_gpp_c13   : ', cpatch%today_gpp_c13(ico)
      ! write(*,*) 'delta gpp       : ', htIsoDelta(cpatch%today_gpp_c13(ico),cpatch%today_gpp(ico))
      ! write(*,*) ''
      ! write(*,*) '!--------------------------------------------------------------------------!'
   ! end if
   
!   call check_patch_c13(cpatch,ico,"end c13_alloc")
end subroutine l_r_diff
!=======================================================================================!

!==========================================================================================!
recursive subroutine nelder_mead(tc13,tC,sol,x,p,pts,c,step,cnt,fn_flg)
   implicit none
   
   !------ Notes: -------------------------------------------------------------------------!
   ! 1) nelder_mead is an implementation of the Nelder-Mead optimization algorithm
   ! 2) For initialization, x just be NaNs.
   ! 3) If you don't want to constrain x,y, set c(1,:) = (/-huge(1.0),+huge(1.0)/) and sim 
   !    for c(2,:)
   !
   
   !------ Arguments. ---------------------------------------------------------------------!
   real                , intent(inout)  :: tc13    ! Total 13C
   real, dimension(3)  , intent(in)     :: tC      ! Total C, by pool
   real, dimension(2)  , intent(inout)  :: sol     ! Best solution found
   real, dimension(4,2), intent(inout)  :: x       ! Used for passing in centroid etc.
   real, dimension(4)  , intent(in)     :: p       ! Params, see below.
   real, dimension(3,2), intent(inout)  :: pts     ! Simplex
   real, dimension(2,2), intent(inout)  :: c       ! Constraints on domains of coordinates
   integer             , intent(inout)  :: step    ! Branch of the algorithm being executed
   integer             , intent(inout)  :: cnt     ! Current recursion depth
   integer             , intent(in)     :: fn_flg  ! Determines error fn.
   
   !------ Local Var. ---------------------------------------------------------------------!
   integer, parameter      :: maxcnt = 490 ! Recursion limit
   real,    dimension(3)   :: labels       ! Used to label points for sort
   real,    dimension(3)   :: fpts         ! iso_err_fn value on pts
   real,    dimension(3,2) :: tpts         ! temporary points storage
   
   real,    dimension(2)   :: x0           ! Centroid
   real,    dimension(2)   :: xr           ! Reflected
   real,    dimension(2)   :: xc           ! Contracted point
   real,    dimension(2)   :: xe           ! Extended point
   
   integer                 :: i            ! Counter
   integer                 :: j            ! Counter

   
   if (step .eq. 1) then
      cnt = cnt+1
      !-------- Order the pts according to their fn values -------------------------------!
      fpts(1) = iso_err_fn(pts(1,:),tc13,tC,fn_flg)
      fpts(2) = iso_err_fn(pts(2,:),tc13,tC,fn_flg)
      fpts(3) = iso_err_fn(pts(3,:),tc13,tC,fn_flg)
      labels = (/1.,2.,3./)
      
      if (isNaN(fpts(1)) .or. isNaN(fpts(2)) .or. isNaN(fpts(3))) then
         write(*,*) "---------- ERROR! ------------"
         write(*,*) "About to call sort, fpts is..."
         write(*,*) fpts 
         write(*,*) "pts is"
         write (*,*) "      ", pts(1,1), pts(1,2)
         write (*,*) "      ", pts(2,1), pts(2,2)
         write (*,*) "      ", pts(3,1), pts(3,2)
         write(*,*) "cnt is ", cnt
      end if

      call sort(fpts,labels)
      tpts(1,:) = pts(nint(labels(1)),:)
      tpts(2,:) = pts(nint(labels(2)),:)
      tpts(3,:) = pts(nint(labels(3)),:)
      pts = tpts

      
      !------- Check to see if we're done -------------------------------------------------!
      if (iso_err_fn(pts(1,:),tc13,tC,fn_flg) < 0.0002 .or. cnt > maxcnt) then
         sol = pts(1,:)
      else
         !Calculate center of gravity of the points excluding the worst.
         !Worst point is last in s_pts, and centroid is analogue of mean:
         x0 = sum(pts(1:2,:),1)/2
         x(1,:) = x0
         step = 3
         call nelder_mead(tc13,tC,sol,x,p,pts,c,step,cnt,fn_flg)
         
      end if
   end if
    
   if (step .eq. 3) then
      cnt = cnt+1
      ! Unpack x as nec. (done for ease of coding / clarity)
      x0 = x(1,:);
      ! Compute domain constrained reflected point:
      xr = x0 + p(1)*(x0 - pts(3,:))

      i = 0
      do i = 1,2
         if (xr(i) .lt. c(i,1)) then
            xr(i) = c(i,1)
         else if (xr(i) .gt. c(i,2)) then
            xr(i) = c(i,2);
         end if
      enddo
      x(2,:) = xr;
      
      ! Compare with best and second worst points.
      if (      iso_err_fn(pts(1,:),tc13,tC,fn_flg) .le. iso_err_fn(xr       ,tc13,tC,fn_flg)             &
          .and. iso_err_fn(xr      ,tc13,tC,fn_flg) .lt. iso_err_fn(pts(2,:) ,tc13,tC,fn_flg) ) then

         pts(3,1) = xr(1)
         pts(3,2) = xr(2)

         ! Then go to step 1
         step = 1
         call nelder_mead(tc13,tC,sol,x,p,pts,c,step,cnt,fn_flg)
        
      else if (iso_err_fn(xr,tc13,tC,fn_flg) .lt. iso_err_fn(pts(1,:),tc13,tC,fn_flg)) then
         ! Then go to step 4
         step = 4
         call nelder_mead(tc13,tC,sol,x,p,pts,c,step,cnt,fn_flg)
        
      else
         ! We have f(xr) >= f(pts(end-1,:))
         ! Go to step 5
         step = 5
         call nelder_mead(tc13,tC,sol,x,p,pts,c,step,cnt,fn_flg)
      end if
   end if

   if (step .eq. 4) then
      cnt = cnt+1
      ! Unpack x as nec. (done for ease of coding / clarity)
      x0 = x(1,:)
      xr = x(2,:)

      ! Compute domain constrained expanded point:
      xe = x0 + p(2)*(x0 - pts(3,:))
      do i = 1,2
        if (xe(i) .lt. c(i,1)) then
            xe(i) = c(i,1)
        else if (xe(i) .gt. c(i,2)) then
            xe(i) = c(i,2)
        end if
      end do
      x(3,:) = xe
      
      ! Compare with reflected point:
      if (iso_err_fn(xe,tc13,tC,fn_flg) .lt. iso_err_fn(xr,tc13,tC,fn_flg)) then
         pts(3,1) = xe(1)
         pts(3,2) = xe(2)
         step = 1
         
         call nelder_mead(tc13,tC,sol,x,p,pts,c,step,cnt,fn_flg)
      else
         pts(3,1) = xr(1)
         pts(3,2) = xr(2)

         step = 1
         call nelder_mead(tc13,tC,sol,x,p,pts,c,step,cnt,fn_flg)
      end if
   end if
   
   if (step .eq. 5) then
      cnt = cnt+1
      ! Unpack x as nec. (done for ease of coding / clarity)
      x0 = x(1,:)

      ! Compute the domain constrained contracted point:
      xc = x0 + p(3)*(x0 - pts(3,:))
      do i = 1,2
        if (xc(i) .lt. c(i,1)) then
            xc(i) = c(i,1)
        else if (xc(i) .gt. c(i,2)) then
            xc(i) = c(i,2)
        end if
      end do
      x(4,:) = xc
      
      if (iso_err_fn(xc,tc13,tC,fn_flg) .lt. iso_err_fn(pts(3,:),tc13,tC,fn_flg)) then
         pts(3,1) = xc(1)
         pts(3,2) = xc(2)
         step = 1
         call nelder_mead(tc13,tC,sol,x,p,pts,c,step,cnt,fn_flg)
      else
         step = 6
         call nelder_mead(tc13,tC,sol,x,p,pts,c,step,cnt,fn_flg)
      end if
   end if

   if (step .eq. 6) then
      cnt = cnt+1;
      ! Reduce...

      i = 0
      do i = 2,3;
        pts(i,:) = pts(1,:) + p(4)*(pts(i,:) - pts(1,:));    
      end do

      ! Enforce domain constraints
      i = 0
      j = 0
      do j = 1,3
        do i = 1,2
            if (pts(j,i) .lt. c(i,1)) then
                pts(j,i) = c(i,1)
            else if (pts(j,i) .gt. c(i,2)) then
                pts(j,i) = c(i,2)
            end if
         end do
      end do
      step = 1
      call nelder_mead(tc13,tC,sol,x,p,pts,c,step,cnt,fn_flg)
   end if
   
end subroutine nelder_mead
!==========================================================================================!


!==========================================================================================!
real function iso_err_fn(m,tc13,tC,fn_flg)
   use isotopes,  only : iso_P1       & ! intent(in)
                           , iso_P2       & ! intent(in)
                           , R_std        ! ! intent(in)
   implicit none
   
   !------ Arguments. ---------------------------------------------------------------------!
   real, dimension(2),   intent(in)  :: m       ! Point from simplex (leaf,root)
   real,                 intent(in)  :: tc13    ! Total c13 being distributed
   real, dimension(3),   intent(in)  :: tC      ! Leaf, root, s-item carbon 
   integer,              intent(in)  :: fn_flg  ! Determines which function to use
   
   !------ Local Vars. --------------------------------------------------------------------!
   real                              :: s, ts   ! s- tissue 13C and total C
   real                              :: dleaf
   real                              :: droot
   real                              :: ds
   
   s    = tc13 - m(1) - m(2)
   ts   = tC(3)

   dleaf = htIsoDelta(m(1),tC(1))
   droot = htIsoDelta(m(2),tC(2))
   ds    = htIsoDelta(s   ,ts   )
   
   
   if (fn_flg == 0) then
   ! We have leaves to compare with.
   iso_err_fn = (dleaf - droot - iso_P1)**2 + (dleaf - ds - iso_P2)**2
                 
   else if (fn_flg == 1) then 
   ! There are roots but no leaves to compute deltas for.
   iso_err_fn = (ds - droot - iso_P1 + iso_P2)**2
   
   else if (fn_flg == 2) then
   ! There are leaves to compute deltas for, but there no roots.
   iso_err_fn = (dleaf - ds - iso_P2)**2

   end if
   
   
end function iso_err_fn
!==========================================================================================!


!==========================================================================================!
recursive subroutine sort(A,labels)
   real,    intent(inout), dimension(:) :: A
   real,    intent(inout), dimension(:) :: labels
   integer                              :: piv_index
 
  
   if(size(A) > 1) then
      call partition(A,piv_index,labels)
      call sort(A(:piv_index-1),labels(:piv_index-1))
      call sort(A(piv_index:),labels(piv_index:))
   endif

end subroutine sort
!==========================================================================================!


!==========================================================================================!
subroutine partition(A, piv_index, labels)
   real,    intent(inout), dimension(:)  :: A            ! vector to sort
   real,                   dimension(:)  :: labels       ! vector of labels for A
   integer, intent(inout)                :: piv_index    ! index of the pivot point
   integer                               :: i, j         ! up stepper, down stepper.
   real                                  :: temp         ! storage variable
   real                                  :: pivot        ! pivot point

   pivot = A(1)

   i  = 0
   j  = size(A) + 1

   do
      j = j-1
      ! Decrement j until A(j) <= pivot
      do
         if (A(j) <= pivot) exit
         j = j-1
      end do
      i = i+1
      
      ! Increment i until A(i) >= pivot
      do
         if (A(i) >= pivot) exit
         i = i+1
      end do

      ! Swap A(i) and A(j), or 
      if (i < j) then
         temp = A(i)
         A(i) = A(j)
         A(j) = temp
         
         temp      = labels(i)
         labels(i) = labels(j)
         labels(j) = temp
         
      elseif (i == j) then
         piv_index = i+1
         return
      else
         piv_index = i
         return
      endif
   end do

end subroutine partition
!==========================================================================================!
   
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
   
   if (total > tiny(1.0)) then
      htIsoDelta = sngl((heavy8/(total8-heavy8))/R_std - 1.d0) * 1000.0
   else
      if (heavy < tiny(1.0)) then
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
subroutine get_lhc_target(ipft,bleaf,bleaf_c13,today_gpp,today_gpp_c13,lh2tc_in,sth2tc_in  &
                         ,bleaf_in,lhc_target)
   use decomp_coms,  only  :  f_labile
   implicit none

   !---- Arguments ------------------------------------------------------------------------!
   integer, intent(in)     :: ipft
   real   , intent(in)     :: bleaf
   real   , intent(in)     :: bleaf_c13
   real   , intent(in)     :: today_gpp
   real   , intent(in)     :: today_gpp_c13
   real   , intent(in)     :: lh2tc_in             ! Pre alloc_c_bal. leaf ratio
   real   , intent(in)     :: sth2tc_in            ! Pre alloc_c_bal. storage ratio
   real   , intent(in)     :: bleaf_in
   real   , intent(inout)  :: lhc_target
   
   !---- Local Vars -----------------------------------------------------------------------!
   real                 :: lh2tc          ! Leaf heavy to total carbon
   real                 :: gpmC           ! Greatest possible amnt. leaf mixing C
   real                 :: gpmc13         ! Greatest possible amnt. leaf mixing c13
   real                 :: lpmc13         ! Least possible amnt. leaf mixing c13, UNDER THE
                                          ! ASSUMPTION of max C mixing.
   real                 :: mix            ! Mixture c13 under case(1)
   real                 :: res            ! Resident (non mixing) c13 under case(1)
   real                 :: max_bleaf_c13  ! Maximum new bleaf c13 under case(1)
   real                 :: min_bleaf_c13  ! Minimum new bleaf c13 under case(1)
   
   real                 :: max_bleaf      ! Max(bleaf,bleaf_in)
   real                 :: init_lhc_targ  ! Initial leaf heavy carbon target (gpp ratio)
   
   select case(0)
   case(0)
   !---------------------------------------------------------------------------------------!
   ! d13C Leaf = d13C Today GPP                                                            !  
   !---------------------------------------------------------------------------------------!
      if (bleaf > tiny(1.0)) then
         if (today_gpp > tiny(1.0)) then
            lh2tc = today_gpp_c13/today_gpp
         else
            if (lh2tc_in > tiny(1.0)) then
               lh2tc = lh2tc_in
            else
               lh2tc = sth2tc_in
            end if
         end if
      else
         lh2tc  = 0.0
      end if
      
      lhc_target = lh2tc *bleaf
      
   case(1) 
   !---------------------------------------------------------------------------------------!
   !     New leaf composition: Total C remains the same, some fraction of the leaf is      !
   ! 'non-mixing' on the diel scale, and some fraction mixes with assimilates before they  !
   ! get exported. Using this fraction daily we would have:                                !
   ! [[   bleaf_c13 = (1-f)*bleaf_c13 + mixed   ]] where f is the mixing fraction, 'mixed' !
   ! has ratio given by (f*bleaf_c13 + assim_c13)/(f*bleaf + assim), and total amnt given  !
   ! by f*bleaf so the mixed c13 added is just the total multiplied by the ratio.          !
   !                                                                                       !
   !     Unfortunately, this would suppose a constant bleaf which is not the case. Hence   !
   ! we instead use this eq. to constrain the maximum rate of change of leaf 13C. At this  !
   ! point bleaf is final for the day.                                                     !
   !---------------------------------------------------------------------------------------!
   ! Suppose the target is higher than the current 13C. The max. increase in leaf 13C we   !
   ! could plausibly obtain would be if somehow none of the 13C presently in the leaf was  !
   ! 'mixing' so that we would only be adding 13C to the leaf. Likewise the maximum dec.   !
   ! in leaf 13C we could obtain, if the target is lower than current 13C, would be if all !
   ! the mixing C were c13.                                                                !
   !                                                                                       !
   ! Neither of these cases may themselves be obtainable for two reasons: The mixing frac. !
   ! of the leaf may be too large to avoid having 13C in it or too small to contain none.  !
   ! However, the max and min will still be found as the 'closest' approx. of these states.!
   !---------------------------------------------------------------------------------------!
   mix = 0.0; min_bleaf_c13 = 0.0
   res = 0.0; max_bleaf_c13 = 0.0
      if (bleaf > tiny(1.0)) then
         max_bleaf = max(bleaf,bleaf_in)
         gpmC      = f_labile(ipft) * max_bleaf
         gpmc13    = min(bleaf_c13,gpmC)
         
         if (bleaf - bleaf_c13 >= gpmC) then
            lpmc13 = 0.0
         else
            lpmc13 = gpmC - (bleaf - bleaf_c13)
         end if
         
         !------------------------------------------------------------------------------------!
         ! Now the maximum plausible increase of leaf c13 is obtained using the mixing eq. w/ !
         ! the mixing C being gpmC and the mixing c13 being lpmc13. The most we can plausibly !
         ! decrease c13 is if we use gpmC as mixing C and gpmc13 as mixing c13.               !
         !------------------------------------------------------------------------------------!
         if (today_gpp > tiny(1.0)) then
            init_lhc_targ = today_gpp_c13 /today_gpp * bleaf
         else
            init_lhc_targ = lh2tc_in * bleaf
         end if
         
         if (init_lhc_targ > bleaf_c13) then
            mix           = (lpmc13 + today_gpp_c13) /(gpmC + today_gpp) * gpmC
            res           = bleaf_c13 - lpmc13
            max_bleaf_c13 = res + mix
            if (init_lhc_targ > max_bleaf_c13) then
               lhc_target = max_bleaf_c13
            else
               lhc_target = init_lhc_targ
            end if
         elseif (init_lhc_targ < bleaf_c13) then
            mix           = (gpmc13 + today_gpp_c13) /(gpmC + today_gpp) * gpmC
            res           = bleaf_c13 - gpmc13
            min_bleaf_c13 = res + mix
            if (init_lhc_targ < min_bleaf_c13) then
               lhc_target = min_bleaf_c13
            else
               lhc_target = init_lhc_targ
            end if
         else
            lhc_target = init_lhc_targ
            return
         end if
           

         
      else
         lhc_target = 0.0
      end if
      !---------------------------------------------------------------------------------------!
      ! Check that if we're adding c13 to the leaf that the addition doesn't exceed the newly !
      ! assimilated GPP c13.
      !---------------------------------------------------------------------------------------!
      if (lhc_target > today_gpp_c13 + bleaf_c13) then
         lhc_target = today_gpp_c13 + bleaf_c13
      end if
         ! t13C_t+1 = r13C + m13C
         ! m13C = (t13C_t - r13C + g13C) * f * tC_t / (tC_t + gC)
         ! hence
         ! t13C_t+1 = (r13C *(tC_t + gC - (f *tC_t)) + (t13C_t + g13C) *f *tC_t) /(tC_t + gC)         
   end select
   
   ! write(*,*) '!--------------------------------------------------------------------------!'
   ! write(*,*) ' bleaf,   bleaf_in,  blc13  :', bleaf, bleaf_in, bleaf_c13
   ! write(*,*) ' gpmC ,   gpmc13  ,  lpmc13 :', gpmC, gpmc13,lpmc13
   ! write(*,*) ''
   ! write(*,*) ' [min, max], init_lhc_t     :', min_bleaf_c13, max_bleaf_c13, init_lhc_targ
   ! write(*,*) ''
   ! write(*,*) ' mix       , res            :', mix, res
   ! write(*,*) ''
   ! write(*,*) ' lhc_target,  dbl, dbl_in   :', lhc_target, htIsoDelta(bleaf_c13,bleaf), htIsoDelta(bleaf_c13,bleaf_in)
   ! write(*,*) ''
   ! write(*,*) ' dlhc_targ ,  dinit_lhc_t   :', htIsoDelta(lhc_target,bleaf)               &
                                             ! , htIsoDelta(today_gpp_c13,today_gpp)

   !(*) NOTE:
   
end subroutine get_lhc_target
!==========================================================================================!





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
!     This sub-routine does a pretty comprehensive sanity check on C-13 variables.         !
!------------------------------------------------------------------------------------------!
subroutine check_patch_c13(cpatch,ico,call_loc,fname,dtlsm_C_gain,dtlsm_c13_gain          &
                           ,assim_h2tc,leaf_h2tc)
   use ed_max_dims    , only : str_len            ! ! intent(in)
   use ed_state_vars  , only : patchtype          ! ! structure
   use isotopes       , only : c13af              & ! intent(in)
                             , c_alloc_flg        ! ! intent(in)
   use isotopes       , only : larprop            ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(patchtype)           , intent(in) :: cpatch      ! Current patch
   integer                   , intent(in) :: ico         ! Current cohort number
   character(*)              , intent(in) :: call_loc    ! What routine called this check?
   character(*)              , intent(in) :: fname       ! What file is that routine in?
   real        ,optional     , intent(in) :: dtlsm_C_gain
   real        ,optional     , intent(in) :: dtlsm_c13_gain
   real        ,optional     , intent(in) :: assim_h2tc  ! Ratio of C-13/C in assimilate
   real        ,optional     , intent(in) :: leaf_h2tc   ! Ratio of C-13/C in leaves
   !----- Local variables. ----------------------------------------------------------------!
   logical        :: error_found = .false.               ! Is there a problem?
   character(60)  :: reason                              ! Error or warning diagnosis reason
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
   !---------------------------------------------------------------------------------------!
   Cfmt = '(11A18)'
   Rfmt = '(11E18.2)'

   !------------------------------------------------------------------------------!
   ! Perform a great big sanity check. Today_XXX vars are included to potentially !
   ! detect problems generated elsewhere in the code. Clearly it is the case in   !
   ! this routine that if a today var is messed up it is because it's sub-daily   !
   ! analogue is to blame.                                                        !
   !------------------------------------------------------------------------------!
   ! Assimilated C-13 should not exceed assimilated C.
   call check_c13(cpatch%gpp_c13(ico),cpatch%gpp(ico),gpp_delta,valid)
   if ( .not. valid) then
      reason = 'GPP C-13 too high.'
      error_found = .true.
   end if
   
   call check_c13(cpatch%leaf_respiration_c13(ico) &
                 ,cpatch%leaf_respiration    (ico) &
                 ,lr_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in leaf respiration...'
      error_found = .true.
   end if
   
   call check_c13(cpatch%root_respiration_c13(ico) &
                 ,cpatch%root_respiration    (ico) &
                 ,rr_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in root respiration...'
      error_found = .true.
   end if
   
   call check_c13(cpatch%leaf_growth_resp_c13(ico) &
                 ,cpatch%leaf_growth_resp    (ico) &
                 ,lg_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in leaf growth respiration...'
      error_found = .true.
   end if
   
   call check_c13(cpatch%root_growth_resp_c13(ico) &
                 ,cpatch%root_growth_resp    (ico) &
                 ,rg_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in root growth respiration...'
      error_found = .true.
   end if
   
   call check_c13(cpatch%sapa_growth_resp_c13(ico) &
                 ,cpatch%sapa_growth_resp    (ico) &
                 ,sag_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in sapa growth respiration...'
      error_found = .true.
   end if
   
   call check_c13(cpatch%sapb_growth_resp_c13(ico) &
                 ,cpatch%sapb_growth_resp    (ico) &
                 ,sbg_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in sapb growth respiration...'
      error_found = .true.
   end if
   
   call check_c13(cpatch%leaf_storage_resp_c13(ico) &
                 ,cpatch%leaf_storage_resp    (ico) &
                 ,ls_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in leaf storage respiration...'
      error_found = .true.
   end if
   
   call check_c13(cpatch%root_storage_resp_c13(ico) &
                 ,cpatch%root_storage_resp    (ico) &
                 ,rs_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in root storage respiration...'
      error_found = .true.
   end if
   
   call check_c13(cpatch%sapa_storage_resp_c13(ico) &
                 ,cpatch%sapa_storage_resp    (ico) &
                 ,sas_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in sapa storage respiration...'
      error_found = .true.
   end if
   
   call check_c13(cpatch%sapb_storage_resp_c13(ico) &
                 ,cpatch%sapb_storage_resp    (ico) &
                 ,sbs_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in sapb storage respiration...'
      error_found = .true.
   end if
   
   call check_c13(cpatch%bleaf_c13(ico) &
                 ,cpatch%bleaf    (ico) &
                 ,leaf_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in leaves...'
      error_found = .true.
   end if
   
   call check_c13(cpatch%broot_c13(ico) &
                 ,cpatch%broot    (ico) &
                 ,root_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in roots...'
      error_found = .true.
   end if
   
   call check_c13(cpatch%bsapwooda_c13(ico) &
                 ,cpatch%bsapwooda    (ico) &
                 ,bsa_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in aboveground sapwood...'
      error_found = .true.
   end if
   
   call check_c13(cpatch%bsapwoodb_c13(ico) &
                 ,cpatch%bsapwoodb    (ico) &
                 ,bsb_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in belowground sapwood...'
      error_found = .true.
   end if
   
   call check_c13(cpatch%bstorage_c13(ico) &
                 ,cpatch%bstorage    (ico) &
                 ,bst_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in storage...'
      error_found = .true.
   end if
   
   non_stor_resp = cpatch%leaf_respiration(ico) + cpatch%root_respiration(ico)        &
                 + cpatch%leaf_growth_resp(ico) + cpatch%root_growth_resp(ico)        &
                 + cpatch%sapa_growth_resp(ico) + cpatch%sapb_growth_resp(ico)
                 
   nsr_c13 = cpatch%leaf_respiration_c13(ico) + cpatch%root_respiration_c13(ico)      &
           + cpatch%leaf_growth_resp_c13(ico) + cpatch%root_growth_resp_c13(ico)      &
           + cpatch%sapa_growth_resp_c13(ico) + cpatch%sapb_growth_resp_c13(ico)
   
   if (error_found) then
      write(*,*) '======================================================================='
      write(*,*) ' C-13 sanity check error in ', call_loc, '!'
      write(*,*) ' ', reason
      write(*,*) '======================================================================='
      write(*,*) ' cohort, pft, larprop : ', ico, cpatch%pft(ico), larprop
      write(*,*) ' '
      write(*,*) '-----------------------------------------------------------------------'
      write(*,*) 'Flux Diagnostics. Row 1: Totals. Row 2: Heavy C. Row 3: d13C.'
      write(*,*) '-----------------------------------------------------------------------'
      write(*,'(3A18)') '               GPP'          &
                         ,'  LEAF_RESPIRATION'        &
                         ,'  ROOT_RESPIRATION'        &
                         ,'        **RESP_SUM'
      write(*,'(3E18.2)')  cpatch%gpp(ico)                &
                          ,cpatch%leaf_respiration(ico)   &
                          ,cpatch%root_respiration(ico)   &
                          ,non_stor_resp
      write(*,'(3E18.2)')  cpatch%gpp_c13(ico)               &
                          ,cpatch%leaf_respiration_c13(ico)  &
                          ,cpatch%root_respiration_c13(ico)  &
                          ,nsr_c13
      write(*,'(3E18.2)') gpp_delta,lr_delta,rr_delta, hotc(nsr_c13,non_stor_resp)
      write(*,*) ' '
      write(*,'(4A18)') '  LEAF_GROWTH_RESP','  ROOT_GROWTH_RESP','  SAPA_GROWTH_RESP' &
                       ,'  SAPB_GROWTH_RESP'
      write(*,'(4E18.2)') cpatch%leaf_growth_resp(ico) ,cpatch%root_growth_resp(ico)       &
                       ,cpatch%sapa_growth_resp(ico) ,cpatch%sapb_growth_resp(ico)
      write(*,'(4E18.2)') cpatch%leaf_growth_resp_c13(ico)    &
                       ,cpatch%root_growth_resp_c13(ico)    &
                       ,cpatch%sapa_growth_resp_c13(ico)    &
                       ,cpatch%sapb_growth_resp_c13(ico)
      write(*,'(4E18.2)') lg_delta,rg_delta,sag_delta,sbg_delta
      write(*,*) ' '
      write(*,'(4A18)') ' LEAF_STORAGE_RESP',' ROOT_STORAGE_RESP',' SAPA_STORAGE_RESP' &
                         ,' SAPB_STORAGE_RESP'
      write(*,'(4E18.2)') cpatch%leaf_storage_resp(ico)       &
                       ,cpatch%root_storage_resp(ico)       &
                       ,cpatch%sapa_storage_resp(ico)       &
                       ,cpatch%sapb_storage_resp(ico)
      write(*,'(4E18.2)') cpatch%leaf_storage_resp_c13(ico)   &
                       ,cpatch%root_storage_resp_c13(ico)   &
                       ,cpatch%sapa_storage_resp_c13(ico)   &
                       ,cpatch%sapb_storage_resp_c13(ico)
      write(*,'(4E18.2)') ls_delta,rs_delta,sas_delta,sbs_delta

      write(*,*) ' **RESP_SUM is the sum of non-storage respiration terms.'
      write(*,*) ' '
      write(*,*) '-----------------------------------------------------------------------'
      write(*,*) 'Pool Diagnostics. Row 1: Totals. Row 2: Heavy C. Row 3: d13C.'
      write(*,*) '-----------------------------------------------------------------------'
      write(*,'(4A10)') '     BLEAF','     BROOT',' BSAPWOODA','BSAPWOODB','  BSTORAGE'
      write(*,'(4E10.2)') cpatch%bleaf(ico),cpatch%broot(ico),cpatch%bsapwooda(ico)      &
                         ,cpatch%bsapwoodb(ico), cpatch%bstorage(ico)
      write(*,'(5E10.2)') cpatch%bleaf_c13(ico),cpatch%broot_c13(ico)                    &
                         ,cpatch%bsapwooda_c13(ico),cpatch%bsapwoodb_c13(ico)            &
                         ,cpatch%bstorage_c13(ico)

      write(*,'(5E10.2)') leaf_delta,root_delta,bsa_delta,bsb_delta,bst_delta

      write(*,*) ' '
      write(*,*) ' '
      write(*,*) '-----------------------------------------------------------------------'
      if (present(assim_h2tc)) then
         write(*,*) '-----------------------------------------------------------------------'
         write(*,*) ' Pieces of leaf C-13 respiration computation'
         write(*,Cfmt) 'R Leaf Tissue', 'R Leaf Assim',                        &
                       'Leaf C-13/C'  , '"leaf_h2tc" '
         write(*,Rfmt) cpatch%leaf_respiration(ico) - cpatch%lassim_resp(ico), &
                       cpatch%lassim_resp     (ico)                          , &
                       cpatch%bleaf_c13       (ico) / cpatch%bleaf      (ico), &
                       leaf_h2tc
         write(*,*)
         write(*,*) ' Leaf assimilate C-13/C ratio and recalculation from GPP and GPP_C13'
         write(*,*) ' If this msg is not from canopy_photosynthesis, ignore them.        '
         write(*,*) ' assim_h2tc : ', assim_h2tc, cpatch%gpp_c13(ico) / cpatch%gpp(ico)
         write(*,*) '-----------------------------------------------------------------------'
      end if
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
subroutine check_site_c13(csite,ipa,call_loc,fname)
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
   !----- Local variables. ----------------------------------------------------------------!
   logical        :: error_found = .false.               ! Is there a problem?
   character(60)  :: reason                              ! Error or warning diagnosis reason
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
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   Cfmt = '(6A18)'
   Rfmt = '(6E18.2)'

   call check_c13(csite%fast_soil_c13(ipa) &
                 ,csite%fast_soil_c  (ipa) &
                 ,fsc_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in fast_soil_c...'
      error_found = .true.
   end if
   
   call check_c13(csite%slow_soil_c13(ipa) &
                 ,csite%slow_soil_c  (ipa) &
                 ,ssc_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in slow_soil_c...'
      error_found = .true.
   end if
   
   call check_c13(csite%structural_soil_c13(ipa) &
                 ,csite%structural_soil_c  (ipa) &
                 ,stsc_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in structural_soil_c...'
      error_found = .true.
   end if
   
   call check_c13(csite%structural_soil_L_c13(ipa) &
                 ,csite%structural_soil_L    (ipa) &
                 ,stsl_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in structural_soil_L...'
      error_found = .true.
   end if
   
   call check_c13(csite%can_co2_c13(ipa) &
                 ,csite%can_co2    (ipa) &
                 ,can_co2_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in can_co2...'
      error_found = .true.
   end if
   
   call check_c13(csite%rh_c13(ipa) &
                 ,csite%rh    (ipa) &
                 ,rh_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in heterotrophic respiration...'
      error_found = .true.
   end if
   
   call check_c13(csite%cwd_rh_c13(ipa) &
                 ,csite%cwd_rh    (ipa) &
                 ,cwd_delta, valid)
   if ( .not. valid) then
      reason = 'There is too much C-13 in heterotrophic respiration...'
      error_found = .true.
   end if
   
   if (error_found) then
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
                   ,' STRUCTURAL_SOIL_L','           CAN_CO2'
      write(*,Rfmt) csite%fast_soil_C(ipa),csite%slow_soil_C(ipa)                        &
                   ,csite%structural_soil_C(ipa),csite%structural_soil_L(ipa)            &
                   ,csite%can_co2(ipa)
      write(*,Rfmt) csite%fast_soil_c13(ipa),csite%slow_soil_c13(ipa)                    &
                   ,csite%structural_soil_c13(ipa),csite%structural_soil_L_c13(ipa)      &
                   ,csite%can_co2_c13(ipa)
      write(*,Rfmt) fsc_delta, ssc_delta, stsc_delta, stsl_delta, can_co2_delta

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
subroutine check_c13(heavy,total,delta,valid)
   use consts_coms, only : tiny_num    ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real   , intent(in)  :: heavy
   real   , intent(in)  :: total
   real   , intent(out) :: delta
   logical, intent(out) :: valid
   !---------------------------------------------------------------------------------------!
   
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
   
   delta = htIsoDelta(heavy,total)

end subroutine check_c13
!==========================================================================================!
!==========================================================================================!



!==========================================================================================!
!==========================================================================================!
! This subroutine computes the leaf and root respiration 13C composition.                  !
!------------------------------------------------------------------------------------------!
subroutine leaf_root_resp_c13(csite,ipa)
   use ed_state_vars  , only : sitetype                  & ! structure
                             , patchtype                 ! ! structure
   use ed_misc_coms   , only : dtlsm                     & ! intent(in)
                             , frqsum                    ! ! intent(in)
   use isotopes       , only : c13af                     ! ! intent(in)
   use consts_coms    , only : tiny_num                  & ! intent(in)
                             , umols_2_kgCyr             ! ! intent(in)
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
         
         !if (non_stor_resp <= cpatch%gpp(ico)) then
         !   ratio_resp = ratio_gpp
         !else
         !   ratio_stor  = hotc(cpatch%bstorage_c13(ico),cpatch%bstorage(ico))
         !   c_from_gpp  = cpatch%gpp(ico)
         !   c_from_stor = min(non_stor_resp - c_from_gpp,cpatch%bstorage(ico))
         !   ratio_resp  = (cpatch%gpp(ico) - non_stor_resp)/non_stor_resp
         !end if
         
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

end subroutine leaf_root_resp_c13
!==========================================================================================!
!==========================================================================================!


end module iso_alloc
