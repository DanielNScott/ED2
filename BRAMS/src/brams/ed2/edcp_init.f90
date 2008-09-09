subroutine node_ed_init
  
  
  ! BRAMS MODULES
  use node_mod,only: master_num,mchnum,mynum,nmachs,machs
  use mem_grid,only: ngrids
  use rpara,only: mainnum
  use soil_coms,only: layer_index,nlon_lyr,nlat_lyr
  use mem_leaf,only:isfcl
  ! ED MODULES
  use ed_state_vars,only:  &
       gdpy,   &
       py_off, &
       allocate_edglobals, &
       allocate_edtype,    &
       edgrid_g

  implicit none
  include 'mpif.h'
  integer :: ifm,ierr
  
  ! First we must transfer over the parallel information
  ! from BRAMS to ED2
  
  if(isfcl.ne.5)return

  call copy_in_bramsmpi(master_num,mchnum,mynum,nmachs,machs)
  

  ! Calculate the polygon list on the current node
  call init_node_work

  ! Recieve the lowest soil layer index array
  if (allocated(layer_index)) deallocate(layer_index)
  allocate(layer_index(nlat_lyr,nlon_lyr))
  call MPI_Barrier(MPI_COMM_WORLD,ierr) ! Safe to receive the data.
  call MPI_Bcast(layer_index,nlat_lyr*nlon_lyr,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  
  ! Allocate the polygons on edgrid
  write (unit=*,fmt='(a,i5,a)') ' + Polygon array allocation, node ',mynum,';'
  call allocate_edglobals(ngrids)
  do ifm=1,ngrids
     call ed_newgrid(ifm)
     call allocate_edtype(edgrid_g(ifm),gdpy(mynum,ifm))
  end do
  
  write (unit=*,fmt='(a,i5,a)') ' + Memory successfully allocated on none ',mynum,';'

  call ed_coup_driver()


  return
end subroutine node_ed_init

! ===========================================================================================

subroutine master_ed_init()

  use rpara,only: mainnum
  use soil_coms,only: layer_index,nlon_lyr,nlat_lyr
  use mem_leaf,only:isfcl
  implicit none
  include 'mpif.h'
  integer :: ifm,ierr

  if(isfcl.ne.5)return

  ! Initialize the work arrays
  
  call init_master_work()

  ! Read the soil depth database
  call read_soil_depth()

  ! Send soil depths to the nodes
  call MPI_Barrier(MPI_COMM_WORLD,ierr) ! Just to wait until the matrix is allocated
  call MPI_Bcast(layer_index,nlat_lyr*nlon_lyr,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  


  return
end subroutine master_ed_init

! ===========================================================================================

subroutine copy_in_bramsmpi(master_num_b,mchnum_b,mynum_b,nmachs_b,machs_b)

  use ed_node_coms,only: master_num,mchnum,mynum,nmachs,machs,nnodetot,recvnum,sendnum
  use ed_para_coms,only: iparallel,mainnum
  use max_dims,only:maxgrds
  use grid_coms, only : ngrids,nnxp,nnyp

  implicit none
  include 'mpif.h'
  integer :: ierr
  integer :: master_num_b
  integer :: mchnum_b
  integer :: mynum_b
  integer :: nmachs_b
  integer :: igr,ngrids_b
  integer,dimension(maxgrds) :: nnxp_b,nnyp_b
  integer,dimension(nmachs_b) :: machs_b

  ! If we are calling this subroutine, it is already a paraellel run
  
  iparallel = 1
  
  master_num = master_num_b
  mchnum     = mchnum_b
  mynum      = mynum_b
  nmachs     = nmachs_b
  machs      = machs_b
  nnodetot   = nmachs_b

  recvnum = mynum-1
  sendnum = mynum+1

  return
end subroutine copy_in_bramsmpi

!==========================================================================================!

subroutine init_node_work
  
  
  use ed_work_vars       , only : work_e,                 & ! intent(out)
       ed_alloc_work,          & ! subroutine
       ed_nullify_work         ! ! subroutine
  
  use grid_coms, only : ngrids

  use ed_node_coms, only: mynum,nmachs,nnodetot,mchnum,machs,master_num,sendnum,recvnum,mmxp,mmyp

  use rpara,only : mainnum
  
  use ed_state_vars,only:  gdpy,   &
       py_off, &
       allocate_edglobals, &
       allocate_edtype,    &
       edgrid_g

  implicit none
  include 'mpif.h'
  integer :: ifm,nm,nm2,ierr,i,j
  integer,       dimension(MPI_STATUS_SIZE) :: status

  allocate(work_e(ngrids))
  do ifm = 1,ngrids
          
     do nm=1,nmachs 

        call MPI_Bcast(nm2,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
        
        if(nm2/=nm) then
           print*,"THE BROADCAST NODE DOES NOT MATCH THE LOOP NODE"
           stop
        endif

        call MPI_Bcast(gdpy(nm,ifm),1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(py_off(nm,ifm),1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)


        if(nm==mynum) then
           
           call MPI_Recv(mmxp(ifm),1,MPI_INTEGER,mainnum,190,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(mmyp(ifm),1,MPI_INTEGER,mainnum,191,MPI_COMM_WORLD,status,ierr)
           
           call ed_nullify_work(work_e(ifm))
           call ed_alloc_work(work_e(ifm),mmxp(ifm),mmyp(ifm))

           call MPI_Recv(work_e(ifm)%glat,mmxp(ifm)*mmyp(ifm),MPI_REAL,mainnum,192,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(work_e(ifm)%glon,mmxp(ifm)*mmyp(ifm),MPI_REAL,mainnum,193,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(work_e(ifm)%work,mmxp(ifm)*mmyp(ifm),MPI_REAL,mainnum,194,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(work_e(ifm)%land,mmxp(ifm)*mmyp(ifm),MPI_LOGICAL,mainnum,195,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(work_e(ifm)%landfrac,mmxp(ifm)*mmyp(ifm),MPI_REAL,mainnum,196,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(work_e(ifm)%ntext,mmxp(ifm)*mmyp(ifm),MPI_INTEGER,mainnum,197,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(work_e(ifm)%xatm,mmxp(ifm)*mmyp(ifm),MPI_INTEGER,mainnum,198,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(work_e(ifm)%yatm,mmxp(ifm)*mmyp(ifm),MPI_INTEGER,mainnum,199,MPI_COMM_WORLD,status,ierr)

        endif

     enddo
  enddo
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
        
end subroutine init_node_work


!==========================================================================================!

subroutine init_master_work()
  
  use mem_grid,only: grid_g,ngrids
  use rpara,only: ixb,ixe,iyb,iye,nmachs,mainnum,machnum
  use ed_work_vars,only : &
       work_e,                 & ! intent(out)
       ed_alloc_work,          & ! subroutine
       ed_nullify_work,        &         ! ! subroutine
       ed_dealloc_work

  use ed_state_vars,only:  &
       gdpy,   &
       py_off, &
       allocate_edglobals, &
       allocate_edtype,    &
       edgrid_g

  implicit none
  include 'mpif.h'
  integer :: ierr
  integer :: nm,ifm
  integer :: offset,npolys
  integer :: nxp,nyp
  integer :: i,j,il,jl
  
  allocate(work_e(ngrids))
  
  do ifm = 1,ngrids
     
     npolys=0
     offset=0
     do nm=1,nmachs 

        nxp = ixe(nm,ifm)-ixb(nm,ifm)+1
        nyp = iye(nm,ifm)-iyb(nm,ifm)+1
        
        call ed_nullify_work(work_e(ifm))
        call ed_alloc_work(work_e(ifm),nxp,nyp)

        il=0
        do i=ixb(nm,ifm),ixe(nm,ifm)
           jl=0
           il=il+1
           do j=iyb(nm,ifm),iye(nm,ifm)
              jl=jl+1
              work_e(ifm)%glon(il,jl) = grid_g(ifm)%glon(i,j)
              work_e(ifm)%glat(il,jl) = grid_g(ifm)%glat(i,j)
              work_e(ifm)%xatm(il,jl)  = i - ixb(nm,ifm) + 2       ! Remember that all tiles have a 
              work_e(ifm)%yatm(il,jl)  = j - iyb(nm,ifm) + 2       ! buffer cell so add one extra 
           end do
        enddo
        
        call get_work(ifm,nxp,nyp)
        
        offset=offset+npolys
        npolys=count(work_e(ifm)%land)
        gdpy(nm,ifm)   = npolys
        py_off(nm,ifm) = offset
        call MPI_Bcast(nm,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(npolys,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(offset,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

        ! Send the work grids to the nodes

        call MPI_Send(nxp,1,MPI_INTEGER,machnum(nm),190,MPI_COMM_WORLD,ierr)
        call MPI_Send(nyp,1,MPI_INTEGER,machnum(nm),191,MPI_COMM_WORLD,ierr)

        call MPI_Send(work_e(ifm)%glat,nxp*nyp,MPI_REAL,machnum(nm),192,MPI_COMM_WORLD,ierr)
        call MPI_Send(work_e(ifm)%glon,nxp*nyp,MPI_REAL,machnum(nm),193,MPI_COMM_WORLD,ierr)
        call MPI_Send(work_e(ifm)%work,nxp*nyp,MPI_REAL,machnum(nm),194,MPI_COMM_WORLD,ierr)
        call MPI_Send(work_e(ifm)%land,nxp*nyp,MPI_LOGICAL,machnum(nm),195,MPI_COMM_WORLD,ierr)
        call MPI_Send(work_e(ifm)%landfrac,nxp*nyp,MPI_REAL,machnum(nm),196,MPI_COMM_WORLD,ierr)
        call MPI_Send(work_e(ifm)%ntext,nxp*nyp,MPI_INTEGER,machnum(nm),197,MPI_COMM_WORLD,ierr)
        call MPI_Send(work_e(ifm)%xatm,nxp*nyp,MPI_INTEGER,machnum(nm),198,MPI_COMM_WORLD,ierr)
        call MPI_Send(work_e(ifm)%yatm,nxp*nyp,MPI_INTEGER,machnum(nm),199,MPI_COMM_WORLD,ierr)
        

        call ed_dealloc_work(work_e(ifm))

     end do
     
      
  enddo
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  return
end subroutine init_master_work