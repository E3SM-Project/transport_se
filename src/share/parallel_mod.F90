#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
module parallel_mod
  ! ---------------------------
  use kinds, only : real_kind, int_kind, iulog
  ! ---------------------------
  use dimensions_mod, only : nmpi_per_node, nlev, qsize_d


  implicit none

  public 

#include <mpif.h>

  integer, parameter, public :: ORDERED = 1
  integer, parameter, public :: FAST = 2
  integer,      public            :: MaxNumberFrames, numframes
  integer,      public            :: useframes 
  logical,      public            :: PartitionForNodes,PartitionForFrames
  integer,      public :: MPIreal_t,MPIinteger_t,MPIChar_t,MPILogical_t
  integer,      public :: iam

  integer,      public, allocatable    :: status(:,:)
  integer,      public, allocatable    :: Rrequest(:)
  integer,      public, allocatable    :: Srequest(:)

  real(kind=4), public,allocatable :: FrameWeight(:)
  integer,      public,allocatable :: FrameIndex(:)
  integer,      public,allocatable :: FrameCount(:)

  ! ==================================================
  ! Define type parallel_t for distributed memory info
  ! ==================================================

  integer, parameter :: ncomponents=1
  integer,public     :: nComPoints,nPackPoints

  type, public :: parallel_t
    integer :: rank                       ! local rank
    integer :: root                       ! local root
    integer :: nprocs                     ! number of processes in group
    integer :: comm                       ! local communicator
    integer :: intercomm(0:ncomponents-1) ! inter communicator list
    logical :: masterproc                 
  end type

  integer, parameter :: nrepro_vars=MAX(10,nlev*qsize_d)
  real(kind=8), public, allocatable :: global_shared_buf(:,:)
  real(kind=8), public :: global_shared_sum(nrepro_vars)

  ! ===================================================
  ! Module Interfaces
  ! ===================================================

  interface assignment ( = )
    module procedure copy_par
  end interface



  public :: initmp
  public :: haltmp
  public :: abortmp
  public :: split
  public :: connect
  public :: syncmp
  public :: psum_1d
  public :: pmax_1d,pmin_1d

contains

! ================================================
!   copy_par: copy constructor for parallel_t type
!
!
!   Overload assignment operator for parallel_t
! ================================================

  subroutine copy_par(par2,par1)
    type(parallel_t), intent(out) :: par2
    type(parallel_t), intent(in)  :: par1

    par2%rank       = par1%rank
    par2%root       = par1%root
    par2%nprocs     = par1%nprocs
    par2%comm       = par1%comm
    par2%intercomm  = par1%intercomm
    par2%masterproc = par1%masterproc

  end subroutine copy_par

! ================================================
!  initmp:
!  Initializes the parallel (message passing)
!  environment, returns a parallel_t structure..
! ================================================
     
  function initmp(npes_in) result(par)

    integer, intent(in), optional ::  npes_in
    type (parallel_t) par

#include <mpif.h>
#ifdef _AIX
    integer(kind=int_kind)                              :: ii         
    character(len=2)                                    :: cfn
#endif

    integer(kind=int_kind)                              :: ierr,tmp
    integer(kind=int_kind)                              :: FrameNumber
    logical :: running   ! state of MPI at beginning of initmp call
    character(len=MPI_MAX_PROCESSOR_NAME)               :: my_name
    character(len=MPI_MAX_PROCESSOR_NAME), allocatable  :: the_names(:)

    integer(kind=int_kind),allocatable                  :: tarray(:)
    integer(kind=int_kind)                              :: namelen,i

    !================================================
    !     Basic MPI initialization
    ! ================================================

    call MPI_initialized(running,ierr)

    if (.not.running) then
       call MPI_init(ierr)
    end if

    par%root      = 0
    par%intercomm = 0
    par%comm      = MPI_COMM_WORLD

    call MPI_comm_rank(par%comm,par%rank,ierr)
    call MPI_comm_size(par%comm,par%nprocs,ierr)

    par%masterproc = .FALSE.
    if(par%rank .eq. par%root) par%masterproc = .TRUE.
    if (par%masterproc) write(iulog,*)'number of MPI processes: ',par%nprocs
           
    if (MPI_DOUBLE_PRECISION==20 .and. MPI_REAL8==18) then
       ! LAM MPI defined MPI_REAL8 differently from MPI_DOUBLE_PRECISION
       ! and LAM MPI's allreduce does not accept on MPI_REAL8
       MPIreal_t    = MPI_DOUBLE_PRECISION
    else
       MPIreal_t    = MPI_REAL8
    endif
    MPIinteger_t = MPI_INTEGER
    MPIchar_t    = MPI_CHARACTER 
    MPILogical_t = MPI_LOGICAL

    ! ================================================ 
    !  Determine where this MPI process is running 
    !   then use this information to determined the 
    !   number of MPI processes per node    
    ! ================================================ 

    my_name(:) = ''
    call MPI_Get_Processor_Name(my_name,namelen,ierr)

    allocate(the_names(par%nprocs))
    do i=1,par%nprocs
       the_names(i)(:) =  ''
    enddo
    ! ================================================ 
    !   Collect all the machine names 
    ! ================================================ 
    call MPI_Allgather(my_name,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER, &
           the_names,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,par%comm,ierr)

    ! ======================================================================
    !   Calculate how many other MPI processes are on my node 
    ! ======================================================================
    nmpi_per_node = 0
    do i=1,par%nprocs
      if( TRIM(ADJUSTL(my_name)) .eq. TRIM(ADJUSTL(the_names(i)))   ) then 
        nmpi_per_node = nmpi_per_node + 1
      endif
    enddo

    ! =======================================================================
    !  Verify that everybody agrees on this number otherwise do not do 
    !  the multi-level partitioning
    ! =======================================================================
    call MPI_Allreduce(nmpi_per_node,tmp,1,MPIinteger_t,MPI_BAND,par%comm,ierr)
    if(tmp .ne. nmpi_per_node) then 
      write(iulog,*)'initmp:  disagrement accross nodes for nmpi_per_node'
      nmpi_per_node = 1
      PartitionForNodes=.FALSE.
    else
      PartitionForNodes=.TRUE.
    endif


#ifdef _AIX
    PartitionForFrames=.FALSE.
    if((my_name(1:2) .eq. 'bv') .or.  &   ! Bluvista
       (my_name(1:2) .eq. 'bl') .or.  &    ! Blueice
       (my_name(1:2) .eq. 'bh') .or.  &    ! Blue Horizon
       (my_name(1:2) .eq. 'bs') .or.  &    ! BlueSky
       (my_name(1:2) .eq. 's0')       &    ! Seaborg
      ) then

       ! ================================================================================
       ! Note: the frame based optimization is only supported on blackforest or babyblue
       ! ================================================================================
       cfn = my_name(3:4)
       read(cfn,'(i2)') FrameNumber

       ! ======================================================
       ! Make sure that the system does not have too may frames 
       ! ======================================================
       call MPI_Allreduce(FrameNumber,MaxNumberFrames,1,MPIinteger_t,MPI_MAX,par%comm,ierr)
       MaxNumberFrames=MaxNumberFrames+1
       allocate(FrameCount(MaxNumberFrames))
       allocate(tarray(MaxNumberFrames))

       call MPI_Allreduce(useframes,tmp,1,MPIinteger_t,MPI_BAND,par%comm,ierr)
       if(tmp .ne. useframes) then 
          write(iulog,*) "initmp:  disagreement accross nodes for useframes"
          PartitionForFrames=.FALSE.
       endif

       if(PartitionForFrames) then 
         tarray(:) = 0
         tarray(FrameNumber+1) = 1
         call MPI_Allreduce(tarray,FrameCount,MaxNumberFrames,MPIinteger_t,MPI_SUM,par%comm,ierr)
         if(par%masterproc)  write(iulog,*)'initmp: FrameCount : ',FrameCount
           numFrames = COUNT(FrameCount .ne. 0)
           allocate(FrameWeight(numFrames))
           allocate(FrameIndex(numFrames))

         ii=1
         do i=1,MaxNumberFrames
           if(FrameCount(i) .ne. 0) then 
             FrameWeight(ii) = real(FrameCount(i),kind=4)/real(par%nprocs,kind=4)
             FrameIndex(ii)  = i
             ii=ii+1
           endif
         enddo
       endif

      FrameCount(:)=FrameCount(:)/nmpi_per_node
      ! ==========================================
      ! We are not running on more than one frame 
      ! ==========================================
      if(numFrames .eq. 1)  PartitionForFrames=.FALSE.

    endif

    write(iulog,*) 'initmp: mpi task ',par%rank,': ',nmpi_per_node,' task(s) on node ',my_name(1:namelen),  &
                   'on frame # ',FrameNumber 
#endif
    if(PartitionForFrames) then 
      if(par%masterproc) write(iulog,*)'initmp: FrameWeight: ',FrameWeight
      if(par%masterproc) write(iulog,*)'initmp: FrameIndex: ',FrameIndex
    endif

    deallocate(the_names)

    !===================================================
    !  Kind of lame but set this variable to be 1 based 
    !===================================================
    iam = par%rank+1

  end function initmp

  ! =========================================================
  ! abortmp:
  !
  ! Tries to abort the parallel (message passing) environment
  ! and prints a message
  ! =========================================================
  subroutine abortmp(string)

    integer info,ierr

    character*(*) string

    write(*,*) iam,' ABORTING WITH ERROR: ',string
#ifdef _AIX
    call xl__trbk()
#endif
    call MPI_Abort(MPI_COMM_WORLD,info,ierr)
    call MPI_finalize(info)

  end subroutine abortmp
       
  ! =========================================================
  ! haltmp:
  !
  !> stops the parallel (message passing) environment 
  !! and prints a message.
  !
  !> Print the message and call MPI_finalize. 
  !! @param[in] string The message to be printed.
  ! =========================================================
  subroutine haltmp(string)
         
  integer info

  character*(*) string
  if(iam .eq. 1) then 
    write(*,*) string
  endif

  call MPI_finalize(info)

  ! This can send a non-zero error code to the shell
  stop
end subroutine haltmp

  ! =========================================================
  ! split:
  !
  ! splits the message passing world into components
  ! and returns a new parallel structure for the
  ! component resident at this process, i.e. lcl_component
  ! =========================================================
  function split(par,leader,lcl_component) result(newpar)

    type (parallel_t)  :: par
    type (parallel_t)  :: newpar
    integer            :: lcl_component
    integer            :: leader(0:ncomponents-1)

    integer ierr
    integer info
    integer            :: key

    lcl_component=ncomponents-1
    do while(leader(lcl_component) > par%rank)
      lcl_component=lcl_component-1
    end do

    key=par%rank   ! simplest key for most cases

    call MPI_comm_split(par%comm, lcl_component, key, newpar%comm,ierr);

    call MPI_comm_rank(newpar%comm,newpar%rank,info)
    call MPI_comm_size(newpar%comm,newpar%nprocs,info)
    newpar%root=0

  end function split

  ! =========================================================
  ! connect:
  !
  ! connects this MPI component to all others by constructing
  ! intercommunicator array and storing it in the local parallel
  ! structure lcl_par. Connect assumes you have called split
  ! to create the lcl_par structure.
  !
  ! =========================================================
  subroutine connect(gbl_par, lcl_par, lcl_component, leader)

    type (parallel_t) :: gbl_par
    type (parallel_t) :: lcl_par
    integer           :: lcl_component
    integer           :: leader(0:ncomponents-1) ! leader rank in bridge group

    integer tag
    integer i
    integer ierr

    do i=0,ncomponents-1

      if (i > lcl_component) then
        tag=ncomponents*lcl_component + i
      else
        tag=ncomponents*i+lcl_component
      end if

      if (i .ne. lcl_component) then
#ifdef _DEBUG
        write(iulog,10) lcl_component,
     &                  gbl_par%rank,
     &                  leader(lcl_component),
     &                  leader(i),
                        tag
10      format("component=",i4, 
     &         " gbl rank =",i4,   
     &         " lcl leader=",i4,  
     &         " rem leader=",i4, 
     &         " tag=",i4)
#endif
        call MPI_Intercomm_create(lcl_par%comm, lcl_par%root, gbl_par%comm, &
                                  leader(i), tag, lcl_par%intercomm(i), ierr)  
      end if
    end do     
  end subroutine connect

! =====================================
! syncmp:
! 
! sychronize message passing domains 
!
! =====================================
  subroutine syncmp(par)

    type (parallel_t) par

#include <mpif.h>
    integer                         :: errorcode,errorlen,ierr
    character(len=MPI_MAX_ERROR_STRING)               :: errorstring

    call MPI_barrier(par%comm,ierr)

    if(ierr.eq.MPI_ERROR) then
      errorcode=ierr
      call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
      call abortmp(errorstring)
    endif

  end subroutine syncmp

#if 0
  ! =============================================
  ! psum_1d:
  ! 1D version of the parallel SUM
  ! =============================================
  function psum_1d(variable,type,par) result(res)
    implicit none
    ! ==========================
    !     Arguments   
    ! ==========================
    real(kind=real_kind),intent(in)  :: variable(:)
    integer,intent(in)               :: type 
    type (parallel_t),intent(in)     :: par

    ! ==========================
    !       Local Variables 
    ! ==========================
    real(kind=real_kind)             :: res
    real(kind=real_kind)             :: local_sum
    !
    ! Note this is a real kludge here since it may be used for 
    !  arrays of size other then nelem
    ! 
    real(kind=real_kind),allocatable :: Global(:),buffer(:)
    integer                          :: ierr,i,ie,ig,disp,nelemr,ip

    local_sum=SUM(variable)

    call MPI_Allreduce(local_sum,res,1,MPIreal_t, MPI_SUM,par%comm,ierr)

  end function psum_1d
#endif
  
  ! =============================================
  ! pmin_1d:
  ! 1D version of the parallel MIN
  ! =============================================
  function pmin_1d(variable,par) result(res)

    implicit none

    real(kind=real_kind),intent(in)  :: variable(:)
    type (parallel_t),intent(in)     :: par
    real(kind=real_kind)             :: res
         
    real(kind=real_kind)             :: local_sum
    integer                          :: ierr

    local_sum=MINVAL(variable)

    call MPI_Allreduce(local_sum,res,1,MPIreal_t, &
                       MPI_MIN,par%comm,ierr)
  end function pmin_1d
  
  ! =============================================
  ! pmax_1d:
  ! 1D version of the parallel MAX
  ! =============================================
  function pmax_1d(variable,par) result(res)
    implicit none

    real(kind=real_kind),intent(in)  :: variable(:)
    type (parallel_t),intent(in)     :: par
    real(kind=real_kind)             :: res
    
    real(kind=real_kind)             :: local_sum
    integer                          :: ierr
    local_sum=MAXVAL(variable)
    call MPI_Allreduce(local_sum,res,1,MPIreal_t, &
                       MPI_MAX,par%comm,ierr)

  end function pmax_1d

  ! =============================================
  ! psum_1d:
  ! 1D version of the parallel MAX
  ! =============================================
  function psum_1d(variable,par) result(res)
    implicit none

    real(kind=real_kind),intent(in)  :: variable(:)
    type (parallel_t),intent(in)     :: par
    real(kind=real_kind)             :: res
     
    real(kind=real_kind)             :: local_sum
    integer                          :: ierr

    local_sum=SUM(variable)

    call MPI_Allreduce(local_sum,res,1,MPIreal_t, &
                       MPI_SUM,par%comm,ierr)

  end function psum_1d
  

end module parallel_mod
