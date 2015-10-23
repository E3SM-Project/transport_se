#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module edge_mod
  use kinds, only : int_kind, log_kind, real_kind
  use dimensions_mod, only : max_neigh_edges, nelemd
  use perf_mod, only: t_startf, t_stopf, t_adj_detailf ! _EXTERNAL
  use thread_mod, only: omp_get_num_threads, omp_get_thread_num
  use coordinate_systems_mod, only : cartesian3D_t
#ifdef MPI_PERSISTENT
  use parallel_mod, only : abortmp, haltmp, MPIreal_t, iam,parallel_t
  use schedtype_mod, only : cycle_t, schedule_t, schedule
#else
  use parallel_mod,   only : abortmp, haltmp, parallel_t
#endif



  implicit none
  private
  save

  type, public :: rotation_t
     integer  :: nbr                                        ! nbr direction: north south east west
     integer  :: reverse                                    ! 0 = do not reverse order
     ! 1 = reverse order
     real (kind=real_kind), dimension(:,:,:), pointer :: R => null()  ! rotation matrix
  end type rotation_t

  type, public :: EdgeDescriptor_t
     integer(kind=int_kind)  :: use_rotation
     integer(kind=int_kind)  :: padding
     integer(kind=int_kind), pointer  :: putmapP(:) => null()
     integer(kind=int_kind), pointer  :: getmapP(:) => null()
     integer(kind=int_kind), pointer  :: globalID(:) => null()
     integer(kind=int_kind), pointer  :: loc2buf(:) => null()
     type (cartesian3D_t)  , pointer  :: neigh_corners(:,:) => null()
     integer                          :: actual_neigh_edges  
     logical(kind=log_kind), pointer  :: reverse(:) => null()
     type (rotation_t), dimension(:), pointer :: rot => null() ! Identifies list of edges
     !  that must be rotated, and how
  end type EdgeDescriptor_t


  type, public :: EdgeBuffer_t
     real (kind=real_kind), dimension(:,:), pointer :: buf => null()
     real (kind=real_kind), dimension(:,:), pointer :: receive => null()
     integer :: nlyr ! Number of layers
     integer :: nbuf ! size of the horizontal dimension of the buffers.
#ifdef MPI_PERSISTENT
     integer, public, pointer :: Rrequest(:)
     integer, public, pointer :: Srequest(:)
#endif
  end type EdgeBuffer_t

  type, public :: LongEdgeBuffer_t
     integer :: nlyr
     integer :: nbuf
     integer (kind=int_kind), dimension(:,:), pointer :: buf => null()
     integer (kind=int_kind), dimension(:,:), pointer :: receive => null()
  end type LongEdgeBuffer_t

  public :: initEdgeBuffer, initLongEdgeBuffer
  public :: FreeEdgeBuffer, FreeLongEdgeBuffer
  public :: edgeVpack,  edgeDGVpack, LongEdgeVpack


  public :: edgeVunpack,edgeDGVunpack, edgeVunpackVert
  public :: edgeVunpackMIN, LongEdgeVunpackMIN
  public :: edgeVunpackMAX

  public :: edgerotate
  public :: buffermap

  real(kind=real_kind), parameter, public :: edgeDefaultVal = 1.11e+100_real_kind

! NOTE ON ELEMENT ORIENTATION
!
! Element orientation:  index V(i,j)
!
!           (1,np) NWEST      (np,np) NEAST 
!
!           (1,1) SWEST       (np,1) SEAST
!
!
! for the edge neighbors:  
!    we set the "reverse" flag if two elements who share an edge use a 
!    reverse orientation.  The data is reversed during the *pack* stage
! For corner neighbors:  
!    for edge buffers, there is no orientation because two corner neighbors
!    only share a single point.
!
! This routine only works for meshes with at most 1 corner element.  It's
! not called and the corner orientation flag is not set for unstructured meshes

! Wrap pointer so we can make an array of them.
  type :: wrap_ptr
     real (kind=real_kind), dimension(:,:), pointer :: ptr => null()
  end type wrap_ptr

  type(wrap_ptr) :: edgebuff_ptrs(0:1)

contains

  ! =========================================
  ! initEdgeBuffer:
  !
  ! create an Real based communication buffer
  ! =========================================
  subroutine initEdgeBuffer(par,edge,nlyr, buf_ptr,receive_ptr)
    use dimensions_mod, only : np, nelemd, max_corner_elem
    implicit none
    type (parallel_t), intent(in) :: par
    integer,intent(in)                :: nlyr
    type (EdgeBuffer_t),intent(out), target :: edge
    real(kind=real_kind), optional, pointer :: buf_ptr(:), receive_ptr(:)

    ! Notes about the buf_ptr/receive_ptr options:
    !
    ! You can pass in 1D pointers to this function. If they are not
    ! associated, they will be allocated and used as buffer space. If they
    ! are associated, their targets will be used as buffer space.
    !
    ! The pointers must not be thread-private.
    !
    ! If an EdgeBuffer_t object is initialized from pre-existing storage
    ! (i.e. buf_ptr is provided and not null), it must *not* be freed,
    ! and must not be used if the underlying storage has been deallocated.
    !
    ! All these restrictions also applied to the old newbuf and newreceive
    ! options.

    ! Workaround for NAG bug.
    ! NAG 5.3.1 dies if you use pointer bounds remapping to set
    ! a pointer that is also a component. So remap to temporary,
    ! then use that to set component pointer.
    real(kind=real_kind), pointer :: tmp_ptr(:,:)

    ! Local variables
    integer :: nbuf,ith
#ifdef MPI_PERSISTENT
    integer :: nSendCycles, nRecvCycles
    integer :: icycle, ierr
    type (Cycle_t), pointer :: pCycle
    type (Schedule_t), pointer :: pSchedule
    integer :: dest, source, length, tag, iptr
#endif

    nbuf=4*(np+max_corner_elem)*nelemd
    edge%nlyr=nlyr
    edge%nbuf=nbuf
    if (nlyr==0) return  ! tracer code might call initedgebuffer() with zero tracers

!$OMP BARRIER

!   only master thread should allocate the buffer
#if (defined HORIZ_OPENMP)
!$OMP MASTER
#endif
    if (present(buf_ptr)) then
       ! If buffer is passed in but not allocated, allocate it.
       if (.not. associated(buf_ptr)) allocate(buf_ptr(nlyr*nbuf))
       ! Verify dimensions
       if (size(buf_ptr) < nlyr*nbuf) then
          print *,'size(buf_ptr),nlyr,nbuf=',size(buf_ptr),nlyr,nbuf
          call abortmp('Error: user provided edge buffer is too small')
       end if
#ifdef HAVE_F2003_PTR_BND_REMAP
       tmp_ptr(1:nlyr,1:nbuf) => buf_ptr
       edge%buf => tmp_ptr
#else
       ! call F77 routine which will reshape array.
       call remap_2D_ptr_buf(edge,nlyr,nbuf,buf_ptr)
#endif
    else
       allocate(edge%buf    (nlyr,nbuf))
    end if

    if (present(receive_ptr)) then
       ! If buffer is passed in but not allocated, allocate it.
       if (.not. associated(receive_ptr)) allocate(receive_ptr(nlyr*nbuf))
       ! Verify dimensions
       if (size(receive_ptr) < nlyr*nbuf) then
          print *,'size(receive_ptr),nlyr,nbuf=',size(receive_ptr),nlyr,nbuf
          call abortmp('Error: user provided edge buffer is too small')
       end if
#ifdef HAVE_F2003_PTR_BND_REMAP
       tmp_ptr(1:nlyr,1:nbuf) => receive_ptr
       edge%receive => tmp_ptr
#else
       ! call F77 routine which will reshape array.
       call remap_2D_ptr_receive(edge,nlyr,nbuf,receive_ptr)
#endif
    else
       allocate(edge%receive(nlyr,nbuf))
    endif
    edge%buf    (:,:)=0.0D0
    edge%receive(:,:)=0.0D0

#ifdef MPI_PERSISTENT

    pSchedule => Schedule(1)
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles
!    print *,'iam: ',iam, ' nSendCycles: ',nSendCycles, ' nRecvCycles: ',
!    nRecvCycles
    allocate(edge%Rrequest(nRecvCycles))
    allocate(edge%Srequest(nSendCycles))
    do icycle=1,nSendCycles
       pCycle => pSchedule%SendCycle(icycle)
       dest   = pCycle%dest -1
       length = nlyr * pCycle%lengthP
       tag    = pCycle%tag
       iptr   = pCycle%ptrP
!       print *,'IAM: ',iam, ' length: ',length,' dest: ',dest,' tag: ',tag
       call MPI_Send_init(edge%buf(1,iptr),length,MPIreal_t,dest,tag,par%comm, edge%Srequest(icycle),ierr)
    enddo
    do icycle=1,nRecvCycles
       pCycle => pSchedule%RecvCycle(icycle)
       source   = pCycle%source -1
       length = nlyr * pCycle%lengthP
       tag    = pCycle%tag
       iptr   = pCycle%ptrP
!       print *,'IAM: ',iam, 'length: ',length,' dest: ',source,' tag: ',tag
       call MPI_Recv_init(edge%receive(1,iptr),length,MPIreal_t,source,tag,par%comm, edge%Rrequest(icycle),ierr)
    enddo
#endif

#if (defined HORIZ_OPENMP)
!$OMP END MASTER
#endif
!   make sure all threads wait until buffer is allocated
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif

    ! sanity check on edge.  edge must NOT be thread-prive, but must be shared by all threads
    ! the calling program needs to instantiate 'edge' outside the threaded region.
    ! if 'edge' is thread private, it creates flaky openMP problems that are difficut to debug
    ! so lets try and detect it here:
    if (omp_get_num_threads()>1) then
       ith=omp_get_thread_num()
       if (ith <= 1 ) then
          edgebuff_ptrs(ith)%ptr => edge%buf
       endif
#if (defined HORIZ_OPENMP)
       !$OMP BARRIER
       !$OMP MASTER
#endif
       if (.not. associated(edgebuff_ptrs(0)%ptr, edgebuff_ptrs(1)%ptr)) then
          call haltmp('ERROR: edge struct appears to be thread-private.')
       endif
#if (defined HORIZ_OPENMP)
       !$OMP END MASTER
#endif
    endif

  end subroutine initEdgeBuffer
  ! =========================================
  ! initLongEdgeBuffer:
  !
  ! create an Integer based communication buffer
  ! =========================================
  subroutine initLongEdgeBuffer(edge,nlyr)
    use dimensions_mod, only : np, nelemd, max_corner_elem
    implicit none
    integer,intent(in)                :: nlyr
    type (LongEdgeBuffer_t),intent(out) :: edge

    ! Local variables

    integer :: nbuf

    ! sanity check for threading
    if (omp_get_num_threads()>1) then
       call haltmp('ERROR: initLongEdgeBuffer must be called before threaded reagion')
    endif

    nbuf=4*(np+max_corner_elem)*nelemd
    edge%nlyr=nlyr
    edge%nbuf=nbuf
    allocate(edge%buf(nlyr,nbuf))
    edge%buf(:,:)=0

    allocate(edge%receive(nlyr,nbuf))
    edge%receive(:,:)=0

  end subroutine initLongEdgeBuffer
  ! =========================================
  ! edgeDGVpack:
  !
  ! Pack edges of v into buf for DG stencil
  ! =========================================
  subroutine edgeDGVpack(edge,v,vlyr,kptr,desc)
    use dimensions_mod, only : np
    type (EdgeBuffer_t)                      :: edge
    integer,              intent(in)   :: vlyr
    real (kind=real_kind),intent(in)   :: v(np,np,vlyr)
    integer,              intent(in)   :: kptr
    type (EdgeDescriptor_t)            :: desc

    ! =========================================
    ! This code is just a wrapper call the 
    !   normal edgeVpack
    ! =========================================
    call edgeVpack(edge,v,vlyr,kptr,desc)

  end subroutine edgeDGVpack

  ! ===========================================
  !  FreeEdgeBuffer:
  !
  !  Freed an edge communication buffer
  ! =========================================
  subroutine FreeEdgeBuffer(edge) 
    implicit none
    type (EdgeBuffer_t),intent(inout) :: edge

#if (defined HORIZ_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif
    edge%nbuf=0
    edge%nlyr=0
    deallocate(edge%buf)
    deallocate(edge%receive)
#if (defined HORIZ_OPENMP)
!$OMP END MASTER
#endif

  end subroutine FreeEdgeBuffer

  ! ===========================================
  !  FreeLongEdgeBuffer:
  !
  !  Freed an edge communication buffer
  ! =========================================
  subroutine FreeLongEdgeBuffer(edge) 
    implicit none
    type (LongEdgeBuffer_t),intent(inout) :: edge

    edge%nbuf=0
    edge%nlyr=0
    deallocate(edge%buf)
    deallocate(edge%receive)

  end subroutine FreeLongEdgeBuffer

  ! =========================================
  !
  !> @brief Pack edges of v into an edge buffer for boundary exchange.
  !
  !> This subroutine packs for one or more vertical layers into an edge 
  !! buffer. If the buffer associated with edge is not large enough to 
  !! hold all vertical layers you intent to pack, the method will 
  !! halt the program with a call to parallel_mod::haltmp().
  !! @param[in] edge Edge Buffer into which the data will be packed.
  !! This buffer must be previously allocated with initEdgeBuffer().
  !! @param[in] v The data to be packed.
  !! @param[in] vlyr Number of vertical level coming into the subroutine
  !! for packing for input v.
  !! @param[in] kptr Vertical pointer to the place in the edge buffer where 
  !! data will be located.
  ! =========================================
  subroutine edgeVpack(edge,v,vlyr,kptr,desc)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod,    only : north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t)                :: edge
    integer,                intent(in) :: vlyr
    real (kind=real_kind),  intent(in) :: v(np,np,vlyr)
    integer,                intent(in) :: kptr
    type (EdgeDescriptor_t),intent(in) :: desc

    ! Local variables
    integer :: i,k,is,ie,in,iw,sw,se,nw,ne,kk

    if (edge%nlyr < (kptr+vlyr) ) then
       call haltmp('edgeVpack: Buffer overflow: size of the vertical dimension must be increased!')
    endif

    call t_adj_detailf(+2)
    call t_startf('edge_pack')

    is = desc%putmapP(south)
    ie = desc%putmapP(east)
    in = desc%putmapP(north)
    iw = desc%putmapP(west)

    ! Copy the 4 corners
    sw = desc%putmapP(swest)
    se = desc%putmapP(seast)
    ne = desc%putmapP(neast)
    nw = desc%putmapP(nwest)

#if (defined COLUMN_OPENMP)
    !$omp parallel private(k,i,kk)
#endif
    if(MODULO(np,4) == 0) then 
#if (defined COLUMN_OPENMP)
       !$omp do
#endif
       do k=1,vlyr
          kk=kptr+k
          do i=1,np,4
             edge%buf(kk,is+i)   = v(i  ,1  ,k)  
             edge%buf(kk,is+i+1) = v(i+1,1  ,k)  
             edge%buf(kk,is+i+2) = v(i+2,1  ,k)  
             edge%buf(kk,is+i+3) = v(i+3,1  ,k)  
             edge%buf(kk,ie+i)   = v(np ,i  ,k)  
             edge%buf(kk,ie+i+1) = v(np ,i+1,k)
             edge%buf(kk,ie+i+2) = v(np ,i+2,k)
             edge%buf(kk,ie+i+3) = v(np ,i+3,k)
             edge%buf(kk,in+i)   = v(i  ,np ,k)
             edge%buf(kk,in+i+1) = v(i+1,np ,k)
             edge%buf(kk,in+i+2) = v(i+2,np ,k)
             edge%buf(kk,in+i+3) = v(i+3,np ,k)
             edge%buf(kk,iw+i)   = v(1  ,i  ,k)  
             edge%buf(kk,iw+i+1) = v(1  ,i+1,k)
             edge%buf(kk,iw+i+2) = v(1  ,i+2,k)
             edge%buf(kk,iw+i+3) = v(1  ,i+3,k)
          enddo
       end do
    else
#if (defined COLUMN_OPENMP)
       !$omp do
#endif
       do k=1,vlyr
          kk=kptr+k
          do i=1,np
             edge%buf(kk,is+i)   = v(i  ,1 ,k) ! South
             edge%buf(kk,ie+i)   = v(np ,i ,k) ! East
             edge%buf(kk,in+i)   = v(i  ,np,k) ! North
             edge%buf(kk,iw+i)   = v(1  ,i ,k) ! West
          enddo
       end do
    endif

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing

    if(desc%reverse(south)) then
#if (defined COLUMN_OPENMP)
       !$omp do
#endif
       do k=1,vlyr
          do i=1,np
             edge%buf(kptr+k,is+np-i+1)=v(i,1,k)
          enddo
       enddo
    endif

    if(desc%reverse(east)) then
#if (defined COLUMN_OPENMP)
       !$omp do
#endif
       do k=1,vlyr
          do i=1,np
             edge%buf(kptr+k,ie+np-i+1)=v(np,i,k)
          enddo
       enddo
    endif

    if(desc%reverse(north)) then
#if (defined COLUMN_OPENMP)
       !$omp do
#endif
       do k=1,vlyr
          do i=1,np
             edge%buf(kptr+k,in+np-i+1)=v(i,np,k)
          enddo
       enddo
    endif

    if(desc%reverse(west)) then
#if (defined COLUMN_OPENMP)
       !$omp do
#endif
       do k=1,vlyr
          do i=1,np
             edge%buf(kptr+k,iw+np-i+1)=v(1,i,k)
          enddo
       enddo
    endif

#if (defined COLUMN_OPENMP)
       !$omp do
#endif
    do k=1,vlyr
       kk=kptr+k
       if (sw /= -1) then
          edge%buf(kk,sw+1)=v(1 ,1 ,k) ! SWEST
       end if
       if (se /= -1) then
          edge%buf(kk,se+1)=v(np,1 ,k) ! SEAST
       end if
       if (ne /= -1) then
          edge%buf(kk,ne+1)=v(np,np,k) ! NEAST
       end if
       if (nw /= -1) then
          edge%buf(kk,nw+1)=v(1 ,np,k) ! NWEST
       end if
    end do
#if (defined COLUMN_OPENMP)
    !$omp end parallel
#endif

    call t_stopf('edge_pack')
    call t_adj_detailf(-2)
  end subroutine edgeVpack


  ! =========================================
  ! LongEdgeVpack:
  !
  ! Pack edges of v into buf...
  ! =========================================
  subroutine LongEdgeVpack(edge,v,vlyr,kptr,desc)
    use control_mod,    only : north, south, east, west, neast, nwest, seast, swest
    use dimensions_mod, only : np, max_corner_elem

    type (LongEdgeBuffer_t)            :: edge
    integer,                intent(in) :: vlyr
    integer (kind=int_kind),intent(in) :: v(np,np,vlyr)
    integer,                intent(in) :: kptr
    type (EdgeDescriptor_t),intent(in) :: desc

    ! Local variables
    logical, parameter :: UseUnroll = .TRUE.
    integer :: i,k,ir,l,is,ie,in,iw

    is = desc%putmapP(south)
    ie = desc%putmapP(east)
    in = desc%putmapP(north)
    iw = desc%putmapP(west)

    if(MODULO(np,2) == 0 .and. UseUnroll) then 
       do k=1,vlyr
          do i=1,np,2
             edge%buf(kptr+k,is+i)   = v(i  ,1  ,k)
             edge%buf(kptr+k,is+i+1) = v(i+1,1  ,k)
             edge%buf(kptr+k,ie+i)   = v(np ,i  ,k)
             edge%buf(kptr+k,ie+i+1) = v(np ,i+1,k)
             edge%buf(kptr+k,in+i)   = v(i  ,np ,k)
             edge%buf(kptr+k,in+i+1) = v(i+1,np ,k)
             edge%buf(kptr+k,iw+i)   = v(1  ,i  ,k)
             edge%buf(kptr+k,iw+i+1) = v(1  ,i+1,k)
          enddo
       end do
    else
       do k=1,vlyr
          do i=1,np
             edge%buf(kptr+k,is+i)   = v(i  ,1 ,k)
             edge%buf(kptr+k,ie+i)   = v(np ,i ,k)
             edge%buf(kptr+k,in+i)   = v(i  ,np,k)
             edge%buf(kptr+k,iw+i)   = v(1  ,i ,k)
          enddo
       end do
    endif

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    if(desc%reverse(south)) then
       is = desc%putmapP(south)
       do k=1,vlyr
          do i=1,np
             ir = np-i+1
             edge%buf(kptr+k,is+ir)=v(i,1,k)
          enddo
       enddo
    endif

    if(desc%reverse(east)) then
       ie = desc%putmapP(east)
       do k=1,vlyr
          do i=1,np
             ir = np-i+1
             edge%buf(kptr+k,ie+ir)=v(np,i,k)
          enddo
       enddo
    endif

    if(desc%reverse(north)) then
       in = desc%putmapP(north)
       do k=1,vlyr
          do i=1,np
             ir = np-i+1
             edge%buf(kptr+k,in+ir)=v(i,np,k)
          enddo
       enddo
    endif

    if(desc%reverse(west)) then
       iw = desc%putmapP(west)
       do k=1,vlyr
          do i=1,np
             ir = np-i+1
             edge%buf(kptr+k,iw+ir)=v(1,i,k)
          enddo
       enddo
    endif

    ! Corners
    ! SWEST
    do l=swest,swest+max_corner_elem-1
        if (desc%putmapP(l) /= -1) then
            do k=1,vlyr
                edge%buf(kptr+k,desc%putmapP(l)+1)=v(1  ,1 ,k)
            end do
        end if
    end do

    ! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
        if (desc%putmapP(l) /= -1) then
            do k=1,vlyr
                edge%buf(kptr+k,desc%putmapP(l)+1)=v(np ,1 ,k)
            end do
        end if
    end do

    ! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        if (desc%putmapP(l) /= -1) then
            do k=1,vlyr
                edge%buf(kptr+k,desc%putmapP(l)+1)=v(np ,np,k)
            end do
        end if
    end do

    ! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        if (desc%putmapP(l) /= -1) then
            do k=1,vlyr
                edge%buf(kptr+k,desc%putmapP(l)+1)=v(1  ,np,k)
            end do
        end if
    end do
  end subroutine LongEdgeVpack


  ! ========================================
  ! edgeVunpack:
  !
  ! Unpack edges from edge buffer into v...
  ! ========================================
  subroutine edgeVunpack(edge,v,vlyr,kptr,desc)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod,    only : north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t),   intent(in)    :: edge
    integer,               intent(in)    :: vlyr
    real (kind=real_kind), intent(inout) :: v(np,np,vlyr)
    integer,               intent(in)    :: kptr
    type (EdgeDescriptor_t)              :: desc

    ! Local
    integer :: i,k,is,ie,in,iw,sw,se,ne,nw,kk

    call t_adj_detailf(+2)
    call t_startf('edge_unpack')

    is=desc%getmapP(south)
    ie=desc%getmapP(east)
    in=desc%getmapP(north)
    iw=desc%getmapP(west)

    ! Copy the 4 corners
    sw = desc%getmapP(swest)
    se = desc%getmapP(seast)
    ne = desc%getmapP(neast)
    nw = desc%getmapP(nwest)

#if (defined COLUMN_OPENMP)
    !$omp parallel private(k,i,kk)
#endif
    if(MODULO(np,4) == 0) then 
#if (defined COLUMN_OPENMP)
       !$omp do
#endif
       do k=1,vlyr
          kk=kptr+k
          do i=1,np,4
             v(i  ,1  ,k) = v(i  ,1  ,k)+edge%buf(kk,is+i  ) ! South
             v(i+1,1  ,k) = v(i+1,1  ,k)+edge%buf(kk,is+i+1)
             v(i+2,1  ,k) = v(i+2,1  ,k)+edge%buf(kk,is+i+2)
             v(i+3,1  ,k) = v(i+3,1  ,k)+edge%buf(kk,is+i+3)
             v(np ,i  ,k) = v(np ,i  ,k)+edge%buf(kk,ie+i  ) ! East
             v(np ,i+1,k) = v(np ,i+1,k)+edge%buf(kk,ie+i+1)
             v(np ,i+2,k) = v(np ,i+2,k)+edge%buf(kk,ie+i+2)
             v(np ,i+3,k) = v(np ,i+3,k)+edge%buf(kk,ie+i+3)
             v(i  ,np ,k) = v(i  ,np ,k)+edge%buf(kk,in+i  ) ! North
             v(i+1,np ,k) = v(i+1,np ,k)+edge%buf(kk,in+i+1)
             v(i+2,np ,k) = v(i+2,np ,k)+edge%buf(kk,in+i+2)
             v(i+3,np ,k) = v(i+3,np ,k)+edge%buf(kk,in+i+3)
             v(1  ,i  ,k) = v(1  ,i  ,k)+edge%buf(kk,iw+i  ) ! West
             v(1  ,i+1,k) = v(1  ,i+1,k)+edge%buf(kk,iw+i+1)
             v(1  ,i+2,k) = v(1  ,i+2,k)+edge%buf(kk,iw+i+2)
             v(1  ,i+3,k) = v(1  ,i+3,k)+edge%buf(kk,iw+i+3)
          enddo
       end do
    else
#if (defined COLUMN_OPENMP)
       !$omp do
#endif
       do k=1,vlyr
         kk=kptr+k
         do i=1,np
           v(i ,1 ,k) = v(i ,1 ,k)+edge%buf(kk,is+i) ! South
           v(np,i ,k) = v(np,i ,k)+edge%buf(kk,ie+i) ! East
           v(i ,np,k) = v(i ,np,k)+edge%buf(kk,in+i) ! North
           v(1 ,i ,k) = v(1 ,i ,k)+edge%buf(kk,iw+i) ! West
         end do
       end do
    endif

#if (defined COLUMN_OPENMP)
    !$omp do
#endif
    do k=1,vlyr
       kk=kptr+k
       if (sw /= -1) then
          v(1 ,1 ,k)=v(1 ,1 ,k)+edge%buf(kk,sw+1) ! SWEST
       end if
       if (se /= -1) then
          v(np,1 ,k)=v(np,1 ,k)+edge%buf(kk,se+1) ! SEAST
       end if
       if (ne /= -1) then
          v(np,np,k)=v(np,np,k)+edge%buf(kk,ne+1) ! NEAST
       end if
       if (nw /= -1) then
          v(1 ,np,k)=v(1 ,np,k)+edge%buf(kk,nw+1) ! NWEST
       end if
    end do
#if (defined COLUMN_OPENMP)
    !$omp end parallel
#endif

    call t_stopf('edge_unpack')
    call t_adj_detailf(-2)
  end subroutine edgeVunpack



  subroutine edgeVunpackVert(edge,v,desc)
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use dimensions_mod, only : np, max_corner_elem, ne
    use coordinate_systems_mod, only: cartesian3D_t

    type (EdgeBuffer_t),   intent(inout)  :: edge
    type (cartesian3D_t), intent(out) :: v(:,:,:)
    type (EdgeDescriptor_t)            :: desc

    ! Local
    logical, parameter :: UseUnroll = .TRUE.
    integer :: i,k,l
    integer :: is,ie,in,iw

    if (max_corner_elem.ne.1 .and. ne==0) then
        ! MNL: this is used to construct the dual grid on the cube,
        !      currently only supported for the uniform grid. If
        !      this is desired on a refined grid, a little bit of
        !      work will be required.
        call haltmp("edgeVunpackVert should not be called with unstructured meshes")
    end if

    is=desc%getmapP(south)
    ie=desc%getmapP(east)
    in=desc%getmapP(north)
    iw=desc%getmapP(west)

    ! N+S
    do i=1,np/2
       ! North
       v(3,i ,np)%x = edge%buf(1,in+i) 
       v(3,i ,np)%y = edge%buf(2,in+i) 
       v(3,i ,np)%z = edge%buf(3,in+i) 
       ! South
       v(2,i ,1)%x  = edge%buf(1,is+i) 
       v(2,i ,1)%y  = edge%buf(2,is+i) 
       v(2,i ,1)%z  = edge%buf(3,is+i) 
    enddo

    do i=np/2+1,np
       ! North
       v(4,i ,np)%x = edge%buf(1,in+i) 
       v(4,i ,np)%y = edge%buf(2,in+i) 
       v(4,i ,np)%z = edge%buf(3,in+i) 
       ! South
       v(1,i ,1)%x  = edge%buf(1,is+i) 
       v(1,i ,1)%y  = edge%buf(2,is+i) 
       v(1,i ,1)%z  = edge%buf(3,is+i)        
    enddo

    do i=1,np/2
       ! East
       v(3,np,i)%x = edge%buf(1,ie+i)
       v(3,np,i)%y = edge%buf(2,ie+i)
       v(3,np,i)%z = edge%buf(3,ie+i)       
       ! West
       v(4,1,i)%x  = edge%buf(1,iw+i)
       v(4,1,i)%y  = edge%buf(2,iw+i)
       v(4,1,i)%z  = edge%buf(3,iw+i)
    end do

    do i=np/2+1,np
       ! East
       v(2,np,i)%x = edge%buf(1,ie+i)
       v(2,np,i)%y = edge%buf(2,ie+i)
       v(2,np,i)%z = edge%buf(3,ie+i)       
       ! West
       v(1,1,i)%x  = edge%buf(1,iw+i)
       v(1,1,i)%y  = edge%buf(2,iw+i)
       v(1,1,i)%z  = edge%buf(3,iw+i)
    end do

! SWEST
    do l=swest,swest+max_corner_elem-1
       ! find the one active corner, then exist
        if(desc%getmapP(l) /= -1) then 
            v(1,1,1)%x=edge%buf(1,desc%getmapP(l)+1)
            v(1,1,1)%y=edge%buf(2,desc%getmapP(l)+1)
            v(1,1,1)%z=edge%buf(3,desc%getmapP(l)+1)
            exit 
        else
            v(1,1,1)%x=0_real_kind
            v(1,1,1)%y=0_real_kind
            v(1,1,1)%z=0_real_kind
        endif
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
       ! find the one active corner, then exist
        if(desc%getmapP(l) /= -1) then 
            v(2,np,1)%x=edge%buf(1,desc%getmapP(l)+1)
            v(2,np,1)%y=edge%buf(2,desc%getmapP(l)+1)
            v(2,np,1)%z=edge%buf(3,desc%getmapP(l)+1)
            exit
        else
            v(2,np,1)%x=0_real_kind
            v(2,np,1)%y=0_real_kind
            v(2,np,1)%z=0_real_kind
        endif
    end do

! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
       ! find the one active corner, then exist
        if(desc%getmapP(l) /= -1) then 
            v(3,np,np)%x=edge%buf(1,desc%getmapP(l)+1)
            v(3,np,np)%y=edge%buf(2,desc%getmapP(l)+1)
            v(3,np,np)%z=edge%buf(3,desc%getmapP(l)+1)
            exit
        else
            v(3,np,np)%x=0_real_kind
            v(3,np,np)%y=0_real_kind
            v(3,np,np)%z=0_real_kind
        endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
       ! find the one active corner, then exist
        if(desc%getmapP(l) /= -1) then 
            v(4,1,np)%x=edge%buf(1,desc%getmapP(l)+1)
            v(4,1,np)%y=edge%buf(2,desc%getmapP(l)+1)
            v(4,1,np)%z=edge%buf(3,desc%getmapP(l)+1)
            exit
        else
            v(4,1,np)%x=0_real_kind
            v(4,1,np)%y=0_real_kind
            v(4,1,np)%z=0_real_kind
        endif
    end do

    ! Fill the missing vertex info

    do i=2,np/2
       ! North
       v(4,i ,np)%x = v(3,i-1 ,np)%x 
       v(4,i ,np)%y = v(3,i-1 ,np)%y
       v(4,i ,np)%z = v(3,i-1 ,np)%z
       ! South
       v(1,i ,1)%x  = v(2,i-1 ,1)%x 
       v(1,i ,1)%y  = v(2,i-1 ,1)%y 
       v(1,i ,1)%z  = v(2,i-1 ,1)%z 
    enddo

    do i=np/2+1,np-1
       ! North
       v(3,i ,np)%x = v(4,i+1 ,np)%x 
       v(3,i ,np)%y = v(4,i+1 ,np)%y
       v(3,i ,np)%z = v(4,i+1 ,np)%z
       ! South
       v(2,i ,1)%x  = v(1,i+1 ,1)%x 
       v(2,i ,1)%y  = v(1,i+1 ,1)%y
       v(2,i ,1)%z  = v(1,i+1 ,1)%z
    enddo

    do i=2,np/2
       ! East
       v(2,np,i)%x = v(3,np,i-1)%x
       v(2,np,i)%y = v(3,np,i-1)%y
       v(2,np,i)%z = v(3,np,i-1)%z
       ! West
       v(1,1,i)%x  = v(4,1,i-1)%x
       v(1,1,i)%y  = v(4,1,i-1)%y
       v(1,1,i)%z  = v(4,1,i-1)%z
    end do

    do i=np/2+1,np-1
       ! East
       v(3,np,i)%x = v(2,np,i+1)%x 
       v(3,np,i)%y = v(2,np,i+1)%y
       v(3,np,i)%z = v(2,np,i+1)%z
       ! West
       v(4,1,i)%x  = v(1,1,i+1)%x 
       v(4,1,i)%y  = v(1,1,i+1)%y
       v(4,1,i)%z  = v(1,1,i+1)%z
    end do
  end subroutine edgeVunpackVert


  ! ========================================
  ! edgeDGVunpack:
  !
  ! Unpack edges from edge buffer into v...
  ! ========================================
  subroutine edgeDGVunpack(edge,v,vlyr,kptr,desc)
    use dimensions_mod, only : np
    use control_mod,    only : north, south, east, west

    type (EdgeBuffer_t),   intent(in)  :: edge
    integer,               intent(in)  :: vlyr
    real (kind=real_kind), intent(out) :: v(0:np+1,0:np+1,vlyr)
    integer,               intent(in)  :: kptr
    type (EdgeDescriptor_t)            :: desc

    ! Local
    integer :: i,k
    integer :: is,ie,in,iw

    is=desc%getmapP(south)
    ie=desc%getmapP(east)
    in=desc%getmapP(north)
    iw=desc%getmapP(west)
    do k=1,vlyr
       do i=1,np
          v(i   ,0   ,k)=edge%buf(kptr+k,is+i)
          v(np+1,i   ,k)=edge%buf(kptr+k,ie+i)
          v(i   ,np+1,k)=edge%buf(kptr+k,in+i)
          v(0   ,i   ,k)=edge%buf(kptr+k,iw+i)
       end do
    end do
  end subroutine edgeDGVunpack


  ! ========================================
  ! edgeVunpackMIN/MAX:
  !
  ! Finds the Min/Max edges from edge buffer into v...
  ! ========================================
  subroutine edgeVunpackMAX(edge,v,vlyr,kptr,desc)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod,    only : north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t),    intent(in)    :: edge
    integer,                intent(in)    :: vlyr
    real (kind=real_kind),  intent(inout) :: v(np,np,vlyr)
    integer,                intent(in)    :: kptr
    type (EdgeDescriptor_t),intent(in)    :: desc

    ! Local
    integer :: i,k,kk
    integer :: is,ie,in,iw,sw,se,ne,nw

    call t_startf('edgeVunpackMAX')

    is=desc%getmapP(south)
    ie=desc%getmapP(east)
    in=desc%getmapP(north)
    iw=desc%getmapP(west)

    sw = desc%getmapP(swest)
    se = desc%getmapP(seast)
    ne = desc%getmapP(neast)
    nw = desc%getmapP(nwest)

#if (defined COLUMN_OPENMP)
    !$omp parallel private(k,i,kk)
    !$omp do
#endif
    do k=1,vlyr
       kk=kptr+k
       do i=1,np
          v(i ,1 ,k) = MAX(v(i ,1 ,k),edge%buf(kk,is+i))
          v(np,i ,k) = MAX(v(np,i ,k),edge%buf(kk,ie+i))
          v(i ,np,k) = MAX(v(i ,np,k),edge%buf(kk,in+i))
          v(1 ,i ,k) = MAX(v(1 ,i ,k),edge%buf(kk,iw+i))
       end do
    end do

#if (defined COLUMN_OPENMP)
    !$omp do
#endif
    do k=1,vlyr
       kk=kptr+k
       if (sw /= -1) then
          v(1 ,1 ,k)=MAX(v(1 ,1 ,k),edge%buf(kk,sw+1))
       end if
       if (se /= -1) then
          v(np,1 ,k)=MAX(v(np,1 ,k),edge%buf(kk,se+1))
       end if
       if (ne /= -1) then
          v(np,np,k)=MAX(v(np,np,k),edge%buf(kk,ne+1))
       end if
       if (nw /= -1) then
          v(1 ,np,k)=MAX(v(1 ,np,k),edge%buf(kk,nw+1))
       end if
    end do
#if (defined COLUMN_OPENMP)
    !$omp end parallel
#endif
    
    call t_stopf('edgeVunpackMAX')
  end subroutine edgeVunpackMAX


  subroutine edgeVunpackMIN(edge,v,vlyr,kptr,desc)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod,    only : north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t),    intent(in)    :: edge
    integer,                intent(in)    :: vlyr
    real (kind=real_kind),  intent(inout) :: v(np,np,vlyr)
    integer,                intent(in)    :: kptr
    type (EdgeDescriptor_t),intent(in)    :: desc

    ! Local
    integer :: i,k,kk
    integer :: is,ie,in,iw,sw,se,ne,nw

    call t_startf('edgeVunpackMIN')

    is=desc%getmapP(south)
    ie=desc%getmapP(east)
    in=desc%getmapP(north)
    iw=desc%getmapP(west)
#if (defined COLUMN_OPENMP)
    !$omp parallel private(k,i,kk)
    !$omp do
#endif
    do k=1,vlyr
       kk=kptr+k
       do i=1,np
          v(i ,1 ,k) = MIN(v(i ,1 ,k),edge%buf(kk,is+i))
          v(np,i ,k) = MIN(v(np,i ,k),edge%buf(kk,ie+i))
          v(i ,np,k) = MIN(v(i ,np,k),edge%buf(kk,in+i))
          v(1 ,i ,k) = MIN(v(1 ,i ,k),edge%buf(kk,iw+i))
       end do
    end do

    sw = desc%getmapP(swest)
    se = desc%getmapP(seast)
    ne = desc%getmapP(neast)
    nw = desc%getmapP(nwest)
#if (defined COLUMN_OPENMP)
    !$omp do
#endif
    do k=1,vlyr
       kk=kptr+k
       if (sw /= -1) then
          v(1 ,1 ,k)=MIN(v(1 ,1 ,k),edge%buf(kk,sw+1))
       end if
       if (se /= -1) then
          v(np,1 ,k)=MIN(v(np,1 ,k),edge%buf(kk,se+1))
       end if
       if (ne /= -1) then
          v(np,np,k)=MIN(v(np,np,k),edge%buf(kk,ne+1))
       end if
       if (nw /= -1) then
          v(1 ,np,k)=MIN(v(1 ,np,k),edge%buf(kk,nw+1))
       end if
    end do
#if (defined COLUMN_OPENMP)
    !$omp end parallel
#endif
    
    call t_stopf('edgeVunpackMIN')
  end subroutine edgeVunpackMIN


  ! ========================================
  ! LongEdgeVunpackMIN:
  !
  ! Finds the Min edges from edge buffer into v...
  ! ========================================
  subroutine LongEdgeVunpackMIN(edge,v,vlyr,kptr,desc)
    use control_mod,    only : north, south, east, west, neast, nwest, seast, swest
    use dimensions_mod, only : np, max_corner_elem

    type (LongEdgeBuffer_t), intent(in)    :: edge
    integer,                 intent(in)    :: vlyr
    integer (kind=int_kind), intent(inout) :: v(np,np,vlyr)
    integer,                 intent(in)    :: kptr
    type (EdgeDescriptor_t), intent(in)    :: desc

    ! Local
    integer :: i,k,kk
    integer :: is,ie,in,iw,sw,se,ne,nw

    call t_startf('LongEdgeVunpackMIN')

    is=desc%getmapP(south)
    ie=desc%getmapP(east)
    in=desc%getmapP(north)
    iw=desc%getmapP(west)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,kk)
#endif
    do k=1,vlyr
       kk=kptr+k
       do i=1,np
          v(i ,1 ,k) = MIN(v(i ,1 ,k),edge%buf(kk,is+i))
          v(np,i ,k) = MIN(v(np,i ,k),edge%buf(kk,ie+i))
          v(i ,np,k) = MIN(v(i ,np,k),edge%buf(kk,in+i))
          v(1 ,i ,k) = MIN(v(1 ,i ,k),edge%buf(kk,iw+i))
       end do
    end do

    sw = desc%getmapP(swest)
    se = desc%getmapP(seast)
    ne = desc%getmapP(neast)
    nw = desc%getmapP(nwest)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,kk)
#endif
    do k=1,vlyr
       kk=kptr+k
       if (sw /= -1) then
          v(1 ,1 ,k)=MIN(v(1 ,1 ,k),edge%buf(kk,sw+1))
       end if
       if (se /= -1) then
          v(np,1 ,k)=MIN(v(np,1 ,k),edge%buf(kk,se+1))
       end if
       if (ne /= -1) then
          v(np,np,k)=MIN(v(np,np,k),edge%buf(kk,ne+1))
       end if
       if (nw /= -1) then
          v(1 ,np,k)=MIN(v(1 ,np,k),edge%buf(kk,nw+1))
       end if
    end do

    call t_stopf('LongEdgeVunpackMIN')
  end subroutine LongEdgeVunpackMIN


  ! =============================
  ! edgerotate:
  !
  ! Rotate edges in buffer...
  ! =============================

  subroutine edgerotate(edge,vlyr,kptr,desc)
    use dimensions_mod, only : np
    type (EdgeBuffer_t)           :: edge         ! edge struct
    integer, intent(in)           :: vlyr         ! number of 2d vector fields to rotate
    integer, intent(in)           :: kptr         ! layer pointer into edge buffer
    type (EdgeDescriptor_t)       :: desc

    ! Local variables

    integer :: i,k,k1,k2
    integer :: irot,ia,nbr

    real(kind=real_kind), dimension(:,:,:), pointer :: R
    real(kind=real_kind)  :: tmp1,tmp2

#ifdef _USEASSOCIATED
    if (associated(rot)) then
#else
       if (desc%use_rotation == 1) then
#endif

          do irot=1,SIZE(desc%rot)

             nbr  =  desc%rot(irot)%nbr
             R    => desc%rot(irot)%R

             ia=desc%putmapP(nbr)

             ! ========================================
             ! If nbr direction is (1-4) => is an edge
             ! ========================================

             if (nbr <= 4) then

                ! ========================================================
                !  Is an edge. Rotate it in place
                ! ========================================================

                do i=1,np
                   do k=1,vlyr,2
                      k1 = kptr + k
                      k2 = k1 + 1
                      tmp1=R(1,1,i)*edge%buf(k1,ia+i) + R(1,2,i)*edge%buf(k2,ia+i)
                      tmp2=R(2,1,i)*edge%buf(k1,ia+i) + R(2,2,i)*edge%buf(k2,ia+i)
                      edge%buf(k1,ia+i)=tmp1
                      edge%buf(k2,ia+i)=tmp2
                   end do
                end do

             else

                ! ===================================================
                ! Is an element corner point, but not a cube corner
                ! point, just rotate it in place.
                ! ===================================================

                if (ia /= -1) then
                   do k=1,vlyr,2
                      k1 = kptr + k
                      k2 = k1+1
                      tmp1=R(1,1,1)*edge%buf(k1,ia+1) + R(1,2,1)*edge%buf(k2,ia+1)
                      tmp2=R(2,1,1)*edge%buf(k1,ia+1) + R(2,2,1)*edge%buf(k2,ia+1)
                      edge%buf(k1,ia+1)=tmp1
                      edge%buf(k2,ia+1)=tmp2
                   end do
                end if

             end if

          end do

       endif
  end subroutine edgerotate

     ! =============================================
     ! buffermap:
     !
     ! buffermap translates element number, inum and
     ! element edge/corner, facet, into an edge buffer 
     ! memory location, loc.
     ! =============================================

     function buffermap(inum,facet) result(loc)
       use dimensions_mod, only : np
       integer, intent(in) :: inum   
       integer, intent(in) :: facet
       integer :: loc

       if (facet>4) then
          if (inum == -1) then
             loc = inum
          else
             loc=(inum-1)*(4*np+4)+4*np+(facet-5)
          end if
       else
          loc=(inum-1)*(4*np+4)+np*(facet-1)
       end if

     end function buffermap

End module edge_mod



#ifndef HAVE_F2003_PTR_BND_REMAP
!
! subroutine to allow sharing edge buffers
! this has to be outside a module to allow us to (F77 style) access the same chunk 
! of memory with a different shape
!
! some compilers dont allow the 'target' attribute to be used in a F77 subroutine
! such as cray.  if that is the case, try compiling with -DHAVE_F2003_PTR_BND_REMAP
!
subroutine remap_2D_ptr_buf(edge,nlyr,nbuf,src_array)
  use edge_mod, only : EdgeBuffer_t ! _EXTERNAL
  use kinds,    only : real_kind
  ! input
  type (EdgeBuffer_t)           :: edge
  integer                       :: nlyr,nbuf
  real(kind=real_kind) , target :: src_array(nlyr,nbuf)

  edge%buf  => src_array

end subroutine remap_2D_ptr_buf


subroutine remap_2D_ptr_receive(edge,nlyr,nbuf,src_array)
  use edge_mod, only : EdgeBuffer_t ! _EXTERNAL
  use kinds,    only : real_kind
  ! input
  type (EdgeBuffer_t)           :: edge
  integer                       :: nlyr,nbuf
  real(kind=real_kind) , target :: src_array(nlyr,nbuf)

  edge%receive => src_array

end subroutine remap_2D_ptr_receive
#endif
