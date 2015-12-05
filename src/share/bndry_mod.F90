#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module bndry_mod
  implicit none
  private
  public :: bndry_exchangeV

  interface bndry_exchangeV
     module procedure bndry_exchangeV_nonth
     module procedure long_bndry_exchangeV_nonth
     module procedure bndry_exchangeV_thsave 
  end interface

contains 

  !********************************************************************************
  ! Single-threaded version of bndry_exchangeV: for calls by horizontal master thread only
  !********************************************************************************
  subroutine bndry_exchangeV_nonth(par,buffer)
    use kinds,         only : log_kind
    use edge_mod,      only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod,    only : omp_get_thread_num
    use parallel_mod,  only : parallel_t, status, srequest, rrequest, &
                                mpireal_t, mpiinteger_t, mpi_success

    type (parallel_t)                :: par
    type (EdgeBuffer_t)              :: buffer
    type (Schedule_t),pointer        :: pSchedule
    type (Cycle_t),pointer           :: pCycle
    integer                          :: dest,length,tag
    integer                          :: i,icycle,ierr
    integer                          :: iptr,source,nlyr
    integer                          :: nSendCycles,nRecvCycles
    integer                          :: errorcode,errorlen
    character*(80) errorstring
    logical(kind=log_kind),parameter :: Debug=.FALSE.

    if(omp_get_thread_num() > 0) then
       print *,'bndry_exchangeV: Warning you are calling a non-thread safe'
       print *,'		 routine inside a threaded region....     '
       print *,'                 Results are not predictable!!            '
    endif

    nlyr = buffer%nlyr

#ifdef MPI_PERSISTENT
    pSchedule => Schedule(1)
    nSendCycles = SIZE(buffer%Srequest)
    nRecvCycles = SIZE(buffer%Rrequest)

    call MPI_startall(nRecvCycles,buffer%Rrequest,ierr)
    call MPI_startall(nSendCycles,buffer%Srequest,ierr)

    call MPI_Waitall(nSendCycles,buffer%Srequest,status,ierr)
    call MPI_Waitall(nRecvCycles,buffer%Rrequest,status,ierr)
#else

    ! Setup the pointer to proper Schedule
#ifdef _PREDICT
    pSchedule => Schedule(iam)
#else
    pSchedule => Schedule(1)
#endif

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    !==================================================
    !  Fire off the sends
    !==================================================
#if (defined COLUMN_OPENMP)
!$omp parallel do private(icycle,pCycle,ierr,dest,length,tag,iptr,errorcode)
#endif
    do icycle=1,nSendCycles
       pCycle         => pSchedule%SendCycle(icycle)
       dest           = pCycle%dest - 1
       length         = nlyr * pCycle%lengthP
       tag            = pCycle%tag
       iptr           = pCycle%ptrP
       !DBG if(Debug) print *,'bndry_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call MPI_Isend(buffer%buf(1,iptr),length,MPIreal_t,dest,tag,par%comm,Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Isend: ',errorstring
       endif
    end do

    !==================================================
    !  Post the Receives 
    !==================================================
#if (defined COLUMN_OPENMP)
!$omp parallel do private(icycle,pCycle,ierr,source,length,tag,iptr,errorcode)
#endif
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source         = pCycle%source - 1
       length         = nlyr * pCycle%lengthP
       tag            = pCycle%tag
       iptr           = pCycle%ptrP
       !DBG if(Debug) print *,'bndry_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call MPI_Irecv(buffer%receive(1,iptr),length,MPIreal_t, &
            source,tag,par%comm,Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
       endif
    end do

    !==================================================
    !  Wait for all the receives to complete
    !==================================================
    call MPI_Waitall(nSendCycles,Srequest,status,ierr)
    call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)

#endif

#if (defined COLUMN_OPENMP)
!$omp parallel do private(icycle,pCycle,i)
#endif
    do icycle=1,nRecvCycles
       pCycle => pSchedule%RecvCycle(icycle)
       do i=0,pCycle%lengthP-1
          buffer%buf(1:nlyr,pCycle%ptrP+i) = buffer%receive(1:nlyr,pCycle%ptrP+i)
       enddo
    end do

  end subroutine bndry_exchangeV_nonth


  !********************************************************************************
  ! Pure-MPI version of bndry_exchangeV
  !********************************************************************************
  subroutine long_bndry_exchangeV_nonth(par,buffer)
    use kinds,         only : log_kind
    use edge_mod,      only : LongEdgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod,    only : omp_in_parallel
    use parallel_mod,  only : parallel_t, status, srequest, rrequest, &
                                mpireal_t, mpiinteger_t, mpi_success

    type (parallel_t)                :: par
    type (LongEdgeBuffer_t)          :: buffer
    type (Schedule_t),pointer        :: pSchedule
    type (Cycle_t),pointer           :: pCycle
    integer                          :: dest,length,tag
    integer                          :: i,icycle,ierr
    integer                          :: iptr,source,nlyr
    integer                          :: nSendCycles,nRecvCycles
    integer                          :: errorcode,errorlen
    character*(80) errorstring
    logical(kind=log_kind),parameter :: Debug=.FALSE.

    if(omp_in_parallel()) then
       print *,'bndry_exchangeV: Warning you are calling a non-thread safe'
       print *,'		 routine inside a threaded region....     '
       print *,'                 Results are not predictable!!            '
    endif

    nlyr = buffer%nlyr
    ! Setup the pointer to proper Schedule
#ifdef _PREDICT
    pSchedule => Schedule(iam)
#else
    pSchedule => Schedule(1)
#endif
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    !==================================================
    !  Fire off the sends
    !==================================================
    do icycle=1,nSendCycles
       pCycle         => pSchedule%SendCycle(icycle)
       dest           = pCycle%dest - 1
       length         = nlyr * pCycle%lengthP
       tag            = pCycle%tag
       iptr           = pCycle%ptrP
       !DBG if(Debug) print *,'bndry_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call MPI_Isend(buffer%buf(1,iptr),length,MPIinteger_t,dest,tag,par%comm,Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Isend: ',errorstring
       endif
    end do

    !==================================================
    !  Post the Receives 
    !==================================================
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source         = pCycle%source - 1
       length         = nlyr * pCycle%lengthP
       tag            = pCycle%tag
       iptr           = pCycle%ptrP
       !DBG if(Debug) print *,'bndry_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call MPI_Irecv(buffer%receive(1,iptr),length,MPIinteger_t, &
            source,tag,par%comm,Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
       endif
    end do

    !==================================================
    !  Wait for all the receives to complete
    !==================================================
    call MPI_Waitall(nSendCycles,Srequest,status,ierr)
    call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)

    do icycle=1,nRecvCycles
       pCycle => pSchedule%RecvCycle(icycle)
       do i=0,pCycle%lengthP-1
          buffer%buf(1:nlyr,pCycle%ptrP+i) = buffer%receive(1:nlyr,pCycle%ptrP+i)
       enddo
    end do

  end subroutine long_bndry_exchangeV_nonth


  !********************************************************************************
  !  Thread-safe version of hybrid MPI+OMP bndry_exchangeV
  !********************************************************************************
  subroutine bndry_exchangeV_thsave(hybrid,buffer)
    use hybrid_mod, only : hybrid_t
    use edge_mod,   only : Edgebuffer_t
    use perf_mod,   only : t_startf, t_stopf, t_adj_detailf
    implicit none

    type (hybrid_t)      :: hybrid
    type (EdgeBuffer_t)  :: buffer

    call t_adj_detailf(+2)
    call t_startf('bndry_exchange')
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
    !$OMP MASTER
    ! Only the root thread calls the bndry_exchangeV_nonth
#endif
    call bndry_exchangeV_nonth(hybrid%par,buffer)
#if (defined HORIZ_OPENMP)
    !$OMP END MASTER
    !$OMP BARRIER
#endif
    call t_stopf('bndry_exchange')
    call t_adj_detailf(-2)

  end subroutine bndry_exchangeV_thsave


end module bndry_mod
