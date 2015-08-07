#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module domain_mod
  use perf_mod, only: t_startf, t_stopf ! _EXTERNAL
  implicit none
  private

  type, public :: domain1D_t
     integer :: start             ! start index
     integer :: end               ! end index
  end type domain1D_t

  public :: decompose

contains           

  function decompose(start,end,nthr,i) result(domain)

    ! distribute elements evenly amongst threads

    integer, intent(in) :: start  ! start index
    integer, intent(in) :: end    ! end   index
    integer, intent(in) :: nthr   ! number of threads
    integer, intent(in) :: i      ! thread index

    type (domain1D_t) :: domain
    integer :: n(nthr),l,q,r      ! n_elements,length,quotient,remainder
    call t_startf('decompose')

    l     = end-start+1           ! get length
    q     = l / nthr              ! get quotient
    r     = mod(l,nthr)           ! get remainder
    n(:)  = q                     ! set q elements per thread
    n(1:r)= q+1                   ! add extra element to first r threads

    domain%start= start + sum( n(1:i) )
    domain%end  = domain%start + n(i+1) -1

    call t_stopf('decompose')
  end function decompose

end module domain_mod





