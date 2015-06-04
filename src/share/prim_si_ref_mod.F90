#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_si_ref_mod
  use kinds, only: r8 => real_kind, iulog
  use dimensions_mod, only: plev => nlev, plevp => nlevp
  implicit none
  private

  public  :: prim_set_mass

contains

  subroutine prim_set_mass(elem, tl,hybrid,hvcoord,nets,nete)
  use kinds, only : real_kind
  use control_mod, only : initial_total_mass
  use physical_constants, only : g
  use element_mod, only : element_t
  use time_mod, only : timelevel_t 
  use hybvcoord_mod, only : hvcoord_t 
  use hybrid_mod, only : hybrid_t
  use dimensions_mod, only : np
  use global_norms_mod, only : global_integral 

  type (element_t), intent(inout) :: elem(:)
  type (TimeLevel_t), target, intent(in) :: tl
  type (hybrid_t),intent(in)     :: hybrid
  type (hvcoord_t), intent(in)   :: hvcoord
  integer,intent(in)             :: nets,nete
  
  ! local 
  real (kind=real_kind)  :: tmp(np,np,nets:nete)
  real (kind=real_kind)  :: scale,mass0
  integer :: n0,nm1,np1,ie

  if (initial_total_mass == 0) return;
  
  n0=tl%n0
  nm1=tl%nm1
  np1=tl%np1
  
  scale=1/g                                  ! assume code is using Pa
  if (hvcoord%ps0 <  2000 ) scale=100*scale  ! code is using mb
  ! after scaling, Energy is in J/m**2,  Mass kg/m**2
  
  do ie=nets,nete
     tmp(:,:,ie)=elem(ie)%state%ps_v(:,:,n0)
  enddo
  mass0 = global_integral(elem, tmp(:,:,nets:nete),hybrid,np,nets,nete)
  mass0 = mass0*scale;  
  
  do ie=nets,nete
     elem(ie)%state%ps_v(:,:,n0)=elem(ie)%state%ps_v(:,:,n0)*(initial_total_mass/mass0)
     elem(ie)%state%ps_v(:,:,np1)=elem(ie)%state%ps_v(:,:,n0)
     elem(ie)%state%ps_v(:,:,nm1)=elem(ie)%state%ps_v(:,:,n0)
  enddo
  if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
     write (*,'(a,e24.15)') "Initializing Total Mass (kg/m^2) = ",initial_total_mass
  endif
  end subroutine prim_set_mass


end module prim_si_ref_mod
