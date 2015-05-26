;
; measure error norms for test 1-2
;_______________________________________________________________________

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   

begin

  ; parameters

  pi  = acos(-1.)
  g   = 9.80616d0       ; gravity (m/s^2)
  rad = pi/180.         ; radians per degree
  R   = 6.37122e6       ; radius of the earth

  ; open netcdf file

	f1  = addfile("./movies/dcmip1-21.nc","r")
	;f1  = addfile("../../../run/movies/dcmip1-21.nc","r")

  ; read grid information

	lat    = f1->lat * rad
	lon    = f1->lon * rad
	lev    = f1->lev
  time   = f1->time
  
	nlat   = getfilevardimsizes(f1, "lat" )
	nlon   = getfilevardimsizes(f1, "lon" )
	nlev   = getfilevardimsizes(f1, "lev" )
  nt     = dimsizes(time)

  ; read initial and final states (time,lev,lat,long)

  q_i = f1->Q2(0,:,:,:)                 ; initial state
	q_f = f1->Q2(nt-1,:,:,:)              ; final state
  dq = q_f - q_i

  ; get vertical extent from geopotential
 
  h   = f1->geo(0,:,0,0)/g              ; get level height
  dh  = new(nlev,typeof(h))
  base= 0.0*h(0)                        ; assume bottom is at z=0

  do i = 1,nlev
    dh(nlev-i) = 2.0*(h(nlev-i)-base)   ; get level width
    base       = base + dh(nlev-i)      ; move interface up
  end do

  ; get sample volumes

  dlat  = lat(1)-lat(0)
  dlon  = lon(1)-lon(0)
  dx = R*dlon
  dy = conform(q_i, R*cos(lat)*dlat, 1)
  dz = conform(q_i, dh, 0)
  dV = dx * dy * dz

  ; get error norms

  L1  =      sum(abs(dq)^1*dV) /      sum(abs(q_i - avg(q_i))^1*dV )
  L2  = sqrt(sum(abs(dq)^2*dV))/ sqrt(sum(abs(q_i - avg(q_i))^2*dV ))
  Linf=      max(abs(dq)  *dV) /      max(abs(q_i - avg(q_i))*dV   )
  q_max = max(q_f)
  q_min = min(q_f)
  
  print("DCMIP 1-2: L1="+L1+" L2= "+L2+" Linf= "+Linf+" q_max= "+q_max+" q_min= "+q_min)

end
