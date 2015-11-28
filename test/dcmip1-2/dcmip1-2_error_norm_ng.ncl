;
; measure error norms for test 1-2 from native grid data
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

  f1 = addfile("./movies/dcmip1-21.nc","r")

  ne      = f1@ne                           ; get num elem per cube face
  np      = f1@np                           ; get num pts per elem edge

  dsizes  = getfiledimsizes(f1)             ; get dimension sizes
  dnames  = getvardims(f1)                  ; get dimension names

  ncol    = dsizes(ind(dnames.eq."ncol" ))  ; get num columns
  nlev    = dsizes(ind(dnames.eq."lev"  ))  ; get num midpoint levels

  lat     = f1->lat *rad ; lat(ncol)        ; get latitudes of columns
  lon     = f1->lon *rad ; lon(ncol)        ; get longitudes of columns
  lev     = f1->lev      ; lev(nlev)        ; get eta levels within column

  time    = f1->time; time(UNLIMITED)       ; get output time indices
  nt      = dimsizes(time)                  ; count num output times

  ; read initial and final states (time,lev,lat,long)

  ; Q2(time,lev,ncol)
  q_i = f1->Q2(0,:,:)                       ; get initial state
	q_f = f1->Q2(nt-1,:,:)                    ; get final state
  dq = q_f - q_i

  ; get vertical extent from geopotential
 
  h   = f1->geo(0,:,0)/g                    ; get level height
  dh  = new(nlev,typeof(h))
  base= 0.0*h(0)                            ; assume bottom is at z=0

  do i = 1,nlev
    dh(nlev-i) = 2.0*(h(nlev-i)-base)       ; get level thickness
    base       = base + dh(nlev-i)          ; move interface up
  end do

  ; approximate element volumes (assuming constant angular extent)

  dlat = 0.5*pi/(ne*(np-1))                 ; estimate angular extent
  dlon = dlat                               ; assume square elements

  dx = conform(q_i, R*cos(lat)*dlon, 1)
  dy = R*dlat
  dz = conform(q_i, dh, 0)
  dV = dx * dy * dz

  ; get error norms

  L1  =      sum(abs(dq)^1*dV) /      sum(abs(q_i - avg(q_i))^1*dV )
  L2  = sqrt(sum(abs(dq)^2*dV))/ sqrt(sum(abs(q_i - avg(q_i))^2*dV ))
  Linf=      max(abs(dq)  *dV) /      max(abs(q_i - avg(q_i))*dV   )
  q_max = max(q_f)
  q_min = min(q_f)
  
  print("DCMIP 1-2: L1="+sprintf("%8.6f",L1)+\
" L2="+sprintf("%8.6f",L2)+" Linf="+sprintf("%8.6f",Linf)+\
" q_max="+sprintf("%8.6f",q_max)+" q_min="+sprintf("%14.6e",q_min))

end
