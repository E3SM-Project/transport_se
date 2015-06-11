;
;  plot lat-height tracer cross-sections using native-grid data
;
;  usage: ncl dcmip1-2_lat_height.ncl
;_______________________________________________________________________

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 

begin

  ; load native-grid data

  f1 = addfile("./movies/dcmip1-21.nc","r")

  ne      = f1@ne                           ; get num elem per cube face
  np      = f1@np                           ; get num pts per elem edge

  dsizes  = getfiledimsizes(f1)             ; get dimension sizes
  dnames  = getvardims(f1)                  ; get dimension names

  ncol    = dsizes(ind(dnames.eq."ncol" ))  ; get num columns
  nlev    = dsizes(ind(dnames.eq."lev"  ))  ; get num midpoint levels
  nelem   = dsizes(ind(dnames.eq."nelem"))  ; get num elemens

  lat     = f1->lat ; lat(ncol)             ; get latitudes of columns
  lon     = f1->lon ; lon(ncol)             ; get longitudes of columns
  lev     = f1->lev ; lev(nlev)             ; get eta levels within column
 
  time    = f1->time; time(UNLIMITED)       ; get output time indices
  nt      = dimsizes(time)                  ; count num output times

  ; get slice at lon = 180 dg

  midlon= 180
  inds  = ind(lon .eq. midlon)              ; find indices where lon=180
  lat1d = lat(inds)                         ; get lats along slice
  pind  = dim_pqsort(lat1d,2)               ; sort lats1d. get permutation indices

  ; read q2(time,lev,ncol)

  q2_t0 = f1->Q2(0,:,inds)                  ; get q2 data along slice
  q2_t1 = f1->Q2(1,:,inds)
  q2_t2 = f1->Q2(2,:,inds)

  q2_t0 = q2_t0(:,pind)                     ; permute data to match lats
  q2_t1 = q2_t1(:,pind)
  q2_t2 = q2_t2(:,pind)

  ; get height from geopotential

  g     = 9.80616d0; gravity (m s^2)
  z     = f1->geo(0,:,inds)/g               ; get z data along slice
  z     = z(:,pind)                         ; permute z data to match lats
  z1d   = z(:,0)                            ; assume all columns have same z levels

  ; delete units and long names

  q2_t0@units     = ""
  q2_t1@units     = ""
  q2_t2@units     = ""
  q2_t0@long_name = ""
  q2_t1@long_name = ""
  q2_t2@long_name = ""

  ; setup the plot

  res1                      = True            ; Plot resource
  res1@gsnDraw              = False          	; panel plot
  res1@gsnFrame             = False           ; don't draw yet
  res1@lbLabelBarOn         = False           ; label bar
  res1@vpWidthF             = 2               ; set aspect ratio 2:1
  res1@vpHeightF            = 1
  res1@sfXArray             = lat1d          	; uses lat as plot x-axis
  res1@sfYArray             = z1d/1000.0        ; vertical axis = z in km

  res1@cnInfoLabelOn        = False
  res1@cnFillOn             = True            ; fill contours
  res1@cnLinesOn            = True            ; draw contour lines
  res1@cnLineLabelsOn       = False           ; don't draw labels
  res1@cnMonoLineColor      = True            ; single color lines
  res1@cnLineColor          = "black"         ; black lines
  res1@cnLineThicknessF     = 0.01             ; thin lines
  res1@cnLineDrawOrder      = "PreDraw"
  res1@cnLevelSelectionMode = "ManualLevels"
  res1@cnMaxLevelValF       = 1.0             ; max contour level
  res1@cnMinLevelValF       = 0.0             ; min contour level
  res1@cnLevelSpacingF      = 0.1             ; contour level spacing

  res1@tiYAxisString        = "height (km)"   ; y-axis title
  res1@tiXAxisString        = "latitude"      ; x-axis title
  res1@txFont               = 22              ; helvetica-bold
  res1@tiMainFont           = 22              ; helvetica-bold
  res1@tiMainFontHeightF    = 0.06            ; set font sizes
  res1@tiXAxisFontHeightF   = 0.05
  res1@tiYAxisFontHeightF   = 0.05
  res1@tmXBLabelFontHeightF = 0.04
  res1@gsnYAxisIrregular2Linear = True

  ; panel properties

  pres                      = True            ; Panel resource
  pres@gsnMaximize          = True            ; maximize panel size
  pres@gsnFrame             = False           ; don't advance yet
  pres@gsnPaperOrientation  = "portrait"      ; paper orientation
  pres@gsnPanelLabelBar     = True            ; add common label bar
  pres@lbLabelFontHeightF   = 0.012           ; set labelbar font height
  pres@txFontHeightF        = 0.016           ; set text font height
  pres@txString             = "DCMIP 1-2   q2    NE="+ne+"    lon="+midlon+" dg"

  ; Create workstation for file output

  name  = "image_dcmip1-2_lat_height_NE"+ne
  print("writing "+name)

  wks_type                  = "pdf"
  wks_type@wkPaperHeightF   = 8
  wks_type@wkPaperWidthF    = 4.5
  wks = gsn_open_wks(wks_type,name)
  wks@wkAntiAlias           = True

  ; define custom colormap

  gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")  ; choose default colormap
  colors       = gsn_retrieve_colormap(wks)           ; retrieve color map for editing
  colors(25 ,:)= (/ 1, 1, 1/)                         ; white for tracer-free region
  colors(2  ,:)= (/ 0.95, 0.95, 0.95 /)               ; light-gray for undershoots
  colors(255,:)= (/ 0.0, 0.0, 0.0 /)                  ; black for overshoots
  gsn_define_colormap(wks,colors)
	  
  ; plot countours and label bar

  plot    = new (3, graphic)
  
  res1@tiMainString = "t=0 hr "
  plot(0) = gsn_csm_contour(wks,q2_t0(:,:),res1)
  
  res1@tiMainString = "t=12 hr "
  plot(1) = gsn_csm_contour(wks,q2_t1(:,:),res1)
  
  res1@tiMainString = "t=24 hr "
  plot(2) = gsn_csm_contour(wks,q2_t2(:,:),res1)

	gsn_panel(wks,plot,(/3,1/),pres)

	frame(wks)

end
