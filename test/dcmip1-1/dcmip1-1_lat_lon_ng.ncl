;
;  plot lattitude-longitude tracer cross-sections
;
;  usage: ncl dcmip1-1_lat_lon.ncl NE=<number of elements per edge>
;_______________________________________________________________________

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   

begin

  interval = 3.0; number of days per output

  ; load native-grid data

  f = addfile("./movies/dcmip1-11.nc","r")
  NE      = f@ne                            ; get num elem per cube face
  np      = f@np                            ; get num pts per elem edge
  dsizes  = getfiledimsizes(f)              ; get dimension sizes
  dnames  = getvardims(f)                   ; get dimension names
  ncol    = dsizes(ind(dnames.eq."ncol" ))  ; get num columns
  nlev    = dsizes(ind(dnames.eq."lev"  ))  ; get num midpoint levels
  nelem   = dsizes(ind(dnames.eq."nelem"))  ; get num elemens
  lat     = f->lat ; lat(ncol)              ; get latitudes of columns
  lon     = f->lon ; lon(ncol)              ; get longitudes of columns
  lev     = f->lev ; lev(nlev)              ; get eta levels within column
  time    = f->time; time(UNLIMITED)        ; get output time indices
  nt      = dimsizes(time)                  ; count num output times

  ; locate vertical level closest to z=5100 m

  g       = 9.80616d0; Gravity (m s^2)
  heights = f->geo(0,:,0)/g
  level   = (nlev-1) - closest_val(5100,heights(::-1))
  z       = round( heights(level), 3)

  ; load tracer data

  q1_t0 = f->Q (0,level,:)
  q1_t1 = f->Q (1,level,:)
  q1_t2 = f->Q (2,level,:)

  q3_t0 = f->Q3(0,level,:)
  q3_t1 = f->Q3(1,level,:)
  q3_t2 = f->Q3(2,level,:)

  ; delete units and long names

  q1_t0@units = " "
  q1_t1@units = " "
  q1_t2@units = " "
  
  q3_t0@units = " "
  q3_t1@units = " "
  q3_t2@units = " "

  q1_t0@long_name = " "
  q1_t1@long_name = " "
  q1_t2@long_name = " "

  q3_t0@long_name = " "
  q3_t1@long_name = " "
  q3_t2@long_name = " "

  ; setup the plot

  res1                      = True
  res1@gsnDraw              = False        		; panel plot
  res1@gsnFrame             = False           ; don't draw yet
  res1@gsnCenterString      = ""
  res1@vpWidthF             = 2
  res1@vpHeightF            = 1
  res1@sfXArray             = lon           	; uses lon as plot x-axis
  res1@sfYArray             = lat         		; uses lat for y axis
  res1@tiYAxisString        = "latitiude"     ; y-axis title
  res1@tiXAxisString        = "longitude"     ; x-axis title
  res1@lbLabelBarOn         = False           ; label bar

  res1@cnLinesOn            = False ;True
  res1@cnFillOn             = True
  res1@cnInfoLabelOn        = False           ; don't give each
  res1@cnLineLabelsOn       = False           ; panel its own
  res1@cnMonoLineColor      = True            ; single color lines
  res1@cnLineColor          = "black"         ; black lines
  res1@cnLineDrawOrder      = "PreDraw"
  res1@cnLineThicknessF     = 0.1
  res1@cnFillOpacityF       = 1.0
  res1@cnLevelSelectionMode = "ManualLevels"
  res1@cnMaxLevelValF       = 1.00            ; max contour color label
  res1@cnMinLevelValF       = 0.00            ; min contour color label
  res1@cnLevelSpacingF      = 0.10            ; contour color spacing

  ; text properties

  res1@tmXBMajorLengthF     = 0.03
  res1@txFont               = 22              ; helvetica-bold
  res1@tiMainFont           = 22              ; helvetica-bold
  res1@tiMainFontHeightF    = 0.06            ; set font sizes
  res1@tiXAxisFontHeightF   = 0.05
  res1@tiYAxisFontHeightF   = 0.05
  res1@tmXBLabelFontHeightF = 0.04
  ;res1@tiXAxisOffsetYF      = 0.010           ; shift x label
  ;res1@tiYAxisOffsetXF      = 0.010           ; shift y label

  ; panel properties

  pres                      = True            ; Panel resource
  pres@gsnMaximize          = True            ; maximize panel size
  pres@gsnFrame             = False           ; don't advance yet
  pres@gsnPaperOrientation  = "portrait"      ; paper orientation
  pres@gsnPanelLabelBar     = True            ; add common label bar
  pres@lbLabelFontHeightF   = 0.012            ; set labelbar font height
  pres@txFontHeightF        = 0.018            ; set text font height
  pres@txString             = "DCMIP 1-1    NE="+ NE +"    z="+z+" m"

  ; set output type and name

  filename ="image_dcmip1-1_lat_lon_NE"+NE
  print("writing "+filename)

  wks_type                  = "pdf"
  wks_type@wkPaperHeightF   = 8
  wks_type@wkPaperWidthF    = 8
	wks = gsn_open_wks(wks_type,filename)

  ; define custom colormap

  gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")      ; choose default colormap
  colors       = gsn_retrieve_colormap(wks)               ; retrieve color map for editing
  colors(25 ,:)= (/ 1, 1, 1/)                             ; white for tracer-free region
  colors(2  ,:)= (/ 0.95, 0.95, 0.95 /)                   ; light-gray for undershoots
  colors(255,:)= (/ 0.0, 0.0, 0.0 /)                      ; black for overshoots
  gsn_define_colormap(wks,colors)

  ; plot countours and label bar

  plot    = new (6, graphic)
  res1@tiMainString = "q1 t=0 days min/max="+sprintf("%7.4f", min(q1_t0))+" "+sprintf("%7.4f", max(q1_t0))
  plot(0) = gsn_csm_contour(wks,q1_t0(:),res1)
  
  res1@tiMainString = "q1 t=6 days min/max="+sprintf("%7.4f", min(q1_t1))+" "+sprintf("%7.4f", max(q1_t1))
  plot(2) = gsn_csm_contour(wks,q1_t1(:),res1)
  
  res1@tiMainString = "q1 t=12 days min/max="+sprintf("%7.4f", min(q1_t2))+" "+sprintf("%7.4f", max(q1_t2))
  plot(4) = gsn_csm_contour(wks,q1_t2(:),res1)

  
  res1@tiMainString = "q3 t=0 days min/max="+sprintf("%7.4f", min(q3_t0))+" "+sprintf("%7.4f", max(q3_t0))
  plot(1) = gsn_csm_contour(wks,q3_t0(:),res1)
    
  res1@tiMainString = "q3 t=6 days min/max="+sprintf("%7.4f", min(q3_t1))+" "+sprintf("%7.4f", max(q3_t1))
  plot(3) = gsn_csm_contour(wks,q3_t1(:),res1)

  res1@tiMainString = "q3 t=12 days min/max="+sprintf("%7.4f", min(q3_t2))+" "+sprintf("%7.4f", max(q3_t2))
  plot(5) = gsn_csm_contour(wks,q3_t2(:),res1)

  gsn_panel(wks,plot,(/3,2/),pres)

  frame(wks)

end




