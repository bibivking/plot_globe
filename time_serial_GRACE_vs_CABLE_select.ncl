;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

begin

case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs-standard/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para_redistrb-on/Outputs-standard/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-standard/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs-standard/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-Haverd2013/gw_on"\
             /)

case_name_brf = (/\
              "NoGW",\
              "NoGW_redstrb",\
              "GW",\
              "NoAQ",\
              "GW_Hvrd"\
               /)

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)

; ______________________________ Pre-define ________________________
year_s   = 2003
year_e   = 2010
year_sum = year_e-year_s+1
total    = year_sum*12
case_sum = dimsizes(case_name)

; ____________________________ Options _______________________________
use_month    = "May" ; "Jul"
no_greenland = True ;False ;True
time_scale   = "monthly" ; "annual"; "monthly"
grace_version= "Gap-fill" ;"Gap-fill"; "Original"
scale        = "45S-N";"30S-N";"Global" ; "50S-N" ; "SE-AU" ; "AU"

if scale .eq. "Global" then
  A = -90.
  B = 90.
  C = 0.
  D = 360.
else if scale .eq. "30S-N" then
  A = -30.
  B = 30.
  C = 0.
  D = 360.
else if scale .eq. "45S-N" then
  A = -45.
  B = 45.
  C = 0.
  D = 360.
else if scale .eq. "50S-N" then
  A = -50.
  B = 50.
  C = 0.
  D = 360.
else if scale .eq. "SE-AU" then
  A = -40.
  B = -27.
  C = 140.
  D = 155.
else if scale .eq. "AU" then
  A = -44.
  B = -10.
  C = 112.
  D = 154.
end if
end if
end if
end if
end if
end if
lat_num = tointeger((B - A)/0.5)
print(lat_num)
; _____________________________ Check ______________________________
pic = "check"
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"ViBlGrWhYeOrRe")
res = True
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                            ;???
res@tmBorderThicknessF = 3.0
res@cnLinesOn          = False                       ; ???????????
res@gsnDraw            = True  ; Don't draw plots
res@gsnFrame           = True  ; ?????

; _______________________ Read No Greenland Mask _____________________
if no_greenland .eq. True then
   filename   = "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/plots/land_mask_0.5_no_Greenland.nc"
   f          = addfile (filename,"r")
   landmask   = lonFlip(f->mask(0,:,:))
end if

; _______________________ read GRACE and fill gap ____________________________
GRACE      = new((/96,360,720/),float)
uncertainty= new((/96,360,720/),float)
GRACE      = 0.
uncertainty= 0.

if grace_version .eq. "Original" then
filename   = "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/plots/GRACE_JPL_RL06/GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01CRIv01.nc"
f          = addfile (filename,"r")

GRACE(0:4,:,:)       = dble2flt(f->lwe_thickness(7:11,:,:))*10.
GRACE(6:95,:,:)      = dble2flt(f->lwe_thickness(12:101,:,:))*10.

uncertainty(0:4,:,:) = dble2flt(f->uncertainty(7:11,:,:))*10.
uncertainty(6:95,:,:)= dble2flt(f->uncertainty(12:101,:,:))*10.

printVarSummary(GRACE)
printVarSummary(uncertainty)

May_ave = new((/360,720/),float)
Jun_ave = new((/360,720/),float)
Jul_ave = new((/360,720/),float)
May_ave = 0.
Jun_ave = 0.
Jul_ave = 0.

GRACE@_FillValue       = dble2flt(f->lwe_thickness@_FillValue)
uncertainty@_FillValue = dble2flt(f->lwe_thickness@_FillValue)

printVarSummary(GRACE)
printVarSummary(uncertainty)

plot = gsn_csm_contour(wks,GRACE(0,:,:),res)

do year = year_s+1,year_e
print("add year :" + year)
May_ave = May_ave + GRACE((year-year_s)*12+4,:,:)
Jun_ave = Jun_ave + GRACE((year-year_s)*12+5,:,:)
Jul_ave = Jul_ave + GRACE((year-year_s)*12+6,:,:)
end do

May_ave = May_ave/(year_e-year_s)
Jun_ave = Jun_ave/(year_e-year_s)
Jul_ave = Jul_ave/(year_e-year_s)

plot = gsn_csm_contour(wks,May_ave,res)
plot = gsn_csm_contour(wks,Jun_ave,res)
plot = gsn_csm_contour(wks,Jul_ave,res)

if use_month .eq. "Jul" then
   GRACE(5,:,:) = GRACE(6,:,:) * (Jun_ave/Jul_ave)
else if use_month .eq. "May" then
   GRACE(5,:,:) = GRACE(4,:,:) * (Jun_ave/May_ave)
end if
end if

uncertainty(5,:,:) = (uncertainty(4,:,:) + uncertainty(6,:,:))/2.

delete(May_ave)
delete(Jun_ave)
delete(Jul_ave)

plot = gsn_csm_contour(wks,GRACE(5,:,:),res)
plot = gsn_csm_contour(wks,uncertainty(5,:,:),res)

else if grace_version .eq. "Gap-fill" then

do year = year_s,year_e

filename   = "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/plots/GRACE_JPL_RL06/fill_gap_filter_singular_Sanaa/GRACERel06_lwe_"+ year+".nc"
f          = addfile (filename,"r")

a = (year-year_s)*12
b = (year-year_s)*12+11

GRACE(a:b,:,:)       = lonFlip(f->lwe)*10.
uncertainty(a:b,:,:) = lonFlip(f->uncertainty)*10.

end do

printVarSummary(GRACE)
printVarSummary(uncertainty)

GRACE@_FillValue       = -1.e+33
uncertainty@_FillValue = -1.e+33

replace_ieeenan (GRACE, GRACE@_FillValue, 0)
replace_ieeenan (uncertainty, uncertainty@_FillValue, 0)

printVarSummary(GRACE)
printVarSummary(uncertainty)

plot = gsn_csm_contour(wks,GRACE(0,:,:),res)
plot = gsn_csm_contour(wks,uncertainty(0,:,:),res)


end if
end if

;______________________________ INPUT __________________________________
TWS     = new((/case_sum,96,360,720/),float)
TWS_ave = new((/case_sum,360,720/),float)
TWS     = 0.0
TWS_ave = 0.0

do case_num = 0,case_sum-1

   print(case_name(case_num))

   counter = 0

   do year = year_s,year_e
      print("year = " +year)

      filename = case_name(case_num)+"/cable_out_"+year+".nc"
      f = addfile (filename,"r")

      do month = 0,11
         if case_num .le. 3 then
            TWS(case_num,counter,:,:) = f->CanopInt(month,:,:) + f->SWE(month,:,:)
         else
            TWS(case_num,counter,:,:) = f->GWMoist(month,:,:)*f->GWdz*1000. + f->CanopInt(month,:,:) + f->SWE(month,:,:)
         end if
         do layer = 0,5
            TWS(case_num,counter,:,:) = TWS(case_num,counter,:,:) + f->SoilMoist(month,layer,:,:)*Soil_thickness(layer)*1000.
         end do
         counter = counter + 1
      end do ;month
   end do ;year

end do ;case

; Calculate 2004-2009 TWS average
TWS_ave =  dim_avg_n_Wrap(TWS(:,12:83,:,:),1)

; Calculate TWS anormaly
do time = 0,95
  TWS(:,time,:,:) = TWS(:,time,:,:) - TWS_ave
end do

delete(TWS_ave)

print("point 1 ")

if no_greenland .eq. True then
do case_num = 0,case_sum-1
   do time = 0,total-1
      TWS(case_num,time,:,:) = where( landmask .eq. 1, TWS(case_num,time,:,:), TWS@_FillValue)
   end do
end do
delete(landmask)
end if

GRACE       = where( .not. ismissing(TWS(0,:,:,:)), GRACE, GRACE@_FillValue)
uncertainty = where( .not. ismissing(TWS(0,:,:,:)), uncertainty, uncertainty@_FillValue)

printVarSummary(GRACE)
printVarSummary(uncertainty)

plot = gsn_csm_contour(wks,GRACE(5,:,:),res)
plot = gsn_csm_contour(wks,uncertainty(5,:,:),res)
plot = gsn_csm_contour(wks,TWS(0,5,:,:),res)

TWS!0 = "case"
TWS!1 = "time"
TWS!2 = "lat"
TWS!3 = "lon"
TWS&lat = f->Qle&y
TWS&lon = f->Qle&x
TWS@_FillValue   = f->Qle@_FillValue

GRACE!0 = "time"
GRACE!1 = "lat"
GRACE!2 = "lon"
GRACE&lat = f->Qle&y
GRACE&lon = f->Qle&x

uncertainty!0 = "time"
uncertainty!1 = "lat"
uncertainty!2 = "lon"
uncertainty&lat = f->Qle&y
uncertainty&lon = f->Qle&x


; _________________________________  Zonal AVE  ______________________________________
print("zonal average")
GRACE_zonal_ave       = new((/lat_num,96/),float)
TWS_zonal_ave         = new((/case_sum,lat_num,96/),float)
GRACE_zonal_ave       = 0.0
TWS_zonal_ave         = 0.0

printVarSummary(GRACE(0,0,{A:B}))
do time = 0,95
   GRACE_zonal_ave(:,time) = dim_avg_n_Wrap(GRACE(time,{A:B},:),1)
   do case_num = 0,case_sum-1
      TWS_zonal_ave(case_num,:,time) = dim_avg_n_Wrap(TWS(case_num,time,{A:B},:),1)
   end do
end do

printVarSummary(f->y({A:B}))
GRACE_zonal_ave!0 = "lat"
GRACE_zonal_ave!1 = "time"
GRACE_zonal_ave&lat  = f->y({A:B})
GRACE_zonal_ave&time = ispan(1,total,1)

TWS_zonal_ave!0 = "case"
TWS_zonal_ave!1 = "lat"
TWS_zonal_ave!2 = "time"
TWS_zonal_ave&lat  = f->y({A:B})
TWS_zonal_ave&time = ispan(1,total,1)

; multi-year zonal average
;GRACE_zonal_ave_multiyear = new(lat_num,float)
TWS_zonal_ave_multiyear   = new((/case_sum+1,lat_num/),float)
;GRACE_zonal_ave_multiyear = 0.0
TWS_zonal_ave_multiyear   = 0.0

;GRACE_zonal_ave_multiyear = dim_avg_n_Wrap(GRACE_zonal_ave({A:B},:),1)
do case_num = 0,case_sum-1
   TWS_zonal_ave_multiyear(case_num,:) = dim_avg_n_Wrap(TWS_zonal_ave(case_num,{A:B},:),1)
end do
TWS_zonal_ave_multiyear(case_sum,:) = dim_avg_n_Wrap(GRACE_zonal_ave({A:B},:),1)
print("finish zonal average")

; _________________________________  Global AVE  ______________________________________
print("global average")
GRACE_global_ave       = new((/3,96/),float)
uncertainty_global_ave = new(96,float)
TWS_global_ave         = new((/case_sum,96/),float)
GRACE_global_ave       = 0.0
uncertainty_global_ave = 0.0
TWS_global_ave         = 0.0

GRACE_TWS              = new(case_sum,float)
GRACE_TWS              = 0.0

if scale .eq. "Global" then
  A = -90.
  B = 90.
  C = 0.
  D = 360.
else if scale .eq. "30S-N" then
  A = -30.
  B = 30.
  C = 0.
  D = 360.
else if scale .eq. "45S-N" then
  A = -45.
  B = 45.
  C = 0.
  D = 360.
else if scale .eq. "50S-N" then
  A = -50.
  B = 50.
  C = 0.
  D = 360.
else if scale .eq. "SE-AU" then
  A = -40.
  B = -27.
  C = 140.
  D = 155.
else if scale .eq. "AU" then
  A = -44.
  B = -10.
  C = 112.
  D = 154.
end if
end if
end if
end if
end if
end if

do time = 0,95
   GRACE_global_ave(0,time) = avg(GRACE(time,{A:B},{C:D}))
   print("GRACE_global_ave: "+GRACE_global_ave(0,time))
   uncertainty_global_ave(time) = avg(uncertainty(time,{A:B},{C:D}))
   do case_num = 0,case_sum-1
      TWS_global_ave(case_num,time) = avg(TWS(case_num,time,{A:B},{C:D}))
   end do
end do

GRACE_global_ave(1,:) = GRACE_global_ave(0,:) - uncertainty_global_ave
GRACE_global_ave(2,:) = GRACE_global_ave(0,:) + uncertainty_global_ave


do case_num = 0,case_sum-1
   GRACE_TWS(case_num) = escorc(GRACE_global_ave(0,:), TWS_global_ave(case_num,:))
   print(case_name_brf(case_num)+" "+GRACE_TWS(case_num))
end do


GRACE_global_ave_year       = new((/3,8/),float)
uncertainty_global_ave_year = new(8,float)
TWS_global_ave_year         = new((/case_sum,8/),float)
GRACE_global_ave_year       = 0.0
uncertainty_global_ave_year = 0.0
TWS_global_ave_year         = 0.0

do year = year_s,year_e
GRACE_global_ave_year(:,year-year_s) = dim_avg_n_Wrap(GRACE_global_ave(:,((year-year_s)*12):((year-year_s)*12+11)),1)
TWS_global_ave_year(:,year-year_s)   = dim_avg_n_Wrap(TWS_global_ave(:,((year-year_s)*12):((year-year_s)*12+11)),1)
uncertainty_global_ave_year(year-year_s) = dim_avg_n_Wrap(uncertainty_global_ave(((year-year_s)*12):((year-year_s)*12+11)),0)
end do

delete(GRACE)
delete(TWS)
delete(uncertainty)

;GRACE & TWS is anormaly here

;_____________________________ PLOT ___________________________________
if no_greenland .eq. True then
   pic = "time_serial_GRACE_vs_CABLE_no_Greenland_"+time_scale+"_GRACE06-"+grace_version+"_"+scale
else
   pic = "time_serial_GRACE_vs_CABLE_"+time_scale+"_GRACE06-"+grace_version+"_"+scale
end if

wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")

pres                   = True
pres@gsnMaximize        = False
res                    = True              ; plot mods desired
res@gsnDraw            = False             ; don't draw yet
res@gsnFrame           = False             ; don't advance frame yet
res@gsnMaximize        = False

res@vpHeightF 	 = 0.4               ; change aspect ratio of plot
res@vpWidthF 	   = 0.7

res@tiMainString   = "(a) "
res@tiYAxisString	 = "mm"              ; add an axis title
res@txFontHeightF	 = 0.01            ; change title font heights

res@trYMaxF = 70.
res@trYMinF = -60.

LineColors       = (/"darkgreen","green","red","pink","orange"/) ; colors chosen
LineThicknesses	 = (/2.,2.,2.,2.,2./)      ; line thicknesses
DashPatterns	   = (/0.,0.,0.,0.,0./)      ; make all lines solid

if time_scale .eq. "monthly" then
  x = ispan(1,total,1)
  res@tmXBMode  = "Explicit"
  res@tmXBValues=(/7,19,31,43,55,67,79,91/)
  res@tmXBLabels=(/"2003","2004","2005","2006","2007","2008","2009","2010"/)

  res@xyMonoLineColor    = False             ; want colored lines
  res@xyLineColor        = "Blue"
  res@xyLineThicknessF   = 2.
  top_plot = gsn_csm_xy (wks,x,GRACE_global_ave(0,:),res)       ; create line plot

; Create a plot with the area between both curves filled in blue.
;  delete(res@xyLineColor)
;  res@gsnXYFillColors = "LightBlue"
;  res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
;  bot_plot  = gsn_csm_xy (wks,x,GRACE_global_ave(1:2,:),res)  ; Create filled XY plot.

; Create a plot with the area between both curves filled in pink.
  ;delete(res@gsnXYFillColors)
  res@xyLineColors       = LineColors
  res@xyLineThicknesses  = LineThicknesses
  res@xyDashPatterns	   = DashPatterns
  mid_plot  = gsn_csm_xy (wks,x,TWS_global_ave,res)  ; Create another filled XY plot.

  overlay(mid_plot,top_plot)
;  overlay(bot_plot,mid_plot)
;  overlay(bot_plot,top_plot)

;  maximize_output(wks,pres)

else if time_scale .eq. "annual" then

  x = ispan(1,8,1)
  res@tmXBValues=(/1,2,3,4,5,6,7,8/)
  res@tmXBLabels=(/"2003","2004","2005","2006","2007","2008","2009","2010"/)

  res@xyMonoLineColor    = False             ; want colored lines
  res@xyLineColor        = "Blue"
  res@xyLineThicknessF   = 2.
  top_plot = gsn_csm_xy (wks,x,GRACE_global_ave_year(0,:),res)       ; create line plot

; Create a plot with the area between both curves filled in blue.
  delete(res@xyLineColor)
  res@gsnXYFillColors    = "LightBlue"
  res@xyLineColor        = -1                           ; We don't want the line, so make it transparent.
  bot_plot  = gsn_csm_xy (wks,x,GRACE_global_ave_year(1:2,:),res)  ; Create filled XY plot.

; Create a plot with the area between both curves filled in pink.
  delete(res@gsnXYFillColors)
  res@xyLineColors       = LineColors
  res@xyLineThicknesses	 = LineThicknesses
  res@xyDashPatterns	   = DashPatterns
  mid_plot  = gsn_csm_xy (wks,x,TWS_global_ave_year,res)  ; Create another filled XY plot.

  overlay(bot_plot,mid_plot)
  overlay(bot_plot,top_plot)
;  maximize_output(wks,pres)

end if
end if

;*****************************************************
; Manually create and attach legend
;*****************************************************
dum1 = new(6,graphic)
dum2 = new(6,graphic)

res_text                = True                  ; text mods desired
res_text@txFontHeightF  = 0.013                ; change text size
res_text@txJust         = "CenterLeft"          ; text justification
res_lines               = True                  ; polyline mods desired

do case_num = 0,case_sum
  if case_num .lt. case_sum then
     res_lines@gsLineDashPattern = DashPatterns(case_num)
     res_lines@gsLineThicknessF  = LineThicknesses(case_num)
     res_lines@gsLineColor       = LineColors(case_num)
  else
     res_lines@gsLineDashPattern = 0.
     res_lines@gsLineThicknessF  = 2.
     res_lines@gsLineColor       = "blue"                                 ; change to blue
  end if
  
  if case_num .eq. 5 then
     point_xl = 1
     point_xr = 5
     point_xt = 5
     point_y  = 65
  else if case_num .ge. 0 .and. case_num .le. 2 then
     point_xl = (case_num+1)*24+1
     point_xr = (case_num+1)*24+6
     point_xt = (case_num+1)*24+6
     point_y = 65
  else
     point_xl = (case_num-3)*24+1
     point_xr = (case_num-3)*24+6
     point_xt = (case_num-3)*24+6
     point_y = 58
  end if
  end if
  xx = (/point_xl,point_xr/)
  yy = (/point_y,point_y/)
  dum1(case_num) = gsn_add_polyline(wks,mid_plot,xx,yy,res_lines)              ; add polyline
  if case_num .lt. case_sum then
     dum2(case_num) = gsn_add_text(wks,mid_plot,case_name_brf(case_num),point_xt,point_y,res_text); add text
  else
     dum2(case_num) = gsn_add_text(wks,mid_plot,"GRACE",point_xt,point_y,res_text); add text
  end if
end do
maximize_output(wks,pres)

res1 = True
res1@gsnMaximize        = False

res1@vpHeightF          = 0.4
res1@vpWidthF           = 0.3
;res1@vpXF               = 0.2
;res1@vpYF               = 0.9
;res@pmLegendWidthF     = 0.30                      ; Resize legend width
;res@pmLegendHeightF    = 0.15                     ; and height
res1@pmLegendSide       = "Right"
;res1@tmYUseLeft         = False
res1@tmXBLabelFontThicknessF = 0.015
res1@tmYRLabelFontThicknessF = 0.015
res1@tmXBLabelFontHeightF = 0.015
res1@tmYLLabelFontHeightF = 0.015
res1@tiMainFontHeightF = 0.025
res1@txFontHeightF	 = 0.025            ; change title font heights
res1@tmYLMode  = "Explicit"
res1@tmYLValues=(/45.,30,15,0,-15,-30,-45/)
res1@tmYLLabels=(/"45~S~o~N~N","30~S~o~N~N","15~S~o~N~N","0~S~o~N~","15~S~o~N~S","30~S~o~N~S","45~S~o~N~S"/)
res1@xyMonoLineColor    = False   
res1@xyLineColors = (/"darkgreen","green","red","pink","orange","blue"/) ; colors chosen
res1@xyLineThicknesses  = (/1.,1.,1.,1.,1.,1./) 
res1@xyDashPatterns       = (/0.,0.,0.,0.,0.,0./)

res1@trYMinF = -45.
res1@trYMaxF = 45.

res1@tiMainString   = "(b)"
plot = gsn_csm_xy(wks,TWS_zonal_ave_multiyear,f->y({A:B}),res1) ;_anomaly


;;*****************************************************
;; Manually create and attach titles
;;*****************************************************
;
;;
;; Attach some titles at the top.
;;
;  res_text               = True
;  res_text@txFontHeightF = 0.03                       ; change font size
;  txid_top = gsn_create_text(wks, "Parallel Climate Model Ensembles", res_text)
;
;  amres                  = True
;  amres@amJust           = "BottomCenter"
;  amres@amParallelPosF   =  0.0    ; This is the center of the plot.
;  amres@amOrthogonalPosF = -0.72   ; This is above the top edge of the plot.
;  annoid_top = gsn_add_annotation(bot_plot, txid_top, amres)
;
;  res_text@txFontHeightF = 0.02                       ; change font size
;  txid_mid = gsn_create_text(wks, "Global Temperature Anomalies",res_text)
;
;  amres@amOrthogonalPosF = -0.62  ; This is just below the previous title.
;  annoid_mid = gsn_add_annotation(bot_plot, txid_mid, amres)
;
;  res_text@txFontHeightF = 0.015                      ; change font size
;  txid_bot = gsn_create_text(wks,"from 1890-1919 average",res_text)
;
;  amres@amOrthogonalPosF = -0.55  ; This is just below the previous title.
;  annoid_bot = gsn_add_annotation(bot_plot, txid_bot, amres)
;
;  pres = True
;  maximize_output(wks,pres)

end
