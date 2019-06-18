;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs-standard/gw_off/"\
             /)
             
var_name= (/\
           "Evap","TVeg","ESoil",\ ;"kg/m^2/s"
           "Qs","Qsb",\      ;"kg/m^2/s"
           "SoilMoist",\     ;"m^3/m^3"  Average layer soil moisture
           "GWMoist"\       ;"mm3/mm3" Aquifer mositure content, Groundwater_off doesn't have
           /)
           ;"LAI" ; the same in Groundwater_on and Groundwater_off

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers

; ______________________________ Pre-define ________________________
year_s   = 2003
year_e   = 2010
year_sum = year_e-year_s+1
total    = year_sum*12
case_sum = dimsizes(case_name)

; ____________________________ Options _______________________________
use_month    = "May" ; "Jul"
no_greenland = True ;False ;True
;time_scale   = "monthly" ; "annual"; "monthly"
grace_version= "Original";"Gap-fill"; "Original"
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
GRACE      = new((/total,360,720/),float)
uncertainty= new((/total,360,720/),float)
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
TWS     = new((/case_sum,total,360,720/),float)
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
do time = 0,total-1
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

; ____________________________ Deseasonalize & Detrend ________________________________
print("deseasonlize and detrend")
GRACE = rmMonAnnCycTLL(GRACE) ; dimensioned time x lat x lon
do case_num = 0,case_sum-1
   TWS(case_num,:,:,:) = rmMonAnnCycTLL(TWS(case_num,:,:,:))
end do

plot = gsn_csm_contour(wks,GRACE(0,:,:),res)
plot = gsn_csm_contour(wks,TWS(0,0,:,:),res)

GRACE = dtrend_n(GRACE, True, 0)*1.0
TWS   = dtrend_n(TWS,   True, 1)*1.0
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

delete(GRACE)
delete(TWS)
delete(uncertainty)

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
GRACE_zonal_ave_multiyear = new(lat_num,float)
TWS_zonal_ave_multiyear   = new((/case_sum,lat_num/),float)
GRACE_zonal_ave_multiyear = 0.0
TWS_zonal_ave_multiyear   = 0.0

GRACE_zonal_ave_multiyear = dim_avg_n_Wrap(GRACE_zonal_ave({A:B},:),1)
do case_num = 0,case_sum-1
   TWS_zonal_ave_multiyear(case_num,:) = dim_avg_n_Wrap(TWS_zonal_ave(case_num,{A:B},:),1)
end do

print("finish zonal average")

;_____________________________ PLOT ___________________________________
if no_greenland .eq. True then
   pic = "time_serial_zonal_avg_GRACE_vs_CABLE_no_Greenland_GRACE06-"+grace_version+"_"+scale
else
   pic = "time_serial_zonal_avg_GRACE_vs_CABLE_GRACE06-"+grace_version+"_"+scale
end if

wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"MPL_PuOr");"MPL_RdBu")

res                    = True

res@gsnDraw            = True  ; Don't draw plots
res@gsnFrame           = True  ; ?????
res@gsnMaximize        = True

res@cnFillOn           = True
res@cnFillMode         = "RasterFill"            ; Raster Mode                           
res@tmBorderThicknessF = 3.0
res@lbLabelBarOn       = True

;************** ????labelbar ***************
res@lbBoxLinesOn       = False  
res@lbTitleFont        = 25
res@lbLabelFont        = 25
res@lbTitleFontHeightF = 0.013
res@lbLabelFontHeightF = 0.013
res@tmXBLabelFont      = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont      = 25

;*************** ???????? *****************
res@cnLinesOn          = False;True                       ; ???????????
res@cnLineColor        = "black"
res@cnLineThicknessF   = 1.0
res@cnLineLabelsOn     = False
;res@cnExplicitLabelBarLabelsOn = True   ;?
res@cnLevelSelectionMode = "ExplicitLevels"
;res@tmXUseBottom        = False
;************ ?????????????¦¶ **************
res@tmXTLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXTLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015
res@tmYLMode  = "Explicit"
res@tmXBMode  = "Explicit"

res@trYMinF   = A
res@trYMaxF   = B
res@trXMinF   = 1
res@trXMaxF   = 96

res@tmYLValues=(/45.,30,15,0,-15,-30,-45/)
res@tmYLLabels=(/"45~S~o~N~N","30~S~o~N~N","15~S~o~N~N","0~S~o~N~","15~S~o~N~S","30~S~o~N~S","45~S~o~N~S"/)
res@tmXBValues=(/1,13,25,37,49,61,73,85/)
res@tmXBLabels=(/"2003","2004","2005","2006","2007","2008","2009","2010"/)

pres                    = True                                      ; ???? panel resources.
pres@gsnMaximize        = True;False
pres@gsnPanelLabelBar   = False                    ; ???? panel labelbar.
pres@lbBoxLinesOn       = False                       ;???lbar box ????
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
pres@lbLabelBarOn       = False ; True                    ; ??????plot??labelbar
pres@lbTitleFontHeightF = 0.5
pres@lbLabelFontHeightF = 0.5
pres@lbAutoManage       = False
pres@lbTitleDirection   = "Across"             ; title direction
pres@lbTitlePosition    = "Right";"Left";"Bottom";Right"


print(GRACE_zonal_ave&time);VarSummary
print(GRACE_zonal_ave&lat)
plots = new(2,graphic)

res@tiYAxisString	 = "Latitude (~S~o~N~)"              ; add an axis title
;res@tiXAxisString	 = "Year"              ; add an axis title
res@txFontHeightF	 = 0.015            ; change title font heights
levels = fspan(-10.,10.,101)
res@cnLevels = levels
;pres@cnLevels = levels
res@tiMainString   = "(c)"; Total Water Storage Anormaly"
plots(0) = gsn_csm_contour(wks,GRACE_zonal_ave,res)
print("finish printing plot 1")

res@tiMainString   = "(d)"; Total Water Storage Anormaly"
plots(1) = gsn_csm_contour(wks,TWS_zonal_ave(0,:,:),res)

gsn_panel(wks,(/plots/),(/1,2/),pres)


;res1 = True
;res1@gsnMaximize        = False
;
;res1@vpHeightF          = 0.7
;res1@vpWidthF           = 0.3
;res1@vpXF               = 0.2
;res1@vpYF               = 0.9
;;res@pmLegendWidthF     = 0.30                      ; Resize legend width
;;res@pmLegendHeightF    = 0.15                     ; and height
;res1@pmLegendSide       = "Right"
;res1@tmYUseLeft         = False
;res1@tmXBLabelFontThicknessF = 0.015
;res1@tmYRLabelFontThicknessF = 0.015
;res1@tmXBLabelFontHeightF = 0.015
;res1@tmYRLabelFontHeightF = 0.015
;res1@tmYRMode  = "Explicit"
;res1@tmYRValues=(/45.,30,15,0,-15,-30,-45/)
;res1@tmYRLabels=(/"45~S~o~N~N","30~S~o~N~N","15~S~o~N~N","0~S~o~N~","15~S~o~N~S","30~S~o~N~S","45~S~o~N~S"/)
;
;res1@tiMainString   = "(c) Total Water Storage Anormaly"
;plots(0) = gsn_csm_xy(wks,GRACE_zonal_ave_multiyear,f->y({A:B}),res1) ;_anomaly
;
;res1@tiMainString   = "(e) Total Water Storage Anormaly"
;plots(1) = gsn_csm_xy(wks,TWS_zonal_ave_multiyear(0,:),f->y({A:B}),res1) 
;
;gsn_panel(wks,(/plots/),(/1,2/),pres)





end


