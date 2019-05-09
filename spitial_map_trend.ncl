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
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off/Outputs/gw_off/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs-standard/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para_redistrb-on/Outputs-standard/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs-Haverd2013/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on/Outputs/gw_on/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-standard/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs-standard/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-Haverd2013/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs-Haverd2013/gw_on"\
             /)

case_name_brf = (/\
              "off-default-std",\
              "off-std",\
              "off-redstrb-std",\
              "off-Hvrd",\
              "on-default-std",\
              "on-std",\
              "on_norchrg-std",\
              "on-Hvrd",\
              "on_norchrg-Hvrd"\
               /)
               
Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)

; ______________________________ Pre-define ________________________
year_s   = 2003
year_e   = 2010
year_sum = year_e-year_s+1
total    = year_sum*12
case_sum = dimsizes(case_name)

; ____________________________ Options _______________________________
scale        = "World" ; "AU"
use_month    = "May" ; "Jul"
no_greenland = False ;True
grace_version= "Gap-fill" ;"Gap-fill"; "Original"

; _____________________________ Check ______________________________
pic = "check"
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"ViBlGrWhYeOrRe")
res = True
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                            ;???
res@tmBorderThicknessF = 3.0
res@cnLinesOn          = False                       ; ???????????
res@gsnDraw            = False  ; Don't draw plots
res@gsnFrame           = False  ; ?????


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

; Calculate 2003-2009 TWS average
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
      TWS(case_num,time,:,:) = where(landmask .eq. 1, TWS(case_num,time,:,:), TWS@_FillValue)
   end do
end do
delete(landmask)
end if

GRACE       = where( .not. ismissing(TWS(0,:,:,:)), GRACE, GRACE@_FillValue)
uncertainty = where( .not. ismissing(TWS(0,:,:,:)), uncertainty, uncertainty@_FillValue)

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


;GRACE & TWS is anormaly here

; _________________________ GRACE CABLE Trend _________________________
x = ispan(1,total,1)
GRACE_trend = regCoef(x,GRACE(lat|:,lon|:,time|:))
TWS_trend   = regCoef(x,TWS(case|:,lat|:,lon|:,time|:))            
 
GRACE_trend!0 = "lat"
GRACE_trend!1 = "lon"
GRACE_trend&lat = f->Qle&y
GRACE_trend&lon = f->Qle&x

TWS_trend!0 = "case"
TWS_trend!1 = "lat"
TWS_trend!2 = "lon"
TWS_trend&lat = f->Qle&y
TWS_trend&lon = f->Qle&x


;_____________________________ PLOT ___________________________________
if no_greenland .eq. True .and. scale .eq. "World" then
   pic = "spitial_map_trend_GRACE_vs_CABLE_no_Greenland_GRACE06-"+grace_version+"_"+scale
else 
   pic = "spitial_map_trend_GRACE_vs_CABLE_GRACE06-"+grace_version+"_"+scale
end if
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"ViBlGrWhYeOrRe")
;"WhiteBlueGreenYellowRed")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;???????

res                    = True
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                            ;???
res@tmBorderThicknessF = 3.0

res@gsnDraw            = True  ; Don't draw plots
res@gsnFrame           = True  ; ?????
res@lbLabelBarOn       = False ;True

;************** ????labelbar ***************
res@lbBoxLinesOn       = True                       ;???lbar box ????
res@lbTitleFont        = 25
res@lbLabelFont        = 25
res@lbTitleFontHeightF = 0.013
res@lbLabelFontHeightF = 0.013
res@txString           = ""
res@tmXBLabelFont      = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont      = 25

;*************** ???????? *****************
res@cnLinesOn          = False                       ; ???????????
res@cnLineColor        = "black"
res@cnLineThicknessF   = 1.5
res@cnLineLabelsOn     = False
res@gsnMaximize        = True
res@cnExplicitLabelBarLabelsOn = True   ;?
res@cnLevelSelectionMode = "ExplicitLevels"

;************ ?????????????¦¶ **************
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015

if scale .eq. "World" then
   res@trYMinF            = -60
   res@trYMaxF            = 80
   res@mpMinLatF          = -60
   res@mpMaxLatF          = 80
   res@mpCenterLonF       = 160

   res@tmYLMode  = "Explicit"             ;??Y????????????????,????????????
   res@tmYLValues=(/60,30,0,-30,-60/)
   res@tmYLLabels=(/"60~S~o~N~N","30~S~o~N~N","0~S~o~N~","30~S~o~N~S","60~S~o~N~S"/)
   res@tmXBMode  = "Explicit"
   res@tmXBValues=(/240,300,0,60,120,180/)
   res@tmXBLabels=(/"120~S~o~N~W","60~S~o~N~W","0~S~o~N~","60~S~o~N~E","120~S~o~N~E","180~S~o~N~"/)
else if scale .eq. "AU" then
   res@gsnAddCyclic       = False
   res@trYMinF            = -40
   res@trYMaxF            = -27
   res@mpMinLatF          = -40
   res@mpMaxLatF          = -27
   res@trXMinF            = 140
   res@trXMaxF            = 155
   res@mpMinLonF          = 140
   res@mpMaxLonF          = 155

   res@tmYLMode  = "Explicit"             ;??Y????????????????,????????????
   res@tmYLValues=(/-30,-40/)
   res@tmYLLabels=(/"30~S~o~N~S","40~S~o~N~S"/)
   res@tmXBMode  = "Explicit"
   res@tmXBValues=(/140,150/)
   res@tmXBLabels=(/"140~S~o~N~E","150~S~o~N~E"/)
end if
end if


;**************** ????????? *****************
pres                    = True                                      ; ???? panel resources.
pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False                    ; ???? panel labelbar.
pres@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
pres@lbBoxLinesOn       = True                       ;???lbar box ????
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
pres@lbLabelBarOn       = False ; True                    ; ??????plot??labelbar
pres@lbTitleFontHeightF = 0.013
pres@lbLabelFontHeightF = 0.013
pres@lbAutoManage       = True
pres@lbTitleDirection   = "Across"             ; title direction
pres@lbTitlePosition    = "Bottom";Right"

plots = new(10,graphic)

levels = fspan(-20.,20.,21)
res@cnLevels = levels

res@tiMainString = "GRACE"
if scale .eq. "World" then
   plots(0) = gsn_csm_contour_map(wks,GRACE_trend,res)
else if scale .eq. "AU" then
   plots(0) = gsn_csm_contour_map(wks,GRACE_trend({-40:-27},{140:155}),res)
end if
end if  
  
do case_num = 0,case_sum-1
  res@tiMainString = "CABLE "+case_name_brf(case_num)
  if scale .eq. "World" then
    print(case_name_brf(case_num))
    plots(case_num+1) = gsn_csm_contour_map(wks,TWS_trend(case_num,:,:),res)
  else if scale .eq. "AU" then
    print(case_name_brf(case_num))
    plots(case_num+1) = gsn_csm_contour_map(wks,TWS_trend(case_num,{-40:-27},{140:155}),res)
  end if
  end if
end do

gsn_panel(wks,(/plots/),(/5,2/),pres)

end
