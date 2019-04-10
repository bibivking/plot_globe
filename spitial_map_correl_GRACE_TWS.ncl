;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-Haverd2013/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs-Haverd2013/gw_on"\
             /)

             ;"/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs/gw_off",\
             ;"/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para_redistrb-on/Outputs/gw_off"\

             ;"/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs/gw_on",\
             ;"/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs/gw_on",\

var_name= (/\
           "Evap","TVeg","ESoil",\ ;"kg/m^2/s"
           "Qs","Qsb",\      ;"kg/m^2/s"
           "Qle","Qh",\      ;"W/m^2"
           "GPP","NPP",\     ;"umol/m^2/s"
           "RadT","VegT",\   ; "K"
           "SoilMoist",\     ;"m^3/m^3"  Average layer soil moisture
           "GWMoist",\       ;"mm3/mm3" Aquifer mositure content, Groundwater_off doesn't have
           "WatTable"\       ;"m"
           /)
           ;"LAI" ; the same in Groundwater_on and Groundwater_off

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872, 20./) ; 6 soil layers and 1 aquifer
case_compare = "gw_on_Haverd2013_vs_gw_norecharge_Haverd2013"
;month_name= (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
;season_name= (/"Spr","Sum","Aut","Win"/)
time_scale = "all" ;"annual" ; "seasonal", "monthly", "all"

;print(time_scale .eq. "seasonal")

;______________________________ Data _________________________________
year_s = 2001
year_e = 2010
year_sum = year_e-year_s+1
var_sum  = dimsizes(var_name)


var = new((/3,17,11,360,720/),float)
var = 0.0
Moist = new((/3,17,7,360,720/),float)
Moist = 0.0
TotSM   = new((/3,17,360,720/),float)
TotSM   = 0.0

do case_num = 0,1
print(case_name(case_num))

do year = year_s,year_e

    print("year = " +year)
    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if

filename = case_name(case_num)+"/cable_out_"+year+".nc"

f = addfile (filename,"r")

var!0 = "case"
var!1 = "time"
var!2 = "var"
var!3 = "lat"
var!4 = "lon"
var&lat = f->Qle&y
var&lon = f->Qle&x

Moist!0 = "case"
Moist!1 = "time"
Moist!2 = "layer"
Moist!3 = "lat"
Moist!4 = "lon"
Moist&lat = f->Qle&y
Moist&lon = f->Qle&x

TotSM!0   = "case"
TotSM!1   = "time"
TotSM!2   = "lat"
TotSM!3   = "lon"
TotSM&lat = f->Qle&y
TotSM&lon = f->Qle&x

var@_FillValue   = f->Qle@_FillValue
Moist@_FillValue = f->Qle@_FillValue
TotSM@_FillValue   = f->Qle@_FillValue

do var_num = 0,10
   if (var_num .le. 4) then
     do month = 0,11
        var(case_num,month,var_num,:,:) = var(case_num,month,var_num,:,:) + f->$var_name(var_num)$(month,:,:)*3600*24*dom(month)
     end do
   else if (var_num .eq. 9) .or. (var_num .eq. 10) then
      var(case_num,0:11,var_num,:,:) = var(case_num,0:11,var_num,:,:) + (f->$var_name(var_num)$ - 273.16)
   else
      var(case_num,0:11,var_num,:,:) = var(case_num,0:11,var_num,:,:) + f->$var_name(var_num)$
   end if
   end if
end do ;var_num

Moist(case_num,0:11,0:5,:,:) = Moist(case_num,0:11,0:5,:,:) + f->SoilMoist
Moist(case_num,0:11,6,:,:)   = Moist(case_num,0:11,6,:,:)   + f->GWMoist


do layers = 0,5
;  m3/m3
TotSM(case_num,0:11, :, :) = TotSM(case_num,0:11, :, :)+ (Moist(case_num,0:11,layers,:,:)*Soil_thickness(layers))/sum(Soil_thickness)
end do

end do ;year

end do ;case

; seasonal
;spring
var(0:1,12,:,:,:)   = (var(0:1,2,:,:,:) + var(0:1,3,:,:,:) + var(0:1,4,:,:,:))/3.
Moist(0:1,12,:,:,:) = (Moist(0:1,2,:,:,:) + Moist(0:1,3,:,:,:) + Moist(0:1,4,:,:,:))/3.
TotSM(0:1,12,:,:)   = (TotSM(0:1,2,:,:) + TotSM(0:1,3,:,:) + TotSM(0:1,4,:,:))/3.
; summer
var(0:1,13,:,:,:)   = (var(0:1,5,:,:,:) + var(0:1,6,:,:,:) + var(0:1,7,:,:,:))/3.
Moist(0:1,13,:,:,:) = (Moist(0:1,5,:,:,:) + Moist(0:1,6,:,:,:) + Moist(0:1,7,:,:,:))/3.
TotSM(0:1,13,:,:)   = (TotSM(0:1,5,:,:) + TotSM(0:1,6,:,:) + TotSM(0:1,7,:,:))/3.
; autumn
var(0:1,14,:,:,:)   = (var(0:1,8,:,:,:) + var(0:1,9,:,:,:) + var(0:1,10,:,:,:))/3.
Moist(0:1,14,:,:,:) = (Moist(0:1,8,:,:,:) + Moist(0:1,9,:,:,:) + Moist(0:1,10,:,:,:))/3.
TotSM(0:1,14,:,:)   = (TotSM(0:1,8,:,:) + TotSM(0:1,9,:,:) + TotSM(0:1,10,:,:))/3.
; winter
var(0:1,15,:,:,:)   = (var(0:1,0,:,:,:) + var(0:1,1,:,:,:) + var(0:1,11,:,:,:))/3.
Moist(0:1,15,:,:,:) = (Moist(0:1,0,:,:,:) + Moist(0:1,1,:,:,:) + Moist(0:1,11,:,:,:))/3.
TotSM(0:1,15,:,:)   = (TotSM(0:1,0,:,:) + TotSM(0:1,1,:,:) + TotSM(0:1,11,:,:))/3.

;annual
var(0:1,16,:,:,:)   = dim_avg_n_Wrap(var(0:1,0:11,:,:,:),1)
Moist(0:1,16,:,:,:) = dim_avg_n_Wrap(Moist(0:1,0:11,:,:,:),1)
TotSM(0:1,16,:,:)   = dim_avg_n_Wrap(TotSM(0:1,0:11,:,:),1)

var(2,:,:,:,:)   = var(1,:,:,:,:) - var(0,:,:,:,:)
Moist(2,:,:,:,:) = Moist(1,:,:,:,:) - Moist(0,:,:,:,:)
TotSM(2,:,:,:)   = TotSM(1,:,:,:) - TotSM(0,:,:,:)

var   = var/year_sum
Moist = Moist/year_sum
TotSM = TotSM/year_sum

;_____________________________ PLOT ___________________________________
pic = "GW_on_norecharge-"+time_scale+"-"+case_compare
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

res                    = True
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                            ;��ɫ
res@tmBorderThicknessF = 3.0

res@gsnDraw            = False  ; Don't draw plots
res@gsnFrame           = False  ; ����ҳ
res@lbLabelBarOn       = True

;************** ����labelbar ***************
res@lbBoxLinesOn       = True                       ;�ر�lbar box ����
res@lbTitleFont        = 25
res@lbLabelFont        = 25
res@lbTitleFontHeightF = 0.013
res@lbLabelFontHeightF = 0.013
res@txString           = ""
res@tmXBLabelFont      = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont      = 25

;*************** ���õ�ֵ�� *****************
res@cnLinesOn          = False                       ; �رյ�ֵ������
res@cnLineColor        = "black"
res@cnLineThicknessF   = 1.5
res@cnLineLabelsOn     = False
res@gsnMaximize        = True
res@cnExplicitLabelBarLabelsOn = True   ;?
res@cnLevelSelectionMode = "ExplicitLevels"

;************ ����������ֵ��Χ **************
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015

res@trYMinF            = -60
res@trYMaxF            = 80
res@mpMinLatF          = -60
res@mpMaxLatF          = 80
res@mpCenterLonF       = 160

res@tmYLMode  = "Explicit"             ;��Y����������������,��ʾ�Լ���Ҫ��ֵ
res@tmYLValues=(/60,30,0,-30,-60/)
res@tmYLLabels=(/"60~S~o~N~N","30~S~o~N~N","0~S~o~N~","30~S~o~N~S","60~S~o~N~S"/)
res@tmXBMode  = "Explicit"
res@tmXBValues=(/240,300,0,60,120,180/)
res@tmXBLabels=(/"120~S~o~N~W","60~S~o~N~W","0~S~o~N~","60~S~o~N~E","120~S~o~N~E","180~S~o~N~"/)

;**************** ͼ�������� *****************
pres                    = True                                      ; ���� panel resources.
pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
pres@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
pres@lbBoxLinesOn       = True                       ;�ر�lbar box ����
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
pres@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
pres@lbTitleFontHeightF = 0.013
pres@lbLabelFontHeightF = 0.013
pres@lbAutoManage       = False
pres@lbTitleDirection   = "Across"             ; title direction
pres@lbTitlePosition    = "Bottom";Right"
;pres@pmLabelBarSide     = "Bottom";"Right"
;pres@pmLabelBarOrthogonalPosF = -0.021
;pres@pmLabelBarParallelPosF = 0.

plots = new(3,graphic)

; "Evap","TVeg","ESoil",\ ;"kg/m^2/s"
do var_num = 0,2
  if (time_scale .eq. "monthly") .or. (time_scale .eq. "all") then
    do month = 0,11
      res@cnFillPalette = "WhiteBlueGreenYellowRed"
      pres@txString = f->$var_name(var_num)$@long_name+"-Month "+month+" (mm/month)"
      levels = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120./)
      res@cnLevels = levels
      res@tiMainString = "(a) GW_off"
      plots(0) = gsn_csm_contour_map(wks,var(0,month,var_num,:,:),res)
      res@tiMainString = "(b) GW-off"
      plots(1) = gsn_csm_contour_map(wks,var(1,month,var_num,:,:),res)
      delete(levels)
      delete(res@cnLevels)
      res@cnFillPalette = "ViBlGrWhYeOrRe"
      levels = (/-25.,-20.,-15.,-10.,-5.,-2.,-1.,1.,2.,5.,10.,15.,20.,25./)
      ;levels = (/-20.,-18.,-16.,-14.,-12.,-10.,-8.,-6.,-4.,-2.,2.,4.,6.,8.,10.,12.,14.,16.,18.,20./)
      res@cnLevels = levels
      res@tiMainString = "(c) GW_on - GW_off"
      plots(2) = gsn_csm_contour_map(wks,var(2,month,var_num,:,:),res)
      gsn_panel(wks,(/plots/),(/2,2/),pres)
      delete(levels)
      delete(res@cnLevels)
    end do
  end if
  if (time_scale .eq. "seasonal") .or. (time_scale .eq. "all") then
    do season = 0,3
      res@cnFillPalette = "WhiteBlueGreenYellowRed"
      pres@txString = f->$var_name(var_num)$@long_name+"-Season "+season+" (mm/month)"
      levels = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120./)
      res@cnLevels = levels
      res@tiMainString = "(a) GW_off"
      plots(0) = gsn_csm_contour_map(wks,var(0,11+season,var_num,:,:),res)
      res@tiMainString = "(b) GW-off"
      plots(1) = gsn_csm_contour_map(wks,var(1,11+season,var_num,:,:),res)
      delete(levels)
      delete(res@cnLevels)
      res@cnFillPalette = "ViBlGrWhYeOrRe"
      levels = (/-25.,-20.,-15.,-10.,-5.,-2.,-1.,1.,2.,5.,10.,15.,20.,25./)
      res@cnLevels = levels
      res@tiMainString = "(c) GW_on - GW_off"
      plots(2) = gsn_csm_contour_map(wks,var(2,11+season,var_num,:,:),res)
      gsn_panel(wks,(/plots/),(/2,2/),pres)
      delete(levels)
      delete(res@cnLevels)
    end do
  end if
  if (time_scale .eq. "annual") .or. (time_scale .eq. "all") then
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = f->$var_name(var_num)$@long_name+"-Annual (mm/month)"
    levels = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120./)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,var(0,16,var_num,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,var(1,16,var_num,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-25.,-20.,-15.,-10.,-5.,-2.,-1.,1.,2.,5.,10.,15.,20.,25./)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,var(2,16,var_num,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end if
end do

; "Qs","Qsb" ;"kg/m^2/s"
do var_num = 3,4
  if (time_scale .eq. "monthly") .or. (time_scale .eq. "all") then
    do month = 0,11
      res@cnFillPalette = "WhiteBlueGreenYellowRed"
      pres@txString = f->$var_name(var_num)$@long_name+"-Month "+month+" (mm/month)"
      levels = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120./)
      res@cnLevels = levels
      res@tiMainString = "(a) GW_off"
      plots(0) = gsn_csm_contour_map(wks,var(0,month,var_num,:,:),res)
      res@tiMainString = "(b) GW-off"
      plots(1) = gsn_csm_contour_map(wks,var(1,month,var_num,:,:),res)
      delete(levels)
      delete(res@cnLevels)
      res@cnFillPalette = "ViBlGrWhYeOrRe"
      levels = (/-25.,-20.,-15.,-10.,-5.,-2.,-1.,1.,2.,5.,10.,15.,20.,25./)
      res@cnLevels = levels
      res@tiMainString = "(c) GW_on - GW_off"
      plots(2) = gsn_csm_contour_map(wks,var(2,month,var_num,:,:),res)
      gsn_panel(wks,(/plots/),(/2,2/),pres)
      delete(levels)
      delete(res@cnLevels)
    end do
  end if
  if (time_scale .eq. "seasonal") .or. (time_scale .eq. "all") then
    do season = 0,3
      res@cnFillPalette = "WhiteBlueGreenYellowRed"
      pres@txString = f->$var_name(var_num)$@long_name+"-Season "+season+" (mm/month)"
      levels = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120./)
      res@cnLevels = levels
      res@tiMainString = "(a) GW_off"
      plots(0) = gsn_csm_contour_map(wks,var(0,11+season,var_num,:,:),res)
      res@tiMainString = "(b) GW-off"
      plots(1) = gsn_csm_contour_map(wks,var(1,11+season,var_num,:,:),res)
      delete(levels)
      delete(res@cnLevels)
      res@cnFillPalette = "ViBlGrWhYeOrRe"
      levels = (/-25.,-20.,-15.,-10.,-5.,-2.,-1.,1.,2.,5.,10.,15.,20.,25./)
      res@cnLevels = levels
      res@tiMainString = "(c) GW_on - GW_off"
      plots(2) = gsn_csm_contour_map(wks,var(2,11+season,var_num,:,:),res)
      gsn_panel(wks,(/plots/),(/2,2/),pres)
      delete(levels)
      delete(res@cnLevels)
    end do
  end if
  if (time_scale .eq. "annual") .or. (time_scale .eq. "all") then
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = f->$var_name(var_num)$@long_name+"-Annual (mm/month)"
    levels = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120./)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,var(0,16,var_num,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,var(1,16,var_num,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-25.,-20.,-15.,-10.,-5.,-2.,-1.,1.,2.,5.,10.,15.,20.,25./)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,var(2,16,var_num,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end if
end do

; "Qle" "W/m^2"
if (time_scale .eq. "monthly") .or. (time_scale .eq. "all") then
  do month = 0,11
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = f->$var_name(5)$@long_name+"-Month "+month+" ("+f->$var_name(5)$@units+")"
    levels = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120./)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,var(0,month,5,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,var(1,month,5,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-25.,-20.,-15.,-10.,-5.,-2.,-1.,1.,2.,5.,10.,15.,20.,25./)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,var(2,month,5,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end do
end if
if (time_scale .eq. "seasonal") .or. (time_scale .eq. "all") then
  do season = 0,3
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = f->$var_name(5)$@long_name+"-Season "+season+" ("+f->$var_name(5)$@units+")"
    levels = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120./)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,var(0,11+season,5,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,var(1,11+season,5,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-25.,-20.,-15.,-10.,-5.,-2.,-1.,1.,2.,5.,10.,15.,20.,25./)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,var(2,11+season,5,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end do
end if
if (time_scale .eq. "annual") .or. (time_scale .eq. "all") then
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  pres@txString = f->$var_name(5)$@long_name+"-Annual"+" ("+f->$var_name(5)$@units+")"
  levels = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW_off"
  plots(0) = gsn_csm_contour_map(wks,var(0,16,5,:,:),res)
  res@tiMainString = "(b) GW-off"
  plots(1) = gsn_csm_contour_map(wks,var(1,16,5,:,:),res)
  delete(levels)
  delete(res@cnLevels)
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = (/-25.,-20.,-15.,-10.,-5.,-2.,-1.,1.,2.,5.,10.,15.,20.,25./)
  res@cnLevels = levels
  res@tiMainString = "(c) GW_on - GW_off"
  plots(2) = gsn_csm_contour_map(wks,var(2,16,5,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)
end if

; "Qh" ;"W/m^2"
if (time_scale .eq. "monthly") .or. (time_scale .eq. "all") then
  do month = 0,11
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = f->$var_name(6)$@long_name+"-Month "+month+" ("+f->$var_name(6)$@units+")"
    levels = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,110.,120.,130.,140.,150.,160,170.,180.,190.,200.,210.,220./)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,var(0,month,6,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,var(1,month,6,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-25.,-20.,-15.,-10.,-5.,-2.,-1.,1.,2.,5.,10.,15.,20.,25./)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,var(2,month,6,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end do
end if
if (time_scale .eq. "seasonal") .or. (time_scale .eq. "all") then
  do season = 0,3
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = f->$var_name(6)$@long_name+"-Season "+season+" ("+f->$var_name(6)$@units+")"
    levels = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,110.,120.,130.,140.,150.,160,170.,180.,190.,200.,210.,220./)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,var(0,11+season,6,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,var(1,11+season,6,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-25.,-20.,-15.,-10.,-5.,-2.,-1.,1.,2.,5.,10.,15.,20.,25./)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,var(2,11+season,6,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end do
end if
if (time_scale .eq. "annual") .or. (time_scale .eq. "all") then
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  pres@txString = f->$var_name(6)$@long_name+"-Annual"+" ("+f->$var_name(6)$@units+")"
  levels = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,110.,120.,130.,140.,150.,160,170.,180.,190.,200.,210.,220./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW_off"
  plots(0) = gsn_csm_contour_map(wks,var(0,16,6,:,:),res)
  res@tiMainString = "(b) GW-off"
  plots(1) = gsn_csm_contour_map(wks,var(1,16,6,:,:),res)
  delete(levels)
  delete(res@cnLevels)
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = (/-25.,-20.,-15.,-10.,-5.,-2.,-1.,1.,2.,5.,10.,15.,20.,25./)
  res@cnLevels = levels
  res@tiMainString = "(c) GW_on - GW_off"
  plots(2) = gsn_csm_contour_map(wks,var(2,16,6,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)
end if

; "GPP"
if (time_scale .eq. "monthly") .or. (time_scale .eq. "all") then
  do month = 0,11
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = f->GPP@long_name+"-Month "+month+" ("+f->GPP@units+")"
    levels = (/1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12./)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,var(0,month,7,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,var(1,month,7,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1./)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,var(2,month,7,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end do
end if
if (time_scale .eq. "seasonal") .or. (time_scale .eq. "all") then
  do season = 0,3
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = f->GPP@long_name+"-Season "+season+" ("+f->GPP@units+")"
    levels = (/1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12./)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,var(0,11+season,7,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,var(1,11+season,7,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1./)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,var(2,11+season,7,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end do
end if
if (time_scale .eq. "annual") .or. (time_scale .eq. "all") then
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  pres@txString = f->GPP@long_name+"-Annual"+" ("+f->GPP@units+")"
  levels = (/1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW_off"
  plots(0) = gsn_csm_contour_map(wks,var(0,16,7,:,:),res)
  res@tiMainString = "(b) GW-off"
  plots(1) = gsn_csm_contour_map(wks,var(1,16,7,:,:),res)
  delete(levels)
  delete(res@cnLevels)
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1./)
  res@cnLevels = levels
  res@tiMainString = "(c) GW_on - GW_off"
  plots(2) = gsn_csm_contour_map(wks,var(2,16,7,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)
end if

;"NPP"
if (time_scale .eq. "monthly") .or. (time_scale .eq. "all") then
  do month = 0,11
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = f->NPP@long_name+"-Month "+month+" ("+f->NPP@units+")"
    levels = (/-5.,-4.,-3.,-2.,-1.,1.,2.,3.,4.,5./)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,var(0,month,8,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,var(1,month,8,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1./)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,var(2,month,8,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end do
end if
if (time_scale .eq. "seasonal") .or. (time_scale .eq. "all") then
  do season = 0,3
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = f->NPP@long_name+"-Season "+season+" ("+f->NPP@units+")"
    levels = (/-5.,-4.,-3.,-2.,-1.,1.,2.,3.,4.,5./)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,var(0,11+season,8,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,var(1,11+season,8,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1./)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,var(2,11+season,8,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end do
end if
if (time_scale .eq. "annual") .or. (time_scale .eq. "all") then
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  pres@txString = f->NPP@long_name+"-Annual"+" ("+f->NPP@units+")"
  levels = (/-5.,-4.,-3.,-2.,-1.,1.,2.,3.,4.,5./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW_off"
  plots(0) = gsn_csm_contour_map(wks,var(0,16,8,:,:),res)
  res@tiMainString = "(b) GW-off"
  plots(1) = gsn_csm_contour_map(wks,var(1,16,8,:,:),res)
  delete(levels)
  delete(res@cnLevels)
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1./)
  res@cnLevels = levels
  res@tiMainString = "(c) GW_on - GW_off"
  plots(2) = gsn_csm_contour_map(wks,var(2,16,8,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)
end if

; "RadT"
if (time_scale .eq. "monthly") .or. (time_scale .eq. "all") then
  do month = 0,11
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = f->RadT@long_name+"-Month "+month+" (~S~o~N~C)"
    levels = (/5.,10.,15.,20.,25.,30.,35.,40./)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,var(0,month,9,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,var(1,month,9,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6/)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,var(2,month,9,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end do
end if
if (time_scale .eq. "seasonal") .or. (time_scale .eq. "all") then
  do season = 0,3
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = f->RadT@long_name+"-Season "+season+" (~S~o~N~C)"
    levels = (/5.,10.,15.,20.,25.,30.,35.,40./)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,var(0,11+season,9,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,var(1,11+season,9,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6/)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,var(2,11+season,9,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end do
end if
if (time_scale .eq. "annual") .or. (time_scale .eq. "all") then
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  pres@txString = f->RadT@long_name+"-Annual"+" (~S~o~N~C)"
  levels = (/5.,10.,15.,20.,25.,30.,35.,40./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW_off"
  plots(0) = gsn_csm_contour_map(wks,var(0,16,9,:,:),res)
  res@tiMainString = "(b) GW-off"
  plots(1) = gsn_csm_contour_map(wks,var(1,16,9,:,:),res)
  delete(levels)
  delete(res@cnLevels)
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = (/-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6/)
  res@cnLevels = levels
  res@tiMainString = "(c) GW_on - GW_off"
  plots(2) = gsn_csm_contour_map(wks,var(2,16,9,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)
end if


;"VegT"
if (time_scale .eq. "monthly") .or. (time_scale .eq. "all") then
  do month = 0,11
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = f->VegT@long_name+"-Month "+month+" (~S~o~N~C)"
    levels = (/5.,10.,15.,20.,25.,30.,35.,40./)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,var(0,month,10,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,var(1,month,10,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6/)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,var(2,month,10,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end do
end if
if (time_scale .eq. "seasonal") .or. (time_scale .eq. "all") then
  do season = 0,3
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = f->VegT@long_name+"-Season "+season+" (~S~o~N~C)"
    levels = (/5.,10.,15.,20.,25.,30.,35.,40./)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,var(0,11+season,10,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,var(1,11+season,10,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6/)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,var(2,11+season,10,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end do
end if
if (time_scale .eq. "annual") .or. (time_scale .eq. "all") then
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  pres@txString = f->VegT@long_name+"-Annual"+" (~S~o~N~C)"
  levels = (/5.,10.,15.,20.,25.,30.,35.,40./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW_off"
  plots(0) = gsn_csm_contour_map(wks,var(0,16,10,:,:),res)
  res@tiMainString = "(b) GW-off"
  plots(1) = gsn_csm_contour_map(wks,var(1,16,10,:,:),res)
  delete(levels)
  delete(res@cnLevels)
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = (/-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6/)
  res@cnLevels = levels
  res@tiMainString = "(c) GW_on - GW_off"
  plots(2) = gsn_csm_contour_map(wks,var(2,16,10,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)
end if

do layers = 0,6
  if (time_scale .eq. "monthly") .or. (time_scale .eq. "all") then
    do month = 0,11
      res@cnFillPalette = "WhiteBlueGreenYellowRed"
      pres@txString = "The soil moisture of the "+(layers+1)+ " layer-Month "+month+" (m~S~3~N~/m~S~3~N~)"
      levels = (/0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.60,0.65,0.70/)
      res@cnLevels = levels
      res@tiMainString = "(a) GW_off"
      plots(0) = gsn_csm_contour_map(wks,Moist(0,month,layers,:,:),res)
      res@tiMainString = "(b) GW-off"
      plots(1) = gsn_csm_contour_map(wks,Moist(1,month,layers,:,:),res)
      delete(levels)
      delete(res@cnLevels)
      res@cnFillPalette = "ViBlGrWhYeOrRe"
      levels = (/-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,-0.02,0.02,0.05,0.1,0.15,0.2,0.25,0.3/)
      res@cnLevels = levels
      res@tiMainString = "(c) GW_on - GW_off"
      plots(2) = gsn_csm_contour_map(wks,Moist(2,month,layers,:,:),res)
      gsn_panel(wks,(/plots/),(/2,2/),pres)
      delete(levels)
      delete(res@cnLevels)
    end do
  end if
  if (time_scale .eq. "seasonal") .or. (time_scale .eq. "all") then
    do season = 0,3
      res@cnFillPalette = "WhiteBlueGreenYellowRed"
      pres@txString = "The soil moisture of the "+(layers+1)+ " layer-Season "+season+" (m~S~3~N~/m~S~3~N~)"
      levels = (/0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.60,0.65,0.70/)
      res@cnLevels = levels
      res@tiMainString = "(a) GW_off"
      plots(0) = gsn_csm_contour_map(wks,Moist(0,11+season,layers,:,:),res)
      res@tiMainString = "(b) GW-off"
      plots(1) = gsn_csm_contour_map(wks,Moist(1,11+season,layers,:,:),res)
      delete(levels)
      delete(res@cnLevels)
      res@cnFillPalette = "ViBlGrWhYeOrRe"
      levels = (/-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,-0.02,0.02,0.05,0.1,0.15,0.2,0.25,0.3/)
      res@cnLevels = levels
      res@tiMainString = "(c) GW_on - GW_off"
      plots(2) = gsn_csm_contour_map(wks,Moist(2,11+season,layers,:,:),res)
      gsn_panel(wks,(/plots/),(/2,2/),pres)
      delete(levels)
      delete(res@cnLevels)
    end do
  end if
  if (time_scale .eq. "annual") .or. (time_scale .eq. "all") then
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = "The soil moisture of the "+(layers+1)+ " layer-Annual (m~S~3~N~/m~S~3~N~)"
    levels = (/0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.60,0.65,0.70/)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,Moist(0,16,layers,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,Moist(1,16,layers,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,-0.02,0.02,0.05,0.1,0.15,0.2,0.25,0.3/)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,Moist(2,16,layers,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end if
end do


if (time_scale .eq. "monthly") .or. (time_scale .eq. "all") then
  do month = 0,11
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = "Total Soil Moist-Month "+month+" (m~S~3~N~/m~S~3~N~)"
    levels = (/0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.60,0.65,0.70/)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,TotSM(0,month,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,TotSM(1,month,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,-0.02,0.02,0.05,0.1,0.15,0.2,0.25,0.3/)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,TotSM(2,month,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end do
end if
if (time_scale .eq. "seasonal") .or. (time_scale .eq. "all") then
  do season = 0,3
    res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = "Total Soil Moist-Season "+season+" (m~S~3~N~/m~S~3~N~)"
    levels = (/0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.60,0.65,0.70/)
    res@cnLevels = levels
    res@tiMainString = "(a) GW_off"
    plots(0) = gsn_csm_contour_map(wks,TotSM(0,11+season,:,:),res)
    res@tiMainString = "(b) GW-off"
    plots(1) = gsn_csm_contour_map(wks,TotSM(1,11+season,:,:),res)
    delete(levels)
    delete(res@cnLevels)
    res@cnFillPalette = "ViBlGrWhYeOrRe"
    levels = (/-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,-0.02,0.02,0.05,0.1,0.15,0.2,0.25,0.3/)
    res@cnLevels = levels
    res@tiMainString = "(c) GW_on - GW_off"
    plots(2) = gsn_csm_contour_map(wks,TotSM(2,11+season,:,:),res)
    gsn_panel(wks,(/plots/),(/2,2/),pres)
    delete(levels)
    delete(res@cnLevels)
  end do
end if
if (time_scale .eq. "annual") .or. (time_scale .eq. "all") then
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
    pres@txString = "Total Soil Moist-Annual (m~S~3~N~/m~S~3~N~)"
  levels = (/0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.60,0.65,0.70/)
  res@cnLevels = levels
  res@tiMainString = "(a) GW_off"
  plots(0) = gsn_csm_contour_map(wks,TotSM(0,16,:,:),res)
  res@tiMainString = "(b) GW-off"
  plots(1) = gsn_csm_contour_map(wks,TotSM(1,16,:,:),res)
  delete(levels)
  delete(res@cnLevels)
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = (/-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,-0.02,0.02,0.05,0.1,0.15,0.2,0.25,0.3/)
  res@cnLevels = levels
  res@tiMainString = "(c) GW_on - GW_off"
  plots(2) = gsn_csm_contour_map(wks,TotSM(2,16,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)
end if

end
