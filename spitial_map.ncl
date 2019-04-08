;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_new_run_on/Outputs/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_new_run_off/Outputs/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_new_run_on_norecharge/Outputs/gw_on"\
             /)

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


;______________________________ Data _________________________________
year_s = 1951;2001
year_e = 1952;2010
year_sum = year_e-year_s+1
var_sum  = dimsizes(var_name)

var = new((/3,11,360,720/),float)
var = 0.0
Moist = new((/3,7,360,720/),float)
Moist = 0.0
TotSM   = new((/3,360,720/),float)
TotSM   = 0.0


do case_num = 0,2
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
var!1 = "var"
var!2 = "lat"
var!3 = "lon"
var&lat = f->Qle&y
var&lon = f->Qle&x

Moist!0 = "case"
Moist!1 = "layer"
Moist!2 = "lat"
Moist!3 = "lon"
Moist&lat = f->Qle&y
Moist&lon = f->Qle&x


TotSM!0   = "case"
TotSM!1   = "lat"
TotSM!2   = "lon"
TotSM&lat = f->Qle&y
TotSM&lon = f->Qle&x

var@_FillValue   = f->Qle@_FillValue
Moist@_FillValue = f->Qle@_FillValue
TotSM@_FillValue   = f->Qle@_FillValue

do var_num = 0,10
   if (var_num .le. 4) then
     do month = 0,11
        var(case_num,var_num,:,:) = var(case_num,var_num,:,:) + f->$var_name(var_num)$(month,:,:)*3600*24*dom(month)
     end do
   else
      var(case_num,var_num,:,:) = var(case_num,var_num,:,:) + dim_avg_n_Wrap(f->$var_name(var_num)$,0)
   end if
end do ;var_num

Moist(case_num,0:5,:,:) = Moist(case_num,0:5,:,:) + dim_avg_n_Wrap(f->SoilMoist,0)
Moist(case_num,6,:,:)   = Moist(case_num,6,:,:)   + dim_avg_n_Wrap(f->GWMoist,0)


do layers = 0,5
;  m3/m2
TotSM(case_num, :, :) = TotSM(case_num, :, :)+ Moist(case_num,layers,:,:)*Soil_thickness(layers)
end do

end do ;year

end do ;case

var   = var/year_sum
var(:,9:10,:,:) = var(:,9:10,:,:) -273.16
Moist = Moist/year_sum
TotSM   = TotSM/year_sum

;_____________________________ PLOT ___________________________________
pic = "GW_Spitial_Map"
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
res@txString      = ""
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
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  pres@txString = f->$var_name(var_num)$@long_name+" (mm/year)"
  levels = (/100.,200.,300.,400.,500.,600.,700.,800.,900.,1000.,1100.,1200,1300,1400,1500,1600,1700/)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,var(0,var_num,:,:),res)
  res@tiMainString = "(b) NoGW"
  plots(1) = gsn_csm_contour_map(wks,var(1,var_num,:,:),res)
  res@tiMainString = "(c) GW-norecharge"
  plots(2) = gsn_csm_contour_map(wks,var(2,var_num,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)
end do

; "Qs","Qsb" ;"kg/m^2/s"
do var_num = 3,4
  pres@txString = f->$var_name(var_num)$@long_name+" (mm/year)"
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  levels = (/100.,200.,300.,400.,500.,600.,700.,800.,900.,1000.,1100.,1200,1300,1400.,1500.,1600./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,var(0,var_num,:,:),res)
  res@tiMainString = "(b) NoGW"
  plots(1) = gsn_csm_contour_map(wks,var(1,var_num,:,:),res)
  res@tiMainString = "(c) GW-norecharge"
  plots(2) = gsn_csm_contour_map(wks,var(2,var_num,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)
end do

; "Qle" "W/m^2"
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  pres@txString = f->$var_name(5)$@long_name+" ("+f->$var_name(5)$@units+")"
  levels = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,var(0,5,:,:),res)
  res@tiMainString = "(b) NoGW"
  plots(1) = gsn_csm_contour_map(wks,var(1,5,:,:),res)
  res@tiMainString = "(c) GW-norecharge"
  plots(2) = gsn_csm_contour_map(wks,var(2,5,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)

; "Qh" ;"W/m^2"
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  pres@txString = f->$var_name(6)$@long_name+" ("+f->$var_name(6)$@units+")"
  levels = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,110.,120.,130.,140.,150.,160,170.,180.,190.,200.,210.,220./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,var(0,6,:,:),res)
  res@tiMainString = "(b) NoGW"
  plots(1) = gsn_csm_contour_map(wks,var(1,6,:,:),res)
  res@tiMainString = "(c) GW-norecharge"
  plots(2) = gsn_csm_contour_map(wks,var(2,6,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)

; "GPP"
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  pres@txString = f->GPP@long_name+" ("+f->GPP@units+")"
  levels = (/1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,var(0,7,:,:),res)
  res@tiMainString = "(b) NoGW"
  plots(1) = gsn_csm_contour_map(wks,var(1,7,:,:),res)
  res@tiMainString = "(c) GW-norecharge"
  plots(2) = gsn_csm_contour_map(wks,var(2,7,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)

;"NPP"
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  pres@txString = f->NPP@long_name+" ("+f->NPP@units+")"
  levels = (/-5.,-4.,-3.,-2.,-1.,1.,2.,3.,4.,5./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,var(0,8,:,:),res)
  res@tiMainString = "(b) NoGW"
  plots(1) = gsn_csm_contour_map(wks,var(1,8,:,:),res)
  res@tiMainString = "(c) GW-norecharge"
  plots(2) = gsn_csm_contour_map(wks,var(2,8,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)

; "RadT"
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  pres@txString = f->RadT@long_name+" (~S~o~N~C)"
  levels = (/5.,10.,15.,20.,25.,30.,35.,40./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,var(0,9,:,:),res)
  res@tiMainString = "(b) NoGW"
  plots(1) = gsn_csm_contour_map(wks,var(1,9,:,:),res)
  res@tiMainString = "(c) GW-norecharge"
  plots(2) = gsn_csm_contour_map(wks,var(2,9,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)

;"VegT"
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  pres@txString = f->VegT@long_name+" (~S~o~N~C)"
  levels = (/5.,10.,15.,20.,25.,30.,35.,40./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,var(0,10,:,:),res)
  res@tiMainString = "(b) NoGW"
  plots(1) = gsn_csm_contour_map(wks,var(1,10,:,:),res)
  res@tiMainString = "(c) GW-norecharge"
  plots(2) = gsn_csm_contour_map(wks,var(2,10,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)

do layers = 0,6
  pres@txString = "The soil moisture of the "+(layers+1)+ " layer (m~S~3~N~/m~S~3~N~)"
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  levels = (/0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.60,0.65,0.70/)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,Moist(0,layers,:,:),res)
  res@tiMainString = "(b) NoGW"
  plots(1) = gsn_csm_contour_map(wks,Moist(1,layers,:,:),res)
  res@tiMainString = "(c) GW-norecharge"
  plots(2) = gsn_csm_contour_map(wks,Moist(2,layers,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)
end do

  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  levels = (/10.,20.,30.,40.,50.,60.,70.,80.,90.,100./)
  res@cnLevels = levels
  pres@txString = "Total Soil Moist (m~S~3~N~/m~S~2~N~)"
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,TotSM(0,:,:),res)
  res@tiMainString = "(b) NoGW"
  plots(1) = gsn_csm_contour_map(wks,TotSM(1,:,:),res)
  res@tiMainString = "(c) GW-norecharge"
  plots(2) = gsn_csm_contour_map(wks,TotSM(2,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)

end
