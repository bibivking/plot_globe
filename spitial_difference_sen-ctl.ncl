;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para_-precp_+temp/Outputs/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para_-precp_+temp/Outputs/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para_-precp_+temp/Outputs/gw_off"\
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

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers


;______________________________ Data _________________________________
year_s = 1971;2001
year_e = 1980;2010
year_sum = year_e-year_s+1
var_sum  = dimsizes(var_name)

var_diff = new((/3,11,360,720/),float)
var_diff = 0.0
Moist_diff = new((/3,7,360,720/),float)
Moist_diff = 0.0
TotSM_diff = new((/3,360,720/),float)
TotSM_diff = 0.0

do case_num = 0,2
print(case_name(case_num))

do year = year_s,year_e

    print("year = " +year)
    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if

filename1 = case_name(case_num)+"/cable_out_"+year+".nc"
filename2 = case_name(case_num+3)+"/cable_out_"+year+".nc"

f1 = addfile (filename1,"r")
f2 = addfile (filename2,"r")

var_diff!0 = "case"
var_diff!1 = "var"
var_diff!2 = "lat"
var_diff!3 = "lon"
var_diff&lat = f1->Qle&y
var_diff&lon = f1->Qle&x

Moist_diff!0 = "case"
Moist_diff!1 = "layer"
Moist_diff!2 = "lat"
Moist_diff!3 = "lon"
Moist_diff&lat = f1->Qle&y
Moist_diff&lon = f1->Qle&x


TotSM_diff!0   = "case"
TotSM_diff!1   = "lat"
TotSM_diff!2   = "lon"
TotSM_diff&lat = f1->Qle&y
TotSM_diff&lon = f1->Qle&x

var_diff@_FillValue   = f1->Qle@_FillValue
Moist_diff@_FillValue = f1->Qle@_FillValue
TotSM_diff@_FillValue = f1->Qle@_FillValue

do var_num = 0,10
   if (var_num .le. 4) then
     do month = 0,11
        var_diff(case_num,var_num,:,:) = var_diff(case_num,var_num,:,:) \
                                      + (f2->$var_name(var_num)$(month,:,:)-f1->$var_name(var_num)$(month,:,:)) \
                                      * 3600*24*dom(month)
     end do
   else
      var_diff(case_num,var_num,:,:) = var_diff(case_num,var_num,:,:) \
                                     + (dim_avg_n_Wrap(f2->$var_name(var_num)$,0) -dim_avg_n_Wrap(f1->$var_name(var_num)$,0))
   end if
end do ;var_num

Moist_diff(case_num,0:5,:,:) = Moist_diff(case_num,0:5,:,:) \
                               + ( dim_avg_n_Wrap(f2->SoilMoist,0) - dim_avg_n_Wrap(f1->SoilMoist,0))
Moist_diff(case_num,6,:,:)   = Moist_diff(case_num,6,:,:)  \
                               + ( dim_avg_n_Wrap(f2->GWMoist,0) - dim_avg_n_Wrap(f1->GWMoist,0))


do layers = 0,5
;  m3/m2
TotSM_diff(case_num, :, :) = TotSM_diff(case_num, :, :) \
                            + Moist_diff(case_num,layers,:,:)*Soil_thickness(layers)/sum(Soil_thickness)
end do

end do ;year

end do ;case

var_diff   = var_diff/year_sum
Moist_diff = Moist_diff/year_sum
TotSM_diff = TotSM_diff/year_sum

;_____________________________ PLOT ___________________________________
pic = "GW_Spitial_Map_difference_sen-ctl"
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
  pres@txString = f1->$var_name(var_num)$@long_name+" (mm/year)"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = (/-300.,-250.,-200.,-150.,-100.,-50.,-30.,-10.,10.,30.,50.,100.,150.,200.,250.,300./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,var_diff(0,var_num,:,:),res)
  res@tiMainString = "(b) GW-norecharge"
  plots(1) = gsn_csm_contour_map(wks,var_diff(1,var_num,:,:),res)
;  delete(levels)
;  delete(res@cnLevels)
;  levels = (/-80.,-70.,-60.,-50.,-40.,-30.,-20.,-10.,-5.,5.,10.,20.,30.,40.,50.,60.,70.,80./)
;  res@cnLevels = levels
  res@tiMainString = "(c) NoGW"
  plots(2) = gsn_csm_contour_map(wks,var_diff(2,var_num,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)
end do

; "Qs","Qsb" ;"kg/m^2/s"
do var_num = 3,4
  pres@txString = f1->$var_name(var_num)$@long_name+" (mm/year)"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = (/-800.,-600.,-400,-200.,-100.,-50.,50,100.,200.,400,600,800/)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,var_diff(0,var_num,:,:),res)
  res@tiMainString = "(b) GW-norecharge"
  plots(1) = gsn_csm_contour_map(wks,var_diff(1,var_num,:,:),res)
;  delete(levels)
;  delete(res@cnLevels)
;  levels = (/-80.,-70.,-60.,-50.,-40.,-30.,-20.,-10.,-5.,5.,10.,20.,30.,40.,50.,60.,70.,80./)
;  res@cnLevels = levels
  res@tiMainString = "(c) NoGW"
  plots(2) = gsn_csm_contour_map(wks,var_diff(2,var_num,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)
end do

; "Qle" "W/m^2"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  pres@txString = f1->$var_name(5)$@long_name+" ("+f1->$var_name(5)$@units+")"
  levels = (/-30.,-25.,-20.,-15.,-10.,-5.,5.,10.,15.,20.,25.,30./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,var_diff(0,5,:,:),res)
  res@tiMainString = "(b) GW-norecharge"
  plots(1) = gsn_csm_contour_map(wks,var_diff(1,5,:,:),res)
;  delete(levels)
;  delete(res@cnLevels)
;  levels = (/-10.,-8.,-6.,-4.,-2.,2.,4.,6.,8.,10./)
;  res@cnLevels = levels
  res@tiMainString = "(c) NoGW"
  plots(2) = gsn_csm_contour_map(wks,var_diff(2,5,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)

; "Qh" ;"W/m^2"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  pres@txString = f1->$var_name(6)$@long_name+" ("+f1->$var_name(6)$@units+")"
  levels = (/-40.,-35.,-30.,-25.,-20.,-15.,-10.,-5.,5.,10.,15.,20.,25.,30.,35.,40./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,var_diff(0,6,:,:),res)
  res@tiMainString = "(b) GW-norecharge"
  plots(1) = gsn_csm_contour_map(wks,var_diff(1,6,:,:),res)
;  delete(levels)
;  delete(res@cnLevels)
;  levels = (/-10.,-8.,-6.,-4.,-2.,2.,4.,6.,8.,10./)
;  res@cnLevels = levels
  res@tiMainString = "(c) NoGW"
  plots(2) = gsn_csm_contour_map(wks,var_diff(2,6,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)



; "GPP"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  pres@txString = f1->GPP@long_name+" ("+f1->GPP@units+")"
  levels = (/-2.5,-2.,-1.5,-1.,-0.5,-0.2,0.2,0.5,1.,1.5,2.,2.5/)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,var_diff(0,7,:,:),res)
  res@tiMainString = "(b) GW-norecharge"
  plots(1) = gsn_csm_contour_map(wks,var_diff(1,7,:,:),res)
;  delete(levels)
;  delete(res@cnLevels)
;  levels = (/-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5/)
;  res@cnLevels = levels
  res@tiMainString = "(c) NoGW"
  plots(2) = gsn_csm_contour_map(wks,var_diff(2,7,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)

;"NPP"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  pres@txString = f1->NPP@long_name+" ("+f1->NPP@units+")"
  levels = (/-2.5,-2.,-1.5,-1.,-0.5,-0.2,0.2,0.5,1.,1.5,2.,2.5/)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,var_diff(0,8,:,:),res)
  res@tiMainString = "(b) GW-norecharge"
  plots(1) = gsn_csm_contour_map(wks,var_diff(1,8,:,:),res)
;  delete(levels)
;  delete(res@cnLevels)
;  levels = (/-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5/)
;  res@cnLevels = levels
  res@tiMainString = "(c) NoGW"
  plots(2) = gsn_csm_contour_map(wks,var_diff(2,8,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)




; "RadT"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  pres@txString = f1->RadT@long_name+" (~S~o~N~C)"
  levels = (/-5.,-4,-3,-2,-1.,-0.5,-0.2,0.2,0.5,1.,2.,3.,4.,5./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,var_diff(0,9,:,:),res)
  res@tiMainString = "(b) GW-norecharge"
  plots(1) = gsn_csm_contour_map(wks,var_diff(1,9,:,:),res)
;  delete(levels)
;  delete(res@cnLevels)
;  levels = (/-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5/)
;  res@cnLevels = levels
  res@tiMainString = "(c) NoGW"
  plots(2) = gsn_csm_contour_map(wks,var_diff(2,9,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)


;"VegT"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  pres@txString = f1->VegT@long_name+" (~S~o~N~C)"
  levels = (/-3,-2.,-1.,-0.5,-0.2,0.2,0.5,1.,2.,3./)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,var_diff(0,10,:,:),res)
  res@tiMainString = "(b) GW-norecharge"
  plots(1) = gsn_csm_contour_map(wks,var_diff(1,10,:,:),res)
;  delete(levels)
;  delete(res@cnLevels)
;  levels = (/-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5/)
;  res@cnLevels = levels
  res@tiMainString = "(c) NoGW"
  plots(2) = gsn_csm_contour_map(wks,var_diff(2,10,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)

do layers = 0,6
  pres@txString = "The soil moisture of the "+(layers+1)+ " layer (m~S~3~N~/m~S~3~N~)"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = (/-0.25,-0.2,-0.15,-0.1,-0.05,-0.02,0.02,0.05,0.1,0.15,0.2,0.25/)
  res@cnLevels = levels
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,Moist_diff(0,layers,:,:),res)
  res@tiMainString = "(b) GW-norecharge"
  plots(1) = gsn_csm_contour_map(wks,Moist_diff(1,layers,:,:),res)
;  delete(levels)
;  delete(res@cnLevels)
;  levels = (/-0.1,-0.08,-0.06,-0.04,-0.02,0.02,0.04,0.06,0.08,0.1/)
;  res@cnLevels = levels
  res@tiMainString = "(c) NoGW"
  plots(2) = gsn_csm_contour_map(wks,Moist_diff(2,layers,:,:),res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)
end do

  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = (/-0.25,-0.2,-0.15,-0.1,-0.05,-0.02,0.02,0.05,0.1,0.15,0.2,0.25/)
  res@cnLevels = levels
  pres@txString = "Total Soil Moist (m~S~3~N~/m~S~3~N~)"
  res@tiMainString = "(a) GW"
  plots(0) = gsn_csm_contour_map(wks,TotSM_diff(0,:,:) ,res)
  res@tiMainString = "(b) GW-norecharge"
  plots(1) = gsn_csm_contour_map(wks,TotSM_diff(1,:,:),res)
;  delete(levels)
;  delete(res@cnLevels)
;  levels = (/-0.1,-0.08,-0.06,-0.04,-0.02,0.02,0.04,0.06,0.08,0.1/)
;  res@cnLevels = levels
  res@tiMainString = "(c) NoGW"
  plots(2) = gsn_csm_contour_map(wks,TotSM_diff(2,:,:) ,res)
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  delete(levels)
  delete(res@cnLevels)

end
