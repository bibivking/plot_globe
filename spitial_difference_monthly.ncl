;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs-standard/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para_redistrb-on/Outputs-standard/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs-Haverd2013/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-standard/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs-standard/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-Haverd2013/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs-Haverd2013/gw_on"\
             /)

case_name_brf = (/\
              "NoGW",\
              "NoGW_redstrb",\
              "NoGW_Hvrd",\
              "GW",\
              "NoAQ",\
              "GW_Hvrd",\
              "NoAQ_Hvrd"\
               /)

plot_name = (/\
               "(a) on-std - off-std",\
               "(b) on-std - on_norchrg-std",\
               "(c) off-redstrb-std - off-std",\
               "(d) off-Hvrd - off-std",\
               "(e) on-Hvrd - on-std",\
               "(f) on_norchrg-Hvrd - on_norchrg-std"\
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
year_s = 2001
year_e = 2010
year_sum = year_e-year_s+1
var_sum  = dimsizes(var_name)
case_sum = dimsizes(case_name)

scale = "Global" ;"Global" ; "50S-N" ; "SE-AU" ; "AU"

var   = new((/case_sum,12,11,360,720/),float)
var   = 0.0
Moist = new((/case_sum,12,7,360,720/),float)
Moist = 0.0
TotSM = new((/case_sum,12,360,720/),float)
TotSM = 0.0

do case_num = 0,case_sum-1

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

    lat = f->Qle&y
    lon = f->Qle&x

    var!0 = "case"
    var!1 = "month"
    var!2 = "var"
    var!3 = "lat"
    var!4 = "lon"
    var&lat = lat
    var&lon = lon

    Moist!0 = "case"
    Moist!1 = "month"
    Moist!2 = "layer"
    Moist!3 = "lat"
    Moist!4 = "lon"
    Moist&lat = lat
    Moist&lon = lon

    TotSM!0   = "case"
    TotSM!1   = "month"
    TotSM!2   = "lat"
    TotSM!3   = "lon"
    TotSM&lat = lat
    TotSM&lon = lon

    var@_FillValue   = f->Qle@_FillValue
    Moist@_FillValue = f->Qle@_FillValue
    TotSM@_FillValue = f->Qle@_FillValue

    do var_num = 0,10

      if (var_num .le. 4) then
         do month = 0,11
            var(case_num,month,var_num,:,:) = var(case_num,month,var_num,:,:) + f->$var_name(var_num)$(month,:,:)*3600*24*dom(month)
         end do
      else
         var(case_num,:,var_num,:,:) = var(case_num,:,var_num,:,:) + f->$var_name(var_num)$
      end if
    end do ;var_num

    Moist(case_num,:,0:5,:,:) = Moist(case_num,:,0:5,:,:) + f->SoilMoist
    Moist(case_num,:,6,:,:)   = Moist(case_num,:,6,:,:)   + f->GWMoist

    do layers = 0,5
    ;  m3/m2
        TotSM(case_num, :, :, :) = TotSM(case_num,:, :, :)+ Moist(case_num,:,layers,:,:)*Soil_thickness(layers)/sum(Soil_thickness)
    end do

end do ;year

end do ;case

var   = var/year_sum
var(:,:,9:10,:,:) = var(:,:,9:10,:,:) -273.16
Moist = Moist/year_sum
TotSM = TotSM/year_sum

var_diff = new((/6,12,11,360,720/),float)
var_diff = 0.0

var_diff!0 = "case"
var_diff!1 = "month"
var_diff!2 = "var"
var_diff!3 = "lat"
var_diff!4 = "lon"
var_diff&lat = lat
var_diff&lon = lon
var_diff@_FillValue   = f->Qle@_FillValue

var_diff(0,:,:,:,:) = var(3,:,:,:,:) - var(0,:,:,:,:) ; 3 on-std - 0 off-std
var_diff(1,:,:,:,:) = var(3,:,:,:,:) - var(4,:,:,:,:) ; 3 on-std - 4 on_norchrg-std
var_diff(2,:,:,:,:) = var(1,:,:,:,:) - var(0,:,:,:,:) ; 1 off-redstrb-std - 0 off-std
var_diff(3,:,:,:,:) = var(2,:,:,:,:) - var(0,:,:,:,:) ; 2 off-Hvrd - 0 off-std
var_diff(4,:,:,:,:) = var(5,:,:,:,:) - var(3,:,:,:,:) ; 5 on-Hvrd - 3 on-std
var_diff(5,:,:,:,:) = var(6,:,:,:,:) - var(4,:,:,:,:) ; 6 on_norchrg-Hvrd - 4 on_norchrg-std
delete(var)

Moist_diff = new((/6,12,7,360,720/),float)
Moist_diff = 0.0
Moist_diff!0 = "case"
Moist_diff!1 = "month"
Moist_diff!2 = "layer"
Moist_diff!3 = "lat"
Moist_diff!4 = "lon"
Moist_diff&lat = lat
Moist_diff&lon = lon
Moist_diff@_FillValue = f->Qle@_FillValue

Moist_diff(0,:,:,:,:) = Moist(3,:,:,:,:) - Moist(0,:,:,:,:)
Moist_diff(1,:,:,:,:) = Moist(3,:,:,:,:) - Moist(4,:,:,:,:)
Moist_diff(2,:,:,:,:) = Moist(1,:,:,:,:) - Moist(0,:,:,:,:)
Moist_diff(3,:,:,:,:) = Moist(2,:,:,:,:) - Moist(0,:,:,:,:)
Moist_diff(4,:,:,:,:) = Moist(5,:,:,:,:) - Moist(3,:,:,:,:)
Moist_diff(5,:,:,:,:) = Moist(6,:,:,:,:) - Moist(4,:,:,:,:)
delete(Moist)

TotSM_diff = new((/6,12,360,720/),float)
TotSM_diff = 0.0
TotSM_diff!0   = "case"
TotSM_diff!1   = "month"
TotSM_diff!2   = "lat"
TotSM_diff!3   = "lon"
TotSM_diff&lat = lat
TotSM_diff&lon = lon
TotSM_diff@_FillValue = f->Qle@_FillValue

TotSM_diff(0,:,:,:) = TotSM(3,:,:,:) - TotSM(0,:,:,:)
TotSM_diff(1,:,:,:) = TotSM(3,:,:,:) - TotSM(4,:,:,:)
TotSM_diff(2,:,:,:) = TotSM(1,:,:,:) - TotSM(0,:,:,:)
TotSM_diff(3,:,:,:) = TotSM(2,:,:,:) - TotSM(0,:,:,:)
TotSM_diff(4,:,:,:) = TotSM(5,:,:,:) - TotSM(3,:,:,:)
TotSM_diff(5,:,:,:) = TotSM(6,:,:,:) - TotSM(4,:,:,:)
delete(TotSM)

;_____________________________ PLOT ___________________________________
res                    = True
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                            ;��ɫ
res@tmBorderThicknessF = 3.0

res@gsnDraw            = False  ; Don't draw plots
res@gsnFrame           = False  ; ����ҳ
res@lbLabelBarOn       = False

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
res@tmYLMode  = "Explicit"
res@tmXBMode  = "Explicit"

if scale .eq. "Global" then
  A = -60.
  B = 90.
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
res@trYMinF   = A
res@trYMaxF   = B
res@mpMinLatF = A
res@mpMaxLatF = B
res@trXMinF   = C
res@trXMaxF   = D
res@mpMinLonF = C
res@mpMaxLonF = D

if scale .eq. "Global" then
   res@mpCenterLonF       = 160
   res@tmYLValues=(/60,30,0,-30,-60/)
   res@tmYLLabels=(/"60~S~o~N~N","30~S~o~N~N","0~S~o~N~","30~S~o~N~S","60~S~o~N~S"/)
   res@tmXBValues=(/240,300,0,60,120,180/)
   res@tmXBLabels=(/"120~S~o~N~W","60~S~o~N~W","0~S~o~N~","60~S~o~N~E","120~S~o~N~E","180~S~o~N~"/)
else if scale .eq. "50S-N" then
  res@mpCenterLonF       = 160
  res@tmYLValues=(/50,25,0,-25,-50/)
  res@tmYLLabels=(/"50~S~o~N~N","25~S~o~N~N","0~S~o~N~","25~S~o~N~S","50~S~o~N~S"/)
  res@tmXBValues=(/240,300,0,60,120,180/)
  res@tmXBLabels=(/"120~S~o~N~W","60~S~o~N~W","0~S~o~N~","60~S~o~N~E","120~S~o~N~E","180~S~o~N~"/)
else if scale .eq. "SE-AU" then
   res@gsnAddCyclic = False
   res@tmYLValues=(/-30,-40/)
   res@tmYLLabels=(/"30~S~o~N~S","40~S~o~N~S"/)
   res@tmXBValues=(/140,150/)
   res@tmXBLabels=(/"140~S~o~N~E","150~S~o~N~E"/)
else if scale .eq. "AU" then
  res@gsnAddCyclic       = False
  res@tmYLValues=(/-10,-20,-30,-40/)
  res@tmYLLabels=(/"10~S~o~N~S","20~S~o~N~S","30~S~o~N~S","40~S~o~N~S"/)
  res@tmXBValues=(/120,130,140,150/)
  res@tmXBLabels=(/"120~S~o~N~E","130~S~o~N~E","140~S~o~N~E","150~S~o~N~E"/)
end if
end if
end if
end if

;**************** ͼ�������� *****************
pres                    = True                                      ; ���� panel resources.
pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = True                    ; ���� panel labelbar.
pres@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
pres@lbBoxLinesOn       = True                       ;�ر�lbar box ����
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
;pres@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
pres@lbTitleFontHeightF = 0.01
pres@lbLabelFontHeightF = 0.01
;pres@lbAutoManage       = False
;pres@lbTitleDirection   = "Across"             ; title direction
;pres@lbTitlePosition    = "Bottom";Right"
;pres@lbLabelBarOn       = True

plots = new(6,graphic)

do month = 0,11

pic = "spitial_difference_monthly-"+(month+1)+"_"+scale
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

; "Evap","TVeg","ESoil",\ ;"kg/m^2/s"
do var_num = 0,2

  pres@txString = f->$var_name(var_num)$@long_name+" (mm/mon)"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = (/-50,-40,-30,-20.,-10.,0.,10,20.,30.,40,50/)
  res@cnLevels = levels
  pres@cnLevels = levels
  do plot_num = 0,5
     res@tiMainString = plot_name(plot_num)
     plots(plot_num) = gsn_csm_contour_map(wks,var_diff(plot_num,month,var_num,:,:),res)
  end do
  gsn_panel(wks,(/plots/),(/3,2/),pres)
  delete(levels)
  delete(res@cnLevels)
  delete(pres@cnLevels)
end do

; "Qs","Qsb" ;"kg/m^2/s"
;do var_num = 3,4
  pres@txString = "Runoff (mm/mon)"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = (/-50,-40,-30,-20.,-10.,0.,10,20.,30.,40,50/)
  res@cnLevels = levels
  pres@cnLevels = levels

  do plot_num = 0,5
     tmp = var_diff(plot_num,month,3,:,:)
     tmp = tmp + var_diff(plot_num,month,4,:,:)
     res@tiMainString = plot_name(plot_num)
     plots(plot_num) = gsn_csm_contour_map(wks,tmp,res)
     delete(tmp)
  end do
  gsn_panel(wks,(/plots/),(/3,2/),pres)
  delete(levels)
  delete(res@cnLevels)
  delete(pres@cnLevels)
;end do

; "Qle" "W/m^2"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  pres@txString = f->$var_name(5)$@long_name+" (W/m~S~2~N~)"
  levels = (/-30.,-25.,-20.,-15.,-10.,-5.,5.,10.,15.,20.,25.,30./)
  ;  levels = (/-10.,-8.,-6.,-4.,-2.,2.,4.,6.,8.,10./)
  res@cnLevels = levels
  pres@cnLevels = levels
  do plot_num = 0,5
     res@tiMainString = plot_name(plot_num)
     plots(plot_num) = gsn_csm_contour_map(wks,var_diff(plot_num,month,5,:,:),res)
  end do
  gsn_panel(wks,(/plots/),(/3,2/),pres)
  delete(levels)
  delete(res@cnLevels)
  delete(pres@cnLevels)

; "Qh" ;"W/m^2"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  pres@txString = f->$var_name(6)$@long_name+" (W/m~S~2~N~)"
  levels = (/-30.,-25.,-20.,-15.,-10.,-5.,5.,10.,15.,20.,25.,30./)
  ;  levels = (/-10.,-8.,-6.,-4.,-2.,2.,4.,6.,8.,10./)
  res@cnLevels = levels
  pres@cnLevels = levels
  do plot_num = 0,5
     res@tiMainString = plot_name(plot_num)
     plots(plot_num) = gsn_csm_contour_map(wks,var_diff(plot_num,month,6,:,:),res)
  end do
  gsn_panel(wks,(/plots/),(/3,2/),pres)
  delete(levels)
  delete(res@cnLevels)
  delete(pres@cnLevels)

; "GPP"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  pres@txString = f->GPP@long_name+" (umol/m~S~2~N~/s)"
  levels = (/-2.5,-2.,-1.5,-1.,-0.5,-0.2,0.2,0.5,1.,1.5,2.,2.5/)
;  levels = (/-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5/)
  res@cnLevels = levels
  pres@cnLevels= levels
  do plot_num = 0,5
     res@tiMainString = plot_name(plot_num)
     plots(plot_num) = gsn_csm_contour_map(wks,var_diff(plot_num,month,7,:,:),res)
  end do
  gsn_panel(wks,(/plots/),(/3,2/),pres)
  delete(levels)
  delete(res@cnLevels)
  delete(pres@cnLevels)

; "RadT"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  pres@txString = f->RadT@long_name+" (~S~o~N~C)"
  levels = (/-5.,-4,-3,-2,-1.,-0.5,-0.2,0.2,0.5,1.,2.,3.,4.,5./)
  ;levels = (/-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5/)
  res@cnLevels = levels
  pres@cnLevels = levels
  do plot_num = 0,5
     res@tiMainString = plot_name(plot_num)
     plots(plot_num) = gsn_csm_contour_map(wks,var_diff(plot_num,month,9,:,:),res)
  end do
  gsn_panel(wks,(/plots/),(/3,2/),pres)
  delete(levels)
  delete(res@cnLevels)
  delete(pres@cnLevels)

;"VegT"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  pres@txString = f->VegT@long_name+" (~S~o~N~C)"
  levels = (/-3,-2.,-1.,-0.5,-0.2,0.2,0.5,1.,2.,3./)
  ;levels = (/-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5/)
  res@cnLevels = levels
  pres@cnLevels = levels
  do plot_num = 0,5
     res@tiMainString = plot_name(plot_num)
     plots(plot_num) = gsn_csm_contour_map(wks,var_diff(plot_num,month,10,:,:),res)
  end do
  gsn_panel(wks,(/plots/),(/3,2/),pres)
  delete(levels)
  delete(res@cnLevels)
  delete(pres@cnLevels)

do layers = 0,6
  pres@txString = "The soil moisture of the "+(layers+1)+ " layer (m~S~3~N~/m~S~3~N~)"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = (/-0.25,-0.2,-0.15,-0.1,-0.05,-0.02,0.02,0.05,0.1,0.15,0.2,0.25/)
  ;levels = (/-0.1,-0.08,-0.06,-0.04,-0.02,0.02,0.04,0.06,0.08,0.1/)
  res@cnLevels = levels
  pres@cnLevels = levels
  do plot_num = 0,5
     res@tiMainString = plot_name(plot_num)
     plots(plot_num) = gsn_csm_contour_map(wks,Moist_diff(plot_num,month,layers,:,:),res)
  end do
  gsn_panel(wks,(/plots/),(/3,2/),pres)
  delete(levels)
  delete(res@cnLevels)
  delete(pres@cnLevels)
end do

  pres@txString = "Soil moisture in 6 layers (m~S~3~N~/m~S~3~N~)"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = (/-0.25,-0.2,-0.15,-0.1,-0.05,-0.02,0.02,0.05,0.1,0.15,0.2,0.25/)
  res@cnLevels = levels
  pres@cnLevels = levels
  do plot_num = 0,5
     res@tiMainString = plot_name(plot_num)
     plots(plot_num) = gsn_csm_contour_map(wks,TotSM_diff(plot_num,month,:,:),res)
  end do
  gsn_panel(wks,(/plots/),(/3,2/),pres)
  delete(levels)
  delete(res@cnLevels)
  delete(pres@cnLevels)
end do

end
