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

plot_name = (/\
               "(a) NoGW - GW",\
               "(b) NoAQ - GW",\
               "(c) NoGW_redstrb - NoGW",\
               "(d) GW_Hvrd - GW"\
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

var   = new((/case_sum,11,360,720/),float)
var   = 0.0
Moist = new((/case_sum,7,360,720/),float)
Moist = 0.0
TotSM = new((/case_sum,360,720/),float)
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
    var!1 = "var"
    var!2 = "lat"
    var!3 = "lon"
    var&lat = lat
    var&lon = lon

    Moist!0 = "case"
    Moist!1 = "layer"
    Moist!2 = "lat"
    Moist!3 = "lon"
    Moist&lat = lat
    Moist&lon = lon

    TotSM!0   = "case"
    TotSM!1   = "lat"
    TotSM!2   = "lon"
    TotSM&lat = lat
    TotSM&lon = lon 
    
    var@_FillValue   = f->Qle@_FillValue
    Moist@_FillValue = f->Qle@_FillValue
    TotSM@_FillValue = f->Qle@_FillValue

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
        TotSM(case_num, :, :) = TotSM(case_num, :, :)+ Moist(case_num,layers,:,:)*Soil_thickness(layers)/sum(Soil_thickness)
    end do

end do ;year

end do ;case

var   = var/year_sum
var(:,9:10,:,:) = var(:,9:10,:,:) -273.16
Moist = Moist/year_sum
TotSM = TotSM/year_sum

var_diff = new((/4,11,360,720/),float)
var_diff = 0.0
Moist_diff = new((/4,7,360,720/),float)
Moist_diff = 0.0
TotSM_diff = new((/4,360,720/),float)
TotSM_diff = 0.0

var_diff!0 = "case"
var_diff!1 = "var"
var_diff!2 = "lat"
var_diff!3 = "lon"
var_diff&lat = lat
var_diff&lon = lon

Moist_diff!0 = "case"
Moist_diff!1 = "layer"
Moist_diff!2 = "lat"
Moist_diff!3 = "lon"
Moist_diff&lat = lat
Moist_diff&lon = lon

TotSM_diff!0   = "case"
TotSM_diff!1   = "lat"
TotSM_diff!2   = "lon"
TotSM_diff&lat = lat
TotSM_diff&lon = lon 
    
var_diff@_FillValue   = f->Qle@_FillValue
Moist_diff@_FillValue = f->Qle@_FillValue
TotSM_diff@_FillValue = f->Qle@_FillValue
               
var_diff(0,:,:,:) = var(0,:,:,:) - var(2,:,:,:) ; 2 on-std - 0 off-std
var_diff(1,:,:,:) = var(3,:,:,:) - var(2,:,:,:) ; 2 on-std - 3 on_norchrg-std
var_diff(2,:,:,:) = var(1,:,:,:) - var(0,:,:,:) ; 1 off-redstrb-std - 0 off-std
var_diff(3,:,:,:) = var(4,:,:,:) - var(2,:,:,:) ; 4 on-Hvrd - 2 on-std

Moist_diff(0,:,:,:) = Moist(0,:,:,:) - Moist(2,:,:,:)
Moist_diff(1,:,:,:) = Moist(3,:,:,:) - Moist(2,:,:,:)
Moist_diff(2,:,:,:) = Moist(1,:,:,:) - Moist(0,:,:,:)
Moist_diff(3,:,:,:) = Moist(4,:,:,:) - Moist(2,:,:,:)

TotSM_diff(0,:,:) = TotSM(0,:,:) - TotSM(2,:,:)
TotSM_diff(1,:,:) = TotSM(3,:,:) - TotSM(2,:,:)
TotSM_diff(2,:,:) = TotSM(1,:,:) - TotSM(0,:,:)
TotSM_diff(3,:,:) = TotSM(4,:,:) - TotSM(2,:,:)
;_____________________________ PLOT ___________________________________
pic = "spitial_difference_Qle_"+scale
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

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
res@gsnLeftString      = ""
res@gsnRightString     = ""
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

plots = new(4,graphic)

; "Qle" "W/m^2"
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  pres@txString = f->$var_name(5)$@long_name+" (W/m~S~2~N~)"
  levels = fspan(-20,20,21)
  res@cnLevels = levels
  pres@cnLevels = levels
  do plot_num = 0,3;1
     res@tiMainString = plot_name(plot_num)
     plots(plot_num) = gsn_csm_contour_map(wks,var_diff(plot_num,5,:,:),res)
  end do
  gsn_panel(wks,(/plots/),(/2,2/),pres)
  ;do plot_num = 2,3
  ;   res@tiMainString = plot_name(plot_num)
  ;   plots(plot_num-2) = gsn_csm_contour_map(wks,var_diff(plot_num,5,:,:),res)
  ;end do
  ;gsn_panel(wks,(/plots/),(/2,1/),pres)
  delete(levels)
  delete(res@cnLevels)
  delete(pres@cnLevels)

end
