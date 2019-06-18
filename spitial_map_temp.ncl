;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-standard/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs-standard/gw_off"\
              /)

case_name_brf = (/\
              "GW",\
              "CLASS",\
              "CERES"\
               /)

var_name= (/\
           "RadT","VegT",\   ; "K"
           "SnowT","BaresoilT",\
           "SoilTemp",\
           "SoilMoist",\     ;"m^3/m^3"  Average layer soil moisture
           "SoilMoistIce",\
           "GWMoist",\       ;"mm3/mm3" Aquifer mositure content, Groundwater_off doesn't have
           "WatTable"\       ;"m"
           /)
           ;"LAI" ; the same in Groundwater_on and Groundwater_off

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers

;______________________________ Data _________________________________
year_s = 2003
year_e = 2009
year_sum = year_e-year_s+1
var_sum  = dimsizes(var_name)
case_sum = dimsizes(case_name)

scale = "Global" ;"Global" ; "50S-N" ; "SE-AU" ; "AU"

CABLE = new((/4,360,720/),float)
CABLE = 0.0
SoilTemp = new((/6,360,720/),float)
SoilTemp = 0.0

CABLE_dif = new((/4,360,720/),float)
CABLE_dif = 0.0
SoilTemp_dif = new((/6,360,720/),float)
SoilTemp_dif = 0.0


do year = year_s,year_e

    print("year = " +year)
    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if

    filename = case_name(0)+"/cable_out_"+year+".nc"
    f = addfile (filename,"r")

    filename1 = case_name(1)+"/cable_out_"+year+".nc"
    f1 = addfile (filename1,"r")

    lat = f->Qle&y
    lon = f->Qle&x

    CABLE!0 = "var"
    CABLE!1 = "lat"
    CABLE!2 = "lon"
    CABLE&lat = lat
    CABLE&lon = lon

    CABLE@_FillValue   = f->Qle@_FillValue

    SoilTemp!0 = "layer"
    SoilTemp!1 = "lat"
    SoilTemp!2 = "lon"
    SoilTemp&lat = lat
    SoilTemp&lon = lon

    SoilTemp@_FillValue   = f->Qle@_FillValue
    
    CABLE_dif!0 = "var"
    CABLE_dif!1 = "lat"
    CABLE_dif!2 = "lon"
    CABLE_dif&lat = lat
    CABLE_dif&lon = lon

    CABLE_dif@_FillValue   = f->Qle@_FillValue

    SoilTemp_dif!0 = "layer"
    SoilTemp_dif!1 = "lat"
    SoilTemp_dif!2 = "lon"
    SoilTemp_dif&lat = lat
    SoilTemp_dif&lon = lon

    SoilTemp_dif@_FillValue   = f->Qle@_FillValue
    
    
    do var_num = 0,3
       CABLE(var_num,:,:) = CABLE(var_num,:,:) + dim_avg_n_Wrap(f->$var_name(var_num)$,0)
       CABLE_dif(var_num,:,:) = CABLE_dif(var_num,:,:) + \
                 (dim_avg_n_Wrap(f->$var_name(var_num)$,0) - dim_avg_n_Wrap(f1->$var_name(var_num)$,0))
    end do ;var_num
    SoilTemp = SoilTemp + dim_avg_n_Wrap(f->SoilTemp,0)
    SoilTemp_dif = SoilTemp_dif + (dim_avg_n_Wrap(f->SoilTemp,0) - dim_avg_n_Wrap(f1->SoilTemp,0))
end do ;year


CABLE = CABLE/year_sum
CABLE_dif = CABLE_dif/year_sum
CABLE = CABLE -273.16
SoilTemp = SoilTemp/year_sum
SoilTemp_dif = SoilTemp_dif/year_sum
SoilTemp = SoilTemp - 273.16

;_____________________________ PLOT ___________________________________
pic = "spitial_map_temp_on_"+scale
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
pres@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
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
  res@cnFillPalette = "WhiteBlueGreenYellowRed";"ViBlGrWhYeOrRe"
  levels = fspan(0,50,21)
  res@cnLevels = levels
  res@tiMainString = "(a) CABLE RadT (~S~o~N~C)"
  plots(0) = gsn_csm_contour_map(wks,CABLE(0,:,:),res)
  res@tiMainString = "(b) CABLE VegT (~S~o~N~C)"
  plots(1) = gsn_csm_contour_map(wks,CABLE(1,:,:),res)
  res@tiMainString = "(c) CABLE SnowT (~S~o~N~C)"
  plots(2) = gsn_csm_contour_map(wks,CABLE(2,:,:),res)
  res@tiMainString = "(d) CABLE BaresoilT (~S~o~N~C)"
  plots(3) = gsn_csm_contour_map(wks,CABLE(3,:,:),res)

  gsn_panel(wks,(/plots/),(/2,2/),pres)   
    
  delete(levels)
  delete(res@cnLevels)
  delete(plots)

  plots = new(6,graphic)
  res@cnFillPalette = "WhiteBlueGreenYellowRed";"ViBlGrWhYeOrRe"
  levels = fspan(0,50,21)
  res@cnLevels = levels
  res@tiMainString = "(a) CABLE Soil Temp Layer1 (~S~o~N~C)"
  plots(0) = gsn_csm_contour_map(wks,SoilTemp(0,:,:),res)
  res@tiMainString = "(b) CABLE Soil Temp Layer2 (~S~o~N~C)"
  plots(1) = gsn_csm_contour_map(wks,SoilTemp(1,:,:),res)
  res@tiMainString = "(c) CABLE Soil Temp Layer3 (~S~o~N~C)"
  plots(2) = gsn_csm_contour_map(wks,SoilTemp(2,:,:),res)
  res@tiMainString = "(d) CABLE Soil Temp Layer4 (~S~o~N~C)"
  plots(3) = gsn_csm_contour_map(wks,SoilTemp(3,:,:),res)
  res@tiMainString = "(e) CABLE Soil Temp Layer5 (~S~o~N~C)"
  plots(4) = gsn_csm_contour_map(wks,SoilTemp(4,:,:),res)
  res@tiMainString = "(f) CABLE Soil Temp Layer6 (~S~o~N~C)"
  plots(5) = gsn_csm_contour_map(wks,SoilTemp(5,:,:),res)
  
  delete(levels)
  delete(res@cnLevels)
  
  gsn_panel(wks,(/plots/),(/3,2/),pres)  
  
  delete(plots)
  
print(max(CABLE_dif))
print(min(CABLE_dif))
  plots = new(4,graphic)
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = fspan(-10,10,21)
  res@cnLevels = levels
  res@tiMainString = "(a) CABLE RadT (~S~o~N~C)"
  plots(0) = gsn_csm_contour_map(wks,CABLE_dif(0,:,:),res)
  res@tiMainString = "(b) CABLE VegT (~S~o~N~C)"
  plots(1) = gsn_csm_contour_map(wks,CABLE_dif(1,:,:),res)
  res@tiMainString = "(c) CABLE SnowT (~S~o~N~C)"
  plots(2) = gsn_csm_contour_map(wks,CABLE_dif(2,:,:),res)
  res@tiMainString = "(d) CABLE BaresoilT (~S~o~N~C)"
  plots(3) = gsn_csm_contour_map(wks,CABLE_dif(3,:,:),res)

  gsn_panel(wks,(/plots/),(/2,2/),pres)   
    
  delete(levels)
  delete(res@cnLevels)
  delete(plots)
  
print(min(SoilTemp_dif))
print(max(SoilTemp_dif))
  plots = new(6,graphic)
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = fspan(-10,10,21)
  res@cnLevels = levels
  res@tiMainString = "(a) CABLE Soil Temp Layer1 (~S~o~N~C)"
  plots(0) = gsn_csm_contour_map(wks,SoilTemp_dif(0,:,:),res)
  res@tiMainString = "(b) CABLE Soil Temp Layer2 (~S~o~N~C)"
  plots(1) = gsn_csm_contour_map(wks,SoilTemp_dif(1,:,:),res)
  res@tiMainString = "(c) CABLE Soil Temp Layer3 (~S~o~N~C)"
  plots(2) = gsn_csm_contour_map(wks,SoilTemp_dif(2,:,:),res)
  res@tiMainString = "(d) CABLE Soil Temp Layer4 (~S~o~N~C)"
  plots(3) = gsn_csm_contour_map(wks,SoilTemp_dif(3,:,:),res)
  res@tiMainString = "(e) CABLE Soil Temp Layer5 (~S~o~N~C)"
  plots(4) = gsn_csm_contour_map(wks,SoilTemp_dif(4,:,:),res)
  res@tiMainString = "(f) CABLE Soil Temp Layer6 (~S~o~N~C)"
  plots(5) = gsn_csm_contour_map(wks,SoilTemp_dif(5,:,:),res)
  
  delete(levels)
  delete(res@cnLevels)
  
  gsn_panel(wks,(/plots/),(/3,2/),pres)  
  
  delete(plots)
  

end
