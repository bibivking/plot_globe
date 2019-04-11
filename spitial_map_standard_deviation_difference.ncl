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
               "gw_off-standard",\
               "gw_off-redistrb_on-standard",\
               "gw_off-Haverd2013",\
               "gw_on-standard",\
               "gw_on_norecharge-standard",\
               "gw_on-Haverd2013",\
               "gw_on_norecharge-Haverd2013"\
                /)

var_name= (/\
             "TVeg",  \
             "ESoil", \ ;"kg/m^2/s"
             "Evap",  \
             "Rainf"  \ ; "kg/m^2/s" ; "Rainfall+snowfall"
         	 /)

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers and 1 aquifer

; ____________________________ Options _______________________________
deseason = True
scale    = "World";"AU" ; World
;______________________________ Data _________________________________
year_s = 2001
year_e = 2010
year_sum = year_e-year_s+1
var_sum  = dimsizes(var_name)
total    = year_sum*12

pic1 = "check"
wks1 = gsn_open_wks("pdf",pic1)
gsn_define_colormap(wks1,"ViBlGrWhYeOrRe")
res1 = True
res1@cnFillMode         = "RasterFill"            ; Raster Mode
res1@cnFillOn           = True                            ;��ɫ
res1@tmBorderThicknessF = 3.0
res1@cnLinesOn          = False                       ; �رյ�ֵ������

;______________________________ INPUT __________________________________
std   = new((/7,360,720/),float)
std   = 0.0

do case_num = 0,6
  print(case_name(case_num))

  var   = new((/total,360,720/),float)
  var   = 0.0
  counter = 0

  do year = year_s,year_e
    print("year = " +year)
    filename = case_name(case_num)+"/cable_out_"+year+".nc"
    f = addfile (filename,"r")

    do month = 0,11
      do layer = 0,5
        var(counter,:,:) = var(counter,:,:) + (f->SoilMoist(month,layer,:,:)*Soil_thickness(layer))*1000.
      end do
      if case_num .ge. 3 then
        var(counter,:,:) = var(counter,:,:) + (f->GWMoist(month,:,:)*f->GWdz)*1000.
      end if
      counter = counter + 1
    end do ;month

  end do ;year

  var!0 = "time"
  var!1 = "lat"
  var!2 = "lon"
  var&lat = f->Qle&y
  var&lon = f->Qle&x
  var@_FillValue   = f->Qle@_FillValue

  plot = gsn_csm_contour_map(wks1,dim_avg_n_Wrap(var(:,:,:),0),res1)
  plot = gsn_csm_contour(wks1,dim_avg_n_Wrap(var(:,:,:),0)*1.0,res1)

  ; __________________________ Deseasonalize _____________________________
  ; rmMonAnnCycTLL : Removes the annual cycle from "monthly" data.
  if deseason .eq. True then
    var = rmMonAnnCycTLL(var)
    plot = gsn_csm_contour_map(wks1,dim_avg_n_Wrap(var,0),res1)
    plot = gsn_csm_contour(wks1,dim_avg_n_Wrap(var,0)*1.0,res1)
  end if

  ; _____________________  CALCULATE stddev ________________________
  print(case_name_brf(case_num))
  std(case_num,:,:) = dim_stddev_n_Wrap(var,0)
  delete(var)
  std!0 = "case"
  std!1 = "lat"
  std!2 = "lon"
  std&lat = f->Qle&y
  std&lon = f->Qle&x
  std@_FillValue   = f->Qle@_FillValue


end do ;case

;_____________________________ PLOT ___________________________________
pic = "Standard_deviation_difference_"+var_name(0)+"_"+scale
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"ViBlGrWhYeOrRe")
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

plots = new(7,graphic)

res@cnFillPalette = "ViBlGrWhYeOrRe"
levels = (/-30.,-25.,-20,-15.,-10.,-5.,5.,10.,15.,20.,25.,30./)
;levels = (/20.,40.,60.,80.,100.,120.,140.,160.,180.,200./)
res@cnLevels = levels

std(0,:,:) = std(0,:,:)-std(3,:,:)
std(1,:,:) = std(1,:,:)-std(3,:,:)
std(2,:,:) = std(2,:,:)-std(3,:,:)
std(4,:,:) = std(4,:,:)-std(3,:,:)
std(5,:,:) = std(5,:,:)-std(3,:,:)
std(6,:,:) = std(6,:,:)-std(3,:,:)
counter = 0
do case_num = 0,6
if case_num .ne. 3 then
  res@tiMainString = case_name_brf(case_num)+" - "+ case_name_brf(3)
  if scale .eq. "AU" then
    print(case_name_brf(case_num)+" "+avg(std(case_num,{-40:-27},{140:155})))
    plots(counter) = gsn_csm_contour_map(wks,std(case_num,{-40:-27},{140:155}),res)
  else if scale .eq. "World" then
    print(case_name_brf(case_num)+" "+avg(std(case_num,:,:)))
    plots(counter) = gsn_csm_contour_map(wks,std(case_num,:,:),res)
  end if
  end if
  counter = counter + 1
end if

end do
gsn_panel(wks,(/plots(0:5)/),(/3,2/),pres)

end
