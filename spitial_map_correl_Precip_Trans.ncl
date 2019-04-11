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

Precip_data = "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs-standard/gw_off"

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


;_____________________________ PLOT ___________________________________
pic = "Precip_"+var_name(0)+"_correl_"+scale
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
;levels = (/-1.,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
levels = (/-1.,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)
res@cnLevels = levels

;______________________________ Precip __________________________________
Precip  = new((/total,360,720/),float)
Precip  = 0.0
counter = 0

do year = year_s,year_e
  filename = Precip_data+"/cable_out_"+year+".nc"
  f = addfile (filename,"r")
  do month = 0,11
    Precip(counter,:,:) = f->Rainf(month,:,:)*1.0
    counter = counter +1
  end do ;month
end do ;year

Precip!0 = "time"
Precip!1 = "lat"
Precip!2 = "lon"
Precip&lat = f->Qle&y
Precip&lon = f->Qle&x
Precip@_FillValue   = f->Qle@_FillValue
plot = gsn_csm_contour_map(wks1,dim_avg_n_Wrap(Precip(:,:,:),0),res1)
plot = gsn_csm_contour(wks1,dim_avg_n_Wrap(Precip(:,:,:),0)*1.0,res1)

if deseason .eq. True then
  Precip = rmMonAnnCycTLL(Precip) ; dimensioned time x lat x lon
  plot = gsn_csm_contour_map(wks1,dim_avg_n_Wrap(Precip(:,:,:),0),res1)
  plot = gsn_csm_contour(wks1,dim_avg_n_Wrap(Precip(:,:,:),0)*1.0,res1)
end if

delete(counter)
delete(f)


do case_num = 0,6
  print(case_name(case_num))

  ;______________________________ INPUT __________________________________
  var          = new((/total,360,720/),float)
  Precip_var   = new((/360,720/),float)
  var          = 0.0
  Precip_var   = 0.0
  counter      = 0

  do year = year_s,year_e
    print("year = " +year)

    filename = case_name(case_num)+"/cable_out_"+year+".nc"
    f = addfile (filename,"r")

    do month = 0,11
      var(counter,:,:) = f->$var_name(0)$(month,:,:)*1.0
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

  ; _____________________  CALCULATE CORRELATION ________________________
  ;if case_num .eq. 6 then
  ;  do time = 0,total-1
  ;    var(time,:,:) = where( dim_stddev_n_Wrap(var(:,:,:),0) .eq. 0.0, var@_FillValue,  var(time,:,:))
  ;  end do
  ;end if

  print(case_name_brf(case_num))
  Precip_var = escorc(var(lat|:,lon|:,time|:), Precip(lat|:,lon|:,time|:))
; escorc: Computes the (Pearson) sample linear cross-correlations at lag 0 only.
; The rightmost dimension is usually time.
  delete(var)
  Precip_var!0 = "lat"
  Precip_var!1 = "lon"
  Precip_var&lat = f->Qle&y
  Precip_var&lon = f->Qle&x
  Precip_var@_FillValue   = f->Qle@_FillValue

  ; _______________________________ PLOT _________________________________
  res@tiMainString = case_name_brf(case_num)
  if scale .eq. "AU" then
    print(case_name_brf(case_num)+" "+avg(Precip_var({-40:-27},{140:155})))
    plots(case_num) = gsn_csm_contour_map(wks,Precip_var({-40:-27},{140:155}),res)
  else if scale .eq. "World" then
    print(case_name_brf(case_num)+" "+avg(Precip_var))
    plots(case_num) = gsn_csm_contour_map(wks,Precip_var,res)
  end if
  end if

  delete(Precip_var)

end do ;case

gsn_panel(wks,(/plots/),(/4,2/),pres)

end
