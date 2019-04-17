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
              "off-std",\
              "off-redstrb-std",\
              "off-Hvrd",\
              "on-std",\
              "on_norchrg-std",\
              "on-Hvrd",\
              "on_norchrg-Hvrd"\
               /)
; ____________________________ Options _______________________________

scale = "World" ; "AU"


; fill gap for 2003 June
year_s = 2003
year_e = 2010
year_sum = year_e-year_s+1
var_sum  = dimsizes(var_name)

; _____________________________ Check ______________________________
pic = "check"
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"ViBlGrWhYeOrRe")
res = True
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                            ;��ɫ
res@tmBorderThicknessF = 3.0
res@cnLinesOn          = False                       ; �رյ�ֵ������
res@gsnDraw            = True  ; Don't draw plots
res@gsnFrame           = True  ; ����ҳ

; _______________________ read GRACE and fill gap ____________________________
GRACE      = new((/96,360,720/),float)
uncertainty= new((/96,360,720/),float)

GRACE      = 0.
uncertainty= 0.

filename = "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/plots/GRACE_JPL_RL06/GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01CRIv01.nc"
f = addfile (filename,"r")

GRACE(0:4,:,:)       = f->lwe_thickness(7:11,:,:)*10.
GRACE(6:95,:,:)      = f->lwe_thickness(12:101,:,:)*10.

uncertainty(0:4,:,:) = f->uncertainty(7:11,:,:)*10.
uncertainty(6:95,:,:)= f->uncertainty(12:101,:,:)*10.

Jun_ave = new((/360,720/),float)
Jul_ave = new((/360,720/),float)
Jun_ave = 0.
Jul_ave = 0.

GRACE@_FillValue       = f->lwe_thickness@_FillValue
uncertainty@_FillValue = f->lwe_thickness@_FillValue


plot = gsn_csm_contour(wks,dim_avg_n_Wrap(GRACE,0),res)

do year = year_s+1,year_e

Jun_ave = Jun_ave + GRACE((year-year_s)*12+5,:,:)
Jul_ave = Jul_ave + GRACE((year-year_s)*12+6,:,:)

end do

Jun_ave = Jun_ave/(year_e-year_s)
Jul_ave = Jul_ave/(year_e-year_s)

plot = gsn_csm_contour(wks,Jun_ave,res)
plot = gsn_csm_contour(wks,Jul_ave,res)

GRACE(5,:,:) = GRACE(6,:,:) * (Jun_ave/Jul_ave)
uncertainty(5,:,:) = (uncertainty(4,:,:) + uncertainty(6,:,:))/2.

delete(Jun_ave)
delete(Jul_ave)

plot = gsn_csm_contour(wks,GRACE(5,:,:),res)
plot = gsn_csm_contour(wks,uncertainty(5,:,:),res)




;______________________________ INPUT __________________________________
TWS     = new((/7,96,360,720/),float)
TWS_ave = new((/7,360,720/),float)
TWS     = 0.0
TWS_ave = 0.0

do case_num = 0,6

   print(case_name(case_num))

   counter = 0

   do year = year_s,year_e
      print("year = " +year)
    
      filename = case_name(case_num)+"/cable_out_"+year+".nc"
      f = addfile (filename,"r")

      do month = 0,11
         if case_num .le. 2 then
            TWS(case_num,counter,:,:) = f->SWE(month,:,:) + f->CanopInt(month,:,:)
         else
            TWS(case_num,counter,:,:) = f->GWMoist(month,:,:)*f->GWdz*1000. + f->SWE(month,:,:) + f->CanopInt(month,:,:)
         end if
         do layer = 0,5
            TWS(case_num,counter,:,:) = TWS(case_num,counter,:,:) + f->SoilMoist(month,layer,:,:)*Soil_thickness(layer)*1000.
         end do
         counter = counter + 1
      end do ;month
   end do ;year

end do ;case

; Calculate TWS average
TWS_ave =  dim_avg_n_Wrap(TWS(:,12:83,:,:),1)

do time = 0,95
  TWS(:,time,:,:) = TWS(:,time,:,:) - TWS_ave
end if

delete(TWS_ave)

print("point 1 ")

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

; _________________________________  Global AVE  ______________________________________
GRACE_global_ave       = new(96,float)
uncertainty_global_ave = new(96,float)
TWS_global_ave         = new((/7,96/),float)
GRACE_global_ave       = 0.0
uncertainty_global_ave = 0.0
TWS_global_ave         = 0.0

GRACE_TWS              = new(7,float)
GRACE_TWS              = 0.0

do time = 0,95
   GRACE_global_ave(time) = avg(GRACE(time,:,:))
   uncertainty_global_ave(time) = avg(uncertainty(time,:,:))
   do case_num = 0,6
      TWS_global_ave(case_num,time) = avg(TWS(case_num,time,:,:)))
   end do
end do

do case_num = 0,6
   GRACE_TWS(case_num) = escorc(GRACE_global_ave, TWS_global_ave(case_num,:))
   print(case_name_brf(case_num)+" "+GRACE_TWS(case_num))
end do

;GRACE & TWS is anormaly here

; ______________________________ Calculate Variance __________________________________
GRACE_variance = dim_variance_n_Wrap(GRACE,0)
TWS_variance   = dim_variance_n_Wrap(TWS,1)

variance_diff  = new((/7,360,720/),float)
variance_diff  = 0.0

do case_num = 0,6
   variance_diff(case_num,:,:) = TWS_variance(case_num,:,:) - GRACE_variance
end do

; ____________________________ Deseasonalize & Detrend ________________________________
GRACE = rmMonAnnCycTLL(GRACE) ; dimensioned time x lat x lon
do case_num = 0,6
   TWS(case_num,:,:,:) = rmMonAnnCycTLL(TWS(case_num,:,:,:))
end do

plot = gsn_csm_contour(wks,GRACE(0,:,:),res)
plot = gsn_csm_contour(wks,TWS(0,:,:),res)

GRACE = dtrend_n(GRACE, True, 0)
TWS   = dtrend_n(TWS,   True, 1)

GRACE_year = new((/8,360,720/),float)
TWS_year   = new((/7,8,360,720/),float)

do year = year_s,year_e
   a = (year-year_s)*12
   b = (year-year_s)*12 +11
   GRACE_year(year-year_s,:,:) = dim_avg_n_Wrap(GRACE(a:b,:,:),0)
   TWS_year(:,year-year_s,:,:) = dim_avg_n_Wrap(TWS(:,a:b,:,:),1)
end do

plot = gsn_csm_contour(wks,GRACE_year(0,:,:),res)
plot = gsn_csm_contour(wks,TWS_year(0,:,:),res)

;GRACE & TWS is deseasonalized and detrended here

; ____________________________ Calculate fraction of IAV ________________________________
GRACE_FIAV = new((/360,720/),float)
TWS_FIAV   = new((/7,360,720/),float)

GRACE_FIAV = dim_variance_n_Wrap(GRACE_year,0)/dim_variance_n_Wrap(GRACE,0)
TWS_FIAV   = dim_variance_n_Wrap(TWS_year,1)/dim_variance_n_Wrap(TWS,1)


delete(GRACE_year)
delete(TWS_year)

delete(GRACE)
delete(TWS)


printVarSummary(GRACE_FIAV)
printVarSummary(TWS_FIAV)

GRACE_FIAV!0 = "lat"
GRACE_FIAV!1 = "lon"
GRACE_FIAV&lat = f->Qle&y
GRACE_FIAV&lon = f->Qle&x

TWS_FIAV!0 = "case"
TWS_FIAV!1 = "lat"
TWS_FIAV!2 = "lon"
TWS_FIAV&lat = f->Qle&y
TWS_FIAV&lon = f->Qle&x

plot = gsn_csm_contour(wks,GRACE_FIAV(0,:,:),res)
plot = gsn_csm_contour(wks,TWS_FIAV(0,:,:),res)

delete(res)
delete(pic)
delete(wks)


;_____________________________ PLOT ___________________________________
pic = "GRACE_vs_CABLE_"+scale

wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")

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

if scale .eq. "World" then
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

   res@tmYLMode  = "Explicit"             ;��Y����������������,��ʾ�Լ���Ҫ��ֵ
   res@tmYLValues=(/-30,-40/)
   res@tmYLLabels=(/"30~S~o~N~S","40~S~o~N~S"/)
   res@tmXBMode  = "Explicit"
   res@tmXBValues=(/140,150/)
   res@tmXBLabels=(/"140~S~o~N~E","150~S~o~N~E"/)
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

; ____________________ all PFT __________________________
print("for all PFTs")
do case_num = 0,6

res@tiMainString = case_name_brf(case_num)

if scale .eq. "World" then
   print(case_name_brf(case_num)+" "+avg(GRACE_TWS(case_num,:,:)))
   levels = (/0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
   res@cnLevels = levels
   plots(case_num) = gsn_csm_contour_map(wks,GRACE_TWS(case_num,:,:),res)
else if scale .eq. "AU" then
   print(case_name_brf(case_num)+" "+avg(GRACE_TWS(case_num,{-40:-27},{140:155})))
   levels = (/0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6/)
   res@cnLevels = levels
   plots(case_num) = gsn_csm_contour_map(wks,GRACE_TWS(case_num,{-40:-27},{140:155}),res)
end if
end if

end do
gsn_panel(wks,(/plots/),(/4,2/),pres)

; tmp var
tmp   = new((/7,360,720/),float)
tmp   = 0.0
tmp!0 = "case"
tmp!1 = "lat"
tmp!2 = "lon"
tmp&lat = f->Qle&y
tmp&lon = f->Qle&x
tmp@_FillValue   = f->Qle@_FillValue

; ____________________ forest __________________________
print("for forest")
pres@txString = "forest"
do case_num = 0,6
  tmp(case_num,:,:) = where(iveg .ge. 1 .and. iveg .le. 4, GRACE_TWS(case_num,:,:), GRACE_TWS@_FillValue )
  res@tiMainString = case_name_brf(case_num)
  if scale .eq. "World" then
    print(case_name_brf(case_num)+" "+avg(tmp(case_num,:,:)))
    levels = (/0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
    res@cnLevels = levels
    plots(case_num) = gsn_csm_contour_map(wks,tmp(case_num,:,:),res)
  else if scale .eq. "AU" then
    print(case_name_brf(case_num)+" "+avg(tmp(case_num,{-40:-27},{140:155})))
    levels = (/0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6/)
    res@cnLevels = levels
    plots(case_num) = gsn_csm_contour_map(wks,tmp(case_num,{-40:-27},{140:155}),res)
  end if
  end if
end do
gsn_panel(wks,(/plots/),(/4,2/),pres)


; ____________________ grass+shrub+crop __________________________
print("for grass+shrub+crop")
pres@txString = "grass+shrub+crop"
do case_num = 0,6
  tmp(case_num,:,:) = where(iveg .ge. 5 .and. iveg .le. 11, GRACE_TWS(case_num,:,:), GRACE_TWS@_FillValue )
  res@tiMainString = case_name_brf(case_num)
  if scale .eq. "World" then
    print(case_name_brf(case_num)+" "+avg(tmp(case_num,:,:)))
    levels = (/0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
    res@cnLevels = levels
    plots(case_num) = gsn_csm_contour_map(wks,tmp(case_num,:,:),res)
  else if scale .eq. "AU" then
    print(case_name_brf(case_num)+" "+avg(tmp(case_num,{-40:-27},{140:155})))
    levels = (/0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6/)
    res@cnLevels = levels
    plots(case_num) = gsn_csm_contour_map(wks,tmp(case_num,{-40:-27},{140:155}),res)
  end if
  end if
end do
gsn_panel(wks,(/plots/),(/4,2/),pres)

end
