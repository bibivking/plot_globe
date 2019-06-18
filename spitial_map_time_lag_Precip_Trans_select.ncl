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
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-standard/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs-standard/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs-standard/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-Haverd2013/gw_on"\
             /)

case_name_brf = (/\
              "GW",\
              "NoGW",\
              "NoAQ",\
              "GW_Hvrd"\
               /)

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)

order = (/"(a) ","(b) ","(c) ","(d) ","(e) ","(f) ","(g) "/)

; ______________________________ Pre-define ________________________
year_s   = 2001
year_e   = 2010
year_sum = year_e-year_s+1
total    = year_sum*12
case_sum = dimsizes(case_name)

scale = "SE-AU" ;"Global" ; "50S-N" ; "SE-AU" ; "AU"

; _____________________________ For Check ___________________________
pic1 = "check"
wks1 = gsn_open_wks("pdf",pic1)
gsn_define_colormap(wks1,"ViBlGrWhYeOrRe")
res1 = True
res1@cnFillMode         = "RasterFill"            ; Raster Mode
res1@cnFillOn           = True                            ;��ɫ
res1@tmBorderThicknessF = 3.0
res1@cnLinesOn          = False                       ; �رյ�ֵ������

;_______________________________ PFT ___________________________________
filename = case_name(0)+"/cable_out_"+year_s+".nc"
f = addfile (filename,"r")
print(filename)
iveg = f->iveg

iveg!0 = "lat"
iveg!1 = "lon"
iveg&lat = f->Qle&y
iveg&lon = f->Qle&x
delete(f)
delete(filename)

;______________________________ Data read in __________________________________
Precip  = new((/total,360,720/),float)
Precip  = 0.0
counter = 0

do year = year_s,year_e
  print("year = " +year)
  if ( mod(year,4) .eq. 0) then
     dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
  else
     dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  end if
     
  filename = case_name(0)+"/cable_out_"+year+".nc"
  f = addfile (filename,"r")
  do month = 0,11
    Precip(counter,:,:) = f->Rainf(month,:,:)*3600*24*dom(month)
    counter = counter +1
  end do ;month
end do ;year

lat = f->Qle&y
lon = f->Qle&x

Precip!0 = "time"
Precip!1 = "lat"
Precip!2 = "lon"
Precip&lat = lat
Precip&lon = lon

Precip@_FillValue   = f->Qle@_FillValue
Precip = rmMonAnnCycTLL(Precip) ; dimensioned time x lat x lon
delete(counter)
delete(f)


time_lag = new((/case_sum,360,720/),integer)
time_lag = 0
time_lag@_FillValue = -999

time_lag!0 = "case" 
time_lag!1 = "lat"
time_lag!2 = "lon"
time_lag&lat = lat
time_lag&lon = lon

do case_num = 0,case_sum-1
  print(case_name(case_num))

  TVeg    = new((/total,360,720/),float)
  Correl  = new((/13,360,720/),float)
  TVeg    = 0.0
  Correl  = 0.0
  counter = 0

  if ( mod(year,4) .eq. 0) then
     dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
  else
     dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  end if

  do year = year_s,year_e

    filename = case_name(case_num)+"/cable_out_"+year+".nc"
    f = addfile (filename,"r")

    do month = 0,11
      TVeg(counter,:,:) = f->TVeg(month,:,:)*3600*24*dom(month)
      counter = counter + 1
    end do ;month

  end do ;year

  TVeg!0 = "time"
  TVeg!1 = "lat"
  TVeg!2 = "lon"
  TVeg&lat = lat
  TVeg&lon = lon
  TVeg@_FillValue   = Precip@_FillValue
  plot = gsn_csm_contour_map(wks1,TVeg(6,:,:),res1)

  ; __________________________ Deseasonalize _____________________________
  ; rmMonAnnCycTLL : Removes the annual cycle from "monthly" data.
  TVeg = rmMonAnnCycTLL(TVeg)
  plot = gsn_csm_contour_map(wks1,TVeg(6,:,:),res1)

  ; _____________________  CALCULATE CORRELATION ________________________
  print(case_name_brf(case_num))
  Correl(0,:,:)  = escorc(TVeg(lat|:,lon|:,time|:), Precip(lat|:,lon|:,time|:))
  plot = gsn_csm_contour(wks1,abs(Correl(0,:,:)),res1)
  Correl(1,:,:)  = escorc(TVeg(lat|:,lon|:,time|1:), Precip(lat|:,lon|:,time|:(total-2)))
  plot = gsn_csm_contour(wks1,abs(Correl(1,:,:)),res1)
  Correl(2,:,:)  = escorc(TVeg(lat|:,lon|:,time|2:), Precip(lat|:,lon|:,time|:(total-3)))
;  plot = gsn_csm_contour(wks1,Correl(2,:,:),res1)
  Correl(3,:,:)  = escorc(TVeg(lat|:,lon|:,time|3:), Precip(lat|:,lon|:,time|:(total-4)))
;  plot = gsn_csm_contour(wks1,Correl(3,:,:),res1)
  Correl(4,:,:)  = escorc(TVeg(lat|:,lon|:,time|4:), Precip(lat|:,lon|:,time|:(total-5)))
;  plot = gsn_csm_contour(wks1,Correl(4,:,:),res1)
  Correl(5,:,:)  = escorc(TVeg(lat|:,lon|:,time|5:), Precip(lat|:,lon|:,time|:(total-6)))
  Correl(6,:,:)  = escorc(TVeg(lat|:,lon|:,time|6:), Precip(lat|:,lon|:,time|:(total-7)))
  Correl(7,:,:)  = escorc(TVeg(lat|:,lon|:,time|7:), Precip(lat|:,lon|:,time|:(total-8)))
  Correl(8,:,:)  = escorc(TVeg(lat|:,lon|:,time|8:), Precip(lat|:,lon|:,time|:(total-9)))
  Correl(9,:,:)  = escorc(TVeg(lat|:,lon|:,time|9:), Precip(lat|:,lon|:,time|:(total-10)))
  Correl(10,:,:) = escorc(TVeg(lat|:,lon|:,time|10:), Precip(lat|:,lon|:,time|:(total-11)))
  Correl(11,:,:) = escorc(TVeg(lat|:,lon|:,time|11:), Precip(lat|:,lon|:,time|:(total-12)))
  Correl(12,:,:) = escorc(TVeg(lat|:,lon|:,time|12:), Precip(lat|:,lon|:,time|:(total-13)))
  ; Note that the time periods are different. This may cause correl's different.
  ; escorc: Computes the (Pearson) sample linear cross-correlations at lag 0 only.
  ; the rightmost dimension is usually time.
  delete(TVeg)
  
  Correl!0 = "time_lag"
  Correl!1 = "lat"
  Correl!2 = "lon"
  Correl&lat = lat
  Correl&lon = lon
  Correl@_FillValue   = Precip@_FillValue
  Correl_max = abs(Correl(1,:,:))
  Correl_num = new((/360,720/),integer)
  Correl_num = 1
  
  do time = 0,12
     Correl_num = where(Correl_max .ge. abs(Correl(time,:,:)), Correl_num, time )
     Correl_max = where(Correl_max .ge. abs(Correl(time,:,:)), Correl_max, abs(Correl(time,:,:)) )
  end do
  time_lag(case_num,:,:) = Correl_num*1
  
  plot = gsn_csm_contour(wks1,abs(Correl_max),res1)
  plot = gsn_csm_contour(wks1,Correl_num,res1)
  plot = gsn_csm_contour(wks1,time_lag(case_num,:,:),res1)
  delete(Correl)
  delete(Correl_num)
  delete(Correl_max)
end do

;_____________________________ PLOT ___________________________________
pic = "spitial_map_time_lag_Precip_Trans_select_"+scale

wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"ViBlGrWhYeOrRe")
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
res@txString           = ""
res@tmXBLabelFont      = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont      = 25

;*************** ���õ�ֵ�� *****************
res@cnLinesOn          = False                       ; �رյ�ֵ������
res@cnLineColor        = "black"
res@cnLineThicknessF   = 1.5
res@cnLineLabelsOn     = False
res@gsnMaximize        = False
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
pres@lbTitleFontHeightF = 0.013
pres@lbLabelFontHeightF = 0.013

plots = new(8,graphic)


; _______________________________ PLOT _________________________________
res@cnFillPalette = "BlGrYeOrReVi200"
levels = (/0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5/)
res@cnLevels = levels

; for all PFTs
printVarSummary(time_lag)

tmp = where(iveg .ge. 1 .and. iveg .le. 4, time_lag(0,:,:),time_lag@_FillValue) ; for forest
tmp!0 = "lat"
tmp!1 = "lon"
tmp&lat = lat
tmp&lon = lon
res@tiMainString = "(a) GW forest"
plots(0) = gsn_csm_contour_map(wks,tmp({A:B},{C:D}),res)
delete(tmp)
tmp = where(iveg .ge. 5 .and. iveg .le. 12, time_lag(0,:,:),time_lag@_FillValue) ; for grass+shrub+crop
tmp!0 = "lat"
tmp!1 = "lon"
tmp&lat = lat
tmp&lon = lon
res@tiMainString = "(b) GW other"
plots(1) = gsn_csm_contour_map(wks,tmp({A:B},{C:D}),res)
delete(tmp)

tmp = where(iveg .ge. 1 .and. iveg .le. 4, time_lag(1,:,:),time_lag@_FillValue) ; for forest
tmp!0 = "lat"
tmp!1 = "lon"
tmp&lat = lat
tmp&lon = lon
res@tiMainString = "(c) NoGW forest"
plots(2) = gsn_csm_contour_map(wks,tmp({A:B},{C:D}),res)
delete(tmp)
tmp = where(iveg .ge. 5 .and. iveg .le. 12, time_lag(1,:,:),time_lag@_FillValue) ; for grass+shrub+crop
tmp!0 = "lat"
tmp!1 = "lon"
tmp&lat = lat
tmp&lon = lon
res@tiMainString = "(d) NoGW other"
plots(3) = gsn_csm_contour_map(wks,tmp({A:B},{C:D}),res)
delete(tmp)


tmp = where(iveg .ge. 1 .and. iveg .le. 4, time_lag(2,:,:),time_lag@_FillValue) ; for forest
tmp!0 = "lat"
tmp!1 = "lon"
tmp&lat = lat
tmp&lon = lon
res@tiMainString = "(e) NoAQ forest"
plots(4) = gsn_csm_contour_map(wks,tmp({A:B},{C:D}),res)
delete(tmp)
tmp = where(iveg .ge. 5 .and. iveg .le. 12, time_lag(2,:,:),time_lag@_FillValue) ; for grass+shrub+crop
tmp!0 = "lat"
tmp!1 = "lon"
tmp&lat = lat
tmp&lon = lon
res@tiMainString = "(f) NoAQ other"
plots(5) = gsn_csm_contour_map(wks,tmp({A:B},{C:D}),res)
delete(tmp)

tmp = where(iveg .ge. 1 .and. iveg .le. 4, time_lag(3,:,:),time_lag@_FillValue) ; for forest
tmp!0 = "lat"
tmp!1 = "lon"
tmp&lat = lat
tmp&lon = lon
res@tiMainString = "(g) GW_Hvrd forest"
plots(6) = gsn_csm_contour_map(wks,tmp({A:B},{C:D}),res)
delete(tmp)
tmp = where(iveg .ge. 5 .and. iveg .le. 12, time_lag(3,:,:),time_lag@_FillValue) ; for grass+shrub+crop
tmp!0 = "lat"
tmp!1 = "lon"
tmp&lat = lat
tmp&lon = lon
res@tiMainString = "(h) GW_Hvrd other"
plots(7) = gsn_csm_contour_map(wks,tmp({A:B},{C:D}),res)
delete(tmp)

gsn_panel(wks,(/plots/),(/4,2/),pres)


end
