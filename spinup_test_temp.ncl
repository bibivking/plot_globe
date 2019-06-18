;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin


case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on/Outputs/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_spinup/Outputs_61-120/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_spinup/Outputs_121-180/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_spinup/Outputs_181-240/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off/Outputs/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs-standard/gw_on"\
            /)

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers


;______________________________ Data _________________________________
year_s = 2001;1951
year_e = 2010
year_sum = year_e-year_s+1
case_sum = dimsizes(case_name)
total = year_sum * 12

scale = "AU";"Global" ;"Sahara"

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
else if scale .eq. "Sahara" then
  A = 14.
  B = 25.
  C = 0.
  D = 30.
end if
end if
end if
end if
end if


SoilTemp = new((/case_sum,6,total/),float)
SoilTemp = 0.0

;______________________________________ process data _________________________________
do case_num = 0,case_sum-1
   print(case_name(case_num))

   do year = year_s,year_e

      filename = case_name(case_num)+"/cable_out_"+year+".nc"
      f = addfile (filename,"r")
      SoilTemp@_FillValue = f->Qle@_FillValue

      time_s = (year-year_s)*12
      do month = 0,11
         do layers = 0,5
            SoilTemp(case_num,layers,time_s+month) = avg(f->SoilTemp(month,layers,{A:B},{C:D}))
         end do
      end do
   end do ;year
end do ;case

SoilTemp =  SoilTemp - 273.16

;_____________________________ PLOT ___________________________________
x = ispan(1,total,1)
print(x)
pic = "Spinup_check_temp_"+scale
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

res                    = True
res@gsnMaximize        = True
res@vpHeightF          = 0.6
res@vpWidthF           = 0.6
res@vpXF               = 0.3
res@vpYF               = 0.85
res@pmLegendWidthF     = 0.14                     ; Resize legend width
res@pmLegendHeightF    = 0.11                     ; and height

res@xyLineThicknesses      = (/1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)          ; make second line thicker
res@pmLegendDisplayMode    = "Always"
res@xyLineColors           = (/"black","blue","green","orange","pink","red"/)          ; change line color
res@xyExplicitLegendLabels = (/"layer1","layer2","layer3","layer4","layer5","layer6"/)
;res@trYMinF = 0

plots = new(6,graphic)
do case_num = 0,case_sum-1

if case_num .le. 3 then
   spinup_length = (case_num+1)*60
   res@tiMainString = "Spinup "+spinup_length+" years"
else if case_num .eq. 4 then
   res@tiMainString = "NoGW"
else if case_num .eq. 5 then
   res@tiMainString = "NoAQ"
end if
end if
end if

res@tiYAxisString  = " Soil Temperature (~S~o~N~C)"
res@tiXAxisString  = "month"

plots(case_num) = gsn_csm_xy(wks,x,SoilTemp(case_num,:,:),res)

end do

pres             = True
pres@gsnMaximize = True
gsn_panel(wks,(/plots/),(/3,2/),pres)

end
