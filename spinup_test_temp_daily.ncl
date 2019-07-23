;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin


case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_daily/Outputs/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_daily/Outputs/gw_off"\
            /)

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers


;______________________________ Data _________________________________
year_s = 1951
year_e = 1951
year_sum = year_e-year_s+1
case_sum = dimsizes(case_name)
total = 365

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


SoilTemp = new((/case_sum,6,50/),float)
SoilTemp = 0.0

;______________________________________ process data _________________________________
do case_num = 0,case_sum-1
   print(case_name(case_num))

   do year = year_s,year_e

      filename = case_name(case_num)+"/cable_out_"+year+".nc"
      f = addfile (filename,"r")
      SoilTemp@_FillValue = f->Qle@_FillValue

      do day = 0,49
         do layers = 0,5
            ;SoilTemp(case_num,layers,day) = avg(f->SoilTemp(day,layers,{A:B},{C:D}))
            SoilTemp(case_num,layers,day) = avg(f->SoilMoist(day,layers,{A:B},{C:D}))
         end do
      end do
   end do ;year
end do ;case

;SoilTemp =  SoilTemp - 273.16

;_____________________________ PLOT ___________________________________
x = ispan(1,50,1)
print(x)
;pic = "Spinup_check_temp_daily_"+scale
pic = "Spinup_check_moist_daily_"+scale
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
res@xyLineColors           = (/"red","blue","green","orange","pink","red"/)          ; change line color
res@xyExplicitLegendLabels = (/"on","off","layer3","layer4","layer5","layer6"/)
;res@trYMinF = 0

plots = new(6,graphic)

do layer = 0,5

;if (case_num .eq. 0) then
   res@tiMainString = "layer "+(layer+1)
;else 
;   res@tiMainString = "off"
;end if

;res@tiYAxisString  = " Soil Temperature (~S~o~N~C)"
res@tiYAxisString  = " Soil Moisture (mm~S~3~N~/mm~S~3~N~)"
res@tiXAxisString  = "month"

plots(layer) = gsn_csm_xy(wks,x,SoilTemp(:,layer,:),res)

end do

pres             = True
pres@gsnMaximize = True
gsn_panel(wks,(/plots/),(/3,2/),pres)

end
