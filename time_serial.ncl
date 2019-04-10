;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs/gw_on/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs/gw_on/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs/gw_off/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/plots/PDSI/pdsisc.monthly.maps.1850-2014.nc"/)

var_name= (/\
           "Evap","TVeg","ESoil",\ ;"kg/m^2/s"
           "Qs","Qsb",\      ;"kg/m^2/s"
           "Qle","Qh",\      ;"W/m^2"
           "WatTable",\       ;"m"
           "GPP","NPP",\     ;"umol/m^2/s"
           "SoilMoist",\     ;"m^3/m^3"  Average layer soil moisture
           "GWMoist"\       ;"mm3/mm3" Aquifer mositure content, Groundwater_off doesn't have
           /)
           ;"LAI" ; the same in Groundwater_on and Groundwater_off

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872, 20./) ; 6 soil layers and 1 aquifer


;______________________________ Data _________________________________
year_s = 2001
year_e = 2010
year_sum = year_e-year_s+1
var_sum  = dimsizes(var_name)
total = year_sum * 12

var = new((/3,10,total/),float)
var = 0.0
Moist = new((/3,7,total/),float)
Moist = 0.0
TWS   = new((/3,total/),float)
TWS   = 0.0


;______________________________________ process data _________________________________
do case_num = 0,2
   print(case_name(case_num))
   counter = 0
   do year = year_s,year_e
      print("year = " +year)
      if ( mod(year,4) .eq. 0) then
         dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
      else
         dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      end if

      filename = case_name(case_num)+"/cable_out_"+year+".nc"
      f = addfile (filename,"r")

      var@_FillValue   = f->Qle@_FillValue
      Moist@_FillValue = f->Qle@_FillValue
      TWS@_FillValue   = f->Qle@_FillValue
      ;print(f->Qle@_FillValue)

      do month = 0,11
         do var_num = 0,9
            if (var_num .le. 4) then
               var(case_num, var_num, counter) = avg(f->$var_name(var_num)$(month,{-38:-28},{140.5:154}))*3600*24*dom(month)
            else if (var_num .le. 8) then
               var(case_num, var_num, counter) = avg(f->$var_name(var_num)$(month,{-38:-28},{140.5:154}))
            end if
            end if
         end do ;var_num

         do layers = 0,5
            Moist(case_num,layers,counter) = avg(f->SoilMoist(month,layers,{-38:-28},{140.5:154}))
         end do
         if (case_num .eq. 1)
            Moist(case_num,6,counter)      = 0.0
         else
            Moist(case_num,6,counter)      = avg(f->GWMoist(month,{-38:-28},{140.5:154}))
         end if

         do layers = 0,5
         ;  m3/m3
            TWS(case_num,counter) = Moist(case_num,layers,counter)*Soil_thickness(layers)/sum(Soil_thickness(0:5))
         end do ;varend
      counter = counter + 1
     end do ;month
   end do ;year

end do ;case
; _________________________ Calculate Anomaly ___________________________
var_avg =  dim_avg_n_Wrap(var,2)
Moist_avg =  dim_avg_n_Wrap(Moist,2)
TWS_avg =  dim_avg_n_Wrap(TWS,1)

var_anomaly = new((/3,10,total/),float)
var_anomaly = 0.0
Moist_anomaly = new((/3,7,total/),float)
Moist_anomaly = 0.0
TWS_anomaly   = new((/3,total/),float)
TWS_anomaly   = 0.0

do case_num = 0,2

do var_num = 0,9
   var_anomaly(case_num,var_num,:)= var(case_num,var_num,:)-var_avg(case_num,var_num)
end do

do layers = 0,6
   Moist_anomaly(case_num,layers,:)= Moist(case_num,layers,:)-Moist_avg(case_num,layers)
end do

   TWS_anomaly(case_num,:)= TWS(case_num,:)-TWS_avg(case_num)
end do

; ___________________________ PDSI ______________________________
PDSI = new(120,float)
PDSI = 0.0
filename1 = case_name(3)
f1 = addfile (filename1,"r")
do cnt = 1812,1931
   PDSI(cnt-1812) = avg(f1->sc_PDSI_pm(cnt,{-38:-28},{140.5:154}))
end do

;_____________________________ PLOT ___________________________________
x = ispan(1,total,1)
print(x)
pic = "GW_time_serial_anomaly";
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

res                    = True
res@gsnMaximize        = False
res@gsnDraw            = False  ; Don't draw plots
res@gsnFrame           = False  ; ����ҳ

res@vpHeightF          = 0.35
res@vpWidthF           = 0.7
;res@vpXF               = 0.0
;res@vpYF               = 0.8
res@pmLegendWidthF     = 0.15                      ; Resize legend width
res@pmLegendHeightF    = 0.15                     ; and height
res@pmLegendSide       = "Right"
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015


;res@trYMinF = 0

pres                    = True                                      ; ���� panel resources.
pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
pres@lbBoxLinesOn       = False                       ;�ر�lbar box ����
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
pres@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
pres@lbTitleFontHeightF = 0.5
pres@lbLabelFontHeightF = 0.5
pres@lbAutoManage       = True
pres@lbTitleDirection   = "Across"             ; title direction
pres@lbTitlePosition    = "Bottom";Right"


plots = new(5,graphic)

res@tmXBMode  = "Explicit"
res@tmXBValues=(/6,12,18,24,30,36,42,48,54,60,66,72,78,84,90,96,102,108,114,120/)
res@tmXBLabels=(/"J","D","J","D","J","D","J","D","J","D","J","D","J","D","J","D",\
                 "J","D","J","D","J","D","J","D"/)

res@pmLegendDisplayMode = "Always"

res@xyLineThicknesses  = (/  2.0, 2.0, 2.0/)          ; make second line thicker
res@xyLineColor           = "black"          ; change line color
res@xyExplicitLegendLabel = "PDSI"
res@xyDashPattern         = 0

;res@trYMaxF = 1.3*max(var(:,0,:))
;res@tiMainString = f->Evap@long_name+" (mm/month)"
res@tiYAxisString  = "PDSI (-)"
;res@tiXAxisString  = "month since Jan 2001"
plots(0) = gsn_csm_xy(wks,x,PDSI,res)

delete(res@xyLineColor)
delete(res@xyExplicitLegendLabel)
delete(res@xyDashPatterns)

res@xyLineThicknesses  = (/  2.0, 2.0, 2.0/)          ; make second line thicker
res@xyLineColors       = (/"blue","red","darkgreen"/)          ; change line color
res@xyExplicitLegendLabels = (/"GW","Norecharge","NoGW"/)
res@xyDashPatterns    = (/0,0,0/)

res@pmLabelBarDisplayMode = "Never"
;res@trYMaxF = 1.3*max(var(:,0,:))
;res@tiMainString = f->Evap@long_name+" (mm/month)"
res@tiYAxisString  = "Evap_anomaly (mm/month)"
;res@tiXAxisString  = "month since Jan 2001"
plots(1) = gsn_csm_xy(wks,x,var_anomaly(:,0,:),res) ;_anomaly

;res@trYMaxF = 1.3*max(var(:,1,:))
;res@tiMainString = f->TVeg@long_name+" (mm/month)"
res@tiYAxisString  = "TVeg_anomaly (mm/month)"
;res@tiXAxisString  = "month since Jan 2001"
plots(2) = gsn_csm_xy(wks,x,var_anomaly(:,1,:),res) ;_anomaly

;res@trYMaxF = 1.3*max(var(:,2,:))
;res@tiMainString = f->ESoil@long_name+" (mm/month)"
res@tiYAxisString  = "ESoil_anomaly (mm/month)"
;res@tiXAxisString  = "month since Jan 2001"
plots(3) = gsn_csm_xy(wks,x,var_anomaly(:,2,:),res) ;_anomaly

;delete(res@tmXBMode)
;delete(res@tmXBValues)
;delete(res@tmXBLabels)

res@pmLabelBarDisplayMode = "Always"
;res@trYMaxF = 1.3*max(TWS(:,:))
;res@tiMainString = "Soil Moisture (m~S~3~N~/m~S~3~N~)"
res@tiYAxisString  = "Soil Moisture_anomaly (m~S~3~N~/m~S~3~N~)"
res@tiXAxisString  = "month since Jan 2001"
plots(4) = gsn_csm_xy(wks,x,TWS_anomaly(:,:),res) ;_anomaly

gsn_panel(wks,(/plots/),(/5,1/),pres)

;res@tmXBMode  = "Explicit"
;res@tmXBValues=(/1,13,25,37,49,61,73,85,97,109/)
;res@tmXBLabels=(/"2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010"/)

;do var_num = 0,4
;res@trYMaxF = 1.3*max(var(:,var_num,:))
;res@tiYAxisString  = f->$var_name(var_num)$@long_name+" (mm/month)"
;res@tiXAxisString  = "month since Jan 2001"
;plot = gsn_csm_xy(wks,x,var(:,var_num,:),res)
;end do

;do var_num = 5,8
;res@trYMaxF = 1.3*max(var(:,var_num,:))
;res@tiYAxisString  = f->$var_name(var_num)$@long_name+" ("+f->$var_name(var_num)$@units+")"
;res@tiXAxisString  = "month since Jan 2001"
;plot = gsn_csm_xy(wks,x,var(:,var_num,:),res)
;end do

;res@trYMaxF = 1.3*max(var(:,9,:))
;res@tiYAxisString  = f->$var_name(9)$@long_name+" ("+f->$var_name(9)$@units+")"
;res@tiXAxisString  = "month since Jan 2001"
;plot = gsn_csm_xy(wks,x,var(0:1,9,:),res)


;do layers = 0,6
;  res@trYMaxF = 1.3*max(Moist(:,layers,:))
;  res@tiYAxisString  = "The soil moisture of the "+(layers+1)+ " layer (m~S~3~N~/m~S~3~N~)"
;  res@tiXAxisString  = "month since Jan 2001"
;  plot = gsn_csm_xy(wks,x,Moist(:,layers,:),res)
;end do


;res@trYMaxF = 1.3*max(TWS)
;res@tiYAxisString  = "Soil Moisture (m~S~3~N~/m~S~3~N~)"
;res@tiXAxisString  = "month since Jan 2001"
;plot = gsn_csm_xy(wks,x,TWS,res)

end
