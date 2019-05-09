;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

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

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers
order = (/"(a) ","(b) ","(c) ","(d) ","(e) ","(f) ","(g) "/)


;______________________________ Predefined __________________________________

flag = "Haverd2013" ;"standard" ;  Haverd2013

if flag .eq. "standard" then

case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-standard/gw_on/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs-standard/gw_on/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs-standard/gw_off/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/plots/PDSI/pdsisc.monthly.maps.1850-2014.nc"/)

else if flag .eq. "Haverd2013" then

case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-Haverd2013/gw_on/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs-Haverd2013/gw_on/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs-Haverd2013/gw_off/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/plots/PDSI/pdsisc.monthly.maps.1850-2014.nc"/)

end if
end if

;______________________________ Data _________________________________

year_s = 2001
year_e = 2010
year_sum = year_e-year_s+1
var_sum  = dimsizes(var_name)
total = year_sum * 12

; ________ 1 forest; 2 no forest _____________
var = new((/6,10,total/),float)
var = 0.0
Moist = new((/6,7,total/),float)
Moist = 0.0
TWS   = new((/6,total/),float)
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
      iveg = f->iveg
      do month = 0,11
         do var_num = 0,9
            if (var_num .le. 4) then
               var(case_num, var_num, counter) = sum( where(iveg({-40:-27},{140:155}) .ge. 1 .and. iveg({-40:-27},{140:155}) .le. 4,\
                                                  f->$var_name(var_num)$(month,{-40:-27},{140:155}),0.))*3600*24 \
                                                  / sum( where(iveg({-40:-27},{140:155}) .ge. 1 .and. iveg({-40:-27},{140:155}) .le. 4,\
                                                  1.,0.))

               var(case_num+3, var_num, counter) = sum( where(iveg({-40:-27},{140:155}) .lt. 1 .or. iveg({-40:-27},{140:155}) .gt. 4,\
                                                  f->$var_name(var_num)$(month,{-40:-27},{140:155}),0.))*3600*24 \
                                                  / sum( where(iveg({-40:-27},{140:155}) .lt. 1 .or. iveg({-40:-27},{140:155}) .gt. 4,\
                                                  1.,0.))
            else if (var_num .le. 8) then
               var(case_num, var_num, counter) = sum( where(iveg({-40:-27},{140:155}) .ge. 1 .and. iveg({-40:-27},{140:155}) .le. 4,\
                                                  f->$var_name(var_num)$(month,{-40:-27},{140:155}),0.)) \
                                                  / sum( where(iveg({-40:-27},{140:155}) .ge. 1 .and. iveg({-40:-27},{140:155}) .le. 4,\
                                                  1.,0.))
               var(case_num+3, var_num, counter) = sum( where(iveg({-40:-27},{140:155}) .lt. 1 .or. iveg({-40:-27},{140:155}) .gt. 4,\
                                                  f->$var_name(var_num)$(month,{-40:-27},{140:155}),0.)) \
                                                  / sum( where(iveg({-40:-27},{140:155}) .lt. 1 .or. iveg({-40:-27},{140:155}) .gt. 4,\
                                                  1.,0.))
            end if
            end if
         end do ;var_num

         do layers = 0,5
            Moist(case_num,layers,counter) = sum(where(iveg({-40:-27},{140:155}) .ge. 1 .and. iveg({-40:-27},{140:155}) .le. 4,\
                                              f->SoilMoist(month,layers,{-40:-27},{140:155}),0.)) \
                                              / sum( where(iveg({-40:-27},{140:155}) .ge. 1 .and. iveg({-40:-27},{140:155}) .le. 4,\
                                              1.,0.))
            Moist(case_num+3,layers,counter) = sum(where(iveg({-40:-27},{140:155}) .lt. 1 .or. iveg({-40:-27},{140:155}) .gt. 4,\
                                              f->SoilMoist(month,layers,{-40:-27},{140:155}),0.)) \
                                              / sum( where(iveg({-40:-27},{140:155}) .lt. 1 .or. iveg({-40:-27},{140:155}) .gt. 4,\
                                              1.,0.))
         end do
         if (case_num .eq. 1)
            Moist(case_num,6,counter)      = 0.0
            Moist(case_num+3,6,counter)    = 0.0
         else
            Moist(case_num,6,counter)      = sum(where(iveg({-40:-27},{140:155}) .ge. 1 .and. iveg({-40:-27},{140:155}) .le. 4,\
                                              f->GWMoist(month,{-40:-27},{140:155}),0.)) \
                                              / sum( where(iveg({-40:-27},{140:155}) .ge. 1 .and. iveg({-40:-27},{140:155}) .le. 4,\
                                              1.,0.))
            Moist(case_num+3,6,counter)    = sum(where(iveg({-40:-27},{140:155}) .lt. 1 .or. iveg({-40:-27},{140:155}) .gt. 4,\
                                              f->GWMoist(month,{-40:-27},{140:155}),0.)) \
                                              / sum( where(iveg({-40:-27},{140:155}) .lt. 1 .or. iveg({-40:-27},{140:155}) .gt. 4,\
                                              1.,0.))
         end if

         do layers = 0,5
         ;  m3/m3
            TWS(case_num,counter)   = Moist(case_num,layers,counter)*Soil_thickness(layers)/sum(Soil_thickness(0:5))
            TWS(case_num+3,counter) = Moist(case_num+3,layers,counter)*Soil_thickness(layers)/sum(Soil_thickness(0:5))
         end do ;varend
      counter = counter + 1
     end do ;month
   end do ;year

end do ;case
; _________________________ Calculate Anomaly ___________________________
var_avg =  dim_avg_n_Wrap(var,2)
Moist_avg =  dim_avg_n_Wrap(Moist,2)
TWS_avg =  dim_avg_n_Wrap(TWS,1)

var_anomaly = new((/6,10,total/),float)
var_anomaly = 0.0
Moist_anomaly = new((/6,7,total/),float)
Moist_anomaly = 0.0
TWS_anomaly   = new((/6,total/),float)
TWS_anomaly   = 0.0

do case_num = 0,5

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
   PDSI(cnt-1812) = avg(f1->sc_PDSI_pm(cnt,{-40:-27},{140:155}))
end do

;_____________________________ PLOT ___________________________________
x = ispan(1,total,1)
print(x)

pic = "GW_time_serial_pft_" + flag ; _anomaly
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

res                    = True
res@gsnMaximize        = True
res@gsnDraw            = False  ; Don't draw plots
res@gsnFrame           = False  ; ����ҳ

res@vpHeightF          = 0.35
res@vpWidthF           = 0.7
res@vpXF               = 0.05
res@vpYF               = 0.8
res@tiMainFontHeightF  = 0.02
res@tiYAxisFontHeightF = 0.02
res@tmXBLabelFontThicknessF = 0.02
res@tmYLLabelFontThicknessF = 0.02
res@tmXBLabelFontHeightF    = 0.02
res@tmYLLabelFontHeightF    = 0.02


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


plots = new(6,graphic)


res@tmXBLabelsOn = True

res@tmXBMode  = "Explicit"
res@tmXBValues=(/6,12,18,24,30,36,42,48,54,60,66,72,78,84,90,96,102,108,114,120/)
res@tmXBLabels=(/"J","D","J","D","J","D","J","D","J","D","J","D","J","D","J","D",\
                 "J","D","J","D","J","D","J","D"/)

res@xyLineThicknesses  = (/  2.0, 2.0, 2.0, 2.0, 2.0, 2.0/)           ; make second line thicker
res@xyLineColor           = "black"          ; change line color
;res@xyExplicitLegendLabel = "PDSI"
res@xyDashPattern         = 0

res@tiMainString   = "(a) PDSI"
res@tiYAxisString  = ""
plots(0) = gsn_csm_xy(wks,x,PDSI,res)
delete(res@xyLineColor)

res@xyLineThicknesses  = (/  2.0, 2.0, 2.0, 2.0, 2.0, 2.0/)          ; make second line thicker
res@xyLineColors       = (/"red","orange","darkgreen", "violet","yellow2","green"/)          ; change line color
res@xyExplicitLegendLabels = (/"on_forest","norchrg_forest","off_forest","on_other","norchrg_other","off_other"/)
res@xyDashPatterns    = (/0,0,0,0,0,0/)


res@tiMainString  = "(b) Soil Moisture"
res@tiYAxisString  ="m~S~3~N~/m~S~3~N~"
;res@tiXAxisString  = "month since Jan 2001"
plots(1) = gsn_csm_xy(wks,x,TWS(:,:),res) ;_anomaly

res@tiMainString   = "(c) Transpiration"
res@tiYAxisString  = "mm/day"
plots(2) = gsn_csm_xy(wks,x,var(:,1,:),res) ;_anomaly

res@tiMainString   = "(d) GPP"
res@tiYAxisString  = "umol/m~S~2~N~/s"
plots(3) = gsn_csm_xy(wks,x,var(:,8,:),res) ;_anomaly

res@tiMainString   = "(e) Latent heat"
res@tiYAxisString  = "W/m~S~2~N~"
plots(4) = gsn_csm_xy(wks,x,var(:,5,:),res) ;_anomaly

res@tiMainString   = "(f) Sensible heat"
res@tiYAxisString  = "W/m~S~2~N~"
plots(5) = gsn_csm_xy(wks,x,var(:,6,:),res) ;_anomaly



gsn_panel(wks,(/plots/),(/3,2/),pres)

delete(plots)
plots = new(7,graphic)

do layer = 0,6
if layer .lt. 6 then
   res@tiMainString = order(layer)+"Moisture of "+(layer+1)+" soil layer"
else
   res@tiMainString = order(layer)+"Aquifer Moisture"
end if
res@tiYAxisString= "m~S~3~N~/m~S~3~N~"
if layer .eq. 6 then
   res@pmLegendWidthF     = 0.15                     ; Resize legend width
   res@pmLegendHeightF    = 0.2                      ; and height
   res@lgAutoManage = False
   res@lgBoxMinorExtentF  =  0.3 ; contral the length of legend line
   res@lgLabelFontHeightF =  0.015; contral the size of legend font
   res@pmLegendDisplayMode = "Always" ; add legend
   res@pmLegendSide        = "Right"
end if
plots(layer) = gsn_csm_xy(wks,x,Moist(:,layer,:),res) ;_anomaly
end do

gsn_panel(wks,(/plots/),(/4,2/),pres)


end
