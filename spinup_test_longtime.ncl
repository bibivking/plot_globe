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
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_spinup/Outputs_181-240/gw_on"\
            /)

var_name= (/\
           "Evap","TVeg","ESoil",\ ;"kg/m^2/s"
           "Qs","Qsb",\      ;"kg/m^2/s"
           "Qle","Qh","Qg",\      ;"W/m^2"
           "Rnet","Albedo",\
           "SWE","CanopInt",\ ;"GPP","NPP",\     ;"umol/m^2/s"
           "RadT","VegT",\   ; "K"
           "SnowT","BaresoilT",\
           "SoilTemp",\
           "SoilMoist",\     ;"m^3/m^3"  Average layer soil moisture
           "SoilMoistIce",\
           "GWMoist",\       ;"mm3/mm3" Aquifer mositure content, Groundwater_off doesn't have
           "WatTable"\       ;"m"
           /)
           ;"LAI" ; the same in Groundwater_on and Groundwater_off

Eq_var_name = (/\
               "wbeq","GWwbeq","wb","GWwb","zq","GWzq","SMP"\
               /)

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers


;______________________________ Data _________________________________
year_s = 1951
year_e = 2010
year_sum = year_e-year_s+1
var_sum  = dimsizes(var_name)
case_sum = dimsizes(case_name)

spinup_time = 4

total = year_sum * 12 * spinup_time
year_sum_total = year_sum * spinup_time

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


var   = new((/16,year_sum_total/),float)
var   = 0.0
Moist = new((/7,year_sum_total/),float)
Moist = 0.0
SoilTemp = new((/6,year_sum_total/),float)
SoilTemp = 0.0
TWS      = new((/year_sum_total/),float)
TWS      = 0.0

;______________________________________ process data _________________________________
do sp_time = 0,spinup_time-1
   print(case_name(sp_time))

   do year = year_s,year_e
      print("year = " +year)
      if ( mod(year,4) .eq. 0) then
         dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
      else
         dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      end if

      filename = case_name(sp_time)+"/cable_out_"+year+".nc"
      f = addfile (filename,"r")

      var@_FillValue   = f->Qle@_FillValue
      Moist@_FillValue = f->Qle@_FillValue
      SoilTemp@_FillValue = f->Qle@_FillValue
      TWS@_FillValue   = f->Qle@_FillValue
      ;print(f->Qle@_FillValue)

      do month = 0,11
         do var_num = 0,15
            if (var_num .le. 4) then
               var(var_num, year-year_s+sp_time*60) = var(var_num, year-year_s+sp_time*60) \
                                                  + avg(f->$var_name(var_num)$(month,{A:B},{C:D}))*3600*24*dom(month)
            else
               var(var_num, year-year_s+sp_time*60) = var(var_num, year-year_s+sp_time*60) \
                                                  + avg(f->$var_name(var_num)$(month,{A:B},{C:D}))/12.
            end if
         end do ;var_num

         do layers = 0,5
            Moist(layers,year-year_s+sp_time*60) = Moist(layers,year-year_s+sp_time*60) \
                                               + avg(f->SoilMoistIce(month,layers,{A:B},{C:D}))/12.
            SoilTemp(layers,year-year_s+sp_time*60) = SoilTemp(layers,year-year_s+sp_time*60) \
                                               + avg(f->SoilTemp(month,layers,{A:B},{C:D}))/12.

         end do
         Moist(6,year-year_s+sp_time*60) = 0.0
     end do ;month

      do layers = 0,5
         ;  m3/m2
         TWS(year-year_s+sp_time*60) = TWS(year-year_s+sp_time*60)+ \
                                     Moist(layers,year-year_s+sp_time*60)*Soil_thickness(layers)
      end do ;varend
      TWS(year-year_s+sp_time*60) = TWS(year-year_s+sp_time*60)+ \
                                    Moist(6,year-year_s+sp_time*60)*avg(f->GWdz)
      print("GWdz "+avg(f->GWdz))
   end do ;year

end do ;spinup_time

;________________________ Check Soil Equilibrium ______________________
wb = new((/2,total,7/),float)
wb = 0.0
SMP  = new((/2,total,7/),float)
SMP  = 0.0
do sp_time = 0,spinup_time-1
do year = year_s,year_e
   print("year = " +year)
   filename = case_name(sp_time)+"/cable_out_"+year+".nc"
   f = addfile (filename,"r")

   time_s = (year-year_s)*12 + sp_time*60*12
   time_e = (year-year_s)*12 + 11 + sp_time*60*12
   wb(0,time_s:time_e,0:5) = dim_avg_n_Wrap(dim_avg_n_Wrap(f->SoilMoist(:,:,{A:B},{C:D}),3),2)  ; Actual
   wb(0,time_s:time_e,6)   = dim_avg_n_Wrap(dim_avg_n_Wrap(f->GWMoist(:,{A:B},{C:D}),2),1)        ; Equilibrium
;   wb(1,time_s:time_e,0:5) = dim_avg_n_Wrap(dim_avg_n_Wrap(f->wbeq(:,:,{A:B},{C:D}),3),2)       ; Actual
;   wb(1,time_s:time_e,6) = dim_avg_n_Wrap(dim_avg_n_Wrap(f->GWwbeq(:,{A:B},{C:D}),2),1)         ; Equilibrium
   SMP(0,time_s:time_e,0:5) = dim_avg_n_Wrap(dim_avg_n_Wrap(f->SMP(:,:,{A:B},{C:D}),3),2)       ; Actual
;   SMP(1,time_s:time_e,0:5) = dim_avg_n_Wrap(dim_avg_n_Wrap(f->zq(:,:,{A:B},{C:D}),3),2)        ; Equilibrium
;   SMP(1,time_s:time_e,6) = dim_avg_n_Wrap(dim_avg_n_Wrap(f->GWzq(:,{A:B},{C:D}),2),1)          ; Equilibrium
end do ;year
end do

printVarSummary(wb)
printVarSummary(SMP)

SMP = SMP*9.8*0.001*0.001
wb@_FillValue   = f->Qle@_FillValue
SMP@_FillValue  = f->Qle@_FillValue


;_____________________________ PLOT ___________________________________
x = ispan(1,year_sum_total,1)
print(x)
pic = "Spinup_test_longtime_"+scale
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

res@xyLineThicknessF    = 1.0        ; make second line thicker
res@pmLegendDisplayMode = "Always"
res@xyLineColor         = "blue"     ; change line color
;res@xyExplicitLegendLabel = "GW"
;res@trYMinF = 0



do var_num = 0,15

res@tiYAxisString  = f->$var_name(var_num)$@long_name+" ("+f->$var_name(var_num)$@units+")"
res@tiXAxisString  = "Year"

if var_num .le. 11 then
;  res@trYMaxF = 1.3*max(var(:,var_num,:))
  plot = gsn_csm_xy(wks,x,var(var_num,:),res)
else
;  res@trYMaxF = 1.3*max(var(:,var_num,:)-273.16)
  plot = gsn_csm_xy(wks,x,var(var_num,:)-273.16,res)
end if

end do


do layers = 0,6
;  res@trYMaxF = 1.3*max(Moist(:,layers,:))
  res@tiYAxisString  = "The soil ice moisture of the "+(layers+1)+ " layer (m~S~3~N~/m~S~3~N~)"
  res@tiXAxisString  = "Year"
  plot = gsn_csm_xy(wks,x,Moist(layers,:),res)
end do
do layers = 0,5
;  res@trYMaxF = 1.3*max(SoilTemp(:,layers,:)-273.16)
  res@tiYAxisString  = "The soil temperature of the "+(layers+1)+ " layer (~S~o~N~C)"
  res@tiXAxisString  = "Year"
  plot = gsn_csm_xy(wks,x,SoilTemp(layers,:)-273.16,res)
end do

;res@trYMaxF = 1.3*max(TWS*1000.)
res@tiYAxisString  = "Total soil water (mm)";"Soil Moisture (m~S~3~N~/m~S~3~N~)"
res@tiXAxisString  = "Year"
plot = gsn_csm_xy(wks,x,TWS*1000.,res)

delete(x)
delete(res)
delete(pic)
delete(wks)

;____________________________ PLOT Soil Equilibrium ____________________________
x = ispan(1,total,1)
print(x)
pic = "Spinup_check_SMP_longtime_"+scale
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow")

res                    = True
res@gsnMaximize        = True
res@vpHeightF          = 0.6
res@vpWidthF           = 0.6
res@vpXF               = 0.3
res@vpYF               = 0.85
res@pmLegendWidthF     = 0.14                     ; Resize legend width
res@pmLegendHeightF    = 0.11                     ; and height

res@pmLegendDisplayMode   = "Always"
res@xyLineThicknessF      = 1.0         ; make second line thicker
res@xyLineColor           = "blue"
;res@xyExplicitLegendLabel = "Spinup"
;res@trYMinF = 0

pres             = True
pres@gsnMaximize = True

plots = new(7,graphic)

do layer = 0,6
   res@tiYAxisString  = "Soil Moisture in layer-"+(layer+1)+" (mm~S~3~N~/mm~S~3~N~)"
   res@tiXAxisString  = "Month since 1951 Jan"
   plots(layer) = gsn_csm_xy(wks,x,wb(0,:,layer),res)
end do
gsn_panel(wks,(/plots/),(/4,2/),pres)

do layer = 0,6
   res@tiYAxisString  = "Soil Matric Potential in layer-"+(layer+1)+" (MPa)"
   res@tiXAxisString  = "Month since 1951 Jan"
   plots(layer) = gsn_csm_xy(wks,x,SMP(0,:,layer),res)
end do

gsn_panel(wks,(/plots/),(/4,2/),pres)

end
