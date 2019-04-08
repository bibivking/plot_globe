;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on/Outputs/gw_on/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off/Outputs/gw_off/"\
              /)

var_name= (/\
           "Evap","TVeg","ESoil",\ ;"kg/m^2/s"
           "Qs","Qsb",\      ;"kg/m^2/s"
           "Qle","Qh",\      ;"W/m^2"
           "GPP","NPP",\     ;"umol/m^2/s"
           "SoilMoist",\     ;"m^3/m^3"  Average layer soil moisture
           "GWMoist",\       ;"mm3/mm3" Aquifer mositure content, Groundwater_off doesn't have
           "WatTable"\       ;"m"
           /)
           ;"LAI" ; the same in Groundwater_on and Groundwater_off

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872, 20./) ; 6 soil layers and 1 aquifer


;______________________________ Data _________________________________
year_s = 1951
year_e = 2010
year_sum = year_e-year_s+1
var_sum  = dimsizes(var_name)
total = year_sum * 12

var = new((/2,9,year_sum/),float)
var = 0.0
Moist = new((/2,7,year_sum/),float)
Moist = 0.0
TWS   = new((/2,year_sum/),float)
TWS   = 0.0

;var!0 = "case"
;var!1 = "var"
;var!2 = "time"

;Moist!0 = "case"
;Moist!1 = "layer"
;Moist!2 = "time"

;TWS!0   = "case"
;TWS!1   = "time"

;____________________________ produce mask_data _______________________________
filename = "/short/w35/mm3972/cable/src/CABLE-AUX/offline/gridinfo_mmy_MD_elev_orig_landmask.nc"
f = addfile (filename,"r")

mask_data = f->landsea ;-9999

;do case_num = 0,1
;   print(case_name(case_num))

;   do year = year_s,year_e
;      print("year = " +year)
;      if ( mod(year,4) .eq. 0) then
;         dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
;      else
;         dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
;      end if

;      filename = case_name(case_num)+"/cable_out_"+year+".nc"
;      f = addfile (filename,"r")

      ;do month = 0,11
      ;   do var_num = 0,8
      ;      temp = f->$var_name(var_num)$(month,:,:)
            ;if (any(isnan_ieee(temp))) then
            ;   replace_ieeenan (temp, temp@_FillValue, 0)
            ;end if
            ;mask_data = where(.not. ismissing(temp) .and. abs(temp) .lt. 10000., mask_data , mask_data@_FillValue)
            ;delete(temp)
      ;   end do ;var_num

      ;   do layers = 0,5
      ;      temp = f->SoilMoist(month,layers,:,:)
      ;      if (any(isnan_ieee(temp))) then
      ;         replace_ieeenan (temp, temp@_FillValue, 0)
      ;      end if
      ;      mask_data = where(.not. ismissing(temp) .and. abs(temp) .lt. 10000., mask_data , mask_data@_FillValue)
      ;      delete(temp)
      ;   end do
      ;   if (case_num .ne. 1)
      ;      temp = f->GWMoist(month,:,:)
      ;      if (any(isnan_ieee(temp))) then
      ;         replace_ieeenan (temp, temp@_FillValue, 0)
      ;      end if
      ;      mask_data = where(.not. ismissing(temp) .and. abs(temp) .lt. 10000., mask_data , mask_data@_FillValue)
      ;      delete(temp)
      ;  end if
     ;end do ;month
;   end do ;year
;end do ;case

;______________________________________ process data _________________________________
do case_num = 0,1
   print(case_name(case_num))

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
         do var_num = 0,8
            if (var_num .le. 4) then
               var(case_num, var_num, year-year_s) = var(case_num, var_num, year-year_s) \
                                           + avg(f->$var_name(var_num)$(month,:,:))*3600*24*dom(month)
            else
               var(case_num, var_num, year-year_s) = var(case_num, var_num, year-year_s) \
                                                  + avg(f->$var_name(var_num)$(month,:,:))/12.
            end if
         end do ;var_num

         do layers = 0,5
            Moist(case_num,layers,year-year_s) = Moist(case_num,layers,year-year_s) \
                                               + avg(f->SoilMoist(month,layers,:,:))/12.
         end do
         if (case_num .ne. 1)
            Moist(case_num,6,year-year_s)      = Moist(case_num,6,year-year_s)   \
                                               + avg(f->GWMoist(month,:,:))/12.
         else
            Moist(case_num,6,year-year_s)      = 0.0
         end if
     end do ;month

      if ( case_num .ne. 1) then
         varend = 5
      else
         varend = 5
      end if

      do layers = 0,varend
         ;  m3/m2
         TWS(case_num,year-year_s) = TWS(case_num,year-year_s)+ \
                                     Moist(case_num,layers,year-year_s)*Soil_thickness(layers)/sum(Soil_thickness(0:5))
      end do ;varend

   end do ;year

end do ;case

;_____________________________ PLOT ___________________________________
x = ispan(year_s,year_e,1)
print(x)
pic = "GW_Spinup_check"
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

res                    = True
res@gsnMaximize        = False
res@vpHeightF          = 0.6
res@vpWidthF           = 0.6
res@vpXF               = 0.3
res@vpYF               = 0.85
res@pmLegendWidthF     = 0.14                     ; Resize legend width
res@pmLegendHeightF    = 0.11                     ; and height

res@xyLineThicknesses  = (/  2.0, 2.0, 2.0/)          ; make second line thicker
res@pmLegendDisplayMode = "Always"
res@xyLineColors       = (/"blue","red","green"/)          ; change line color
res@xyExplicitLegendLabels = (/"GW","NoGW","MD"/)
res@trYMinF = 0



do var_num = 0,8
res@trYMaxF = 1.3*max(var(:,var_num,:))
res@tiYAxisString  = f->$var_name(var_num)$@long_name+" ("+f->$var_name(var_num)$@units+")"
res@tiXAxisString  = "Year"
plot = gsn_csm_xy(wks,x,var(:,var_num,:),res)
end do


do layers = 0,6
  res@trYMaxF = 1.3*max(Moist(:,layers,:))
  res@tiYAxisString  = "The soil moisture of the "+(layers+1)+ " layer (m~S~3~N~/m~S~3~N~)"
  res@tiXAxisString  = "Year"
  plot = gsn_csm_xy(wks,x,Moist(:,layers,:),res)
end do


res@trYMaxF = 1.3*max(TWS)
res@tiYAxisString  = "Soil Moisture (m~S~3~N~/m~S~3~N~)"
res@tiXAxisString  = "Year"
plot = gsn_csm_xy(wks,x,TWS,res)

end
