;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-standard/gw_on/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs-standard/gw_on/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs-standard/gw_off/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on/Outputs/gw_on/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off/Outputs/gw_off/"\
             /)
             
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


;______________________________ Data _________________________________
year_s = 2001
year_e = 2010
year_sum = year_e-year_s+1
var_sum  = dimsizes(var_name)
total = year_sum * 12

var = new((/5,10,360/),float)
var = 0.0
Moist = new((/5,7,360/),float)
Moist = 0.0
TWS   = new((/5,360/),float)
TWS   = 0.0


;______________________________________ process data _________________________________
do case_num = 0,4
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

         do var_num = 0,9
            if (var_num .le. 4) then
               var(case_num, var_num,:) = var(case_num, var_num,:) + dim_avg_n_Wrap(dim_avg_n_Wrap(f->$var_name(var_num)$,0),1)*3600*24
            else if (var_num .le. 8) then
               var(case_num, var_num,:) = var(case_num, var_num,:) + dim_avg_n_Wrap(dim_avg_n_Wrap(f->$var_name(var_num)$,0),1)
            end if
            end if
         end do ;var_num

         do layers = 0,5
            Moist(case_num,layers,:) = Moist(case_num,layers,:) + dim_avg_n_Wrap(dim_avg_n_Wrap(f->SoilMoist(:,layers,:,:),0),1) 
         end do
         if (case_num .eq. 1)
            Moist(case_num,6,:)      = 0.0
         else
            Moist(case_num,6,:)      = Moist(case_num,layers,:) + dim_avg_n_Wrap(dim_avg_n_Wrap(f->GWMoist,0),1)  
         end if

         do layers = 0,5
         ;  m3/m3
            TWS(case_num,:) = Moist(case_num,layers,:)*Soil_thickness(layers)/sum(Soil_thickness(0:5))
         end do ;varend
   end do ;year

end do ;case

var   = var/year_sum
Moist = Moist/year_sum
TWS   = TWS/year_sum

;_____________________________ PLOT ___________________________________
pic = "GW_zonal_avg";
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"rainbow")

res                    = True
res@gsnMaximize        = False
res@gsnDraw            = False  ; Don't draw plots
res@gsnFrame           = False  ; ?????

res@vpHeightF          = 0.7
res@vpWidthF           = 0.3
;res@vpXF               = 0.0
;res@vpYF               = 0.8
res@pmLegendWidthF     = 0.15                      ; Resize legend width
res@pmLegendHeightF    = 0.15                     ; and height
res@pmLegendSide       = "Bottom"
res@tmXBLabelFontThicknessF = 0.015
res@tmYLLabelFontThicknessF = 0.015
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015

pres                    = True                                      ; ???? panel resources.
pres@gsnMaximize        = True
pres@gsnPanelLabelBar   = False                    ; ???? panel labelbar.
pres@lbBoxLinesOn       = False                       ;???lbar box ????
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
pres@lbLabelBarOn       = False ; True                    ; ??????plot??labelbar
pres@lbTitleFontHeightF = 0.5
pres@lbLabelFontHeightF = 0.5
pres@lbAutoManage       = True
pres@lbTitleDirection   = "Across"             ; title direction
pres@lbTitlePosition    = "Bottom";Right"


plots = new(5,graphic)

res@xyLineThicknesses  = (/  2.0, 2.0, 2.0, 2.0, 2.0/)          ; make second line thicker
res@xyLineColors       = (/"blue","red","darkgreen","cyan","green"/)          ; change line color
res@xyExplicitLegendLabels = (/"GW","Norecharge","NoGW","GW-CTL","NoGW-CTL"/)
res@xyDashPatterns    = (/0,0,0,0,0/)

;res@trYMinF = 0
;res@trYMaxF = 1.3*max(var(:,0,:))
res@tmYLLabelsOn = True
res@tiXAxisString  = "Evap (mm/day)"
plots(0) = gsn_csm_xy(wks,var(:,0,:),f->y,res) ;_anomaly

res@tmYLLabelsOn = False
res@tiXAxisString  = "TVeg (mm/day)"
plots(1) = gsn_csm_xy(wks,var(:,1,:),f->y,res) ;_anomaly

res@tiXAxisString  = "ESoil (mm/day)"
plots(2) = gsn_csm_xy(wks,var(:,2,:),f->y,res) ;_anomaly

res@tiXAxisString  = "Runoff (mm/day)"
plots(3) = gsn_csm_xy(wks,var(:,3,:)+var(:,4,:),f->y,res) ;_anomaly

res@tiXAxisString  = "Soil Moistur (m~S~3~N~/m~S~3~N~)"
plots(4) = gsn_csm_xy(wks,TWS(:,:),f->y,res) ;_anomaly

gsn_panel(wks,(/plots/),(/1,5/),pres)
res@tmYLLabelsOn = True
res@tiXAxisString  = "Qle (W/m~S~2~N~)"
plots(0) = gsn_csm_xy(wks,var(:,5,:),f->y,res) ;_anomaly

res@tmYLLabelsOn = False
res@tiXAxisString  = "Qh (W/m~S~2~N~)"
plots(1) = gsn_csm_xy(wks,var(:,6,:),f->y,res) ;_anomaly

res@tiXAxisString  = "WatTable (m)"
plots(2) = gsn_csm_xy(wks,var(:,7,:),f->y,res) ;_anomaly

res@tiXAxisString  = "GPP (umol/m~S~2~N~/s)"
plots(3) = gsn_csm_xy(wks,var(:,8,:),f->y,res) ;_anomaly

gsn_panel(wks,(/plots(0:3)/),(/1,4/),pres)


  
;  plot=gsn_csm_xy(wks,zave(0,:),ts&lat,False)       ; plot zonal ave

end


