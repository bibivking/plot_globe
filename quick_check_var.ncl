;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-standard/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs-standard/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-Haverd2013/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs-Haverd2013/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs-standard/gw_on"\
             /)

var_name = "Fwsoil"
var_name1= "sand_vec"
mon_name= (/\
              "Jan","Feb","Mar","Apr","May",\
              "Jun","Jul","Aug","Sep","Oct",\
              "Nov","Dec"\
              /)
; Note that: Cannot use the name of "month_name", while I have no idea about the reason...
;______________________________ Data _________________________________
year_s = 2001
year_e = 2010
year_sum = year_e-year_s+1
var_sum  = dimsizes(var_name)
;year = 1951

filename1 = case_name(0)+"/cable_out_"+year_e+".nc"
filename2 = case_name(1)+"/cable_out_"+year_e+".nc"
filename3 = case_name(2)+"/cable_out_"+year_e+".nc"
filename4 = case_name(3)+"/cable_out_"+year_e+".nc"
;filename3 = case_name(2)+"/cable_out_"+year+".nc"

;filename1 = case_name(0)+"/restart_1951.nc"
;filename2 = case_name(1)+"/restart_1951.nc"
;filename3 = case_name(2)+"/restart_1951.nc"

;filename4 = "/short/w35/mm3972/cable/src/CABLE-AUX/offline/gridinfo_mmy_MD_elev_orig_std_avg-sand_mask.nc"
;filename2 = "/short/w35/mm3972/cable/src/CABLE-AUX/offline/gridinfo_mmy_MD_elev_orig_std_mask.nc"


f1 = addfile (filename1,"r")
f2 = addfile (filename2,"r")
f3 = addfile (filename3,"r")
f4 = addfile (filename4,"r")

var1 = f1->$var_name$
var2 = f2->$var_name$
var3 = f3->$var_name$
var4 = f4->$var_name$

;print(sum(var1-var2))
;print(sum(var1-var3))
;print(sum(var2-var3))

;var1 = f1->$var_name$
;var2 = f2->$var_name$
;print(sum(var1-var2))

;var1 = dim_avg_n_Wrap(f1->$var_name$,0)
;var2 = dim_avg_n_Wrap(f2->$var_name$,0)
;var3 = dim_avg_n_Wrap(f3->$var_name$,0)
;var4 = dim_avg_n_Wrap(f4->$var_name1$,0)

;var4 = dim_avg_n_Wrap(where( (f1->sand+f1->clay+f1->silt) .gt. 1. , 1., 0.), 0)
;var5 = dim_avg_n_Wrap(where( (f2->sand+f2->clay+f2->silt) .gt. 1. , 1., 0.), 0)
;var6 = where( (f3->sand+f3->clay+f3->silt) .gt. 1., 1., 0.)

pic = "fwsoil_2010_on-off"
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed");"ViBlGrWhYeOrRe")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

res                    = True
res@gsnDraw            = False  ; Don't draw plots
res@gsnFrame           = False  ; ����ҳ
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                            ;��ɫ
res@tmBorderThicknessF = 3.0
res@lbLabelBarOn       = True
res@cnLinesOn          = False

res@cnExplicitLabelBarLabelsOn = True   
res@cnLevelSelectionMode = "ExplicitLevels"

pres = True
plots = new(6,graphic)

do month = 0,11
res@cnFillPalette = "WhiteBlueGreenYellowRed"
levels = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)
res@cnLevels = levels
res@tiMainString = "on_std"
plots(0) = gsn_csm_contour(wks,var1(month,:,:),res)
res@tiMainString = "off_std"
plots(1) = gsn_csm_contour(wks,var2(month,:,:),res)
delete(res@cnLevels)
delete(levels)

levels = (/-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)
res@cnLevels = levels
res@cnFillPalette = "ViBlGrWhYeOrRe"
res@tiMainString = "difference_std (on - off)"
plots(2) = gsn_csm_contour(wks,var1(month,:,:)-var2(month,:,:),res)
delete(levels)
delete(res@cnLevels)

res@cnFillPalette = "WhiteBlueGreenYellowRed"
levels = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)
res@cnLevels = levels
res@tiMainString = "on_Hvrd"
plots(3) = gsn_csm_contour(wks,var3(month,:,:),res)
res@tiMainString = "off_Hvrd"
plots(4) = gsn_csm_contour(wks,var4(month,:,:),res)
delete(res@cnLevels)
delete(levels)

levels = (/-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)
res@cnLevels = levels
res@cnFillPalette = "ViBlGrWhYeOrRe"
res@tiMainString = "difference_Hvrd (on - off)"
plots(5) = gsn_csm_contour(wks,var3(month,:,:)-var4(month,:,:),res)
delete(levels)
delete(res@cnLevels)

pres@txString = mon_name(month)
gsn_panel(wks,(/plots/),(/2,3/),pres)

end do
;levels = (/0,1,2,3,4,5/)
;res@cnLevels = levels

;forest = var1
;forest = 0
;forest = where(var1 .ge. 1 .and. var1 .le. 4, var1, 0)
;res@tiMainString = "forest"
;plot = gsn_csm_contour(wks,forest,res)

;plot = gsn_csm_contour(wks,forest({-40:-27},{140:155}),res)
;plot = gsn_csm_contour(wks,var1, res)

;res@cnExplicitLabelBarLabelsOn = True   ;?
;res@cnLevelSelectionMode = "ExplicitLevels"

;levels = (/-35,-30,-25,-20,-15,-10,-5,-4,-3,-2,-1,1,2,3,4,5,10,15,20,25,30,35/)
;res@cnLevels = levels

;layer = 5
; res@tiMainString = "gw_on"
; plot = gsn_csm_contour(wks,var1(layer,:,:),res)
; res@tiMainString = "gw_off"
; plot = gsn_csm_contour(wks,var2(layer,:,:),res)
; res@tiMainString = "gw_on_no_recharge"
; plot = gsn_csm_contour(wks,var3(layer,:,:),res)
; res@tiMainString = "gridinfo"
; plot = gsn_csm_contour(wks,var4,res)
;
; plot = gsn_csm_contour(wks,var4,res)
; plot = gsn_csm_contour(wks,var5,res)
; plot = gsn_csm_contour(wks,var6,res)

; res@tiMainString = "gw_on - gridinfo"
; plot = gsn_csm_contour(wks,(abs(var1)-abs(var4)),res)
; res@tiMainString = "gw_off -gridinfo"
; plot = gsn_csm_contour(wks,(abs(var2)-abs(var4)),res)
; res@tiMainString = "gw_off -gridinfo"
; plot = gsn_csm_contour(wks,(abs(var2)-abs(var4)),res)
; res@tiMainString = "gw_on_no_recharge -gridinfo"
; plot = gsn_csm_contour(wks,(abs(var3)-abs(var4)),res)


; res@tiMainString = "gw_on - gw_off"
; plot = gsn_csm_contour(wks,(var1(layer,:,:)-var2(layer,:,:)),res)
; res@tiMainString = "gw_on_no_recharge - gw_on"
; plot = gsn_csm_contour(wks,(var3(layer,:,:)-var1(layer,:,:)),res)


end
