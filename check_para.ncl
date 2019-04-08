;*****************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_norecharge_unify_para/Outputs/gw_on"\
             /)

var_name = (/\
             "bch","clay","sand","silt","ssat","sfc","swilt","hyds","sucs",\
	           "css","rhosoil"\
             /)
             ;,"albsoil","zse","froot" ;3D
             ;"iveg","isoil","vbeta","slope","elev","slope_std","GWdz","Qhmax",\ ;2D
             ;"QhmaxEfold","QhmaxEfold","SatFracmax","HKefold","HKdepth",\
var_name1= (/\
             "bch_vec","clay_vec","sand_vec","silt_vec","ssat_vec","sfc_vec","swilt_vec",\
             "hyds_vec","sucs_vec","css_vec","rhosoil_vec"\
             /)
;______________________________ Data _________________________________
year_s = 2001
year_e = 2010
year_sum = year_e-year_s+1
var_sum  = dimsizes(var_name)
year = 1952

pic = "var_check"
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

res                    = True
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                            ;��ɫ
res@tmBorderThicknessF = 3.0
res@lbLabelBarOn       = True
res@cnLinesOn          = False

filename1 = case_name(0)+"/cable_out_"+year+".nc"
filename2 = case_name(1)+"/cable_out_"+year+".nc"
filename3 = case_name(2)+"/cable_out_"+year+".nc"
filename4 = "/short/w35/mm3972/cable/src/CABLE-AUX/offline/gridinfo_mmy_MD_elev_orig_std_avg-sand_mask.nc"

f1 = addfile (filename1,"r")
f2 = addfile (filename2,"r")
f3 = addfile (filename3,"r")
f4 = addfile (filename4,"r")

do var_num = 0,10

var1 = dim_avg_n_Wrap(f1->$var_name(var_num)$,0)
var2 = dim_avg_n_Wrap(f2->$var_name(var_num)$,0)
var3 = dim_avg_n_Wrap(f3->$var_name(var_num)$,0)

if var_num .eq. 7 then
var4 = dim_avg_n_Wrap(f4->$var_name1(var_num)$,0)/1000.
else
var4 = dim_avg_n_Wrap(f4->$var_name1(var_num)$,0)
end if

 res@tiMainString = var_name(var_num)+" gw_on"
 plot = gsn_csm_contour(wks,var1,res)
 res@tiMainString = "gw_off"
 plot = gsn_csm_contour(wks,var2,res)
 res@tiMainString = "gw_on_no_recharge"
 plot = gsn_csm_contour(wks,var3,res)
 res@tiMainString = "gridinfo"
 plot = gsn_csm_contour(wks,var4,res)

 res@tiMainString = "gw_on - gridinfo"
 plot = gsn_csm_contour(wks,(abs(var1)-abs(var4)),res)
 res@tiMainString = "gw_off -gridinfo"
 plot = gsn_csm_contour(wks,(abs(var2)-abs(var4)),res)
 res@tiMainString = "gw_on_no_recharge -gridinfo"
 plot = gsn_csm_contour(wks,(abs(var3)-abs(var4)),res)

 res@tiMainString = "gw_on - gw_off"
 plot = gsn_csm_contour(wks,(abs(var1)-abs(var2)),res)
 res@tiMainString = "gw_on_no_recharge - gw_on"
 plot = gsn_csm_contour(wks,(abs(var3)-abs(var1)),res)

end do
end
