;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

begin

case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_debug/Outputs/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on/Outputs/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_debug/Outputs/gw_on"\
              /)
              
case_name_brf = (/\
              "GW",\
              "CLASS",\
              "CERES"\
               /)

var_name= (/\
           "Rnet","Qh","Qg","Qle",\
           "Fwsoil",\
           "WatTable",\       ;"m"
           "Qrecharge",\
           "Evap","TVeg","ESoil",\
           "Qs","Qsb",\
           "SoilTemp",\
           "SoilMoist","GWMoist",\    ;"m^3/m^3"  Average layer soil moisture
           "wbeq","GWwbeq",\
           "zq","GWzq","SMP",\
           "SoilMoistIce",\
           "RadT","VegT",\   ; "K"
           "SnowT","BaresoilT"\
           /)
           ;"LAI" ; the same in Groundwater_on and Groundwater_off

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers

;______________________________ Data _________________________________
year_s   = 1991
year_e   = 2010
year_sum = year_e-year_s+1
var_sum  = dimsizes(var_name)
case_sum = dimsizes(case_name)

scale = "Global" ;"Global" ; "50S-N" ; "SE-AU" ; "AU"
time_scale = "annual" ;"DJF";"MAM";"JJA";"SON" ;"annual"

if time_scale .eq. "DJF" then
   month_s = 0
   month_e = 1
else if time_scale .eq. "MAM" then
   month_s = 2
   month_e = 4
else if time_scale .eq. "JJA" then
   month_s = 5
   month_e = 7
else if time_scale .eq. "SON" then
   month_s = 8
   month_e = 10
else 
   month_s = 0
   month_e = 11
end if 
end if 
end if
end if


CABLE     = new((/12,360,720/),float)
CABLE     = 0.0
SoilMoist = new((/2,7,360,720/),float)
SoilMoist = 0.0
SMP       = new((/2,6,360,720/),float)
SMP       = 0.0
SoilTemp  = new((/6,360,720/),float)
SoilTemp  = 0.0

CABLE_dif     = new((/12,360,720/),float)
CABLE_dif     = 0.0
SoilMoist_dif = new((/7,360,720/),float)
SoilMoist_dif = 0.0
SoilTemp_dif  = new((/6,360,720/),float)
SoilTemp_dif  = 0.0

filename0 = case_name(2)+"/cable_out_2000.nc"
f0 = addfile (filename0,"r")
lat = f0->Qle&y
lon = f0->Qle&x

CABLE!0 = "var"
CABLE!1 = "lat"
CABLE!2 = "lon"
CABLE&lat = lat
CABLE&lon = lon
CABLE@_FillValue   = f0->Qle@_FillValue

SoilMoist!0 = "actualorequilibrium"
SoilMoist!1 = "layer"
SoilMoist!2 = "lat"
SoilMoist!3 = "lon"
SoilMoist&lat = lat
SoilMoist&lon = lon
SoilMoist@_FillValue   = f0->Qle@_FillValue

SMP!0 = "actualorequilibrium"
SMP!1 = "layer"
SMP!2 = "lat"
SMP!3 = "lon"
SMP&lat = lat
SMP&lon = lon
SMP@_FillValue   = f0->Qle@_FillValue

SoilTemp!0 = "layer"
SoilTemp!1 = "lat"
SoilTemp!2 = "lon"
SoilTemp&lat = lat
SoilTemp&lon = lon
SoilTemp@_FillValue   = f0->Qle@_FillValue

CABLE_dif!0 = "var"
CABLE_dif!1 = "lat"
CABLE_dif!2 = "lon"
CABLE_dif&lat = lat
CABLE_dif&lon = lon
CABLE_dif@_FillValue   = f0->Qle@_FillValue

SoilMoist_dif!0 = "layer"
SoilMoist_dif!1 = "lat"
SoilMoist_dif!2 = "lon"
SoilMoist_dif&lat = lat
SoilMoist_dif&lon = lon
SoilMoist_dif@_FillValue   = f0->Qle@_FillValue

SoilTemp_dif!0 = "layer"
SoilTemp_dif!1 = "lat"
SoilTemp_dif!2 = "lon"
SoilTemp_dif&lat = lat
SoilTemp_dif&lon = lon
SoilTemp_dif@_FillValue   = f0->Qle@_FillValue

do year = year_s,year_e

    print("year = " +year)
    if ( mod(year,4) .eq. 0) then
       dom = (/31,29,31,30,31,30,31,31,30,31,30,31/)
    else
       dom = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    end if

    filename  = case_name(0)+"/cable_out_"+year+".nc"
    f  = addfile (filename,"r")

    filename1 = case_name(1)+"/cable_out_"+year+".nc"
    f1 = addfile (filename1,"r")

    filename2 = case_name(2)+"/cable_out_"+year+".nc"
    f2 = addfile (filename2,"r")


  if time_scale .eq. "DJF" then
    do var_num = 0,5
       CABLE(var_num,:,:) = CABLE(var_num,:,:) + (f->$var_name(var_num)$(11,:,:)+f->$var_name(var_num)$(0,:,:)+f->$var_name(var_num)$(1,:,:))/3.
       CABLE_dif(var_num,:,:) = CABLE_dif(var_num,:,:) \
                       + ( f1->$var_name(var_num)$(11,:,:)+f1->$var_name(var_num)$(0,:,:)+f1->$var_name(var_num)$(1,:,:) )/3.\
                       - ( f->$var_name(var_num)$(11,:,:) +f->$var_name(var_num)$(0,:,:) +f->$var_name(var_num)$(1,:,:)  )/3.
    end do ;var_num
    ; Qrecharge
       CABLE(6,:,:) = CABLE(6,:,:) + f->$var_name(6)$(11,:,:)*3600.*24.*dom(11) \
                                   + f->$var_name(6)$(0,:,:)*3600.*24.*dom(0)   \
                                   + f->$var_name(6)$(1,:,:)*3600.*24.*dom(1)
       CABLE_dif(6,:,:) = CABLE_dif(6,:,:) + (f1->$var_name(6)$(11,:,:) - f->$var_name(6)$(11,:,:))*3600.*24.*dom(11) \
                                           + (f1->$var_name(6)$(0,:,:) - f->$var_name(6)$(0,:,:))*3600.*24.*dom(0)    \  
                                           + (f1->$var_name(6)$(1,:,:) - f->$var_name(6)$(1,:,:))*3600.*24.*dom(1)
    ; Evap             
       CABLE(7,:,:) = CABLE(7,:,:) + f->$var_name(7)$(11,:,:)*3600.*24.*dom(11) \
                                   + f->$var_name(7)$(0,:,:)*3600.*24.*dom(0)   \
                                   + f->$var_name(7)$(1,:,:)*3600.*24.*dom(1)
       CABLE_dif(7,:,:) = CABLE_dif(7,:,:) + (f1->$var_name(7)$(11,:,:) - f->$var_name(7)$(11,:,:))*3600.*24.*dom(11) \
                                           + (f1->$var_name(7)$(0,:,:) - f->$var_name(7)$(0,:,:))*3600.*24.*dom(0)    \
                                           + (f1->$var_name(7)$(1,:,:) - f->$var_name(7)$(1,:,:))*3600.*24.*dom(1)
    ; TVeg
       CABLE(8,:,:) = CABLE(8,:,:) + f->$var_name(8)$(11,:,:)*3600.*24.*dom(11) \
                                   + f->$var_name(8)$(0,:,:)*3600.*24.*dom(0)   \
                                   + f->$var_name(8)$(1,:,:)*3600.*24.*dom(1)
       CABLE_dif(8,:,:) = CABLE_dif(8,:,:) + (f1->$var_name(8)$(11,:,:) - f->$var_name(8)$(11,:,:))*3600.*24.*dom(11) \
                                           + (f1->$var_name(8)$(0,:,:) - f->$var_name(8)$(0,:,:))*3600.*24.*dom(0)    \
                                           + (f1->$var_name(8)$(1,:,:) - f->$var_name(8)$(1,:,:))*3600.*24.*dom(1)
    ; ESoil
       CABLE(9,:,:) = CABLE(9,:,:) + f->$var_name(9)$(11,:,:)*3600.*24.*dom(11) \
                                   + f->$var_name(9)$(0,:,:)*3600.*24.*dom(0)   \
                                   + f->$var_name(9)$(1,:,:)*3600.*24.*dom(1)
       CABLE_dif(9,:,:) = CABLE_dif(9,:,:) + (f1->$var_name(9)$(11,:,:) - f->$var_name(9)$(11,:,:))*3600.*24.*dom(11) \
                                           + (f1->$var_name(9)$(0,:,:) - f->$var_name(9)$(0,:,:))*3600.*24.*dom(0)    \
                                           + (f1->$var_name(9)$(1,:,:) - f->$var_name(9)$(1,:,:))*3600.*24.*dom(1)
    ; Qs
       CABLE(10,:,:) = CABLE(10,:,:) + f->$var_name(10)$(11,:,:)*3600.*24.*dom(11) \
                                     + f->$var_name(10)$(0,:,:)*3600.*24.*dom(0)   \
                                     + f->$var_name(10)$(1,:,:)*3600.*24.*dom(1)
       CABLE_dif(10,:,:) = CABLE_dif(10,:,:) + (f1->$var_name(10)$(11,:,:) - f->$var_name(10)$(11,:,:))*3600.*24.*dom(11) \
                                             + (f1->$var_name(10)$(0,:,:) - f->$var_name(10)$(0,:,:))*3600.*24.*dom(0)    \
                                             + (f1->$var_name(10)$(1,:,:) - f->$var_name(10)$(1,:,:))*3600.*24.*dom(1)
    ; Qsb
       CABLE(11,:,:) = CABLE(11,:,:) + f->$var_name(11)$(11,:,:)*3600.*24.*dom(11) \
                                     + f->$var_name(11)$(0,:,:)*3600.*24.*dom(0)   \
                                     + f->$var_name(11)$(1,:,:)*3600.*24.*dom(1)
       CABLE_dif(11,:,:) = CABLE_dif(11,:,:) + (f1->$var_name(11)$(11,:,:) - f->$var_name(11)$(11,:,:))*3600.*24.*dom(11) \
                                             + (f1->$var_name(11)$(0,:,:) - f->$var_name(11)$(0,:,:))*3600.*24.*dom(0)    \
                                             + (f1->$var_name(11)$(1,:,:) - f->$var_name(11)$(1,:,:))*3600.*24.*dom(1)
    SoilTemp = SoilTemp     + (f->SoilTemp(11,:,:,:)+f->SoilTemp(0,:,:,:)+f->SoilTemp(1,:,:,:))/3.

    SoilMoist(0,0:5,:,:) = SoilMoist(0,0:5,:,:) + ( f->SoilMoist(11,:,:,:)\
                                                +   f->SoilMoist(0,:,:,:) \
                                                +   f->SoilMoist(1,:,:,:) ) /3. ; Actual
    SoilMoist(1,0:5,:,:) = SoilMoist(1,0:5,:,:) + ( f1->SoilMoist(11,:,:,:) - f->SoilMoist(11,:,:,:)    \
                                                +   f1->SoilMoist(0,:,:,:)  - f->SoilMoist(0,:,:,:)     \
                                                +   f1->SoilMoist(1,:,:,:)  - f->SoilMoist(1,:,:,:)     ) /3.      ; Equilibrium
    SoilMoist(0,6,:,:)   = SoilMoist(0,6,:,:)   + ( f->GWMoist(11,:,:)    \
                                                +   f->GWMoist(0,:,:)     \
                                                +   f->GWMoist(1,:,:)     ) /3. ; Actual
    SoilMoist(1,6,:,:)   = SoilMoist(1,6,:,:)   + ( f1->GWMoist(11,:,:) - f->GWMoist(11,:,:)    \
                                                +   f1->GWMoist(0,:,:)  - f->GWMoist(0,:,:)     \
                                                +   f1->GWMoist(1,:,:)  - f->GWMoist(1,:,:)     ) /3.      ; Equilibrium
    SMP(0,:,:,:)       = SMP(0,:,:,:)       + ( f->SMP(11,:,:,:)      \
                                            +   f->SMP(0,:,:,:)       \
                                            +   f->SMP(1,:,:,:)       ) /3. ; Actual
    SMP(1,:,:,:)       = SMP(1,:,:,:)       + ( f1->SMP(11,:,:,:)   - f->SMP(11,:,:,:)   \
                                            +   f1->SMP(0,:,:,:)    - f->SMP(0,:,:,:)    \
                                            +   f1->SMP(1,:,:,:)    - f->SMP(1,:,:,:)    ) /3.     ; Equilibrium


    SoilTemp_dif  = SoilTemp_dif + ( f1->SoilTemp(11,:,:,:) + f1->SoilTemp(0,:,:,:) + f1->SoilTemp(1,:,:,:) \
                                 -   f->SoilTemp(11,:,:,:)  - f->SoilTemp(0,:,:,:)  - f->SoilTemp(1,:,:,:)  ) /3.
    SoilMoist_dif(0:5,:,:) = SoilMoist_dif(0:5,:,:) + (  f1->SoilMoist(11,:,:,:) + f1->SoilMoist(0,:,:,:) + f1->SoilMoist(1,:,:,:)\
                                                    -    f->SoilMoist(11,:,:,:)  - f->SoilMoist(0,:,:,:)  - f->SoilMoist(1,:,:,:) ) /3.
    SoilMoist_dif(6,:,:)   = SoilMoist_dif(6,:,:) + ( f1->GWMoist(11,:,:) + f1->GWMoist(0,:,:) + f1->GWMoist(1,:,:) \
                                                  -   f->GWMoist(11,:,:)  - f->GWMoist(0,:,:)  - f->GWMoist(1,:,:)  ) /3.
  else 
    do var_num = 0,5
       CABLE(var_num,:,:) = CABLE(var_num,:,:) + dim_avg_n_Wrap(f->$var_name(var_num)$(month_s:month_e,:,:),0)
       CABLE_dif(var_num,:,:) = CABLE_dif(var_num,:,:) + \
                 (dim_avg_n_Wrap(f1->$var_name(var_num)$(month_s:month_e,:,:),0) - dim_avg_n_Wrap(f->$var_name(var_num)$(month_s:month_e,:,:),0))
    end do ;var_num
    do month = month_s,month_e
    ; Qrecharge
       CABLE(6,:,:) = CABLE(6,:,:) + f->$var_name(6)$(month,:,:)*3600.*24.*dom(month)
       CABLE_dif(6,:,:) = CABLE_dif(6,:,:) + (f1->$var_name(6)$(month,:,:) - f->$var_name(6)$(month,:,:))*3600.*24.*dom(month)
    ; Evap
       CABLE(7,:,:) = CABLE(7,:,:) + f->$var_name(7)$(month,:,:)*3600.*24.*dom(month)
       CABLE_dif(7,:,:) = CABLE_dif(7,:,:) + (f1->$var_name(7)$(month,:,:) - f->$var_name(7)$(month,:,:))*3600.*24.*dom(month)
    ; TVeg
       CABLE(8,:,:) = CABLE(8,:,:) + f->$var_name(8)$(month,:,:)*3600.*24.*dom(month)
       CABLE_dif(8,:,:) = CABLE_dif(8,:,:) + (f1->$var_name(8)$(month,:,:) - f->$var_name(8)$(month,:,:))*3600.*24.*dom(month)
    ; ESoil
       CABLE(9,:,:) = CABLE(9,:,:) + f->$var_name(9)$(month,:,:)*3600.*24.*dom(month)
       CABLE_dif(9,:,:) = CABLE_dif(9,:,:) + (f1->$var_name(9)$(month,:,:) - f->$var_name(9)$(month,:,:))*3600.*24.*dom(month)
    ; Qs
       CABLE(10,:,:) = CABLE(10,:,:) + f->$var_name(10)$(month,:,:)*3600.*24.*dom(month)
       CABLE_dif(10,:,:) = CABLE_dif(10,:,:) + (f1->$var_name(10)$(month,:,:) - f->$var_name(10)$(month,:,:))*3600.*24.*dom(month)
    ; Qsb
       CABLE(11,:,:) = CABLE(11,:,:) + f->$var_name(11)$(month,:,:)*3600.*24.*dom(month)
       CABLE_dif(11,:,:) = CABLE_dif(11,:,:) + (f1->$var_name(11)$(month,:,:) - f->$var_name(11)$(month,:,:))*3600.*24.*dom(month)
    end do
    SoilTemp = SoilTemp     + dim_avg_n_Wrap(f->SoilTemp,0)

    SoilMoist(0,0:5,:,:) = SoilMoist(0,0:5,:,:) + dim_avg_n_Wrap(f->SoilMoist(month_s:month_e,:,:,:),0) ; Actual
    SoilMoist(1,0:5,:,:) = SoilMoist(1,0:5,:,:) + dim_avg_n_Wrap(f1->SoilMoist(month_s:month_e,:,:,:),0) - dim_avg_n_Wrap(f->SoilMoist(month_s:month_e,:,:,:),0)      ; Equilibrium
    SoilMoist(0,6,:,:)   = SoilMoist(0,6,:,:)   + dim_avg_n_Wrap(f->GWMoist(month_s:month_e,:,:),0) ; Actual
    SoilMoist(1,6,:,:)   = SoilMoist(1,6,:,:)   + dim_avg_n_Wrap(f1->GWMoist(month_s:month_e,:,:),0) - dim_avg_n_Wrap(f->GWMoist(month_s:month_e,:,:),0)      ; Equilibrium
    SMP(0,:,:,:)       = SMP(0,:,:,:) + dim_avg_n_Wrap(f->SMP(month_s:month_e,:,:,:),0) ; Actual
    SMP(1,:,:,:)       = SMP(1,:,:,:) + dim_avg_n_Wrap(f1->SMP(month_s:month_e,:,:,:),0) - dim_avg_n_Wrap(f->SMP(month_s:month_e,:,:,:),0)     ; dif

    SoilTemp_dif  = SoilTemp_dif + (dim_avg_n_Wrap(f1->SoilTemp(month_s:month_e,:,:,:),0) - dim_avg_n_Wrap(f->SoilTemp(month_s:month_e,:,:,:),0))
    SoilMoist_dif(0:5,:,:) = SoilMoist_dif(0:5,:,:)+ (dim_avg_n_Wrap(f1->SoilMoist(month_s:month_e,:,:,:),0) - dim_avg_n_Wrap(f->SoilMoist(month_s:month_e,:,:,:),0))
    SoilMoist_dif(6,:,:)   = SoilMoist_dif(6,:,:)+ (dim_avg_n_Wrap(f1->GWMoist(month_s:month_e,:,:),0) - dim_avg_n_Wrap(f->GWMoist(month_s:month_e,:,:),0))
  end if 
    
end do ;year

CABLE     = CABLE/year_sum
SoilMoist = SoilMoist/year_sum
SMP       = SMP/year_sum * 9.8*0.001*0.001
SoilTemp  = SoilTemp/year_sum - 273.16

CABLE_dif     = CABLE_dif/year_sum
SoilTemp_dif  = SoilTemp_dif/year_sum
SoilMoist_dif = SoilMoist_dif/year_sum


;_____________________________ PLOT ___________________________________
pic = "spitial_map_debug_LH_vs_debug-tk-tvrad_"+scale+"_"+time_scale
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
;"BlueWhiteOrangeRed") ;"BlueYellowRed")             ;����ɫ��

res                    = True
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                            ;��ɫ
res@tmBorderThicknessF = 3.0

res@gsnDraw            = False  ; Don't draw plots
res@gsnFrame           = False  ; ����ҳ
res@lbLabelBarOn       = True

;************** ����labelbar ***************
res@lbBoxLinesOn       = True                       ;�ر�lbar box ����
res@lbTitleFont        = 25
res@lbLabelFont        = 25
res@lbTitleFontHeightF = 0.013
res@lbLabelFontHeightF = 0.013
res@txString      = ""
res@tmXBLabelFont      = 25 ;Sets the font index for the bottom X-Axis labels.
res@tmYLLabelFont      = 25

;*************** ���õ�ֵ�� *****************
res@cnLinesOn          = False                       ; �رյ�ֵ������
res@cnLineColor        = "black"
res@cnLineThicknessF   = 1.5
res@cnLineLabelsOn     = False
res@gsnMaximize        = True
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
pres@gsnPanelLabelBar   = False                    ; ���� panel labelbar.
pres@cnLevelSelectionMode = "ExplicitLevels" ;"ManualLevels"
pres@lbBoxLinesOn       = True                       ;�ر�lbar box ����
pres@lbTitleFont        = 25
pres@lbLabelFont        = 25
;pres@lbLabelBarOn       = False ; True                    ; �رո���plot��labelbar
pres@lbTitleFontHeightF = 0.01
pres@lbLabelFontHeightF = 0.01
;pres@lbAutoManage       = False
;pres@lbTitleDirection   = "Across"             ; title direction
;pres@lbTitlePosition    = "Bottom";Right"
;pres@lbLabelBarOn       = True

  plots = new(10,graphic)

; _______ CABLE _______
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  levels = fspan(-100,250,36)
  res@cnLevels = levels
  res@tiMainString = "(a) Debug_tk Rnet (W/m~S~2~N~)"
  plots(0) = gsn_csm_contour_map(wks,CABLE(0,{A:B},{C:D}),res)
  res@tiMainString = "(c) Debug_tk Qh (W/m~S~2~N~)"
  plots(2) = gsn_csm_contour_map(wks,CABLE(1,{A:B},{C:D}),res)
  res@tiMainString = "(e) Debug_tk Qg (W/m~S~2~N~)"
  plots(4) = gsn_csm_contour_map(wks,CABLE(2,{A:B},{C:D}),res)
  res@tiMainString = "(g) Debug_tk Qle (W/m~S~2~N~)"
  plots(6) = gsn_csm_contour_map(wks,CABLE(3,{A:B},{C:D}),res)

  delete(levels)
  delete(res@cnLevels)

  levels = fspan(0,1,21)
  res@cnLevels = levels
  res@tiMainString = "(i) Debug_tk Fsoil (-)"
  plots(8) = gsn_csm_contour_map(wks,CABLE(4,{A:B},{C:D}),res)

  delete(levels)
  delete(res@cnLevels)
  delete(res@cnFillPalette)
   
; _______ CABLE_dif _______
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = fspan(-40.,40.,41)
  res@cnLevels = levels
  res@tiMainString = "(b) Rnet_dif (Debug_tk - GW)(W/m~S~2~N~)"
  plots(1) = gsn_csm_contour_map(wks,CABLE_dif(0,{A:B},{C:D}),res)
  res@tiMainString = "(d) Qh_dif (W/m~S~2~N~)"
  plots(3) = gsn_csm_contour_map(wks,CABLE_dif(1,{A:B},{C:D}),res)
  res@tiMainString = "(f) Qg_dif (W/m~S~2~N~)"
  plots(5) = gsn_csm_contour_map(wks,CABLE_dif(2,{A:B},{C:D}),res)
  res@tiMainString = "(h) Qle_dif (W/m~S~2~N~)"
  plots(7) = gsn_csm_contour_map(wks,CABLE_dif(3,{A:B},{C:D}),res)

  delete(levels)
  delete(res@cnLevels)

  levels = fspan(-0.5,0.5,21)
  res@cnLevels = levels
  res@tiMainString = "(j) Fsoil_dif (-)"
  plots(9) = gsn_csm_contour_map(wks,CABLE_dif(4,{A:B},{C:D}),res)

  delete(levels)
  delete(res@cnLevels)

  gsn_panel(wks,(/plots/),(/5,2/),pres)

  delete(plots)
  

  plots = new(10,graphic)  

; _______ CABLE _______
if time_scale .eq. "annual" then
  unit = "(mm/y)"
  levels = fspan(0,800,41)
else
  unit = "(mm/mon)"
  levels = fspan(0,400,41)
end if
  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  res@cnLevels = levels
  res@tiMainString = "(a) Debug_tk Evap " + unit
  plots(0) = gsn_csm_contour_map(wks,CABLE(7,{A:B},{C:D}),res)
  res@tiMainString = "(c) Debug_tk TVeg " + unit
  plots(2) = gsn_csm_contour_map(wks,CABLE(8,{A:B},{C:D}),res)
  res@tiMainString = "(e) Debug_tk ESoil " + unit
  plots(4) = gsn_csm_contour_map(wks,CABLE(9,{A:B},{C:D}),res)
  res@tiMainString = "(g) Debug_tk Qs " + unit
  plots(6) = gsn_csm_contour_map(wks,CABLE(10,{A:B},{C:D}),res)
  res@tiMainString = "(i) Debug_tk Qsb " + unit
  plots(8) = gsn_csm_contour_map(wks,CABLE(11,{A:B},{C:D}),res)


  delete(levels)
  delete(res@cnLevels)
  delete(res@cnFillPalette)

; _______ CABLE_dif _______
if time_scale .eq. "annual" then
  levels = fspan(-200.,200.,41)
else
  levels = fspan(-50.,50.,21)
end if
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  res@cnLevels = levels
  res@tiMainString = "(b) Evap_dif (Debug_tk - GW)"+unit
  plots(1) = gsn_csm_contour_map(wks,CABLE_dif(7,{A:B},{C:D}),res)
  res@tiMainString = "(d) TVeg_dif "+unit
  plots(3) = gsn_csm_contour_map(wks,CABLE_dif(8,{A:B},{C:D}),res)
  res@tiMainString = "(f) Esoil_dif "+unit
  plots(5) = gsn_csm_contour_map(wks,CABLE_dif(9,{A:B},{C:D}),res)
  res@tiMainString = "(h) Qs_dif "+unit
  plots(7) = gsn_csm_contour_map(wks,CABLE_dif(10,{A:B},{C:D}),res)
  res@tiMainString = "(j) Qsb_dif "+unit
  plots(9) = gsn_csm_contour_map(wks,CABLE_dif(11,{A:B},{C:D}),res)
  delete(levels)
  delete(res@cnLevels)

  gsn_panel(wks,(/plots/),(/5,2/),pres)

  delete(plots)

; ______________________ Soil Layer _______________________
  plots = new(6,graphic)

do layer = 0,5

  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  levels = fspan(0.,0.6,31)
  res@cnLevels = levels
  res@tiMainString = "(a) Actual SoilMoist Layer-"+(layer+1)+" (mm~S~3~N~/mm~S~3~N~)"
  plots(0) = gsn_csm_contour_map(wks,SoilMoist(0,layer,{A:B},{C:D}),res)
  delete(levels)
  delete(res@cnLevels)

  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = fspan(-0.3,0.3,31)
  res@cnLevels = levels
  res@tiMainString = "(b) SoilMoist_dif Layer-"+(layer+1)+" (mm~S~3~N~/mm~S~3~N~)"
  plots(1) = gsn_csm_contour_map(wks,SoilMoist_dif(layer,{A:B},{C:D}),res)
  delete(levels)
  delete(res@cnLevels)

  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  levels = fspan(0.,40.,21)
  res@cnLevels = levels
  res@tiMainString = "(c) SoilTemp Layer-"+(layer+1)+" (~S~o~N~C)"
  plots(2) = gsn_csm_contour_map(wks,SoilTemp(layer,{A:B},{C:D}),res)
  delete(levels)
  delete(res@cnLevels)

  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = fspan(-2.,2.,21)
  res@cnLevels = levels
  res@tiMainString = "(d) SoilTemp_dif Layer-"+(layer+1)+" (~S~o~N~C)"
  plots(3) = gsn_csm_contour_map(wks,SoilTemp_dif(layer,{A:B},{C:D}),res)
  delete(levels)
  delete(res@cnLevels)

  res@cnFillPalette = "MPL_RdYlGn"
  levels = fspan(-1000.,0.,21)
  res@cnLevels = levels
  res@tiMainString = "(e) Actual SMP Layer-"+(layer+1)+" (MPa)"
  plots(4) = gsn_csm_contour_map(wks,SMP(0,layer,{A:B},{C:D}),res)
  delete(levels)
  delete(res@cnLevels)

  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = fspan(-500.,500.,21)
  res@cnLevels = levels
  res@tiMainString = "(f) Actual SMP_dif Layer-"+(layer+1)+" (MPa)"
  plots(5) = gsn_csm_contour_map(wks,SMP(1,layer,{A:B},{C:D}),res)
  delete(levels)
  delete(res@cnLevels)

  gsn_panel(wks,(/plots/),(/3,2/),pres)
end do
  delete(plots)
  plots = new(6,graphic)

  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  levels = fspan(0.,1.,21)
  res@cnLevels = levels
  res@tiMainString = "(a) Actual SoilMoist Layer-7 (mm~S~3~N~/mm~S~3~N~)"
  plots(0) = gsn_csm_contour_map(wks,SoilMoist(0,6,{A:B},{C:D}),res)
  delete(levels)
  delete(res@cnLevels)

  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = fspan(-0.3,0.3,31)
  res@cnLevels = levels
  res@tiMainString = "(b) SoilMoist_dif Layer-7 (mm~S~3~N~/mm~S~3~N~)"
  plots(1) = gsn_csm_contour_map(wks,SoilMoist_dif(6,{A:B},{C:D}),res)
  delete(levels)
  delete(res@cnLevels)


  res@cnFillPalette = "WhiteBlueGreenYellowRed"
  levels = fspan(0.,30.,31)
  res@cnLevels = levels
  res@tiMainString = "(c) WatTable (m)"
  plots(2) = gsn_csm_contour_map(wks,CABLE(5,{A:B},{C:D}),res)
  delete(levels)
  delete(res@cnLevels)
  
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = fspan(-10.,10,41)
  res@cnLevels = levels
  res@tiMainString = "(d) WatTable_dif (m)"
  plots(3) = gsn_csm_contour_map(wks,CABLE_dif(5,{A:B},{C:D}),res)
  delete(levels)
  delete(res@cnLevels)

  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = fspan(-20.,20.,41)
  res@cnLevels = levels
  res@tiMainString = "(e) Qrecharge "+unit
  plots(4) = gsn_csm_contour_map(wks,CABLE(6,{A:B},{C:D}),res)
  delete(levels)
  delete(res@cnLevels)
  
  res@cnFillPalette = "ViBlGrWhYeOrRe"
  levels = fspan(-10.,10,41)
  res@cnLevels = levels
  res@tiMainString = "(f) Qrecharge_dif "+unit
  plots(5) = gsn_csm_contour_map(wks,CABLE_dif(6,{A:B},{C:D}),res)
  delete(levels)
  delete(res@cnLevels)
  gsn_panel(wks,(/plots/),(/3,2/),pres)
end
