;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

begin

case_name = (/\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off/Outputs/gw_off/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs-standard/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para_redistrb-on/Outputs-standard/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_off_unify_para/Outputs-Haverd2013/gw_off",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on/Outputs/gw_on/",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-standard/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs-standard/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_unify_para/Outputs-Haverd2013/gw_on",\
              "/short/w35/mm3972/cable/runs/Marks_latest_branch_with_fixes/GSWP3_gw_on_no_aquifer_influence_unify_para/Outputs-Haverd2013/gw_on"\
             /)

case_name_brf = (/\
              "off-def-std",\
              "off-std",\
              "off-redstrb-std",\
              "off-Hvrd",\
              "on-def-std",\
              "on-std",\
              "on_norchrg-std",\
              "on-Hvrd",\
              "on_norchrg-Hvrd"\
               /)

Soil_thickness = (/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/) ; 6 soil layers and 1 aquifer

; ____________________________ Options _______________________________
scale = "SE-AU" ;"Global" ; "50S-N" ; "SE-AU" ; "AU"

if scale .eq. "Global" then
  A = -90.
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

; _____________________________ Check ______________________________
pic = "check"
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"ViBlGrWhYeOrRe")
res = True
res@cnFillMode         = "RasterFill"            ; Raster Mode
res@cnFillOn           = True                            ;???
res@tmBorderThicknessF = 3.0
res@cnLinesOn          = False                       ; ???????????
res@gsnDraw            = True  ; Don't draw plots
res@gsnFrame           = True  ; ?????

;______________________________ Data _________________________________
year_s = 2001
year_e = 2010
year_sum = year_e-year_s+1
total = year_sum*12
case_sum  = dimsizes(case_name)


LAI      = new((/3,total/),float)
Precip   = new((/3,total/),float)
Trans    = new((/case_sum,3,total/),float)

LAI      = 0.
Precip   = 0.
Trans    = 0.

filename = case_name(0)+"/cable_out_2000.nc"
f = addfile (filename,"r")

iveg = f->iveg

plot = gsn_csm_contour(wks,iveg,res)

LAI@_FillValue    = f->TVeg@_FillValue
Precip@_FillValue = f->TVeg@_FillValue
Trans@_FillValue  = f->TVeg@_FillValue

do case_num = 0, case_sum-1
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
       if case_num .eq. 0 then
          plot = gsn_csm_contour(wks,f->Rainf(1,{A:B},{C:D}),res)
          plot = gsn_csm_contour(wks,f->Rainf(1,{A:B},{C:D})*1.0,res)
          tmp = where(iveg({A:B},{C:D}) .ge. 1 .and. iveg({A:B},{C:D}) .le. 4 ,\
                                     f->LAI(1,{A:B},{C:D}),LAI@_FillValue )
          plot = gsn_csm_contour(wks,tmp,res)
       end if
       do month = 0,11
          if case_num .eq. 0 then
             Precip(0,counter) = avg(f->Rainf(month,{A:B},{C:D}))*3600*24*dom(month)
             Precip(1,counter) = avg(where(iveg({A:B},{C:D}) .ge. 1 .and. iveg({A:B},{C:D}) .le. 4 , \
                                     f->Rainf(month,{A:B},{C:D}),Precip@_FillValue ))\
                                     *3600*24*dom(month)
             Precip(2,counter) = avg(where(iveg({A:B},{C:D}) .ge. 5 .and. iveg({A:B},{C:D}) .le. 12 , \
                                     f->Rainf(month,{A:B},{C:D}),Precip@_FillValue ))\
                                     *3600*24*dom(month)
                                     
             LAI(0,counter)    = avg(f->LAI(month,{A:B},{C:D}))
             LAI(1,counter)    = avg(where(iveg({A:B},{C:D}) .ge. 1 .and. iveg({A:B},{C:D}) .le. 4 ,\
                                     f->LAI(month,{A:B},{C:D}),LAI@_FillValue ))
             LAI(2,counter)    = avg(where(iveg({A:B},{C:D}) .ge. 5 .and. iveg({A:B},{C:D}) .le. 12 ,\
                                     f->LAI(month,{A:B},{C:D}),LAI@_FillValue ))
             print(avg(f->Rainf(month,{A:B},{C:D}))*3600*24*dom(month)+" "+avg(f->Rainf(month,{A:B},{C:D})*1.0)*3600*24*dom(month))
          end if
          Trans(case_num,0,counter)    = avg(f->TVeg(month,{A:B},{C:D}))*3600*24*dom(month)
          Trans(case_num,1,counter)    = avg(where(iveg({A:B},{C:D}) .ge. 1 .and. iveg({A:B},{C:D}) .le. 4 , \
                                             f->TVeg(month,{A:B},{C:D}),Trans@_FillValue ))\
                                             *3600*24*dom(month)
          Trans(case_num,2,counter)    = avg(where(iveg({A:B},{C:D}) .ge. 5 .and. iveg({A:B},{C:D}) .le. 12 , \
                                             f->TVeg(month,{A:B},{C:D}),Trans@_FillValue ))\
                                             *3600*24*dom(month)
          counter = counter+1
       end do 
    end do ;year
end do ;case

LAI_avg    = dim_avg_n_Wrap(LAI,1)
Precip_avg = dim_avg_n_Wrap(Precip,1)
Trans_avg  = dim_avg_n_Wrap(Trans,2)

print(LAI_avg+" "+Precip_avg+" "+Trans_avg(0,:))

do pft = 0,2 
   LAI(pft,:)    = LAI(pft,:)/LAI_avg(pft)
   Precip(pft,:) = Precip(pft,:)/Precip_avg(pft)
   do case_num = 0,case_sum-1
      Trans(case_num,pft,:) = Trans(case_num,pft,:)/Trans_avg(case_num,pft)
   end do
end do

; __________________________ Deseasonalize _____________________________
; rmMonAnnCycTLL : Removes the annual cycle from "monthly" data.
do pft = 0,2 
   LAI(pft,:)   =rmAnnCycle1D(LAI(pft,:))
   Precip(pft,:)=rmAnnCycle1D(Precip(pft,:))
end do
Trans = rmMonAnnCycLLT(Trans) ; dimensioned time x lat x lon


; ________________________ Running average _____________________________

LAI   =runave_n_Wrap(LAI,3,0,1)
Precip=runave_n_Wrap(Precip,3,0,1)
Trans = runave_n_Wrap(Trans,3,0,2)



delete(pic)
delete(wks)
delete(res)
; _____________________________ PLOT ___________________________________
pic = "time_serial_Precip_LAI_Trans_"+scale
wks = gsn_open_wks("pdf",pic)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")

res                    = True              ; plot mods desired
res@gsnDraw            = False             ; don't draw yet
res@gsnFrame           = False             ; don't advance frame yet

res@vpHeightF 	 = 0.4               ; change aspect ratio of plot
res@vpWidthF 	   = 0.7

res@tiYAxisString	 = "Scale Anormaly (mm)"              ; add an axis title
res@txFontHeightF	 = 0.015            ; change title font heights

res@trYMaxF = 2.
res@trYMinF = -2.

LineColors       = (/"darkgreen","darkgreen","darkgreen","green","red","red","red","orange","orange"/) ; colors chosen
LineThicknesses	 = (/1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5/)      ; line thicknesses
DashPatterns	   = (/0.,0.,0.,0.,0.,0.,1.,0.,1./)      ; make all lines solid

x = ispan(1,total,1)
res@tmXBMode  = "Explicit"
res@tmXBValues=(/7,19,31,43,55,67,79,91,103,115/)
res@tmXBLabels=(/"2001","2002","2003","2004","2005","2006","2007","2008","2009","2010"/)
res@xyLineThicknessF   = 1.5

do pft = 0,2
   res@xyMonoLineColor    = False             ; want colored lines
   res@xyLineColor        = "Blue"
   top_plot = gsn_csm_xy (wks,x,Precip(pft,:),res)       ; create line plot

   res@xyLineColor        = "Purple"
   mid_plot  = gsn_csm_xy (wks,x,LAI(pft,:),res)  ; Create filled XY plot.
 
   res@xyLineColors       = LineColors
   res@xyLineThicknesses  = LineThicknesses
   res@xyDashPatterns	    = DashPatterns
   bot_plot  = gsn_csm_xy (wks,x,Trans(:,pft,:),res)  ; Create another filled XY plot.
   overlay(bot_plot,mid_plot)
   overlay(bot_plot,top_plot)
   pres = True
   maximize_output(wks,pres)
   delete(res@xyLineColors)
   delete(res@xyLineThicknesses)
   delete(res@xyDashPatterns)
end do

end 