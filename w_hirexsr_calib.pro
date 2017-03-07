
; function for checking old calib

FUNCTION check_calib_exist, u
;PRO check_prev_calib, u
shot = u.shot
module = u.module
path = '\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+strcompress(string(module), /remove_all)+'.ELLIPSE:COEFS'
;print, path
mdsopen, 'spectroscopy', shot
val = mdsvalue(path, /quiet, status = status)
mdsclose

if status EQ 1 then return, 1 else return, 0
return, 0
end

; check if any calibrations have been done
PRO check_prev_calib, u
shot = u.shot
module = u.module
path = '\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+strcompress(string(module), /remove_all)+'.ELLIPSE:COEFS'
;print, path
mdsopen, 'spectroscopy', shot
val = mdsvalue(path, /quiet, status = status)
mdsclose
;help, val

if status EQ 1 then begin
    ;   load data
    	WIDGET_CONTROL,/hourglass
        load_info,u
        load_fits,u
        load_ellipse,u
        plot_image,u
        plot_fits,u
        plot_ellipse,u	
        print_ellipse,u
        calib_load_preparer_name, u
endif

END

; place the user's name in the widget_box
PRO calib_load_user_name, u

user = REFORM(logname())
user = user(0)           ; recast

;widget_control, u.id.fitname, set_value = STRCOMPRESS(STRING(user), /remove_all)

END

;procedure for placing a previous preparers name in the widget box
PRO calib_load_preparer_name,u
shot = u.shot
module = u.module
path = '\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+strcompress(string(module), /remove_all)+':COMMENT'

mdsopen, 'spectroscopy', shot
string_name  = mdsvalue(path, status = status, /quiet)
mdsclose
; data exists
if status EQ 1 then begin
    ;substring = strsplit(string_name, /extract)
    ; 0- name
    ; 1- time
    ; etc. not yet defined
    widget_control, u.id.fitname, set_value = string_name
endif

END

; saves user name 
PRO calib_save_user_name, u

widget_control, u.id.fitname, get_value = user_name
date_time = systime()
output = strcompress(string(user_name), /remove_all)+','+date_time
shot = u.shot
module = u.module

path = '\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+strcompress(string(module), /remove_all)+':COMMENT'

mdsopen, 'spectroscopy', shot
;check if node exists
;265388144 - error for no node, intuitive really
;265388258 - error for no data in node
;mdsput, path, '$', output, status = status
test = mdsvalue(path, status = status, /quiet)
if status EQ 265388144L then begin
	;create the node
	mdstcl, "set verify"
	mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
	mdstcl, 'add node /usage=text '+path
	mdstcl, 'write'
	mdstcl, 'close'

	; load the node
	mdsopen, 'spectroscopy', shot
	mdsput, path, '$', output
	mdsclose
endif else begin
	mdsopen, 'spectroscopy', shot
	mdsput, path, '$', output
	mdsclose
endelse

;commented by C. Gao 03/17/2014
;mdsclose
	
END


PRO fit_ellipse,u
	widget_control,u.id.out_slider,get_value=out_slider
	x=size(u.image)
	index=u.ellipse.i
	good=u.ellipse.good
	bad=u.ellipse.bad
	xpt=u.fits.peaks[*,index]
	xerr=u.fits.err[*,index]
	ypt=findgen(x[1])
	out=nonlin_fit_ellipse(xpt,ypt,xerr,tilt=u.ellipse.tilt,good=good-bad)	;fits ellipse with all values not flagged as bad
	IF NOT u.ellipse.tilt THEN out=[out,0.0]				
	IF out_slider GT 0 THEN BEGIN						
		xplot=eq_ellipse(ypt,out)
		resid=xplot-xpt
		tmp=where(finite(xerr) EQ 0 OR bad EQ 1 OR good EQ 0)
		resid[tmp]=0							;set outliers of unused points to zero so they are not selected 
		out=reverse(sort(abs(resid)))					;make list of largest outliers
		u.ellipse.out[*,index]=0
		u.ellipse.out[out[0:out_slider-1],index]=1
		temp_good=good
		temp_good[out[0:out_slider-1]]=0
		out=nonlin_fit_ellipse(xpt,ypt,xerr,tilt=u.ellipse.tilt,good=temp_good-bad)	;fits ellipse with bad and outliers removed
		IF NOT u.ellipse.tilt THEN out=[out,0.0]
	ENDIF
	u.ellipse.param[*,index]=out
	tmp=where(u.ellipse.bad NE 1 AND u.ellipse.out[*,index] NE 1 AND finite(xerr) EQ 1)	;avoid either set as bad spectra of flagged as outliers
	xpt=u.fits.peaks[tmp,index]
	ypt=ypt[tmp]
	xellipse=eq_ellipse(ypt,u.ellipse.param[*,index])
	resid=xellipse-xpt
	mu=mean(resid)
	sigma=stdev(resid)
	u.ellipse.mu[index]=mu
	u.ellipse.sigma[index]=sigma
	widget_control,u.id.base, set_uvalue=u
END

PRO fit_all_ellipse,u
	index=u.ellipse.i
	npeaks=n(u.fits.peaks[0,*])+1
	FOR i=0,npeaks-1 DO BEGIN
		u.ellipse.i=i
		widget_control,u.id.base, set_uvalue=u
		fit_ellipse,u	
	ENDFOR
	u.ellipse.i=index
	widget_control,u.id.base, set_uvalue=u
END

PRO plot_ellipse,u	
	widget_control,u.id.ilow_slider,get_value=ilow
	widget_control,u.id.ihigh_slider,get_value=ihigh
	widget_control,u.id.vert_slider,get_value=gain
	gain/=100.0
	widget_control,u.id.draw3,get_value=draw_win
	IF NOT u.stat.ps THEN window,0,xsize=u.stat.xsize,ysize=u.stat.ysize,/pixmap
	x=size(u.image)
	ypt=findgen(x[1])
	index=u.ellipse.i
	tmp=where(u.ellipse.bad NE 1 AND u.ellipse.out[*,index] NE 1 AND finite(u.fits.err[*,index]) EQ 1)	;either set as bad spectra of flagged as outliers)
	xpt=u.fits.peaks[tmp,index]
	mu=u.ellipse.mu[index]
	sigma=u.ellipse.sigma[index]
	xerr=u.fits.err[tmp,index]
	ypt=ypt[tmp]
	xplot=eq_ellipse(ypt,u.ellipse.param[*,index])
	resid=xplot-xpt
	plot,ypt,resid,psym=8,yr=[-0.1,0.1]*10.0^(gain),symsize=0.5,xr=[0,x[1]-1],/xsty,ytit='RESIDUAL [PIX]',tit=n2g('mu')+'='+num2str(mu,dp=4)+'  '$
		+n2g('sigma')+'='+num2str(sigma,dp=3),/ysty
	oploterror,ypt,resid,xerr,psym=8,symsize=0.5
	oplot,[0,x[1]-1],[0,0],linestyle=2.0,thick=2.0
	IF ilow GT 0 THEN oplot,[ilow,ilow],[-0.1,0.1]*10^gain,linestyle=2,color=30,thick=4.0
	IF ihigh LT x[1]-1 THEN oplot,[ihigh,ihigh],[-0.1,0.1]*10^gain,linestyle=2,color=30,thick=4.0
	
	IF u.stat.iplot THEN BEGIN
		resid=u.info[*,index]-xpt
		oploterror,ypt,resid,xerr,psym=8,symsize=0.5,color=30,errcolor=30	
	ENDIF
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.stat.xsize,u.stat.ysize,0,0,0]	
	ENDIF	
END

PRO print_ellipse,u
	index=u.ellipse.i
	widget_control,u.id.xo,set_value=num2str(u.ellipse.param[0,index],dp=2)
	widget_control,u.id.yo,set_value=num2str(u.ellipse.param[1,index],dp=2)
	widget_control,u.id.a,set_value=num2str(u.ellipse.param[2,index],dp=2)
	widget_control,u.id.b,set_value=num2str(u.ellipse.param[3,index],dp=2)
	widget_control,u.id.phi,set_value=num2str(u.ellipse.param[4,index],dp=2)
END

PRO save_ellipse,u
	widget_control,u.id.ilow_slider,get_value=ilow
	widget_control,u.id.ihigh_slider,get_value=ihigh
	widget_control,u.id.out_slider,get_value=out_slider
	module=u.module
	shot=u.shot
	hirexsr_write_ellipse,shot,module,u.ellipse.param,int([ilow,ihigh]),u.ellipse.out,u.fits.lambda,u.ellipse.mu,u.ellipse.sigma,$
		u.ellipse.bad
	widget_control,u.id.message,set_value='ELLIPSE TO TREE FOR MOD'+num2str(module)+' ON SHOT: '+num2str(shot,1),/append
 	hirexsr_write_lambda,shot,module,good=u.stat.lamgood,order=u.stat.order
	widget_control,u.id.message,set_value='LAMBDA TO TREE FOR MOD'+num2str(module)+' ON SHOT: '+num2str(shot,1),/append

        ; now save the user name as well
        calib_save_user_name, u
END

PRO load_ellipse,u
	module=u.module
	shot=u.shot
	hirexsr_load_ellipse,shot,module,coefs,ifit_e,outl,lambda,mu,sigma,bad
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+'.ELLIPSE'
	u.ellipse.param=coefs
	u.fits.lambda=lambda
	i_fit=ifit_e
	u.ellipse.out=outl
	u.ellipse.mu=mu
	u.ellipse.sigma=sigma
	widget_control,u.id.ilow_slider,set_value=i_fit[0]
	widget_control,u.id.ihigh_slider,set_value=i_fit[1]
	widget_control,u.id.out_slider,set_value=total(u.ellipse.out[*,0])
	u.ellipse.bad=bad
	widget_control,u.id.base, set_uvalue=u
	widget_control,u.id.message,set_value='ELLIPSE FROM TREE FOR MOD'+num2str(module)+' ON SHOT: '+num2str(shot,1),/append
        ; now get the old user
        ;calib_load_preparer_name,u
END

PRO fit_image,u
	image=u.image
	rawimage=u.rawimage
	module=u.module
	shot=u.shot
	t1=u.t1
	t2=u.t2
	path=u.path.save
	double=u.stat.dp
	break=*u.brkpts
	splot=u.stat.splot
	x=size(image)
	resid=fltarr(x[1],x[2])
	widget_control,u.id.flow_slider,get_value=istart
	widget_control,u.id.fhigh_slider,get_value=istop
	FOR i=istart,istop DO BEGIN
		widget_control,u.id.message,set_value='Fitting '+num2str(i,1)+' out of '+num2str(x[1]-1,1),/append
		CASE u.stat.fitcase OF 
			0 : out=hirexsr_fit_he(reform(image[i,*]),double=double,break=break)
			1 : out=hirexsr_fit_h(reform(image[i,*]),double=double,break=break)
			2 : BEGIN
				IF shot GT 1140101000 THEN BEGIN
					out=hirexsr_fit_he_ca_mod(reform(image[i,*]),double=double,break=break,/focus)
                                ENDIF ELSE BEGIN
					out=hirexsr_fit_he_ca(reform(image[i,*]),double=double,break=break,focus=u.stat.focus)
                                ENDELSE
			END
			3 : out=hirexsr_fit_h_ca(reform(image[i,*]),double=double,break=break)
		ENDCASE
		resid[i,*]=out.residual
		IF i EQ istart THEN BEGIN
			IF keyword_set(double) THEN BEGIN
				n_lines=n(out.lab)+1
				nphot=dblarr(x[1],n_lines)
				peaks=dblarr(x[1],n_lines)
				width=dblarr(x[1],n_lines)
				off=dblarr(x[1])
				err=dblarr(x[1],n_lines)
			ENDIF ELSE BEGIN
				n_lines=n(out.lab)+1
				nphot=fltarr(x[1],n_lines)
				peaks=fltarr(x[1],n_lines)
				width=fltarr(x[1],n_lines)
				off=fltarr(x[1])
				err=fltarr(x[1],n_lines)
			ENDELSE
		ENDIF
		nphot[i,*]=2.0*sqrt(!pi/2.0)*out.c0
		peaks[i,*]=out.c1
		width[i,*]=out.c2
		off[i]=out.off
		IF splot THEN BEGIN
			p=fltarr(n_lines*3+1)
			p[n_lines*3]=off[i]
			p[indgen(n_lines)*3+0]=nphot[i,*]/(2.0*sqrt(!pi/2.0))
			p[indgen(n_lines)*3+1]=peaks[i,*]
			p[indgen(n_lines)*3+2]=width[i,*]
			plot_spec_whilefitting,u,p,i
		ENDIF
	ENDFOR
	label=out.lab
	IF u.stat.stree THEN BEGIN
		hirexsr_write_calibfits,shot,module,[t1,t2],image,label,nphot,off,peaks,width,resid,break,double,[istart,istop]
		widget_control,u.id.message,set_value='FIT TO TREE FOR MOD'+num2str(module)+' ON SHOT: '+num2str(shot,1),/append
	ENDIF ELSE BEGIN
		save,image,rawimage,resid,nphot,peaks,width,off,label,shot,module,t1,t2,double,break,istart,istop,filename=path
		widget_control,u.id.message,set_value='FIT TO SAVEFILE FOR MOD'+num2str(module)+' ON SHOT: '+num2str(shot,1),/append
	ENDELSE
END

PRO plot_spec_whilefitting,u,p,index
	x=size(u.image)
	pixel=indgen(x[2])
	spec=u.image[index,*]
	pixfit=make(0,x[2],2000)
	spec_fit = gaussian_fits(pixfit, p)
	resid=spec-gaussian_fits(pixel,p)
	n_lines=n(p)/3.0

	;plot total fit and individual lines
	widget_control,u.id.draw2,get_value=draw_win
	wset,draw_win
	plot, pixel,spec, psym=8,xtitle = 'Pixel #', ytitle = 'Counts',/xsty
	oplot, pixfit, spec_fit, color = 100
    	FOR i = 0, n_lines-1 DO oplot, pixfit, p[3*i]*exp(-(pixfit-p[3*i+1])^2/(2.*p[3*i+2]^2)), color=200, line=2,symsize=0.75
		
	;plot residual
	widget_control,u.id.draw3,get_value=draw_win
	wset,draw_win
	plot,pixel,resid,psym=8,ytit=n2g('Delta')+'Counts',/nodata,/xsty
	oplot, [0,x[2]],[0,0],linestyle=2.0
	oplot,pixel,resid,psym=8,color=200,symsize=0.75
END

PRO plot_calib_spec,u
	x=size(u.image)
	pixel=findgen(x[2])
	widget_control,u.id.spec_slider,get_value=index

	spec=u.image[index,*]
	resid=u.fits.resid[index,*]
	spec_fit = spec+resid
	
	;plot total fit and individual lines
	widget_control,u.id.draw2,get_value=draw_win
	IF u.stat.ps THEN BEGIN
		lthick=5.0
		ssize=0.85
        ENDIF ELSE BEGIN
		window,0,xsize=u.stat.xsize,ysize=u.stat.ysize,/pixmap
		lthick=1.0
		ssize=1.0
	ENDELSE
	plot, pixel,spec, psym=8,xtitle = 'Pixel #', ytitle = 'Counts',/xsty,symsize=ssize
	tmp=where(u.spec.label NE 'noline')
	totspec=fltarr(x[2])
	FOR i=0,n(tmp) DO BEGIN
		ispec=gaussian_fits(pixel,[u.spec.a0[index,i],u.spec.a1[index,i],u.spec.a2[index,i]])
		oplot,pixel,ispec,color=80,thick=lthick
		totspec+=ispec
	ENDFOR
	totspec+=u.spec.a3[index]
	oplot, pixel,totspec, color = 50,thick=lthick
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.stat.xsize,u.stat.ysize,0,0,0]	
	ENDIF	

	;plot residual
	widget_control,u.id.draw3,get_value=draw_win
	IF NOT u.stat.ps THEN window,0,xsize=u.stat.xsize,ysize=u.stat.ysize,/pixmap
	plot,pixel,resid,psym=8,ytit=n2g('Delta')+'Counts',/nodata,/xsty
	oplot, [0,x[2]],[0,0],linestyle=2.0
	oplot,pixel,resid,psym=8,color=200,symsize=0.75
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.stat.xsize,u.stat.ysize,0,0,0]	
	ENDIF
END

PRO load_fits,u
	path=u.path.save
	module=u.module
	shot=u.shot
	IF u.stat.stree THEN BEGIN
		hirexsr_load_calibfits,shot,module,times,image,label,nphot,offset,peaks,width,resid,break,double,ifit_f
		rawimage=image
		t1=times[0]
		t2=times[1]
		istart=ifit_f[0]
		istop=ifit_f[1]
	ENDIF ELSE restore, u.path.save		;save,image,rawimage,resid,nphot,peaks,width,off,label,shot,module,t1,t2,double,break,istart,istop,filename=path
	
        ;determine fitcase by looking at wavelength range if not stored (poorly hard-coded)
	x=size(label)
	IF label[1] EQ 'n5' THEN  u.stat.fitcase=0
	IF label[0] EQ 'lyas1' THEN  u.stat.fitcase=1
	IF x[1] LT 6 THEN BEGIN
		IF label[0] EQ 'w' THEN u.stat.fitcase=2
        ENDIF ELSE IF label[4] EQ 'w4n' OR label[5] EQ 'w4n' THEN  u.stat.fitcase=2
	IF last(label) EQ 'w7n'  THEN u.stat.fitcase=3
	reset_fitcase_buttons,u
	CASE u.stat.fitcase OF 
		0 : widget_control,u.id.ar_helikefit,set_button=1
		1 : widget_control,u.id.ar_hlikefit,set_button=1
		2 : widget_control,u.id.ca_helikefit,set_button=1
		3 : widget_control,u.id.ca_hlikefit,set_button=1
	ENDCASE	
	CASE u.stat.fitcase OF 
		0 : BEGIN					;he-like Ar
			findlab=['w','x','y','z']
			zlab=intarr(4)+18
			lamgood=[1,1,1,1]
			order=u.stat.order
                END
		1 : BEGIN					;h-like Ar
			findlab=['lya1','lya2','T','J']
			zlab=intarr(4)+18
			lamgood=[1,1,0,1]	
			order=u.stat.order
		END
		2 : BEGIN					;he-like Ca
			IF shot GT 1140101000 THEN BEGIN
				findlab=['w','wn3','x','y']
				zlab=intarr(4)+20
				lamgood=[1,0,0,1]
				order=1
                        ENDIF ELSE BEGIN
				findlab=['w','y','w4n','z']
				zlab=[20,20,18,20]
				lamgood=[1,1,1,1]
				order=u.stat.order
			ENDELSE
		END
		3 : BEGIN					;h-like Ca
			findlab=['w11n','w9n','w8n','w7n']
			zlab=[18,18,18,18]
			lamgood=[0,0,1,1]
			order=1
		END
	ENDCASE
	
	u.stat.lamgood=lamgood
	u.stat.order=order
	hirexsr_load_wavelengths,-1,lam_o,z_o,label_o	;wavelengths should be loaded from the hardcoded shot
	u.fits.resid=resid
	FOR i=0,n(findlab) DO BEGIN
		index=where(label EQ findlab[i])
		u.fits.peaks[*,i]=peaks[*,index]
		u.fits.err[*,i]=width[*,index]/sqrt(2.0*(nphot[*,index]-1.0))
		u.fits.lambda[i]=lam_o[where(z_o EQ zlab[i] AND label_o EQ findlab[i])]
	ENDFOR
	nlines=n(label)+1
	u.spec.a0[*,0:nlines-1]=nphot/(sqrt(!pi/2.0)*2.0)
	u.spec.a1[*,0:nlines-1]=peaks
	u.spec.a2[*,0:nlines-1]=width
	u.spec.a3=offset
	u.spec.label[0:nlines-1]=label
	u.spec.label[nlines:*]='noline'

	u.image=image
	u.rawimage=rawimage
	;name=''
	widget_control,u.id.shotid,set_value=num2str(shot,1)
	widget_control,u.id.module,set_value=num2str(module,1)
	widget_control,u.id.t1,set_value=num2str(t1,dp=2)
	widget_control,u.id.t2,set_value=num2str(t2,dp=2)
	widget_control,u.id.ilow_slider,set_value=istart
	widget_control,u.id.ihigh_slider,set_value=istop
	widget_control,u.id.flow_slider,set_value=istart
	widget_control,u.id.fhigh_slider,set_value=istop
	;widget_control,u.id.fitname,set_value=name
	*u.brkpts=break
	IF double THEN BEGIN
		u.stat.dp=1
		widget_control,u.id.precision,set_button=1
	ENDIF ELSE widget_control,u.id.precision,set_button=0
	widget_control,u.id.base, set_uvalue=u
	setbrk_text,u
	widget_control,u.id.message,set_value='FIT FROM TREE FOR'+num2str(module)+' ON SHOT: '+num2str(shot,1),/append
END

PRO plot_fits,u
	widget_control,u.id.draw2,get_value=draw_win
	IF NOT u.stat.ps THEN window,0,xsize=u.stat.xsize,ysize=u.stat.ysize,/pixmap
	x=size(u.image)
	plot, [0],[0],xr=[0,x[1]-1],yr=[0,x[2]-1],/xsty,/ysty
	pix=findgen(x[1])
	npeaks=n(u.fits.peaks[0,*])+1
	FOR i=0,npeaks-1 DO BEGIN
		oplot, pix,u.fits.peaks[*,i],psym=3
		IF total(u.ellipse.param[*,i]) NE 0 AND u.ellipse.i EQ i THEN BEGIN
			xplot=eq_ellipse(pix,u.ellipse.param[*,i])
			oplot,pix,xplot,color=100
			IF u.stat.iplot THEN oplot,pix,u.info[*,i],color=50
		ENDIF		
	ENDFOR
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.stat.xsize,u.stat.ysize,0,0,0]	
	ENDIF
END

PRO load_info,u
	IF u.stat.itree THEN info=hirexsr_read_treeinfo(u.shot,u.module) ELSE info=hirexsr_load_info(u.module)
	CASE u.stat.fitcase OF 
		0 : BEGIN					;he-like Ar
			findlab=['w','x','y','z']
			zlab=intarr(4)+18
                END
		1 : BEGIN					;h-like Ar
			findlab=['lya1','lya2','T','J']
			zlab=intarr(4)+18
		END
		2 : BEGIN					;he-like Ca
			findlab=['w','y','w4n','z']
			zlab=[20,20,18,20]
		END
		3 : BEGIN					;h-like Ca

		END
	ENDCASE
	hirexsr_load_wavelengths,-1,lam_o,z_o,label_o
	x=size(u.image)
	pix=findgen(x[1])
	FOR i=0,n(findlab) DO BEGIN
		lambda=lam_o[where(z_o EQ zlab[i] AND label_o EQ findlab[i])]
		param=genpos_spherical2quadcurve(info,lambda)
		u.info[*,i]=reverse(x[2]-1.0-ellipse_xpt(param,pix))		;change from (xi,zeta) to (i,j) of the image
	ENDFOR
	widget_control,u.id.base, set_uvalue=u
END


PRO get_image,u
	widget_control,u.id.shotid,get_value=shot
	widget_control,u.id.module,get_value=module
	widget_control,u.id.t1,get_value=t1
	widget_control,u.id.t2,get_value=t2
	

	u.module=module
	u.shot = shot
	u.t1 = t1
	u.t2 = t2
	widget_control,u.id.message,set_value='Loading Module '+num2str(module)+' for SHOT: '+num2str(shot,1),/append
	hirexsr_load_image,long(shot[0]),int(module[0]),image,t,break=brkpts
	ilow=ipt(t,float(t1[0]))
	ihigh=ipt(t,float(t2[0]))
	image=sum_array(float(image[*,*,ilow:ihigh]),/k)
	*u.brkpts=brkpts
	u.image=image
	u.rawimage=image
	widget_control,u.id.base, set_uvalue=u
	setbrk_text,u
END

PRO setbrk_text,u
	brkpts=*u.brkpts
	brkstr=num2str(int(brkpts[0]),1)
	FOR i=1,n(brkpts) DO brkstr+=','+num2str(int(brkpts[i]),1)
	widget_control,u.id.brkpts, set_value=brkstr
END

PRO plot_image,u
	widget_control,u.id.flow_slider,get_value=istart
	widget_control,u.id.fhigh_slider,get_value=istop
	widget_control,u.id.gain_slider,get_value=gain
	widget_control,u.id.draw1,get_value=draw_win
	IF NOT keyword_set(gain) THEN gain= 1.0
	IF NOT u.stat.ps THEN window,0,xsize=u.stat.xsize,ysize=u.stat.ysize,/pixmap
	image=u.image
	x=size(image)
	nx=x[1]
	ny=x[2]
	zoom=1.5
	i_pt=indgen(nx*zoom)/zoom
	j_pt=indgen(ny*zoom)/zoom
	new_pic=interpolate(image,i_pt,j_pt,/grid)
	loadct,39,/silent
	order=sort(new_pic)
	norm=max(new_pic[order[0:int(nx*zoom)*int(ny*zoom)-20]])
	tv,new_pic/norm*256.0*gain
	loadct,12,/silent
	IF istart NE 0 THEN plots,istart/486.0*[1,1],[0,1],linestyle=2.0,thick=3.0,color=255,/norm
	IF istop NE 486 THEN plots,istop/486.0*[1,1],[0,1],linestyle=2.0,thick=3.0,color=255,/norm
	brkpts=*u.brkpts
	IF n(brkpts) GT 0 THEN FOR i=0,n(brkpts) DO plots,[0,1],brkpts[i]/195.0,thick=2.0,color=30,/norm,linestyle=4
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.stat.xsize,u.stat.ysize,0,0,0]	
	ENDIF
END

FUNCTION spathstr,u
	user=logname()
	module=u.module
	shot=u.shot
	t1=u.t1
	t2=u.t2
	IF u.stat.stree THEN path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1) ELSE $
		path='/home/'+user+'/idl/hirexsr/calib/ellipse_'+num2str(module,1)+'_'+num2str(shot,1)+'_'+num2str(int(t1*1000),1)+'_'+num2str(int(t2*1000),1)+'.dat'
	RETURN,path
END

FUNCTION ipathstr,u
	user=logname()
	module=u.module
	shot=u.shot
	t1=u.t1
	t2=u.t2
	IF u.stat.itree THEN path='\SPECTROSCOPY::TOP.HIREXSR.INFO.MOD'+num2str(u.module,1) ELSE $
		path='/home/'+user+'/idl/genie/data/info/hirexsr/hirexsr_'+num2str(module,1)+'.info'
	RETURN,path
END
PRO reset_fit_buttons,u
	widget_control,u.id.fit_lya1,set_button=0
	widget_control,u.id.fit_lya2,set_button=0
	widget_control,u.id.fit_T,set_button=0
	widget_control,u.id.fit_J,set_button=0
	widget_control,u.id.fit_w,set_button=0
	widget_control,u.id.fit_x,set_button=0
	widget_control,u.id.fit_y,set_button=0
	widget_control,u.id.fit_z,set_button=0
END

PRO reset_fitcase_buttons,u
	widget_control,u.id.ar_helikefit,set_button=0
	widget_control,u.id.ar_hlikefit,set_button=0
	widget_control,u.id.ca_helikefit,set_button=0
	widget_control,u.id.ca_hlikefit,set_button=0
END

PRO reset_ellipse_text,u

END

PRO w_hirexsr_calib_event,event

	widget_control,event.top,get_uvalue=u
	id = u.id
	tag = tag_names(event,/st)
	button=' '
	idtags=tag_names(id)
	FOR i=0,n(idtags) DO IF id.(i) EQ event.id THEN ename=idtags[i]
	CASE tag OF
		"WIDGET_BUTTON": BEGIN
			widget_control,event.id,get_value=button,get_uvalue=uvalue
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF
				
	 			"QUIT":  BEGIN
					widget_control,event.top,/destroy
					heap_free,u.brkpts
					heap_gc
				END
				"PRINT": BEGIN
					u.stat.ps=1
					d_old=!d
					;!p.multi=[0,0,3]
					set_plot,'ps'
					!p.font=00
					device, /color
					device, /portrait
					device, xsize=7.0,ysize=7.0*195/487.0,/in
					device, encapsulated=0
					device, preview=0
					device, yoffset=3.0
					plot_image,u
					plot_calib_spec,u
					plot_fits,u
					plot_ellipse,u	
					;!p.multi=0
					psc
					device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
					xwplot
					set_plot,'x'
					u.stat.ps=0
				END
				"STOPBUTTON" : stop
				"GET":  BEGIN
					get_image,u
					load_info,u
					plot_image,u
                                        ;check_prev_calib, u
                                        
                                	;reset FITS/SPEC/ELLIPSE structures
					npeaks=4
					mxp=40
					u.fits={resid:fltarr(487,195)+1,peaks:fltarr(487,npeaks),err:fltarr(487,npeaks),lambda:fltarr(npeaks)}
					u.spec={a0:fltarr(487,mxp),a1:fltarr(487,mxp),a2:fltarr(487,mxp),a3:fltarr(487),label:strarr(mxp)}
					u.ellipse={i:0,tilt:0,param:fltarr(5,npeaks),good:intarr(487)+1,bad:intarr(487),out:intarr(487,npeaks),mu:fltarr(npeaks),sigma:fltarr(npeaks)}
					print_ellipse,u
					widget_control,u.id.ilow_slider,set_value=0
					widget_control,u.id.ihigh_slider,set_value=486
					widget_control,u.id.out_slider,set_value=0

					;clear screens
					blank=fltarr(u.stat.xsize,u.stat.ysize)
					widget_control,u.id.draw2,get_value=draw_win
					wset,draw_win
					tv,blank
					widget_control,u.id.draw3,get_value=draw_win
					wset,draw_win
					tv,blank
				END
				"LOADDATA" : BEGIN
					WIDGET_CONTROL,/hourglass
					load_info,u
					load_fits,u
					load_ellipse,u
					plot_image,u
					plot_fits,u
					plot_ellipse,u	
					print_ellipse,u
                                        calib_load_preparer_name, u
                                END

				"LOADINFO" : BEGIN
					load_info,u
					IF u.stat.iplot THEN BEGIN
						plot_fits,u
						plot_ellipse,u
					ENDIF	
				END
				"FITSPEC" : BEGIN
					; now check to make sure user actually wants 
					  cont = 0
					  exist = check_calib_exist(u)
					  if exist EQ 1 then begin
					     inquiry = dialog_message('Previos calibration exists. Overwrite?', /default_no, /question)
					     if inquiry EQ 'Yes' then cont =1 else cont = 0
					  endif else begin
					  	cont = 1
					  endelse
					;
					if cont EQ 1 then begin
					  WIDGET_CONTROL,/hourglass
					  calib_load_user_name, u
					  fit_image,u
					  load_fits,u
					  plot_fits,u
					  fit_all_ellipse,u
					  save_ellipse,u
					  plot_fits,u
					  plot_ellipse,u	
					  print_ellipse,u
                                          calib_load_user_name, u   
					endif else begin
					    calib_load_preparer_name, u
					endelse
				END
				"FITELLIPSE" : BEGIN
					WIDGET_CONTROL,/hourglass
					fit_all_ellipse,u
					plot_fits,u
					plot_ellipse,u	
					print_ellipse,u
                                        calib_load_user_name, u
				END
				"WRTELLIPSE" : BEGIN
					WIDGET_CONTROL,/hourglass
					save_ellipse,u	
				END
				"TREESAVE" : BEGIN
					IF event.select EQ 1 THEN u.stat.stree=1 ELSE u.stat.stree=0
					u.path.save=spathstr(u)
					widget_control,id.savepath,set_value=u.path.save
				END
				"TREEINFO" : BEGIN
					IF event.select EQ 1 THEN u.stat.itree=1 ELSE u.stat.itree=0
					u.path.info=ipathstr(u)
					widget_control,id.infopath,set_value=u.path.info
					load_info,u
					IF u.stat.iplot THEN BEGIN
						plot_fits,u
						plot_ellipse,u
					ENDIF	
				END

				"PRECISION" : IF event.select EQ 1 THEN u.stat.dp=1 ELSE u.stat.dp=0
				"SPLOT" : IF event.select EQ 1 THEN u.stat.splot=1 ELSE u.stat.splot=0
				"PHIFIT" : BEGIN
					IF event.select EQ 1 THEN u.ellipse.tilt=1 ELSE u.ellipse.tilt=0
				END
				"INFOPLOT" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						u.stat.iplot=1 
						plot_fits,u
						plot_ellipse,u
					ENDIF ELSE BEGIN
						u.stat.iplot=0
						plot_fits,u
						plot_ellipse,u
					ENDELSE
                                     END
				"HLIKEFIT": IF event.select EQ 1 THEN u.stat.hlikefit=1 ELSE u.stat.hlikefit=0
				"AR_HELIKEFIT" : BEGIN
						reset_fitcase_buttons,u
						widget_control,u.id.ar_helikefit,set_button=1
						u.stat.fitcase=0
						widget_control,u.id.message,set_value='Spectra will be fitted assuming He-like Ar',/append
                                 END
				"AR_HLIKEFIT" : BEGIN
						reset_fitcase_buttons,u
						widget_control,u.id.ar_hlikefit,set_button=1
						u.stat.fitcase=1
						widget_control,u.id.message,set_value='Spectra will be fitted assuming H-like Ar',/append
                                            END
				"CA_HELIKEFIT" : BEGIN
						reset_fitcase_buttons,u
						widget_control,u.id.ca_helikefit,set_button=1
						u.stat.fitcase=2
						widget_control,u.id.message,set_value='Spectra will be fitted assuming He-like Ca',/append
                                 END
				"CA_HLIKEFIT" : BEGIN
						reset_fitcase_buttons,u
						widget_control,u.id.ca_hlikefit,set_button=1
						u.stat.fitcase=3
						widget_control,u.id.message,set_value='Spectra will be fitted assuming H-like Ca',/append
				END
				"FIT_LYA1" : BEGIN		
					IF u.stat.fitcase EQ 1 OR u.stat.fitcase EQ 3 THEN BEGIN
						u.ellipse.i=0
						reset_fit_buttons,u
						widget_control,u.id.fit_lya1,set_button=1
						plot_fits,u
						print_ellipse,u
						plot_ellipse,u
					ENDIF ELSE widget_control,u.id.fit_lya1,set_button=0
				END
				"FIT_LYA2" : BEGIN
					IF u.stat.fitcase EQ 1 OR u.stat.fitcase EQ 3 THEN BEGIN
						u.ellipse.i=1
						reset_fit_buttons,u
						widget_control,u.id.fit_lya2,set_button=1
						plot_fits,u
						print_ellipse,u
						plot_ellipse,u
					ENDIF ELSE widget_control,u.id.fit_lya2,set_button=0
				END
				"FIT_T" : BEGIN
					IF u.stat.fitcase EQ 1 OR u.stat.fitcase EQ 3 THEN BEGIN
						u.ellipse.i=2
						reset_fit_buttons,u
						widget_control,u.id.fit_T,set_button=1
						plot_fits,u
						print_ellipse,u
						plot_ellipse,u
					ENDIF ELSE widget_control,u.id.fit_T,set_button=0
				END
				"FIT_J" : BEGIN
					IF u.stat.fitcase EQ 1 OR u.stat.fitcase EQ 3 THEN BEGIN
						u.ellipse.i=3
						reset_fit_buttons,u
						widget_control,u.id.fit_J,set_button=1
						plot_fits,u
						print_ellipse,u
						plot_ellipse,u	
					ENDIF ELSE widget_control,u.id.fit_J,set_button=0
				END
				"FIT_W" : BEGIN
					IF u.stat.fitcase EQ 0 OR u.stat.fitcase EQ 2 THEN BEGIN
						u.ellipse.i=0
						reset_fit_buttons,u
						widget_control,u.id.fit_w,set_button=1
						plot_fits,u
						print_ellipse,u
						plot_ellipse,u	
					ENDIF ELSE widget_control,u.id.fit_w,set_button=0
				END
				"FIT_X" : BEGIN
					IF u.stat.fitcase EQ 0 OR u.stat.fitcase EQ 2 THEN BEGIN
						u.ellipse.i=1
						reset_fit_buttons,u
						widget_control,u.id.fit_x,set_button=1
						plot_fits,u
						print_ellipse,u
						plot_ellipse,u	
					ENDIF ELSE widget_control,u.id.fit_x,set_button=0
				END
				"FIT_Y" : BEGIN
					IF u.stat.fitcase EQ 0 OR u.stat.fitcase EQ 2 THEN BEGIN
						u.ellipse.i=2
						reset_fit_buttons,u
						widget_control,u.id.fit_y,set_button=1
						plot_fits,u
						print_ellipse,u
						plot_ellipse,u	
					ENDIF ELSE widget_control,u.id.fit_y,set_button=0
				END
				"FIT_Z" : BEGIN
					IF u.stat.fitcase EQ 0 OR u.stat.fitcase EQ 2 THEN BEGIN
						u.ellipse.i=3
						reset_fit_buttons,u
						widget_control,u.id.fit_z,set_button=1
						plot_fits,u
						print_ellipse,u
						plot_ellipse,u	
					ENDIF ELSE widget_control,u.id.fit_z,set_button=0
				END
				"IBAD" :  BEGIN
					widget_control,u.id.spec_slider,get_value=index
					IF event.select EQ 1 THEN u.ellipse.bad[index]=1 ELSE u.ellipse.bad[index]=0
				END
                              	"LEFT_SLIDER_BUTTON": BEGIN
					widget_control, u.id.spec_slider, get_value = index
 	                               	IF index GT 0 THEN BEGIN
                                       		widget_control, u.id.spec_slider, set_value = index-1
                                       		index = index-1
                                       		IF u.fits.resid[0,0] NE 1.0 THEN plot_calib_spec,u
                                       		IF u.ellipse.bad[index] EQ 1 THEN widget_control,u.id.ibad,set_button=1 ELSE widget_control,u.id.ibad,set_button=0 
                                 	 ENDIF
                                END
                                "RIGHT_SLIDER_BUTTON": BEGIN
                                	widget_control, u.id.spec_slider, get_value = index
                                   	sliderminmax = widget_info(u.id.spec_slider, /slider_min_max)
                                   	slidermax = sliderminmax(1);  two element array
                                   	IF index LT slidermax THEN BEGIN
                                       		widget_control, u.id.spec_slider, set_value = index+1
                                       		index = index +1
                                       		IF u.fits.resid[0,0] NE 1.0 THEN plot_calib_spec,u
                                      		IF u.ellipse.bad[index] EQ 1 THEN widget_control,u.id.ibad,set_button=1 ELSE widget_control,u.id.ibad,set_button=0 
                                   	ENDIF
                                END
				ELSE:
			ENDCASE
		END
  		"WIDGET_SLIDER": BEGIN
			widget_control,event.id,get_value=slider
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF 
				'GAIN_SLIDER' : IF total(u.image NE 0) THEN plot_image,u

				'BIN_SLIDER' : BEGIN
					IF total(u.image NE 0) THEN BEGIN
						IF slider NE 1 THEN BEGIN
							num_bins=floor(487.0/slider)
							FOR i=0,num_bins-2 DO BEGIN
								line_sum=total(u.rawimage[i*slider:(i+1)*slider-1,*],1)
								FOR j=0,slider-1 DO u.image[i*slider+j,*]=line_sum
							ENDFOR
						ENDIF ELSE u.image=u.rawimage
						plot_image,u
					ENDIF
				END

				'VERT_SLIDER' : IF total(u.ellipse.param[*,u.ellipse.i]) NE 0 THEN plot_ellipse,u

				'FLOW_SLIDER' : IF total(u.image NE 0) THEN plot_image,u
				
				'FHIGH_SLIDER' : IF total(u.image NE 0) THEN plot_image,u	
			
				'ILOW_SLIDER' : BEGIN
					widget_control,u.id.ihigh_slider,get_value=ihigh
					u.ellipse.good[slider+1:ihigh]=1
					u.ellipse.good[0:slider]=0
					IF total(u.ellipse.param[*,u.ellipse.i]) THEN BEGIN
						fit_ellipse,u
						plot_fits,u
						plot_ellipse,u
						print_ellipse,u
					ENDIF		
				END
				'IHIGH_SLIDER' : BEGIN
					widget_control,u.id.ilow_slider,get_value=ilow
					u.ellipse.good[ilow:slider-1]=1
					u.ellipse.good[slider:*]=0
					IF total(u.ellipse.param[*,u.ellipse.i]) THEN BEGIN
						fit_ellipse,u
						plot_fits,u
						plot_ellipse,u
						print_ellipse,u
					ENDIF		

				END
				'OUT_SLIDER' : BEGIN
					IF total(u.ellipse.param[*,u.ellipse.i]) THEN BEGIN
						IF slider EQ 0 THEN u.ellipse.out[*]=0
						fit_ellipse,u
						plot_fits,u
						plot_ellipse,u
						print_ellipse,u
					ENDIF	
				END
				'SPEC_SLIDER' : BEGIN
					IF u.fits.resid[0,0] NE 1.0 THEN plot_calib_spec,u
					IF u.ellipse.bad[slider] EQ 1 THEN widget_control,u.id.ibad,set_button=1 ELSE widget_control,u.id.ibad,set_button=0
				END
				ELSE:
			ENDCASE
		END
   		"WIDGET_TEXT_CH": BEGIN
			CASE event.id OF 
				u.id.shotid : BEGIN
					u.path.save=spathstr(u)
					widget_control,u.id.shotid,get_value=shot
					u.shot=shot
					u.path.save=spathstr(u)
					widget_control,id.savepath,set_value=u.path.save
				END
				u.id.module : BEGIN
					widget_control,u.id.module,get_value=module
					u.module=module
					u.path.save=spathstr(u)
					u.path.info=ipathstr(u)
					widget_control,id.savepath,set_value=u.path.save
					widget_control,id.infopath,set_value=u.path.info
				END
				u.id.t1 : BEGIN
					widget_control,u.id.t1,get_value=t1
					u.t1=t1
					u.path.save=spathstr(u)
					widget_control,id.savepath,set_value=u.path.save
				END
				u.id.t2 : BEGIN
					widget_control,u.id.t2,get_value=t2
					u.t2=t2
					u.path.save=spathstr(u)
					widget_control,id.savepath,set_value=u.path.save
				END
				u.id.savepath : BEGIN
					widget_control,u.id.savepath,get_value=path
					u.path.save=path
					widget_control,id.savepath,set_value=u.path.save
				END
				u.id.infopath : BEGIN
					widget_control,u.id.info,get_value=path
					u.path.info=path
					widget_control,id.info,set_value=u.path.info
				END
				u.id.brkpts : BEGIN
					widget_control,u.id.brkpts,get_value=brkpts
					brkpts=int(strsplit(brkpts,',',/extract))
					*u.brkpts=brkpts
					plot_image,u	
				END

				ELSE: 
			ENDCASE
                END
		"WIDGET_DRAW" : BEGIN
			CASE event.id OF 
				u.id.draw3 : BEGIN
					IF event.release NE 0 THEN BEGIN
						click_loc = convert_coord(event.x, event.y, /device, /to_data)
						widget_control,u.id.spec_slider,set_value=int(click_loc[0])
						IF u.ellipse.bad[int(click_loc[0])] EQ 1 THEN widget_control,u.id.ibad,set_button=1 ELSE widget_control,u.id.ibad,set_button=0
						IF u.fits.resid[0,0] NE 1.0 THEN plot_calib_spec,u
					ENDIF
				END
				ELSE :
			ENDCASE	
		END
		ELSE: 
	ENDCASE
	IF button NE 'QUIT' THEN widget_control,event.top,set_uvalue=u
      	
	
END

;+
;NAME:
;	W_HIREXSR_CALIB
;
;PURPOSE:
;	This launches the HIREXSR spectral calibration widget which is used to general the wavelength map
;	used by the W_HIREXSR_HE_MOMENTS and W_HIREXSR_H_MOMENTS.  It performs the spectral fits which are used by W_HIREXSR_DET_ALIGN.
;
;CALLING SEQUENCE:
;	W_HIREXSR_CALIB can be called after running @W_HIREXSR_CALIB.BAT
;
;PROCEDURE:
;	A shot number, module and time range (t1,t2) are entered for one of the locked mode calibrations.
;	The fit low and fit high sliders are used to identify a subset of the image on which to perform the non-linear fitting
;	Spectra are fit, the coefficients saved to the tree or the selected .DAT file.  These are used to calculate the (xo,yo,a,b,phi) values for an
;	ellipse that is fit to a set of the bright lines in the spectrum.
;	Fits can be inspected and removed from the ellipse fitting process in order to reduce the residual of ellipse fit to the peaks of various lines.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 2/10
;	9/9/10		M.L. Reinke - removed zlow slider and changed break to use a vector
;
;-

PRO w_hirexsr_calib,shot=shot
	user=logname()
	
	loadct,12,/silent
	mdim=get_screen_size()
	IF mdim[0] NE 1600 AND mdim[1] NE 1200 THEN base=widget_base(title='HIREXSR SPECTRAL CALIBRATION',/row,/scroll,$
		x_scroll_size=mdim[0]*0.95,y_scroll_size=mdim[1]*0.85) ELSE base=widget_base(title='HIREXSR SPECTRAL CALIBRATION',/row)
	left=widget_base(base,/column)
	center=widget_base(base,/column)
	right=widget_base(base,/column)
	
	;setup plotting windows
	l1=widget_base(left,frame=5)
	draw1=widget_draw(l1,xsize=487*1.5,ysize=195*1.5)
	l2=widget_base(left,frame=5)
	draw2=widget_draw(l2,xsize=487*1.5,ysize=195*1.5)
	l3=widget_base(left,frame=5)
	draw3=widget_draw(l3,xsize=487*1.5,ysize=195*1.5,/button_events)

	c1space=widget_base(center,ysize=195*3+75)
	vert_slider=widget_slider(center,ysize=195,min=1.0,max=300.0,value=10.0,/drag,/vert,/suppress)

	;setup information
	rtop=widget_base(right,/row)
	ra=widget_base(rtop,/column,/frame)
	rb=widget_base(rtop,/column,/nonexcl,/frame)
	ar_helikefit=widget_button(rb,value='Use He-like Ar Fitting')
	ar_hlikefit=widget_button(rb,value='Use H-like Ar Fitting')
	ca_helikefit=widget_button(rb,value='Use He-like Ca Fitting')
	ca_hlikefit=widget_button(rb,value='Use H-like Ca Fitting')

	r1=widget_base(ra,/row)
	dum = widget_label(r1,value='SHOT: ')
	shotid = widget_text(r1,xsize=10,ysize=1,/edit)
	dum = widget_label(r1,value='   MODULE: ')
	module = widget_text(r1,xsize=2,ysize=1,/edit)
	dum = widget_label(r1,value='  ')
	get = widget_button(r1,value='LOAD')
	dum = widget_label(r1,value='  ')
	quit= widget_button(r1,value='QUIT')
	dum = widget_label(r1,value='  ')
	print= widget_button(r1,value='PRINT')
	stopbutton=widget_button(r1,value='STOP')

	r2=widget_base(ra,/row)
	dum = widget_label(r2,value='t1=')
	t1 = widget_text(r2,xsize=5,ysize=1,/edit)
	dum = widget_label(r2,value='     t2=')
	t2 = widget_text(r2,xsize=5,ysize=1,/edit)
	dum = widget_label(r2,value='  BIN')
	bin_slider=widget_slider(r2,xsize=100,min=1.0,max=40.0,value=1.0,/drag)
	dum=widget_label(r2,value=' BREAK PTS:  ')
	brkpts=widget_text(r2,xsize=15,ysize=1,/edit)
	r2p1=widget_base(ra,/row)
	dum = widget_label(r2p1,value='GAIN  ')
	gain_slider=widget_slider(r2p1,xsize=400,min=1.0,max=80.0,value=1.0,/drag)
	r2p2=widget_base(right,/row)
	message = widget_text(r2p2,xsize=70,ysize=3,/scroll)

	r2p3=widget_base(right,/row)
	dum = widget_label(r2p3,value='FIT LOW ')
	flow_slider = widget_slider(r2p3,xsize=486,min=0,max=486,value=0,/drag)
	r2p4=widget_base(right,/row)
	dum = widget_label(r2p4,value='FIT HIGH')
	fhigh_slider = widget_slider(r2p4,xsize=486,min=0,max=486,value=486,/drag)
	r2space=widget_base(right,/row,ysize=25)
	r5=widget_base(right,/row)
	r6=widget_base(r5,/row,/nonexcl)

	precision = widget_button(r6,value='DOUBLE PRECISION')
	splot = widget_button(r6,value='PLOT SPECTRAL FITS')
	fitspec = widget_button(r5,value='FIT/SAVE SPECTRA')
	dum = widget_label(r5,value='  CALIB BY: ')
	fitname=widget_text(r5,xsize=27,ysize=1)
	
	r3=widget_base(right,/row)
	r3p1=widget_base(r3,/row,/nonexclu)
	treesave=widget_button(r3p1,value='USE TREE')
	dum = widget_label(r3,value='SAVE PATH')
	savepath = widget_text(r3,xsize=70,ysize=1,/edit)
	dum = widget_label(r3,value='  ')
	loaddata = widget_button(r3,value='LOAD')
	r7p1=widget_base(right,/row)
	r7p1p1=widget_base(r7p1,/row)
	dum = widget_label(r7p1p1,value='SPEC ')
	spec_slider = widget_slider(r7p1p1,xsize=486,min=0,max=486,value=0,/drag)
        left_slider_button  = widget_button(r7p1p1, VALUE = '<-', xsize = 30)        
        right_slider_button = widget_button(r7p1p1, VALUE = '->', xsize = 30)
	r7p1p2=widget_base(r7p1,/row,/nonexcl)
	ibad=widget_button(r7p1p2,value='BAD')

	r2space=widget_base(right,/row,ysize=10)

	r4=widget_base(right,/row)
	r4p1=widget_base(r4,/row,/nonexclu)
	treeinfo=widget_button(r4p1,value='USE TREE')
	dum = widget_label(r4,value='INFO PATH')
	infopath = widget_text(r4,xsize=70,ysize=1,/edit)
	dum = widget_label(r4,value='  ')
	loadinfo = widget_button(r4,value='LOAD')
	r7space=widget_base(right,/row,ysize=55)



	r8p2=widget_base(right,/row)
	dum = widget_label(r8p2,value='                         SELECT LINE FOR ELLIPSE RESIDUAL PLOTTING            ')
	r8p1=widget_base(right,/row,/nonexcl)
	infoplot=widget_button(r8p1,value='PLOT INFO  ')
	phifit=widget_button(r8p1,value='ALLOW TILT        ')
	fit_lya1=widget_button(r8p1,value='lya1')
	fit_lya2=widget_button(r8p1,value='lya2')
	fit_T=widget_button(r8p1,value='T')
	fit_J=widget_button(r8p1,value='J     ')
	fit_w=widget_button(r8p1,value='w')
	fit_x=widget_button(r8p1,value='x')
	fit_y=widget_button(r8p1,value='y')
	fit_z=widget_button(r8p1,value='z')
	r8=widget_base(right,/row)
	dum = widget_label(r8,value='xo=')
	xo = widget_text(r8,xsize=7)
	dum = widget_label(r8,value=' yo=')
	yo = widget_text(r8,xsize=7)
	dum = widget_label(r8,value=' a=')
	a = widget_text(r8,xsize=7)
	dum = widget_label(r8,value=' b=')
	b = widget_text(r8,xsize=7)
	dum = widget_label(r8,value=' phi=')
	phi = widget_text(r8,xsize=7)
	dum = widget_label(r8,value=' ')
	fitellipse = widget_button(r8,value='FIT ELLIPSES')
	wrtellipse = widget_button(r8,value='SAVE ELLIPSES')

	r9=widget_base(right,/row)
	dum = widget_label(r9,value='LOW ')
	ilow_slider = widget_slider(r9,xsize=486,min=0,max=486,value=0,/drag)
	r10=widget_base(right,/row)
	dum = widget_label(r10,value='HIGH')
	ihigh_slider = widget_slider(r10,xsize=486,min=0,max=486,value=486,/drag)
	r11=widget_base(right,/row)
	dum = widget_label(r11,value='OUTL ')
	out_slider = widget_slider(r11,xsize=486,min=0,max=486,value=0,/drag)

	
	;build u structure
	id={base:base,draw1:draw1,draw2:draw2,draw3:draw3,$
		shotid:shotid,module:module,get:get,quit:quit,print:print,stopbutton:stopbutton,t1:t1,t2:t2,$
		gain_slider:gain_slider,bin_slider:bin_slider,brkpts:brkpts,$
		ar_helikefit:ar_helikefit,ar_hlikefit:ar_hlikefit,ca_helikefit:ca_helikefit,ca_hlikefit:ca_hlikefit,$
		savepath:savepath,treesave:treesave,loaddata:loaddata,$
		infopath:infopath,treeinfo:treeinfo,loadinfo:loadinfo,$
		precision:precision,splot:splot,fitspec:fitspec,fitname:fitname,ibad:ibad,message:message,$
		infoplot:infoplot,$
		phifit:phifit,fit_lya1:fit_lya1,fit_lya2:fit_lya2,fit_T:fit_T,fit_J:fit_J,fit_w:fit_w,fit_x:fit_x,fit_y:fit_y,fit_z:fit_z,$
		xo:xo,yo:yo,a:a,b:b,phi:phi,fitellipse:fitellipse,wrtellipse:wrtellipse,$
		ilow_slider:ilow_slider,ihigh_slider:ihigh_slider,out_slider:out_slider,vert_slider:vert_slider,spec_slider:spec_slider,flow_slider:flow_slider,$
			fhigh_slider:fhigh_slider,left_slider_button:left_slider_button, right_slider_button:right_slider_button}
	stat={ps:0,dp:0,bs:1,splot:0,iplot:0,itree:1,stree:1,lya1:0,lya2:0,J:0,T:0,w:0,x:0,y:0,z:0,focus:0,name:'',lamgood:[0,0,0,0],order:0,fitcase:0,xsize:487*1.5,ysize:195*1.5}

	;set defaults	
	IF NOT keyword_set(shot) THEN shot=1070830020
	module=2
	IF module EQ 4 then stat.fitcase=1 ELSE stat.fitcase=0
	t1=1.05000
	t2=1.45000
	
	IF stat.stree EQ 0 THEN $
		psave='/home/'+user+'/idl/hirexsr/calib/ellipse_'+num2str(module,1)+'_'+num2str(shot,1)+'_'+num2str(int(t1*1000),1)+'_'+num2str(int(t2*1000),1)+'.dat' ELSE $
		psave='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)
	IF stat.itree EQ 0 THEN $
		pinfo='/home/'+user+'/idl/genie/data/info/hirexsr/hirexsr_'+num2str(module,1)+'.info' ELSE  $
		pinfo='\SPECTROSCOPY::TOP.HIREXSR.INFO.MOD'+num2str(module,1)
	path={save:psave,info:pinfo}

	npeaks=4
	mxp=40
	fits={resid:fltarr(487,195)+1,peaks:fltarr(487,npeaks),err:fltarr(487,npeaks),lambda:fltarr(npeaks)}
	spec={a0:fltarr(487,mxp),a1:fltarr(487,mxp),a2:fltarr(487,mxp),a3:fltarr(487),label:strarr(mxp)}
	ellipse={i:0,tilt:0,param:fltarr(5,npeaks),good:intarr(487)+1,bad:intarr(487),out:intarr(487,npeaks),mu:fltarr(npeaks),sigma:fltarr(npeaks)}
	info=fltarr(487,npeaks)

	u={id:id,stat:stat,path:path,shot:shot,module:module,t1:t1,t2:t2,image:fltarr(487,195)+1,rawimage:fltarr(487,195)+1,brkpts:ptr_new([0],/allocate),fits:fits,ellipse:ellipse,info:info,spec:spec}
	widget_control,base,/realize
	widget_control,base,set_uvalue=u
	widget_control,id.shotid,set_value=strtrim(u.shot,2)

	widget_control,id.module,set_value=strtrim(u.module,2)
	widget_control,id.t1,set_value=num2str(t1,dp=2)
	widget_control,id.t2,set_value=num2str(t2,dp=2)
	widget_control,id.savepath,set_value=u.path.save
	widget_control,id.infopath,set_value=u.path.info
	index=0
	widget_control,u.id.fit_w,set_button=1
	widget_control,u.id.treeinfo,set_button=u.stat.itree
	widget_control,u.id.treesave,set_button=u.stat.stree
	widget_control,u.id.infoplot,set_button=u.stat.iplot
	CASE u.stat.fitcase OF 
		0 : widget_control,u.id.ar_helikefit,set_button=1
		1 : widget_control,u.id.ar_hlikefit,set_button=1
		2 : widget_control,u.id.ca_helikefit,set_button=1
		3 : widget_control,u.id.ca_hlikefit,set_button=1
	ENDCASE
	widget_control,u.id.xo,set_value=num2str(u.ellipse.param[0,index],dp=2)
	widget_control,u.id.yo,set_value=num2str(u.ellipse.param[1,index],dp=2)
	widget_control,u.id.a,set_value=num2str(u.ellipse.param[2,index],dp=2)
	widget_control,u.id.b,set_value=num2str(u.ellipse.param[3,index],dp=2)
	widget_control,u.id.phi,set_value=num2str(u.ellipse.param[4,index],dp=2)
        calib_load_user_name, u ; write the current user name

	!except=0
	xmanager,'w_hirexsr_calib',base
END
