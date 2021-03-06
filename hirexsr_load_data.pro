;+
;NAME:
;	HIREXSR_LOAD_WAVELENGTHS
;PURPOSE:
;	This procedure load the wavelength table from the tree
;
;CALLING SEQUENCE
;	HIREXSR_LOAD_WAVELENGTHS,shot,lam_o,z,label
;
;INPUTS:
;	shot	LONG	shot number
;	
;OUPUTS:
;	lam_o	FLTARR	[n_lines] rest wavelength of major transition [Ang]
;	z	INTARR	[n_lines] atomic number of element
;	label	STRARR	[n_lines] indexing label used to identify transition (usually Gabriel notation)
;
;RESTRICTIONS:
;	Literature searches have been done for lines in this table to provide the most accurate
;	EXPERIMENTALLY measured rest wavelengths
;
;MODIFICATION HISTORY:
;	Written by: 	M.L. Reinke 6/8/2010
;
;-


PRO hirexsr_load_wavelengths,shot,lam_o,z,label
	mdsopen,'spectroscopy',shot
	lam_o=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB:LAM_O')
	z=mdsvalue('dim_of(\SPECTROSCOPY::TOP.HIREXSR.CALIB:LAM_O,0)')
	label=mdsvalue('dim_of(\SPECTROSCOPY::TOP.HIREXSR.CALIB:LAM_O,1)')
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_WRITE_WAVELENGTHS
;PURPOSE:
;	This procedure writes the wavelength table to the tree from a IDL save file
;
;CALLING SEQUENCE
;	HIREXSR_WRITE_WAVELENGTHS,shot
;
;INPUTS:
;	shot	LONG	shot number
;
;KEYWORD PARAMETERS:
;	convert		/convert overwrites the DAT file with data from the CSV file
;	
;OUPUTS:
;	Data is written to \SPECTROSCOPY::TOP.HIREXSR.CALIB:LAM_O node in the tree
;
;RESTRICTIONS:
;	This currently loads the lam_o,z,label arrays from /usr/local/cmod/idl/HIREXSR/hirexsr_wavelengths.dat
;	or converts ~/HIREXSR/hirexsr_wavelengths.csv to ~/HIREXSR/hirexsr_wavelengths.dat and writes those values to the tree.
;	The distributed .dat file should then be updated.
;
;MODIFICATION HISTORY:
;	Written by: 	M.L. Reinke 6/8/2010
;	12/20/10:	M.L. Reinke - changed the CSV and DAT path away from /home/mlreinke/...
;
;-


PRO hirexsr_write_wavelengths,shot,convert=convert,local=local
	csv_path='/home/'+logname()+'/HIREXSR/hirexsr_wavelengths.csv'
	IF keyword_set(local) OR keyword_set(convert) THEN dat_path='/home/'+logname()+'/HIREXSR/hirexsr_wavelengths.dat' ELSE $
		dat_path='/usr/local/cmod/idl/HIREXSR/hirexsr_wavelengths.dat'
	IF keyword_set(convert) THEN csv_convert,csv_path,dat_path
	restore,dat_path
	mdsopen,'spectroscopy',shot
	mdsput,'\SPECTROSCOPY::TOP.HIREXSR.CALIB:LAM_O','build_signal(build_with_units($1,"Ang"),*,build_with_units($2,""),build_with_units($3,""))',lam_o,z,label
	mdsclose,'spectroscopy',shot
END

; +
;NAME: 
;	HIREXSR_GET_TIME
;
;MODIFICATION HISTORY:
;	Written by	A. Ince-Cushman (2007)
;
; -

FUNCTION hirexsr_get_time, shot_num, quiet = quiet

mdsopen, "spectroscopy", shot_num, quiet=quiet
  n_frames   = mdsvalue ('\spectroscopy::top.hirex_sr.n_frames'  , status =status1, /quiet)
  exp_period = mdsvalue ('\spectroscopy::top.hirex_sr.exp_period', status =status2, /quiet)
  exp_time   = mdsvalue ('\spectroscopy::top.hirex_sr.exp_time'  , status =status3, /quiet)
  t_offset   = mdsvalue ('\SPECTROSCOPY::TOP.X_RAY_PHA:GAS_DECODER.CHANNEL_3:P1' , /quiet, status =status4)
  ; this t_offset is from the amount of time before 'pulse that the data acq begins: default is -0.1s
mdsclose, quiet=quiet

if status1*status2*status3*status4 then begin
    time = findgen(n_frames)*exp_period + 0.5*exp_time+t_offset
    return, time
endif else begin
    if not keyword_set(quiet) then print, 'Problem getting timing information'
    return, -1
endelse

END

;+
;NAME:
;	HIREXSR_LOAD_IMAGE
;	
;PURPOSE:
;	This procedure loads the raw image file from the tree
;
;CALLING SEQUENCE:
;	HIREXSR_LOAD_IMAGE,shot,module,image,t
;
;INPUTS:
;	shot:	LONG	shot number
;	module	FLOAT	module number (1-4)
;
;KEYWORD PARAMETERS:
;	noimage	/noimage will not load the image, allowing just the time scale to be loaded
;
;OUTPUTS:
;	image	FLTARR	[487,195,n_frame] of the PILATUS camera images
;	t	FLTARR	[n_frame] of the time points
;
;OPTIONAL OUTPUTS:
;	break	FLTARR	of break points where the spectral units should be found (adjust in W_HIREXSR_CALIB)
;
;MODFICATION HISTORY:
;	Written by:	M.L. Reinke 1/10
;	3/25/11		M.L. Reinke - added the /noimage keyword
;	9/4/11		M.L. Reinke - modified the loading of the time vector to use HIREXSR_GET_TIMES rather then use the dimension of
;				      of the raw data.  Should speed up load times.
;
;-

PRO hirexsr_load_image,shot,module,image,t,break=break,noimage=noimage
	mdsopen,'spectroscopy', shot
	IF NOT keyword_set(noimage) THEN image=float(mdsvalue('\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD'+num2str(module,1)))
	break=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+'.FITS.SETUP:BREAK',/quiet,status=status)
	IF NOT status THEN break=0
	mdsclose,'spectroscopy',shot
	t=hirexsr_get_time(shot)

END

;+
;NAME:
;	HIREXSR_WRITE_BREAK
;
;-

PRO hirexsr_write_break,shot,module,break
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+'.FITS.SETUP:BREAK'
	mdsput,path,'build_with_units($,"pix")',int(break)
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_BREAK
;
;-

PRO hirexsr_load_break,shot,module,break
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+'.FITS.SETUP:BREAK'
	break=mdsvalue(path,/quiet,status=status)
	IF NOT status THEN break=0
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_DT
;
;-

PRO hirexsr_load_dt,shot,dt
	mdsopen,'spectroscopy',shot
	dt=mdsvalue('\SPECTROSCOPY::TOP.HIREX_SR:EXP_TIME')
	mdsclose,'spectroscopy',shot
END


;+
;NAME:
;	HIREXSR_WRITE_CALIBFITS
;
;-

PRO hirexsr_write_calibfits,shot,module,times,image,label,nphot,offset,peaks,width,resid,break,double,ifit_f
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)
	mdsopen,'spectroscopy',shot
	mdsput,'\SPECTROSCOPY::TOP.HIREXSR.CALIB:SHOT','build_with_units($,"")',shot
	mdsput,'\SPECTROSCOPY::TOP.HIREXSR.CALIB:TIMES','build_with_units($,"sec")',times
	mdsput,path+':IMAGE','build_with_units($,"counts")',image
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+'.FITS'
	mdsput,path+':LABEL','build_with_units($,"")',label
	mdsput,path+':NPHOT','build_with_units($,"counts")',nphot
	mdsput,path+':OFFSET','build_with_units($,"counts")',offset
	mdsput,path+':PEAKS','build_with_units($,"pixel")',peaks
	mdsput,path+':WIDTH','build_with_units($,"pixel")',width
	mdsput,path+':RESID','build_with_units($,"pixel")',resid
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+'.FITS.SETUP:'
	mdsput,path+'BREAK','build_with_units($,"")',int(break)
	mdsput,path+'DOUBLE','build_with_units($,"")',int(double)
	mdsput,path+'IFIT_F','build_with_units($,"")',ifit_f
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_CALIBFITS
;
;-

PRO hirexsr_load_calibfits,shot,module,times,image,label,nphot,offset,peaks,width,resid,break,double,ifit_f,name
	mdsopen,'spectroscopy',shot
	times=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB:TIMES')
	image=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module)+':IMAGE')
	label=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module)+'.FITS:LABEL')
	nphot=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module)+'.FITS:NPHOT')
	offset=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module)+'.FITS:OFFSET')
	peaks=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module)+'.FITS:PEAKS')
	width=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module)+'.FITS:WIDTH')
	resid=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module)+'.FITS:RESID')
	break=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module)+'.FITS.SETUP:BREAK')
	double=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module)+'.FITS.SETUP:DOUBLE')
	ifit_f=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module)+'.FITS.SETUP:IFIT_F')
	mdsclose,'spectroscopy',shot
END


;+
;NAME:
;	HIREXSR_WRITE_ELLIPSE
;
;-

PRO hirexsr_write_ellipse,shot,module,coefs,ifit_e,outl,lambda,mu,sigma,bad
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+'.ELLIPSE'
	mdsput,path+':COEFS','build_with_units($,"")',coefs
	mdsput,path+':IFIT_E','build_with_units($,"pixel")',ifit_e
	mdsput,path+':OUTL','build_with_units($,"pixel")',outl
	mdsput,path+':LAMBDA','build_with_units($,"Ang")',lambda
	mdsput,path+':MU','build_with_units($,"Ang")',mu
	mdsput,path+':SIGMA','build_with_units($,"Ang")',sigma
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+'.FITS.SETUP'
	mdsput,path+':BAD','build_with_units($,"")',bad
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_ELLIPSE
;
;-

PRO hirexsr_load_ellipse,shot,module,coefs,ifit_e,outl,lambda,mu,sigma,bad
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+'.ELLIPSE'
	coefs=mdsvalue(path+':COEFS')
	ifit_e=mdsvalue(path+':IFIT_E')
	outl=mdsvalue(path+':OUTL')
	lambda=mdsvalue(path+':LAMBDA')
	mu=mdsvalue(path+':MU')
	sigma=mdsvalue(path+':SIGMA')
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+'.FITS.SETUP'
	bad=mdsvalue(path+':BAD')
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_ECOEFS
;
;PURPOSE:
;	This procedure loads the elliptical fit coefficients from a locked mode calibration 
;	and the corresponding center wavelengths.  These can be used to generate a wavelength map
;
;CALLING SEQUENCE:
;	HIREXSR_LOAD_ECOEFS,shot,module,coefs,lam_o
;
;INPUTS:
;	shot:	LONG	shot number
;	module	FLOAT	module number (1-4)
;
;OUTPUTS:
;	coefs:	FLTARR	[5,n_lines] of the elliptic curve coefficients
;	lam_o:	FLTARR	[n_lines] of the line centers 
;
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 3/10
;
;-

PRO hirexsr_load_ecoefs,shot,module,coefs,lam_o
	mdsopen,'spectroscopy',shot
	coefs=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+'.ELLIPSE:COEFS')
	lam_o=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+'.ELLIPSE:LAMBDA')
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_COEFS2LAM
;
;PURPOSE:
;	This function takes the elliptic curve coefficients and line centers and generates
;	a 2D map of wavelength values
;
;CALLING SEQUENCE:
;	result=HIREXSR_COEFS2LAM(coefs,lam_o)
;
;INPUTS:
;	coefs:	FLTARR	[5,n_lines] of the elliptic curve coefficients
;	lam_o:	FLTARR	[n_lines] of the line centers 
;
;OPTIONAL INPUTS:
;	order:	INT	order of the polynomial fit to calculate wavelength over the row DEFAULT: 1
;	nx:	INT	number of pixels in the non-wavelenth direction DEFAULT: 487
;	ny:	INT	number of pixesl in the wavelength direction DEFAULT: 195
;
;KEYWORD PARAMETERS:
;	plot:	/plot will plot the results of the polynomial fit
;	debug:	/debug will stop the code before the RETURN
;
;OUTPUS:
;	result:	FLTARR 	[nx,ny] of the wavelength values at each pixel [units of lam_o]
;
;PROCEDURE:
;	This uses the data loaded from HIREXSR_LOAD_COEFS and the function EQ_ELLIPSE
;	to generate the data.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 3/10
;	7/28/16		M.L. Reinke - fixed a bug in the /plot feature
;                                     for dispersion plot
;-

FUNCTION hirexsr_coefs2lam,coefs,lam_o,order=order,nx=nx,ny=ny,plot=plot,debug=debug,good=good
	IF NOT keyword_set(order) THEN order=2
	IF NOT keyword_set(nx) THEN nx=487
	IF NOT keyword_set(ny) THEN ny=195
	lam=fltarr(nx,ny)
	nlines=n(lam_o)+1
	pix_o=fltarr(nlines)
	pix=findgen(ny)
	IF NOT keyword_set(good) THEN good=intarr(nlines)+1
	tmp=where(good EQ 1)
	FOR i=0,nx-1 DO BEGIN
		FOR j=0,nlines-1 DO pix_o[j]=eq_ellipse(float(i),coefs[*,j])
		poly_coefs=poly_fit(pix_o[tmp],lam_o[tmp],order)
		lam[i,*]=poly(pix,poly_coefs)
		IF keyword_set(plot) THEN BEGIN
			openwin,0			;plot the residual of the poly fit to the points in units of apparent velocity
			IF i EQ 0 THEN plot,[0],[0],xr=[0,nlines-1],yr=[-10.0,10.0],ytit=n2g('Delta')+'V [km/s]',/ysty,/xsty
			oplot,(lam_o-poly(pix_o,poly_coefs))/lam_o*3.0e8/1.0e3,psym=8
		ENDIF
	ENDFOR
	IF keyword_set(plot) THEN BEGIN
		openwin,1				;plot the dispersion Ang/Pix
		plot,[0],[0],xr=[0,ny-1],yr=[2.0,4.0]*1.0e-4,ytit='Dispersion [Ang/Pix]',/xsty,/ysty,xtit='Pixel #'
		FOR i=0,nx-1 DO oplot, pix,deriv(pix,lam[i,*])
	ENDIF
	IF keyword_set(debug) THEN stop
	RETURN,lam
	
END

;+
;NAME:
;	HIREXSR_WRITE_LAMBDA
;
;-

PRO hirexsr_write_lambda,shot,module,order=order,nx=nx,ny=ny,lambda=lambda,good=good
	IF NOT keyword_set(lambda) THEN BEGIN
		hirexsr_load_ecoefs,shot,module,coefs,lam_o
		lambda=hirexsr_coefs2lam(coefs,lam_o,order=order,nx=nx,ny=ny,good=good)
	ENDIF
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)
	mdsput,path+':LAMBDA','build_with_units($,"angstroms")',lambda
	mdsopen,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_LAMBDA
;
;-

PRO hirexsr_load_lambda,shot,module,lambda
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)
	lambda=mdsvalue(path+':LAMBDA')
	mdsclose,'spectroscopy',shot	
END

;+
;NAME:
;	HIREXSR_WRITE_CALIB_INST
;
;PURPOSE:
;	This procedure writes the instrumental shift and width [Ang] calculated for
;	each pixel on a given detector module.  Tree nodes are created for the given shot if not present
;
;WRITTEN BY:
;	M.L. Reinke - 9/4/11
;
;-

PRO hirexsr_write_calib_inst,shot,module,iwidth,ishift
	hirexsr_add_inst_nodes,shot,/quiet				;add the instrumental nodes if not already there
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)
	mdsput,path+':INST','build_signal(build_with_units($1,"Ang"),*,build_with_units($2,"Ang"))',iwidth,ishift
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_CALIB_INST
;
;PURPOSE:
;	This procedure loads the instrumental shift and width [Ang] calculated for
;	each pixel on a given detector module
;
;WRITTEN BY:
;	M.L. Reinke - 9/4/11
;
;-

PRO hirexsr_load_calib_inst,shot,module,iwidth,ishift,status=status
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)
	iwidth=mdsvalue(path+':INST',/quiet,status=status)
	ishift=mdsvalue('dim_of('+path+':INST,0)',/quiet,status=status)
	mdsclose,'spectroscopy',shot
END	

;+
;NAME:
;	HIREXSR_READ_TREEINFO
;
;PURPOSE:
;	This function reads INFO file data from the HIREXSR tree
;
;CALLING SEQUENCE
;	result=HIREXSR_READ_TREEINFO(shot,moldule)
;
;INPUTS:
;	shot:	LONG	shot number
;	module	FLOAT	module number (1-4)
;
;OUTPUTS:
;	result:	STRUC	same as the output of HIREXSR_LOAD_INFO
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 3/10
;
;-

FUNCTION hirexsr_read_treeinfo,shot,module
	path='\SPECTROSCOPY::TOP.HIREXSR.INFO.MOD'+num2str(module,1)
	mdsopen,'spectroscopy',shot
	m_rad=mdsvalue(path+'.MIRROR:RAD')
	m_rot=mdsvalue(path+'.MIRROR:ROT')
	m_vec=mdsvalue(path+'.MIRROR:VEC')
	m_size=mdsvalue(path+'.MIRROR:SIZE')
	m_iref=mdsvalue(path+'.MIRROR.BRAGG:IREF')
	m_rwid=mdsvalue(path+'.MIRROR.BRAGG:RWID')
	m_twod=mdsvalue(path+'.MIRROR.BRAGG:TWOD')
	d_xi=mdsvalue(path+'.DET:DET_XI')
	d_zeta=mdsvalue(path+'.DET:DET_ZETA')
	n_xi=mdsvalue(path+'.DET:N_XI')
	n_zeta=mdsvalue(path+'.DET:N_ZETA')
	x0=mdsvalue(path+'.DET:X0')
	x1=mdsvalue(path+'.DET:X1')
	x2=mdsvalue(path+'.DET:X2')
	xi_ch=mdsvalue(path+'.DET:XI_CH')
	xi_o=mdsvalue(path+'.DET:XI_O')
	zeta_ch=mdsvalue(path+'.DET:ZETA_CH')
	zeta_o=mdsvalue(path+'.DET:ZETA_O')
	mdsclose,'spectroscopy',shot

	;setup mirror structure
	bragg=create_struct('twod',m_twod,'iref',m_iref,'rwid',m_rwid)
	m_struc=create_struct('vec',m_vec,'rot',m_rot,'size',m_size,'bragg',bragg,'rad',m_rad)

	;setup detector structure
	det_size=[d_xi,d_zeta]
	cntr=0L
	num=n_xi*n_zeta
	IF n_xi EQ 0 THEN num=n_zeta
	IF n_zeta EQ 0 THEN num=n_xi
	xi=fltarr(num)
	zeta=fltarr(num)
	IF n_zeta EQ 0 THEN n_zeta=1
	IF n_xi EQ 0 THEN n_xi = 1
	FOR i=0,n_zeta-1 DO BEGIN
		FOR j=0,n_xi-1 DO BEGIN
			xi[cntr]=xi_ch*(j-xi_o)
			zeta[cntr]=zeta_ch*(i-zeta_o)
			cntr+=1
		ENDFOR
	ENDFOR
	det_struc=create_struct('x0',x0,'x1',x1,'x2',x2,'xi',xi, 'zeta',zeta,'size',det_size,'n_xi',int(n_xi),'n_zeta',int(n_zeta))	
	name='hirexsr_0'+num2str(module,1)
	type='spherical'
	author=logname()
	output={name:name, m:m_struc, det:det_struc, type:type, author:author[0]}

	RETURN,output
END

;+
;NAME:
;	HIREXSR_CONSTRAINED_INFO
;
;-

PRO hirexsr_constrained_info,shot,i1,i3,debug=debug
	hirexsr_load_morder,shot,morder
	i2=hirexsr_read_treeinfo(shot,morder[1])
	i1=hirexsr_read_treeinfo(shot,morder[0])
	i3=hirexsr_read_treeinfo(shot,morder[2])

	;use CAD drawing constraints to adjust the i1 and i3 x0,x1,x2 vectors based on the position of i2
	x0=i2.det.x0
	x1=i2.det.x1
	x2=i2.det.x2
	zhat=(x1-x0)/sqrt(total((x1-x0)*(x1-x0)))
	xhat=(x2-x0)/sqrt(total((x2-x0)*(x2-x0)))
	nhat=crossp(xhat,zhat)
	i1.det.x0=i2.det.x0-0.34e-3*xhat+88.82e-3*zhat+0.2e-3*nhat
	i1.det.x2=i2.det.x0+32.89e-3*xhat+93.06e-3*zhat+0.35e-3*nhat
	i1.det.x1=i2.det.x0-10.90e-3*xhat+171.33e-3*zhat+6.29e-3*nhat

	IF keyword_set(debug) THEN stop

	x0=i2.det.x0
	x1=i2.det.x1
	x2=i2.det.x2
	zhat=(x1-x0)/sqrt(total((x1-x0)*(x1-x0)))
	xhat=(x2-x0)/sqrt(total((x2-x0)*(x2-x0)))
	nhat=crossp(xhat,zhat)
	i3.det.x0=i2.det.x0-10.90e-3*xhat-87.93e-3*zhat+6.29e-3*nhat
	i3.det.x2=i2.det.x0+22.33e-3*xhat-92.17e-3*zhat+6.45e-3*nhat
	i3.det.x1=i2.det.x0-0.34e-3*xhat-5.42e-3*zhat+0.2e-3*nhat

	IF keyword_set(debug) THEN stop
END


;+
;NAME:
;	HIREXSR_WRITE_POS
;
;-

PRO hirexsr_write_pos,shot,module,info=info,pos=pos
	
	IF NOT keyword_set(pos) THEN BEGIN
		IF NOT keyword_set(info) THEN info=hirexsr_read_treeinfo(shot,module)
		pos=genpos_spherical2pos(info,etendue=u)
		newpos=fltarr(4,info.det.n_zeta,info.det.n_xi)
		FOR i=0,info.det.n_xi-1 DO newpos[*,*,i]=pos[*,indgen(info.det.n_zeta)*float(info.det.n_xi)+i]
		FOR i=0,3 DO newpos[i,*,*]=rotate(reform(newpos[i,*,*]),2) 	;this converts from the (xi,zeta) coordinates to the customary (i,j) of the collected image
	ENDIF ELSE newpos=pos
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+':POS'
	mdsput,path,'build_with_units($,"")',newpos
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_POS
;
;-

PRO hirexsr_load_pos,shot,module,pos
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+':POS'
	pos=mdsvalue(path)
	mdsclose,'spectroscopy',shot
END


;+
;NAME:
;	HIREXSR_WRITE_ETENDUE
;
;-

PRO hirexsr_write_etendue,shot,module,u=u							
	;eventually replace w/ ray tracing output
	IF NOT keyword_set(u) THEN BEGIN
		IF NOT keyword_set(info) THEN info=hirexsr_read_treeinfo(shot,module)
		pos=genpos_spherical2pos(info,etendue=u)
		newu=fltarr(info.det.n_zeta,info.det.n_xi)
		FOR i=0,info.det.n_xi-1 DO newu[*,i]=u[indgen(info.det.n_zeta)*float(info.det.n_xi)+i]
		newu=rotate(newu,2)						;this converts from the (xi,zeta) coordinates to the customary (i,j) of the collected image
	ENDIF ELSE newu=u
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+':U'
	mdsput,path,'build_with_units($,"m^2-str")',newu
	mdsclose,'spectroscopy',shot
END

;+
;NAME: 
;	HIREXSR_LOAD_ETENDUE
;
;-

PRO hirexsr_load_etendue,shot,module,u
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+':U'
	u=mdsvalue(path)
	mdsclose,'spectroscopy',shot
END


;+
;NAME:
;	HIREXSR_WRITE_TRANS
;
;-

PRO hirexsr_write_trans,shot,module,trans=trans
	IF NOT keyword_set(trans) THEN BEGIN
		trans=fltarr(487,195)+1.0
	ENDIF
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+':TRANS'
	mdsput,path,'build_with_units($,"")',trans					;eventually replace w/ ray tracing output
	mdsclose,'spectroscopy',shot

END

;+
;NAME:
;	HIREXSR_LOAD_TRANS
;
;-

PRO hirexsr_load_trans,shot,module,trans
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+':TRANS'
	trans=mdsvalue(path)
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_WRITE_WHITE
;
;-

PRO hirexsr_write_white,shot,module,white=white
	IF NOT keyword_set(white) THEN BEGIN
		white=fltarr(487,195)+1.0
	ENDIF
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+':WHTFLD'
	mdsput,path,'build_with_units($,"")',white					;eventually replace w/ measured data
	mdsclose,'spectroscopy',shot

END

;+
;NAME:
;	HIREXSR_LOAD_WHITE
;
;-

PRO hirexsr_load_white,shot,module,white
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+':WHTFLD'
	white=mdsvalue(path)
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_CONST
;
;MODIFICATION HISTORY:
;	8/10/11		M.L. Reinke - noticed that const did NOT include the dt, included this now
;	11/5/11		M.L. Reinke - added the /wf keyword.  Default will be to use wf=1.0, use /wf to force loading from tree
;-

PRO hirexsr_load_const,shot,module,const,wf=wf
	hirexsr_load_dt,shot,dt
	hirexsr_load_lambda,shot,module,lambda
	hirexsr_load_trans,shot,module,trans
	hirexsr_load_etendue,shot,module,u
	IF keyword_set(wf) THEN hirexsr_load_white,shot,module,white ELSE white=trans*0.0+1.0
	x=size(lambda)
	dlam=fltarr(x[1],x[2])
	FOR i=0,x[1]-1 DO dlam[i,*]=fltarr(x[2])+(lambda[i,x[2]-1]-lambda[i,0])/(x[2]-1.0)
	const=4.0*!pi/(trans*u*white*dlam*dt)*1.0e-17						;adjust absolute value
END


;+
;NAME:
;	HIREXSR_PTR_IMAGE
;
;PURPOSE:
;	This procedure loads each frame of the HIREXSR data into an element of a PTRARR 
;	and calculates the wavelength map.  For the 3 module system, the images are combined
;	into one image.
;
;CALLING SEQUENCE:
;	HIREXSR_PTR_IMAGE,shot,cnts,lambda,t
;
;INPUTS:
;	shot:	LONG	shot number
;
;OPTIONAL INPUTS:
;	morder:	INTARR	[3] ([1] if /h) of the order to assemble the images DEFAULT: [1,2,3] or [4] when /h
;
;KEYWORD PARAMETERS:
;	h	/h will load a single module with the default being module 4
;	noimage	/noimage will avoid skip loading the image data and cnts will be an empty pointer.
;		This is useful in FITSPEC2TREE when the POS/CONST/U/INST data is useful but no the raw data
;	wf	/wf will load white-field stored in tree else wf=1.0 is assumed 
;OUTPUTS:
;	cnts:	PTRARR	[n_time] of pointers for each frame
;	lambda:	FLTARR	[3*nx,ny] of the wavelength values [Ang] ([nx,ny] for /h)
;	t:	FLTARR	[n_time] of the time points [sec]
;
;OPTIONAL OUTPUTS:
;	pos	FLTARR	[4,3*nx,ny] of the pos vectors for each pixel ([nx,ny] for /h)
;	u	FLTARR	[4,3*nx,ny] of the etendue for each pixel ([nx,ny] for /h)
;	const	FLTARR	[4,3*nx,ny] of the brightness constant for each pixel ([nx,ny] for /h)
;	iwidth  FLTARR	[4,3*nx,ny] of the instrumental width [Ang] for each pixel ([nx,ny] for /h) DEFAULTS=0.0 if not stored
;	ishift	FLTARR	[4,3*nx,ny] of the instrumental shift [Ang] for each pixel ([nx,ny] for /h) DEFAULTS=0.0 if not stored
;
;MODIFICATION HISTORY:
;	Written by: 	M.L. Reinke 4/5/2010
;	8/17/10		M.L. Reinke - modified to allowing loading of u, pos and const when using /h
;	9/4/11		M.L. Reinke - modified to allow loading of instrumental data & updated documentation for optional outputs	
;-

PRO hirexsr_ptr_image,shot,cnts,lambda,t,pos=pos,u=u,const=const,iwidth=iwidth,ishift=ishift,h=h,morder=morder,noimage=noimage,wf=wf
	IF NOT keyword_set(morder) THEN BEGIN
		mdsopen,'spectroscopy',shot
		IF keyword_set(h) THEN morder=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE:MORDER') ELSE $
			morder=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE:MORDER')
	ENDIF
	IF NOT keyword_set(h) THEN BEGIN
		;load raw image data and ellipse coefs fits
		;module1
		hirexsr_load_image,shot,morder[0],image1,t,noimage=noimage
 		hirexsr_load_lambda,shot,morder[0],lambda1
		hirexsr_load_etendue,shot,morder[0],u1
		hirexsr_load_pos,shot,morder[0],pos1
		hirexsr_load_const,shot,morder[0],const1,wf=wf
		hirexsr_load_calib_inst,shot,morder[0],iw1,is1,status=status1
		;module 2
		hirexsr_load_image,shot,morder[1],image2,t,noimage=noimage
 		hirexsr_load_lambda,shot,morder[1],lambda2
		hirexsr_load_etendue,shot,morder[1],u2
		hirexsr_load_pos,shot,morder[1],pos2
		hirexsr_load_const,shot,morder[1],const2,wf=wf
		hirexsr_load_calib_inst,shot,morder[1],iw2,is2,status=status2
		;module 3
		hirexsr_load_image,shot,morder[2],image3,t,noimage=noimage
 		hirexsr_load_lambda,shot,morder[2],lambda3
		hirexsr_load_etendue,shot,morder[2],u3
		hirexsr_load_pos,shot,morder[2],pos3
		hirexsr_load_const,shot,morder[2],const3,wf=wf
		hirexsr_load_calib_inst,shot,morder[2],iw3,is3,status=status3
		;combine
		lambda=[lambda1,lambda2,lambda3]
		u=[u1,u2,u3]
		pos=[[pos1],[pos2],[pos3]]
		const=[const1,const2,const3]
		IF status1 AND status2 AND status3 THEN BEGIN
			iwidth=[iw1,iw2,iw3]
			ishift=[is1,is2,is3]
		ENDIF ELSE BEGIN
			iwidth=u*0.0		;default both to 0.0
			ishift=u*0.0
		ENDELSE
	ENDIF ELSE BEGIN
		hirexsr_load_image,shot,morder[0],image,t,noimage=noimage
		hirexsr_load_lambda,shot,morder[0],lambda
		hirexsr_load_etendue,shot,morder[0],u
		hirexsr_load_pos,shot,morder[0],pos
		hirexsr_load_const,shot,morder[0],const,wf=wf
		hirexsr_load_calib_inst,shot,morder[0],iwidth,ishift,status=status
		IF NOT status THEN BEGIN
			iwidth=u*0.0		;default both to 0.0
			ishift=u*0.0
		ENDIF	
	ENDELSE

	;form PTRARR for the image data
	n_time=n(t)+1
	cnts=ptrarr(n_time,/allocate_heap)
	IF NOT keyword_set(noimage) THEN BEGIN
		IF keyword_set(h) THEN BEGIN
			FOR i=0,n_time-1 DO *cnts[i]=image[*,*,i]
		ENDIF ELSE BEGIN
			FOR i=0,n_time-1 DO *cnts[i]=[image1[*,*,i],image2[*,*,i],image3[*,*,i]]
		ENDELSE
	ENDIF

	
END




;+
;NAME:
;	HIREXSR_LOAD_INFO
;
;PURPOSE:
;	This function loads the INFO file for each of the HIREXSR modules
;
;CALLING SEQUENCE:
;	result=HIREXSR_LOAD_INFO(det)
;
;INPUTS:
;	module	INT 	detector module number (1-4)
;
;OPTIONAL INPUTS:
;	shot:	LONG	shot number of the tree to load with /tree is invoked
;
;KEYWORD PARAMETERS:
;	ca:	/ca (and not /tree) will load the info file for viewing He-like Ca
;	tree:	/tree will load the info from the tre usinG HRIEXSR_READ_TREEINFO
;
;OUTPUTS:
;	result:	STRUC	INFO structure for the module (see GENPOS_SPHERICAL INFO)
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 1/10
;	4/7/10:		M.L. Reinke - added optional inputs to allow this function to call data from the tree
;
;-

FUNCTION hirexsr_load_info,module,ca=ca,shot=shot,tree=tree
	path='/home/mlreinke/idl/genie/data/info/hirexsr/'
	IF keyword_set(ca) THEN file='hirexsr_4_ca.info' ELSE file='hirexsr_'+num2str(module,1)+'.info'
	IF NOT keyword_set(tree) THEN info=genpos_spherical_info(path+file) ELSE info=hirexsr_read_treeinfo(shot,module)
	RETURN,info
END


;+
;NAME:
;	HIREXSR_LOAD_INFO2TREE
;
;PURPOSE:
;	This procedure loads data from an spherical INFO file and writes it to the tree'
;
;CALLING SEQUENCE:
;	HIREXSR_LOAD_INFO2TREE,shot,module
;
;INPUTS:
;	shot:	LONG	shot number
;	module	FLOAT	module number (1-4)
;
;OPTIONAL INPUTS:
;	info:	STRUC	info structure DEFAULT: HIREXSR_LOAD_INFO(module)
;	
;KEYWORD PARAMETERS:
;	ca:	/ca will load the INFO file to look at He-like Ca
;
;OUTPUTS:
;	Data written to the \SPECTROSCOPY::TOP.HIREXSR.INFO.MOD# nodes
;
;MODIFICATION HISTORY:
;	Written by: 	M.L. Reinke 3/10
;
;-	

PRO hirexsr_load_info2tree,shot,module,ca=ca,info=info
	IF NOT keyword_set(info) THEN info=hirexsr_load_info(module,ca=ca)
	xi_ch=info.det.xi[1]-info.det.xi[0]
	xi_o=-1.0*info.det.xi[0]/xi_ch
	zeta_ch=info.det.zeta[info.det.n_xi]-info.det.zeta[0]
	zeta_o=-1.0*info.det.zeta[0]/zeta_ch

	path='\SPECTROSCOPY::TOP.HIREXSR.INFO.MOD'+num2str(module,1)

	mdsopen,'spectroscopy',shot
	;load mirror parameters
	mdsput,path+'.MIRROR:RAD','build_with_units($,"m")',info.m.rad
	mdsput,path+'.MIRROR:ROT','build_with_units($,"radians")',info.m.rot
	mdsput,path+'.MIRROR:SIZE','build_with_units($,"m")',info.m.size
	mdsput,path+'.MIRROR:VEC','build_with_units($,"m")',info.m.vec
	mdsput,path+'.MIRROR.BRAGG:IREF','build_with_units($,"")',info.m.bragg.iref
	mdsput,path+'.MIRROR.BRAGG:RWID','build_with_units($,"mRad")',info.m.bragg.rwid
	mdsput,path+'.MIRROR.BRAGG:TWOD','build_with_units($,"Angstroms")',info.m.bragg.twod

	;load detector parameters
	mdsput,path+'.DET:DET_XI','build_with_units($,"m")',info.det.size[0]
	mdsput,path+'.DET:DET_ZETA','build_with_units($,"m")',info.det.size[1]
	mdsput,path+'.DET:N_XI','build_with_units($,"")',info.det.n_xi
	mdsput,path+'.DET:N_ZETA','build_with_units($,"")',info.det.n_zeta
	mdsput,path+'.DET:X0','build_with_units($,"m")',info.det.x0
	mdsput,path+'.DET:X1','build_with_units($,"m")',info.det.x1
	mdsput,path+'.DET:X2','build_with_units($,"m")',info.det.x2
	mdsput,path+'.DET:XI_CH','build_with_units($,"m")',xi_ch
	mdsput,path+'.DET:XI_O','build_with_units($,"pixel")',xi_o
	mdsput,path+'.DET:ZETA_CH','build_with_units($,"m")',zeta_ch
	mdsput,path+'.DET:ZETA_O','build_with_units($,"pixel")',zeta_o
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_WRITE_DETALIGN
;
;-

PRO hirexsr_write_detalign,shot,module,align
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+'.ALIGN:DETALIGN'
	mdsput,path,'build_with_units($,"")',align
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_DETALIGN
;
;-

PRO hirexsr_load_detalign,shot,module,align,status=status
	mdsopen,'spectroscopy',shot
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+num2str(module,1)+'.ALIGN:DETALIGN'
	align=mdsvalue(path,/quiet,status=status)
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_BINNING
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees
;-

PRO hirexsr_load_binning,shot,chmap,tch,tmap,good,chmax,h=h,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	IF NOT keyword_set(h) THEN path='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.BINNING:' ELSE path='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.BINNING:'
	mdsopen,'spectroscopy',shot
	chmap=mdsvalue(path+'CHMAP')
	tch=mdsvalue(path+'TCH')
	tmap=mdsvalue(path+'TMAP')
	good=mdsvalue(path+'GOOD')
	chmax=mdsvalue(path+'CHMAX')
	mdsclose,'spectroscopy',shot

END

;+
;NAME:
;	HIREXSR_WRITE_BINNING
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees
;	
;-

PRO hirexsr_write_binning,shot,chmap=chmap,tch=tch,tmap=tmap,good=good,chmax=chmax,h=h,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	IF NOT keyword_set(h) THEN path='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.BINNING:' ELSE path='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.BINNING:'	
	mdsopen,'spectroscopy',shot
	IF keyword_set(chmap) THEN mdsput,path+'CHMAP','build_with_units($,"")',chmap
	IF keyword_set(tch) THEN BEGIN
		IF tch[0] EQ -1 THEN tch=[0]
		mdsput,path+'TCH','build_with_units($,"")',float(tch)
	ENDIF
	IF keyword_set(tmap) THEN mdsput,path+'TMAP','build_with_units($,"")',tmap
	IF keyword_set(good) THEN mdsput,path+'GOOD','build_with_units($,"")',int(good)
	IF keyword_set(chmax) THEN mdsput,path+'CHMAX','build_with_units($,"")',int(chmax)
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_COEFS_PTR2ARR
;	
;-

FUNCTION hirexsr_coefs_ptr2arr,coefs,chmax,ntime
	cntr=0
	length=1
	WHILE length EQ 1 AND cntr LE chmax*ntime DO BEGIN
		icoefs=*coefs[cntr]
		length=n(icoefs)+1
		cntr+=1
	ENDWHILE
	arr=fltarr(chmax,ntime,length)
	FOR i=0,chmax-1 DO BEGIN
		FOR j=0,ntime-1 DO BEGIN
			icoefs=*coefs[i,j]
			x=size(icoefs)
			IF x[1] EQ 1 THEN arr[i,j,*]=fltarr(length)+icoefs[0] ELSE arr[i,j,*]=icoefs
		ENDFOR
	ENDFOR

	RETURN,arr	 
END

;+
;NAME:
;	HIREXSR_COEFS_ARR2PTR
;	
;-

FUNCTION hirexsr_coefs_arr2ptr,coefs
	x=size(coefs)
	chmax=x[1]
	ntime=x[2]
	ptr=ptrarr(chmax,ntime,/allocate)
	FOR i=0,chmax-1 DO BEGIN
		FOR j=0,ntime-1 DO BEGIN
			icoefs=reform(coefs[i,j,*])
			x=size(icoefs)
			IF x[1] EQ 1 THEN *ptr[i,j]=icoefs[0] ELSE *ptr[i,j]=icoefs
		ENDFOR
	ENDFOR

	RETURN,ptr
END

;+
;NAME:
;	HIREXSR_LOAD_FITS
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;-

PRO hirexsr_load_fits,shot,line,coefs,ave,double,labels,ptr=ptr,ch=ch,tau=tau,quiet=quiet,status=status,tht=tht
	IF keyword_set(tht) THEN hepath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS'+num2str(tht,1)+'.HELIKE.FITS.' ELSE hepath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.'
	IF keyword_set(tht) THEN hpath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS'+num2str(tht,1)+'.HLIKE.FITS.' ELSE hpath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.FITS.'
	CASE line OF
		0 : path=hepath+'WN3:'	
		1 : path=hepath+'XY:'
		2 : path=hepath+'ZJK:'
		3 : path=hpath+'LYA:'
		4 : path=hpath+'LYA:'
		5 : path=hpath+'TJ:'
		6 : path=hpath+'LYA:'
		7 : path=hepath+'ZJK:'
		8 : path=hepath+'XY:'
		9 : path=hpath+'LYA:'
	ENDCASE
	mdsopen,'spectroscopy',shot
	ave=mdsvalue(path+'AVE',quiet=quiet)
	coefs=mdsvalue(path+'COEFS',quiet=quiet,status=status)
	ch=mdsvalue('dim_of('+path+'COEFS,0)',quiet=quiet)
	tau=mdsvalue('dim_of('+path+'COEFS,1)',quiet=quiet)
	double=mdsvalue(path+'DOUBLE',quiet=quiet)
	labels=mdsvalue(path+'LABELS',quiet=quiet)
	mdsclose,'spectroscopy',shot
	
	IF keyword_set(ptr) THEN coefs=hirexsr_coefs_arr2ptr(coefs)
END

;+
;NAME:
;	HIREXSR_WRITE_FITS
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees	
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;
;-

PRO hirexsr_write_fits,shot,line,coefs,tau,ave,double,labels,tht=tht
	IF keyword_set(tht) THEN hepath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS'+num2str(tht,1)+'.HELIKE.FITS.' ELSE hepath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.'
	IF keyword_set(tht) THEN hpath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS'+num2str(tht,1)+'.HLIKE.FITS.' ELSE hpath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.FITS.'
	CASE line OF
		0 : path=hepath+'WN3:'	
		1 : path=hepath+'XY:'
		2 : path=hepath+'ZJK:'
		3 : path=hpath+'LYA:'
		4 : path=hpath+'LYA:'
		5 : path=hpath+'TJ:'
		6 : path=hpath+'LYA:'		;write wn3 for he-like ca to LYA1
		7 : path=hepath+'ZJK:'		;write high-Te lya1 in ZJK
		8 : path=hepath+'XY:'		;write high-Te Mo4d in XY
		9 : path=hpath+'LYA:'		;write high-Te He-like Ca in LYA1
	ENDCASE
	x=size(coefs)
	IF NOT keyword_set(rhotang) THEN ch=indgen(x[1])+1 ELSE ch=rhotang

	mdsopen,'spectroscopy',shot
	mdsput,path+'AVE', 'build_with_units($,"")',ave
	mdsput,path+'COEFS','build_signal(build_with_units($1,""),*,build_with_units($2,"CH"),build_with_units($3,"seconds"))',coefs,ch,tau
	mdsput,path+'DOUBLE', 'build_with_units($,"")',double
	mdsput,path+'LABELS', 'build_with_units($,"")',labels
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_MLINTPTR
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees	
;	9/18/11		M.L. Reinke - added a check on presence of valid EFIT before filling pointer
;	11/12		M.L. Reinke - added the poloidal velocity [*,6,7]
;	12/6/12		M.L. Reinke - added the brightness [*,8,9]
;-

PRO hirexsr_load_mlintptr,shot,line,mlint,tau,tht=tht,time=time,good=good,clear=clear
	IF line LT 3 THEN modcheck=2 ELSE modcheck=4

	hirexsr_load_image,shot,modcheck,image,rawtime,/noimage
	ntau=n(rawtime)+1
	mlint=ptrarr(ntau,/allocate)
	c=2.998e8 			;speed of light
	e=1.602e-19			;conversion for eV -> J
	mconv=1.661e-27			;conversion for amu -> kg
	
	mdsopen,'analysis',shot
	rmid=mdsvalue('\efit_rmid')
	etime=mdsvalue('dim_of(\efit_rmid,0)')
	epsin=mdsvalue('dim_of(\efit_rmid,1)')
	mdsclose,'analysis',shot

	hirexsr_load_momentptr,shot,line,mom,tau,pos,tpos,lam_o,z,tht=tht,good=good,clear=clear
	mass=read_atomic_mass(z)
	conv_factor=(lam_o/c)^2*(e*1.0e3/(mass*mconv))		;conversion factor for 
	FOR j=0,n(tau)-1 DO BEGIN
		imom=*mom[j]
		eindex=ipt(etime,tau[j])
		IF imom[0] NE -1 AND eindex[0] NE -1 THEN BEGIN
			pindex=last(where(tpos LE tau[j]))
			jpos=pos[*,*,pindex]
			br=imom[*,0]
			brerr=imom[*,3]
			v=-1.0*(lam_o-imom[*,11])*c/lam_o/(2.0*!pi*jpos[2,*]*cos(jpos[3,*]))*1.0e-3
			verr=imom[*,14]*c/lam_o/(2.0*!pi*jpos[2,*]*cos(jpos[3,*]))*1.0e-3
			ti=(imom[*,12])^2/conv_factor
			tierr=2.0*imom[*,15]*sqrt(imom[*,12]^2)/conv_factor
			psintang=imom[*,16]
			rmidtang=interpol(reform(rmid[eindex,*]),epsin,psintang)
			rhotang=(rmidtang-rmid[eindex,0])/(last(rmid[eindex,*])-rmid[eindex,0])

			;compute estimate of poloidal velocity
			vobs=-1.0*(lam_o-imom[*,11])*c/lam_o*1.0e-3		;km/s of observed Doppler shift
			voerr=imom[*,14]*c/lam_o*1.0e-3
			k=minloc(psintang)
			vlow=vobs[0:k]
			vlerr=voerr[0:k]
			vhigh=vobs[k:*]
			vherr=voerr[k:*]
			rlow=psintang[0:k]
			rhigh=psintang[k:*]
			vp=v*0.0
			vperr=v*0.0
			vp[0:k]=(vlow-interpol(vhigh,rhigh,rlow))/2.0		;vp = 1/2 the differnce between below-above
			vperr[0:k]=sqrt(vlerr^2+interpol(vherr,rhigh,rlow)^2)/2.0

			ilint=[[v],[verr],[ti],[tierr],[psintang],[rmidtang],[vp],[vperr],[br],[brerr],[rhotang]]
			*mlint[j]=ilint
		ENDIF ELSE *mlint[j]=-1
	ENDFOR			
END

;+
;NAME:
;	HIREXSR_LOAD_TLINTPTR
;
;PURPOSE:
;	This procedure loads the MLINTPTR data, but organizes it in
;	time rather than space
;
;MODIFICATION HISOTRY:
;	Written by:	M.L. Reinke - adataped from W_HIREXSR_MOMENTS (1/31/13)
;
;-

PRO hirexsr_load_tlintptr,shot,line,tlint,tau,tht=tht
	hirexsr_load_momentptr,shot,line,mom,tau,pos,tpos,lam_o,z,status=status,tht=tht,dlam=dlam,/clear
	hirexsr_load_mlintptr,shot,line,lint,tau,tht=tht,/clear

	;make time-evolving lint pointer
	ntau=n(tau)+1
	ilint=*lint[0]
	nch=n(ilint[*,0])+1		;assume a non time-evolving CHMAP

	tlint=ptrarr(nch,/allocate)
	FOR i=0,nch-1 DO BEGIN
		IF i EQ 0 THEN BEGIN
			ilint=*lint[0]
			nlint=n(ilint[0,*])+1
		ENDIF
		itlint=fltarr(ntau,nlint+1)
		FOR j=0,ntau-1 DO BEGIN
			IF tau[j] NE -1 THEN BEGIN
				ilint=*lint[j]
				imom=*mom[j]
				itlint[j,0:nlint-1]=ilint[i,*]
				itlint[j,nlint]=imom[i,24]	;add the fitcase
                        ENDIF
                ENDFOR
		*tlint[i]=itlint
 	ENDFOR
	heap_free,mom
	heap_free,lint
END

;+
;NAME:
;	HIREXSR_LOAD_LINE_POS
;
;MODIFICATION HISTORY:
;	4/17/12		M.L. Reinke - adopted from HIREXSR_LOAD_LINTPTR
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;
;-

PRO hirexsr_load_line_pos,shot,line,pos,tpos=tpos,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	CASE line OF
		0 : pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.W:POS'
		1 : pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.X:POS'
		2 : pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.Z:POS'
		3 : pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.LYA1:POS'		
		4 : pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.MO4D:POS'
		5 : pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.J:POS'
		6 : pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.LYA1:POS'
		7 : pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.Z:POS'
		8 : pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.X:POS'
		9 : pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.LYA1:POS'
	ENDCASE
	mdsopen,'spectroscopy',shot
	pos=mdsvalue(pospath)
	tpos=mdsvalue('dim_of('+pospath+')')
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_LINTPTR
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees	
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;	4/7/14		M.L. Reinke - hardcoded load_wavelengths to -1
;
;-

PRO hirexsr_load_lintptr,shot,line,lint,tau,tht=tht
	hirexsr_load_wavelengths,-1,lam,z_o,label
	hirexsr_load_fits,shot,line,coefs,ave,double,labels,tau=tau,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	CASE line OF
		0 : BEGIN
			z=18
			lam_o=last(lam[where(z_o EQ z AND label EQ 'w')])
			index=where(labels EQ 'w')
			rhopath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.W:MOM'
			pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.W:POS'
		END
		1 : BEGIN
			z=18
			lam_o=last(lam[where(z_o EQ z AND label EQ 'x')])
			index=where(labels EQ 'x')
			rhopath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.X:MOM'
			pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.X:POS'
		END
		2 : BEGIN
			z=18
			lam_o=last(lam[where(z_o EQ z AND label EQ 'z')])
			index=where(labels EQ 'z')
			rhopath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.Z:MOM'
			pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.Z:POS'
		END
		3 : BEGIN
			z=18
			lam_o=last(lam[where(z_o EQ z AND label EQ 'lya1')])
			index=where(labels EQ 'lya1')
			rhopath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.LYA1:MOM'
			pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.LYA1:POS'
		END		
		4 : BEGIN
			z=42
			lam_o=last(lam[where(z_o EQ z AND label EQ '4d')])
			index=where(labels EQ '4d')
			rhopath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.MO4D:MOM'
			pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.MO4D:POS'
		END
		5 : BEGIN
			z=18
			lam_o=last(lam[where(z_o EQ z AND label EQ 'J')])
			index=where(labels EQ 'J')
			rhopath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.J:MOM'
			pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.J:POS'
		END
		6 : BEGIN
			z=20
			lam_o=last(lam[where(z_o EQ z AND label EQ 'w')])
			index=where(labels EQ 'w')
			rhopath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.LYA1:MOM'
			pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.LYA1:POS'
                    END
		7 : BEGIN									;write high-Te lya1 in Z
			z=18
			lam_o=last(lam[where(z_o EQ z AND label EQ 'lya1')])
			index=where(labels EQ 'lya1')
			rhopath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.Z:MOM'
			pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.Z:POS'
		END		
		8 : BEGIN									;write high-Te 4d in X
			z=42
			lam_o=last(lam[where(z_o EQ z AND label EQ '4d')])
			index=where(labels EQ '4d')
			rhopath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.X:MOM'
			pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.X:POS'
		END
		9 : BEGIN
			z=20
			lam_o=last(lam[where(z_o EQ z AND label EQ 'lya1')])
			index=where(labels EQ 'lya1')
			rhopath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.LYA1:MOM'
			pospath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.LYA1:POS'
                    END
	ENDCASE
	mdsopen,'spectroscopy',shot
	rhotang=mdsvalue('dim_of('+rhopath+',0)')
	pos=mdsvalue(pospath)
	tpos=mdsvalue('dim_of('+pospath+')')
	mdsclose,'spectroscopy',shot
	index=index[0]
	c=2.998e8 			;speed of light
	e=1.602e-19			;conversion for eV -> J
	mconv=1.661e-27			;conversion for amu -> kg
	mass=read_atomic_mass(z)
	conv_factor=(lam_o/c)^2*(e*1.0e3/(mass*mconv))		;conversion factor for 
	
	ntime=n(tau)+1

	lint=ptrarr(ntime,/allocate)
	FOR i=0,ntime-1 DO BEGIN
		IF tau[i] NE -1 THEN BEGIN
			pindex=last(where(tpos LE tau[i]))
			ipos=pos[*,*,pindex]
			icoefs=reform(coefs[*,i,*])
			tmp=where(icoefs[*,0] NE -1)	;number of channels
			icoefs=icoefs[tmp,*]
			irho=rhotang[tmp,i]
			ipos=ipos[*,tmp]		
			nch=n(tmp)+1
			v=fltarr(nch)
			ti=fltarr(nch)
			tmp=where(icoefs[*,0] GE 0)	;indicates a valid fit
			v[tmp]=-1.0*(lam_o-icoefs[tmp,index*3+1])*c/lam_o*1.0e-3/(2.0*!pi*ipos[2,tmp]*cos(ipos[3,tmp]))		;velocity in kHz
			ti[tmp]=icoefs[tmp,index*3+2]^2*mass*mconv*c^2/(e*1.0e3*lam_o^2)					;iontemp in [keV]
			*lint[i]=[[v],[ti],[irho]]
		ENDIF ELSE *lint[i]=-1
	ENDFOR
	;heap_free,rhotang
END
;+
;NAME:
;	HIREXSR_AVESPEC_PTR2ARR
;	
;-

FUNCTION hirexsr_avespec_ptr2arr,avespec,chmax,ntime,maxave	
	arr=fltarr(chmax,ntime,maxave,4)-1
	FOR i=0,chmax-1 DO BEGIN
		FOR j=0,ntime-1 DO BEGIN
			iave=*avespec[i,j]
			x=size(iave)
			IF x[1] NE 1 THEN arr[i,j,0:x[1]-1,*]=iave
		ENDFOR
	ENDFOR

	RETURN,arr
END

;+
;NAME:
;	HIREXSR_AVESPEC_ARR2PTR
;	
;-

FUNCTION hirexsr_avespec_arr2ptr,specbr,lam,sig,resid,maxave=maxave
	x=size(specbr)
	chmax=x[1]
	ntime=x[2]
	maxave=x[3]
	ptr=ptrarr(chmax,ntime,/allocate)
	FOR i=0,chmax-1 DO BEGIN
		FOR j=0,ntime-1 DO BEGIN
			tmp=where(specbr[i,j,*] NE -1)
			IF tmp[0] NE -1 THEN BEGIN
				iptr=fltarr(n(tmp)+1,4)
				iptr[*,0]=specbr[i,j,tmp]
				iptr[*,1]=lam[i,j,tmp]
				iptr[*,2]=sig[i,j,tmp]
				iptr[*,3]=resid[i,j,tmp]
				*ptr[i,j]=iptr
			ENDIF ELSE *ptr[i,j]=-1
		ENDFOR
	ENDFOR

	RETURN,ptr
END

;+
;NAME:
;	HIREXSR_WRITE_AVESPEC
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees	
;	
;-

PRO hirexsr_write_avespec,shot,specbr,lam,sig,resid,tau,nave,h=h,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	IF NOT keyword_set(h) THEN path='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.SPEC:' ELSE path='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.SPEC:'
	IF NOT keyword_set(h) THEN npath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.FITS.WN3:AVE' ELSE npath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.FITS.LYA:AVE'

	x=size(specbr)
	ch=indgen(x[1])+1
	mdsopen,'spectroscopy',shot
	mdsput,path+'LAM','build_signal(build_with_units($1,"Ang"),*,build_with_units($2,"CH"),build_with_units($3,"seconds"))',lam,ch,tau
	mdsput,path+'RESID','build_signal(build_with_units($1,"ph/s/m^2/ang"),*,build_with_units($2,"CH"),build_with_units($3,"seconds"))',resid,ch,tau
	mdsput,path+'SIG','build_signal(build_with_units($1,"ph/s/m^2/ang"),*,build_with_units($2,"CH"),build_with_units($3,"seconds"))',sig,ch,tau
	mdsput,path+'SPECBR','build_signal(build_with_units($1,"ph/s/m^2/ang"),*,build_with_units($2,"CH"),build_with_units($3,"seconds"))',specbr,ch,tau
	mdsput,npath,'build_with_units($,"")',nave
	mdsclose,'spectroscopy',shot
	
END

;+
;NAME:
;	HIREXSR_LOAD_AVESPEC
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees	
;	
;-

PRO hirexsr_load_avespec,shot,specbr,lam,sig,resid,tau,nave,ptr=ptr,h=h,tht=tht,status=status
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	IF NOT keyword_set(h) THEN path='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.SPEC:' ELSE path='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.SPEC:'
	IF NOT keyword_set(h) THEN npath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.FITS.WN3:AVE' ELSE npath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.FITS.LYA:AVE'

	mdsopen,'spectroscopy',shot
	lam=mdsvalue(path+'LAM')
	resid=mdsvalue(path+'RESID')
	sig=mdsvalue(path+'SIG')
	specbr=mdsvalue(path+'SPECBR',status=status)
	tau=mdsvalue('dim_of('+path+'SPECBR,1)')
	nave=mdsvalue(npath)
	mdsclose,'spectroscopy',shot

	IF keyword_set(ptr) THEN specbr=hirexsr_avespec_arr2ptr(specbr,lam,sig,resid)
END



;+
;NAME:
;	HIREXSR_FILL_RESID
;
;PURPOSE:
;	This is a higher level function which fills in the residual in the spec PTRARR once all the fits are completed.
;
;-	

PRO hirexsr_fill_resid,shot
	hirexsr_load_avespec,shot,specbr,lam,sig,resid,tau,/ptr,h=h,tht=tht
	IF keyword_set(h) THEN BEGIN
	ENDIF ELSE BEGIN

	ENDELSE
	
END

;+
;NAME:
;	HIREXSR_LOAD_INST
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 9/4/11
;	4/24/12		M.L. Reinke - added THT optional input
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;
;-

PRO hirexsr_load_inst,shot,line,iwidth,ishift,tpos,status=status,quiet=quiet,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	hepath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS'
	CASE line OF 
		0 : path=hepath+'.W:'
		1 : path=hepath+'.X:'
		2 : path=hepath+'.Z:'
		3 : path=hpath+'.LYA1:'
		4 : path=hpath+'.MO4D:'
		5 : path=hpath+'.J:'
		6 : path=hpath+'.LYA1:'		;write wn3 for he-like ca to LYA1
		7 : path=hepath+'.Z:'		;write high-Te lya1 in Z
		8 : path=hepath+'.X:'		;write high-Te 4d in X
		9 : path=hpath+'.LYA1:'		;write lya1 for h-like ca to LYA1
	ENDCASE
	mdsopen,'spectroscopy',shot
	iwidth=mdsvalue(path+'INST',quiet=quiet,status=status)
	ishift=mdsvalue('dim_of('+path+'INST,0)',quiet=quiet)
	tpos=mdsvalue('dim_of('+path+'INST,1)',quiet=quiet)
	mdsclose,'spectroscopy',shot
END


;+
;NAME:
;	HIREXSR_WRITE_INST
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 9/4/11
;	4/24/12		M.L. Reinke - added THT optional input
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout	
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;
;-

PRO hirexsr_write_inst,shot,line,iwidth,ishift,tpos,status=status,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	hepath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS'
	CASE line OF 
		0 : path=hepath+'.W:'
		1 : path=hepath+'.X:'
		2 : path=hepath+'.Z:'
		3 : path=hpath+'.LYA1:'
		4 : path=hpath+'.MO4D:'
		5 : path=hpath+'.J:'
		6 : path=hpath+'.LYA1:'		;write wn3 for he-like ca to LYA1
		7 : path=hepath+'.Z:'		;write high-Te lya1 in Z
		8 : path=hepath+'.X:'		;write high-Te 4d in X
		9 : path=hpath+'.LYA1:'		;write lya1 for h-like ca to LYA1
	ENDCASE
	mdsopen,'spectroscopy',shot
	chk=mdsvalue('getnci('+path+'INST,"usage")',status=status,/quiet)
	IF status THEN mdsput,path+'INST','build_signal(build_with_units($1,"Ang"),*,build_with_units($2,"Ang"),build_with_units($3,"sec"))',iwidth,ishift,tpos
	mdsclose,'spectroscopy',shot
END


;+
;NAME:
;	HIREXSR_LOAD_MOMENTS
;
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees	
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;
;-

PRO hirexsr_load_moments,shot,line,mom,err,pmom,perr,tau,pos,rhotang,bfrac,scale,u,tpos,dlam,double,fitcase,tree,quiet=quiet,status=status,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	hepath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS'
	CASE line OF 
		0 : path=hepath+'.W:'
		1 : path=hepath+'.X:'
		2 : path=hepath+'.Z:'
		3 : path=hpath+'.LYA1:'
		4 : path=hpath+'.MO4D:'
		5 : path=hpath+'.J:'
		6 : path=hpath+'.LYA1:'	
		7 : path=hepath+'.Z:'		;write high-Te lya1 in Z
		8 : path=hepath+'.X:'		;write high-Te 4d in X
		9 : path=hpath+'.LYA1:'		;write lya1 for h-like ca to LYA1
	ENDCASE
	mdsopen,'spectroscopy',shot
	mom=mdsvalue(path+'MOM',quiet=quiet,status=status)
	IF NOT status THEN RETURN
	rhotang=mdsvalue('dim_of('+path+'MOM,0)')
	pmom=mdsvalue('dim_of('+path+'MOM,2)')
	bfrac=mdsvalue('dim_of('+path+'MOM,3)')
	fitcase=mdsvalue('dim_of('+path+'MOM,4)',/quiet,status=fitstat)
	IF NOT fitstat THEN BEGIN
		x=size(bfrac)
		fitcase=fltarr(x[1],x[2])
		tmp=where(bfrac EQ 0)
		fitcase[tmp]=-1
	ENDIF
	err=mdsvalue(path+'ERR',quiet=quiet)
	perr=mdsvalue('dim_of('+path+'ERR,2)')
	scale=mdsvalue('dim_of('+path+'ERR,3)')
	tau=mdsvalue('dim_of('+path+'MOM,1)',quiet=quiet)
	pos=mdsvalue(path+'POS',quiet=quiet)
	u=mdsvalue(path+'U',quiet=quiet)
	tpos=mdsvalue('dim_of('+path+'POS,0)',quiet=quiet)
	tree=mdsvalue('dim_of('+path+'POS,1)',quiet=quiet,status=treestatus)
	IF NOT treestatus THEN tree='analysis'
	dlam=mdsvalue(path+'DLAM',quiet=quiet)
	double=mdsvalue(path+'DOUBLE',quiet=quiet)
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_WRITE_MOMENTS
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees	
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;
;-

PRO hirexsr_write_moments,shot,line,mom,err,pmom,perr,tau,pos,rhotang,bfrac,scale,u,tpos,dlam,double,fitcase,tree,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	hepath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS'
	CASE line OF 
		0 : path=hepath+'.W:'
		1 : path=hepath+'.X:'
		2 : path=hepath+'.Z:'
		3 : path=hpath+'.LYA1:'
		4 : path=hpath+'.MO4D:'
		5 : path=hpath+'.J:'
		6 : path=hpath+'.LYA1:'		;write wn3 for he-like ca to LYA1
		7 : path=hepath+'.Z:'		;write high-Te lya1 in Z
		8 : path=hepath+'.X:'		;write high-Te 4d in X
		9 : path=hpath+'.LYA1:'		;write lya1 for h-like ca to LYA1
	ENDCASE
	x=size(mom)

	mdsopen,'spectroscopy',shot
	mdsput,path+'MOM','build_signal(build_with_units($1,""),*,build_with_units($2,"PSIN"),build_with_units($3,"seconds"),build_with_units($4,""),build_with_units($5,""),build_with_units($6,""))',$
		mom,rhotang,tau,pmom,bfrac,fitcase
	mdsput,path+'ERR','build_signal(build_with_units($1,""),*,build_with_units($2,"PSIN"),build_with_units($3,"seconds"),build_with_units($4,""),build_with_units($5,""))',err,rhotang,tau,perr,scale
	mdsput,path+'POS','build_signal(build_with_units($1,"POS"),*,build_with_units($2,"seconds"),build_with_units($3," "))',pos,tpos,tree
	mdsput,path+'U','build_signal(build_with_units($1,"m^2-str"),*,build_with_units($2,"seconds"))',u,tpos
	mdsput,path+'DOUBLE', 'build_with_units($,"")',double
	mdsput,path+'DLAM', 'build_with_units($,"mAng")',dlam
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIXERSR_LOAD_MOMENTPTR,shot,line
;
;OUTPUTS:
;	momptr		PTRARR [ntime] which has the inputs to the inversion code in format that's used in HIREXSER_CALC_PROFILES
;			arr=[[reform(mom[tmp,i,*])],$		;0:2 are the 0th, 1st and 2nd moments
;			     [reform(err[tmp,i,*])],$		;3:5 are the 0th, 1st and 2nd moment errors
;		             [igood],$				;6 is the good vector
;			     [icheck],$				;7:9 is the moment check
;			     [reform(pmom[tmp,i,*])],$		;10-12 is the profile moment
;			     [reform(perr[tmp,i,*])],$		;13-15 is profile moment error
;			     [rhotang[tmp,i]],$			;16 is the rhotang
;			     [bfrac[tmp,i]],$			;17 is the background fraction
;			     [scale[tmp,i]],$			;18 is the MOM0=scale*NPHOT constant
;			     [iw],$				;19 is the instrumental width [Ang]
;			     [is],$				;20 is the instrumental shift [Ang]
;			     [isubcheck]]			;21-23 are the sub vectors (inst_v_vec, inst_ti_vec, sub_vec)	
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees	
;	8/25/11		M.L. Reinke - fixed a rather important bug in loading BRFAC & SCALE.  It only used the 0th time slice data for all points :(
;	9/4/11		M.L. Reinke - added even more crap to the momptr, now has width/shift instrumentals [19,20] and subvec checks [21-23]
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout
;	11/20/12	M.L. Reinke - removed keyword NEW and added GOOD optional input
;	11/27/12	M.L. Reinke - added the keyowrd CLEAR to toggle CHECK/SUBCHECK loading
;	12/6/12		M.L. Reinke - added the dlam optional output
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca	
;	4/7/14		M.L. Reinke - hard coded load_wavelengths to use -1
;
;-

PRO hirexsr_load_momentptr,shot,line,momptr,tau,pos,tpos,lam_o,z,tree=tree,status=status,good=good,tht=tht,clear=clear,dlam=dlam
	hirexsr_load_wavelengths,-1,lam,z_o,label
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	CASE line OF
		0 : BEGIN
			z=18
			lam_o=last(lam[where(z_o EQ z AND label EQ 'w')])
			gpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.BINNING:GOOD'
			fpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.FITS.WN3:'
		END
		1 : BEGIN
			z=18
			lam_o=last(lam[where(z_o EQ z AND label EQ 'x')])
			gpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.BINNING:GOOD'
			fpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.FITS.XY:'
		END
		2 : BEGIN
			z=18
			lam_o=last(lam[where(z_o EQ z AND label EQ 'z')])
			gpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.BINNING:GOOD'
			fpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.FITS.ZJK:'
		END
		3 : BEGIN
			z=18
			lam_o=last(lam[where(z_o EQ z AND label EQ 'lya1')])
			gpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.BINNING:GOOD'
			fpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.FITS:'
		END		
		4 : BEGIN
			z=42
			lam_o=last(lam[where(z_o EQ z AND label EQ '4d')])
			gpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.BINNING:GOOD'
			fpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.FITS.LYA:'
		END
		5 : BEGIN
			z=18
			lam_o=last(lam[where(z_o EQ z AND label EQ 'J')])
			gpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.BINNING:GOOD'
			fpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.FITS.TJ:'
		END
		6 : BEGIN
			z=20
			lam_o=last(lam[where(z_o EQ z AND label EQ 'w')])
			gpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.BINNING:GOOD'
			fpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.FITS:'
                END
		7 : BEGIN									;write high-Te lya1 in Z
			z=18
			lam_o=last(lam[where(z_o EQ z AND label EQ 'lya1')])
			gpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.BINNING:GOOD'
			fpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.FITS.ZJK:'
		END		
		8 : BEGIN									;write high-Te 4d in X
			z=42
			lam_o=last(lam[where(z_o EQ z AND label EQ '4d')])
			gpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.BINNING:GOOD'
			fpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.FITS.XY:'
		END
		9 : BEGIN		
			z=20
			lam_o=last(lam[where(z_o EQ z AND label EQ 'lya1')])
			gpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.BINNING:GOOD'
			fpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.FITS:'
                END
	ENDCASE
	mdsopen,'spectroscopy',shot
	bingood=mdsvalue(gpath)				;load good from BINNING to initialize the inversion if not present
	mdsclose,'spectroscopy',shot
	hirexsr_load_moments,shot,line,mom,err,pmom,perr,tau,pos,rhotang,bfrac,scale,u,tpos,dlam,double,fitcase,tree,status=status,tht=tht,/quiet
	hirexsr_load_inst,shot,line,iwidth,ishift,tht=tht,status=istatus,/quiet
	IF NOT status THEN RETURN
	hirexsr_load_check,shot,line,check,chkgood,subcheck,gs=gs,cs=cs,tht=tht						
	
	;select appropriate GOOD vector
	xgood=bingood								;use binning as a default if not externally or internally available
      	IF keyword_set(good) THEN BEGIN
		IF good[0] EQ -1 AND gs THEN xgood=chkgood 			;use check good for good=-1 if available
		IF good[0] NE -1 THEN BEGIN						
			ntau=n(tau)+1
			ngood=n(good)+1
			xgood=intarr(ngood,ntau)							
			FOR i=0,ntau-1 DO IF tau[i] NE -1 THEN xgood[*,i]=good		;copy specified good vector for all times
		ENDIF
        ENDIF

	ntime=n(tau)+1
	momptr=ptrarr(ntime,/allocate)
	FOR i=0,ntime-1 DO BEGIN			;at each index write a [nch,24] array to the pointer array
		IF tau[i] NE -1 THEN BEGIN
			tmp=where(mom[*,i,0] NE -1)
			pindex=last(where(tpos LE tau[i]))
			IF istatus THEN BEGIN
				iw=iwidth[tmp,pindex]		 
				is=ishift[tmp,pindex]
			ENDIF ELSE BEGIN
				iw=fltarr(n(tmp)+1)		;if no instrumental width stored, set to 0.0
				is=fltarr(n(tmp)+1)		;if no instrumental shift stored, set to 0.0
			ENDELSE
			IF NOT gs THEN igood=xgood[tmp,i] ELSE igood=xgood[0:n(tmp),i]
			IF NOT cs OR keyword_set(clear) THEN BEGIN
				icheck=fltarr(n(tmp)+1,3)	;set check to zero
				isubcheck=fltarr(n(tmp)+1,3)	;set subcheck to zero
			ENDIF ELSE BEGIN
				icheck=reform(check[0:n(tmp),i,*])
				isubcheck=reform(subcheck[0:n(tmp),i,*])
			ENDELSE
			arr=[[reform(mom[tmp,i,*])],$		;0:2 are the 0th, 1st and 2nd moments
			     [reform(err[tmp,i,*])],$		;3:5 are the 0th, 1st and 2nd moment errors
		             [igood],$				;6 is the good vector
			     [icheck],$				;7:9 is the moment check
			     [reform(pmom[tmp,i,*])],$		;10-12 is the profile moment
			     [reform(perr[tmp,i,*])],$		;13-15 is profile moment error
			     [rhotang[tmp,i]],$			;16 is the rhotang
			     [bfrac[tmp,i]],$			;17 is the background fraction
			     [scale[tmp,i]],$			;18 is the MOM0=scale*NPHOT constant
			     [iw],$				;19 is the instrumental width [Ang]
			     [is],$				;20 is the instrumental shift [Ang]
			     [isubcheck],$			;21-23 are the sub vectors (inst_v_vec, inst_ti_vec, sub_vec)
			     [fitcase[tmp,i]]]			;24 is the fitcase
			*momptr[i]=arr	
		ENDIF ELSE *momptr[i]=-1
	ENDFOR
END

;+
;NAME:
;	HIREXSR_EXTRACT_CHECK
;
;MODIFICATION HISTORY:
;	9/4/11		M.L. Reinke - added ability to extract the subcheck data from the moment ptr
;	
;-

PRO hirexsr_extract_check,mom,check,good,subcheck
	x=size(mom)
	ntime=x[1]
	nch=0
	FOR i=0,ntime-1 DO BEGIN
		x=size(*mom[i])
		IF x[1] GT nch THEN nch=x[1]
	ENDFOR
	
	check=fltarr(nch,ntime,3)
	subcheck=fltarr(nch,ntime,3)
	good=intarr(nch,ntime)
	FOR i=0,ntime-1 DO BEGIN
		imom=*mom[i]
		x=size(imom)
		IF x[0] NE 0 THEN BEGIN
			check[0:x[1]-1,i,*]=imom[*,7:9]
			good[0:x[1]-1,i]=imom[*,6]
		ENDIF
		IF x[2] GT 19 THEN subcheck[0:x[1]-1,i,*]=imom[*,21:23]
	ENDFOR
END


;+
;NAME:
;	HIXEXSR_VOXEL_PTR2ARR
;
;-

FUNCTION hirexsr_voxel_ptr2arr,voxel
	x=size(voxel)
	ntime=x[1]
	nch=0
	nrho=0
	FOR i=0,ntime-1 DO BEGIN
		x=size(*voxel[i])
		IF x[1] GT nch THEN nch=x[1]
		IF x[2] GT nrho THEN nrho=x[2]
	ENDFOR
	x=size(*voxel[0])			;determine if it's voxel or velvoxel based on size
	IF x[0] EQ 3 THEN BEGIN
		voxarr=fltarr(nch,nrho,ntime,2)-1.0
		FOR i=0,ntime-1 DO BEGIN
			ivox=*voxel[i]
			x=size(ivox)
			IF x[0] EQ 3 THEN BEGIN
				voxarr[0:x[1]-1,0:x[2]-1,i,0]=ivox[*,*,0]
				voxarr[0:x[1]-1,0:x[2]-1,i,1]=ivox[*,*,1]
			ENDIF
		ENDFOR
	ENDIF ELSE BEGIN
		voxarr=fltarr(nch,nrho,ntime)-1.0
		FOR i=0,ntime-1 DO BEGIN
			ivox=*voxel[i]
			x=size(ivox)
			IF x[0] EQ 2 THEN voxarr[0:x[1]-1,0:x[2]-1,i]=ivox
		ENDFOR
	ENDELSE
	RETURN,voxarr
END

;+
;NAME:
;	HIXEXSR_VOXEL_ARR2PTR
;
;-

FUNCTION hirexsr_voxel_arr2ptr,voxarr
	x=size(voxarr)
	ntime=x[3]
	voxel=ptrarr(ntime,/allocate)
	FOR i=0,ntime-1 DO BEGIN
		ivoxarr=reform(voxarr[*,*,i,*])
		IF ivoxarr[0] NE -1.0 THEN BEGIN
			tmp=where(ivoxarr[*,0,0] NE -1.0)
			nch=last(tmp)+1
			tmp=where(ivoxarr[0,*,0] NE -1.0)
			nrho=last(tmp)+1
			*voxel[i]=ivoxarr[0:nch-1,0:nrho-1,*]
		ENDIF ELSE *voxel[i]=-1.0
	ENDFOR
	RETURN,voxel
END

;+
;NAME:
;	HIREXSR_RHO_PTR2ARR
;
;-

FUNCTION hirexsr_rho_ptr2arr,rho
	x=size(rho)
	ntime=x[1]
	nrho=0
	FOR i=0,ntime-1 DO BEGIN
		x=size(*rho[i])
		IF x[1] GT nrho THEN nrho=x[1]
	ENDFOR
	rhoarr=fltarr(nrho,ntime)-1.0
	FOR i=0,ntime-1 DO BEGIN
		irho=*rho[i]
		IF irho[0] NE -1 THEN rhoarr[0:n(irho),i]=irho
	ENDFOR
	RETURN,rhoarr
END

;+
;NAME:
;	HIREXSR_RHO_ARR2PTR
;
;-

FUNCTION hirexsr_rho_arr2ptr,rhoarr
	x=size(rhoarr)
	ntime=x[2]
	rho=ptrarr(ntime,/allocate)
	FOR i=0,ntime-1 DO BEGIN
		tmp=where(rhoarr[*,i] NE -1.0)
		IF tmp[0] EQ -1 THEN *rho[i]=-1.0 ELSE *rho[i]=rhoarr[tmp,i]
	ENDFOR
	RETURN,rho
END

;+
;NAME:
;	HIREXSR_WRITE_VOXEL
;
;
;MODIFICATION HISTORY:
;	4/2/11		M.L. Reinke - modified to store the m1s voxel matrix with the radial only voxel matrix
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;
;-

PRO hirexsr_write_voxel,shot,line,voxel,velvoxel,rho,tau,tree,m1svoxel=m1svoxel,m1svelvoxel=m1svelvoxel,tht=tht
	x=size(voxel,/type)
	IF x EQ 10 THEN varr=hirexsr_voxel_ptr2arr(voxel) ELSE varr=voxel
	x=size(velvoxel,/type)
	IF x EQ 10 THEN vvarr=hirexsr_voxel_ptr2arr(velvoxel) ELSE vvarr=velvoxel
	x=size(rho,/type)
	IF x EQ 10 THEN rhoarr=hirexsr_rho_ptr2arr(rho) ELSE rhoarr=rho
	IF keyword_set(m1svoxel) THEN BEGIN
		x=size(m1svoxel,/type)
		IF x EQ 10 THEN v1arr=hirexsr_voxel_ptr2arr(m1svoxel) ELSE v1arr=m1svoxel
		varr=[[varr],[v1arr]]
	ENDIF
	IF keyword_set(m1svelvoxel) THEN BEGIN
		x=size(m1svelvoxel,/type)
		IF x EQ 10 THEN vv1arr=hirexsr_voxel_ptr2arr(m1svelvoxel) ELSE vv1arr=m1svelvoxel
		vvarr=[[vvarr],[vv1arr]]
	ENDIF
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	hepath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.PROFILES'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.PROFILES'
	CASE line OF 
		0 : path=hepath+'.W'
		1 : path=hepath+'.X'
		2 : path=hepath+'.Z'
		3 : path=hpath+'.LYA1'
		4 : path=hpath+'.MO4D'
		5 : path=hpath+'.J'
		6 : path=hpath+'.LYA1' 		;write wn3 for he-like ca to LYA1
		7 : path=hepath+'.Z'		;write high-Te lya1 in Z
		8 : path=hepath+'.X'		;write high-Te 4d in X
   		9 : path=hpath+'.LYA1' 		;write lya1 for h-like ca to LYA1
     	ENDCASE	
	mdsopen,'spectroscopy',shot
	mdsput,path+'.CONFIG:VOXEL', 'build_with_units($,"m^3")',varr
	mdsput,path+'.CONFIG:VELVOXEL', 'build_with_units($,"m^3")',vvarr
	mdsput,path+':RHO','build_signal(build_with_units($1,"psin"),*,build_with_units($2,"seconds"),build_with_units($3," "))',rhoarr,tau,tree
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_VOXEL
;
;MODIFICATION HISTORY:
;	4/2/11		M.L. Reinke - modified to check and if present, load the m1s voxel matrix stored w/ the radial only voxel matrix
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees
;	8/22/11		M.L. Reinke - fixed a bug that when NOT use /ptr the m1svelvoxel wouldn't get filled
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;
;-

PRO hirexsr_load_voxel,shot,line,voxel,velvoxel,rho,tau,tree,ptr=ptr,m1svoxel=m1svoxel,m1svelvoxel=m1svelvoxel,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	hepath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.PROFILES'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.PROFILES'
	CASE line OF 
		0 : path=hepath+'.W'
		1 : path=hepath+'.X'
		2 : path=hepath+'.Z'
		3 : path=hpath+'.LYA1'
		4 : path=hpath+'.MO4D'
		5 : path=hpath+'.J'
		6 : path=hpath+'.LYA1'		;write wn3 for he-like ca to LYA1
		7 : path=hepath+'.Z'		;write high-Te lya1 in Z
		8 : path=hepath+'.X'		;write high-Te 4d in X
		9 : path=hpath+'.LYA1'		;write lya1 for h-like ca to LYA1
	ENDCASE
	mdsopen,'spectroscopy',shot
	varr=mdsvalue(path+'.CONFIG:VOXEL')
	vvarr=mdsvalue(path+'.CONFIG:VELVOXEL')
	rhoarr=mdsvalue(path+':RHO')
	tau=mdsvalue('dim_of('+path+':RHO,0)')
	tree=mdsvalue('dim_of('+path+':RHO,1)',/quiet,status=treestatus)
	IF NOT treestatus THEN tree='analysis'
	mdsclose,'spectroscopy',shot

	;check to see if voxel arrays have m1s with it
	x=size(rhoarr) 
	y=size(varr)
	z=size(vvarr)
	IF y[2] EQ 2*x[1] THEN BEGIN		;indicates that a m1svox has been stored
		v1arr=varr[*,x[1]:*,*]
		varr=varr[*,0:x[1]-1,*]
	ENDIF ELSE v1arr=-1
	IF z[2] EQ 2*x[1] THEN BEGIN		;indicates that a m1svox has been stored
		vv1arr=vvarr[*,x[1]:*,*,*]
		vvarr=vvarr[*,0:x[1]-1,*,*]
	ENDIF ELSE vv1arr=-1
	
	IF keyword_set(ptr) THEN BEGIN
		voxel=hirexsr_voxel_arr2ptr(varr)
		velvoxel=hirexsr_voxel_arr2ptr(vvarr)
		rho=hirexsr_rho_arr2ptr(rhoarr)
		IF v1arr[0] NE -1 THEN m1svoxel=hirexsr_voxel_arr2ptr(v1arr)
		IF vv1arr[0] NE -1 THEN m1svelvoxel=hirexsr_voxel_arr2ptr(vv1arr)
	ENDIF ELSE BEGIN
		voxel=varr
		velvoxel=vvarr
		rho=rhoarr
		IF v1arr[0] NE -1 THEN m1svoxel=v1arr
		IF vv1arr[0] NE -1 THEN m1svelvoxel=vv1arr
	ENDELSE
END

;+
;NAME:
;	HIREXSR_PROFILE_PTR2ARR
;
;-

PRO hirexsr_profile_ptr2arr,profiles,pro_arr,proerr_arr,rho_arr
	x=size(profiles)
	ntime=x[1]
	nrho=0
	FOR i=0,ntime-1 DO BEGIN
		ipro=*profiles[i]
		x=size(ipro)
		IF x[1] GT nrho THEN nrho=x[1]
	ENDFOR
	pro_arr=fltarr(nrho,ntime,4)-1.0
	proerr_arr=fltarr(nrho,ntime,4)-1.0
	rho_arr=fltarr(nrho,ntime)-1.0
	FOR i=0,ntime-1 DO BEGIN
		ipro=*profiles[i]
		x=size(ipro)
		IF x[0] NE 0 THEN BEGIN
			pro_arr[0:x[1]-1,i,*]=ipro[*,0:3]
			proerr_arr[0:x[1]-1,i,*]=ipro[*,4:7]
			rho_arr[0:x[1]-1,i]=ipro[*,8]
		ENDIF
	ENDFOR
END

;+
;NAME:
;	HIREXSR_PROFILE_ARR2PTR
;
;-

FUNCTION hirexsr_profile_arr2ptr,pro_arr,proerr_arr,rho_arr
	x=size(rho_arr)
	ntime=x[2]
	profiles=ptrarr(ntime,/allocate_heap)
	FOR i=0,ntime-1 DO BEGIN
		tmp=where(rho_arr[*,i] NE -1.0)
		IF tmp[0] NE -1 THEN BEGIN
			arr=[[reform(pro_arr[tmp,i,*])],$
			      [reform(proerr_arr[tmp,i,*])],$
			      [rho_arr[tmp,i]]]
			*profiles[i]=arr
		ENDIF ELSE *profiles[i]=-1.0
	ENDFOR	

	RETURN,profiles
END

;+
;NAME:
;	HIREXSR_WRITE_PROFILE
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees	
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;
;-

PRO hirexsr_write_profile,shot,line,pro_arr,proerr_arr,rho_arr,tau,tinst=tinst,tgood=tgood,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	hepath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.PROFILES'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.PROFILES'
	CASE line OF 
		0 : path=hepath+'.W:'
		1 : path=hepath+'.X:'
		2 : path=hepath+'.Z:'
		3 : path=hpath+'.LYA1:'
		4 : path=hpath+'.MO4D:'
		5 : path=hpath+'.J:'
		6 : path=hpath+'.LYA1:'		;write wn3 for he-like ca to LYA1
		7 : path=hepath+'.Z:'		;write high-Te lya1 in Z
		8 : path=hepath+'.X:'		;write high-Te 4d in X
		9 : path=hpath+'.LYA1:'		;write lya1 for h-like ca to LYA1
	ENDCASE
	IF NOT keyword_set(tinst) THEN tinst=0.0
	IF NOT keyword_set(tgood) THEN BEGIN
		tgood=tau
		tmp=where(tau NE -1)
		tgood[tmp]=1
	ENDIF
	mdsopen,'spectroscopy',shot
	mdsput,path+'PRO','build_signal(build_with_units($1,""),*,build_with_units($2,"psin"),build_with_units($3,"seconds"),build_with_units($4,"keV"),build_with_units($5,""))',$
		pro_arr,rho_arr,tau,tinst,tgood
	mdsput,path+'PROERR','build_signal(build_with_units($1,""),*,build_with_units($2,"psin"),build_with_units($3,"seconds"))',proerr_arr,rho_arr,tau
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_PROFILE
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees	
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;
;-

PRO hirexsr_load_profile,shot,line,pro_arr,proerr_arr,rho_arr,tau,tinst=tinst,tgood=tgood,tht=tht,status=status
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	hepath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.PROFILES'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.PROFILES'
	CASE line OF 
		0 : path=hepath+'.W:'
		1 : path=hepath+'.X:'
		2 : path=hepath+'.Z:'
		3 : path=hpath+'.LYA1:'
		4 : path=hpath+'.MO4D:'
		5 : path=hpath+'.J:'
		6 : path=hpath+'.LYA1:'		;write wn3 for he-like ca to LYA1
		7 : path=hepath+'.Z:'		;write high-Te lya1 in Z
		8 : path=hepath+'.X:'		;write high-Te 4d in X
		9 : path=hpath+'.LYA1:'		;write lya1 for h-like ca to LYA1
	ENDCASE	

	mdsopen,'spectroscopy',shot
	pro_arr=mdsvalue(path+'PRO',/quiet,status=status)
	proerr_arr=mdsvalue(path+'PROERR',/quiet)
	rho_arr=mdsvalue('dim_of('+path+'PRO,0)',/quiet)
	tau=mdsvalue('dim_of('+path+'PRO,1)',/quiet)
	tinst=mdsvalue('dim_of('+path+'PRO,2)',/quiet,status=tstatus)
	IF NOT tstatus THEN tinst=0.0
	tgood=mdsvalue('dim_of('+path+'PRO,3)',/quiet,status=gstatus)
	IF gstatus AND status THEN BEGIN	
		tgood=tau
		tmp=where(tau NE -1)
		tgood[tmp]=1
	ENDIF
	mdsclose,'spectroscopy',shot	
END

;+
;NAME:
;	HIREXSR_BSFIT_PTR2ARR
;
;PURPOSE:
;	Converts the BSFIT pointer array to float arrays for storage in the tree
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 6/24/13
;	7/11/13		M.L. Reinke - modified to 'BSFIT' from 'BSTI' to handle generic structure
;
;-

PRO hirexsr_bsfit_ptr2arr,bsfit,fit,data,config,rlab=rlab
	ndat=0
	nfit=0
	ntime=n(bsfit)+1
	FOR i=0,ntime-1 DO BEGIN
		ibs=*bsfit[i]
		nfit = nfit > ibs.config.nrho
		ndat = ndat > (n(ibs.dat.(0))+1)
        ENDFOR
	
	fit=fltarr(nfit,ntime,5)-1.0		;ti, rho, err, dti, derr
	data=fltarr(ndat,ntime,5)-1.0		;ti,rho,err,good,type
	config=fltarr(ntime,21)			;in order in the config nodes
	FOR i=0,ntime-1 DO BEGIN
		ibs=*bsfit[i]
		fit[0:ibs.config.nrho-1,i,0]=ibs.fit.(0)
		fit[0:ibs.config.nrho-1,i,1]=ibs.fit.rho
		fit[0:ibs.config.nrho-1,i,2]=ibs.fit.err
		fit[0:ibs.config.nrho-1,i,3]=ibs.fit.(3)
		fit[0:ibs.config.nrho-1,i,4]=ibs.fit.derr
		idat=n(ibs.dat.(0))+1
		data[0:idat-1,i,0]=ibs.dat.(0)
		data[0:idat-1,i,1]=ibs.dat.rho
		data[0:idat-1,i,2]=ibs.dat.err
		data[0:idat-1,i,3]=ibs.fit.good
		data[0:idat-1,i,4]=ibs.dat.type
		IF i EQ 0 THEN rlab=ibs.dat.rlab
		config[i,0]=ibs.error.a
		config[i,1]=ibs.error.al
		config[i,2:3]=ibs.limits.alrho
		config[i,4:5]=ibs.limits.arho
		config[i,6]=ibs.error.b
		config[i,7]=ibs.error.bl
		config[i,8:9]=ibs.limits.blrho
		config[i,10:11]=ibs.limits.brho
		config[i,12]=ibs.limits.errmax
		config[i,13]=ibs.config.fitcase
		config[i,14]=ibs.config.nknots
		config[i,15]=ibs.config.nrho	
		config[i,16]=ibs.config.nsigma
		config[i,17]=ibs.config.ntrials
		config[i,18]=ibs.config.order
		config[i,19]=ibs.limits.ymax
		config[i,20]=ibs.limits.ymin
	ENDFOR
END

;+
;NAME:
;	HIREXSR_BSTI_ARR2PTR
;
;PURPOSE:
;	Converts float arrays loaded from the tree using HIREXSR_LOAD_BSTI into the PTRARR format used by other codes
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 6/24/13
;
;-

PRO hirexsr_bsti_arr2ptr,fit,data,config,shot,tht,time,bsti,rlab=rlab
	IF NOT keyword_set(rlab) THEN rlab='r/a'
	ntime=n(time)+1
	bsti=ptrarr(ntime,/allocate_heap)

	FOR i=0,ntime-1 DO BEGIN
		ftmp=where(fit[*,i,1] NE -1)		;identify non-empty fit indices
		dtmp=where(data[*,i,1] NE -1)		;identify non-empty data indices
		ifit={prof:fit[ftmp,i,0],rho:fit[ftmp,i,1],err:fit[ftmp,i,2],dprof:fit[ftmp,i,3],derr:fit[ftmp,i,4],good:data[dtmp,i,3]}
		idat={ti:data[dtmp,i,0],rho:data[dtmp,i,1],err:data[dtmp,i,2],type:data[dtmp,i,4],rlab:rlab}
		ierror={a:config[i,0],b:config[i,6],al:config[i,1],bl:config[i,7]}
		ilimits={arho:reform(config[i,4:5]),brho:reform(config[i,10:11]),alrho:reform(config[i,2:3]),blrho:reform(config[i,8:9]),$
			ymin:config[i,20],ymax:config[i,19s],errmax:config[i,12]}
		iconfig={order:config[i,18],nknots:config[i,14],ntrials:config[i,17],nsigma:config[i,16],nrho:config[i,15],fitcase:config[i,13]}
		ibs={fit:ifit,dat:idat,limits:ilimits,error:ierror,config:iconfig,shot:shot,tht:tht,time:time[i]}
		*bsti[i]=ibs
        ENDFOR
END

;+
;NAME:
;	HIREXSR_BSOM_ARR2PTR
;
;PURPOSE:
;	Converts float arrays loaded from the tree using HIREXSR_LOAD_BSOM into the PTRARR format used by other codes
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 7/11/13 (based on HIREXSR_BSTI_ARR2PTR)
;
;-

PRO hirexsr_bsom_arr2ptr,fit,data,config,shot,tht,time,bsom,rlab=rlab
	IF NOT keyword_set(rlab) THEN rlab='r/a'
	ntime=n(time)+1
	bsom=ptrarr(ntime,/allocate_heap)

	FOR i=0,ntime-1 DO BEGIN
		ftmp=where(fit[*,i,1] NE -1)		;identify non-empty fit indices
		dtmp=where(data[*,i,1] NE -1)		;identify non-empty data indices
		ifit={prof:fit[ftmp,i,0],rho:fit[ftmp,i,1],err:fit[ftmp,i,2],dprof:fit[ftmp,i,3],derr:fit[ftmp,i,4],good:data[dtmp,i,3]}
		idat={om:data[dtmp,i,0],rho:data[dtmp,i,1],err:data[dtmp,i,2],type:data[dtmp,i,4],rlab:rlab}
		ierror={a:config[i,0],b:config[i,6],al:config[i,1],bl:config[i,7]}
		ilimits={arho:reform(config[i,4:5]),brho:reform(config[i,10:11]),alrho:reform(config[i,2:3]),blrho:reform(config[i,8:9]),$
			ymin:config[i,20],ymax:config[i,19s],errmax:config[i,12]}
		iconfig={order:config[i,18],nknots:config[i,14],ntrials:config[i,17],nsigma:config[i,16],nrho:config[i,15],fitcase:config[i,13]}
		ibs={fit:ifit,dat:idat,limits:ilimits,error:ierror,config:iconfig,shot:shot,tht:tht,time:time[i]}
		*bsom[i]=ibs
        ENDFOR
END

;+
;NAME:
;	HIREXSR_WRITE_BSTI
;
;PURPOSE:
;	This procedure writes the bsti array information to the tree
;
;OPTIONAL INPUTS:
;	dc	FLOAT	of the instrumental ti [keV] to set for INST for all channels/times.  Uncertainty is set to 25% of dc value
;
;MODFIICATION HISTORY:
;	Written by:	M.L. Reinke - 6/24/13
;	9/23/14		M.L. Reinke - modified INST writing to allow over writing back to no instrumental
;
;-

PRO hirexsr_write_bsti,shot,fit,data,inst,config,time,rlab=rlab,tht=tht,dc=dc
	IF keyword_set(dc) THEN BEGIN
		inst=[[[fit[*,*,0]*0.0+dc]],[[fit[*,*,0]*0.0+0.25*dc]]]
	ENDIF
	IF NOT keyword_set(rlab) THEN rlab=''
	IF keyword_set(tht) THEN rstr='RESULTS'+num2str(tht) ELSE rstr='RESULTS'
	path='\SPECTROSCOPY::TOP.HIREXSR.'+rstr+'.BSTI:'
	mdsopen,'spectroscopy',shot
	mdsput,path+'FIT','build_signal(build_with_units($,"kev"),*,build_with_units($2,"'+rlab+'"), build_with_units($3,"seconds"), build_with_units($4,"keV"),build_with_units($5,"kev/rho"), build_with_units($6,"keV/rho"))',$
		fit[*,*,0],fit[*,*,1],time,fit[*,*,2],fit[*,*,3],fit[*,*,4]
	mdsput,path+'DATA','build_signal(build_with_units($,"kev"),*,build_with_units($2,"'+rlab+'"), build_with_units($3,"seconds"), build_with_units($4,"keV"),build_with_units($5,""), build_with_units($6,""))',$
		data[*,*,0],data[*,*,1],time,data[*,*,2],data[*,*,3],data[*,*,4]
	IF inst[0] NE -1 THEN mdsput,path+'INST','build_signal(build_with_units($,"kev"),*,build_with_units($2,"'+rlab+'"), build_with_units($3,"seconds"), build_with_units($4,"keV"))',$
		inst[*,*,0],fit[*,*,1],time,inst[*,*,1] ELSE mdsput,path+'INST','*'
	mdsput,path+'CONFIG', 'build_with_units($,"")',config
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_WRITE_BSOM
;
;PURPOSE:
;	This procedure writes the bsom array information to the tree
;
;OPTIONAL INPUTS:
;	dc	FLOAT	of the shift in rotation [kHz] to set for INST for all channels/times.  Uncertainty is set to 25% of dc value
;
;MODFIICATION HISTORY:
;	Written by:	M.L. Reinke - 7/11/13 (based on HIREXSR_WRITE_BSOM)
;
;-

PRO hirexsr_write_bsom,shot,fit,data,inst,config,time,rlab=rlab,tht=tht,dc=dc
	IF keyword_set(dc) THEN BEGIN
		inst=[[[fit[*,*,0]*0.0+dc]],[[fit[*,*,0]*0.0+0.25*dc]]]
	ENDIF
	IF NOT keyword_set(rlab) THEN rlab=''
	IF keyword_set(tht) THEN rstr='RESULTS'+num2str(tht) ELSE rstr='RESULTS'
	path='\SPECTROSCOPY::TOP.HIREXSR.'+rstr+'.BSOM:'
	mdsopen,'spectroscopy',shot
	mdsput,path+'FIT','build_signal(build_with_units($,"kHz"),*,build_with_units($2,"'+rlab+'"), build_with_units($3,"seconds"), build_with_units($4,"kHz"),build_with_units($5,"kHz/rho"), build_with_units($6,"kHz/rho"))',$
		fit[*,*,0],fit[*,*,1],time,fit[*,*,2],fit[*,*,3],fit[*,*,4]
	mdsput,path+'DATA','build_signal(build_with_units($,"kHz"),*,build_with_units($2,"'+rlab+'"), build_with_units($3,"seconds"), build_with_units($4,"kHz"),build_with_units($5,""), build_with_units($6,""))',$
		data[*,*,0],data[*,*,1],time,data[*,*,2],data[*,*,3],data[*,*,4]
	IF keyword_set(inst) THEN mdsput,path+'INST','build_signal(build_with_units($,"kHz"),*,build_with_units($2,"'+rlab+'"), build_with_units($3,"seconds"), build_with_units($4,"kHz"))',$
		inst[*,*,0],fit[*,*,1],time,inst[*,*,1]
	mdsput,path+'CONFIG', 'build_with_units($,"")',config
	mdsclose,'spectroscopy',shot
END
;+
;NAME
;	HIREXSR_READ_BSTI
;
;PURPOSE:
;	This procedure reads the bspline data from the tree into
;	arrays and can be used to form the BSTI pointer
;
;KEYWORD PARAMETERS:
;	/ptr	this will fill the BSTI optional output
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 6/24/13
;
;-

PRO hirexsr_read_bsti,shot,fit,data,inst,config,time,tht=tht,ptr=ptr,bsti=bsti
	IF NOT keyword_set(tht) THEN tht=0
	IF keyword_set(tht) THEN rstr='RESULTS'+num2str(tht) ELSE rstr='RESULTS'
	path='\SPECTROSCOPY::TOP.HIREXSR.'+rstr+'.BSTI:'
	mdsopen,'spectroscopy',shot
	chk=mdsvalue('_dat='+path+'FIT',/quiet,status=status)
	IF status THEN BEGIN
		x=size(chk)
		fit=fltarr(x[1],x[2],5)
		fit[*,*,0]=chk
		fit[*,*,1]=mdsvalue('dim_of(_dat,0)')
		rlab=mdsvalue('units_of(dim_of(_dat,0))')
		time=mdsvalue('dim_of(_dat,1)')
		fit[*,*,2]=mdsvalue('dim_of(_dat,2)')
		fit[*,*,3]=mdsvalue('dim_of(_dat,3)')
		fit[*,*,4]=mdsvalue('dim_of(_dat,4)')
	ENDIF
	chk=mdsvalue('_dat='+path+'DATA',/quiet,status=status)
	IF status THEN BEGIN
		x=size(chk)
		data=fltarr(x[1],x[2],5)
		data[*,*,0]=chk
		data[*,*,1]=mdsvalue('dim_of(_dat,0)')
		data[*,*,2]=mdsvalue('dim_of(_dat,2)')
		data[*,*,3]=mdsvalue('dim_of(_dat,3)')
		data[*,*,4]=mdsvalue('dim_of(_dat,4)')		
	ENDIF
	config=mdsvalue(path+'CONFIG',/quiet,status=status)
	chk=mdsvalue('_dat='+path+'INST',/quiet,status=status)
	IF status THEN BEGIN
		x=size(chk)
		inst=fltarr(x[1],x[2],2)
		inst[*,*,0]=chk
		inst[*,*,1]=mdsvalue('dim_of(_dat,2)')
        ENDIF
	mdsclose,'spectroscopy',shot
	IF keyword_set(ptr) THEN hirexsr_bsti_arr2ptr,fit,data,config,shot,tht,time,bsti,rlab=rlab	
END

;+
;NAME
;	HIREXSR_READ_BSOM
;
;PURPOSE:
;	This procedure reads the bspline data from the tree into
;	arrays and can be used to form the BSOM pointer
;
;KEYWORD PARAMETERS:
;	/ptr	this will fill the BSOM optional output
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 7/11/13 (based on HIREXSR_READ_BSTI)
;
;-

PRO hirexsr_read_bsom,shot,fit,data,inst,config,time,tht=tht,ptr=ptr,bsom=bsom
	IF NOT keyword_set(tht) THEN tht=0
	IF keyword_set(tht) THEN rstr='RESULTS'+num2str(tht) ELSE rstr='RESULTS'
	path='\SPECTROSCOPY::TOP.HIREXSR.'+rstr+'.BSOM:'
	mdsopen,'spectroscopy',shot
	chk=mdsvalue('_dat='+path+'FIT',/quiet,status=status)
	IF status THEN BEGIN
		x=size(chk)
		fit=fltarr(x[1],x[2],5)
		fit[*,*,0]=chk
		fit[*,*,1]=mdsvalue('dim_of(_dat,0)')
		rlab=mdsvalue('units_of(dim_of(_dat,0))')
		time=mdsvalue('dim_of(_dat,1)')
		fit[*,*,2]=mdsvalue('dim_of(_dat,2)')
		fit[*,*,3]=mdsvalue('dim_of(_dat,3)')
		fit[*,*,4]=mdsvalue('dim_of(_dat,4)')
	ENDIF
	chk=mdsvalue('_dat='+path+'DATA',/quiet,status=status)
	IF status THEN BEGIN
		x=size(chk)
		data=fltarr(x[1],x[2],5)
		data[*,*,0]=chk
		data[*,*,1]=mdsvalue('dim_of(_dat,0)')
		data[*,*,2]=mdsvalue('dim_of(_dat,2)')
		data[*,*,3]=mdsvalue('dim_of(_dat,3)')
		data[*,*,4]=mdsvalue('dim_of(_dat,4)')		
	ENDIF
	config=mdsvalue(path+'CONFIG',/quiet,status=status)
	chk=mdsvalue('_dat='+path+'INST',/quiet,status=status)
	IF status THEN BEGIN
		x=size(chk)
		inst=fltarr(x[1],x[2],2)
		inst[*,*,0]=chk
		inst[*,*,1]=mdsvalue('dim_of(_dat,2)')
        ENDIF
	mdsclose,'spectroscopy',shot
	IF keyword_set(ptr) THEN hirexsr_bsom_arr2ptr,fit,data,config,shot,tht,time,bsom,rlab=rlab	
END

;+
;NAME:
;	HIREXSR_LOAD_RHOGRID
;
;PURPOSE:
;	This function loads the PSINORM grid used in a previous shot
;	in order to use it in analyzing another shot
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 11/20/12
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;
;-

FUNCTION hirexsr_load_rhogrid,shot,line,time,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	hepath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.PROFILES'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.PROFILES'
	CASE line OF 
		0 : path=hepath+'.W:'
		1 : path=hepath+'.X:'
		2 : path=hepath+'.Z:'
		3 : path=hpath+'.LYA1:'
		4 : path=hpath+'.MO4D:'
		5 : path=hpath+'.J:'
		6 : path=hpath+'.LYA1:'		;write wn3 for he-like ca to LYA1
		7 : path=hepath+'.Z:'		;write high-Te lya1 in Z
		8 : path=hepath+'.X:'		;write high-Te 4d in X
		9 : path=hpath+'.LYA1:'		;write lya1 for he-like ca to LYA1
	ENDCASE	
	mdsopen,'spectroscopy',shot
	rho=mdsvalue('dim_of('+path+'PRO,0)',/quiet,status=status)
	tau=mdsvalue('dim_of('+path+'PRO,1)',/quiet,status=status)
	mdsclose,'spectroscopy',shot
	IF status THEN BEGIN
		tmp=where(tau GT 0)
		index=ipt(tau[tmp],time) 
        ENDIF  ELSE index=-1
	IF index[0] NE -1 THEN BEGIN
		irho=rho[*,index]
		nrho=n(irho)+1
		IF irho[0] EQ irho[nrho/2] THEN irho=irho[0:nrho/2-1]
        ENDIF ELSE irho=-1
	RETURN,irho
END
;+
;NAME:
;	HIREXSR_WRITE_PROCONFIG
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees	
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;
;-

PRO hirexsr_write_proconfig,shot,line,eps,eta,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	hepath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.PROFILES'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.PROFILES'
	CASE line OF 
		0 : path=hepath+'.W.CONFIG:'
		1 : path=hepath+'.X.CONFIG:'
		2 : path=hepath+'.Z.CONFIG:'
		3 : path=hpath+'.LYA1.CONFIG:'
		4 : path=hpath+'.MO4D.CONFIG:'
		5 : path=hpath+'.J.CONFIG:'
		6 : path=hpath+'.LYA1.CONFIG:'		;write wn3 for he-like ca to LYA1
		7 : path=hepath+'.Z.CONFIG:'		;write high-Te lya1 in Z
		8 : path=hepath+'.X.CONFIG:'		;write high-Te 4d in X
		9 : path=hpath+'.LYA1.CONFIG:'		;write lya1 for h-like ca to LYA1
	ENDCASE	

	mdsopen,'spectroscopy',shot
	mdsput,path+'EPS', 'build_with_units($,"")',eps
	mdsput,path+'ETA', 'build_with_units($,"")',eta
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_PROCONFIG
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees	
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;
;-

PRO hirexsr_load_proconfig,shot,line,eps,eta,quiet=quiet,status=status,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	hepath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.PROFILES'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.PROFILES'
	CASE line OF 
		0 : path=hepath+'.W.CONFIG:'
		1 : path=hepath+'.X.CONFIG:'
		2 : path=hepath+'.Z.CONFIG:'
		3 : path=hpath+'.LYA1.CONFIG:'
		4 : path=hpath+'.MO4D.CONFIG:'
		5 : path=hpath+'.J.CONFIG:'
		6 : path=hpath+'.LYA1.CONFIG:'		;write wn3 for he-like ca to LYA1
		7 : path=hepath+'.Z.CONFIG:'		;write high-Te lya1 in Z
		8 : path=hepath+'.X.CONFIG:'		;write high-Te 4d in X
		9 : path=hpath+'.LYA1.CONFIG:'		;write lya1 for h-like ca to LYA1
	ENDCASE	
	mdsopen,'spectroscopy',shot
	eps=mdsvalue(path+'EPS',quiet=quiet,status=status)
	eta=mdsvalue(path+'ETA',quiet=quiet,status=status)
	mdsclose,'spectroscopy',shot
END


;+
;NAME:
;	HIREXSR_WRITE_CHECK
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees	
;	9/4/11		M.L. Reinke - added the ability to store subcheck data
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;
;-

PRO hirexsr_write_check,shot,line,check,good,subcheck,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	hepath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.PROFILES'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.PROFILES'
	CASE line OF 
		0 : path=hepath+'.W'
		1 : path=hepath+'.X'
		2 : path=hepath+'.Z'
		3 : path=hpath+'.LYA1'
		4 : path=hpath+'.MO4D'
		5 : path=hpath+'.J'
		6 : path=hpath+'.LYA1'		;write wn3 for he-like ca to LYA1
		7 : path=hepath+'.Z'		;write high-Te lya1 in Z
		8 : path=hepath+'.X'		;write high-Te 4d in X
		9 : path=hpath+'.LYA1'		;write lya1 for h-like ca to LYA1
	ENDCASE	
	totcheck=[[[check]],[[subcheck]]]		;combine check and subcheck into array for storage
	mdsopen,'spectroscopy',shot
	mdsput,path+':CHECK','build_with_units($,"")',totcheck
	mdsput,path+'.CONFIG:GOOD', 'build_with_units($,"")',good
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_CHECK
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees	
;	9/4/11		M.L. Reinke - added the ability to load subcheck data
;	8/25/12		M.L. Reinke - added the line=7,8 paths for  high-Te layout
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;
;-

PRO hirexsr_load_check,shot,line,check,good,subcheck,cs=cs,gs=gs,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	hepath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.PROFILES'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.PROFILES'
	CASE line OF 
		0 : path=hepath+'.W'
		1 : path=hepath+'.X'
		2 : path=hepath+'.Z'
		3 : path=hpath+'.LYA1'
		4 : path=hpath+'.MO4D'
		5 : path=hpath+'.J'
		6 : path=hpath+'.LYA1'		;write wn3 for he-like ca to LYA1
		7 : path=hepath+'.Z'		;write high-Te lya1 in Z
		8 : path=hepath+'.X'		;write high-Te 4d in X
		9 : path=hpath+'.LYA1'		;write lya1 for h-like ca to LYA1
	ENDCASE	
	mdsopen,'spectroscopy',shot
	totcheck=mdsvalue(path+':CHECK',/quiet,status=cs)
	IF cs THEN BEGIN
		check=totcheck[*,*,0:2]
		x=size(totcheck)		
		IF x[3] EQ 3 THEN subcheck=check*0.0 ELSE subcheck=totcheck[*,*,3:5]  	;default subcheck to 0.0 if check has already been stored	
	ENDIF	
	good=mdsvalue(path+'.CONFIG:GOOD',/quiet,status=gs)
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_INVGOOD
;
;PURPOSE:
;	This function loads the GOOD vector used to invert a specic
;	time point for a given SHOT/LINE/THT in order to copy it to
;	another shot
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 11/20/12
;	3/28/14		M.L. Reinke - added the line=9 for H-like Ca
;
;-

FUNCTION hirexsr_load_invgood,shot,line,time,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	hepath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.PROFILES'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.PROFILES'
	CASE line OF 
		0 : path=hepath+'.W'
		1 : path=hepath+'.X'
		2 : path=hepath+'.Z'
		3 : path=hpath+'.LYA1'
		4 : path=hpath+'.MO4D'
		5 : path=hpath+'.J'
		6 : path=hpath+'.LYA1'		;write wn3 for he-like ca to LYA1
		7 : path=hepath+'.Z'		;write high-Te lya1 in Z
		8 : path=hepath+'.X'		;write high-Te 4d in X
		9 : path=hpath+'.LYA1'		;write lya1 for h-like ca to LYA1
	ENDCASE	
	mdsopen,'spectroscopy',shot
	good=mdsvalue(path+'.CONFIG:GOOD',/quiet,status=gs)
	tau=mdsvalue('dim_of('+path+':PRO,1)',/quiet,status=status)
	mdsclose,'spectroscopy',shot
	IF status THEN BEGIN
		tmp=where(tau GT 0)
		index=ipt(tau[tmp],time) 
        ENDIF  ELSE index=-1
	IF index[0] NE -1 THEN BEGIN
		igood=good[*,index]
        ENDIF ELSE igood=-1
	RETURN,igood
END

;+
;NAME:
;	HIREXR_WRITE_RAYTRC
;
;-

PRO hirexsr_write_raytrc,shot,bport,vacves,reduc,hecrys,h=h
	hepath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.RAYTRC:'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.RAYTRC:'
	mdsopen,'spectroscopy',shot
	mdsput,path+'BPORT','build_with_units($,"m")',bport
	mdsput,path+'VACVES','build_with_units($,"m")',vacves
	mdsput,path+'REDUC','build_with_units($,"m")',reduc
	IF keyword_set(h) THEN mdsput,path+'HECRYS','build_with_units($,"m")',hecrys
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_RAYTRC
;
;-

PRO hirexsr_load_raytrc,shot,bport,vacves,reduc,hecrys,h=h
	hepath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.RAYTRC:'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.RAYTRC:'
	mdsopen,'spectroscopy',shot
	bport=mdsvalue(path+'BPORT')
	vacves=mdsvalue(path+'VACVES')
	reduc=mdsvalue(path+'REDUC')
	IF keyword_set(h) THEN hecrys=mdsvalue(path+'HECRYS')
	mdsclose,'spectroscopy',shot
END


;+
;NAME:
;	HIREXSR_UPDATE_TREE
;
;PURPOSE:
;	This procedure adds the HIREXSR tree to a shot which THACO has previously not been deployed
;
;-

PRO hirexsr_update_tree,shot,force=force

	IF keyword_set(force) THEN BEGIN
		mdstcl, "set verify"
		mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
		mdstcl, 'delete node \SPECTROSCOPY::TOP.HIREXSR /noconfirm'
		mdstcl, 'write'
		mdstcl, 'close'
	ENDIF
	
	mdstcl, "set verify"
	mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
	mdstcl, 'add node \SPECTROSCOPY::TOP.HIREXSR'
	mdstcl, 'add node \SPECTROSCOPY::TOP.HIREXSR.ANALYSIS'
	mdstcl, '@/usr/local/cmod/idl/HIREXSR/HIREXSR_WRITE_ANALYSISTREE.TCL'
	mdstcl, 'add node \SPECTROSCOPY::TOP.HIREXSR.INFO'
	mdstcl, '@/usr/local/cmod/idl/HIREXSR/HIREXSR_WRITE_INFOTREE.TCL'
	mdstcl, 'add node \SPECTROSCOPY::TOP.HIREXSR.CALIB'
	mdstcl, '@/usr/local/cmod/idl/HIREXSR/HIREXSR_WRITE_CALIBTREE.TCL'
	mdstcl, 'write'
	mdstcl, 'close'
END

;+
;NAME:
;	HIREXSR_IS_ANALYSIS
;
;PURPOSE:
;	This function is used to check for the presence of ANALYSIS# trees by checking
;	for a populated MORDER node in the ANALYSIS#.HELIKE.MORDER
;
;CALLING SEQUENCE:
;	result=HIREXSR_IS_ANALYSIS(shot,tht)
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 8/19/11
;
;-

FUNCTION hirexsr_is_analysis,shot,tht
	mdsopen,'spectroscopy',shot
	IF tht EQ 0 THEN path='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE:MORDER' ELSE $
		path='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS'+num2str(tht,1)+'.HELIKE:MORDER'
	chk=mdsvalue('getnci('+path+',"usage")',status=status,/quiet)
	mdsclose,'spectroscopy',shot
	IF status THEN output=1 ELSE output=0
	RETURN,output
END

;+
;NAME:
;	HIREXSR_IS_TIRESULTS
;
;PURPOSE:
;	This function is used to check for the presence of RESULTS# trees by checking
;	for the RESULTS#:BSTI:CONFIG node
;
;CALLING SEQUENCE:
;	result=HIREXSR_IS_TIRESULTS(shot,tht)
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 6/24/13
;	7/11/13 	M.L. Reinke - modified naming from '_RESULTS' to '_TIRESULTS'
;       9/22/14		M.L. Reinke - added a check to CONFIG node see if it is filled                             
;
;-

FUNCTION hirexsr_is_tiresults,shot,tht,filled=filled
	mdsopen,'spectroscopy',shot
	IF tht EQ 0 THEN path='\SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSTI:CONFIG' ELSE $
		path='\SPECTROSCOPY::TOP.HIREXSR.RESULTS'+num2str(tht,1)+'.BSTI:CONFIG'
	chk=mdsvalue('getnci('+path+',"usage")',status=status,/quiet)
	IF status THEN config=mdsvalue(path,status=filled,/quiet) ELSE filled=status			;if node exists, check if data
	mdsclose,'spectroscopy',shot
	IF status THEN output=1 ELSE output=0
	RETURN,output
END

;+
;NAME:
;	HIREXSR_IS_OMRESULTS
;
;PURPOSE:
;	This function is used to check for the presence of RESULTS# trees by checking
;	for the RESULTS#:BSOM:CONFIG node
;
;CALLING SEQUENCE:
;	result=HIREXSR_IS_OMRESULTS(shot,tht)
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 7/11/13
;       9/22/14		M.L. Reinke - added a check to CONFIG node see if it is filled                             
;
;-

FUNCTION hirexsr_is_omresults,shot,tht,filled=filled
	mdsopen,'spectroscopy',shot
	IF tht EQ 0 THEN path='\SPECTROSCOPY::TOP.HIREXSR.RESULTS.BSOM:CONFIG' ELSE $
		path='\SPECTROSCOPY::TOP.HIREXSR.RESULTS'+num2str(tht,1)+'.BSOM:CONFIG'
	chk=mdsvalue('getnci('+path+',"usage")',status=status,/quiet)
	IF status THEN config=mdsvalue(path,status=filled,/quiet) ELSE filled=status			;if node exists, check if data
	mdsclose,'spectroscopy',shot
	IF status THEN output=1 ELSE output=0
	RETURN,output
END

;+
;NAME:
;	HIREXSR_REMOVE_ANALYSIS
;
;PURPOSE:
;	This procedure deletes a numbered analysis tree (\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS#).
;	The code will stop and require a .c input from the user to cofirm deletion.
;
;CALLING SEQUENCE:
;	HIREXSR_REMOVE_ANALYSIS,shot,tht
;
;INPUTS:
;	shot	LONG	shot number
;	tht	INT	THACO tree number
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 8/19/11
;
;-

PRO hirexsr_remove_analysis,shot,tht
	status=hirexsr_is_analysis(shot,tht)
	IF NOT status THEN BEGIN
		print,'NO ANALYSIS'+num2str(tht,1)+' PRESENT'
		RETURN
	ENDIF ELSE BEGIN
		print, 'CONTINUE TO DELETE '+num2str(tht,1)+' TREE'
		stop
		mdstcl, "set verify"
		mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
		mdstcl, 'delete node \SPECTROSCOPY::TOP.HIREXSR.ANALYSIS'+num2str(tht,1)+' /noconfirm'
		mdstcl, 'write'
		mdstcl, 'close'	
	ENDELSE

END

;+
;NAME:
;	HIREXSR_REMOVE_RESULTS
;
;PURPOSE:
;	This procedure deletes a numbered results tree (\SPECTROSCOPY::TOP.HIREXSR.RESULTS#).
;	The code will stop and require a .c input from the user to cofirm deletion.
;
;CALLING SEQUENCE:
;	HIREXSR_REMOVE_RESULTS,shot,tht
;
;INPUTS:
;	shot	LONG	shot number
;	tht	INT	THACO tree number
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 6/24/13
;	7/11/13		M.L. Reinke - added a check to see if either Ti or Om RESULTS existed
;-

PRO hirexsr_remove_results,shot,tht
	tistat=hirexsr_is_tiresults(shot,tht)
	omstat=hirexsr_is_omresults(shot,tht)
	IF NOT tistat AND NOT omstat THEN BEGIN
		print,'NO RESULTS'+num2str(tht,1)+' PRESENT'
		RETURN
	ENDIF ELSE BEGIN
		print, 'CONTINUE TO DELETE '+num2str(tht,1)+' TREE'
		stop
		mdstcl, "set verify"
		mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
		mdstcl, 'delete node \SPECTROSCOPY::TOP.HIREXSR.RESULTS'+num2str(tht,1)+' /noconfirm'
		mdstcl, 'write'
		mdstcl, 'close'	
	ENDELSE

END

;+
;NAME:
;	HIREXSR_ADDNEW_ANALYSIS
;
;PURPOSE:
;	This procedure adds a new numbered analysis tree (\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS#)
;	for use in parallel runs of THACO
;
;CALLING SEQUENCE
;	shot	LONG	shot number
;
;OPTIONAL INPUTS:
;	tht	INT	THACO tree number DEFAULT: the next available number
;
;KEYWORD PARAMETERS:
;	/force 	this will force an overwrite of a THACO tree that already exists
;
;OPTIONAL OUTPUTS:
;	chk	INT	indicating result of process:
;			  0 - THT existed andnothing done
;			 -1 THT existed and overwritten
;			  1 THT did not exist and was added 
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 8/19/11
;	M.L. Reinke 	4/28/16 - updated the /force command to use
;                                 the /confirm instead of /noconfirm to automate delition and added removal of RESULTS nodes.
;-	

PRO hirexsr_addnew_analysis,shot,tht=tht,force=force,chk=chk,ftht=ftht,setup=setup
	lownum=0
	status=1
	WHILE lownum LT 20 AND status EQ 1 DO BEGIN
		lownum+=1
		status=hirexsr_is_analysis(shot,lownum)
		;print, lownum,status
	ENDWHILE
	IF keyword_set(tht) THEN BEGIN
		IF tht LT lownum THEN BEGIN			
			IF keyword_set(force) THEN BEGIN	;if specified and /force, delete existint ANALYSISX tree
				mdstcl, "set verify"
				mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
				mdstcl, 'delete node \SPECTROSCOPY::TOP.HIREXSR.ANALYSIS'+num2str(tht,1)+' /confirm'
				mdstcl, 'write'
				mdstcl, 'close'
				tistat=hirexsr_is_tiresults(shot,tht)
				omstat=hirexsr_is_omresults(shot,tht)
				IF tistat OR omstat THEN BEGIN
					mdstcl, "set verify"
					mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
					mdstcl, 'delete node \SPECTROSCOPY::TOP.HIREXSR.RESULTS'+num2str(tht,1)+' /confirm'
					mdstcl, 'write'
					mdstcl, 'close'
				ENDIF
			ENDIF ELSE BEGIN
				print, 'HIREXSR ANALYSIS'+num2str(tht,1)+' ALREADY EXISTS'
				chk=0
				RETURN
			ENDELSE
		ENDIF 				;if specified and available then use 
	ENDIF ELSE tht=lownum			;if not specified default to next available
	hstr='ANALYSIS'+num2str(tht,1)

	mdstcl, "set verify"
	mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
	mdstcl, 'add node \SPECTROSCOPY::TOP.HIREXSR.'+hstr
	filename='/usr/local/cmod/idl/HIREXSR/HIREXSR_WRITE_ANALYSISTREE.TCL'
	openr,lun,filename,/get_lun
	line=strarr(1)
	readf,lun,line
	WHILE line NE 'WRITE' DO BEGIN
		command=strsplit(line,'ANALYSIS',/extract,/regex)
		newline=command[0]+hstr+command[1]
		mdstcl,newline
		readf,lun,line
	ENDWHILE
	mdstcl, 'write'
	mdstcl, 'close'
	IF NOT keyword_set(ftht) THEN ftht=0			;fill from default tree
	hirexsr_load_morder,shot,morder,tht=ftht				
	hirexsr_write_morder,shot,morder=morder,tht=tht
	hirexsr_load_morder,shot,morder,/h,tht=ftht			
	hirexsr_write_morder,shot,morder=morder,/h,tht=tht
	IF keyword_set(setup) THEN BEGIN
		hirexsr_load_binning,shot,chmap,tch,tmap,good,chmax,tht=ftht
		IF tch[0] EQ 0 THEN tch[0]= -1
		hirexsr_write_binning,shot,chmap=chmap,tch=tch,tmap=tmap,good=good,chmax=chmax,tht=tht
		hirexsr_load_binning,shot,chmap,tch,tmap,good,chmax,tht=ftht,/h
		IF tch[0] EQ 0 THEN tch[0]= -1
		hirexsr_write_binning,shot,chmap=chmap,tch=tch,tmap=tmap,good=good,chmax=chmax,tht=tht,/h
		print,'morder and binning information copied from THT='+num2str(ftht,1)
	ENDIF
	free_lun,lun
	print, hstr+' TREE added to SHOT: '+num2str(shot,1)
	IF keyword_set(force) THEN chk=-1 ELSE chk=1
END

;+
;NAME:
;	HIREXSR_ADDNEW_TIRESULTS
;
;PURPOSE:
;	This procedure adds a new default or numbered results tree (\SPECTROSCOPY::TOP.HIREXSR.RESULTS#:BSTI)
;	for use in parallel runs of THACO
;
;CALLING SEQUENCE:
;	HIREXSR_ADDNEW_TIRESULTS,shot
;
;INPUTS:
;	shot	LONG	shot number
;
;OPTIONAL INPUTS:
;	tht	INT	THACO tree number DEFAULT: the next available number
;
;OPTIONAL OUTPUTS:
;	chk	INT	indicating result of process:
;			  0 - THT existed andnothing done
;			  1 THT did not exist and was added 
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 6/24/13
;	7/11/13:	M.L. Reinke - modified to check for BSOM
;
;-

PRO hirexsr_addnew_tiresults,shot,tht=tht,chk=chk
	IF NOT keyword_set(tht) THEN tht=0
	tistat=hirexsr_is_tiresults(shot,tht)
	omstat=hirexsr_is_omresults(shot,tht)

	IF tistat THEN BEGIN
		print, 'HIREXSR RESULTS'+num2str(tht,1)+':BSTI ALREADY EXISTS'
		chk=0
		RETURN
	ENDIF 
	IF tht NE 0 THEN hstr='RESULTS'+num2str(tht,1) ELSE hstr='RESULTS'
	
	mdstcl, "set verify"
	mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
	IF NOT omstat THEN mdstcl, 'add node \SPECTROSCOPY::TOP.HIREXSR.'+hstr
         
	filename='/usr/local/cmod/idl/HIREXSR/HIREXSR_WRITE_TIRESULTSTREE.TCL'
	openr,lun,filename,/get_lun
	line=strarr(1)
	readf,lun,line
	WHILE line NE 'WRITE' DO BEGIN
		command=strsplit(line,'RESULTS',/extract,/regex)
		newline=command[0]+hstr+command[1]
		mdstcl,newline
		readf,lun,line
	ENDWHILE
	mdstcl, 'write'
	mdstcl, 'close'
	free_lun,lun
	print, hstr+' TREE added to SHOT: '+num2str(shot,1)
	chk=1
END

;+
;NAME:
;	HIREXSR_ADDNEW_OMRESULTS
;
;PURPOSE:
;	This procedure adds a new default or numbered results tree (\SPECTROSCOPY::TOP.HIREXSR.RESULTS#:BSOM)
;	for use in parallel runs of THACO
;
;CALLING SEQUENCE:
;	HIREXSR_ADDNEW_OMRESULTS,shot
;
;INPUTS:
;	shot	LONG	shot number
;
;OPTIONAL INPUTS:
;	tht	INT	THACO tree number DEFAULT: the next available number
;
;OPTIONAL OUTPUTS:
;	chk	INT	indicating result of process:
;			  0 - THT existed andnothing done
;			  1 THT did not exist and was added 
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 7/11/13
;
;-

PRO hirexsr_addnew_omresults,shot,tht=tht,chk=chk
	IF NOT keyword_set(tht) THEN tht=0
	tistat=hirexsr_is_tiresults(shot,tht)
	omstat=hirexsr_is_omresults(shot,tht)

	IF omstat THEN BEGIN
		print, 'HIREXSR RESULTS'+num2str(tht,1)+':BSOM ALREADY EXISTS'
		chk=0
		RETURN
	ENDIF 
	IF tht NE 0 THEN hstr='RESULTS'+num2str(tht,1) ELSE hstr='RESULTS'
	
	mdstcl, "set verify"
	mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
	IF NOT tistat THEN mdstcl, 'add node \SPECTROSCOPY::TOP.HIREXSR.'+hstr
         
	filename='/usr/local/cmod/idl/HIREXSR/HIREXSR_WRITE_OMRESULTSTREE.TCL'
	openr,lun,filename,/get_lun
	line=strarr(1)
	readf,lun,line
	WHILE line NE 'WRITE' DO BEGIN
		command=strsplit(line,'RESULTS',/extract,/regex)
		newline=command[0]+hstr+command[1]
		mdstcl,newline
		readf,lun,line
	ENDWHILE
	mdstcl, 'write'
	mdstcl, 'close'
	free_lun,lun
	print, hstr+' TREE added to SHOT: '+num2str(shot,1)
	chk=1
END

;+
;NAME:
;	HIREXSR_COPY_INFO
;
;-

PRO hirexsr_copy_info,fromshot,toshot,module=module
	IF NOT keyword_set(module) THEN module=[1,2,3,4]
	FOR j=0,n(module) DO BEGIN
		i=module[j]
		info=hirexsr_load_info(i,shot=fromshot,/tree)
		hirexsr_load_info2tree,toshot,i,info=info
		hirexsr_load_break,fromshot,i,break
		hirexsr_write_break,toshot,i,break
	ENDFOR
	print, 'INFO files and BREAK copied from '+num2str(fromshot,1)+' to '+num2str(toshot,1)
END

;+
;NAME:
;	HIREXSR_COPY_COMMENT
;
;-

PRO HIREXSR_COPY_COMMENT, fromshot, toshot, module

path = '\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD'+strcompress(string(module), /remove_all)+':COMMENT'
;mdsopen, 'spectroscopy', fromshot
;val = mdsvalue(path, /quiet, status = status)
;mdsclose
; check that comment exists
;if status EQ 1 then begin
;	subs = strsplit(val, ',', /extract)
;	if n_elements(subs) GE 3 then begin
;		subs(2) = STRCOMPRESS(STRING(fromshot), /remove_all)
;		new_output = STRJOIN(subs,',', /single)
;	endif else begin
;		if n_elements(subs) LT 2 then begin
;			user = logname()
;			date_time = systime()
;			new_output = strcompress(string(user_name), /remove_all)+','+date_time+','+STRCOMPRESS(STRING(fromshot), /remove_all)
;		endif else begin
;			new_output = val + ',' +STRCOMPRESS(STRING(fromshot), /remove_all)
;		endelse	
;	endelse
;endif else begin
;	user = logname()
;	date_time = systime()
;	new_output = strcompress(string(user_name), /remove_all)+','+date_time+'

user = logname()
date_time = systime()
new_output = strcompress(string(user_name), /remove_all)+','+date_time+','+STRCOMPRESS(STRING(fromshot), /remove_all)

mdsopen, 'spectroscopy', toshot
test = mdsvalue(path, /quiet, status = status)		
;265388144 - error for no node
;265388258 - error for no data in node

if status EQ 265388144 then begin
	mdstcl, "set verify"
	mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
	mdstcl, 'add node /usage=text '+path
	mdstcl, 'write'
	mdstcl, 'close'
endif

	mdsopen, 'spectroscopy', shot
	mdsput, path, '$', output
	mdsclose
END


;+
;NAME:
;	HIREXSR_COPY_CALIB
;	
;MODIFICATION HISTORY:
;	9/4/11		M.L. Reinke - added INST to be copied from target shot if available
;
;-

PRO hirexsr_copy_calib,fromshot,toshot,module=module

	;write from shot information for the calibration to be copied
	mdsopen,'spectroscopy',fromshot
	times=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB:TIMES')
	mdsclose,'spectroscopy',fromshot
	mdsopen,'spectroscopy',toshot
	mdsput,'\SPECTROSCOPY::TOP.HIREXSR.CALIB:TIMES','build_with_units($,"sec")',times
	mdsput,'\SPECTROSCOPY::TOP.HIREXSR.CALIB:SHOT','build_with_units($,"")',fromshot
	mdsclose,'spectroscopy',toshot	

	IF NOT keyword_set(module) THEN module=[1,2,3,4]
	FOR i=0,n(module) DO BEGIN
		hirexsr_load_lambda,fromshot,module[i],lambda
		hirexsr_write_lambda,toshot,module[i],lambda=lambda
		hirexsr_load_pos,fromshot,module[i],pos
		hirexsr_write_pos,toshot,module[i],pos=pos
		hirexsr_load_etendue,fromshot,module[i],u
		hirexsr_write_etendue,toshot,module[i],u=u
		hirexsr_load_trans,fromshot,module[i],trans
		hirexsr_write_trans,toshot,module[i],trans=trans
		hirexsr_load_white,fromshot,module[i],white
		hirexsr_write_white,toshot,module[i],white=white
		hirexsr_load_calib_inst,fromshot,module[i],iwidth,ishift,status=status
		IF status THEN hirexsr_write_calib_inst,toshot,module[i],iwidth,ishift
		IF status THEN mvstr='LAMBDA,POS,U,TRANS,WHTFLD,INST' ELSE mvstr='LAMBDA,POS,U,TRANS,WHTFLD'
		print, 'MOD'+num2str(module[i],1)+' '+mvstr+' copied from '+num2str(fromshot,1)+' to '+num2str(toshot,1)
	ENDFOR

END

;+
;NAME:
;	HIREXSR_COPY_BINNING
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees	
;
;-

PRO hirexsr_copy_binning,fromshot,toshot,chmap=chmap,tmap=tmap,good=good,h=h,fromtht=fromtht,totht=totht
	hirexsr_load_binning,fromshot,fchmap,ftch,ftmap,fgood,fchmax,h=h,tht=fromtht
	IF ftch[0] EQ 0 THEN ftch[0]= -1
	IF keyword_set(chmap) THEN hirexsr_write_binning,toshot,chmap=fchmap,tch=ftch,chmax=fchmax,h=h,tht=totht
	IF keyword_set(tmap) THEN hirexsr_write_binning,toshot,tmap=ftmap,h=h,tht=totht
	IF keyword_set(good) THEN hirexsr_write_binning,toshot,good=fgood,h=h,tht=totht
END

;+
;NAME:
;	HIREXSR_COPY_SETUP
;
;MODIFICATION HISTORY:
;	8/19/11		M.L. Reinke - added optional input THT to work with numbered analysis trees	
;
;-
PRO hirexsr_copy_setup,fshot,tshot,h=h,fromtht=fromtht,totht=totht
	hirexsr_load_morder,fshot,module,h=h						;load the target morder and copy morder 
	FOR i=0,n(module) DO hirexsr_copy_info,fshot,tshot,module=module[i]
	FOR i=0,n(module) DO hirexsr_copy_calib,fshot,tshot,module=module[i]
	hirexsr_write_morder,tshot,h=h,morder=morder
	hirexsr_copy_binning,fshot,tshot,/chmap,/tmap,/good,h=h,fromtht=fromtht,totht=totht
END

;+
;NAME:
;	HIREXSR_WRITE_MORDER
;
;-

PRO hirexsr_write_morder,shot,h=h,morder=morder,tht=tht
	IF NOT keyword_set(morder) THEN BEGIN
		IF keyword_set(h) THEN morder=[4] ELSE morder=[1,2,3]
	ENDIF
	mdsopen,'spectroscopy',shot
	IF keyword_set(tht) THEN hstr='ANALYSIS'+num2str(tht,1) ELSE hstr='ANALYSIS'
	hepath='\SPECTROSCOPY::TOP.HIREXSR.'+hstr+'.HELIKE.MORDER'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.'+hstr+'.HLIKE:MORDER'
	IF keyword_set(h) THEN path=hpath ELSE path=hepath
	mdsput,path,'build_with_units($,"")',morder
	mdsclose,'spectroscopy',shot
END

;NAME:
;	HIREXSR_LOAD_MORDER
;
;-

PRO hirexsr_load_morder,shot,morder,h=h,tht=tht
	IF keyword_set(tht) THEN hstr='ANALYSIS'+num2str(tht,1) ELSE hstr='ANALYSIS'
	mdsopen,'spectroscopy',shot
	hepath='\SPECTROSCOPY::TOP.HIREXSR.'+hstr+'.HELIKE.MORDER'
	hpath='\SPECTROSCOPY::TOP.HIREXSR.'+hstr+'.HLIKE:MORDER'
	IF keyword_set(h) THEN path=hpath ELSE path=hepath
	morder=mdsvalue(path)
	mdsclose,'spectroscopy',shot	
END

;+
;NAME:
;	HIREXSR_WRITE_ANALYSIS_COMMENT
;
;-

PRO hirexsr_write_analysis_comment,shot,comment,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht,1) ELSE astr='ANALYSIS'
	path='\SPECTROSCOPY::TOP.HIREXSR.'+astr+':COMMENT'
	hirexsr_add_analysis_comment,shot,tht=tht,/quiet			;will add COMMENT node if it's not already there
	mdsopen,'spectroscopy',shot
	mdsput,path,'$',comment	
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_LOAD_ANALYSIS_COMMENT
;
;-

PRO hirexsr_load_analysis_comment,shot,comment,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht,1) ELSE astr='ANALYSIS'
	path='\SPECTROSCOPY::TOP.HIREXSR.'+astr+':COMMENT'
	mdsopen,'spectroscopy',shot
	comment=mdsvalue(path,/quiet,status=status)
	mdsclose,'spectroscopy',shot
	IF NOT status THEN comment='error'
END



;+
;NAME:
; HIREXSR_GET_LINT_PROFILE
;-
function hirexsr_get_lint_profile, shot, line = line
c = 2.99e8
e = 1.602e-19
mconv=1.661e-27	
NUMBER_CHANNELS = 48 ; bad hardcode but good for expediency
if not keyword_set(line) then line =2 ; take z line
	
mdsopen,'analysis',shot
; grab the efit stuff
MdsOpen,'analysis',shot
MdsSetDefault,'\efit_aeqdsk'
rmid = mdsvalue ('\ANALYSIS::EFIT_RMID')
efit_time = mdsvalue('dim_of(\ANALYSIS::EFIT_RMID,0)') ; time of rmid measurement
rpsi = mdsvalue('dim_of(\ANALYSIS::EFIT_RMID,1)') ; psi values of different rpsi values
mdsclose

;r_proj = multi_interpol(rmid, rpsi,efit_time, subpsinorm, subtime); get rmid values at sr times


hirexsr_load_momentptr,shot,line,mom,tau,pos,tpos,lam_o,z
valindices = where(tau NE -1)
minindex = valindices(0)
maxindex = valindices(n_elements(valindices)-1)
subtime = tau(where(tau NE -1))

mass=read_atomic_mass(z)
conv_factor=(lam_o/c)^2*(e*1.0e3/(mass*mconv)) ;conversion factor for 
v = fltarr(NUMBER_CHANNELS, n_elements(valindices))
verr = fltarr(NUMBER_CHANNELS, n_elements(valindices))
ti = fltarr(NUMBER_CHANNELS, n_elements(valindices))
tierr = fltarr(NUMBER_CHANNELS, n_elements(valindices))
rhotang = fltarr(NUMBER_CHANNELS, n_elements(valindices))
FOR j=minindex, maxindex DO BEGIN
    imom=*mom[j]
    IF imom[0] NE -1 THEN BEGIN
        pindex=last(where(tpos LE tau[j]))
        jpos=pos[*,*,pindex]
        v[*, j]=-1.0*(lam_o-imom[*,11])*c/lam_o/(2.0*!pi*jpos[2,*]*cos(jpos[3,*]))*1.0e-3
        verr[*, j]=imom[*,14]*c/lam_o/(2.0*!pi*jpos[2,*]*cos(jpos[3,*]))*1.0e-3
        ti[*, j]=(imom[*,12])^2/conv_factor
        tierr[*,j]=2.0*imom[*,15]*sqrt(imom[*,12]^2)/conv_factor
        rhotang[*, j]=imom[*,16]
        ;ilint=[[v],[verr],[ti],[tierr],[rhotang]]
        ;*lint[i,j]=ilint
    ENDIF
ENDFOR			
heap_free,mom

r_proj = multi_interpol(rmid, rpsi,efit_time, rhotang, subtime); get rmid values at sr times
r_ave =  total(r_proj, 2)/n_elements(r_proj(0, *))
retstr = {shot:shot, line:line, v:v, verr:verr, ti:ti, tierr:tierr, rhotang:rhotang, r_proj:r_proj, r_ave:r_ave, time:subtime}
return, retstr
end


;+
;NAME:
;	multi_interpol projects ne_tsc (with r_tsc, t_tsc) onto the r_major, times grid
;-
function multi_interpol, rmid, rpsi, efit_time, psinorm, times

xdim = n_elements(psinorm(*, 0))
ydim = n_elements(psinorm(0, *))
r_proj = fltarr(xdim, ydim)

for index = 0, n_elements(times) - 1 do begin
    temp_sr_r = psinorm (*,index) ; select the psi values of sr at this time point
    temp_sr_t = times(index)
    index_t = where(temp_sr_t LT efit_time) ; find where the shot lies in time
    index_t = index_t(0) ; select the lower index
    time_frac = (temp_sr_t-efit_time(index_t-1))/(efit_time(index_t) - efit_time(index_t-1)) ; fraction of the way to next index
    
    ; now get the averaged r, and ne values on the grid
    ;rn_low = rpsi(index_t, *)
    ;rn_high = rpsi(index_t+1, *)

    ;project the ne values onto the time point of the velocity measurement
    ;r_project = rn_low + time_frac*(rn_high-rn_low) 
    ne_project = rmid(index_t, *) + time_frac* (rmid( index_t, *)-rmid(index_t-1, *))
    
    ; interpolate the ne values on the r values of the velocity measurements
    r_proj(*, index) = interpol (Reform(ne_project), rpsi, temp_sr_r)
    
endfor


return, r_proj
end



;+
;NAME:
;	HIREXSR_GET_PROFILE
;PURPOSE:
;	This function returns the inverted profile from the tree and
;	transforms it into r/a space
;
;CALLING SEQUENCE
;	INVERSION = HIREXSR_GET_PROFILE(shot)
;
;INPUTS:
;	shot        LONG  shot number
;
;OUPUTS:
;       Returns struct of inversion or -1
;NOTES:
;	quiet keyword has not been implemented yet
;
;
;MODIFICATION HISTORY:
;	Written by: 	Y. Podpaly  02/03/2011
;                       Y.P 05/25/2011 to deal with m=1
;			Y.P 06/06/2011 to add moly, lya1
;                       Y.P 08/30/2011 adding dc_shift and adding shot
;                       to outputted struct
;			Y.P. 08/30 added tht capability
;			Y.P 09/14 dc_shift tied to omega
;			Y.P 09/16 override command added to remove safety
;			Y.P. 1/16 added time keyword
;-

FUNCTION HIREXSR_GET_PROFILE, shot, w = w, moly = moly, lya1 = lya1, quiet = quiet, dc_shift = dc_shift, tht = tht, override = override, seltime= seltime

if keyword_set(dc_shift) then dc_shift = dc_shift else dc_shift= 0.0

if keyword_set(tht) then begin
      if tht EQ 0 then begin
      	 initstring = '\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.'
	 endif else begin
	       initstring = '\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS' + STRCOMPRESS(STRING(tht), /remove_all) + '.'
	 endelse 
endif else begin
      initstring = '\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.'
endelse


mdsopen, 'spectroscopy', shot
mod1=mdsvalue('getnci(\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD1, "length")')
mod2=mdsvalue('getnci(\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD2, "length")')
mod3=mdsvalue('getnci(\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD3, "length")')
mod4=mdsvalue('getnci(\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD4, "length")')

print,'image sizes: '
print,mod1,mod2,mod3,mod4
IF mod1-mod2 NE 0 OR mod1-mod3 NE 0 OR mod1-mod4 NE 0 THEN BEGIN
    print, 'Likely module malfunction on this shot'
    if not keyword_set(override) then return, -1
ENDIF

; get the data
if not keyword_set(lya1) AND not keyword_set(moly) then begin
	if not keyword_set(w) then begin
	    invstruc = mdsvalue(initstring+'HELIKE.PROFILES.Z:PRO')
	    invstrucerr = mdsvalue(initstring + 'HELIKE.PROFILES.Z:PROERR')
	    psinorm = mdsvalue('dim_of('+initstring+'HELIKE.PROFILES.Z:PRO, 0)')
	    time = mdsvalue('dim_of('+initstring+'HELIKE.PROFILES.Z:PRO, 1)')
	    lineid = 'z'
	endif else begin
	    invstruc = mdsvalue(initstring+'HELIKE.PROFILES.W:PRO')
	    invstrucerr = mdsvalue(initstring+'HELIKE.PROFILES.W:PROERR')
	    psinorm = mdsvalue('dim_of('+initstring+'HELIKE.PROFILES.W:PRO, 0)')
	    time = mdsvalue('dim_of('+initstring+'HELIKE.PROFILES.W:PRO, 1)')
	    lineid = 'w'
	endelse
endif else begin
if keyword_set(lya1) then begin
    invstruc = mdsvalue(initstring+'HLIKE.PROFILES.LYA1:PRO')
    invstrucerr = mdsvalue(initstring+'HLIKE.PROFILES.LYA1:PROERR')
    psinorm = mdsvalue('dim_of('+initstring+'HLIKE.PROFILES.LYA1:PRO, 0)')
    time = mdsvalue('dim_of('+initstring+'HLIKE.PROFILES.LYA1:PRO, 1)')
    lineid = 'lya1'
endif

if keyword_set(moly) then begin
    invstruc = mdsvalue(initstring+'HLIKE.PROFILES.MO4D:PRO')
    invstrucerr = mdsvalue(initstring+'HLIKE.PROFILES.MO4D:PROERR')
    psinorm = mdsvalue('dim_of('+initstring+'HLIKE.PROFILES.MO4D:PRO, 0)')
    time = mdsvalue('dim_of('+initstring+'HLIKE.PROFILES.MO4D:PRO, 1)')
    lineid = 'moly4d'
endif
endelse

mdsclose

emissstr = invstruc(*, where(time NE -1), 0)
emissstrerr = invstruc(*, where(time NE -1), 0)
tistr = invstruc(*, where(time NE -1), 3)
tistrerr = invstrucerr(*, where(time NE -1), 3)
omgstr = invstruc(*, where(time NE -1), 1)+dc_shift
omgstrerr = invstrucerr(*, where(time NE -1), 1)
subtime = time(where(time NE -1))
subpsinorm = psinorm(*, where(time NE -1))

; select only times of interest
if keyword_set(seltime) then begin
   time_low = seltime(0)
   time_high = seltime(1)
   indlow = where(subtime GT time_low)
   indlow = indlow(0)
   indhigh = where(subtime GT time_high)
   indhigh = indhigh(0)

   emissstr = emissstr(*, indlow:indhigh)
   emissstrerr = emissstrerr(*, indlow:indhigh)
   omgstr = omgstr(*, indlow:indhigh)
   omgstrerr = omgstrerr(*, indlow:indhigh)
   tistr = tistr(*,indlow:indhigh )
   tistrerr = tistrerr(*, indlow:indhigh)
   subtime = subtime(indlow:indhigh)
   subpsinorm = subpsinorm(*,indlow:indhigh)
endif

; grab the efit stuff
MdsOpen,'analysis',shot
MdsSetDefault,'\efit_aeqdsk'
rmid = mdsvalue ('\ANALYSIS::EFIT_RMID')
efit_time = mdsvalue('dim_of(\ANALYSIS::EFIT_RMID,0)') ; time of rmid measurement
rpsi = mdsvalue('dim_of(\ANALYSIS::EFIT_RMID,1)') ; psi values of different rpsi values
mdsclose
;stop
r_proj = multi_interpol(rmid, rpsi,efit_time, subpsinorm, subtime); get rmid values at sr times

; now create km/s

rotstr = omgstr
rotstrerr = omgstrerr
for index = 0, n_elements(subtime) -1 do begin 
    rotstr(*, index) = 2.*!Pi * r_proj(*, index)*(omgstr(*, index) )
    rotstrerr(*, index) = rotstr(*, index)*(omgstrerr(*, index)/omgstr(*, index)) 
endfor
r_ave =  total(r_proj, 2)/n_elements(r_proj(0, *))

; we now have all of the data in proper form
; m=1 needs to be checked here
; m=1 causes the spatial array to be doubled

numspatial = n_elements(subpsinorm(*, 0))

; not very elegant
logictest1 =  subpsinorm(0, *) EQ subpsinorm(numspatial/2, *) ; check element 0
logictest2 =  subpsinorm(1, *) EQ subpsinorm(1+numspatial/2, *) ; check element 1

if (sum(logictest1) EQ n_elements(logictest1)) AND(sum(logictest2) EQ n_elements(logictest2)) then begin
	print, 'Found m=1 components'
	finindex = numspatial/2 -1
	startindex = numspatial/2
	endindex = numspatial-1

	emiss = emissstr(0:finindex, *)
	emissm1 = emissstr(startindex:endindex, *)
	emisserr = emissstrerr(0:finindex, *)
	emisserrm1 = emissstrerr(startindex:endindex,*)

	ti = tistr(0:finindex, *)
	tim1 = tistr(startindex:endindex, *)
	tierr = tistrerr(0:finindex, *)
	tierrm1 = tistrerr(startindex:endindex,*)	

	omg = omgstr(0:finindex, *)
	omgm1 = omgstr(startindex:endindex, *)
	omgerr = omgstrerr(0:finindex, *)
	omgerrm1 = omgstrerr(startindex:endindex,*)

	rot = rotstr(0:finindex, *)
	rotm1 = rotstr(startindex:endindex, *)
	roterr = rotstrerr(0:finindex, *)
	roterrm1 = rotstrerr(startindex:endindex,*)

	r_maj = r_proj(0:finindex, *)
	psi = subpsinorm(0:finindex, *)
	r_ave = r_ave(0:finindex, *)

	inversion = {emiss:emiss, emissm1:emissm1, emisserr:emisserr,emisserrm1:emisserrm1, ti:ti, tim1:tim1, tierr:tierr,tierrm1:tierrm1, rot:rot, rotm1:rotm1, roterr:roterr, roterrm1:roterrm1, omg:omg, omgm1:omgm1, omgerr:omgerr, omgerrm1:omgerrm1, r_maj:r_maj, time:subtime, psi:psi, r_ave:r_ave, lineid:lineid, shot:shot}
	return, inversion
endif else begin
	
	inversion = {emiss:emissstr, emisserr:emissstrerr, ti:tistr, tierr:tistrerr, rot:rotstr, roterr:rotstrerr, omg:omgstr, omgerr:omgstrerr, r_maj:r_proj, time:subtime, psi:subpsinorm, r_ave:r_ave, lineid:lineid, shot:shot}
	return, inversion
endelse

	return, -1; should never get here

end

