;+
;NAME:
;	GAUSSIAN_FITS
;
;PURPOSE:
;	This function calculates the gaussian line profile for the sum of 
;	an arbitrary number of gaussians plus an optional DC offset.  The format
;	is made to be used with MPFIT
;
;CALLING_SEQUENCE:
;	result=GAUSSIAN_FITS(x,p)
;
;INPUTS:
;	x	FLTARR	[n_x] of the points for which to calculate the spectra
;	p	FLTARR 	[n_gauss*3+nb] where nb will determine the order of the baseline fit
;
;OUTPUTS:
;	result:	FLTARR  [n_x] of the sum of all the gaussians at each point
;
;PROCEDURE:
;	Gaussians are specified as a*exp(-(x-b)^2/(2*c^2)).  If the parameter values of p
;	are a double precision array, then the result is returned as a double.
;	
;	IF n_elements(p) MOD 3 = 0 then the last three are assume to be a quadraditc baseline
;	IF n_elements(p) MOD 3 = 1 then the last three are assume to be a constant baseline
;	IF n_elements(p) MOD 3 = 2 then the last three are assume to be a linear baseline
;
;MOFICATION HISTORY:
;	Written by 	A. Ince-Cushman 2007
;	1/10		M.L. Reinke - modified to auto detect for DP or SP
;	5/13/11		M.L. Reinke - modified to allow for higher order baseline terms in the fits
;				      currently hard coded for use with the zjk fits
;	4/17/12		M.L. Reinke - noticed that n_elements(p) was called multiple times and 
;				      replaced with a single call to be more efficient
;
;-

FUNCTION gaussian_fits, x, p,base=base
	np=n_elements(p)
	CASE np MOD 3 OF
		0 : BEGIN
			IF np EQ 3 THEN BEGIN
				n_line=1
				basen=3
			ENDIF ELSE BEGIN		;quadradic baseline
				n_line=np/3-1
				basen=2
			ENDELSE
		END
		1 : BEGIN		;constant baseline
			n_line=np/3
			basen=0
		END

		2 : BEGIN		;linear baseline
			n_line=np/3
			basen=1	
		END
	ENDCASE
			
	type=size(p,/type)
	IF type EQ 5 THEN L = dblarr(3,n_line) ELSE L=fltarr(3,n_line)

	FOR i = 0, n_line-1 DO FOR j = 0,2 DO L[j,i] = p[3*i+j]
	nx=n(x)+1
	dx=0.01		;hard coded for zjk
	x0=3.94		;hard coded for zjk
	y=fltarr(nx)
	FOR i = 0, n_line-1 DO y = y + L[0,i]*exp(-(x-L[1,i])^2/(2.*L[2,i]^2))
	CASE basen OF 
		3 : base=0
		2 : base=p[n(p)]+p[n(p)-1]*(x-x0)/dx+p[n(p)-2]*(x-x0)^2/dx^2	;
		1 : base=p[n(p)]+p[n(p)-1]*(x-x0)/dx				; 
		0 : base=p[n(p)]
		ELSE :	
	ENDCASE
	y+=base

	RETURN, y
END

;+
;NAME:
;	HIREXSR_IV2IS
;
;PURPOSE:
;	This function converts a velocity in [km/s] into a line shift in [Ang]
;
;MODIFICATION HISTORY:
;	Written by: 	M.L. Reinke - 9/5/11
;
;-

FUNCTION hirexsr_iv2is,iv,lam_o						;input instrumental velocity [km/s], get instrumental shift [Ang]
	c=3.0e8  			;speed of light
	vconv=lam_o/(c*1.0e-3)
	output=iv*vconv
	RETURN,output
END

;+
;NAME:
;	HIREXSR_IS2IV
;
;PURPOSE:
;	This function converts a line shift in [Ang] to a velocity in [km/s] 
;
;MODIFICATION HISTORY:
;	Written by: 	M.L. Reinke - 9/5/11
;
;-

FUNCTION hirexsr_is2iv,is,lam_o						;input instrumental shift [Ang], get instrumental velocity [km/s]
	c=3.0e8  			;speed of light
	vconv=lam_o/(c*1.0e-3)
	output=is/vconv
	RETURN,output
END

;+
;NAME:
;	HIREXSR_ITI2IW
;
;PURPOSE:
;	This function converts a temperature [keV] to a line width [Ang]
;
;MODIFICATION HISTORY:
;	Written by: 	M.L. Reinke - 9/4/11
;
;-


FUNCTION hirexsr_iti2iw,iti,z,lam_o					;input instrumental ti [keV], get instrumental width [Ang]
	c=3.0e8  			;speed of light
	e=1.60e-19			;conversion for eV -> J
	mconv=1.66e-27			;conversion for amu -> kg
	mass=read_atomic_mass(z)
	ticonv=(lam_o/c)^2*(e*1.0e3/(mass*mconv))
	output=sqrt(iti*ticonv)
	RETURN,output
END

;+
;NAME:
;	HIREXSR_IW2ITI
;
;PURPOSE:
;	This function converts a line width [Ang] to temperature [keV]
;
;MODIFICATION HISTORY:
;	Written by: 	M.L. Reinke - 9/4/11
;
;-

FUNCTION hirexsr_iw2iti,iw,z,lam_o					;input instrumental width [Ang], get instrumental ti [keV]
	c=3.0e8  			;speed of light
	e=1.60e-19			;conversion for eV -> J
	mconv=1.66e-27			;conversion for amu -> kg
	mass=read_atomic_mass(z)
	ticonv=(lam_o/c)^2*(e*1.0e3/(mass*mconv))
	output=iw^2/ticonv
	RETURN,output
END

;+
;NAME:
;	HIREXSR_FIT_WN3
;
;PURPOSE:
;	This function performs a multi-gaussian fit on the He-like Ar w n>=3 spectrum.  Currently
;	this is done using a 3 gaussian fit for the w and n=3 and n=4 satellite lines.
;
;CALLING SEQUENCE:
;	result=HIREXSR_FIT_WN3(spec,lam)
;
;INPUTS:
;	spec	FLTARR [n] of the spectral brightness at each wavelength
;	lam	FLTARR [n] of the wavelengths [Angstroms]
;
;OPTIONAL INPUTS:
;	sig	FLTARR [n] of the uncertainty in the spec data DEFAULT: sqrt(spec)
;	limit	FLOAT	of a lower limit on the signal or signal to noise ratio required to proceed with fitting DEFAULT: 0
;	win	INT	of the window number to plot the estimate and fit+residual DEFAULT: (no plot)
;
;KEYWORD PARAMETERS:
;	/verb sets verbose mode and outputs messages to the terminal
;
;OUTPUTS:
;	result	FLTARR	[3*nlines+1] of the gaussian fit coefficients a*exp(-(x-b)^2/(2*c^2)) along with the baseline as the last element
;
;OPTIONAL OUTPUTS:
;	label	STRARR 	[nlines] of the line labels for the corresponding fit coefficients
;	
;
;PROCEDURE:
;	This function calls GAUSSIAN_FITS using MPFITFUN to perfom the fit.
;	Rest wavelengths are called using HIREXSR_LOAD_WAVELENGTHS
;	The following restrictions are used to constrain the fit
;		- all intensities are >= 0
;		- widths of all lines are restricted to 0.5 -> 1.5 the seeded width of the w-line
;
;	Two checks are made to ensure proper signal and seeding before MPFITFUN is called.
;		1a) if limit is > 0 then total(spec-mean(spec)) > limit
;		1b) if limit is < 0 and sig is given then max(spec/sig) < |limit|
;		2) the seed lambda for the w-line must be such that velocity < 350 km/s 
;
;	If 1a or 1b is not satisfied result=-2 is return and if 2 is not satisfied then result=-3 is returned.
;
;MODIFICATION HISTORY;
;	Written by	M.L. Reinke (adapted from HIREXSR_FIT_HE) 3/10
;	8/10		M.L. Reinke - added the SNR check before fitting
;	8/27/10		M.L. Reinke - modifided the SNR limit check to look only around the w-line and do a basline subtraction
;	8/30/10		M.L. Reinke - added the n_lines optional output so that HIREXSR_FIT_SPECTRA can form a valid baseline only coefs vector
;	9/14/10		M.L. Reinke - added a RETURN=-4 if trying to seed off the edge of the spectral region
;	12/17/10	M.L. Reinke - added the n=5 line to the fit, and changed the n=4,n=3 constraints on the fit.  Adjusted the w1 seed -3.0e-4
;-

FUNCTION hirexsr_fit_wn3,spec,lam,sig=sig,win=win,double=double,label=label,limit=limit,verb=verb,n_lines=n_lines

	n_lines = 4
	IF keyword_set(double) THEN L = dblarr(3,n_lines) ELSE L=fltarr(3,n_lines)

	;basic check if spectra is worth fitting	
	IF NOT keyword_set(limit) THEN limit=0.0
	IF limit GT 0 THEN BEGIN
		meanoff=total(spec-min(spec))
		IF meanoff LE limit THEN BEGIN
			IF keyword_set(verb) THEN print, 'Too Few Photons: COEFS=-2'
			RETURN,-2		;fail w/ too few photons
		ENDIF
	ENDIF
	ilow=ipt(lam,3.946)
	iup=ipt(lam,3.954)
	IF limit LT 0 AND keyword_set(sig) THEN BEGIN
		snr_limit=abs(limit)
		snr_max=max((spec[ilow:iup]-min(spec))/sig[ilow:iup])
		IF snr_max LE snr_limit THEN BEGIN
			IF keyword_set(verb) THEN print, 'SNR too low: COEFS=-2'
			RETURN,-2		;fail w/ too low of signal to noise ratio
		ENDIF
	ENDIF

	;make w-line a sub-array
	w_line = spec[ilow:iup]
	w_lam=lam[ilow:iup]

	;fit the w-line to get the average line width
	dummy = gaussfit(w_lam,w_line,coefs, nterm=4)
 	W0 = coefs[0]                   ; max(dummy)
	W1 = coefs[1]-3.0e-4            ; where(dummy Eq max(dummy))
	W2 = coefs[2]
	width = W2

	lam_o_shot=-1		;change to model eventually
	hirexsr_load_wavelengths,lam_o_shot,lam_o,z,label
	L_n5=last(lam_o[where(z EQ 18 AND label EQ 'wn5')])
	L_n4=last(lam_o[where(z EQ 18 AND label EQ 'wn4')])
	L_n3=last(lam_o[where(z EQ 18 AND label EQ 'wn3')])
	L_w=last(lam_o[where(z EQ 18 AND label EQ 'w')])


	shift=L_w-w1	;estimate of shift for seeding
	IF abs(shift) GT 0.005 THEN RETURN, -3	;fail with bad seed (.005 > 350 km/s)
	
	;set the initial guess
	;----------------------------

	i = 0
	cent=L_w-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4
	L[0,i] = spec[ipt(lam,cent)]	; w height
	L[1,i] = cent		   	; w pos
	L[2,i] = width                  ; w wid
	i++

	cent=L_n5-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4		;return if seed would push other lines off the spectral range
	L[0,i] = spec[ipt(lam,cent)]*0.0; n5 height
	L[1,i] = cent		   	; n5 pos
	L[2,i] = width                  ; n5 wid
	i++

	cent=L_n4-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4		;return if seed would push other lines off the spectral range
	L[0,i] = spec[ipt(lam,cent)]	; n4 height
	L[1,i] = cent		   	; n4 pos
	L[2,i] = width                  ; n4 wid
	i++

	cent=L_n3-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4		;return if seed would push other lines off the spectral range
	L[0,i] = spec[ipt(lam,cent)]	; n3 height
	L[1,i] = cent		   	; n3 pos
	L[2,i] = width                  ; n3 wid
	i++

	label=['w','wn5','wn4','wn3']
	
	base_line = min(spec)
	IF keyword_set(double) THEN estimate = dblarr(n_lines*3+1) ELSE estimate=fltarr(n_lines*3+1)
	FOR i = 0,n_lines-1 DO FOR j = 0,2 DO estimate[3*i+j] = L[j,i]
	estimate[n_lines*3] = base_line 
	
	;define MPFITFUN constraints
	IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(estimate)) ELSE $
	 	parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(estimate))

	;write starting values
	parinfo[*].value=estimate


	; Tie n4,n5 shift to n3 shift
	i = 1
	parinfo(3*i+1).tied = 'P('+strtrim(3*(i+2)+1,1)+')-'+strtrim((L_n3-L_n5),1)
	i = 2
	parinfo(3*i+1).tied = 'P('+strtrim(3*(i+1)+1,1)+')-'+strtrim((L_n3-L_n4),1)

	; tie n4 width to n3 width
	i=2
	parinfo(3*i+2).tied = 'P('+strtrim(3*(i+1)+2,1)+')' ;

	FOR i=0,n_lines-1 DO parinfo[3*i].limited[0]=1		;states that all intensity coefs will have a lower bound
	FOR i=0,n_lines-1 DO parinfo[3*i].limits[0]=0.0      	;sets the lower bound to be 0.0

	;width
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limited[0] = 1  ; states that width parameters will have an lower bound
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limited[1] = 1  ; states that width parameters will have an upper bound
	FOR i = 0,n_lines-1 DO parinfo(3*i+2).limits[0] = 0.5*width ; sets the lower bound to be 50% of w line width 
	FOR i = 0,n_lines-1 DO parinfo(3*i+2).limits[1] = 1.5*width ; sets the upper bound to be 150% of w line width

	;perfom fit
	IF NOT keyword_set(sig) THEN sig=sqrt(spec)
	coefs = mpfitfun('gaussian_fits', lam,spec,sig+1, estimate,/quiet, parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)

	IF keyword_set(win) THEN BEGIN
		resid=spec-gaussian_fits(lam,coefs)
		spec_est=gaussian_fits(lam,estimate)
		spec_fit=gaussian_fits(lam,coefs)
		makesym,10
		IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
		d_old=!d
		ls=1.0
		IF keyword_set(ps) THEN BEGIN
			device, xsize=7.0,ysize=7.0*700.0/1200.0,/inches
			ls=0.6
		ENDIF ELSE openwin,win+1,xsize=1200,ysize=700
		pos=[0.1,0.35,0.95,0.95]
 		plot, lam,spec, psym=8,title='wn3 Spectral Fit', xtitle = 'Pixel #', ytitle = 'Counts',pos=pos,chars=1.5*ls,/xsty,symsize=1.0*ls
	    	oplot, lam, spec_fit, color = 100
    		FOR i = 0, n_lines-1 DO oplot, lam, coefs[3*i]*exp(-(lam-coefs[3*i+1])^2/(2.*coefs[3*i+2]^2)), color = 200, line =2
		pos=[0.1,0.05,0.95,0.25]
		plot,lam,resid,psym=8,ytit=n2g('Delta')+'Counts',pos=pos,/noerase,/nodata,chars=1.5*ls,/xsty,symsize=1.0*ls
		oplot, [0,max(lam)],[0,0],linestyle=2.0
		oplot,lam,resid,psym=8,color=200,symsize=1.0*ls
		IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm,/port	
	ENDIF

	RETURN,coefs
END


;+
;NAME:
;	HIREXSR_FIT_XY
;
;PURPOSE:
;	This function performs a multi-gaussian fit on the He-like Ar x+y spectrum.  Currently
;	this is done using a 6 gaussian fit for the n,x,st,y and two n=3/4 satellite groups.
;
;CALLING SEQUENCE:
;	result=HIREXSR_FIT_XY(spec,lam)
;
;INPUTS:
;	spec	FLTARR [n] of the spectral brightness at each wavelength
;	lam	FLTARR [n] of the wavelengths [Angstroms]
;
;OPTIONAL INPUTS:
;	sig	FLTARR [n] of the uncertainty in the spec data DEFAULT: sqrt(spec)
;	limit	FLOAT	of a lower limit on the signal or signal to noise ratio required to proceed with fitting DEFAULT: 0
;	win	INT	of the window number to plot the estimate and fit+residual DEFAULT: (no plot)
;
;KEYWORD PARAMETERS:
;	/verb sets verbose mode and outputs messages to the terminal
;
;OUTPUTS:
;	result	FLTARR	[3*nlines+1] of the gaussian fit coefficients a*exp(-(x-b)^2/(2*c^2)) along with the baseline as the last element
;
;OPTIONAL OUTPUTS:
;	label	STRARR 	[nlines] of the line labels for the corresponding fit coefficients
;	
;
;PROCEDURE:
;	This function calls GAUSSIAN_FITS using MPFITFUN to perfom the fit.
;	Rest wavelengths are called using HIREXSR_LOAD_WAVELENGTHS
;	The following restrictions are used to constrain the fit
;		- all intensities are >= 0
;		- widths of all lines are restricted to 0.5 -> 1.5 the seeded width of the y-line
;		- shift and width of the st line is fixed to that of the y-line
;
;	Two checks are made to ensure proper signal and seeding before MPFITFUN is called.
;		1a) if limit is > 0 then total(spec-mean(spec)) > limit
;		1b) if limit is < 0 and sig is given then max(spec/sig) < |limit|
;		2) the seed lambda for the w-line must be such that velocity < 350 km/s 
;
;	If 1a or 1b is not satisfied result=-2 is return and if 2 is not satisfied then result=-3 is returned.
;
;MODIFICATION HISTORY;
;	Written by	M.L. Reinke (adapted from HIREXSR_FIT_HE) 3/10
;	8/10		M.L. Reinke - added the SNR check before fitting
;	8/30/10		M.L. Reinke - added the n_lines optional output so that HIREXSR_FIT_SPECTRA can form a valid baseline only coefs vector
;	9/14/10		M.L. Reinke - added a RETURN=-4 if trying to seed off the edge of the spectral region
;
;-


FUNCTION hirexsr_fit_xy,spec,lam,sig=sig,win=win,double=double,label=label,limit=limit,verb=verb,n_lines=n_lines

	n_lines = 6
	IF keyword_set(double) THEN L = dblarr(3,n_lines) ELSE L=fltarr(3,n_lines)

	;basic check if spectra is worth fitting
	IF NOT keyword_set(limit) THEN limit=0.0
	IF limit GT 0 THEN BEGIN
		meanoff=total(spec-min(spec))
		IF meanoff LE limit THEN BEGIN
			IF keyword_set(verb) THEN print, 'Too Few Photons: COEFS=-2'
			RETURN,-2		;fail w/ too few photons
		ENDIF
	ENDIF
	ilow=ipt(lam,3.9675)
	iup=ipt(lam,3.973)
	IF limit LT 0 AND keyword_set(sig) THEN BEGIN
		snr_limit=abs(limit)
		snr_max=max((spec[ilow:iup]-min(spec))/sig[ilow:iup])
		IF snr_max LE snr_limit THEN BEGIN
			IF keyword_set(verb) THEN print, 'SNR too low: COEFS=-2'
			RETURN,-2		;fail w/ too low of signal to noise ratio
		ENDIF
	ENDIF

	;make xy-line a sub-array
	y_line = spec[ilow:iup]
	y_lam=lam[ilow:iup]

	;fit the y-line to get the average line width
	dummy = gaussfit(y_lam,y_line,coefs, nterm=4)
 	W0 = coefs[0]                   ; max(dummy)
	W1 = coefs[1]                   ; where(dummy Eq max(dummy))
	W2 = coefs[2]
	width = W2

	lam_o_shot=-1		;change to model eventually
	hirexsr_load_wavelengths,lam_o_shot,lam_o,z,label_o
	L_n=last(lam_o[where(z EQ 18 AND label_o EQ 'n')])
	L_x=last(lam_o[where(z EQ 18 AND label_o EQ 'x')])
	L_y=last(lam_o[where(z EQ 18 AND label_o EQ 'y')])
	L_st=last(lam_o[where(z EQ 18 AND label_o EQ 'st')])
	L_yn4=last(lam_o[where(z EQ 18 AND label_o EQ 'yn4')])
	L_yn3=last(lam_o[where(z EQ 18 AND label_o  EQ 'yn3')])

	shift=L_y-w1	;estimate of shift for seeding
	IF abs(shift) GT 0.005 THEN RETURN, -3	;fail with bad seed (.005 > 350 km/s)
	
	;set the initial guess
	;----------------------------

	i = 0
	cent=L_n-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4		;return if seed would push other lines off the spectral range
	L[0,i] = spec[ipt(lam,cent)]	; n height
	L[1,i] = cent		   	; n pos
	L[2,i] = width                  ; n wid
	i++

	cent=L_x-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4		
	L[0,i] = spec[ipt(lam,cent)]	; x height
	L[1,i] = cent		   	; x pos
	L[2,i] = width                  ; x wid
	i++

	cent=L_y-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4		
	L[0,i] = spec[ipt(lam,cent)]	; y height
	L[1,i] = cent		   	; y pos
	L[2,i] = width                  ; y wid
	i++

	cent=L_st-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4		
	L[0,i] = spec[ipt(lam,cent)]*0.0; st height
	L[1,i] = cent		   	; st pos
	L[2,i] = width                  ; st wid
	i++

	cent=L_yn4-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4		
	L[0,i] = spec[ipt(lam,cent)]	; yn4 height
	L[1,i] = cent		   	; yn4 pos
	L[2,i] = width                  ; yn4 wid
	i++

	cent=L_yn3-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4		
	L[0,i] = spec[ipt(lam,cent)]	; yn3 height
	L[1,i] = cent		   	; yn3 pos
	L[2,i] = width                  ; yn3 wid
	i++

	label=['n','x','y','st','yn4','yn3']
	
	base_line = min(spec)
	IF keyword_set(double) THEN estimate = dblarr(n_lines*3+1) ELSE estimate=fltarr(n_lines*3+1)
	FOR i = 0,n_lines-1 DO FOR j = 0,2 DO estimate[3*i+j] = L[j,i]
	estimate[n_lines*3] = base_line 
	
	;define MPFITFUN constraints
	IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(estimate)) ELSE $
	 	parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(estimate))

	;write starting values
	parinfo[*].value=estimate

	i=3
	parinfo(3*i+1).tied = 'P('+strtrim(3*(i-1)+1,1)+')-'+strtrim((L_y-L_st),1) ;tie y and st shift
	parinfo(3*i+2).tied = 'P('+strtrim(3*(i-1)+2,1)+')' 	;tie y and st width

	;intensity
	FOR i=0,n_lines-1 DO parinfo[3*i].limited[0]=1		;states that all intensity coefs will have a lower bound
	FOR i=0,n_lines-1 DO parinfo[3*i].limits[0]=0.0      	;sets the lower bound to be 0.0

	;width
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limited[0]=1  	; states that width parameters will have an lower bound
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limited[1]= 1  	; states that width parameters will have an upper bound
	FOR i=0,n_lines-1 DO parinfo(3*i+2).limits[0]=0.5*width ;sets the lower bound to be 50% of y line width 
	FOR i=0,n_lines-1 DO parinfo(3*i+2).limits[1]=1.5*width ; sets the upper bound to be 150% of y line width

	;perfom fit
	IF NOT keyword_set(sig) THEN sig=sqrt(spec)
	coefs = mpfitfun('gaussian_fits', lam,spec,sig+1, estimate, parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)

	IF keyword_set(win) THEN BEGIN
		resid=spec-gaussian_fits(lam,coefs)
		spec_est=gaussian_fits(lam,estimate)
		spec_fit=gaussian_fits(lam,coefs)
		makesym,10
		IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
		d_old=!d
		ls=1.0
		IF keyword_set(ps) THEN BEGIN
			device, xsize=7.0,ysize=7.0*700.0/1200.0,/inches
			ls=0.6
		ENDIF ELSE openwin,win+1,xsize=1200,ysize=700
		pos=[0.1,0.35,0.95,0.95]
 		plot, lam,spec, psym=8,title='xy Spectral Fit', xtitle = 'Pixel #', ytitle = 'Counts',pos=pos,chars=1.5*ls,/xsty,symsize=1.0*ls
	    	oplot, lam, spec_fit, color = 100
    		FOR i = 0, n_lines-1 DO oplot, lam, coefs[3*i]*exp(-(lam-coefs[3*i+1])^2/(2.*coefs[3*i+2]^2)), color = 200, line =2
		pos=[0.1,0.05,0.95,0.25]
		plot,lam,resid,psym=8,ytit=n2g('Delta')+'Counts',pos=pos,/noerase,/nodata,chars=1.5*ls,/xsty,symsize=1.0*ls
		oplot, [0,max(lam)],[0,0],linestyle=2.0
		oplot,lam,resid,psym=8,color=200,symsize=1.0*ls
		IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm,/port	
	ENDIF
	RETURN,coefs
END

;+
;NAME:
;	HIREXSR_FIT_ZJK
;
;PURPOSE:
;	This function performs a multi-gaussian fit on the He-like Ar qra+zjk spectrum.  Currently
;	this is done using a 5 gaussian fit for the q,r,a,k and z lines.  The j line is degenerate w/ z
;	to within the HIREXSR resolution.
;
;CALLING SEQUENCE:
;	result=HIREXSR_FIT_ZJK(spec,lam)
;
;INPUTS:
;	spec	FLTARR [n] of the spectral brightness at each wavelength
;	lam	FLTARR [n] of the wavelengths [Angstroms]
;
;OPTIONAL INPUTS:
;	sig	FLTARR [n] of the uncertainty in the spec data DEFAULT: sqrt(spec)
;	limit	FLOAT	of a lower limit on the signal or signal to noise ratio required to proceed with fitting DEFAULT: 0
;	win	INT	of the window number to plot the estimate and fit+residual DEFAULT: (no plot)
;	nback	INT	of the polynomial order to use for the background
;
;KEYWORD PARAMETERS:
;	/verb sets verbose mode and outputs messages to the terminal
;
;OUTPUTS:
;	result	FLTARR	[3*nlines+1] of the gaussian fit coefficients a*exp(-(x-b)^2/(2*c^2)) along with the baseline as the last element
;
;OPTIONAL OUTPUTS:
;	label	STRARR 	[nlines] of the line labels for the corresponding fit coefficients
;	
;
;PROCEDURE:
;	This function calls GAUSSIAN_FITS using MPFITFUN to perfom the fit.
;	Rest wavelengths are called using HIREXSR_LOAD_WAVELENGTHS
;	The following restrictions are used to constrain the fit
;		- all intensities are >= 0
;		- width of the r and a line is fixed to that of the q-line
;
;	Two checks are made to ensure proper signal and seeding before MPFITFUN is called.
;		1a) if limit is > 0 then total(spec-mean(spec)) > limit
;		1b) if limit is < 0 and sig is given then max(spec/sig) < |limit|
;		2) the seed lambda for the w-line must be such that velocity < 350 km/s 
;
;	If 1a or 1b is not satisfied result=-2 is return and if 2 is not satisfied then result=-3 is returned.
;
;MODIFICATION HISTORY;
;	Written by	M.L. Reinke (adapted from HIREXSR_FIT_HE) 3/10
;	8/10		M.L. Reinke - added the SNR check before fitting
;	8/27/10		M.L. Reinke - modifided the SNR limit check to look only around the z-line and do a basline subtraction
;	8/30/10		M.L. Reinke - added the n_lines optional output so that HIREXSR_FIT_SPECTRA can form a valid baseline only coefs vector
;	9/14/10		M.L. Reinke - added a RETURN=-4 if trying to seed off the edge of the spectral region
;	5/13/11		M.L. Reinke - adjusted fitting constraints for qra and kj to prevent them from wandering around 
;       3/13/14         C. Gao - added the sanity check of ilow and iup. If ilow=-1, then ilow=0; if iup=-1, then iup=n_elements(lam)-1
;-


FUNCTION hirexsr_fit_zjk,spec,lam,sig=sig,win=win,double=double,label=label,limit=limit,verb=verb,n_lines=n_lines,nback=nback

	n_lines = 6
	IF keyword_set(double) THEN L = dblarr(3,n_lines) ELSE L=fltarr(3,n_lines)

	;basic check if spectra is worth fitting
	IF NOT keyword_set(limit) THEN limit=0.0
	IF limit GT 0 THEN BEGIN
		meanoff=total(spec-min(spec))
		IF meanoff LE limit THEN BEGIN
			IF keyword_set(verb) THEN print, 'Too Few Photons: COEFS=-2'
			RETURN,-2		;fail w/ too few photons
		ENDIF
	ENDIF
	ilow=ipt(lam,3.9915)
	iup=ipt(lam,3.9975)
        if (ilow eq -1) then ilow = 0
        if (iup eq -1) then iup = n_elements(lam)-1

	IF limit LT 0 AND keyword_set(sig) THEN BEGIN
		snr_limit=abs(limit)
		snr_max=max((spec[ilow:iup]-min(spec))/sig[ilow:iup])
		IF snr_max LE snr_limit THEN BEGIN
			IF keyword_set(verb) THEN print, 'SNR too low: COEFS=-2'
			RETURN,-2		;fail w/ too low of signal to noise ratio
		ENDIF
	ENDIF ELSE snr_maxloc=-1


	;make z-line a sub-array
	z_line = spec[ilow:iup]
	z_lam=lam[ilow:iup]

	;fit the z-line to get the average line width
	dummy = gaussfit(z_lam,z_line,coefs, nterm=4)
 	W0 = coefs[0]                   ; max(dummy)
	W1 = coefs[1]                   ; where(dummy Eq max(dummy))
	W2 = coefs[2]
	width = W2

	lam_o_shot=-1		;change to model eventually
	hirexsr_load_wavelengths,lam_o_shot,lam_o,z,label_o
	L_z=last(lam_o[where(z EQ 18 AND label_o EQ 'z')])
	L_q=last(lam_o[where(z EQ 18 AND label_o EQ 'q')])
	L_r=last(lam_o[where(z EQ 18 AND label_o EQ 'r')])
	L_a=last(lam_o[where(z EQ 18 AND label_o EQ 'a')])
	L_k=last(lam_o[where(z EQ 18 AND label_o EQ 'k')])
	L_j=last(lam_o[where(z EQ 18 AND label_o EQ 'j')])

	shift=L_z-w1	;estimate of shift for seeding
	IF abs(shift) GT 0.005 THEN RETURN, -3	;fail with bad seed (.005 > 350 km/s)
	
	;set the initial guess
	;----------------------------
	bl = min(spec)

	i = 0
	cent=L_q-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4		;return if seed would push other lines off the spectral range
	L[0,i] = spec[ipt(lam,cent)]-bl	; q height
	L[1,i] = cent		   	; q pos
	L[2,i] = width                  ; q wid
	i++

	cent=L_r-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4
	L[0,i] = spec[ipt(lam,cent)]-bl	; r height
	L[1,i] = cent		   	; r pos
	L[2,i] = width                  ; r wid
	i++

	cent=L_a-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4
	L[0,i] = spec[ipt(lam,cent)]-bl	; a height
	L[1,i] = cent		   	; a pos
	L[2,i] = width                  ; a wid
	i++

	cent=L_k-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4
	L[0,i] = spec[ipt(lam,cent)]-bl	; k height
	L[1,i] = cent		   	; k pos
	L[2,i] = width                  ; k wid
	i++

	cent=L_j-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4
	L[0,i] = L[0,i-1]*1.3576	; j height
	L[1,i] = cent		   	; j pos
	L[2,i] = width                  ; j wid
	i++

	cent=L_z-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4
	L[0,i] = spec[ipt(lam,cent)]-bl	; z height
	L[1,i] = cent		   	; z pos
	L[2,i] = width                  ; z wid
	i++

	label=['q','r','a','k','j','z']
	
	IF keyword_set(double) THEN estimate = dblarr(n_lines*3+1+nback) ELSE estimate=fltarr(n_lines*3+1+nback)
	FOR i = 0,n_lines-1 DO FOR j = 0,2 DO estimate[3*i+j] = L[j,i]
	estimate[n_lines*3+nback] = bl 	;if lin/quad background set to zero
	
	;define MPFITFUN constraints
	IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(estimate)) ELSE $
	 	parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(estimate))

	; Tie r & a  width to q
	i = 1
	parinfo(3*i+2).tied = 'P('+strtrim(3*(i-1)+2,1)+')' ;
	i = 2
	parinfo(3*i+2).tied = 'P('+strtrim(3*(i-2)+2,1)+')' ;

	;tie j and k together 
	i=4
	parinfo(3*i+0).tied = '1.3576*P('+strtrim(3*(i-1),1)+')'
	parinfo(3*i+1).tied = 'P('+strtrim(3*(i-1)+1,1)+')+'+strtrim((L_j-L_k),1) ;
	parinfo(3*i+2).tied = 'P('+strtrim(3*(i-1)+2,1)+')' ;

	FOR i=0,n_lines-1 DO parinfo[3*i].limited[0]=1		;states that all intensity coefs will have a lower bound
	FOR i=0,n_lines-1 DO parinfo[3*i].limits[0]=0.0      	;sets the lower bound to be 0.0

	;shift
	xline=[0,4]	;limit shift of q and k
	nx=n(xline)
	FOR i=0,nx DO parinfo[3*xline[i]+1].limited[0]=1  ; states that shift parameters will have an lower bound
	FOR i=0,nx DO parinfo[3*xline[i]+1].limited[1]=1  ; states that shift parameters will have an upper bound
	FOR i=0,nx DO parinfo[3*xline[i]+1].limits[0]=last(lam_o[where(z EQ 18 AND label_o EQ label[xline[i]])])-0.0005
	FOR i=0,nx DO parinfo[3*xline[i]+1].limits[1]=last(lam_o[where(z EQ 18 AND label_o EQ label[xline[i]])])+0.0005

	;width
	xline=[0,4]	;limit shift of q and k
	nx=n(xline)
	FOR i=0,nx DO parinfo[3*xline[i]+2].limited[0]=1  ; states that shift parameters will have an lower bound
	FOR i=0,nx DO parinfo[3*xline[i]+2].limited[1]=1  ; states that shift parameters will have an upper bound
	FOR i=0,nx DO parinfo[3*xline[i]+2].limits[0]=width*0.75
	FOR i=0,nx DO parinfo[3*xline[i]+2].limits[1]=width*1.25

	;perfom fit
	IF NOT keyword_set(sig) THEN sig=sqrt(spec)
	coefs = mpfitfun('gaussian_fits', lam,spec,sig+1, estimate,/quiet,parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)
	IF keyword_set(win) THEN BEGIN
		resid=spec-gaussian_fits(lam,coefs)
		spec_est=gaussian_fits(lam,estimate)
		spec_fit=gaussian_fits(lam,coefs)
		makesym,10
		IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
		d_old=!d
		ls=1.0
		IF keyword_set(ps) THEN BEGIN
			device, xsize=7.0,ysize=7.0*700.0/1200.0,/inches
			ls=0.6
		ENDIF ELSE openwin,win+1,xsize=1200,ysize=700
		pos=[0.1,0.35,0.95,0.95]
 		plot, lam,spec, psym=8,title='zjk Spectral Fit', xtitle = 'Pixel #', ytitle = 'Counts',pos=pos,chars=1.5*ls,/xsty,symsize=1.0*ls
		IF keyword_set(sig) THEN oploterror,lam,spec,sig,symsize=1.0*ls,psym=8
	    	oplot, lam, spec_fit, color = 100
    		FOR i = 0, n_lines-1 DO oplot, lam, coefs[3*i]*exp(-(lam-coefs[3*i+1])^2/(2.*coefs[3*i+2]^2)), color = 200, line =2
		pos=[0.1,0.05,0.95,0.25]
		plot,lam,resid,psym=8,ytit=n2g('Delta')+'Counts',pos=pos,/noerase,/nodata,chars=1.5*ls,/xsty,symsize=1.0*ls
		oplot, [0,max(lam)],[0,0],linestyle=2.0
		IF keyword_set(sig) THEN oploterror,lam,resid,sig,psym=8,color=200,symsize=1.0*ls,errcolor=200 ELSE oplot,lam,resid,psym=8,color=200,symsize=1.0*ls
		IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm,/port	
	ENDIF

	RETURN,coefs
END

;+
;NAME:
;	HIREXSR_FIT_HE
;
;PURPOSE:
;	This function is used to perform multigaussian fitting for the He-like Ar spectra in pixel space for use with W_HIREXSR_CALIB
;
;CALLING SEQUENCE:
;	result=HIREXSR_FIT_HE(spec)
;
;MODFICATION HISTORY
;	Written by	M.L. Reinke (adapted from HE_RAW by A. Ince-Cushman) - 2/10
;	9/9/10		M.L. Reinke - modified to have break on as a default and use a break vector
;	12/17/10	M.L. Reinke - added the j satellite to the fit, fully constrained by the k fit, added another high-n w satellite constrained by n=4,n=3
;
;-
;adapted from HE_RAW  Written by A. Ince-Cushman 12/07
FUNCTION hirexsr_fit_he,spec,win=win,double=double,break=break,verb=verb
	
	;setup 
	IF NOT keyword_set(break) THEN break=[10,55,105,194] 

	lam_o_shot=-1		;change to model eventually
	hirexsr_load_wavelengths,lam_o_shot,lam_o,z,label
	L_n5=lam_o[where(z EQ 18 AND label EQ 'wn5')]*1.0e3
	L_n4=lam_o[where(z EQ 18 AND label EQ 'wn4')]*1.0e3
	L_n3=lam_o[where(z EQ 18 AND label EQ 'wn3')]*1.0e3
	L_w=lam_o[where(z EQ 18 AND label EQ 'w')]*1.0e3
	L_x=lam_o[where(z EQ 18 AND label EQ 'x')]*1.0e3
	L_y=lam_o[where(z EQ 18 AND label EQ 'y')]*1.0e3
	L_st=lam_o[where(z EQ 18 AND label EQ 'st')]*1.0e3
	L_n=lam_o[where(z EQ 18 AND label EQ 'n')]*1.0e3
	L_yn4=lam_o[where(z EQ 18 AND label EQ 'yn4')]*1.0e3
	L_yn3=lam_o[where(z EQ 18 AND label EQ 'yn3')]*1.0e3
	L_q=lam_o[where(z EQ 18 AND label EQ 'q')]*1.0e3
	L_r=lam_o[where(z EQ 18 AND label EQ 'r')]*1.0e3
	L_a=lam_o[where(z EQ 18 AND label EQ 'a')]*1.0e3
	L_k=lam_o[where(z EQ 18 AND label EQ 'k')]*1.0e3
	L_j=lam_o[where(z EQ 18 AND label EQ 'j')]*1.0e3
	L_z=lam_o[where(z EQ 18 AND label EQ 'z')]*1.0e3
	
	;make w-line a sub-array
	n_lines = 16
	IF keyword_set(double) THEN L = dblarr(3,n_lines) ELSE l=fltarr(3,n_lines)

	;fit the w-line to get the average line width
	center = maxloc(spec[break[0]:break[1]])+break[0]	;force w to be between break[0] and break[1]
	range=7
	nspec=n(spec)
	lowpix=0 > (center-range)			;prevent from indexing off the array if w isn't the brightest
	highpix=nspec < (center+range)		
	w_line = spec[lowpix:highpix]	
	x_w=indgen(n(w_line)+1)
	dummy = gaussfit(x_w,w_line, coefs, nterm=4)
	W0 = coefs[0]                    
	P1 = coefs[1]+lowpix      	
	W2 = coefs[2]
	width = W2

	;estimate the dispersion by finding the z-line
	pz=maxloc(spec[break[2]:break[3]])+break[2]
   	lam2p = (L_z - L_w)/(pz-p1)
	pORlam = 1.0                	; mAngst per pixel
    	x = findgen(n_elements(spec))
    	W1 = P1
  	xrange = [ 0, n_elements(spec)]

	;set the initial guess
	;----------------------------

	i = 0
	L[0,i] = W0                     ; w height
	L[1,i] = W1                     ; w pos
	L[2,i] = width                  ; w width
	i++

	n_wn5 = (L_n5-L_w)/lam2p
	L[0,i] = spec[P1+n_wn5]*0.0	; n=5 Sats height
	L[1,i] = W1 + n_wn5*pORLam      ; n=5 Sats pos
	L[2,i] = width                  ; n=5 Sats wid
	i++

	n_wn4 = (L_n4-L_w)/lam2p
	L[0,i] = spec[P1+n_wn4]         ; n=4 Sats height
	L[1,i] = W1 + n_wn4*pORLam      ; n=4 Sats pos
	L[2,i] = width		        ; n=4 Sats wid
	i++

	n_wn3 = (L_n3-L_w)/lam2p
	L[0,i] = spec[P1+n_wn3]         ; n=3 Sats height
	L[1,i] = W1 + n_wn3*pORLam      ; n=3 Sats pos
	L[2,i] = width                  ; n=3 Sats wid
	i++

	n_n = (L_n-L_w)/lam2p
	L[0,i] = spec[P1+n_n]        	; n height
	L[1,i] = W1 + n_n*pORLam     	; n pos
	L[2,i] = width                  ; n wid
	i++

	n_wx =  (L_x-L_w)/lam2p 
	L[0,i] = spec[P1+n_wx]          ; x height
	L[1,i] = W1+n_wx*pORLam         ; x pos
	L[2,i] = width                  ; x wid
	i++

	n_wy = (L_y-L_w)/lam2p 
	L[0,i] = spec[P1+n_wy]          ; y height
	L[1,i] = W1 + n_wy*pORLam       ; y pos
	L[2,i] = width                  ; y wid
	i++

	n_st=(L_st-L_w)/lam2p
	L[0,i] = spec[P1+n_st]*0.0   	; st height
	L[1,i] = W1 + n_st*pORLam     	; st pos
	L[2,i] = width                  ; st wid
	i++

	n_yn4=(L_yn4-L_w)/lam2p
	L[0,i] = spec[P1+n_yn4]        	; yn4 height
	L[1,i] = W1 + n_yn4*pORLam     	; yn4 pos
	L[2,i] = width                  ; yn4 wid
	i++

	n_yn3=(L_yn3-L_w)/lam2p
	L[0,i] = spec[P1+n_yn3]        	; yn3 height
	L[1,i] = W1 + n_yn3*pORLam     	; yn3 pos
	L[2,i] = width                  ; yn3 wid
	i++

	n_wq = (L_q - l_w)/lam2p 
	L[0,i] = spec[P1+n_wq]          ; q height
	L[1,i] = W1+n_wq*pORLam         ; q pos
	L[2,i] = width                  ; q wid
	i++

	n_wr = (L_r - l_w)/lam2p 
	L[0,i] = spec[P1+n_wr]          ; r height
	L[1,i] = W1 + n_wr*pORLam       ; r pos
	L[2,i] = width                  ; r wid
	i++

	n_wa = (L_a - l_w)/lam2p 
	L[0,i] = spec[P1+n_wa]          ; a height
	L[1,i] = W1 + n_wa*pORLam       ; a pos
	L[2,i] =  width                 ; a wid
	i++

	n_wk = (L_k - l_w)/lam2p    
	L[0,i] = spec[P1+n_wk]          ; k height
	L[1,i] = W1 + n_wk*pORLam       ; k pos
	L[2,i] = width                  ; k wid
	i++

	n_wj = (L_j - l_w)/lam2p    
	L[0,i] = spec[P1+n_wk]*1.3	; j height 
	L[1,i] = W1 + n_wj*pORLam       ; j pos
	L[2,i] = width                  ; j wid
	i++

	n_wz = (L_z - l_w)/lam2p      
	L[0,i] = spec[P1+n_wz]          ; z height
	L[1,i] = W1 + n_wz*pORLam   	; z pos; 5 adjustment for 2009
	L[2,i] = width                  ; z wid
	i++
	
	label=['w','n5','n4','n3','n','x','y','st','yn4','yn3','q','r','a','k','j','z']
	
	base_line = min(spec)
	IF keyword_set(double) THEN estimate = dblarr(n_lines*3+1) ELSE estimate=fltarr(n_lines*3+1)
	FOR i = 0,n_lines-1 DO FOR j = 0,2 DO estimate[3*i+j] = L[j,i]
	estimate[n_lines*3] = base_line 
	spec_est = gaussian_fits(x, estimate)
	
	IF keyword_set(win) THEN BEGIN
    		openwin, win
    		plot, x, spec, psym = 4, xrange = xrange, thick = 2, title = 'Guesses for Ar'
    		oplot, x, spec_est, color = 100
    		p = estimate
    		FOR i = 0, n_lines-1 DO oplot, x, p[3*i]*exp(-(x-p[3*i+1])^2/(2.*p[3*i+2]^2))+p[n_lines*3], color = 200, line =2, thick =2
	ENDIF

	IF n(break) EQ 0 THEN BEGIN			;fit full spectrum
		;define MPFITFUN constraints
		IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(estimate)) ELSE $
		 	parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(estimate))
		parinfo(*).value = estimate

		; Tie n4 shift to n3 shift
		i = 1
		parinfo(3*i+1).tied = 'P('+strtrim(3*(i+1)+1,1)+')-'+strtrim((L_n3-L_n4)/lam2p*pORlam,1)

		; tie n4 width to n3 width
		parinfo(3*i+2).tied = 'P('+strtrim(3*(i+1)+2,1)+')' ;

		; Tie y shift & width to x shift & width
		i = 5
		parinfo(3*i+1).tied = 'P('+strtrim(3*(i-1)+1,1)+')+'+strtrim((L_y-L_x)/lam2p*pORlam,1)

		; Tie r width to x width
		i = 8
		parinfo(3*i+2).tied = 'P('+strtrim(3*(i-4)+2,1)+')' ;

		; Tie a & k  width to r
		i = 9
		parinfo(3*i+2).tied = 'P('+strtrim(3*(i-1)+2,1)+')' ;
		i = 10
		parinfo(3*i+2).tied = 'P('+strtrim(3*(i-2)+2,1)+')' ;


		parinfo[*].limited(0) = 1  ; states that all parameters will have an lower bound
		parinfo[*].limited(1) = 1  ; states that all parameters will have an upper bound

		; intens
		FOR i = 0, n_lines-1 DO parinfo[3*i+0].limits[0] = 0.0      		;sets the lower bound to be 0.0
		FOR i = 0, n_lines-1 DO parinfo[3*i+0].limits[1] = max(spec) 		;1.1 sets the upper bound to be less than max(spec) 

		;pos
		FOR i = 0, n_lines-1 DO parinfo(3*i+1).limits[0] = L[1,i]-3.0  ;sets the lower bound to be -3 pix of the estimate
		FOR i = 0, n_lines-1 DO parinfo(3*i+1).limits[1] = L[1,i]+3.0  ;sets the upper bound to be +3 pix of the estimate 

		;width
		FOR i = 0, n_lines-1 DO parinfo(3*i+2).limits[0] = 0.2*width ; sets the lower bound to be 80% of w line width
		FOR i = 0, n_lines-1 DO parinfo(3*i+2).limits[1] = 1.5*width ; sets the upper bound to be 150% of w line width

		; baseline
		parinfo[n_lines*3].limited(*) = 1
		parinfo[n_lines*3].limits[0] = 0 
		parinfo[n_lines*3].limits[1] = mean(spec)

		;perfom fit 
		p = mpfitfun('gaussian_fits', x,spec,spec^0.5+1, estimate,/quiet, parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)

		IF keyword_set(verb) THEN print,'MPFITFUN returned '+num2str(status,1)+' after '+num2str(niter,1)+' iterations.'
	ENDIF ELSE BEGIN
		;fit w+n>=3 portion
		;--------------------------------------
		start=0
		sub_lines=4 	;w, n=5,n=3, n=4
		sub_est=[estimate[start*3:(start+sub_lines)*3-1],estimate[n_lines*3]]	;take lines + baseline offset

		;define MPFITFUN constraints
		IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(sub_est)) ELSE $
		 	parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(sub_est))
		
		;write starting values
		parinfo(*).value = sub_est

		; Tie n4,n5 shift to n3 shift
		i = 1
		parinfo(3*i+1).tied = 'P('+strtrim(3*(i+2)+1,1)+')-'+strtrim((L_n3-L_n5)/lam2p*pORlam,1)
		i = 2
		parinfo(3*i+1).tied = 'P('+strtrim(3*(i+1)+1,1)+')-'+strtrim((L_n3-L_n4)/lam2p*pORlam,1)

		; tie n4 width to n3 width
		i=2
		parinfo(3*i+2).tied = 'P('+strtrim(3*(i+1)+2,1)+')' ;

		parinfo[*].limited(0) = 1  ; states that all parameters will have an lower bound
		parinfo[*].limited(1) = 1  ; states that all parameters will have an upper bound

		; intens
		FOR i = 0, sub_lines-1 DO parinfo[3*i+0].limits[0] = 0.0      		;sets the lower bound to be 0.0
		FOR i = 0, sub_lines-1 DO parinfo[3*i+0].limits[1] = max(spec)*1.2	;1.1 sets the upper bound to be less than max(spec) 

		;pos
		FOR i = 0, sub_lines-1 DO parinfo(3*i+1).limits[0] = L[1,i]-3  ;sets the lower bound to be -3 pix of the estimate
		FOR i = 0, sub_lines-1 DO parinfo(3*i+1).limits[1] = L[1,i]+3  ;sets the upper bound to be +3 pix of the estimate 

		;width
		FOR i = 0, sub_lines-1 DO parinfo(3*i+2).limits[0] = 0.5*width ; sets the lower bound to be 80% of w line width
		FOR i = 0, sub_lines-1 DO parinfo(3*i+2).limits[1] = 1.5*width ; sets the upper bound to be 150% of w line width

		; baseline
		parinfo[sub_lines*3].limited(*) = 1
		parinfo[sub_lines*3].limits[0] = 0 
		parinfo[sub_lines*3].limits[1] = mean(spec)

		;perfom fit
		sub_spec=spec[break[0]:break[1]]
		sub_x=x[break[0]:break[1]]
		p_wn = mpfitfun('gaussian_fits', sub_x,sub_spec,sub_spec^0.5+1, sub_est,/quiet, parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)
		IF keyword_set(verb) THEN print,'MPFITFUN returned '+num2str(status,1)+' after '+num2str(niter,1)+' iterations.
		;fit x+y portion
		;--------------------------------------
		start=4
		sub_lines=6	;n, x , y, st, yn3
		sub_est=[estimate[start*3:(start+sub_lines)*3-1],estimate[n_lines*3]]	;take lines + baseline offset
		sub_spec=spec[break[1]:break[2]]
		sub_x=x[break[1]:break[2]]


		;define MPFITFUN constraints
		IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(sub_est)) ELSE $
		 	parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(sub_est))

		;write starting values
		parinfo[*].value = sub_est
		
		parinfo[*].limited(0) = 1  ; states that all parameters will have an lower bound
		parinfo[*].limited(1) = 1  ; states that all parameters will have an upper bound

		; intens
		FOR i = 0, sub_lines-1 DO parinfo[3*i+0].limits[0] = 0.0      		;sets the lower bound to be 0.0
		FOR i = 0, sub_lines-1 DO parinfo[3*i+0].limits[1] = max(sub_spec) 	;1.1 sets the upper bound to be less than max(spec) 

		;pos
		FOR i = 0, sub_lines-1 DO parinfo(3*i+1).limits[0] = L[1,i+start]-3  ;sets the lower bound to be -3 pix of the estimate
		FOR i = 0, sub_lines-1 DO parinfo(3*i+1).limits[1] = L[1,i+start]+3  ;sets the upper bound to be +3 pix of the estimate 

		; Tie y shift & width to st shift & width
		i = 3
		parinfo(3*i+1).tied = 'P('+strtrim(3*(i-1)+1,1)+')-'+strtrim((L_y-L_st)/lam2p*pORlam,1)
		parinfo(3*i+2).tied = 'P('+strtrim(3*(i-1)+2,1)+')'

		;width
		FOR i = 0, sub_lines-1 DO parinfo(3*i+2).limits[0] = 0.8*width ; sets the lower bound to be 80% of w line width
		FOR i = 0, sub_lines-1 DO parinfo(3*i+2).limits[1] = 1.5*width ; sets the upper bound to be 150% of w line width

		; baseline
		parinfo[sub_lines*3].limited(*) = 1
		parinfo[sub_lines*3].limits[0] = 0 
		parinfo[sub_lines*3].limits[1] = mean(spec)

		;perfom fit
		p_xy = mpfitfun('gaussian_fits', sub_x,sub_spec,sub_spec^0.5+1, sub_est,parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)
		IF keyword_set(verb) THEN print,'MPFITFUN returned '+num2str(status,1)+' after '+num2str(niter,1)+' iterations.'

		;fit qra+zjk portion
		;--------------------------------------
		start=10
		sub_lines=6	;q,r,a,k,j,z
		sub_est=[estimate[start*3:(start+sub_lines)*3-1],estimate[n_lines*3]]	;take lines + baseline offset
		sub_spec=spec[break[2]:break[3]]
		sub_x=x[break[2]:break[3]]

		;define MPFITFUN constraints
		IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(sub_est)) ELSE $
		 	parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(sub_est))

		;write starting values
		parinfo[*].value = sub_est

		; Tie r & a  width to q
		i = 1
		parinfo(3*i+2).tied = 'P('+strtrim(3*(i-1)+2,1)+')' ;
		i = 2
		parinfo(3*i+2).tied = 'P('+strtrim(3*(i-2)+2,1)+')' ;

		i=4
		parinfo(3*i+0).tied = '1.3576*P('+strtrim(3*(i-1),1)+')'
		parinfo(3*i+1).tied = 'P('+strtrim(3*(i-1)+1,1)+')+'+strtrim((L_j-L_k)/lam2p*pORlam,1) ;
		parinfo(3*i+2).tied = 'P('+strtrim(3*(i-1)+2,1)+')' ;

		parinfo[*].limited(0) = 1  ; states that all parameters will have an lower bound
		parinfo[*].limited(1) = 1  ; states that all parameters will have an upper bound

		; intens
		FOR i = 0, sub_lines-1 DO parinfo[3*i+0].limits[0] = 0.0      		;sets the lower bound to be 0.0
		FOR i = 0, sub_lines-1 DO parinfo[3*i+0].limits[1] = max(sub_spec) 		;1.1 sets the upper bound to be less than max(spec) 

		;pos
		FOR i = 0, sub_lines-1 DO parinfo(3*i+1).limits[0] = L[1,i+start]-10  ;sets the lower bound to be -3 pix of the estimate
		FOR i = 0, sub_lines-1 DO parinfo(3*i+1).limits[1] = L[1,i+start]+10  ;sets the upper bound to be +3 pix of the estimate 

		;width
		FOR i = 0, sub_lines-1 DO parinfo(3*i+2).limits[0] = 0.8*width ; sets the lower bound to be 80% of w line width
		FOR i = 0, sub_lines-1 DO parinfo(3*i+2).limits[1] = 1.5*width ; sets the upper bound to be 150% of w line width

		; baseline
		parinfo[sub_lines*3].limited(*) = 1
		parinfo[sub_lines*3].limits[0] = 0 
		parinfo[sub_lines*3].limits[1] = mean(spec)
	
		;perfom fit
		p_zjk = mpfitfun('gaussian_fits', sub_x,sub_spec,sub_spec^0.5+1, sub_est,/quiet, parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)
		IF keyword_set(verb) THEN print, 'MPFITFUN returned '+num2str(status,1)+' after '+num2str(niter,1)+' iterations'
		p=[p_wn[0:n(p_wn)-1],p_xy[0:n(p_xy)-1],p_zjk]
	ENDELSE

	xfit=make(min(x),max(x),2000)
	spec_fit = gaussian_fits(xfit, p)
	resid=spec-gaussian_fits(x,p)

	IF keyword_set(win) THEN BEGIN
		IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
		d_old=!d
		ls=1.0
		IF keyword_set(ps) THEN BEGIN
			device, xsize=7.0,ysize=7.0*700.0/1200.0,/inches
			ls=0.6
		ENDIF ELSE openwin,win+1,xsize=1200,ysize=700
		pos=[0.1,0.35,0.95,0.95]
 		plot, x,spec, psym=8,title='He-like Ar Spectra', xtitle = 'Pixel #', ytitle = 'Counts',pos=pos,chars=1.5*ls,/xsty,symsize=1.0*ls
	    	oplot, xfit, spec_fit, color = 100
    		FOR i = 0, n_lines-1 DO oplot, xfit, p[3*i]*exp(-(xfit-p[3*i+1])^2/(2.*p[3*i+2]^2)), color = 200, line =2
		pos=[0.1,0.05,0.95,0.25]
		plot,x,resid,psym=8,ytit=n2g('Delta')+'Counts',pos=pos,/noerase,/nodata,chars=1.5*ls,/xsty,symsize=1.0*ls
		oplot, [0,max(x)],[0,0],linestyle=2.0
		oplot,x,resid,psym=8,color=200,symsize=1.0*ls
		IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm,/port	
	ENDIF
	output={spec:spec, residual:resid, xfit:xfit, yfit:spec_fit, c0:p[indgen(n_lines)*3+0], c1:p[indgen(n_lines)*3+1] , c2:p[indgen(n_lines)*3+2], off:last(p), lab:label}

	RETURN,output
END

;+
;NAME:
;	HIREXSR_FIT_LYA
;
;PURPOSE:
;	This function performs a multi-gaussian fit on the H-like Ar Lyman-a spectrum.  Currently
;	this is done using a 7 gaussian fit for the Lya1, Lya2, Mo4d and 4 satellite "groups".
;
;CALLING SEQUENCE:
;	result=HIREXSR_FIT_LYA(spec,lam)
;
;INPUTS:
;	spec	FLTARR [n] of the spectral brightness at each wavelength
;	lam	FLTARR [n] of the wavelengths [Angstroms]
;
;OPTIONAL INPUTS:
;	sig	FLTARR [n] of the uncertainty in the spec data DEFAULT: sqrt(spec)
;	limit	FLOAT	of a lower limit on the signal or signal to noise ratio required to proceed with fitting DEFAULT: 0
;	win	INT	of the window number to plot the estimate and fit+residual DEFAULT: (no plot)
;
;KEYWORD PARAMETERS:
;	/verb sets verbose mode and outputs messages to the terminal
;
;OUTPUTS:
;	result	FLTARR	[3*nlines+1] of the gaussian fit coefficients a*exp(-(x-b)^2/(2*c^2)) along with the baseline as the last element
;
;OPTIONAL OUTPUTS:
;	label	STRARR 	[nlines] of the line labels for the corresponding fit coefficients
;	
;
;PROCEDURE:
;	This function calls GAUSSIAN_FITS using MPFITFUN to perfom the fit.
;	Rest wavelengths are called using HIREXSR_LOAD_WAVELENGTHS
;	The following restrictions are used to constrain the fit
;		- all intensities are >= 0
;		- the #1, #2 and #3 satellites are tied to the lya2 line
;		- the #4 satellite is restricted to 1 mA shift from the seeded position (based on a shifted lya1)
;		- widths of lya1 and lya2 lines are restricted to 0.5 -> 1.25 the seeded width of the lya1
;		- widths of satellite lines are restricted to 0.5 -> 1.5 the seeded width of the lya1
;		- width of the mo4d line is restricted to 0.2 -> 0.9 of the seeded lya1
;		- if Mo 4d line is brighter than Ar lya1 then lya2 position is restricted
;
;	Two checks are made to ensure proper signal and seeding before MPFITFUN is called.
;		1a) if limit is > 0 then total(spec-mean(spec)) > limit
;		1b) if limit is < 0 and sig is given then max(spec/sig) < |limit|
;		2) the seed lambda for the lya1-line must be such that velocity < 350 km/s 
;
;	If 1a or 1b is not satisfied result=-2 is return and if 2 is not satisfied then result=-3 is returned.
;
;MODIFICATION HISTORY;
;	Written by	M.L. Reinke (adapted from HIREXSR_FIT_WN3) 8/18/10
;	8/30/10		M.L. Reinke - added the n_lines optional output so that HIREXSR_FIT_SPECTRA can form a valid baseline only coefs vector
;	6/14/11		M.L. Reinke - added a conditional statement that constrains Lya2 when Mo line is bright and reduced max Lya1,2 width to 1.25*seed
;					subtraced baseline from seeding heights, and reduced seed and limit of Mo line width.
;	9/2/11		M.L. Reinke - increased the height of lyas2 and tied lyas1 to the lya2 position
;	3/5/12		M.L. Reinke - modified the esitmate for the background to be the mean of the
;                                     first and last pixel
;-

FUNCTION hirexsr_fit_lya,spec,lam,sig=sig,win=win,double=double,label=label,limit=limit,verb=verb,n_lines=n_lines

	n_lines = 7
	IF keyword_set(double) THEN L = dblarr(3,n_lines) ELSE L=fltarr(3,n_lines)

	;basic check if spectra is worth fitting	
	IF NOT keyword_set(limit) THEN limit=0.0
	IF limit GT 0 THEN BEGIN
		meanoff=total(spec-min(spec))
		IF meanoff LE limit THEN BEGIN
			IF keyword_set(verb) THEN print, 'Too Few Photons: COEFS=-2'
			RETURN,-2		;fail w/ too few photons
		ENDIF
	ENDIF
	IF limit LT 0 AND keyword_set(sig) THEN BEGIN
		snr_limit=abs(limit)
		snr_max=max(spec/sig)
		IF snr_max LE snr_limit THEN BEGIN
			IF keyword_set(verb) THEN print, 'SNR too low: COEFS=-2'
			RETURN,-2		;fail w/ too low of signal to noise ratio
		ENDIF
	ENDIF

	;make lya1-line a sub-array
	ilow=ipt(lam,3.729)
	iup=ipt(lam,3.734)
	lya1_line = spec[ilow:iup]
	lya1_lam=lam[ilow:iup]

	;fit the lya1-line to get the average line width
	dummy = gaussfit(lya1_lam,lya1_line,coefs, nterm=4)
 	W0 = coefs[0]                   ; max(dummy)
	W1 = coefs[1]                   ; where(dummy Eq max(dummy))
	W2 = coefs[2]
	width = W2

	lam_o_shot=-1		;change to model eventually
	hirexsr_load_wavelengths,lam_o_shot,lam_o,z,label
	L_sat1=last(lam_o[where(z EQ 18 AND label EQ 'lyas1')])
   	L_lya1=last(lam_o[where(z EQ 18 AND label EQ 'lya1')])
   	L_sat2=last(lam_o[where(z EQ 18 AND label EQ 'lyas2')])
	L_lya2=last(lam_o[where(z EQ 18 AND label EQ 'lya2')])
	L_sat3=last(lam_o[where(z EQ 18 AND label EQ 'lyas3')])
	L_mo4d=last(lam_o[where(z EQ 42 AND label EQ '4d')])
   	L_sat4=last(lam_o[where(z EQ 18 AND label EQ 'lyas4')])

	shift=L_lya1-w1	;estimate of shift for seeding
	IF abs(shift) GT 0.005 THEN RETURN, -3	;fail with bad seed (.005 > 350 km/s)
	
	;set the initial guess
	;----------------------------

	i = 0
	cent=L_sat1-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4			;return if seed would push other lines off the spectral range
	L[0,i] = 0.0			; lya sat#1 height
	L[1,i] = cent		   	; lya sat#1  pos
	L[2,i] = width                  ; lya sat#1 wid
	i++

	cent=L_lya1-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4	
	L[0,i] = spec[ipt(lam,cent)]	; lya1 height
	L[1,i] = cent		   	; lya1  pos
	L[2,i] = width                  ; lya1 wid
	i++

	cent=L_sat2-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4	
	L[0,i] = 0.0			; lya sat#2 height
	L[1,i] = cent		   	; lya sat#2  pos
	L[2,i] = width                  ; lya sat#2 wid
	i++
	
	cent=L_lya2-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4	
	L[0,i] = spec[ipt(lam,cent)]	; lya2 height
	L[1,i] = cent		   	; lya2  pos
	L[2,i] = width                  ; lya2 wid
	i++

	cent=L_sat3-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4	
	L[0,i] = 0.0			; lya sat#3 height
	L[1,i] = cent		   	; lya sat#3  pos
	L[2,i] = width                  ; lya sat#3 wid
	i++

	cent=L_mo4d-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4	
	L[0,i] = spec[ipt(lam,cent)]	; mo4d height
	L[1,i] = cent		   	; mo4d  pos
	L[2,i] = width*0.66                  ; mo4d wid
	i++

	cent=L_sat4-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4	
	L[0,i] = 0.0			; lya sat#4 height
	L[1,i] = cent		   	; lya sat#4  pos
	L[2,i] = width                  ; lya sat#4 wid
	i++

	label=['lyas1','lya1','lyas2','lya2','lyas3','4d','lyas4']
	satlines=[0,2,4,6]

	base_line = 0.5*(spec[0]+last(spec))
	L[0,[1,3,5]]-=base_line
	IF keyword_set(double) THEN estimate = dblarr(n_lines*3+1) ELSE estimate=fltarr(n_lines*3+1)
	FOR i = 0,n_lines-1 DO FOR j = 0,2 DO estimate[3*i+j] = L[j,i]
	estimate[n_lines*3] = base_line 
	
	;define MPFITFUN constraints
	IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(estimate)) ELSE $
	 	parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(estimate))

	;write starting values
	parinfo[*].value=estimate
	
	;intens
	FOR i=0,n_lines-1 DO parinfo[3*i].limited[0]=1		;states that all intensity coefs will have a lower bound
	FOR i=0,n_lines-1 DO parinfo[3*i].limits[0]=0.0      	;sets the lower bound to be 0.0
	satlines=[0,2,4,6]
	FOR j = 0,3 DO BEGIN 						
		i=satlines[j]
		parinfo[3*i].limited[1]=1				;satellites will have an upper bound
		parinfo[3*i].limits[1]=0.05*(max(spec)-min(spec))	;satellites < 5% of the peak line brightness in spectral region adjusted for background
	ENDFOR
	parinfo[3*2].limits[1]=0.10*(max(spec)-min(spec))		;lya2 satellites < 10% of the peak line brightness

	;if Mo line is brighter than Ar lya1 then prevent lya2 from wandering
	IF L[0,1] LT L[0,5] THEN BEGIN
		i=3
		parinfo[3*i+1].limited[0]=1			;lya2 line will have lower bound on shift
		parinfo[3*i+1].limited[1]=1			;lya2 line will have an upper bound on shift						
		parinfo(3*i+1).limits[0] = L[1,i]-0.002		;prevent wander of lya2
		parinfo(3*i+1).limits[1] = L[1,i]+0.001			
	ENDIF

	; pos
	i=0
	parinfo[3*i+1].limited[0]=1			;satellite line will have lower bound on shift
	parinfo[3*i+1].limited[1]=1			;satellite line will have an upper bound on shift						
	parinfo(3*i+1).limits[0] = L[1,i]-0.001		;prevent wander of lya sat#1
	parinfo(3*i+1).limits[1] = L[1,i]+0.001

	i=6
	parinfo[3*i+1].limited[0]=1			;satellite line will have lower bound on shift
	parinfo[3*i+1].limited[1]=1			;satellite line will have an upper bound on shift
	parinfo(3*i+1).limits[0] = L[1,i]-0.001		;prevent wander of lya sat#4
	parinfo(3*i+1).limits[1] = L[1,i]+0.001

	i=0
	parinfo(3*i+1).tied = 'P('+strtrim(3*(i+1)+1,1)+')-'+strtrim((L_lya2-L_sat1),1) 	;tie lya sat#1 pos to lya2
	i=2
	parinfo(3*i+1).tied = 'P('+strtrim(3*(i+1)+1,1)+')-'+strtrim((L_lya2-L_sat2),1) 	;tie lya sat#2 pos to lya2
	i=4
	parinfo(3*i+1).tied = 'P('+strtrim(3*(i-1)+1,1)+')+'+strtrim((L_sat3-L_lya2),1) 	;tie lya sat#3 pos to lya2

	;width
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limited[0]=1	 ;all width coefs will have a lower bound	
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limited[1]=1	 ;all width coefs will have an upper bound
	i=1 
	parinfo(3*i+2).limits[0] = 0.5*width 		;set width for lya1 line	
	parinfo(3*i+2).limits[1] = 1.25*width 
	i=3
	parinfo(3*i+2).limits[0] = 0.5*width 		;set width for lya2 line
	parinfo(3*i+2).limits[1] = 1.25*width
	satlines=[0,2,4,6]
	FOR j = 0,3 DO BEGIN 				;set satellite width
		i=satlines[j]
		parinfo(3*i+2).limits[0] = 0.5*width 
		parinfo(3*i+2).limits[1] = 1.5*width
	ENDFOR

	i=5
	parinfo(3*i+2).limits[0] = 0.2*width		;set width for moly line
	parinfo(3*i+2).limits[1] = 0.9*width

	; baseline
	parinfo[n_lines*3].limited(*) = 1
	parinfo[n_lines*3].limits[0] = 0 
	parinfo[n_lines*3].limits[1] = mean(spec)

	;perfom fit
	IF NOT keyword_set(sig) THEN sig=sqrt(spec)
	coefs = mpfitfun('gaussian_fits', lam,spec,sig+1, estimate, parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)

	;stop
	IF keyword_set(win) THEN BEGIN
		resid=spec-gaussian_fits(lam,coefs)
		spec_est=gaussian_fits(lam,estimate)
		spec_fit=gaussian_fits(lam,coefs)
		makesym,10
		IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
		d_old=!d
		ls=1.0
		IF keyword_set(ps) THEN BEGIN
			device, xsize=7.0,ysize=7.0*700.0/1200.0,/inches
			ls=0.6
		ENDIF ELSE openwin,win+1,xsize=1200,ysize=700
		pos=[0.1,0.35,0.95,0.95]
 		plot, lam,spec, psym=8,title='Lya Spectral Fit', xtitle = 'Pixel #', ytitle = 'Counts',pos=pos,chars=1.5*ls,/xsty,symsize=1.0*ls
	    	oplot, lam, spec_fit, color = 100
    		FOR i = 0, n_lines-1 DO oplot, lam, coefs[3*i]*exp(-(lam-coefs[3*i+1])^2/(2.*coefs[3*i+2]^2)), color = 200, line =2
		pos=[0.1,0.05,0.95,0.25]
		plot,lam,resid,psym=8,ytit=n2g('Delta')+'Counts',pos=pos,/noerase,/nodata,chars=1.5*ls,/xsty,symsize=1.0*ls
		oplot, [0,max(lam)],[0,0],linestyle=2.0
		oplot,lam,resid,psym=8,color=200,symsize=1.0*ls
		IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm,/port
	ENDIF

	RETURN,coefs
END

;+
;NAME:
;	HIREXSR_FIT_TJ
;
;PURPOSE:
;	This function performs a multi-gaussian fit on the H-like Ar TJ satellite spectrum.  Currently
;	this is done using a 7 gaussian fit for the multiple satellite lines
;
;CALLING SEQUENCE:
;	result=HIREXSR_FIT_TJ(spec,lam)
;
;INPUTS:
;	spec	FLTARR [n] of the spectral brightness at each wavelength
;	lam	FLTARR [n] of the wavelengths [Angstroms]
;
;OPTIONAL INPUTS:
;	sig	FLTARR [n] of the uncertainty in the spec data DEFAULT: sqrt(spec)
;	limit	FLOAT	of a lower limit on the signal or signal to noise ratio required to proceed with fitting DEFAULT: 0
;	win	INT	of the window number to plot the estimate and fit+residual DEFAULT: (no plot)
;
;KEYWORD PARAMETERS:
;	/verb sets verbose mode and outputs messages to the terminal
;
;OUTPUTS:
;	result	FLTARR	[3*nlines+1] of the gaussian fit coefficients a*exp(-(x-b)^2/(2*c^2)) along with the baseline as the last element
;
;OPTIONAL OUTPUTS:
;	label	STRARR 	[nlines] of the line labels for the corresponding fit coefficients
;	
;
;PROCEDURE:
;	This function calls GAUSSIAN_FITS using MPFITFUN to perfom the fit.
;	Rest wavelengths are called using HIREXSR_LOAD_WAVELENGTHS
;	The following restrictions are used to constrain the fit
;		- all intensities are >= 0
;
;	Two checks are made to ensure proper signal and seeding before MPFITFUN is called.
;		1a) if limit is > 0 then total(spec-mean(spec)) > limit
;		1b) if limit is < 0 and sig is given then max(spec/sig) < |limit|
;		2) the seed lambda for the lya1-line must be such that velocity < 350 km/s 
;
;	If 1a or 1b is not satisfied result=-2 is return and if 2 is not satisfied then result=-3 is returned.
;
;MODIFICATION HISTORY;
;	Written by	M.L. Reinke (adapted from HIREXSR_FIT_LYA) 8/9/11
;
;-

FUNCTION hirexsr_fit_tj,spec,lam,sig=sig,win=win,double=double,label=label,limit=limit,verb=verb,n_lines=n_lines

	n_lines = 7
	IF keyword_set(double) THEN L = dblarr(3,n_lines) ELSE L=fltarr(3,n_lines)

	;basic check if spectra is worth fitting	
	IF NOT keyword_set(limit) THEN limit=0.0
	IF limit GT 0 THEN BEGIN
		meanoff=total(spec-min(spec))
		IF meanoff LE limit THEN BEGIN
			IF keyword_set(verb) THEN print, 'Too Few Photons: COEFS=-2'
			RETURN,-2		;fail w/ too few photons
		ENDIF
	ENDIF
	IF limit LT 0 AND keyword_set(sig) THEN BEGIN
		snr_limit=abs(limit)
		snr_max=max(spec/sig)
		IF snr_max LE snr_limit THEN BEGIN
			IF keyword_set(verb) THEN print, 'SNR too low: COEFS=-2'
			RETURN,-2		;fail w/ too low of signal to noise ratio
		ENDIF
	ENDIF

	;make J-line a sub-array
	ilow=ipt(lam,3.770)
	iup=ipt(lam,3.775)
	J_line = spec[ilow:iup]
	J_lam=lam[ilow:iup]

	;fit the lya1-line to get the average line width
	dummy = gaussfit(J_lam,J_line,coefs, nterm=4)
 	W0 = coefs[0]                   ; max(dummy)
	W1 = coefs[1]                   ; where(dummy Eq max(dummy))
	W2 = coefs[2]
	width = W2

	lam_o_shot=-1		;change to model eventually
	hirexsr_load_wavelengths,lam_o_shot,lam_o,z,label
	L_T=last(lam_o[where(z EQ 18 AND label EQ 'T')])
   	L_K=last(lam_o[where(z EQ 18 AND label EQ 'K')])
   	L_Q=last(lam_o[where(z EQ 18 AND label EQ 'Q')])
	L_B=last(lam_o[where(z EQ 18 AND label EQ 'B')])
	L_R=last(lam_o[where(z EQ 18 AND label EQ 'R')])
	L_A=last(lam_o[where(z EQ 18 AND label EQ 'A')])
   	L_J=last(lam_o[where(z EQ 18 AND label EQ 'J')])

	shift=L_J-w1	;estimate of shift for seeding
	IF abs(shift) GT 0.005 THEN RETURN, -3	;fail with bad seed (.005 > 350 km/s)
	
	;set the initial guess
	;----------------------------

	i = 0
	cent=L_T-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4			;return if seed would push other lines off the spectral range
	L[0,i] = spec[ipt(lam,cent)]			; T height
	L[1,i] = cent		   	; T pos
	L[2,i] = width                  ; T wid
	i++

	cent=L_K-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4	
	L[0,i] = 0.0			; K height
	L[1,i] = cent		   	; K  pos
	L[2,i] = width                  ; K wid
	i++

	cent=L_Q-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4	
	L[0,i] = 0.0			; Q height
	L[1,i] = cent		   	; Q  pos
	L[2,i] = width                  ; Q wid
	i++
	
	cent=L_B-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4	
	L[0,i] = 0.0			; B height
	L[1,i] = cent		   	; B  pos
	L[2,i] = width                  ; B wid
	i++

	cent=L_R-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4	
	L[0,i] = 0.0			; R height
	L[1,i] = cent		   	; R  pos
	L[2,i] = width                  ; R wid
	i++

	cent=L_A-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4	
	L[0,i] = 0.0			; A height
	L[1,i] = cent		   	; A  pos
	L[2,i] = width                  ; A wid
	i++

	cent=L_J-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4	
	L[0,i] = spec[ipt(lam,cent)]	; J height
	L[1,i] = cent		   	; J  pos
	L[2,i] = width                  ; J wid
	i++

	label=['T','K','Q','B','R','A','J']

	base_line = min(spec)
	L[0,[0,6]]-=base_line
	IF keyword_set(double) THEN estimate = dblarr(n_lines*3+1) ELSE estimate=fltarr(n_lines*3+1)
	FOR i = 0,n_lines-1 DO FOR j = 0,2 DO estimate[3*i+j] = L[j,i]
	estimate[n_lines*3] = base_line 
	
	;define MPFITFUN constraints
	IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(estimate)) ELSE $
	 	parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(estimate))

	;write starting values
	parinfo[*].value=estimate

	;intens
	FOR i=0,n_lines-1 DO parinfo[3*i].limited[0]=1		;states that all intensity coefs will have a lower bound
	FOR i=0,n_lines-1 DO parinfo[3*i].limits[0]=0.0      	;sets the lower bound to be 0.0

	; pos
	FOR i=0,n_lines-1 DO parinfo[3*i+1].limited[0]=1			;satellite line will have lower bound on shift
	FOR i=0,n_lines-1 DO parinfo[3*i+1].limited[1]=1			;satellite line will have an upper bound on shift						
	FOR i=0,n_lines-1 DO parinfo(3*i+1).limits[0] = L[1,i]-0.001		;prevent wander of lines 
	FOR i=0,n_lines-1 DO parinfo(3*i+1).limits[1] = L[1,i]+0.001

	;width
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limited[0]=1	 ;all width coefs will have a lower bound	
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limited[1]=1	 ;all width coefs will have an upper bound

	FOR i=0,n_lines-1 DO parinfo(3*i+2).limits[0] = 0.5*width 		
	FOR i=0,n_lines-1 DO parinfo(3*i+2).limits[1] = 1.25*width 

	;baseline
	parinfo[n_lines*3].limited(*) = 1
	parinfo[n_lines*3].limits[0] = 0 
	parinfo[n_lines*3].limits[1] = mean(spec)	

	;perfom fit
	IF NOT keyword_set(sig) THEN sig=sqrt(spec)
	coefs = mpfitfun('gaussian_fits', lam,spec,sig+1, estimate, parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)

	IF keyword_set(win) THEN BEGIN
		resid=spec-gaussian_fits(lam,coefs)
		spec_est=gaussian_fits(lam,estimate)
		spec_fit=gaussian_fits(lam,coefs)
		makesym,10
		IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
		d_old=!d
		ls=1.0
		IF keyword_set(ps) THEN BEGIN
			device, xsize=7.0,ysize=7.0*700.0/1200.0,/inches
			ls=0.6
		ENDIF ELSE openwin,win+1,xsize=1200,ysize=700
		pos=[0.1,0.35,0.95,0.95]
 		plot, lam,spec, psym=8,title='TJ Spectral Fit', xtitle = 'Pixel #', ytitle = 'Counts',pos=pos,chars=1.5*ls,/xsty,symsize=1.0*ls
	    	oplot, lam, spec_fit, color = 100
    		FOR i = 0, n_lines-1 DO oplot, lam, coefs[3*i]*exp(-(lam-coefs[3*i+1])^2/(2.*coefs[3*i+2]^2)), color = 200, line =2
		pos=[0.1,0.05,0.95,0.25]
		plot,lam,resid,psym=8,ytit=n2g('Delta')+'Counts',pos=pos,/noerase,/nodata,chars=1.5*ls,/xsty,symsize=1.0*ls
		oplot, [0,max(lam)],[0,0],linestyle=2.0
		oplot,lam,resid,psym=8,color=200,symsize=1.0*ls
		IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm,/port
	ENDIF

	RETURN,coefs
END


;+
;NAME:
;	HIREXSR_FIT_H
;
;-

;8/25/12 - increase RANGE from 4 to 6 to help fitting out of focus or
;          hot lines, also set satellites to < 10% of main peak height

;adapted from GET_H_LINES  Written byA. Ince-Cushman 12/07
FUNCTION hirexsr_fit_h, spec, win=win,double=double,break=break,debug=debug
	x=indgen(n(spec)+1)

	IF NOT keyword_set(break) THEN  break=[25,90,95,170] 

	;fit the lya1-line to get the average line width
	center = maxloc(spec[break[0]:break[1]])+break[0]	;force lya to be between break[0] and break[1]
	range=6
	nspec=n(spec)
	lowpix=0 > (center-range)			;prevent from indexing off the array if LYA1 isn't the brightest
	highpix=nspec < (center+range)		
	lya1_line = spec[lowpix:highpix]	
	x_lya1=indgen(n(lya1_line)+1)
	dummy = gaussfit(x_lya1,lya1_line, coefs, nterm=4)
	W0 = coefs[0]                    
	W1 = coefs[1]+lowpix      	
	W2 = coefs[2]
	width = W2
	
	IF w1 GT last(x) OR w1 LT x[0] THEN BEGIN	;indicate a crazy fit and set to seom defaults
		w1=center
		width=2.0
		w0=spec[center]/(width*sqrt(2.0*!pi))
	ENDIF
	
	n_lines=14
	
	;load rest wavelengths
	lam_o_shot=-1		;change to model eventually
	hirexsr_load_wavelengths,lam_o_shot,lam_o,z,label
	L_sat1=lam_o[where(z EQ 18 AND label EQ 'lyas1')]*1.0e3
   	L_lya1=lam_o[where(z EQ 18 AND label EQ 'lya1')]*1.0e3
   	L_sat2=lam_o[where(z EQ 18 AND label EQ 'lyas2')]*1.0e3
	L_lya2=lam_o[where(z EQ 18 AND label EQ 'lya2')]*1.0e3
	L_sat3=lam_o[where(z EQ 18 AND label EQ 'lyas3')]*1.0e3
	L_mo4d=lam_o[where(z EQ 42 AND label EQ '4d')]*1.0e3
   	L_sat4=lam_o[where(z EQ 18 AND label EQ 'lyas4')]*1.0e3
  	L_T=lam_o[where(z EQ 18 AND label EQ 'T')]*1.0e3
	L_K=lam_o[where(z EQ 18 AND label EQ 'K')]*1.0e3
	L_Q=lam_o[where(z EQ 18 AND label EQ 'Q')]*1.0e3
	L_B=lam_o[where(z EQ 18 AND label EQ 'B')]*1.0e3
	L_R=lam_o[where(z EQ 18 AND label EQ 'R')]*1.0e3
	L_A=lam_o[where(z EQ 18 AND label EQ 'A')]*1.0e3
  	L_J=lam_o[where(z EQ 18 AND label EQ 'J')]*1.0e3

	IF keyword_set(double) THEN L = dblarr(3,n_lines) ELSE L=fltarr(3,n_lines)
	pJ=maxloc(spec[break[2]:break[3]])+break[2]		;assume brightest line in 2nd spectra is J
	IF pJ-w1 LT 70 THEN dpix=104.0 ELSE dpix=pJ-w1			;set to nominal value if pJ is picking up noise
	lam2p=(L_J-L_lya1)/dpix

	i = 0
	n = (L_sat1-L_lya1)/lam2p
	L[0,i] = 0.0	         	; lya sat#1 height
	L[1,i] = W1 + n		      	; lya sat#1 pos
	L[2,i] = width                  ; lya sat#1 wid
	i++

	L[0,i] = spec[W1]		; lya1 height
	L[1,i] = W1 			; lya1 pos
	L[2,i] = width			; lya1 width
	i++

	n = (L_sat2-L_lya1)/lam2p
	L[0,i] = 0.0        		; lya sat#2 height
	L[1,i] = W1 + n		      	; lya sat#2 pos
	L[2,i] = width                  ; lya sat#2 wid
	i++
	
	n = (L_lya2-L_lya1)/lam2p
	L[0,i] = spec[W1+n]         	; lya2 height
	L[1,i] = W1 + n		      	; lya2 pos
	L[2,i] = width                  ; lya2 wid
	i++

	n = (L_sat3-L_lya1)/lam2p
	L[0,i] = 0.0       	 	; lya sat#3 height
	L[1,i] = W1 + n		      	; lya sat#3 pos
	L[2,i] = width                  ; lya sat#3 wid
	i++

	n = (L_mo4d-L_lya1)/lam2p
	L[0,i] = spec[W1+n]         	; mo4d height
	L[1,i] = W1 + n		      	; mo4d pos
	L[2,i] = width/2.0              ; mo4d wid
	i++

	n = (L_sat4-L_lya1)/lam2p
	L[0,i] = 0.0	         	; lya sat#4 height
	L[1,i] = W1 + n		      	; lya sat#4 pos
	L[2,i] = width                  ; lya sat#4 wid
	i++

	n = (L_T-L_lya1)/lam2p
	L[0,i] = spec[W1+n]         	; T height
	L[1,i] = W1 + n		      	; T pos
	L[2,i] = width                  ; T wid
	i++

	n = (L_K-L_lya1)/lam2p
	L[0,i] = spec[W1+n]         	; K height
	L[1,i] = W1 + n		      	; K pos
	L[2,i] = width                  ; K wid
	i++

	n = (L_Q-L_lya1)/lam2p
	L[0,i] = spec[W1+n]         	; Q height
	L[1,i] = W1 + n		      	; Q pos
	L[2,i] = width                  ; Q wid
	i++

	n = (L_B-L_lya1)/lam2p
	L[0,i] = spec[W1+n]         	; B height
	L[1,i] = W1 + n		      	; B pos
	L[2,i] = width                  ; B wid
	i++

	n = (L_R-L_lya1)/lam2p
	L[0,i] = spec[W1+n]         	; R height
	L[1,i] = W1 + n		      	; R pos
	L[2,i] = width                  ; R wid
	i++

	n = (L_A-L_lya1)/lam2p
	L[0,i] = spec[W1+n]         	; A height
	L[1,i] = W1 + n		      	; A pos
	L[2,i] = width                  ; A wid
	i++

	n = (L_J-L_lya1)/lam2p
	L[0,i] = spec[W1+n]         	; J height
	L[1,i] = W1 + n		      	; J pos
	L[2,i] = width                  ; J wid

	label=['lyas1','lya1','lyas2','lya2','lyas3','4d','lyas4','T','K','Q','B','R','A','J']
	
	base_line = min(spec)
	IF keyword_set(double) THEN estimate = dblarr(n_lines*3+1) ELSE estimate = fltarr(n_lines*3+1)
	FOR i = 0,n_lines-1 DO FOR j = 0,2 DO estimate[3*i+j] = L[j,i]
	estimate[n_lines*3] = base_line 
	spec_est = gaussian_fits(x, estimate)

	IF keyword_set(win) THEN BEGIN
    		openwin, win
    		plot, x, spec, psym = 4, xrange = xrange, thick = 2, title = 'Guesses for Ar'
    		oplot, x, spec_est, color = 100
    		p = estimate
    		FOR i = 0, n_lines-1 DO oplot, x, p[3*i]*exp(-(x-p[3*i+1])^2/(2.*p[3*i+2]^2))+p[n_lines*3], color = 200, line =2, thick =2
	ENDIF

	;fit lya+mo w/ satellites
	start=0
	sub_lines=7
	satlines=[0,2,4,6]
	sub_est=[estimate[start*3:(start+sub_lines)*3-1],estimate[n_lines*3]]	;take lines + baseline offset
	sub_x=x[break[0]:break[1]]
	sub_spec=spec[break[0]:break[1]]

	IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(sub_est)) ELSE $
		parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(sub_est)) 

	
	parinfo(*).value = sub_est	;write starting values
	parinfo(*).limited(0) = 1  	;states that all parameters will have a lower bound
	parinfo(*).limited(1) = 1  	;states that all parameters will have an upper bound

	; intens
	FOR i = 0, sub_lines-1 DO parinfo(3*i+0).limits[0] = 0.0  	;sets the lower bound to zero
	FOR i = 0, sub_lines-1 DO parinfo(3*i+0).limits[1] = max(spec) 	;sets the upper bound to be less than max(spec) 
	FOR j = 0,n(satlines) DO BEGIN 						
		i=satlines[j]
		parinfo[3*i].limits[1]=0.1*(max(spec)-min(spec))	;satellites < 5% of the peak line brightness in spectral region adjusted for background
	ENDFOR

	; pos
	FOR i = 0, sub_lines-1 DO parinfo(3*i+1).limits[0] = L[1,i]-3.0  
	FOR i = 0, sub_lines-1 DO parinfo(3*i+1).limits[1] = L[1,i]+3.0  
	i=2
	parinfo(3*i+1).tied = 'P('+strtrim(3*(i+1)+1,1)+')-'+strtrim((L_lya2-L_sat2)/lam2p,1) 	;tie sat2 pos to lya2
	i=4
	parinfo(3*i+1).tied = 'P('+strtrim(3*(i-1)+1,1)+')+'+strtrim((L_sat3-L_lya2)/lam2p,1) 	;tie sat3 pos to lya2

	;width
	i=1 ;set width for lya1 line
	parinfo(3*i+2).limits[0] = 0.5*width 
	parinfo(3*i+2).limits[1] = 1.25*width 
	i=3 ;set width for lya2 line
	parinfo(3*i+2).limits[0] = 0.5*width 
	parinfo(3*i+2).limits[1] = 1.25*width
	;set satellite width
	FOR j = 0,3 DO BEGIN
		i=satlines[j]
		parinfo(3*i+2).limits[0] = 0.5*width 
		parinfo(3*i+2).limits[1] = 1.5*width
	ENDFOR
	i=5 ;set width for moly line
	parinfo(3*i+2).limits[0] = 0.2*width
	parinfo(3*i+2).limits[1] = 0.9*width


	; baseline
	parinfo[sub_lines*3].limited(*) = 1
	parinfo[sub_lines*3].limits[0] = 0 
	parinfo[sub_lines*3].limits[1] = mean(spec)	

    	p_lya= mpfitfun('gaussian_fits', sub_x,sub_spec,sub_spec^0.5+1, sub_est,parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)
	
	;fit the n=2 satellite collection from T->J
	start=7
	sub_lines=7
	sub_est=[estimate[start*3:(start+sub_lines)*3-1],estimate[n_lines*3]]	;take lines + baseline offset
	sub_x=x[break[2]:break[3]]
	sub_spec=spec[break[2]:break[3]]

	IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(sub_est)) ELSE $
		parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(sub_est)) 

	parinfo(*).value = sub_est	;write starting values
	parinfo(*).limited(0) = 1  	;states that all parameters will have a lower bound
	parinfo(*).limited(1) = 1  	;states that all parameters will have an upper bound

	; intens
	FOR i = 0, sub_lines-1 DO parinfo(3*i+0).limits[0] = 0.0  	;sets the lower bound to zero
	FOR i = 0, sub_lines-1 DO parinfo(3*i+0).limits[1] = max(spec) 	;sets the upper bound to be less than max(spec) 

	; pos
	FOR i = 0, sub_lines-1 DO parinfo(3*i+1).limits[0] = L[1,i+start]-3.0  
	FOR i = 0, sub_lines-1 DO parinfo(3*i+1).limits[1] = L[1,i+start]+3.0  
	i=1
	parinfo(3*i+1).tied = 'P('+strtrim(3*(i-1)+1,1)+')+'+strtrim((L_K-L_T)/lam2p,1) 	;ties T pos to K

	;width
	FOR i = 0, sub_lines-1 DO parinfo(3*i+2).limits[0] = 0.5*width
	FOR i = 0, sub_lines-1 DO parinfo(3*i+2).limits[1] = 1.5*width  
	i=0
	parinfo(3*i+2).tied = 'P('+strtrim(3*(i+6)+2,1)+')'	;ties T and J width

	; baseline
	parinfo[sub_lines*3].limited(*) = 1
	parinfo[sub_lines*3].limits[0] = 0 
	parinfo[sub_lines*3].limits[1] = mean(sub_spec)	

    	p_TJ= mpfitfun('gaussian_fits', sub_x,sub_spec,sub_spec^0.5+1, sub_est,parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)

	p=[p_lya[0:n(p_lya)-1],p_TJ[0:n(p_TJ)-1],last(p_lya)]

	xfit=make(min(x),max(x),2000)
	spec_fit = gaussian_fits(xfit, p)
	resid=spec-gaussian_fits(x,p)
	IF keyword_set(win) THEN BEGIN
		IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
		d_old=!d
		ls=1.0
		IF keyword_set(ps) THEN BEGIN
			device, xsize=7.0,ysize=7.0*700.0/1200.0,/inches
			ls=0.6
		ENDIF ELSE openwin,win+1,xsize=1200,ysize=700
		pos=[0.1,0.35,0.95,0.95]
 		plot, x,spec, psym=8,title='H-like Ar Spectra', xtitle = 'Pixel #', ytitle = 'Counts',pos=pos,chars=1.5*ls,/xsty,symsize=1.0*ls
	    	oplot, xfit, spec_fit, color = 100
    		FOR i = 0, n_lines-1 DO oplot, xfit, p[3*i]*exp(-(xfit-p[3*i+1])^2/(2.*p[3*i+2]^2)), color = 200, line =2
		pos=[0.1,0.05,0.95,0.25]
		plot,x,resid,psym=8,ytit=n2g('Delta')+'Counts',pos=pos,/noerase,/nodata,chars=1.5*ls,/xsty,symsize=1.0*ls
		oplot, [0,max(x)],[0,0],linestyle=2.0
		oplot,x,resid,psym=8,color=200,symsize=1.0*ls
		IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm,/port	
	ENDIF
	IF keyword_set(debug) THEN stop
	output={spec:spec, residual:resid, xfit:xfit, yfit:spec_fit, c0:p[indgen(n_lines)*3+0], c1:p[indgen(n_lines)*3+1] , c2:p[indgen(n_lines)*3+2], off:last(p), lab:label}
	RETURN,output

END

;+
;NAME:
;	HIREXSR_FIT_CA_WN3
;
;MODIFIATION HISTORY:
;	10/6/2014 - M.L. Reinke - applied ILOW/IUP fix
;-

FUNCTION hirexsr_fit_ca_wn3,spec,lam,sig=sig,win=win,double=double,label=label,limit=limit,verb=verb,n_lines=n_lines

	n_lines = 3
	IF keyword_set(double) THEN L = dblarr(3,n_lines) ELSE L=fltarr(3,n_lines)

	;basic check if spectra is worth fitting	
	IF NOT keyword_set(limit) THEN limit=0.0
	IF limit GT 0 THEN BEGIN
		meanoff=total(spec-min(spec))
		IF meanoff LE limit THEN BEGIN
			IF keyword_set(verb) THEN print, 'Too Few Photons: COEFS=-2'
			RETURN,-2		;fail w/ too few photons
		ENDIF
	ENDIF
	ilow=ipt(lam,3.171)
	iup=ipt(lam,3.184)
     	if (ilow eq -1) then ilow = 0
        if (iup eq -1) then iup = n_elements(lam)-1

	IF limit LT 0 AND keyword_set(sig) THEN BEGIN
		snr_limit=abs(limit)
		snr_max=max((spec[ilow:iup]-min(spec))/sig[ilow:iup])
		IF snr_max LE snr_limit THEN BEGIN
			IF keyword_set(verb) THEN print, 'SNR too low: COEFS=-2'
			RETURN,-2		;fail w/ too low of signal to noise ratio
		ENDIF
	ENDIF

	;make w-line a sub-array
	w_line = spec[ilow:iup]
	w_lam=lam[ilow:iup]

	;fit the w-line to get the average line width
	dummy = gaussfit(w_lam,w_line,coefs, nterm=4)
 	W0 = coefs[0]                   ; max(dummy)
	W1 = coefs[1]                   ; where(dummy Eq max(dummy))
	W2 = coefs[2]
	width = W2

	lam_o_shot=-1		;change to model eventually
	hirexsr_load_wavelengths,lam_o_shot,lam_o,z,label
	L_n4=last(lam_o[where(z EQ 20 AND label EQ 'wn4')])
	L_n3=last(lam_o[where(z EQ 20 AND label EQ 'wn3')])
	L_w=last(lam_o[where(z EQ 20 AND label EQ 'w')])


	shift=L_w-w1	;estimate of shift for seeding
	IF abs(shift) GT 0.005 THEN RETURN, -3	;fail with bad seed (.005 > 350 km/s)
	
	;set the initial guess
	;----------------------------

	i = 0
	cent=L_w-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4
	L[0,i] = spec[ipt(lam,cent)]	; w height
	L[1,i] = cent		   	; w pos
	L[2,i] = width                  ; w wid
	i++

	cent=L_n3-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4		;return if seed would push other lines off the spectral range
	L[0,i] = spec[ipt(lam,cent)]	; n3 height
	L[1,i] = cent		   	; n3 pos
	L[2,i] = width                  ; n3 wid
	i++

	cent=L_n4-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4		;return if seed would push other lines off the spectral range
	L[0,i] = spec[ipt(lam,cent)]	; n4 height
	L[1,i] = cent		   	; n4 pos
	L[2,i] = width                  ; n4 wid
	i++

	label=['w','wn3','wn4']
	
	base_line = min(spec)
	IF keyword_set(double) THEN estimate = dblarr(n_lines*3+1) ELSE estimate=fltarr(n_lines*3+1)
	FOR i = 0,n_lines-1 DO FOR j = 0,2 DO estimate[3*i+j] = L[j,i]
	estimate[n_lines*3] = base_line 
	
	;define MPFITFUN constraints
	IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(estimate)) ELSE $
	 	parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(estimate))

	;write starting values
	parinfo[*].value=estimate

	i=1
	;parinfo(3*i+2).tied = 'P('+strtrim(3*(i+1)+2,1)+')' 				;tie n4 width to n3 width
	;parinfo(3*i+1).tied = 'P('+strtrim(3*(i+1)+1,1)+')+'+strtrim((L_n3-L_n4),1)	; Tie n4 shift to n3 shift

	FOR i=0,n_lines-1 DO parinfo[3*i].limited[0]=1		;states that all intensity coefs will have a lower bound
	FOR i=0,n_lines-1 DO parinfo[3*i].limits[0]=0.0      	;sets the lower bound to be 0.0

	;width
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limited[0] = 1  ; states that width parameters will have an lower bound
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limited[1] = 1  ; states that width parameters will have an upper bound
	FOR i = 0,n_lines-1 DO parinfo(3*i+2).limits[0] = 0.5*width ; sets the lower bound to be 50% of w line width 
	FOR i = 0,n_lines-1 DO parinfo(3*i+2).limits[1] = 1.5*width ; sets the upper bound to be 150% of w line width

	;perfom fit
	IF NOT keyword_set(sig) THEN sig=sqrt(spec)
	coefs = mpfitfun('gaussian_fits', lam,spec,sig+1, estimate,/quiet, parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)

	IF keyword_set(win) THEN BEGIN
		resid=spec-gaussian_fits(lam,coefs)
		spec_est=gaussian_fits(lam,estimate)
		spec_fit=gaussian_fits(lam,coefs)
		makesym,10
		IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
		d_old=!d
		ls=1.0
		IF keyword_set(ps) THEN BEGIN
			device, xsize=7.0,ysize=7.0*700.0/1200.0,/inches
			ls=0.6
		ENDIF ELSE openwin,win+1,xsize=1200,ysize=700
		pos=[0.1,0.35,0.95,0.95]
 		plot, lam,spec, psym=8,title='Ca wn3 Spectral Fit', xtitle = 'Pixel #', ytitle = 'Counts',pos=pos,chars=1.5*ls,/xsty,symsize=1.0*ls
	    	oplot, lam, spec_fit, color = 100
    		FOR i = 0, n_lines-1 DO oplot, lam, coefs[3*i]*exp(-(lam-coefs[3*i+1])^2/(2.*coefs[3*i+2]^2)), color = 200, line =2
		pos=[0.1,0.05,0.95,0.25]
		plot,lam,resid,psym=8,ytit=n2g('Delta')+'Counts',pos=pos,/noerase,/nodata,chars=1.5*ls,/xsty,symsize=1.0*ls
		oplot, [0,max(lam)],[0,0],linestyle=2.0
		oplot,lam,resid,psym=8,color=200,symsize=1.0*ls
		IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm,/port	
	ENDIF

	RETURN,coefs
END

;+
;NAME:
;	HIREXSR_FIT_CA_LYA
;
;PURPOSE:
;	This function performs a multi-gaussian fit on the H-like Ca Lyman-a spectrum.  Currently
;	this is done using a 2 gaussian fit for the Lya1, Lya2.
;
;CALLING SEQUENCE:
;	result=HIREXSR_FIT_CA_LYA(spec,lam)
;
;INPUTS:
;	spec	FLTARR [n] of the spectral brightness at each wavelength
;	lam	FLTARR [n] of the wavelengths [Angstroms]
;
;OPTIONAL INPUTS:
;	sig	FLTARR [n] of the uncertainty in the spec data DEFAULT: sqrt(spec)
;	limit	FLOAT	of a lower limit on the signal or signal to noise ratio required to proceed with fitting DEFAULT: 0
;	win	INT	of the window number to plot the estimate and fit+residual DEFAULT: (no plot)
;
;KEYWORD PARAMETERS:
;	/verb sets verbose mode and outputs messages to the terminal
;
;OUTPUTS:
;	result	FLTARR	[3*nlines+1] of the gaussian fit coefficients a*exp(-(x-b)^2/(2*c^2)) along with the baseline as the last element
;
;OPTIONAL OUTPUTS:
;	label	STRARR 	[nlines] of the line labels for the corresponding fit coefficients
;	
;
;PROCEDURE:
;	This function calls GAUSSIAN_FITS using MPFITFUN to perfom the fit.
;	Rest wavelengths are called using HIREXSR_LOAD_WAVELENGTHS
;	The following restrictions are used to constrain the fit
;			 - Lya1 and Lya2 are fully tied together w/ Lya2 at 0.53*Lya1
;			 - 0.1 < Ti < 4.0 keV used as constraint
;			 - -0.001 < dlam < +0.001 used to avoid line wandering
;
;	Two checks are made to ensure proper signal and seeding before MPFITFUN is called.
;		1a) if limit is > 0 then total(spec-mean(spec)) > limit
;		1b) if limit is < 0 and sig is given then max(spec/sig) < |limit|
;
;	If 1a or 1b is not satisfied result=-2 is returned
;
;MODIFICATION HISTORY;
;	Written by	M.L. Reinke (adapted from HIREXSR_FIT_LYA) 3/28/2014
;	4/1/2014	M.L. Reinke - updated the fitting to use a fully constrainted Lya2 line
;
;-

FUNCTION hirexsr_fit_ca_lya,spec,lam,sig=sig,win=win,double=double,label=label,limit=limit,verb=verb,n_lines=n_lines

	n_lines = 2
	IF keyword_set(double) THEN L = dblarr(3,n_lines) ELSE L=fltarr(3,n_lines)

	;basic check if spectra is worth fitting	
	IF NOT keyword_set(limit) THEN limit=0.0
	IF limit GT 0 THEN BEGIN
		meanoff=total(spec-min(spec))
		IF meanoff LE limit THEN BEGIN
			IF keyword_set(verb) THEN print, 'Too Few Photons: COEFS=-2'
			RETURN,-2		;fail w/ too few photons
		ENDIF
	ENDIF
	IF limit LT 0 AND keyword_set(sig) THEN BEGIN
		snr_limit=abs(limit)
		snr_max=max(spec/sig)
		IF snr_max LE snr_limit THEN BEGIN
			IF keyword_set(verb) THEN print, 'SNR too low: COEFS=-2'
			RETURN,-2		;fail w/ too low of signal to noise ratio
		ENDIF
	ENDIF

	;make lya1-line a sub-array
	ilow=ipt(lam,3.015)
	iup=ipt(lam,3.020)
	lya1_line = spec[ilow:iup]
	lya1_lam=lam[ilow:iup]

	;fit the lya1-line to get the average line width
	dummy = gaussfit(lya1_lam,lya1_line,coefs, nterm=4)
 	W0 = coefs[0]                   ; max(dummy)
	W1 = coefs[1]                   ; where(dummy Eq max(dummy))
	W2 = coefs[2]
	width = w2

	lam_o_shot=-1		;change to model eventually
	hirexsr_load_wavelengths,lam_o_shot,lam_o,z,label
	;L_sat1=last(lam_o[where(z EQ 18 AND label EQ 'lyas1')])
   	L_lya1=last(lam_o[where(z EQ 20 AND label EQ 'lya1')])
   	;L_sat2=last(lam_o[where(z EQ 18 AND label EQ 'lyas2')])
	L_lya2=last(lam_o[where(z EQ 20 AND label EQ 'lya2')])
	;L_sat3=last(lam_o[where(z EQ 18 AND label EQ 'lyas3')])
	;L_mo4d=last(lam_o[where(z EQ 42 AND label EQ '4d')])
   	;L_sat4=last(lam_o[where(z EQ 18 AND label EQ 'lyas4')])

	shift=L_lya1-w1	;estimate of shift for seeding
	wmin=hirexsr_iti2iw(0.1,20,L_lya1)
	wmax=hirexsr_iti2iw(4.0,20,L_lya1)
	IF width LT wmin THEN width=1.01*wmin
	IF width GT wmax THEN width=0.99*wmax
	IF abs(shift) GT 0.005 THEN shift=0.0	;fail with bad seed (.005 > 350 km/s)
	
	;set the initial guess
	;----------------------------

	i = 0
	cent=L_lya1-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4	
	L[0,i] = spec[ipt(lam,cent)]	; lya1 height
	L[1,i] = cent		   	; lya1  pos
	L[2,i] = width                  ; lya1 wid
	i++

	cent=L_lya2-shift
	IF cent GT max(lam) OR cent LT min(lam) THEN RETURN,-4	
	L[0,i] = spec[ipt(lam,cent)]	; lya2 height
	L[1,i] = cent		   	; lya2  pos
	L[2,i] = width                  ; lya2 wid
	i++

	label=['lya1','lya2']

	base_line = 0.5*(spec[0]+last(spec))
	L[0,*]=L[0,*]-base_line > 0.0
	IF keyword_set(double) THEN estimate = dblarr(n_lines*3+1) ELSE estimate=fltarr(n_lines*3+1)
	FOR i = 0,n_lines-1 DO FOR j = 0,2 DO estimate[3*i+j] = L[j,i]
	estimate[n_lines*3] = base_line 
	
	;define MPFITFUN constraints
	IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(estimate)) ELSE $
	 	parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(estimate))

	;write starting values
	parinfo[*].value=estimate
	
	;intens
	FOR i=0,n_lines-1 DO parinfo[3*i].limited[0]=1		;states that all intensity coefs will have a lower bound
	FOR i=0,n_lines-1 DO parinfo[3*i].limits[0]=0.0      	;sets the lower bound to be 0.0

	i=1
	parinfo(3*i+0).tied = '0.53*P('+strtrim(3*(i-1),1)+')'
	parinfo(3*i+1).tied = 'P('+strtrim(3*(i-1)+1,1)+')+'+strtrim((L_lya2-L_lya1),1) ;
	parinfo(3*i+2).tied = 'P('+strtrim(3*(i-1)+2,1)+')' ;

	; pos
	FOR i=0,n_lines-1 DO parinfo[3*i+1].limited[0]=1			;line will have lower bound on shift
	FOR i=0,n_lines-1 DO parinfo[3*i+1].limited[1]=1			;line will have an upper bound on shift					
	FOR i=0,n_lines-1 DO parinfo(3*i+1).limits[0] = L[1,i]-0.001		;prevent wander of lines 
	FOR i=0,n_lines-1 DO parinfo(3*i+1).limits[1] = L[1,i]+0.001

	;width
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limited[0]=1	;states that all width coefs will have a lower bound
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limits[0]=wmin  	;sets the lower bound to be ~100 eV
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limited[1]=1	;states that all width coefs will have a lower bound
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limits[1]=wmax  	;sets the lower bound to be ~4.0 keV eV


	; baseline
;	parinfo[n_lines*3].limited(*) = 1
;	parinfo[n_lines*3].limits[0] = 0 
;	parinfo[n_lines*3].limits[1] = mean(spec)

	;perfom fit
	IF NOT keyword_set(sig) THEN sig=sqrt(spec)
	coefs = mpfitfun('gaussian_fits', lam,spec,sig+1, estimate, parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)

	;stop
	;win=10
	IF keyword_set(win) THEN BEGIN
		resid=spec-gaussian_fits(lam,coefs)
		spec_est=gaussian_fits(lam,estimate)
		spec_fit=gaussian_fits(lam,coefs)
		makesym,10
		IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
		d_old=!d
		ls=1.0
		IF keyword_set(ps) THEN BEGIN
			device, xsize=7.0,ysize=7.0*700.0/1200.0,/inches
			ls=0.6
		ENDIF ELSE openwin,win+1,xsize=1200,ysize=700
		pos=[0.1,0.35,0.95,0.95]
 		plot, lam,spec, psym=8,title='Ca Lya Spectral Fit', xtitle = 'Pixel #', ytitle = 'Counts',pos=pos,chars=1.5*ls,/xsty,symsize=1.0*ls
	    	oplot, lam, spec_fit, color = 100
    		FOR i = 0, n_lines-1 DO oplot, lam, coefs[3*i]*exp(-(lam-coefs[3*i+1])^2/(2.*coefs[3*i+2]^2)), color = 200, line =2
		pos=[0.1,0.05,0.95,0.25]
		plot,lam,resid,psym=8,ytit=n2g('Delta')+'Counts',pos=pos,/noerase,/nodata,chars=1.5*ls,/xsty,symsize=1.0*ls
		oplot, [0,max(lam)],[0,0],linestyle=2.0
		oplot,lam,resid,psym=8,color=200,symsize=1.0*ls
		IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm,/port
	ENDIF
	;stop

	RETURN,coefs
END

;+
;NAME:
;	HIREXSR_FIT_HE_CA
;
;PURPOSE:
;	This function is used to perform multigaussian fitting for the He-like Ca spectra in pixel space for use with W_HIREXSR_CALIB
;
;CALLING SEQUENCE:
;	result=HIREXSR_FIT_HE_CA(spec)
;
;MODFICATION HISTORY
;	Written by	Y. Podpaly (adapted from HIREXSR_FIT_HE) - 8/10
;	9/16/10		M.L. Reinke - increased range to 9 pixels	
;	1/14/10		M.L. Reinke - tweaked the the /focus fitting parameters based on data from 1101014
;	12/3/12		M.L. Reinke - changed name to HE_CA to allow for H_CA naming scheme
;-

FUNCTION hirexsr_fit_he_ca,spec,win=win,double=double,focus=focus,break=break
	IF NOT keyword_set(break) THEN break=[40,65,70,170]

	lam_o_shot=-1		;change to model eventually
	hirexsr_load_wavelengths,lam_o_shot,lam_o,z,label
	L_n4=lam_o[where(z EQ 20 AND label EQ 'wn4')]*1.0e3
	L_n3=lam_o[where(z EQ 20 AND label EQ 'wn3')]*1.0e3
	L_w=lam_o[where(z EQ 20 AND label EQ 'w')]*1.0e3
	L_x=lam_o[where(z EQ 20 AND label EQ 'x')]*1.0e3
	L_y=lam_o[where(z EQ 20 AND label EQ 'y')]*1.0e3
	L_z=lam_o[where(z EQ 20 AND label EQ 'z')]*1.0e3
	L_j=lam_o[where(z EQ 20 AND label EQ 'j')]*1.0e3
	L_q=lam_o[where(z EQ 20 AND label EQ 'q')]*1.0e3
	L_r=lam_o[where(z EQ 20 AND label EQ 'r')]*1.0e3
	L_a=lam_o[where(z EQ 20 AND label EQ 'a')]*1.0e3
	L_k=lam_o[where(z EQ 20 AND label EQ 'k')]*1.0e3
	L_4n=lam_o[where(z EQ 18 AND label EQ 'w4n')]*1.0e3
	
	;make w-line sub-array
	center = maxloc(spec[break[0]:break[1]])+break[0]	;maximum location of pixel
	range=9
	nspec=n(spec)
	lowpix=0 > (center-range)				;prevent from indexing off the array if line is close to the edge
	highpix=nspec < (center+range)		
	w_line = spec[lowpix:highpix]	
	x_w=indgen(n(w_line)+1)
	dummy = gaussfit(x_w,w_line, coefs, nterm=4)
	W0 = coefs[0]                    
	P1 = coefs[1]+lowpix      	
	W2 = coefs[2]
	width = W2

	p4n=maxloc(spec[break[2]:break[3]])+break[2]
    	lam2p = (L_4n - L_w)/(p4n-p1)				;use He-like Ar 4n and He-like Ca w line to determine average dispersion
	pORlam = 1.0                				;no idea what this does, but vestigial from Alex's POS code
    	x = findgen(n_elements(spec))
    	W1 = P1
    	xrange = [ 0, n_elements(spec)]

	IF keyword_set(focus) THEN n_lines = 12 ELSE n_lines=9
	IF keyword_set(double) THEN L = dblarr(3,n_lines) ELSE l=fltarr(3,n_lines)

	;set the initial guess
	;----------------------------

	i = 0
	L[0,i] = W0                     ; w height
	L[1,i] = W1                     ; w pos
	L[2,i] = width                  ; w width
	i++

	IF keyword_set(focus) THEN BEGIN
		n_wn4 = (L_n4-L_w)/lam2p
		L[0,i] = spec[P1+n_wn4]*0.0         ; n=4 Sats height
		L[1,i] = W1 + n_wn4*pORLam      ; n=4 Sats pos
		L[2,i] = width                  ; n=4 Sats wid
		i++
	ENDIF

	n_wn3 = (L_n3-L_w)/lam2p
	L[0,i] = spec[P1+n_wn3]*0.0     ; n=3 Sats height
	IF keyword_set(focus) THEN L[0,i] = spec[P1+n_wn3]
	L[1,i] = W1 + n_wn3*pORLam      ; n=3 Sats pos
	L[2,i] = width                  ; n=3 Sats wid
	i++

	n_wx =  (L_x-L_w)/lam2p 
	L[0,i] = spec[P1+n_wx]          ; x height
	L[1,i] = W1+n_wx*pORLam         ; x pos
	L[2,i] = width                  ; x wid
	i++

	n_wy = (L_y-L_w)/lam2p 
	L[0,i] = spec[P1+n_wy]          ; y height
	L[1,i] = W1 + n_wy*pORLam       ; y pos
	L[2,i] = width                  ; y wid
	i++

	n_w4n = (L_4n-L_w)/lam2p
	L[0,i] = spec[P1+n_w4n]		; He-like Ar n=4 height
	L[1,i] = W1 + n_w4n*pORLam      ; He-like Ar n=4 Sats pos
	L[2,i] = width                  ; He-like Ar n=4 Sats wid
	i++

	IF keyword_set(focus) THEN BEGIN
		n_wq = (L_q - l_w)/lam2p 
		L[0,i] = 0.0          		; q height
		L[1,i] = W1+n_wq*pORLam         ; q pos
		L[2,i] = width                  ; q wid
		i++
	ENDIF

	n_wr = (L_r - l_w)/lam2p 
	L[0,i] = spec[P1+n_wr]*0.0      ; r height
	L[1,i] = W1 + n_wr*pORLam       ; r pos
	L[2,i] = width                  ; r wid
	i++

	n_wa = (L_a - l_w)/lam2p 
	L[0,i] = spec[P1+n_wa]          ; a height
	L[1,i] = W1 + n_wa*pORLam       ; a pos
	L[2,i] =  width                 ; a wid
	i++

	n_wk = (L_k - l_w)/lam2p    
	L[0,i] = spec[P1+n_wk]          ; k height
	L[1,i] = W1 + n_wk*pORLam       ; k pos
	L[2,i] = width                  ; k wid
	i++

	IF keyword_set(focus) THEN BEGIN
		n_wj = (L_j - l_w)/lam2p    
		L[0,i] = spec[P1+n_wk]	        ; j height
		L[1,i] = W1 + n_wj*pORLam       ; j pos
		L[2,i] = width                  ; j wid
		i++

	ENDIF

	n_wz = (L_z - l_w)/lam2p      
	L[0,i] = spec[P1+n_wz]          ; z height	
	L[1,i] = W1 + n_wz*pORLam   	; z pos; 5 adjustment for 2009
	L[2,i] = width                  ; z wid
	i++

	
       	IF keyword_set(focus) THEN BEGIN
            label=['w','n4','n3','x','y', 'w4n', 'q', 'r','a','k','j', 'z']
            ;       0    1    2   3   4    5      6    7   8   9   10   11 
        ENDIF ELSE BEGIN
            label=['w','n3','x','y','w4n','r','a','k','z']
            ;       0   1    2   3    4    5   6   7   8    
        ENDELSE
	
	base_line = min(spec)
	IF keyword_set(double) THEN estimate = dblarr(n_lines*3+1) ELSE estimate=fltarr(n_lines*3+1)
	FOR i = 0,n_lines-1 DO FOR j = 0,2 DO estimate[3*i+j] = L[j,i]
	estimate[n_lines*3] = base_line 
	spec_est = gaussian_fits(x, estimate)
	
	IF keyword_set(win) THEN BEGIN
    		openwin, win
    		plot, x, spec, psym = 4, xrange = xrange, thick = 2, title = 'Guesses for Ar'
    		oplot, x, spec_est, color = 100
    		p = estimate
    		FOR i = 0, n_lines-1 DO oplot, x, p[3*i]*exp(-(x-p[3*i+1])^2/(2.*p[3*i+2]^2))+p[n_lines*3], color = 200, line =2, thick =2
	ENDIF

	;fit full spectrum (spectrum is compressed relative to He-like Ar and I don't feel like coming up with a breaking algorithm)

	;define MPFITFUN constraints
	IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(estimate)) ELSE $
	 	parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(estimate))
	parinfo(*).value = estimate

	IF keyword_set(focus) THEN BEGIN
		; Tie n4 shift to n3 shift
		i = 1
		parinfo(3*i+1).tied = 'P('+strtrim(3*(i+1)+1,1)+')-'+strtrim((L_n3-L_n4)/lam2p,1)

		; tie n4 width to n3 width
		parinfo(3*i+2).tied = 'P('+strtrim(3*(i+1)+2,1)+')' ;
	ENDIF

	IF keyword_set(focus) THEN BEGIN
		;tie q shift to r a shift
		i=6
		parinfo(3*i+1).tied = 'P('+strtrim(3*(i+2)+1,1)+')-'+strtrim((L_a-L_q)/lam2p,1)

		; tie q width to a width
		parinfo(3*i+2).tied = 'P('+strtrim(3*(i+2)+2,1)+')' ;
	ENDIF

	IF keyword_set(focus) THEN i=7 ELSE i=5
	;Tie r shift to a shift
	parinfo(3*i+1).tied = 'P('+strtrim(3*(i+1)+1,1)+')-'+strtrim((L_a-L_r)/lam2p,1)

	; tie r width to a width
	parinfo(3*i+2).tied = 'P('+strtrim(3*(i+1)+2,1)+')' ;


	IF keyword_set(focus) THEN BEGIN
		i=10
		parinfo(3*i+1).tied = 'P('+strtrim(3*(i+1)+1,1)+')-'+strtrim((L_z-L_j)/lam2p,1)
	ENDIF

	parinfo[*].limited(0) = 1  ; states that all parameters will have an lower bound
	parinfo[*].limited(1) = 1  ; states that all parameters will have an upper bound


	; intens
	FOR i = 0, n_lines-1 DO parinfo[3*i+0].limits[0] = 0.0      		;sets the lower bound to be 0.0
	FOR i = 0, n_lines-1 DO parinfo[3*i+0].limits[1] = max(spec) 		;1.1 sets the upper bound to be less than max(spec) 

	;pos
	FOR i = 0, n_lines-1 DO parinfo(3*i+1).limits[0] = L[1,i]-3.0  ;sets the lower bound to be -3 pix of the estimate
	FOR i = 0, n_lines-1 DO parinfo(3*i+1).limits[1] = L[1,i]+3.0  ;sets the upper bound to be +3 pix of the estimate 

	;width
	FOR i = 0, n_lines-1 DO parinfo(3*i+2).limits[0] = 0.2*width ; sets the lower bound to be 80% of w line width
	FOR i = 0, n_lines-1 DO parinfo(3*i+2).limits[1] = 1.5*width ; sets the upper bound to be 150% of w line width

	; baseline
	parinfo[n_lines*3].limited(*) = 1
	parinfo[n_lines*3].limits[0] = 0 
	parinfo[n_lines*3].limits[1] = mean(spec)

	;perfom fit 
	p = mpfitfun('gaussian_fits', x,spec,spec^0.5+1, estimate, parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)
	;print, 'MPFITFUN returned '+num2str(status,1)+' after '+num2str(niter,1)+' iterations.


	xfit=make(min(x),max(x),2000)
	spec_fit = gaussian_fits(xfit, p)
	resid=spec-gaussian_fits(x,p)

	IF keyword_set(win) THEN BEGIN
		IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
		d_old=!d
		ls=1.0
		IF keyword_set(ps) THEN BEGIN
			device, xsize=7.0,ysize=7.0*700.0/1200.0,/inches
			ls=0.6
		ENDIF ELSE openwin,win+1,xsize=1200,ysize=700
		pos=[0.1,0.35,0.95,0.95]
 		plot, x,spec, psym=8,title='He-like Ca Spectra', xtitle = 'Pixel #', ytitle = 'Counts',pos=pos,chars=1.5*ls,/xsty,symsize=1.0*ls
	    	oplot, xfit, spec_fit, color = 100
    		FOR i = 0, n_lines-1 DO oplot, xfit, p[3*i]*exp(-(xfit-p[3*i+1])^2/(2.*p[3*i+2]^2)), color = 200, line =2
		pos=[0.1,0.05,0.95,0.25]
		plot,x,resid,psym=8,ytit=n2g('Delta')+'Counts',pos=pos,/noerase,/nodata,chars=1.5*ls,/xsty,symsize=1.0*ls
		oplot, [0,max(x)],[0,0],linestyle=2.0
		oplot,x,resid,psym=8,color=200,symsize=1.0*ls
		IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm,/port	
	ENDIF
	output={spec:spec, residual:resid, xfit:xfit, yfit:spec_fit, c0:p[indgen(n_lines)*3+0], c1:p[indgen(n_lines)*3+1] , c2:p[indgen(n_lines)*3+2], off:last(p), lab:label}

	RETURN,output
END

;+
;NAME:
;	HIREXSR_FIT_HE_CA_SUB
;
;PURPOSE:
;	This function is used to perform multigaussian fitting for the He-like Ca spectra in pixel space for use with W_HIREXSR_CALIB.
;	Unlike HIREXSR_FIT_HE_CA, this uses only the w and x y lines to allow for fitting FY14 data which had the spectrum shifted to avoid a bad sub module.
;
;CALLING SEQUENCE:
;	result=HIREXSR_FIT_HE_CA_MOD(spec)
;
;MODFICATION HISTORY
;	Written by	M.L. Reinke 8-30-2014 (adapted from HIREXSR_FIT_HE_CA)
;
;-

FUNCTION hirexsr_fit_he_ca_mod,spec,win=win,double=double,focus=focus,break=break
	IF NOT keyword_set(break) THEN break=[110,146,146,178]

	lam_o_shot=-1		;change to model eventually
	hirexsr_load_wavelengths,lam_o_shot,lam_o,z,label
	L_n4=lam_o[where(z EQ 20 AND label EQ 'wn4')]*1.0e3
	L_n3=lam_o[where(z EQ 20 AND label EQ 'wn3')]*1.0e3
	L_w=lam_o[where(z EQ 20 AND label EQ 'w')]*1.0e3
	L_x=lam_o[where(z EQ 20 AND label EQ 'x')]*1.0e3
	L_y=lam_o[where(z EQ 20 AND label EQ 'y')]*1.0e3
	L_4n=lam_o[where(z EQ 18 AND label EQ 'w4n')]*1.0e3
	
	;make w-line sub-array
	center = maxloc(spec[break[0]:break[1]])+break[0]	;maximum location of pixel
	range=9
	nspec=n(spec)
	lowpix=0 > (center-range)				;prevent from indexing off the array if line is close to the edge
	highpix=nspec < (center+range)		
	w_line = spec[lowpix:highpix]	
	x_w=indgen(n(w_line)+1)
	dummy = gaussfit(x_w,w_line, coefs, nterm=4)
	W0 = coefs[0]                    
	P1 = coefs[1]+lowpix      	
	W2 = coefs[2]
	width = W2

	p4n=maxloc(spec[break[2]:*])+break[2]
    	lam2p = (L_4n - L_w)/(p4n-p1)				;use He-like Ar 4n and He-like Ca w line to determine average dispersion
	pORlam = 1.0                				;no idea what this does, but vestigial from Alex's POS code
    	x = findgen(n_elements(spec))
    	W1 = P1
    	xrange = [ 0, n_elements(spec)]

	IF keyword_set(focus) THEN n_lines = 5 ELSE n_lines=4
	IF keyword_set(double) THEN L = dblarr(3,n_lines) ELSE l=fltarr(3,n_lines)

	;set the initial guess
	;----------------------------

	i = 0
	L[0,i] = W0                     ; w height
	L[1,i] = W1                     ; w pos
	L[2,i] = width                  ; w width
	i++

	IF keyword_set(focus) THEN BEGIN
		n_wn4 = (L_n4-L_w)/lam2p
		L[0,i] = spec[P1+n_wn4]*0.0         ; n=4 Sats height
		L[1,i] = W1 + n_wn4*pORLam      ; n=4 Sats pos
		L[2,i] = width                  ; n=4 Sats wid
		i++
	ENDIF

	n_wn3 = (L_n3-L_w)/lam2p
	L[0,i] = spec[P1+n_wn3]*0.0     ; n=3 Sats height
	IF keyword_set(focus) THEN L[0,i] = spec[P1+n_wn3]
	L[1,i] = W1 + n_wn3*pORLam      ; n=3 Sats pos
	L[2,i] = width                  ; n=3 Sats wid
	i++

	n_wx =  (L_x-L_w)/lam2p 
	L[0,i] = spec[P1+n_wx]          ; x height
	L[1,i] = W1+n_wx*pORLam         ; x pos
	L[2,i] = width                  ; x wid
	i++

	n_wy = (L_y-L_w)/lam2p 
	L[0,i] = spec[P1+n_wy]          ; y height
	L[1,i] = W1 + n_wy*pORLam       ; y pos
	L[2,i] = width                  ; y wid
	i++
	
       	IF keyword_set(focus) THEN BEGIN
            label=['w','wn4','wn3','x','y']
            ;       0    1    2   3   4 
        ENDIF ELSE BEGIN
            label=['w','wn3','x','y']
            ;       0   1    2   3 
        ENDELSE
	
	base_line = min(spec)
	IF keyword_set(double) THEN estimate = dblarr(n_lines*3+1) ELSE estimate=fltarr(n_lines*3+1)
	FOR i = 0,n_lines-1 DO FOR j = 0,2 DO estimate[3*i+j] = L[j,i]
	estimate[n_lines*3] = base_line 
	spec_est = gaussian_fits(x, estimate)
	
	IF keyword_set(win) THEN BEGIN
    		openwin, win
    		plot, x, spec, psym = 4, xrange = xrange, thick = 2, title = 'Guesses for Ar'
    		oplot, x, spec_est, color = 100
    		p = estimate
    		FOR i = 0, n_lines-1 DO oplot, x, p[3*i]*exp(-(x-p[3*i+1])^2/(2.*p[3*i+2]^2))+p[n_lines*3], color = 200, line =2, thick =2
	ENDIF

	;fit spectrum from break[0]:break[3]
	sub_x=x[break[0]:break[3]]
	sub_spec=spec[break[0]:break[3]]

	;define MPFITFUN constraints
	IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(estimate)) ELSE $
	 	parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(estimate))
	parinfo(*).value = estimate

	IF keyword_set(focus) THEN BEGIN
		; Tie n4 shift to n3 shift
		i = 1
		parinfo(3*i+1).tied = 'P('+strtrim(3*(i+1)+1,1)+')-'+strtrim((L_n3-L_n4)/lam2p,1)

		; tie n4 width to n3 width
		parinfo(3*i+2).tied = 'P('+strtrim(3*(i+1)+2,1)+')' ;
	ENDIF

	; intens
	FOR i = 0, n_lines-1 DO parinfo[3*i+0].limits[0] = 0.0      		;sets the lower bound to be 0.0
	FOR i = 0, n_lines-1 DO parinfo[3*i+0].limits[1] = max(spec) 		;1.1 sets the upper bound to be less than max(spec) 

	;pos
	FOR i = 0, n_lines-1 DO parinfo(3*i+1).limits[0] = L[1,i]-3.0  ;sets the lower bound to be -3 pix of the estimate
	FOR i = 0, n_lines-1 DO parinfo(3*i+1).limits[1] = L[1,i]+3.0  ;sets the upper bound to be +3 pix of the estimate 

	;width
	FOR i = 0, n_lines-1 DO parinfo(3*i+2).limits[0] = 0.2*width ; sets the lower bound to be 80% of w line width
	FOR i = 0, n_lines-1 DO parinfo(3*i+2).limits[1] = 1.5*width ; sets the upper bound to be 150% of w line width

	; baseline
	parinfo[n_lines*3].limited(*) = 1
	parinfo[n_lines*3].limits[0] = 0 
	parinfo[n_lines*3].limits[1] = mean(spec)

	;perfom fit 
	p = mpfitfun('gaussian_fits', sub_x,sub_spec,sub_spec^0.5+1, estimate, parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)
	;print, 'MPFITFUN returned '+num2str(status,1)+' after '+num2str(niter,1)+' iterations.


	xfit=make(min(x),max(x),2000)
	spec_fit = gaussian_fits(xfit, p)
	resid=spec-gaussian_fits(x,p)

	IF keyword_set(win) THEN BEGIN
		IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
		d_old=!d
		ls=1.0
		IF keyword_set(ps) THEN BEGIN
			device, xsize=7.0,ysize=7.0*700.0/1200.0,/inches
			ls=0.6
		ENDIF ELSE openwin,win+1,xsize=1200,ysize=700
		pos=[0.1,0.35,0.95,0.95]
 		plot, x,spec, psym=8,title='He-like Ca Spectra', xtitle = 'Pixel #', ytitle = 'Counts',pos=pos,chars=1.5*ls,/xsty,symsize=1.0*ls
	    	oplot, xfit, spec_fit, color = 100
    		FOR i = 0, n_lines-1 DO oplot, xfit, p[3*i]*exp(-(xfit-p[3*i+1])^2/(2.*p[3*i+2]^2)), color = 200, line =2
		pos=[0.1,0.05,0.95,0.25]
		plot,x,resid,psym=8,ytit=n2g('Delta')+'Counts',pos=pos,/noerase,/nodata,chars=1.5*ls,/xsty,symsize=1.0*ls
		oplot, [0,max(x)],[0,0],linestyle=2.0
		oplot,x,resid,psym=8,color=200,symsize=1.0*ls
		IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm,/port	
	ENDIF
	output={spec:spec, residual:resid, xfit:xfit, yfit:spec_fit, c0:p[indgen(n_lines)*3+0], c1:p[indgen(n_lines)*3+1] , c2:p[indgen(n_lines)*3+2], off:last(p), lab:label}

	RETURN,output
END

;+
;NAME:
;	HIREXSR_FIT_H_CA
;
;PURPOSE:
;	This function is used to perform multigaussian fitting for the H-like Ca spectra in pixel space for use with W_HIREXSR_CALIB
;
;CALLING SEQUENCE:
;	result=HIREXSR_FIT_H_CA(spec)
;
;MODFICATION HISTORY
;	Written by	M.L. Reinke  (adapted from HIREXSR_FIT_HE_CA) - 12/12
;
;-

FUNCTION hirexsr_fit_h_ca,spec,win=win,double=double,focus=focus,break=break
	IF NOT keyword_set(break) THEN break=[50,100,160,190]

	lam_o_shot=-1		;change to model eventually
	hirexsr_load_wavelengths,lam_o_shot,lam_o,z,label
	L_7n=lam_o[where(z EQ 18 AND label EQ 'w7n')]*1.0e3
	L_8n=lam_o[where(z EQ 18 AND label EQ 'w8n')]*1.0e3
	L_9n=lam_o[where(z EQ 18 AND label EQ 'w9n')]*1.0e3
	L_10n=lam_o[where(z EQ 18 AND label EQ 'w10n')]*1.0e3
	L_11n=lam_o[where(z EQ 18 AND label EQ 'w11n')]*1.0e3
	L_12n=lam_o[where(z EQ 18 AND label EQ 'w12n')]*1.0e3

	
	;make 7n-line sub-array
	center = maxloc(spec[break[2]:break[3]])+break[2]	;maximum location of pixel
	range=9
	nspec=n(spec)
	lowpix=0 > (center-range)				;prevent from indexing off the array if line is close to the edge
	highpix=nspec < (center+range)		
	w_line = spec[lowpix:highpix]	
	x_w=indgen(n(w_line)+1)
	dummy = gaussfit(x_w,w_line, coefs, nterm=4)
	W0 = coefs[0]                    
	P1 = coefs[1]+lowpix      	
	W2 = coefs[2]
	width = W2
	bl=coefs[3]

	p10n=maxloc(spec[break[0]:break[1]])+break[0]
    	lam2p = (L_7n - L_10n)/(p10n-p1)			;use He-like Ar 7n and He-like Ar 10n line to determine average dispersion
	pORlam = 1.0                				;no idea what this does, but vestigial from Alex's POS code
    	x = findgen(n_elements(spec))
    	W1 = P1
    	xrange = [ 0, n_elements(spec)]

	IF keyword_set(focus) THEN n_lines = 12 ELSE n_lines=9
	IF keyword_set(double) THEN L = dblarr(3,n_lines) ELSE l=fltarr(3,n_lines)

	n_lines=6

	;set the initial guess
	;----------------------------

	i = 0
	n_w12 =  (L_7n-L_12n)/lam2p 
	L[0,i] = spec[P1+n_w12]-bl          ; 12n height
	L[1,i] = W1+n_w12*pORLam         ; 12n pos
	L[2,i] = width                  ; 12n wid
	i++

	n_w11 =  (L_7n-L_11n)/lam2p 
	L[0,i] = spec[P1+n_w11]-bl          ; 11n height
	L[1,i] = W1+n_w11*pORLam         ; 11n pos
	L[2,i] = width                  ; 11n wid
	i++

	n_w10 =  (L_7n-L_10n)/lam2p 
	L[0,i] = spec[P1+n_w10]-bl          ; 10n height
	L[1,i] = W1+n_w10*pORLam         ; 10n pos
	L[2,i] = width                  ; 10n wid
	i++

	n_w9 =  (L_7n-L_9n)/lam2p 
	L[0,i] = spec[P1+n_w9]-bl          ; 9n height
	L[1,i] = W1+n_w9*pORLam         ; 9n pos
	L[2,i] = width                  ; 9n wid
	i++

	n_w8 =  (L_7n-L_8n)/lam2p 
	L[0,i] = spec[P1+n_w8]-bl          ; 8n height
	L[1,i] = W1+n_w8*pORLam         ; 8n pos
	L[2,i] = width                  ; 8n wid
	i++

	L[0,i] = W0                     ; w7n height
	L[1,i] = W1                     ; w7n pos
	L[2,i] = width                  ; w7n width
	i++

        label=['w12n','w11n','w10n','w9n','w8n','w7n']
        ;       0      1      2      3     4     5
	
	base_line = min(spec)
	IF keyword_set(double) THEN estimate = dblarr(n_lines*3+1) ELSE estimate=fltarr(n_lines*3+1)
	FOR i = 0,n_lines-1 DO FOR j = 0,2 DO estimate[3*i+j] = L[j,i]
	estimate[n_lines*3] = base_line 
	spec_est = gaussian_fits(x, estimate)
	
	IF keyword_set(win) THEN BEGIN
    		openwin, win
    		plot, x, spec, psym = 4, xrange = xrange, thick = 2, title = 'Guesses for High-n Ar'
    		oplot, x, spec_est, color = 100
    		p = estimate
    		FOR i = 0, n_lines-1 DO oplot, x, p[3*i]*exp(-(x-p[3*i+1])^2/(2.*p[3*i+2]^2))+p[n_lines*3], color = 200, line =2, thick =2
	ENDIF

	;fit lines <= n=10
	start=0
	sub_lines=3
	sub_est=[estimate[start*3:(start+sub_lines)*3-1],estimate[n_lines*3]]	;take lines + baseline offset
	sub_x=x[break[0]:break[1]]
	sub_spec=spec[break[0]:break[1]]

	IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(sub_est)) ELSE $
		parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(sub_est)) 

	parinfo(*).value = sub_est	;write starting values
	parinfo(*).limited(0) = 1  	;states that all parameters will have a lower bound
	parinfo(*).limited(1) = 1  	;states that all parameters will have an upper bound

	; intens
	FOR i = 0, sub_lines-1 DO parinfo(3*i+0).limits[0] = 0.0  	;sets the lower bound to zero
	FOR i = 0, sub_lines-1 DO parinfo(3*i+0).limits[1] = max(spec) 	;sets the upper bound to be less than max(spec) 

	; pos
	FOR i = 0, sub_lines-1 DO parinfo(3*i+1).limits[0] = L[1,i]-3.0  
	FOR i = 0, sub_lines-1 DO parinfo(3*i+1).limits[1] = L[1,i]+3.0  

	; width
	FOR i = 0, sub_lines-1 DO parinfo(3*i+2).limits[0] = 0.5*width
	FOR i = 0, sub_lines-1 DO parinfo(3*i+2).limits[1] = 1.25*width

	; baseline
	parinfo[sub_lines*3].limited(*) = 1
	parinfo[sub_lines*3].limits[0] = 0 
	parinfo[sub_lines*3].limits[1] = mean(spec[break[1]:*])	

    	p_highn= mpfitfun('gaussian_fits', sub_x,sub_spec,sub_spec^0.5+1, sub_est,parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)

	;fit lines > 10
	start=3
	sub_lines=3
	sub_est=[estimate[start*3:(start+sub_lines)*3-1],estimate[n_lines*3]]	;take lines + baseline offset
	sub_x=x[break[1]:break[3]]
	sub_spec=spec[break[1]:break[3]]

	IF keyword_set(double) THEN parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0], tied:''}, n_elements(sub_est)) ELSE $
		parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(sub_est)) 

	parinfo(*).value = sub_est	;write starting values
	parinfo(*).limited(0) = 1  	;states that all parameters will have a lower bound
	parinfo(*).limited(1) = 1  	;states that all parameters will have an upper bound

	; intens
	FOR i = 0, sub_lines-1 DO parinfo(3*i+0).limits[0] = 0.0  	;sets the lower bound to zero
	FOR i = 0, sub_lines-1 DO parinfo(3*i+0).limits[1] = max(spec) 	;sets the upper bound to be less than max(spec) 

	; pos
	FOR i = 0, sub_lines-1 DO parinfo(3*i+1).limits[0] = L[1,i+start]-3.0  
	FOR i = 0, sub_lines-1 DO parinfo(3*i+1).limits[1] = L[1,i+start]+3.0  

	; width
	FOR i = 0, sub_lines-1 DO parinfo(3*i+2).limits[0] = 0.5*width
	FOR i = 0, sub_lines-1 DO parinfo(3*i+2).limits[1] = 1.25*width

	; baseline
	parinfo[sub_lines*3].limited(*) = 1
	parinfo[sub_lines*3].limits[0] = 0 
	parinfo[sub_lines*3].limits[1] = mean(spec[0:break[1]])	
	
    	p_lown= mpfitfun('gaussian_fits', sub_x,sub_spec,sub_spec^0.5+1, sub_est,parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30)

	p=[p_highn[0:n(p_highn)-1],p_lown[0:n(p_lown)-1],last(p_lown)]

	xfit=make(min(x),max(x),2000)
	spec_fit = gaussian_fits(xfit, p)
	resid=spec-gaussian_fits(x,p)
	IF keyword_set(win) THEN BEGIN
		IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
		d_old=!d
		ls=1.0
		IF keyword_set(ps) THEN BEGIN
			device, xsize=7.0,ysize=7.0*700.0/1200.0,/inches
			ls=0.6
		ENDIF ELSE openwin,win+1,xsize=1200,ysize=700
		pos=[0.1,0.35,0.95,0.95]
 		plot, x,spec, psym=8,title='H-like Ca Spectra', xtitle = 'Pixel #', ytitle = 'Counts',pos=pos,chars=1.5*ls,/xsty,symsize=1.0*ls
	    	oplot, xfit, spec_fit, color = 100
    		FOR i = 0, n_lines-1 DO oplot, xfit, p[3*i]*exp(-(xfit-p[3*i+1])^2/(2.*p[3*i+2]^2)), color = 200, line =2
		pos=[0.1,0.05,0.95,0.25]
		plot,x,resid,psym=8,ytit=n2g('Delta')+'Counts',pos=pos,/noerase,/nodata,chars=1.5*ls,/xsty,symsize=1.0*ls
		oplot, [0,max(x)],[0,0],linestyle=2.0
		oplot,x,resid,psym=8,color=200,symsize=1.0*ls
		IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm,/port	
	ENDIF
	output={spec:spec, residual:resid, xfit:xfit, yfit:spec_fit, c0:p[indgen(n_lines)*3+0], c1:p[indgen(n_lines)*3+1] , c2:p[indgen(n_lines)*3+2], off:last(p), lab:label}

	RETURN,output
END



;+
;NAME:
;	HIREXSR_FIT_IMAGE
;
;-

PRO hirexsr_fit_image,shot,module,t1,t2,plot=plot,path=path,double=double,break=break,save=save,istart=istart,istop=istop
	IF keyword_set(double) THEN extra_str='dp' ELSE extra_str='sp'
	IF NOT keyword_set(path) THEN path='/home/mlreinke/idl/hirexsr/'
	IF NOT keyword_set(save) THEN save=path+'hirexsr_fit_image_'+num2str(shot,1)+'_'+num2str(module,1)+'_'+extra_str+'.dat'
	hirexsr_load_image,shot,module,image,t
	ilow=ipt(t,t1)
	ihigh=ipt(t,t2)
	image=sum_array(float(image[*,*,ilow:ihigh]),/k)
	x=size(image)
	resid=image*0.0
	IF NOT keyword_set(istart) THEN istart=x[1]-1
	IF NOT keyword_set(istop) THEN istop=x[1]-1


	FOR i=istart,istop DO BEGIN
		IF i MOD 10 EQ 0 THEN print, 'Fitting '+num2str(i,1)+' out of '+num2str(x[1]-1,1)
		IF module EQ 4 THEN out=hirexsr_fit_h(reform(image[i,*]),win=plot,double=double) ELSE out=hirexsr_fit_he(reform(image[i,*]),win=plot,double=double,break=break)
		resid[i,*]=out.residual
		IF i EQ istart THEN BEGIN
			IF keyword_set(double) THEN BEGIN
				nphot=dblarr(x[1],n(out.lab)+1)
				peaks=dblarr(x[1],n(out.lab)+1)
				width=dblarr(x[1],n(out.lab)+1)
				err=dblarr(x[1],n(out.lab)+1)
			ENDIF ELSE BEGIN
				nphot=fltarr(x[1],n(out.lab)+1)
				peaks=fltarr(x[1],n(out.lab)+1)
				width=fltarr(x[1],n(out.lab)+1)
				err=fltarr(x[1],n(out.lab)+1)
			ENDELSE
		ENDIF
		nphot[i,*]=2.0*sqrt(!pi/2.0)*out.c0
		peaks[i,*]=out.c1
		width[i,*]=out.c2
		err[i,*]=out.c2/sqrt(2.0*(nphot-1.0))
	ENDFOR
	label=out.lab
	save,image,resid,nphot,peaks,width,err,label,shot,module,t1,t2,filename=save
END



;+
;NAME:
;	HIREXSR_FIT_SPECTRA
;
;PURPOSE:
;	This is a higher level function which performs the multigaussian fits on a spec or avespec PTRARR
;	
;CALLING SEQUENCE:
;	result=HIREXSR_FIT_SPECTRA(spec)
;
;INPUTS:
;	spec	PTRARR	[chmax,ntime] of pointers to the spectral data (see HIREXSR_AVESPEC_ARR2PTR)
;	line	INT	of the index to the line of interest (0 - w, 1 - x, 2 - z, 3 - lya1, 4 - mo4d, 5 -J)
;
;OPTIONAL INPUTS:
;	limit	FLOAT	value sent to HIREXSR_FIT_###
;	wl	FLTARR	[2] of the wavelength bounds over which to perform fit (see PROCEDURE for DEFAULTS)
;	fits	PTRARR	[chmanx,ntime] of the previous call to HIREXSR_FIT_SPECTRA.  To be used to update fits.
;	fix	INTARR	[chmax,ntime] of 1's (yes) or 0's (no) whether to perform the fit.  This can be used to save time after rebinning a few channels.
;	nback	INT	of the polynomial order for the background [DEFAULT: 0] (currently only linked into ZJK fitting) 
;
;KEYWORD PARAMETERS:
;	/plot sends a win=1 do the HIREXSR_FIT_### program to display plots
;	/verb and /double are sent to HIREXSR_FIT_###

;OUTPUTS:
;	result	PTRARR	[chmax,ntime] where if a spectra exists the then [i,j] pointer is filled with the output of HIREXSR_FIT_###
;			If ispec[0]=-1 THEN [i,j] is set to -1 if no spectra exists.
;	
;OPTIONAL OUTPUTS:
;	label 	STRARR 	[nlines] of the line labels for the corresponding fit coefficients
;	using /resid will fill in the residual after a successful fit in the [*,3] elements of the spec PTRARR
;
;PROCEDURE:
;	The DEFAULTS for the wl and limit are set seperately for each line
;		line = 0	'w'	wl=[3.944,3.9607]	limit=-3.0
;		line = 1	'x'	wl=[3.960,3.9760]	limit=-3.0
;		line = 2	'z'	wl=[3.977,4.0000]	limit=-3.0
;		line = 3	'lya1'	wl=[3.725,4.747]	limit=-3.0
;		line = 4	'mo4d'	wl=[3.725,3.747]	limit=-3.0
;		line = 5	'J'	wl=[3.725,3.747]	limit=-3.0
;		line = 6	'w'(Ca)	wl=[3.170,3.185]	limit=-3.0
;		line = 7	'lya1'	wl=[3.725,4.747] 	limit=-3.0	(high-Te layout)
;		line = 8	'mo4d'	wl=[3.725,3.747]	limit=-3.0	(high-Te layout)
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 7/10
;	8/18/10		M.L. Reinke - added the ability to work on H-like lines (line=3,4)
;	8/27/10		M.L. Reinke - modified the filling of label to do it only once after a successful fit
;	8/30/10		M.L. Reinke - on the failure of a fit for too few photons, COEFS is filled with the baseline offset
;	9/14/10		M.L. Reinke - increased all SNR failures to -3.0 after revising the AVESPEC uncertainty calculation
;	5/13/11		M.L. Reinke - added the nback optional input for use in HIREXSR_FIT_XXX
;	8/25/12		M.L. Reinke - added line=7 and line=8 for the high-Te layout	
;
;-

FUNCTION hirexsr_fit_spectra,spec,line,double=double,label=label,limit=limit,wl=wl,verb=verb,plot=plot,resid=resid,fix=fix,fits=fits,nback=nback
	IF NOT keyword_set(nback) THEN nback=0
	x=size(spec)
	chmax=x[1]
	ntime=x[2]
	CASE line OF 
		0 : BEGIN	;wn3
			wl=[3.944,3.9607]
			snr_limit=-3.0
		END
		1 : BEGIN 	;xy
			wl=[3.960,3.976]
			snr_limit=-3.0
		END
		2 : BEGIN	;zjk
			wl=[3.978,3.9995]
			snr_limit=-3.0
		END
		3 : BEGIN	;lya
			wl=[3.725,3.747]
			snr_limit=-3.0
		END
		4 : BEGIN	;mo4d
			wl=[3.725,3.747]
			snr_limit=-3.0			
		END
		5 : BEGIN	;J
			wl=[3.750,3.780]
			snr_limit=-3.0			
		END
		6 : BEGIN	;wn3 for He-like Ca
			wl=[3.170,3.185]
			snr_limit=-3.0			
                END
		7 : BEGIN	;lya		;lyman on the x3 modules
			wl=[3.725,3.747]
			snr_limit=-3.0
		END
		8 : BEGIN	;mo4d		;Mo on the x3 modules
			wl=[3.725,3.747]
			snr_limit=-3.0			
		END
		9 : BEGIN	;ca lya		;on the x1 module when h-like ar is on x3
			wl=[3.010,3.030]
			snr_limit=-3.0		
		END
		ELSE: RETURN,-1
	ENDCASE	
	IF NOT keyword_set(fits) THEN fits=ptrarr(chmax,ntime,/allocate_heap)		;allow fits to be updated and not created
	IF NOT keyword_set(fix) THEN fix=intarr(chmax,ntime)+1				;look for all ch's and times to try and fit spectra
	label=0

	FOR i=0,chmax-1 DO BEGIN
		IF keyword_set(verb) THEN print, 'Fitting CH '+num2str(i,1)+' OF '+num2str(chmax-1,1)
		FOR j=0,ntime-1 DO BEGIN
			IF fix[i,j] EQ 1 THEN BEGIN		;check if this spectra requires refitting or is new
				ispec=*spec[i,j]
				IF ispec[0] NE -1 THEN BEGIN
					tmp=where(ispec[*,1] GE wl[0] AND ispec[*,1] LE wl[1])
                                        CASE line OF 
						0 : coefs=hirexsr_fit_wn3(ispec[tmp,0],ispec[tmp,1],sig=ispec[tmp,2],double=double,label=ilabel,win=plot,limit=snr_limit,n_lines=n_lines,verb=verb)
						1 : coefs=hirexsr_fit_xy(ispec[tmp,0],ispec[tmp,1],sig=ispec[tmp,2],double=double,label=ilabel,win=plot,limit=snr_limit,n_lines=n_lines,verb=verb)
						2 : coefs=hirexsr_fit_zjk(ispec[tmp,0],ispec[tmp,1],sig=ispec[tmp,2],double=double,label=ilabel,win=plot,limit=snr_limit,n_lines=n_lines,verb=verb,nback=nback)
						3 : coefs=hirexsr_fit_lya(ispec[tmp,0],ispec[tmp,1],sig=ispec[tmp,2],double=double,label=ilabel,win=plot,limit=snr_limit,n_lines=n_lines,verb=verb)
						4 : coefs=hirexsr_fit_lya(ispec[tmp,0],ispec[tmp,1],sig=ispec[tmp,2],double=double,label=ilabel,win=plot,limit=snr_limit,n_lines=n_lines,verb=verb)
						5 : coefs=hirexsr_fit_TJ(ispec[tmp,0],ispec[tmp,1],sig=ispec[tmp,2],double=double,label=ilabel,win=plot,limit=snr_limit,n_lines=n_lines,verb=verb)
						6 : coefs=hirexsr_fit_ca_wn3(ispec[tmp,0],ispec[tmp,1],sig=ispec[tmp,2],double=double,label=ilabel,win=plot,limit=snr_limit,n_lines=n_lines,verb=verb)
						7 : coefs=hirexsr_fit_lya(ispec[tmp,0],ispec[tmp,1],sig=ispec[tmp,2],double=double,label=ilabel,win=plot,limit=snr_limit,n_lines=n_lines,verb=verb)
						8 : coefs=hirexsr_fit_lya(ispec[tmp,0],ispec[tmp,1],sig=ispec[tmp,2],double=double,label=ilabel,win=plot,limit=snr_limit,n_lines=n_lines,verb=verb)
						9 : coefs=hirexsr_fit_ca_lya(ispec[tmp,0],ispec[tmp,1],sig=ispec[tmp,2],double=double,label=ilabel,win=plot,limit=snr_limit,n_lines=n_lines,verb=verb)

					ENDCASE
					IF n(coefs) NE 0 AND keyword_set(resid) THEN BEGIN
						ispec[tmp,3]=ispec[tmp,0]-gaussian_fits(ispec[tmp,1],coefs)
						*spec[i,j]=ispec
					ENDIF
					IF n(coefs) NE 0 AND size(label,/type) NE 7 THEN label=ilabel
					IF coefs[0] EQ -2 THEN BEGIN			;fit failed w/ too few photons - do a background subtraction of the DC offset
						coefs=fltarr(3*n_lines+1+nback)
						coefs[3*n_lines+nback]=median(ispec[tmp,0])	;fill baseline with median
						coefs[indgen(n_lines)*3+2]=1.0			;make width non-zero so you don't get divide by zero errors, leave shift and height = 0
					ENDIF
					*fits[i,j]=coefs
				ENDIF ELSE *fits[i,j]=-1
			ENDIF
		ENDFOR
	ENDFOR
	RETURN,fits
END

;+
;NAME:
;	HIREXSR_CALC_RESID
;
;-

PRO hirexsr_calc_resid,spec,coefs
	x=size(spec)
	chmax=x[1]
	ntime=x[2]
	FOR i=0,chmax-1 DO BEGIN
		FOR j=0,ntime-1 DO BEGIN
			icoefs=*coefs[i,j]
			IF icoefs[0] GT 0 THEN BEGIN		;spectra > 0 means there is a valid fit coefficient for 
				ispec=*spec[i,j]
				ispec[*,3]=ispec[*,0]-gaussian_fits(ispec[*,1],icoefs)
				*spec[i,j]=ispec
			ENDIF
		ENDFOR
	ENDFOR
END
