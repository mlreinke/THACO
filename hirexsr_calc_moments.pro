;+
;NAME:
;	HIREXSR_LINE_SUBTRACT
;
;PURPOSE:
;	This function removes neighboring lines using coeficients from multigaussian fits leaving the line
;	of interest for use in moment generation
;
;CALLING SEQUENCE:
;	result=HIREXSR_LINE_SUBTRACT(lam,specbr,coefs,labels,line)
;
;INPUTS:
;	LAM	FLTARR	[n] of the wavelength values
;	SPECBR	FLTARR	[n] of the spectral brightness at each value of lam
;	COEFS	FLTARR	[3*n_lines+nb] of the fit coefficients (nb specifies the order of the background)
;	LABELS	STRARR	[n_lines] of the labels for each line
;	LINE	INT	of the line indentification (w:0, x:1, z:2, lya1:3, mo4d: 4)
;
;OUTPUTS:
;	result	FLTARR	[n] of the spectral brightness at each lam with various lines removed
;
;OPTIONAL OUTPUTS:
;	back	FLTARR	[n] of the background subtraction done at each lam
;
;PROCEDURE:
;	LINE=0  - wn3 and wn4 are removed
;	LINE=1	- n, st and y are removed
;	LINE=2	- k and a are removed
;	LINE=3 	- lyas1, lyas2,lya2 are removed
;	LINE=4	- lya2, lyas3, lyas4 are removed
;	LINE=5
;	LINE=6	- wn3 and wn4 are removed
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 7/10
;	8/18/10		M.L. Reinke - added the subtraction setup for line=3,4
;	9/27/10		M.L. Reinke - added the back optional output
;	12/17/10	M.L. Reinke - added the subtraction for the n=5 on w and the j on z
;	12/20/10	M.L. Reinke - added an exisistence check on the n=5, and j lines so it won't crash
;				      if it doesn't find the fits for these lines
;	5/13/11		M.L. Reinke - added the ability to subtract the  n=0,1 or 2 order polynomial background 
;	8/9/11		M.L. Reinke - added the ability to use line=5 (TJ) in moment analysis
;	8/25/12		M.L. Reinke - added the line=7,8 configurations  high-Te layout
;	12/10/12	M.L. Reinke - added the lyas3 to the subtraction for line=3,7 (lya1)
;	3/28/14		M.L. Reinke - added line=9 for H-like Ca
;-

FUNCTION hirexsr_line_subtract,lam,specbr,coefs,labels,line,back=back
	CASE line OF
		0 : BEGIN
			i=where(labels EQ 'wn3')
			j=where(labels EQ 'wn4')
			k=where(labels EQ 'wn5')
			IF k[0] EQ -1 THEN subcoefs=[coefs[3*i:3*(i+1)-1],coefs[3*j:3*(j+1)-1]] ELSE $
				subcoefs=[coefs[3*i:3*(i+1)-1],coefs[3*j:3*(j+1)-1],coefs[3*k:3*(k+1)-1]]
		END
		1 : BEGIN
			i=where(labels EQ 'n')
			j=where(labels EQ 'st')
			k=where(labels EQ 'y')
			subcoefs=[coefs[3*i:3*(i+1)-1],coefs[3*j:3*(j+1)-1],coefs[3*k:3*(k+1)-1]]
		END
		2 : BEGIN
			i=where(labels EQ 'k')
			j=where(labels EQ 'a')
			k=where(labels EQ 'j')
			IF k[0] EQ -1 THEN subcoefs=[coefs[3*i:3*(i+1)-1],coefs[3*j:3*(j+1)-1]]
			subcoefs=[coefs[3*i:3*(i+1)-1],coefs[3*j:3*(j+1)-1],coefs[3*k:3*(k+1)-1]]
		END
		3 : BEGIN
			i=where(labels EQ 'lyas1')
			j=where(labels EQ 'lyas2')
			k=where(labels EQ 'lya2')
			l=where(labels EQ 'lyas3')
			subcoefs=[coefs[3*i:3*(i+1)-1],coefs[3*j:3*(j+1)-1],coefs[3*k:3*(k+1)-1],coefs[3*l:3*(l+1)-1]]
		END
		4 : BEGIN
			i=where(labels EQ 'lya2')
			j=where(labels EQ 'lyas3')
			k=where(labels EQ 'lyas4')
			subcoefs=[coefs[3*i:3*(i+1)-1],coefs[3*j:3*(j+1)-1],coefs[3*k:3*(k+1)-1]]
		END
		5 : BEGIN
			i=where(labels EQ 'R')
			j=where(labels EQ 'A')
			subcoefs=[coefs[3*i:3*(i+1)-1],coefs[3*j:3*(j+1)-1]]
		END
		6 : BEGIN
			i=where(labels EQ 'wn3')
			j=where(labels EQ 'wn4')
			subcoefs=[coefs[3*i:3*(i+1)-1],coefs[3*j:3*(j+1)-1]]
                 END
		7 : BEGIN
			i=where(labels EQ 'lyas1')
			j=where(labels EQ 'lyas2')
			k=where(labels EQ 'lya2')
			l=where(labels EQ 'lyas3')
			subcoefs=[coefs[3*i:3*(i+1)-1],coefs[3*j:3*(j+1)-1],coefs[3*k:3*(k+1)-1],coefs[3*l:3*(l+1)-1]]
		END
		8 : BEGIN
			i=where(labels EQ 'lya2')
			j=where(labels EQ 'lyas3')
			k=where(labels EQ 'lyas4')
			subcoefs=[coefs[3*i:3*(i+1)-1],coefs[3*j:3*(j+1)-1],coefs[3*k:3*(k+1)-1]]
		END
		9 : BEGIN
			i=where(labels EQ 'lya2')
			subcoefs=[coefs[3*i:3*(i+1)-1]]
		END


	ENDCASE
	
	ncoefs=n(coefs)+1
	CASE ncoefs MOD 3 OF 
		0 : subcoefs=[subcoefs,coefs[ncoefs-3],coefs[ncoefs-2],coefs[ncoefs-1]]		;linear+sqrt baseline plus the peak(s) of interest	
		1 : subcoefs=[subcoefs,coefs[ncoefs-1]]						;DC baseline plus the peak(s) of interest	
		2 : subcoefs=[subcoefs,coefs[ncoefs-2],coefs[ncoefs-1]]				;linear baseline plus the peak(s) of interest
	ENDCASE		

	back=gaussian_fits(lam,subcoefs)
	subspec=specbr-back
	RETURN,subspec
END

;+
;NAME:
;	HIREXSR_CALC_MOMENTS
;
;PURPOSE:
;	This is a lower level function with calculates the spectral moments and their uncertainties from
;	the spectral brightness and wavelength data.  This should be run on the output of HIREXSR_LINE_SUBTRACT.
;
;CALLING SEQUENCE:
;	HIREXSR_CALC_MOMENTS,specbr,lam,lam_o,sig,mom,err,scale,bfrac
;
;INPUTS:
;	specbr	FLTARR	[n] of spectral brightness data 
;	lam	FLTARR	[n] of the wavelength values
;	lam_o	FLOAT	wavelength about which to take the 1st and 2nd moments (rest wavelength)
;	sig	FLTARR	[n] of the uncertainties in the spectral brightness
;
;KEYWORD PARAMETERS:
;	/double 	is linked through to INT_TABULATED and will force the output to be DBLARR [3] instead of FLTARR [3]
;	/trap_int	uses the function TRAP_INT to perform the numerical integration
;	/gauss		uses the GAUSSFIT function (nterms=4) to fit line and calculate the moments
;	/mpfit		uses MPFIT (constrained, a[0] and a[4] > 0)  to fit line and calculate the moments
;
;OUTPUTS:
;	mom	FLTARR	[3] of the 0th, 1st and 2nd moments
;	err	FLTARR	[3] of the uncertainty in those moments
;	pmom	FLTARR	[3] of the #phot, mu and w^2 of the distribution
;	perr	FLTARR	[3] of the uncertainty in those values (via Hutch)
;	scale	FLT	of the constant between MOM[0] and PMOM[0] w/o background removed
;	bfrac	FLT	fraction of the # of background photons to background+signal photons		
;
;OPTIONAL OUTPUTS:
;	afit	FLTARR 	[4] of the single guassian + DC offset fitting parameters if /gauss or /mpfit is used
;	asig	FLTARR 	[4] of the uncertainties in the afit parameters
;
;PROCEDURE:
;	The methods for calculating the moments and errors are detailed in the XCIS equations paper.
;	
;	The various mathematical techniques have their benifits and drawbacks.  In an exact sense, the numerical moment
;	from INT_TABULATED or TRAP_INT is correct moment since the instrument function and line-integral make the
;	spectral line non-gaussian.  But, at low signal-to-background, the noise associated with taking moments
;	dramatically increases making the data useless.  Thus the best thing to do would be do estimate using 
;	HIREXSR_MOMENTS_COMPARE to see for a given width and number of #photons where the cutoff is.  Then feedback on this
;	in HIREXSR_SPEC2MOMENTS.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 7/10
;	9/14/10		M.L. Reinke - finally got the better moment/error calculation equations impliments
;	12/7/10		M.L. Reinke - large overhaul to include multiple techniques to find the mom, err, pmom and perr
;	12/8/10		M.L. Reinke - fixed the INT_TAB and TRAP_INT methods of calculating moment error.  1st moment error looks reasonable now
;				      added, the afit and asig optional output
;	12/9/10		M.L. Reinke - fixed a bug when using GAUSSFIT or MPFIT to find the moments that caused a crash when sig=0
;	1/5/11		M.L. Reinke - fixed a case where a crash would occur if median(sig)=0.  This case is not treated like total(specbr)=0,
;				      and thrown out as a bad frame
;	5/13/11		M.L. Reinke - adjusted the MPFIT and GAUSSFIT options to remove the baseline (should be removed by HIREXSR_LINE_SUBTRACT)
;-

PRO hirexsr_calc_moments,specbr,lam,lam_o,sig,mom,err,pmom,perr,scale,bfrac,afit=afit,asig=asig,double=double,back=back,debug=debug,trap=trap,gauss=gauss,mpfit=mpfit
	fitcase=0
	IF keyword_set(trap) THEN fitcase=1
	IF keyword_set(gauss) THEN fitcase=2
	IF keyword_set(mpfit) THEN fitcase=3
	IF keyword_set(double) THEN BEGIN
		mom=dblarr(3) 
		err=dblarr(3)
		pmom=dblarr(3)
		perr=dblarr(3)
	ENDIF ELSE BEGIN
		mom=fltarr(3)
		err=fltarr(3)
		pmom=fltarr(3)
		perr=fltarr(3)

	ENDELSE
	IF total(specbr) EQ 0 OR median(sig) EQ 0 THEN RETURN		;if no signal or median of the sig=0 (disruption frame?) then RETURN empty moments
	del=lam-lam_o
	order=sort(lam)
	tmp=uniq(lam[order])				;remove possible "doubles" in wavelength
	n=n(tmp)+1
	IF NOT keyword_set(back) THEN back=fltarr(n)	;if no background subtraction given set = 0

	IF NOT keyword_set(wseed) THEN wseed=4.6e-4	;roughly 500 eV 
	CASE fitcase OF
		2 : BEGIN
			afit=[max(specbr),lam_o,wseed]
			tmpzero=where(sig EQ 0)
			IF tmpzero[0] NE -1 THEN sig[tmpzero]=median(sig)/2.0
			out=gaussfit(lam,specbr,afit,nterms=n(afit)+1,measure_errors=sig,sigma=asig)
		END
		3 : BEGIN
			estimate=[abs(max(specbr)),lam_o,wseed]		
			parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(estimate))
			parinfo[*].value=estimate
			parinfo[0].limited[0]=1
			parinfo[0].limits[0]=0.0
			parinfo[2].limited[0]=1
			parinfo[2].limits[0]=0.5*(last(lam)-lam[0])/(n(lam)+1.0)	;minimum ~1/2 pixel width (Ti_Ar ~ 100 eV)
			tmpzero=where(sig EQ 0)
			IF tmpzero[0] NE -1 THEN sig[tmpzero]=median(sig)/2.0
			afit=mpfitfun('gaussian_fits',lam,specbr,sig, estimate,parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30,perror=asig)
			IF afit[0] EQ 0 THEN afit[0]=abs(median(specbr))		;rough fix for channels that fail with zero brightness, crashing error analysis in INVERT
		END
		ELSE : 
	ENDCASE

	x=lam[order[tmp]]					;sort the the wavelength
	
	;zeroth moment				
	y=specbr[order[tmp]]					;sort the pectral brightness and uncertainty values
	ysig=sig[order[tmp]]
	CASE fitcase OF
		0 : BEGIN
			mom[0]=int_tabulated(x,y,double=double)
			FOR i=0,n-2 do err[0]+=(0.5*(x[i+1]-x[i]))^2*ysig[i]^2+(0.5*(x[i+1]-x[i]))^2*ysig[i+1]^2
		END
		1 : BEGIN
			mom[0]=trap_int(x,y)
			FOR i=0,n-2 do err[0]+=(0.5*(x[i+1]-x[i]))^2*ysig[i]^2+(0.5*(x[i+1]-x[i]))^2*ysig[i+1]^2
		END
		2 : BEGIN
			mom[0]=sqrt(2.0*!pi)*afit[0]*afit[2]
			err[0]=2.0*!pi*(afit[2]^2*asig[0]^2+afit[0]^2*asig[2]^2)
		END
		3 :BEGIN
			mom[0]=sqrt(2.0*!pi)*afit[0]*afit[2]
			err[0]=2.0*!pi*(afit[2]^2*asig[0]^2+afit[0]^2*asig[2]^2)
		END
	ENDCASE 

	;first moment	
	y=specbr[order[tmp]]*del[order[tmp]]
	ysig=sig[order[tmp]]*del[order[tmp]]
	CASE fitcase OF

		0 : BEGIN
			mom[1]=int_tabulated(x,y,double=double)
			FOR i=0,n-2 do err[1]+=(0.5*(x[i+1]-x[i]))^2*ysig[i]^2+(0.5*(x[i+1]-x[i]))^2*ysig[i+1]^2
		END
		1 : BEGIN
			mom[1]=trap_int(x,y)
			FOR i=0,n-2 do err[1]+=(0.5*(x[i+1]-x[i]))^2*ysig[i]^2+(0.5*(x[i+1]-x[i]))^2*ysig[i+1]^2
		END
		2 : BEGIN
			mom[1]=mom[0]*(afit[1]-lam_o)
			err[1]=mom[0]^2*asig[1]^2+err[0]^2*(afit[1]-lam_o)^2
		END
		 
		3 :BEGIN
			mom[1]=mom[0]*(afit[1]-lam_o)
			err[1]=mom[0]^2*asig[1]^2+err[0]^2*(afit[1]-lam_o)^2
		END
	ENDCASE

	;second moment
	y=specbr[order[tmp]]*del[order[tmp]]^2
	ysig=sig[order[tmp]]*del[order[tmp]]^2
	CASE fitcase OF

		0 : BEGIN
			mom[2]=int_tabulated(x,y,double=double)
			FOR i=0,n-2 do err[2]+=(0.5*(x[i+1]-x[i]))^2*ysig[i]^2+(0.5*(x[i+1]-x[i]))^2*ysig[i+1]^2
		END
		1 : BEGIN
			mom[2]=trap_int(x,y)
			FOR i=0,n-2 do err[2]+=(0.5*(x[i+1]-x[i]))^2*ysig[i]^2+(0.5*(x[i+1]-x[i]))^2*ysig[i+1]^2
		END
		2 : BEGIN
			mom[2]=mom[0]*(afit[2]^2+(lam_o-afit[1])^2)
			err[2]=err[0]^2*(afit[2]^2+(lam_o-afit[1])^2)^2+(2.0*mom[0]*afit[2])^2*asig[2]^2+(2.0*mom[0]*(afit[1]-lam_o))^2*asig[1]^2
		END
		3 :BEGIN
			mom[2]=mom[0]*(afit[2]^2+(lam_o-afit[1])^2)
			err[2]=err[0]^2*(afit[2]^2+(lam_o-afit[1])^2)^2+(2.0*mom[0]*afit[2])^2*asig[2]^2+(2.0*mom[0]*(afit[1]-lam_o))^2*asig[1]^2
		END
	ENDCASE
	err=sqrt(err)	
	
	tmp=uniq(lam[order])
	x=lam[order[tmp]]
	y=specbr[order[tmp]]	
	yback=back[order[tmp]]
	CASE fitcase OF
		0 : BEGIN
			mtot=int_tabulated(x,y+yback,double=double)
			bfrac=int_tabulated(x,yback,double=double)/mtot		;# of background photons to back+signal
		END
		1 :BEGIN
			mtot=trap_int(x,y+yback)
			bfrac=trap_int(x,yback)/mtot
		END
		2 : BEGIN
			btot=int_tabulated(x,yback,double=double)
			mtot=sqrt(2.0*!pi)*afit[0]*afit[2]+btot
			bfrac=btot/mtot
		END
		3 : BEGIN
			btot=int_tabulated(x,yback,double=double)
			mtot=sqrt(2.0*!pi)*afit[0]*afit[2]+btot
			bfrac=btot/mtot
		END
		ELSE :
	END
	nphot=((specbr+back)/sig)^2				;this recovers the number of photons from units of spectral brightness
	fine=where(sig NE 0)					;avoid divide by zero errors
	scale=mtot/total(nphot[fine])				;M0=scale*Nph	
	IF fitcase EQ 0 OR fitcase EQ 1 THEN BEGIN			
		pmom[0]=mom[0]/scale					
		pmom[1]=mom[1]/mom[0]+lam_o				;mu
		pmom[2]=sqrt(mom[2]/mom[0]-(mom[1]/mom[0])^2)		;w
		perr[0]=sqrt(pmom[0])
		perr[1]=sqrt(pmom[2]^2/pmom[0])				;uncertainty in N, mu, w
		perr[2]=sqrt(pmom[2]^2/(2.0*(pmom[0]-1)))
	ENDIF
	IF fitcase EQ 2 OR fitcase EQ 3 THEN BEGIN
		pmom[0]=mom[0]/scale		
		pmom[1]=afit[1]
		pmom[2]=afit[2]
		perr[0]=sqrt(pmom[0])
		perr[1]=asig[1]
		perr[2]=asig[2]
	ENDIF
	IF keyword_set(debug) THEN stop
END

;+
;NAME:
;	HIREXSR_SPEC2MOMENTS
;
;PURPOSE:
;	This is a higher level function with takes in a spec PTRARR and a coefs PTRARR  and generates the
;	temporal and spatial moment profiles.
;
;CALLING SEQUENCE:
;	result=HIREXSR_SPEC2MOMENTS(spec,coefs,labels,lam_o,line,chmax,time)
;
;INPUTS:
;	spec	PTRARR	[chmax,ntime] of pointers to the spectral data (see HIREXSR_AVESPEC_ARR2PTR)
;	coefs	PTRARR	[chmax,ntime] of pointers to the fit coefficient data (see HIREXSR_COEFS_ARR2PTR)
;	labels	STRARR	[nlines] of the line labels
;	lam_o	FLOAT	of the rest wavelength of the line of interest
;	line	INT	of the index of the line of interest (0 - w, 1 - x, 2 - z, 3 - lya1, 4 - mo4d, 5 - J)
;	
;OPTIONAL INPUTS:
;	dlam	FLOAT	of the wavelength interval over which the moment is calculated (lam_o-dlam,lam_o+dlam)
;	mom	FLTARR	[chmax,ntime,3,4] of the moments from a previous call to HIREXSR_SPEC2MOMENTS
;	fix	INTARR	[chmax,ntime] of 1's (yes) and 0's (no) to update the moments
;
;KEYWORD_PARAMETERS:
;	/double calculates the moments in double preceision and is sent to HIREXSR_CALC_MOMENTS
;	/fitmom	forces the moment calculation to use MPFIT, usuafully used only when count rates are low
;
;OUTPUTS:
;	result	FLTARR	[chmax,ntime,3,4] where [*,*,i,0] are the ith moments and [*,*,i,1] are the uncertainty in those moments
;			[*,*,i,2] are the "profile" moments and [*,*,i,3] are their errors. If there is no spectral data for a given (ch,time)
;			THEN [ch,time,*,*] is set to -1
;
;OPTIONAL OUTPUTS:
;	bfrac	FLTARR	[chmax,ntime] of the fraction of background to signal+background photons for use in noise analysis
;	scale	FLTARR	[chmax,ntime] of the constant between MOM[0] and PMOM[0] w/o background removed
;	fitcase	INTARR	[chmax,ntime] indicating the technique used to calculate the moments (0 - INT_TAB, 3 - MPFIT) 
;	
;PROCEDURE:
;	For each (ch,time) point, a check is made to see if there is a valid spectral array (ispec[0] NE -1).
;	If there is a valid fit (n(icoefs) GT 0) then HIREXSR_LINE_SUBTRACT is called.  Then HIREXSR_CALC_MOMENTS
;	is called to calculate the moment over the (lam_o-dlam, lam_o+dlam) interval.
;
;MODFICATION HISTORY:
;	Written by:	M.L.Reinke - 7/10
;	9/12/10		M.L.Reinke - added the mom,fix optional inputs so that the function can be used to updated a selection of moments.
;				     fixed a bug that would allow COEFS < 0 (fails) to be sent to moment calculation routines
;	9/14/10		M.L.Reinke - modified mom
;	8/25/12		M.L. Reinke - added the line=7,8 configurations  high-Te layout
;	12/10/12	M.L. Reinke - changed default line=7 dlam=3.6 since it is for higher-T plasmas
;	3/28/14		M.L. Reinke - line=9
;	7/27/16		M.L. Reinke - added /fitmom to force using MPFIT to compute moments
;-

FUNCTION hirexsr_spec2moments,spec,coefs,labels,lam_o,line,dlam=dlam,double=double,mom=mom,bfrac=bfrac,scale=scale,fitcase=fitcase,fix=fix,fitmom=fitmom

	;get maximum size of moments
	x=size(spec)
	chmax=x[1]
	ntime=x[2]

	CASE line OF 
		0 : BEGIN
			IF NOT keyword_set(dlam) THEN dlam=3.1
		END
		1 : BEGIN
			IF NOT keyword_set(dlam) THEN dlam=3.1
		END
		2 : BEGIN
			IF NOT keyword_set(dlam) THEN dlam=3.1	
		END
		3 : BEGIN
			IF NOT keyword_set(dlam) THEN dlam=3.0
		END
		4 : BEGIN
			IF NOT keyword_set(dlam) THEN dlam=3.0
		END
		5 : BEGIN
			IF NOT keyword_set(dlam) THEN dlam=3.0 
		END
		6 : BEGIN
			IF NOT keyword_set(dlam) THEN dlam=3.0
                END
		7 : BEGIN
			IF NOT keyword_set(dlam) THEN dlam=3.6
		END
		8 : BEGIN
			IF NOT keyword_set(dlam) THEN dlam=3.0
		END
		9 : BEGIN
			IF NOT keyword_set(dlam) THEN dlam=3.0
		END	
	ENDCASE

	IF NOT keyword_set(mom) THEN BEGIN
		IF keyword_set(double) THEN BEGIN 
			mom=dblarr(chmax,ntime,3,4)
			bfrac=dblarr(chmax,ntime)
			scale=dblarr(chmax,ntime)
		ENDIF ELSE BEGIN
			mom=fltarr(chmax,ntime,3,4)
			bfrac=dblarr(chmax,ntime)
			scale=dblarr(chmax,ntime)		
		ENDELSE
		fitcase=intarr(chmax,ntime)-1
	ENDIF
	IF NOT keyword_set(fix) THEN fix=intarr(chmax,ntime)+1

	FOR i=0,chmax-1 DO BEGIN
		FOR j=0,ntime-1 DO BEGIN
			IF fix[i,j] THEN BEGIN
				ispec=*spec[i,j]
				icoefs=*coefs[i,j]
				IF ispec[0] NE -1 THEN BEGIN
					tmp=where(ispec[*,1] GE lam_o-dlam/1.0e3 AND ispec[*,1] LE lam_o+dlam/1.0e3)
					lam=ispec[tmp,1]
					specbr=ispec[tmp,0]
					sig=ispec[tmp,2]
					IF icoefs[0] GE 0 THEN specbr=hirexsr_line_subtract(lam,specbr,icoefs,labels,line,back=back)
					hirexsr_calc_moments,specbr,lam,lam_o,sig,imom,ierr,ipmom,iperr,iscale,ibfrac,double=double,back=back,debug=debug
					fitcase[i,j]=0
					IF 1.0/ibfrac LE 4.0 OR keyword_set(fitmom) THEN BEGIN	       ;if low signal to background then use MPFIT 
						hirexsr_calc_moments,specbr,lam,lam_o,sig,imom,ierr,ipmom,iperr,iscale,ibfrac,double=double,back=back,debug=debug,/mpfit
						fitcase[i,j]=3
					ENDIF
					mom[i,j,*,0]=imom
					mom[i,j,*,1]=ierr
					mom[i,j,*,2]=ipmom
					mom[i,j,*,3]=iperr
					bfrac[i,j]=ibfrac
					scale[i,j]=iscale
				ENDIF ELSE mom[i,j,*,*]=-1.0	;fill w/ -1 if no spectral data available
			ENDIF	
		ENDFOR
	ENDFOR
	RETURN,mom
END

;+
;NAME:
;	HIREXSR_MOMENT_CONV
;
;PURPOSE:
;	This function calculates the spectral moments for a range of wavelength ranges in order to show
;	that the moments have converged and that error hasn't begun to creap in.
;
;CALLING SEQUENCE
;	result=HIREXSR_MOMENT_CONV(ispec,icoefs,labels,lam_o,line,dlam)
;
;INPUTS:
;	ispec	FLTARR	[n,4] element of a spec or avespec PTRARR that holds the spectral information
;	icoefs	FLTARR	[3*nlines+1] element of a coef PTRARR that holds the fit coefficients
;	labels	STRARR	[nlines] of the line labels
;	lam_o	FLOAT	the rest wavelength
;	line	INT	of the line index (0 - w, 1 - x, 2 - z, 3 - lya1, 4 - mo4d, 6 - w (Ca))
;	dlam	FLTARR	[m] of the various ranges (lam_o-dlam to lam_o+dlam) to take the moments 
;			    [units are 10^3 of those in lam_o and ispec[*,1]]
;
;KEYWORD PARAMETERS:
;	/double calculates the moments in double precision and is passed to HIREXSR_CALC_MOMENTS
;
;OUTPUTS:
;	result	FLTARR	[m,3,4] where [*,*,0] are the moments and [*,*,1] are the errors and [*,*,2],[*,*,3] are the profile
;			moments and errors
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 7/10
;	9/14/20		M.L. Reinke - added the use of profile moments
;
;-

FUNCTION hirexsr_moment_conv,ispec,icoefs,labels,lam_o,line,dlam,double=double

	ndlam=n(dlam)+1
	IF keyword_set(double) THEN mom=dblarr(ndlam,3,4) ELSE mom=fltarr(ndlam,3,4)
	FOR i=0,ndlam-1 DO BEGIN
		tmp=where(ispec[*,1] GE lam_o-dlam[i]/1.0e3 AND ispec[*,1] LE lam_o+dlam[i]/1.0e3)	
		lam=ispec[tmp,1]
		specbr=ispec[tmp,0]
		sig=ispec[tmp,2]
		specbr=hirexsr_line_subtract(lam,specbr,icoefs,labels,line,back=back)
		hirexsr_calc_moments,specbr,lam,lam_o,sig,imom,ierr,ipmom,iperr,iscale,ibfrac,double=double,back=back
		mom[i,*,0]=imom
		mom[i,*,1]=ierr
		mom[i,*,2]=ipmom
		mom[i,*,3]=iperr
	ENDFOR
	RETURN,mom
END

;+
;NAME:
;	HIREXSR_CALC_RHOTANG
;
;PURPOSE:
;	This function calculates the rho of closests approach for an array of POS vectors
;	which can be used for plotting line-integrated moments/fit coefficients.  This uses
;	EFIT_RZ2RHO which outputs "rho" as normalized poloidal flux.
;
;CALLING SEQUENCE:
;	result=HIREXSR_CALC_RHOTANG(shot,pos,tpos,tau)
;
;OPTIONAL INPUTS:
;	tree:		STRING	of the EFIT tree to use for calculating tangency radii [DEFAULT: 'analysis']
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 9/24/10
;	6-13-11:	ML Reinke - added the tree optional input for use with EFIT_RZ2XXX codes in GENPOS_POS2RMIDTANG
;
;-

FUNCTION hirexsr_calc_rhotang,shot,pos,tpos,tau,ptr=ptr,tree=tree
	ntime=n(tau)+1
	x=size(pos)
	nch=x[2]
	
	IF keyword_set(ptr) THEN rhotang=ptrarr(ntime,/allocate) ELSE rhotang=fltarr(nch,ntime)-1.0
	FOR i=0,ntime-1 DO BEGIN
		IF tau[i] NE -1 THEN BEGIN
			pindex=last(where(tpos LE tau[i]))
			ipos=pos[*,*,pindex]			
			tmp=where(ipos[0,*] NE -1)
			ipos=ipos[*,tmp]
			genpos_pos_reform,ipos,[0.44,1.0,-0.6,0.6]
			irho=genpos_pos2rmidtang(ipos,shot,tau[i],/psin,tree=tree)
			IF keyword_set(ptr) THEN *rhotang[i]=irho ELSE rhotang[0:n(tmp),i]=irho
		ENDIF ELSE IF keyword_set(ptr) THEN *rhotang[i]=-1
	ENDFOR
	RETURN,rhotang
END

;+
;NAME:
;	HIREXSR_MOMENTS_COMPARE
;
;PURPOSE:
;	This procedure is used to generate noisy single gaussians to emulate experimental spectral lines so
;	that the different technqiues of HIREXSR_CALC_MOMENTS can be compared.
;
;MODIFICATION HISTORY:
;	Written by	M.L. Reinke - 12/7/2010
;
;-

PRO hirexsr_moments_compare,fitplot=fitplot,back=back,istop=istop,dlam=dlam,ti=ti,v=v,z=z,seed=seed,debug=debug
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0

	lam_o=3.994
	c=3.0e8  			;speed of light
	e=1.60e-19			;conversion for eV -> J
	mconv=1.66e-27			;conversion for amu -> kg
	IF NOT keyword_set(z) THEN z=18
	mass=read_atomic_mass(z)	;mass of Ar
	
	IF NOT keyword_set(back) THEN back=0.0
	nterms=4
	IF NOT keyword_set(dlam) THEN dlam=4.0

	;number of photons scan
	nspec=25
	IF NOT keyword_set(v) THEN v=0.0	;[km/s]
	IF NOT keyword_set(ti) THEN ti=1.0e3 	;[eV]
	w=sqrt(lam_o^2*e*ti/(mass*mconv*c^2))
	mu=lam_o-lam_o/c*v*1.0e3
	x=double(make(lam_o-0.004,lam_o+0.004,nspec))
	x_o=lam_o					;unshifted line_center
	nphot=float(int(10.0^make(1.6,3.5,50)))
	num=n(nphot)+1
	counts=fltarr(num,8)
	center=fltarr(num,8)
	width=fltarr(num,8)
	xplot=make(min(x),max(x),100)
	FOR i=0,n(nphot) DO BEGIN
		a=[10.0,mu,w]
		y=a[0]*exp(-(x-a[1])^2/(2.0*a[2]^2))
		scale=nphot[i]/total(y)
		y*=scale
		ynoise=fltarr(nspec)
		FOR j=0,nspec-1 DO ynoise[j]=randomn(seed,1,poisson=y[j])			;generate spectrum
		IF back NE 0 THEN ynoise+=randomn(seed,nspec,poisson=back)			;add background	
		ysig=sqrt(ynoise)								;calc photon statistics
		ynoise-=fltarr(nspec)+back							;remove DC offset
		tmp=where(ysig LT 1.0)
		IF tmp[0] NE -1 THEN ysig[tmp]=1.0
		
		tmp=where(x GE lam_o-dlam*1.0e-3 AND x LE lam_o+dlam*1.0e-3)
		hirexsr_calc_moments,ynoise[tmp],x[tmp],x_o,ysig[tmp],mom,err,pmom,perr,/gauss,back=fltarr(n(tmp)+1)+back,afit=afit_gf
		counts[i,0]=pmom[0]
		counts[i,1]=perr[0]
		center[i,0]=pmom[1]
		center[i,1]=perr[1]
		width[i,0]=pmom[2]
		width[i,1]=perr[2]		
		
		tmp=where(x GE lam_o-dlam*1.0e-3 AND x LE lam_o+dlam*1.0e-3)
		hirexsr_calc_moments,ynoise[tmp],x[tmp],x_o,ysig[tmp],mom,err,pmom,perr,back=fltarr(n(tmp)+1)+back
		counts[i,2]=pmom[0]
		counts[i,3]=perr[0]
		center[i,2]=pmom[1]
		center[i,3]=perr[1]
		width[i,2]=pmom[2]
		width[i,3]=perr[2]

		tmp=where(x GE lam_o-dlam*1.0e-3 AND x LE lam_o+dlam*1.0e-3)
		hirexsr_calc_moments,ynoise[tmp],x[tmp],x_o,ysig[tmp],mom,err,pmom,perr,/trap,back=fltarr(n(tmp)+1)+back
		counts[i,4]=pmom[0]
		counts[i,5]=perr[0]
		center[i,4]=pmom[1]
		center[i,5]=perr[1]
		width[i,4]=pmom[2]
		width[i,5]=perr[2]

		estimate=[max(ynoise),x_o,w/2.0,0.0]
		parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(estimate))
		parinfo[*].value=estimate
		parinfo[0].limited[0]=1
		parinfo[0].limits[0]=0.0	
		parinfo[1].limited[0]=1
		parinfo[1].limited[1]=1
		parinfo[1].limits[0]=lam_o*(1.0-100.0e3/c)
		parinfo[1].limits[1]=lam_o*(1.0+100.0e3/c)
		parinfo[2].limited[0]=1
		parinfo[2].limited[1]=1
		parinfo[2].limits[0]=sqrt(lam_o^2*e*100.0/(mass*mconv*c^2))
		parinfo[2].limits[1]=sqrt(lam_o^2*e*5000.0/(mass*mconv*c^2))
		parinfo[3].limited[0]=1
		parinfo[3].limits[0]=0.0

		;afit=mpfitfun('gaussian_fits', x,ynoise,ysig, estimate,/quiet, parinfo=parinfo,status=status,niter=niter,xtol=1.0e-30,perror=asig)
		;counts[i,6]=total(gaussian_fits(x,afit))
		;counts[i,7]=sqrt(counts[i,6])
		;center[i,6]=afit[1]
		;center[i,7]=asig[1]
		;width[i,6]=afit[2]
		;width[i,7]=asig[2]

		tmp=where(x GE lam_o-dlam*1.0e-3 AND x LE lam_o+dlam*1.0e-3)
		hirexsr_calc_moments,ynoise[tmp],x[tmp],x_o,ysig[tmp],mom,err,pmom,perr,/mpfit,back=fltarr(n(tmp)+1)+back,afit=afit_mp
		counts[i,6]=pmom[0]
		counts[i,7]=perr[0]
		center[i,6]=pmom[1]
		center[i,7]=perr[1]
		width[i,6]=pmom[2]
		width[i,7]=perr[2]

		IF keyword_set(fitplot) THEN BEGIN
			yr=[min(ynoise-ysig)*1.05 < 0,max(ynoise+ysig)+1.05]
			plot,x,ynoise,psym=8,xtit=n2g('lambda')+' [Ang]',ytit='# Counts',tit='Total # = '+num2str(nphot[i],dp=1),/xsty,/ysty,yr=yr,/nodata
			yplot=a[0]*exp(-(xplot-a[1])^2/(2.0*a[2]^2))
			oploterror,x,ynoise,ysig,psym=7
			oplot,xplot,yplot*scale,linestyle=2.0
			yplot=afit_gf[0]*exp(-(xplot-afit_gf[1])^2/(2.0*afit_gf[2]^2))
			;IF nterms EQ 4 THEN yplot+=fltarr(100)+afit_gf[3]
			oplot,xplot,yplot,color=200
			yplot=afit_mp[0]*exp(-(xplot-afit_mp[1])^2/(2.0*afit_mp[2]^2))
			;IF nterms EQ 4 THEN yplot+=fltarr(100)+afit_mp[3]
			oplot,xplot,yplot,color=120

			oplot,x_o*[1,1],[-max(yr),max(yr)],linestyle=1.0
			oplot,[min(x),max(x)],[0,0],linestyle=1.0
		ENDIF
		IF keyword_set(istop) THEN stop
	ENDFOR
	IF NOT keyword_set(back) THEN tit='No Background, d'+n2g('lambda')+'/dpix = '+num2str((x[1]-x[0])*1.0e3,dp=3)+' mAng' ELSE $
			tit='Back = '+num2str(int(back),1)+', d'+n2g('lambda')+'/dpix = '+num2str((x[1]-x[0])*1.0e3,dp=3)+' mAng'

	!p.multi=[0,0,3]
	ls=1.0
	IF keyword_set(ps) THEN BEGIN
			d_old=!d
			xsize=7.0
			ysize=10.0
			device, xsize=xsize, ysize=ysize, /inches
			ls=0.8
	ENDIF ELSE openwin,0,ysize=1000

	;counts
	xr=[min(nphot)/2,2.0*max(nphot)]
	yr=[1,max(counts[*,0]+counts[*,1]) > max(counts[*,2]+counts[*,3])*2.0]
	plot,[0],[0],xr=xr,/xsty,yr=yr,xtit='Number Photons',ytit='Calc. Number of Photons',/xlog,tit=tit,/ysty,/ylog,chars=2.0
	oploterror,nphot,counts[*,0],counts[*,1],psym=-8,color=200,errcolor=200
	oploterror,nphot,counts[*,2],counts[*,3],psym=-8,color=100,errcolor=100
	makesym,9	
	oploterror,nphot,counts[*,4],counts[*,5],psym=-8,color=120,errcolor=120
	oploterror,nphot,counts[*,6],counts[*,7],psym=-8,color=150,errcolor=150
	makesym,10
	oplot,nphot,nphot,linestyle=2.0
	oplot,nphot,nphot+sqrt(nphot+back*nspec),color=30
	oplot,nphot,nphot-sqrt(nphot+back*nspec),color=30
	xyouts,1000,8,'dlam='+num2str(dlam,dp=1)
	xyouts,1000,20,'v!lI!n='+num2str(v,dp=1)+' [km/s]'
	xyouts,1000,40,'T!lI!n='+num2str(ti/1.0e3,dp=1)+' [keV]'


	;shift
	xr=[min(nphot)/2,2.0*max(nphot)]
	tmp=where(finite(center[*,0]) EQ 1 AND finite(center[*,1]) EQ 1 AND finite(center[*,2]) EQ 1 AND  finite(center[*,3]) EQ 1)
	yr=[(min(center[tmp,0]-center[tmp,1]) < min(center[tmp,2]-center[tmp,3])) > a[1]-1.0e-3,(max(center[tmp,0]+center[tmp,1]) > max(center[tmp,2]+center[tmp,3])) < a[1]+1.0e-3]
	yr=[yr[0]-0.05*(yr[1]-yr[0]),yr[1]+0.05*(yr[1]-yr[0])]
	plot,[0],[0],xr=xr,/xsty,yr=yr,xtit='Number Photons',ytit='Calc. Line Center',/xlog,/ysty,chars=2.0
	oploterror,nphot,center[*,0],center[*,1],psym=-8,color=200,errcolor=200
	oploterror,nphot,center[*,2],center[*,3],psym=-8,color=100,errcolor=100
	makesym,9
	oploterror,nphot,center[*,4],center[*,5],psym=-8,color=120,errcolor=120
	oploterror,nphot,center[*,6],center[*,7],psym=-8,color=150,errcolor=150
	makesym,10
	oplot,xr,a[1]*[1,1],linestyle=2.0
	oplot,nphot,a[1]+a[2]/sqrt(nphot)*sqrt(1+(dlam*1.0e-3)^2/a[2]^2*back/nphot),color=30
	oplot,nphot,a[1]-a[2]/sqrt(nphot)*sqrt(1+(dlam*1.0e-3)^2/a[2]^2*back/nphot),color=30
	

	;width
	xr=[min(nphot)/2,2.0*max(nphot)]
	tmp=where(finite(width[*,0]) EQ 1 AND finite(width[*,1]) EQ 1 AND finite(width[*,2]) EQ 1 AND  finite(width[*,3]) EQ 1)
	yr=[(min(width[tmp,0]-width[tmp,1]) < min(width[tmp,2]-width[tmp,3])) > a[2]*0.5,max(width[tmp,0]+width[tmp,1]) > max(width[tmp,2]+width[tmp,3])  < a[2]*2.0]
	yr=[yr[0]-0.05*(yr[1]-yr[0]),yr[1]+0.05*(yr[1]-yr[0])]
	plot,[0],[0],xr=xr,/xsty,yr=yr,xtit='Number Photons',ytit='Calc. Line Width',/xlog,/ysty,chars=2.0
	oploterror,nphot,width[*,0],width[*,1],psym=-8,color=200,errcolor=200
	oploterror,nphot,width[*,2],width[*,3],psym=-8,color=100,errcolor=100
	makesym,9
	oploterror,nphot,width[*,4],width[*,5],psym=-8,color=120,errcolor=120
	oploterror,nphot,width[*,6],width[*,7],psym=-8,color=150,errcolor=150
	makesym,10
	oplot,xr,a[2]*[1,1],linestyle=2.0
	xyouts,11,yr[0]+1.04*(yr[1]-yr[0]),'GAUSSFIT',color=200,chars=ls
	xyouts,30,yr[0]+1.04*(yr[1]-yr[0]),'MPFIT',color=150,chars=ls
	xyouts,60,yr[0]+1.04*(yr[1]-yr[0]),'INT_TABULATED',color=100,chars=ls
	xyouts,250,yr[0]+1.04*(yr[1]-yr[0]),'TRAP_INT',color=120,chars=ls
	xyouts,800,yr[0]+1.04*(yr[1]-yr[0]),'HUTCH',color=30,chars=ls
	oplot,nphot,a[2]+a[2]/sqrt(2*(nphot-1))*sqrt(1+(dlam*1.0e-3)^4/a[2]^4*back/nphot),color=30
	oplot,nphot,a[2]-a[2]/sqrt(2.0*(nphot-1))*sqrt(1+(dlam*1.0e-3)^4/a[2]^4*back/nphot),color=30

	!p.multi=0
	IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
	IF keyword_set(debug) THEN stop

END
