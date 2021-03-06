;+
;NAME:
;	EQ_ELLIPSE
;
;PURPOSE:
;	This function calculates the x values of an ellipse given the y values
;	and the parameters.  A tilted ellipse can be specified.  This is made to work with MPFIT
;
;CALLING SEQUENCE:
;	result=EQ_ELLIPSE(y,param)
;	
;INPUTS:
;	y	FLTARR [n_y] of y values
;	param	FLTARR [4] or [5] of the parameters (xo,yo,a,b,phi)
;
;OUTPUTS:
;	result:	FLTARR [n_y] of the x values of the ellipse
;
;PROCEDURE:
;	This function re-organizes the equation (x-xo)^2/a^2+(y-yo)^2/b^2=1 to
;	solve for x w/r/t y and the parameters a,b,xo,yo.  For a
;	rotated ellipse, the quadradic curve is used
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 1/10
;	6/29/12		M.L. Reinke - fixed bug in fitting function when using a titled ellipse
;
;-

FUNCTION eq_ellipse,y,param
	l=size(param)
	tilt=1
	IF l[1] EQ 4 THEN tilt=0
	IF l[1] EQ 5 THEN BEGIN
		IF  param[4] EQ 0.0 THEN tilt=0
	ENDIF
	IF NOT tilt THEN BEGIN
		xo=param[0]
		yo=param[1]
		a=param[2]
		b=param[3]
		x=-1.0*a*sqrt(1.0-(y-yo)^2/b^2)+xo
	ENDIF ELSE BEGIN
		xo=param[0]
		yo=param[1]
		a=param[2]
		b=param[3]
		phi=param[4]

		;alpha=y*sin(phi)-xo
		;beta=y*cos(phi)-yo
	
		;quad_a=b^2*cos(phi)^2+a^2*sin(phi)^2
		;quad_b=2.0*alpha*b^2*cos(phi)-2.0*beta*a^2*sin(phi)
		;quad_c=alpha^2*b^2+beta^2*a^2-a^2*b^2
                
                ;x=(-1.0*quad_b-sqrt(quad_b^2-4.0*quad_a*quad_c))/(2.0*quad_a)

		bigA=(cos(phi)/a)^2+(sin(phi)/b^2)
		bigB=-2.0*cos(phi)*sin(phi)*(1/a^2-1/b^2)
		bigC=(sin(phi)/a)^2+(cos(phi)/b)^2
		quad_a=bigA
		quad_b=bigB*y-(2*bigA*xo+yo*bigB)
		quad_c=-1.0*(2*bigC*yo+bigB*xo)*y+(bigA*xo^2+bigB*xo*yo+bigC*yo^2-1)+bigC*y^2

 		x=(-1.0*quad_b-sqrt(quad_b^2-4.0*quad_a*quad_c))/(2.0*quad_a)
	ENDELSE
	RETURN,x
END

;+
;NAME:
;	NONLIN_FIT_ELLIPSE
;
;PURPOSE:
;	This function performs a non-linear least squares fit of an ellipse to a set of (x,y) points.
;	It calles EQ_ELLIPSE which fits to the form (x*sin(phi)-xo)^2/a^2+(y*cos(phi)-yo)^2/b^2=1
;
;CALLING SEQUENCE:
;	result=NONLIN_FIT_ELLIPSE(xpt,ypt,xerr)
;
;INPUTS:
;	xpt	FLTARR	[n] of the x values
;	ypt	FLTARR	[n] of the y values
;	xerr	FLTARR	[n] of the uncertainty in the x point (usually from the spectral fit)
;	
;OPTIONAL INPUTS:
;	good	INTARR 	[n] of 1's or 0's where the	
;
;KEYWORD PARAMETERS
;	/tilt will allow the ellipse to be tilted so that the major/minor axes are not aligned with the x/y axes.
;
;PROCEDURE:
;	The values of xpt and xerr are checked to make sure they are both finite.  If not then those points
;	are set to GOOD=0.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 3/10
;	
;-

FUNCTION nonlin_fit_ellipse,xpt,ypt,xerr,tilt=tilt,good=good,plot=plot,bins=bins,bad=bad,yr=yr

	npts=n(xpt)+1.0
	IF NOT keyword_set(good) THEN good=intarr(npts)+1
	tmp=where(finite(xpt) EQ 0)
	IF tmp[0] NE -1 THEN good[tmp]=0
	tmp=where(finite(xerr) EQ 0)
	IF tmp[0] NE -1 THEN good[tmp]=0


	estimate=double([1.0e3,250.0,1.0e3,1.0e3])
	IF keyword_set(tilt) THEN estimate=[estimate,0.01]
	tmp=where(good EQ 1)
	p=mpfitfun('eq_ellipse',ypt[tmp],xpt[tmp],xerr[tmp],estimate)
	
	IF keyword_set(plot) THEN BEGIN
		openwin,0
		plot,ypt,xpt,psym=3
		xplot=eq_ellipse(ypt,p)
		oplot,ypt,xplot,color=200
		resid=xplot-xpt
		openwin,1
		IF keyword_set(yr) THEN yrange=[-1.0*yr,yr] ELSE yr=[-0.5,0.5]
		plot,ypt,resid,psym=8,yr=yrange,symsize=0.5
		
		IF keyword_set(xerr) THEN oploterror,ypt,resid,xerr,psym=8,symsize=0.5
		IF NOT keyword_set(bins) THEN bins=20
		num_in_bin=floor(npts/bins)
		FOR i=0,bins-1 DO BEGIN
			sub_array=resid[i*num_in_bin:(i+1)*num_in_bin-1]
			mean=mean(sub_array)
			stdev=stdev(sub_array)
			oploterror, [num_in_bin*(i+0.5)],[mean],[stdev],color=200,errcolor=200,psym=8
		ENDFOR
		bad=reverse(sort(abs(resid)))
	ENDIF			
	RETURN,p
END

;+
;NAME:
;	LINLSQ_FIT_ELLIPSE
;
;PURPOSE:
;	This function uses a linear least squares fit to find the quadradic curve coefficients.
;-

FUNCTION linlsq_fit_ellipse,xpt,ypt,debug=debug,plot=plot
	
	tmp=where(xpt NE 0)
	xfit=xpt[tmp]
	yfit=ypt[tmp]
	npts=n(tmp)+1
	D=fltarr(npts,6)
	FOR i=0,npts-1 DO D[i,*]=[xfit[i]^2,xfit[i]*yfit[i],yfit[i]^2,xfit[i],yfit[i],1.0]

	C=fltarr(6,6)
	C[*,0]=[0,0,2.,0,0,0]
	C[*,1]=[0,-1.,0,0,0,0]
	C[*,2]=[2.,0,0,0,0,0]

	S=transpose(D)#D

	evals=la_eigenproblem(S,C,eigenvectors=evecs,/double)
	index=where(real_part(evals) GT 0 AND finite(real_part(evals)) EQ 1)
	evecs=real_part(evecs)
	u=evecs[*,index]
	mu=sqrt(1.0/(transpose(u)#C#u))
	a=mu[0]*u


	IF keyword_set(plot) THEN BEGIN
		openwin,0
		plot,ypt,xpt,psym=3
		xplot=ellipse_xpt(a,ypt)
		oplot,ypt,xplot,color=200
		openwin,1
		plot,ypt,xplot-xpt,psym=8
	ENDIF

	IF keyword_set(debug) THEN stop
	RETURN,a
	
END

;+
;NAME:
;	ELLIPSE_XPT
;
;PURPOSE:
;	This function calculates the x values of an ellipse given the y values the parameters of a quadratic curve [a,b,c,d,f,g]
;	according to the function a*x^2+2*b*x*y+c*y^2+2*d*x+2*f*y+g=0
;
;CALLING SEQUENCE:
;	result=ellipse_xpt(vec,ypt)
;
;INPUTS:
;	vec	FLTARR [5] of the quadradic curve coefficients [a,b,c,d,f,g]
;	ypt	FLTARR [n] of the y values
;
;OPTIONAL INPUTS:
;	noise	FLOAT that is used to add noise to the xpoint (xpt=xpt+noise) using randomn
;
;KEYWORD PARAMETERS:
;	/neg uses the negative root to solve the equation
;
;OUTPUTS:
;	result	FLTARR [n] of the x values on the ellipse
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 3/10
;
;-

FUNCTION ellipse_xpt,vec,ypt,neg=neg,noise=noise
	IF keyword_set(neg) THEN ff=-1.0 ELSE ff=1.0
	IF NOT keyword_set(noise) THEN noise=0.0
		
	npts=n(ypt)+1
	xpt=fltarr(npts)
	
	FOR i=0,npts-1 DO BEGIN
		alpha=vec[0]
		beta=vec[1]*ypt[i]+vec[3]
		gamma=vec[2]*ypt[i]^2+vec[4]*ypt[i]+vec[5]
		xpt[i]=(-1.0*beta+ff*sqrt(beta^2-4.0*alpha*gamma))/(2.0*alpha)+noise*randomn(seed)
	ENDFOR
	RETURN,xpt
END

