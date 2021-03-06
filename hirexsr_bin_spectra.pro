;+
;NAME:
;	OPTIMIZE_BINNING
;
;PURPOSE:
;	This procedure is used to create single or two gaussian spectra that are distributed on an ellipse
;	for demonstration of binning techniques
;
;CALLING SEQUENCE:
;	OPTIMIZE_BINNING
;
;OPTIONAL INPUTS:
;	nrow	INT	number of rows to define bin [DEFAULT=1]
;	nphot	INT	number of photons in the main spectral line [DEFAULT=1000]
;	nback	INT	number of DC offset counts in each pixel [DEFAULT=50]
;	nbin	INT	number of measurements to bin in wavelength [DEFAULT=1 (no binning)]
;	ngauss	INT	number of gaussians to include (1 or 2 for now) [DEFAULT=1]
;	r0	INT	lower bound of rows (r0 < rows < r0+nrow) [DEFAULT=242] (stay within 1 < 487 to be a PILATUS)
;	ti	FLT	ion temperature [eV] DEFAULT=1.5e3
;	v	FLT	flow velocity [m/s] DEFAULT=10.0e3
;	
;KEYWORD PARAMETERS:
;	plot	/plot	will display a plot of the data and the fits
;
;OUTPUTS:
;	dv	FLT	difference between the fit and defined velocity [m/s]
;	verr	FLT	uncertainty in the velocity [m/s]
;
;PROCEDURE:
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke (2010)
;	4/20/12		M.L. Reinke - added second gaussian capability 
;
;-

PRO optimize_binning,nrow=nrow,ti=ti,v=v,nphot=nphot,nback=nback,nbin=nbin,r0=r0,plot=plot,dv=dv,verr=verr,ngauss=ngauss
	IF NOT keyword_set(nrow) THEN nrow=1
	IF NOT keyword_set(nphot) THEN nphot=1.0e3	;number of photons
	IF NOT keyword_set(ti) THEN ti=1.5e3		;ion temperature [eV]
	IF NOT keyword_set(v) THEN v=10.0e3		;observed doppler velocity [m/s]
	IF NOT keyword_set(nbin) THEN nbin=1		;number of pixels to bin together
	IF NOT keyword_set(nback) THEN nback=50		;number of background photons
	IF NOT keyword_set(ngauss) THEN ngauss=1	;number of gaussians to include
	IF NOT keyword_set(r0) THEN r0=242		;row center 

	npix=46
	lam=fltarr(npix,nrow)
	elam=fltarr(npix,nrow)
	sig=fltarr(npix,nrow)
	esig=fltarr(npix,nrow)

	lam_w=[3.99417,3.98999]				;z-line in Angstroms, k-line
	e=1.61e-19					;conversion for Joules/eV
	m=40.0*1.66e-27					;argon mass in kg
	c=2.998e8					;speed on light m/s

	c2=sqrt(lam_w^2*e*ti/(m*c^2))*[1.0,1.0]		;linewidth in Angstroms
	c0=nphot/(2.0*sqrt(!pi/2.0))/c2*[1.0,0.2]	;constant normalized so integral over line has nphot 
	c1=(lam_w-lam_w/c*v)*[1.0,1.0]			;lineshift in Angstroms

	;load ellipse coefs from core He-like module on 1070830020
	;load_ellipse_coefs,1070830020,2,coefs,lam_o
	hirexsr_load_ellipse,1101209008,2,coefs,ifit_e,outl,lam_o,mu,sigma,bad

	nlines=n(lam_o)+1
	pix_o=fltarr(nlines)
	pix=make(140,185,npix)
	FOR i=0,nrow-1 DO BEGIN
		FOR j=0,nlines-1 DO pix_o[j]=eq_ellipse(r0+i,coefs[*,j])
		poly_coefs=poly_fit(pix_o,lam_o,1)
		edgeu=poly(pix+0.5,poly_coefs)
		edgel=poly(pix-0.5,poly_coefs)
		FOR j=0,npix-1 DO BEGIN
			lam[j,i]=0.5*(edgeu[j]+edgel[j])
			elam[j,i]=0.5*(edgeu[j]-edgel[j])
			FOR k=0,ngauss-1 DO sig[j,i]+=fix(sqrt(!pi/2.0)*c0[k]*c2[k]*(erf((c1[k]-edgel[j])/(sqrt(2.0)*c2[k]))-erf((c1[k]-edgeu[j])/(sqrt(2.0)*c2[k]))))
			sig[j,i]+=nback
			sig[j,i]+=fix(randomn(seed,/normal)*sqrt(sig[j,i]))
			sig[j,i]=sig[j,i] > 0		;make sure that signal is > 0
			esig[j,i]=sqrt(sig[j,i]) > 1
		ENDFOR
	ENDFOR
	order=sort(lam)
	lam=lam[order]
	elam=elam[order]
	sig=sig[order]
	esig=esig[order]

	IF nbin NE 1 AND nrow NE 1 THEN BEGIN
		nnew=fix(npix*nrow/nbin)
		new_lam=fltarr(nnew)
		new_elam=fltarr(nnew)
		new_sig=fltarr(nnew)
		new_esig=fltarr(nnew)
			
		FOR i=0,nnew-1 DO BEGIN
			new_elam[i]=mean(elam[i*nbin:(i+1)*nbin-1])
			new_lam[i]=total(lam[i*nbin:(i+1)*nbin-1]*sig[i*nbin:(i+1)*nbin-1])/total(sig[i*nbin:(i+1)*nbin-1])
			new_sig[i]=total(sig[i*nbin:(i+1)*nbin-1])/nbin
			new_esig[i]=sqrt(total(esig[i*nbin:(i+1)*nbin-1]^2))/nbin
		ENDFOR
		lam=new_lam
		elam=new_elam
		sig=new_sig
		esig=new_esig
	ENDIF

	IF ngauss EQ 1 THEN BEGIN
		estimate=[c0[0]*0.9,c1[0]-0.001,c2[0]*0.8,min(sig)]
	ENDIF ELSE BEGIN
		estimate=[c0[0]*0.9,c1[0]-0.001,c2[0]*0.8,c0[1]*0.9,c1[1]-0.001,c2[1]*0.8,min(sig)]	
	ENDELSE
	p = mpfitfun('gaussian_fits', lam,sig,esig, estimate,/quiet,niter=niter,perror=perror,xtol=1.0e-20)

	vfit=(lam_w[0]-p[1])/lam_w[0]*c
	dv=(vfit-v)
	verr1=(lam_w[0]-(p[1]-perror[1]))/lam_w[0]*c
	verr2=(lam_w[0]-(p[1]+perror[1]))/lam_w[0]*c
	verr=0.5*(abs(vfit-verr1)+abs(vfit-verr2))

	IF keyword_set(plot) THEN BEGIN
		lamfit=make(min(lam),max(lam),2000)
		sigfit = gaussian_fits(lamfit, p)
		tit='|VFIT-VREAL|='+num2str(dv/1.0e3,dp=2)+' +/- '+num2str(verr/1.0e3,dp=2)+' [km/s]'
		xr=[min(lam),max(lam)]
		yr=[0,max(sig)*1.05]
		plot,[0],[0],yr=yr,xr=xr,/xsty,/ysty,xtit=n2g('lambda')+' [Ang]',ytit='# Photons',chars=1.2,tit=tit
		oploterror,lam,sig,elam,esig,psym=8
		oplot,lamfit,sigfit,color=200
		oplot,lam_w[0]*[1.0,1.0],[0,2*max(sig)],linestyle=2,color=100
		xyouts,0.1*(xr[1]-xr[0])+xr[0],0.90*yr[1],'#PHOT='+num2str(int(nphot),1)
		xyouts,0.1*(xr[1]-xr[0])+xr[0],0.82*yr[1],'#ROW='+num2str(nrow,1)
		xyouts,0.1*(xr[1]-xr[0])+xr[0],0.74*yr[1],'#BIN='+num2str(nbin,1)
		xyouts,0.1*(xr[1]-xr[0])+xr[0],0.66*yr[1],'R!i0!n='+num2str(r0,1)
	ENDIF
END

PRO plot_binning,load=load
	
	path='/home/mlreinke/idl/hirexsr/binning_test.dat'

	IF NOT keyword_set(load) THEN BEGIN
		nphot=[1.0e2,5.0e2,1.0e3,5.0e3,1.0e4,5.0e4]
		r0=[5,115,242]
		row=[1,5,10,20]
		nfits=100
		num_phot=n(nphot)+1
		num_r0=n(r0)+1
		num_rows=n(row)+1
		
		
		dv=fltarr(num_phot,num_rows,num_r0,nfits)
		verr=fltarr(num_phot,num_rows,num_r0,nfits)
		bin_dv=fltarr(num_phot,num_rows,num_r0,nfits)
		bin_verr=fltarr(num_phot,num_rows,num_r0,nfits)
		
		FOR i=0,num_phot-1 DO BEGIN
			FOR j=0,num_rows-1 DO BEGIN
				FOR k=0,num_r0-1 DO BEGIN
					FOR m=0,nfits-1 DO BEGIN
						OPTIMIZE_BINNING,nrow=row[j],r0=r0[k],nphot=nphot[i],dv=idv,verr=iverr,ngauss=2
						dv[i,j,k,m]=idv
						verr[i,j,k,m]=iverr
						OPTIMIZE_BINNING,nrow=row[j],r0=r0[k],nphot=nphot[i],dv=idv,verr=iverr,nbin=nrow,ngauss=2
						bin_dv[i,j,k,m]=idv
						bin_verr[i,j,k,m]=iverr

					ENDFOR
				ENDFOR
			ENDFOR
		ENDFOR

		save,dv,verr,bin_dv,bin_verr,nphot,r0,row,nfits,filename=path
	ENDIF ELSE restore, path
	dv/=1.0e3
	verr/=1.0e3
	bin_dv/=1.0e3
	bin_verr/=1.0e3	
	color=[0,30,100,120,150,200]
	sym=[4,5,6,7]
	xlab1=10^(make(4.3,5.6,3))
	xlab1=[xlab1,xlab1]
	ylab1=[-10,-10,-10,-15,-15,-15]+2
	xlab2=10^(make(5.0,5.6,2))
	xlab2=[xlab2,xlab2]
	ylab2=[5,5,10,10]
	
	;stop

	FOR k=0,n(r0) DO BEGIN
		openwin,k
		!p.multi=[0,0,2]
		plot,[1],[1],xr=[0.5e2,2e6],/xlog,yr=[-15,15],/ysty,xtit='# Photons in Bin',ytit=n2g('Delta')+'V [km/s]',tit='ROW: '+num2str(r0[k],1)+' NBIN=1',/xsty
		nfits=nfits
		num_rows=n(row)+1
		num_phot=n(nphot)+1
		makesym,10
		xyouts,2e3,mean(ylab1),'N!lPHOT!n in LINE:'
		xyouts,8e3,mean(ylab2),'N!lROWS!n in BIN:'
		FOR i=0,num_phot-1 DO BEGIN
			xyouts,xlab1[i],ylab1[i],num2str(nphot[i],dp=1,/sn),color=color[i]
			FOR j=0,num_rows-1 DO BEGIN
				IF i EQ 0 THEN BEGIN
					oplot,[xlab2[j]],[ylab2[j]],psym=sym[j]
					xyouts,xlab2[j],ylab2[j]-1,'   : '+num2str(row[j],1)
				ENDIF
				oploterror,[nphot[i]*row[j]],mean(dv[i,j,k,*]),stdev(dv[i,j,k,*]),psym=sym[j],color=color[i],errcolor=color[i]
			ENDFOR
		ENDFOR
		oplot,[1,1.0e8],[0,0],linestyle=1
		xplot=10^(make(2,6,100))
		oplot,xplot,1/sqrt(xplot)*2.0e2,linestyle=2
		oplot,xplot,-1/sqrt(xplot)*2.0e2,linestyle=2
	
		makesym,9
		plot,[1],[1],xr=[0.5e2,2e6],/xlog,yr=[-15,15],/ysty,xtit='# Photons in Bin',ytit=n2g('Delta')+'V [km/s]',tit='ROW: '+num2str(r0[k],1)+' NBIN=NROW',/xsty
		FOR i=0,num_phot-1 DO BEGIN
			FOR j=0,num_rows-1 DO BEGIN
				oploterror,[nphot[i]*row[j]],mean(bin_dv[i,j,k,*]),stdev(bin_dv[i,j,k,*]),psym=sym[j],color=color[i],errcolor=color[i]
			ENDFOR
		ENDFOR
		oplot,[1,1.0e8],[0,0],linestyle=1
		oplot,xplot,1/sqrt(xplot)*2.0e2,linestyle=2
		oplot,xplot,-1/sqrt(xplot)*2.0e2,linestyle=2
	ENDFOR
	!p.multi=0


END

;+
;NAME:
;	HIREXSR_AUTOCHMAP
;
;PURPOSE:
;	This function proceduces a channel map by assuming a fixed number of channels per submodule that evenly
;	use the area offset by a given amount from the submodule boundaries.
;
;CALLING SEQUENCE:
;	result=HIREXSR_AUTOCHMAP(nsub,noff)
;	
;INPUTS:
;	nsub	INT 	number of channels per sub module
;	noff	INT	number of pixels to offset the nsub channels from the submodule boundaries
;	
;OPTIONAL INPUTS:
;	nx	INT	image size in the spatial direction DEFAULT: 3*487
;	ny	INT	image size in the spectral direction DEFAULT: 195
;	submod	INT	number of total submodules in the spatial direction DEFAULT: 24
;
;OUTPUTS:
;	result:	INTARR	[nx,ny] of the chmap to be used in HIREXSR_BIN_SPEC.  Values will be -1 for pixels
;			that are not used and 0->chmax for pixels that correspond to various channels
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 4/10
;
;-

FUNCTION hirexsr_autochmap,nsub,noff,nx=nx,ny=ny,submod=submod
	IF NOT keyword_set(nx) THEN nx=3*487
	IF NOT keyword_set(ny) THEN ny=195
	IF NOT keyword_set(submode) THEN submod=24

	rpersub=61			;487/8	for some reason is not an INT
	chmap=intarr(nx,ny)-1		;set values to -1 by default and set desired pixels to CH#
	nch=nsub*submod
	cntr=0
	delta=int((rpersub-2.0*noff)/nsub)	;width of channel
	FOR i=0,submod-1 DO BEGIN
		FOR j=0,nsub-1 DO BEGIN
			start=int(i*rpersub+noff)+j*delta
			endpt=start+delta-1
			IF endpt GE nx THEN delta-=endpt-nx+1
			FOR k=start,start+delta-1 DO chmap[k,*]=fltarr(ny)+cntr
			cntr+=1
		ENDFOR
	ENDFOR

	RETURN,chmap
END

;+
;NAME:
;	HIREXSR_AUTOTMAP
;
;PURPOSE:
;	This function creates a time binning map by binning together a fixed number of frames 
;	for the whole time base.
;
;CALLING SEQUENCE:
;	result=HIREXSR_AUTOTMAP(nt,nframe)
;
;INPUTS:
;	nt	INT	number of raw images collected
;	nframe	INT	number of frames to bin together to form a data frame
;
;OPTIONAL INPUTS:
;	tr	FLTARR	[2] of the lower and upper time points to truncate the map to [sec]
;	shot	LONG	of the shot number to check the time base for
;
;OUTPUTS:
;	result	INTARR	[nt] of the data frame each raw image corresponds to.  This is used
;			as an input in HIREXSR_BIN_IMAGE and HIREXSR_BIN_TIME
;
;MODFICIATION HISTORY:
;	Written by:	M.L. Reinke 4/10
;	4/27/11		M.L. Reinke - added the ability to truncate the autotmap between tr[0] and tr[1] for a given shot
;
;-

FUNCTION hirexsr_autotmap,nt,nframe,tr=tr,shot=shot
	tmap=intarr(nt)
	nbins=int(nt/nframe)
	FOR i=0,nbins-1 DO tmap[i*nframe:(i+1)*nframe-1]=i	
	ex=nt-int(nt/nframe)*nframe	;number left over
	IF ex NE 0 THEN tmap[nt-1-(ex-1):nt-1]=nbins	;set the rest into the last time bin
	IF keyword_set(tr) AND keyword_set(shot) THEN BEGIN 	
		hirexsr_load_image,shot,2,image,time,/noimage
		ilow=ipt(time,tr[0])
		tmp=where(tmap EQ tmap[ilow]-1)
		tmap[0:last(tmp)]=-1
		tmap[last(tmp)+1:*]-=tmap[ilow]
		iup=ipt(time,tr[1])
		tmp=where(tmap EQ tmap[iup]+1)
		tmap[tmp[0]:*]=-1
	ENDIF
	RETURN,tmap
END

;+
;NAME:
;	HIREXSR_BIN_TIME
;
;PURPOSE:
;	This function calculates a new time base given a mapping of collected frames to data frames
;
;CALLING SEQUENCE:
;	result=HIREXSR_BIN_TIME(t,tmap)
;
;INPUTS:
;	t	FLTARR	[n_time] of the time points at which the raw images were collected
;	tmap	INTARR	[n_time] of what data frame to assign each raw image to (-1 if frame is not to be used)
;
;OUTPUTS:
;	result:	FLTARR 	[n_time] of the time values of the data frames.
;
;PROCEDURE:
;	If there is any binning specified by TMAP [tmap NE indgen(n_time)] then values of the output (tau) will be equal to -1.
;	This allows the values of tau to be modified to reflect dynamic binning adjustments w/o modifying the size of the array.
;	You can find the number of data bins by using the length of where(tau NE -1).
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke 4/7/10
;		
;-


FUNCTION hirexsr_bin_time,t,tmap
	ntbins=n(tmap)+1	;less than or equal to number of frames collected which is variable depending on integration time
	tau=fltarr(ntbins)-1
	FOR j=0,ntbins-1 DO BEGIN
		tmp=where(tmap EQ j)				;find all frames to be included in time bin j
		IF tmp[0] NE -1 THEN tau[j]=mean(t[tmp])		;average to find new time
	ENDFOR

	RETURN,tau
END

;+
;NAME:
;	HIREXSR_BIN_IMAGE
;
;PURPOSE:
;	This function applies time binning to the ptrarr array of raw data, creating a new ptrarr of binned images (data frames)
;
;CALLING SEQUENCE:
;	result=HIREXSR_BIN_IMAGE(raw,tmap)
;
;INPUTS:
;	raw	PTRARR	[n_time] of pointers to the raw images
;	tmap	FLTARR	[n_time] of the data frames
;
;OUTPUTS:
;	result	PTRARR	[n_time] of pointers to the binned data frames.  
;
;PROCEDURE:
;	Like with HIREXSR_BIN_TIME, the length of the output is a fixed length but the data values will change.
;	If there is a valid data frame then *result[i]=FLTARR instead of -1.  Using tmap will be a more efficient
;	way to see how many elements of output will be images.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 4/8/10
;
;-

FUNCTION hirexsr_bin_image,raw,tmap
	x=size(*raw[0])
	ntbins=n(tmap)+1			
	cnts=ptrarr(ntbins,/allocate)
	FOR j=0,ntbins-1 DO BEGIN
		tmp=where(tmap EQ j)		;find all frames to be included in time bin j
		iraw=fltarr(x[1],x[2])
		IF tmp[0] NE -1 THEN BEGIN		
			FOR k=0,n(tmp) DO iraw+=*raw[tmp[k]]	;sum frames over time bin assuming raw is a pointer array of the frames
			*cnts[j]=iraw
		ENDIF ELSE *cnts[j]=-1
	ENDFOR

	RETURN,cnts
END
	
;+
;NAME:
;	HIREXSR_BIN_SPEC
;
;PURPOSE:
;	This function takes a set of data frames and uses a channel mapping to create a set of
;	(cnts,lambda) points that define a spectral channel.
;
;CALLING SEQUENCE
;	result=HIREXSR_BIN_SPEC(cnts,tau,chmap,tch,lambda)
;	
;INPUTS:
;	cnts	PTRARR	[n_time] where each element calls a FLTARR [nx,ny] or a -1 (output of HIREXSR_BIN_IMAGE)
;	tau	FLTARR	[n_time] where each element is a time point or -1 (output of HIREXSR_BIN_TIME)
;	chmap	INTARR	[nx,ny,nmaps] where elements specificy the channel number that the pixel belongs to
;			Can make it time evolving (nmaps > 1) and make the number of channels change with time
;	tch	FLTARR	[nmaps] of the time points of when to change the chmap
;	lambda	FLTARR	[nx,ny] of the wavelength values for each pixel
;	tmap	INTARR	[n_time] indicating which frames are included in each time point
;
;OPTIONAL INPUTS:
;	nchbins	INT	max number of of channels DEFAULT: 96	
;	const	FLTARR	[nx,ny] of the multiplicative constant to turn data from counts to spectral brightness [ph/s/m^2/Ang]
;
;OUTPUTS:
;	result:	PTRARR	[nchbins,n_time] of pointers to the FLTARR [npts,3] or -1 values that define the (specbr,lambda,sigma) data for each
;			channel at each point in time.  Both specbr and sigma are in units of [ph/s/m^2/Ang] if the optional input CONST is used properly.
;	
;PROCEDURE:
;	As with other binning programs, the size of the pointer array is fixed at some maxium number of space and time points.
;	Where there is no data for that channel there will be -1 value.  The time-evolving channel mapping allows for each time point
;	to have a variable number of channels and each channel to be made up of a variable number of pixels.
;
;	The sigma is simply the sqrt(cnts)*const.  Doing it explicitely saves time when browsing through the data in the widget.  Also
;	when HIREXSR_AVE_SPEC is used, it won't be that sigma = sqrt(cnts).
;
;MODFICATION HISTORY:
;	Written by:	M.L. Reinke 4/10
;	4/9/10		M.L. Reinke - make each ptr [npts,3] with [*,2]=uncertainty
;	8/9/10		M.L. Reinke - added the use of the CONST optional input
;	3/30/11		M.L. Reinke - added the use of TMAP input to adjust spectra brightness when binning
;	7/26/16		M.L. Reinke - added a sorting call for the  elements w/r/t lambda, avoiding
;                                     an AVESPEC issue when there's high curvature
;-

FUNCTION hirexsr_bin_spec,cnts,tau,chmap,tch,lambda,tmap,const=const,nchbins=nchbins,coefs=coefs

	x=size(*cnts[0])
	IF NOT keyword_set(const) THEN const=fltarr(x[1],x[2])+1.0
	ntbins=n(tau)+1							;less than or equal to number of frames collected which is variable depending on integration time
	IF NOT keyword_set(nchbins) THEN nchbins=96			;unless specified, have hard limit at 96 channels
	nchmap=n(tch)+1							;number of different channel maps
	ichmap=chmap[*,*,0]						;initialize map index which will save time if it's the only one
	mapindex=0

	spec=ptrarr(nchbins,ntbins, /allocate_heap)
	FOR j=0,ntbins-1 DO BEGIN
		icnts=*cnts[j]
		IF icnts[0] NE -1 THEN BEGIN				;if a data frame exsists, apply chmap to get spectra
			nframes=n(where(tmap EQ j))+1			;number of frames binned to get spectra			
			IF nchmap NE 1 THEN BEGIN			;select proper channel map
				tmp=where(tch LE tau[j])	
				mapindex=last(tmp)
				ichmap=chmap[*,*,mapindex]
			ENDIF	
			FOR i=0,nchbins-1 DO BEGIN
				tmp=where(ichmap EQ i)
				IF tmp[0] NE -1 THEN BEGIN
					ispec=fltarr(n(tmp)+1,4)			;(cnts,lam,sig_cnts,resid) data 
					ispec[*,0]=icnts[tmp]*const[tmp]/nframes	;const is calculated assuming a single frames integration time
					ispec[*,1]=lambda[tmp]
					ispec[*,2]=sqrt(icnts[tmp])*const[tmp]/nframes
					IF keyword_set(coefs) THEN BEGIN		;if fit coefs are given then apply them
						;eventually just add all the fits together 
						p=*coefs[i,j,2]
						ispec[*,3]=ispec[*,0]-gaussian_fits(ispec[*,1],p)
					ENDIF ELSE ispec[*,3]=fltarr(n(tmp)+1)-1	;if no coefs then set residual to -1
					order=sort(ispec[*,1])		;reorder pixels to ensure increasing order to aid AVE_SPEC
					ispec[*,0]=ispec[order,0]
					ispec[*,1]=ispec[order,1]
					ispec[*,2]=ispec[order,2]
					ispec[*,3]=ispec[order,3]
					*spec[i,j]=ispec				;write the (cnts,lam,sigma) data to the [i,j]th pointer
				ENDIF ELSE *spec[i,j]=-1
			ENDFOR
		ENDIF ELSE FOR i=0,nchbins-1 DO *spec[i,j]=-1.0		;set nodata in all channels for time bin j

	ENDFOR
	
	RETURN,spec		;returns the pointer array containing the (cnts,lam) data for each ch at each time bin
END

;+
;NAME:
;	HIREXSR_REBIN_AFTERKILL
;
;PURPOSE:
;	This procedure is used to rebin the spectra using a modified CHMAP such as after a small number of pixels
;	have been labeled bad using HIREXSR_ZOOM_KILL in W_HIREXSR_HE_MOMENTS
;
;CALLING SEQUENCE:
;	HIREXSR_REBIN_AFTERKILL,cnts,tau,chmap,tch,lambda,spec,newmap
;
;INPUTS:
;	cnts		PTRARR	[n_time] where each element calls a FLTARR [nx,ny] or a -1 (output of HIREXSR_BIN_IMAGE)
;	tau		FLTARR	[n_time] where each element is a time point or -1 (output of HIREXSR_BIN_TIME)
;	chmap		INTARR	[nx,ny,nmaps] where elements specificy the channel number that the pixel belongs to
;				Can make it time evolving (nmaps > 1) and make the number of channels change with time
;	tch		FLTARR	[nmaps] of the time points of when to change the chmap
;	lambda		FLTARR	[nx,ny] of the wavelength values for each pixel
;	tmap	INTARR	[n_time] indicating which frames are included in each time point
;	spec		PTRARR	[chmax,ntime] of the binned spectra, output from HIREXSR_BIN_SPEC that is to be updated
;	newmap		INTARR	[nx,ny,nmaps] of the new CHMAP 
;
;OPTIONAL INPUTS:
;	const		FLTARR 	[nx,ny] of the multiplicative constant to turn data from counts to spectral brightness [ph/s/m^2/Ang]
;
;OPTIONAL OUTPUTS:
;	fix		INTARR 	[chmax,ntime] of 1's or 0's indicating the spectra has changed and requires updating
;
;PROCEDURE:
;	This procedure looks for changes in the length of the # of points in a specific channel so removing one pixel and adding another
;	will NOT result in that spectrum being updated.
;
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 9/12/10
;	3/30/11		M.L. Reinke - added the use of TMAP input
;
;-


PRO hirexsr_rebin_afterkill,cnts,tau,chmap,tch,lambda,tmap,spec,newmap,const=const,fix=fix


	x=size(*cnts[0])
	IF NOT keyword_set(const) THEN const=fltarr(x[1],x[2])+1.0
	ntbins=n(tau)+1							;less than or equal to number of frames collected which is variable depending on integration time
	x=size(spec)
	nchbins=x[1]
	nchmap=n(tch)+1							;number of different channel maps
	ichmap=chmap[*,*,0]						;initialize map index which will save time if it's the only one
	inewmap=newmap[*,*,0]
	mapindex=0
	fix=intarr(nchbins,ntbins)			;initialize an array of 0's or 1's to indicate which spectra require recalculation

	chk=chmap-newmap
	IF total(chk) EQ 0 THEN RETURN			;do nothing if reference and newmap are the same


	FOR j=0,ntbins-1 DO BEGIN
		icnts=*cnts[j]
		IF icnts[0] NE -1 THEN BEGIN				;if a data frame exsists, apply chmap to get spectra
			nframes=n(where(tmap EQ j))+1			;number of frames binned to get spectra		
			IF nchmap NE 1 THEN BEGIN			;select proper channel map
				tmp=where(tch LE tau[j])	
				mapindex=last(tmp)
				ichmap=champ[*,*,mapindex]
				inewmap=inewmap[*,*,mapindex]
			ENDIF	
			FOR i=0,nchbins-1 DO BEGIN
				ich=where(ichmap EQ i)
				inew=where(inewmap EQ i)
				IF n(ich) NE n(inew) THEN BEGIN				;only change *spec[i,j] data if the two maps are different size
					tmp=inew
					ispec=fltarr(n(tmp)+1,4)			;(cnts,lam,sig_cnts,resid) data 
					ispec[*,0]=icnts[tmp]*const[tmp]/nframes
					ispec[*,1]=lambda[tmp]
					ispec[*,2]=sqrt(icnts[tmp])*const[tmp]/nframes
					ispec[*,3]=-1					
					*spec[i,j]=ispec				;write the (cnts,lam,sigma) data to the [i,j]th pointer
					fix[i,j]=1					;set fix=1 to indicate activity needed
				ENDIF
			ENDFOR
		ENDIF 
	ENDFOR	
END


;+
;NAME:
;	HIREXSR_BIN_BOUNDS
;
;PURPOSE:
;	This function calculates the row boundaries of the bins describe by an input channel map
;
;CALLING SEQUENCE:
;	result=HIREXSR_BIN_BOUNDS(chmap,chmax)
;	
;INPUTS:
;	chmap	INTARR	[nx,ny,nmaps] where elements specificy the channel number that the pixel belongs to
;			Can make it time evolving (nmaps > 1) and make the number of channels change with time
;	chmax	INT	max number of channels
;
;OUTPUTS:
;	result	INTARR	[2,chmax,nmaps] of the upper [0,*,*] and lower [1,*,*] bounds of a channel
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 4/8/10
;
;-

FUNCTION hirexsr_bin_bounds,chmap,chmax
	x=size(chmap)
	nx=x[1]
	ny=x[2]
	IF x[0] EQ 3 THEN nmaps=x[3] ELSE nmaps=1
	bounds=intarr(2,chmax,nmaps)-1			;for specified channels < chmax set bounds to -1
	FOR i=0,nmaps-1 DO BEGIN
		mapslice=chmap[*,0,i]
		n_ch=max(chmap[*,*,i])+1
		FOR j=0,n_ch-1 DO BEGIN
			tmp=where(mapslice EQ j)
			bounds[0,j,i]=tmp[0]
			bounds[1,j,i]=max(tmp)
		ENDFOR
	ENDFOR

	RETURN,bounds
END

;+
;NAME:
;	HIREXSR_AVE_SPEC
;
;PURPOSE:
;	This function is used to average raw spectra over a given number of points to increase SNR and reduce
;	the number of points prior to the nonlinear gaussian fit.  This increases speed without a reduction in accuracy.
;
;CALLING SEQUENCE:
;	result=HIREXSR_AVE_SPEC(spec,nave)
;	
;INPUTS:
;	spec		PTRARR	[chmax,ntime] where each [i,j] is a pointer reference to an [n,4] array (output of HIREXSR_BIN_SPEC)
;	nave		INT	number of points to average over
;
;OPTIONAL INPUTS:
;	avespec		PTRARR	[chmax,ntime] of a previous call to HIREXSR_AVE_SPEC to be used to update isolated spectra
;	fix		INTARR	[chmax,ntime] of 1's or 0's indicating which spectra need to be reaveraged
;	chmap		INTARR	[nx,ny] including chmap enables adaptive NAVE 
;
;OUTPUTS:
;	result		PTRARR	[chmax,ntime] where each [i,j] is a pointer to an [m,4] array w/ m=n/nave.
;
;OPTIONAL OUTPUTS:
;	maxave		INT	of the maximum number of points in any of the [i,j] spectral pointers and
;				is used when writing the AVESPEC to the tree
;
;PROCEDURE:
;	For each spectra the spectral brightness is averaged, the wavelength is a weighted averaged of the lambda values
;	and the noise is added in quadrature and divided by the # of averaging points.  If n/nave is not a whole number then
;	the final point will be averaged over a number < nave.  For the best results on future steps where the moments
;	are taken using a Newton-Cotes integration it's best to average over the # of rows used to form the spectra which 
;	eliminates near degenercy of wavelengths in some regions of the image.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 3/10
;	9/12/10		M.L. Reinke - fixed problem where the uncertainty was being averaged incorrectly
;	9/21/10		M.L. Reinke - modifed the lambda weighted average calculation to avoid NaN when zero photons
;	1/25/13		M.L. Reinke - added the chmap optional input which, if specified sets nave
;                                     for each CH based on #rows in CHMAP
;
;-

FUNCTION hirexsr_ave_spec,spec,nave,maxave=maxave,avespec=avespec,fix=fix,chmap=chmap
	x=size(spec)
	chmax=x[1]
	ntime=x[2]	
	maxave=0
	IF NOT keyword_set(avespec) THEN avespec=ptrarr(chmax,ntime,/allocate_heap)
	IF NOT keyword_set(fix) THEN fix=intarr(chmax,ntime)+1
	FOR i=0,chmax-1 DO BEGIN
		FOR j=0,ntime-1 DO BEGIN
			IF fix[i,j] EQ 1 THEN BEGIN
				ispec=*spec[i,j]
				IF ispec[0] NE -1 THEN BEGIN
					IF keyword_set(chmap) THEN inave=n(where(chmap[*,0,0] EQ i))+1 ELSE inave=nave
					x=size(ispec)
					npts=x[1]
					new=int(npts/inave)
					ex=npts-new*inave
					IF ex EQ 0 THEN BEGIN
						iave=fltarr(new,4) 
						IF maxave LT new THEN maxave=new
					ENDIF ELSE BEGIN
						iave=fltarr(new+1,4)
						IF maxave LT new+1 THEN maxave=new+1
					ENDELSE
				
					FOR k=0,new-1 DO BEGIN
						k1=k*inave
						k2=(k+1)*inave-1
						sum=total(ispec[k1:k2,0])
						iave[k,0]=sum/inave						;average spectra brightness
						;weighted average lambda value (allow for no photons)
						IF sum NE 0 THEN iave[k,1]=total(ispec[k1:k2,1]*ispec[k1:k2,0])/sum ELSE iave[k,1]=mean(ispec[k1:k2,1])
						iave[k,2]=sqrt(total(ispec[k1:k2,2]*ispec[k1:k2,2]))/inave	;new noise
						iave[k,3]=total(ispec[k1:k2,3])/inave

					ENDFOR
					IF ex NE 0 THEN BEGIN
						k1=new*inave
						k2=npts-1
						sum=total(ispec[k1:k2,0])
						iave[new,0]=sum/ex
						IF sum NE 0 THEN iave[new,1]=total(ispec[k1:k2,1]*ispec[k1:k2,0])/sum ELSE iave[k,1]=mean(ispec[k1:k2,1])
						iave[new,2]=sqrt(total(ispec[k1:k2,2]*ispec[k1:k2,2]))/ex
						iave[new,3]=total(ispec[k1:k2,3])/ex		
					ENDIF
					order=sort(iave[*,1])
					FOR k=0,3 DO iave[*,k]=iave[order,k]
					*avespec[i,j]=iave
				ENDIF ELSE *avespec[i,j]=-1.0
			ENDIF
		ENDFOR
	ENDFOR

	RETURN,avespec
END

;+
;NAME:
;	HIREXSR_BIN_POS
;
;PURPOSE:
;	This function forms POS vectors based on the CHMAP and DLAM used to form the moments.
;	It can also bin the etendue as well.
;
;CALLING SEQUENCE:
;	result=HIREXSR_BIN_POS(pos,dlam,chmap,lambda,chmax)
;
;INPUTS:
;	pos		FLTARR	[4,nx,ny] of the POS vector for each point on the detector
;	dlam		FLTARR	[2] of the lower [0] and upper[1] wavelength values used to form the moment
;	chmap		INTARR	[nx,ny,nmaps]	that describes the binning used to form the spectra (see HIREXSR_BIN_SPEC)
;	chmax		INT	of the maximum number of spatial channels allowed by the software
;	lambda		FLTARR	[nx,ny] of the wavelength values at each pixel
;
;OPTIONAL INPUTS:
;	u		FLTARR	[nx,ny,nmaps] of the etendue values at each pixel
;	const		FLTARR	[nx,ny] of the calibration constant to turn #photons into brightness
;
;OUTPUTS:
;	result		FLTARR	[4,chmax,nmaps] of the POS vector for each channel.  If no pixels selected for a channel
;				then [*,i,j] is filled with -1
;
;OPTIONAL OUTPUTS:
;	newu		FLTARR	[chmax,nmaps] of the etendue values for each each channel.  If no pixels selected for that
;				channel then [i,j]=-1
;	newc		FLTARR 	[chmax,nmaps] of the constant brightness constant (used to recover #photons from brightness)
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 8/10
;
;-


FUNCTION hirexsr_bin_pos,pos,dlam,chmap,lambda,chmax,u=u,newu=newu,const=const,newc=newc
	x=size(chmap)
	IF x[0] EQ 3 THEN nmaps=x[3] ELSE nmaps=1
	newpos=fltarr(4,chmax,nmaps)
	IF keyword_set(u) THEN newu=fltarr(chmax,nmaps)-1.0
	IF keyword_set(const) THEN newc=fltarr(chmax,nmaps)-1.0
	p0=reform(pos[0,*,*])
	p1=reform(pos[1,*,*])
	p2=reform(pos[2,*,*])
	p3=reform(pos[3,*,*])
	FOR i=0,nmaps-1 DO BEGIN
		FOR j=0,chmax-1 DO BEGIN
			tmp=where(chmap[*,*,i] EQ j AND lambda GE dlam[0] AND lambda LE dlam[1])
			IF tmp[0] NE -1 THEN BEGIN
				newpos[0,j,i]=mean(p0[tmp])
				newpos[1,j,i]=mean(p1[tmp])
				newpos[2,j,i]=mean(p2[tmp])
				newpos[3,j,i]=mean(p3[tmp])
				IF keyword_set(u) THEN newu[j,i]=mean(u[tmp])
				IF keyword_set(const) THEN newc[j,i]=mean(const[tmp])
			ENDIF ELSE newpos[*,j,i]=-1
		ENDFOR
	ENDFOR
	RETURN,newpos
END

;+
;NAME:
;	HIREXSR_BIN_INST
;
;PURPOSE:
;	This function forms instrumental data based on the CHMAP and DLAM used to form the moments.
;
;CALLING SEQUENCE:
;	result=HIREXSR_BIN_INST(inst,dlam,chmap,lambda,chmax)
;
;INPUTS:
;	inst		FLTARR	[nx,ny] of the instrumental for each point on the detector
;	dlam		FLTARR	[2] of the lower [0] and upper[1] wavelength values used to form the moment
;	chmap		INTARR	[nx,ny,nmaps]	that describes the binning used to form the spectra (see HIREXSR_BIN_SPEC)
;	chmax		INT	of the maximum number of spatial channels allowed by the software
;	lambda		FLTARR	[nx,ny] of the wavelength values at each pixel
;;
;OUTPUTS:
;	result		FLTARR	[chmax,nmaps] of the inst vector for each channel.  If no pixels selected for a channel
;				then [i,j] is filled with -1
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 9/4/11 (adapted from HIREXSR_BIN_POS)
;
;-


FUNCTION hirexsr_bin_inst,inst,dlam,chmap,lambda,chmax
	x=size(chmap)
	IF x[0] EQ 3 THEN nmaps=x[3] ELSE nmaps=1
	newinst=fltarr(chmax,nmaps)
	FOR i=0,nmaps-1 DO BEGIN
		FOR j=0,chmax-1 DO BEGIN
			tmp=where(chmap[*,*,i] EQ j AND lambda GE dlam[0] AND lambda LE dlam[1])
			IF tmp[0] NE -1 THEN BEGIN
				newinst[j,i]=mean(inst[tmp])
			ENDIF ELSE newinst[j,i]=-1
		ENDFOR
	ENDFOR
	RETURN,newinst
END

