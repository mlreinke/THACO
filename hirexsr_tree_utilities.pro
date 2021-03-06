;+
;NAME:
;	HIREXSR_CALC_DATA
;
;-

PRO hirexsr_calc_data,shot,line,verb=verb,debug=debug,nofit=nofit,tree=tree
	nave=20
	double=0	
	CASE line OF 
		0 : gline='w'
		1 : gline='x'
		2 : gline='z'
		3 : gline='lya1'
		4 : gline='4d'
		5 : gline='J'
		6 : gline='w'
	ENDCASE
	z=18
	IF line EQ 4 THEN z=42
	IF line LE 2 THEN h=0 ELSE h=1
	IF line EQ 6 THEN z=20

	hirexsr_load_wavelengths,-1,lam_o,z_o,label_o
	tmp=where(z_o EQ z AND label_o EQ gline)
	lam_o=lam_o[tmp[0]]

	hirexsr_ptr_image,shot,raw,lambda,t,const=const,pos=pos,u=u,h=h
	ntime=n(t)+1
	print, 'raw data loaded to PTRARR'

	;bin images in time
	hirexsr_load_binning,shot,chmap,tch,tmap,good,chmax,h=h
	cnts=hirexsr_bin_image(raw,tmap)
	tau=hirexsr_bin_time(t,tmap)
	tmp=where(tmap EQ -1)
	IF n(tmap) NE n(t) THEN BEGIN
		print, 'TMAP, T are unequal length'
		heap_free,raw
		heap_free,cnts
		heap_gc
		RETURN
	ENDIF
	IF tmp[0] EQ -1 THEN nused=ntime ELSE nused=ntime-(n(tmp)+1)
	print, 'tmap applied: using '+num2str(nused,1)+' of '+num2str(ntime,1)+' frames'

	;form spectra according to chmap
	spec=hirexsr_bin_spec(cnts,tau,chmap,tch,lambda,tmap,const=const,nchbins=chmax,coefs=coefs)
	print, 'chmap applied'

	;average over nave # of points
	avespec= hirexsr_ave_spec(spec,nave,maxave=maxave,chmap=chmap)
	print, 'spectra averaged'

	;fit all spectra for all times
	IF keyword_set(nofit) THEN BEGIN
		hirexsr_load_fits,shot,line,coefs,nave,double,labels,/ptr
		print, 'FITTING SKIPPED - fits loaded'
	ENDIF ELSE BEGIN
		coefs=hirexsr_fit_spectra(avespec,line,double=double,label=labels,verb=verb)
		print, 'spectra fitted'
	ENDELSE

	;calculate moments
	mom=hirexsr_spec2moments(avespec,coefs,labels,lam_o,line,dlam=dlam,bfrac=bfrac,scale=scale,fitcase=fitcase)
	print, 'moments generated'
	
	;calculate pos for each channel
	chpos=hirexsr_bin_pos(pos,lam_o+dlam*[-1.0e-3,1.0e3],chmap,lambda,chmax,u=u,newu=chu,const=const,newc=chc)

	;calculate rhotang for each channel
	rhotang=hirexsr_calc_rhotang(shot,chpos,tch,tau,tree=tree)

	;covert averaged spectra pointer array and write to tree
	avespec_arr=hirexsr_avespec_ptr2arr(avespec,chmax,ntime,maxave)
	specbr=avespec_arr[*,*,*,0]
	lam=avespec_arr[*,*,*,1]
	sig=avespec_arr[*,*,*,2]
	resid=avespec_arr[*,*,*,3]
	hirexsr_write_avespec,shot,specbr,lam,sig,resid,tau,nave,h=h
	print, 'avespec written to tree SHOT: '+num2str(shot,1)

	;convert fit coefs pointer array and write to tree
	coefs_arr=hirexsr_coefs_ptr2arr(coefs,chmax,ntime)
	hirexsr_write_fits,shot,line,coefs_arr,tau,nave,double,labels
	print, 'fit coefficients written to tree SHOT: '+num2str(shot,1)

	;write moments to tree
	hirexsr_write_moments,shot,line,mom[*,*,*,0],mom[*,*,*,1],mom[*,*,*,2],mom[*,*,*,3],tau,chpos,rhotang,bfrac,scale,chu,tch,dlam,double,fitcase,tree
	print, 'moments written to tree SHOT: '+num2str(shot,1)

	IF keyword_set(debug) THEN stop

	;clean up memory
	heap_free,raw
	heap_free,cnts
	heap_free,coefs
	heap_free,spec
	heap_free,avespec
	heap_gc
END

;+
;NAME:
;	HIREXSR_INVERT_DATA
;
;-


PRO hirexsr_invert_data,shot,line,treevox=treevox,nosub=nosub

	IF NOT keyword_set(treevox) THEN BEGIN
		hirexsr_load_moments,shot,line,mom,err,pmom,perr,tau,pos,rhotang,bfrac,scale,u,tpos,dlam,double,fitcase
		ntime=n(tau)+1
		tmp=where(tau GT 0)
		nrho=intarr(ntime)-1
		nrho[tmp]=25
		rhomin=fltarr(ntime)-1.0
		rhomin[tmp]=0
		rhomax=fltarr(ntime)-1.0
		rhomax[tmp]=1.05
		rho=hirexsr_make_rhovec(nrho,rhomin,rhomax)
		hirexsr_pos2voxel,shot,tau,pos,tpos,rho,voxel,velvoxel,/verb
		hirexsr_write_voxel,shot,line,voxel,velvoxel,rho,tau
		heap_free,voxel
		heap_free,velvoxel
		heap_free,rho
	ENDIF
	IF keyword_set(treevox) THEN new=0 ELSE new=1
	hirexsr_load_momentptr,shot,line,mom,tau,pos,tpos,lam_o,z,new=new				;load moment data
	hirexsr_load_voxel,shot,line,voxel,velvoxel,rho,/ptr						;load voxel data
	hirexsr_load_proconfig,shot,line,eps,eta,/quiet,status=status					;load config data
	IF NOT status THEN BEGIN
		ntime=n(tau)+1
		eps=fltarr(3,ntime)+1.0
		eta=fltarr(3,ntime)+0.1
	ENDIF
	hirexsr_calc_profiles,mom,voxel,velvoxel,eps,eta,lam_o,z,rho,profiles,/solidbody,nosub=nosub	;perform inversions
	hirexsr_profile_ptr2arr,profiles,pro_arr,proerr_arr,rho_arr		;convert from PTR to array
	hirexsr_write_profile,shot,line,pro_arr,proerr_arr,rho_arr,tau		;write profiles to the tree

	hirexsr_extract_check,mom,check,good					;extract check and good from momentptr
	hirexsr_write_check,shot,line,check,good				;write check/good to tree
	hirexsr_write_proconfig,shot,line,eps,eta				;write eps/eta to tree

	heap_free,mom
	heap_free,voxel
	heap_free,velvoxel
	heap_free,profiles
	heap_free,rho
	heap_gc


END

PRO test_mom
	shot=1070830020
	line=2
	nave=20
	double=0
	gline='z'
	chmax=96
	ntime=11

	hirexsr_load_wavelengths,-1,lam_o,z,label_o
	tmp=where(z EQ 18 AND label_o EQ gline)
	lam_o=lam_o[tmp[0]]

	hirexsr_load_fits,shot,coefs,ave,double,labels,/zjk,/ptr
	hirexsr_load_avespec,shot,avespec,lam,sig,resid,tau,/ptr

	mom=hirexsr_spec2moments(avespec,coefs,labels,lam_o,line,dlam=dlam)
	stop
END

PRO hirexsr_setup_realtime,shot
	HIREXSR_ADDNEW_ANALYSIS,shot,tht=99
	HIREXSR_COPY_BINNING,-1,shot,totht=99,/chmap,/tmap,/good
	HIREXSR_COPY_BINNING,-1,shot,totht=99,/chmap,/tmap,/good,/h
END

PRO hirexsr_setup_variation,shot,nfits
	FOR i=0,nfits-1 DO BEGIN
		HIREXSR_ADDNEW_ANALYSIS,shot,tht=100+i
		HIREXSR_COPY_BINNING,shot,shot,fromtht=0,totht=100+i,/chmap,/tmap,/good
		HIREXSR_AVESPEC2TREE,shot,tht=100+i,/random
	ENDFOR
END


;+
;NAME:
;	HIREXSR_AUTOBIN
;
;PURPOSE:
;	This procedure streamlines making spatially/temporally uniform
;	CHMAP and TMAPs and writes them to the tree
;
;MODIFICATION HISOTRY:
;	Written by:	M.L. Reinke - 11/2012
;       1/14/13		M.L. Reinke - added the tch=-1 and (hardcoded) chmax optional inputs into the
;				      HIREXSR_WRITE_BINNING CALL
;-

PRO hirexsr_autobin,shot,nt,tr,nb,noff,tht=tht,h=h
	time=hirexsr_get_time(shot)
	tmap=hirexsr_autotmap(n(time)+1,nt,tr=tr,shot=shot)
	ntmap=max(tmap)
	chmap=hirexsr_autochmap(nb,noff)
	IF keyword_set(h) THEN BEGIN
		chmap=chmap[0:486,*]		;truncate to one module
		nchmap=max(chmap[*,0])
		good=intarr(32,n(time)+1)
		good[0:nchmap,0:ntmap]=1
		chmax=32
        ENDIF ELSE BEGIN
		nchmap=max(chmap[*,0])
		good=intarr(96,n(time)+1)
		good[0:nchmap,0:ntmap]=1
		chmax=96
        ENDELSE
	HIREXSR_WRITE_BINNING,shot,chmap=chmap,tch=-1,tmap=tmap,good=good,tht=tht,chmax=chmax,h=h
END

;+
;NAME:
;	HIREXSR_RUN_THACO
;
;PURPOSE:
;	This is a high level procedure as a fire & forget way to run THACO for a given shot
;	it just runs AVESPEC/FITSPEC/INVERT in one call
;-

PRO hirexsr_run_thaco,shot,line,tht=tht,nave=nave,nback=nback,h=h,newtht=newtht,note=note,rmax=rmax,nr=nr,rhovec=rhovec,tree=tree,sine=sine,nofit=nofit,wf=wf
	CASE line OF
		3 : h=1
		4 : h=1
		5 : h=1
		6 : h=1
		ELSE : h=0
	ENDCASE
	IF line EQ 2 THEN rhovec=[0.00,0.01,0.03,0.05,0.07,0.10,0.14,0.18,0.22,0.26,0.30,0.34,0.38,0.42,0.46,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95]
	IF line EQ 3 THEN rhovec=[0.00,0.01,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.45,0.50,0.55]
	IF line EQ 7 THEN rhovec=[0.00,0.01,0.03,0.05,0.07,0.10,0.14,0.18,0.22,0.26,0.30,0.34,0.38,0.42,0.46,0.50,0.55,0.60,0.65,0.70,0.75,0.80]
	hirexsr_avespec2tree,shot,debug=debug,h=h,nave=nave,double=double,tht=tht,note=note,newtht=newtht,wf=wf
	hirexsr_fitspec2tree,shot,line,verb=verb,debug=debug,nofit=nofit,nback=nback,tree=tree,tht=tht
	hirexsr_invert2tree,shot,line,treevox=treevox,nosub=nosub,rmax=rmax,nr=nr,rhovec=rhovec,sine=sine,tht=tht
END

;+
;NAME:
;	HIREXSR_CLONE_ANALYSIS
;	
;PURPOSE:
;	This function will clone an analysis run from an existing THACO# and copy
;	(for a given line) the analysis to a new THACO#
;
;CALLING SEQUENCE:
;	result=HIREXSR_CLONE_ANALYSIS(shot,line,thta,thtb)
;
;INPUTS:
;	shot	LONG	shot number
;	line	INT	line # for the analysis	
;	thta	INT	THACO# to copy analysis from
;	thtb	INT	THACO# to copy analysis to - set to an unassigned variable chose the next available#
;
;KEYWORD PARAMETERS:
;	nobin	/nobin does not copy any of the binning and AVESPEC
;	nomom	/nomom does not copy any of the fitting/moment results
;	nopro	/nopro does not copy any of the inversion results
;	quiet	/quiet supresses text updates to the terminal
;	force	/force will force delete an existing THACO# node and rewrite with new data
;
;OUTPUTS:
;	result	INT 	indicating the result of the analysis
;			-2 : invalid THACO# to copy from (thta)
;			-1 : THACO# to copy to (thtb) already exists
;			 1 : copy of thta -> thtb completed
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 8/19/11
;	11/20/12	M.L. Reinke - added /nobin and moved fitting  to be specified by /nomom
;
;-
	
FUNCTION hirexsr_clone_analysis,shot,line,thta,thtb,nobin=nobin,nomom=nomom,nopro=nopro,quiet=quiet,force=force,add=add
	status=hirexsr_is_analysis(shot,thta)
	IF NOT status THEN BEGIN
		print, 	'THT '+num2str(thta,1)+' does not exist'
		RETURN,-2
	ENDIF

	hirexsr_addnew_analysis,shot,tht=thtb,chk=chk,force=force
	IF chk EQ 0 AND NOT keyword_set(add) THEN BEGIN
		print, 'THT '+num2str(thtb,1)+' already exists, use /force to delete or /add to apply new data'
		RETURN,-1
	ENDIF
	CASE line OF
		3 : h=1
		4 : h=1
		5 : h=1
		6 : h=1
		ELSE : h=0
	ENDCASE	

	;copy binning/avespec
	IF NOT keyword_set(nobin) THEN BEGIN
		hirexsr_load_binning,shot,chmap,tch,tmap,good,chmax,h=h,tht=thta
		IF tch[0] EQ 0 THEN tch[0]= -1
		hirexsr_write_binning,shot,chmap=chmap,tch=tch,tmap=tmap,good=good,chmax=chmax,h=h,tht=thtb
		IF NOT keyword_set(quiet) THEN print, 'binning copied from THT '+num2str(thta,1)+' to '+num2str(thtb,1)
		hirexsr_load_avespec,shot,specbr,lam,sig,resid,tau,nave,h=h,tht=thta
		hirexsr_write_avespec,shot,specbr,lam,sig,resid,tau,nave,h=h,tht=thtb
		IF NOT keyword_set(quiet) THEN print, 'avespec copied from THT '+num2str(thta,1)+' to '+num2str(thtb,1)
	ENDIF

	;copy fitting & moments
	IF NOT keyword_set(nofit) THEN BEGIN
		hirexsr_load_fits,shot,line,coefs,ave,double,labels,ch=ch,tau=tau,tht=thta
		hirexsr_write_fits,shot,line,coefs,tau,ave,double,labels,tht=thtb
		IF NOT keyword_set(quiet) THEN print, 'line='+num2str(line,1)+ ' fitting copied from THT '+num2str(thta,1)+' to '+num2str(thtb,1)
		hirexsr_load_moments,shot,line,mom,err,pmom,perr,tau,pos,rhotang,bfrac,scale,u,tpos,dlam,double,fitcase,tree,tht=thta
		hirexsr_write_moments,shot,line,mom,err,pmom,perr,tau,pos,rhotang,bfrac,scale,u,tpos,dlam,double,fitcase,tree,tht=thtb
		IF NOT keyword_set(quiet) THEN print, 'line='+num2str(line,1)+ ' moments copied from THT '+num2str(thta,1)+' to '+num2str(thtb,1)
	ENDIF

	IF NOT keyword_set(nopro) THEN BEGIN
		;copy voxels
		hirexsr_load_voxel,shot,line,voxel,velvoxel,rho,tau,tree,m1svoxel=m1svoxel,m1svelvoxel=m1svelvoxel,tht=thta
		hirexsr_write_voxel,shot,line,voxel,velvoxel,rho,tau,tree,m1svoxel=m1svoxel,m1svelvoxel=m1svelvoxel,tht=thtb

		;copy profiles
		hirexsr_load_profile,shot,line,pro_arr,proerr_arr,rho_arr,tau,tinst=tinst,tgood=tgood,tht=thta
		hirexsr_write_profile,shot,line,pro_arr,proerr_arr,rho_arr,tau,tinst=tinst,tgood=tgood,tht=thtb
		hirexsr_load_proconfig,shot,line,eps,eta,tht=thta
		hirexsr_write_proconfig,shot,line,eps,eta,tht=thtb
		hirexsr_load_check,shot,line,check,good,subcheck,tht=thta
		hirexsr_write_check,shot,line,check,good,subcheck,tht=thtb
		IF NOT keyword_set(quiet) THEN print, 'line='+num2str(line,1)+ ' inversions copied from THT '+num2str(thta,1)+' to '+num2str(thtb,1)
	ENDIF

	;write comment
	note='analysis cloned from THT='+num2str(thta,1)
	comment=logname()+','+systime()+','+note
	hirexsr_write_analysis_comment,shot,comment,tht=thtb
	RETURN,1
END

;+
;NAME:
;	HIREXSR_CHECK_CALIB
;
;PURPOSE:
;	This procedure looks into the tree to see what shot/times the
;	wavelength/spatial calibration was based on for the given shot
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 11/19/12
;
;-
PRO hirexsr_check_calib,shot,cshot,ctimes,quiet=quiet,status=status
	mdsopen,'spectroscopy',shot
	cshot=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB:SHOT',/quiet,status=status)
	ctimes=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.CALIB:TIMES',/quiet,status=status)
	mdsclose,'spectroscopy',shot
	IF NOT status THEN BEGIN
		cshot=-1
		ctimes=-1	
		IF NOT keyword_set(quiet) THEN print, 'No Calibration Data Stored for SHOT: '+num2str(shot,1)
        ENDIF ELSE BEGIN
		IF NOT keyword_set(quiet) THEN print, 'SHOT '+num2str(shot,1)+' calibrated using '+num2str(ctimes[0],dp=2)+' < t < '+num2str(ctimes[1],dp=2)+' of SHOT '+num2str(cshot,1)
	ENDELSE
END

;+
;NAME:
;	HIREXSR_CHECK_ANALYSIS
;
;PURPOSE:
;	This procedure checks and sees what THACO ANALYSIS# nodes
;	have been populated, by whom, when and for what purpose
;
;CALLING SEQUENCE:
;	HIREXSR_CHECK_ANALYSIS,shot,user,time,note
;	
;INPUTS:
;	shot	LONG	shot number
;
;OPTIONAL INPUTS:
;	max	INT	maximum THACHO# to check DEFAULT: 10
;
;KEYWORD PARAMETERS:
;	quiet	/quiet will supress COMMENT node data being printed to the terminal
;
;OUTPUTS:
;	user	STRARR	[max+2] of the username of the author
;	time	STRARR	[max+2] of the systime() stamp when the analysis was done
;	note	STRARR	[max+2] of note input by the author (quality, purpose...etc)
;	thts	INTARR	[max+2] of the tht#	
;
;OPTIONAL OUTPUTS:
;	tstatus	INTARR	[max+1] indicating the existance of a tree at THACO#
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 8/19/11
;	M.L. Reinke	11/19/2012 - added the automatic check of THT=99 (real-time) and included thts ooutputs
;
;-	

PRO hirexsr_check_analysis,shot,user,time,note,thts,max=max,quiet=quiet,tstatus=tstatus
	IF NOT keyword_set(max) THEN max=10
	tstatus=fltarr(max+2)
	tstatus[0]=1
	thts=int([make(0,max,max+1),99])
	mdsopen,'spectroscopy',shot
	FOR j=1,n(thts) DO BEGIN
		i=int(thts[j])
		path='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS'+num2str(i,1)+'.HELIKE:MORDER'
		chk=mdsvalue('getnci('+path+',"usage")',status=status,/quiet)
		IF status THEN tstatus[j]=1
	ENDFOR
	user=strarr(max+2)
	time=strarr(max+2)
	note=strarr(max+2)
	path='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS:COMMENT'
	comment=mdsvalue(path,/quiet,status=status)
	IF status THEN BEGIN
		IF NOT keyword_set(quiet) THEN print, comment
		split=strsplit(comment,',',/extract)
		user[0]=split[0]
		time[0]=split[1]
		note[0]=split[2]
	ENDIF ELSE BEGIN
		IF NOT keyword_set(quiet) THEN print, 'no COMMENT data'
		user[0]='unknown'
		time[0]='unknown'
		note[0]='unknown'
	ENDELSE
	FOR j=1,n(thts) DO BEGIN
		i=int(thts[j])
		astr='ANALYSIS'+num2str(i,1)
		path='\SPECTROSCOPY::TOP.HIREXSR.'+astr+':COMMENT'
		comment=mdsvalue(path,/quiet,status=status)
		user[j]='none'
		time[j]='none'
		note[j]='none'
		IF tstatus[j] THEN BEGIN
			IF status THEN BEGIN
				IF NOT keyword_set(quiet) THEN print, comment
				split=strsplit(comment,',',/extract)
				user[j]=split[0]
				time[j]=split[1]
				note[j]=split[2]
			ENDIF ELSE IF NOT keyword_set(quiet) THEN print, 'no COMMENT data'
		ENDIF ELSE IF NOT keyword_set(quiet) THEN print, 'no '+astr+' node'
	ENDFOR
	mdsclose,'spectroscopy',shot
END

;+
;NAME:
;	HIREXSR_AVESPEC2TREE
;
;PURPOSE:
;	This procedure loads raw HIREXSR images and applies the TMAP and CHMAP
;	to form averaged spectra for use in FITSPEC2TREE
;
;CALLING SEQUENCE:
;	HIREXSR_AVESPEC2TREE,shot
;	
;INPUTS:
;	shot	LONG	shot number
;
;OPTIONAL INPUTS:
;	nave	INT	number of lambda points over which to average DEFAULT: # of chmap[*,0,0] pts
;	tht	INT	THACO tree # DEFAULT:0 (just ANALYSIS tree)
;
;KEYWORD PARMATERS:
;	h	/h will run the the code for the H-like module rather then for He-like
;	double	/double will run the code in double precisions
;	debug	/debug stops before the end of the code (and the HEAP cleanup)
;	wf	/wf will load white-field stored in tree else wf=1.0 is assumed 
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 11/10
;	8/19/11		M.L. Reinke - adapted to use THT optional input and set the default from 20 to read CHMAP
;				      added the functionality of writing a logname and systime stamp into a comment node
;	1/25/13		M.L. Reinke - added the chmap optional input in call to HIREXSR_AVE_SPEC to
;                                     improve results for non-uniform binning
; 	7/26/13		M.L. Reinke - added the /random keyword for storing data w/ variation for
;				      spectral fitting studies                              			    
;-

PRO hirexsr_avespec2tree,shot,debug=debug,h=h,nave=nave,double=double,tht=tht,note=note,newtht=newtht,wf=wf,random=random
	IF NOT keyword_set(double) THEN double=0
	IF keyword_set(newtht) THEN hirexsr_addnew_analysis,shot,tht=tht,/setup
	IF NOT keyword_set(tht) THEN tht=0
	status=hirexsr_is_analysis(shot,tht)
	IF NOT status THEN BEGIN
		print, 'ANALYSIS'+num2str(tht,1)+' TREE DOES NOT EXIST'
		RETURN
	ENDIF	

	hirexsr_ptr_image,shot,raw,lambda,t,const=const,pos=pos,u=u,h=h,wf=wf
	ntime=n(t)+1
	print, 'raw data loaded to PTRARR'

	;bin images in time
	hirexsr_load_binning,shot,chmap,tch,tmap,good,chmax,h=h,tht=tht
	cnts=hirexsr_bin_image(raw,tmap)
	tau=hirexsr_bin_time(t,tmap)
	tmp=where(tmap EQ -1)
	IF n(tmap) NE n(t) THEN BEGIN
		print, 'TMAP, T are unequal length'
		heap_free,raw
		heap_free,cnts
		heap_gc
		RETURN
	ENDIF
	IF tmp[0] EQ -1 THEN nused=ntime ELSE nused=ntime-(n(tmp)+1)
	print, 'tmap applied: using '+num2str(nused,1)+' of '+num2str(ntime,1)+' frames'
	heap_free,raw

	;form spectra according to chmap
	spec=hirexsr_bin_spec(cnts,tau,chmap,tch,lambda,tmap,const=const,nchbins=chmax,coefs=coefs)
	print, 'chmap applied'
	heap_free,cnts

	;average over nave # of points
	IF NOT keyword_set(nave) THEN BEGIN
		tmp=where(chmap[*,0,0] EQ 0)
		nave=n(tmp)+1
	ENDIF
	avespec= hirexsr_ave_spec(spec,nave,maxave=maxave,chmap=chmap)
	print, 'spectra averaged - nave taken from CHMAP'

	;covert averaged spectra pointer array and write to tree
	avespec_arr=hirexsr_avespec_ptr2arr(avespec,chmax,ntime,maxave)
	specbr=avespec_arr[*,*,*,0]
	lam=avespec_arr[*,*,*,1]
	sig=avespec_arr[*,*,*,2]
	resid=avespec_arr[*,*,*,3]
	IF keyword_set(random) THEN BEGIN
		tmp=where(lam NE -1)
		specbr[tmp]=specbr[tmp]+randomn(seed,n(tmp)+1)*sig[tmp]
	ENDIF
	hirexsr_write_avespec,shot,specbr,lam,sig,resid,tau,nave,h=h,tht=tht
	IF keyword_set(tht) THEN tstr='ANALYSIS'+num2str(tht) ELSE tstr='ANALYSIS'
	print, 'avespec written to '+tstr+' tree SHOT: '+num2str(shot,1)

	IF NOT keyword_set(note) THEN note='x'
	IF keyword_set(random) THEN note='random error added'
	IF strlen(logname()) GT 0 THEN lname=logname() ELSE lname='mdsplus'
	comment=lname+','+systime()+','+note
	print, comment
	hirexsr_write_analysis_comment,shot,comment,tht=tht
	IF keyword_set(debug) THEN stop

	;clean up memory
	heap_free,spec
	heap_free,avespec
	heap_gc
END

;+
;NAME:
;	HIREXSR_FITSPEC2TREE
;
;PURPOSE:
;	This procedure runs the spectral fitting codes to remove neighboring lines
;	and then calculates the moment profiles to be used in the inversions
;
;CALLING SEQUENCE:
;	HIREXSR_FITSPEC2TREE,shot,line
;
;INPUTS:
;	shot	LONG	shot number
;	line	INT	line # for the analysis 0->6 are Ar w,x,z,lya1,Mo 4d, Ar J, Ca w
;
;OPTIONAL INPUTS:
;	nback	INT	0,1,2 of the order of the polynomial to use in fitting the background (currently only works for line=2)
;	tree	STRING	of the EFIT tree to use for calculating the tangency radii DEFAULT: 'ANALYSIS'
;	tht	INT	of the THACO tree # DEFAULT:0 (just \SPECTROSCOPY::TOP.HIREXSR.ANALYSIS)
;
;KEYWORD PARAMETERS:
;	nofit	/nofit will skip the HIREXSR_FIT_SPECTRA step and load from tree (assuming it's been done before)
;		in order to save time when recalculating moments or POS vectors
;	verb	/verb will be sent to HIREXSR_FIT_SPECTRA
;	debug	/debug stops before the end of the code (and the HEAP cleanup)
;
;MODIFICATION HISTORY:	
;	Written by:	M.L. Reinke - 11/10
;	8/19/11		M.L. Reinke - adapted to use THT optional input
;	8/25/12		M.L. Reinke - adjusted CASE statement selection of h & z, added line=7-9 for high-Te layout
;
;-

PRO hirexsr_fitspec2tree,shot,line,verb=verb,debug=debug,nofit=nofit,nback=nback,tree=tree,tht=tht,noinst=noinst
	IF NOT keyword_set(tht) THEN tht=0
	status=hirexsr_is_analysis(shot,tht)
	IF NOT status THEN BEGIN
		print, 'ANALYSIS'+num2str(tht,1)+' TREE DOES NOT EXIST'
		RETURN
	ENDIF	
	double=0	
	CASE line OF 
		0 : BEGIN
			gline='w'
			h=0
			z=18
		END
		1 : BEGIN
			gline='x'
			h=0
			z=18
		END
		2 : BEGIN
			gline='z'
			h=0
			z=18
		END
		3 : BEGIN
			gline='lya1'
			h=1
			z=18
		END
		4 : BEGIN
			gline='4d'
			h=1
			z=42
		END
		5 : BEGIN
			gline='J' 
			h=1
			z=18
		END
		6 : BEGIN
			gline='w'
			h=1
			z=20
		END
		7 : BEGIN
			gline='lya1'	;h-like Ar on the x3 modules
			h=0
			z=18
		END
		8 : BEGIN
			gline='4d'	;Mo 4d on the x3 modules
			h=0
			z=42
		END
		9 : BEGIN
			gline='lya1'	;h-like Ca on module 4
			h=1
			z=20
		END
	ENDCASE
	IF keyword_set(tht) THEN tstr='ANALYSIS'+num2str(tht) ELSE tstr='ANALYSIS'

	hirexsr_load_wavelengths,-1,lam_o,z_o,label_o
	tmp=where(z_o EQ z AND label_o EQ gline)
	lam_o=lam_o[tmp[0]]

	hirexsr_ptr_image,shot,raw,lambda,t,const=const,pos=pos,u=u,iwidth=iwidth,ishift=ishift,h=h,/noimage
	ntime=n(t)+1
	print, 'raw data loaded to PTRARR'

	;load binning
	hirexsr_load_binning,shot,chmap,tch,tmap,good,chmax,h=h,tht=tht
	print, 'binning data loaded'

	;load avespec
	hirexsr_load_avespec,shot,avespec,lam,sig,resid,tau,nave,/ptr,h=h,tht=tht
	print, 'AVESPEC loaded from '+tstr+' tree'

	;fit all spectra for all times
	IF keyword_set(nofit) THEN BEGIN
		IF line EQ 8 OR line EQ 4 THEN hirexsr_load_fits,shot,line-1,coefs,nave,double,labels,/ptr,tht=tht  ELSE hirexsr_load_fits,shot,line,coefs,nave,double,labels,/ptr,tht=tht
		print, 'FITTING SKIPPED - fits loaded from '+tstr+' tree'
	ENDIF ELSE BEGIN
		IF NOT keyword_set(nback) THEN nback=0
		coefs=hirexsr_fit_spectra(avespec,line,double=double,label=labels,verb=verb,nback=nback)
		print, 'spectra fitted - nback='+num2str(nback,1)
	ENDELSE

	;calculate moments
	mom=hirexsr_spec2moments(avespec,coefs,labels,lam_o,line,dlam=dlam,bfrac=bfrac,scale=scale,fitcase=fitcase)
	print, 'moments generated'
	
	;calculate pos for each channel
	chpos=hirexsr_bin_pos(pos,lam_o+dlam*[-1.0e-3,1.0e3],chmap,lambda,chmax,u=u,newu=chu,const=const,newc=chc)

	;calculate instrumentals for each channel
	IF keyword_set(noinst) THEN BEGIN
		iwidth*=0.0
		ishift*=0.0
	ENDIF
	chwidth=hirexsr_bin_inst(iwidth,lam_o+dlam*[-1.0e-3,1.0e3],chmap,lambda,chmax)
	chshift=hirexsr_bin_inst(ishift,lam_o+dlam*[-1.0e-3,1.0e3],chmap,lambda,chmax)
	hirexsr_write_inst,shot,line,chwidth,chshift,tch,tht=tht,status=status			;write ANALYSIS instrumentals to tree
	IF status THEN BEGIN
		tmp=where(chwidth NE -1)
		print, 'average iwidth stored = '+num2str(mean(chwidth[tmp])*1.0e3,dp=1)+' [mAng]'
		print, 'average ishift stored = '+num2str(mean(chshift[tmp])*1.0e3,dp=1)+' [mAng]'
	ENDIF ELSE print, 'no instrumental nodes - run HIREXSR_ADD_INST_NODES to update tree'

	;calculate rhotang for each channel
	rhotang=hirexsr_calc_rhotang(shot,chpos,tch,tau,tree=tree)

	;convert fit coefs pointer array and write to tree
	coefs_arr=hirexsr_coefs_ptr2arr(coefs,chmax,ntime)
	hirexsr_write_fits,shot,line,coefs_arr,tau,nave,double,labels,tht=tht
	print, 'fit coefficients written to '+tstr+' tree SHOT: '+num2str(shot,1)

	;write moments to tree
	hirexsr_write_moments,shot,line,mom[*,*,*,0],mom[*,*,*,1],mom[*,*,*,2],mom[*,*,*,3],tau,chpos,rhotang,bfrac,scale,chu,tch,dlam,double,fitcase,tree,tht=tht
	print, 'moments written to '+tstr+' tree SHOT: '+num2str(shot,1)

	IF keyword_set(debug) THEN stop

	;clean up memory
	heap_free,raw
	heap_free,coefs
	heap_free,avespec
	heap_gc
END

;+
;NAME:
;	HIREXSR_INVERT2TREE
;
;PURPOSE:
;	This procedure loads the moments calculated using FITSPEC2TREE, calculates the spatial and velocity
;	weighting matrices and then does the inversion to find emiss, rotation and temperature profiles.
;
;CALLING SEQUENCE:
;	HIREXSR_INVERT2TREE,shot,line
;
;INPUTS:
;	shot	LONG	shot number
;	line	INT	line # for the analysis 0->6 are Ar w,x,z,lya1,Mo 4d, Ar J, Ca w
;
;OPTIONAL INPUTS:
;	tht	INT	of the THACO tree # DEFAULT:0 (just \SPECTROSCOPY::TOP.HIREXSR.ANALYSIS)
;	good	INTARR	of the moment channels to use for inversion DEFAULT : see HIREXSR_LOAD_MOMENTPTR
;	rmax	FLOAT	of the maximum psinorm for inversion DEFAULT: line specific
;	nr	INT	of the number of (evenly spaced) psinorm points to use DEFAULT: line specific
;	rhovec	FLTARR	of the psinorm points to use for the inversion DEFAULT: uses HIREXSER_MAKE_RHOVEC
;
;KEYWORD PARAMETERS:
;	sine	/sine with run the inversion calculating m=1 terms (up/down asymmetry) in all profiles (use only w/ He-like lines)
;	treevox	/treevox will load voxels from the tree
;	nosub	/nosub will be sent to HIREXSR_CALC_PROFILES and can be used to avoid issues of a bad velocity calibration on Ti profiles
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 11/10
;	1/17/10		M.L. Reinke - added the line specific rhomax and nrho
;	11/20/12	M.L. Reinke - added the GOOD optional input
;    
;-

PRO hirexsr_invert2tree,shot,line,treevox=treevox,nosub=nosub,rmax=rmax,nr=nr,rhovec=rhovec,sine=sine,good=good,tht=tht
	IF NOT keyword_set(tht) THEN tht=0
	status=hirexsr_is_analysis(shot,tht)
	IF NOT status THEN BEGIN
		print, 'ANALYSIS'+num2str(tht,1)+' TREE DOES NOT EXIST'
		RETURN
	ENDIF
	IF keyword_set(tht) THEN tstr='ANALYSIS'+num2str(tht) ELSE tstr='ANALYSIS'
	IF NOT keyword_set(treevox) THEN BEGIN
		hirexsr_load_moments,shot,line,mom,err,pmom,perr,tau,pos,rhotang,bfrac,scale,u,tpos,dlam,double,fitcase,tree,tht=tht
		ntime=n(tau)+1
		tmp=where(tau GT 0)

		;setup default psinorm gridding - make constant for all time points -> use W_HIREXSR_PROFILES to manage time-evolving gridding
		nrho=intarr(ntime)-1
		rhomin=fltarr(ntime)-1.0
		rhomax=fltarr(ntime)-1.0
		rhomin[tmp]=0										;all have default minimum rho of 0.0
		CASE line OF 
			0 : BEGIN
				rhomax[tmp]=1.05 	;'w'
				nrho[tmp]=25
			END
			1 : BEGIN
				rhomax[tmp]=1.05 	;'x'
				nrho[tmp]=25
			END
			2 : BEGIN
				rhomax[tmp]=1.05 	;'z'
				nrho[tmp]=25
			END
			3 : BEGIN
				rhomax[tmp]=0.55 	;'lya1'
				nrho[tmp]=15
			END
			4 : BEGIN
				rhomax[tmp]=0.55 	;'4d'
				nrho[tmp]=15
			END
			5 : BEGIN
				rhomax[tmp]=0.55 	;'J'
				nrho[tmp]=15
			END
			6 : BEGIN
				rhomax[tmp]=0.55 	;'w' for he-like ca
				nrho[tmp]=15
                        END
			7 : BEGIN
				rhomax[tmp]=0.8 	;'lya1' (high-Te)
				nrho[tmp]=20
                        END
			8 : BEGIN
				rhomax[tmp]=0.8 	;'4d' (high Te)
				nrho[tmp]=20
			END
			9 : BEGIN
				rhomax[tmp]=0.55 	;'lya1' for he-like ca
				nrho[tmp]=15
                        END
		ENDCASE
		IF keyword_set(rmax) THEN rhomax[tmp]=rmax						;over ride defaults to set rhomax for all
		IF keyword_set(nr) THEN nrho[tmp]=nr							;over ride defaults to set #rho for all
		IF keyword_set(rhovec) THEN BEGIN					
			x=size(nrho)
			rho=ptrarr(x[1],/allocate_heap)
			FOR i=0,x[1]-1 DO IF nrho[i] GT 0 THEN *rho[i]=rhovec ELSE *rho[i]=-1		;over ride defaults to set custom rhovec for all
		ENDIF ELSE rho=hirexsr_make_rhovec(nrho,rhomin,rhomax)

		hirexsr_pos2voxel,shot,tau,pos,tpos,rho,voxel,velvoxel,/verb,sine=sine,m1svoxel=m1svoxel,m1svelvoxel=m1svelvoxel,tree=tree
		hirexsr_write_voxel,shot,line,voxel,velvoxel,rho,tau,tree,m1svoxel=m1svoxel,m1svelvoxel=m1svelvoxel,tht=tht
		print, 'weighting matrices written to '+tstr+' SHOT: '+num2str(shot,1)
		heap_free,voxel
		heap_free,velvoxel
		heap_free,rho
		IF keyword_set(m1svoxel) THEN heap_free,m1svoxel
		IF keyword_set(m1svelvoxel) THEN heap_free,m1svelvoxel
	ENDIF
	hirexsr_load_momentptr,shot,line,mom,tau,pos,tpos,lam_o,z,good=good,tht=tht,/clear					;load moment data
	hirexsr_load_voxel,shot,line,voxel,velvoxel,rho,tau,tree,/ptr,m1svoxel=m1svoxel,m1svelvoxel=m1svelvoxel,tht=tht		;load voxel data
	hirexsr_load_proconfig,shot,line,eps,eta,/quiet,status=status,tht=tht							;load config data
	IF NOT status THEN BEGIN	
		ntime=n(tau)+1
		eps=fltarr(3,ntime)+1.0
		eta=fltarr(3,ntime)+0.1
	ENDIF
	IF NOT keyword_set(sine) AND keyword_set(m1svoxel) THEN  BEGIN
		heap_free,m1svoxel
		m1svoxel=0
	ENDIF
	IF NOT keyword_set(sine) AND keyword_set(m1svelvoxel) THEN  BEGIN
		heap_free,m1svelvoxel
		m1svelvoxel=0
	ENDIF

	hirexsr_calc_profiles,mom,voxel,velvoxel,eps,eta,lam_o,z,rho,profiles,/solidbody,$
		m1svoxel=m1svoxel,m1svelvoxel=m1svelvoxel,nosub=nosub,sine=sine			;perform inversions
	hirexsr_profile_ptr2arr,profiles,pro_arr,proerr_arr,rho_arr				;convert from PTR to array
	hirexsr_write_profile,shot,line,pro_arr,proerr_arr,rho_arr,tau,tht=tht			;write profiles to the tree
	print, 'profiles written to '+tstr+' SHOT: '+num2str(shot,1)

	hirexsr_extract_check,mom,check,good,subcheck				;extract check and good from momentptr
	hirexsr_write_check,shot,line,check,good,subcheck,tht=tht		;write check/good to tree
	hirexsr_write_proconfig,shot,line,eps,eta,tht=tht			;write eps/eta to tree

	heap_free,mom
	heap_free,voxel
	heap_free,velvoxel
	heap_free,profiles
	heap_free,rho
	IF keyword_set(m1svoxel) THEN heap_free,m1svoxel
	IF keyword_set(m1svelvoxel) THEN heap_free,m1svelvoxel
	heap_gc


END

;+
;NAME:
;	HIREXSR_TREE_ACTIONS
;
;Written by:
;	Written by:	M.L. Reinke - 10/10
;	1/17/11 	M.L. Reinke - added functionality for lya1 and mo4d	
;-

PRO hirexsr_tree_actions,on=on,off=off,x=x,z=z,w=w,lya1=lya1,mo4d=mo4d
	IF NOT keyword_set(on) AND NOT keyword_set(off) THEN BEGIN
		print, 'set /on or /off'
		RETURN
	ENDIF
	IF NOT keyword_set(w) AND NOT keyword_set(x) AND NOT keyword_set(z) AND NOT keyword_set(lya1) AND NOT keyword_set(mo4d) THEN BEGIN
		w=1
		x=1
		z=1
		lya1=1
		mo4d=1
	ENDIF
	
	wnodes=['\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.WN3:FITSPEC',$
	       '\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.PROFILES.W:INVERT']
	xnodes=['\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.XY:FITSPEC',$
	       '\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.PROFILES.X:INVERT']
	znodes=['\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.ZJK:FITSPEC',$
               '\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.PROFILES.Z:INVERT']
	lya1nodes=['\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.FITS.LYA:FITSPEC',$
               '\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.PROFILES.LYA1:INVERT']
	mo4dnodes=['\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.FITS.LYA:FITSPEC',$
               '\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.PROFILES.MO4D:INVERT']

	mdstcl, 'set tree spectroscopy/shot=-1'
	IF keyword_set(w) THEN BEGIN
		nodes=wnodes
		FOR i=0,n(nodes) DO BEGIN
			IF keyword_set(on) THEN mdstcl, 'set node '+nodes[i]+' /on'
			IF keyword_set(off) THEN mdstcl, 'set node '+nodes[i]+' /off'
		ENDFOR
	ENDIF
	IF keyword_set(x) THEN BEGIN
		nodes=xnodes
		FOR i=0,n(nodes) DO BEGIN
			IF keyword_set(on) THEN mdstcl, 'set node '+nodes[i]+' /on'
			IF keyword_set(off) THEN mdstcl, 'set node '+nodes[i]+' /off'
		ENDFOR
	ENDIF
	IF keyword_set(z) THEN BEGIN
		nodes=znodes
		FOR i=0,n(nodes) DO BEGIN
			IF keyword_set(on) THEN mdstcl, 'set node '+nodes[i]+' /on'
			IF keyword_set(off) THEN mdstcl, 'set node '+nodes[i]+' /off'
		ENDFOR
	ENDIF
	IF keyword_set(lya1) THEN BEGIN
		nodes=lya1nodes
		FOR i=0,n(nodes) DO BEGIN
			IF keyword_set(on) THEN mdstcl, 'set node '+nodes[i]+' /on'
			IF keyword_set(off) THEN mdstcl, 'set node '+nodes[i]+' /off'
		ENDFOR
	ENDIF
	IF keyword_set(mo4d) THEN BEGIN
		nodes=mo4dnodes
		FOR i=0,n(nodes) DO BEGIN
			IF keyword_set(on) THEN mdstcl, 'set node '+nodes[i]+' /on'
			IF keyword_set(off) THEN mdstcl, 'set node '+nodes[i]+' /off'
		ENDFOR
	ENDIF
	
	mdstcl, 'close'
END

;
; HIREXSR_SAFE_AVESPEC2TREE
;

PRO hirexsr_safe_avespec2tree,shot,verb=verb,debug=debug,h=h,nave=nave,double=double, threshold = threshold

shot = LONG(shot)

if hirexsr_is_plasma(shot, threshold = threshold) then begin
	hirexsr_avespec2tree,shot,debug=debug,h=h,nave=nave,double=double
	print, 'hirexsr_avespec2tree run'
endif else begin
	print, 'NO PLASMA'
endelse

end

;
; HIREXSR_SAFE_FITSPEC2TREE
;


PRO hirexsr_safe_fitspec2tree,shot,line,verb=verb,debug=debug,nofit=nofit, time_thresh= time_thresh, lev_thresh = lev_thresh, volt_thresh = volt_thresh, threshold = threshold

shot = LONG(shot)

if hirexsr_is_plasma(shot,threshold = threshold ) then begin
	if hirexsr_ar_exist(shot, time_thresh = time_thresh, lev_thresh = lev_thresh, volt_thresh = volt_thresh) then begin
		hirexsr_fitspec2tree,shot,line,verb=verb,debug=debug,nofit=nofit
	endif else begin
		print, 'NO ARGON'
	endelse
endif else begin
	print, 'NO PLASMA'
endelse

end

;
; HIREXSR_SAFE_INVERT2TREE
;


PRO hirexsr_safe_invert2tree,shot,line,treevox=treevox,nosub=nosub, time_thresh= time_thresh, lev_thresh = lev_thresh, volt_thresh = volt_thresh, threshold = threshold

shot = LONG(shot)

if hirexsr_is_plasma(shot, threshold= threshold) then begin
	if hirexsr_ar_exist(shot, time_thresh = time_thresh, lev_thresh = lev_thresh, volt_thresh = volt_thresh) then begin
		hirexsr_invert2tree, shot,line,treevox=treevox,nosub=nosub
		print, 'hirexsr_avespec2tree run'
	endif else begin
		print, 'NO ARGON'
	endelse
endif else begin
	print, 'NO PLASMA'
endelse


end

;+
;NAME:
;	HIREXSR_CHANGE_CRYSTAL_BETA
;
;-

PRO hirexsr_change_crystal_beta,shot,beta,h=h,nomom=nomom,ca=ca,invline=invline,morder=morder
	IF NOT keyword_set(morder) THEN hirexsr_load_morder,shot,morder,h=h
	IF keyword_set(h) OR keyword_set(ca) THEN BEGIN
		IF keyword_set(ca) THEN lines=[6,7,8] ELSE lines=[3,4,5]
	ENDIF ELSE BEGIN
		lines=[0,1,2]
	ENDELSE
	FOR i=0,n(morder) DO BEGIN
		info=hirexsr_load_info(morder[i],shot=shot,/tree)
		info.m.rot[1]=beta						;set new declination angle
		hirexsr_load_info2tree,shot,morder[i],info=info
		hirexsr_write_pos,shot,morder[i]				;refill pos w/ modified info file
		print, 'beta='+num2str(beta,dp=3)+' [rad] written to MOD'+num2str(morder[i],1)+' INFO file and POS updated'
	ENDFOR

	IF NOT keyword_set(nomom) THEN BEGIN
		IF NOT keyword_set(h) THEN BEGIN						;load wavelengths
			hirexsr_load_lambda,shot,morder[0],lambda1
			hirexsr_load_lambda,shot,morder[1],lambda2
 			hirexsr_load_lambda,shot,morder[2],lambda3
			lambda=[lambda1,lambda2,lambda3]
			hirexsr_load_pos,shot,morder[0],pos1
			hirexsr_load_pos,shot,morder[1],pos2
			hirexsr_load_pos,shot,morder[2],pos3
			pos=[[pos1],[pos2],[pos3]]
		ENDIF ELSE BEGIN
			hirexsr_load_lambda,shot,morder[0],lambda
			hirexsr_load_pos,shot,morder[0],pos
		ENDELSE
		hirexsr_load_wavelengths,-1,lam_o,z_o,label_o					;get lam_o
		hirexsr_load_binning,shot,chmap,tch,tmap,good,chmax,h=h
		FOR i=0,n(lines) DO BEGIN
			line=lines[i]
			CASE line OF 
				0 : gline='w'
				1 : gline='x'
				2 : gline='z'
				3 : gline='lya1'
				4 : gline='4d'
				5 : gline='J'
			ENDCASE
			z=18
			IF line EQ 4 THEN z=42
			tmp=where(z_o EQ z AND label_o EQ gline)
			lam0=lam_o[tmp[0]]
			hirexsr_load_moments,shot,line,mom,err,pmom,perr,tau,chpos,rhotang,bfrac,scale,u,tpos,dlam,double,fitcase	;get dlam and moment data for rewritting
	
			;calculate new pos for each channel
			chpos=hirexsr_bin_pos(pos,lam0+dlam*[-1.0e-3,1.0e3],chmap,lambda,chmax)
			rhotang=hirexsr_calc_rhotang(shot,chpos,tch,tau)
			hirexsr_write_moments,shot,line,mom,err,pmom,perr,tau,chpos,rhotang,bfrac,scale,u,tpos,dlam,double,fitcase
			print, 'POS vectors updated for LINE='+num2str(line,1)
		ENDFOR
	ENDIF
	IF keyword_set(invline) THEN FOR i=0,n(invline) DO hirexsr_invert_data,shot,invline[i]
END

;+
;NAME:
;	HIREXSR_CHANGE_CRYSTAL_GAMMA
;
;PURPOSE:
;	This procedure writes a new toroidal viewing angle for the
;	crystal to the INFO files and updates the POS
;
;-

PRO hirexsr_change_crystal_gamma,shot,gamma,h=h,morder=morder
	IF NOT keyword_set(morder) THEN hirexsr_load_morder,shot,morder,h=h
	FOR i=0,n(morder) DO BEGIN
		info=hirexsr_load_info(morder[i],shot=shot,/tree)
		info.m.rot[2]=gamma						;set new toroidal angle
		hirexsr_load_info2tree,shot,morder[i],info=info
		hirexsr_write_pos,shot,morder[i]				;refill pos w/ modified info file
		print, 'gamma='+num2str(gamma,dp=3)+' [rad] written to MOD'+num2str(morder[i],1)+' INFO file and POS updated'
	ENDFOR
END

;+
;NAME:
;	HIREXSR_CHANGE_CRYSTAL_DECL
;
;-

PRO hirexsr_change_crystal_decl,shot,dpsi,h=h,nomom=nomom,ca=ca,invline=invline,morder=morder
	IF NOT keyword_set(morder) THEN hirexsr_load_morder,shot,morder,h=h
	IF keyword_set(h) OR keyword_set(ca) THEN BEGIN
		IF keyword_set(ca) THEN lines=[6,7,8] ELSE lines=[3,4,5]
	ENDIF ELSE BEGIN
		lines=[0,1,2]
	ENDELSE
	FOR i=0,n(morder) DO BEGIN
		info=hirexsr_load_info(morder[i],shot=shot,/tree)
		info.m.rot[1]=0.0						;set beta=0.0
		hirexsr_load_info2tree,shot,morder[i],info=info
		hirexsr_write_pos,shot,morder[i]				;write original POS to tree
		hirexsr_load_pos,shot,morder[i],pos
		pos[3,*,*]+=dpsi						;shift view by small angle
		hirexsr_write_pos,shot,morder[i],pos=pos			;refill pos w/ modified info file
		print, 'beta='+num2str(0.0,dp=3)+' [rad] written to MOD'+num2str(morder[i],1)+' INFO file and POS updated with dpsi='+num2str(dpsi,dp=3)+' [rad]'
	ENDFOR

	IF NOT keyword_set(nomom) THEN BEGIN
		IF NOT keyword_set(h) THEN BEGIN						;load wavelengths
			hirexsr_load_lambda,shot,morder[0],lambda1
			hirexsr_load_lambda,shot,morder[1],lambda2
 			hirexsr_load_lambda,shot,morder[2],lambda3
			lambda=[lambda1,lambda2,lambda3]
			hirexsr_load_pos,shot,morder[0],pos1
			hirexsr_load_pos,shot,morder[1],pos2
			hirexsr_load_pos,shot,morder[2],pos3
			pos=[[pos1],[pos2],[pos3]]
		ENDIF ELSE BEGIN
			hirexsr_load_lambda,shot,morder[0],lambda
			hirexsr_load_pos,shot,morder[0],pos
		ENDELSE
		hirexsr_load_wavelengths,-1,lam_o,z_o,label_o					;get lam_o
		hirexsr_load_binning,shot,chmap,tch,tmap,good,chmax,h=h
		FOR i=0,n(lines) DO BEGIN
			line=lines[i]
			CASE line OF 
				0 : gline='w'
				1 : gline='x'
				2 : gline='z'
				3 : gline='lya1'
				4 : gline='4d'
				5 : gline='J'
			ENDCASE
			z=18
			IF line EQ 4 THEN z=42
			tmp=where(z_o EQ z AND label_o EQ gline)
			lam0=lam_o[tmp[0]]
			hirexsr_load_moments,shot,line,mom,err,pmom,perr,tau,chpos,rhotang,bfrac,scale,u,tpos,dlam,double,fitcase	;get dlam and moment data for rewritting
	
			;calculate new pos for each channel
			chpos=hirexsr_bin_pos(pos,lam0+dlam*[-1.0e-3,1.0e3],chmap,lambda,chmax)
			rhotang=hirexsr_calc_rhotang(shot,chpos,tch,tau)
			hirexsr_write_moments,shot,line,mom,err,pmom,perr,tau,chpos,rhotang,bfrac,scale,u,tpos,dlam,double,fitcase
			print, 'POS vectors updated for LINE='+num2str(line,1)
		ENDFOR
	ENDIF
	IF keyword_set(invline) THEN FOR i=0,n(invline) DO hirexsr_invert_data,shot,invline[i]
END

;fixes a bug in the original deployment that forgot a MO4D VOXEL node
PRO hirexsr_fixmo4d_tree,shot,quiet=quiet,tht=tht
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht) ELSE astr='ANALYSIS'
	mdsopen,'spectroscopy',shot
	chk=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.PROFILES.MO4D.CONFIG:VOXEL',/quiet,status=status)
	mdsclose,'spectroscopy',shot	
	IF NOT keyword_set(tht) THEN ntht=0 ELSE ntht=tht

	IF NOT status THEN BEGIN
		mdstcl, "set verify"
		mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
		mdstcl, 'add node /usage=numeric \SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.PROFILES.MO4D.CONFIG:VOXEL'
		mdstcl, 'write'
		mdstcl, 'close'
	ENDIF ELSE IF NOT keyword_set(quiet) THEN print, 'SHOT: '+num2str(shot,1)+' THT='+num2str(ntht)+' has MO4D VOXEL node already'
END

;adds a COMMENT node for the default ANALYSIS or tht
PRO hirexsr_add_analysis_comment,shot,quiet=quiet,tht=tht
	IF NOT keyword_set(tht) THEN astr='ANALYSIS' ELSE astr='ANALYSIS'+num2str(tht,1)
	path='\SPECTROSCOPY::TOP.HIREXSR.'+astr+':COMMENT'
	mdsopen,'spectroscopy',shot
	chk=mdsvalue('getnci('+path+',"usage")',status=status,/quiet)
	mdsclose,'spectroscopy',shot
	
	IF NOT status THEN BEGIN
		mdstcl, "set verify"
		mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
		mdstcl, 'add node /usage=text \SPECTROSCOPY::TOP.HIREXSR.'+astr+':COMMENT'
		mdstcl, 'write'
		mdstcl, 'close'
	ENDIF ELSE IF NOT keyword_set(quiet) THEN print, 'SHOT: '+num2str(shot,1)+' has ANALYSIS:COMMENT node already'
END

;adds INST nodes to CALIB and the selected ANALYSIS (DEFAULT: tht=0)
PRO hirexsr_add_inst_nodes,shot,quiet=quiet,tht=tht
	;write the INST nodes into the CALIB tree
	path='\SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD1:INST'
	mdsopen,'spectroscopy',shot
	chk=mdsvalue('getnci('+path+',"usage")',status=status,/quiet)
	mdsclose,'spectroscopy',shot
	IF NOT status THEN BEGIN
		mdstcl, "set verify"
		mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
		mdstcl, 'add node /usage=signal \SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD1:INST'
		mdstcl, 'add node /usage=signal \SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD2:INST'
		mdstcl, 'add node /usage=signal \SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD3:INST'
		mdstcl, 'add node /usage=signal \SPECTROSCOPY::TOP.HIREXSR.CALIB.MOD4:INST'
		mdstcl, 'write'
		mdstcl, 'close'
	ENDIF ELSE IF NOT keyword_set(quiet) THEN print, 'SHOT: '+num2str(shot,1)+' has CALIB INST nodes already'

	;write the INST nodes into the ANALYSIS tree
	IF NOT keyword_set(tht) THEN astr='ANALYSIS' ELSE astr='ANALYSIS'+num2str(tht,1)
	path='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.W:INST'
	mdsopen,'spectroscopy',shot
	chk=mdsvalue('getnci('+path+',"usage")',status=status,/quiet)
	mdsclose,'spectroscopy',shot
	IF NOT status THEN BEGIN
		mdstcl, "set verify"
		mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
		mdstcl, 'add node /usage=signal \SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.W:INST'
		mdstcl, 'add node /usage=signal \SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.X:INST'
		mdstcl, 'add node /usage=signal \SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.MOMENTS.Z:INST'
		mdstcl, 'add node /usage=signal \SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.J:INST'
		mdstcl, 'add node /usage=signal \SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.LYA1:INST'
		mdstcl, 'add node /usage=signal \SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.MOMENTS.MO4D:INST'
		mdstcl, 'write'
		mdstcl, 'close'
	ENDIF ELSE IF NOT keyword_set(quiet) THEN print, 'SHOT: '+num2str(shot,1)+' has '+astr+' MOMENTS INST nodes already'
END
