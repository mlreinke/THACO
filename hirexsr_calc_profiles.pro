;+
;NAME:
;	HIREXSR_MAKE_RHOVEC
;
;PURPOSE:
;	
;
;-

FUNCTION hirexsr_make_rhovec,nrho,rhomin,rhomax
	x=size(nrho)
	rho=ptrarr(x[1],/allocate_heap)
	FOR i=0,x[1]-1 DO IF nrho[i] GT 0 THEN *rho[i]=make(rhomin[i],rhomax[i],nrho[i]) ELSE *rho[i]=-1
	RETURN,rho
END

;+
;NAME:
;	HIREXSR_UPDATE_VOXEL
;
;PURPOSE:
;	This procedure can be used to update specific voxel weightings
;
;-

PRO hirexsr_update_voxel,voxel,velvoxel,index,shot,tau,pos,tpos,rho,n_s=n_s,r_ap=r_ap,tree=tree
	IF NOT keyword_set(r_ap) THEN r_ap=0.95
	IF NOT keyword_set(n_s) THEN n_s=500
	FOR j=0,n(index) DO BEGIN
		i=index[j]
		irho=*rho[i]
		pindex=last(where(tpos LE tau[i]))
		ipos=pos[*,*,pindex]
		tmp=where(pos[2,*] GE 0)
		iu=fltarr(n(tmp)+1)+4.0*!pi			;this puts voxel in units of length so brightness can be inverted rather than power
		ivoxel=genpos_pos2voxel_matrix(ipos[*,tmp],iu,shot,tpts=tau[i],rho_vec=irho,r_ap=r_ap,n_s=n_s,rhopts=rhopts,/psinorm,tree=tree)
		*voxel[i]=ivoxel
		ivelvoxel=genspec_pos2voxel_vel_matrix(ipos[*,tmp],iu,shot,tpts=tau[i],rho_vec=irho,r_ap=r_ap,n_s=n_s,/psinorm,tree=tree)
		*velvoxel[i]=reform(ivelvoxel)
	ENDFOR
END

;+
;NAME:
;	HIREXSR_POS2VOXEL
;
;PURPOSE:
;	This function creates the voxel PTRARR
;
;MODIFICATION HISTORY:
;	6/13/11 - added the ability to specific EFIT tree
;-

PRO hirexsr_pos2voxel,shot,tau,pos,tpos,rho,voxel,velvoxel,n_s=n_s,r_ap=r_ap,verb=verb,index=index,m1svoxel=m1svoxel,m1svelvoxel=m1svelvoxel,sine=sine,tree=tree
	IF NOT keyword_set(r_ap) THEN r_ap=0.95
	IF NOT keyword_set(n_s) THEN n_s=500

	ntime=n(tau)+1
	voxel=ptrarr(ntime,/allocate_heap)
	velvoxel=ptrarr(ntime,/allocate_heap)
	IF keyword_set(sine) THEN BEGIN
		m1svoxel=ptrarr(ntime,/allocate_heap)
		m1svelvoxel=ptrarr(ntime,/allocate_heap)
	ENDIF
	npos=n(tpos)+1
	start_time=systime(/seconds)
	IF keyword_set(index) THEN BEGIN
		start=index
		stop=index
	ENDIF ELSE BEGIN
		start=0
		stop=ntime-1
	ENDELSE	

	;runs p2v for each time slice seperately
	FOR i=start,stop DO BEGIN
		irho=*rho[i]
		IF irho[0] NE -1 THEN BEGIN
			pindex=last(where(tpos LE tau[i]))
			ipos=pos[*,*,pindex]
			tmp=where(ipos[2,*] GE 0)
			iu=fltarr(n(tmp)+1)+4.0*!pi			;this puts voxel in units of length so brightness can be inverted rather than power
			ivoxel=genpos_pos2voxel_matrix(ipos[*,tmp],iu,shot,tpts=tau[i],rho_vec=irho,r_ap=r_ap,n_s=n_s,/psinorm,tree=tree)
			*voxel[i]=ivoxel
			ivelvoxel=genspec_pos2voxel_vel_matrix(ipos[*,tmp],iu,shot,tpts=tau[i],rho_vec=irho,r_ap=r_ap,n_s=n_s,/psinorm,tree=tree)
			*velvoxel[i]=reform(ivelvoxel)
			IF keyword_set(sine) THEN BEGIN
				im1svoxel=genpos_pos2voxel_matrix(ipos[*,tmp],iu,shot,tpts=tau[i],rho_vec=irho,r_ap=r_ap,n_s=n_s,/psinorm,m=1,/sine,tree=tree)
				*m1svoxel[i]=im1svoxel
				im1svelvox=genspec_pos2voxel_vel_matrix(ipos[*,tmp],iu,shot,tpts=tau[i],rho_vec=irho,r_ap=r_ap,n_s=n_s,/psinorm,m=1,/sine,tree=tree)
				*m1svelvoxel[i]=reform(im1svelvox)
			ENDIF
		ENDIF ELSE BEGIN
			*voxel[i]=-1
			*velvoxel[i]=-1
			IF keyword_set(sine) THEN BEGIN
				*m1svoxel[i]=-1
				*m1svelvoxel[i]=-1
			ENDIF
		ENDELSE
	ENDFOR
	ctime=systime(/seconds)-start_time
	IF tree EQ 'ANALYSIS' THEN treestr='ANALYSIS EFIT' ELSE treestr=strupcase(tree)
	IF keyword_set(verb) THEN print, 'VOXELs calculated: '+num2str(ctime,dp=2)+' using '+treestr
END



;+
;NAME:
;	HIREXSR_INVERT_MOMENTS
;
;
;MODIFICATION HISTORY:
;	4/4/11		M.L. Reinke - included the m=1 terms in the inversions but haven't included noise propigation
;	4/6/11		M.L. Reinke - included the m=1 error propigation
;	5/12/11		M.L. Reinke - modified to work with newer version of GENPOS_PROFILE_INVERT that accepts inv_matrix as an optional input
;	9/2/11		M.L. Reinke - added the ability to insert a non-uniform instrumental Ti [keV] to subtract off prior to inversion
;				      can see the inst_ti_vec and sub_vec as optional outputs
;	9/4/11		M.L. Reinke - modified handling of instrumentals to work with INST nodes in tree
;-

PRO hirexsr_invert_moments,imom,ivox,ivelvox,eps,eta,imomerr,igood,lam_o,z,ipro,iproerr,icheck,isubcheck,nofirst=nofirst,solidbody=solidbody,ipmom=ipmom,nosub=nosub,$
		im1svox=im1svox,im1svelvox=im1svelvox,iwidth=iwidth,ishift=ishift

	c=3.0e8  			;speed of light
	e=1.60e-19			;conversion for eV -> J
	mconv=1.66e-27			;conversion for amu -> kg
	mass=read_atomic_mass(z)
	
	x=size(ivox)
	nch=x[1]
	nrho=x[2]
	IF keyword_set(im1svox) THEN BEGIN
		ipro=fltarr(2*nrho,4)
		iproerr=fltarr(2*nrho,4)
		nprof=2
	ENDIF ELSE BEGIN
		ipro=fltarr(nrho,4)
		iproerr=fltarr(nrho,4)
		nprof=1
	ENDELSE
	icheck=fltarr(nch,3)				;store the comparison of the derived profiles with the input moment data
	isubcheck=fltarr(nch,3)				;store the vectors subtracted off the moment profiles to instrumentals and non-zero velocity on Ti
	tmp=where(igood EQ 1)
	IF keyword_set(ipmom) THEN BEGIN
		totphot=ipmom[*,0]
		mu=ipmom[*,1]
		width=ipmom[*,2]
		bfrac=ipmom[*,3]
		bad=where(finite(width) EQ 0)
		IF bad[0] NE -1 THEN width[bad]=sqrt(abs(imom[bad,2]/imom[bad,0]))	;
		scale=ipmom[*,4]
		bckphot=totphot*bfrac
		bad=where(totphot-bckphot LT 2)
		IF bad[0] NE -1 THEN bckphot[bad]=totphot[bad]-2
		sigphot=totphot-bckphot
	ENDIF

	;inversion of 0th moment
	a=max(ivox)
	IF keyword_set(im1svox) THEN ivoxtot=[[ivox[tmp,*]],[im1svox[tmp,*]]] ELSE ivoxtot=ivox[tmp,*]
	ipro[*,0]=genpos_profile_invert(imom[tmp,0],ivoxtot,double(eps[0]*a^2),brchk=br0,nofirst=nofirst,eta=double(eta[0]*a^2),err=imomerr[tmp,0],inv_matrix=inv_matrix0,nprof=nprof)
	icheck[tmp,0]=br0.mom
	iproerr[*,0]=br0.inverr
	dmat=inv_matrix0#transpose(ivoxtot)

	;inversion of 1st moment
	IF keyword_set(solidbody) THEN BEGIN			;invert the 1st moment profile assuming only solidbody toroidal rotation
		IF keyword_set(ishift) THEN BEGIN
			ivox_inst=ivox
			FOR i=0,nch-1 DO ivox_inst[i,*]*=ishift[i]
			inst_v_vec=ivox_inst#ipro[0:nrho-1,0]			;derive the instrument vector to subtract off before inversion
		ENDIF ELSE inst_v_vec=fltarr(nch)
		a=max(ivelvox[*,*,0])
		IF keyword_set(im1svelvox) THEN ivelvox_tot=[[ivelvox[tmp,*,0]],[im1svelvox[tmp,*,0]]] ELSE ivelvox_tot=ivelvox[tmp,*,0]
		w_vel_inv=genpos_profile_invert(reform(imom[tmp,1])-inst_v_vec[tmp],ivelvox_tot,double(eps[1]*a^2),brchk=br1,nofirst=nofirst,eta=double(eta[1]*a^2),err=imomerr[tmp,1],inv_matrix=inv_matrix1)
		emat=inv_matrix1#transpose(ivelvox_tot)
		w_err=br1.inverr             
		u_vel_inv=fltarr(nrho)		;set poloidal term to zero
		u_err=fltarr(nrho)
	ENDIF ELSE BEGIN					;invert assuming both toroidal and poloidal
		inst_v_vec=fltarr(nch)
		ivelvox_tot=[[ivelvox[*,*,0]],[ivelvox[*,*,1]]]
		a=max(abs(ivelvox_tot))
		vel_inv=genpos_profile_invert(imom[tmp,1]-inst_v_vec[tmp],ivelvox_tot[tmp,*],double(eps[1]*a^2),brchk=br1,nofirst=nofirst,eta=double(eta[1]*a^2),err=imomerr[tmp,1],inv_matrix=inv_matrix1)
		emat=inv_matrix1#transpose(ivelvox_tot[tmp,*])
		w_vel_inv=vel_inv[0:nrho-1]
		w_err=br1.inverr[0:nrho-1]		;w profile and error
		u_vel_inv=vel_inv[nrho:*]
		u_err=br1.inverr[nrho:*]		;u profile and error                          
	ENDELSE
	vtconv=lam_o/c*2.0*!pi*1.0e3
	IF keyword_set(im1svelvox) THEN BEGIN
		ipro[0:nrho-1,1]=w_vel_inv/ipro[0:nrho-1,0]/vtconv	
		ipro[nrho:*,1]=(w_vel_inv[nrho:*]-ipro[nrho:*,0]*w_vel_inv[0:nrho-1]/ipro[0:nrho-1,0])/ipro[0:nrho-1,0]/vtconv
		ipro[*,2]=0.0
	ENDIF ELSE BEGIN
		ipro[0:nrho-1,1]=w_vel_inv/(lam_o/c)/ipro[*,0]/(2.0*!pi)/1.0e3					;rotation frequency in [kHz]		
		ipro[0:nrho-1,2]=u_vel_inv/(lam_o/c)/ipro[*,0]/1.0e3						;rotation in km/s/T (u*B = vpol)
	ENDELSE
	icheck[tmp,1]=br1.mom
	isubcheck[*,0]=inst_v_vec
	IF keyword_set(nosub) THEN sub_vec=fltarr(n(imom[*,2])+1) ELSE sub_vec=genspec_matrix_sub_vector(ivox,ivelvox[*,*,0],w_vel_inv[0:nrho-1],ipro[0:nrho-1,0])
	IF keyword_set(ipmom) THEN BEGIN
		IF NOT keyword_set(solidbody) THEN BEGIN
			;get to this eventually
		ENDIF ELSE BEGIN
			v_var_n=fltarr(nrho)
			v_var_mu=fltarr(nrho)
			v_var_n1=fltarr(nrho)
			v_var_mu1=fltarr(nrho)	
			FOR i=0,nrho-1 DO BEGIN
				k=i+nrho	;index for matrices for the m=1 profiles
				v_var_n[i]+=total((scale[tmp,*]*emat[i,*]*imom[tmp,1]/imom[tmp,0]/ipro[i,0]-w_vel_inv[i]/ipro[i,0]^2*dmat[i,*]*scale[tmp,*])^2*(sigphot[tmp]))
				v_var_mu[i]+=total((emat[i,*]*imom[tmp,0]/ipro[i,0])^2*width[tmp]^2/sigphot[tmp]*(1.0+bckphot[tmp]/sigphot[tmp]))
				IF keyword_set(im1svelvox) THEN BEGIN
					v_var_n1[i]+=total((scale[tmp,*]*emat[k,*]*imom[tmp,1]/imom[tmp,0]/ipro[i,0]-w_vel_inv[i]/ipro[i,0]^2*dmat[k,*]*scale[tmp,*]-$
						ipro[k,0]/ipro[i,0]*(scale[tmp,*]*emat[i,*]*imom[tmp,1]/imom[tmp,0]/ipro[i,0]-w_vel_inv[i]/ipro[i,0]^2*dmat[i,*]*scale[tmp,*])-$
						(w_vel_inv[k]-ipro[k,0]*w_vel_inv[i]/ipro[i,0])/(ipro[i,0]^2)*dmat[i,*]*scale[tmp,*])^2*sigphot[tmp])
					v_var_mu1[i]+=total((emat[k,*]*imom[tmp,0]/ipro[i,0]-ipro[k,0]/ipro[i,0]*emat[i,*]*imom[tmp,0]/ipro[i,0])^2*$
						width[tmp]^2/sigphot[tmp]*(1.0+bckphot[tmp]/sigphot[tmp]))
				ENDIF
			ENDFOR
			v_err=sqrt(v_var_n+v_var_mu)
			iproerr[0:nrho-1,1]=v_err/vtconv
			IF keyword_set(im1svelvox) THEN BEGIN
			
				v_err1=sqrt(v_var_n1+v_var_mu1)
				iproerr[nrho:*,1]=v_err1/vtconv
			ENDIF
		ENDELSE
	ENDIF ELSE BEGIN
		iproerr[0:nrho-1,1]=ipro[*,1]*sqrt(w_err^2/w_vel_inv^2+iproerr[*,0]^2/ipro[*,0]^2)
		IF NOT keyword_set(solidbody) THEN iproerr[*,2]=ipro[*,2]*sqrt(u_err^2/u_vel_inv^2+iproerr[*,0]^2/ipro[*,0]^2)
	ENDELSE
	
	;inversion of 2nd moment
	a=max(ivox)
	ticonv=(lam_o/c)^2*(e*1.0e3/(mass*mconv))	
	IF keyword_set(im1svox) THEN ivoxtot=[[ivox[tmp,*]],[im1svox[tmp,*]]] ELSE ivoxtot=ivox[tmp,*]
	IF keyword_set(iwidth) THEN BEGIN
		ivox_inst=ivox
		FOR i=0,nch-1 DO ivox_inst[i,*]*=iwidth[i]^2
		inst_ti_vec=ivox_inst#ipro[0:nrho-1,0]			;derive the instrument vector to subtract off before inversion
	ENDIF ELSE inst_ti_vec=fltarr(nch)

	ti_inv=genpos_profile_invert(imom[tmp,2]-sub_vec[tmp]-inst_ti_vec[tmp],ivoxtot,double(eps[2]*a^2),brchk=br2,nofirst=nofirst,eta=double(eta[2]*a^2),err=imomerr[tmp,2],$
		inv_matrix=inv_matrix2,nprof=nprof)
	fmat=inv_matrix2#transpose(ivoxtot)
	icheck[tmp,2]=br2.mom+sub_vec[tmp]+inst_ti_vec[tmp]
	isubcheck[*,1]=inst_ti_vec
	isubcheck[*,2]=sub_vec
	IF keyword_set(im1svox) THEN BEGIN
		ipro[0:nrho-1,3]=ti_inv[0:nrho-1]/ipro[0:nrho-1,0]/ticonv
		ipro[nrho:*,3]=(ti_inv[nrho:*]-ipro[nrho:*,0]*ti_inv[0:nrho-1]/ipro[0:nrho-1,0])/ipro[0:nrho-1,0]/ticonv
	ENDIF ELSE ipro[*,3]=ti_inv/((lam_o/c)^2*(e*1.0e3/(mass*mconv)))/ipro[*,0]			;Ti in [keV]
	IF keyword_set(ipmom) THEN BEGIN
		ti_var_n=fltarr(nrho)
		ti_var_mu=fltarr(nrho)
		ti_var_w=fltarr(nrho)
		ti_var_n1=fltarr(nrho)
		ti_var_mu1=fltarr(nrho)
		ti_var_w1=fltarr(nrho)
		FOR i=0,nrho-1 DO BEGIN
			ti_var_n[i]+=total((scale[tmp,*]*fmat[i,*]*imom[tmp,2]/imom[tmp,0]/ipro[i,0]-ti_inv[i]/ipro[i,0]^2*dmat[i,*]*scale[tmp,*])^2*sigphot[tmp])		;dTi/dN
			ti_var_mu[i]+=total((2.0*fmat[i,*]*imom[tmp,1]/ipro[i,0])^2*width[tmp]^2/sigphot[tmp]*(1.0+bckphot[tmp]/sigphot[tmp]))					;dTi/dmu	
			ti_var_w[i]+=total((2.0*fmat[i,*]*width[tmp]*imom[tmp,0]/ipro[i,0])^2*width[tmp]^2/(2.0*(sigphot[tmp]))*(1.0+bckphot[tmp]/sigphot[tmp]))		;dTi/dw
			IF keyword_set(im1svox) THEN BEGIN
				k=i+nrho	;index for matrices for the m=1 profiles
				ti_var_n1[i]+=total((scale[tmp,*]*fmat[k,*]*imom[tmp,2]/imom[tmp,0]/ipro[i,0]-ti_inv[i]/ipro[i,0]^2*dmat[k,*]*scale[tmp,*]-$
					ipro[k,0]/ipro[i,0]*(scale[tmp,*]*fmat[i,*]*imom[tmp,2]/imom[tmp,0]/ipro[i,0]-ti_inv[i]/ipro[i,0]^2*dmat[i,*]*scale[tmp,*])-$
					(ti_inv[k]-ipro[k,0]*ti_inv[i]/ipro[i,0])/(ipro[i,0]^2)*dmat[i,*]*scale[tmp,*])^2*sigphot[tmp])
				ti_var_mu1[i]+=total((2.0*fmat[k,*]*imom[tmp,1]/ipro[i,0]-ipro[k,0]/ipro[i,0]*2.0*fmat[i,*]*imom[tmp,1]/ipro[i,0])^2*$
					width[tmp]^2/sigphot[tmp]*(1.0+bckphot[tmp]/sigphot[tmp]))
				ti_var_w1[i]+=total((2.0*fmat[k,*]*width[tmp]*imom[tmp,0]/ipro[i,0]-ipro[k,0]/ipro[i,0]*2.0*fmat[i,*]*width[tmp]*imom[tmp,0]/ipro[i,0])^2*$
					width[tmp]^2/(2.0*(sigphot[tmp]))*(1.0+bckphot[tmp]/sigphot[tmp]))
			ENDIF
		ENDFOR
		ti_err=sqrt(ti_var_n+ti_var_mu+ti_var_w)
		iproerr[0:nrho-1,3]=ti_err/ticonv
		IF keyword_set(im1svox) THEN BEGIN
			ti_err1=sqrt(ti_var_n1+ti_var_mu1+ti_var_w1)
			iproerr[nrho:*,3]=ti_err1/ticonv
		ENDIF
	ENDIF ELSE BEGIN
		ti_err=br2.inverr
		iproerr[0:nrho-1,3]=ipro[*,3]*sqrt(ti_err^2/ti_inv^2+iproerr[*,0]^2/ipro[*,0]^2)
	ENDELSE
END

	     
;+
;NAME:
;	HIREXSR_CALC_PROFILES
;
;MODIFICATION HISTORY:
;	9/4/11		M.L. Reinke - added the ability to use the iwidth,ishift terms in the moment pointer and store the subcheck
;-

PRO hirexsr_calc_profiles,mom,voxel,velvoxel,eps,eta,lam_o,z,rho,profiles,nofirst=nofirst,solidbody=solidbody,nosub=nosub,fix=fix,m1svoxel=m1svoxel,m1svelvoxel=m1svelvoxel,sine=sine
	x=size(mom)
	ntime=x[1]
	IF NOT keyword_set(fix) THEN BEGIN
		profiles=ptrarr(ntime,/allocate)
		fix=intarr(ntime)+1
	ENDIF
	
	FOR i=0,ntime-1 DO BEGIN
		irho=*rho[i]
		IF irho[0] NE -1 AND fix[i] THEN BEGIN
			xmom=*mom[i]
			imom=xmom[*,0:2]
			imomerr=xmom[*,3:5]
			ipmom=xmom[*,[10,11,12,17,18]]		;N, mu, w, const
			igood=xmom[*,6]
			iwidth=xmom[*,19]			;instrumental width [Ang]
			ishift=xmom[*,20]			;instrumental shift [Ang]
			ieps=eps[*,i]
			ieta=eta[*,i]
			ivox=*voxel[i]
			ivelvox=*velvoxel[i]
			IF keyword_set(m1svoxel) AND keyword_set(sine) THEN im1svox=*m1svoxel[i] ELSE im1svox=0
			IF keyword_set(m1svelvoxel) AND keyword_set(sine) THEN im1svelvox=*m1svelvoxel[i] ELSE im1svelvoxel=0
			IF total(ivox) NE 0 THEN BEGIN		;if empty voxel matrix (TMAP includes a disruption frame) then avoid doing the inversion
				hirexsr_invert_moments,imom,ivox,ivelvox,ieps,ieta,imomerr,igood,lam_o,z,ipro,iproerr,icheck,isubcheck,solidbody=solidbody,ipmom=ipmom,$
					nosub=nosub,im1svox=im1svox,im1svelvox=im1svelvox,iwidth=iwidth,ishift=ishift
			ENDIF ELSE BEGIN
				ipro=fltarr(n(irho)+1,4)
				iproerr=ipro
			ENDELSE
			IF keyword_set(sine) THEN BEGIN
				arr=[[ipro],$
				     [iproerr],$
				     [irho,irho]]
			ENDIF ELSE BEGIN
				arr=[[ipro],$
				     [iproerr],$
				     [irho]]
			ENDELSE
			*profiles[i]=arr
			xmom[*,7:9]=icheck
			xmom[*,21:23]=isubcheck
			*mom[i]=xmom
		ENDIF ELSE IF irho[0] EQ -1 THEN *profiles[i]=-1
	ENDFOR
END