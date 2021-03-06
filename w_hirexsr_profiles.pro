FUNCTION is_dat,u
	taglist=tag_names(u)
	tmp=where(taglist EQ 'DAT')
	IF tmp[0] NE -1 THEN output=1 ELSE output=0
	RETURN,output
END

PRO wprof_write_linelist,u
	line=[0,1,2,3,4,5,6,7,8,9]
	line_label=['w','x','z','lya1','4d','J','w','lya1','4d','lya1']
	line_z=[18,18,18,18,42,18,20,18,42,20]
	hirexsr_load_wavelengths,-1,lam_o,z,label
	FOR i=0,n(line) DO BEGIN
		ilam=lam_o[where(z EQ line_z[i] AND label EQ line_label[i])]
		IF line[i] EQ 7 or line[i] EQ 8 THEN extra=' (on Branch A)' ELSE extra=''
		str='LINE='+num2str(line[i],1)+' '+read_atomic_name(line_z[i])+' '+line_label[i]+' @ '+num2str(ilam,dp=5)+' [Ang]'+extra
		widget_control,u.id.linelist,set_value=str,/app
	ENDFOR
END

PRO wprofiles2ps,u
	psplot
	u.stat.ps=1
	plot_moments,u
	plot_profiles,u
	psc
	xwplot
	set_plot,'x'
	u.stat.ps=0
END

PRO update_current_voxel,u
	IF u.stat.m1s AND u.stat.m1scalc THEN sine=1 ELSE sine=0	;if m1s data already exists and you've selected to 
	hirexsr_pos2voxel,u.shot,u.dat.tau,u.dat.pos,u.dat.tpos,u.dat.rho,new_voxel,new_velvoxel,index=u.index,m1svoxel=new_m1sv,m1svelvoxel=new_m1svv,sine=sine
	i=u.index
	inew=*new_voxel[i]
	*u.dat.voxel[i]=inew
	inew=*new_velvoxel[i]
	*u.dat.velvoxel[i]=inew
	heap_free,new_voxel
	heap_free,new_velvoxel
	IF keyword_set(sine) THEN BEGIN
		inew=*new_m1sv[i]
		*u.dat.m1svoxel[i]=inew
		inew=*new_m1svv[i]
		*u.dat.m1svelvoxel[i]=inew	
		heap_free,new_m1sv
		heap_free,new_m1svv
	ENDIF
END

PRO update_all_voxel,u	
	IF u.stat.m1scalc THEN sine=1 ELSE sine=0	;select desire to calculate m1s terms
	hirexsr_pos2voxel,u.shot,u.dat.tau,u.dat.pos,u.dat.tpos,u.dat.rho,new_voxel,new_velvoxel,sine=sine,m1svoxel=new_m1sv,m1svelvoxel=new_m1svv
	FOR i=0,u.stat.ntime-1 DO BEGIN
		inew=*new_voxel[i]
		*u.dat.voxel[i]=inew
		inew=*new_velvoxel[i]
		*u.dat.velvoxel[i]=inew
		IF keyword_set(sine) THEN BEGIN
			inew=*new_m1sv[i]
			*u.dat.m1svoxel[i]=inew
			inew=*new_m1svv[i]
			*u.dat.m1svelvoxel[i]=inew
			u.stat.m1s=1
		ENDIF ELSE BEGIN		
			*u.dat.m1svoxel[i]=0
			*u.dat.m1svelvoxel[i]=0
			u.stat.m1s=0
		ENDELSE
	ENDFOR
	heap_free,new_voxel
	heap_free,new_velvoxel
	IF u.stat.m1s THEN BEGIN
		heap_free,new_m1sv
		heap_free,new_m1svv
	ENDIF
	
END


PRO invert_all_profiles,u
	fix=intarr(u.stat.ntime)+1
	hirexsr_calc_profiles,u.dat.mom,u.dat.voxel,u.dat.velvoxel,u.dat.eps,u.dat.eta,u.dat.lam_o,u.dat.z,u.dat.rho,u.dat.profiles,/solidbody,fix=fix,nosub=u.dat.nosub,$
		m1svoxel=u.dat.m1svoxel,m1svelvoxel=u.dat.m1svelvoxel,sine=u.inv.m1s
	widget_control,u.id.base, set_uvalue=u
END

PRO invert_current_profile,u
	fix=intarr(u.stat.ntime)
	fix[u.index]=1
	hirexsr_calc_profiles,u.dat.mom,u.dat.voxel,u.dat.velvoxel,u.dat.eps,u.dat.eta,u.dat.lam_o,u.dat.z,u.dat.rho,u.dat.profiles,/solidbody,fix=fix,nosub=u.dat.nosub,$
		m1svoxel=u.dat.m1svoxel,m1svelvoxel=u.dat.m1svelvoxel,sine=u.inv.m1s
	widget_control,u.id.base, set_uvalue=u
END

PRO save_profiles,u
	IF u.stat.m1s THEN BEGIN
		m1svoxel=u.dat.m1svoxel
		m1svelvoxel=u.dat.m1svelvoxel
	ENDIF ELSE BEGIN
		m1svoxel=0
		m1svelvoxel=0
	ENDELSE
	hirexsr_write_proconfig,u.shot,u.stat.line,u.dat.eps,u.dat.eta,tht=u.tht
	hirexsr_profile_ptr2arr,u.dat.profiles,pro_arr,proerr_arr,rho_arr
	hirexsr_write_profile,u.shot,u.stat.line,pro_arr,proerr_arr,rho_arr,u.dat.tau,tinst=u.dat.tinst,tgood=u.dat.tgood,tht=u.tht
	hirexsr_write_voxel,u.shot,u.stat.line,u.dat.voxel,u.dat.velvoxel,u.dat.rho,u.dat.tau,u.dat.tree,m1svoxel=m1svoxel,m1svelvoxel=m1svelvoxel,tht=u.tht
	hirexsr_extract_check,u.dat.mom,check,good,subcheck
	hirexsr_write_check,u.shot,u.stat.line,check,good,subcheck,tht=u.tht
	comment=logname()+','+systime()+','+u.dat.note
	hirexsr_write_analysis_comment,u.shot,comment,tht=u.tht
END


PRO plot_moments,u
	
	imom=*u.dat.mom[u.index]
	x=size(imom)
	imom[*,1]*=1.0e3
	imom[*,4]*=1.0e3	
	imom[*,8]*=1.0e3
	imom[*,2]*=1.0e6
	imom[*,5]*=1.0e6	
	imom[*,9]*=1.0e6
	imom[*,21]*=1.0e3
	imom[*,22]*=1.0e6
	imom[*,23]*=1.0e6

	nch=x[1]
	ch=indgen(nch)+1
	xr=[1,nch] 
	good=imom[*,6]
	tmp=where(good EQ 1)
	IF u.stat.ps THEN col=u.stat.pscol ELSE col=u.stat.col

	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw1,get_value=draw_win
		window,0,xsize=u.stat.plot.msize[0],ysize=u.stat.plot.msize[1],/pixmap
	ENDIF ELSE BEGIN
		xsize=4.0
		ysize=4.0*u.stat.plot.msize[1]/u.stat.plot.msize[0]
		ls=0.7
		!p.thick=6
		!p.charthick=2.0
		!x.thick=5.0
		!y.thick=5.0	
		device, xsize=xsize, ysize=ysize, /inches
	ENDELSE
	nplot=3.0

	;plot 0th moment and residual
	i=2
	pos=[0.125,0.3*1/nplot+i/nplot,0.97,0.975*1/nplot+i/nplot]
	IF u.stat.auto.m0 THEN BEGIN
		yr=[0,max(imom[tmp,0]+imom[tmp,3])]*1.05
		u.stat.plot.m0=yr
	ENDIF ELSE yr=u.stat.plot.m0
	plot,[0],[0],xr=xr,yr=yr,ytit='0!uth!n Moment',pos=pos,/xsty,/ysty
	makesym,9
	oploterror,ch,imom[*,0],imom[*,3],psym=8,color=col[0]
	makesym,10
	oploterror,ch[tmp],imom[tmp,0],imom[tmp,3],psym=8,color=col[0]
	oplot,ch[tmp],imom[tmp,7],color=col[1]
	IF NOT u.stat.ps THEN BEGIN
		oplot, ch[u.ch-1]*[1,1],yr,color=col[2],linestyle=2.0
		oplot, xr,imom[u.ch-1,0]*[1,1],color=col[2],linestyle=2.0
	ENDIF
	resid=imom[tmp,0]-imom[tmp,7]
	pos=[0.125,0.1*1/nplot+i/nplot,0.97,0.3*1/nplot+i/nplot]
	yr=(abs(max(resid+imom[tmp,3])) > abs(max(resid-imom[tmp,3])))*[-1.0,1.0]
	plot,[0],[0],xr=xr,yr=yr,ytit=n2g('Delta')+'!l0!n',xtit='CH #',pos=pos,/noerase,/xsty,/ysty
	oploterror,ch[tmp],resid,imom[tmp,3],psym=8,color=col[1],errcolor=col[1]
	oplot,xr,[0,0],linestyle=2.0,col=col[0]
	IF NOT u.stat.ps THEN BEGIN
		oplot, ch[u.ch-1]*[1,1],yr,color=col[2],linestyle=2.0
;		wset,draw_win
;		device,copy=[0,0,u.stat.plot.msize[0],u.stat.plot.msize[1],0,0,0]
	ENDIF		

	;plot 1st moment and residual
	i=1
	IF NOT u.stat.ps THEN BEGIN
;		widget_control,u.id.draw2,get_value=draw_win
;		window,0,xsize=u.stat.plot.msize[0],ysize=u.stat.plot.msize[1],/pixmap
	ENDIF
	pos=[0.125,0.3*1/nplot+i/nplot,0.97,0.975*1/nplot+i/nplot]
	IF u.stat.auto.m1 THEN BEGIN
		yr=[min(imom[tmp,1]-imom[tmp,4]),max(imom[tmp,1]+imom[tmp,4])]
		IF yr[0] LT 0 THEN yr[0]=yr[0]*1.05 ELSE yr[0]=yr[0]*0.95
		IF yr[1] LT 0 THEN yr[1]=yr[1]*0.95 ELSE yr[1]=yr[1]*1.05
		u.stat.plot.m1=yr
	ENDIF ELSE yr=u.stat.plot.m1
	plot,[0],[0],xr=xr,yr=yr,ytit='10!u-3!n 1!ust!n Moment',pos=pos,/xsty,/ysty,/noerase
	makesym,9
	oploterror,ch,imom[*,1],imom[*,4],psym=8,color=col[0]
	makesym,10
	oploterror,ch[tmp],imom[tmp,1],imom[tmp,4],psym=8,color=col[0]
	oplot,ch[tmp],imom[tmp,8],color=col[1]					;total check
	oplot,ch,imom[*,21],color=col[3],linestyle=4				;instrumental sub_vec
	oplot,xr,[0,0],linestyle=2.0,col=col[0]
	IF NOT u.stat.ps THEN BEGIN
		oplot, ch[u.ch-1]*[1,1],yr,color=col[2],linestyle=2.0
		oplot, xr,imom[u.ch-1,1]*[1,1],color=col[2],linestyle=2.0
	ENDIF
	resid=imom[tmp,1]-imom[tmp,8]
	pos=[0.125,0.1*1/nplot+i/nplot,0.97,0.3*1/nplot+i/nplot]
	yr=(abs(max(resid+imom[tmp,4])) > abs(max(resid-imom[tmp,4])))*[-1.0,1.0]
	plot,[0],[0],xr=xr,yr=yr,ytit=n2g('Delta')+'!l0!n',xtit='CH #',pos=pos,/noerase,/xsty,/ysty
	oploterror,ch[tmp],resid,imom[tmp,4],psym=8,color=col[1],errcolor=col[1]
	oplot,xr,[0,0],linestyle=2.0,col=col[0]
	IF NOT u.stat.ps THEN BEGIN
		oplot, ch[u.ch-1]*[1,1],yr,color=col[2],linestyle=2.0
;		wset,draw_win
;		device,copy=[0,0,u.stat.plot.msize[0],u.stat.plot.msize[1],0,0,0]
	ENDIF

	;plot 2nd moment and residual
	i=0
	IF NOT u.stat.ps THEN BEGIN
;		widget_control,u.id.draw3,get_value=draw_win
;		window,0,xsize=u.stat.plot.msize[0],ysize=u.stat.plot.msize[1],/pixmap
	ENDIF
	pos=[0.125,0.3*1/nplot+i/nplot,0.97,0.975*1/nplot+i/nplot]
	IF u.stat.auto.m2 THEN BEGIN
		ymin=min(imom[*,23]) < 0
		yr=[ymin,max(imom[tmp,2]+imom[tmp,5])]*1.05
		u.stat.plot.m2=yr
	ENDIF ELSE yr=u.stat.plot.m2
	plot,[0],[0],xr=xr,yr=yr,ytit='10!u-6!n 2!und!n Moment',pos=pos,/xsty,/ysty,/noerase
	makesym,9
	oploterror,ch,imom[*,2],imom[*,5],psym=8,color=col[0]
	makesym,10
	oploterror,ch[tmp],imom[tmp,2],imom[tmp,5],psym=8,color=col[0]
	oplot,ch[tmp],imom[tmp,9],color=col[1]					;total check
	oplot,ch,imom[*,22],color=col[3],linestyle=4				;instrumental sub_vec
	oplot,ch,imom[*,23],color=col[1],linestyle=3				;rotation sub_vec
	oplot,xr,[0,0],linestyle=2.0,col=col[0]
	IF NOT u.stat.ps THEN BEGIN
		oplot, ch[u.ch-1]*[1,1],yr,color=col[2],linestyle=2.0
		oplot, xr,imom[u.ch-1,2]*[1,1],color=col[2],linestyle=2.0
	ENDIF
	resid=imom[tmp,2]-imom[tmp,9]
	pos=[0.125,0.1*1/nplot+i/nplot,0.97,0.3*1/nplot+i/nplot]
	yr=(abs(max(resid+imom[tmp,5])) > abs(max(resid-imom[tmp,5])))*[-1.0,1.0]
	plot,[0],[0],xr=xr,yr=yr,ytit=n2g('Delta')+'!l2!n',xtit='CH #',pos=pos,/noerase,/xsty,/ysty
	oploterror,ch[tmp],resid,imom[tmp,5],psym=8,color=col[1],errcolor=col[1]
	oplot,xr,[0,0],linestyle=2.0,col=col[0]
	IF NOT u.stat.ps THEN BEGIN
		oplot, ch[u.ch-1]*[1,1],yr,color=col[2],linestyle=2.0
		wset,draw_win
		device,copy=[0,0,u.stat.plot.msize[0],u.stat.plot.msize[1],0,0,0]
	ENDIF
	
	update_mom_plotting_text,u
END

PRO plot_profiles,u
	c=2.998e8 				;speed of light
	e=1.602e-19				;conversion for eV -> J
	mconv=1.661e-27				;conversion for amu -> kg
	mass=read_atomic_mass(u.dat.z)		;mass of argon in amu
	nrho=n(*u.dat.rho[u.index])+1		;number of radial values
	lam_o=u.dat.lam_o
	
	ipro=*u.dat.profiles[u.index]
	imom=*u.dat.mom[u.index]
	IF n(ipro[*,0]) NE nrho-1 THEN m1splot=1 ELSE m1splot=0
	IF u.stat.auto.rho THEN BEGIN
		xr=[min(ipro[*,8]),max(ipro[*,8])] 
		u.stat.plot.rho=xr
	ENDIF ELSE xr=u.stat.plot.rho
	IF u.stat.ps THEN col=u.stat.pscol ELSE col=u.stat.col

	;plot emissivity
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw4,get_value=draw_win
		window,0,xsize=u.stat.plot.psize[0],ysize=u.stat.plot.psize[1],/pixmap
	ENDIF ELSE BEGIN
		xsize=4.5
		ysize=4.5*u.stat.plot.psize[1]/u.stat.plot.psize[0]
		ls=0.7
		!p.thick=6
		!p.charthick=2.0
		!x.thick=5.0
		!y.thick=5.0	
		device, xsize=xsize, ysize=ysize, /inches
	ENDELSE
	nplot=4.0

	;plot emissivity
	i=3
	pos=[0.125,0.12*1/nplot+i/nplot,0.97,0.975*1/nplot+i/nplot]
	IF u.stat.auto.em THEN BEGIN
		yr=[0 < min(ipro[*,0]-ipro[*,4]),max(ipro[*,0]+ipro[*,4])]*1.05
		u.stat.plot.em=yr
	ENDIF ELSE yr=u.stat.plot.em
	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit=n2g('psi')+'!ln!n',ytit='EMISS',pos=pos
	oploterror,ipro[0:nrho-1,8],ipro[0:nrho-1,0],ipro[0:nrho-1,4],color=col[0]
	IF m1splot THEN oploterror,ipro[nrho:*,8],ipro[nrho:*,0],ipro[nrho:*,4],color=col[0],linestyle=3.0
	oplot,xr,[0,0],linestyle=2.0,color=col[0]
	IF NOT u.stat.ps THEN BEGIN
;		wset,draw_win
;		device,copy=[0,0,u.stat.plot.psize[0],u.stat.plot.psize[1],0,0,0]
	ENDIF

	;plot toroidal rotation freqency
	i=2
	pos=[0.125,0.12*1/nplot+i/nplot,0.97,0.975*1/nplot+i/nplot]
	winst_dc=0.0
	winst_ch=hirexsr_is2iv(imom[*,20],lam_o)
	widget_control,u.id.ishift_text,set_value=num2str(winst_ch[u.ch-1]*1.0e3,dp=1)
	good=imom[*,6]
	IF NOT u.stat.ps THEN BEGIN
;		widget_control,u.id.draw5,get_value=draw_win
;		window,0,xsize=u.stat.plot.psize[0],ysize=u.stat.plot.psize[1],/pixmap
	ENDIF
	IF u.stat.auto.vt THEN BEGIN
		yr=[min(ipro[*,1]-ipro[*,5]),max(ipro[*,1]+ipro[*,5])]
		IF yr[0] LT 0 THEN yr[0]=yr[0]*1.05 ELSE yr[0]=yr[0]*0.95
		IF yr[1] LT 0 THEN yr[1]=yr[1]*0.95 ELSE yr[1]=yr[1]*1.05
		u.stat.plot.vt=yr
	ENDIF ELSE yr=u.stat.plot.vt
	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit=n2g('psi')+'!ln!n',ytit=n2g('omega')+' [kHz]',pos=pos,/noerase
	oploterror,ipro[0:nrho-1,8],ipro[0:nrho-1,1],ipro[0:nrho-1,5],color=col[0]
	IF m1splot THEN oploterror,ipro[nrho:*,8],ipro[nrho:*,1],ipro[nrho:*,5],color=col[0],linestyle=3.0
	oplot,xr,[0,0],linestyle=2.0,color=col[0]
	IF u.stat.plot.lifit THEN BEGIN
		imlint=*u.dat.mlint[u.index]
		tmp=where(good EQ 1)
		makesym,9
		oploterror,imlint[tmp,4],imlint[tmp,0],imlint[tmp,1],psym=8,color=col[4],errcolor=col[4]
		makesym,10
	ENDIF
	IF u.stat.plot.limom THEN BEGIN
		tmp=where(good EQ 1)
		pindex=last(where(u.dat.tpos LE u.stat.time))
		jpos=u.dat.pos[*,*,pindex]
		v=-1.0*(lam_o-imom[tmp,11])*c/lam_o/(2.0*!pi*jpos[2,tmp]*cos(jpos[3,tmp]))*1.0e-3
		verr=imom[tmp,14]*c/lam_o/(2.0*!pi*jpos[2,tmp]*cos(jpos[3,tmp]))*1.0e-3
		oploterror,imom[tmp,16],v,verr,psym=5,color=col[5],errcolor=col[5]
		IF NOT u.stat.ps THEN BEGIN
			oplot, imom[u.ch-1,16]*[1,1],yr,color=col[2],linestyle=2.0
			IF good[u.ch-1] EQ 1 THEN BEGIN
				ypt=v[where(tmp EQ u.ch-1)]
				oplot, xr,ypt[0]*[1,1],color=col[2],linestyle=2.0
			ENDIF
		ENDIF
	ENDIF
	IF NOT u.stat.ps THEN BEGIN
;		wset,draw_win
;		device,copy=[0,0,u.stat.plot.psize[0],u.stat.plot.psize[1],0,0,0]
	ENDIF

	;plot poloidal 
	i=1
	pos=[0.125,0.12*1/nplot+i/nplot,0.97,0.975*1/nplot+i/nplot]
	IF NOT u.stat.ps THEN BEGIN
;		widget_control,u.id.draw6,get_value=draw_win
;		window,0,xsize=u.stat.plot.psize[0],ysize=u.stat.plot.psize[1],/pixmap
	ENDIF
	IF u.stat.auto.vp THEN BEGIN
		yr=[min(ipro[*,2]-ipro[*,6]),max(ipro[*,2]+ipro[*,6])]
		IF yr[0] LT 0 THEN yr[0]=yr[0]*1.05 ELSE yr[0]=yr[0]*0.95
		IF yr[1] LT 0 THEN yr[1]=yr[1]*0.95 ELSE yr[1]=yr[1]*1.05
		IF yr[0] EQ yr[1] THEN yr=[-1,1]
		u.stat.plot.vp=yr
	ENDIF ELSE yr=u.stat.plot.vp
	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit=n2g('psi')+'!ln!n',ytit='vpol [km/s]',pos=pos,/noerase
	oploterror,ipro[0:nrho-1,8],ipro[0:nrho-1,2],ipro[0:nrho-1,6],color=col[0]
	IF m1splot THEN oploterror,ipro[nrho:*,8],ipro[nrho:*,2],ipro[nrho:*,6],color=col[0]
	oplot,xr,[0,0],linestyle=2.0,color=col[0]
	good=imom[*,6]
	IF u.stat.plot.lifit THEN BEGIN
		imlint=*u.dat.mlint[u.index]
		tmp=where(imlint[*,6] NE 0 AND good EQ 1)
		makesym,9
		if tmp[0] ne -1 then oploterror,imlint[tmp,4],imlint[tmp,6],imlint[tmp,7],psym=8,color=col[4],errcolor=col[4]
		makesym,10
	ENDIF
	IF NOT u.stat.ps THEN BEGIN
;		wset,draw_win
;		device,copy=[0,0,u.stat.plot.psize[0],u.stat.plot.psize[1],0,0,0]
        ENDIF


	;plot ion temperature 
	i=0
	pos=[0.125,0.12*1/nplot+i/nplot,0.97,0.975*1/nplot+i/nplot]
	tinst_dc=u.dat.tinst
	tinst_ch=hirexsr_iw2iti(imom[*,19],u.dat.z,lam_o)
	widget_control,u.id.iwidth_text,set_value=num2str(tinst_ch[u.ch-1]*1.0e3,dp=1)
	IF NOT u.stat.ps THEN BEGIN
;		widget_control,u.id.draw7,get_value=draw_win
;		window,0,xsize=u.stat.plot.psize[0],ysize=u.stat.plot.psize[1],/pixmap
	ENDIF
	IF u.stat.auto.ti THEN BEGIN
		yr=[min(ipro[*,3]-ipro[*,7]),max(ipro[*,3]+ipro[*,7])]
		u.stat.plot.ti=yr
	ENDIF ELSE yr=u.stat.plot.ti
	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit=n2g('psi')+'!ln!n',ytit='T!lI!n [keV]',pos=pos,/noerase
	oploterror,ipro[0:nrho-1,8],ipro[0:nrho-1,3]-tinst_dc,ipro[0:nrho-1,7],color=col[0]
	IF m1splot THEN oploterror,ipro[nrho:*,8],ipro[nrho:*,3],ipro[nrho:*,7],color=col[0],linestyle=3.0
	oplot,xr,[0,0],linestyle=2.0,color=col[0]
	IF u.stat.plot.ts THEN BEGIN
		its=ipt(u.dat.ts.time,u.stat.time)
		IF its[0] NE -1 THEN oploterror,u.dat.ts.rho[*,its],u.dat.ts.temp[*,its],u.dat.ts.err[*,its],color=col[3],errcolor=col[3],psym=u.dat.ts.usym
	ENDIF

	good=imom[*,6]
	IF u.stat.plot.lifit THEN BEGIN
		imlint=*u.dat.mlint[u.index]
		tmp=where(imlint[*,1] NE 0 AND good EQ 1)
		makesym,9
		oploterror,imlint[tmp,4],imlint[tmp,2]-tinst_ch[tmp]-tinst_dc,imlint[tmp,3],psym=8,color=col[4]
		makesym,10
	ENDIF
	IF u.stat.plot.limom THEN BEGIN
		tmp=where(good EQ 1)
		jmom=(imom[*,12])^2			;w^2
		jerr=imom[*,15]*2.0*sqrt(jmom)		;uncertainty prop through the square
		jmom*=mass*mconv*c^2/(e*1.0e3*u.dat.lam_o^2)
		jerr*=mass*mconv*c^2/(e*1.0e3*u.dat.lam_o^2)
		oploterror,imom[tmp,16],jmom[tmp]-tinst_ch[tmp]-tinst_dc,jerr[tmp],psym=5,color=col[5],errcolor=col[5]
		IF NOT u.stat.ps THEN BEGIN
			oplot, imom[u.ch-1,16]*[1,1],yr,color=col[2],linestyle=2.0
			IF good[u.ch-1] EQ 1 THEN BEGIN
				ypt=jmom[u.ch-1]-tinst_ch[u.ch-1]-tinst_dc
				oplot, xr,ypt[0]*[1,1],color=col[2],linestyle=2.0
			ENDIF
		ENDIF
	ENDIF
	update_pro_plotting_text,u

	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.stat.plot.psize[0],u.stat.plot.psize[1],0,0,0]
	ENDIF
		
END

PRO load_profile_data,u
	shot=u.shot
	tht=u.tht
	line=u.stat.line
	IF line GT 2 THEN u.ch=2
	IF u.stat.dat THEN BEGIN					;if LOAD_DATA already ran, clean out the heap
		cleanheap_w_hirexsr_profiles,u
		heap_gc
	ENDIF

	IF u.stat.plot.ts THEN BEGIN
		IF NOT u.stat.dat THEN BEGIN		;load TS profile only once
			cmod_ts,u.shot,ts
			rho=efit_rz2rho(ts.r_ts,ts.z_ts,ts.time,shot=u.shot,/psinorm)
			ts={temp:rotate(ts.te,4),err:rotate(ts.dte,4),time:ts.time,rho:rho,usym:7,shot:shot}
		ENDIF ELSE BEGIN
			IF shot NE u.dat.ts.shot THEN BEGIN
				cmod_ts,u.shot,ts
				rho=efit_rz2rho(ts.r_ts,ts.z_ts,ts.time,shot=u.shot,/psinorm)
				ts={temp:rotate(ts.te,4),err:rotate(ts.dte,4),time:ts.time,rho:rho,usym:7,shot:shot}
			ENDIF ELSE ts=u.dat.ts
		ENDELSE
	ENDIF ELSE ts=-1

	;load moments and profile data from tree
	hirexsr_load_momentptr,shot,line,mom,tau,pos,tpos,lam_o,z,status=status,tht=tht,good=-1
	IF NOT status THEN BEGIN
		widget_control,u.id.message,set_value='No LINE='+num2str(line)+' available for'+num2str(shot,1)+' THT-'+num2str(tht,1),/app
		RETURN
	ENDIF
	hirexsr_load_mlintptr,shot,line,mlint,tau,tht=tht
	hirexsr_load_voxel,shot,line,voxel,velvoxel,rho,tau,tree,/ptr,m1svoxel=m1svoxel,m1svelvoxel=m1svelvoxel,tht=tht
	hirexsr_load_proconfig,shot,line,eps,eta,tht=tht
	hirexsr_load_profile,shot,line,pro_arr,proerr_arr,rho_arr,tau,tinst=tinst,tgood=tgood,tht=tht
	hirexsr_load_analysis_comment,shot,comment,tht=tht
	profiles=hirexsr_profile_arr2ptr(pro_arr,proerr_arr,rho_arr)
	widget_control,u.id.tinst_text,set_value=num2str(tinst*1.0e3,dp=1)
	widget_control,u.id.tinst_slider,set_value=tinst*1.0e3
	widget_control,u.id.message,set_value='LINE='+num2str(line)+' MOM/PRO Loaded: '+num2str(shot,1)+' THT-'+num2str(tht,1),/app
	IF comment NE 'error' THEN BEGIN
		split=strsplit(comment,',',/extract)
		widget_control,u.id.message,set_value=split[0]+','+split[1],/app
		note=split[2]
	ENDIF ELSE note=comment
	widget_control,u.id.note,set_value=note

	IF u.stat.time GT max(tau) THEN u.stat.time=max(tau)
	u.index=ipt(tau,u.stat.time)	;initialize index
	u.stat.ntime=n(tau)+1
	imom=*mom[u.index]
	u.stat.nch=n(imom[*,0])+1
	u.stat.good=[0,n(imom[*,0])]					;set no large good filtering
	widget_control,u.id.glow_text,set_value=num2str(u.stat.good[0]+1,1)
	widget_control,u.id.glow_slider,set_value=u.stat.good[0]+1
	widget_control,u.id.ghigh_text,set_value=num2str(u.stat.good[1]+1,1)
	widget_control,u.id.ghigh_slider,set_value=u.stat.good[1]+1
	IF keyword_set(m1svoxel) AND keyword_set(m1svelvoxel) THEN BEGIN
		u.stat.m1s=1
		u.stat.m1scalc=1
	ENDIF ELSE BEGIN
		m1svoxel=ptrarr(u.stat.ntime,/allocate_heap)
		FOR i=0,u.stat.ntime-1 DO *m1svoxel[i]=0
		m1svelvoxel=ptrarr(u.stat.ntime,/allocate_heap)
		FOR i=0,u.stat.ntime-1 DO *m1svelvoxel[i]=0
		u.stat.m1s=0
		u.stat.m1scalc=0
	ENDELSE
	widget_control,u.id.m1scalc,set_button=u.stat.m1scalc

	;make a copy of the good
	good=ptrarr(u.stat.ntime,/allocate)
	FOR i=0,u.stat.ntime-1 DO BEGIN
		imom=*mom[i]
		IF n(imom) EQ 0 THEN *good[i]=-1 ELSE *good[i]=imom[*,6]
	ENDFOR

	dat={mom:mom,mlint:mlint,pos:pos,tpos:tpos,voxel:voxel,velvoxel:velvoxel,m1svoxel:m1svoxel,m1svelvoxel:m1svelvoxel,rho:rho,good:good,eps:eps,eta:eta,nosub:0,tinst:tinst,tgood:tgood,$
		profiles:profiles,tau:tau,tree:tree,lam_o:lam_o,z:z,ts:ts,note:note}
	u={id:u.id,shot:u.shot,tht:u.tht,ch:u.ch,index:u.index,inv:u.inv,stat:u.stat,dat:dat}
	u.stat.dat=1
	IF u.dat.tgood[u.index] THEN widget_control,u.id.tgood,set_button=1 ELSE widget_control,u.id.tgood,set_button=0
	igood=*u.dat.good[u.index]
	IF igood[u.ch-1] THEN widget_control,u.id.chgood,set_button=1 ELSE widget_control,u.id.chgood,set_button=0
	update_ch_text,u
	widget_control,u.id.base, set_uvalue=u
END

PRO cleanheap_w_hirexsr_profiles,u
	taglist=tag_names(u)
	tmp=where(taglist EQ 'DAT')
	IF tmp[0] NE -1 THEN BEGIN
		heap_free,u.dat.mom
		heap_free,u.dat.mlint
		heap_free,u.dat.voxel
		heap_free,u.dat.velvoxel
		heap_free,u.dat.profiles
		heap_free,u.dat.rho
		heap_free,u.dat.good
		IF size(u.dat.m1svoxel,/type) EQ 10 THEN heap_free,u.dat.m1svoxel
		IF size(u.dat.m1svelvoxel,/type) EQ 10 THEN heap_free,u.dat.m1svelvoxel
	ENDIF
END

PRO reset_invline_button,u
	id=[u.id.inv_w,u.id.inv_x,u.id.inv_z,u.id.inv_lya,u.id.inv_m32,u.id.inv_caw]
	widget_control,id[u.stat.line],set_button=0
END

PRO update_rho_text,u,str=str
	index=u.index
	irho=*u.dat.rho[index]
	min=min(irho)
	max=max(irho)
	npts=n(irho)+1
	IF NOT keyword_set(str) THEN BEGIN
		widget_control,u.id.rmin_text,set_value=num2str(min,dp=2)
		widget_control,u.id.rmin_slider,set_value=int(min*100)
		widget_control,u.id.rmax_text,set_value=num2str(max,dp=2)
		widget_control,u.id.rmax_slider,set_value=int(max*100)
		widget_control,u.id.npts_text,set_value=num2str(npts,1)
		widget_control,u.id.npts_slider,set_value=npts
	ENDIF
	rhostr=num2str(irho[0],dp=2)
	FOR i=1,npts-1 DO rhostr=rhostr+','+num2str(irho[i],dp=2)
	widget_control,u.id.rho_text,set_value=rhostr
END

PRO update_temporal_text,u
	index=u.index
	imom=*u.dat.mom[index]
	seta_id=[u.id.em_eta_slider,u.id.vt_eta_slider,u.id.ti_eta_slider]
	seps_id=[u.id.em_eps_slider,u.id.vt_eps_slider,u.id.ti_eps_slider]
	teta_id=[u.id.em_eta_text,u.id.vt_eta_text,u.id.ti_eta_text]
	teps_id=[u.id.em_eps_text,u.id.vt_eps_text,u.id.ti_eps_text]
	inv=[u.inv.em,u.inv.vt,u.inv.ti]
	FOR i=0,2 DO BEGIN
		IF inv[i] THEN BEGIN
			eps=u.dat.eps[i,index]
			widget_control,seps_id[i],set_value=int(eps*100.0)		;eps=slider/100.0 (0.01,10)
			widget_control,teps_id[i],set_value=num2str(eps,dp=2)
			eta=u.dat.eta[i,index]
			widget_control,seta_id[i],set_value=int(alog10(eta)*10)		;eta=10^(slider/10) (1,10^5)
			widget_control,teta_id[i],set_value=num2str(alog10(eta),dp=2)
		ENDIF
	ENDFOR
	widget_control,u.id.itpt,set_value=num2str(index,1)
	u.stat.nch=n(imom[*,0])+1
	update_rho_text,u
	
END

PRO update_ch_text,u
	widget_control,u.id.ch_text,set_value=num2str(u.ch,1)
	widget_control,u.id.ch, set_value=num2str(u.ch,1)
	widget_control,u.id.chtot, set_value=num2str(u.stat.nch,1)	
END

PRO update_mom_plotting_text,u
	widget_control,u.id.m0min,set_value=num2str(u.stat.plot.m0[0],dp=2)
	widget_control,u.id.m0max,set_value=num2str(u.stat.plot.m0[1],dp=2)
	widget_control,u.id.m1min,set_value=num2str(u.stat.plot.m1[0],dp=2)
	widget_control,u.id.m1max,set_value=num2str(u.stat.plot.m1[1],dp=2)
	widget_control,u.id.m2min,set_value=num2str(u.stat.plot.m2[0],dp=2)
	widget_control,u.id.m2max,set_value=num2str(u.stat.plot.m2[1],dp=2)
END

PRO update_pro_plotting_text,u
	widget_control,u.id.xmin,set_value=num2str(u.stat.plot.rho[0],dp=2)
	widget_control,u.id.xmax,set_value=num2str(u.stat.plot.rho[1],dp=2)
	widget_control,u.id.emmin,set_value=num2str(u.stat.plot.em[0],dp=2)
	widget_control,u.id.emmax,set_value=num2str(u.stat.plot.em[1],dp=2)
	widget_control,u.id.vtmin,set_value=num2str(u.stat.plot.vt[0],dp=2)
	widget_control,u.id.vtmax,set_value=num2str(u.stat.plot.vt[1],dp=2)
	widget_control,u.id.vpmin,set_value=num2str(u.stat.plot.vp[0],dp=2)
	widget_control,u.id.vpmax,set_value=num2str(u.stat.plot.vp[1],dp=2)
	widget_control,u.id.timin,set_value=num2str(u.stat.plot.ti[0],dp=2)
	widget_control,u.id.timax,set_value=num2str(u.stat.plot.ti[1],dp=2)
END

PRO w_hirexsr_profiles_event,event
	widget_control,event.top,get_uvalue=u
	id = u.id
	tag = tag_names(event,/st)
	button=' '
	idtags=tag_names(id)
	FOR i=0,n(idtags) DO IF id.(i) EQ event.id THEN ename=idtags[i]
	CASE tag OF
		"WIDGET_BASE" : BEGIN

		END
		"WIDGET_BUTTON": BEGIN
			widget_control,event.id,get_value=button,get_uvalue=uvalue
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF
	 			"QUIT": BEGIN
					widget_control,event.top,/destroy
					cleanheap_w_hirexsr_profiles,u
					heap_gc
					!except=1
				
				END
				"SAVE" : BEGIN
					save_profiles,u
					widget_control,u.id.message,set_value='LINE='+num2str(u.stat.line)+' PROFILES and VOXELS saved to THT-'+num2str(u.tht,1),/app
				END
				"LOAD" : BEGIN
					load_profile_data,u
					IF u.stat.dat THEN BEGIN
						update_temporal_text,u
						plot_profiles,u
						plot_moments,u
						tmp=where(u.dat.tau GT 0)
						widget_control,u.id.ntpt,set_value=num2str(n(tmp),1)
					ENDIF
				END
				"PRINT" : BEGIN
					IF u.stat.dat THEN wprofiles2ps,u
				END
				"STOP" : BEGIN
					stop
				END
				"INV_W" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						reset_invline_button,u
						u.stat.line=0
						load_profile_data,u
						IF u.stat.dat THEN BEGIN
							update_temporal_text,u
							plot_profiles,u
							plot_moments,u
						ENDIF
					ENDIF ELSE widget_control,u.id.inv_w,set_button=1
				END
 				"INV_X" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						reset_invline_button,u
						u.stat.line=1
						load_profile_data,u
						IF u.stat.dat THEN BEGIN
							update_temporal_text,u
							plot_profiles,u
							plot_moments,u
						ENDIF
					ENDIF ELSE widget_control,u.id.inv_x,set_button=1
				END
 				"INV_Z" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						reset_invline_button,u
						u.stat.line=2
						load_profile_data,u
						IF u.stat.dat THEN BEGIN
							update_temporal_text,u
							plot_profiles,u
							plot_moments,u
						ENDIF
					ENDIF ELSE widget_control,u.id.inv_z,set_button=1
				END
 				"INV_LYA" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						reset_invline_button,u
						u.stat.line=3
						load_profile_data,u
						IF u.stat.dat THEN BEGIN
							update_temporal_text,u
							plot_profiles,u
							plot_moments,u
						ENDIF
					ENDIF ELSE widget_control,u.id.inv_lya,set_button=1
				END
				"INV_M32" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						reset_invline_button,u
						u.stat.line=4
						load_profile_data,u
						IF u.stat.dat THEN BEGIN
							update_temporal_text,u
							plot_profiles,u
							plot_moments,u
						ENDIF
					ENDIF ELSE widget_control,u.id.inv_m32,set_button=1
				END
				"INV_CAW" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						reset_invline_button,u
						u.stat.line=6
						load_profile_data,u
						IF u.stat.dat THEN BEGIN
							update_temporal_text,u
							plot_profiles,u
							plot_moments,u
						ENDIF
					ENDIF ELSE widget_control,u.id.inv_caw,set_button=1
				END
				'MT' : BEGIN
					IF u.index NE 0 AND u.stat.dat THEN BEGIN
						IF u.dat.tau[u.index-1] NE -1 THEN BEGIN
							u.index=u.index-1
							u.stat.time=u.dat.tau[u.index]
							widget_control,u.id.t_slider,set_value=u.stat.time*1.0e3
							update_temporal_text,u
							update_ch_text,u
							widget_control,id.time,set_value=num2str(u.stat.time,dp=2)
							plot_profiles,u
							plot_moments,u
							IF u.dat.tgood[u.index] THEN widget_control,u.id.tgood,set_button=1 ELSE widget_control,u.id.tgood,set_button=0
							igood=*u.dat.good[u.index]
							IF igood[u.ch-1] THEN widget_control,u.id.chgood,set_button=1 ELSE widget_control,u.id.chgood,set_button=0
						ENDIF
					ENDIF
				END
				'PT' : BEGIN
					IF u.index NE u.stat.ntime-1 AND u.stat.dat THEN BEGIN
						IF u.dat.tau[u.index+1] NE -1 THEN BEGIN
							u.index=u.index+1
							u.stat.time=u.dat.tau[u.index]
							widget_control,u.id.t_slider,set_value=u.stat.time*1.0e3
							update_temporal_text,u
							update_ch_text,u
							widget_control,id.time,set_value=num2str(u.stat.time,dp=2)
							plot_profiles,u
							plot_moments,u
							IF u.dat.tgood[u.index] THEN widget_control,u.id.tgood,set_button=1 ELSE widget_control,u.id.tgood,set_button=0
							igood=*u.dat.good[u.index]
							IF igood[u.ch-1] THEN widget_control,u.id.chgood,set_button=1 ELSE widget_control,u.id.chgood,set_button=0
						ENDIF
					ENDIF
				END
				'TGOOD' : BEGIN
					IF u.stat.dat THEN BEGIN
						IF event.select EQ 1 THEN u.dat.tgood[u.index]=1 ELSE u.dat.tgood[u.index]=0
					ENDIF
				END
				'MCH' : BEGIN
					IF u.ch NE 1 THEN BEGIN
						u.ch=u.ch-1
						widget_control,u.id.ch_slider,set_value=u.ch
						update_ch_text,u
						plot_moments,u
						plot_profiles,u
						igood=*u.dat.good[u.index]
						IF igood[u.ch-1] THEN widget_control,u.id.chgood,set_button=1 ELSE widget_control,u.id.chgood,set_button=0
					ENDIF
				END
				'PCH' : BEGIN
					IF u.ch NE u.stat.nch THEN BEGIN
						u.ch=u.ch+1
						widget_control,u.id.ch_slider,set_value=u.ch
						update_ch_text,u
						plot_moments,u
						plot_profiles,u
						igood=*u.dat.good[u.index]
						IF igood[u.ch-1] THEN widget_control,u.id.chgood,set_button=1 ELSE widget_control,u.id.chgood,set_button=0
					ENDIF
				END
				'CHGOOD' : BEGIN
					IF u.stat.dat THEN BEGIN
						imom=*u.dat.mom[u.index]
						igood=*u.dat.good[u.index]
						IF imom[0] NE -1 THEN BEGIN
							IF event.select EQ 1 THEN imom[u.ch-1,6]=1 ELSE imom[u.ch-1,6]=0
							*u.dat.mom[u.index]=imom
							IF event.select EQ 1 THEN igood[u.ch-1]=1 ELSE igood[u.ch-1]=0
							*u.dat.good[u.index]=igood
						ENDIF
					ENDIF
				END

				'REGOOD' : BEGIN
					IF u.stat.dat THEN BEGIN
						FOR i=0,u.stat.ntime-1 DO BEGIN
							imom=*u.dat.mom[i]
							igood=*u.dat.good[i]
							IF imom[0] NE -1 THEN BEGIN
								imom[u.ch-1,6]=0
								*u.dat.mom[i]=imom
								igood[u.ch-1]=0
								*u.dat.good[i]=igood
							ENDIF	
						ENDFOR	
						widget_control,u.id.message,set_value='CH '+num2str(u.ch,1)+' REMOVED FOR ALL TIMES',/app
					ENDIF
				END
				'ADDGOOD' : BEGIN
					IF u.stat.dat THEN BEGIN
						FOR i=0,u.stat.ntime-1 DO BEGIN
							imom=*u.dat.mom[i]
							igood=*u.dat.good[i]
							IF imom[0] NE -1 THEN BEGIN
								imom[u.ch-1,6]=1
								*u.dat.mom[i]=imom
								igood[u.ch-1]=1
								*u.dat.good[i]=igood
							ENDIF	
						ENDFOR	
						widget_control,u.id.message,set_value='CH '+num2str(u.ch,1)+' ADDED FOR ALL TIMES',/app
					ENDIF
				END
				'COPYGOOD' : BEGIN
					IF u.stat.dat THEN BEGIN
						imom=*u.dat.mom[u.index]
						good=imom[*,6]
						FOR i=0,u.stat.ntime-1 DO BEGIN
							imom=*u.dat.mom[i]
							igood=*u.dat.good[i]
							IF imom[0] NE -1 THEN BEGIN
								imom[*,6]=good
								*u.dat.mom[i]=imom
								igood[u.ch-1]=0
								*u.dat.good[i]=igood
							ENDIF	
						ENDFOR	
						widget_control,u.id.message,set_value='CH '+num2str(u.ch,1)+' FRAME '+num2str(u.index,1)+' GOOD SET FOR ALL TIME',/app
					ENDIF				
				END
				'GSLIDE' : BEGIN
					IF u.stat.dat THEN BEGIN
						IF event.select EQ 1 THEN BEGIN
							FOR i=0,u.stat.ntime-1 DO BEGIN
								imom=*u.dat.mom[i]
								IF imom[0] NE -1 THEN BEGIN
									igood=intarr(n(imom[*,0])+1)
									igood[u.stat.good[0]:u.stat.good[1]]=1
									imom[*,6]=igood
									*u.dat.mom[i]=imom
								ENDIF	
							ENDFOR
							widget_control,u.id.message,set_value='GOOD SET BY SLIDERS',/app
							u.stat.gslide=1
						ENDIF ELSE BEGIN
							FOR i=0,u.stat.ntime-1 DO BEGIN
								imom=*u.dat.mom[i]
								IF imom[0] NE -1 THEN BEGIN
									igood=*u.dat.good[i]
									imom[*,6]=igood
									*u.dat.mom[i]=imom
								ENDIF
							ENDFOR
							widget_control,u.id.message,set_value='GOOD ARRAY RESTORED',/app
							u.stat.gslide=0
						ENDELSE
					ENDIF
				END
				'M1SINV' : BEGIN
					IF u.stat.m1s THEN BEGIN
						IF event.select EQ 1 THEN u.inv.m1s=1 ELSE u.inv.m1s=0
					ENDIF ELSE BEGIN
						widget_control,u.id.message,set_value='NO m=1 VOXELS STORED',/app
						widget_control,u.id.m1sinv,set_button=0
						u.inv.m1s=0
					ENDELSE
				END
				'M1SCALC' : BEGIN
					IF event.select EQ 1 THEN u.stat.m1scalc=1 ELSE u.stat.m1scalc=0
				END
				'COPYRHO' : BEGIN
					IF u.stat.dat THEN BEGIN
						rhocopy=*u.dat.rho[u.index]
						FOR i=0,u.stat.ntime-1 DO BEGIN
							IF u.dat.tau[i] NE -1 THEN *u.dat.rho[i]=rhocopy
						ENDFOR	
						widget_control,u.id.message,set_value='FRAME '+num2str(u.index,1)+' RHO SET FOR ALL TIME',/app
					ENDIF				
				END

				"VOXCURR" : BEGIN
					WIDGET_CONTROL,/hourglass
					update_current_voxel,u
				 	invert_current_profile,u
					plot_profiles,u
					plot_moments,u
					widget_control,u.id.message,set_value='CURRENT VOXELS CALCULATED',/app
					widget_control,u.id.message,set_value='CURRENT PROFILE INVERTED',/app
				END
				"VOXALL" : BEGIN
					WIDGET_CONTROL,/hourglass
					update_all_voxel,u
					invert_all_profiles,u
					plot_profiles,u
					plot_moments,u
					widget_control,u.id.message,set_value='ALL VOXELS CALCULATED',/app
					widget_control,u.id.message,set_value='ALL PROFILES INVERTED',/app
				END
				"INVCURR" : BEGIN
				 	invert_current_profile,u
					plot_profiles,u
					plot_moments,u
					widget_control,u.id.message,set_value='CURRENT PROFILE INVERTED',/app
				END
				"INVALL" : BEGIN
					WIDGET_CONTROL,/hourglass
					invert_all_profiles,u
					plot_profiles,u
					plot_moments,u
					widget_control,u.id.message,set_value='ALL PROFILES INVERTED',/app
				END
				"TSPLOT" : BEGIN
					IF event.select EQ 1 THEN u.stat.plot.ts=1 ELSE u.stat.plot.ts=0
					IF u.stat.dat THEN plot_profiles,u
				END
				"ECEPLOT" : BEGIN
					IF event.select EQ 1 THEN u.stat.plot.ece=1 ELSE u.stat.plot.ece=0
					IF u.stat.dat THEN plot_profiles,u
				END
				"LIMOM" : BEGIN
					IF event.select EQ 1 THEN u.stat.plot.limom=1 ELSE u.stat.plot.limom=0
					IF u.stat.dat THEN plot_profiles,u
				END
				"LIFIT" : BEGIN
					IF event.select EQ 1 THEN u.stat.plot.lifit=1 ELSE u.stat.plot.lifit=0
					IF u.stat.dat THEN plot_profiles,u
				END
				"XAUTO" : IF event.select EQ 1 THEN u.stat.auto.rho=1 ELSE u.stat.auto.rho=0
				"EMAUTO" : IF event.select EQ 1 THEN u.stat.auto.em=1 ELSE u.stat.auto.em=0
				"VPAUTO" : IF event.select EQ 1 THEN u.stat.auto.vp=1 ELSE u.stat.auto.vp=0
				"VTAUTO" : IF event.select EQ 1 THEN u.stat.auto.vt=1 ELSE u.stat.auto.vt=0
				"TIAUTO" : IF event.select EQ 1 THEN u.stat.auto.ti=1 ELSE u.stat.auto.ti=0
				"M0AUTO" : IF event.select EQ 1 THEN u.stat.auto.m0=1 ELSE u.stat.auto.m0=0
				"M1AUTO" : IF event.select EQ 1 THEN u.stat.auto.m1=1 ELSE u.stat.auto.m1=0
				"M2AUTO" : IF event.select EQ 1 THEN u.stat.auto.m2=1 ELSE u.stat.auto.m2=0		
				ELSE:
			ENDCASE
		END
		"WIDGET_SLIDER": BEGIN
			widget_control,event.id,get_value=slider
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF 
				'CH_SLIDER' : BEGIN
					u.ch=slider
					IF u.ch LE u.stat.nch THEN BEGIN
						IF u.stat.dat THEN BEGIN
							igood=*u.dat.good[u.index]
							IF igood[u.ch-1] THEN widget_control,u.id.chgood,set_button=1 ELSE widget_control,u.id.chgood,set_button=0
							plot_moments,u
							plot_profiles,u
						ENDIF
						update_ch_text,u
					ENDIF ELSE u.ch=u.stat.nch			
				END
				'EM_EPS_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						u.dat.eps[0,u.index]=float(slider/100.0)
						widget_control,u.id.em_eps_text,set_value=num2str(u.dat.eps[0,u.index],dp=2)
						invert_current_profile,u
						plot_moments,u
						plot_profiles,u	
					ENDIF	
				END
				'EM_ETA_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						u.dat.eta[0,u.index]=10^(slider/10.0)
						widget_control,u.id.em_eta_text,set_value=num2str(slider/10.0,dp=2)
						invert_current_profile,u
						plot_moments,u
						plot_profiles,u	
					ENDIF	
				END
				'VT_EPS_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						u.dat.eps[1,u.index]=float(slider/100.0)
						widget_control,u.id.vt_eps_text,set_value=num2str(u.dat.eps[1,u.index],dp=2)
						invert_current_profile,u
						plot_moments,u
						plot_profiles,u	
					ENDIF	
				END
				'VT_ETA_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						u.dat.eta[1,u.index]=10^(slider/10.0)
						widget_control,u.id.vt_eta_text,set_value=num2str(slider/10.0,dp=2)
						invert_current_profile,u
						plot_moments,u
						plot_profiles,u	
					ENDIF	
				END
				'TI_EPS_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						u.dat.eps[2,u.index]=float(slider/100.0)
						widget_control,u.id.ti_eps_text,set_value=num2str(u.dat.eps[2,u.index],dp=2)
						invert_current_profile,u
						plot_moments,u
						plot_profiles,u	
					ENDIF	
				END
				'TI_ETA_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						u.dat.eta[2,u.index]=10^(slider/10.0)
						widget_control,u.id.ti_eta_text,set_value=num2str(slider/10.0,dp=2)
						invert_current_profile,u
						plot_moments,u
						plot_profiles,u	
					ENDIF	
				END
				'GLOW_SLIDER' : BEGIN
					IF u.stat.dat AND slider LT u.stat.good[1] THEN BEGIN
						imom=*u.dat.mom[u.index]
						IF slider GT u.stat.good[0] THEN BEGIN
							imom[0:slider-1,6]=0
						ENDIF ELSE BEGIN
							up=u.stat.good[1] < u.stat.nch-1
							IF u.stat.gslide THEN BEGIN
								imom[slider:up,6]=1
							ENDIF ELSE BEGIN
								igood=*u.dat.good[u.index]
								imom[slider:up,6]=igood[slider:up]
							ENDELSE
						ENDELSE
						widget_control,u.id.glow_text,set_value=num2str(u.stat.good[0]+1,1)
						u.stat.good[0]=slider
						*u.dat.mom[u.index]=imom
						invert_current_profile,u
						plot_moments,u
						plot_profiles,u		
					ENDIF
				END
				'GHIGH_SLIDER' : BEGIN
					IF u.stat.dat AND slider GT u.stat.good[0] AND slider LT u.stat.nch THEN BEGIN
						imom=*u.dat.mom[u.index]
						IF slider LT u.stat.good[1] THEN BEGIN
							imom[slider:u.stat.nch-1,6]=0
						ENDIF ELSE BEGIN
							low=u.stat.good[0] > 0
							IF u.stat.gslide THEN BEGIN
								imom[low:slider,6]=1
							ENDIF ELSE BEGIN
								igood=*u.dat.good[u.index]
								imom[low:slider,6]=igood[low:slider]
							ENDELSE
						ENDELSE
						widget_control,u.id.ghigh_text,set_value=num2str(u.stat.good[1]+1,1)
						u.stat.good[1]=slider
						*u.dat.mom[u.index]=imom
						invert_current_profile,u
						plot_moments,u
						plot_profiles,u		
					ENDIF ELSE widget_control,u.id.ghigh_slider,set_value=u.stat.nch
				END

				'T_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						index=ipt(u.dat.tau,slider/1.0e3)
						tmp=where(u.dat.tau GT 0)
						IF slider/1.0e3 GE max(u.dat.tau[tmp]) THEN index=n(tmp)
						IF slider/1.0e3 LE min(u.dat.tau[tmp]) THEN index=0
						IF u.stat.time NE u.dat.tau[index] THEN BEGIN
							u.stat.time=u.dat.tau[index]
							u.index=index
							update_temporal_text,u
							update_ch_text,u
							widget_control,u.id.time,set_value=num2str(u.stat.time,dp=2)
							plot_profiles,u
							plot_moments,u
							IF u.dat.tgood[u.index] THEN widget_control,u.id.tgood,set_button=1 ELSE widget_control,u.id.tgood,set_button=0
							igood=*u.dat.good[u.index]
							IF igood[u.ch-1] THEN widget_control,u.id.chgood,set_button=1 ELSE widget_control,u.id.chgood,set_button=0
						ENDIF
					ENDIF
				END
				'TINST_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						u.dat.tinst=slider/1.0e3
						widget_control,u.id.tinst_text,set_value=num2str(int(slider),1)
						plot_profiles,u
					ENDIF
				END
				'WINST_SLIDER' : BEGIN
	
				END
				'NPTS_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						npts=slider
						widget_control,u.id.rmin_slider,get_value=rmin
						widget_control,u.id.rmax_slider,get_value=rmax
						irho=make(rmin/100.0,rmax/100.0,npts)
						*u.dat.rho[u.index]=irho
						update_current_voxel,u
				 		invert_current_profile,u
						plot_profiles,u
						plot_moments,u	
						widget_control,u.id.npts_text,set_value=num2str(npts,1)
						update_rho_text,u,/str
					ENDIF		
				END
				'RMIN_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						rmin=slider
						widget_control,u.id.npts_slider,get_value=npts
						widget_control,u.id.rmax_slider,get_value=rmax
						irho=make(rmin/100.0,rmax/100.0,npts)
						*u.dat.rho[u.index]=irho
						update_current_voxel,u
				 		invert_current_profile,u
						plot_profiles,u
						plot_moments,u
						widget_control,u.id.rmin_text,set_value=num2str(rmin/100.0,dp=2)
						update_rho_text,u,/str
					ENDIF		
				END
				'RMAX_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						rmax=slider
						widget_control,u.id.npts_slider,get_value=npts
						widget_control,u.id.rmin_slider,get_value=rmin
						irho=make(rmin/100.0,rmax/100.0,npts)
						*u.dat.rho[u.index]=irho
						update_current_voxel,u
				 		invert_current_profile,u
						plot_profiles,u
						plot_moments,u
						widget_control,u.id.rmax_text,set_value=num2str(rmax/100.0,dp=2)
						update_rho_text,u,/str
					ENDIF		
				END
				ELSE:
			ENDCASE
		END
   		"WIDGET_TEXT_CH": BEGIN
			CASE event.id OF 
			u.id.shotid : BEGIN
				widget_control,u.id.shotid,get_value=shot
				u.shot=shot
                        END
			u.id.lineid : BEGIN
				widget_control,u.id.lineid,get_value=line
				u.stat.line=int(line)
			END
			u.id.thtid : BEGIN
				widget_control,u.id.thtid,get_value=tht
				IF hirexsr_is_analysis(u.shot,int(tht[0])) THEN BEGIN
					u.tht=int(tht[0])
				ENDIF ELSE BEGIN
					widget_control,u.id.message,set_value='invalid THT',/app
					widget_control,u.id.thtid,set_value=num2str(u.tht,1)
				ENDELSE
			END
			u.id.note : BEGIN
				widget_control,u.id.note,get_value=note	
				IF is_dat(u) THEN u.dat.note=note
			END
			u.id.rmax_text : BEGIN
				widget_control,u.id.rmax_text,get_value=rmax
				rmax*=100.0
				widget_control,u.id.rmax_slider,set_value=rmax
				IF u.stat.dat THEN BEGIN
					widget_control,u.id.npts_slider,get_value=npts
					widget_control,u.id.rmin_slider,get_value=rmin
					irho=make(rmin/100.0,rmax/100.0,npts)
					*u.dat.rho[u.index]=irho
					update_current_voxel,u
			 		invert_current_profile,u
					plot_profiles,u
					plot_moments,u
					update_rho_text,u,/str
				ENDIF
			END
			u.id.rmin_text : BEGIN
				widget_control,u.id.rmin_text,get_value=rmin
				rmin*=100.0
				widget_control,u.id.rmin_slider,set_value=rmin
				IF u.stat.dat THEN BEGIN
					widget_control,u.id.npts_slider,get_value=npts
					widget_control,u.id.rmax_slider,get_value=rmax
					irho=make(rmin/100.0,rmax/100.0,npts)
					*u.dat.rho[u.index]=irho
					update_current_voxel,u
			 		invert_current_profile,u
					plot_profiles,u
					plot_moments,u
					update_rho_text,u,/str
				ENDIF
			END
			u.id.rho_text : BEGIN
				widget_control,u.id.rho_text,get_value=rhostr
				irho=float(strsplit(rhostr,',',/extract))
				IF u.stat.dat THEN BEGIN	
					*u.dat.rho[u.index]=irho
					update_current_voxel,u
			 		invert_current_profile,u
					plot_profiles,u
					plot_moments,u
					update_rho_text,u
				ENDIF
				
			END
			u.id.xmin : BEGIN
				widget_control,u.id.xmin,get_value=value
				u.stat.plot.rho[0]=float(value)
				plot_profiles,u
			END
			u.id.xmax : BEGIN
				widget_control,u.id.xmax,get_value=value
				u.stat.plot.rho[1]=float(value)
				plot_profiles,u
			END
			u.id.emmin : BEGIN
				widget_control,u.id.emmin,get_value=value
				u.stat.plot.em[0]=float(value)
				plot_profiles,u
			END
			u.id.emmax : BEGIN
				widget_control,u.id.emmax,get_value=value
				u.stat.plot.em[1]=float(value)
				plot_profiles,u
			END
			u.id.vtmin : BEGIN
				widget_control,u.id.vtmin,get_value=value
				u.stat.plot.vt[0]=float(value)
				plot_profiles,u
			END
			u.id.vtmax : BEGIN
				widget_control,u.id.vtmax,get_value=value
				u.stat.plot.vt[1]=float(value)
				plot_profiles,u
			END
			u.id.vpmin : BEGIN
				widget_control,u.id.vpmin,get_value=value
				u.stat.plot.vp[0]=float(value)
				plot_profiles,u
			END
			u.id.vpmax : BEGIN
				widget_control,u.id.vpmax,get_value=value
				u.stat.plot.vp[1]=float(value)
				plot_profiles,u
			END
			u.id.timin : BEGIN
				widget_control,u.id.timin,get_value=value
				u.stat.plot.ti[0]=float(value)
				plot_profiles,u
			END
			u.id.timax : BEGIN
				widget_control,u.id.timax,get_value=value
				u.stat.plot.ti[1]=float(value)
				plot_profiles,u
			END
			u.id.m0min : BEGIN
				widget_control,u.id.m0min,get_value=value
				u.stat.plot.m0[0]=float(value)
				plot_moments,u
			END
			u.id.m0max : BEGIN
				widget_control,u.id.m0max,get_value=value
				u.stat.plot.m0[1]=float(value)
				plot_moments,u
			END
			u.id.m1min : BEGIN
				widget_control,u.id.m1min,get_value=value
				u.stat.plot.m1[0]=float(value)
				plot_moments,u
			END
			u.id.m1max : BEGIN
				widget_control,u.id.m1max,get_value=value
				u.stat.plot.m1[1]=float(value)
				plot_moments,u
			END
			u.id.m2min : BEGIN
				widget_control,u.id.m2min,get_value=value
				u.stat.plot.m2[0]=float(value)
				plot_moments,u
			END
			u.id.m2max : BEGIN
				widget_control,u.id.m2max,get_value=value
				u.stat.plot.m2[1]=float(value)
				plot_moments,u
			END
			u.id.em_eps_text : BEGIN
				widget_control,u.id.em_eps_text,get_value=value
				widget_control,u.id.em_eps_slider,set_value=int(value)*100
				IF u.stat.dat THEN BEGIN
						u.dat.eps[0,u.index]=float(value)
						widget_control,u.id.em_eps_text,set_value=num2str(u.dat.eps[0,u.index],dp=2)
						invert_current_profile,u
						plot_moments,u
						plot_profiles,u	
				ENDIF
			END
			u.id.vt_eps_text : BEGIN
				widget_control,u.id.vt_eps_text,get_value=value
				widget_control,u.id.vt_eps_slider,set_value=int(value)*100
				IF u.stat.dat THEN BEGIN
						u.dat.eps[1,u.index]=float(value)
						widget_control,u.id.vt_eps_text,set_value=num2str(u.dat.eps[1,u.index],dp=2)
						invert_current_profile,u
						plot_moments,u
						plot_profiles,u	
				ENDIF
			END
			u.id.ti_eps_text : BEGIN
				widget_control,u.id.ti_eps_text,get_value=value
				widget_control,u.id.ti_eps_slider,set_value=int(value)*100
				IF u.stat.dat THEN BEGIN
						u.dat.eps[2,u.index]=float(value)
						widget_control,u.id.ti_eps_text,set_value=num2str(u.dat.eps[2,u.index],dp=2)
						invert_current_profile,u
						plot_moments,u
						plot_profiles,u	
				ENDIF
			END
			ENDCASE
		END
		ELSE:
	ENDCASE
	IF button NE 'QUIT' THEN widget_control,event.top,set_uvalue=u		
END

;+
;NAME:
;	W_HIREXSR_PROFILES
;
;-

PRO w_hirexsr_profiles,shot=shot,tht=tht,line=line,loaddata=loaddata
	user=logname()
	loadct,12,/silent
	mdim=get_screen_size()
	IF mdim[0] LT 1600 OR mdim[1] LT 1100 THEN base=widget_base(title='HIREXSR GENSPEC ANALYSIS',/row,tlb_size_events=1,/scroll,$
		x_scroll_size=mdim[0]*0.95,y_scroll_size=mdim[1]*0.85) ELSE base=widget_base(title='HIREXSR GENSPEC ANALYSIS',/row,tlb_size_events=1)
	A=widget_base(base,/column)
	B=widget_base(base,/row)
	space=widget_base(base,/row,xsize=2)	
	C=widget_base(base,/column)
	space=widget_base(base,/row,xsize=2)	
	D=widget_base(base,/column)

	dum = widget_label(A,value='MOMENT PROFIELS')
	A1=widget_base(A,frame=5)
	draw1=widget_draw(A1,xsize=500,ysize=1050)
	;A2=widget_base(A,frame=5)
	;draw2=widget_draw(A2,xsize=500,ysize=348)
	;A3=widget_base(A,frame=5)
	;draw3=widget_draw(A3,xsize=500,ysize=348)

	B1=widget_base(B,/column)
	dum = widget_label(B1,value='CH #')
	ch_slider=widget_slider(B1,ysize=330*3.14,min=1,max=24*4,value=24.,/drag,/vert,/supp)
	ch_text=widget_text(B1,xsize=2,ysize=1)
	B2=widget_base(B,/column)
	dum = widget_label(B2,value='GLOW')
	glow_slider=widget_slider(B2,ysize=330*3.14,min=0,max=97,value=0.,/drag,/vert,/supp)
	glow_text=widget_text(B2,xsize=2,ysize=1)
	B3=widget_base(B,/column)
	dum = widget_label(B3,value='GHIGH')
	ghigh_slider=widget_slider(B3,ysize=330*3.14,min=0,max=97,value=97.,/drag,/vert,/supp)
	ghigh_text=widget_text(B3,xsize=2,ysize=1)

	

	dum = widget_label(C,value='INPUT/OUTPUT SETUP')
	xsize=325
	;setup tabs
	Ctb=widget_tab(C,location=3)
	Ctb1=widget_base(Ctb,title=' MOMENTS ',/column,group_leader=base,/frame)
	Ctb2=widget_base(Ctb,title=' PLOTTING ',/column,group_leader=base,/frame)

	;moment setup options
	C1=widget_base(Ctb1,xsize=xsize,ysize=400,/column)
	linelist = widget_text(C1,xsize=45,ysize=12)

; 	dum=widget_label(C1,value='Moment Setup')
; 	C1p1=widget_base(C1,/row)
; 	dum=widget_label(C1p1,value='LINE: ')
; 	dum=widget_label(C1p1,value='  Ar XVII:')
; 	C1p1p1=widget_base(C1p1,/row,/nonexcl)
; 	inv_w=widget_button(C1p1p1,value=' w ')
; 	inv_x=widget_button(C1p1p1,value=' x ')
; 	inv_z=widget_button(C1p1p1,value=' z ')
; 	C1p2=widget_base(C1,/row)
; 	dum=widget_label(C1p2,value='     ')
; 	dum=widget_label(C1p2,value='  Ar XVIII:')
; 	C1p2p1=widget_base(C1p2,/row,/nonexcl)
; 	inv_lya=widget_button(C1p2p1,value=' Lya ')
; 	dum=widget_label(C1p2,value='     Mo XXXIII:')
; 	C1p2p2=widget_base(C1p2,/row,/nonexcl)
; 	inv_m32=widget_button(C1p2p2,value=' 4c ')
; 	C1p3=widget_base(C1,/row)
; 	dum=widget_label(C1p3,value='     ')
; 	dum=widget_label(C1p3,value='    Ca XIX:')
; 	C1p3p1=widget_base(C1p3,/row,/nonexcl)
; 	inv_caw=widget_button(C1p3p1,value=' w ')
	

	;time spatial selection
	C1x1=widget_base(C1,/row)
	dum = widget_label(C1x1,value='TIME PT ')
	itpt=widget_text(C1x1,xsize=3,ysize=1)
	dum = widget_label(C1x1,value=' OF ')
	ntpt=widget_text(C1x1,xsize=3,ysize=1)
	dum = widget_label(C1x1,value='  ')
	mt=widget_button(C1x1,value=' - ')
	pt=widget_button(C1x1,value=' + ')
	C1x1p1=widget_base(C1x1,/row,/nonexcl)
	tgood=widget_button(C1x1p1,value=' GOOD')
	C1x2=widget_base(C1,/row)
	dum = widget_label(C1x2,value='    CH# ')
	ch=widget_text(C1x2,xsize=2,ysize=1)
	dum = widget_label(C1x2,value='  OF  ')
	chtot=widget_text(C1x2,xsize=2,ysize=1)
	dum = widget_label(C1x2,value='  ')
	mch=widget_button(C1x2,value=' - ')
	pch=widget_button(C1x2,value=' + ')
	C1x2p1=widget_base(C1x2,/row,/nonexcl)
	chgood=widget_button(C1x2p1,value=' GOOD')
	C1x3=widget_base(C1,/row)
	regood=widget_button(C1x3,value='REMOVE')
	dum=widget_label(C1x3, value='CH FOR ALL TIMES  ')
	copygood=widget_button(C1x3,value='COPY')
	dum=widget_label(C1x3, value='GOOD TO ALL TIMES')
	C1x4=widget_base(C1,/row)
	addgood=widget_button(C1x4,value='ADD')
	dum=widget_label(C1x4, value='CH FOR ALL TIMES  ')
	C1x5=widget_base(C1,/row,/nonexcl)
	gslide=widget_button(C1x5,value='USE GOOD FROM SLIDERS')
	C1x6=widget_base(C1,/row,/nonexcl)
	m1sinv=widget_button(C1x6,value='INCLUDE m=1 SINE TERMS IN INVERSION')


	;plotting setup options
	Px=widget_base(Ctb2,xsize=xsize,ysize=400,/column)
	P1=widget_base(Px, /row)
	P1a=widget_base(P1,/row,/nonexcl)
	xauto=widget_button(P1a,value='AUTOSCALE   ')
	xmin=widget_text(P1,xsize=5,ysize=1,/edit)
	dum=widget_label(P1,value=' < r/a < ')
	xmax=widget_text(P1,xsize=5,ysize=1,/edit)
	P2=widget_base(Px, /row)
	P2a=widget_base(P2,/row,/nonexcl)
	emauto=widget_button(P2a,value='AUTOSCALE   ')
	emmin=widget_text(P2,xsize=5,ysize=1,/edit)
	dum=widget_label(P2,value=' < EMISS < ')
	emmax=widget_text(P2,xsize=5,ysize=1,/edit)
	P3=widget_base(Px, /row)
	P3a=widget_base(P3,/row,/nonexcl)
	vtauto=widget_button(P3a,value='AUTOSCALE   ')
	vtmin=widget_text(P3,xsize=5,ysize=1,/edit)
	dum=widget_label(P3,value=' < VTOR < ')
	vtmax=widget_text(P3,xsize=5,ysize=1,/edit)	
	P4=widget_base(Px, /row)
	P4a=widget_base(P4,/row,/nonexcl)
	vpauto=widget_button(P4a,value='AUTOSCALE   ')
	vpmin=widget_text(P4,xsize=5,ysize=1,/edit)
	dum=widget_label(P4,value=' < VPOL < ')
	vpmax=widget_text(P4,xsize=5,ysize=1,/edit)
	P41=widget_base(Px, /row)
	P41a=widget_base(P41,/row,/nonexcl)
	tiauto=widget_button(P41a,value='AUTOSCALE   ')
	timin=widget_text(P41,xsize=5,ysize=1,/edit)
	dum=widget_label(P41,value=' < TI < ')
	timax=widget_text(P41,xsize=5,ysize=1,/edit)
	P5=widget_base(Px, /row)
	P5a=widget_base(P5,/row,/nonexcl)
	m0auto=widget_button(P5a,value='AUTOSCALE   ')
	m0min=widget_text(P5,xsize=5,ysize=1,/edit)
	dum=widget_label(P5,value=' < M0 < ')
	m0max=widget_text(P5,xsize=5,ysize=1,/edit)	
	P6=widget_base(Px, /row)
	P6a=widget_base(P6,/row,/nonexcl)
	m1auto=widget_button(P6a,value='AUTOSCALE   ')
	m1min=widget_text(P6,xsize=5,ysize=1,/edit)
	dum=widget_label(P6,value=' < M1 < ')
	m1max=widget_text(P6,xsize=5,ysize=1,/edit)	
	P7=widget_base(Px, /row)
	P7a=widget_base(P7,/row,/nonexcl)
	m2auto=widget_button(P7a,value='AUTOSCALE   ')
	m2min=widget_text(P7,xsize=5,ysize=1,/edit)
	dum=widget_label(P7,value=' < M2 < ')
	m2max=widget_text(P7,xsize=5,ysize=1,/edit)
	P8=widget_base(Px, /row,/nonexcl)
	tsplot=widget_button(P8,value='OPLOT Te (TS)')
	eceplot=widget_button(P8,value='  OPLOT Te (GPC)')
	P9a=widget_base(Px,/row)
	dum=widget_label(P9a,value='OPLOT LINE INTEGRATED ')
	P9=widget_base(P9a, /row,/nonexcl)
	limom=widget_button(P9,value='MOMENTS')
	lifit=widget_button(P9,value='FITS')

	;inversion setup options
	Ctc=widget_tab(C,location=3)
	Ctc1=widget_base(Ctc,title=' RHO ',/column,group_leader=base,/frame)
	Ctc2=widget_base(Ctc,title=' EPS/ETA ',/column,group_leader=base,/frame)

	C2=widget_base(Ctc1,xsize=xsize,ysize=415,/column)
	dum=widget_label(C2,value='Inversion Setup')
	C2p1=widget_base(C2,/row)
	dum=widget_label(C2p1,value='RHOMIN: ')
	rmin_slider=widget_slider(C2p1,xsize=200,min=0,max=120,value=0,/drag,/suppress)
	rmin_text=widget_text(C2p1,xsize=4,ysize=1,/edit)
	C2p2=widget_base(C2,/row)
	dum=widget_label(C2p2,value='RHOMAX: ')
	rmax_slider=widget_slider(C2p2,xsize=200,min=0,max=120,value=100,/drag,/suppress)
	rmax_text=widget_text(C2p2,xsize=4,ysize=1,/edit)
	C2p3=widget_base(C2,/row)
	dum=widget_label(C2p3,value='NPTS:   ')
	npts_slider=widget_slider(C2p3,xsize=200,min=0,max=96,value=25,/drag,/suppress)
	npts_text=widget_text(C2p3,xsize=4,ysize=1,/edit)
	C2p17=widget_base(C2,/row)
	dum=widget_label(C2p17,value='RHO: ')
	rho_text=widget_text(C2p17,xsize=41,ysize=1,/edit)
	C2p16=widget_base(C2,/row)
	copyrho=widget_button(C2p16,value='COPY')
	dum=widget_label(C2p16, value='RHO TO ALL TIMES ')
	C2p16x=widget_base(C2p16,/row,/nonexcl)
	m1scalc=widget_button(C2p16x,value='CALCULATE m=1 TERMS')
	C2p13=widget_base(C2,/row)
	dum=widget_label(C2p13,value='CALCULATE VOXELS FOR: ')
	voxcurr=widget_button(C2p13,value='CURRENT')
	voxall=widget_button(C2p13,value=' ALL ')
	space=widget_base(C2,/row,ysize=10)
	C2p12=widget_base(C2,/row)
	dum=widget_label(C2p12,value='DO INVERSIONS FOR: ')
	invcurr=widget_button(C2p12,value='CURRENT')
	invall=widget_button(C2p12,value=' ALL ')
	C2p16=widget_base(C2,/row)
	dum=widget_label(C2p16,value='                               DC      CALIB')
	C2p15=widget_base(C2,/row)
	dum=widget_label(C2p15,value='V_INST: ')
	winst_slider=widget_slider(C2p15,xsize=115,min=-1000,max=1000,value=0,/drag,/suppress)
	winst_text=widget_text(C2p15,xsize=5,ysize=1,/edit)
	ishift_text=widget_text(C2p15,xsize=5,ysize=1)
	dum=widget_label(C2p15,value='[km/s]')
	C2p14=widget_base(C2,/row)
	dum=widget_label(C2p14,value='T_INST: ')
	tinst_slider=widget_slider(C2p14,xsize=115,min=0,max=1000,value=0,/drag,/suppress)
	tinst_text=widget_text(C2p14,xsize=5,ysize=1,/edit)
	iwidth_text=widget_text(C2p14,xsize=5,ysize=1)
	dum=widget_label(C2p14,value='[eV]')
	
	C2x=widget_base(Ctc2,xsize=xsize,ysize=415,/column)
	C2p4=widget_base(C2x,/row)
	dum=widget_label(C2p4,value='EMISS: ')
	C2p4p1=widget_base(C2p4,/row,/nonexcl)
	invem=widget_button(C2p4p1,value=' ')
	dum=widget_label(C2p4,value='ETA: ')
	em_eta_slider=widget_slider(C2p4,xsize=100,min=-20,max=50,value=10,/drag,/suppress)
	dum=widget_label(C2p4,value=' 10^')
	em_eta_text=widget_text(C2p4,xsize=4,ysize=1,/edit)
	C2p5=widget_base(C2x,/row)
	dum=widget_label(C2p5,value='              ')
	dum=widget_label(C2p5,value='EPS: ')
	em_eps_slider=widget_slider(C2p5,xsize=100,min=1,max=1000,value=10,/drag,/suppress)
	dum=widget_label(C2p5,value='    ')
	em_eps_text=widget_text(C2p5,xsize=4,ysize=1,/edit)
	space=widget_base(C2x,/row,ysize=10)
	C2p6=widget_base(C2x,/row)
	dum=widget_label(C2p6,value='VTOR:  ')
	C2p6p1=widget_base(C2p6,/row,/nonexcl)
	invvt=widget_button(C2p6p1,value=' ')
	dum=widget_label(C2p6,value='ETA: ')
	vt_eta_slider=widget_slider(C2p6,xsize=100,min=-20,max=50,value=10,/drag,/suppress)
	dum=widget_label(C2p6,value=' 10^')
	vt_eta_text=widget_text(C2p6,xsize=4,ysize=1,/edit)
	C2p7=widget_base(C2x,/row)
	dum=widget_label(C2p7,value='VPOL:  ')
	C2p7p1=widget_base(C2p7,/row,/nonexcl)
	invvp=widget_button(C2p7p1,value=' ')
	dum=widget_label(C2p7,value='EPS: ')
	vt_eps_slider=widget_slider(C2p7,xsize=100,min=0,max=1000,value=10,/drag,/suppress)
	dum=widget_label(C2p7,value='    ')
	vt_eps_text=widget_text(C2p7,xsize=4,ysize=1,/edit)
	space=widget_base(C2x,/row,ysize=10)
	C2p10=widget_base(C2x,/row)
	dum=widget_label(C2p10,value='TI:    ')
	C2p10p1=widget_base(C2p10,/row,/nonexcl)
	invti=widget_button(C2p10p1,value=' ')
	dum=widget_label(C2p10,value='ETA: ')
	ti_eta_slider=widget_slider(C2p10,xsize=100,min=-20,max=50,value=10,/drag,/suppress)
	dum=widget_label(C2p10,value=' 10^')
	ti_eta_text=widget_text(C2p10,xsize=4,ysize=1,/edit)
	C2p11=widget_base(C2x,/row)
	dum=widget_label(C2p11,value='              ')
	dum=widget_label(C2p11,value='EPS: ')
	ti_eps_slider=widget_slider(C2p11,xsize=100,min=1,max=1000,value=10,/drag,/suppress)
	dum=widget_label(C2p11,value='    ')
	ti_eps_text=widget_text(C2p11,xsize=4,ysize=1,/edit)

	;tree input output
	C3=widget_base(C,frame=2,xsize=xsize,ysize=185,/column)
	dum=widget_label(C3,value='Tree I/O')
	C3p1=widget_base(C3,/row)
	dum = widget_label(C3p1,value='SHOT:')
	shotid = widget_text(C3p1,xsize=10,ysize=1,/edit)
	dum = widget_label(C3p1,value='L#:')
	lineid = widget_text(C3p1,xsize=2,ysize=1,/edit)
	dum = widget_label(C3p1,value='THT:')
	thtid = widget_text(C3p1,xsize=2,ysize=1,/edit)
	dum = widget_label(C3p1,value='')
	load= widget_button(C3p1,value='LOAD')
	save= widget_button(C3p1,value='SAVE')
	C3p1x=widget_base(C3,/row)
	dum = widget_label(C3p1x,value='NOTE: ')
	note = widget_text(C3p1x,xsize=17,ysize=1,/edit)
	quit= widget_button(C3p1x,value='QUIT')
	print= widget_button(C3p1x,value='PRINT')
	stop= widget_button(C3p1x,value='STOP')
	C3p2=widget_base(C3,/row)
	message = widget_text(C3p2,xsize=45,ysize=3,/scroll)

	C4=widget_base(C,/row)
	dum=widget_label(C4,value='T: ')
	time=widget_text(C4,xsize=5,ysize=1,/edit)
	dum=widget_label(C4,value='[msec]')
	t_slider=widget_slider(C4,xsize=230,min=0,max=2000,value=1000,/drag,/suppress)

	dum = widget_label(D,value='INVERTED PROFILES')
	D1=widget_base(D,frame=5)
	draw4=widget_draw(D1,xsize=550,ysize=1050)
	;D2=widget_base(D,frame=5)
	;draw5=widget_draw(D2,xsize=550,ysize=258)
	;D3=widget_base(D,frame=5)
	;draw6=widget_draw(D3,xsize=550,ysize=258)
	;D4=widget_base(D,frame=5)
	;draw7=widget_draw(D4,xsize=550,ysize=258)

	;build u structure
	id={base:base,draw1:draw1,draw4:draw4,$
		ch_slider:ch_slider,glow_slider:glow_slider,ghigh_slider:ghigh_slider,$
		ch_text:ch_text,glow_text:glow_text,ghigh_text:ghigh_text,$
		;inv_w:inv_w,inv_x:inv_x,inv_z:inv_z,inv_lya:inv_lya,inv_m32:inv_m32,inv_caw:inv_caw,$
		linelist:linelist,itpt:itpt,ntpt:ntpt,mt:mt,pt:pt,tgood:tgood,ch:ch,chtot:chtot,mch:mch,pch:pch,chgood:chgood,$
		regood:regood,copygood:copygood,addgood:addgood,gslide:gslide,m1sinv:m1sinv,$
		xauto:xauto,emauto:emauto,vtauto:vtauto,vpauto:vpauto,tiauto:tiauto,m0auto:m0auto,m1auto:m1auto,m2auto:m2auto,$
		xmin:xmin,emmin:emmin,vtmin:vtmin,vpmin:vpmin,timin:timin,m0min:m0min,m1min:m1min,m2min:m2min,$
		xmax:xmax,emmax:emmax,vtmax:vtmax,vpmax:vpmax,timax:timax,m0max:m0max,m1max:m1max,m2max:m2max,$
		tsplot:tsplot,eceplot:eceplot,limom:limom,lifit:lifit,$
		rmin_slider:rmin_slider,rmin_text:rmin_text,rmax_slider:rmax_slider,rmax_text:rmax_text,rho_text:rho_text,copyrho:copyrho,m1scalc:m1scalc,$
		npts_slider:npts_slider,npts_text:npts_text,$
		invem:invem,invvt:invvt,invvp:invvp,invti:invti,$
		em_eta_slider:em_eta_slider,em_eta_text:em_eta_text,em_eps_slider:em_eps_slider,em_eps_text:em_eps_text,$
		vt_eta_slider:vt_eta_slider,vt_eta_text:vt_eta_text,vt_eps_slider:vt_eps_slider,vt_eps_text:vt_eps_text,$
		ti_eta_slider:ti_eta_slider,ti_eta_text:ti_eta_text,ti_eps_slider:ti_eps_slider,ti_eps_text:ti_eps_text,$
		tinst_slider:tinst_slider,tinst_text:tinst_text,iwidth_text:iwidth_text,winst_slider:winst_slider,winst_text:winst_text,ishift_text:ishift_text,$
		voxcurr:voxcurr,voxall:voxall,invcurr:invcurr,invall:invall,$
		shotid:shotid,lineid:lineid,thtid:thtid,save:save,load:load,quit:quit,print:print,stop:stop,message:message,note:note,$
		time:time,t_slider:t_slider}

	;set defaults
	IF NOT keyword_set(shot) THEN shot=1101209012
	IF NOT keyword_set(line) THEN line=0
	IF NOT keyword_set(tht) THEN tht=0
	time=1.0
	ch=24
	index=0
	rplt=[0.0,0.0]
	emplt=[0.0,0.0]
	vtplt=[-25.0,25.0]
	vpplt=[-5.0,5.0]
	tiplt=[0.0,2.5]
	m0plt=[0.0,0.0]
	m1plt=[0.0,0.0]
	m2plt=[0.0,0.0]
	plot={rho:rplt,em:emplt,vt:vtplt,vp:vpplt,ti:tiplt,m0:m0plt,m1:m1plt,m2:m2plt,ts:1,ece:0,limom:1,lifit:1,msize:[500,1050],psize:[550,1050]}
	auto={rho:1,em:1,vt:0,vp:0,ti:0,m0:1,m1:1,m2:1}
	inv={em:1,vt:1,vp:0,ti:1,m1s:0}
	stat={plot:plot,auto:auto,ps:0,line:line,time:time,ntime:0,nch:0,dat:0,good:[0,0],gslide:0,col:[255,50,200,80,150,200],pscol:[0,30,200,100,150,200],m1s:0,m1scalc:0}

	u={id:id,shot:shot,tht:tht,ch:ch,index:index,stat:stat,inv:inv}
	widget_control,base,set_uvalue=u
	widget_control,u.id.shotid,set_value=num2str(u.shot,1)
	widget_control,u.id.lineid,set_value=num2str(u.stat.line,1)
	widget_control,u.id.thtid,set_value=num2str(u.tht,1)
	widget_control,u.id.t_slider,set_value=u.stat.time*1.0e3
	widget_control,u.id.time,set_value=num2str(u.stat.time,dp=2)
	widget_control,u.id.ch_slider,set_value=u.ch
	widget_control,u.id.invem,set_button=u.inv.em
	widget_control,u.id.invvt,set_button=u.inv.vt
	widget_control,u.id.invvp,set_button=u.inv.vp
	widget_control,u.id.invti,set_button=u.inv.ti
;	invid=[u.id.inv_w,u.id.inv_x,u.id.inv_z,u.id.inv_lya,u.id.inv_m32]
;	widget_control,invid[u.stat.line],set_button=1
	widget_control,u.id.ch_slider,get_value=slider
	widget_control,u.id.ch_text,set_value=num2str(slider,1)
	widget_control,u.id.glow_slider,get_value=slider
	widget_control,u.id.glow_text,set_value=num2str(slider,1)
	widget_control,u.id.ghigh_slider,get_value=slider
	widget_control,u.id.ghigh_text,set_value=num2str(slider,1)
	widget_control,u.id.tsplot, set_button=u.stat.plot.ts
	widget_control,u.id.limom, set_button=u.stat.plot.limom
	widget_control,u.id.lifit, set_button=u.stat.plot.lifit
	IF u.stat.auto.rho THEN widget_control,u.id.xauto, set_button=1 ELSE BEGIN
		widget_control,u.id.xmin,set_value=num2str(u.stat.plot.rho[0],dp=2)
		widget_control,u.id.xmax,set_value=num2str(u.stat.plot.rho[1],dp=2)
	ENDELSE
	IF u.stat.auto.em  THEN widget_control,u.id.emauto, set_button=1 ELSE BEGIN
		widget_control,u.id.emmin,set_value=num2str(u.stat.plot.em[0],dp=2)
		widget_control,u.id.emmax,set_value=num2str(u.stat.plot.em[1],dp=2)
	ENDELSE
	IF u.stat.auto.vt THEN widget_control,u.id.vtauto, set_button=1 ELSE BEGIN
		widget_control,u.id.vtmin,set_value=num2str(u.stat.plot.vt[0],dp=2)
		widget_control,u.id.vtmax,set_value=num2str(u.stat.plot.vt[1],dp=2)
	ENDELSE
	IF u.stat.auto.vp THEN widget_control,u.id.vpauto, set_button=1 ELSE BEGIN
		widget_control,u.id.vpmin,set_value=num2str(u.stat.plot.vp[0],dp=2)
		widget_control,u.id.vpmax,set_value=num2str(u.stat.plot.vp[1],dp=2)
	ENDELSE
	IF u.stat.auto.ti THEN widget_control,u.id.tiauto, set_button=1 ELSE BEGIN
		widget_control,u.id.timin,set_value=num2str(u.stat.plot.ti[0],dp=2)
		widget_control,u.id.timax,set_value=num2str(u.stat.plot.ti[1],dp=2)
	ENDELSE
	IF u.stat.auto.m0 THEN widget_control,u.id.m0auto, set_button=1 ELSE BEGIN
		widget_control,u.id.m0min,set_value=num2str(u.stat.plot.m0[0],dp=2)
		widget_control,u.id.m0max,set_value=num2str(u.stat.plot.m0[1],dp=2)
	ENDELSE
	IF u.stat.auto.m1 THEN widget_control,u.id.m1auto, set_button=1 ELSE BEGIN
		widget_control,u.id.m1min,set_value=num2str(u.stat.plot.m1[0],dp=2)
		widget_control,u.id.m1max,set_value=num2str(u.stat.plot.m1[1],dp=2)
	ENDELSE
	IF u.stat.auto.m2 THEN widget_control,u.id.m2auto, set_button=1 ELSE BEGIN
		widget_control,u.id.m2min,set_value=num2str(u.stat.plot.m2[0],dp=2)
		widget_control,u.id.m2max,set_value=num2str(u.stat.plot.m2[1],dp=2)
	ENDELSE

	!except=0
	widget_control,base,/realize
	IF keyword_set(loaddata) THEN BEGIN
		load_profile_data,u
		IF u.stat.dat THEN BEGIN
			update_temporal_text,u
			plot_profiles,u
			plot_moments,u
			tmp=where(u.dat.tau GT 0)
			widget_control,u.id.ntpt,set_value=num2str(n(tmp),1)
		ENDIF
	END
	wprof_write_linelist,u
	xmanager,'w_hirexsr_profiles',base

END
