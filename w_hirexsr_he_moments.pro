FUNCTION is_dat,u
	taglist=tag_names(u)
	tmp=where(taglist EQ 'DAT')
	IF tmp[0] NE -1 THEN output=1 ELSE output=0
	RETURN,output
END

PRO save_bins,u
	hirexsr_write_binning,u.shot,chmap=*u.dat.chmap,tch=*u.dat.tch,tmap=u.dat.tmap,good=u.dat.good,tht=u.tht
	widget_control,u.id.message,set_value='CHMAP, TCH and TMAP saved: '+num2str(u.shot,1)+' THT-'+num2str(u.tht,1),/app
END

PRO save_mom,u
	line=u.stat.line
	pos=*u.dat.pos
	etendue=*u.dat.etendue
	hirexsr_write_moments,u.shot,line,u.dat.mom[*,*,0:2,line],u.dat.mom[*,*,3:5,line],u.dat.mom[*,*,6:8,line],u.dat.mom[*,*,9:11,line],$
		u.dat.tau,pos[*,*,line,*],u.dat.mom[*,*,12,line],u.dat.mom[*,*,13,line],u.dat.mom[*,*,14,line],etendue[*,line,*],$
		*u.dat.tch,u.dat.dlam[line],u.stat.double,u.dat.mom[*,*,15,line],u.dat.tree,tht=u.tht
	widget_control,u.id.message,set_value='LINE = '+num2str(u.stat.line,1)+' MOM Saved - '+num2str(u.shot,1)+' THT-'+num2str(u.tht,1),/app
END

PRO save_fits,u
	line=u.stat.line
	coefs_arr=hirexsr_coefs_ptr2arr(u.dat.coefs[*,*,line],u.stat.chmax,u.stat.ntime)
	hirexsr_write_fits,u.shot,line,coefs_arr,u.dat.tau,u.dat.nave,u.stat.double,*u.dat.labels[line],tht=u.tht
	widget_control,u.id.message,set_value='LINE = '+num2str(u.stat.line,1)+' COEFS Saved - '+num2str(u.shot,1)+' THT-'+num2str(u.tht,1),/app
END

PRO save_spec,u
	avespec_arr=hirexsr_avespec_ptr2arr(u.dat.avespec,u.stat.chmax,u.stat.ntime,u.stat.maxave)
	specbr=avespec_arr[*,*,*,0]
	lam=avespec_arr[*,*,*,1]
	sig=avespec_arr[*,*,*,2]
	resid=avespec_arr[*,*,*,3]
	hirexsr_write_avespec,u.shot,specbr,lam,sig,resid,u.dat.tau,u.dat.nave,tht=u.tht
	widget_control,u.id.message,set_value='AVESPEC Saved - '+num2str(u.shot,1)+' THT-'+num2str(u.tht,1),/app
END

PRO plot_vessel,u
	;IF u.stat.ves EQ 1 THEN plot_detailed_vessel,u ELSE plot_rough_vessel,u
	plot_rough_vessel,u	
END

PRO plot_rough_vessel,u
	IF NOT u.stat.ps THEN BEGIN
		window,0,xsize=u.plot.xsize[0],ysize=u.plot.ysize[0],/pixmap
	ENDIF
	yrange=[-0.54,0.54]
	xrange=[0.40,0.98]
	plot, [0],[0],title=title,chars=1.0,xrange=xrange,yrange=yrange,xtit='R (m)',ytit='Z (m)',/xsty, /ysty
	oplot,u.dat.fs.rw,u.dat.fs.zw
END

PRO plot_detailed_vessel,u
	x_offset=0.0
	y_offset=0.0
	widget_control,u.id.draw1,get_value=draw_win
	wset,draw_win
	yrange=[-0.67,0.67]
	xrange=[0.38,1.1]
	IF NOT u.stat.ps THEN BEGIN
		window,0,xsize=u.plot.xsize[0],ysize=u.plot.ysize[0],/pixmap
	ENDIF
	plot, [0],[0],title=title,chars=1.0,xrange=xrange,yrange=yrange,xtit='R (m)',ytit='Z (m)',/xsty, /ysty
	shot=u.shot
	IF shot GT 1070101000 THEN BEGIN
		restore, "/home/labombard/minicad/vv_tiles_cryo_2007_s.vctr"
		for i=0,nvctr-1 do oplot,xvctr(0:lvctr(i)-1,i)+x_offset,yvctr(0:lvctr(i)-1,i)+y_offset,color=color,thick=thick
	ENDIF		

	IF shot GT 1020101000 AND shot LT 1070101000 THEN BEGIN
		restore, "/home/labombard/minicad/vv_and_tiles_2002_s.vctr"
		for i=0,nvctr-1 do oplot,xvctr(0:lvctr(i)-1,i)+x_offset,yvctr(0:lvctr(i)-1,i)+y_offset,color=color,thick=thick
	ENDIF
	
	IF shot LT 1020101000 THEN BEGIN
		restore, '/home/labombard/minicad/vacuumvessel.vctr'
		for i=0,nvctr-1 do oplot,xvctr(0:lvctr(i)-1,i)+x_offset,yvctr(0:lvctr(i)-1,i)+y_offset,color=color,thick=thick
		restore, '/home/labombard/minicad/tiles_2002_s.vctr'
		for i=0,nvctr-1 do oplot,xvctr(0:lvctr(i)-1,i)+x_offset,yvctr(0:lvctr(i)-1,i)+y_offset,color=color,thick=thick
	ENDIF
END

PRO plot_fs,u
;	widget_control,u.id.draw1,get_value=draw_win
;	wset,draw_win
	index=ipt(u.dat.fs.t,u.time)
	IF index[0] EQ -1 THEN RETURN
	z_xpt=u.dat.fs.zxpt[index]
	rho=[make(min(u.dat.fs.rhofine),-.01,8),-0.004,0.0,0.004,0.008]
        FOR i=0,n(rho) DO BEGIN
        	contour,u.dat.fs.rhofine[*,*,index],u.dat.fs.rfine,u.dat.fs.zfine,levels=[rho[i]],/overplot,path_xy=path_fs,path_info=info_fs,/path_data
        	rfs=reform(path_fs[0,*])
        	zfs=reform(path_fs[1,*])
                good=inside(rfs,zfs,u.dat.fs.rw,u.dat.fs.zw)
        	;IF rho[i] LT 0 THEN BEGIN
                ;	tmp=where(good EQ 1 AND zfs GT u.dat.fs.zxpt[index]) 
                ;        color=80
                ;ENDIF ELSE BEGIN
                ;        tmp=where(good EQ 1)
                ;        IF rho[i] EQ 0 THEN color=200 ELSE color=50
                ;ENDELSE
        	;IF tmp[0] NE -1 THEN oplot,rfs[tmp],zfs[tmp],color=color
                ;IF rho[i] LT 0 THEN BEGIN
              	;	tmp=where(good EQ 1 AND zfs LT u.dat.fs.zxpt[index]) 
                ;        color=50
                ;        IF tmp[0] NE -1 THEN oplot,rfs[tmp],zfs[tmp],color=color
                ;ENDIF
                IF z_xpt LT 0 THEN BEGIN
                	IF rho[i] LT 0 THEN BEGIN
                        	tmp=where(good EQ 1 AND zfs GT z_xpt) 
	                        color=100
        	        ENDIF ELSE BEGIN
                	        tmp=where(good EQ 1)
                        	IF rho[i] EQ 0 THEN color=200 ELSE color=30
	                ENDELSE
        		IF tmp[0] NE -1 THEN oplot,rfs[tmp],zfs[tmp],color=color
                	IF rho[i] LT 0 THEN BEGIN
              			tmp=where(good EQ 1 AND zfs LT z_xpt) 
                                color=30
	                        IF tmp[0] NE -1 THEN oplot,rfs[tmp],zfs[tmp],color=color
                        ENDIF
                ENDIF ELSE BEGIN
                	IF rho[i] LT 0 THEN BEGIN
                        	tmp=where(good EQ 1 AND zfs LT z_xpt) 
	                        color=100
        	        ENDIF ELSE BEGIN
                	        tmp=where(good EQ 1 AND zfs GT -0.4)
                        	IF rho[i] EQ 0 THEN color=200 ELSE color=30
	                ENDELSE
        		IF tmp[0] NE -1 THEN oplot,rfs[tmp],zfs[tmp],color=color
                	IF rho[i] LT 0 THEN BEGIN
              			tmp=where(good EQ 1 AND zfs GT z_xpt)
	                        IF tmp[0] NE -1 THEN oplot,rfs[tmp],zfs[tmp],color=color
                        ENDIF
                
                ENDELSE
        ENDFOR        
END

PRO plot_pos,u
	ns=100
;	widget_control,u.id.draw1,get_value=draw_win
;	wset,draw_win
	pos=*u.dat.pos
	pos=reform(pos[*,*,u.stat.line,*])
	tpos=*u.dat.tch
	time=u.dat.tau[u.index]
	IF n(tpos) NE 0 THEN BEGIN
		tmp=where(tpos LE time)
		pos=pos[*,*,last(tmp)]
	ENDIF
	
	IF u.stat.bins THEN BEGIN
		start=0
		tmp=where(pos[2,*] NE -1)
		stop=last(tmp)
	ENDIF ELSE BEGIN
		start=u.ch-1
		stop=u.ch-1
	ENDELSE
	FOR i=start,stop DO BEGIN
		IF i EQ u.ch-1 THEN col=200 ELSE col=255
		ipos=pos[*,i]
		starts=line_s(ipos,r=1.0)	;where it enters the plot range
		starts=min(starts)		;take first intersection
		ends=line_s(ipos,r=0.44)	;where it hits the inner wall
		ends=min(ends)			;take first intersection
		splot=make(starts,ends,ns)
		rplot=line_r(splot,ipos)
		zplot=line_z(splot,ipos)
		oplot,rplot,zplot,color=col
	ENDFOR		
END

PRO plot_cx,u
	widget_control,u.id.draw1,get_value=draw_win
	plot_vessel,u
	IF u.stat.fs THEN plot_fs,u
	IF u.stat.pos THEN plot_pos,u
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize[0],u.plot.ysize[0],0,0,0]
	ENDIF
END


PRO plot_all_moments,u
	plot_spec_moments,u
	IF NOT u.stat.zoom THEN plot_time_moments,u
	IF u.stat.mconv THEN plot_conv_moments,u ELSE plot_spatial_moments,u
END

PRO plot_time_moments,u
	c=2.998e8 			;speed of light
	e=1.602e-19			;conversion for eV -> J
	mconv=1.661e-27			;conversion for amu -> kg
	CASE u.stat.line OF
		4 : mass = 95.94			;mass of molybdenum in amu
		ELSE : mass=39.9477			;mass of argon in amu
	ENDCASE
	line=u.stat.line
	lam_o=u.wl.lam_o[line]
	conv_factor=(lam_o/c)^2*(e*1.0e3/(mass*mconv))		;conversion factor for 
	makesym,9

	ch=u.ch-1
	index=u.index
	ntau=n(where(u.dat.tau GE 0))+1
	ifitc=u.dat.mom[ch,0:ntau-1,15,line]

	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw5,get_value=draw_win
		window,0,xsize=u.plot.xsize[4],ysize=u.plot.ysize[4],/pixmap
	ENDIF
	jgood=u.dat.good[ch,*]
	tmpg=where(jgood EQ 1)
	xr=[0,2.0]
	
	;zeroth moment
	IF u.stat.m0 NE 0 AND tmpg[0] NE -1 THEN BEGIN
		IF u.stat.mplot THEN BEGIN
			jmom=reform(u.dat.mom[ch,0:ntau-1,0,line])
			jerr=reform(u.dat.mom[ch,0:ntau-1,3,line])
			yr=[0,max(jmom+jerr)*1.05]
			plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='0!uth!n Moment',/xsty,/ysty
		ENDIF ELSE BEGIN
			jmom=reform(u.dat.mom[ch,0:ntau-1,6,line])
			jerr=reform(u.dat.mom[ch,0:ntau-1,9,line])
			yr=[0,max(jmom+jerr)*1.05]
			plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='# of Photons',/xsty,/ysty		
		ENDELSE
		oplot,u.dat.tau[0:ntau-1],jmom[0:ntau-1]
		tmp=where(ifitc EQ 0)
		IF tmp[0] NE -1 THEN oploterror,u.dat.tau[tmp],jmom[tmp],jerr[tmp],psym=8
		tmp=where(ifitc EQ 3)
		IF tmp[0] NE -1 THEN oploterror,u.dat.tau[tmp],jmom[tmp],jerr[tmp],psym=7
		oplot,u.dat.tau[index]*[1,1],yr,linestyle=2.0,color=60
	ENDIF
	
	;first moment
	IF u.stat.m1 NE 0 AND tmpg[0] NE -1 THEN BEGIN
		IF u.stat.mplot THEN BEGIN			;plot moment
			jmom=reform(u.dat.mom[ch,0:ntau-1,1,line])
			jerr=reform(u.dat.mom[ch,0:ntau-1,4,line])
			IF min(jmom-jerr) LT 0 THEN ymin=min(jmom-jerr)*1.05 ELSE ymin=min(jmom-jerr)*0.95
 			IF max(jmom+jerr) GT 0 THEN ymax=max(jmom+jerr)*1.05 ELSE ymax=min(jmom+jerr)*0.95
			yr=[ymin,ymax]
			plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='1!ust!n Moment',/xsty,/ysty
		ENDIF ELSE BEGIN				;plot profile moment
			mu=reform(u.dat.mom[ch,0:ntau-1,7,line])
			sig_mu=reform(u.dat.mom[ch,0:ntau-1,10,line])
			jmom=(mu-lam_o)*c/lam_o*1.0e-3		;velocity
			jerr=sig_mu*c/lam_o*1.0e-3		;velocity error
			IF min(jmom-jerr) LT 0 THEN ymin=min(jmom-jerr)*1.05 ELSE ymin=min(jmom-jerr)*0.95
 			IF max(jmom+jerr) GT 0 THEN ymax=max(jmom+jerr)*1.05 ELSE ymax=min(jmom+jerr)*0.95
			yr=[ymin,ymax]
			plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='Doppler Shift [km/s]',/xsty,/ysty
		ENDELSE
		oplot,u.dat.tau[0:ntau-1],jmom
		tmp=where(ifitc EQ 0)
		IF tmp[0] NE -1 THEN oploterror,u.dat.tau[tmp],jmom[tmp],jerr[tmp],psym=8
		tmp=where(ifitc EQ 3)
		IF tmp[0] NE -1 THEN oploterror,u.dat.tau[tmp],jmom[tmp],jerr[tmp],psym=7
		oplot,u.dat.tau[index]*[1,1],yr,linestyle=2.0,color=60
		oplot,xr,[0,0],linestyle=2.0

	ENDIF

	;second moment
	IF u.stat.m2 NE 0 AND tmpg[0] NE -1 THEN BEGIN
		IF u.stat.mplot THEN BEGIN
			jmom=reform(u.dat.mom[ch,0:ntau-1,2,line])
			jerr=reform(u.dat.mom[ch,0:ntau-1,5,line])
			yr=[0,max(jmom+jerr)*1.05]
			plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='2!und!n Moment',/xsty,/ysty
		ENDIF ELSE BEGIN
			jmom=(reform(u.dat.mom[ch,0:ntau-1,8,line]))^2			;w^2
			jerr=reform(u.dat.mom[ch,0:ntau-1,11,line])*2.0*sqrt(jmom)	;uncertainty prop through the square sig_ti=2*a*w*sig_w
			jmom*=mass*mconv*c^2/(e*1.0e3*lam_o^2)
			jerr*=mass*mconv*c^2/(e*1.0e3*lam_o^2)
			yr=[0,max(jmom+jerr)*1.05]
			plot,[0],[0],xr=xr,yr=yr,xtit='Time [sec]',ytit='Ion Temp [keV]',/xsty,/ysty
		ENDELSE
		oplot,u.dat.tau[0:ntau-1],jmom
		tmp=where(ifitc EQ 0)
		IF tmp[0] NE -1 THEN oploterror,u.dat.tau[tmp],jmom[tmp],jerr[tmp],psym=8
		tmp=where(ifitc EQ 3)
		IF tmp[0] NE -1 THEN oploterror,u.dat.tau[tmp],jmom[tmp],jerr[tmp],psym=7
		oplot,u.dat.tau[index]*[1,1],yr,linestyle=2.0,color=60
	ENDIF
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize[4],u.plot.ysize[4],0,0,0]
	ENDIF
END

PRO plot_spatial_moments,u
	c=2.998e8 			;speed of light
	e=1.602e-19			;conversion for eV -> J
	mconv=1.661e-27			;conversion for amu -> kg
	CASE u.stat.line OF
		4 : mass = 95.94			;mass of molybdenum in amu
		ELSE : mass=39.9477			;mass of argon in amu
	ENDCASE
	line=u.stat.line
	lam_o=u.wl.lam_o[line]
	conv_factor=(lam_o/c)^2*(e*1.0e3/(mass*mconv))		;conversion factor for 
	
	ch=u.ch-1
	index=u.index
	ifitc=u.dat.mom[*,index,15,line]

	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw6,get_value=draw_win
		window,0,xsize=u.plot.xsize[5],ysize=u.plot.ysize[5],/pixmap
	ENDIF
	chars=2.4
	!p.multi=[0,0,3]
	IF u.stat.mplot THEN BEGIN
		imom=reform(u.dat.mom[*,index,0:2,line])
		imom[*,1]*=1.0e3
		imom[*,2]*=1.0e6
		ierr=reform(u.dat.mom[*,index,3:5,line])
		ierr[*,1]*=1.0e3
		ierr[*,2]*=1.0e6
		tit0='0!uth!n Moment'
		tit1='10!u-3!n 1!ust!n Moment'
		tit2='10!u-6!n 2!und!n Moment'
	ENDIF ELSE BEGIN
		imom=reform(u.dat.mom[*,index,6:8,line])
		imom[*,1]=(imom[*,1]-lam_o)*c/lam_o*1.0e-3
		w=imom[*,2]
		imom[*,2]=w^2*mass*mconv*c^2/(e*1.0e3*lam_o^2)
		ierr=reform(u.dat.mom[*,index,9:11,line])
		ierr[*,1]=ierr[*,1]*c/lam_o*1.0e-3
		ierr[*,2]=ierr[*,2]*2.0*w*mass*mconv*c^2/(e*1.0e3*lam_o^2)
		tit0='# of Photons'
		tit1='Doppler Shift [km/s]'
		tit2='Ion Temp [keV]'
	ENDELSE
	igood=u.dat.good[*,index]
	
	
	tmp=where(igood EQ 1)
	chplot=tmp+1
	maxpt=max(imom[tmp,0]+ierr[tmp,0])*1.05
	yr=[0,maxpt]
	plot,[0],[0],xr=[0,max(tmp)+1],yr=yr,ytit=tit0,chars=chars,/xsty,/ysty
	makesym,9
	tmpc=where(ifitc EQ 0)
	IF tmpc[0] NE -1 THEN oploterror,chplot[tmpc],imom[tmpc,0],ierr[tmpc,0],psym=8
	tmpc=where(ifitc EQ 3)
	IF tmpc[0] NE -1 THEN oploterror,chplot[tmpc],imom[tmpc,0],ierr[tmpc,0],psym=7
	makesym,10
	;oploterror,ch,imom[ch,0],ierr[ch,0],psym=8,color=150,errcolor=150
	oplot,[ch,ch]+1,yr,linestyle=2.0,color=150

	tmp=where(igood EQ 1 AND finite(imom[*,1]) EQ 1 AND finite(ierr[*,1]) EQ 1)
	chplot=tmp+1
	IF u.stat.mplot THEN BEGIN
		minpt=min(imom[tmp,1]-ierr[tmp,1])
		maxpt=max(imom[tmp,1]+ierr[tmp,1])
	ENDIF ELSE BEGIN
		maxpt=u.stat.vplt[1]
		minpt=u.stat.vplt[0]
	ENDELSE
	IF minpt LT 0 THEN minpt*=1.05 ELSE minpt*=0.95
	yr=[minpt,maxpt]
	xr=[0,max(tmp)+1]
	plot,[0],[0],xr=xr,yr=yr,ytit=tit1,chars=chars,/xsty,/ysty
	makesym,9
	tmpc=where(ifitc EQ 0)
	IF tmpc[0] NE -1 THEN oploterror,chplot[tmpc],imom[tmpc,1],ierr[tmpc,1],psym=8
	tmpc=where(ifitc EQ 3)
	IF tmpc[0] NE -1 THEN oploterror,chplot[tmpc],imom[tmpc,1],ierr[tmpc,1],psym=7
	makesym,10
	;IF finite(imom[ch,1]) AND finite(ierr[ch,1]) THEN oploterror,ch,imom[ch,1],ierr[ch,1],psym=8,color=150,errcolor=150
	oplot,[ch,ch]+1,yr,linestyle=2.0,color=150
	oplot,xr,[0,0],linestyle=2.0

	tmp=where(igood EQ 1 AND finite(imom[*,2]) EQ 1 AND finite(ierr[*,2]) EQ 1)
	chplot=tmp+1
	IF u.stat.mplot THEN BEGIN
		minpt=0
		maxpt=max(imom[tmp,2]+ierr[tmp,2])
	ENDIF ELSE BEGIN
		maxpt=u.stat.tplt[1]
		minpt=u.stat.tplt[0]
	ENDELSE
	yr=[minpt,maxpt]
	plot,[0],[0],xr=[0,max(tmp)+1],yr=yr,ytit=tit2,chars=chars,/xsty,/ysty
	makesym,9
	tmpc=where(ifitc EQ 0)
	IF tmpc[0] NE -1 THEN oploterror,chplot[tmpc],imom[tmpc,2],ierr[tmpc,2],psym=8
	tmpc=where(ifitc EQ 3)
	IF tmpc[0] NE -1 THEN oploterror,chplot[tmpc],imom[tmpc,2],ierr[tmpc,2],psym=7
	makesym,10
	;IF finite(imom[ch,2]) AND finite(ierr[ch,2]) THEN oploterror,ch,imom[ch,2],ierr[ch,2],psym=8,color=150,errcolor=150
	oplot,[ch,ch]+1,yr,linestyle=2.0,color=150

	!p.multi=0
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize[5],u.plot.ysize[5],0,0,0]
	ENDIF
END

PRO plot_conv_moments,u
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw6,get_value=draw_win
		window,0,xsize=u.plot.xsize[5],ysize=u.plot.ysize[5],/pixmap
	ENDIF
	line=u.stat.line
	lam_o=u.wl.lam_o[line]
	ch=u.ch-1
	index=u.index
	IF u.stat.ave THEN ispec=*u.dat.avespec[ch,index] ELSE ispec=*u.dat.spec[ch,index]
	icoefs=*u.dat.coefs[ch,index,u.stat.line]
	ilabels=*u.dat.labels[u.stat.line]
	dlam=make(0.4,6.0,25)
	convmom=hirexsr_moment_conv(ispec,icoefs,ilabels,lam_o,line,dlam)
	!p.multi=[0,0,3]
	widget_control,u.id.dlampt,get_value=dlamstr
	dlampt=last(float(dlamstr))
	xr=[0,max(dlam)]

	;0th moment
	maxpt=max(convmom[*,0,0]+convmom[*,0,1])
	yr=[0,maxpt*1.05]
	plot,dlam,convmom[*,0,0],xtit='d'+n2g('lambda')+' [mAng]',ytit='0!uth!n Moment',chars=2.4,xr=xr,/xsty,yr=yr,/ysty
	oploterror,dlam,convmom[*,0,0],convmom[*,0,1]
	loadct,39,/silent
	oplot,dlampt*[1,1],yr,color=190,linestyle=3.0
	m0pt=interpol(convmom[*,0,0],dlam,dlampt)
	oplot,xr,m0pt*[1,1],color=190,linestyle=3.0	
	loadct,12,/silent

	;1st moment
	maxpt=max(convmom[*,1,0]+convmom[*,1,1])
	minpt=min(convmom[*,1,0]-convmom[*,1,1])
	IF maxpt GT 0 THEN maxpt*=1.05 ELSE maxpt*=0.95
	IF minpt LT 0 THEN minpt*=1.05 ELSE minpt*=0.95
	yr=[minpt,maxpt]
	plot,dlam,convmom[*,1,0],xtit='d'+n2g('lambda')+' [mAng]',ytit='1!ust!n Moment',chars=2.4,xr=xr,/xsty,yr=yr,/ysty
	oploterror,dlam,convmom[*,1,0],convmom[*,1,1]
	loadct,39,/silent
	oplot,dlampt*[1.0,1.0],yr,color=190,linestyle=3.0
	m1pt=interpol(convmom[*,1,0],dlam,dlampt)
	oplot,xr,m1pt*[1,1],color=190,linestyle=3.0	
	loadct,12,/silent
	oplot,xr,[0,0],linestyle=2.0

	;2nd moment
	maxpt=max(convmom[*,2,0]+convmom[*,2,1])
	yr=[0,maxpt*1.05]
	plot,dlam,convmom[*,2,0],xtit='d'+n2g('lambda')+' [mAng]',ytit='2!und!n Moment',chars=2.4,xr=xr,/xsty,yr=yr,/ysty
	oploterror,dlam,convmom[*,2,0],convmom[*,2,1]
	loadct,39,/silent
	oplot,dlampt*[1,1],yr,color=190,linestyle=3.0
	m2pt=interpol(convmom[*,2,0],dlam,dlampt)
	oplot,xr,m2pt*[1,1],color=190,linestyle=3.0	
	loadct,12,/silent

	!p.multi=0
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize[5],u.plot.ysize[5],0,0,0]
	ENDIF
END

PRO plot_spec_moments,u	 
	;plot background spectra
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw4,get_value=draw_win
		window,0,xsize=u.plot.xsize[3],ysize=u.plot.ysize[3],/pixmap
	ENDIF
	line=u.stat.line
	CASE line OF
		0 : BEGIN
			wlplot=[3.943,3.955]
		END
		1 : BEGIN
			wlplot=[3.962,3.972]
		END
		2 : BEGIN
			wlplot=[3.989,4.0]
                END
		3 : BEGIN
			wlplot=[3.723,3.739]
                END
		4 : BEGIN
			wlplot=[3.735,3.744]
                END
	ENDCASE
	lam_o=u.wl.lam_o[line]
	
	ch=u.ch-1
	index=u.index
	IF u.stat.ave THEN ispec=*u.dat.avespec[ch,index] ELSE ispec=*u.dat.spec[ch,index]
	lam=ispec[*,1]
	cnts=ispec[*,0]
	sig=ispec[*,2]
	icoefs=*u.dat.coefs[ch,index,u.stat.line]
	ilabels=*u.dat.labels[u.stat.line]
	cnts=hirexsr_line_subtract(lam,cnts,icoefs,ilabels,line)
	tmp=where(lam GE wlplot[0] AND lam LE wlplot[1])
	maxpt=max(cnts[tmp]+sig[tmp])*1.05
	minpt=(min(cnts[tmp]-sig[tmp]) < 0)*1.05
	plot,[0],[0],xr=wlplot,yr=[minpt,maxpt],ytit='B!l'+n2g('lambda')+'!n [AU]',/xsty,/ysty
	oploterror,lam[tmp],cnts[tmp],sig[tmp],psym=8,symsize=0.25
	oplot,wlplot,[0,0],linestyle=2.0
	oplot,lam_o*[1,1],[minpt,maxpt],color=40,linestyle=3.0
	idlam=u.dat.dlam[line]
	loadct,39,/silent
	oplot, lam_o*[1,1]-idlam/1.0e3,[minpt,maxpt],linestyle=3.0,color=190
	oplot, lam_o*[1,1]+idlam/1.0e3,[minpt,maxpt],linestyle=3.0,color=190
	loadct,12,/silent
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize[3],u.plot.ysize[3],0,0,0]
	ENDIF

END

PRO plot_spec,u,debug=debug
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw3,get_value=draw_win
		window,0,xsize=u.plot.xsize[2],ysize=u.plot.ysize[2],/pixmap
	ENDIF
	wl=u.wl.(u.stat.line)
	ch=u.ch-1
	index=u.index
	IF u.stat.ave THEN ispec=*u.dat.avespec[ch,index] ELSE ispec=*u.dat.spec[ch,index]
	lam=ispec[*,1]
	cnts=ispec[*,0]
	sig=ispec[*,2]
	resid=ispec[*,3]
	icoefs=*u.dat.coefs[ch,index,u.stat.line]
	tmp=where(lam GE wl[0] AND lam LE wl[1])
	lamfit=make(wl[0],wl[1],100)
	maxpt=max(cnts[tmp])*1.05
	pos=[0.05,0.3,0.95,0.98]
	plot,[0],[0],xr=wl,yr=[0,maxpt],ytit='B!l'+n2g('lambda')+'!n [AU]',/xsty,/ysty,pos=pos
	oploterror,lam[tmp],cnts[tmp],sig[tmp],psym=8,symsize=0.25
	x=size(icoefs)
	IF x[1] NE 1 THEN BEGIN
		ncoefs=x[1]
		CASE ncoefs MOD 3 OF 
			0 : BEGIN
				nline=ncoefs/3-1
				nb=2
			END
			1 : BEGIN
				nline=ncoefs/3
				nb=0
			END
			2 : BEGIN
				nline=ncoefs/3
				nb=1
                        END
		ENDCASE
		FOR i=0,nline-1 DO oplot,lamfit,gaussian_fits(lamfit,icoefs[i*3:(i+1)*3-1]),color=80
		oplot,lamfit,gaussian_fits(lamfit,icoefs,base=base),color=200,thick=2
		IF n(base) NE 0 THEN oplot,lamfit,base,color=60 ELSE oplot,[min(lamfit),max(lamfit)],base*[1.,1],color=60
		;CASE nb OF 
		;	0 : oplot,[lamfit[0],last(lamfit)],icoefs[ncoefs-1]*[1.0,1.0],color=60
		;	1 : oplot,lamfit,lamfit*icoefs[ncoefs-2]+icoefs[ncoefs-1],color=60
		;	2 : oplot,lamfit,lamfit^2*icoefs[ncoefs-3]+lamfit*icoefs[ncoefs-2]+icoefs[ncoefs-1],color=60
		;ENDCASE
		oplot,lamfit,gaussian_fits(lamfit,icoefs),color=200,thick=2
        ENDIF
	IF keyword_set(debug) THEN stop
	pos=[0.05,0.05,0.95,0.3]
	maxpt=max(abs(resid[tmp]+sig[tmp]))
	plot,[0],[0],xr=wl,yr=[-1.0*maxpt,maxpt],xtit='Wavelength [Ang]',ytit='RESID',/xsty,/ysty,pos=pos,/noerase
	oploterror,lam[tmp],resid[tmp],sig[tmp],psym=8,symsize=0.25,color=200,errcolor=200
	oplot,wl,[0,0],linestyle=2.0
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize[2],u.plot.ysize[2],0,0,0]
	ENDIF
END

PRO plot_hemom_image,u
	IF NOT u.stat.image THEN ff=0.0 ELSE ff=1.0
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw2,get_value=draw_win
		window,0,xsize=u.plot.xsize[1],ysize=u.plot.ysize[1],/pixmap
	ENDIF
	image=rotate(*u.dat.cnts[u.index],4)
	x=size(image)
	widget_control,u.id.gain_slider,get_value=gain
	gain*=ff
	nx=x[1]
	ny=x[2]
	zoom=0.7

	;mask killed pixels
	killval=0.0
	chmap=*u.dat.chmap
	tch=*u.dat.tch
	imap=rotate(chmap[*,*,last(where(tch LE u.dat.tau[u.index]))],4)
	tmp=where(imap EQ -2)
	IF tmp[0] NE -1 THEN image[tmp]=killval	;set killed pixels to neighbor
	
	i_pt=indgen(nx*zoom)/zoom
	j_pt=indgen(ny*zoom)/zoom
	new_pic=interpolate(image,i_pt,j_pt,/grid)
	loadct,39,/silent
	norm=max(new_pic)
	tv,new_pic/norm*256.0*gain
	loadct,12,/silent
	IF NOT u.stat.image THEN RETURN

	bnds=*u.dat.bnds
	tmp=where(*u.dat.tch LE u.time)
	r0=reform(bnds[0,*,last(tmp)])
	r1=reform(bnds[1,*,last(tmp)])
	ch=u.ch-1
	
	IF u.stat.bins THEN BEGIN
		tmp=where(r0 NE -1)
		FOR i=0,n(tmp) DO BEGIN
			IF i EQ ch THEN color=200  ELSE color=240
			plots, [0,1,1,0,0],[r0[i],r0[i],r1[i],r1[i],r0[i]]/float(ny),color=color,/norm
		ENDFOR
	ENDIF ELSE plots, [0,1,1,0,0],[r0[ch],r0[ch],r1[ch],r1[ch],r0[ch]]/float(ny),color=200,/norm
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize[1],u.plot.ysize[1],0,0,0]
	ENDIF		

	IF u.stat.zoom THEN BEGIN
		zoom=3.076
		height=int(250/zoom)
		ex=int((height-(r1[ch]-r0[ch]))/2.0)
		cent=int((r0[ch]+r1[ch])/2.0)
                IF r0[ch]-ex LT 0 THEN BEGIN
                	lowbnd=0 
	                zoomx = zoom
        	        zoomy = zoom*(2.*ex+r1[ch]-r0[ch])/((2.*ex+r1[ch]-r0[ch])+2.*(r0[ch]-ex)) ;expand the y channel until it fills the screen
                    	height = int(250/zoomy)	
                    	ex=int((height-(r1[ch]-r0[ch]))/2.0)
                    	IF r1[ch]+ex GT ny THEN highbnd=ny-1 ELSE highbnd=r1[ch]+ex
                ENDIF ELSE BEGIN 
                    	lowbnd=r0[ch]-ex
                    	zoomx = zoom
                    	zoomy = zoom 
			IF r1[ch]+ex GT ny THEN BEGIN
                        	highbnd = ny-1
                        	zoomx = zoom
                        	zoomy = zoom * (2.*ex+r1[ch]-r0[ch])/((2.*ex+r1[ch]-r0[ch]) -2.*(r1[ch]+ex-ny))
                        	height = int(250/zoomy)
                        	ex=int((height-(r1[ch]-r0[ch]))/2.0)
                        	lowbnd = r0[ch]-ex; overwrite with new ex
                    	ENDIF ELSE BEGIN
                        	highbnd = r1[ch]+ex
                    	ENDELSE
                ENDELSE
		sub=image[*,lowbnd:highbnd]
		submap=imap[*,lowbnd:highbnd]
		tmp=where(submap EQ -2)
		IF tmp[0] NE -1 THEN sub[tmp]=killval	 ;set killed pixels to zero
		x=size(sub)
		nx=x[1]
		ny=x[2]
		i_pt=indgen(nx*zoomx)/zoomx
		j_pt=indgen(ny*zoomy)/zoomy
		new_pic=interpolate(sub,i_pt,j_pt,/grid)
		IF NOT u.stat.ps THEN BEGIN
			widget_control,u.id.draw5,get_value=draw_win
			window,0,xsize=u.plot.xsize[4],ysize=u.plot.ysize[4],/pixmap
		ENDIF
		loadct,39,/silent
		tv,new_pic/norm*256.0*gain
		loadct,12,/silent
		color=200
		frac=(r1[ch]-r0[ch])/float(height)/2.0
		plots, [0,1,1,0,0],0.5+frac*[1.0,1.0,-1.0,-1.0,1.0],color=color,/norm
		IF NOT u.stat.ps THEN BEGIN
			wset,draw_win
			device,copy=[0,0,u.plot.xsize[4],u.plot.ysize[4],0,0,0]
		ENDIF	
	ENDIF
END



PRO update_hemom_temporal_text,u
	widget_control,u.id.itpt,set_value=num2str(u.index,1)
	widget_control,u.id.ntpt,set_value=num2str(max(u.dat.tmap),1)
	tmp=where(u.dat.tmap EQ u.index)
	widget_control,u.id.flow,set_value=num2str(tmp[0],1)
	widget_control,u.id.fhigh,set_value=num2str(last(tmp),1)
	widget_control,u.id.tpt,set_value=num2str(u.time,dp=2)
	widget_control,u.id.t_slider,set_value=u.time*1.0e3
END

PRO update_spatial_text,u
	imap=*u.dat.chmap
	imap=imap[*,*,u.mapindex]
	ibnds=*u.dat.bnds
	ibnds=ibnds[*,*,u.mapindex]
	widget_control,u.id.ch,set_value=num2str(u.ch,1)
	widget_control,u.id.chtot,set_value=num2str(max(imap)+1,1)
	widget_control,u.id.rlow,set_value=num2str(ibnds[0,u.ch-1],1)
	widget_control,u.id.rhigh,set_value=num2str(ibnds[1,u.ch-1],1)
	widget_control,u.id.ch_slider,set_value=u.ch
END

PRO update_moment_text,u
	line=u.stat.line
	widget_control,u.id.dlampt,set_value=num2str(u.dat.dlam[line],1)
	widget_control,u.id.lowpt,set_value=num2str(u.wl.lam_o[line]-u.dat.dlam[line]/1.0e3,1)
	widget_control,u.id.highpt,set_value=num2str(u.wl.lam_o[line]+u.dat.dlam[line]/1.0e3,1)
END	

PRO load_hemom_data,u
	shot=u.shot
	tht=u.tht
	IF is_dat(u) THEN BEGIN					;if LOAD_DATA already ran, clean out the heap
		cleanheap_w_hirexsr_he_moments,u.id.base
		heap_gc
	ENDIF

	;load raw counts data
	hirexsr_ptr_image,shot,raw,lambda,t,const=const,morder=morder
	widget_control,u.id.message,set_value='Loaded Raw Data: '+num2str(shot,1),/app
	ntime=n(t)+1
	u.stat.ntime=ntime
	widget_control,u.id.nframes,set_value=num2str(u.stat.ntime,1)

	;load binning data
	hirexsr_load_binning,shot,chmap,tch,tmap,good,chmax,tht=tht
	u.stat.chmax=chmax
	widget_control,u.id.chmax,set_value=num2str(u.stat.chmax,1)

	;load INFO files
	info1=hirexsr_load_info(morder[0],shot=shot,/tree)
	info2=hirexsr_load_info(morder[1],shot=shot,/tree)
	info3=hirexsr_load_info(morder[2],shot=shot,/tree)
	info=[info1,info2,info3]
	widget_control,u.id.message,set_value='Loaded INFO Files: '+num2str(shot,1),/app

	nmaps=n(tch)+1
	u.stat.nmaps=nmaps
	chmap=ptr_new(chmap,/allocate)		;write champ, tch to PTRs so the size can be changed
	chmap_chk=ptr_new(*chmap,/allocate)	;make a backup copy of CHMAP to compare to when using ZOOMKILL
	tch=ptr_new(tch,/allocate)
	widget_control,u.id.message,set_value='Loaded Binning Setup: '+num2str(shot,1)+' THT-'+num2str(tht,1),/app

	;load fluxsurface data
	fs=make_fs_struc(shot)
	widget_control,u.id.message,set_value='Loaded FS Data: '+num2str(shot,1),/app

	;load fit coefficients
	coefs=ptrarr(u.stat.chmax,u.stat.ntime,5,/allocate_heap)
	labels=ptrarr(5,/allocate_heap)
	FOR i=0,n(u.stat.isline) DO BEGIN
		hirexsr_load_fits,shot,u.stat.linenum[i],icoefs,nave,double,ilabels,/quiet,status=status,tht=tht
		IF status THEN BEGIN
			coefsptr=hirexsr_coefs_arr2ptr(icoefs)
			coefs[*,*,i]=coefsptr
			*labels[i]=ilabels
			u.stat.isline[i]=1
		ENDIF ELSE u.stat.isline[i]=0
	ENDFOR
	widget_control,u.id.nave,set_value=num2str(nave,1)
	widget_control,u.id.message,set_value='Loaded Fit Coefs: '+num2str(shot,1)+' THT-'+num2str(tht,1),/app

	;load averaged spectra
	hirexsr_load_avespec,shot,specbr,lam,sig,resid,tht=tht
	avespec=hirexsr_avespec_arr2ptr(specbr,lam,sig,resid,maxave=maxave)
	u.stat.maxave=maxave
	widget_control,u.id.message,set_value='Loaded Ave Spec: '+num2str(shot,1)+' THT-'+num2str(tht,1),/app

	;load moments
	mom=fltarr(chmax,ntime,16,5)
	dlam=fltarr(5)
	pos=fltarr(4,u.stat.chmax,5,u.stat.nmaps)	
	etendue=fltarr(u.stat.chmax,5,u.stat.nmaps)
	FOR i=0,n(u.stat.isline) DO BEGIN
		IF u.stat.isline[i] THEN hirexsr_load_moments,shot,u.stat.linenum[i],imom,ierr,ipmom,iperr,tau,ipos,irhot,ifrac,iscale,iu,tpos,idlam,double,ifitcase,tree,tht=tht,/quiet,status=status ELSE status=0
		IF status THEN BEGIN
			tmp=[[[imom]],[[ierr]],[[ipmom]],[[iperr]],[[irhot]],[[ifrac]],[[iscale]],[[ifitcase]]]
			mom[*,*,*,i]=tmp
			dlam[i]=idlam
			IF n(tch)+1 EQ u.stat.nmaps THEN FOR j=0,u.stat.nmaps-1 DO BEGIN
				pos[*,*,i,j]=ipos
				etendue[*,i,j]=iu
			ENDFOR
		ENDIF
	ENDFOR
	pos=ptr_new(pos,/allocate)		;write POS to PTR to allow easy size change along with CHMAP
	etendue=ptr_new(etendue,/allocate)
	widget_control,u.id.message,set_value='Loaded Moments: '+num2str(shot,1)+' THT-'+num2str(tht,1),/app
	
	;apply tmap to raw images and time
	cnts=hirexsr_bin_image(raw,tmap)
	tau=hirexsr_bin_time(t,tmap)
	bnds=hirexsr_bin_bounds(*chmap,u.stat.chmax)
	bnds=ptr_new(bnds,/allocate)		;write to PTR so that size can be changed	
	u.index=ipt(tau,u.time)
	IF u.index EQ -1 THEN BEGIN
		u.time=max(tau)
		u.index=ipt(tau,u.time)
		widget_control,u.id.tpt,set_value=num2str(u.time,dp=2)
		widget_control,u.id.t_slider,set_value=u.time*1.0e3
	ENDIF
	widget_control,u.id.message,set_value='Raw Data Binned to Data Frames',/app
	
	;fill in temporal mapping text boxes
	widget_control,u.id.itpt,set_value=num2str(u.index,1)
	widget_control,u.id.ntpt,set_value=num2str(max(tmap),1)
	tmp=where(tmap EQ u.index)
	widget_control,u.id.flow,set_value=num2str(tmp[0],1)
	widget_control,u.id.fhigh,set_value=num2str(last(tmp),1)

	;fill in spatial mapping text boxes
	imap=*chmap
	imap=imap[*,*,u.mapindex]
	ibnds=*bnds
	ibnds=ibnds[*,*,u.mapindex]
	widget_control,u.id.ch,set_value=num2str(u.ch,1)
	widget_control,u.id.chtot,set_value=num2str(max(imap)+1,1)
	widget_control,u.id.rlow,set_value=num2str(ibnds[0,u.ch-1],1)
	widget_control,u.id.rhigh,set_value=num2str(ibnds[1,u.ch-1],1)

	;fill in moment text boxes
	line=u.stat.line
	widget_control,u.id.dlampt,set_value=num2str(dlam[line],1)
	widget_control,u.id.lowpt,set_value=num2str(u.wl.lam_o[line]-dlam[line]/1.0e3,1)
	widget_control,u.id.highpt,set_value=num2str(u.wl.lam_o[line]+dlam[line]/1.0e3,1)
	widget_control,u.id.lam_slider,set_value=dlam[line]*10.0

	IF u.stat.display EQ 0 THEN BEGIN		;form spectral bins
		spec=hirexsr_bin_spec(cnts,tau,*chmap,*tch,lambda,tmap,nchbins=u.stat.chmax,coefs=coefs,const=const)
		widget_control,u.id.message,set_value='Data Frames Arranged in Spectra',/app
	ENDIF ELSE spec=-1

	dat={fs:fs,raw:raw,t:t,info:info,pos:pos,etendue:etendue,lambda:lambda,chmap:chmap,chmap_chk:chmap_chk,tch:tch,bnds:bnds,good:good,tmap:tmap,$
		cnts:cnts,const:const,nave:nave,spec:spec,avespec:avespec,mom:mom,dlam:dlam,tau:tau,coefs:coefs,labels:labels,tree:tree}
	
	u={id:u.id,shot:u.shot,tht:u.tht,time:u.time,index:u.index,mapindex:u.mapindex,ch:u.ch,stat:u.stat,wl:u.wl,plot:u.plot,dat:dat}
	u.stat.dat=1
	widget_control,u.id.base, set_uvalue=u
	
END

PRO reset_spec_buttons,u
	widget_control,u.id.wn3mom,set_button=0
	widget_control,u.id.xymom,set_button=0
	widget_control,u.id.zjkmom,set_button=0
	widget_control,u.id.lyamom,set_button=0
	widget_control,u.id.momom,set_button=0
END

PRO cleanheap_w_hirexsr_he_moments,widget
	WIDGET_CONTROL, widget, get_uvalue=u
	taglist=tag_names(u)
	tmp=where(taglist EQ 'DAT')
	IF tmp[0] NE -1 THEN BEGIN
		heap_free,u.dat.raw
		heap_free,u.dat.cnts
		heap_free,u.dat.spec
		heap_free,u.dat.chmap
		heap_free,u.dat.tch
		heap_free,u.dat.good
		heap_free,u.dat.pos
		heap_free,u.dat.etendue
		heap_free,u.dat.const
		heap_free,u.dat.bnds
		heap_free,u.dat.avespec
		heap_free,u.dat.coefs
		heap_free,u.dat.labels
	ENDIF
END
	

PRO w_hirexsr_he_moments_event,event

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
					cleanheap_w_hirexsr_he_moments,event.top
					widget_control,event.top,/destroy
					heap_gc
					!except=1
				
				END
				"SAVE" :BEGIN
					IF u.stat.dat THEN BEGIN
						WIDGET_CONTROL,/hourglass
						line_in=u.stat.line
						save_spec,u
						save_bins,u
						FOR line=0,2 DO BEGIN
							IF u.stat.isline[line] THEN BEGIN
								u.stat.line=line
								save_fits,u
								save_mom,u
							ENDIF
						ENDFOR
						u.stat.line=line_in
					ENDIF
				END
				"LOAD": BEGIN
                   			WIDGET_CONTROL, /hourglass
					load_hemom_data,u
					plot_cx,u
					plot_hemom_image,u
					plot_spec,u
					plot_all_moments,u
				END
      				"PBIN": BEGIN
					IF event.select EQ 1 THEN u.stat.bins=1 ELSE u.stat.bins=0
					plot_hemom_image,u
					IF u.stat.pos THEN plot_cx,u
				END
				"PZOOM" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						u.stat.zoom=1 
						u.stat.m0=0.0
						u.stat.m1=0.0
						u.stat.m2=0.0
						widget_control,u.id.m0plot,set_button=u.stat.m0
						widget_control,u.id.m1plot,set_button=u.stat.m1
						widget_control,u.id.m2plot,set_button=u.stat.m2
					ENDIF ELSE u.stat.zoom=0
					plot_hemom_image,u
				END
				"PIMAGE" : BEGIN
					IF event.select EQ 1 THEN u.stat.image=1 ELSE u.stat.image=0
					plot_hemom_image,u
				END
				"PFS" : BEGIN
					IF event.select EQ 1 THEN u.stat.fs=1 ELSE u.stat.fs=0
					plot_cx,u
				END
				"PPOS" : BEGIN
					IF event.select EQ 1 THEN u.stat.pos=1 ELSE u.stat.pos=0
					plot_cx,u
				END
				"PAVE" : BEGIN
					IF event.select EQ 1 THEN u.stat.ave=1 ELSE u.stat.ave=0
					IF is_dat(u) THEN BEGIN
						type=size(u.dat.spec,/type)
						IF type EQ 10 THEN BEGIN
							plot_spec,u
							plot_all_moments,u
						ENDIF ELSE widget_control,u.id.message,set_value='NOT AVAILABLE IN DISPLAY MODE',/app
					ENDIF
				END
				"PCONV" : BEGIN
					IF event.select EQ 1 THEN u.stat.mconv=1 ELSE u.stat.mconv=0
					IF is_dat(u) THEN BEGIN
						plot_all_moments,u
					ENDIF
				END
				"M0PLOT" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						u.stat.zoom=0
						widget_control,u.id.pzoom,set_button=u.stat.zoom
						u.stat.m0=1.0
						u.stat.m1=0.0
						u.stat.m2=0.0
						widget_control,u.id.m0plot,set_button=u.stat.m0
						widget_control,u.id.m1plot,set_button=u.stat.m1
						widget_control,u.id.m2plot,set_button=u.stat.m2
						plot_time_moments,u
					ENDIF ELSE u.stat.m0=0.0
				END

				"M1PLOT" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						u.stat.zoom=0
						widget_control,u.id.pzoom,set_button=u.stat.zoom
						u.stat.m0=0.0
						u.stat.m1=1.0
						u.stat.m2=0.0
						widget_control,u.id.m0plot,set_button=u.stat.m0
						widget_control,u.id.m1plot,set_button=u.stat.m1
						widget_control,u.id.m2plot,set_button=u.stat.m2
						plot_time_moments,u
					ENDIF ELSE u.stat.m1=0.0
				END

				"M2PLOT" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						u.stat.zoom=0
						widget_control,u.id.pzoom,set_button=u.stat.zoom
						u.stat.m0=0.0
						u.stat.m1=0.0
						u.stat.m2=2.0
						widget_control,u.id.m0plot,set_button=u.stat.m0
						widget_control,u.id.m1plot,set_button=u.stat.m1
						widget_control,u.id.m2plot,set_button=u.stat.m2
						plot_time_moments,u
					ENDIF ELSE u.stat.m2=0.0				
				END
				"AUTOBIN_R" : BEGIN
					widget_control,u.id.nsub,get_value=nsub
					nsub=int(nsub[0])
					IF nsub LE 4 and is_dat(u) THEN BEGIN
						widget_control,u.id.noff,get_value=noff
						noff=int(noff[0])
						chmap=hirexsr_autochmap(nsub,noff)
						bnds=hirexsr_bin_bounds(chmap,u.stat.chmax)
						*u.dat.chmap=chmap
						*u.dat.bnds=bnds
						update_spatial_text,u
						IF u.stat.pos THEN plot_cx,u
						IF u.stat.image THEN plot_hemom_image,u

					ENDIF
				END
				"AUTOBIN_T" : BEGIN
					widget_control,u.id.ntsub,get_value=numframes
					numframes=int(numframes[0])
					IF is_dat(u) THEN BEGIN
						tmap=hirexsr_autotmap(n(u.dat.t)+1,numframes)
						cnts=hirexsr_bin_image(u.dat.raw,tmap)
						tau=hirexsr_bin_time(u.dat.t,tmap)
						widget_control,u.id.message,set_value='Re-binned @'+num2str(numframes,1)+' Images/Data Frame',/app
						u.index=ipt(tau,u.time)
						u.time=tau[u.index]
						u.dat.cnts=cnts
						u.dat.tau=tau
						u.dat.tmap=tmap	
						update_hemom_temporal_text,u
						IF u.stat.image THEN plot_hemom_image,u
					ENDIF
				END
				"REBIN_ALL" : BEGIN

				END
				"REBIN_KILL" : BEGIN
					IF u.stat.dat THEN BEGIN
						;update bins
						WIDGET_CONTROL, /hourglass
						hirexsr_rebin_afterkill,u.dat.cnts,u.dat.tau,*u.dat.chmap,*u.dat.tch,u.dat.lambda,u.dat.tmap,u.dat.spec,*u.dat.chmap_chk,$
							const=u.dat.const,fix=fix
						widget_control,u.id.message,set_value='Killed Pix Removed In '+num2str(total(fix),1)+' Spectra',/app
						*u.dat.chmap=*u.dat.chmap_chk
						IF total(fix) NE 0 THEN BEGIN
							x=size(fix)
							line=u.stat.line
							plot_hemom_image,u

							;update avespec
							avespec=hirexsr_ave_spec(u.dat.spec,u.dat.nave,avespec=u.dat.avespec,fix=fix)
							FOR i=0,x[1]-1 DO FOR j=0,x[2]-1 DO IF fix[i,j] THEN *u.dat.avespec[i,j]=*avespec[i,j]
							widget_control,u.id.message,set_value='Killed Spectra Averaged',/app
							plot_spec,u

							;update fit coefficients
							coefs=u.dat.coefs[*,*,line]
							icoefs=*coefs[u.ch-1,u.index]
							ncoefs=size(icoefs)
							CASE ncoefs[1] MOD 3 OF 
								0 : nback=2
								1 : nback=0
								2 : nback=1
							ENDCASE
							coefs=hirexsr_fit_spectra(avespec,line,double=u.stat.double,verb=verb,fits=coefs,fix=fix,nback=nback)
							FOR i=0,x[1]-1 DO FOR j=0,x[2]-1 DO IF fix[i,j] THEN *u.dat.coefs[i,j,line]=*coefs[i,j]
							widget_control,u.id.message,set_value='Killed Spectra Fit',/app

							;update moments
							imom=u.dat.mom[*,*,0:2,line]
							ierr=u.dat.mom[*,*,3:5,line]
							ipmom=u.dat.mom[*,*,6:8,line]
							iperr=u.dat.mom[*,*,9:11,line]
							irhot=u.dat.mom[*,*,12,line]
							ifrac=u.dat.mom[*,*,13,line]
							iscale=u.dat.mom[*,*,14,line]
							ifitcase=u.dat.mom[*,*,15,line]
							type=size(imom,/type)
							IF type EQ 4 THEN mom=fltarr(x[1],x[2],3,4) ELSE mom=dblarr(x[1],x[2],3,4)
							mom[*,*,*,0]=imom
							mom[*,*,*,1]=ierr
							mom[*,*,*,2]=ipmom
							mom[*,*,*,3]=iperr								
							m=hirexsr_spec2moments(u.dat.avespec,u.dat.coefs[*,*,line],*u.dat.labels[line],u.wl.lam_o[line],line,dlam=u.dat.dlam[line],$
								mom=mom,scale=iscale,bfrac=ifrac,fix=fix,fitcase=ifitcase)
							u.dat.mom[*,*,*,line]=[[[mom[*,*,*,0]]],[[mom[*,*,*,1]]],[[mom[*,*,*,2]]],[[mom[*,*,*,3]]],[[irhot]],[[ifrac]],[[iscale]],[[ifitcase]]]
							widget_control,u.id.message,set_value='Killed Moments Calculated',/app
							plot_all_moments,u
						ENDIF
					ENDIF
				END
				"SAVE_BINS" : BEGIN
					IF u.stat.dat THEN BEGIN
						save_bins,u
					ENDIF
				END
				"LOAD_BINS" : BEGIN
					IF u.stat.dat THEN BEGIN
						hirexsr_load_binning,u.shot,chmap,tch,tmap,good,chmax
						*u.dat.chmap_chk=chmap
						*u.dat.tch=tch
						u.dat.tmap=tmap
						widget_control,u.id.message,set_value='CHMAP, TCH and TMAP loaded: '+num2str(u.shot,1),/app
					ENDIF
				END
				"MT" : BEGIN
					IF is_dat(u) AND u.index NE 0 THEN BEGIN
						u.index-=1
						u.time=u.dat.tau[u.index]	
						update_hemom_temporal_text,u
						IF u.stat.pos OR u.stat.fs THEN plot_cx,u
						IF u.stat.image THEN plot_hemom_image,u
						plot_spec,u
						plot_all_moments,u
					ENDIF
				END
				"PT" : BEGIN
					IF is_dat(u) THEN BEGIN
						IF u.index NE max(u.dat.tmap) THEN BEGIN
							u.index+=1
							u.time=u.dat.tau[u.index]					
							update_hemom_temporal_text,u
							IF u.stat.pos OR u.stat.fs THEN plot_cx,u
							IF u.stat.image THEN plot_hemom_image,u
							plot_spec,u
							plot_all_moments,u
						ENDIF
					ENDIF
				END
				"MCH" : BEGIN
					IF is_dat(u) AND u.ch NE 1 THEN BEGIN
						u.ch-=1
						update_spatial_text,u
						widget_control,u.id.good,set_button=u.dat.good[u.ch-1,u.index]
						IF u.stat.pos OR u.stat.fs THEN plot_cx,u
						IF u.stat.image THEN plot_hemom_image,u
						plot_spec,u
						plot_all_moments,u
					ENDIF
				END
				"PCH" : BEGIN
					IF is_dat(u) THEN BEGIN
						chmap=*u.dat.chmap
						IF u.ch NE max(chmap[*,*,u.mapindex])+1 THEN BEGIN
							widget_control,u.id.good,set_button=u.dat.good[u.ch-1,u.index]
							u.ch+=1					
							update_spatial_text,u
							IF u.stat.pos OR u.stat.fs THEN plot_cx,u
							IF u.stat.image THEN plot_hemom_image,u
							plot_spec,u
							plot_all_moments,u
						ENDIF
					ENDIF
				END
				"WN3MOM" : BEGIN
					IF is_dat(u) AND u.stat.isline[0] THEN BEGIN
						IF event.select EQ 1 THEN u.stat.line=0
						reset_spec_buttons,u
						widget_control,u.id.wn3mom,set_button=1
						plot_spec,u
						plot_all_moments,u
						update_moment_text,u
					ENDIF ELSE BEGIN
						widget_control,u.id.message,set_value='No WN3 Data Available',/app
						widget_control,u.id.wn3mom,set_button=0
					ENDELSE
				END
				"XYMOM" : BEGIN
					IF is_dat(u) AND u.stat.isline[1] THEN BEGIN
						IF event.select EQ 1 THEN u.stat.line=1
						reset_spec_buttons,u
						widget_control,u.id.xymom,set_button=1
						plot_spec,u
						plot_all_moments,u
						update_moment_text,u
					ENDIF ELSE BEGIN
						widget_control,u.id.message,set_value='No XY Data Available',/app
						widget_control,u.id.xymom,set_button=0
					ENDELSE
				END
				"ZJKMOM" : BEGIN
					IF is_dat(u) AND u.stat.isline[2] THEN BEGIN
						IF event.select EQ 1 THEN u.stat.line=2
						reset_spec_buttons,u
						widget_control,u.id.zjkmom,set_button=1
						plot_spec,u
						plot_all_moments,u
						update_moment_text,u
					ENDIF ELSE BEGIN
						widget_control,u.id.message,set_value='No ZJK Data Available',/app
						widget_control,u.id.zjkmom,set_button=0
					ENDELSE
                                 END
				"LYAMOM" : BEGIN
					IF is_dat(u) AND u.stat.isline[3] THEN BEGIN
						IF event.select EQ 1 THEN u.stat.line=3
						reset_spec_buttons,u
						widget_control,u.id.lyamom,set_button=1
						plot_spec,u
						plot_all_moments,u
						update_moment_text,u
					ENDIF ELSE BEGIN
						widget_control,u.id.message,set_value='No LYA Data Available',/app
						widget_control,u.id.lyamom,set_button=0
                                            ENDELSE
				END
  				"MOMOM" : BEGIN
					IF is_dat(u) AND u.stat.isline[3] THEN BEGIN
						IF event.select EQ 1 THEN u.stat.line=4
						reset_spec_buttons,u
						widget_control,u.id.momom,set_button=1
						plot_spec,u
						plot_all_moments,u
						update_moment_text,u
					ENDIF ELSE BEGIN
						widget_control,u.id.message,set_value='No LYA Data Available',/app
						widget_control,u.id.lyamom,set_button=0
					ENDELSE
				END                                  

				"M_MOMPLOT" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						u.stat.mplot=1
						widget_control,u.id.m_intplot,set_button=0
						plot_all_moments,u	
					ENDIF
				END
				"M_INTPLOT" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						u.stat.mplot=0
						widget_control,u.id.m_momplot,set_button=0
						plot_all_moments,u	
					ENDIF
				END
				"GOOD" : BEGIN
					IF is_dat(u) THEN BEGIN
						IF event.select EQ 1 THEN u.dat.good[u.ch-1,u.index]=1 ELSE u.dat.good[u.ch-1,u.index]=0
						plot_all_moments,u
					ENDIF
				END
				"REAVE" : BEGIN
					IF is_dat(u) THEN BEGIN
						type=size(u.dat.spec,/type)
						IF type EQ 10 THEN BEGIN
							avespec=hirexsr_ave_spec(u.dat.spec,u.dat.nave,maxave=maxave)
							u.stat.maxave=maxave	
							u.dat.avespec=avespec
							widget_control,u.id.message,set_value='Averaged Spectra Over '+num2str(u.dat.nave,1)+' Points',/app
							plot_spec,u
							plot_all_moments,u
						ENDIF ELSE widget_control,u.id.message,set_value='NOT AVAILABLE IN DISPLAY MODE',/app
					ENDIF
				END
				"SAVE_AVE" : BEGIN
					IF u.stat.dat THEN BEGIN
						save_spec,u
					ENDIF
				END
				"LOAD_AVE" : BEGIN
					widget_control,u.id.message,set_value='NOT FUNCTIONAL: USE LOAD',/app
				END

				"FIT_SPEC" : BEGIN
					IF u.stat.dat THEN BEGIN
						fix=intarr(u.stat.chmax,u.stat.ntime)
						fix[u.ch-1,u.index]=1
						x=size(fix)
						line=u.stat.line
						coefs=u.dat.coefs[*,*,line]
						coefs=hirexsr_fit_spectra(avespec,line,double=u.stat.double,verb=verb,fits=coefs,fix=fix)
						FOR i=0,x[1]-1 DO FOR j=0,x[2]-1 DO IF fix[i,j] THEN *u.dat.coefs[i,j,line]=*coefs[i,j]
						widget_control,u.id.message,set_value='Current LINE='+num2str(line,1)+' Spectra Refit',/app
						plot_spec,u
					ENDIF
				END			
				"FIT_ALLSPEC" : BEGIN
					IF u.stat.dat THEN BEGIN
						WIDGET_CONTROL, /hourglass
						line=u.stat.line
						coefs=hirexsr_fit_spectra(u.dat.avespec,line,double=u.stat.double,verb=verb,fix=fix)	;rerun all the fits for the current line
						x=size(fix)
						FOR i=0,x[1]-1 DO FOR j=0,x[2]-1 DO IF fix[i,j] THEN *u.dat.coefs[i,j,line]=*coefs[i,j]
						widget_control,u.id.message,set_value='All LINE='+num2str(line,1)+' Spectra Refit',/app
						plot_spec,u
					ENDIF
				END
				"SAVE_FITS" : BEGIN
					IF u.stat.dat THEN BEGIN
						save_fits,u	
					ENDIF
				END
				"LOAD_FITS" : BEGIN
					widget_control,u.id.message,set_value='NOT FUNCTIONAL: USE LOAD',/app
				END
				"CALC_MOM" : BEGIN
					IF u.stat.dat THEN BEGIN
						fix=intarr(u.stat.chmax,u.stat.ntime)
						fix[u.ch-1,u.index]=1
						x=size(fix)
						line=u.stat.line
						imom=u.dat.mom[*,*,0:2,line]
						ierr=u.dat.mom[*,*,3:5,line]
						ipmom=u.dat.mom[*,*,6:8,line]
						iperr=u.dat.mom[*,*,9:11,line]
						irhot=u.dat.mom[*,*,12,line]
						ifrac=u.dat.mom[*,*,13,line]
						iscale=u.dat.mom[*,*,14,line]
						ifitcase=u.dat.mom[*,*,15,line]
						type=size(imom,/type)
						IF type EQ 4 THEN mom=fltarr(x[1],x[2],3,4) ELSE mom=dblarr(x[1],x[2],3,4)
						mom[*,*,*,0]=imom
						mom[*,*,*,1]=ierr
						mom[*,*,*,2]=ipmom
						mom[*,*,*,3]=iperr
						mom=hirexsr_spec2moments(u.dat.avespec,u.dat.coefs[*,*,line],*u.dat.labels[line],u.wl.lam_o[line],line,dlam=u.dat.dlam[line],$
							mom=mom,scale=iscale,bfrac=ifrac,fix=fix,fitcase=ifitcase)
						u.dat.mom[*,*,*,line]=[[[mom[*,*,*,0]]],[[mom[*,*,*,1]]],[[mom[*,*,*,2]]],[[mom[*,*,*,3]]],[[irhot]],[[ifrac]],[[iscale]],[[ifitcase]]]
						widget_control,u.id.message,set_value='Current LINE='+num2str(line,1)+' Momments Recalculated',/app
						plot_all_moments,u
					ENDIF
				END
				"CALC_ALLMOM" : BEGIN
					IF u.stat.dat THEN BEGIN
						WIDGET_CONTROL, /hourglass
						line=u.stat.line
						irhot=u.dat.mom[*,*,12,line]
						mom=hirexsr_spec2moments(u.dat.avespec,u.dat.coefs[*,*,line],*u.dat.labels[line],u.wl.lam_o[line],line,$
							dlam=u.dat.dlam[line],bfrac=ifrac,scale=iscale,fitcase=ifitcase)
						u.dat.mom[*,*,*,line]=[[[mom[*,*,*,0]]],[[mom[*,*,*,1]]],[[mom[*,*,*,2]]],[[mom[*,*,*,3]]],[[irhot]],[[ifrac]],[[iscale]],[[ifitcase]]]
						widget_control,u.id.message,set_value='All LINE='+num2str(line,1)+' Moments Recalculated',/app
						plot_all_moments,u
					ENDIF
				END
				"SAVE_MOM" : BEGIN
					IF u.stat.dat THEN BEGIN
						save_mom,u
					ENDIF
				END

				"LOAD_MOM" : BEGIN
					widget_control,u.id.message,set_value='NOT FUNCTIONAL: USE LOAD',/app
				END			


				"PRINT" : BEGIN
					stop
				END 
				"DISPLAY" : IF event.select EQ 1 THEN u.stat.display=1 ELSE u.stat.display=0
			
			ELSE:
		ENDCASE
		END
  		"WIDGET_SLIDER": BEGIN
			widget_control,event.id,get_value=slider
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF 
				'CH_SLIDER' : BEGIN
					u.ch=slider
                                        IF is_dat(u) THEN BEGIN
						tch=*u.dat.tch
						time=u.dat.tau[u.index]
						tempbnds = *u.dat.bnds
						IF n(tch) NE 0 THEN BEGIN
							tmp=where(tch LE time)
							tempbnds=tempbnds[*,*,last(tmp)]
						ENDIF
                                        	
                                            	IF tempbnds[0,u.ch-1] NE -1 THEN BEGIN
                                              		update_spatial_text,u
                                                	widget_control,u.id.good,set_button=u.dat.good[u.ch-1,u.index]
	                                                plot_hemom_image,u
        	                                        plot_spec,u
                	                                plot_all_moments,u
							IF u.stat.pos THEN plot_cx,u                                                
                                	        ENDIF ELSE BEGIN
                                                	max = where(tempbnds(0, *) GT 0)
	                                                max = FIX(max(n_elements(max)-1))
        	                                        WIDGET_CONTROL,u.id.ch_slider, set_value=max+1
                	                       	ENDELSE
                                        ENDIF
				END
				'GAIN_SLIDER' : IF is_dat(u) THEN plot_hemom_image,u
				'T_SLIDER' : BEGIN
					IF is_dat(u) THEN BEGIN
						index=ipt(u.dat.tau,slider/1.0e3)
						tmp=where(u.dat.tau GT 0)
						IF slider/1.0e3 GE max(u.dat.tau[tmp]) THEN index=n(tmp)
						IF slider/1.0e3 LE min(u.dat.tau[tmp]) THEN index=0
						IF u.time NE u.index THEN BEGIN
							update_hemom_temporal_text,u
							u.time=u.dat.tau[index]
							u.index=index
							widget_control,id.tpt,set_value=num2str(u.time,dp=2)
							IF u.stat.pos OR u.stat.fs THEN plot_cx,u
							IF u.stat.image THEN plot_hemom_image,u
							widget_control,u.id.good,set_button=u.dat.good[u.ch-1,u.index]
							plot_spec,u
							plot_all_moments,u
						ENDIF
					ENDIF
				END
				'LAM_SLIDER' : BEGIN
					IF is_dat(u) THEN BEGIN
						newdlam=slider/10.0
						u.dat.dlam[u.stat.line]=newdlam
						update_moment_text,u
						plot_spec_moments,u
						IF u.stat.mconv THEN plot_conv_moments,u
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
				u.id.thtid : BEGIN
					widget_control,u.id.thtid,get_value=tht
					IF hirexsr_is_analysis(u.shot,int(tht[0])) THEN BEGIN
						u.tht=int(tht[0])
					ENDIF ELSE BEGIN
						widget_control,u.id.message,set_value='invalid THT',/app
						widget_control,u.id.thtid,set_value=num2str(u.tht,1)
					ENDELSE
				END
				u.id.nave : BEGIN
					widget_control,u.id.nave,get_value=nave
					u.dat.nave=int(nave)
				END
				u.id.vmaxplt : BEGIN
					widget_control,u.id.vmaxplt,get_value=vmax
					u.stat.vplt[1]=vmax
					IF u.stat.dat THEN plot_spatial_moments,u
				END
				u.id.vminplt : BEGIN
					widget_control,u.id.vminplt,get_value=vmin
					u.stat.vplt[0]=vmin
					IF u.stat.dat THEN plot_spatial_moments,u
				END
				u.id.tmaxplt : BEGIN
					widget_control,u.id.tmaxplt,get_value=tmax
					u.stat.tplt[1]=tmax
					IF u.stat.dat THEN plot_spatial_moments,u
				END
				u.id.tminplt : BEGIN
					widget_control,u.id.tminplt,get_value=tmin
					u.stat.tplt[0]=tmin
					IF u.stat.dat THEN plot_spatial_moments,u
				END

			ELSE: 
			ENDCASE
		END

               "WIDGET_DRAW": BEGIN
                    CASE event.id OF
                        u.id.draw5 : BEGIN
                            IF WIDGET_INFO(u.id.pzoom, /button_set) EQ 1 AND u.stat.display EQ 0 THEN BEGIN
                                IF event.press EQ 1 THEN BEGIN
                                    ptru = PTR_NEW(u) ; ptr is killed in hirexsr_zoom_kill
                                    hirexsr_zoom_kill, ptru 
                                ENDIF 
                            ENDIF ELSE IF u.stat.display THEN widget_control,u.id.message,set_value='NOT AVAILABLE IN DISPLAY MODE',/app
                        END

                        ELSE: BEGIN
                             ; ODD?
                        END
                    ENDCASE
                END
		ELSE:
	ENDCASE
	IF button NE 'QUIT' THEN widget_control,event.top,set_uvalue=u
END  		

;+
;NAME:
;	W_HIREXSR_HE_MOMENTS
;
;-


PRO w_hirexsr_he_moments,shot=shot,tht=tht,load=load
	user=logname()

	loadct,12,/silent
	base=widget_base(title='HIREXSR He-Like Moment Analysis',/row,tlb_size_events=1, kill_notify= 'cleanheap_w_hirexsr_he_moments')
	A=widget_base(base,/column)
	B=widget_base(base,/column)
	space=widget_base(base,/row,xsize=5)	
	C=widget_base(base,/column)
	space=widget_base(base,/row,xsize=5)	
	D=widget_base(base,/column)

	vzoom=0.6
	dum = widget_label(A,value='POS VECTORS')
	A1=widget_base(A,frame=5)
	draw1=widget_draw(A1,xsize=500*vzoom,ysize=818*vzoom)
	Aprime=widget_base(A,/column,/frame)
	Ax=widget_base(Aprime,/row)
	dum = widget_label(Ax,value='PLOT: ')
	Axp1=widget_base(Ax,/row,/nonexclusive)
	pimage=widget_button(Axp1,value='IMAGE')
	pfs=widget_button(Axp1,value='FS')
	ppos=widget_button(Axp1,value='POS')
	pbin=widget_button(Axp1,value='BINS')
	pzoom=widget_button(Axp1,value='ZOOM')
	Ay=widget_base(Aprime,/row)
	dum = widget_label(Ay,value='      ')
	Ayp1=widget_base(Ay,/row,/nonexclusive)
	pave=widget_button(Ayp1,value='AVE')
	pconv=widget_button(Ayp1,value='MCONV')
	

	A2=widget_base(A,/column,/frame)
	dum = widget_label(A2,value='IMAGE BINNING SETUP')
	A2tb=widget_tab(A2,location=3)

	A2t1=widget_base(A2tb,title=' SPATIAL ',/column,group_leader=base,/frame)
	dum = widget_label(A2t1,value='SPATIAL SETUP')
	A2p1=widget_base(A2t1,/row)
	dum = widget_label(A2p1,value='CH# ')
	ch=widget_text(A2p1,xsize=2,ysize=1)
	dum = widget_label(A2p1,value=' OF ')
	chtot=widget_text(A2p1,xsize=2,ysize=1)
	dum = widget_label(A2p1,value='  ')
	mch=widget_button(A2p1,value=' - ')
	pch=widget_button(A2p1,value=' + ')
	A2p1x=widget_base(A2t1,/row)	
	dum = widget_label(A2p1x,value='USING  ROWS: ')
	rlow=widget_text(A2p1x,xsize=4,ysize=1,/edit)
	dum = widget_label(A2p1x,value='   TO   ')
	rhigh=widget_text(A2p1x,xsize=4,ysize=1,/edit)
	A2p2=widget_base(A2t1,/row)
	dum = widget_label(A2p2,value='           ')
	mlow=widget_button(A2p2,value=' - ')
	plow=widget_button(A2p2,value=' + ')
	dum = widget_label(A2p2,value='    ')
	mhigh=widget_button(A2p2,value=' - ')
	phigh=widget_button(A2p2,value=' + ')
	A2p5=widget_base(A2t1,/row,/nonexcl)
	good=widget_button(A2p5,value='INCLUDE CHANNEL IN ANALYSIS (GOOD)')
	A2p3=widget_base(A2t1,/row)
	ch_del=widget_button(A2p3,value='DELETE')
	ch_ins=widget_button(A2p3,value='INSERT')
	mcomb=widget_button(A2p3,value=' - ')
	dum = widget_label(A2p3,value='COMBINE w/ CH')
	pcomb=widget_button(A2p3,value=' + ')
	A2p6=widget_base(A2t1,/row)
	autobin_r = widget_button(A2p6,value='AUTOBIN' )
	dum=widget_label(A2p6,value=' @ ')
	nsub=widget_text(A2p6,xsize=2,ysize=1,/edit)
	dum = widget_label(A2p6,value='/SUB w/ ')
	noff=widget_text(A2p6,xsize=2,ysize=1,/edit)
	dum = widget_label(A2p6,value=' OFFSET ')	

	A2t2=widget_base(A2tb,title=' TEMPORAL ',/column,group_leader=base,/frame)
	dum = widget_label(A2t2,value='TEMPORTAL SETUP')
	A2x2=widget_base(A2t2,/row)
	dum = widget_label(A2x2,value='TIME BIN ')
	itpt=widget_text(A2x2,xsize=3,ysize=1)
	dum = widget_label(A2x2,value=' OF ')
	ntpt=widget_text(A2x2,xsize=3,ysize=1)
	dum = widget_label(A2x2,value='  ')
	mt=widget_button(A2x2,value=' - ')
	pt=widget_button(A2x2,value=' + ')
	A2x2=widget_base(A2t2,/row)
	dum = widget_label(A2x2,value='USING FRAMES: ')
	flow=widget_text(A2x2,xsize=3,ysize=1,/edit)
	dum = widget_label(A2x2,value='   TO   ')
	fhigh=widget_text(A2x2,xsize=3,ysize=1,/edit)
	A2x4=widget_base(A2t2,/row)
	dum = widget_label(A2x4,value='            ')
	mtlow=widget_button(A2x4,value=' - ')
	ptlow=widget_button(A2x4,value=' + ')
	dum = widget_label(A2x4,value='    ')
	mthigh=widget_button(A2x4,value=' - ')
	pthigh=widget_button(A2x4,value=' + ')

	A2x3=widget_base(A2t2,/row)
	t_del=widget_button(A2x3,value='DELETE')
	t_ins=widget_button(A2x3,value='INSERT')
	mtcomb=widget_button(A2x3,value=' - ')
	dum = widget_label(A2x3,value='COMBINE w/ BIN')
	ptcomb=widget_button(A2x3,value=' + ')
	A2x5=widget_base(A2t2,/row)
	autobin_t = widget_button(A2x5,value='AUTOBIN')
	dum=widget_label(A2x5,value=' @ ')
	ntsub=widget_text(A2x5,xsize=2,ysize=1,/edit)
	dum = widget_label(A2x5,value='FRAMES/TIME BIN')
	A2x6=widget_base(A2t2,/row)
	dum = widget_label(A2x6,value='T_LOW')
	tbinlow=widget_text(A2x6,xsize=4,ysize=1,/edit)
	dum = widget_label(A2x6,value='   T_HIGH')
	tbinhigh=widget_text(A2x6,xsize=4,ysize=1,/edit)

	A2y=widget_base(A2,/row)
	dum = widget_label(A2y,value='BINNING: ')
	load_bins=widget_button(A2y,value='LOAD')
	save_bins=widget_button(A2y,value='SAVE')
	dum = widget_label(A2y,value='  REFILL: ')
	rebin_all=widget_button(A2y,value='ALL')
	rebin_kill=widget_button(A2y,value='KILLED')


	A3=widget_base(A,/column,/frame)
	A3p1=widget_base(A3,/row)
	dum = widget_label(A3p1,value='SHOT:')
	shotid = widget_text(A3p1,xsize=10,ysize=1,/edit)
	dum = widget_label(A3p1,value='THT')
	thtid=widget_text(A3p1,xsize=2,ysize=1,/edit)
	save= widget_button(A3p1,value='SAVE')
	load= widget_button(A3p1,value='LOAD')
	quit= widget_button(A3p1,value='QUIT')

	A3p4=widget_base(A3,/row)
	message = widget_text(A3p4,xsize=40,ysize=6,/scroll)

	imzoom=0.7
	dum = widget_label(B,value='     RAW IMAGE')
	B1=widget_base(B,/row)
	ch_slider=widget_slider(B1,ysize=3*487*imzoom,min=1,max=24*4,value=24.,/drag,/vert)
	B1p1=widget_base(B1,frame=5)
	draw2=widget_draw(B1p1,xsize=195*imzoom,ysize=3*487*imzoom)
	gain_slider=widget_slider(B,xsize=195*imzoom,min=1,max=40,value=1,/drag)

	dum = widget_label(C,value='SPECTRAL LINE FITTING')
	C1=widget_base(C,frame=5)
	draw3=widget_draw(C1,xsize=600,ysize=300)
	space=widget_base(C,/row,ysize=2)	
	C2=widget_base(C,/row)
	C2p1=widget_base(C2,frame=5)
	draw4=widget_draw(C2p1,xsize=300,ysize=360)
	space=widget_base(C2,/row,xsize=4)
	C2p2=widget_base(C2,/column)

	C3=widget_base(C2p2,/column,/frame)
	C3x=widget_base(C3,/row)
	C3p1=widget_base(C3x,/row,/nonexcl)
	double=widget_button(C3p1,value='DBL')
	avespec=widget_button(C3p1,value='AVE')
	dum = widget_label(C3x,value='#AVE ')
	nave=widget_text(C3x,xsize=2,ysize=1,/edit)
	reave=widget_button(C3x,value='RECALC')
	C3y=widget_base(C3,/row)
	dum = widget_label(C3y,value='   AVESPEC TO TREE: ')
	save_ave=widget_button(C3y,value='SAVE')
	load_ave=widget_button(C3y,value='LOAD')

	C3p2=widget_base(C3,/row)
	dum = widget_label(C3p2,value='   FIT SPECTRA FOR: ')
	fit_spec=widget_button(C3p2,value='CURRENT')
	fit_allspec=widget_button(C3p2,value=' ALL ')
	C3p3=widget_base(C3,/row)
	dum = widget_label(C3p3,value='      FITS TO TREE: ')
	save_fits=widget_button(C3p3,value='SAVE')
	load_fits=widget_button(C3p3,value='LOAD')

	space=widget_base(C2p2,/row,ysize=2)	

	C4=widget_base(C2p2,/column,/frame)
	dum = widget_label(C4,value='MOMENT CALCULATION SETUP')
	C4p1=widget_base(C4,/row)
	dum = widget_label(C4p1,value=' DLAM ')
	lam_slider = widget_slider(C4p1,xsize=100,min=1,max=100,value=40,/drag,/suppress)
	dlampt=widget_text(C4p1,xsize=3,ysize=1)
	dum = widget_label(C4p1,value=' [mAng] ')
	C4p2=widget_base(C4,/row)
	dum = widget_label(C4p2,value=' LOW ')

	lowpt=widget_text(C4p2,xsize=6,ysize=1)
	dum = widget_label(C4p2,value=' HIGH ')
	highpt=widget_text(C4p2,xsize=6,ysize=1)
	dum = widget_label(C4p2,value=' [Ang] ')
	C4p3=widget_base(C4,/row,/nonexcl)
	wn3mom=widget_button(C4p3,value='wn3')
	xymom=widget_button(C4p3,value='xy')
	zjkmom=widget_button(C4p3,value='zjk')
	lyamom=widget_button(C4p3,value='lya')
	momom=widget_button(C4p3,value='mo4d')
	
	C4p4=widget_base(C4,/row)
	dum = widget_label(C4p4,value='CALCULATE MOMENTS FOR: ')
	calc_mom=widget_button(C4p4,value='CURRENT')
	calc_allmom=widget_button(C4p4,value=' ALL ')
	C4p5=widget_base(C4,/row)
	dum = widget_label(C4p5,value='      MOMENTS TO TREE: ')
	save_mom=widget_button(C4p5,value='SAVE')
	load_mom=widget_button(C4p5,value='LOAD')



	
	C5=widget_base(C,/row,/frame)
	dum = widget_label(C5,value='MOMENT TIME HISTORY: ')
	C5p1=widget_base(C5,/row,/nonexcl)
	m0plot=widget_button(C5p1,value=' BRIGHT ')
	m1plot=widget_button(C5p1,value=' VELOCITY ')
	m2plot=widget_button(C5p1,value=' TEMPERATURE ')
	C5p2=widget_base(C,frame=5)
	draw5=widget_draw(C5p2,xsize=600,ysize=250, /button_events)
	    	
	C6=widget_base(C,/row)
	dum = widget_label(C6,value='TIME: ')
	tpt=widget_text(C6,xsize=5,ysize=1)
	t_slider=widget_slider(C6,xsize=500,min=0,max=2000,value=1000,/drag,/suppress)
	C7=widget_base(C,/row)
	C7a=widget_base(C7,/row)
	dum = widget_label(C7a,value='# OF FRAMES COLLECTED: ')
	nframes=widget_text(C7a,xsize=4,ysize=1)
	dum = widget_label(C7a,value='   ')
	dum = widget_label(C7a,value='MAX # OF CHANNELS: ')
	chmax=widget_text(C7a,xsize=4,ysize=1)
	print= widget_button(C7,value='PRINT')
	dum = widget_label(C7a,value='   ')
	C7b=widget_base(C7,/row,/nonexcl)
	display=widget_button(c7b, value=' DISPLAY MODE ')


	dum = widget_label(D,value='MOMENT/LINT-INT PROFILES')
	D1=widget_base(D,frame=5)
	draw6=widget_draw(D1,xsize=370,ysize=900)
	D2=widget_base(D,/row)
	dum = widget_label(D2,value='X-AXIS SELECTION: ')
	D2p1=widget_base(D2,/row,/nonexcl)
	m_chplot=widget_button(D2p1,value='CH #')
	m_rmidplot=widget_button(D2p1,value='RMID')
	D3=widget_base(D,/row)
	dum = widget_label(D3,value='Y-AXIS SELECTION: ')
	D3p1=widget_base(D3,/row,/nonexcl)
	m_momplot=widget_button(D3p1,value='MOMENTS')
	m_intplot=widget_button(D3p1,value='LINE-INTEGRATED')
	D4=widget_base(D,/row)
	dum = widget_label(D4,value='V_MIN ')
	vminplt=widget_text(D4,xsize=6,ysize=1,/edit)	
	dum = widget_label(D4,value='V_MAX ')
	vmaxplt=widget_text(D4,xsize=6,ysize=1,/edit)	
	dum = widget_label(D4,value=' [km/s] ')
	D5=widget_base(D,/row)
	dum = widget_label(D5,value='T_MIN ')
	tminplt=widget_text(D5,xsize=6,ysize=1,/edit)	
	dum = widget_label(D5,value='T_MAX ')
	tmaxplt=widget_text(D5,xsize=6,ysize=1,/edit)	
	dum = widget_label(D5,value=' [keV] ')

	
	
	;build u structuru
	id={base:base,draw1:draw1,draw2:draw2,draw3:draw3,draw4:draw4,draw5:draw5,draw6:draw6,$
		pimage:pimage,pfs:pfs,ppos:ppos,pbin:pbin,pzoom:pzoom,pave:pave,pconv:pconv,$
		ch:ch,chtot:chtot,rlow:rlow,rhigh:rhigh,$
		mch:mch,pch:pch,mlow:mlow,plow:plow,mhigh:mhigh,phigh:phigh,$
		good:good,ch_del:ch_del,ch_ins:ch_ins,mcomb:mcomb,pcomb:pcomb,$
		autobin_r:autobin_r,nsub:nsub,noff:noff,$
		itpt:itpt,ntpt:ntpt,mt:mt,pt:pt,flow:flow,fhigh:fhigh,$
		mtlow:mtlow,ptlow:ptlow,mthigh:mthigh,pthigh:pthigh,$
		t_del:t_del,t_ins:t_ins,mtcomb:mtcomb,ptcomb:ptcomb,$
		autobin_t:autobin_t,ntsub:ntsub,tbinlow:tbinlow,tbinhigh:tbinhigh,$
		save_bins:save_bins,load_bins:load_bins,rebin_all:rebin_all,rebin_kill:rebin_kill,$
		shotid:shotid,thtid:thtid,save:save,load:load,quit:quit,print:print,message:message,$
		ch_slider:ch_slider,gain_slider:gain_slider,$
		double:double,avespec:avespec,nave:nave,reave:reave,save_ave:save_ave,load_ave:load_ave,$
		fit_spec:fit_spec,fit_allspec:fit_allspec,save_fits:save_fits,load_fits:load_fits,$
		lowpt:lowpt,lam_slider:lam_slider,highpt:highpt,dlampt:dlampt,$
			wn3mom:wn3mom,xymom:xymom,zjkmom:zjkmom,lyamom:lyamom,momom:momom,$
			calc_mom:calc_mom,calc_allmom:calc_allmom,save_mom:save_mom,load_mom:load_mom,$
		m0plot:m0plot,m1plot:m1plot,m2plot:m2plot,$
		tpt:tpt,t_slider:t_slider,nframes:nframes,chmax:chmax,display:display,$
		m_chplot:m_chplot,m_rmidplot:m_rmidplot,m_momplot:m_momplot,m_intplot:m_intplot,$
		vminplt:vminplt,vmaxplt:vmaxplt,tminplt:tminplt,tmaxplt:tmaxplt}

	;set defaults
	IF NOT keyword_set(shot) THEN shot=1070830020
	IF NOT keyword_set(tht) THEN tht=0
	time=1.0
	ch=24
	ntsub=1
	nsub=2
	noff=10
	nave=0
	tplt=[0,2.0]
	vplt=[-10,10]
	stat={dat:0,double:0,image:1,fs:1,pos:1,bins:0,zoom:0,m0:0,m1:1,m2:0,chmax:96,ntime:0,nmaps:1,maxave:0,line:2,ave:1,display:1,mplot:1,chplot:1,mconv:0,$
		tplt:tplt,vplt:vplt,linenum:[0,1,2,7,8],isline:[0,0,0,0,0],ps:0}
	plot={xsize:[500*vzoom,195*imzoom,600,300,600,370],ysize:[818*vzoom,487*3*imzoom,300,360,250,900]}
	hirexsr_load_wavelengths,shot,lam_o,z,label
	tmp=where(z EQ 18)
	lam_o=[lam_o[where(label EQ 'w' AND z EQ 18)],lam_o[where(label EQ 'x' AND z EQ 18)],lam_o[where(label EQ 'z' AND z EQ 18)],lam_o[where(label EQ 'lya1' AND z EQ 18)],lam_o[where(label EQ '4d' AND z EQ 42)]]
	label=['w','x','z','lya1','4d']
	wl={wn3:[3.943,3.96],xy:[3.96,3.975],zjk:[3.977,4.0],lya:[3.723,3.749],mo4d:[3.723,3.749],lam_o:lam_o,label:label}

	u={id:id,shot:shot,tht:tht,time:time,ch:ch,index:0,mapindex:0,stat:stat,wl:wl,plot:plot}
	widget_control,base,set_uvalue=u
	widget_control,id.shotid,set_value=num2str(u.shot,1)
	widget_control,id.thtid,set_value=num2str(u.tht,1)
	widget_control,id.tpt,set_value=num2str(u.time,dp=2)
	widget_control,id.t_slider,set_value=u.time*1.0e3
	widget_control,id.ch_slider,set_value=u.ch
	widget_control,u.id.double,set_button=u.stat.double
	widget_control,u.id.pimage,set_button=u.stat.image
	widget_control,u.id.pfs,set_button=u.stat.fs
	widget_control,u.id.ppos,set_button=u.stat.pos
	widget_control,u.id.pbin,set_button=u.stat.bins
	widget_control,u.id.pzoom,set_button=u.stat.zoom
	widget_control,u.id.pave,set_button=u.stat.ave
	widget_control,u.id.pconv,set_button=u.stat.mconv
	widget_control,u.id.nsub,set_value=num2str(nsub,1)
	widget_control,u.id.noff,set_value=num2str(noff,1)
	widget_control,u.id.ntsub,set_value=num2str(ntsub,1)
	widget_control,u.id.nave,set_value=num2str(nave,1)
	widget_control,u.id.m0plot,set_button=stat.m0
	widget_control,u.id.m1plot,set_button=stat.m1
	widget_control,u.id.m2plot,set_button=stat.m2
	widget_control,u.id.m_momplot,set_button=stat.mplot
	widget_control,u.id.m_chplot,set_button=stat.chplot
	widget_control,u.id.display,set_button=u.stat.display
	widget_control,u.id.vminplt,set_value=num2str(u.stat.vplt[0],dp=1)
	widget_control,u.id.vmaxplt,set_value=num2str(u.stat.vplt[1],dp=1)
	widget_control,u.id.tminplt,set_value=num2str(u.stat.tplt[0],dp=1)
	widget_control,u.id.tmaxplt,set_value=num2str(u.stat.tplt[1],dp=1)


	CASE u.stat.line OF
		0 : widget_control,u.id.wn3mom,set_button=1
		1 : widget_control,u.id.xymom,set_button=1
		2 : widget_control,u.id.zjkmom,set_button=1
		7 : widget_control,u.id.lyamom,set_button=1
		8 : widget_control,u.id.momom,set_button=1
		ELSE : 
	END


	!except=0			
	widget_control,base,/realize
	xmanager,'w_hirexsr_he_moments',base
END
