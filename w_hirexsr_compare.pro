PRO comp_plot_emiss,u
	IF u.stat.ps THEN BEGIN
		xsize=5.0
		ysize=5.0*900/700.0
		ls=0.65
		col=u.plot.pscol
	ENDIF ELSE BEGIN
		ls=1.0
		col=u.plot.col
		widget_control,u.id.draw1,get_value=draw_win
		window,0,xsize=u.plot.size[0],ysize=u.plot.size[1],/pixmap
	ENDELSE
	index=u.index
	xr=u.plot.xr
	yr=u.plot.em
	IF u.plot.ratio THEN ytit='EMISS RATIO m=1/m=0' ELSE ytit='EMISS'
	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit=u.plot.xtit,ytit=ytit
	FOR i=0,n(u.stat.line) DO BEGIN
		IF u.stat.inv[i] THEN BEGIN
			ipro=*u.dat.inv[i,index[i]]
			npts=n(ipro[*,8])+1
			IF ipro[npts-1,8] EQ ipro[npts/2-1,8] THEN nrho=npts/2 ELSE nrho=npts
			CASE u.plot.x OF
				0 : xplot=ipro[0:nrho-1,8]
				1 : BEGIN
					eindex=ipt(u.dat.etime,u.stat.time)
					rmid=reform(u.dat.rmid[eindex,*])
					xplot=interpol(rmid,u.dat.epsin,ipro[0:nrho-1,8])
				END
				2 : BEGIN
					eindex=ipt(u.dat.etime,u.stat.time)
					rmid=reform(u.dat.rmid[eindex,*])
					xplot=(interpol(rmid,u.dat.epsin,ipro[0:nrho-1,8])-rmid[0])/(last(rmid)-rmid[0])
				END
			ENDCASE
			IF u.plot.ratio AND npts NE nrho THEN BEGIN
				ratio=ipro[nrho:*,0]/ipro[0:nrho-1,0]
				raterr=ratio*sqrt(ipro[0:nrho-1,4]^2/ipro[nrho:*,0]^2+ipro[nrho:*,4]^2/ipro[nrho:*,0]^2)
				oploterror,xplot,ratio,raterr,color=col[i],errcolor=col[i]
			ENDIF ELSE BEGIN	
				IF NOT u.plot.ratio THEN oploterror,xplot,ipro[0:nrho-1,0],ipro[0:nrho-1,4],color=col[i],errcolor=col[i]
				IF npts NE nrho THEN oploterror,xplot,ipro[nrho:*,0],ipro[nrho:*,4],color=col[i],errcolor=col[i],linestyle=3.0
			ENDELSE
		ENDIF
	ENDFOR
	oplot,xr,[0,0],linestyle=2,color=col[4]

	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.size[0],u.plot.size[1],0,0,0]
	ENDIF
END

PRO comp_plot_rot,u
	IF u.stat.ps THEN BEGIN
		xsize=5.0
		ysize=5.0*900/700.0
		ls=0.65
		col=u.plot.pscol
	ENDIF ELSE BEGIN
		ls=1.0
		col=u.plot.col
		widget_control,u.id.draw2,get_value=draw_win
		window,0,xsize=u.plot.size[0],ysize=u.plot.size[1],/pixmap
	ENDELSE
	index=u.index
	xr=u.plot.xr
	yr=u.plot.vt
	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit=u.plot.xtit,ytit=n2g('omega')+' [kHz]'
	FOR i=0,n(u.stat.line) DO BEGIN
		IF u.stat.inv[i] THEN BEGIN
			ipro=*u.dat.inv[i,index[i]]
			npts=n(ipro[*,8])+1
			IF ipro[npts-1,8] EQ ipro[npts/2-1,8] THEN nrho=npts/2 ELSE nrho=npts
			CASE u.plot.x OF
				0 : xplot=ipro[0:nrho-1,8]
				1 : BEGIN
					eindex=ipt(u.dat.etime,u.stat.time)
					rmid=reform(u.dat.rmid[eindex,*])
					xplot=interpol(rmid,u.dat.epsin,ipro[0:nrho-1,8])
				END
				2 : BEGIN
					eindex=ipt(u.dat.etime,u.stat.time)
					rmid=reform(u.dat.rmid[eindex,*])
					xplot=(interpol(rmid,u.dat.epsin,ipro[0:nrho-1,8])-rmid[0])/(last(rmid)-rmid[0])
				END
			ENDCASE	
			oploterror,xplot,ipro[0:nrho-1,1],ipro[0:nrho-1,5],color=col[i],errcolor=col[i]
			IF npts NE nrho THEN oploterror,xplot,ipro[nrho:*,1],ipro[nrho:*,5],color=col[i],errcolor=col[i],linestyle=3.0
		ENDIF
		IF u.stat.lint[i] THEN BEGIN
			IF u.plot.sym[i] EQ 8 THEN makesym,u.plot.psym[i]
			ilint=*u.dat.lint[i,index[i]]
			tmp=where(finite(ilint[*,0]) EQ 1 AND finite(ilint[*,1]) EQ 1)
			CASE u.plot.x OF
				0 : xplot=ilint[*,4]
				1 : BEGIN
					eindex=ipt(u.dat.etime,u.stat.time)
					rmid=reform(u.dat.rmid[eindex,*])
					xplot=interpol(rmid,u.dat.epsin,ilint[*,4])
				END
				2 : BEGIN
					eindex=ipt(u.dat.etime,u.stat.time)
					rmid=reform(u.dat.rmid[eindex,*])
					xplot=(interpol(rmid,u.dat.epsin,ilint[*,4])-rmid[0])/(last(rmid)-rmid[0])
				END
			ENDCASE	
			oploterror,xplot[tmp],ilint[tmp,0],ilint[tmp,1],color=col[i],errcolor=col[i],psym=u.plot.sym[i]
		ENDIF
	ENDFOR
	oplot,xr,[0,0],linestyle=2,color=col[4]

	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.size[0],u.plot.size[1],0,0,0]
	ENDIF	
END

PRO comp_plot_temp,u
	IF u.stat.ps THEN BEGIN
		xsize=5.0
		ysize=5.0*900/700.0
		ls=0.65
		col=u.plot.pscol
	ENDIF ELSE BEGIN
		ls=1.0
		col=u.plot.col
		widget_control,u.id.draw3,get_value=draw_win
		window,0,xsize=u.plot.size[0],ysize=u.plot.size[1],/pixmap
	ENDELSE
	index=u.index
	xr=u.plot.xr
	yr=u.plot.ti
	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit=u.plot.xtit,ytit='T!li!n [keV]'
	FOR i=0,n(u.stat.line) DO BEGIN
		IF u.stat.inv[i] THEN BEGIN
			ipro=*u.dat.inv[i,index[i]]
			npts=n(ipro[*,8])+1
			IF ipro[npts-1,8] EQ ipro[npts/2-1,8] THEN nrho=npts/2 ELSE nrho=npts
			CASE u.plot.x OF
				0 : xplot=ipro[0:nrho-1,8]
				1 : BEGIN
					eindex=ipt(u.dat.etime,u.stat.time)
					rmid=reform(u.dat.rmid[eindex,*])
					xplot=interpol(rmid,u.dat.epsin,ipro[0:nrho-1,8])
				END
				2 : BEGIN
					eindex=ipt(u.dat.etime,u.stat.time)
					rmid=reform(u.dat.rmid[eindex,*])
					xplot=(interpol(rmid,u.dat.epsin,ipro[*,8])-rmid[0])/(last(rmid)-rmid[0])
				END
			ENDCASE	
			oploterror,xplot,ipro[0:nrho-1,3]-u.dat.inst[i],ipro[0:nrho-1,7],color=col[i],errcolor=col[i]
			IF npts NE nrho THEN oploterror,xplot,ipro[nrho:*,3],ipro[nrho:*,7],color=col[i],errcolor=col[i],linestyle=3.0
		ENDIF
		IF u.stat.lint[i] THEN BEGIN
			IF u.plot.sym[i] EQ 8 THEN makesym,u.plot.psym[i]
			ilint=*u.dat.lint[i,index[i]]
			tmp=where(finite(ilint[*,2]) EQ 1 AND finite(ilint[*,3]) EQ 1)
			CASE u.plot.x OF
				0 : xplot=ilint[*,4]
				1 : BEGIN
					eindex=ipt(u.dat.etime,u.stat.time)
					rmid=reform(u.dat.rmid[eindex,*])
					xplot=interpol(rmid,u.dat.epsin,ilint[*,4])
				END
				2 : BEGIN
					eindex=ipt(u.dat.etime,u.stat.time)
					rmid=reform(u.dat.rmid[eindex,*])
					xplot=(interpol(rmid,u.dat.epsin,ilint[*,4])-rmid[0])/(last(rmid)-rmid[0])
				END
			ENDCASE	
			oploterror,xplot[tmp],ilint[tmp,2]-u.dat.inst[i],ilint[tmp,3],color=col[i],errcolor=col[i],psym=u.plot.sym[i]
		ENDIF
	ENDFOR
	IF u.plot.cxrs AND u.stat.cxrs THEN BEGIN
		icxrs=ipt(u.stat.time,u.dat.cxrs_tor.time)
		CASE u.plot.x OF
			0 : BEGIN
				eindex=ipt(u.dat.etime,u.stat.time)
				rmid=reform(u.dat.rmid[eindex,*])
				xplot=interpol(u.dat.epsin,rmid,u.dat.cxrs_tor.rmid)
			END
			1 : BEGIN
				xplot=u.dat.cxrs_tor.rmid
			END
			2 : BEGIN
				eindex=ipt(u.dat.etime,u.stat.time)
				rmid=reform(u.dat.rmid[eindex,*])
				xplot=(u.dat.cxrs_tor.rmid-rmid[0])/(last(rmid)-rmid[0])
			END
                ENDCASE
		;oploterror,xplot,u.dat.cxrs_tor.itemp[*,icxrs]/1.0e3,u.dat.cxrs_tor.iterr[*,icxrs]/1.0e3,color=col[0],errcolor=col[0],psym=u.plot.sym[0]
		icxrs=ipt(u.stat.time,u.dat.cxrs_pol.time)
		CASE u.plot.x OF
			0 : BEGIN
				eindex=ipt(u.dat.etime,u.stat.time)
				rmid=reform(u.dat.rmid[eindex,*])
				xplot=interpol(u.dat.epsin,rmid,u.dat.cxrs_pol.rmid)
			END
			1 : BEGIN
				xplot=u.dat.cxrs_tor.rmid
			END
			2 : BEGIN
				eindex=ipt(u.dat.etime,u.stat.time)
				rmid=reform(u.dat.rmid[eindex,*])
				xplot=(u.dat.cxrs_pol.rmid-rmid[0])/(last(rmid)-rmid[0])
			END
                ENDCASE
		oploterror,xplot,u.dat.cxrs_pol.itemp[*,icxrs]/1.0e3,u.dat.cxrs_pol.iterr[*,icxrs]/1.0e3,color=col[0],errcolor=col[0],psym=u.plot.sym[0]

	ENDIF
	oplot,xr,[0,0],linestyle=2,color=col[4]

	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.size[0],u.plot.size[1],0,0,0]
	ENDIF
END

PRO comp_plot_legend,u
	IF u.stat.ps THEN BEGIN
		xsize=5.0
		ysize=5.0*900/700.0
		ls=0.65
		col=u.plot.pscol
	ENDIF ELSE BEGIN
		ls=1.15
		col=u.plot.col
		widget_control,u.id.draw4,get_value=draw_win
		window,0,xsize=u.plot.lsize[0],ysize=u.plot.lsize[1],/pixmap
	ENDELSE
	label=['Ar XVI - w', 'Ar XVI - x', 'Ar XVI - z', 'Ar XVII - Lya', 'Mo XXXIII - 4d' ,'Ar XVII - J','Ca XVIII - w', 'Ar XVII - Lya', 'Mo XXXIII - 4d', 'Ca XIX - Lya']
	cntr=0
	plot,[0],[0],xr=[0,1],/ysty,yr=[0,1],/xsty,tit='LEGEND',xticks=1,xtickname=replicate(' ',2+1),yticks=1,ytickname=replicate(' ',2+1),$
		pos=[0.05,0.05,0.95,0.89]
	xyouts,0.05,0.9,'LINE',chars=ls
	xyouts,0.35,0.9,'INV(0)',chars=ls
	xyouts,0.55,0.9,'INV(1)',chars=ls
	xyouts,0.75,0.9,'LINE-INT',chars=ls
	FOR i=0,n(label) DO BEGIN
		IF u.stat.inv[i] THEN BEGIN
			oplot,[0.35,0.45],(0.82-0.15*cntr)*[1,1],color=col[i]
			oplot,[0.55,0.65],(0.82-0.15*cntr)*[1,1],color=col[i],linestyle=3.0

		ENDIF
		IF u.stat.lint[i] THEN BEGIN
			IF u.plot.sym[i] EQ 8 THEN makesym,u.plot.psym[i]
			oplot,[0.83],(0.815-0.15*cntr)*[1],psym=u.plot.sym[i],color=col[i]
		END
		IF u.stat.inv[i] OR u.stat.lint[i] THEN BEGIN
			xyouts,0.05,0.80-0.15*cntr,label[i],color=col[i],chars=ls
			cntr+=1
		END
	ENDFOR
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.lsize[0],u.plot.lsize[1],0,0,0]
	ENDIF
END



PRO load_compare_data,u
	shot=u.shot
	tht=u.tht
	nlines=n(u.stat.line)+1
	IF u.stat.dat THEN BEGIN					;if LOAD_DATA already ran, clean out the heap
		cleanup_w_compare,u
		heap_gc
	ENDIF

	mdsopen,'dnb',shot
	vtor=mdsvalue('\DNB::TOP.MIT_CXRS.RESULTS.ACTIVE.TOR_OUT:VEL',quiet=quiet,status=status)
	IF status THEN BEGIN
		sigvt=mdsvalue('\DNB::TOP.MIT_CXRS.RESULTS.ACTIVE.TOR_OUT:VEL_SIGMA')
		titor=mdsvalue('\DNB::TOP.MIT_CXRS.RESULTS.ACTIVE.TOR_OUT:TEMP')
		sigti=mdsvalue('\DNB::TOP.MIT_CXRS.RESULTS.ACTIVE.TOR_OUT:TEMP_SIGMA')
		rmid=mdsvalue('dim_of(\DNB::TOP.MIT_CXRS.RESULTS.ACTIVE.TOR_OUT:VEL,0)')/100.0
		time=mdsvalue('dim_of(\DNB::TOP.MIT_CXRS.RESULTS.ACTIVE.TOR_OUT:VEL,1)')
		omega=vtor
		FOR i=0,n(time) DO omega[*,i]/=2.0*!pi*rmid
		omerr=sigvt
		FOR i=0,n(time) DO sigvt[*,i]/=2.0*!pi*rmid
		cxrs_tor={omega:omega,omerr:omerr,itemp:titor,iterr:sigti,rmid:rmid,time:time}
	ENDIF ELSE BEGIN
		cxrs_tor={omega:[0],omerr:[0],itemp:[0],iterr:[0],rmid:[0],time:[0]}
        ENDELSE
	vpol=mdsvalue('\DNB::TOP.MIT_CXRS.RESULTS.ACTIVE.POLOIDAL:VEL',quiet=quiet,status=status)
	IF status THEN BEGIN
		sigvp=mdsvalue('\DNB::TOP.MIT_CXRS.RESULTS.ACTIVE.POLOIDAL:VEL_SIGMA')
		tipol=mdsvalue('\DNB::TOP.MIT_CXRS.RESULTS.ACTIVE.POLOIDAL:TEMP')
		sigti=mdsvalue('\DNB::TOP.MIT_CXRS.RESULTS.ACTIVE.POLOIDAL:TEMP_SIGMA')
		rmid=mdsvalue('dim_of(\DNB::TOP.MIT_CXRS.RESULTS.ACTIVE.POLOIDAL:VEL,0)')/100.0
		time=mdsvalue('dim_of(\DNB::TOP.MIT_CXRS.RESULTS.ACTIVE.POLOIDAL:VEL,1)')
		cxrs_pol={vpol:vpol,vperr:sigvp,itemp:tipol,iterr:sigti,rmid:rmid,time:time}
		u.stat.cxrs=1
	ENDIF ELSE BEGIN
		cxrs_pol={vpol:[0],vperr:[0],itemp:[0],iterr:[0],rmid:[0],time:[0]}
		u.stat.cxrs=0
        ENDELSE
	mdsclose,'dnb',shot
	hirexsr_load_image,shot,2,image,rawtime,/noimage
	ntime=n(rawtime)+1
	u.stat.ntime=ntime
	inv=ptrarr(nlines,ntime,/allocate)
	lint=ptrarr(nlines,ntime,/allocate)
	tau=ptrarr(nlines,/allocate)
	c=2.998e8 			;speed of light
	e=1.602e-19			;conversion for eV -> J
	mconv=1.661e-27			;conversion for amu -> kg
	
	mdsopen,'analysis',shot
	rmid=mdsvalue('\efit_rmid')
	etime=mdsvalue('dim_of(\efit_rmid,0)')
	epsin=mdsvalue('dim_of(\efit_rmid,1)')
	mdsclose,'analysis',shot

	lintid=[u.id.lint0,u.id.lint1,u.id.lint2,u.id.lint3,u.id.lint4,u.id.lint5,u.id.lint6,u.id.lint7,u.id.lint8,u.id.lint9,u.id.lint10,$
		u.id.lint11,u.id.lint12,u.id.lint13]
	invid=lintid-1
	inst=fltarr(nlines)
	FOR i=0,n(u.stat.line) DO BEGIN
		IF u.stat.inv[i] EQ 1 THEN BEGIN
			hirexsr_load_profile,shot,u.stat.line[i],pro_arr,err_arr,rho_arr,itau,status=status,tht=tht,tinst=tinst
			inst[i]=tinst
			IF status THEN BEGIN
				inv[i,*]=hirexsr_profile_arr2ptr(pro_arr,err_arr,rho_arr)
				widget_control,u.id.message,set_value='LINE='+num2str(u.stat.line[i])+' INV Loaded: '+num2str(shot,1),/app
				*tau[i]=itau
				u.index[i]=ipt(itau,u.stat.time)
			ENDIF ELSE BEGIN
				u.stat.inv[i]=0
				widget_control,u.id.message,set_value='LINE='+num2str(u.stat.line[i])+' INV not available',/app
				widget_control,invid[i],set_button=0
			ENDELSE
		ENDIF ELSE FOR j=0,ntime-1 DO *inv[i,j]=-1
		IF u.stat.lint[i] EQ 1 THEN BEGIN
			hirexsr_load_momentptr,shot,u.stat.line[i],mom,itau,pos,tpos,lam_o,z,status=status,tht=tht
			IF status THEN BEGIN
				mass=read_atomic_mass(z)
				conv_factor=(lam_o/c)^2*(e*1.0e3/(mass*mconv))		;conversion factor for 
				FOR j=0,ntime-1 DO BEGIN
					imom=*mom[j]
					IF imom[0] NE -1 THEN BEGIN
						pindex=last(where(tpos LE itau[j]))
						jpos=pos[*,*,pindex]
						v=-1.0*(lam_o-imom[*,11])*c/lam_o/(2.0*!pi*jpos[2,*]*cos(jpos[3,*]))*1.0e-3
						verr=imom[*,14]*c/lam_o/(2.0*!pi*jpos[2,*]*cos(jpos[3,*]))*1.0e-3
						ti=(imom[*,12])^2/conv_factor
						tierr=2.0*imom[*,15]*sqrt(imom[*,12]^2)/conv_factor
						rhotang=imom[*,16]
						ilint=[[v],[verr],[ti],[tierr],[rhotang]]
						*lint[i,j]=ilint
					ENDIF ELSE *lint[i,j]=-1
				ENDFOR	
			
				*tau[i]=itau
				u.index[i]=ipt(itau,u.stat.time)
				heap_free,mom
				widget_control,u.id.message,set_value='LINE='+num2str(u.stat.line[i])+' LINT Loaded: '+num2str(shot,1),/app
			ENDIF ELSE BEGIN
				u.stat.lint[i]=0
				widget_control,u.id.message,set_value='LINE='+num2str(u.stat.line[i])+' LINT not available',/app
				widget_control,lintid[i],set_button=0
			ENDELSE
		ENDIF ELSE FOR j=0,ntime-1 DO *lint[i,j]=-1
	ENDFOR	
	dat={inv:inv,lint:lint,tau:tau,rmid:rmid,inst:inst,etime:etime,epsin:epsin,time:rawtime,cxrs_tor:cxrs_tor,cxrs_pol:cxrs_pol}
	u={id:u.id,shot:u.shot,tht:u.tht,index:u.index,stat:u.stat,plot:u.plot,dat:dat}
	u.stat.dat=1
	widget_control,u.id.base, set_uvalue=u
END



PRO cleanup_w_compare,u
	IF u.stat.dat THEN BEGIN
		heap_free,u.dat.inv
		heap_free,u.dat.lint
		heap_free,u.dat.tau
	ENDIF
END

PRO reset_xplot_buttons,u
	id=[u.id.rpsin,u.id.rrmid,u.id.rra]
	FOR i=0,n(id) DO widget_control,id[i],set_button=0
END

PRO comp_plot_text,u
	widget_control,u.id.xmin,set_value=num2str(u.plot.xr[0],dp=2)
	widget_control,u.id.xmax,set_value=num2str(u.plot.xr[1],dp=2)
	widget_control,u.id.emmin,set_value=num2str(u.plot.em[0],dp=2)
	widget_control,u.id.emmax,set_value=num2str(u.plot.em[1],dp=2)
	widget_control,u.id.vtmin,set_value=num2str(u.plot.vt[0],dp=2)
	widget_control,u.id.vtmax,set_value=num2str(u.plot.vt[1],dp=2)
	widget_control,u.id.timin,set_value=num2str(u.plot.ti[0],dp=2)
	widget_control,u.id.timax,set_value=num2str(u.plot.ti[1],dp=2)	
END

PRO comp_plot_all,u
	comp_plot_emiss,u
	comp_plot_rot,u
	comp_plot_temp,u	
END

PRO w_hirexsr_compare_event,event
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
					cleanup_w_compare,u
					heap_gc
					!except=1
				
				END
				"LOAD" : BEGIN
					load_compare_data,u
					comp_plot_text,u
					comp_plot_all,u	
					comp_plot_legend,u	
				END
				"PRINT" : BEGIN
					IF u.stat.dat THEN stop
				END
				"STOP" : BEGIN
					stop
				END
				"RPSIN" : BEGIN
					reset_xplot_buttons,u
					widget_control,u.id.rpsin,set_button=1
					u.plot.x=0
					u.plot.xtit=n2g('psi')+'!ln!n'
					u.plot.xr=[0.0,1.0]
					widget_control,u.id.xmin,set_value=num2str(u.plot.xr[0],dp=2)
					widget_control,u.id.xmax,set_value=num2str(u.plot.xr[1],dp=2)
					comp_plot_all,u	
				END
				"RRMID" : BEGIN
					reset_xplot_buttons,u
					widget_control,u.id.rrmid,set_button=1
					u.plot.x=1
					u.plot.xtit='RMID [m]'
					u.plot.xr=[0.67,0.91]
					widget_control,u.id.xmin,set_value=num2str(u.plot.xr[0],dp=2)
					widget_control,u.id.xmax,set_value=num2str(u.plot.xr[1],dp=2)
					comp_plot_all,u	
				END
				"RRA" : BEGIN
					reset_xplot_buttons,u
					widget_control,u.id.rra,set_button=1
					u.plot.x=2
					u.plot.xtit='r/a'
					u.plot.xr=[0.0,1.0]
					widget_control,u.id.xmin,set_value=num2str(u.plot.xr[0],dp=2)
					widget_control,u.id.xmax,set_value=num2str(u.plot.xr[1],dp=2)
					comp_plot_all,u	
				END
				"MT" : BEGIN
					IF ipt(u.dat.time,u.stat.time) NE 0 AND u.stat.dat THEN BEGIN
						newtime=u.dat.time[ipt(u.dat.time,u.stat.time)-1]
						chkindex=u.index
						index=u.index
						FOR i=0,n(u.stat.line) DO BEGIN
							IF u.stat.inv[i] OR u.stat.lint[i] THEN BEGIN
								itau=*u.dat.tau[i]
								index[i]=ipt(itau,newtime)
								tmp=where(itau GT 0)
								IF newtime GE max(itau[tmp]) THEN index[i]=n(tmp)
								IF newtime LE min(itau[tmp]) THEN index[i]=0
							ENDIF
						ENDFOR
						u.stat.time=newtime
						widget_control,u.id.time,set_value=num2str(u.stat.time,dp=2)
						u.index=index
						IF total(index-chkindex) NE 0 THEN comp_plot_all,u
						widget_control,u.id.t_slider,set_value=newtime*1.0e3
					ENDIF
				END
				"PT" : BEGIN
					IF ipt(u.dat.time,u.stat.time) NE u.stat.ntime-1 AND u.stat.dat THEN BEGIN
						newtime=u.dat.time[ipt(u.dat.time,u.stat.time)+1]
						chkindex=u.index
						index=u.index
						FOR i=0,n(u.stat.line) DO BEGIN
							IF u.stat.inv[i] OR u.stat.lint[i] THEN BEGIN
								itau=*u.dat.tau[i]
								index[i]=ipt(itau,newtime)
								tmp=where(itau GT 0)
								IF newtime GE max(itau[tmp]) THEN index[i]=n(tmp)
								IF newtime LE min(itau[tmp]) THEN index[i]=0
							ENDIF
						ENDFOR
						u.stat.time=newtime
						widget_control,u.id.time,set_value=num2str(u.stat.time,dp=2)
						u.index=index
						IF total(index-chkindex) NE 0 THEN comp_plot_all,u
						widget_control,u.id.t_slider,set_value=newtime*1.0e3
					ENDIF
				END
				"INV0" : IF event.select EQ 1 THEN u.stat.inv[0]=1 ELSE u.stat.inv[0]=0
				"LINT0" : IF event.select EQ 1 THEN u.stat.lint[0]=1 ELSE u.stat.lint[0]=0
				"INV1" : IF event.select EQ 1 THEN u.stat.inv[1]=1 ELSE u.stat.inv[1]=0
				"LINT1" : IF event.select EQ 1 THEN u.stat.lint[1]=1 ELSE u.stat.lint[1]=0
				"INV2" : IF event.select EQ 1 THEN u.stat.inv[2]=1 ELSE u.stat.inv[2]=0
				"LINT2" : IF event.select EQ 1 THEN u.stat.lint[2]=1 ELSE u.stat.lint[2]=0
				"INV3" : IF event.select EQ 1 THEN u.stat.inv[3]=1 ELSE u.stat.inv[3]=0
				"LINT3" : IF event.select EQ 1 THEN u.stat.lint[3]=1 ELSE u.stat.lint[3]=0
				"INV4" : IF event.select EQ 1 THEN u.stat.inv[4]=1 ELSE u.stat.inv[4]=0
				"LINT4" : IF event.select EQ 1 THEN u.stat.lint[4]=1 ELSE u.stat.lint[4]=0
				"INV5" : IF event.select EQ 1 THEN u.stat.inv[5]=1 ELSE u.stat.inv[5]=0
				"LINT5" : IF event.select EQ 1 THEN u.stat.lint[5]=1 ELSE u.stat.lint[5]=0
				"INV6" : IF event.select EQ 1 THEN u.stat.inv[6]=1 ELSE u.stat.inv[6]=0
				"LINT6" : IF event.select EQ 1 THEN u.stat.lint[6]=1 ELSE u.stat.lint[6]=0
				"INV7" : IF event.select EQ 1 THEN u.stat.inv[7]=1 ELSE u.stat.inv[7]=0
				"LINT7" : IF event.select EQ 1 THEN u.stat.lint[7]=1 ELSE u.stat.lint[7]=0
				"INV8" : IF event.select EQ 1 THEN u.stat.inv[8]=1 ELSE u.stat.inv[8]=0
				"LINT8" : IF event.select EQ 1 THEN u.stat.lint[8]=1 ELSE u.stat.lint[8]=0
				"INV9" : IF event.select EQ 1 THEN u.stat.inv[9]=1 ELSE u.stat.inv[9]=0
				"LINT9" : IF event.select EQ 1 THEN u.stat.lint[9]=1 ELSE u.stat.lint[9]=0
				"INV10" : IF event.select EQ 1 THEN u.stat.inv[10]=1 ELSE u.stat.inv[10]=0
				"LINT10" : IF event.select EQ 1 THEN u.stat.lint[10]=1 ELSE u.stat.lint[10]=0
				"INV11" : IF event.select EQ 1 THEN u.stat.inv[11]=1 ELSE u.stat.inv[11]=0
				"LINT11" : IF event.select EQ 1 THEN u.stat.lint[11]=1 ELSE u.stat.lint[11]=0
				"INV12" : IF event.select EQ 1 THEN u.stat.inv[12]=1 ELSE u.stat.inv[12]=0
				"LINT12" : IF event.select EQ 1 THEN u.stat.lint[12]=1 ELSE u.stat.lint[12]=0
				"INV13" : IF event.select EQ 1 THEN u.stat.inv[13]=1 ELSE u.stat.inv[13]=0
				"LINT13" : IF event.select EQ 1 THEN u.stat.lint[13]=1 ELSE u.stat.lint[13]=0
				"PRAT" : BEGIN
					IF event.select EQ 1 THEN u.plot.ratio=1 ELSE u.plot.ratio=0
					IF u.stat.dat THEN comp_plot_all,u	
                                END
				"PCXRS" : BEGIN
					IF event.select EQ 1 THEN u.plot.cxrs=1 ELSE u.plot.cxrs=0
					IF u.stat.dat THEN comp_plot_all,u
				END
				ELSE:
			ENDCASE
		END
		"WIDGET_SLIDER": BEGIN
			widget_control,event.id,get_value=slider
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF 
				'T_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						chkindex=u.index
						index=u.index
						FOR i=0,n(u.stat.line) DO BEGIN
							IF u.stat.inv[i] OR u.stat.lint[i] THEN BEGIN
								itau=*u.dat.tau[i]
								index[i]=ipt(itau,slider/1.0e3)
								tmp=where(itau GT 0)
								IF slider/1.0e3 GE max(itau[tmp]) THEN index[i]=n(tmp)
								IF slider/1.0e3 LE min(itau[tmp]) THEN index[i]=0
							ENDIF
						ENDFOR
						u.stat.time=slider/1.0e3
						widget_control,u.id.time,set_value=num2str(u.stat.time,dp=2)
						u.index=index
						IF total(index-chkindex) NE 0 THEN comp_plot_all,u	
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
			u.id.xmin : BEGIN
				widget_control,u.id.xmin,get_value=value
				u.plot.xr[0]=float(value)
				comp_plot_all,u
			END
			u.id.xmax : BEGIN
				widget_control,u.id.xmax,get_value=value
				u.plot.xr[1]=float(value)
				comp_plot_all,u
			END
			u.id.emmin : BEGIN
				widget_control,u.id.emmin,get_value=value
				u.plot.em[0]=float(value)
				comp_plot_emiss,u
			END
			u.id.emmax : BEGIN
				widget_control,u.id.emmax,get_value=value
				u.plot.em[1]=float(value)
				comp_plot_emiss,u
			END
			u.id.vtmin : BEGIN
				widget_control,u.id.vtmin,get_value=value
				u.plot.vt[0]=float(value)
				comp_plot_rot,u
			END
			u.id.vtmax : BEGIN
				widget_control,u.id.vtmax,get_value=value
				u.plot.vt[1]=float(value)
				comp_plot_rot,u
			END
			u.id.timin : BEGIN
				widget_control,u.id.timin,get_value=value
				u.plot.ti[0]=float(value)
				comp_plot_temp,u
			END
			u.id.timax : BEGIN
				widget_control,u.id.timax,get_value=value
				u.plot.ti[1]=float(value)
				comp_plot_temp,u
			END

			ELSE:
			ENDCASE
		END
		ELSE:
	ENDCASE
	IF button NE 'QUIT' THEN widget_control,event.top,set_uvalue=u	
END

;+
;NAME:
;	W_HIREXSR_COMPARE
;
;-

PRO w_hirexsr_compare,shot=shot,tht=tht
	user=logname()
	loadct,12,/silent
	mdim=get_screen_size()
	IF mdim[0] LT 1600 OR mdim[1] LT 1100 THEN base=widget_base(title='HIREXSR PROFILE COMPARISON',/row,tlb_size_events=1,/scroll,$
		x_scroll_size=mdim[0]*0.95,y_scroll_size=mdim[1]*0.85) ELSE base=widget_base(title='HIREXSR PROFILE COMPARISON',/row,tlb_size_events=1)
	A=widget_base(base,/column)
	B=widget_base(base,/column)

	dum = widget_label(A,value='PROFILES')
	A1=widget_base(A,frame=5)
	draw1=widget_draw(A1,xsize=600,ysize=338)		
	A2=widget_base(A,frame=5)
	draw2=widget_draw(A2,xsize=600,ysize=338)
	A3=widget_base(A,frame=5)
	draw3=widget_draw(A3,xsize=600,ysize=338)


	xsize=420
	dum = widget_label(B,value='PLOTTING OPTIONS')
	B1=widget_base(B,xsize=xsize,ysize=270,/column,/frame)
	dum = widget_label(B,value='DATA SELECTION')
	B2=widget_base(B,xsize=xsize,ysize=330,/row)
	B2a=widget_base(B2,xsize=xsize*0.95/2.0,ysize=330,/column,/frame)
	B2b=widget_base(B2,xsize=xsize*0.95/2.0,ysize=330,/column,/frame)

	dum = widget_label(B,value='INPUT/OUTPUT')
	B3=widget_base(B,xsize=xsize,ysize=160,/column,/frame)
	B4=widget_base(B,frame=5)
	draw4=widget_draw(B4,xsize=xsize,ysize=200)			;legend

	;plotting options
	B1p1=widget_base(B1,/row)
	dum=widget_label(B1p1,value='t: ')
	time=widget_text(B1p1,xsize=5,ysize=1,/edit)
	dum=widget_label(B1p1,value='[sec]')
	B1p1x=widget_base(B1p1,/nonexcl,/row)
	rpsin=widget_button(B1p1x,value='PSIN  ')
	rrmid=widget_button(B1p1x,value='RMID  ')
	rra=widget_button(B1p1x,value='r/a')
	B1p2=widget_base(B1,/row)
	t_slider=widget_slider(B1p2,xsize=350,min=0,max=2000,value=1000,/drag,/suppress)
	mt=widget_button(B1p2,value=' - ')
	pt=widget_button(B1p2,value=' + ')
	P1=widget_base(B1, /row)
	P1a=widget_base(P1,/row,/nonexcl)
	xauto=widget_button(P1a,value='AUTOSCALE   ')
	xmin=widget_text(P1,xsize=5,ysize=1,/edit)
	dum=widget_label(P1,value=' < r/a < ')
	xmax=widget_text(P1,xsize=5,ysize=1,/edit)
	P2=widget_base(B1, /row)
	P2a=widget_base(P2,/row,/nonexcl)
	emauto=widget_button(P2a,value='AUTOSCALE   ')
	emmin=widget_text(P2,xsize=5,ysize=1,/edit)
	dum=widget_label(P2,value=' < EMISS < ')
	emmax=widget_text(P2,xsize=5,ysize=1,/edit)
	dum=widget_label(P2,value=' [AU] ')
	P3=widget_base(B1, /row)
	P3a=widget_base(P3,/row,/nonexcl)
	vtauto=widget_button(P3a,value='AUTOSCALE   ')
	vtmin=widget_text(P3,xsize=5,ysize=1,/edit)
	dum=widget_label(P3,value=' < VTOR < ')
	vtmax=widget_text(P3,xsize=5,ysize=1,/edit)
	dum=widget_label(P3,value=' [kHz] ')	
	P41=widget_base(B1, /row)
	P41a=widget_base(P41,/row,/nonexcl)
	tiauto=widget_button(P41a,value='AUTOSCALE   ')
	timin=widget_text(P41,xsize=5,ysize=1,/edit)
	dum=widget_label(P41,value=' < TI < ')
	timax=widget_text(P41,xsize=5,ysize=1,/edit)
	dum=widget_label(P41,value=' [keV] ')	
	P5=widget_base(B1, /row)
	P5a=widget_base(P5,/row,/nonexcl)
	prat=widget_button(P5a,value='PLOT m=1/m=0 RATIO ')
	pcxrs=widget_button(P5a,value='PLOT CXRS DATA')
	
	;data selection
	B2p01=widget_base(B2a,/row)
	dum = widget_label(B2p01,value='LINE                 INV    LINT')
	B2p02=widget_base(B2a,/row)
	dum = widget_label(B2p02,value='--------------------------------')
	B2p1=widget_base(B2a,/row)
	dum=widget_label(B2p1,value='(0) Ar XVI - w      ')
	B2p1x=widget_base(B2p1,/row,/nonexcl)
	inv0=widget_button(B2p1x, value='  ')
	lint0=widget_button(B2p1x, value=' ')
	B2p2=widget_base(B2a,/row)
	dum=widget_label(B2p2,value='(1) Ar XVI - x      ')
	B2p2x=widget_base(B2p2,/row,/nonexcl)
	inv1=widget_button(B2p2x, value='  ')
	lint1=widget_button(B2p2x, value=' ')
	B2p3=widget_base(B2a,/row)
	dum=widget_label(B2p3,value='(2) Ar XVI - z      ')
	B2p3x=widget_base(B2p3,/row,/nonexcl)
	inv2=widget_button(B2p3x, value='  ')
	lint2=widget_button(B2p3x, value=' ')
	B2p4=widget_base(B2a,/row)
	dum=widget_label(B2p4,value='(3) Ar XVII - Lya   ')
	B2p4x=widget_base(B2p4,/row,/nonexcl)
	inv3=widget_button(B2p4x, value='  ')
	lint3=widget_button(B2p4x, value=' ')
	B2p5=widget_base(B2a,/row)
	dum=widget_label(B2p5,value='(4) Mo XXXIII - 4d  ')
	B2p5x=widget_base(B2p5,/row,/nonexcl)
	inv4=widget_button(B2p5x, value='  ')
	lint4=widget_button(B2p5x, value=' ')
	B2p6=widget_base(B2a,/row)
	dum=widget_label(B2p6,value='(5) Ar XVII - J     ')
	B2p6x=widget_base(B2p6,/row,/nonexcl)
	inv5=widget_button(B2p6x, value='  ')
	lint5=widget_button(B2p6x, value=' ')
	B2p7=widget_base(B2a,/row)
	dum=widget_label(B2p7,value='(6) Ca XVIII - w      ')
	B2p7x=widget_base(B2p7,/row,/nonexcl)
	inv6=widget_button(B2p7x, value='  ')
	lint6=widget_button(B2p7x, value=' ')

	B2q01=widget_base(B2b,/row)
	dum = widget_label(B2q01,value='LINE                 INV    LINT')
	B2q02=widget_base(B2b,/row)
	dum = widget_label(B2q02,value='--------------------------------')
	B2q1=widget_base(B2b,/row)
	dum=widget_label(B2q1,value='(7) Ar XVII - Lya    ')
	B2q1x=widget_base(B2q1,/row,/nonexcl)
	inv7=widget_button(B2q1x, value='  ')
	lint7=widget_button(B2q1x, value=' ')
	B2q2=widget_base(B2b,/row)
	dum=widget_label(B2q2,value='(8) Mo XXXIII - 4d   ')
	B2q2x=widget_base(B2q2,/row,/nonexcl)
	inv8=widget_button(B2q2x, value='  ')
	lint8=widget_button(B2q2x, value=' ')
	B2q3=widget_base(B2b,/row)
	dum=widget_label(B2q3,value='(9) Ca XIX - Lya     ')
	B2q3x=widget_base(B2q3,/row,/nonexcl)
	inv9=widget_button(B2q3x, value='  ')
	lint9=widget_button(B2q3x, value=' ')
	B2q4=widget_base(B2b,/row)
	dum=widget_label(B2q4,value='(10) open            ')
	B2q4x=widget_base(B2q4,/row,/nonexcl)
	inv10=widget_button(B2q4x, value='  ')
	lint10=widget_button(B2q4x, value=' ')
	B2q5=widget_base(B2b,/row)
	dum=widget_label(B2q5,value='(11) open            ')
	B2q5x=widget_base(B2q5,/row,/nonexcl)
	inv11=widget_button(B2q5x, value='  ')
	lint11=widget_button(B2q5x, value=' ')
	B2q6=widget_base(B2b,/row)
	dum=widget_label(B2q6,value='(12) open            ')
	B2q6x=widget_base(B2q6,/row,/nonexcl)
	inv12=widget_button(B2q6x, value='  ')
	lint12=widget_button(B2q6x, value=' ')
	B2q7=widget_base(B2b,/row)
	dum=widget_label(B2q7,value='(13) open            ')
	B2q7x=widget_base(B2q7,/row,/nonexcl)
	inv13=widget_button(B2q7x, value='  ')
	lint13=widget_button(B2q7x, value=' ')

	;input/output	
	B3p1=widget_base(B3,/row)
	dum = widget_label(B3p1,value='SHOT: ')
	shotid = widget_text(B3p1,xsize=10,ysize=1,/edit)
	dum = widget_label(B3p1,value='  THT:')
	thtid=widget_text(B3p1,xsize=2,ysize=1,/edit)
	dum = widget_label(B3p1,value=' ')
	load= widget_button(B3p1,value='LOAD')
	dum = widget_label(B3p1,value=' ')
	quit= widget_button(B3p1,value='QUIT')
	dum = widget_label(B3p1,value=' ')
	print= widget_button(B3p1,value='PS')
	dum = widget_label(B3p1,value=' ')
	stop= widget_button(B3p1,value='STOP')
	B3p2=widget_base(B3,/row)
	message = widget_text(B3p2,xsize=45,ysize=5,/scroll)

	id={base:base,draw1:draw1,draw2:draw2,draw3:draw3,draw4:draw4,$
		time:time,rpsin:rpsin,rrmid:rrmid,rra:rra,t_slider:t_slider,pt:pt,mt:mt,$
		xauto:xauto,emauto:emauto,vtauto:vtauto,tiauto:tiauto,$
		xmin:xmin,emmin:emmin,vtmin:vtmin,timin:timin,$
		xmax:xmax,emmax:emmax,vtmax:vtmax,timax:timax,$
		prat:prat,pcxrs:pcxrs,$
		inv0:inv0,inv1:inv1,inv2:inv2,inv3:inv3,inv4:inv4,inv5:inv5,inv6:inv6,inv7:inv7,inv8:inv8,inv9:inv9,inv10:inv10,inv11:inv11,inv12:inv12,inv13:inv13,$
		lint0:lint0,lint1:lint1,lint2:lint2,lint3:lint3,lint4:lint4,lint5:lint5,lint6:lint6,lint7:lint7,lint8:lint8,lint9:lint9,lint10:lint10,lint11:lint11,$
			lint12:lint12,lint13:lint13,$
		shotid:shotid,thtid:thtid,load:load,quit:quit,print:print,stop:stop,message:message}

	;0-w, 1-x, 2-z,3-lya,4-mo4d,

	IF NOT keyword_set(shot) THEN shot=1110316014
	IF NOT keyword_set(tht) THEN tht=0
	time=1.0
	stat={line:[0,1,2,3,4,5,6,7,8,9,10,11,12],ps:0,inv:intarr(14),lint:intarr(14),time:time,dat:0,ntime:125,cxrs:0}
	stat.inv[2]=1
	stat.lint[2]=1
	index=intarr(n(stat.line)+1)
	emplt=[0,5.0]
	vtplt=[-10.0,10.0]
	tiplt=[0.0,3.0]
	plot={sym:[4,5,6,7,8,2,7,7,8,2,4,5,6],psym:[0,0,0,0,9,0,0,0,8,0,0,0,0],col:[40,80,130,200,255,100,200,200,255,100,40,80,130],pscol:[30,80,130,200,0,200],size:[600,338],$
		lsize:[320,200],x:0,xtit:'',xr:[0.0,1.0],em:emplt,vt:vtplt,ti:tiplt,ratio:0,cxrs:0}
		
	u={id:id,shot:shot,tht:tht,index:index,stat:stat,plot:plot}

	widget_control,u.id.shotid,set_value=num2str(u.shot,1)
	widget_control,u.id.thtid,set_value=num2str(u.tht,1)
	lintid=[u.id.lint0,u.id.lint1,u.id.lint2,u.id.lint3,u.id.lint4,u.id.lint5,u.id.lint6,u.id.lint7,u.id.lint8,u.id.lint9,u.id.lint10,$
		u.id.lint11,u.id.lint12,u.id.lint13]
	invid=lintid-1
	FOR i=0,n(lintid) DO IF u.stat.lint[i] EQ 1 THEN widget_control,lintid[i],set_button=1
	FOR i=0,n(invid) DO IF u.stat.inv[i] EQ 1 THEN widget_control,invid[i],set_button=1
	CASE u.plot.x OF
		0 : BEGIN
			widget_control,u.id.rpsin,set_button=1
			u.plot.xtit=n2g('psi')+'!ln!n'
			u.plot.xr=[0.0,1.0]
		END
		1 : BEGIN
			widget_control,u.id.rpsin,set_button=1
			u.plot.xtit='RMID [m]'
			u.plot.xr=[0.68,0.91]
		END
		2 : BEGIN
			widget_control,u.id.rpsin,set_button=1
			u.plot.xtit='r/a'
			u.plot.xr=[0.0,1.0]
		END
	ENDCASE
	widget_control,u.id.time,set_value=num2str(u.stat.time,dp=2)
	widget_control,base,set_uvalue=u

	!except=0
	widget_control,base,/realize
	xmanager,'w_hirexsr_compare',base
END
