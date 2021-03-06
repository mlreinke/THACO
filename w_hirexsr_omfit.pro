PRO omfit_update_text,u
	omfit_write_atext,u
	IF NOT u.stat.nob THEN omfit_write_btext,u
	omfit_write_otext,u
	IF u.dat.chk[u.index] THEN widget_control,u.id.chk,set_value='1' ELSE widget_control,u.id.chk,set_value='0'
END

PRO omfit_write_otext,u
	ibs=*u.dat.bsom[u.index]
	widget_control,u.id.ommin,set_value=num2str(ibs.limits.(4),dp=2)
	widget_control,u.id.ommax,set_value=num2str(ibs.limits.(5),dp=2)
	widget_control,u.id.emax,set_value=num2str(ibs.limits.errmax,dp=2)
	widget_control,u.id.nsigma,set_value=num2str(ibs.config.nsigma,1)
	widget_control,u.id.order,set_value=num2str(ibs.config.order,1)
	widget_control,u.id.nknots,set_value=num2str(ibs.config.nknots,1)
	widget_control,u.id.ntrials,set_value=num2str(ibs.config.ntrials,1)
	widget_control,u.id.nrho,set_value=num2str(ibs.config.nrho,1)
	xv=make(ibs.fit.rho[0],last(ibs.fit.rho),ibs.config.nknots+ibs.config.order)
	bsnak,xv,ibs.config.order,xknot
	low=where(xknot EQ xknot[0])
	high=where(xknot EQ last(xknot))
	xknot=xknot[last(low):high[0]]
	xkstr=num2str(xknot[0],dp=2)
	FOR i=1,n(xknot) DO xkstr+=','+num2str(xknot[i],dp=2)
	widget_control,u.id.xknots,set_value=xkstr
	widget_control,u.id.fit_iter,set_button=0
	widget_control,u.id.fit_bs,set_button=0
	widget_control,u.id.fit_con,set_button=0	
	CASE ibs.config.fitcase OF
		0 : widget_control,u.id.fit_iter,set_button=1
		1 : widget_control,u.id.fit_bs,set_button=1
		2 : widget_control,u.id.fit_con,set_button=1
	ENDCASE	
END

PRO omfit_write_atext,u
	ibs=*u.dat.bsom[u.index]
	ainvb=[u.id.ainve_tree,u.id.ainve_pct,u.id.ainve_ev]
	FOR i=0,2 DO widget_control,ainvb[i],set_button=0
	IF ibs.error.a EQ 0 THEN BEGIN
		widget_control,u.id.ainve_tree,set_button=1
		widget_control,u.id.ainve,set_value='x'
		u.stat.estat[0]=0
        ENDIF
	IF ibs.error.a GT 0 THEN BEGIN
		widget_control,u.id.ainve_pct,set_button=1
		widget_control,u.id.ainve,set_value=num2str(ibs.error.a*100.0,dp=1)
		u.stat.estat[0]=1
   
        ENDIF	
	IF ibs.error.a LT 0 THEN BEGIN
		widget_control,u.id.ainve_ev,set_button=1
		widget_control,u.id.ainve,set_value=num2str(-1.0*ibs.error.a,dp=2)
		u.stat.estat[0]=2
        ENDIF
	widget_control,u.id.ainvl,set_value=num2str(ibs.limits.arho[0],dp=2)
	widget_control,u.id.ainvls,set_value=ibs.limits.arho[0]*1.0e3
	widget_control,u.id.ainvh,set_value=num2str(ibs.limits.arho[1],dp=2)
	widget_control,u.id.ainvhs,set_value=ibs.limits.arho[1]*1.0e3

	alintb=[u.id.alinte_tree,u.id.alinte_pct,u.id.alinte_ev]
	FOR i=0,2 DO widget_control,alintb[i],set_button=0
	IF ibs.error.al EQ 0 THEN BEGIN
		widget_control,u.id.alinte_tree,set_button=1
		widget_control,u.id.alinte,set_value='x'
		u.stat.estat[1]=0
        ENDIF
	IF ibs.error.al GT 0 THEN BEGIN
		widget_control,u.id.alinte_pct,set_button=1
		widget_control,u.id.alinte,set_value=num2str(ibs.error.al*100.0,dp=1)
		u.stat.estat[1]=1
        ENDIF	
	IF ibs.error.al LT 0 THEN BEGIN
		widget_control,u.id.alinte_ev,set_button=1
		widget_control,u.id.alinte,set_value=num2str(-1.0*ibs.error.al,dp=2)
		u.stat.estat[1]=2
        ENDIF
	widget_control,u.id.alintl,set_value=num2str(ibs.limits.alrho[0],dp=2)
	widget_control,u.id.alintls,set_value=ibs.limits.alrho[0]*1.0e3
	widget_control,u.id.alinth,set_value=num2str(ibs.limits.alrho[1],dp=2)
	widget_control,u.id.alinths,set_value=ibs.limits.alrho[1]*1.0e3

END

PRO omfit_write_btext,u
	ibs=*u.dat.bsom[u.index]
	binvb=[u.id.binve_tree,u.id.binve_pct,u.id.binve_ev]
	FOR i=0,2 DO widget_control,binvb[i],set_button=0
	IF ibs.error.b EQ 0 THEN BEGIN
		widget_control,u.id.binve_tree,set_button=1
		widget_control,u.id.binve,set_value='x'
		u.stat.estat[2]=0
        ENDIF
	IF ibs.error.b GT 0 THEN BEGIN
		widget_control,u.id.binve_pct,set_button=1
		widget_control,u.id.binve,set_value=num2str(ibs.error.b*100.0,dp=1)
		u.stat.estat[2]=1
        ENDIF	
	IF ibs.error.b LT 0 THEN BEGIN
		widget_control,u.id.binve_ev,set_button=1
		widget_control,u.id.binve,set_value=num2str(-1.0*ibs.error.b,dp=2)
 		u.stat.estat[2]=2
        ENDIF
	widget_control,u.id.binvl,set_value=num2str(ibs.limits.brho[0],dp=2)
	widget_control,u.id.binvls,set_value=ibs.limits.brho[0]*1.0e3
	widget_control,u.id.binvh,set_value=num2str(ibs.limits.brho[1],dp=2)
	widget_control,u.id.binvhs,set_value=ibs.limits.brho[1]*1.0e3

	blintb=[u.id.blinte_tree,u.id.blinte_pct,u.id.blinte_ev]
	FOR i=0,2 DO widget_control,blintb[i],set_button=0
	IF ibs.error.bl EQ 0 THEN BEGIN
		widget_control,u.id.blinte_tree,set_button=1
		widget_control,u.id.blinte,set_value='x'
 		u.stat.estat[3]=0
        ENDIF
	IF ibs.error.bl GT 0 THEN BEGIN
		widget_control,u.id.blinte_pct,set_button=1
		widget_control,u.id.blinte,set_value=num2str(ibs.error.bl*100.0,dp=1)
  		u.stat.estat[3]=1
       ENDIF	
	IF ibs.error.bl LT 0 THEN BEGIN
		widget_control,u.id.blinte_ev,set_button=1
		widget_control,u.id.blinte,set_value=num2str(-1.0*ibs.error.bl,dp=2)
  		u.stat.estat[3]=2
       ENDIF
	widget_control,u.id.blintl,set_value=num2str(ibs.limits.blrho[0],dp=2)
	widget_control,u.id.blintls,set_value=ibs.limits.blrho[0]*1.0e3
	widget_control,u.id.blinth,set_value=num2str(ibs.limits.blrho[1],dp=2)
	widget_control,u.id.blinths,set_value=ibs.limits.blrho[1]*1.0e3

END

PRO omfit_bsfit_index,u
	ibs=*u.dat.bsom[u.index]
	ia=*u.dat.adat[u.index]
	hirexsr_bsfit_adjerror,ia,ibs.error.a,ibs.error.al
	IF NOT u.stat.nob THEN BEGIN
		ib=*u.dat.bdat[u.index]
		hirexsr_bsfit_adjerror,ib,ibs.error.b,ibs.error.bl
        ENDIF ELSE ib=-1
	hirexsr_omfit_input,ia,ib,idat,limits=ibs.limits
	iterfit=0
	conft=0
	CASE ibs.config.fitcase OF
		0 : iterfit=1
		1 : 
		2 : conft=1
	ENDCASE	
	
	hirexsr_run_bsfit,idat,ibsom,conft=conft,iterfit=iterfit,order=ibs.config.order,nknots=ibs.config.nknots,nsigma=ibs.config.nsigma,$
		nrho=ibs.config.nrho,ntrials=ibs.config.ntrials,fitcase=fitcase
	 *u.dat.bsom[u.index]={fit:ibsom,dat:idat,limits:ibs.limits,error:ibs.error,config:ibs.config,shot:u.shot,tht:u.tht,time:ibs.time}	
END


PRO omfit_plot,u
	omfit_plot_profiles,u
	omfit_plot_ranges,u
END

PRO omfit_plot_ranges,u
	IF u.stat.ps THEN BEGIN
		xsize=7.0
		ysize=7.0*600/800.0
		ls=0.65
		col=u.plot.pscol
	ENDIF ELSE BEGIN
		ls=0.8
		col=u.plot.col
		widget_control,u.id.draw3,get_value=draw_win
		window,0,xsize=u.plot.rsize[0],ysize=u.plot.rsize[1],/pixmap
	ENDELSE
	makesym,9
	ibs=*u.dat.bsom[u.index]
	limits=[[ibs.limits.arho],[ibs.limits.alrho],[ibs.limits.brho],[ibs.limits.blrho]]
	sym=[4,5,6,8]
	chk=abs([1,1,u.stat.nob-1,u.stat.nob-1])

	plot,[0],[0],xr=[0,1],yr=[0,4.1],/ysty,xtit='r/a',ytit='TYPE'
	FOR i=0,3 DO BEGIN
		IF chk[i] THEN BEGIN
			x=[0,limits[0,i],limits[*,i],limits[1,i],1.0]
			y=[0,0,i+1,i+1,0,0]
			oplot,x,y,psym=-sym[i],color=col[i+1]
		ENDIF
        ENDFOR
	makesym,10

	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.rsize[0],u.plot.rsize[1],0,0,0]
        ENDIF
END

PRO omfit_plot_profiles,u
	IF u.stat.ps THEN BEGIN
		xsize=7.0
		ysize=7.0*500/800.0
		ls=1.15
		col=u.plot.pscol
		device, xsize=xsize, ysize=ysize, /inches
	ENDIF ELSE BEGIN
		ls=1.2
		col=u.plot.col
		widget_control,u.id.draw1,get_value=draw_win
		window,0,xsize=u.plot.psize[0],ysize=u.plot.psize[1],/pixmap
	ENDELSE
	
	ibs=*u.dat.bsom[u.index]
	IF u.plot.pauto[0] THEN BEGIN
		u.plot.x=[0,1]
		widget_control,u.id.xl,set_value=num2str(u.plot.x[0],dp=2)
		widget_control,u.id.xh,set_value=num2str(u.plot.x[1],dp=2)
	ENDIF
	IF u.plot.pauto[1] THEN BEGIN
		tmp=where(ibs.fit.good EQ 1)
		ymax=max(ibs.dat.om[tmp]+ibs.dat.err[tmp]) > max(ibs.fit.prof+ibs.fit.err) > 0
		ymin=min(ibs.dat.om[tmp]-ibs.dat.err[tmp]) < min(ibs.fit.prof-ibs.fit.err) < 0
		u.plot.yom=[ymin,ymax]*1.05
		widget_control,u.id.yoml,set_value=num2str(u.plot.yom[0],dp=2)
		widget_control,u.id.yomh,set_value=num2str(u.plot.yom[1],dp=2)
	ENDIF

	IF u.plot.pauto[2] THEN BEGIN
		tmp=where(ibs.fit.rho LT 0.8)
		ymin=min(ibs.fit.dprof[tmp]-ibs.fit.derr[tmp]) < 0
		ymax=max(ibs.fit.dprof[tmp]+ibs.fit.derr[tmp]) > 0
		u.plot.ydom=[ymin,ymax]*1.05
		widget_control,u.id.ydoml,set_value=num2str(u.plot.ydom[0],dp=2)
		widget_control,u.id.ydomh,set_value=num2str(u.plot.ydom[1],dp=2)
	ENDIF
	
	xr=u.plot.x
	yr=u.plot.yom
	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit=ibs.dat.rlab,ytit='Tor. Rotation [kHz]',chars=ls
	type=[1,2,3,4]
	IF u.plot.alt THEN BEGIN
		dthick=2.0
		fthick=1.0
        ENDIF ELSE BEGIN
		fthick=2.0
		dthick=1.0
        ENDELSE
	FOR i=0,n(type) DO BEGIN
		tmp=where(ibs.fit.good EQ 1 AND ibs.dat.type EQ type[i])
		IF tmp[0] NE -1 THEN oploterror,ibs.dat.rho[tmp],ibs.dat.om[tmp],ibs.dat.err[tmp],psym=3,color=col[i+1],errcolor=col[i+1],thick=dthick
		tmp=where(ibs.fit.good EQ 0 AND ibs.dat.type EQ type[i])
		IF tmp[0] NE -1 THEN oplot,ibs.dat.rho[tmp],ibs.dat.om[tmp],psym=7,color=col[i+1],thick=dthick
	ENDFOR
	oploterror,ibs.fit.rho,ibs.fit.prof,ibs.fit.err,color=col[0],errcolor=col[0],thick=fthick
	oplot,xr,[0,0],linestyle=1
	IF u.stat.ps THEN xyouts,xr[1]+0.025*(xr[1]-xr[0]),yr[0]+0.02*(yr[1]-yr[0]),num2str(u.shot,1)+'   THT='+num2str(u.tht,1)+'   t='+num2str(u.stat.time,dp=2),$
		orient=90,chars=0.5*ls
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.psize[0],u.plot.psize[1],0,0,0]
		widget_control,u.id.draw2,get_value=draw_win
		window,0,xsize=u.plot.psize[0],ysize=u.plot.psize[1],/pixmap
        ENDIF

	xr=u.plot.x
	yr=u.plot.ydom
	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit=ibs.dat.rlab,ytit='Rotation Gradient [kHz/(r/a)]',chars=ls
	oploterror,ibs.fit.rho,ibs.fit.dprof,ibs.fit.derr,color=col[0],errcolor=col[0],thick=fthick
	oplot,xr,[0,0],linestyle=1
	IF u.stat.ps THEN xyouts,xr[1]+0.025*(xr[1]-xr[0]),yr[0]+0.02*(yr[1]-yr[0]),num2str(u.shot,1)+'   THT='+num2str(u.tht,1)+'   t='+num2str(u.stat.time,dp=2),$
		orient=90,chars=0.5*ls
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.psize[0],u.plot.psize[1],0,0,0]
        ENDIF
END

PRO omfit_load,u
	IF u.stat.dat THEN cleanup_w_omfit,u
	isdat=0
	runomfit=0
	CASE u.stat.ltree OF 
		0 : BEGIN
			path=u.stat.lpath
			chk=is_file(path)
			IF NOT chk THEN BEGIN
				ask = dialog_message(path+' not found, run new analysis?', /question)
				IF ask EQ 'Yes' THEN runomfit=1
                        ENDIF ELSE BEGIN
				restore,path
				isdat=1
				widget_control,u.id.message,set_value='BSOM loaded from: '+path,/app
			ENDELSE
                END

		1 : BEGIN
			chk=hirexsr_is_omresults(u.shot,u.tht,filled=filled)
			IF NOT filled THEN BEGIN
				ask = dialog_message('BSOM data not found, run new analysis?', /question)
				IF ask EQ 'Yes' THEN runomfit=1
                        ENDIF ELSE BEGIN
				hirexsr_read_bsom,u.shot,a,b,c,d,e,tht=u.tht,/ptr,bsom=bsom
				isdat=1
				widget_control,u.id.message,set_value='BSOM loaded from RESULTS node for THT='+num2str(u.tht,1),/app
			ENDELSE	
                END
		ELSE : 
	ENDCASE	
	IF NOT isdat AND NOT runomfit THEN RETURN			;no data and user doesn't want to make new
	IF NOT isdat AND runomfit THEN BEGIN
		hirexsr_omfit,u.shot,omfit,rho,time,err,bsom,tht=u.tht,/ptr,nobranchb=u.stat.nob
		widget_control,u.id.message,set_value='new BSOM computed for '+num2str(u.shot,1)+', THT='+num2str(u.tht,1),/app
        ENDIF

	;at this point, either by loading file, loading from tree or running new the bsom ptrarr should be formed
 	;load the raw ADAT and BDAT data which can be used to reformulate the inputs

	hirexsr_omfit_load,u.shot,u.stat.line,adat,atime,tht=u.tht,/rho
	widget_control,u.id.message,set_value='ADAT loaded for '+num2str(u.shot,1)+', THT='+num2str(u.tht,1),/app
 	IF NOT u.stat.nob THEN BEGIN
		hirexsr_omfit_load,u.shot,u.stat.line+1,bdat,btime,tht=u.tht,/rho
		widget_control,u.id.message,set_value='BDAT loaded for '+num2str(u.shot,1)+', THT='+num2str(u.tht,1),/app
		IF total(atime-btime) NE 0 THEN BEGIN
			print, 'ERROR: Branch A and Branch B timebase differet'
			heap_free,bsom
			heap_free,adat
			heap_free,bdat
			RETURN
		ENDIF
        ENDIF ELSE BEGIN
		bdat=-1
		widget_control,u.id.message,set_value='BDAT loading skipped',/app
	ENDELSE
	time=atime

	u.index=ipt(time,u.stat.time)
	u.stat.nframe=n(time)+1
	IF u.index EQ -1 THEN u.index=0
	u.stat.time=time[u.index]
	widget_control,u.id.t_slider,set_value=u.stat.time*1.0e3
	widget_control,u.id.time,set_value=num2str(u.stat.time,dp=2)
	widget_control,u.id.iframe,set_value=num2str(u.index,1)
	widget_control,u.id.nframe,set_value=num2str(u.stat.nframe-1,1)
	

	dat={bsom:bsom,adat:adat,bdat:bdat,time:time,chk:intarr(u.stat.nframe)}
	u={id:u.id,shot:u.shot,tht:u.tht,index:u.index,stat:u.stat,plot:u.plot,dat:dat}
	u.stat.dat=1
	widget_control,u.id.base, set_uvalue=u
END

PRO omfit_save,u
	IF total(u.dat.chk) NE 0 THEN BEGIN
		widget_control,u.id.message,set_value='CANNOT SAVE: Resolve Chnages in RANGE/SETUP',/app	
		RETURN
	ENDIF
	CASE u.stat.stree OF

		0 : BEGIN 
			bsom=u.dat.bsom
			time=u.dat.time
			path=u.stat.spath
			chk=is_file(path)
			IF chk THEN BEGIN
				ask = dialog_message(path+' exists, overwrite?', /question)
				IF ask EQ 'No' THEN BEGIN
        	                	widget_control,u.id.message,set_value='BSOM write cancelled',/app
					RETURN
				ENDIF
                        ENDIF 
			save,bsom,time,filename=path
			widget_control,u.id.message,set_value='BSOM written to: '+path,/app
                END

		1 : BEGIN
			chk=0
			IF NOT hirexsr_is_omresults(u.shot,u.tht) THEN BEGIN
				ask = dialog_message('No RESULTS node for THT='+num2str(u.tht,1)+'  exists. Create new node?', /question)
				IF ask EQ 'Yes' THEN hirexsr_addnew_omresults,u.shot,tht=u.tht,chk=chk						
                        ENDIF ELSE chk=1
			IF NOT chk THEN BEGIN
				widget_control,u.id.message,set_value='BSOM write cancelled',/app
				RETURN
			ENDIF
			hirexsr_bsfit_ptr2arr,u.dat.bsom,fit,data,config,rlab=rlab
			hirexsr_write_bsom,u.shot,fit,data,inst,config,u.dat.time,rlab=rlab,tht=u.tht
  			widget_control,u.id.message,set_value='BSOM written to RESULTS node for THT='+num2str(u.tht,1),/app
              	END	
	ENDCASE
	
END

PRO omfit_setpath,u
	path='/home/'+logname()+'/fits/omfit_'+num2str(u.shot,1)+'_THT'+num2str(u.tht,1)+'.dat'
	u.stat.spath=path[0]
	u.stat.lpath=path[0]
	IF u.stat.stree THEN widget_control,u.id.stree,set_button=1
	parts=strsplit(u.stat.spath,'/',/extract)
	widget_control,u.id.spath,set_value=last(parts)
	IF u.stat.ltree THEN widget_control,u.id.ltree,set_button=1
	parts=strsplit(u.stat.lpath,'/',/extract)
	widget_control,u.id.lpath,set_value=last(parts)
END

PRO cleanup_w_omfit,u
	IF u.stat.dat THEN BEGIN
		heap_free,u.dat.adat
		heap_free,u.dat.bsom
		IF size(u.dat.bdat,/type) EQ 10 THEN heap_free,u.dat.bdat
	ENDIF
END

PRO w_hirexsr_omfit_event,event
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
					cleanup_w_omfit,u
					heap_gc
					!except=1
				
				END
				"LOAD" : BEGIN
					WIDGET_CONTROL,/hourglass
					omfit_load,u
					IF u.stat.dat THEN BEGIN
						omfit_update_text,u
						omfit_plot,u
					ENDIF
				END
				"SAVE" : BEGIN
					WIDGET_CONTROL,/hourglass
					IF u.stat.dat THEN omfit_save,u
				END

				"PRINT" : BEGIN
					IF u.stat.dat THEN BEGIN
						psplot
						u.stat.ps=1
						omfit_plot_profiles,u
						psc
						xwplot
						set_plot,'x'
						u.stat.ps=0
						widget_control,u.id.message,set_value='profiles plotted to idl.ps',/app
					ENDIF
				END
				"STOP" : BEGIN
					stop
				END
				"T_M" : BEGIN
					IF u.index EQ 0 THEN RETURN
					u.index=u.index-1
					u.stat.time=u.dat.time[u.index]
					widget_control,u.id.t_slider,set_value=u.stat.time*1.0e3
					widget_control,u.id.time,set_value=num2str(u.stat.time,dp=2)
					widget_control,u.id.iframe,set_value=num2str(u.index,1)
					omfit_update_text,u
					omfit_plot,u					
				END
				"T_P" : BEGIN
					IF u.index EQ u.stat.nframe-1 THEN RETURN
					u.index=u.index+1
					u.stat.time=u.dat.time[u.index]
					widget_control,u.id.t_slider,set_value=u.stat.time*1.0e3
					widget_control,u.id.time,set_value=num2str(u.stat.time,dp=2)
					widget_control,u.id.iframe,set_value=num2str(u.index,1)
					omfit_update_text,u
					omfit_plot,u	
                                END
				"OPTFRM" : BEGIN
					widget_control,u.id.iopt,get_value=i
					IF i[0] EQ 'x' THEN RETURN
					i=int(i[0])
					IF i LT 0 OR i GT u.stat.nframe-1 THEN RETURN
					ibs=*u.dat.bsom[u.index]	
					limits=ibs.limits
					config=ibs.config
					ibs=*u.dat.bsom[i]
					ibs.limits.(4)=limits.(4)
					ibs.limits.(5)=limits.(5)
					ibs.limits.errmax=limits.errmax
					ibs.config.order=config.order
					ibs.config.nknots=config.nknots
					ibs.config.ntrials=config.ntrials
					ibs.config.nsigma=config.nsigma
					ibs.config.nrho=config.nrho
					ibs.config.fitcase=config.fitcase
					*u.dat.bsom[i]=ibs
 				 	widget_control,u.id.message,set_value='Current Fitting Options Set for Frame '+num2str(i,1),/app
					widget_control,u.id.message,set_value='***MUST SPLINE FIT FRAME '+num2str(i,1)+' TO UPDATE***',/app
					u.dat.chk[i]=1
				END
				"OPTALL" : BEGIN
					ibs=*u.dat.bsom[u.index]	
					limits=ibs.limits
					config=ibs.config
					FOR i=0,u.stat.nframe-1 DO BEGIN
						ibs=*u.dat.bsom[i]
						ibs.limits.(4)=limits.(4)
						ibs.limits.(5)=limits.(5)
						ibs.limits.errmax=limits.errmax
						ibs.config.order=config.order
						ibs.config.nknots=config.nknots
						ibs.config.ntrials=config.ntrials
						ibs.config.nsigma=config.nsigma
						ibs.config.nrho=config.nrho
						ibs.config.fitcase=config.fitcase
						*u.dat.bsom[i]=ibs
                                        ENDFOR
					widget_control,u.id.message,set_value='Current Fitting Options Set for All Frames',/app
					widget_control,u.id.message,set_value='***MUST SPLINE FIT ALL TO UPDATE***',/app
					u.dat.chk[*]=1			
				END
				"RNGFRM" : BEGIN
					widget_control,u.id.irng,get_value=i
					IF i[0] EQ 'x' THEN RETURN
					i=int(i[0])
					IF i LT 0 OR i GT u.stat.nframe-1 THEN RETURN
					ibs=*u.dat.bsom[u.index]	
					limits=ibs.limits
					error=ibs.error
					ibs=*u.dat.bsom[i]
					ibs.limits.arho=limits.arho
					ibs.limits.alrho=limits.alrho
					ibs.limits.brho=limits.brho
					ibs.limits.blrho=limits.blrho
					ibs.error.a=error.a
					ibs.error.b=error.b
					ibs.error.al=error.al
					ibs.error.bl=error.bl
					*u.dat.bsom[i]=ibs
 					widget_control,u.id.message,set_value='Current Data Range Set for Frame '+num2str(i,1),/app
					widget_control,u.id.message,set_value='***MUST SPLINE FIT FRAME '+num2str(i,1)+' TO UPDATE***',/app					
					u.dat.chk[i]=1
                                END
				"RNGALL" : BEGIN
					ibs=*u.dat.bsom[u.index]	
					limits=ibs.limits
					error=ibs.error
					FOR i=0,u.stat.nframe-1 DO BEGIN
						ibs=*u.dat.bsom[i]
						ibs.limits.arho=limits.arho
						ibs.limits.alrho=limits.alrho
						ibs.limits.brho=limits.brho
						ibs.limits.blrho=limits.blrho
						ibs.error.a=error.a
						ibs.error.b=error.b
						ibs.error.al=error.al
						ibs.error.bl=error.bl
						*u.dat.bsom[i]=ibs
                                        ENDFOR
					widget_control,u.id.message,set_value='Current Data Range Set for All Frames',/app
					widget_control,u.id.message,set_value='***MUST SPLINE FIT ALL TO UPDATE***',/app
					u.dat.chk[*]=1
				END
				"BSFRM" : BEGIN
					omfit_bsfit_index,u
					u.dat.chk[u.index]=0
					omfit_update_text,u
					omfit_plot,u
					widget_control,u.id.message,set_value='BS Fit Current Frame Complete',/app
				END
				"BSALL" : BEGIN
					WIDGET_CONTROL,/hourglass
					j=u.index
					FOR i=0,u.stat.nframe-1 DO BEGIN
						u.index=i
						omfit_bsfit_index,u
						u.dat.chk[i]=0
                                        ENDFOR
					u.index=j
					omfit_update_text,u
					omfit_plot,u
					widget_control,u.id.message,set_value='BS Fit All Frames Complete',/app
				END
				"FIT_ITER" : BEGIN
					widget_control,u.id.fit_iter,set_button=0
					widget_control,u.id.fit_bs,set_button=0
					widget_control,u.id.fit_con,set_button=0	
					ibs=*u.dat.bsom[u.index]
					CASE ibs.config.fitcase OF 
						0 : 
						1 : widget_control,u.id.fit_bs,set_button=1
						2 : widget_control,u.id.fit_con,set_button=1				
					ENDCASE
                                END
				"FIT_BS" : BEGIN
					widget_control,u.id.fit_iter,set_button=0
					widget_control,u.id.fit_bs,set_button=1
					widget_control,u.id.fit_con,set_button=0
					ibs=*u.dat.bsom[u.index]
					IF ibs.config.fitcase NE 1 THEN BEGIN
						ibs.config.fitcase=1
						*u.dat.bsom[u.index]=ibs
						u.dat.chk[u.index]=1
					ENDIF
				END					
				"FIT_CON" : BEGIN
					widget_control,u.id.fit_iter,set_button=0
					widget_control,u.id.fit_bs,set_button=0
					widget_control,u.id.fit_con,set_button=1				
					ibs=*u.dat.bsom[u.index]
					IF ibs.config.fitcase NE 2 THEN BEGIN
						ibs.config.fitcase=2
						*u.dat.bsom[u.index]=ibs
						u.dat.chk[u.index]=1
					ENDIF
				END
				"XAUTO" : IF event.select EQ 1 THEN u.plot.pauto[0]=1 ELSE u.plot.pauto[0]=0
				"YOMAUTO" : IF event.select EQ 1 THEN u.plot.pauto[1]=1 ELSE u.plot.pauto[1]=0
				"YDOMAUTO" : IF event.select EQ 1 THEN u.plot.pauto[2]=1 ELSE u.plot.pauto[2]=0
				"STREE" : IF event.select EQ 1 THEN u.stat.stree=1 ELSE u.stat.stree=0
				"SFILE" : BEGIN
					out=dialog_pickfile(path='/home/'+logname()+'/fits/',filter='omfit*',/write)
					IF out NE '' THEN BEGIN
						u.stat.spath=out
						parts=strsplit(u.stat.spath,'/',/extract)
						widget_control,u.id.spath,set_value=last(parts)
					ENDIF
				END
				"LTREE" : IF event.select EQ 1 THEN u.stat.ltree=1 ELSE u.stat.ltree=0
				"LFILE" : BEGIN
					out=dialog_pickfile(path='/home/'+logname()+'/fits/',filter='omfit*',/read,/must_exist)
					IF out NE '' THEN BEGIN
						u.stat.lpath=out
						parts=strsplit(u.stat.lpath,'/',/extract)
						widget_control,u.id.lpath,set_value=last(parts)
					ENDIF
				END
				"AINVE_TREE" : BEGIN
					widget_control,u.id.ainve_pct,set_button=0
					widget_control,u.id.ainve_ev,set_button=0
					widget_control,u.id.ainve,set_value='x'
					ibs=*u.dat.bsom[u.index]
					ibs.error.a=0.0
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					u.stat.estat[0]=0
				END
				"AINVE_PCT" : BEGIN
					widget_control,u.id.ainve_tree,set_button=0
					widget_control,u.id.ainve_ev,set_button=0
					widget_control,u.id.ainve,set_value=num2str(u.stat.pctdef*100.0,dp=1)
					ibs=*u.dat.bsom[u.index]
					ibs.error.a=u.stat.pctdef
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					u.stat.estat[0]=1
				END
				"AINVE_EV" : BEGIN
					widget_control,u.id.ainve_tree,set_button=0
					widget_control,u.id.ainve_pct,set_button=0
					widget_control,u.id.ainve,set_value=num2str(-1.0*u.stat.evdef,dp=2)
					ibs=*u.dat.bsom[u.index]
					ibs.error.a=u.stat.evdef
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					u.stat.estat[0]=2
				END
				"ALINTE_TREE" : BEGIN
					widget_control,u.id.alinte_pct,set_button=0
					widget_control,u.id.alinte_ev,set_button=0
					widget_control,u.id.alinte,set_value='x'
					ibs=*u.dat.bsom[u.index]
					ibs.error.al=0.0
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					u.stat.estat[1]=0
				END
				"ALINTE_PCT" : BEGIN
					widget_control,u.id.alinte_tree,set_button=0
					widget_control,u.id.alinte_ev,set_button=0
					widget_control,u.id.alinte,set_value=num2str(u.stat.pctdef*100.0,dp=1)
					ibs=*u.dat.bsom[u.index]
					ibs.error.al=u.stat.pctdef
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					u.stat.estat[1]=1
				END
				"ALINTE_EV" : BEGIN
					widget_control,u.id.alinte_tree,set_button=0
					widget_control,u.id.alinte_pct,set_button=0
					widget_control,u.id.alinte,set_value=num2str(-1.0*u.stat.evdef,dp=2)
					ibs=*u.dat.bsom[u.index]
					ibs.error.al=u.stat.evdef
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					u.stat.estat[1]=2
				END
				"BINVE_TREE" : BEGIN
					widget_control,u.id.binve_pct,set_button=0
					widget_control,u.id.binve_ev,set_button=0
					widget_control,u.id.binve,set_value='x'
					ibs=*u.dat.bsom[u.index]
					ibs.error.b=0.0
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					u.stat.estat[2]=0
				END
				"BINVE_PCT" : BEGIN
					widget_control,u.id.binve_tree,set_button=0
					widget_control,u.id.binve_ev,set_button=0
					widget_control,u.id.binve,set_value=num2str(u.stat.pctdef*100.0,dp=1)
					ibs=*u.dat.bsom[u.index]
					ibs.error.b=u.stat.pctdef
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					u.stat.estat[2]=1
				END
				"BINVE_EV" : BEGIN
					widget_control,u.id.binve_tree,set_button=0
					widget_control,u.id.binve_pct,set_button=0
					widget_control,u.id.binve,set_value=num2str(-1.0*u.stat.evdef,dp=2)
					ibs=*u.dat.bsom[u.index]
					ibs.error.b=u.stat.evdef
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					u.stat.estat[2]=2
				END
				"BLINTE_TREE" : BEGIN
					widget_control,u.id.blinte_pct,set_button=0
					widget_control,u.id.blinte_ev,set_button=0
					widget_control,u.id.blinte,set_value='x'
					ibs=*u.dat.bsom[u.index]
					ibs.error.bl=0.0
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					u.stat.estat[3]=0
				END
				"BLINTE_PCT" : BEGIN
					widget_control,u.id.blinte_tree,set_button=0
					widget_control,u.id.blinte_ev,set_button=0
					widget_control,u.id.blinte,set_value=num2str(u.stat.pctdef*100.0,dp=1)
					ibs=*u.dat.bsom[u.index]
					ibs.error.bl=u.stat.pctdef
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					u.stat.estat[3]=1
				END
				"BLINTE_EV" : BEGIN
					widget_control,u.id.blinte_tree,set_button=0
					widget_control,u.id.blinte_pct,set_button=0
					widget_control,u.id.blinte,set_value=num2str(-1.0*u.stat.evdef,dp=2)
					ibs=*u.dat.bsom[u.index]
					ibs.error.bl=u.stat.evdef
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					u.stat.estat[3]=2
				END
				'NOB' : IF event.select EQ 1 THEN u.stat.nob=1 ELSE u.stat.nob=0
				'ALTPLOT' : BEGIN
					IF event.select EQ 1 THEN  u.plot.alt=1 ELSE u.plot.alt=0
					omfit_plot,u
                                END
				ELSE :
			ENDCASE
		END
		"WIDGET_SLIDER": BEGIN
			widget_control,event.id,get_value=slider
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF 
				'T_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						chkindex=u.index
						index=ipt(u.dat.time,slider/1.0e3)
						IF index EQ -1 THEN BEGIN
							IF slider/1.0e3 GE max(u.dat.time) THEN index=u.stat.nframe-1
							IF slider/1.0e3 LE min(u.dat.time) THEN index=0
							widget_control,u.id.t_slider,set_value=u.dat.time[index]*1.0e3
                                                ENDIF
						
						IF total(index-chkindex) NE 0 THEN BEGIN
							u.stat.time=u.dat.time[index]
							widget_control,u.id.time,set_value=num2str(u.stat.time,dp=2)
							u.index=index
							widget_control,u.id.iframe,set_value=num2str(u.index,1)
							omfit_update_text,u
							omfit_plot,u	
						ENDIF
					ENDIF
				END
				'AINVLS' : BEGIN
					IF u.stat.dat THEN BEGIN
						ibs=*u.dat.bsom[u.index]
						ibs.limits.arho[0]=slider/1.0e3
						*u.dat.bsom[u.index]=ibs
						omfit_bsfit_index,u
						omfit_write_atext,u
						omfit_plot,u
					ENDIF
				END
				'AINVHS' : BEGIN
					IF u.stat.dat THEN BEGIN
						ibs=*u.dat.bsom[u.index]
						ibs.limits.arho[1]=slider/1.0e3
						*u.dat.bsom[u.index]=ibs
						omfit_bsfit_index,u
						omfit_write_atext,u
						omfit_plot,u
					ENDIF
				END
				'ALINTLS' : BEGIN
					IF u.stat.dat THEN BEGIN
						ibs=*u.dat.bsom[u.index]
						ibs.limits.alrho[0]=slider/1.0e3
						*u.dat.bsom[u.index]=ibs
						omfit_bsfit_index,u
						omfit_write_atext,u
						omfit_plot,u
					ENDIF
				END
				'ALINTHS' : BEGIN
					IF u.stat.dat THEN BEGIN
						ibs=*u.dat.bsom[u.index]
						ibs.limits.alrho[1]=slider/1.0e3
						*u.dat.bsom[u.index]=ibs
						omfit_bsfit_index,u
						omfit_write_atext,u
						omfit_plot,u
					ENDIF
				END
				'BINVLS' : BEGIN
					IF u.stat.dat THEN BEGIN
						ibs=*u.dat.bsom[u.index]
						ibs.limits.brho[0]=slider/1.0e3
						*u.dat.bsom[u.index]=ibs
						omfit_bsfit_index,u
						omfit_write_btext,u
						omfit_plot,u
					ENDIF
				END
				'BINVHS' : BEGIN
					IF u.stat.dat THEN BEGIN
						ibs=*u.dat.bsom[u.index]
						ibs.limits.brho[1]=slider/1.0e3
						*u.dat.bsom[u.index]=ibs
						omfit_bsfit_index,u
						omfit_write_btext,u
						omfit_plot,u
					ENDIF
				END
				'BLINTLS' : BEGIN
					IF u.stat.dat THEN BEGIN
						ibs=*u.dat.bsom[u.index]
						ibs.limits.blrho[0]=slider/1.0e3
						*u.dat.bsom[u.index]=ibs
						omfit_bsfit_index,u
						omfit_write_btext,u
						omfit_plot,u
					ENDIF
				END
				'BLINTHS' : BEGIN
					IF u.stat.dat THEN BEGIN
						ibs=*u.dat.bsom[u.index]
						ibs.limits.blrho[1]=slider/1.0e3
						*u.dat.bsom[u.index]=ibs
						omfit_bsfit_index,u
						omfit_write_btext,u
						omfit_plot,u
					ENDIF
				END
	
			ELSE:
			ENDCASE
		END
		"WIDGET_TEXT_CH": BEGIN
			widget_control,event.id,get_value=text
			CASE event.id OF 
			u.id.shotid : BEGIN
				widget_control,u.id.shotid,get_value=shot
				u.shot=shot
				IF u.shot GT 1120821000 AND u.shot LT 1121002000 THEN u.stat.line=7 ELSE u.stat.line=2
				IF u.stat.line EQ 2 THEN u.stat.nob=0 ELSE u.stat.nob=1
				widget_control,u.id.nob,set_button=u.stat.nob
				omfit_setpath,u
			END
			u.id.thtid : BEGIN
				widget_control,u.id.thtid,get_value=tht
				IF hirexsr_is_analysis(u.shot,int(tht[0])) THEN BEGIN
					u.tht=int(tht[0])
				ENDIF ELSE BEGIN
					widget_control,u.id.message,set_value='invalid THT',/app
					widget_control,u.id.thtid,set_value=num2str(u.tht,1)
				ENDELSE
				omfit_setpath,u
			END
			u.id.time : BEGIN
				widget_control,u.id.time,get_value=time
				u.stat.time=time
				IF u.stat.dat THEN BEGIN
					u.index=ipt(u.dat.time,u.stat.time)
					IF u.index EQ -1 THEN u.index=0
					u.stat.time=u.dat.time[u.index]
					widget_control,u.id.t_slider,set_value=u.stat.time*1.0e3
					widget_control,u.id.time,set_value=num2str(u.stat.time,dp=2)
					widget_control,u.id.iframe,set_value=num2str(u.index,1)
					omfit_update_text,u
					omfit_plot,u
				ENDIF
			END
			u.id.iframe : BEGIN
				widget_control,u.id.iframe,get_value=index
				IF u.stat.dat THEN BEGIN
					IF index LT 0 OR index GT u.stat.nframe -1 THEN BEGIN
						widget_control,u.id.iframe,set_value=num2str(u.index,1)
						RETURN
					ENDIF
					u.index=index
					u.stat.time=u.dat.time[u.index]
					widget_control,u.id.t_slider,set_value=u.stat.time*1.0e3
					widget_control,u.id.time,set_value=num2str(u.stat.time,dp=2)
					widget_control,u.id.iframe,set_value=num2str(u.index,1)
					omfit_update_text,u
					omfit_plot,u
                                ENDIF ELSE widget_control,u.id.iframe,set_value=num2str(u.index,1)
			END
			u.id.yoml : BEGIN
				u.plot.yom[0]=float(text)
				IF u.stat.dat THEN omfit_plot_profiles,u
			END
			u.id.yomh : BEGIN
				u.plot.yom[1]=float(text)
				IF u.stat.dat THEN omfit_plot_profiles,u
			END
			u.id.ydoml : BEGIN
				u.plot.ydom[0]=float(text)
				IF u.stat.dat THEN omfit_plot_profiles,u
			END
			u.id.ydomh : BEGIN
				u.plot.ydom[1]=float(text)
				IF u.stat.dat THEN omfit_plot_profiles,u
                             END
			u.id.ainvl : BEGIN
				IF u.stat.dat THEN BEGIN
					ibs=*u.dat.bsom[u.index]
					ibs.limits.arho[0]=float(text[0])
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					widget_control,u.id.ainvls,set_value=float(text[0])*1.0e3
				ENDIF
			END
			u.id.ainvh : BEGIN
				IF u.stat.dat THEN BEGIN
					ibs=*u.dat.bsom[u.index]
					ibs.limits.arho[1]=float(text[0])
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					widget_control,u.id.ainvhs,set_value=float(text[0])*1.0e3
				ENDIF
                        END
			u.id.alintl : BEGIN
				IF u.stat.dat THEN BEGIN
					ibs=*u.dat.bsom[u.index]
					ibs.limits.alrho[0]=float(text[0])
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					widget_control,u.id.alintls,set_value=float(text[0])*1.0e3
				ENDIF
			END
			u.id.alinth : BEGIN
				IF u.stat.dat THEN BEGIN
					ibs=*u.dat.bsom[u.index]
					ibs.limits.alrho[1]=float(text[0])
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					widget_control,u.id.alinths,set_value=float(text[0])*1.0e3
				ENDIF
			END
			u.id.binvl : BEGIN
				IF u.stat.dat THEN BEGIN
					ibs=*u.dat.bsom[u.index]
					ibs.limits.brho[0]=float(text[0])
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					widget_control,u.id.binvls,set_value=float(text[0])*1.0e3
				ENDIF
			END
			u.id.binvh : BEGIN
				IF u.stat.dat THEN BEGIN
					ibs=*u.dat.bsom[u.index]
					ibs.limits.brho[1]=float(text[0])
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					widget_control,u.id.binvhs,set_value=float(text[0])*1.0e3
				ENDIF
                        END
			u.id.blintl : BEGIN
				IF u.stat.dat THEN BEGIN
					ibs=*u.dat.bsom[u.index]
					ibs.limits.blrho[0]=float(text[0])
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					widget_control,u.id.blintls,set_value=float(text[0])*1.0e3
				ENDIF
			END
			u.id.blinth : BEGIN
				IF u.stat.dat THEN BEGIN
					ibs=*u.dat.bsom[u.index]
					ibs.limits.blrho[1]=float(text[0])
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
					widget_control,u.id.blinths,set_value=float(text[0])*1.0e3
				ENDIF
			END
			u.id.ainve : BEGIN
				IF u.stat.estat[0] EQ 0 THEN BEGIN
					widget_control,u.id.ainve,set_value='x'
                               	ENDIF ELSE BEGIN
					err=float(text[0])
					ibs=*u.dat.bsom[u.index]
					IF u.stat.estat[0] EQ 1 THEN ibs.error.a=err/100.0
					IF u.stat.estat[0] EQ 2 THEN ibs.error.a=-1.0*err
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
				ENDELSE
			END
			u.id.alinte : BEGIN
				IF u.stat.estat[1] EQ 0 THEN BEGIN
					widget_control,u.id.alinte,set_value='x'
                               	ENDIF ELSE BEGIN
					err=float(text[0])
					ibs=*u.dat.bsom[u.index]
					IF u.stat.estat[1] EQ 1 THEN ibs.error.al=err/100.0
					IF u.stat.estat[1] EQ 2 THEN ibs.error.al=-1.0*err
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
				ENDELSE
			END
			u.id.binve : BEGIN
				IF u.stat.estat[2] EQ 0 THEN BEGIN
					widget_control,u.id.binve,set_value='x'
                               	ENDIF ELSE BEGIN
					err=float(text[0])
					ibs=*u.dat.bsom[u.index]
					IF u.stat.estat[2] EQ 1 THEN ibs.error.b=err/100.0
					IF u.stat.estat[2] EQ 2 THEN ibs.error.b=-1.0*err
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
				ENDELSE
			END
			u.id.blinte : BEGIN
				IF u.stat.estat[3] EQ 0 THEN BEGIN
					widget_control,u.id.alinte,set_value='x'
                               	ENDIF ELSE BEGIN
					err=float(text[0])
					ibs=*u.dat.bsom[u.index]
					IF u.stat.estat[3] EQ 1 THEN ibs.error.bl=err/100.0
					IF u.stat.estat[3] EQ 2 THEN ibs.error.bl=-1.0*err
					*u.dat.bsom[u.index]=ibs
					omfit_bsfit_index,u
					omfit_plot,u
				ENDELSE
			END
			u.id.ommin : BEGIN
				ibs=*u.dat.bsom[u.index]
				ibs.limits.(4)=float(text[0])
				*u.dat.bsom[u.index]=ibs
				omfit_bsfit_index,u
				omfit_plot,u				
			END
			u.id.ommax : BEGIN
				ibs=*u.dat.bsom[u.index]
				ibs.limits.(5)=float(text[0])
				*u.dat.bsom[u.index]=ibs
				omfit_bsfit_index,u
				omfit_plot,u				
			END
			u.id.emax : BEGIN
				ibs=*u.dat.bsom[u.index]
				ibs.limits.errmax=float(text[0])
				*u.dat.bsom[u.index]=ibs
				omfit_bsfit_index,u
				omfit_plot,u				
			END
			u.id.order : BEGIN
				IF int(text[0]) GT 20 THEN BEGIN
					text='20'
					widget_control,u.id.order,set_value=text	
				ENDIF
				IF int(text[0]) LE 1 THEN BEGIN
					text='2'
					widget_control,u.id.order,set_value=text	
				ENDIF				
				ibs=*u.dat.bsom[u.index]
				ibs.config.order=int(text[0])
				*u.dat.bsom[u.index]=ibs
				omfit_bsfit_index,u
				omfit_write_otext,u
				omfit_plot,u				
                        END
			u.id.nknots : BEGIN
				ibs=*u.dat.bsom[u.index]
				ibs.config.nknots=int(text[0])
				*u.dat.bsom[u.index]=ibs
				omfit_bsfit_index,u
				omfit_write_otext,u
				omfit_plot,u				
                        END
			u.id.ntrials : BEGIN
				ibs=*u.dat.bsom[u.index]
				ibs.config.ntrials=int(text[0])
				*u.dat.bsom[u.index]=ibs
				omfit_bsfit_index,u
				omfit_plot,u				
                        END
			u.id.nsigma : BEGIN
				ibs=*u.dat.bsom[u.index]
				ibs.config.nsigma=float(text[0])
				*u.dat.bsom[u.index]=ibs
				omfit_bsfit_index,u
				omfit_plot,u				
                        END
			u.id.nrho : BEGIN
				ibs=*u.dat.bsom[u.index]
				ibs.config.nrho=int(text[0])
				*u.dat.bsom[u.index]=ibs
				omfit_bsfit_index,u
				omfit_plot,u				
                        END
			ELSE:
			ENDCASE
		END
		ELSE:
	ENDCASE
	IF ename NE 'QUIT' THEN BEGIN
		IF u.stat.dat THEN widget_control,u.id.chk,set_value=num2str(u.dat.chk[u.index],1)
		widget_control,event.top,set_uvalue=u
	ENDIF	
END


PRO w_hirexsr_omfit,shot=shot,tht=tht,line=line
	font='-adobe-symbol-medium-r-normal--12-120-75-75-p-74-adobe-fontspecific'
	loadct,12,/silent
	mdim=get_screen_size()
	IF mdim[0] NE 1600 LT  mdim[1] NE 1000 THEN base=widget_base(title='HIREXSR Omega Fit',/row,tlb_size_events=1,/scroll,$
		x_scroll_size=mdim[0]*0.95,y_scroll_size=mdim[1]*0.85) ELSE base=widget_base(title='HIREXSR Omega Fit',/row,tlb_size_events=1)
	A=widget_base(base,/column,/frame)
	B=widget_base(base,/column,/frame)	
	C=widget_base(base,/column,/frame)

	dum = widget_label(A,value='PROFILES')
	A1=widget_base(A,frame=5)
	draw1=widget_draw(A1,xsize=700,ysize=415)		
	A2=widget_base(A,frame=5)
	draw2=widget_draw(A2,xsize=700,ysize=415)
	A3=widget_base(A,/row)
	dum=widget_label(A3,value='Time: ')
	time=widget_text(A3,xsize=4,ysize=1,/edit)
	t_slider=widget_slider(A3,xsize=350,min=0,max=2000,value=1000,/drag,/suppress)
	dum=widget_label(A3,value=' ')
	t_m=widget_button(A3,value=' - ')
	t_p=widget_button(A3,value=' + ')
	dum=widget_label(A3,value=' FRAME #')
	iframe=widget_text(A3,xsize=3,ysize=1,/edit)
	dum=widget_label(A3,value=' of ')
	nframe=widget_text(A3,xsize=3,ysize=1)
	
	dum = widget_label(B,value='DATA RANGE OPTIONS')
	draw3=widget_draw(B,xsize=309,ysize=296,/frame)

	;inverted branch A
	B1=widget_base(B,xsize=xsize,ysize=135,/column,/frame)
	dum = widget_label(B1,value='BRANCH A: INVERTED')
	B1p1=widget_base(B1,/row)
	ainvl=widget_text(B1p1,xsize=6,ysize=1,/edit)
	ainvls=widget_slider(B1p1,xsize=200,min=0.0,max=1000,value=200,/drag,/suppress)
	B1p2=widget_base(B1,/row)
	ainvh=widget_text(B1p2,xsize=6,ysize=1,/edit)
	ainvhs=widget_slider(B1p2,xsize=200,min=0.0,max=1000,value=900,/drag,/suppress)
	B1p3=widget_base(B1,/row)
	dum=widget_label(B1p3,value='Set Unc.')
	ainve=widget_text(B1p3,xsize=5,ysize=1,/edit)
	B1p3a=widget_base(B1p3,/row,/nonexcl)
	ainve_tree=widget_button(B1p3a, value=' TREE ')
	ainve_pct=widget_button(B1p3a, value=' % ')
	ainve_ev=widget_button(B1p3a, value=' [kHz] ')

	;line averaged branch A
	B2=widget_base(B,xsize=xsize,ysize=135,/column,/frame)
	dum = widget_label(B2,value='BRANCH A: LINE. AVE.')
	B2p1=widget_base(B2,/row)
	alintl=widget_text(B2p1,xsize=6,ysize=1,/edit)
	alintls=widget_slider(B2p1,xsize=200,min=0.0,max=1000,value=200,/drag,/suppress)
	B2p2=widget_base(B2,/row)
	alinth=widget_text(B2p2,xsize=6,ysize=1,/edit)
	alinths=widget_slider(B2p2,xsize=200,min=0.0,max=1000,value=900,/drag,/suppress)
	B2p3=widget_base(B2,/row)
	dum=widget_label(B2p3,value='Set Unc.')
	alinte=widget_text(B2p3,xsize=5,ysize=1,/edit)
	B2p3a=widget_base(B2p3,/row,/nonexcl)
	alinte_tree=widget_button(B2p3a,value=' TREE ')
	alinte_pct=widget_button(B2p3a, value=' % ')
	alinte_ev=widget_button(B2p3a, value=' [kHz] ')

	;inverted branch B
	B3=widget_base(B,xsize=xsize,ysize=135,/column,/frame)
	dum = widget_label(B3,value='BRANCH B: INVERTED')
	B3p1=widget_base(B3,/row)
	binvl=widget_text(B3p1,xsize=6,ysize=1,/edit)
	binvls=widget_slider(B3p1,xsize=200,min=0.0,max=1000,value=200,/drag,/suppress)
	B3p2=widget_base(B3,/row)
	binvh=widget_text(B3p2,xsize=6,ysize=1,/edit)
	binvhs=widget_slider(B3p2,xsize=200,min=0.0,max=1000,value=900,/drag,/suppress)
	B3p3=widget_base(B3,/row)
	dum=widget_label(B3p3,value='Set Unc.')
	binve=widget_text(B3p3,xsize=5,ysize=1,/edit)
	B3p3a=widget_base(B3p3,/row,/nonexcl)
	binve_tree=widget_button(B3p3a, value=' TREE ')
	binve_pct=widget_button(B3p3a, value=' % ')
	binve_ev=widget_button(B3p3a, value=' [kHz] ')

	;line averaged branch B
	B4=widget_base(B,xsize=xsize,ysize=135,/column,/frame)
	dum = widget_label(B4,value='BRANCH A: LINE. AVE.')
	B4p1=widget_base(B4,/row)
	blintl=widget_text(B4p1,xsize=6,ysize=1,/edit)
	blintls=widget_slider(B4p1,xsize=200,min=0.0,max=1000,value=200,/drag,/suppress)
	B4p2=widget_base(B4,/row)
	blinth=widget_text(B4p2,xsize=6,ysize=1,/edit)
	blinths=widget_slider(B4p2,xsize=200,min=0.0,max=1000,value=900,/drag,/suppress)
	B4p3=widget_base(B4,/row)
	dum=widget_label(B4p3,value='Set Unc.')
	blinte=widget_text(B4p3,xsize=5,ysize=1,/edit)
	B4p3a=widget_base(B4p3,/row,/nonexcl)
	blinte_tree=widget_button(B4p3a,value=' TREE ')
	blinte_pct=widget_button(B4p3a, value=' % ')
	blinte_ev=widget_button(B4p3a, value=' [kHz] ')	

	C1=widget_base(C,/column)
	C1p1=widget_base(C1,/row)
	dum = widget_label(C1p1,value='SHOT: ')
	shotid = widget_text(C1p1,xsize=10,ysize=1,/edit)
	dum = widget_label(C1p1,value=' THT:')
	thtid=widget_text(C1p1,xsize=2,ysize=1,/edit)
	load= widget_button(C1p1,value=' LOAD ')
	save= widget_button(C1p1,value=' SAVE ')
	quit= widget_button(C1p1,value=' QUIT ')
	C1p2=widget_base(C1,/row)	
	dum = widget_label(C1p2,value='                                ')
	print= widget_button(C1p2,value=' PS ')
	stop= widget_button(C1p2,value=' STOP ')
	C1p3=widget_base(C1,/row)
	dum = widget_label(C1p3,value='SAVE: ')
	C1p3a=widget_base(C1p3,/row,/nonexcl)
	stree=widget_button(C1p3a,value='TREE')
	C1p3b=widget_base(C1p3,/row)
	dum=widget_label(C1p3b,value='FILE:')
	spath=widget_text(C1p3b,xsize=26,ysize=1)
	sfile=widget_button(C1p3b,value=' ? ')
	C1p4=widget_base(C1,/row)
	dum = widget_label(C1p4,value='LOAD: ')
	C1p4a=widget_base(C1p4,/row,/nonexcl)
	ltree=widget_button(C1p4a,value='TREE')
	C1p4b=widget_base(C1p4,/row)
	dum=widget_label(C1p4b,value='FILE:')
	lpath=widget_text(C1p4b,xsize=26,ysize=1)
	C1p5=widget_base(C1,/row)
	lfile=widget_button(C1p4b,value=' ? ')
	defpath='         DEFAULT PATH:  /home/'+logname()+'/fits/'
	dum=widget_label(C1p5, value=defpath[0])
	

	C1p5=widget_base(C1,/row)
	message = widget_text(C1p5,xsize=50,ysize=4,/scroll,/wrap)

	C2=widget_base(C,/column,/frame)
	dum = widget_label(C2,value='FITTING OPTIONS')
	C2p1=widget_base(C2,/row)
	dum = widget_label(C2p1,value='FITTING ROUTINE: ')	
	C2p1a=widget_base(C2p1,/row,/nonexcl)
	fit_iter=widget_button(C2p1a,value='ITERFIT  ')
	fit_bs=widget_button(C2p1a,value='BSLSQ  ')
	fit_con=widget_button(C2p1a,value='CONFT  ')
	C2p2=widget_base(C2,/row)
	dum = widget_label(C2p2,value='ACCEPTANCE: ')
	ommin=widget_text(C2p2,xsize=5,ysize=1,/edit)
	dum = widget_label(C2p2,value=' < Om < ')
	ommax=widget_text(C2p2,xsize=5,ysize=1,/edit)
	dum = widget_label(C2p2,value=' [kHz]  ')
	C2p3=widget_base(C2,/row)
	dum = widget_label(C2p3,value='           ')
	dum = widget_label(C2p3,value=' UNC(Om) < ')
	emax=widget_text(C2p3,xsize=5,ysize=1,/edit)
	dum = widget_label(C2p3,value='    NSIGMA < ')
	nsigma=widget_text(C2p3,xsize=3,ysize=1,/edit)
	C2p3x=widget_base(C2,/row)
	dum = widget_label(C2p3x,value='           ')
	C2p3xa=widget_base(C2p3x,/row,/nonexcl)
	nob=widget_button(C2p3xa,value='REMOVE BRANCH B FROM FIT')
	C2p4=widget_base(C2,/row)
	dum = widget_label(C2p4,value='ORDER:  ')
	order=widget_text(C2p4,xsize=4,ysize=1,/edit)
	dum = widget_label(C2p4,value='  #RADIAL POINTS:  ')
	nrho=widget_text(C2p4,xsize=4,ysize=1,/edit)
	C2p5=widget_base(C2,/row)
	dum = widget_label(C2p5,value='#KNOTS: ')
	nknots=widget_text(C2p5,xsize=4,ysize=1,/edit)
	dum = widget_label(C2p5,value='  X=')
	xknots=widget_text(C2p5,xsize=30,ysize=1)
	C2p6=widget_base(C2,/row)
	dum = widget_label(C2p6,value='#TRIALS')
	ntrials=widget_text(C2p6,xsize=4,ysize=1,/edit)
	C2p6a=widget_base(C2p6,/row,/nonexcl)
	dotrials=widget_button(C2p6a,value='ENABLE ERROR ESTIMATION')

	C3=widget_base(C,/column,/frame)
	C3p1=widget_base(C3,/row)
	optall=widget_button(C3p1,value=' COPY FITTING SETUP TO ALL FRAMES ')	
	C3p2=widget_base(C3,/row)
	optfrm=widget_button(C3p2,value=' COPY FITTING SETUP TO FRAME # ')
	dum = widget_label(C3p2,value=' ')
	iopt=widget_text(C3p2,xsize=4,ysize=1,/edit)
	C3p3=widget_base(C3,/row)
	rngall=widget_button(C3p3,value=' COPY DATA RANGE TO ALL FRAMES ')	
	C3p4=widget_base(C3,/row)
	rngfrm=widget_button(C3p4,value=' COPY DATA RANGE TO FRAME # ')
	dum = widget_label(C3p4,value=' ')
	irng=widget_text(C3p4,xsize=4,ysize=1,/edit)
	C3p5=widget_base(C3,/row)
	bsall=widget_button(C3p5,value=' SPLINE FIT ALL FRAMES ')	
	bsfrm=widget_button(C3p5,value=' SPLINE FIT CURRENT FRAME ')
	chk=widget_text(C3p5,xsize=1,ysize=1)

	C5=widget_base(C,/column,/frame)
	dum = widget_label(C5,value='PROFILE PLOTTING RANGES')
	C5p1=widget_base(C5,/row)
	yoml=widget_text(C5p1,xsize=5,ysize=1,/edit)
	dum=widget_label(C5p1,value='  <  Om  <  ')
	yomh=widget_text(C5p1,xsize=5,ysize=1,/edit)
	dum=widget_label(C5p1,value=' [kHz]  ')
	C5p1a=widget_base(C5p1,/row,/nonexcl)
	yomauto=widget_button(C5p1a,value='AUTO')
	C5p2=widget_base(C5,/row)
	ydoml=widget_text(C5p2,xsize=5,ysize=1,/edit)
	dum=widget_label(C5p2,value='  <  dOm/drho  <  ')
	ydomh=widget_text(C5p2,xsize=5,ysize=1,/edit)
	dum=widget_label(C5p2,value=' [kHz/(r/a)]  ')
	C5p2a=widget_base(C5p2,/row,/nonexcl)
	ydomauto=widget_button(C5p2a,value='AUTO')
	C5p3=widget_base(C5,/row)
	xl=widget_text(C5p3,xsize=5,ysize=1,/edit)
	dum=widget_label(C5p3,value='  <  r/a  <  ')
	xh=widget_text(C5p3,xsize=5,ysize=1,/edit)
	C5p3a=widget_base(C5p3,/row,/nonexcl)
	xauto=widget_button(C5p3a,value='AUTO')
	C5p4=widget_base(C5,/row,/nonexcl)
	altplot=widget_button(C5p4,value='ALT. PLOTTING STYLE')


	id={base:base,draw1:draw1,draw2:draw2,draw3:draw3,$
		time:time,t_slider:t_slider,t_m:t_m,t_p:t_p,iframe:iframe,nframe:nframe,$
		ainvl:ainvl,ainvls:ainvls,ainvh:ainvh,ainvhs:ainvhs,ainve:ainve,ainve_tree:ainve_tree,ainve_pct:ainve_pct,ainve_ev:ainve_ev,$
		alintl:alintl,alintls:alintls,alinth:alinth,alinths:alinths,alinte:alinte,alinte_tree:alinte_tree,alinte_pct:alinte_pct,alinte_ev:alinte_ev,$
		binvl:binvl,binvls:binvls,binvh:binvh,binvhs:binvhs,binve:binve,binve_tree:binve_tree,binve_pct:binve_pct,binve_ev:binve_ev,$
		blintl:blintl,blintls:blintls,blinth:blinth,blinths:blinths,blinte:blinte,blinte_tree:blinte_tree,blinte_pct:blinte_pct,blinte_ev:blinte_ev,$
		shotid:shotid,thtid:thtid,load:load,save:save,quit:quit,stop:stop,print:print,message:message,$
		stree:stree,spath:spath,sfile:sfile,ltree:ltree,lpath:lpath,lfile:lfile,$
		fit_iter:fit_iter,fit_bs:fit_bs,fit_con:fit_con,$
		ommin:ommin,ommax:ommax,emax:emax,nsigma:nsigma,nob:nob,$
		order:order,nrho:nrho,nknots:nknots,xknots:xknots,ntrials:ntrials,dotrials:dotrials,$
		optall:optall,optfrm:optfrm,iopt:iopt,rngall:rngall,rngfrm:rngfrm,irng:irng,bsall:bsall,bsfrm:bsfrm,chk:chk,$
		yoml:yoml,yomh:yomh,yomauto:yomauto,ydoml:ydoml,ydomh:ydomh,ydomauto:ydomauto,xl:xl,xh:xh,xauto:xauto,altplot:altplot}
	

	IF NOT keyword_set(shot) THEN shot=1120815009
	IF NOT keyword_set(tht) THEN tht=0
	IF NOT keyword_set(line) THEN line=2
	IF line EQ 2 THEN nob=0 ELSE nob=1	

	time=1.0
	path='/home/'+logname()+'/fits/omfit_'+num2str(shot,1)+'_THT'+num2str(tht,1)+'.dat'
	stat={time:time,dat:0,ps:0,stree:0,ltree:1,spath:path[0],lpath:path[0],nob:nob,line:line,nframe:0,pctdef:0.05,evdef:-0.1,estat:[0,0,0,0]}
	plot={col:[50,80,120,150,200],pscol:[30,100,120,150,200],psize:[700,415],rsize:[309,296],pauto:[0,1,1],x:[0.0,1.0],yom:[-20,20.0],ydom:[-1.0,1.0],alt:0}

	u={id:id,shot:shot,tht:tht,index:0,stat:stat,plot:plot}
	widget_control,u.id.shotid,set_value=num2str(u.shot,1)
	widget_control,u.id.thtid,set_value=num2str(u.tht,1)
	widget_control,u.id.time,set_value=num2str(u.stat.time,dp=2)
	widget_control,u.id.xl,set_value=num2str(u.plot.x[0],dp=2)
	widget_control,u.id.xh,set_value=num2str(u.plot.x[1],dp=2)
	widget_control,u.id.xauto,set_button=u.plot.pauto[0]
	widget_control,u.id.yoml,set_value=num2str(u.plot.yom[0],dp=2)
	widget_control,u.id.yomh,set_value=num2str(u.plot.yom[1],dp=2)
	widget_control,u.id.yomauto,set_button=u.plot.pauto[1]
	widget_control,u.id.ydoml,set_value=num2str(u.plot.ydom[0],dp=2)
	widget_control,u.id.ydomh,set_value=num2str(u.plot.ydom[1],dp=2)
	widget_control,u.id.ydomauto,set_button=u.plot.pauto[2]
	widget_control,u.id.altplot,set_button=u.plot.alt
	IF u.stat.stree THEN widget_control,u.id.stree,set_button=1
	parts=strsplit(u.stat.spath,'/',/extract)
	widget_control,u.id.spath,set_value=last(parts)
	IF u.stat.ltree THEN widget_control,u.id.ltree,set_button=1
	parts=strsplit(u.stat.lpath,'/',/extract)
	widget_control,u.id.lpath,set_value=last(parts)
	widget_control,base,set_uvalue=u
	widget_control,u.id.iopt,set_value='x'
	widget_control,u.id.irng,set_value='x'
	widget_control,u.id.dotrials, set_button=1
	
	!except=0
	widget_control,base,/realize
	xmanager,'w_hirexsr_omfit',base
	
END
