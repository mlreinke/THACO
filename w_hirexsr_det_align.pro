FUNCTION is_dat,u
	taglist=tag_names(u)
	tmp=where(taglist EQ 'DAT')
	IF tmp[0] NE -1 THEN output=1 ELSE output=0
	RETURN,output	
END

PRO plot_all,u
	plot_xy,u
	plot_yz,u
	plot_xz,u
	plot_xizeta,u
END

PRO plot_xy,u
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw1,get_value=draw_win
		window,0,xsize=u.stat.xsize[0],ysize=u.stat.ysize[0],/pixmap
	ENDIF ELSE BEGIN
		xsize=5.0
		device,xsize=xsize,ysize=xsize*u.stat.ysize[0]/u.stat.xsize[0],/in
	ENDELSE
	origin=u.dat.plotpts[*,0]
	IF u.stat.zoom THEN BEGIN
		xr=origin[1]-u.stat.dxzoom*[-1.0,1.0]
		yr=origin[0]+u.stat.dxzoom*[-1.0,1.0]
	ENDIF ELSE BEGIn
		xr=[0,-1.0*u.stat.dxfull]
		yr=[0,u.stat.dxfull]
	ENDELSE
	rad=u.dat.info.m.rad
	l=sqrt(origin[0]^2+origin[1]^2)
	h=origin[2]
	d=sqrt(rad^2/4.0+l^2-rad*l*sin(u.dat.align[2]))
	phi=atan(h/d)
	nrow=100
	xrow=make(0.0,rad,nrow)
	row=sqrt((rad/2.0)^2-(xrow-rad/2.0)^2)
	yrow=-1.0*cos(phi)*row

	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit='y [m]',ytit='x [m]',/iso
	IF u.stat.sym[0] EQ 8 THEN makesym,u.stat.usym[0]
	oplot,[0,origin[1]],[0,origin[0]],color=u.stat.color[0]
	oplot,[u.dat.info.det.x0[1]],[u.dat.info.det.x0[0]],color=u.stat.color[0],psym=u.stat.sym[0]
	IF u.stat.sym[1] EQ 8 THEN makesym,u.stat.usym[1]
	oplot,[u.dat.info.det.x0[1],u.dat.info.det.x2[1]],[u.dat.info.det.x0[0],u.dat.info.det.x2[0]],color=u.stat.color[1]
	oplot,[u.dat.info.det.x2[1]],[u.dat.info.det.x2[0]],color=u.stat.color[1],psym=u.stat.sym[1]
	IF u.stat.sym[2] EQ 8 THEN makesym,u.stat.usym[2]
	oplot,[u.dat.info.det.x0[1],u.dat.info.det.x1[1]],[u.dat.info.det.x0[0],u.dat.info.det.x1[0]],color=u.stat.color[2]
	oplot,[u.dat.info.det.x1[1]],[u.dat.info.det.x1[0]],color=u.stat.color[2],psym=u.stat.sym[2]
	oplot,yrow,xrow,thick=2.0,color=70,linestyle=2.0
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.stat.xsize[0],u.stat.ysize[0],0,0,0]
	ENDIF
END

PRO plot_yz,u
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw2,get_value=draw_win
		window,0,xsize=u.stat.xsize[0],ysize=u.stat.ysize[0],/pixmap
	ENDIF ELSE BEGIN
		xsize=5.0
		device,xsize=xsize,ysize=xsize*u.stat.ysize[0]/u.stat.xsize[0],/in
	ENDELSE
	origin=[u.dat.align[0],u.dat.align[1]]
	IF u.stat.zoom THEN BEGIN
		xr=origin[0]+u.stat.dxzoom*[-1.0,1.0]
		yr=origin[1]+u.stat.dxzoom*[-1.0,1.0]
	ENDIF ELSE BEGIN
		xr=[0,u.stat.dxfull]
		yr=u.stat.dxfull*[-0.5,0.5]
	ENDELSE
	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit='L [m]',ytit='z [m]',/iso
	IF u.stat.sym[0] EQ 8 THEN makesym,u.stat.usym[0]
	oplot,[0,origin[0]],[0,origin[1]],color=u.stat.color[0]
	oplot,[origin[0]-u.dat.plotpts[1,1]],[origin[1]+u.dat.plotpts[2,1]],color=u.stat.color[0],psym=u.stat.sym[0]
	IF u.stat.sym[1] EQ 8 THEN makesym,u.stat.usym[1]
	oplot,origin[0]-[u.dat.plotpts[1,1],u.dat.plotpts[1,3]],origin[1]+[u.dat.plotpts[2,1],u.dat.plotpts[2,3]],color=u.stat.color[1]
	oplot,[origin[0]-u.dat.plotpts[1,3]],[origin[1]+u.dat.plotpts[2,3]],color=u.stat.color[1],psym=u.stat.sym[1]
	IF u.stat.sym[2] EQ 8 THEN makesym,u.stat.usym[2]
	oplot,origin[0]-[u.dat.plotpts[1,1],u.dat.plotpts[1,2]],origin[1]+[u.dat.plotpts[2,1],u.dat.plotpts[2,2]],color=u.stat.color[2]
	oplot,[origin[0]-u.dat.plotpts[1,2]],[origin[1]+u.dat.plotpts[2,2]],color=u.stat.color[2],psym=u.stat.sym[2]
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.stat.xsize[0],u.stat.ysize[0],0,0,0]
	ENDIF
END


PRO plot_xz,u
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw3,get_value=draw_win
		window,0,xsize=u.stat.xsize[0],ysize=u.stat.ysize[0],/pixmap
	ENDIF ELSE BEGIN
		xsize=5.0
		device,xsize=xsize,ysize=xsize*u.stat.ysize[0]/u.stat.xsize[0],/in
	ENDELSE
	origin=[0.0,u.dat.align[1]]
	IF u.stat.zoom THEN BEGIN
		xr=origin[0]+u.stat.dxzoom*[-1.0,1.0]
		yr=origin[1]+u.stat.dxzoom*[-1.0,1.0]
	ENDIF ELSE BEGIN
		xr=u.stat.dxfull*[-0.5,0.5]
		yr=u.stat.dxfull*[-0.5,0.5]
	ENDELSE
	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit='d!lperp!n [m]',ytit='z [m]',/iso
	IF u.stat.sym[0] EQ 8 THEN makesym,u.stat.usym[0]
	oplot,[0,origin[0]],[0,origin[1]],color=u.stat.color[0]
	oplot,[origin[0]-u.dat.plotpts[0,1]],[origin[1]+u.dat.plotpts[2,1]],color=u.stat.color[0],psym=u.stat.sym[0]
	IF u.stat.sym[1] EQ 8 THEN makesym,u.stat.usym[1]
	oplot,origin[0]-[u.dat.plotpts[0,1],u.dat.plotpts[0,3]],origin[1]+[u.dat.plotpts[2,1],u.dat.plotpts[2,3]],color=u.stat.color[1]
	oplot,[origin[0]-u.dat.plotpts[0,3]],[origin[1]+u.dat.plotpts[2,3]],color=u.stat.color[1],psym=u.stat.sym[1]
	IF u.stat.sym[2] EQ 8 THEN makesym,u.stat.usym[2]
	oplot,origin[0]-[u.dat.plotpts[0,1],u.dat.plotpts[0,2]],origin[1]+[u.dat.plotpts[2,1],u.dat.plotpts[2,2]],color=u.stat.color[2]
	oplot,[origin[0]-u.dat.plotpts[0,2]],[origin[1]+u.dat.plotpts[2,2]],color=u.stat.color[2],psym=u.stat.sym[2]
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.stat.xsize[0],u.stat.ysize[0],0,0,0]
	ENDIF
END

PRO plot_xizeta,u
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw4,get_value=draw_win
		window,0,xsize=u.stat.xsize[1],ysize=u.stat.ysize[1],/pixmap
	ENDIF ELSE BEGIN
		xsize=3.0
		device,xsize=xsize,ysize=xsize*u.stat.ysize[1]/u.stat.xsize[1],/in
	ENDELSE
	n_xi=u.dat.info.det.n_xi
	n_zeta=u.dat.info.det.n_zeta
	x=size(u.dat.peaks.fpeaks)
	zeta=findgen(n_zeta)
	plot,[0],[0],xr=[0,n_xi],yr=[0,n_zeta],/xsty,/ysty,pos=[0.15,0.05,0.975,0.95]
	FOR i=0,u.nlines-1 DO BEGIN
		tmp=where(u.dat.peaks.fpeaks[*,i] LT n_xi)
		oplot,u.dat.peaks.fpeaks[tmp,i],zeta[tmp],color=u.stat.color[0]
		oplot,u.dat.peaks.ipeaks[tmp,i],zeta[tmp],color=u.stat.color[3]
	ENDFOR
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.stat.xsize[1],u.stat.ysize[1],0,0,0]
	ENDIF
	makesym,u.stat.usym[0]
	pos=[0.05,0.05,0.975,0.95]
	dx=0.2

	i=3
	IF u.nlines GE 4 THEN BEGIN
		IF NOT u.stat.ps THEN BEGIN
			widget_control,u.id.draw5,get_value=draw_win
			window,0,xsize=u.stat.xsize[2],ysize=u.stat.ysize[2],/pixmap
		ENDIF ELSE BEGIN
			xsize=2.0
			device,xsize=xsize,ysize=xsize*u.stat.ysize[2]/u.stat.xsize[2],/in
		ENDELSE
		tmp=where(u.dat.peaks.fpeaks[*,i] LT n_xi-1)
		delta=u.dat.peaks.fpeaks[*,i]-u.dat.peaks.ipeaks[*,i]
		xr=median(abs(delta[tmp]))*[-1.5,1.5]
		IF max(xr) LT dx THEN xr=[-dx,dx]
		plot,[0],[0],xr=xr,yr=[0,n_zeta],/xsty,/ysty,pos=pos,tit=n2g('lambda')+'='+num2str(u.dat.lambda[i],dp=5)+' [Ang]'
		loadct,39,/silent
		polyfill,[-dx,dx,dx,-dx,-dx],[0,0,n_zeta,n_zeta,0],color=160,linestyle=3.0
		loadct,12,/silent
		oplot,delta,zeta,psym=8
		IF NOT u.stat.ps THEN BEGIN
			wset,draw_win
			device,copy=[0,0,u.stat.xsize[2],u.stat.ysize[2],0,0,0]
		ENDIF
	ENDIF

	i=2
	IF u.nlines GE 3 THEN BEGIN
		IF NOT u.stat.ps THEN BEGIN
			widget_control,u.id.draw6,get_value=draw_win
			window,0,xsize=u.stat.xsize[2],ysize=u.stat.ysize[2],/pixmap
		ENDIF ELSE BEGIN
			xsize=2.0
			device,xsize=xsize,ysize=xsize*u.stat.ysize[2]/u.stat.xsize[2],/in
		ENDELSE
		tmp=where(u.dat.peaks.fpeaks[*,i] LT n_xi-1)
		delta=u.dat.peaks.fpeaks[*,i]-u.dat.peaks.ipeaks[*,i]
		xr=median(abs(delta[tmp]))*[-1.5,1.5]
		IF max(xr) LT dx THEN xr=[-dx,dx]
		plot,[0],[0],xr=xr,yr=[0,n_zeta],/xsty,/ysty,pos=pos,tit=n2g('lambda')+'='+num2str(u.dat.lambda[i],dp=5)+' [Ang]'
		loadct,39,/silent
		polyfill,[-dx,dx,dx,-dx,-dx],[0,0,n_zeta,n_zeta,0],color=160,linestyle=3.0
		loadct,12,/silent
		oplot,delta,zeta,psym=8
		IF NOT u.stat.ps THEN BEGIN
			wset,draw_win
			device,copy=[0,0,u.stat.xsize[2],u.stat.ysize[2],0,0,0]
		ENDIF
	ENDIF

	i=1
	IF u.nlines GE 2 THEN BEGIN
		IF NOT u.stat.ps THEN BEGIN
			widget_control,u.id.draw7,get_value=draw_win
			window,0,xsize=u.stat.xsize[2],ysize=u.stat.ysize[2],/pixmap
		ENDIF ELSE BEGIN
			xsize=2.0
			device,xsize=xsize,ysize=xsize*u.stat.ysize[2]/u.stat.xsize[2],/in
		ENDELSE
		tmp=where(u.dat.peaks.fpeaks[*,i] LT n_xi-1)
		delta=u.dat.peaks.fpeaks[*,i]-u.dat.peaks.ipeaks[*,i]
		xr=median(abs(delta[tmp]))*[-1.5,1.5]
		IF max(xr) LT dx THEN xr=[-dx,dx]
		plot,[0],[0],xr=xr,yr=[0,n_zeta],/xsty,/ysty,pos=pos,tit=n2g('lambda')+'='+num2str(u.dat.lambda[i],dp=5)+' [Ang]'
		loadct,39,/silent
		polyfill,[-dx,dx,dx,-dx,-dx],[0,0,n_zeta,n_zeta,0],color=160,linestyle=3.0
		loadct,12,/silent
		oplot,delta,zeta,psym=8
		IF NOT u.stat.ps THEN BEGIN
			wset,draw_win
			device,copy=[0,0,u.stat.xsize[2],u.stat.ysize[2],0,0,0]
		ENDIF
	ENDIF

	i=0
	IF u.nlines GE 1 THEN BEGIN
		IF NOT u.stat.ps THEN BEGIN
			widget_control,u.id.draw8,get_value=draw_win
			window,0,xsize=u.stat.xsize[2],ysize=u.stat.ysize[2],/pixmap
		ENDIF ELSE BEGIN
			xsize=2.0
			device,xsize=xsize,ysize=xsize*u.stat.ysize[2]/u.stat.xsize[2],/in
		ENDELSE
		tmp=where(u.dat.peaks.fpeaks[*,i] LT n_xi-1)
		delta=u.dat.peaks.fpeaks[*,i]-u.dat.peaks.ipeaks[*,i]
		xr=median(abs(delta[tmp]))*[-1.5,1.5]
		IF max(xr) LT dx THEN xr=[-dx,dx]
		plot,[0],[0],xr=xr,yr=[0,n_zeta],/xsty,/ysty,pos=pos,tit=n2g('lambda')+'='+num2str(u.dat.lambda[i],dp=5)+' [Ang]'
		loadct,39,/silent
		polyfill,[-dx,dx,dx,-dx,-dx],[0,0,n_zeta,n_zeta,0],color=160,linestyle=3.0
		loadct,12,/silent
		oplot,delta,zeta,psym=8
		IF NOT u.stat.ps THEN BEGIN
			wset,draw_win
			device,copy=[0,0,u.stat.xsize[2],u.stat.ysize[2],0,0,0]
		ENDIF
	ENDIF

END

PRO calc_ipeaks,u
	pix=findgen(u.dat.info.det.n_zeta)
	IF u.stat.finite THEN BEGIN
		FOR i=0,u.nlines-1 DO BEGIN
			u.dat.peaks.ipeaks[*,i]=0
			iphi=u.dat.phi[*,i]
			FOR j=0,n(iphi) DO BEGIN
				param=genpos_spherical2quadcurve(u.dat.info,u.dat.lambda[i],phi=iphi[j])
				u.dat.peaks.ipeaks[*,i]+=ellipse_xpt(param,pix)
			ENDFOR
			u.dat.peaks.ipeaks[*,i]/=n(iphi)+1
		ENDFOR	
	ENDIF ELSE BEGIN
		FOR i=0,u.nlines-1 DO BEGIN
			param=genpos_spherical2quadcurve(u.dat.info,u.dat.lambda[i],phi=0.0)
			u.dat.peaks.ipeaks[*,i]=ellipse_xpt(param,pix)
		ENDFOR	
	ENDELSE
END

PRO fill_tree_table,u
	tree_table=fltarr(3,3)
	tree_table[*,0]=u.dat.infotree.det.x0*100.0
	tree_table[*,1]=u.dat.infotree.det.x1*100.0
	tree_table[*,2]=u.dat.infotree.det.x2*100.0
	widget_control,u.id.tree_table,set_value=tree_table	
END


PRO fill_curr_table,u
	tree_table=fltarr(3,3)
	tree_table[*,0]=u.dat.info.det.x0*100.0
	tree_table[*,1]=u.dat.info.det.x1*100.0
	tree_table[*,2]=u.dat.info.det.x2*100.0
	widget_control,u.id.curr_table,set_value=tree_table
	lam=u.dat.info.m.bragg.twod*sin(u.dat.align[2])
	widget_control,u.id.lampt,set_value=num2str(lam,dp=5)
		
END

PRO fill_align_text,u
	widget_control,u.id.l_text,set_value=num2str(u.dat.align[0],dp=4)
	widget_control,u.id.l_slider,set_value=int(u.dat.align[0]*1000)
	widget_control,u.id.h_text,set_value=num2str(u.dat.align[1],dp=4)
	widget_control,u.id.h_slider,set_value=int(u.dat.align[1]*1000)
	widget_control,u.id.q_text,set_value=num2str(u.dat.align[2]*180/!pi,dp=2)
	widget_control,u.id.q_slider,set_value=int(u.dat.align[2]*180/!pi*1000)
	widget_control,u.id.x_text,set_value=num2str(u.dat.align[3],dp=2)
	widget_control,u.id.x_slider,set_value=int(u.dat.align[3]*10)
	widget_control,u.id.z_text,set_value=num2str(u.dat.align[4],dp=2)
	widget_control,u.id.z_slider,set_value=int(u.dat.align[4]*10)
	widget_control,u.id.a_text,set_value=num2str(u.dat.align[5]*180/!pi,dp=2)
	widget_control,u.id.a_slider,set_value=int(u.dat.align[5]*180/!pi*100)
	widget_control,u.id.b_text,set_value=num2str(u.dat.align[6]*180/!pi,dp=2)
	widget_control,u.id.b_slider,set_value=int(u.dat.align[6]*180/!pi*100)
	widget_control,u.id.c_text,set_value=num2str(u.dat.align[7]*180/!pi,dp=2)
	widget_control,u.id.c_slider,set_value=int(u.dat.align[7]*180/!pi*100)	

END

PRO write_data,u
	shot=u.shot
	module=u.module
	hirexsr_write_detalign,shot,module,u.dat.align
	hirexsr_load_info2tree,shot,module,info=u.dat.info
	widget_control,u.id.message,set_value='WROTE ALIGN/INFO Data: '+num2str(shot,1)+' MOD: '+num2str(module,1),/app
	hirexsr_write_pos,shot,module	
	widget_control,u.id.message,set_value='WROTE POS Data: '+num2str(shot,1)+' MOD: '+num2str(module,1),/app
	hirexsr_write_etendue,shot,module	
	widget_control,u.id.message,set_value='WROTE U Data: '+num2str(shot,1)+' MOD: '+num2str(module,1),/app

END

PRO load_detalign_data,u
	shot=u.shot
	module=u.module
	hirexsr_load_wavelengths,-1,lam_o,z_o,label_o								;load wavelengths table
	hirexsr_load_calibfits,shot,module,times,image,label,nphot,offset,peaks,width,resid,break,double,ifit_f		;load multigaussian fits
	hirexsr_load_ellipse,shot,module,coefs,ifit_e,outl,lambda,mu,sigma,bad						;load elliptical fit coefficients	
	IF u.stat.cad AND module NE 2 AND module NE 4 THEN BEGIN
		hirexsr_constrained_info,shot,i1,i3
		IF module EQ 1 THEN info=i1 ELSE info=i3
		align=genpos_info2align(info,x=0.0,z=0.0)
	ENDIF ELSE BEGIN
		info=hirexsr_read_treeinfo(u.shot,u.module)
		hirexsr_load_detalign,shot,module,align,status=status
		IF NOT status THEN align=genpos_info2align(info,x=0.0,z=0.0)
	ENDELSE
	tmp=genpos_align2info(align,info,plotpts=plotpts)
	widget_control,u.id.message,set_value='Loaded CALIB/ELLPISE Data: '+num2str(shot,1)+' MOD: '+num2str(module,1),/app
	
	;determine size of peaks and fill array
	x=size(coefs)
	nlines=x[2]
	u.nlines=nlines
	x=size(peaks)
	nfits=x[1]
	ipeaks=fltarr(nfits,nlines)
	fpeaks=fltarr(nfits,nlines)
	pix=findgen(nfits)
	FOR i=0,nlines-1 DO BEGIN
		tmp=where(lam_o EQ lambda[i])
		fpeaks[*,i]=reverse(info.det.n_xi-1.0-peaks[*,where(label EQ label_o[tmp[0]])])		;converts from (i,j) of image to (xi,zeta)
		param=genpos_spherical2quadcurve(info,lambda[i])
		ipeaks[*,i]=ellipse_xpt(param,pix)
	ENDFOR
	peaks={fpeaks:fpeaks,ipeaks:ipeaks,itpeaks:ipeaks}

	dphi=atan(info.m.size[0]/(2*info.m.rad))
	aphi=[dphi/2.0,0.0,-dphi/2.0]	;for x,y,z using full crystal
	bphi=[0.4*dphi,0.2*dphi,0.0]	;for w using 3/5 of the c-side of the crystal
	phi=fltarr(3,4)
	phi[*,0]=bphi
	FOR i=1,3 DO phi[*,i]=aphi
	dat={image:image,peaks:peaks,lambda:lambda,infotree:info,info:info,align:align,aligntree:align,plotpts:plotpts,phi:phi}
	u.stat.dat=1
	u={id:u.id,shot:u.shot,module:u.module,nlines:u.nlines,font:u.font,stat:u.stat,dat:dat}
	widget_control,u.id.base,set_uvalue=u		
END


PRO w_hirexsr_det_align_event,event
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
					heap_gc
					!except=1
					makesym,10
				
				END
				"STOP" : BEGIN
					stop
				END 
				"PRINT" : BEGIN
					u.stat.ps=1
					u.stat.color=u.stat.pscolor
					psplot
					plot_all,u
					psc
					xwplot
					set_plot,'x'
					u.stat.ps=0
					u.stat.color=u.stat.xwcolor
				END
				"LOAD" : BEGIN
					WIDGET_CONTROL,/hourglass
					u.stat.color=u.stat.xwcolor
					load_detalign_data,u
					fill_tree_table,u
					fill_curr_table,u
					fill_align_text,u
					plot_all,u
				END
				"WRITE" : BEGIN
					WIDGET_CONTROL,/hourglass
					write_data,u
				END
				"CADBUT" : BEGIN
					IF event.select EQ 1 THEN u.stat.cad=1 ELSE u.stat.cad=0
				END
				"ZOOMBUT" : BEGIN
					IF event.select EQ 1 THEN u.stat.zoom=1 ELSE u.stat.zoom=0
					IF u.stat.dat THEN plot_all,u
				END
				ELSE:
			ENDCASE
		END
		"WIDGET_SLIDER": BEGIN
			widget_control,event.id,get_value=slider
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF
				"L_SLIDER": BEGIN
					u.dat.align[0]=slider/1.0e3
					widget_control,u.id.l_text,set_value=num2str(u.dat.align[0],dp=4)
					u.dat.info=genpos_align2info(u.dat.align,u.dat.info,plotpts=plotpts)
					u.dat.plotpts=plotpts
					fill_curr_table,u
					calc_ipeaks,u
					plot_all,u
				END
				"H_SLIDER": BEGIN
					u.dat.align[1]=slider/1.0e3
					widget_control,u.id.h_text,set_value=num2str(u.dat.align[1],dp=4)
					u.dat.info=genpos_align2info(u.dat.align,u.dat.info,plotpts=plotpts)
					u.dat.plotpts=plotpts
					fill_curr_table,u
					calc_ipeaks,u
					plot_all,u
				END
				"Q_SLIDER": BEGIN
					u.dat.align[2]=slider/1000.0*!pi/180.0
					widget_control,u.id.q_text,set_value=num2str(slider/1000.0,dp=4)
					u.dat.info=genpos_align2info(u.dat.align,u.dat.info,plotpts=plotpts)
					u.dat.plotpts=plotpts
					fill_curr_table,u
					calc_ipeaks,u
					plot_all,u
				END
				"X_SLIDER": BEGIN
					u.dat.align[3]=slider/10.0
					widget_control,u.id.x_text,set_value=num2str(u.dat.align[3],dp=4)
					u.dat.info=genpos_align2info(u.dat.align,u.dat.info,plotpts=plotpts)
					u.dat.plotpts=plotpts
					fill_curr_table,u
					calc_ipeaks,u
					plot_all,u
				END
				"Z_SLIDER": BEGIN
					u.dat.align[4]=slider/10.0
					widget_control,u.id.z_text,set_value=num2str(u.dat.align[4],dp=4)
					u.dat.info=genpos_align2info(u.dat.align,u.dat.info,plotpts=plotpts)
					u.dat.plotpts=plotpts
					fill_curr_table,u
					calc_ipeaks,u
					plot_all,u
				END
				"A_SLIDER": BEGIN
					u.dat.align[5]=slider/100.0*!pi/180.0
					widget_control,u.id.a_text,set_value=num2str(slider/100.0,dp=4)
					u.dat.info=genpos_align2info(u.dat.align,u.dat.info,plotpts=plotpts)
					u.dat.plotpts=plotpts
					fill_curr_table,u
					calc_ipeaks,u
					plot_all,u
				END
				"B_SLIDER": BEGIN
					u.dat.align[6]=slider/100.0*!pi/180.0
					widget_control,u.id.b_text,set_value=num2str(slider/100.0,dp=4)
					u.dat.info=genpos_align2info(u.dat.align,u.dat.info,plotpts=plotpts)
					u.dat.plotpts=plotpts
					fill_curr_table,u
					calc_ipeaks,u
					plot_all,u
				END
				"C_SLIDER": BEGIN
					u.dat.align[7]=slider/100.0*!pi/180.0
					widget_control,u.id.c_text,set_value=num2str(slider/100.0,dp=4)
					u.dat.info=genpos_align2info(u.dat.align,u.dat.info,plotpts=plotpts)
					u.dat.plotpts=plotpts
					calc_ipeaks,u
					fill_curr_table,u
					plot_all,u
				END
				ELSE:
			ENDCASE
		END
   		"WIDGET_TEXT_CH": BEGIN
			widget_control,event.id,get_value=text
			CASE event.id OF 
				u.id.shotid : BEGIN
					u.shot=text
				END
				u.id.module : BEGIN
					u.module=text
				END
				u.id.dxzoom : BEGIN
					u.stat.dxzoom=text
					plot_all,u
				END
				u.id.dxfull : BEGIN
					u.stat.dxfull=text
					plot_all,u
				END
				u.id.l_text : BEGIN
					u.dat.align[0]=text
					widget_control,u.id.l_slider,set_value=int(text*1.0e3)
					u.dat.info=genpos_align2info(u.dat.align,u.dat.info,plotpts=plotpts)
					u.dat.plotpts=plotpts
					calc_ipeaks,u
					fill_curr_table,u
					plot_all,u	
				END
				u.id.h_text : BEGIN
					u.dat.align[1]=text
					widget_control,u.id.h_slider,set_value=int(text*1.0e3)
					u.dat.info=genpos_align2info(u.dat.align,u.dat.info,plotpts=plotpts)
					u.dat.plotpts=plotpts
					calc_ipeaks,u
					fill_curr_table,u
					plot_all,u	
				END
				u.id.q_text : BEGIN
					u.dat.align[2]=text*!pi/180.0
					widget_control,u.id.q_slider,set_value=int(text*1.0e3)
					u.dat.info=genpos_align2info(u.dat.align,u.dat.info,plotpts=plotpts)
					u.dat.plotpts=plotpts
					calc_ipeaks,u
					fill_curr_table,u
					plot_all,u	
				END
				u.id.x_text : BEGIN
					u.dat.align[3]=text
					widget_control,u.id.x_slider,set_value=int(text*1.0e1)
					u.dat.info=genpos_align2info(u.dat.align,u.dat.info,plotpts=plotpts)
					u.dat.plotpts=plotpts
					calc_ipeaks,u
					fill_curr_table,u
					plot_all,u	
				END
				u.id.z_text : BEGIN
					u.dat.align[4]=text
					widget_control,u.id.z_slider,set_value=int(text*1.0e1)
					u.dat.info=genpos_align2info(u.dat.align,u.dat.info,plotpts=plotpts)
					u.dat.plotpts=plotpts
					calc_ipeaks,u
					fill_curr_table,u
					plot_all,u	
				END
				u.id.a_text : BEGIN
					u.dat.align[5]=text*!pi/180.0
					widget_control,u.id.a_slider,set_value=int(text*1.0e2)
					u.dat.info=genpos_align2info(u.dat.align,u.dat.info,plotpts=plotpts)
					u.dat.plotpts=plotpts
					calc_ipeaks,u
					fill_curr_table,u
					plot_all,u	
				END
				u.id.b_text : BEGIN
					u.dat.align[6]=text*!pi/180.0
					widget_control,u.id.b_slider,set_value=int(text*1.0e2)
					u.dat.info=genpos_align2info(u.dat.align,u.dat.info,plotpts=plotpts)
					u.dat.plotpts=plotpts
					calc_ipeaks,u
					fill_curr_table,u
					plot_all,u	
				END
				u.id.c_text : BEGIN
					u.dat.align[7]=text*!pi/180.0
					widget_control,u.id.c_slider,set_value=int(text*1.0e2)
					u.dat.info=genpos_align2info(u.dat.align,u.dat.info,plotpts=plotpts)
					u.dat.plotpts=plotpts
					calc_ipeaks,u
					fill_curr_table,u
					plot_all,u	
				END
				ELSE: 
			ENDCASE
		END
		"WIDGET_TABLE_CH" : BEGIN
			CASE event.id OF 
				u.id.tree_table : BEGIN
				END
				u.id.curr_table : BEGIN
				END
				ELSE :
			ENDCASE
		END
		ELSE:
	ENDCASE
	IF button NE 'QUIT' THEN widget_control,event.top,set_uvalue=u		
END

;+
;NAME:
;	W_HIREXSR_DET_ALIGN
;
;-

PRO w_hirexsr_det_align,shot=shot
	font='-adobe-symbol-medium-r-normal--12-120-75-75-p-74-adobe-fontspecific'
	user=logname()
	loadct,12,/silent
	det=[195,487]
	mdim=get_screen_size()
	IF mdim[0] NE 1600 AND mdim[1] NE 1200 THEN base=widget_base(title='HIREXSR DETECTOR ALIGNMENT',/col,tlb_size_events=1,/scroll,$
		x_scroll_size=mdim[0]*0.95,y_scroll_size=mdim[1]*0.85) ELSE base=widget_base(title='HIREXSR DETECTOR ALIGNMENT',/col,tlb_size_events=1)
	R1=widget_base(base,/row)
	R2=widget_base(base,/row)
	A=widget_base(R2,/column)
	C=widget_base(R2,/column)
	D=widget_base(R2,/column)

	size=510
	DR1=widget_base(R1,frame=5)
	draw1=widget_draw(DR1,xsize=size,ysize=size)
	DR2=widget_base(R1,frame=5)
	draw2=widget_draw(DR2,xsize=size,ysize=size)
	DR3=widget_base(R1,frame=5)
	draw3=widget_draw(DR3,xsize=size,ysize=size)

	dum=widget_label(A,value='SETUP/OPTIONS')
	A3=widget_base(A,frame=2,xsize=400,ysize=160,/column)
	A3p1=widget_base(A3,/row)
	dum = widget_label(A3p1,value='SHOT: ')
	shotid = widget_text(A3p1,xsize=10,ysize=1,/edit)
	dum = widget_label(A3p1,value=' MOD: ')
	module = widget_text(A3p1,xsize=2,ysize=1,/edit)
	dum = widget_label(A3p1,value=' ')
	load= widget_button(A3p1,value='LOAD')
	dum = widget_label(A3p1,value='')
	quit= widget_button(A3p1,value='QUIT')
	dum = widget_label(A3p1,value='')
	print= widget_button(A3p1,value='PRINT')
	dum = widget_label(A3p1,value='')
	stop= widget_button(A3p1,value='STOP')
	A3p2=widget_base(A3,/row)
	message = widget_text(A3p2,xsize=50,ysize=5,/scroll)

	dum = widget_label(A,value='DETECTOR POSITION')
	B1=widget_base(A,frame=2,xsize=420,ysize=330,/column)
	B1p1=widget_base(B1,/row)
	dum=widget_label(B1p1,value=' L ')
	xsize=275
	l_slider=widget_slider(B1p1,xsize=xsize,min=0,max=2000,value=0,/drag,/suppress)
	l_text=widget_text(B1p1,xsize=8,ysize=1,/edit)
	dum=widget_label(B1p1,value='[m]')
	B1p2=widget_base(B1,/row)
	dum=widget_label(B1p2,value=' H ')
	h_slider=widget_slider(B1p2,xsize=xsize,min=-1000,max=1000,value=0,/drag,/suppress)
	h_text=widget_text(B1p2,xsize=8,ysize=1,/edit)
	dum=widget_label(B1p2,value='[m]')
	B1p3=widget_base(B1,/row)
	dum=widget_label(B1p3,value=' q  ',font=font)
	q_slider=widget_slider(B1p3,xsize=xsize,min=0,max=90000,value=0,/drag,/suppress)
	q_text=widget_text(B1p3,xsize=8,ysize=1,/edit)
	dum=widget_label(B1p3,value='[deg]')
	B1p4=widget_base(B1,/row)
	dum=widget_label(B1p4,value=' x ',font=font)
	x_slider=widget_slider(B1p4,xsize=xsize,min=0,max=int(det[0]*10.0),value=0,/drag,/suppress)
	x_text=widget_text(B1p4,xsize=8,ysize=1,/edit)	
	dum=widget_label(B1p4,value='[pix]')
	B1p5=widget_base(B1,/row)
	dum=widget_label(B1p5,value=' z ',font=font)
	z_slider=widget_slider(B1p5,xsize=xsize,min=0,max=int(det[1]*10.0),value=0,/drag,/suppress)
	z_text=widget_text(B1p5,xsize=8,ysize=1,/edit)
	dum=widget_label(B1p5,value='[pix]')
	B1p6=widget_base(B1,/row)
	dum=widget_label(B1p6,value=' a  ',font=font)
	a_slider=widget_slider(B1p6,xsize=xsize,min=-9000,max=9000,value=0,/drag,/suppress)
	a_text=widget_text(B1p6,xsize=8,ysize=1,/edit)	
	dum=widget_label(B1p6,value='[deg]')
	B1p7=widget_base(B1,/row)
	dum=widget_label(B1p7,value=' b  ',font=font)
	b_slider=widget_slider(B1p7,xsize=xsize,min=-9000,max=9000,value=0,/drag,/suppress)
	b_text=widget_text(B1p7,xsize=8,ysize=1,/edit)
	dum=widget_label(B1p7,value='[deg]')
	B1p8=widget_base(B1,/row)
	dum=widget_label(B1p8,value=' g  ',font=font)
	c_slider=widget_slider(B1p8,xsize=xsize,min=-9000,max=9000,value=0,/drag,/suppress)
	c_text=widget_text(B1p8,xsize=8,ysize=1,/edit)	
	dum=widget_label(B1p8,value='[deg]')

	dum = widget_label(C,value='TREE INFO FILE')
	clab=['x [cm]','y [cm]','z [cm]']
	rlab=['x0','x1','x2']
	size=275
	C1=widget_base(C,frame=2,xsize=size,ysize=150,/column)
	tree_table=widget_table(C1,xsize=3,ysize=3,row_lab=rlab,column_lab=clab,value=fltarr(3,3))
	C1a=widget_base(C1,/row)
	dum=widget_label(C1a,value=' x  ',font=font)
	tx=widget_text(C1a,xsize=8,ysize=1,/edit)
	dum=widget_label(C1a,value=' z  ',font=font)
	tz=widget_text(C1a,xsize=8,ysize=1,/edit)

	dum = widget_label(C,value='CURRENT INFO FILE')
	C2=widget_base(C,frame=2,xsize=size,ysize=265,/column)
	curr_table=widget_table(C2,xsize=3,ysize=3,/edit,row_lab=rlab,column_lab=clab,value=fltarr(3,3))
	C2a=widget_base(C2,/row)
	dum=widget_label(C2a,value=' l ',font=font)
	dum=widget_label(C2a,value=' : ')
	lampt=widget_text(C2a,xsize=6,ysize=1)
	dum=widget_label(C2a, value=' [Ang]')
	C2b=widget_base(C2,/row)
	dum=widget_label(C2b,value='FOCUS')
	dum=widget_label(C2b,value=' x ',font=font)
	xfoc=widget_text(C2b,xsize=6,ysize=1)
	dum=widget_label(C2b,value=' z ',font=font)
	zfoc=widget_text(C2b,xsize=6,ysize=1)
	C2c=widget_base(C2,/row)
	dum=widget_label(C2c,value='SAVE TO TREE')
	write= widget_button(C2c,value='WRITE')
	C2d=widget_base(C2,/row)
	dum=widget_label(C2d,value='use to constrain MOD1,MOD3 to MOD2: ') 
	C2dx=widget_base(C2d,/row,/nonexcl)
	cadbut=widget_button(C2dx,value='CAD')

	dum=widget_label(C,value='PLOTTING OPTIONS')
	C3=widget_base(C,frame=2,xsize=size,ysize=45,/column)
	C3a=widget_base(C3,/row)
	dum=widget_label(C3a,value='FULL dx')
	dxfull=widget_text(C3a,xsize=5,ysize=1,/edit)
	dum=widget_label(C3a,value='ZOOM dx')
	dxzoom=widget_text(C3a,xsize=5,ysize=1,/edit)
	C3ax=widget_base(C3a,/row,/nonexcl)
	zoombut=widget_button(C3ax,value='ZOOM')



	zoom=1.1
	xsize=135
	dum = widget_label(D,value='ELLIPSE/RESIDUAL')
	D1=widget_base(D,/row)
	DR4=widget_base(D1,frame=5)
	draw4=widget_draw(DR4,xsize=195*zoom,ysize=487*zoom)
	DR5=widget_base(D1,frame=5)
	draw5=widget_draw(DR5,xsize=xsize,ysize=487*zoom)
	DR6=widget_base(D1,frame=5)
	draw6=widget_draw(DR6,xsize=xsize,ysize=487*zoom)
	DR7=widget_base(D1,frame=5)
	draw7=widget_draw(DR7,xsize=xsize,ysize=487*zoom)
	DR8=widget_base(D1,frame=5)
	draw8=widget_draw(DR8,xsize=xsize,ysize=487*zoom)
	D2=widget_base(D,/row)
	;build u structure
	id={base:base,draw1:draw1,draw2:draw2,draw3:draw3,draw4:draw4,draw5:draw5,draw6:draw6,draw7:draw7,draw8:draw8,$
		shotid:shotid,module:module,load:load,quit:quit,print:print,stop:stop,message:message,$
		l_slider:l_slider,l_text:l_text,h_slider:h_slider,h_text:h_text,q_slider:q_slider,q_text:q_text,x_slider:x_slider,x_text:x_text,$
		z_slider:z_slider,z_text:z_text,a_slider:a_slider,a_text:a_text,b_slider:b_slider,b_text:b_text,c_slider:c_slider,c_text:c_text,$
		tree_table:tree_table,tx:tx,tz:tz,curr_table:curr_table,lampt:lampt,xfoc:xfoc,zfoc:zfoc,write:write,cadbut:cadbut,$
		dxfull:dxfull,dxzoom:dxzoom,zoombut:zoombut}

	;set defaults
	IF NOT keyword_set(shot) THEN shot=1070830020
	det=2
	dxz=0.1
	dxf=2.0

	stat={plot:0,zoom:0,ps:0,color:[0,0,0,0],xwcolor:[255,50,150,150],pscolor:[0,30,200,200],sym:[8,5,2],usym:[9,0,0],dxzoom:dxz,dxfull:dxf,cad:0,dat:0,finite:0,$
		xsize:[510,195*1.1,135],ysize:[510,487*1.1,487*1.1]}
	u={id:id,shot:shot,module:det,nlines:0,font:font,stat:stat}
	widget_control,base,set_uvalue=u
	widget_control,id.shotid,set_value=num2str(u.shot,1)
	widget_control,id.module,set_value=num2str(u.module,1)
	widget_control,id.dxzoom,set_value=num2str(u.stat.dxzoom,dp=2)
	widget_control,id.dxfull,set_value=num2str(u.stat.dxfull,dp=2)
	widget_control,id.tx,set_value='0.0'
	widget_control,id.tz,set_value='0.0'
	
	




	!except=0
	widget_control,base,/realize
	xmanager,'w_hirexsr_det_align',base

END
