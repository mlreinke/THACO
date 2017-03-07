PRO wbin_plot_chmap,u
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw3,get_value=draw_win
		window,0,xsize=u.plot.xsize[2],ysize=u.plot.ysize[2],/pixmap
		col=u.plot.col
        ENDIF
	CASE u.stat.branch OF
		0 : BEGIN
			ymax=max(u.adat.chmap)+1
			xoff=487
			xmax=487*3.0
			chmap=u.adat.chmap
			tit='BRANCH A'
		END
		1 : BEGIN
			ymax=max(u.bdat.chmap)+1
			xoff=0
			xmax=487
			chmap=u.bdat.chmap
			tit='BRANCH B'
		END
        ENDCASE
	plot,[0],[0],xr=[0,xmax],yr=[-2,ymax],/xsty,/ysty,xtit='ROW #',ytit='CHMAP',chars=1.2,tit=tit
	tmp=where(chmap[*,0] EQ u.ch-1)
	oplot,tmp[0]*[1,1],[-2,ymax],color=col[1],thick=2.0,linestyle=2.0
	oplot,last(tmp)*[1,1],[-2,ymax],color=col[1],thick=2.0,linestyle=2.0
	IF u.stat.branch EQ 0 THEN oplot,chmap[*,0],color=col[0]
	IF u.stat.branch EQ 1 THEN oplot,indgen(487)+xoff,chmap[*,0],color=col[0]
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize[2],u.plot.ysize[2],0,0,0]
	ENDIF
END

PRO wbin_plot_tmap,u
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw4,get_value=draw_win
		window,0,xsize=u.plot.xsize[2],ysize=u.plot.ysize[2],/pixmap
		col=u.plot.col
        ENDIF
	widget_control,u.id.nframes,get_value=nfr
	nfr=int(nfr[0])
	xmax=nfr
	CASE u.stat.branch OF
		0 : BEGIN
			tmap=u.adat.tmap
			tit='BRANCH A'
		END
		1 : BEGIN
			tmap=u.bdat.tmap
			tit='BRANCH B'
		END
        ENDCASE
	ymax=max(tmap)+1
	plot,[0],[0],xr=[0,xmax],yr=[-2,ymax],/xsty,/ysty,xtit='ROW #',ytit='TMAP',chars=1.2,tit=tit
	oplot,u.index*[1,1],[-2,ymax],color=col[1],thick=2.0,linestyle=2.0
	oplot,tmap,color=col[0]
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize[2],u.plot.ysize[2],0,0,0]
	ENDIF
END

PRO wbin_plot_vessel,u
	IF NOT u.stat.ps THEN BEGIN
		window,0,xsize=u.plot.xsize[0],ysize=u.plot.ysize[0],/pixmap
	ENDIF
	yrange=[-0.54,0.54]
	xrange=[0.40,0.98]
	plot, [0],[0],title=title,chars=1.0,xrange=xrange,yrange=yrange,xtit='R (m)',ytit='Z (m)',/xsty, /ysty
	oplot,u.dat.fs.rw,u.dat.fs.zw	
END

PRO wbin_plot_fs,u
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

PRO wbin_plot_pos,u
	ns=100
;	widget_control,u.id.draw1,get_value=draw_win
;	wset,draw_win
	CASE u.stat.branch OF
		0 : pos=*u.aspec.pos
		1 : pos=*u.bspec.pos
	ENDCASE
	
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

PRO wbin_plot_cx,u
	widget_control,u.id.draw1,get_value=draw_win
	wbin_plot_vessel,u
	IF u.stat.fs THEN wbin_plot_fs,u
	IF u.stat.pos THEN wbin_plot_pos,u
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize[0],u.plot.ysize[0],0,0,0]
	ENDIF
END

PRO wbin_plot_image,u
	IF NOT u.stat.image THEN ff=0.0 ELSE ff=1.0
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw2,get_value=draw_win
		window,0,xsize=u.plot.xsize[1],ysize=u.plot.ysize[1],/pixmap
	ENDIF
	CASE u.stat.branch OF
		0 : BEGIN
			image=rotate(*u.adat.raw[u.index],4)
			chmap=u.adat.chmap
			bnds=*u.adat.bnds
			offset=0
			chmax=u.adat.chmax
                END
		1 : BEGIN
			image=*u.bdat.raw[u.index]
			image=[image*0.0,image,image*0.0]
			image=rotate(image,4)
			chmap=u.bdat.chmap
			chmap=[chmap*0.0-1,chmap,chmap*0.0-1]
			bnds=*u.bdat.bnds
			offset=487
			chmax=u.bdat.chmax
                END
	ENDCASE
	x=size(image)
	widget_control,u.id.gain_slider,get_value=gain
	gain*=ff
	nx=x[1]
	ny=x[2]
	zoom=0.7

	;mask killed pixels
	killval=0.0
	imap=rotate(chmap,4)
	tmp=where(imap EQ -2)
	IF tmp[0] NE -1 THEN image[tmp]=killval
	
	i_pt=indgen(nx*zoom)/zoom
	j_pt=indgen(ny*zoom)/zoom
	new_pic=interpolate(image,i_pt,j_pt,/grid)
	loadct,39,/silent
	norm=max(new_pic)
	tv,new_pic/norm*256.0*gain
	loadct,12,/silent
	IF NOT u.stat.image THEN RETURN

	bnds+=offset
	r0=reform(bnds[0,*])
	r1=reform(bnds[1,*])
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

	;if zoom selected plot selected spectra and image
	IF u.stat.zoom AND u.itau NE -1 AND u.ch LT chmax THEN BEGIN
		
		;plot SPEC or AVESPEC in draw3
		CASE u.stat.branch OF
			0 : BEGIN
				iave=*u.aspec.ave[u.ch-1,u.itau]
				ispec=*u.aspec.spec[u.ch-1,u.itau]
			END
			1 : BEGIN
				iave=*u.bspec.ave[u.ch-1,u.itau]
				ispec=*u.bspec.spec[u.ch-1,u.itau]
                        END
		ENDCASE		
		IF NOT u.stat.ps THEN BEGIN
			widget_control,u.id.draw3,get_value=draw_win
			window,0,xsize=u.plot.xsize[2],ysize=u.plot.ysize[2],/pixmap
			color=u.plot.col
                ENDIF ELSE color=u.plot.pscol
		xr=u.plot.wl
		yr=[0,max(ispec[where(ispec[*,1] GE xr[0] AND ispec[*,1] LE xr[1]),0])]
		plot,[0],[0],xr=xr,yr=yr,xtit='Wavelength [Ang]',ytit='Spec. Br [AU]',chars=1.2,/xsty,/ysty
		oplot,ispec[*,1],ispec[*,0],psym=4
		oplot,iave[*,1],iave[*,0],color=color[2],thick=3.0

		IF NOT u.stat.ps THEN BEGIN
			wset,draw_win
			device,copy=[0,0,u.plot.xsize[2],u.plot.ysize[2],0,0,0]
                ENDIF ;ELSE stop

		;plot zoomed image in draw4
		zoom=3.076
		height=int(350/zoom)
		ex=int((height-(r1[ch]-r0[ch]))/2.0)
		cent=int((r0[ch]+r1[ch])/2.0)
                IF r0[ch]-ex LT 0 THEN BEGIN
                	lowbnd=0 
	                zoomx = zoom
        	        zoomy = zoom*(2.*ex+r1[ch]-r0[ch])/((2.*ex+r1[ch]-r0[ch])+2.*(r0[ch]-ex)) ;expand the y channel until it fills the screen
                    	height = int(350/zoomy)	
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
                        	height = int(350/zoomy)
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
			widget_control,u.id.draw4,get_value=draw_win
			window,0,xsize=u.plot.xsize[2],ysize=u.plot.ysize[2],/pixmap
		ENDIF
		loadct,39,/silent
		tv,new_pic/norm*256.0*gain
		loadct,12,/silent
		color=200
		frac=(r1[ch]-r0[ch])/float(height)/2.0
		plots, [0,1,1,0,0],0.5+frac*[1.0,1.0,-1.0,-1.0,1.0],color=color,/norm
		IF NOT u.stat.ps THEN BEGIN
			wset,draw_win
			device,copy=[0,0,u.plot.xsize[2],u.plot.ysize[2],0,0,0]
		ENDIF	
	ENDIF
END



PRO wbin_plotall,u
	wbin_spatial_text,u
	wbin_temporal_text,u
	wbin_plot_cx,u	
	wbin_plot_image,u
	IF NOT u.stat.zoom THEN BEGIN
		wbin_plot_chmap,u
		wbin_plot_tmap,u
	ENDIF
     END

;MODIFICATION HISTORY:
;	03/17/14:	C. Gao - added /NaN keyword in max/min. For
;                                some shots (1140227015 for example),
;                                there are NaNs in lambda array, which
;                                might (not always, why?) cause problem when using
;                                max/min

PRO wbin_load_data,u
	IF u.stat.dat THEN wbin_cleanup,u
	shot=u.shot
	tht=u.tht

	;load raw data
	hirexsr_ptr_image,shot,raw,lambda,t,const=const,morder=morder,pos=pos
	hirexsr_load_binning,shot,chmap,tch,tmap,good,chmax,tht=tht
	bnds=hirexsr_bin_bounds(chmap,chmax)
	adat={raw:raw,lambda:lambda,t:t,const:const,morder:morder,chmap:chmap,chk:ptr_new(chmap,/allocate),tch:tch,tmap:tmap,good:good,chmax:chmax,$
		pos:pos,bnds:ptr_new(bnds,/allocate),wl:[min(lambda,/NAN),max(lambda,/NAN)]}
	widget_control,u.id.message,set_value='BRANCH A - raw data and binning loaded '+num2str(shot,1),/app
	morder=0
	const=0
	pos=0
	
	hirexsr_ptr_image,shot,raw,lambda,t,const=const,morder=morder,pos=pos,/h
	hirexsr_load_binning,shot,chmap,tch,tmap,good,chmax,tht=tht,/h
	bnds=hirexsr_bin_bounds(chmap,chmax)
	bdat={raw:raw,lambda:lambda,t:t,const:const,morder:morder,chmap:chmap,chk:ptr_new(chmap,/allocate),tch:tch,tmap:tmap,good:good,chmax:chmax,$
		pos:pos,bnds:ptr_new(bnds,/allocate),wl:[min(lambda,/NAN),max(lambda,/NAN)]}
	widget_control,u.id.message,set_value='BRANCH B - raw data and binning loaded '+num2str(shot,1),/app
	widget_control,u.id.nframes,set_value=num2str(n(t)+1,1)
	widget_control,u.id.chmaxA,set_value=num2str(adat.chmax,1)
	widget_control,u.id.chmaxB,set_value=num2str(bdat.chmax,1)
	

; 	IF u.stat.kill THEN hirexsr_load_avespec,shot,specbr,lam,sig,resid,avetau,naveA,tht=tht,status=statusA ELSE statusA=0
; 	IF statusA THEN BEGIN
; 		avespecA=hirexsr_avespec_arr2ptr(specbr,lam,sig,resid,maxave=maxave)
; 		widget_control,u.id.message,set_value='BRANCH A - AVESPEC loaded '+num2str(shot,1),/app
;         ENDIF ELSE BEGIN
; 		avespecA=-1
; 	ENDELSE

; 	IF u.stat.kill THEN hirexsr_load_avespec,shot,specbr,lam,sig,resid,avetau,naveB,tht=tht,status=statusB,/h ELSE statusB=0
; 	IF statusB THEN BEGIN
; 		avespecB=hirexsr_avespec_arr2ptr(specbr,lam,sig,resid,maxave=maxave)
; 		widget_control,u.id.message,set_value='BRANCH B - AVESPEC loaded '+num2str(shot,1),/app
;         ENDIF ELSE BEGIN
; 		avespecB=-1
;         ENDELSE

	tauA=hirexsr_bin_time(adat.t,adat.tmap)
	naveA=0
	statusA=1
	FOR i=0,adat.chmax-1 DO BEGIN
		tmp=where(adat.chmap[*,0] EQ i)
		IF n(tmp)+1 GT naveA THEN naveA=n(tmp)+1
	ENDFOR
	tauB=hirexsr_bin_time(bdat.t,bdat.tmap)
	naveB=0
	statusB=1
	FOR i=0,bdat.chmax-1 DO BEGIN
		tmp=where(bdat.chmap[*,0] EQ i)
		IF n(tmp)+1 GT naveB THEN naveB=n(tmp)+1
	ENDFOR
	IF u.stat.kill THEN BEGIN
		IF statusA THEN BEGIN
			cnts=hirexsr_bin_image(adat.raw,adat.tmap)
			specA=hirexsr_bin_spec(cnts,tauA,adat.chmap,adat.tch,adat.lambda,adat.tmap,nchbins=adat.chmax,const=adat.const)
			avespecA=hirexsr_ave_spec(specA,naveA,chmap=adat.chmap)
			heap_free,cnts
			widget_control,u.id.message,set_value='BRANCH A - SPEC/AVESPEC computed '+num2str(shot,1),/app

                ENDIF ELSE specA=-1
		IF statusB THEN BEGIN
			cnts=hirexsr_bin_image(bdat.raw,bdat.tmap)
			specB=hirexsr_bin_spec(cnts,tauB,bdat.chmap,bdat.tch,bdat.lambda,bdat.tmap,nchbins=bdat.chmax,const=bdat.const)
			avespecB=hirexsr_ave_spec(specB,naveB,chmap=bdat.chmap)
			heap_free,cnts
			widget_control,u.id.message,set_value='BRANCH B - SPEC/AVE computed '+num2str(shot,1),/app
                ENDIF ELSE specA=-1	
	ENDIF ELSE BEGIN
		widget_control,u.id.message,set_value='SPEC/AVESPEC Loading Skipped '+num2str(shot,1),/app
		specA=-1
		specB=-1
		avespecA=-1
		avespecB=-1
        ENDELSE
	
	aspec={ave:avespecA,spec:specA,tau:tauA,nave:naveA,pos:ptr_new([0],/allocate_heap)}
	*aspec.pos=hirexsr_bin_pos(adat.pos,[min(adat.lambda),max(adat.lambda)],adat.chmap,adat.lambda,adat.chmax)
	bspec={ave:avespecB,spec:specB,tau:tauB,nave:naveB,pos:ptr_new([0],/allocate_heap)}
	*bspec.pos=hirexsr_bin_pos(bdat.pos,[min(bdat.lambda),max(bdat.lambda)],bdat.chmap,bdat.lambda,bdat.chmax)

	;load fluxsurface data
	fs=make_fs_struc(shot)
	widget_control,u.id.message,set_value='Loaded FS Data: '+num2str(shot,1),/app

	dat={fs:fs,t:t}
	u={id:u.id,shot:u.shot,tht:u.tht,time:u.time,index:u.index,itau:u.itau,ch:u.ch,stat:u.stat,plot:u.plot,aspec:aspec,adat:adat,bspec:bspec,bdat:bdat,dat:dat}
	u.stat.dat=1
	u.index=ipt(u.dat.t,u.time)
	CASE u.stat.branch OF
		0 : BEGIN
			u.itau=ipt(u.time,u.aspec.tau)
			u.plot.wl=u.adat.wl
		END
		1 : BEGIN
			u.itau=ipt(u.time,u.bspec.tau)
			u.plot.wl=u.bdat.wl
		END
        ENDCASE
	wbin_wavelength_text,u
	wbin_temporal_text,u
	wbin_spatial_text,u
	widget_control,u.id.base,set_uvalue=u
END

PRO wbin_cleanup,u
	IF u.stat.dat THEN BEGIN
		heap_free,u.aspec.ave
		heap_free,u.aspec.spec
		heap_free,u.aspec.pos
		heap_free,u.adat.raw
		heap_free,u.adat.bnds
		heap_free,u.adat.chk
		heap_free,u.bspec.ave
		heap_free,u.bspec.spec
		heap_free,u.bspec.pos
		heap_free,u.bdat.raw
		heap_free,u.bdat.bnds
		heap_free,u.bdat.chk
	ENDIF
END

PRO wbin_reset_branches,u
	widget_control,u.id.setA,set_button=0
	widget_control,u.id.setB,set_button=0
END

PRO wbin_wavelength_text,u
	widget_control,u.id.wlow,set_value=num2str(u.plot.wl[0],dp=4)
	widget_control,u.id.whigh,set_value=num2str(u.plot.wl[1],dp=4)
END

PRO wbin_temporal_text,u
	CASE u.stat.branch OF 
		0 : tmap=u.adat.tmap
		1 : tmap=u.bdat.tmap
        ENDCASE
	widget_control,u.id.findex,set_value=num2str(u.index,1)
	widget_control,u.id.itpt,set_value=num2str(tmap[u.index],1)
	widget_control,u.id.ntpt,set_value=num2str(max(tmap),1)
	tmp=where(tmap EQ tmap[u.index])
	u.itau=tmap[u.index]
	widget_control,u.id.flow,set_value=num2str(tmp[0],1)
	widget_control,u.id.fhigh,set_value=num2str(last(tmp),1)
	widget_control,u.id.tpt,set_value=num2str(u.time,dp=2)
	widget_control,u.id.t_slider,set_value=u.time*1.0e3
END



;MODIFICATION HISTORY:
;	03/17/14:	C. Gao - implemented the button of chdel
PRO wbin_spatial_text,u
	CASE u.stat.branch OF 
		0 : BEGIN
			chmap=u.adat.chmap
			bnds=*u.adat.bnds
                END
		1 : BEGIN
			chmap=u.bdat.chmap
			bnds=*u.bdat.bnds
                END
	ENDCASE
	widget_control,u.id.ch,set_value=num2str(u.ch,1)
	widget_control,u.id.chtot,set_value=num2str(max(chmap)+1,1)
	widget_control,u.id.rlow,set_value=num2str(bnds[0,u.ch-1],1)
	widget_control,u.id.rhigh,set_value=num2str(bnds[1,u.ch-1],1)
	widget_control,u.id.ch_slider,set_value=u.ch
END

PRO w_hirexsr_binning_event,event

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
					wbin_cleanup,u
					heap_gc
					!except=1
				
				END
				"SAVE" :BEGIN
					IF u.stat.dat THEN BEGIN
						WIDGET_CONTROL,/hourglass
						ask = dialog_message("Save CHAMP and TMAP data to tree?", /default_no, /question)
						IF ask EQ 'Yes' AND u.stat.dat THEN BEGIN
							;save branch A
							tmap=u.adat.tmap
							chmap=u.adat.chmap
							widget_control,u.id.nframes,get_value=nf
							nf=int(nf[0])
							widget_control,u.id.chmaxA,get_value=chmax
							chmax=int(chmax[0])
							good=intarr(chmax,nf)
							good[0:max(chmap[*,0]),0:max(tmap)]=1
							hirexsr_write_binning,u.shot,chmap=chmap,tch=-1,tmap=tmap,good=good,tht=u.tht,chmax=96
							widget_control,u.id.message,set_value='Branch A TMAP/CHMAP Saved',/app

							;save branch B
							tmap=u.bdat.tmap
							chmap=u.bdat.chmap
							widget_control,u.id.nframes,get_value=nf
							nf=int(nf[0])
							widget_control,u.id.chmaxB,get_value=chmax
							chmax=int(chmax[0])
							good=intarr(chmax,nf)
							good[0:max(chmap[*,0]),0:max(tmap)]=1
							hirexsr_write_binning,u.shot,chmap=chmap,tch=-1,tmap=tmap,good=good,tht=u.tht,chmax=32,/h
 							widget_control,u.id.message,set_value='Branch B TMAP/CHMAP Saved',/app
                                              ENDIF

                                        ENDIF	
				END
				"LOAD": BEGIN
                   			WIDGET_CONTROL, /hourglass
					wbin_load_data,u
					wbin_plotall,u	
				END
				"SETA" :BEGIN
					wbin_reset_branches,u
					widget_control,u.id.setA,set_button=1
					u.stat.branch=0
					IF u.stat.dat THEN BEGIN
						u.plot.wl=u.adat.wl
						wbin_wavelength_text,u
						wbin_plotall,u
					ENDIF
				END
				"SETB" :BEGIN
					wbin_reset_branches,u
					widget_control,u.id.setB,set_button=1
					u.stat.branch=1
					u.itau=ipt(u.time,u.aspec.tau)
					IF u.ch GT 16 THEN u.ch=16
					IF u.stat.dat THEN BEGIN
						u.plot.wl=u.bdat.wl
						wbin_wavelength_text,u
						wbin_plotall,u
					END
				END
      				"PBIN": BEGIN
					IF event.select EQ 1 THEN u.stat.bins=1 ELSE u.stat.bins=0
					plot_hemom_image,u
					IF u.stat.pos THEN plot_cx,u
				END
				"PZOOM" : BEGIN
					IF u.stat.kill THEN BEGIN
						IF event.select EQ 1 THEN u.stat.zoom=1 ELSE u.stat.zoom=0
						wbin_plot_image,u
						IF NOT u.stat.zoom THEN BEGIN
							wbin_plot_chmap,u
							wbin_plot_tmap,u
						ENDIF
					ENDIF ELSE BEGIN
						widget_control,u.id.message,set_value='Select KILL mode and reload to use ZOOM',/app
						widget_control,u.id.pzoom,set_button=0
					ENDELSE
				END
				"PFS" : BEGIN
					IF event.select EQ 1 THEN u.stat.fs=1 ELSE u.stat.fs=0
					wbin_plot_cx,u
				END
				"PPOS" : BEGIN
					IF event.select EQ 1 THEN u.stat.pos=1 ELSE u.stat.pos=0
					wbin_plot_cx,u
				END
				"PAVE" : BEGIN
					IF event.select EQ 1 THEN u.stat.ave=1 ELSE u.stat.ave=0
					IF u.stat.kill THEN wbin_plot_image,u
				END
				"APPLYAUTOCH" : BEGIN
					widget_control,u.id.clow,get_value=clow
					clow=int(clow[0])
					widget_control,u.id.chigh,get_value=chigh
					chigh=int(chigh[0])
					widget_control,u.id.noff,get_value=noff
					noff=int(noff[0])
					widget_control,u.id.nsub,get_value=nsub
					nsub=int(nsub[0])
			
					IF nsub LE 4 AND u.stat.dat THEN BEGIN
                                                CASE u.stat.branch OF
							0 : BEGIN
								chmap=u.adat.chmap
								chmax=u.adat.chmax
							END
							1 : BEGIN
								chmap=u.bdat.chmap
								chmax=u.bdat.chmax
							END
                                                ENDCASE
				
						imap=hirexsr_autochmap(nsub,noff)
						a=where(imap[*,0] EQ (clow-1)*nsub)
						b=where(imap[*,0] EQ chigh*nsub-1)
						FOR i=a[0]-2,last(b)+2 DO chmap[i,*]=imap[i,*]	;include extra -1's to account for residuals
						IF clow NE 1 THEN BEGIN				;shift lower channels
							tmp=where(chmap[0:a[0]-1,0] NE -1)
							diff=max(chmap[tmp,0])-chmap[a[0],0]
							tmp=where(chmap[a[0]:*,0] NE -1)+a[0]
							FOR i=0,n(tmp) DO chmap[tmp[i],*]+=diff+1
                                                ENDIF
						IF chigh NE chmax/4 THEN BEGIN			;shift upper channels
							tmp=where(chmap[last(b)+1:*,0] NE -1)+last(b)+1
							diff=chmap[b[0],0]-min(chmap[tmp,0])
							FOR i=0,n(tmp) DO chmap[tmp[i],*]+=diff+1
						ENDIF
						CASE u.stat.branch OF
							0 : BEGIN
								u.adat.chmap=chmap
								*u.adat.chk=chmap
								*u.adat.bnds=hirexsr_bin_bounds(chmap,u.adat.chmax)
								*u.aspec.pos=hirexsr_bin_pos(u.adat.pos,[min(u.adat.lambda),max(u.adat.lambda)],u.adat.chmap,u.adat.lambda,u.adat.chmax)
							END
							1 : BEGIN
								u.bdat.chmap=chmap
								*u.bdat.chk=chmap
								*u.bdat.bnds=hirexsr_bin_bounds(chmap,u.bdat.chmax)
								*u.bspec.pos=hirexsr_bin_pos(u.bdat.pos,[min(u.bdat.lambda),max(u.bdat.lambda)],u.bdat.chmap,u.bdat.lambda,u.bdat.chmax)
							END
                                                ENDCASE
                                        ENDIF
					wbin_spatial_text,u
                                        wbin_plotall,u
                                     END
                                
                                "CHDEL" : BEGIN
                                       print,'chdel'
                                       widget_control,u.id.ch,get_value=ch
                                       ch=int(ch[0])-1
                                      
                                       IF u.stat.dat THEN BEGIN
                                          CASE u.stat.branch OF
                                             0 : BEGIN
                                                chmap=u.adat.chmap
                                             END
                                             1 : BEGIN
                                                chmap=u.bdat.chmap
                                             END
                                          ENDCASE
                                          tmp=where(chmap EQ ch,count) 
                                          IF(count NE 0) THEN chmap[tmp]=-1
                                          tmp=where(chmap GT ch,count)
                                          IF(count NE 0) THEN chmap[tmp]-=1
                                          CASE u.stat.branch OF
                                             0 : BEGIN
                                                u.adat.chmap=chmap
                                                *u.adat.chk=chmap
                                                *u.adat.bnds=hirexsr_bin_bounds(chmap,u.adat.chmax)
                                                *u.aspec.pos=hirexsr_bin_pos(u.adat.pos,[min(u.adat.lambda),max(u.adat.lambda)],u.adat.chmap,u.adat.lambda,u.adat.chmax)
                                             END
                                             1 : BEGIN
                                                u.bdat.chmap=chmap
                                                *u.bdat.chk=chmap
                                                *u.bdat.bnds=hirexsr_bin_bounds(chmap,u.bdat.chmax)
                                                *u.bspec.pos=hirexsr_bin_pos(u.bdat.pos,[min(u.bdat.lambda),max(u.bdat.lambda)],u.bdat.chmap,u.bdat.lambda,u.bdat.chmax)
                                            END
                                          ENDCASE
                                       ENDIF
                                       wbin_spatial_text,u
                                       wbin_plotall,u                                              
                                    END

                                "CHADD" : BEGIN
                                       print,'chadd currently not available'
                                    END                                

				"APPLYAUTOT" : BEGIN
					widget_control,u.id.ntbin,get_value=ntbin
					ntbin=int(ntbin[0])
					widget_control,u.id.tlow,get_value=tlow
					tlow=float(tlow[0])
					widget_control,u.id.thigh,get_value=thigh
					thigh=float(thigh[0])
					widget_control,u.id.nframes,get_value=nframes
					nframes=int(nframes[0])
					IF u.stat.dat THEN BEGIN
						imap=hirexsr_autotmap(nframes,ntbin,shot=u.shot,tr=[tlow,thigh])
						CASE u.stat.branch OF
							0 : tmap=u.adat.tmap
							1 : tmap=u.bdat.tmap
                                                ENDCASE
						tmp=where(imap NE -1)
						low=max(tmap[0:tmp[0]-1])
						IF low NE -1 THEN tmap[tmp]=imap[tmp]+low+1 ELSE tmap[tmp]=imap[tmp]
						chk=where(tmap EQ -1)
						diff=tmap[last(tmp)]-tmap[last(tmp)+1]
						tmap[last(tmp)+1:*]+=diff+1
						tmap[chk]=-1
						CASE u.stat.branch OF
							0 : u.adat.tmap=tmap
							1 : u.bdat.tmap=tmap
						ENDCASE							
						wbin_temporal_text,u
						wbin_plot_tmap,u
					ENDIF
				END
				"COPYTMAP" : BEGIN
					CASE u.stat.branch OF
						0 : question='Confirm copy of Branch A TMAP to Branch B'
						1 : question='Confirm copy of Branch B TMAP to Branch A'
                                        ENDCASE
					ask = dialog_message(question, /default_no, /question)
					IF ask EQ 'Yes' AND u.stat.dat THEN BEGIN
						CASE u.stat.branch OF
							0 : u.bdat.tmap=u.adat.tmap
							1 : u.adat.tmap=u.bdat.tmap
                                                ENDCASE
						wbin_temporal_text,u
						wbin_plot_tmap,u
					ENDIF
				END
				"REBIN" : BEGIN
					IF u.stat.dat AND u.stat.kill THEN BEGIN	
						CASE u.stat.branch OF
							0 : BEGIN
								cnts=hirexsr_bin_image(u.adat.raw,u.adat.tmap)
								tau=hirexsr_bin_time(u.adat.t,u.adat.tmap)
								chmap=u.adat.chmap
								tch=[0]
								tmap=u.adat.tmap
								lambda=u.adat.lambda
								spec=u.aspec.spec
								ave=u.aspec.ave
								newmap=*u.adat.chk
								const=u.adat.const
								nave=u.aspec.nave							
							END
							1 : BEGIN
								cnts=hirexsr_bin_image(u.bdat.raw,u.bdat.tmap)
								tau=hirexsr_bin_time(u.bdat.t,u.bdat.tmap)
								chmap=u.bdat.chmap
								tch=[0]
								tmap=u.bdat.tmap
								lambda=u.bdat.lambda
								spec=u.bspec.spec
								ave=u.bspec.ave
								newmap=*u.bdat.chk
								const=u.bdat.const
								nave=u.bspec.nave
							END
						ENDCASE

						;update bins
						WIDGET_CONTROL, /hourglass
						hirexsr_rebin_afterkill,cnts,tau,chmap,tch,lambda,tmap,spec,newmap,const=const,fix=fix
						widget_control,u.id.message,set_value='Killed Pix Removed In '+num2str(total(fix),1)+' Spectra',/app
						ave=hirexsr_ave_spec(spec,nave,avespec=ave,fix=fix,chmap=chmap)
						widget_control,u.id.message,set_value='Killed Spectra Averaged',/app
						
						CASE u.stat.branch OF
							0 : u.adat.chmap=*u.adat.chk
							1 : u.bdat.chmap=*u.bdat.chk
						ENDCASE
						heap_free,cnts
						wbin_plot_image,u
					ENDIF
				END
				"MT" : BEGIN
					IF u.stat.dat AND u.index NE 0 THEN BEGIN
						u.index-=1
						u.time=u.dat.t[u.index]	
						wbin_temporal_text,u
						IF u.stat.pos THEN wbin_plot_cx,u
						IF u.stat.image THEN wbin_plot_image,u
						IF NOT u.stat.zoom THEN wbin_plot_tmap,u
					ENDIF
				END
				"PT" : BEGIN
					IF u.stat.dat THEN BEGIN
						CASE u.stat.branch OF
							0 : tmap=u.adat.tmap
							1 : tmap=u.bdat.tmap
						ENDCASE
						IF u.index NE max(tmap) THEN BEGIN
							u.index+=1
							u.time=u.dat.t[u.index]					
							wbin_temporal_text,u
							IF u.stat.pos THEN wbin_plot_cx,u
							IF u.stat.image THEN wbin_plot_image,u
							IF NOT u.stat.zoom THEN wbin_plot_tmap,u
						ENDIF
					ENDIF
				END
				"FDEL" : BEGIN
					IF u.stat.dat THEN BEGIN
						CASE u.stat.branch OF
							0 : tmap=u.adat.tmap
							1 : tmap=u.bdat.tmap
						ENDCASE							
						tmap[u.index]=-1
						CASE u.stat.branch OF
							0 : u.adat.tmap=tmap
							1 : u.bdat.tmap=tmap
						ENDCASE	
						wbin_temporal_text,u
						wbin_plot_tmap,u						
					ENDIF
                                END
				"FADD" : BEGIN
					IF u.stat.dat THEN BEGIN
						CASE u.stat.branch OF
							0 : tmap=u.adat.tmap
							1 : tmap=u.bdat.tmap
						ENDCASE							
						tmp=where(tmap EQ -1)
						tmap[u.index+1:*]+=2
						tmap[u.index]+=1
						tmap[tmp]=-1
						CASE u.stat.branch OF
							0 : u.adat.tmap=tmap
							1 : u.bdat.tmap=tmap
						ENDCASE	
						wbin_temporal_text,u
						wbin_plot_tmap,u						
					ENDIF
				END
				"MTLOW" : BEGIN
					IF u.stat.dat THEN BEGIN
						CASE u.stat.branch OF
							0 : tmap=u.adat.tmap
							1 : tmap=u.bdat.tmap
						ENDCASE					
						tmp=where(tmap EQ u.itau)
						tmap[tmp[0]-1]=tmap[tmp[0]]
						IF u.itau NE 0 THEN chk=where(tmap EQ u.itau-1) ELSE chk=0
						IF chk[0] EQ -1 THEN BEGIN
							tmap[tmp[0]-1:*]-=1
							neg=where(tmap LT -1)
							tmap[neg]=-1
						ENDIF
						CASE u.stat.branch OF
							0 : u.adat.tmap=tmap
							1 : u.bdat.tmap=tmap
						ENDCASE	
						wbin_temporal_text,u
						wbin_plot_tmap,u
					ENDIF
				END
				"PTHIGH" : BEGIN
					IF u.stat.dat THEN BEGIN
						CASE u.stat.branch OF
							0 : tmap=u.adat.tmap
							1 : tmap=u.bdat.tmap
						ENDCASE					
						tmp=where(tmap EQ u.itau)
						tmap[last(tmp)+1]=tmap[tmp[0]]
						IF u.itau NE max(tmap) THEN chk=where(tmap EQ u.itau+1) ELSE chk=0
						IF chk[0] EQ -1 THEN BEGIN
							tmap[last(tmp)+2:*]-=1
							neg=where(tmap LT -1)
							tmap[neg]=-1
						ENDIF
						CASE u.stat.branch OF
							0 : u.adat.tmap=tmap
							1 : u.bdat.tmap=tmap
						ENDCASE	
						wbin_temporal_text,u
						wbin_plot_tmap,u
					ENDIF
				END
				"MCH" : BEGIN
					IF u.stat.dat AND u.ch NE 1 THEN BEGIN
						u.ch-=1
						wbin_spatial_text,u
						IF u.stat.pos THEN wbin_plot_cx,u
						IF u.stat.image THEN wbin_plot_image,u
						IF NOT u.stat.zoom THEN wbin_plot_chmap,u
					ENDIF
				END
				"PCH" : BEGIN
					IF u.stat.dat THEN BEGIN
						CASE u.stat.branch OF
							0 : chmap=u.adat.chmap
							1 : chmap=u.bdat.chmap
						ENDCASE
						IF u.ch NE max(chmap)+1 THEN BEGIN
							u.ch+=1					
							wbin_spatial_text,u
							IF u.stat.pos THEN wbin_plot_cx,u
							IF u.stat.image THEN wbin_plot_image,u
							IF NOT u.stat.zoom THEN wbin_plot_chmap,u
						ENDIF
					ENDIF
				END

				"STOP" : BEGIN
					stop
				END 
				"KILLMODE" : IF event.select EQ 1 THEN u.stat.kill=1 ELSE u.stat.kill=0
			
			ELSE:
		ENDCASE
		END
  		"WIDGET_SLIDER": BEGIN
			widget_control,event.id,get_value=slider
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF 
				'CH_SLIDER' : BEGIN
					u.ch=slider
                                        IF u.stat.dat THEN BEGIN
						CASE u.stat.branch OF
							0 : bnds = *u.adat.bnds
							1 : bnds = *u.bdat.bnds
						ENDCASE
                                      		IF bnds[0,u.ch-1] NE -1 THEN BEGIN
                                              		wbin_spatial_text,u
	                                                wbin_plot_image,u
        	                                    	IF u.stat.pos THEN wbin_plot_cx,u
							IF NOT u.stat.zoom THEN wbin_plot_chmap,u                                            
                                	        ENDIF ELSE BEGIN
                                                	max = where(bnds[0, *] GT 0)
	                                                max = FIX(max(n_elements(max)-1))
        	                                        WIDGET_CONTROL,u.id.ch_slider, set_value=max+1
                	                       	ENDELSE
                                        ENDIF
				END
				'GAIN_SLIDER' : IF u.stat.dat THEN wbin_plot_image,u
				'T_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						index=ipt(u.dat.t,slider/1.0e3)
						IF slider/1.0e3 GE max(u.dat.fs.t) THEN index=ipt(u.dat.t,max(u.dat.fs.t))
						IF slider/1.0e3 LE min(u.dat.fs.t) THEN index=ipt(u.dat.t,min(u.dat.fs.t))
						IF index NE u.index THEN BEGIN
							u.time=u.dat.t[index]
							u.index=index
							widget_control,id.tpt,set_value=num2str(u.time,dp=3)
							wbin_temporal_text,u
							IF u.stat.pos OR u.stat.fs THEN wbin_plot_cx,u
							IF u.stat.image THEN wbin_plot_image,u
							IF NOT u.stat.zoom THEN wbin_plot_tmap,u
						ENDIF
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
				u.id.tpt : BEGIN
					IF u.stat.dat THEN BEGIN
						time=float(text[0])
						index=ipt(u.dat.t,time)
						IF time GE max(u.dat.fs.t) THEN index=ipt(u.dat.t,max(u.dat.fs.t))
						IF time LE min(u.dat.fs.t) THEN index=ipt(u.dat.t,min(u.dat.fs.t))
						IF index NE u.index THEN BEGIN
							u.time=u.dat.t[index]
							u.index=index
							widget_control,id.tpt,set_value=num2str(u.time,dp=3)
							wbin_temporal_text,u
							IF u.stat.pos OR u.stat.fs THEN wbin_plot_cx,u
							IF u.stat.image THEN wbin_plot_image,u
							IF NOT u.stat.zoom THEN wbin_plot_tmap,u
						ENDIF
					ENDIF
				END
				u.id.itpt : BEGIN
					IF u.stat.dat THEN BEGIN
						CASE u.stat.branch OF
							0 : BEGIN
								tmap = u.adat.tmap
								tau=u.aspec.tau
							END
							1 : BEGIN
								tmap = u.bdat.tmap
								tau=u.bspec.tau
							END
						ENDCASE
						time=tau[int(text[0])]
						index=ipt(u.dat.t,time)
						IF time GE max(u.dat.fs.t) THEN index=ipt(u.dat.t,max(u.dat.fs.t))
						IF time LE min(u.dat.fs.t) THEN index=ipt(u.dat.t,min(u.dat.fs.t))
						IF index NE u.index THEN BEGIN
							u.time=u.dat.t[index]
							u.index=index
							widget_control,id.tpt,set_value=num2str(u.time,dp=3)
							wbin_temporal_text,u
							IF u.stat.pos OR u.stat.fs THEN wbin_plot_cx,u
							IF u.stat.image THEN wbin_plot_image,u
							IF NOT u.stat.zoom THEN wbin_plot_tmap,u
						ENDIF
					ENDIF
				END
				u.id.ch : BEGIN
					u.ch=int(text[0])
					IF u.stat.dat THEN BEGIN
						CASE u.stat.branch OF
							0 : bnds = *u.adat.bnds
							1 : bnds = *u.bdat.bnds
						ENDCASE
                                      		IF bnds[0,u.ch-1] NE -1 THEN BEGIN
                                              		wbin_spatial_text,u
	                                                wbin_plot_image,u
        	                                    	IF u.stat.pos THEN wbin_plot_cx,u
							IF NOT u.stat.zoom THEN wbin_plot_chmap,u                                            
                                	        ENDIF ELSE BEGIN
                                                	max = where(bnds[0, *] GT 0)
	                                                max = FIX(max(n_elements(max)-1))
        	                                        WIDGET_CONTROL,u.id.ch_slider, set_value=max+1
                	                       	ENDELSE	
                                        ENDIF
				END
				u.id.wlow : BEGIN
					wl=float(text[0])
					u.plot.wl[0]=wl
					IF u.stat.dat THEN wbin_plot_image,u
                                END
				u.id.whigh : BEGIN
					wl=float(text[0])
					u.plot.wl[1]=wl
					IF u.stat.dat THEN wbin_plot_image,u
                                END
			ELSE: 
			ENDCASE
		END

               "WIDGET_DRAW": BEGIN
                    CASE event.id OF
                        u.id.draw4 : BEGIN
				IF u.stat.zoom AND u.stat.kill THEN BEGIN
                                	IF event.press EQ 1 THEN BEGIN
	                                 	ptru = PTR_NEW(u) ;			 ptr is killed in hirexsr_zoom_kill
                                    		wbin_zoom_kill, ptru 
	                                ENDIF 
				ENDIF		
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

PRO w_hirexsr_binning,shot=shot,tht=tht
	font='-adobe-symbol-medium-r-normal--12-120-75-75-p-74-adobe-fontspecific'
	loadct,12,/silent
	mdim=get_screen_size()
	IF mdim[0] NE 1600 AND mdim[1] NE 1200 THEN  base=widget_base(title='HIREXSR Binning Setup',/row,tlb_size_events=1,/scroll,$
		x_scroll_size=mdim[0]*0.95,y_scroll_size=mdim[1]*0.85) ELSE base=widget_base(title='HIREXSR Binning Setup',/row,tlb_size_events=1);, kill_notify= 'cleanheap_w_hirexsr_binning')
	A=widget_base(base,/column)
	B=widget_base(base,/column)
	space=widget_base(base,/row,xsize=5)	
	C=widget_base(base,/column)

	vzoom=0.6
	dum = widget_label(A,value='POS VECTORS')
	A1=widget_base(A,frame=5)
	draw1=widget_draw(A1,xsize=500*vzoom,ysize=818*vzoom)
	Aprime=widget_base(A,/column,/frame)
	Ax=widget_base(Aprime,/row)
	dum = widget_label(Ax,value='PLOT: ')
	Axp1=widget_base(Ax,/row,/nonexclusive)
	pfs=widget_button(Axp1,value='FS')
	ppos=widget_button(Axp1,value='POS')
	pbin=widget_button(Axp1,value='BINS')
	pzoom=widget_button(Axp1,value='ZOOM')
	pave=widget_button(Axp1,value='AVE')
	Ay=widget_base(Aprime,/row)
	dum = widget_label(Ay,value='WAVELENGTH RANGE:')
	wlow=widget_text(Ay,xsize=5,ysize=1,/edit)
	dum = widget_label(Ay,value=' < l < ',font=font)
	whigh=widget_text(Ay,xsize=5,ysize=1,/edit)
	dum=widget_label(Ay,value='[Ang]')


	A2=widget_base(A,/column,/frame)
	space=widget_base(A2,/row,ysize=2)	
	A2p1=widget_base(A2,/row)
	dum = widget_label(A2p1,value='SETUP BINNING FOR: ')
	A2p1x=widget_base(A2p1,/row,/nonexcl)
	setA=widget_button(A2p1x,value='BRANCH A')
	setB=widget_button(A2p1x,value='BRANCH B')
	space=widget_base(A2,/row,ysize=2)	

	A3=widget_base(A,/column,/frame)
	A3p1=widget_base(A3,/row)
	dum = widget_label(A3p1,value='SHOT:')
	shotid = widget_text(A3p1,xsize=10,ysize=1,/edit)
	dum = widget_label(A3p1,value='THT')
	thtid=widget_text(A3p1,xsize=2,ysize=1,/edit)
	A3p1x=widget_base(A3p1,/row,/nonexcl)
	killmode=widget_button(A3p1x, value=' KILL MODE ')
	A3p2=widget_base(A3,/row)
	dum = widget_label(A3p2,value=' ')
	save= widget_button(A3p2,value='SAVE')
	dum = widget_label(A3p2,value=' ')
	load= widget_button(A3p2,value='LOAD')
	dum = widget_label(A3p2,value=' ')
	quit= widget_button(A3p2,value='QUIT')
	dum = widget_label(A3p2,value=' ')
	stop= widget_button(A3p2,value='STOP')
	dum = widget_label(A3p2,value=' ')
	rebin= widget_button(A3p2,value='REBIN')


	A3p4=widget_base(A3,/row)
	message = widget_text(A3p4,xsize=40,ysize=6,/scroll,/wrap)

	A4=widget_base(A,/column,/frame)
	A4p1=widget_base(A4,/row)
	dum = widget_label(A4p1,value='# OF FRAMES COLLECTED: ')
	nframes=widget_text(A4p1,xsize=4,ysize=1)
	A4p2=widget_base(A4,/row)
	dum = widget_label(A4p2,value='MAX # CH:  ')
	dum = widget_label(A4p2,value='BRANCH A ')
	chmaxA=widget_text(A4p2,xsize=4,ysize=1)
	dum = widget_label(A4p2,value='BRANCH B ')
	chmaxB=widget_text(A4p2,xsize=4,ysize=1)

	A5=widget_base(A,/row)
	dum = widget_label(A5,value='TIME: ')
	tpt=widget_text(A5,xsize=5,ysize=1,/edit)
	t_slider=widget_slider(A5,xsize=500*vzoom-90,min=0,max=2000,value=1000,/drag,/suppress)


	imzoom=0.7
	dum = widget_label(B,value='     RAW IMAGE')
	B1=widget_base(B,/row)
	ch_slider=widget_slider(B1,ysize=3*487*imzoom,min=1,max=24*4,value=24.,/drag,/vert)
	B1p1=widget_base(B1,frame=5)
	draw2=widget_draw(B1p1,xsize=195*imzoom,ysize=3*487*imzoom)
	gain_slider=widget_slider(B,xsize=195*imzoom,min=1,max=40,value=1,/drag,/suppress)

	;DISPLAY
	dum = widget_label(C,value='TMAP/CHMAP Display')
	C1=widget_base(C,frame=5)
	draw3=widget_draw(C1,xsize=600,ysize=350)
	space=widget_base(C,/row,ysize=2)	
	C2=widget_base(C,frame=5)
	draw4=widget_draw(C2,xsize=600,ysize=350,/button_events)
	space=widget_base(C2,/row,ysize=2)

	;BINNING CONTROL TABS
	A2tb=widget_tab(C,location=3)

	;SPATIAL BINNING TAB
	A2t1=widget_base(A2tb,title=' SPATIAL ',/column,group_leader=base,/frame)
	A2t1_tit=widget_base(A2t1,/column,xsize=550,/align_center)
	dum = widget_label(A2t1_tit,value='SETUP OF SPATIAL BINNING - CHMAP')

	;channel selection
	A2p1=widget_base(A2t1,/row)
	dum = widget_label(A2p1,value='SELECT SPATIAL CH# ')
	ch=widget_text(A2p1,xsize=3,ysize=1,/edit)
	dum = widget_label(A2p1,value=' OF ')
	chtot=widget_text(A2p1,xsize=3,ysize=1)
	dum = widget_label(A2p1,value='  ')
	mch=widget_button(A2p1,value=' - ')
	pch=widget_button(A2p1,value=' + ')

	;non-uniform CHMAP setup
	A2p2=widget_base(A2t1,/row)
	dum = widget_label(A2p2,value='FOR SUB-MOD ')
	clow=widget_text(A2p2,xsize=3,ysize=1,/edit)
	dum = widget_label(A2p2,value=' TO ')
	chigh=widget_text(A2p2,xsize=3,ysize=1,/edit)
	dum = widget_label(A2p2,value=' OFFSET #: ')
	noff=widget_text(A2p2,xsize=3,ysize=1,/edit)
	dum=widget_label(A2p2,value=' CH/BIN: ')
	nsub=widget_text(A2p2,xsize=3,ysize=1,/edit)
	dum = widget_label(A2p2,value='  ')	
	applyautoch=widget_button(A2p2, value='  APPLY AUTOBIN  ')

	;insert or remove rows from selected CHMAP
	A2p3=widget_base(A2t1,/row)
	dum = widget_label(A2p3,value='ADD ROWS TO CURRENT CHANNEL: ')
	rlow=widget_text(A2p3,xsize=4,ysize=1,/edit)
	dum = widget_label(A2p3,value='   TO   ')
	rhigh=widget_text(A2p3,xsize=4,ysize=1,/edit)
	dum = widget_label(A2p3,value='   ')
	chdel=widget_button(A2p3,value=' REMOVE CH ')
	dum = widget_label(A2p3,value='   ')
	chadd=widget_button(A2p3,value=' ADD CH ')

	A2p4=widget_base(A2t1,/row)
	dum = widget_label(A2p4,value='                              ')
	mlow=widget_button(A2p4,value=' - ')
	dum = widget_label(A2p4,value='          ')
	phigh=widget_button(A2p4,value=' + ')

	;TEMPORAL BINNING TAB
	A2t2=widget_base(A2tb,title=' TEMPORAL ',/column,group_leader=base,/frame)
	A2t2_tit=widget_base(A2t2,/column,xsize=550,/align_center)
	dum = widget_label(A2t2_tit,value='SETUP OF TEMPORAL BINNING - TMAP')

	;time bin selection
	A2x1=widget_base(A2t2,/row)
	dum = widget_label(A2x1,value='TIME BIN ')
	itpt=widget_text(A2x1,xsize=3,ysize=1,/edit)
	dum = widget_label(A2x1,value=' OF ')
	ntpt=widget_text(A2x1,xsize=3,ysize=1)
	dum = widget_label(A2x1,value='      ADJUST ACTIVE FRAME')
	findex=widget_text(A2x1,xsize=4,ysize=1)
	dum = widget_label(A2x1,value=' ')
	mt=widget_button(A2x1,value=' - ')
	pt=widget_button(A2x1,value=' + ')
	space=widget_base(A2t2,/row,ysize=2)

	;non-uniform TMAP setup
	A2x2=widget_base(A2t2,/row)
	dum = widget_label(A2x2,value=' FROM ')
	tlow=widget_text(A2x2,xsize=3,ysize=1,/edit)
	dum = widget_label(A2x2,value=' < t < ')
	thigh=widget_text(A2x2,xsize=3,ysize=1,/edit)
	dum = widget_label(A2x2,value=' [sec]   ')
	dum = widget_label(A2x2,value=' BIN AT #FRAMES/TMAP: ')
	ntbin=widget_text(A2x2,xsize=3,ysize=1,/edit)
	dum = widget_label(A2x2,value='  ')	
	applyautot=widget_button(A2x2, value='  APPLY AUTOBIN  ')

	;insert remove frames from selected TMAP
	A2x3=widget_base(A2t2,/row)
	dum = widget_label(A2x3,value='ADD FRAMES TO CURRENT BIN: ')
	flow=widget_text(A2x3,xsize=3,ysize=1,/edit)
	dum = widget_label(A2x3,value='   TO   ')
	fhigh=widget_text(A2x3,xsize=3,ysize=1,/edit)
	dum = widget_label(A2x3,value='   ')
	fdel=widget_button(A2x3,value=' REMOVE FRAME ')
	dum = widget_label(A2x3,value='   ')
	fadd=widget_button(A2x3,value=' ADD BIN ')
	
	A2x4=widget_base(A2t2,/row)
	dum = widget_label(A2x4,value='                           ')
	mtlow=widget_button(A2x4,value=' - ')
	dum = widget_label(A2x4,value='         ')
	pthigh=widget_button(A2x4,value=' + ')

	A2x5=widget_base(A2t2,/row,/align_center)
	copytmap=widget_button(A2x5,value=' COPY TMAP BETWEEN BRANCHES')
	IF NOT keyword_set(shot) THEN shot=1120224012
	IF NOT keyword_set(tht) THEN tht=0
	IF NOT keyword_set(time) THEN time=1.0
	IF NOT keyword_set(branch) THEN branch=0
	
	id={base:base,shotid:shotid,thtid:thtid,save:save,load:load,quit:quit,stop:stop,rebin:rebin,message:message,$
		draw1:draw1,draw2:draw2,draw3:draw3,draw4:draw4,$
		pfs:pfs,ppos:ppos,pbin:pbin,pzoom:pzoom,pave:pave,setA:setA,setB:setB,killmode:killmode,wlow:wlow,whigh:whigh,$
	   	nframes:nframes,chmaxA:chmaxA,chmaxB:chmaxB,tpt:tpt,t_slider:t_slider,ch_slider:ch_slider,gain_slider:gain_slider,$
		ch:ch,chtot:chtot,mch:mch,pch:pch,clow:clow,chigh:chigh,noff:noff,nsub:nsub,applyautoch:applyautoch,chdel:chdel,chadd:chadd,$
		rlow:rlow,rhigh:rhigh,mlow:mlow,phigh:phigh,$
		itpt:itpt,ntpt:ntpt,findex:findex,mt:mt,pt:pt,$
		tlow:tlow,thigh:thigh,ntbin:ntbin,applyautot:applyautot,$
		flow:flow,fhigh:fhigh,mtlow:mtlow,pthigh:pthigh,fdel:fdel,fadd:fadd,copytmap:copytmap}

	stat={dat:0,kill:0,image:1,fs:1,pos:1,bins:0,zoom:0,ave:0,ps:0,branch:0}
	plot={xsize:[500*vzoom,195*imzoom,600],ysize:[818*vzoom,3*487*imzoom,350],col:[255,50,150],pscol:[0,30,200],wl:[0.0,0.0]}
	u={id:id,shot:shot,tht:tht,time:time,index:0,itau:0,ch:10,stat:stat,plot:plot}
	widget_control,base,set_uvalue=u
	widget_control,id.shotid,set_value=num2str(u.shot,1)
	widget_control,id.thtid,set_value=num2str(u.tht,1)
	widget_control,id.tpt,set_value=num2str(u.time,dp=2)
	widget_control,id.t_slider,set_value=u.time*1.0e3
	widget_control,u.id.pfs,set_button=u.stat.fs
	widget_control,u.id.ppos,set_button=u.stat.pos
	widget_control,u.id.pbin,set_button=u.stat.bins
	widget_control,u.id.pzoom,set_button=u.stat.zoom
	widget_control,u.id.pave,set_button=u.stat.ave
	CASE u.stat.branch OF
		0 : widget_control,u.id.setA,set_button=1
		1 : widget_control,u.id.setB,set_button=1
        ENDCASE
	widget_control,u.id.gain_slider,set_value=10

	!except=0			
	widget_control,base,/realize
	xmanager,'w_hirexsr_binning',base
END
