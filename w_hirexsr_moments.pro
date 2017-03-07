PRO wmoments2ps,u
	psplot
	u.stat.ps=1
	wmom_plot_spatial,u
	wmom_plot_spec,u
	wmom_plot_time,u
	wmom_plot_subspec,u
	psc
	u.stat.ps=0
	xwplot
	set_plot,'x'
END

PRO wmom_plot_mconv,u
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw1,get_value=draw_win
		window,0,xsize=u.plot.xsize[0],ysize=u.plot.ysize[0],/pixmap
		color=u.plot.col
		ls=1.0
        ENDIF ELSE BEGIN
		xsize=4.7
		ysize=xsize*u.plot.ysize[0]/u.plot.xsize[0]
		!p.thick=6
		!p.charthick=2.0
		!x.thick=5.0
		!y.thick=5.0
		device, xsize=xsize, ysize=ysize, /inches
		color=u.plot.pscol
		ls=0.9
        ENDELSE

	ispec=*u.dat.spec[u.ch-1,u.index]
	icoefs=*u.dat.coefs[u.ch-1,u.index]
	dlam=make(0.7,6.0,25)
	convmom=hirexsr_moment_conv(ispec,icoefs,u.dat.labels,u.dat.lam_o,u.line,dlam)
	!p.multi=[0,0,3]
	widget_control,u.id.dlampt,get_value=dlamstr
	dlampt=last(float(dlamstr))
	xr=[0,max(dlam)]

	;0th moment
	maxpt=max(convmom[*,0,0]+convmom[*,0,1])
	yr=[0,maxpt*1.05]
	plot,dlam,convmom[*,0,0],xtit='d'+n2g('lambda')+' [mAng]',ytit='0!uth!n Moment',chars=2.4*ls,xr=xr,/xsty,yr=yr,/ysty
	oploterror,dlam,convmom[*,0,0],convmom[*,0,1]
	loadct,39,/silent
	oplot,dlampt*[1,1],yr,color=color[5],linestyle=3.0
	m0pt=interpol(convmom[*,0,0],dlam,dlampt)
	oplot,xr,m0pt*[1,1],color=color[5],linestyle=3.0	
	loadct,12,/silent

	;1st moment
	maxpt=max(convmom[*,1,0]+convmom[*,1,1])
	minpt=min(convmom[*,1,0]-convmom[*,1,1])
	IF maxpt GT 0 THEN maxpt*=1.05 ELSE maxpt*=0.95
	IF minpt LT 0 THEN minpt*=1.05 ELSE minpt*=0.95
	yr=[minpt,maxpt]*1.0e3
	plot,dlam,convmom[*,1,0],xtit='d'+n2g('lambda')+' [mAng]',ytit='1!ust!n Moment',chars=2.4*ls,xr=xr,/xsty,yr=yr,/ysty
	oploterror,dlam,convmom[*,1,0]*1.0e3,convmom[*,1,1]*1.0e3
	loadct,39,/silent
	oplot,dlampt*[1.0,1.0],yr,color=color[5],linestyle=3.0
	m1pt=interpol(convmom[*,1,0],dlam,dlampt)
	oplot,xr,m1pt*[1,1]*1.0e3,color=color[5],linestyle=3.0	
	loadct,12,/silent
	oplot,xr,[0,0],linestyle=2.0

	;2nd moment
	maxpt=max(convmom[*,2,0]+convmom[*,2,1])
	yr=[0,maxpt*1.05]*1.0e6
	plot,dlam,convmom[*,2,0],xtit='d'+n2g('lambda')+' [mAng]',ytit='2!und!n Moment',chars=2.4*ls,xr=xr,/xsty,yr=yr,/ysty
	oploterror,dlam,convmom[*,2,0]*1.0e6,convmom[*,2,1]*1.0e6
	loadct,39,/silent
	oplot,dlampt*[1,1],yr,color=color[5],linestyle=3.0
	m2pt=interpol(convmom[*,2,0]*1.0e6,dlam,dlampt)
	oplot,xr,m2pt*[1,1],color=color[5],linestyle=3.0	
	loadct,12,/silent
	!p.multi=0
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize[0],u.plot.ysize[0],0,0,0]
	ENDIF

END

PRO wmom_plot_spatial,u
	IF u.stat.my EQ 2 THEN BEGIN
		wmom_plot_mconv,u
		RETURN
	ENDIF
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw1,get_value=draw_win
		window,0,xsize=u.plot.xsize[0],ysize=u.plot.ysize[0],/pixmap
		color=u.plot.col
		ls=1.0
        ENDIF ELSE BEGIN
		xsize=5.0
		ysize=xsize*u.plot.ysize[0]/u.plot.xsize[0]
		!p.thick=6
		!p.charthick=2.0
		!x.thick=5.0
		!y.thick=5.0
		device, xsize=xsize, ysize=ysize, /inches
		color=u.plot.pscol
  		ls=0.9
        ENDELSE
	imom=*u.dat.mom[u.index]
	ilint=*u.dat.lint[u.index]
	nch=n(imom[*,6])+1
	yplot=fltarr(nch,3)
	yerr=fltarr(nch,3)
	CASE u.stat.my OF 
		0 : BEGIN		;moments
			yplot=imom[*,0:2]
			yplot[*,1]*=1.0e3
			yplot[*,2]*=1.0e6
			yerr=imom[*,3:5]
			tit0='0!uth!n Moment'
			tit1='10!u-3!n 1!ust!n Moment'
			tit2='10!u-6!n 2!und!n Moment'
			lab0='M!i0!n'
			lab1='M!i1!n'
			lab2='M!i2!n'
			unit0=''
			unit1=''
			unit2=''									
                END
		1 : BEGIN		;line-integrated profiles
			yplot[*,0]=imom[*,0]
			yerr[*,0]=imom[*,3]
			yplot[*,1]=ilint[*,0]
			yplot[*,2]=ilint[*,2]
			yerr[*,1]=ilint[*,1]
			yerr[*,2]=ilint[*,3]		
			tit0='Brightness'
			tit1='Rotation [kHz]'
			tit2='Ion Temp. [keV]'
			lab0='Br'
			lab1=n2g('omega')
			lab2='T!ii!n'
 			unit0='[AU]'
			unit1='[kHz]'
			unit2='[keV]'
               END
		2 : BEGIN		;convergence tests

                END
	ENDCASE
		
	CASE u.stat.mx OF
		0 : BEGIN		;CH#
			xplot=indgen(nch)+1
			xtit='CH#'
			IF u.plot.mauto[0] THEN u.plot.xr=[0,nch+1]
                END
		1 : BEGIN		;r/a
			xplot=ilint[*,10]
			xtit='r/a'
			IF u.plot.mauto[0] THEN u.plot.xr=[0,max(xplot)]
                END
		2 : BEGIN		;PSIN
			xplot=ilint[*,4]
			xtit=n2g('psi')+'!iN!n'
			IF u.plot.mauto[0] THEN u.plot.xr=[0,max(xplot)]
		END
        ENDCASE
	IF u.plot.mauto[0] THEN BEGIN
		widget_control,u.id.mxlow,set_value=num2str(u.plot.xr[0],dp=2)
		widget_control,u.id.mxhigh,set_value=num2str(u.plot.xr[1],dp=2)
	ENDIF
	!p.multi=[0,0,3]
	tmp=where(xplot GE u.plot.xr[0] AND xplot LE u.plot.xr[1])
	;plot 0th
	IF u.plot.mauto[1] THEN BEGIN
		u.plot.yr[*,0]=[min(yplot[tmp,0]-yerr[tmp,0]) < 0,max(yplot[tmp,0]+yerr[tmp,0]) > 0]*1.1
		widget_control,u.id.my0low,set_value=num2str(u.plot.yr[0,0],dp=2)
		widget_control,u.id.my0high,set_value=num2str(u.plot.yr[1,0],dp=2)
	ENDIF
	plot,[0],[0],xr=u.plot.xr,yr=u.plot.yr[*,0],xtit=xtit,ytit=tit0,/xsty,/ysty,chars=2.4*ls
	IF u.stat.ps THEN xyouts,u.plot.xr[0]+0.8*(u.plot.xr[1]-u.plot.xr[0]),u.plot.yr[1,0]+0.02*(u.plot.yr[1,0]-u.plot.yr[0,0]),'t='+num2str(u.time,dp=2)+' [sec]',$
		chars=0.8*ls,color=color[0]
	fcase3=where(imom[*,24] EQ 3)
	fcase0=where(imom[*,24] EQ 0)
	IF fcase3[0] NE -1 THEN oploterror,xplot[fcase3],yplot[fcase3,0],yerr[fcase3,0],psym=7,color=color[0]
	IF fcase0[0] NE -1 THEN oploterror,xplot[fcase0],yplot[fcase0,0],yerr[fcase0,0],psym=6,color=color[0]
	loadct,2,/silent
	CASE imom[u.ch-1,24] OF
		0 : sym=6
		3 : sym=7
		ELSE : sym=3
	ENDCASE
	oploterror,[xplot[u.ch-1]],[yplot[u.ch-1,0]],[yerr[u.ch-1,0]],psym=sym,color=color[3],errcolor=color[3]
	oplot,u.plot.xr,yplot[u.ch-1,0]*[1,1],linestyle=1,color=color[3]
	oplot,xplot[u.ch-1]*[1,1],u.plot.yr[*,0],linestyle=1,color=color[3]
	xyouts,u.plot.xr[0]+0.1*(u.plot.xr[1]-u.plot.xr[0]),u.plot.yr[1,0]+0.03*(u.plot.yr[1,0]-u.plot.yr[0,0]),lab0+'='+num2str(yplot[u.ch-1,0],dp=2)+' '+unit0,$
		chars=1.25*ls,color=color[3]
	loadct,12,/silent
	IF u.plot.mauto[2] THEN BEGIN
		u.plot.yr[*,1]=[min(yplot[tmp,1]-yerr[tmp,1]) < 0,max(yplot[tmp,1]+yerr[tmp,1]) > 0]*1.1
		widget_control,u.id.my1low,set_value=num2str(u.plot.yr[0,1],dp=2)
		widget_control,u.id.my1high,set_value=num2str(u.plot.yr[1,1],dp=2)
	ENDIF
	plot,[0],[0],xr=u.plot.xr,yr=u.plot.yr[*,1],xtit=xtit,ytit=tit1,/xsty,/ysty,chars=2.4*ls
	IF fcase3[0] NE -1 THEN oploterror,xplot[fcase3],yplot[fcase3,1],yerr[fcase3,1],psym=7,color=color[0],errcolor=color[0]
	IF fcase0[0] NE -1 THEN oploterror,xplot[fcase0],yplot[fcase0,1],yerr[fcase0,1],psym=6,color=color[0],errcolor=color[0]
	loadct,2,/silent
	oploterror,[xplot[u.ch-1]],[yplot[u.ch-1,1]],[yerr[u.ch-1,1]],psym=sym,color=color[3],errcolor=color[3]
	oplot,u.plot.xr,yplot[u.ch-1,1]*[1,1],linestyle=1,color=color[3]
	oplot,xplot[u.ch-1]*[1,1],u.plot.yr[*,1],linestyle=1,color=color[3]
	xyouts,u.plot.xr[0]+0.1*(u.plot.xr[1]-u.plot.xr[0]),u.plot.yr[1,1]+0.03*(u.plot.yr[1,1]-u.plot.yr[0,1]),lab1+'='+num2str(yplot[u.ch-1,1],dp=2)+' '+unit1,$
		chars=1.25*ls,color=color[3]
	loadct,12,/silent
	oplot,u.plot.xr,[0,0],linestyle=2
	IF u.plot.mauto[3] THEN BEGIN
		u.plot.yr[*,2]=[min(yplot[tmp,2]-yerr[tmp,2]) < 0,max(yplot[tmp,2]+yerr[tmp,2]) > 0]*1.1
		widget_control,u.id.my2low,set_value=num2str(u.plot.yr[0,2],dp=2)
		widget_control,u.id.my2high,set_value=num2str(u.plot.yr[1,2],dp=2)
	ENDIF
	plot,[0],[0],xr=u.plot.xr,yr=u.plot.yr[*,2],xtit=xtit,ytit=tit2,/xsty,/ysty,chars=2.4*ls
	IF fcase3[0] NE -1 THEN oploterror,xplot[fcase3],yplot[fcase3,2],yerr[fcase3,2],psym=7,color=color[0]
	IF fcase0[0] NE -1 THEN oploterror,xplot[fcase0],yplot[fcase0,2],yerr[fcase0,2],psym=6,color=color[0]
	loadct,2,/silent
	oploterror,[xplot[u.ch-1]],[yplot[u.ch-1,2]],[yerr[u.ch-1,2]],psym=sym,color=color[3],errcolor=color[3]
	oplot,u.plot.xr,yplot[u.ch-1,2]*[1,1],linestyle=1,color=color[3]
	oplot,xplot[u.ch-1]*[1,1],u.plot.yr[*,2],linestyle=1,color=color[3]
	xyouts,u.plot.xr[0]+0.1*(u.plot.xr[1]-u.plot.xr[0]),u.plot.yr[1,2]+0.03*(u.plot.yr[1,2]-u.plot.yr[0,2]),lab2+'='+num2str(yplot[u.ch-1,2],dp=2)+' '+unit2,$
		chars=1.25*ls,color=color[3]

	loadct,12,/silent
	IF u.stat.ps THEN xyouts,u.plot.xr[1]+0.04*(u.plot.xr[1]-u.plot.xr[0]),0.05*(u.plot.yr[1,2]-u.plot.yr[0,2]),$
		num2str(u.shot,1)+'  LINE='+num2str(u.line,1)+'  THT='+num2str(u.tht,1),orient=90,chars=0.8*ls,color=color[0]
	!p.multi=0
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize[0],u.plot.ysize[0],0,0,0]
	ENDIF
END

PRO wmom_plot_time,u
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw3,get_value=draw_win
		window,0,xsize=u.plot.xsize[2],ysize=u.plot.ysize[2],/pixmap
		color=u.plot.col
		ls=1.0
        ENDIF ELSE BEGIN
		xsize=7.0
		ysize=xsize*u.plot.ysize[2]/u.plot.xsize[2]
		!p.thick=6
		!p.charthick=2.0
		!x.thick=5.0
		!y.thick=5.0
		device, xsize=xsize, ysize=ysize, /inches
		color=u.plot.pscol
		ls=0.85
	ENDELSE

	IF u.stat.tj EQ 0 THEN chlist=[u.ch]
	IF u.stat.tj EQ 1 THEN chlist=*u.plot.list	
	
	tmp=where(u.dat.tau NE -1)
	xplot=u.dat.tau[tmp]
	FOR j=0,n(chlist) DO BEGIN
		itlint=*u.dat.tlint[chlist[j]-1]
		x=size(itlint)
		CASE u.stat.ty OF 
			0 : BEGIN
				yplot=itlint[tmp,8]
				yerr=itlint[tmp,9]
				ytit='Brightness [AU]'		
			END
			1 : BEGIN
				yplot=itlint[tmp,0]
				yerr=itlint[tmp,1]
				ytit=n2g('omega')+' [kHz]'
                	END
			2 : BEGIN
				yplot=itlint[tmp,6]
				yerr=itlint[tmp,7]
				ytit='v!i'+n2g('theta')+'!n [km/s]'
	                END
			3 : BEGIN
				yplot=itlint[tmp,2]
				yerr=itlint[tmp,3]
				ytit='Ion Temperature [keV]'
			END
			ELSE :	
                ENDCASE
		tit=''
		IF j EQ 0 THEN yr=[min(yplot) < 0, max(yplot) > 0]*1.1
		IF j EQ 0 THEN plot,[0],[0],xr=u.plot.tr,xtit='Time [sec]',yr=yr,ytit=ytit,/xsty,/ysty,chars=1.2*ls
		oploterror,xplot,yplot,yerr,col=color[j],errcolor=color[j]
		fcase3=where(itlint[*,x[2]-1] EQ 3)
		fcase0=where(itlint[*,x[2]-1] EQ 0)
		IF fcase3[0] NE -1 THEN oploterror,xplot[fcase3],yplot[fcase3],yerr[fcase3],psym=7,col=color[j],errcolor=color[j]
		IF fcase0[0] NE -1 THEN oploterror,xplot[fcase0],yplot[fcase0],yerr[fcase0],psym=6,col=color[j],errcolor=color[j]
		IF u.plot.jrho THEN BEGIN
			IF j GE u.dat.nch/2 THEN add=' (Z > 0)' ELSE add=' (Z < 0)'
			str='r/a='+num2str(itlint[u.index,10],dp=2)+add
                ENDIF ELSE BEGIN
			str='CH#='+num2str(chlist[j],1)
		ENDELSE
		xyouts,u.plot.tr[0]+0.05*(u.plot.tr[1]-u.plot.tr[0]),yr[0]+(0.90-0.1*j)*yr[1]-yr[0],str,color=color[j],chars=1.0*ls
	ENDFOR
	oplot,u.time*[1,1],yr,linestyle=1
	IF u.stat.ty EQ 1 OR u.stat.ty EQ 2 THEN oplot,u.plot.tr,[0,0],linestyle=2.0
	IF u.stat.ps THEN xyouts,u.plot.tr[1]+0.03*(u.plot.tr[1]-u.plot.tr[0]),yr[0]+0.05*(yr[1]-yr[0]),$
		num2str(u.shot,1)+'  LINE='+num2str(u.line,1)+'  THT='+num2str(u.tht,1),orient=90,chars=0.75*ls,color=color[0]

	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize[2],u.plot.ysize[2],0,0,0]
	ENDIF
END

PRO wmom_plot_spec,u
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw2,get_value=draw_win
		window,0,xsize=u.plot.xsize[1],ysize=u.plot.ysize[1],/pixmap
		color=u.plot.col
		ls=1.0
        ENDIF ELSE BEGIN
		xsize=6.5
		ysize=xsize*u.plot.ysize[1]/u.plot.xsize[1]
		!p.thick=6
		!p.charthick=2.0
		!x.thick=5.0
		!y.thick=5.0
		device, xsize=xsize, ysize=ysize, /inches
		color=u.plot.pscol
		ls=1.0
	ENDELSE
	ispec=*u.dat.spec[u.ch-1,u.index] 
	icoefs=*u.dat.coefs[u.ch-1,u.index]
	tmp=where(ispec[*,1] GE u.plot.wl[0] AND ispec[*,1] LE u.plot.wl[1])
	lamfit=make(u.plot.wl[0],u.plot.wl[1],100)
	maxpt=max(ispec[tmp,0])*1.05
	pos=[0.09,0.4,0.95,0.98]
	plot,[0],[0],xr=u.plot.wl,yr=[0,maxpt],ytit='B!l'+n2g('lambda')+'!n [AU]',/xsty,/ysty,pos=pos,xtit='Wavelength [Ang]'
	oploterror,ispec[tmp,1],ispec[tmp,0],ispec[tmp,2],psym=8,symsize=0.25
	IF u.stat.ps THEN xyouts,u.plot.wl[1]+0.03*(u.plot.wl[1]-u.plot.wl[0]),0.05*maxpt,$
		num2str(u.shot,1)+'  THT='+num2str(u.tht,1)+'  t='+num2str(u.time,dp=2),orient=90,chars=0.75*ls,color=color[0]
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
		FOR i=0,nline-1 DO oplot,lamfit,gaussian_fits(lamfit,icoefs[i*3:(i+1)*3-1]),color=color[2]
		oplot,lamfit,gaussian_fits(lamfit,icoefs,base=base),color=color[4],thick=2
		IF n(base) NE 0 THEN oplot,lamfit,base,color=60 ELSE oplot,[min(lamfit),max(lamfit)],base*[1.,1],color=color[1]
		;CASE nb OF 
		;	0 : oplot,[lamfit[0],last(lamfit)],icoefs[ncoefs-1]*[1.0,1.0],color=60
		;	1 : oplot,lamfit,lamfit*icoefs[ncoefs-2]+icoefs[ncoefs-1],color=60
		;	2 : oplot,lamfit,lamfit^2*icoefs[ncoefs-3]+lamfit*icoefs[ncoefs-2]+icoefs[ncoefs-1],color=60
		;ENDCASE
		oplot,lamfit,gaussian_fits(lamfit,icoefs),color=color[4],thick=2
        ENDIF
	IF keyword_set(debug) THEN stop
	pos=[0.09,0.06,0.95,0.28]
	resid=ispec[tmp,0]-gaussian_fits(ispec[tmp,1],icoefs)
	maxpt=max(abs(resid+ispec[tmp,2])) > max(abs(resid-ispec[tmp,2]))
	plot,[0],[0],xr=u.plot.wl,yr=[-1.0*maxpt,maxpt],ytit='RESID',/xsty,/ysty,pos=pos,/noerase
	oploterror,ispec[tmp,1],resid,ispec[tmp,2],psym=8,symsize=0.25,color=color[4],errcolor=color[4]
	oplot,u.plot.wl,[0,0],linestyle=2.0
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize[1],u.plot.ysize[1],0,0,0]
	ENDIF
END

PRO wmom_plot_subspec,u
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw4,get_value=draw_win
		window,0,xsize=u.plot.xsize[3],ysize=u.plot.ysize[3],/pixmap
 		color=u.plot.col
        ENDIF ELSE BEGIN
		xsize=4.0
		ysize=xsize*u.plot.ysize[3]/u.plot.xsize[3]
		!p.thick=6
		!p.charthick=2.0
		!x.thick=5.0
		!y.thick=5.0
		device, xsize=xsize, ysize=ysize, /inches
		color=u.plot.pscol

	ENDELSE
	ispec=*u.dat.spec[u.ch-1,u.index]
	icoefs=*u.dat.coefs[u.ch-1,u.index]
	cnts=hirexsr_line_subtract(ispec[*,1],ispec[*,0],icoefs,u.dat.labels,u.line)
	tmp=where(ispec[*,1] GE u.plot.swl[0] AND ispec[*,1] LE u.plot.swl[1])
	maxpt=max(cnts[tmp]+ispec[tmp,2])*1.05
	minpt=(min(cnts[tmp]-ispec[tmp,2]) < 0)*1.05
	plot,[0],[0],xr=u.plot.swl,yr=[minpt,maxpt],ytit='B!l'+n2g('lambda')+'!n [AU]',/xsty,/ysty,xtit='Wavelength [Ang]'
	oploterror,ispec[tmp,1],cnts[tmp],ispec[tmp,2],psym=8,symsize=0.25
	oplot,u.plot.swl,[0,0],linestyle=2.0
	oplot,u.dat.lam_o*[1,1],[minpt,maxpt],color=color[1],linestyle=3.0
	loadct,39,/silent
	widget_control,u.id.lam_slider,get_value=slider
	oplot, u.dat.lam_o*[1,1]-slider/1.0e4,[minpt,maxpt],linestyle=3.0,color=color[5]
	oplot, u.dat.lam_o*[1,1]+slider/1.0e4,[minpt,maxpt],linestyle=3.0,color=color[5]
	loadct,12,/silent

	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize[3],u.plot.ysize[3],0,0,0]
	ENDIF
END

PRO wmom_load_data,u
	shot=u.shot
	line=u.line
	tht=u.tht
	IF u.stat.dat THEN BEGIN					;if LOAD_DATA already ran, clean out the heap
		cleanup_w_hirexsr_moments,u
		heap_gc
         ENDIF
	
	;specify subspectral view
	CASE u.line OF
		0 : u.plot.swl=[3.943,3.955]
		1 : u.plot.swl=[3.962,3.972]
		2 : u.plot.swl=[3.989,3.999]
                3 : u.plot.swl=[3.723,3.739]
		4 : u.plot.swl=[3.735,3.744]
		6 : u.plot.swl=[3.173,3.182]
                7 : u.plot.swl=[3.723,3.739]
		8 : u.plot.swl=[3.735,3.744]
  		9 : u.plot.swl=[3.015,3.022]
         ENDCASE
	
	;specify spectral view
	CASE u.line OF
		0 : u.plot.wl=[3.943,3.96]
		1 : u.plot.wl=[3.96,3.975]
		2 : u.plot.wl=[3.977,3.999]
                3 : u.plot.wl=[3.723,3.749]
		4 : u.plot.wl=[3.723,3.749]
		6 : u.plot.wl=[3.172,3.185]
                7 : u.plot.wl=[3.723,3.749]
		8 : u.plot.wl=[3.723,3.749]
		9 : u.plot.wl=[3.010,3.030]
        ENDCASE

	;load spectra
	widget_control,u.id.message,set_value='Loading Data for SHOT/LINE/THT = '+num2str(u.shot,1)+'/'+num2str(u.line,1)+'/'+num2str(u.tht,1),/app
	hlines=[3,4,5,6,9,10,11,12]
	IF min(where(hlines EQ line)) EQ -1 THEN h=0 ELSE h=1
	hirexsr_load_avespec,shot,specbr,lam,sig,resid,tht=tht,h=h
	avespec=hirexsr_avespec_arr2ptr(specbr,lam,sig,resid,maxave=maxave)

	;load fits
	hirexsr_load_fits,shot,line,icoefs,nave,double,labels,/quiet,status=status,tht=tht
	coefs=hirexsr_coefs_arr2ptr(icoefs)

	;load moments
	hirexsr_load_momentptr,shot,line,mom,tau,pos,tpos,lam_o,z,status=status,tht=tht,dlam=dlam,/clear
	hirexsr_load_mlintptr,shot,line,lint,tau,tht=tht,/clear
	
	;make time-evolving lint pointer
	ntau=n(tau)+1
	nch=0
	FOR i=0,n(avespec[*,0]) DO IF total(*avespec[i,0]) NE -1 THEN nch+=1
	tlint=ptrarr(nch,/allocate)
	FOR i=0,nch-1 DO BEGIN
		IF i EQ 0 THEN BEGIN
			ilint=*lint[0]
			nlint=n(ilint[0,*])+1
		ENDIF
		itlint=fltarr(ntau,nlint+1)
		FOR j=0,ntau-1 DO BEGIN
			IF tau[j] NE -1 THEN BEGIN
				ilint=*lint[j]
				imom=*mom[j]
				itlint[j,0:nlint-1]=ilint[i,*]
				itlint[j,nlint]=imom[i,24]	;add the fitcase
                        ENDIF
                ENDFOR
		*tlint[i]=itlint
 	ENDFOR

	u.index=ipt(u.time,tau)
	IF u.index EQ -1 THEN BEGIN
		tmp=where(tau NE -1)
		u.time=tau[last(tmp)]
		u.index=last(tmp)
		widget_control,u.id.t_slider,set_value=u.time*1.0e3
		widget_control,u.id.tpt,set_value=num2str(u.time,dp=2)
	ENDIF
	dat={spec:avespec,coefs:coefs,labels:labels,mom:mom,lint:lint,tlint:tlint,tau:tau,dlam:dlam,lam_o:lam_o,z:z,nch:nch}
	u={id:u.id,shot:u.shot,line:u.line,tht:u.tht,ch:u.ch,index:u.index,time:u.time,stat:u.stat,plot:u.plot,dat:dat}
	u.stat.dat=1
	widget_control,u.id.dlampt,set_value=num2str(u.dat.dlam,dp=2)
	widget_control,u.id.lam_slider,set_value=u.dat.dlam*10
	IF u.ch GT u.dat.nch THEN BEGIN
		u.ch=u.dat.nch
		widget_control,u.id.ch_slider,set_value=u.ch
		widget_control,u.id.ch,set_value=num2str(u.ch,1)
	ENDIF
	widget_control,u.id.base, set_uvalue=u	
END

PRO cleanup_w_hirexsr_moments,u
	IF u.stat.dat THEN BEGIN
		heap_free,u.dat.coefs
		heap_free,u.dat.spec
		heap_free,u.dat.mom
		heap_free,u.dat.lint
		heap_free,u.dat.tlint
	ENDIF
END

PRO wmom_reset_ytplot,u
	id=[u.id.ytplot1,u.id.ytplot2,u.id.ytplot3,u.id.ytplot4]
	FOR i=0,n(id) DO widget_control,id[i],set_button=0
	CASE u.stat.ty OF 
		0 : widget_control,u.id.ytplot1,set_button=1
		1 : widget_control,u.id.ytplot2,set_button=1
		2 : widget_control,u.id.ytplot3,set_button=1
		3 : widget_control,u.id.ytplot4,set_button=1
        ENDCASE
END

PRO wmom_reset_jtplot,u
	id=[u.id.jtplot_ch,u.id.jtplot_list]
	FOR i=0,n(id) DO widget_control,id[i],set_button=0
	CASE u.stat.tj OF
		0 : widget_control,u.id.jtplot_ch,set_button=1
		1 : widget_control,u.id.jtplot_list,set_button=1
	ENDCASE
END

PRO wmom_reset_jrho,u
	id=[u.id.jch,u.id.jrho]
	FOR i=0,n(id) DO widget_control,id[i],set_button=0
	CASE u.plot.jrho OF
		0 : widget_control,u.id.jch,set_button=1
		1 : widget_control,u.id.jrho,set_button=1
	ENDCASE
END

PRO wmom_reset_ymplot,u
	id=[u.id.ymplot1,u.id.ymplot2,u.id.ymplot3]
	FOR i=0,n(id) DO widget_control,id[i],set_button=0
	CASE u.stat.my OF 
		0 : widget_control,u.id.ymplot1,set_button=1
		1 : widget_control,u.id.ymplot2,set_button=1
		2 : widget_control,u.id.ymplot3,set_button=1
	ENDCASE	
END

PRO wmom_reset_xmplot,u
	id=[u.id.xmplot1,u.id.xmplot2,u.id.xmplot3]
	FOR i=0,n(id) DO widget_control,id[i],set_button=0
	CASE u.stat.mx OF 
		0 : widget_control,u.id.xmplot1,set_button=1
		1 : widget_control,u.id.xmplot2,set_button=1
		2 : widget_control,u.id.xmplot3,set_button=1
	ENDCASE	
END

PRO wmom_plotall,u
	wmom_plot_subspec,u
	wmom_plot_spec,u
	wmom_plot_spatial,u
	wmom_plot_time,u
END

PRO w_hirexsr_moments_event,event
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
					cleanup_w_hirexsr_moments,u
					heap_free,u.plot.list
					heap_gc
					!except=1
				
				END
				"SAVE" : BEGIN
					widget_control,u.id.message,set_value='MOMENTS widget currently only used for display',/app
				END
				"LOAD" : BEGIN
					WIDGET_CONTROL,/hourglass
					wmom_load_data,u
					wmom_plot_subspec,u
					wmom_plot_spec,u
					wmom_plot_spatial,u
					wmom_plot_time,u
				END
				"PRINT" : BEGIN
					IF u.stat.dat THEN wmoments2ps,u
					widget_control,u.id.message,set_value='Printing Display Windows to Default PS File',/app
				END
				"STOP" : BEGIN
					stop
				END
				"JCH" : BEGIN
					u.plot.jrho=0
					wmom_reset_jrho,u
					wmom_plot_time,u
				END
				"JRHO" : BEGIN
					u.plot.jrho=1
					wmom_reset_jrho,u
					wmom_plot_time,u
				END
				"JTPLOT_CH" : BEGIN
					u.stat.tj=0
					wmom_reset_jtplot,u
					wmom_plot_time,u
				END
				"JTPLOT_LIST" : BEGIN
					u.stat.tj=1
					wmom_reset_jtplot,u
					wmom_plot_time,u
				END
				"YTPLOT1" : BEGIN
					u.stat.ty=0
					wmom_reset_ytplot,u
					wmom_plot_time,u
				END
				"YTPLOT2" : BEGIN
					u.stat.ty=1
					wmom_reset_ytplot,u
					wmom_plot_time,u
                                END
				"YTPLOT3" : BEGIN
					u.stat.ty=2
					wmom_reset_ytplot,u
					wmom_plot_time,u
				END
				"YTPLOT4" : BEGIN
					u.stat.ty=3
					wmom_reset_ytplot,u
					wmom_plot_time,u
                                 END
				"YMPLOT1" : BEGIN
					u.stat.my=0
					wmom_reset_ymplot,u
					wmom_plot_spatial,u
				END
				"YMPLOT2" : BEGIN
					u.stat.my=1
					wmom_reset_ymplot,u
					wmom_plot_spatial,u
				END
				"YMPLOT3" : BEGIN
					u.stat.my=2
					wmom_reset_ymplot,u
					wmom_plot_spatial,u
				END
				"MXAUTO" : IF event.select EQ 1 THEN u.plot.mauto[0]=1 ELSE u.plot.mauto[0]=0
				"MY0AUTO" : IF event.select EQ 1 THEN u.plot.mauto[1]=1 ELSE u.plot.mauto[1]=0
				"MY1AUTO" : IF event.select EQ 1 THEN u.plot.mauto[2]=1 ELSE u.plot.mauto[2]=0
				"MY2AUTO" : IF event.select EQ 1 THEN u.plot.mauto[3]=1 ELSE u.plot.mauto[3]=0
				"XMPLOT1" : BEGIN
					u.stat.mx=0
					u.plot.mauto[0]=1
					widget_control,u.id.mxauto,set_button=1
					wmom_reset_xmplot,u
					wmom_plot_spatial,u
                                END
				"XMPLOT2" : BEGIN
					u.stat.mx=1
					u.plot.mauto[0]=1
					widget_control,u.id.mxauto,set_button=1
					wmom_reset_xmplot,u
					wmom_plot_spatial,u
				END
				"XMPLOT3" : BEGIN
					u.stat.mx=2
					u.plot.mauto[0]=1
					widget_control,u.id.mxauto,set_button=1
					wmom_reset_xmplot,u
					wmom_plot_spatial,u
                                END
			ENDCASE
		END
		"WIDGET_SLIDER": BEGIN
			widget_control,event.id,get_value=slider
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF 
				'CH_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						IF slider LE u.dat.nch THEN BEGIN
							u.ch=slider
							wmom_plotall,u
							widget_control,u.id.ch,set_value=num2str(u.ch,1)	
                                                ENDIF ELSE widget_control,u.id.ch_slider,set_value=u.ch
					ENDIF			
				END
				'T_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						index=ipt(u.dat.tau,slider/1.0e3)
						tmp=where(u.dat.tau GT 0)
						IF slider/1.0e3 GE max(u.dat.tau[tmp]) THEN index=n(tmp)
						IF slider/1.0e3 LE min(u.dat.tau[tmp]) THEN index=0
						IF u.time NE u.dat.tau[index] THEN BEGIN
							u.time=u.dat.tau[index]
							u.index=index
							widget_control,u.id.tpt,set_value=num2str(u.time,dp=2)
							wmom_plotall,u
						ENDIF
					ENDIF
				END
				"LAM_SLIDER" : BEGIN
					widget_control,u.id.dlampt,set_value=num2str(slider/10.0,dp=2)
					IF u.stat.dat THEN BEGIN
						u.stat.my=2
						wmom_plot_subspec,u
						wmom_plot_spatial,u
                                        ENDIF
				END
			ENDCASE
		END
   		"WIDGET_TEXT_CH": BEGIN
			widget_control,event.id,get_value=text
			CASE event.id OF 
			u.id.shotid : BEGIN
				widget_control,u.id.shotid,get_value=shot
				u.shot=shot
                        END
			u.id.lineid : BEGIN
				widget_control,u.id.lineid,get_value=line
				u.line=int(line)
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
			u.id.mxlow : BEGIN
				u.plot.xr[0]=float(text[0])
				wmom_plot_spatial,u
			END
			u.id.mxhigh : BEGIN
				u.plot.xr[1]=float(text[0])
				wmom_plot_spatial,u
                        END
			u.id.my0low : BEGIN
				u.plot.yr[0,0]=float(text[0])
				wmom_plot_spatial,u
			END
			u.id.my0high : BEGIN
				u.plot.yr[1,0]=float(text[0])
				wmom_plot_spatial,u
                        END
			u.id.my1low : BEGIN
				u.plot.yr[0,1]=float(text[0])
				wmom_plot_spatial,u
			END
			u.id.my1high : BEGIN
				u.plot.yr[1,1]=float(text[0])
				wmom_plot_spatial,u
                        END
			u.id.my2low : BEGIN
				u.plot.yr[0,2]=float(text[0])
				wmom_plot_spatial,u
			END
			u.id.my2high : BEGIN
				u.plot.yr[1,2]=float(text[0])
				wmom_plot_spatial,u
                        END
			u.id.xtlow : BEGIN
				u.plot.tr[0]=float(text[0])
				wmom_plot_time,u
			END
			u.id.xthigh : BEGIN
				u.plot.tr[1]=float(text[0])
				wmom_plot_time,u
			END
			u.id.jlist : BEGIN
				split=strsplit(text,',',/extract)
				*u.plot.list=float(split)
				u.stat.tj=1
				wmom_reset_jtplot,u
				wmom_plot_time,u
			END
			ENDCASE
		END
		ELSE:
	ENDCASE
	IF button NE ' QUIT ' THEN widget_control,event.top,set_uvalue=u		
END

PRO w_hirexsr_moments,shot=shot,line=line,tht=tht,load=load

	loadct,12,/silent
	mdim=get_screen_size()
	IF mdim[0] NE 1600 AND mdim[1] NE 1200 THEN base=widget_base(title='HIREXSR Moment/Line-Int Display',/row,tlb_size_events=1,/scroll,$
		x_scroll_size=mdim[0]*0.95,y_scroll_size=mdim[1]*1.5) ELSE base=widget_base(title='HIREXSR Moment/Line-Int Display',/row,tlb_size_events=1)

	A=widget_base(base,/column)
	B=widget_base(base,/column)
	C=widget_base(base,/column)

	A1=widget_base(A,frame=5)
	draw1=widget_draw(A1,xsize=450,ysize=900)

	B1=widget_base(B,frame=5)
	draw2=widget_draw(B1,xsize=620,ysize=350)
	space=widget_base(B,/row,ysize=2)
	B2=widget_base(B,frame=5)
	draw3=widget_draw(B2,xsize=620,ysize=280)
	B3=widget_base(B,/row)
	dum=widget_label(B3,value='TIME: ')
	tpt=widget_text(B3,xsize=5,ysize=1)
	t_slider=widget_slider(B3,xsize=520,min=0,max=2000,value=1000,/drag,/suppress)
	B4=widget_base(B,/row)
	dum=widget_label(B4,value='CH#:  ')
	ch=widget_text(B4,xsize=5,ysize=1)
	ch_slider=widget_slider(B4,xsize=520,min=1,max=24*4,value=16,/drag,/suppress)
	B5=widget_base(B,frame=2,/column)
	B5a=widget_base(B5,/row)
	dum = widget_label(B5a,value='SHOT:')
	shotid = widget_text(B5a,xsize=10,ysize=1,/edit)
	dum = widget_label(B5a,value=' LINE')
	lineid=widget_text(B5a,xsize=3,ysize=1,/edit)
	dum = widget_label(B5a,value=' THT')
	thtid=widget_text(B5a,xsize=3,ysize=1,/edit)
	dum = widget_label(B5a,value='    ')
	save= widget_button(B5a,value=' SAVE ')
	dum = widget_label(B5a,value=' ')
	load= widget_button(B5a,value=' LOAD ')
	dum = widget_label(B5a,value=' ')
	quit= widget_button(B5a,value=' QUIT ')
	dum = widget_label(B5a,value=' ')
	print= widget_button(B5a,value= ' PRINT ')
	dum = widget_label(B5a,value=' ')
	stop= widget_button(B5a,value=' STOP ')	
	message = widget_text(B5,xsize=97,ysize=7,/scroll,/wrap)
	
	C1=widget_base(C,frame=5)
	draw4=widget_draw(C1,xsize=300,ysize=350)
	C2=widget_base(C,/row)
	dum = widget_label(C2,value=' DLAM ')
	lam_slider = widget_slider(C2,xsize=160,min=1,max=100,value=40,/drag,/suppress)
	dlampt=widget_text(C2,xsize=3,ysize=1)
	dum = widget_label(C2,value=' [mAng] ')

	C3=widget_base(C,frame=2,/column,/align_center)
	dum = widget_label(C3,value='Time History Plotting')
	dum = widget_label(C3,value='==========================================')
	C3a=widget_base(C3,/row)
	dum = widget_label(C3a,value='PLOT: ')
	C3ax=widget_base(C3a,/row,/nonexcl)
	ytplot1=widget_button(C3ax,value='Br ')
	ytplot2=widget_button(C3ax,value='Vtor ')
	ytplot3=widget_button(C3ax,value='Vpol ')
	ytplot4=widget_button(C3ax,value='Ti ')
	C3b=widget_base(C3,/row)
	dum=widget_label(C3b, value=' USE: ')
	C3bx=widget_base(C3b,/row,/nonexclu)
	jtplot_ch=widget_button(C3bx,value='Selected CH')
	jtplot_list=widget_button(C3bx,value='CH List (below)')
	C3c=widget_base(C3,/row)
	dum = widget_label(C3c,value='       ')
	jlist=widget_text(C3c,xsize=16,ysize=1,/edit)
	C3cx=widget_base(C3c,/row,/nonexcl)
	jch=widget_button(C3cx,value='CH')
	jrho=widget_button(C3cx,value='r/a')
	C3d=widget_base(C3,/row)
	dum = widget_label(C3d,value='TRANGE: ')
	xtlow=widget_text(C3d,xsize=4,ysize=1,/edit)
	dum = widget_label(C3d,value='to')
	xthigh=widget_text(C3d,xsize=4,ysize=1,/edit)
	dum = widget_label(C3d,value=' [sec]')

	C4=widget_base(C,frame=2,/column,/align_center)
	dum = widget_label(C4,value='           Moment Profile Plotting          ')
	dum = widget_label(C4,value='==========================================')
	C4a=widget_base(C4,/row)
	dum = widget_label(C4a,value='PLOT: ')
	C4ax=widget_base(C4a,/row,/nonexcl)
	ymplot1=widget_button(C4ax,value='MOM ')
	ymplot2=widget_button(C4ax,value='LINT ')
	ymplot3=widget_button(C4ax,value='MCONV ')
	C4b=widget_base(C4,/row)
	dum = widget_label(C4b,value='  Vs: ')
	C4bx=widget_base(C4b,/row,/nonexcl)
	xmplot1=widget_button(C4bx,value='CH# ')
	xmplot2=widget_button(C4bx,value='r/a ')
	xmplot3=widget_button(C4bx,value='PSIN ')
	C4c=widget_base(C4,/row)
	dum = widget_label(C4c,value='      XRANGE: ')
	mxlow=widget_text(C4c,xsize=4,ysize=1,/edit)
	dum = widget_label(C4c,value='to')
	mxhigh=widget_text(C4c,xsize=4,ysize=1,/edit)
	C4cx=widget_base(C4c,/row,/nonexcl)
	mxauto=widget_button(C4cx,value='AUTO')
	C4d=widget_base(C4,/row)
	dum = widget_label(C4d,value='YRANGE (0th): ')
	my0low=widget_text(C4d,xsize=4,ysize=1,/edit)
	dum = widget_label(C4d,value='to')
	my0high=widget_text(C4d,xsize=4,ysize=1,/edit)
	C4dx=widget_base(C4d,/row,/nonexcl)
	my0auto=widget_button(C4dx,value='AUTO')

	C4e=widget_base(C4,/row)
	dum = widget_label(C4e,value='YRANGE (1st): ')
	my1low=widget_text(C4e,xsize=4,ysize=1,/edit)
	dum = widget_label(C4e,value='to')
	my1high=widget_text(C4e,xsize=4,ysize=1,/edit)	
	C4ex=widget_base(C4e,/row,/nonexcl)
	my1auto=widget_button(C4ex,value='AUTO')

	C4f=widget_base(C4,/row)
	dum = widget_label(C4f,value='YRANGE (2nd): ')
	my2low=widget_text(C4f,xsize=4,ysize=1,/edit)
	dum = widget_label(C4f,value='to')
	my2high=widget_text(C4f,xsize=4,ysize=1,/edit)		
	C4fx=widget_base(C4f,/row,/nonexcl)
	my2auto=widget_button(C4fx,value='AUTO')

	id={base:base,draw1:draw1,draw2:draw2,draw3:draw3,draw4:draw4,$
		tpt:tpt,t_slider:t_slider,ch:ch,ch_slider:ch_slider,$
		shotid:shotid,lineid:lineid,thtid:thtid,$
		save:save,load:load,quit:quit,print:print,stop:stop,message:message,$
		lam_slider:lam_slider,dlampt:dlampt,$
		ytplot1:ytplot1,ytplot2:ytplot2,ytplot3:ytplot3,ytplot4:ytplot4,xtlow:xtlow,xthigh:xthigh,$
		jtplot_ch:jtplot_ch,jtplot_list:jtplot_list,jlist:jlist,jch:jch,jrho:jrho,$
		ymplot1:ymplot1,ymplot2:ymplot2,ymplot3:ymplot3,$
		xmplot1:xmplot1,xmplot2:xmplot2,xmplot3:xmplot3,$
		mxlow:mxlow,mxhigh:mxhigh,mxauto:mxauto,my0low:my0low,my0high:my0high,my0auto:my0auto,$
		my1low:my1low,my1high:my1high,my1auto:my1auto,my2low:my2low,my2high:my2high,my2auto:my2auto}

	IF NOT keyword_set(shot) THEN shot=1120224017
	IF NOT keyword_set(line) THEN line=0
	IF NOT keyword_set(tht) THEN tht=0

	xr=[0,0.9]
	tr=[0.0,2.0]
	yr=fltarr(2,3)
	plot={xsize:[450,620,620,300],ysize:[900,350,280,350],col:[255,50,80,150,200,190],pscol:[0,30,100,200,200,200],jrho:1,mauto:[10,1,1,1],$
		xr:xr,tr:tr,yr:yr,list:ptr_new([0],/allocate_heap),wl:fltarr(2),swl:fltarr(2)}
	stat={dat:0,ps:0,ty:0,tj:0,my:1,mx:1}
	u={id:id,shot:shot,line:line,tht:tht,ch:10,index:0,time:1.0,stat:stat,plot:plot}
	widget_control,base,set_uvalue=u
	widget_control,u.id.shotid,set_value=num2str(u.shot,1)
	widget_control,u.id.lineid,set_value=num2str(u.line,1)
	widget_control,u.id.thtid,set_value=num2str(u.tht,1)
	widget_control,u.id.t_slider,set_value=u.time*1.0e3
	widget_control,u.id.tpt,set_value=num2str(u.time,dp=2)
	widget_control,u.id.ch_slider,set_value=u.ch
	widget_control,u.id.ch,set_value=num2str(u.ch,1)
	widget_control,u.id.xtlow,set_value=num2str(u.plot.tr[0],dp=2)
	widget_control,u.id.xthigh,set_value=num2str(u.plot.tr[1],dp=2)
	widget_control,u.id.mxlow,set_value=num2str(u.plot.xr[0],dp=2)
	widget_control,u.id.mxhigh,set_value=num2str(u.plot.xr[1],dp=2)
	CASE u.stat.ty OF 
		0 : widget_control,u.id.ytplot1,set_button=1
		1 : widget_control,u.id.ytplot2,set_button=1
		2 : widget_control,u.id.ytplot3,set_button=1
		3 : widget_control,u.id.ytplot4,set_button=1
        ENDCASE
	CASE u.stat.tj OF
		0 : widget_control,u.id.jtplot_ch,set_button=1
		1 : widget_control,u.id.jtplot_list,set_button=1
	ENDCASE
	IF u.plot.jrho THEN widget_control,u.id.jrho,set_button=1 ELSE widget_control,u.id.jch,set_button=1
	CASE u.stat.my OF 
		0 : widget_control,u.id.ymplot1,set_button=1
		1 : widget_control,u.id.ymplot2,set_button=1
		2 : widget_control,u.id.ymplot3,set_button=1
	ENDCASE	
	CASE u.stat.mx OF 
		0 : widget_control,u.id.xmplot1,set_button=1
		1 : widget_control,u.id.xmplot2,set_button=1
		2 : widget_control,u.id.xmplot3,set_button=1
	ENDCASE	
	*u.plot.list=[14,12,10,8]+1
	widget_control,u.id.jlist,set_value='14,12,10,8'
	IF u.plot.mauto[0] THEN widget_control,u.id.mxauto,set_button=1
	IF u.plot.mauto[1] THEN widget_control,u.id.my0auto,set_button=1
	IF u.plot.mauto[2] THEN widget_control,u.id.my1auto,set_button=1
	IF u.plot.mauto[3] THEN widget_control,u.id.my2auto,set_button=1

	!except=0			
	widget_control,base,/realize
	IF keyword_set(load) THEN BEGIN
		WIDGET_CONTROL,/hourglass
		wmom_load_data,u
		wmom_plot_subspec,u
		wmom_plot_spec,u
		wmom_plot_spatial,u
		wmom_plot_time,u
	ENDIF
	xmanager,'w_hirexsr_moments',base
END
