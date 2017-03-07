FUNCTION hirexsr_bsfit_bspline,fit,order=order,nknots=nknots,iterfit=iterfit,conft=conft,bfit=bfit,resid=resid,xfit=xfit,yfit=yfit,ypfit=ypfit,$
		good=good,dtifit=dtifit,fitcase=fitcase

	IF NOT keyword_set(order) THEN order=4
	IF NOT keyword_set(nknots) THEN nknots=1
	IF NOT keyword_set(good) THEN good=int(fit.rho*0.0+1.0)
	fitcase=1
	IF keyword_set(iterfit) THEN fitcase=0
	IF keyword_set(conft) THEN fitcase=2
	IF nknots LT 0 THEN fitcase=-1
	tmp=where(good EQ 1)
	CASE fitcase OF
		-1 : BEGIN
			print, 'NKNOTS too small, must be GE 0'
		END
		0 : BEGIN
			bfit=bspline_iterfit(fit.rho[tmp],fit.(0)[tmp],nord=order,bkspace=5)
			bsfit=bspline_valu(fit.rho,bfit)
			IF keyword_set(xfit) THEN BEGIN
				yfit=bspline_valu(xfit,bfit)
				ypfit=deriv(xfit,yfit)
			ENDIF
			dbsfit=deriv(fit.rho,tifit)
                END

		1 : BEGIN
			xv=make(0.0,last(fit.rho[tmp])*1.0001,nknots+order)
			bsnak,xv,order,xknot
			bslsq,float(fit.rho[tmp]),float(fit.(0)[tmp]),1/float(fit.err[tmp]),order,xknot,nknots+order,bscoef
			bsfit=bs1gd(0,float(fit.rho),order,xknot,bscoef)
			dbsfit=bs1gd(1,float(fit.rho),order,xknot,bscoef)
			IF keyword_set(xfit) THEN BEGIN
				yfit=bs1gd(0,xfit,order,xknot,bscoef)
				ypfit=bs1gd(1,xfit,order,xknot,bscoef)
			ENDIF
			bfit={xknot:xknot,bscoef:bscoef}
                     END

		2 : BEGIN
			xv=make(0.0,last(fit.rho),nknots+order)
			bsnak,xv,order,xknot
			conft,float(fit.rho[tmp]),float(fit.(0)[tmp]),1/float(fit.err[tmp]),[0.0],nHard,[1.0],[1],[0.0],[0.0],order,xknot,nknots+order,bscoef
			bsfit=bs1gd(0,float(fit.rho),order,xknot,bscoef)
			dbsfit=bs1gd(1,float(fit.rho),order,xknot,bscoef)
			IF keyword_set(xfit) THEN BEGIN
				yfit=bs1gd(0,xfit,order,xknot,bscoef)
				ypfit=bs1gd(1,xfit,order,xknot,bscoef)
			ENDIF
			bfit={xknot:xknot,bscoef:bscoef}
		END		
        ENDCASE
	resid=fit.(0)-bsfit	

	RETURN,bsfit
END

PRO hirexsr_run_bsfit,dat,bsfit,ntrials=ntrials,nsigma=nsigma,iterfit=iterfit,conft=conft,order=order,nknots=nknots,nrho=nrho,fitcase=fitcase
	IF NOT keyword_set(nsigma) THEN nsigma=3.0
	IF NOT keyword_set(ntrials) THEN ntrials=500
	IF NOT keyword_set(nrho) THEN nrho=50
	bsdat=fltarr(nrho,ntrials)
	dbsdat=fltarr(nrho,ntrials)
	bschk=hirexsr_bsfit_bspline(dat,order=order,nknots=nknots,iterfit=iterfit,conft=conft,resid=resid,good=good,fitcase=fitcase)
	tmp=where(abs(resid/dat.err) GT nsigma)	
	IF tmp[0] NE -1 THEN good[tmp]=0
	xfit=make(0.0,max(dat.rho[where(good EQ 1)]),nrho)
	FOR i=0,ntrials-1 DO BEGIN
		idat=dat
		idat.(0)+=randomn(seed,nrho)*dat.err
		out=hirexsr_bsfit_bspline(idat,order=order,nknots=nknots,iterfit=iterfit,conft=conft,good=good,xfit=xfit,yfit=yfit,ypfit=ypfit,bfit=bfit)
		bsdat[*,i]=yfit
		dbsdat[*,i]=ypfit
	ENDFOR
	prof=fltarr(nrho)
	dprof=fltarr(nrho)
	err=fltarr(nrho)
	derr=fltarr(nrho)
	FOR i=0,nrho-1 DO BEGIN
		prof[i]=mean(bsdat[i,*])
		err[i]=stdev(bsdat[i,*])
		dprof[i]=mean(dbsdat[i,*])
		derr[i]=stdev(dbsdat[i,*])
	ENDFOR
	bsfit={prof:prof,rho:xfit,err:err,dprof:dprof,derr:derr,good:good}
END

PRO hirexsr_bsfit_plotbsti,fit,bsti
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	d_old=!d
	ls=1.0
	IF keyword_set(ps) THEN BEGIN
		device, xsize=7.5,ysize=7.5*6.0/8.0,/inches
		ls=1.0
		!p.thick=6
		!p.charthick=2.0
	ENDIF ELSE openwin,0,xsize=800,ysize=600
	type=[1,2,3,4]
	color=[100,120,150,200]

	!p.multi=[0,0,2]
	xr=[0,1]
	ymax=max(fit.ti+fit.err) > max(bsti.prof+bsti.err)
	yr=[0,ymax]*1.05
	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit=fit.rlab,ytit='Ion Temp. [keV]'
	FOR i=0,n(type) DO BEGIN
		tmp=where(bsti.good EQ 1 AND fit.type EQ type[i])
		IF tmp[0] NE -1 THEN oploterror,fit.rho[tmp],fit.ti[tmp],fit.err[tmp],psym=3,color=color[i],errcolor=color[i]
		tmp=where(bsti.good EQ 0 AND fit.type EQ type[i])
		IF tmp[0] NE -1 THEN oplot,fit.rho[tmp],fit.ti[tmp],psym=7,color=color[i]
	ENDFOR
	oploterror,bsti.rho,bsti.prof,bsti.err,color=30,errcolor=30,thick=2.0

	lti=-1.0*bsti.dprof/bsti.prof
	lerr=lti*sqrt((bsti.derr/bsti.dprof)^2+(bsti.err/bsti.prof)^2)
	xr=[0,1]
	tmp=where(bsti.rho LT 0.8)
	ymax=max(lti[tmp]+lerr[tmp])
	yr=[0,ymax]*1.05
	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit=fit.rlab,ytit='a/LTi'
	oploterror,bsti.rho,lti,lerr,color=30,errcolor=30,thick=2.0
	
	!p.multi=0
	
END

PRO hirexsr_bsfit_plotbsom,fit,bsom
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	d_old=!d
	ls=1.0
	IF keyword_set(ps) THEN BEGIN
		device, xsize=7.5,ysize=7.5*6.0/8.0,/inches
		ls=1.0
		!p.thick=6
		!p.charthick=2.0
	ENDIF ELSE openwin,0,xsize=800,ysize=600
	type=[1,2,3,4]
	color=[100,120,150,200]

	!p.multi=[0,0,2]
	xr=[0,1]
	ymax=max(fit.om+fit.err) > max(bsom.prof+bsom.err) > 0
	ymin=min(fit.om-fit.err) < max(bsom.prof-bsom.err) < 0
	yr=[ymin,ymax]*1.05
	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit=fit.rlab,ytit='Tor. Rotation [kHz]'
	FOR i=0,n(type) DO BEGIN
		tmp=where(bsom.good EQ 1 AND fit.type EQ type[i])
		IF tmp[0] NE -1 THEN oploterror,fit.rho[tmp],fit.om[tmp],fit.err[tmp],psym=3,color=color[i],errcolor=color[i]
		tmp=where(bsom.good EQ 0 AND fit.type EQ type[i])
		IF tmp[0] NE -1 THEN oplot,fit.rho[tmp],fit.om[tmp],psym=7,color=color[i]
	ENDFOR
	oploterror,bsom.rho,bsom.prof,bsom.err,color=30,errcolor=30,thick=2.0
	oplot,xr,[0,0],linestyle=1

	xr=[0,1]
	tmp=where(bsom.rho LT 0.8)
	ymax=max(bsom.dprof[tmp]+bsom.derr[tmp]) > 0
	ymin=min(bsom.dprof[tmp]-bsom.derr[tmp]) < 0

	yr=[ymin,ymax]*1.05
	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit=fit.rlab,ytit='dom/dr'
	oploterror,bsom.rho,bsom.dprof,bsom.derr,color=30,errcolor=30,thick=2.0
	oplot,xr,[0,0],linestyle=1
	
	!p.multi=0
	
END

PRO hirexsr_tifit_input,adat,bdat,fitdat,limits=limits
	IF NOT keyword_set(limits) THEN limits={arho:[0.2,0.9],brho:[0.0,0.6],alrho:[0.8,0.95],blrho:[0.5,0.7],ymin:0.0,ymax:10.0,errmax:0.5}
	tmp=where(adat.rho GE limits.arho[0] AND adat.rho LE limits.arho[1])
	ti=adat.ti[tmp]
	rho=adat.rho[tmp]
 	err=adat.err[tmp]
	type=tmp*0.0+1.0
	tmp=where(adat.lrho GE limits.alrho[0] AND adat.lrho LE limits.alrho[1])
	IF tmp[0] NE -1 THEN BEGIN
		ti=[ti,adat.lti[tmp]]	
		rho=[rho,adat.lrho[tmp]]
		err=[err,adat.lerr[tmp]]
		type=[type,tmp*0.0+2.0]
        ENDIF
	IF size(bdat,/type) EQ 8 THEN BEGIN
		tmp=where(bdat.rho GE limits.brho[0] AND bdat.rho LE limits.brho[1])
		ti=[ti,bdat.ti[tmp]]
		rho=[rho,bdat.rho[tmp]]
		err=[err,bdat.err[tmp]]
		type=[type,tmp*0.0+3.0]
		tmp=where(bdat.lrho GE limits.blrho[0] AND bdat.lrho LE limits.blrho[1])
		IF tmp[0] NE -1 THEN BEGIN
			ti=[ti,bdat.lti[tmp]]
			rho=[rho,bdat.lrho[tmp]]
			err=[err,bdat.lerr[tmp]]
			type=[type,tmp*0.0+4.0]
	        ENDIF	
	ENDIF
	order=sort(rho)
	ti=ti[order]
	err=err[order]
	rho=rho[order]
	type=type[order]
	tmp=where(err LT limits.errmax AND ti LT limits.ymax  AND ti GT limits.ymin)
	fitdat={ti:ti[tmp],err:err[tmp],rho:rho[tmp],type:type[tmp],rlab:adat.rlab}
END

PRO hirexsr_omfit_input,adat,bdat,fitdat,limits=limits
	IF NOT keyword_set(limits) THEN limits={arho:[1.0,1.0],brho:[1.0,1.0],alrho:[0.0,0.9],blrho:[1.0,1.0],ymin:-20.0,ymax:20.0,errmax:4.0}
	tmp=where(adat.lrho GE limits.alrho[0] AND adat.lrho LE limits.alrho[1])
	om=adat.lom[tmp]
	rho=adat.lrho[tmp]
 	err=adat.lerr[tmp]
	type=tmp*0.0+2.0
	tmp=where(adat.rho GE limits.arho[0] AND adat.rho LE limits.arho[1])
	IF tmp[0] NE -1 THEN BEGIN
		om=[om,adat.om[tmp]]
		rho=[rho,adat.rho[tmp]]
		err=[err,adat.err[tmp]]
		type=[type,tmp*0.0+1.0]
        ENDIF
	IF size(bdat,/type) EQ 8 THEN BEGIN
		tmp=where(bdat.rho GE limits.brho[0] AND bdat.rho LE limits.brho[1])
		IF tmp[0] NE -1 THEN BEGIN
			om=[om,bdat.om[tmp]]
			rho=[rho,bdat.rho[tmp]]
			err=[err,bdat.err[tmp]]
			type=[type,tmp*0.0+3.0]
             	ENDIF
		tmp=where(bdat.lrho GE limits.blrho[0] AND bdat.lrho LE limits.blrho[1])
		IF tmp[0] NE -1 THEN BEGIN
			om=[om,bdat.lom[tmp]]
			rho=[rho,bdat.lrho[tmp]]
			err=[err,bdat.lerr[tmp]]
			type=[type,tmp*0.0+4.0]
	        ENDIF	
	ENDIF
	order=sort(rho)
	om=om[order]
	err=err[order]
	rho=rho[order]
	type=type[order]
	tmp=where(err LT limits.errmax AND om LT limits.ymax  AND om GT limits.ymin)
	fitdat={om:om[tmp],err:err[tmp],rho:rho[tmp],type:type[tmp],rlab:adat.rlab}
END

PRO hirexsr_tifit_load,shot,line,tidat,time,tht=tht,rho=rho,etree=etree,m1err=m1merr,inst=inst
	hirexsr_load_mlintptr,shot,line,mlint,tau,tht=tht
	hirexsr_load_profile,shot,line,prof,perr,psin,tau,tht=tht,tinst=tinst
	IF keyword_set(inst) THEN dc=tinst ELSE dc=0.0
	
	tmp=where(tau NE -1)
	time=tau[tmp]
	ntime=n(time)+1
	
	IF NOT keyword_set(etree) THEN etree='ANALYSIS'
	mdsopen,etree,shot
	rmid=mdsvalue('_dat=\efit_rmid')
	etime=mdsvalue('dim_of(_dat,0)')
	epsin=mdsvalue('dim_of(_dat,1)')
	mdsclose,etree,shot

	tidat=ptrarr(ntime,/allocate_heap)
	FOR i=0,ntime-1 DO BEGIN
		tmp=where(psin[*,i] NE -1)
		ipsin=psin[tmp,i]
		nr=n(ipsin)+1
		IF ipsin[0] EQ ipsin[nr/2] THEN BEGIN
			nr/=2
			ipsin=psin[0:nr-1,i]			
			ti=prof[0:nr-1,i,3]
			m1=prof[nr:*,i,3]
			err=perr[0:nr-1,i,3]
                ENDIF ELSE BEGIN
			ti=prof[tmp,i,3]
			m1=ti*0.0
			err=perr[tmp,i,3]
                ENDELSE
		IF keyword_set(m1err) THEN err=sqrt(m1^2+err^2)
		IF keyword_set(rho) THEN BEGIN
			irmid=reform(rmid[ipt(etime,time[i]),*])
			irho=(irmid-irmid[0])/(last(irmid)-irmid[0])
			xrho=interpol(irho,epsin,ipsin)
			xlab='r/a'
                ENDIF ELSE BEGIN
			xrho=ipsin
			xlab='PSIN'
		ENDELSE
		ilint=*mlint[i]
		lti=ilint[*,2]
		lerr=ilint[*,3]
		IF keyword_set(rho) THEN BEGIN
			lrho=interpol(irho,epsin,ilint[*,4])
                ENDIF ELSE lrho=ilint[*,4]
		itidat={ti:ti-dc,err:err,rho:xrho,lti:lti-dc,lerr:lerr,lrho:lrho,shot:shot,tht:tht,time:time[i],rlab:xlab}
		*tidat[i]=itidat
	ENDFOR

	heap_free,mlint
END

PRO hirexsr_omfit_load,shot,line,omdat,time,tht=tht,rho=rho,etree=etree,m1err=m1merr
	hirexsr_load_mlintptr,shot,line,mlint,tau,tht=tht
	hirexsr_load_profile,shot,line,prof,perr,psin,tau,tht=tht
	
	tmp=where(tau NE -1)
	time=tau[tmp]
	ntime=n(time)+1
	
	IF NOT keyword_set(etree) THEN etree='ANALYSIS'
	mdsopen,etree,shot
	rmid=mdsvalue('_dat=\efit_rmid')
	etime=mdsvalue('dim_of(_dat,0)')
	epsin=mdsvalue('dim_of(_dat,1)')
	mdsclose,etree,shot

	omdat=ptrarr(ntime,/allocate_heap)
	FOR i=0,ntime-1 DO BEGIN
		tmp=where(psin[*,i] NE -1)
		ipsin=psin[tmp,i]
		nr=n(ipsin)+1
		IF ipsin[0] EQ ipsin[nr/2] THEN BEGIN
			nr/=2
			ipsin=psin[0:nr-1,i]			
			om=prof[0:nr-1,i,1]
			m1=prof[nr:*,i,1]
			err=perr[0:nr-1,i,1]
                ENDIF ELSE BEGIN
			om=prof[tmp,i,1]
			m1=om*0.0
			err=perr[tmp,i,1]
                ENDELSE
		IF keyword_set(m1err) THEN err=sqrt(m1^2+err^2)
		IF keyword_set(rho) THEN BEGIN
			irmid=reform(rmid[ipt(etime,time[i]),*])
			irho=(irmid-irmid[0])/(last(irmid)-irmid[0])
			xrho=interpol(irho,epsin,ipsin)
			xlab='r/a'
                ENDIF ELSE BEGIN
			xrho=ipsin
			xlab='PSIN'
		ENDELSE
		ilint=*mlint[i]
		lom=ilint[*,0]
		lerr=ilint[*,1]
		IF keyword_set(rho) THEN BEGIN
			lrho=interpol(irho,epsin,ilint[*,4])
                ENDIF ELSE lrho=ilint[*,4]
		iomdat={om:om,err:err,rho:xrho,lom:lom,lerr:lerr,lrho:lrho,shot:shot,tht:tht,time:time[i],rlab:xlab}
		*omdat[i]=iomdat
	ENDFOR

	heap_free,mlint
END

PRO hirexsr_bsfit_adjerror,idat,perr,lerr
	pcase=0
	IF perr LT 0 THEN pcase=1
	IF perr GT 0 THEN pcase=2
	CASE pcase OF
		1 : BEGIN					;set minimum error in absolute temperature
			echk=where(idat.err LT -1.0*perr)
			IF echk[0] NE -1 THEN idat.err[echk]=-1.0*perr			
		END
 		2 : BEGIN					;set minimum error as fractional error
			echk=where(abs(idat.err/idat.(0)) LT perr)
			IF echk[0] NE -1 THEN idat.err[echk]=abs(perr*idat.(0)[echk])
		END
		ELSE :  
        ENDCASE
	lcase=0
	IF lerr LT 0 THEN lcase=1
	IF lerr GT 0 THEN lcase=2
	CASE lcase OF
		1 : BEGIN					;set minimum error in absolute temperature
			echk=where(idat.lerr LT -1.0*lerr)
			IF echk[0] NE -1 THEN idat.lerr[echk]=-1.0*lerr			
		END
 		2 : BEGIN					;set minimum error as fractional error
			echk=where(abs(idat.lerr/idat.(2)) LT lerr)
			IF echk[0] NE -1 THEN idat.lerr[echk]=abs(lerr*idat.(2)[echk])
		END
		ELSE :  
        ENDCASE
END

;IF bin=-1 then use all in tr, if bin=0, fit each time slice,
;otherwise collect every BIN # of frames together
PRO hirexsr_tifit,shot,tifit,rho,time,err,bsti,tht=tht,plot=plot,nrho=nrho,bin=bin,tr=tr,emin=emin,m1err=m1err,conft=conft,iterfit=iterfit,order=order,nknots=nknots,$
		nobranchb=nobranchb,ptr=ptr,limits=limits,error=error,config=config,inst=inst
	IF NOT keyword_set(bin) THEN bin=0
	IF NOT keyword_set(emin) THEN emin=0.05
	IF NOT keyword_set(error) THEN error={a:emin,b:emin,al:0.0,bl:0.0}
	IF NOT keyword_set(inst) THEN inst=0
	IF shot GT 1120821000 AND shot LT 1121002000 THEN line=7 ELSE line=2
	IF line EQ 7 OR keyword_set(nobranchb) THEN nob=1 ELSE nob=0
	IF NOT keyword_set(limits) THEN limits={arho:[0.2,0.9],brho:[0.0,0.6],alrho:[0.8,0.95],blrho:[0.5,0.7],ymin:0.0,ymax:10.0,errmax:0.5}
	hirexsr_tifit_load,shot,line,adat,atime,tht=tht,/rho,m1err=m1err,inst=inst
 	IF NOT nob THEN BEGIN
		hirexsr_tifit_load,shot,line+1,bdat,btime,tht=tht,/rho,inst=inst
		IF total(atime-btime) NE 0 THEN BEGIN
			print, 'ERROR: Branch A and Branch B timebase differet'
			RETURN
		ENDIF
        ENDIF ELSE BEGIN
		bdat=-1
		limits.arho=[0,0.9]
	ENDELSE
	time=atime
	IF NOT keyword_set(tr) THEN tr=[time[0],last(time)]
	bcase=0
	IF bin GT 1 THEN bcase=1
	IF bin LT 0 THEN bcase=2

	CASE bcase OF
		0 : BEGIN
			tmp=where(time GE tr[0] AND time LE tr[1])
			time=time[tmp]
			ntime=n(tmp)+1
			IF keyword_set(ptr) THEN bsti=ptrarr(ntime,/allocate_heap) ELSE bsti=0		
			FOR j=0,n(tmp) DO BEGIN
				i=tmp[j]
				ia=*adat[i]
				hirexsr_bsfit_adjerror,ia,error.a,error.al
				IF NOT nob THEN BEGIN
					ib=*bdat[i]
					hirexsr_bsfit_adjerror,ib,error.b,error.bl
                                ENDIF ELSE BEGIN
					ib=-1
				ENDELSE
				hirexsr_tifit_input,ia,ib,idat,limits=limits
				hirexsr_run_bsfit,idat,ibsti,conft=conft,iterfit=iterfit,order=order,nknots=nknots,nsigma=nsigma,nrho=nrho,ntrials=ntrials,fitcase=fitcase
				IF j EQ 0 THEN BEGIN
					tifit=fltarr(nrho,ntime)
					err=fltarr(nrho,ntime)
					rho=fltarr(nrho,ntime)
				ENDIF
				tifit[*,i]=ibsti.(0)
				rho[*,i]=ibsti.rho
				err[*,i]=ibsti.err
				config={order:order,nknots:nknots,ntrials:ntrials,nsigma:nsigma,nrho:nrho,fitcase:fitcase}
				IF keyword_set(ptr) THEN *bsti[i]={fit:ibsti,dat:idat,limits:limits,error:error,config:config,shot:shot,tht:tht,time:time[i]}
				IF keyword_set(plot) THEN hirexsr_bsfit_plotbsti,idat,ibsti
			ENDFOR
		END
	ENDCASE
	heap_free,adat
	heap_free,bdat

	;stop
END

PRO hirexsr_omfit,shot,omfit,rho,time,err,bsom,tht=tht,plot=plot,nrho=nrho,bin=bin,tr=tr,emin=emin,m1err=m1err,conft=conft,iterfit=iterfit,order=order,nknots=nknots,$
		nobranchb=nobranchb,ptr=ptr,limits=limits,error=error,config=config
	IF NOT keyword_set(bin) THEN bin=0
	IF NOT keyword_set(emin) THEN emin=0.05
	IF NOT keyword_set(error) THEN error={a:emin,b:emin,al:0.0,bl:0.0}
	IF shot GT 1120821000 AND shot LT 1121002000 THEN line=7 ELSE line=2
	IF line EQ 7 OR keyword_set(nobranchb) THEN nob=1 ELSE nob=0
	IF NOT keyword_set(limits) THEN limits={arho:[1.0,1.0],brho:[1.0,1.0],alrho:[0.0,0.9],blrho:[1.0,1.0],ymin:-20.0,ymax:20.0,errmax:10.0}
	hirexsr_omfit_load,shot,line,adat,atime,tht=tht,/rho,m1err=m1err
 	IF NOT nob THEN BEGIN
		hirexsr_omfit_load,shot,line+1,bdat,btime,tht=tht,/rho
		IF total(atime-btime) NE 0 THEN BEGIN
			print, 'ERROR: Branch A and Branch B timebase differet'
			RETURN
		ENDIF
        ENDIF ELSE BEGIN
		bdat=-1
	ENDELSE
	time=atime
	IF NOT keyword_set(tr) THEN tr=[time[0],last(time)]
	bcase=0
	IF bin GT 1 THEN bcase=1
	IF bin LT 0 THEN bcase=2

	CASE bcase OF
		0 : BEGIN
			tmp=where(time GE tr[0] AND time LE tr[1])
			time=time[tmp]
			ntime=n(tmp)+1
			IF keyword_set(ptr) THEN bsom=ptrarr(ntime,/allocate_heap) ELSE bsom=0		
			FOR j=0,n(tmp) DO BEGIN
				i=tmp[j]
				ia=*adat[i]
				hirexsr_bsfit_adjerror,ia,error.a,error.al
				IF NOT nob THEN BEGIN
					ib=*bdat[i]
					hirexsr_bsfit_adjerror,ib,error.b,error.bl
                                ENDIF ELSE BEGIN
					ib=-1
				ENDELSE
				hirexsr_omfit_input,ia,ib,idat,limits=limits
				hirexsr_run_bsfit,idat,ibsom,conft=conft,iterfit=iterfit,order=order,nknots=nknots,nsigma=nsigma,nrho=nrho,ntrials=ntrials,fitcase=fitcase
				IF j EQ 0 THEN BEGIN
					omfit=fltarr(nrho,ntime)
					err=fltarr(nrho,ntime)
					rho=fltarr(nrho,ntime)
				ENDIF
				omfit[*,i]=ibsom.(0)
				rho[*,i]=ibsom.rho
				err[*,i]=ibsom.err
				config={order:order,nknots:nknots,ntrials:ntrials,nsigma:nsigma,nrho:nrho,fitcase:fitcase}
				IF keyword_set(ptr) THEN *bsom[i]={fit:ibsom,dat:idat,limits:limits,error:error,config:config,shot:shot,tht:tht,time:time[i]}
				IF keyword_set(plot) THEN hirexsr_bsfit_plotbsom,idat,ibsom
			ENDFOR
		END
	ENDCASE
	heap_free,adat
	heap_free,bdat

	;stop
END
