PRO hirexsr_vpol_plot,shot,time,line=line,tht=tht,rho=rho,rtmax=rtmax,tir=tir
	IF NOT keyword_set(drmid) THEN drmid=0
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	IF NOT keyword_set(line) THEN line=2

    	hirexsr_load_mlintptr,shot,line,lint,tau,tht=tht		;ilint=[[v],[verr],[ti],[tierr],[rhotang],[rmidtang]]
	mdsopen,'analysis',shot
	rmid=mdsvalue('\efit_rmid')
	etime=mdsvalue('dim_of(\efit_rmid,0)')
	mdsclose,'analysis',shot
	rmid=reform(rmid[ipt(etime,time),*])

	IF keyword_set(ps) THEN BEGIN
		xsize=7.75
		ysize=7.75*1100/1200.0
		ls=0.65
	ENDIF ELSE BEGIN
		xsize=1200.0
		ysize=1100.0
		ls=1.0
	ENDELSE
	IF NOT keyword_set(ps) THEN BEGIN
		device, window_state=var
		IF var[1] EQ 0 THEN window,1,xsize=xsize,ysize=ysize,xpos=2340,ypos=670,$
			title='BRCHK and Residual' ELSE wset,1
	ENDIF ELSE BEGIN
		device, xsize=xsize, ysize=ysize, /inches
	ENDELSE
	IF NOT keyword_set(rtmax) THEN rtmax=0.86
	
	

	index=ipt(tau,time)
	ilint=*lint[index]
	heap_free,lint

	j=minloc(ilint[*,4])
	simp_v_low=ilint[0:j,0]
	v_err_low=ilint[0:j,1]
	simp_v_high=ilint[j:*,0]
	v_err_high=ilint[j:*,1]
	IF NOT keyword_set(rho) THEN BEGIN
		rmid_low=ilint[0:j,5]
		rmid_high=ilint[j:*,5]
		xtit='RMID [m]'
		xr=[0.68,0.91]
	ENDIF ELSE BEGIN
		rmid_low=(ilint[0:j,5]-rmid[0])/(last(rmid)-rmid[0])
		rmid_high=(ilint[j:*,5]-rmid[0])/(last(rmid)-rmid[0])
		xtit='r/a'	
		xr=[0,1]
	ENDELSE

	v_low_interp=interpol(simp_v_low,rmid_low,rmid_high)
	v_err_low_interp=interpol(v_err_low,rmid_low,rmid_high)
	dv=simp_v_high-v_low_interp
	dv_err=sqrt(v_err_low_interp^2+v_err_high^2)

	simp_ti_low=ilint[0:j,2]
	ti_err_low=ilint[0:j,3]
	simp_ti_high=ilint[j:*,2]
	ti_err_high=ilint[j:*,3]
	ti_low_interp=interpol(simp_ti_low,rmid_low,rmid_high)
	ti_err_low_interp=interpol(ti_err_low,rmid_low,rmid_high)
	dti=simp_ti_high-ti_low_interp
	dti_err=sqrt(ti_err_low_interp^2+ti_err_high^2)

	tmp1=where(simp_v_low NE 0 AND rmid_low LT rtmax)
	tmp2=where(simp_v_high NE 0 AND rmid_high LT rtmax)
	maxpt=max(simp_v_low[tmp1]) > max(simp_v_high[tmp2]) 
	minpt=min(simp_v_low[tmp1]) < min(simp_v_high[tmp2]) < 0
	pos=[0.075,0.35,0.475,0.95]
	IF NOT keyword_set(vr) THEN vr=[minpt,maxpt]
	plot,[0],[0],xtit=xtit,ytit='Line-Integrated v!lI!n [km/s]',chars=1.2*ls,yr=vr*1.1,xr=xr,/xsty,/ysty,pos=pos,$
		tit='                                                                                       HIREXSR Poloidal Velocity SHOT: '$
		+num2str(shot,1)+' t='+num2str(time,dp=2)
	makesym,10
	tmp=where(simp_v_low NE 0 AND rmid_low LE rtmax)
	oploterror,rmid_low[tmp],simp_v_low[tmp],fltarr(n(tmp)+1),v_err_low[tmp],psym=-8,color=200,errcolor=200
	makesym,9
	tmp=where(simp_v_high NE 0 AND rmid_high LE rtmax)
	oploterror,rmid_high[tmp],simp_v_high[tmp],fltarr(n(tmp)+1),v_err_high[tmp],psym=-8	

	tmp1=where(rmid_high LT rtmax)
	maxpt=max(dv[tmp1]+dv_err[tmp1]) > 0
	maxpt=maxpt < 3
	minpt=min(dv[tmp1]-dv_err[tmp1]) < 0
	minpt=minpt > 3*(-1.0)
	pos=[0.075,0.075,0.475,0.275]
	IF NOT keyword_set(dvr) THEN dvr=[minpt,maxpt]
	plot, [0],[0],xtit=xtit,ytit=n2g('Delta')+'v!lI!n [km/s]',chars=1.2*ls,yr=dvr*1.1,xr=xr,/xsty,/ysty,pos=pos,/noerase
	oploterror,rmid_high[tmp1],dv[tmp1],fltarr(n(tmp1)+1),dv_err[tmp1],color=30,psym=-8,errcolor=30
	oplot,xr,[0,0],linestyle=1.0

	pos=[0.55,0.35,0.95,0.95]
	tmp1=where(simp_ti_low NE 0 AND rmid_low LT rtmax)
	tmp2=where(simp_ti_high NE 0 AND rmid_high LT rtmax)
	maxpt=max(simp_ti_low[tmp1]) > max(simp_ti_high[tmp2]) 
	IF keyword_set(tir) THEN maxpt=tir[1] ELSE tir=[0,maxpt]
	plot,[0],[0],xtit=xit,ytit='Line-Integrated T!lI!n [keV]',chars=1.2*ls,yr=tir*1.1,xr=xr,/xsty,/ysty,pos=pos,/noerase
	makesym,10
	tmp=where(simp_ti_low NE 0 AND rmid_low LT rtmax)
	oploterror,rmid_low[tmp],simp_ti_low[tmp],fltarr(n(tmp)+1),ti_err_low[tmp],psym=-8,color=200,errcolor=200

	IF keyword_set(rho) THEN BEGIN
		xyouts,0.15,0.39*maxpt,'Views Below Midplane',chars=1.2*ls 
		oplot, [0.10],[0.4*maxpt],psym=8,color=200
		oplot, [0.08,0.12],[1,1]*0.4*maxpt,color=200
	ENDIF ELSE BEGIN
		xyouts,0.72,0.39*maxpt,'Views Below Midplane',chars=1.2*ls
		oplot, [0.7],[0.4*maxpt],psym=8,color=200
		oplot, [0.695,0.705],[1,1]*0.4*maxpt,color=200
	ENDELSE
	makesym,9
	tmp=where(simp_ti_high NE 0 AND rmid_high LT rtmax)
	oploterror,rmid_high[tmp],simp_ti_high[tmp],fltarr(n(tmp)+1),ti_err_high[tmp],psym=-8
	IF keyword_set(rho) THEN BEGIN
		xyouts,0.15,0.29*maxpt,'Views Above Midplane',chars=1.2*ls 
		oplot, [0.10],[0.3*maxpt],psym=8
		oplot, [0.08,0.12],[1,1]*0.3*maxpt
	ENDIF ELSE BEGIN
		xyouts,0.72,0.29*maxpt,'Views Above Midplane',chars=1.2*ls
		oplot, [0.7],[0.3*maxpt],psym=8
		oplot, [0.695,0.705],[1,1]*0.3*maxpt
	ENDELSE

	tmp1=where(rmid_high LT rtmax)
	maxpt=max(dti[tmp1]+dti_err[tmp1]) > 0
	minpt=min(dti[tmp1]-dti_err[tmp1]) < 0
	IF NOT keyword_set(dtir) THEN dtir=[minpt,maxpt]
	pos=[0.55,0.075,0.95,0.275]
	plot, [0],[0],xtit=xtit,ytit=n2g('Delta')+'T!lI!n [eV]',chars=1.2*ls,yr=dtir*1.1,xr=xr,/xsty,/ysty,pos=pos,/noerase
	oploterror,rmid_high[tmp1],dti[tmp1],fltarr(n(tmp1)+1),dti_err[tmp1],color=30,psym=-8,errcolor=30
	oplot,xr,[0,0],linestyle=1.0

	makesym,10
	rmid=rmid_high


END
