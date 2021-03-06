PRO viewbin_plot_chmap,u
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw3,get_value=draw_win
		window,0,xsize=u.plot.xsize,ysize=u.plot.ysize,/pixmap
		col=u.plot.col
        ENDIF
	IF u.stat.a THEN BEGIN
		ymax=u.dat.a.chmax
		xoff=487
		xmax=487*3.0
        ENDIF ELSE BEGIN 
		ymax=u.dat.b.chmax
		xoff=0
		xmax=487
        ENDELSE
	plot,[0],[0],xr=[0,xmax],yr=[-2,ymax],/xsty,/ysty,xtit='ROW #',ytit='CHMAP',chars=1.2
	IF u.stat.a THEN oplot,*u.dat.a.chmap,color=col[0]
	IF u.stat.b THEN oplot,indgen(487)+xoff,*u.dat.b.chmap,color=col[1]
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize,u.plot.ysize,0,0,0]
	ENDIF
END

PRO viewbin_plot_tmap,u
	IF NOT u.stat.ps THEN BEGIN
		widget_control,u.id.draw4,get_value=draw_win
		window,0,xsize=u.plot.xsize,ysize=u.plot.ysize,/pixmap
		col=u.plot.col
        ENDIF
	xmax=u.dat.nfr
	IF u.stat.a THEN ymax=u.dat.a.tmax ELSE ymax=u.dat.b.tmax
        
	plot,[0],[0],xr=[0,xmax],yr=[-2,ymax],/xsty,/ysty,xtit='ROW #',ytit='TMAP',chars=1.2
	IF u.stat.a THEN oplot,*u.dat.a.tmap,color=col[0]
	IF u.stat.b THEN oplot,*u.dat.b.tmap,color=col[1],linestyle=3.0
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.xsize,u.plot.ysize,0,0,0]
	ENDIF
END

PRO viewbin_load_data,u
	hirexsr_load_binning,u.shot,chmap,tch,tmap,good,tht=u.tht
	*u.dat.a.chmap=chmap[*,0]
	*u.dat.a.tmap=tmap
	u.dat.a.chmax=max(chmap[*,0])+1
	u.dat.a.tmax=max(tmap)+1
	u.dat.nfr=n(tmap)+1
	hirexsr_load_binning,u.shot,chmap,tch,tmap,good,tht=u.tht,/h
	*u.dat.b.chmap=chmap[*,0]
	*u.dat.b.tmap=tmap
	u.dat.b.chmax=max(chmap[*,0])+1
	u.dat.b.tmax=max(tmap)+1
	viewbin_plot_chmap,u
	viewbin_plot_tmap,u
END

PRO hirexsr_view_binning_event,event

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
					heap_free,u.dat.a.chmap
					heap_free,u.dat.a.tmap
					heap_free,u.dat.b.chmap
					heap_free,u.dat.b.tmap
					heap_gc
					!except=1
				
				END
				"LOAD": BEGIN
                   			WIDGET_CONTROL, /hourglass
					viewbin_load_data,u
                                    END
				"STOP" : stop
    				"SETA" : BEGIN
					IF event.select EQ 1 THEN u.stat.a=1 ELSE u.stat.a=0
					viewbin_plot_chmap,u
					viewbin_plot_tmap,u					
                                END
   				"SETB" : BEGIN
					IF event.select EQ 1 THEN u.stat.b=1 ELSE u.stat.b=0
					viewbin_plot_chmap,u
					viewbin_plot_tmap,u					
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
					u.tht=int(tht[0])
				END

			ELSE: 
			ENDCASE
		END
		ELSE:
        ENDCASE
	IF button NE 'QUIT' THEN widget_control,event.top,set_uvalue=u  	
END

PRO hirexsr_view_binning,shot=shot,tht=tht
	loadct,12,/silent
	base=widget_base(title='HIREXSR Binning Viewer',/row,tlb_size_events=1)
	C=widget_base(base,/column)
	dum = widget_label(C,value='TMAP/CHMAP Display')
	C1=widget_base(C,frame=5)
	draw3=widget_draw(C1,xsize=700,ysize=350)
	space=widget_base(C,/row,ysize=2)	
	C2=widget_base(C,frame=5)
	draw4=widget_draw(C2,xsize=700,ysize=350)
	space=widget_base(C2,/row,ysize=2)

	A2=widget_base(C,/column,/frame)
	space=widget_base(A2,/row,ysize=2)
	A2p1=widget_base(A2,/row)
	dum = widget_label(A2p1,value='PLOT BINNING FOR: ')
	A2p1x=widget_base(A2p1,/row,/nonexcl)
	setA=widget_button(A2p1x,value='BRANCH A')
	setB=widget_button(A2p1x,value='BRANCH B')
	space=widget_base(A2,/row,ysize=2)

	A3p1=widget_base(A2p1,/row)
	dum = widget_label(A3p1,value='SHOT:')
	shotid = widget_text(A3p1,xsize=10,ysize=1,/edit)
	dum = widget_label(A3p1,value='THT')
	thtid=widget_text(A3p1,xsize=2,ysize=1,/edit)
	dum = widget_label(A3p1,value=' ')
	load= widget_button(A3p1,value='LOAD')
	dum = widget_label(A3p1,value=' ')
	quit= widget_button(A3p1,value='QUIT')
	dum = widget_label(A3p1,value=' ')
	stop= widget_button(A3p1,value='STOP')

	IF NOT keyword_set(shot) THEN shot=1120224012
	IF NOT keyword_set(tht) THEN tht=0

	id={base:base,shotid:shotid,thtid:thtid,quit:quit,load:load,stop:stop,$
	    draw3:draw3,draw4:draw4,setA:setA,setB:setB}

	plot={xsize:700,ysize:350,col:[255,50],pscol:[0,30]}
	stat={a:1,b:1,ps:0}
	a={chmap:ptr_new([0],/allocate_heap),tmap:ptr_new([0],/allocate_heap),chmax:0,tmax:0}
	b={chmap:ptr_new([0],/allocate_heap),tmap:ptr_new([0],/allocate_heap),chmax:0,tmax:0}
	dat={a:a,b:b,nfr:0}
	u={id:id,shot:shot,tht:tht,stat:stat,dat:dat,plot:plot}
	widget_control,base,set_uvalue=u
	widget_control,u.id.shotid,set_value=num2str(u.shot,1)
	widget_control,u.id.thtid,set_value=num2str(u.tht,1)
	widget_control,u.id.setA,set_button=u.stat.a
	widget_control,u.id.setB,set_button=u.stat.b

	!except=0			
	widget_control,base,/realize
	viewbin_load_data,u
	xmanager,'hirexsr_view_binning',base
END
