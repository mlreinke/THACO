; miniprocedure for killing individual pixels

; plot the zoomed view
PRO wbin_zoom_plot, v, debug = debug
u = *v.ptru
imzoom = *v.ptrimzoom

; get the slider positions
WIDGET_CONTROL, u.id.ch_slider, get_value =ch
ch = ch-1
WIDGET_CONTROL, u.id.t_slider, get_value = slider
index=ipt(u.dat.t,slider/1.0e3)
u.index=index

WIDGET_CONTROL, v.drawid, get_value = draw_win
wset,draw_win

; directly from plot_image in w_hirexsr_he_moments

IF NOT u.stat.image THEN ff=0.0 ELSE ff=1.0
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

if not keyword_set(debug) then begin
    widget_control,u.id.gain_slider,get_value=gain
    gain*=ff
endif else begin
    gain = 1
endelse
	x=size(image)
	nx=x[1]
	ny=x[2]
	zoom=0.7;*1.2/0.7
	i_pt=indgen(nx*zoom)/zoom
	j_pt=indgen(ny*zoom)/zoom
	new_pic=interpolate(image,i_pt,j_pt,/grid)
;	loadct,39,/silent
	norm=max(new_pic)

	bnds+=offset
	r0=reform(bnds[0,*])
	r1=reform(bnds[1,*])
	;ch=u.ch-1 deleted to allow changing slider position on main widget

zoom=v.imxorig*imzoom/195.; 600*1.2/195
height=int(250*imzoom/zoom) ; scaling to 1.2*250
ex=int((height-(r1[ch]-r0[ch]))/2.0)
cent=int((r0[ch]+r1[ch])/2.0)
;IF r0[ch]-ex LT 0 THEN lowbnd=0 ELSE lowbnd=r0[ch]-ex
;IF r1[ch]+ex GT ny THEN highbnd=ny-1 ELSE highbnd=r1[ch]+ex
IF r0[ch]-ex LT 0 THEN BEGIN
    lowbnd=0 
    zoomx = zoom
    zoomy = zoom*(2.*ex+r1[ch]-r0[ch])/((2.*ex+r1[ch]-r0[ch])+2.*(r0[ch]-ex)) ;expand the y channel until it fills the screen
    height = int(250*imzoom/zoomy)	
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
         height = int(250*imzoom/zoomy)
         ex=int((height-(r1[ch]-r0[ch]))/2.0)
         lowbnd = r0[ch]-ex     ; overwrite with new ex
     ENDIF ELSE BEGIN
         highbnd = r1[ch]+ex
     ENDELSE
ENDELSE
sub=image[*,lowbnd:highbnd]
x=size(sub)
nx=x[1]
ny=x[2]
;i_pt=indgen(nx*zoom)/zoom
;j_pt=indgen(ny*zoom)/zoom
i_pt=indgen(nx*zoomx)/zoomx
j_pt=indgen(ny*zoomy)/zoomy
new_pic=interpolate(sub,i_pt,j_pt,/grid)
loadct,39,/silent
tv,new_pic/norm*256.0*gain
;contour, new_pic/norm*256*gain, nlevels = 30
loadct,12,/silent
color=200
frac=(r1[ch]-r0[ch])/float(height)/2.0
plots, [0,1,1,0,0],0.5+frac*[1.0,1.0,-1.0,-1.0,1.0],color=color,/norm
;stop
end

; find and return the pixel numbers and channel number for the
; provided x and y fractions
FUNCTION wbin_get_channel, v, xnorm, ynorm

imzoom = *v.ptrimzoom
DETECTOR_WIDTH = 195.
u = *v.ptru


WIDGET_CONTROL, u.id.ch_slider, get_value =ch
ch = ch-1
WIDGET_CONTROL, u.id.t_slider, get_value = slider
index=ipt(u.dat.t,slider/1.0e3)
u.index=index

CASE u.stat.branch OF
	0 : BEGIN
		image=rotate(*u.adat.raw[u.index],4)
		bnds=*u.adat.bnds
		offset=0
        END
	1 : BEGIN
		image=*u.bdat.raw[u.index]
		image=[image*0.0,image,image*0.0]
		image=rotate(image,4)
		bnds=*u.bdat.bnds
		offset=487
        END
ENDCASE
IF NOT u.stat.image THEN ff=0.0 ELSE ff=1.0

x=size(image)
if not keyword_set(debug) then begin
    widget_control,u.id.gain_slider,get_value=gain
    gain*=ff
endif else begin
    gain = 1
endelse

;	widget_control,u.id.gain_slider,get_value=gain
;	gain*=ff
	nx=x[1]
	ny=x[2]
	zoom=0.7;*1.2/0.7
	i_pt=indgen(nx*zoom)/zoom
	j_pt=indgen(ny*zoom)/zoom
	new_pic=interpolate(image,i_pt,j_pt,/grid)
;	loadct,39,/silent
	norm=max(new_pic)

	bnds+=offset
	r0=reform(bnds[0,*])
	r1=reform(bnds[1,*])

	;ch=u.ch-1 deleted to allow changing slider position on main widget
        zoom=v.imxorig*imzoom/195. ; 600*1.2/195
        height=int(250*imzoom/zoom) ; scaling to 1.2*250
        ex=int((height-(r1[ch]-r0[ch]))/2.0)
        cent=int((r0[ch]+r1[ch])/2.0)
        ;IF r0[ch]-ex LT 0 THEN lowbnd=0 ELSE lowbnd=r0[ch]-ex
        ;IF r1[ch]+ex GT ny THEN highbnd=ny-1 ELSE highbnd=r1[ch]+ex
        IF r0[ch]-ex LT 0 THEN BEGIN
            lowbnd=0 
            zoomx = zoom
            zoomy = zoom*(2.*ex+r1[ch]-r0[ch])/((2.*ex+r1[ch]-r0[ch])+2.*(r0[ch]-ex)) ;expand the y channel until it fills the screen
            height = int(250*imzoom/zoomy)	
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
                height = int(250*imzoom/zoomy)
                ex=int((height-(r1[ch]-r0[ch]))/2.0)
                lowbnd = r0[ch]-ex ; overwrite with new ex
            ENDIF ELSE BEGIN
                highbnd = r1[ch]+ex
            ENDELSE
        ENDELSE

        sub=image[*,lowbnd:highbnd]
        x=size(sub)
        nx=x[1]
        ny=x[2]
        ;i_pt=indgen(nx*zoom)/zoom
        ;j_pt=indgen(ny*zoom)/zoom
        i_pt=indgen(nx*zoomx)/zoomx
        j_pt=indgen(ny*zoomy)/zoomy    
        new_pic=interpolate(sub,i_pt,j_pt,/grid)

        channely = FIX(lowbnd+(highbnd-lowbnd)*ynorm)
        channelx = FIX(DETECTOR_WIDTH*xnorm) ; 195 is detector width
        chmap =  *v.ptrinternalchmap
	;print, channely
        chval = chmap(channely-offset, channelx)
        pixval = image(channelx, channely-offset)

        returnstr = {chval:chval, channelx:channelx, channely:channely-offset, pixval:pixval}
        return, returnstr
end

; set the selected pixel to ignored
FUNCTION wbin_kill_pixel, v
u = *v.ptru
chmap = *v.ptrinternalchmap

limitx = n_elements(chmap(0, *))
limity = n_elements(chmap(*, 0))

WIDGET_CONTROL, v.infoxpixtext, get_value = channelx
WIDGET_CONTROL, v.infoypixtext, get_value = channely

; sanity checks
if channelx LT 0 OR channelx GT limitx then begin
    print, "NOT VALID PIXEL"
    return, -1
endif
if channely LT 0 OR channely GT limity then begin
    print, "NOT VALID PIXEL"
    return, -1
endif
chmapval = (*v.ptrinternalchmap)[channely, channelx]

if chmapval LT 0 then begin
    ; this pixel is already being ignored. 
    return, 0
endif

; now store the value for undoing
if n_elements(*v.ptrundoarr) LT 3 then begin
    undoarr = [channelx, channely, chmapval]
    *v.ptrundoarr = undoarr
endif else begin
    temparr = [channelx, channely, chmapval]
    undoarr = *v.ptrundoarr
    undoarr = [temparr, undoarr]
    *v.ptrundoarr = undoarr
endelse

; now kill the pixel
(*v.ptrinternalchmap)[channely, channelx] = -2.

return, 1; success

end

; undo last kill pixel event
FUNCTION wbin_undo_killpixel, v
undoarr = *v.ptrundoarr
chmap = *v.ptrinternalchmap

if n_elements(undoarr) GE 3 then begin
    channelx = undoarr(0)
    channely = undoarr(1)
    chval = undoarr(2)
    
    ;reset the channel map
    (*v.ptrinternalchmap)[channely, channelx] = chval

    if n_elements(undoarr) EQ 3 then begin
        undoarr = fltarr(1)-1; make a dummy array
        *v.ptrundoarr = undoarr
        return, 1
    endif else begin
        nelem = n_elements(undoarr)
        temp_undoarr = undoarr(3:nelem-1)
        undoarr = temp_undoarr
        *v.ptrundoarr = undoarr
        return, 2
    endelse

endif else begin
    return, -1 ; no undo array exists
endelse


end

; overplot killed pixels with black 
PRO wbin_oplot_rempixels, v

DETECTOR_WIDTH = 195.
;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;; standard loading stuff
u = *v.ptru
imzoom = *v.ptrimzoom

; get the slider positions
WIDGET_CONTROL, u.id.ch_slider, get_value =ch
ch = ch-1

WIDGET_CONTROL, u.id.t_slider, get_value = slider
index=ipt(u.dat.t,slider/1.0e3)
u.index=index

WIDGET_CONTROL, v.drawid, get_value = draw_win
wset,draw_win

; directly from plot_image in w_hirexsr_he_moments

IF NOT u.stat.image THEN ff=0.0 ELSE ff=1.0
CASE u.stat.branch OF
	0 : BEGIN
		image=rotate(*u.adat.raw[u.index],4)
		bnds=*u.adat.bnds
		offset=0
        END
	1 : BEGIN
		image=*u.bdat.raw[u.index]
		image=[image*0.0,image,image*0.0]
		image=rotate(image,4)
		bnds=*u.bdat.bnds
		offset=487
        END
ENDCASE
x=size(image)
if not keyword_set(debug) then begin
    widget_control,u.id.gain_slider,get_value=gain
    gain*=ff
endif else begin
    gain = 1
endelse

;	widget_control,u.id.gain_slider,get_value=gain
;	gain*=ff
	nx=x[1]
	ny=x[2]
	zoom=0.7;*1.2/0.7
	i_pt=indgen(nx*zoom)/zoom
	j_pt=indgen(ny*zoom)/zoom
	new_pic=interpolate(image,i_pt,j_pt,/grid)
;	loadct,39,/silent
	norm=max(new_pic)

	bnds+=offset
	r0=reform(bnds[0,*])
	r1=reform(bnds[1,*])
  
	;ch=u.ch-1 deleted to allow changing slider position on main widget

zoom=v.imxorig*imzoom/DETECTOR_WIDTH; 600*1.2/195
height=int(250*imzoom/zoom) ; scaling to 1.2*250
ex=int((height-(r1[ch]-r0[ch]))/2.0)
cent=int((r0[ch]+r1[ch])/2.0)
;IF r0[ch]-ex LT 0 THEN lowbnd=0 ELSE lowbnd=r0[ch]-ex
;IF r1[ch]+ex GT ny THEN highbnd=ny-1 ELSE highbnd=r1[ch]+ex
IF r0[ch]-ex LT 0 THEN BEGIN
    lowbnd=0 
    zoomx = zoom
    zoomy = zoom*(2.*ex+r1[ch]-r0[ch])/((2.*ex+r1[ch]-r0[ch])+2.*(r0[ch]-ex)) ;expand the y channel until it fills the screen
    height = int(250*imzoom/zoomy)	
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
        height = int(250*imzoom/zoomy)
        ex=int((height-(r1[ch]-r0[ch]))/2.0)
        lowbnd = r0[ch]-ex      ; overwrite with new ex
    ENDIF ELSE BEGIN
        highbnd = r1[ch]+ex
    ENDELSE
ENDELSE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; end standard loading stuff
; now zero out the bad pixels


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; YP remove all bad pixels
;undoarr = *v.ptrundoarr
;index = 0
;numundoarr = n_elements(undoarr)
;if numundoarr GT 1 then begin
;    while index*3 LT numundoarr do begin
;        channelx = undoarr(index*3)
;        channely = undoarr(index*3+1)
;        image(channelx, channely) = 0
;        index++
;    endwhile
;endif
internalchmap = *(v.ptrinternalchmap)
killedpixels= where(internalchmap EQ -2)
if n_elements(killedpixels) GT 0 then begin
    ;print, killedpixels
    if killedpixels(0) NE -1 then begin
        arr_index = array_indices(internalchmap, killedpixels)
        ;if n_elements(killedpixels) GT 4 then stop
        image(arr_index(1, *), arr_index(0, *)+offset) =0 ; zero out bad pixels
    endif 
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

sub=image[*,lowbnd:highbnd]
x=size(sub)
nx=x[1]
ny=x[2]
;i_pt=indgen(nx*zoom)/zoom
;j_pt=indgen(ny*zoom)/zoom
i_pt=indgen(nx*zoomx)/zoomx
j_pt=indgen(ny*zoomy)/zoomy

new_pic=interpolate(sub,i_pt,j_pt,/grid)
loadct,39,/silent
tv,new_pic/norm*256.0*gain
;contour, new_pic/norm*256*gain, nlevels = 30
loadct,12,/silent
color=200
frac=(r1[ch]-r0[ch])/float(height)/2.0
plots, [0,1,1,0,0],0.5+frac*[1.0,1.0,-1.0,-1.0,1.0],color=color,/norm

end


; widget end event
PRO wbin_zoom_end, widget
WIDGET_CONTROL, widget, GET_UVALUE = v, /no_copy

if size(v, /type) EQ 8 then begin
    PTR_FREE, v.ptru
    PTR_FREE, v.ptrinternalchmap
    PTR_FREE, v.ptrundoarr
    PTR_FREE, v.ptrimzoom
    widget_control, v.base, /destroy
endif

end

PRO wbin_zoom_kill_event, ev

widget_control, ev.top, get_uvalue = v
tag = tag_names(ev, /st)
case tag OF
    "WIDGET_BUTTON": BEGIN
        WIDGET_CONTROL, ev.id, get_value = button
        case button OF
            'QUIT': BEGIN
                wbin_zoom_end, ev.top
            END

            'STOP': begin
                stop ; debugging only
            end

            
            'AUTO': BEGIN
                WIDGET_CONTROL, v.killbutton, sensitive = 0
                ; turn on motion events on draw screen
                WIDGET_CONTROL, v.drawid, DRAW_MOTION_EVENTS = 1
            END


           'MANUAL': BEGIN
               WIDGET_CONTROL, v.killbutton, sensitive = 1
               ; turn off motion events on draw screen
               WIDGET_CONTROL, v.drawid, DRAW_MOTION_EVENTS = 0
           END

           'REMOVE': BEGIN
               stat = wbin_kill_pixel(v)
               if stat EQ 1 then begin
                   WIDGET_CONTROL, v.infottext, set_value = STRCOMPRESS(STRING(-2), /remove_all)
                   wbin_oplot_rempixels, v
               endif
           END

           'LOAD MAP': BEGIN
               internalchmap = *v.ptrinternalchmap
               u = *v.ptru
	 	CASE u.stat.branch OF
			0 : *u.adat.chk = internalchmap
			1 : *u.bdat.chk = internalchmap
		ENDCASE
           END

           'RAW': BEGIN
               wbin_zoom_plot, v
           END
           
           'IMAGE': BEGIN
               wbin_oplot_rempixels, v
           END


           'UNDO':BEGIN
               stat = wbin_undo_killpixel(v)
               wbin_oplot_rempixels, v
           END


           ELSE: BEGIN
               print, 'HIREX_SR_ZOOM_KILL: UNHANDLED EXCEPTION 2'
           END

        endcase


    END

    "WIDGET_TEXT_CH": BEGIN
        WIDGET_CONTROL, v.zoomtext, get_value = ztext
        imzoom = FLOAT(ztext)
        *v.ptrimzoom = imzoom
        print, *v.ptrimzoom
        wbin_oplot_rempixels, v
    END


    "WIDGET_TRACKING" : BEGIN
        if ev.enter EQ 1 then begin
            wbin_oplot_rempixels, v
        endif
    END


    "WIDGET_DRAW": BEGIN
         ;click_loc = convert_coord(ev.x, ev.y)
        xval = ev.x
        yval = ev.y
        imzoom = *v.ptrimzoom
        xnorm = ev.x/(imzoom*v.imxorig)
        ynorm = ev.y/(imzoom*v.imyorig)

        if WIDGET_INFO(v.buttonauto, /button_set) then begin
            u = *v.ptru
            ;WIDGET_CONTROL, v.infoxtext, set_value = STRCOMPRESS(STRING(xnorm), /remove_all)
            ;WIDGET_CONTROL, v.infoytext, set_value = STRCOMPRESS(STRING(ynorm), /remove_all)
            ; now calculate which 
            WIDGET_CONTROL, u.id.ch_slider, get_value =ch
            str_ch = wbin_get_channel(v, xnorm, ynorm)
            chmapval = str_ch.chval
            chxpix = str_ch.channelx
            chypix = str_ch.channely
            pixval  = str_ch.pixval
            WIDGET_CONTROL, v.infoxpixtext, set_value =STRCOMPRESS(STRING(chxpix),/remove_all)
            WIDGET_CONTROL, v.infoypixtext, set_value =STRCOMPRESS(STRING(chypix),/remove_all)
            WIDGET_CONTROL, v.infozpixtext, set_value =STRCOMPRESS(STRING(pixval),/remove_all)
            WIDGET_CONTROL, v.infottext, set_value =STRCOMPRESS(STRING(chmapval),/remove_all)

        endif       

        
        if (ev.press eq 1) then begin
            u = *v.ptru
            ;WIDGET_CONTROL, v.infoxtext, set_value = STRCOMPRESS(STRING(xnorm), /remove_all)
            ;WIDGET_CONTROL, v.infoytext, set_value = STRCOMPRESS(STRING(ynorm), /remove_all)
            ; now calculate which 
            WIDGET_CONTROL, u.id.ch_slider, get_value =ch
            str_ch = wbin_get_channel(v, xnorm, ynorm)
            chmapval = str_ch.chval
            chxpix = str_ch.channelx
            chypix = str_ch.channely
            pixval  = str_ch.pixval
            WIDGET_CONTROL, v.infoxpixtext, set_value =STRCOMPRESS(STRING(chxpix),/remove_all)
            WIDGET_CONTROL, v.infoypixtext, set_value =STRCOMPRESS(STRING(chypix),/remove_all)
            WIDGET_CONTROL, v.infozpixtext, set_value =STRCOMPRESS(STRING(pixval),/remove_all)
            WIDGET_CONTROL, v.infottext, set_value =STRCOMPRESS(STRING(chmapval),/remove_all)

            if WIDGET_INFO(v.buttonauto, /button_set) then begin
                stat = wbin_kill_pixel(v)
                if stat EQ 1 then begin
                    WIDGET_CONTROL, v.infottext, set_value = STRCOMPRESS(STRING(-2), /remove_all)
                    wbin_oplot_rempixels, v
                endif
            endif
        endif



	   
    END

    ELSE: BEGIN
        print, 'HIREX_SR_ZOOM_KILL: UNHANDLED EXCEPTION 1 ' + STRING(tag)
        ; do nothing else other than the warning
    END
endcase


end


; u is the array of information from the master w_hirexsr_he_moments
; debug activates debug mode, which should only be used during
; development
; written by YP 08/2010
PRO wbin_zoom_kill, ptru, debug=debug
	u= *ptru
	CASE u.stat.branch OF
		0 : internalchmap = *u.adat.chk
		1 : internalchmap = *u.bdat.chk
	ENDCASE
	ptrinternalchmap = PTR_NEW(internalchmap)
	undoarr = findgen(1) -1; init undoarr
	ptrundoarr = PTR_NEW(undoarr)

	user=logname()
	imzoom = 2.25 ; zoomed image level
	plotzoom = 2.25 ; size of plot window
	imxorig = 600. ; zoomed values from parent widget
	imyorig = 250. ; zoomed values from parent widget

	ptrimzoom = PTR_NEW(imzoom)

	base = WIDGET_BASE(title = 'HiReX Sr Pixel Removing', /row,tlb_size_events=1, group_leader = u.id.base, /TRACKING_EVENTS, kill_notify = 'wbin_zoom_end')

	textbase = WIDGET_BASE(base, /column, xsize = 170)
	drawbase = WIDGET_BASE(base, /column, xsize = imzoom*imxorig)

	;text base
	clickbase  = WIDGET_BASE(textbase, /row, /exclusive, frame = 1)
	buttonauto = WIDGET_BUTTON(clickbase, VALUE = 'AUTO')
	buttonman  = WIDGET_BUTTON(clickbase, VALUE = 'MANUAL')

	; information on pixel selected
	infobase = WIDGET_BASE(textbase, /column, frame = 1)

	infoxpixbase  = WIDGET_BASE(infobase, /row)
	infoxpixlabel = WIDGET_LABEL(infoxpixbase, VALUE = 'X pix: ', /align_left, xsize = 40)
	infoxpixtext  = WIDGET_TEXT(infoxpixbase, xsize = 10)
	infoypixbase  = WIDGET_BASE(infobase, /row)
	infoypixlabel = WIDGET_LABEL(infoypixbase, VALUE = 'Y pix: ', /align_left, xsize = 40)
	infoypixtext  = WIDGET_TEXT(infoypixbase, xsize = 10)
	infozpixbase  = WIDGET_BASE(infobase, /row)
	infozpixlabel = WIDGET_LABEL(infozpixbase, VALUE = 'V pix: ', /align_left, xsize = 40)
	infozpixtext  = WIDGET_TEXT(infozpixbase, xsize = 10)
	infotbase  = WIDGET_BASE(infobase, /row)
	infotlabel = WIDGET_LABEL(infotbase, VALUE = 'Chan: ', /align_left, xsize = 40)
	infottext  = WIDGET_TEXT(infotbase, xsize = 10)

	; various buttons
	buttonbase = WIDGET_BASE(textbase, /column)
	topbuttonbase = WIDGET_BASE(buttonbase, /column)
	bottombuttonbase = WIDGET_BASE(buttonbase, /row,  frame = 1)
	topubuttonbase = WIDGET_BASE(topbuttonbase, /row)
	topbbuttonbase = WIDGET_BASE(topbuttonbase, /row)
	killbutton = WIDGET_BUTTON(topubuttonbase, VALUE = 'REMOVE')
	undobutton = WIDGET_BUTTON(topubuttonbase, VALUE = 'UNDO')
	rawbutton = WIDGET_BUTTON(topbbuttonbase, VALUE = 'RAW') ; raw image with killed pixels shown
	reloadbutton = WIDGET_BUTTON(topbbuttonbase, VALUE = 'IMAGE') ; killed pixels not shown


	if keyword_set(debug) then stopbutton = WIDGET_BUTTON(bottombuttonbase, value = 'STOP')
	respondbutton= WIDGET_BUTTON(bottombuttonbase, VALUE = 'LOAD MAP') ; load into original channel map
	quitbutton = WIDGET_BUTTON(bottombuttonbase, VALUE = 'QUIT') ; quit


	;draw base
	drawid = WIDGET_DRAW(drawbase, xsize = 600*plotzoom, ysize = 250*plotzoom, /button_events, /motion_events);, /scroll, x_scroll_size = 600*plotzoom, y_scroll_size = 250*plotzoom)

	v = {base:base, drawid:drawid, buttonauto:buttonauto, buttonman:buttonman,  infozpixtext:infozpixtext, $
	     infottext:infottext, infoxpixtext:infoxpixtext, infoypixtext:infoypixtext, killbutton:killbutton, undobutton:undobutton, $
	     quitbutton:quitbutton, ptrimzoom:ptrimzoom, plotzoom:plotzoom, imxorig:imxorig, imyorig:imyorig, reloadbutton:reloadbutton, $
	    respondbutton:respondbutton, ptrinternalchmap:ptrinternalchmap,  ptrundoarr:ptrundoarr, ptru:ptru}
	WIDGET_CONTROL, base, set_uvalue=v
	WIDGET_CONTROL, buttonauto, set_button=1
	WIDGET_CONTROL, v.killbutton, sensitive = 0

	WIDGET_CONTROL, base, /realize

	if not keyword_set(debug) then begin
	   wbin_zoom_plot, v
	endif

	xmanager, 'wbin_zoom_kill', base

END
