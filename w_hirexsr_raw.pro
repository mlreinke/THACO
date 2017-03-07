; This is a utility for looking at HiReX Sr raw data
; It is a copy of the McPherson utility written by Jim Terry for the
; Hirex Sr
;
; I (YP) deny all culpability for the common block usage and goto commands

;**************************************************************************
; a small utility from Matt Reinke
;**************************************************************************
PRO xwplot
	set_plot,'x'
	device,true_color=24,retain=2,decomposed=0  
END


;**************************************************************************
Pro GET_TIMES
;**************************************************************************
; this procedure sets up a widget that gets the times for curve fit.

common hirex_base, w1, time1_label, time2_label, w2, w3, ok_button, reply1_text, reply2_text 
common hirex_times, time1, time2

base = widget_base(/column, title='User Dialogue', xoffset=200, yoffset=200)
  w1 = widget_base(base, /row)
    time1_label = widget_label(w1, value='Enter start and stop time for the average[times seperated by comma]:')
    reply1_text = widget_text(w1, value='', xsize=15, ysize=1, /editable)
;  w3 = widget_base(base, /row)
;    time2_label = widget_label(w3, value='select the stop time for the average:')
;    reply2_text = widget_text(w3, value='', xsize=5, ysize=1, /editable)
  w2 = widget_base(base, /column)
    ok_button = widget_button(w2, value='OK')

widget_control, base, /realize
xmanager, 'get_times', base, event_handler='get_times_event', /modal

END

;*********************************************************************
PRO get_times_event, ev
;*********************************************************************

common hirex_base, w1, time1_label, time2_label, w2, w3, ok_button, reply1_text, reply2_text
common hirex_times, time1, time2

type = tag_names(ev, /structure)
widget_control, ev.id, get_value=button_name
widget_control, reply1_text, get_value=time_text
;widget_control, reply2_text, get_value=time2
comma=strpos(time_text,',')
time1=float(strmid(time_text,0,comma(0)))
text_len=strlen(time_text)
time2=float(strmid(time_text,comma(0)+1,text_len(0)-comma(0)))

widget_control, /destroy, ev.top

End

;*************************************************************************
FUNCTION yes, message, ALERT=alert
;*************************************************************************
; from chr_view

common hirex_return_value, return_value

base = widget_base(/column, title='User Alert', xoffset=200, yoffset=200)
  w1 = widget_base(base, /row)
    message_text = widget_text(w1, value=message, xsize=strlen(message), ysize=1)
  w2 = widget_base(base, /row, space=25)
    if keyword_set(alert) then ok_button = widget_button(w2, value='OK') $
    else begin
    yes_button = widget_button(w2, value='Yes')
    no_button = widget_button(w2, value='No')
    endelse

widget_control, base, /realize
xmanager, 'alertbox', base, event_handler='yes_event', /modal

return, return_value

END

;******************************************************************************
PRO yes_event, ev
;******************************************************************************
; from chr_view

common hirex_return_value, return_value

type = tag_names(ev, /structure)
widget_control, ev.id, get_value=button_name

;if (type eq 'WIDGET_BUTTON') then begin
  widget_control, /destroy, ev.top
  if (button_name eq 'Yes') then return_value=-1 else return_value=0
;endif

END

;******************************************************************************
FUNCTION ask, message
;******************************************************************************
; This is a widget that emulates a Mac-style dialog box. It asks the user for
; input and returns the reply as a string.
; from chr_view
common hirex_reply, reply
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type
common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type, data_type

base = widget_base(/column, title='User Dialogue', xoffset=200, yoffset=200)
  w1 = widget_base(base, /row)
    message_text = widget_text(w1, value=message, xsize=strlen(message), ysize=1)
    reply_text = widget_text(w1, value='', xsize=32, ysize=1, /editable)
  w2 = widget_base(base, /column)
    ok_button = widget_button(w2, value='OK', uvalue=reply_text)




widget_control, base, /realize
xmanager, 'ask', base, event_handler='ask_event', /modal

return, reply

END

;******************************************************************************
PRO ask_event, ev
;******************************************************************************
; event handler from chr_view

common hirex_reply, reply
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4,
;ppflag, xra
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type
common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type, data_type
type = tag_names(ev, /structure)
widget_control, ev.id, get_value=button_name

if (type eq 'WIDGET_BUTTON') then begin
  widget_control, ev.id, get_uvalue=reply_text
  widget_control, reply_text, get_value=reply
  widget_control, /destroy, ev.top
endif

END

;**************************************************************************
pro gfunct,x,a,f,pder
;**************************************************************************

common hirex_widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, $
        w3c, w3d, w4, draw, shot_label, shot_text, shot_button, chord_slider, $
        group_slider, frame_slider, groupplot_button, frameplot_button, $
        zoom1_button, zoom2_button, timehist_button, movie_button, from_label, from_text, $
        from_label2, to_label, to_text, to_label2, xclick_label, $
        xclick_text, xclick_label2, yclick_label, yclick_text, $
        yclick_label2, print_button, done_button, curvefit_button, $
        ylabel2, di_box, w1w, linewavelength_button, current_button
common hirex_l1, wave0, width, grating
common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
        w_stop, dmaxval, maxval1, maxval, minval, tstring, d_maxval, pos, $
	test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
	 d_fit_time, center_fit, fwhm,print_queue_num, hirex_chord_index, hirex_sz

big=n_elements(x)
pks=(n_elements(a)-1)/2
f=fltarr(big)
pder=fltarr(big,n_elements(a))
f(*)=0.

if(grating eq 600) then fact1=.15 & fact2=3.2
if(grating eq 300) then fact1=.15 & fact2=2.7

for i=0,pks-1 do begin
ad=2*i
n_width=n_elements(a)-1
width=a(n_width)

f=f+a(0+ad)*(exp(-(x-a(1+ad))^2/width^2)+fact1*exp(-(x-a(1+ad))^2/(fact2*width)^2))

pder(*,0+ad)=exp(-(x-a(1+ad))^2/width^2)+fact1*exp(-(x-a(1+ad))^2/(fact2*width)^2)
pder(*,1+ad)=a(0+ad)*(exp(-(x-a(1+ad))^2/width^2)*2*(x-a(1+ad))/width^2+$
fact1*exp(-(x-a(1+ad))^2/(fact2*width)^2)*2*(x-a(1+ad))/(fact2*width)^2)
endfor
pder(*,n_width)=a(0+ad)*(exp(-(x-a(1+ad))^2/width^2)*2*(x-a(1+ad))^2/width^3+$
fact1*exp(-(x-a(1+ad))^2/(fact2*width)^2)*2*(x-a(1+ad))^2/(fact2^2*width^3))

return

end

;**********************************************************************
PRO hirex_get_shot, type = type
;**********************************************************************
; this procedure gets the data for a given shot.

common hirex_widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, $
        w3c, w3d, w4, draw, shot_label, shot_text, shot_button, chord_slider, $
        group_slider, frame_slider, groupplot_button, frameplot_button, $
        zoom1_button, zoom2_button, timehist_button, movie_button, from_label, from_text, $
        from_label2, to_label, to_text, to_label2, xclick_label, $
        xclick_text, xclick_label2, yclick_label, yclick_text, $
        yclick_label2, print_button, done_button, curvefit_button, $
	ylabel2, di_box, w1w, linewavelength_button, current_button
common hirex_l1,wave0,width,grating
;common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
;        w_stop, dmaxval, maxval1, maxval, minval, tstring, d_maxval, pos, $
;	test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
;	 d_fit_time, center_fit,fwhm,print_queue_num, hirex_chord_index
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra
common hirex_print,x,y,tit,xtit,ytit,yra,ox1,oy1,ox2,oy2,ox3,oy3, line
common hirex_zmcp, z_mcp ;for the known wavelengths
common hirex_line_data, w, ident, line_plot, zoomed
common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
        w_stop, dmaxval, maxval1, maxval, minval, tstring, d_maxval, pos, $
	test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
	 d_fit_time, center_fit, print_queue_num, hirex_chord_index, hirex_sz
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type
common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type, data_type
;print, hirex_chord_index
pflag1 = 0
pflag2 = 0
pflag3 = 0
pflag4 = 0
ppflag = 0
; get a new shot

widget_control, shot_text, get_value = shot_str
shot_num = long(shot_str(0))

print, 'accessing database'
widget_control, /hourglass
mdsopen,'spectroscopy',shot_num
;d=mdsvalue('\SPECTROSCOPY::TOP.VUV.ANALYSIS:MCP_BRI_MOD',/quiet,status=mod_status)
;if(mod_status and d(0) ne -1) then begin
;    good_time=mdsvalue('dim_of(\SPECTROSCOPY::TOP.VUV.ANALYSIS:MCP_BRI_MOD)',/quiet,status=status)
;    waves=mdsvalue('dim_of(\SPECTROSCOPY::TOP.VUV.ANALYSIS:MCP_BRI_MOD,1)',/quiet,status=status)
;endif else begin
;    mdsset_def,'\spectroscopy::top.vuv'
;    d= MDSVALUE('mcp_bright')
;    good_time=mdsvalue('dim_of(mcp_bright)')
;    waves=mdsvalue('dim_of(mcp_bright,1)')
;endelse
if type EQ 0 then begin ; Helium
;d = mdsvalue('\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:HE_LIKE:INTEN',
;/quiet, status =mod_status)
ave = mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.SPEC:SPECBR', /quiet, status = mod_status)
d = transpose(ave, [2, 0, 1])
if (mod_status NE -1) then begin
    ; initialize to middle chord
    hirex_chord_index = 23
    good_time = mdsvalue('dim_of(\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.SPEC:SPECBR,1)', /quiet, status = status)
    waves = transpose(mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.SPEC:LAM', /quiet, status = status), [2, 0, 1])
    hirex_sz = mdsvalue('dim_of(\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.SPEC:SPECBR,0)', /quiet, status = status)   
    ; now remove extraneous data points
    valid_times = where(good_time NE -1)
    d = d(*, *, valid_times)
    waves = waves(*, *, valid_times)
    hirexsr_load_binning,shot_num, chmap,  tch, tmap, good, chmax
    maxchan = max(chmap)
    d = d(*, 0:maxchan, *)
    waves = waves(*, 0:maxchan, *)
    hirex_sz = hirex_sz(0:maxchan)
    widget_control, chord_slider, set_slider_max = n_elements(hirex_sz)
    widget_control, chord_slider, set_value = FIX(n_elements(hirex_sz)/2)

    hirex_chord_index = FIX(n_elements(hirex_sz)/2)-1   
    
    z_mcp = hirex_sz(hirex_chord_index)
    pos = ' Line-of-sight at SZ = ' + strtrim(string(z_mcp), 2)
    

endif
endif else begin ;hydrogen
;d = mdsvalue('\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:H_LIKE:INTEN',
;/quiet, status =mod_status)
ave = mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.SPEC:SPECBR', /quiet, status = mod_status)
d = transpose(ave, [2, 0, 1])
if (mod_status NE -1) then begin
    ;hirex_chord_index = 7
    good_time = mdsvalue('dim_of(\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.SPEC:SPECBR,1)', /quiet, status = status)
    waves = transpose(mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.SPEC:LAM', /quiet, status = status), [2, 0, 1])
    hirex_sz = mdsvalue('dim_of(\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.SPEC:SPECBR,0)', /quiet, status = status)    
; now remove extraneous data points
    valid_times = where(good_time NE -1)
    d = d(*, *, valid_times)
    waves = waves(*, *, valid_times)

    hirexsr_load_binning,shot_num, chmap,  tch, tmap, good, chmax
    maxchan = max(chmap)
    d = d(*, 0:maxchan, *)
    waves = waves(*, 0:maxchan, *)
    hirex_sz = hirex_sz(0:maxchan)

    widget_control, chord_slider, set_slider_max = n_elements(hirex_sz)
    ;widget_control, chord_slider, set_value = FIX(n_elements(hirex_sz)/2)
    ;hirex_chord_index = FIX(n_elements(hirex_sz)/2)-1   

    ;h-like should be set off axis
    hirex_chord_index = 2;
    widget_control, chord_slider, set_value=2
    
    z_mcp = hirex_sz(hirex_chord_index)
    pos = ' Line-of-sight at SZ = ' + strtrim(string(z_mcp), 2)

endif

endelse

;if(shot_num eq 940601022) then begin
;        jackpos=5.02
;        goto, jackjump
;endif
;jackpos=mdsvalue('\spectroscopy::top.vuv:jack_pos',/quiet,status=status)
if(not status) then begin
  dummy = yes('Could not access the shot requested. Enter another shot number.',$
  /alert)
return
endif
;grating=mdsvalue('\top.vuv:mcp_grating',/quiet,status=status)
;if(not status) then begin
;        grating=600.
;endif

;if(shot_num gt 970501000) then begin
;       at this time I put the correctly calculated nodes in the tree
;        view_ang=mdsvalue('\SPECTROSCOPY::MCP_VIEW_ANG')
;        pivot=mdsvalue('\SPECTROSCOPY::MCP_PIVOT_PT')
;       z_mcp=pivot(1)-(pivot(0)-.67)*tan(view_ang)
;        z_mcp=pivot(1)-(pivot(0)-.44)*tan(view_ang)
;        pos=' (Z of line-of-sight at R=0.44='+strtrim(string(form='(f7.4)',z_mcp),2)+')'
;endif else begin
;        if(shot_num gt 950315000) then jackpos=jackpos+0.36
;        jackjump:
;        dis=33.926+.8441*jackpos
;        plas_ang=6.98e-2+.3047-asin((59.25-dis)/69.88)+0.02
;        pivot=[3.738,.30607]
;;       z_mcp=.30607-(3.738-.67)*tan(plas_ang)
;        z_mcp=pivot(1)-(pivot(0)-.44)*tan(plas_ang)
;        pos=' (Z of line-of-sight at R=0.44='+strtrim(string(form='(f7.4)',z_mcp),2)+')'
;endelse
mdsclose
frames=n_elements(d(0, hirex_chord_index, *))
if(frames ge 6) then begin
;        w_start=(where(d(*,5) ne 0))(0)
;        w_stop=(where(d(*,5) ne 0))(n_elements(where(d(*,5) ne 0))-1)
        w_start=(where(d(*,hirex_chord_index, frames/2) ne 0))(0)
        w_stop=(where(d(*,hirex_chord_index, frames/2) ne 0))(n_elements(where(d(*,hirex_chord_index, frames/2) ne 0))-1)
endif else begin
        w_start=(where(d(*,hirex_chord_index, frames-1) ne 0))(0)
        w_stop=(where(d(*,hirex_chord_index, frames-1) ne 0))(n_elements(where(d(*,hirex_chord_index, frames-1) ne $
0))-1)
endelse
;endelse
;d_maxval=d(w_start+30:w_stop-90,0:frames-2)
d_maxval=d(w_start:w_stop, hirex_chord_index, 0:frames-1)
maxval=max(d_maxval)
minval=min(d_maxval)

widget_control, frame_slider, set_slider_max=frames-1
;widget_control, frame_slider, set_slider_max=frames-3
;widget_control, group_slider, set_slider_max=frames-20
if(frames lt 20) then widget_control, group_slider, set_slider_max=frames-1
if(frames ge 20) then widget_control, group_slider, set_slider_max=frames-20
widget_control, yclick_text, set_value= strcompress(string(z_mcp), /remove_all)
;reset wavelengths
w=0
line_plot=0

END

;*************************************************************************
PRO plot_group
;*************************************************************************
;this procedure plots 20 frames at a time.

common hirex_widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, $
        w3c, w3d, w4, draw, shot_label, shot_text, shot_button, chord_slider, $
        group_slider, frame_slider, groupplot_button, frameplot_button, $
        zoom1_button, zoom2_button, timehist_button, movie_button, from_label, from_text, $
        from_label2, to_label, to_text, to_label2, xclick_label, $
        xclick_text, xclick_label2, yclick_label, yclick_text, $
        yclick_label2, print_button, done_button, curvefit_button, $
	ylabel2, di_box, w1w, linewavelength_button, current_button
common hirex_l1,wave0,width,grating
;common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
;        w_stop, dmaxval, maxval1, maxval,  minval, tstring, d_maxval, pos, $
;	test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
;	 d_fit_time, center_fit,fwhm,print_queue_num, hirex_chord_index
common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
        w_stop, dmaxval, maxval1, maxval, minval, tstring, d_maxval, pos, $
	test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
	 d_fit_time, center_fit, print_queue_num, hirex_chord_index, hirex_sz
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra
common hirex_print,x,y,tit,xtit,ytit,yra,ox1,oy1,ox2,oy2,ox3,oy3, line
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type
common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type, data_type
pflag1 = 0
pflag2 = 0
pflag3 = 0
pflag4 = 0
ppflag = 0
widget_control, group_slider, sensitive = 1
widget_control, frame_slider, sensitive = 1
widget_control, linewavelength_button, sensitive = 0
widget_control, zoom1_button, sensitive = 0
widget_control, zoom2_button, sensitive = 0
widget_control, timehist_button, sensitive = 0
widget_control, movie_button, sensitive = 0
;widget_control, curvefit_button, sensitive = 0
widget_control, print_button, sensitive = 0
!p.multi=[0,5,4]

widget_control, group_slider, get_value = x
if (x+19 lt frames-1) then begin
  for i=x,x+19 do begin

  tstring='t='+strcompress(string(good_time(i)))
  plot,waves,abs(d(*,hirex_chord_index, i)),yra=[0,maxval],title=tstring,$
  charsize=1.5,xticks=2,/xstyle,xtitle='frame '+strcompress(i)
  endfor
endif else begin
  for i = x, frames-1 do begin

  tstring='t='+strcompress(string(good_time(i)))
  plot,waves,abs(d(*,hirex_chord_index, i)),yra=[0,maxval],title=tstring,$
  charsize=1.5,xticks=2,/xstyle,xtitle='frame '+strcompress(i)
  endfor
endelse
;print, frames, n_elements(good_time)
plot_type = 1
END



;**************************************************************************
PRO plot_frame
;**************************************************************************
; this procedure plots one frame

common hirex_widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, $
        w3c, w3d, w4, draw, shot_label, shot_text, shot_button, chord_slider,  $
        group_slider, frame_slider, groupplot_button, frameplot_button, $
        zoom1_button, zoom2_button, timehist_button, movie_button, from_label, from_text, $
        from_label2, to_label, to_text, to_label2, xclick_label, $
        xclick_text, xclick_label2, yclick_label, yclick_text, $
        yclick_label2, print_button, done_button, curvefit_button, $
	ylabel2, di_box, w1w, linewavelength_button, current_button
common hirex_l1,wave0,width,grating
;common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
;        w_stop, dmaxval, minval, maxval, maxval1,  tstring, d_maxval, pos, $
;        test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
;	 d_fit_time, center_fit, fwhm,print_queue_num,  hirex_chord_index
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra 
common hirex_print,x,y,tit,xtit,ytit,yra,ox1,oy1,ox2,oy2,ox3,oy3, line
common hirex_line_data, w, ident, line_plot, zoomed
common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
        w_stop, dmaxval, maxval1, maxval, minval, tstring, d_maxval, pos, $
	test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
	 d_fit_time, center_fit, print_queue_num, hirex_chord_index, hirex_sz
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type
common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type, data_type
widget_control, group_slider, sensitive = 1
widget_control, frame_slider, sensitive = 1
widget_control, zoom1_button, sensitive = 1
widget_control, zoom2_button, sensitive = 1
widget_control, linewavelength_button, sensitive = 1
widget_control, timehist_button, sensitive = 0, set_value='Time history'
widget_control, movie_button, sensitive = 1
;widget_control, curvefit_button, sensitive = 0
widget_control, frame_slider, get_value=test_fr
widget_control, print_button, sensitive = 1
;widget_control, xclick_label2, set_value='  mA  '
;widget_control, ylabel2, set_value='arb.'

pflag1 = 0
pflag2 = 0
pflag3 = 0
pflag4 = 0
ppflag = 0
again:
!p.multi=[0,1,1]
maxval1=max(d_maxval(*, 0, test_fr))
;maxval1 = max(d_maxval(*, hirex_chord_index, test_fr))
tstring=strcompress(string(good_time(test_fr)))
pflag1 = 1
plot,waves,abs(d(*,hirex_chord_index, test_fr)),title=string(shot_num)+' -- from t='+tstring+'  frame '+strcompress(test_fr)+pos(0),$
ytitle='brightness (arb.)',xtit='wavelength (mA)',yra=[0,maxval1]

; these variables are set for the print_it procedure
x = waves
y = abs(d(*,hirex_chord_index, test_fr))
tit = string(shot_num)+' -- from t='+tstring+'  frame '+strcompress(test_fr)+pos(0)
ytit = 'brightness (arb.)'
xtit = 'wavelength (mA)'
yra = [0,maxval1]

if(n_elements(where(d(*,hirex_chord_index, test_fr) lt -0.05)) gt 1) then begin
  pflag2 = 1
  oplot,waves(where(d(*,hirex_chord_index, test_fr) lt -0.05)),abs(d(where(d(*,hirex_chord_index, test_fr) lt -0.05),test_fr)),psym=2
  ox1 = waves(where(d(*,hirex_chord_index, test_fr) lt -0.05))
  oy1 = abs(d(where(d(*,hirex_chord_index, test_fr) lt -0.05),test_fr))

endif

;checks to see if wavelengths are plotted
zoomed = 0
if(n_elements(line_plot) ne 1) then line_plot = 0 
if(line_plot) then begin
  plot_line_wv
endif
plot_type= 0
;help,d
END

;***********************************************************
PRO line_wavelength
;***********************************************************
; this procedure toggles the annotation of known wavelengths

common hirex_widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, $
        w3c, w3d, w4, draw, shot_label, shot_text, shot_button, chord_slider, $
        group_slider, frame_slider, groupplot_button, frameplot_button, $
        zoom1_button, zoom2_button, timehist_button, movie_button, from_label,$
        from_text, from_label2, to_label, to_text, to_label2, xclick_label, $
        xclick_text, xclick_label2, yclick_label, yclick_text, $
        yclick_label2, print_button, done_button, curvefit_button, $
	ylabel2, di_box, w1w, linewavelength_button, current_button
common hirex_l1, wave0,width,grating
common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
        w_stop, dmaxval, minval, maxval, maxval1, tstring, d_maxval, pos, $
	test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
	d_fit_time, center_fit, fwhm,print_queue_num,  hirex_chord_index, hirex_sz
common hirex_print, x,y,tit,xtit,ytit,yra,ox1,oy1,ox2,oy2,ox3,oy3, line
common hirex_line_data, w, ident, line_plot, zoomed
common hirex_zmcp, z_mcp

if(n_elements(w) lt 2) then begin
     ;Read the wavelength table
  w=fltarr(500)
  ident=strarr(500)
  ss=0.
  sss=' '
;  if (z_mcp le -.26) Then BEGIN
;;     openr, 33, 'cmod$models:[spectroscopy.mcp]XPT.tex'
;;     print, 'opening cmod$models:[spectroscopy.mcp]XPT.tex'
;     openr, 33, '/usr/local/cmod/codes/spectroscopy/mcpher/xpt.tex'
;     print, 'opening /usr/local/cmod/codes/spectroscopy/mcpher/xpt.tex'
;  endif else BEGIN
;;     openr, 33, 'cmod$models:[spectroscopy.mcp]ar_.tex'
;;     print, 'opening cmod$models:[spectroscopy.mcp]ar_.tex'
     openr, 33, '/home/ypodpaly/Senior/WHIREX/senior_.tex'
     print, 'opening /home/ypodpaly/Senior/WHIREX/senior_.tex'
;  endelse

  i=0
  while (not eof(33)) and (i lt 500) do begin
     readf,33,ss,sss,form='(f7.2,a50)'
     w(i)=ss
     ident(i)=sss
     i=i+1
  endwhile
  close,33
  print, 'closing file'
  w=w(0:i)
  ident=ident(0:i)
endif

if(n_elements(line_plot) ne 1) then line_plot=0

if(line_plot) then line_plot=0 else line_plot=1
if(line_plot) then plot_line_wv else Begin
    if(zoomed eq 1) then hirex_zoom_in else plot_frame
endelse
END

;***********************************************************
PRO plot_line_wv
;***********************************************************
;this procedure realizes the plot of known wavelengths

common hirex_widgets,base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, $
        w3c, w3d, w4, draw, shot_label, shot_text, shot_button, chord_slider, $
        group_slider, frame_slider, groupplot_button, frameplot_button, $
        zoom1_button, zoom2_button, timehist_button, movie_button, from_label,$
        from_text, from_label2, to_label, to_text, to_label2, xclick_label, $
        xclick_text, xclick_label2, yclick_label, yclick_text, $
        yclick_label2, print_button, done_button, curvefit_button, $
        ylabel2, di_box, w1w, linewavelength_button, current_button
common hirex_print, x,y,tit,xtit,ytit,yra,ox1,oy1,ox2,oy2,ox3,oy3, line
common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
        w_stop, dmaxval, minval, maxval, maxval1, tstring, d_maxval, pos, $
	test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
	 d_fit_time, center_fit, fwhm,print_queue_num,  hirex_chord_index, hirex_sz
common hirex_l1, wave0,width,grating
common hirex_line_data, w, ident, line_plot, zoomed

 
if (zoomed eq 1) then begin
   myyra = [min([d_blow(*,test_fr),0]), max(d_blow(*,test_fr))]   
endif else begin
 if(n_elements(yra) ne 2) then begin
     myyra=[(min(y)-0.1*(max(y)-min(y))),(max(y)+0.1*(max(y)-min(y)))]
 endif else begin
     myyra=yra
 endelse
endelse
 wline=where(w gt min(x) and w lt max(x))
 if(wline(0) ne -1) then begin
   yfrac=.25
    for  i=0,n_elements(wline)-1 do begin
      oplot, [w(wline(i)),w(wline(i))], $
            [myyra(0),myyra(0)*yfrac+myyra(1)*(1.-yfrac)], lines=1
      xyouts, w(wline(i)), myyra(0)*yfrac+myyra(1)*(1.-yfrac), $
             ident(wline(i)),alignment=0.0,orientation=90.
    endfor
 endif


END

;***********************************************************
PRO hirex_zoom_in
;***********************************************************
; this procedure zooms in on a selected region of a plot.

common hirex_widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, $
        w3c, w3d, w4, draw, shot_label, shot_text, shot_button, chord_slider, $
        group_slider, frame_slider, groupplot_button, frameplot_button, $
        zoom1_button, zoom2_button, timehist_button, movie_button, from_label, from_text, $
        from_label2, to_label, to_text, to_label2, xclick_label, $
        xclick_text, xclick_label2, yclick_label, yclick_text, $
        yclick_label2, print_button, done_button, curvefit_button, $
	ylabel2, di_box, w1w, linewavelength_button, current_button
common hirex_l1,wave0,width,grating
;common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
;        w_stop, dmaxval, minval, maxval, maxval1,  tstring, d_maxval, pos, $
;        test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
;	 d_fit_time, center_fit, fwhm,print_queue_num,  hirex_chord_index, hirex_sz
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra
common hirex_print,x,y,tit,xtit,ytit,yra,ox1,oy1,ox2,oy2,ox3,oy3, line
common hirex_line_data, w, ident, line_plot, zoomed
common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
        w_stop, dmaxval, maxval1, maxval, minval, tstring, d_maxval, pos, $
	test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
	 d_fit_time, center_fit, print_queue_num, hirex_chord_index, hirex_sz
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type
common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type, data_type
pflag1 = 0
pflag2 = 0
pflag3 = 0
pflag4 = 0
ppflag = 0
widget_control, print_button, sensitive = 1
widget_control, from_text, get_value = from
widget_control, to_text, get_value = totext 
;widget_control, curvefit_button, sensitive = 1
widget_control, timehist_button, sensitive = 1

loblow=(where(waves lt from(0)))(n_elements(where(waves lt from(0)))-1)
hiblow= (where(waves gt totext(0)))(0)                       
wave_blow=waves(loblow:hiblow)
d_blow=d(loblow:hiblow, hirex_chord_index, *)
;d_blow = fltarr(hiblow-loblow+1, n_elements(d(0,0, *)))
d_blow = reform(d_blow) ; get rid of annoying one dimension

plot,wave_blow,d_blow(*,test_fr),title=string(shot_num)+' -- from t='+tstring+'  frame '+strcompress(test_fr)+pos(0),$
ytitle='brightness (arb.)',xtit='wavelength (mA)'
pflag1 = 1

x = wave_blow
y = d_blow(*,test_fr)
tit = string(shot_num)+' -- from t='+tstring+'  frame '+strcompress(test_fr)+pos(0)
ytit = 'brightness (arb.)'
xtit = 'wavelength (mA)'
yra = [0,maxval1]

;checks to see if wavelengths are plotted
zoomed = 1
if (line_plot) then plot_line_wv


END

;@spect$root:[chromex.mo]curvefit_old.pro
;*******************************************************************
PRO curve_fit
;*******************************************************************
; this procedure does curve fitting on selected peaks.


common hirex_widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, $
        w3c, w3d, w4, draw, shot_label, shot_text, shot_button, chord_slider, $
        group_slider, frame_slider, groupplot_button, frameplot_button, $
        zoom1_button, zoom2_button, timehist_button, movie_button, from_label, from_text, $
        from_label2, to_label, to_text, to_label2, xclick_label, $
        xclick_text, xclick_label2, yclick_label, yclick_text, $
        yclick_label2, print_button, done_button, curvefit_button,$
	ylabel2, di_box, w1w, linewavelength_button, current_button
common hirex_l1,wave0,width,grating
;common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
;        w_stop, dmaxval, minval, maxval, maxval1,  tstring, d_maxval, pos, $
;        test_fr,d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
;	 d_fit_time, center_fit, fwhm,print_queue_num,  hirex_chord_index, hirex_sz
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra
common hirex_print,x,y,tit,xtit,ytit,yra,ox1,oy1,ox2,oy2,ox3,oy3, line
common hirex_fits, st1, st2, fff, lllast, d_blow_ave
common hirex_times, time1, time2
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type
common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
        w_stop, dmaxval, maxval1, maxval, minval, tstring, d_maxval, pos, $
	test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
	 d_fit_time, center_fit, print_queue_num, hirex_chord_index, hirex_sz
common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type, data_type
widget_control, print_button, sensitive = 1
pflag1 = 0
pflag2 = 0
pflag3 = 0
pflag4 = 0
ppflag = 0

;if(grating eq 600) then fact1=.15 & fact2=3.2
;if(grating eq 300) then fact1=.15 & fact2=2.7

fact1 = 0. ; temporary fix 

get_times

st1 = time1(0)
st2 = time2(0)

fff=locate(good_time,st1)
lllast=locate(good_time,st2)

if (fff ne lllast) then begin
        d_blow_ave=total(d_blow(*,fff:lllast),2)/(lllast-fff+1.)
        d_blow(*,fff)=d_blow_ave
endif else begin
        d_blow_ave=d_blow(*,fff)
endelse
titvar1 = good_time(fff)
titvar2 = good_time(lllast)
        plot,wave_blow,d_blow(*,fff),$
tit='frames'+string(fff)+' thru '+string(lllast)+' from t='+$
string(form='(f7.3)',good_time(fff))+'to '+string(form='(f7.3)', good_time(lllast))+pos(0),$
xtit='wavelength (mA)', ytit='brightness (arb.)'
pflag1 = 1

x = wave_blow
y = d_blow(*,fff)
tit= 'frames'+string(fff)+' thru '+string(lllast)+' from t='+$
string(form='(f7.3)', good_time(fff))+'to '+string(form='(f7.3)',good_time(lllast))+pos(0)
xtit = 'wavelength (mA)'
ytit = 'brightness (arb.)'


jlt=where(d_blow(*,fff) lt 0)

if(jlt(0) ne -1) then begin
oplot,wave_blow(jlt),abs(d_blow(jlt,fff)),psym=2
pflag2 = 1
ox1 = wave_blow(jlt)
oy1 = abs(d_blow(jlt,fff))

endif

g=fltarr(75)
ctr=-1
 
;widget_control, curvefit_button, set_value='Done'
widget_control, di_box, /append, set_value='select peak(s)'
label1:

;event1 = widget_event(curvefit_button, /nowait)
;if (event1.id ne 0) then if event1.select then goto, label2

;event2 = widget_event(draw, /nowait)
;if (event2.id ne 0) then begin
;    if (event2.press eq 1) then begin
;    click_loc1 = convert_coord(event2.x, event2.y, /device, /to_data)
;    goto, label1a
;    endif
;endif

;goto, label1
;label1a:
wait,.25
cursor,aaa,bbb & print,aaa,bbb
wait,.25
xpk = aaa
if xpk gt wave_blow(n_elements(wave_blow)-1) then goto, label2
location = string(xpk)
ctr=ctr+1
add=2*ctr
stop
g(0+add)= d_blow_ave(locate(wave_blow,xpk))/(1.+fact1)
g(1+add)=xpk
widget_control, di_box, /append, set_value='peak selected at '+location+' mA (click to right if DONE)'
goto,label1

label2:
;widget_control, curvefit_button, set_value='Curve fit'
g(2+add)=.8
fitsum=fltarr(ctr+1,frames)
;for frs=fff,lllast do begin
frs=fff
a=g(0:(ctr+1)*2)
pvar=a
hold=d_blow(*,frs)
xhold=wave_blow
nn=float(n_elements(hold))
nnn=wave_blow(n_elements(wave_blow)-1)-wave_blow(0)
lin=(-total(hold(0:2)/3.)+total(hold(nn-3:nn-1)/3.))/nnn
offset=total(hold(0:2)/3.)-lin*wave_blow(0)
wait,.25
widget_control, di_box, /append, set_value='select background pts (both wavelength and intensity)'
cursor,aaa,bbb & print,aaa,bbb
wait,.25
mm1 = aaa
lb1 = bbb
xtemp = string(form='(f7.2)',mm1)
ytemp = string(form='(f7.2)',lb1)
widget_control, di_box, /append, $
set_value='left background pt selected at '+xtemp+' mA   y:'+ytemp+' (arb.)'

cursor,aaa,bbb & print,aaa,bbb
wait,.25
mm2 = aaa
lb2= bbb
xtemp = string(form='(f7.2)',mm2)
ytemp = string(form='(f7.2)',lb2)
widget_control, di_box, /append, $
set_value='right background pt selected at '+xtemp+' mA   y:'+ytemp+' (arb.)'

;labelback1:
;event3 = widget_event(draw, /nowait)
;if (event3.id ne 0) then begin
;    if (event3.press eq 1) then begin
;    click_loc1 = convert_coord(event3.x, event3.y, /device, /to_data)
;    goto, labelback2
;    endif
;endif
;goto, labelback1

;labelback2:
;mmm = click_loc1(1)
;mm1 = click_loc1(0) 
;lb1=mmm
;xtemp = string(form='(f7.2)',mm1)
;ytemp = string(form='(f7.2)',mmm)
;widget_control, di_box, /append, set_value='x:'+xtemp+' (A)   y:'+ytemp+' (10^18 Photons/sec/cm^2/ster/pixel)'
;wait, .5

;labelback3:

;event4 = widget_event(draw, /nowait)
;if (event4.id ne 0) then begin
;    if (event4.press eq 1) then begin
;    click_loc1 = convert_coord(event4.x, event4.y, /device, /to_data)
;    goto, labelback4
;    endif
;endif
;goto, labelback3
;
;labelback4:
;
;mmm = click_loc1(1)
;mm2 = click_loc1(0)
;
;xtemp = string(form='(f7.2)',mm2)
;ytemp = string(form='(f7.2)',mmm)
;widget_control, di_box, /append, set_value='x:'+xtemp+' (A)   y:'+ytemp+' (10^18 Photons/sec/cm^2/ster/pixel)'
;
;lb2=mmm
nnn=wave_blow(locate(wave_blow,mm2))-wave_blow(locate(wave_blow,mm1))
lin=(-lb1+lb2)/nnn
offset=lb1-lin*wave_blow(locate(wave_blow,mm1))

y1=hold-offset-lin*xhold
x1=xhold
;x1=float(indgen(n_elements(y1)))
;w=x1-x1+1.
w=y1
;
label3:
r=curvefit(x1,y1,w,a,function_name='gfunct')
;rsub=r-a(lin)-a(lin+1)*x1
counts=total(r)

pflag1 = 1
pflag2 = 0
pflag3 = 0
pflag4 = 0
ppflag = 0
plot,xhold,y1,thi=2,$
tit=string(shot_num)+' - frames '+strtrim(string(fix(fff)),2)$
+' thru '+strtrim(string(fix(lllast)),2)+' from t='+$
strtrim(string(form='(f7.3)',good_time(fff)),2)+' to '+$
strtrim(string(form='(f7.3)',good_time(lllast)),2)+pos(0), $
 xtit= 'wavelength (mA)', ytit='brightness (arb.)'

x = xhold
y = y1
tit = string(shot_num)+' - frames '+strtrim(string(fix(fff)),2)$
+' thru '+strtrim(string(fix(lllast)),2)+' from t='+$
strtrim(string(form='(f7.3)',good_time(fff)),2)+' to '+$
strtrim(string(form='(f7.2)',good_time(lllast)),2)+pos(0)
xtit = 'wavelength (mA)'
ytit = 'brightness (arb.)'

pflag2 = 1
pflag4 = 2
oplot,xhold,r
ox1 = xhold
oy1 = r

;oplot,xhold,rsub,line=5
p=fltarr(ctr+1,n_elements(hold))
psub=fltarr(ctr+1,n_elements(hold))
pvsub=fltarr(ctr+1,1000)
back=fltarr(n_elements(hold))

v=indgen(1000)
v=v-375.

for i=0,ctr do begin
addd=2*i
psub(i,*)=a(0+addd)*(exp(-(x-a(1+addd))^2/a(n_elements(a)-1)/a(n_elements(a)-1))+$
                fact1*exp(-(x-a(1+addd))^2/(fact2*a(n_elements(a)-1))^2))

oplot,xhold,psub(i,*),line=2
xyouts,a(1+addd),0.,$
'!7k!3!Io!N='+strmid(strtrim(a(1+addd),2),0,8)$
+', B(x10!E14!N)='+strmid(strtrim(total(psub(i,*)),2),0,6),$
charsiz=1.2,orient=90.

endfor
if(grating eq 600) then fwhm=a(n_elements(a)-1)*2.*sqrt(.8286)

; note that the factor of 2*sqrt(.8286) is true only for the line shape
; having the form       exp(-x^2/a^2)+.15*exp(-x^2/(3.2*a)^2)
; i.e. fact1=0.15 and fact2=3.2
if(grating eq 300) then fwhm=a(n_elements(a)-1)*2.*sqrt(.814)
; note that the factor of 2*sqrt(.814) is true only for the line shape
; having the form       exp(-x^2/a^2)+.15*exp(-x^2/(2.7*a)^2)
; i.e. fact1=0.15 and fact2=2.7

xyouts,wave_blow(10),max(y)*.8,'FWHM='+strmid(strtrim(fwhm,2),0,4),charsiz=1.2

back=offset+lin*xhold
;oplot,xhold,back
;oplot,y,psym=4
framex = string('frame ',frs)
widget_control, di_box, /append, set_value=framex
for i=0,ctr do begin
addd=2*i
framex2 = string('line at ', a(1+addd), ' has ' ,total(psub(i,*)),' counts.')
widget_control, di_box, /append, set_value=framex2
fitsum(i,frs)=total(psub(i,*))
;print,fix(frs),a(1+addd),where(psub(i,*) eq max(psub(i,*)))+star,$
;max(psub(i,*)),fwhm,total(psub(i,*)),total(y)
endfor

;print,summing/3.,tsumming/3.
summing=0.
tsumming=0.
;endif
wait,.25
fans='n'
 
END

;**************************************************************************
PRO Movie
;**************************************************************************
; this procedure plays a movie by ploting each of the frames in order.

common hirex_widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, $
        w3c, w3d, w4, draw, shot_label, shot_text, shot_button, chord_slider, $
        group_slider, frame_slider, groupplot_button, frameplot_button, $
        zoom1_button, zoom2_button, timehist_button, movie_button, from_label, from_text, $
        from_label2, to_label, to_text, to_label2, xclick_label, $
        xclick_text, xclick_label2, yclick_label, yclick_text, $
        yclick_label2, print_button, done_button, curvefit_button, $
	ylabel2, di_box, w1w, linewavelength_button, current_button
common hirex_l1,wave0,width,grating
;common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
;        w_stop, dmaxval, maxval1, maxval, minval, tstring, d_maxval, pos, $
;	test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
;	 d_fit_time, center_fit, fwhm,print_queue_num, hirex_chord_index, hirex_sz
;common hirex_flags, from_flag,  pflag1, pflag2, pflag3, pflag4, ppflag, xra
common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
        w_stop, dmaxval, maxval1, maxval, minval, tstring, d_maxval, pos, $
	test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
	 d_fit_time, center_fit, print_queue_num, hirex_chord_index, hirex_sz
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type

common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type, data_type
widget_control, shot_button, sensitive = 0
widget_control, group_slider, sensitive = 0 
widget_control, frame_slider, sensitive = 0
widget_control, groupplot_button, sensitive = 0
widget_control, frameplot_button, sensitive = 0
widget_control, zoom1_button, sensitive = 0
widget_control, zoom2_button, sensitive = 0
widget_control, timehist_button, sensitive = 0
widget_control, linewavelength_button, sensitive = 0
widget_control, print_button, sensitive = 0
widget_control, movie_button, set_value = 'Stop Movie', xsize = 160
widget_control, frame_slider, get_value = start
;widget_control, curvefit_button, sensitive = 0

widget_control, frame_slider, get_value=test_fr
pflag1 = 0
pflag2 = 0
pflag3 = 0
pflag4 = 0
ppflag = 0               
!p.multi=[0,1,1]

plot,waves,abs(d(*, hirex_chord_index, test_fr)),title=string(shot_num),$
ytitle='brightness (arb.)',xtit='wavelength (mA)',yra=[-.1 *(maxval),maxval]
pflag1 = 1

x = waves
y = abs(d(*,hirex_chord_index, test_fr))
tit = string(shot_num)
ytit = 'brightness (arb.)'
xtit = 'wavelength (mA)'
;yra = [0,maxval]
yra = [-.1*(maxval), maxval]

if(n_elements(where(d(*,hirex_chord_index, test_fr) lt -0.05)) gt 1) then begin
  pflag2 = 1
  oplot,waves(where(d(*,hirex_chord_index, test_fr) lt -0.05)),abs(d(where(d(*,hirex_chord_index, test_fr) lt -0.05),test_fr)),psym=2
  ox1 = waves(where(d(*,hirex_chord_index, test_fr) lt -0.05))
  oy1 = abs(d(where(d(*,hirex_chord_index, test_fr) lt -0.05),hirex_chord_index,test_fr))
endif
oplot,waves,abs(d(*,hirex_chord_index, test_fr)), color = 0
if(n_elements(where(d(*,hirex_chord_index, test_fr) lt -0.05)) gt 1) then $
oplot,waves(where(d(*,hirex_chord_index, test_fr) lt -0.05)),abs(d(where(d(*,hirex_chord_index, test_fr) lt -0.05),hirex_chord_index, test_fr)),psym=2, color = 0

repeat begin
for i = start, frames-3  do begin
widget_control, frame_slider, set_value = i
outstr= string(format='(f7.3)', good_time(i)) + ' sec'
oplot,waves,abs(d(*,hirex_chord_index, i))
xyouts, .75, .9, /normal,charsize = 1.5, outstr 
if(n_elements(where(d(*,hirex_chord_index, i) lt -0.05)) gt 1) then $
oplot,waves(where(d(*,hirex_chord_index, i) lt -0.05)),abs(d(where(d(*,hirex_chord_index, i) lt -0.05),hirex_chord_index, i)),psym=2

event = widget_event(movie_button, /nowait)
if (event.id ne 0) then $
    if event.select then goto, finished
wait, .3
oplot,waves,abs(d(*,hirex_chord_index, i)), color = 0
xyouts, .75, .9, /normal, charsize = 1.5, outstr, color = 0
if(n_elements(where(d(*,hirex_chord_index, i) lt -0.05)) gt 1) then $
oplot,waves(where(d(*,hirex_chord_index, i) lt -0.05)),abs(d(where(d(*,hirex_chord_index, i) lt -0.05),hirex_chord_index, i)),psym=2, color = 0

endfor
endrep until 0

finished:
widget_control, shot_button, sensitive = 1
widget_control, group_slider, sensitive = 0
widget_control, frame_slider, sensitive = 1
widget_control, groupplot_button, sensitive = 1
widget_control, frameplot_button, sensitive = 1
widget_control, zoom1_button, sensitive = 1
widget_control, zoom2_button, sensitive = 1
widget_control, linewavelength_button, sensitive = 1
widget_control, timehist_button, sensitive = 1
widget_control, movie_button, set_value = 'Movie', xsize = 160
;widget_control, curvefit_button, sensitive = 1
END

;*******************************************************************************
PRO Time_history
;*******************************************************************************
; this procedure gives a time history under a selected portion of the plot.

common hirex_widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, $
        w3c, w3d, w4, draw, shot_label, shot_text, shot_button, chord_slider, $
        group_slider, frame_slider, groupplot_button, frameplot_button, $
        zoom1_button, zoom2_button, timehist_button, movie_button, from_label, from_text, $
        from_label2, to_label, to_text, to_label2, xclick_label, $
        xclick_text, xclick_label2, yclick_label, yclick_text, $
        yclick_label2, print_button, done_button, curvefit_button, $
        ylabel2, di_box, w1w, linewavelength_button, current_button
common hirex_l1,wave0,width,grating
common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
        w_stop, dmaxval, maxval1, maxval, minval, tstring, d_maxval, pos, $
        test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
	 d_fit_time, center_fit, fwhm,print_queue_num,  hirex_chord_index, hirex_sz
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra
common hirex_print, x,y,tit,xtit,ytit,yra,ox1,oy1,ox2,oy2,ox3,oy3, line
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type
common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type, data_type
pflag1 = 0
pflag2 = 0
pflag3 = 0
pflag4 = 0
ppflag = 0
;widget_control, xclick_label2, set_value='sec'
;widget_control, ylabel2, set_value='arb.'
widget_control, timehist_button, set_value='fit tau'

widget_control, from_text, get_value = lorx
lor = lorx(0)
loblow=(where(wave_blow lt lor))(n_elements(where(wave_blow lt lor))-1)

wait,.5

widget_control, to_text, get_value = hirx
hir = hirx(0)
hiblow=(where(wave_blow gt hir))(0)

wait,.5
;center_fit=(wave_blow(hiblow)+wave_blow(loblow))/2.
wave_fit=wave_blow(loblow:hiblow)
d_fit=d_blow(loblow:hiblow,*)
center_fit=wave_fit((where(d_fit(*,test_fr) eq max(d_fit(*,test_fr))))(0))

d_fit_time=total(d_fit(*,*),1)-(d_fit(0,*)+d_fit(hiblow-loblow,*))/2.*(hiblow-loblow+1.)
d_fit_time_noback=total(abs(d_fit(*,*)),1)

plot,good_time(0:frames-1),d_fit_time,title=string(shot_num)+' -- '+strtrim(center_fit)+' A',$
xtitle='time (s)',ytitle='brightness (arb.)',xra=[good_time(0),good_time(frames-1)],$
/xstyle,yra=[0,max(d_fit_time(0:frames-3))]
pflag1 = 2
x = good_time(0:frames-1)
y = d_fit_time
tit = string(shot_num)+' -- '+strtrim(center_fit)+' mA'
xtit = 'time (s)'
ytit = 'brightness (arb.)'
yra = [0,max(d_fit_time(0:frames-3))]
xra = [good_time(0),good_time(frames-1)]
sat=intarr(frames)
for jtii=0,frames-1 do begin
        if(total(where(d_fit(*,jtii) lt -0.05)) ne -1) then sat(jtii)=1
endfor
if(total(where(sat eq 1)) ge 2 and n_elements(where(sat eq 1)) gt 1) then begin
  pflag2 = 1
  oplot,good_time(where(sat eq 1)),d_fit_time(where(sat eq 1)),psym=2
  ox1 = good_time(where(sat eq 1))
  oy1 = d_fit_time(where(sat eq 1))
endif
if(n_elements(where(sat eq 1)) eq 1 and total(where(sat eq 1)) ne -1) then begin
  pflag3 = 1
  oplot,[good_time(where(sat eq 1)),good_time(where(sat eq 1))],$
  [d_fit_time(where(sat eq 1)),d_fit_time(where(sat eq 1))],psym=2
  ox2 = [good_time(where(sat eq 1)),good_time(where(sat eq 1))]
  oy2 = [d_fit_time(where(sat eq 1)),d_fit_time(where(sat eq 1))]
endif
hrd = yes('add RF signal?')

if (hrd eq -1) then begin
  mdsopen,'rf',shot_num
  rfpower=mdsvalue('\rf_power_net')
  rftime=mdsvalue('dim_of(\rf_power_net)')
  ;rfpower=mdsvalue('\rf::top.antenna.results.d_port:pwr_net')
  ;rftime=mdsvalue('dim_of(\rf::top.antenna.results.d_port:pwr_net)')
  mdsclose

  oplot,rftime,rfpower/max(rfpower)*max(d_fit_time)/3.,line=2
  pflag4 = 1
  ox3 = rftime
  oy3 = rfpower/max(rfpower)*max(d_fit_time)/3
  line = 2
endif

END

;******************************************************************************
PRO tau_fit
;******************************************************************************

common hirex_widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, $
        w3c, w3d, w4, draw, shot_label, shot_text, shot_button, chord_slider, $
        group_slider, frame_slider, groupplot_button, frameplot_button, $
        zoom1_button, zoom2_button, timehist_button, movie_button, from_label, from_text, $
        from_label2, to_label, to_text, to_label2, xclick_label, $
        xclick_text, xclick_label2, yclick_label, yclick_text, $
        yclick_label2, print_button, done_button, curvefit_button, $
        ylabel2, di_box, w1w, linewavelength_button, current_button
common hirex_l1,wave0,width,grating
common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
        w_stop, dmaxval, maxval1, maxval, minval, tstring, d_maxval, pos, $
        test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
	 d_fit_time, center_fit, fwhm,print_queue_num,  hirex_chord_index, hirex_sz
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra
common hirex_print, prntvl, prntvl2, prntvl3, prntvl4
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type
common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type, data_type
pflag1 = 0
pflag2 = 0
pflag3 = 0
pflag4 = 0
ppflag = 0
time_blow=good_time(0:frames-1)
xdata_blow=d_fit_time
plot,time_blow,xdata_blow
pflag1 = 1
x = time_blow
y = xdata_blow
tit = ''
xtit = ''
ytit = ''

wait,.5

widget_control, di_box, /append, set_value='blow up range...'

label1tau:

event2 = widget_event(draw, /nowait)
if (event2.id ne 0) then begin
    if (event2.press eq 1) then begin
    click_loc1 = convert_coord(event2.x, event2.y, /device, /to_data)
    goto, label1aa
    endif
endif

goto, label1tau
label1aa:

tlow = click_loc1(0)
yyy = click_loc1(1)
wait, .5
widget_control, di_box, /append, set_value= 'OK left side...'
;if (tlow gt max(time_blow)) then goto, no_more_blow

label2tau:

event1 = widget_event(draw, /nowait)
if (event1.id ne 0) then begin
    if (event1.press eq 1) then begin
    click_loc1 = convert_coord(event1.x, event1.y, /device, /to_data);
    goto, label2aa
    endif
endif

goto, label2tau

label2aa:
thigh = click_loc1(0)
yyy = click_loc1(1)
wait,.5
widget_control, di_box, /append, set_value= 'OK right side...'
jt=min(abs(time_blow-tlow),plow)
jt1=min(abs(time_blow-thigh),phigh)

xdata_blow=xdata_blow(plow:phigh)
time_blow=time_blow(plow:phigh)
plot,time_blow,xdata_blow
pflag1 = 1 
x = time_blow
y = xdata_blow
tit = ''
xtit = ''
ytit = ''

widget_control, di_box, /append , set_value= 'channel 1'
widget_control, di_box, /append, set_value='select background end points...'

label3tau:

event3 = widget_event(draw, /nowait)
if (event3.id ne 0) then begin
    if (event3.press eq 1) then begin
    click_loc1 = convert_coord(event3.x, event3.y, /device, /to_data)
    goto, label3aa
    endif
endif

goto, label3tau
label3aa:

bklow = click_loc1(0)
yyy1 = click_loc1(1)
wait,.5
if(bklow gt max(time_blow)) then goto,skipping
widget_control, di_box,/append, set_value= 'OK left side...'

label4tau:

event4 = widget_event(draw, /nowait)
if (event4.id ne 0) then begin
    if (event4.press eq 1) then begin
    click_loc1 = convert_coord(event4.x, event4.y, /device, /to_data)
    goto, label4aa
    endif
endif

goto, label4tau
label4aa:

bkhigh = click_loc1(0)
yyy2 = click_loc1(1)
wait,.5
widget_control, di_box, /append, set_value= 'OK right side...'
jt=min(abs(time_blow-bklow),low_back)
jt1=min(abs(time_blow-bkhigh),high_back)
left_ave=yyy1
right_ave=yyy2
x_slope=(right_ave-left_ave)/(high_back-low_back)
x_bkgnd=x_slope*findgen(n_elements(xdata_blow))+left_ave-x_slope*(low_back-2.)
oplot,time_blow,x_bkgnd
pflag2 = 1
ox1 = time_blow
oy1 = x_bkgnd

label5tau:

event5 = widget_event(draw, /nowait)
if (event5.id ne 0) then begin
    if (event5.press eq 1) then begin
    click_loc1 = convert_coord(event5.x, event5.y, /device, /to_data)
    goto, label5aa
    endif
endif

goto, label5tau
label5aa:

mm = click_loc1(0)
nn = click_loc1(1)

wait,.5
xdata_corr=(xdata_blow-x_bkgnd)(low_back:high_back)
time_blow=time_blow(low_back:high_back)
plot,time_blow,xdata_corr
pflag1 = 1
x = time_blow
y = xdata_corr
tit = ''
xtit = ''
ytit = ''

pflag2 = 0
pflag3 = 0
pflag4 = 0

xdata_smoo=xdata_corr
;xdata_smoo=smooth(xdata_corr,25)
oplot,time_blow,xdata_smoo
pflag2 = 1
ox1 = time_blow
oy1 = xdata_smoo
widget_control, di_box, /append, set_value='select fit end points ...'

label6tau:


event6 = widget_event(draw, /nowait)
if (event6.id ne 0) then begin
    if (event6.press eq 1) then begin
    click_loc1 = convert_coord(event6.x, event6.y, /device, /to_data);
    goto, label6aa
    endif
endif

goto, label6tau
label6aa:
lf = click_loc1(0)
yyy = click_loc1(1)

wait,.5
widget_control, di_box, /append, set_value='OK left side...'

label7tau:

event7 = widget_event(draw, /nowait)
if (event7.id ne 0) then begin
    if (event7.press eq 1) then begin
    click_loc1 = convert_coord(event7.x, event7.y, /device, /to_data)
    goto, label7aa
    endif
endif

goto, label7tau
label7aa:

hf = click_loc1(0)
yyy = click_loc1(1)


wait,.5
widget_control, di_box, /append, set_value='OK right side...'
jt=min(abs(time_blow-lf),low_fit)
jt1=min(abs(time_blow-hf),high_fit)
plot_io,time_blow(low_fit:high_fit),xdata_corr(low_fit:high_fit)
pflag1 = 1
x = time_blow(low_fit:high_fit)
y = xdata_corr(low_fit:high_fit)

oplot,time_blow(low_fit:high_fit),xdata_smoo(low_fit:high_fit)
pflag2 = 1
x = time_blow(low_fit:high_fit)
y = xdata_smoo(low_fit:high_fit)
xfitable=alog(xdata_smoo(low_fit:high_fit))
time_fitable=time_blow(low_fit:high_fit)
xcoeffs=poly_fit(time_fitable,xfitable,1,yfitable)
plot,time_blow,xdata_smoo,title=string(shot_num)+' -- '+strtrim(center_fit)+' mA'+pos(0),$
        xtit=' time (s)'
pflag1 = 1
x = time_blow
y = xdata_smoo
tit = string(shot_num)+' -- '+strtrim(center_fit)+' mA'+pos(0)
xtit = ' time (s)'
ytit = ''


oplot, time_fitable,exp(yfitable)
pflag2 = 1
ox1 = time_fitable
oy1 = exp(yfitable)

xyouts,time_blow(10),exp(yfitable(0)),'!4s!3!dimp!N='+strtrim(fix((-1./xcoeffs(1)*1000.)),2)+' ms',charsiz=2.
widget_control, di_box, /append, set_value='measured tau for channel 1 = '+strtrim(-1./xcoeffs(1)*1000.,2)
meas_tau=-1./xcoeffs(1)*1000.

skipping:
end

;**************************************************************************
PRO hirex_print_it
;**************************************************************************
; this function should print whatever is in the draw widget. However, all cases
; have not been tested.

common hirex_widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, $
        w3c, w3d, w4, draw, shot_label, shot_text, shot_button, chord_slider,  $
        group_slider, frame_slider, groupplot_button, frameplot_button, $
        zoom1_button, zoom2_button, timehist_button, movie_button, from_label, from_text, $
        from_label2, to_label, to_text, to_label2, xclick_label, $
        xclick_text, xclick_label2, yclick_label, yclick_text, $
        yclick_label2, print_button, done_button, curvefit_button, $
        ylabel2, di_box, w1w, linewavelength_button, current_button
common hirex_l1,wave0,width,grating
common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
        w_stop, dmaxval, minval, maxval, maxval1,  tstring, d_maxval, pos, $
        test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
	 d_fit_time, center_fit,fwhm,print_queue_num,  hirex_chord_index, hirex_sz
;common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra
common hirex_print,x,y,tit,xtit,ytit,yra,ox1,oy1,ox2,oy2,ox3,oy3, line
common hirex_line_data, w, ident, line_plot, zoomed

common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type, data_type



ppflag = 1

tmp = ask('Which printer number?[114 is cmod control room]')
tmp = strtrim(tmp,2)
if (tmp(0) eq '') then begin
    if (print_queue_num eq '') then print_queue_num = 114
    tmp(0) = print_queue_num
  endif else begin
    print_queue_num = tmp(0)
  endelse
set_plot, 'ps'
device, /landscape

if (pflag1 eq 1) then begin
  plot, x, y, tit= tit, xtit=xtit, ytit=ytit, yrange=[0,max(y)]
endif

if (pflag1 eq 2) then begin
 plot, x, y, tit= tit, xtit=xtit, ytit=ytit, yra= yra, xrange=xra, /xstyle
endif

if (pflag2 eq 1) then begin
  oplot, ox1, oy1, psym=2
endif

if (pflag3 eq 1) then begin
  oplot, ox2, oy2, psym = 2
endif

if (pflag4 eq 1) then begin
  oplot, ox3, oy3, line = line
endif

if (pflag4 eq 2) then begin
  a = pvar
  for i=0,ctr do begin
    addd=2*i
    oplot,xhold,psub(i,*),line=2
    xyouts,a(1+addd),0.,$
    '!7k!3!Io!N='+strmid(strtrim(a(1+addd),2),0,8)$
    +', B(x10!E14!N)='+strmid(strtrim(total(psub(i,*)),2),0,6),$
    charsiz=1.2,orient=90.
  endfor
  xyouts,wave_blow(10),max(psub)*.8,'FWHM='+strmid(strtrim(fwhm,2),0,4),charsiz=1.2
endif

if (line_plot) then plot_line_wv

device, /close

spawn, 'lpr -Pp'+tmp(0)+' idl.ps'
; /param=(data=ps, page_o=landscape) /que=p'+tmp(0)
;set_plot, 'x'
xwplot
end

;**************************************************************************
PRO hirex_look_event, ev
;**************************************************************************
; event handler for the main widget.

common hirex_widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, $
        w3c, w3d, w4, draw, shot_label, shot_text, shot_button, chord_slider, $
        group_slider, frame_slider, groupplot_button, frameplot_button, $
        zoom1_button, zoom2_button, timehist_button, movie_button, from_label, from_text, $
        from_label2, to_label, to_text, to_label2, xclick_label, $
        xclick_text, xclick_label2, yclick_label, yclick_text, $
        yclick_label2, print_button, done_button, curvefit_button, $
	ylabel2, di_box, w1w, linewavelength_button, current_button
common hirex_l1,wave0,width,grating
common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
        w_stop, dmaxval, maxval1, maxval, minval, tstring, d_maxval, pos, $
	test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
	 d_fit_time, center_fit, print_queue_num, hirex_chord_index, hirex_sz

common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type, data_type
type = tag_names(ev, /structure)
;help,ev,/str

case type of
	'WIDGET_SLIDER': begin
            widget_control, ev.id, get_uvalue = slide
       	    case slide of
               'sz': begin
                        widget_control, chord_slider, get_value = temp_slider
                        hirex_chord_index = temp_slider - 1
                        z_mcp = hirex_sz(hirex_chord_index)
                        if z_mcp GT 0 then begin
                            pos = ' Line-of-sight at SZ = ' + strtrim(string(z_mcp), 2) 
                            endif else begin
                            pos = ' Line-of-sight at SZ = -' + strtrim(string( abs(z_mcp)), 2) 
                        endelse
                        widget_control, yclick_text, set_value= strcompress(string(z_mcp), /remove_all)
                        if plot_type EQ 0 then plot_frame
                        if plot_type EQ 1 then plot_group
                        end
	       'group': plot_group
               'frame': plot_frame
               ELSE: message, 'UNHANDLED SLIDER'
            endcase
                         end

	'WIDGET_BUTTON': begin
            widget_control, ev.id, get_value = button
	    case button of
	        'Quit': widget_control, /destroy, ev.top
                'Debug': stop ; a debug option
	       'Apply': begin
                        if data_type EQ 0 then  hirex_get_shot, type = 0
                        if data_type EQ 1 then hirex_get_shot, type = 1
	                plot_group
                        end
      'Set To Current': begin
                        cshot = mdsvalue('current_shot("cmod")')
                        widget_control, shot_text, set_value = strtrim(cshot,2)
                        end
          'Plot frame': plot_frame
          'Plot group': plot_group
      'Identify Lines': line_wavelength
        'Time history': Time_history
               'Movie': Movie
             'Zoom In': hirex_zoom_in
            'Zoom Out': plot_frame
               'Print': hirex_print_it
           ;'Curve fit': curve_fit
             'fit tau': tau_fit 
                  'He': begin
                      data_type = 0
                      hirex_get_shot, type = 0
                      if plot_type EQ 0 then plot_frame
                      if plot_type EQ 1 then plot_group
                  end
                   'H': begin
                       data_type = 1
                       hirex_get_shot, type = 1
                       if plot_type EQ 0 then plot_frame
                       if plot_type EQ 1 then plot_group
                   end

             ELSE: message, 'UNHANDLED BUTTON'
                       endcase         
                          end

        'WIDGET_DRAW': begin
           format_str = '(f7.2)'
           click_loc = convert_coord(ev.x, ev.y, /device, /to_data)

       	 if (ev.press eq 1) then $
	    if not from_flag then begin
		widget_control, from_text, set_value=string(form=format_str, click_loc(0))
		from_flag = 1
	    endif else begin
		widget_control, to_text, set_value=string(form=format_str, click_loc(0))
		from_flag = 0
	    endelse
	   
         if ((not ev.press) and (not ev.release)) then begin
	   ;widget_control, xclick_text, set_value=string(form = format_str, click_loc(0))
           ;widget_control, yclick_text, set_value=string(click_loc(1))
         endif
end

	'WIDGET_TEXT_CH': begin
              if data_type EQ 0 then  hirex_get_shot, type = 0
              if data_type EQ 1 then hirex_get_shot, type = 1
              plot_group
          end

          ELSE: message, 'UNHANDLED EXCEPTION'
endcase

end

;***************************************************************************
;                    M A I N      P R O G R A M
;***************************************************************************
; This is the part of the program that gets called from the idl prompt.
; This is where the widgets are built and realized.
;------------
PRO hirex_look, shot = givenshot, debug = debug
;------------


common hirex_widgets, base, lcol, rcol, w1, w1a, w1b, w1c, w2, w3, w3a, w3b, $
	w3c, w3d, w4, draw, shot_label, shot_text, shot_button, chord_slider, $
	group_slider, frame_slider, groupplot_button, frameplot_button, $
	zoom1_button, zoom2_button, timehist_button, movie_button, from_label, from_text, $
	from_label2, to_label, to_text, to_label2, xclick_label, $
	xclick_text, xclick_label2, yclick_label, yclick_text, $
	yclick_label2, print_button, done_button, curvefit_button, $
	ylabel2, di_box, w1w, linewavelength_button, current_button
common hirex_l1,wave0,width,grating
common hirex_data_block, shot_num, d, good_time, waves, jackpos, frames, w_start, $
        w_stop, dmaxval, maxval1, maxval, minval, tstring, d_maxval, pos, $
	test_fr, d_blow, wave_blow, fact1,  ctr, fact2, psub, pvar, xhold,$
	 d_fit_time, center_fit, print_queue_num, hirex_chord_index, hirex_sz

common hirex_flags, from_flag, pflag1, pflag2, pflag3, pflag4, ppflag, xra, plot_type, data_type

xwplot
hirex_chord_index = 23 ; initialize the chord to the central
plot_type = 1 ; initialize to group
data_type = 0; initialize to He

shotstr = string(form='(i9)', mdscur_shot("cmod"))
if keyword_set(debug) then shotstr = string(1070803011) ; working shot
if n_elements(givenshot) GT 0 then shotstr = strcompress(string(givenshot), /remove_all)
;window, color=2
;wdelete

base = widget_base(title='HIREX_LOOK', /row)

; the left column is for buttons, and cursors.
  lcol = widget_base(base, /column, space = 28)

; the right column has the draw window.
  rcol = widget_base(base, /frame, /column)
;    draw = widget_draw(rcol, xsize=800, ysize = 700, retain=2, /motion_events, $
;                   /button_events)
    draw = widget_draw(rcol, xsize=800, ysize = 700, retain=2, $
                   /button_events)
    w1w = widget_base(rcol, /column)
      di_box = widget_text(w1w, xsize = 95, ysize = 3, /wrap, /scroll, value='This is a dialogue box')
    w1 = widget_base(lcol, /column, /frame)
      w1a = widget_base(w1, /row)
        shot_label = widget_label(w1a, value='Shot:')
        shot_text = widget_text(w1a, /editable, xsize = 10, ysize = 1, value=shotstr)
    current_button  = widget_button(w1, value='Set To Current', xsize = 180)

; adding new data source widget buttons
    use_mcp_bri=intarr(2)
    w11=widget_base(lcol, /frame,/col,/EXCLUSIVE)
;;       c=widget_label(w11,value='use MOD spectral data')
       use_mcp_bri(0) = widget_button(w11,frame=1 ,value='He',uvalue=1, sensitive = 1)
;       c=widget_label(w11,value='use VIRGIN spectral data')
       use_mcp_bri(1) = widget_button(w11,frame=1 ,value='H',uvalue=0, sensitive = 1)
widget_control, use_mcp_bri(0), set_button = 1 
;;;;;;;;;;;;;;;

       shot_button = widget_button(w1a, value='Apply', xsize = 80)
      w1chord = widget_base(w1, /row)
         chord_slider = widget_slider(w1chord, min = 1, max = 48, title = '(SZ)', $
                                      xsize = 230, value = 24, uvalue = 'sz')
      w1b = widget_base(w1, /row)
        group_slider = widget_slider(w1b, min=1, max=150, title='(Group)',$
                            xsize= 230, uvalue = 'group')
      w1c = widget_base(w1, /row)
        frame_slider = widget_slider(w1c, min=1, max=150, title='(Frame)', $
                             xsize=230, uvalue = 'frame') 

    w2 = widget_base(lcol, /column, /frame)
      groupplot_button = widget_button(w2, value='Plot group', xsize=230)
      frameplot_button = widget_button(w2, value='Plot frame', xsize=230)
      zoom1_button = widget_button(w2, value='Zoom In', xsize=230)
      zoom2_button = widget_button(w2, value='Zoom Out', xsize=230)
      linewavelength_button = widget_button(w2, value='Identify Lines',$
                              xsize=230)
      ;curvefit_button = widget_button(w2, value='Curve fit', xsize=230)
      timehist_button = widget_button(w2, value='Time history',xsize=230)
      movie_button = widget_button(w2, value='Movie', xsize=230)

    w3 = widget_base(lcol, /column, /frame)
      w3a = widget_base(w3, /row)
        from_label = widget_label(w3a, value='   From:')
        from_text = widget_text(w3a, value='', xsize=7, ysize=1, /editable)
        from_label2 = widget_label(w3a, value='A')
      w3b = widget_base(w3, /row)
        to_label = widget_label(w3b, value='       To:')
        to_text = widget_text(w3b, value='', xsize=7, ysize=1, /editable)
        to_label2 = widget_label(w3b, value='mA')
      ;w3c = widget_base(w3, /row)
        ;xclick_label = widget_label(w3c, value='Cursor:')
        ;xclick_text = widget_text(w3c, value='', xsize=7, ysize=1)
        ;xclick_label2 = widget_label(w3c, value='  mA  ')
      w3d = widget_base(w3, /row)
        yclick_label = widget_label(w3d, value='     ')    
        yclick_label2 = widget_label(w3d, value= 'SZ')
        yclick_text = widget_text(w3d, value='', xsize= 7, ysize=1)
        ylabel2 = widget_label(w3d, value='cm.')

    w4 = widget_base(lcol, /column)
    if keyword_set(debug) then begin
        debug_button = widget_button(w4, value = 'Debug', xsize = 230)
    endif

      print_button = widget_button(w4, value='Print', xsize=230)
      done_button = widget_button(w4, value='Quit', xsize = 230)

widget_control, base, /realize
widget_control, get_value=draw_window, draw

; fix the widget_sliders
widget_control, frame_slider, set_slider_min = 0, set_value = 0
widget_control, group_slider, set_slider_min = 0, set_value = 0

wset, draw_window

from_flag = 0

if keyword_set(debug) then begin
hirex_get_shot, type = 0
plot_group
endif

xmanager, 'hirex_look', base, event_handler='hirex_look_event'
end

hirex_look

end

