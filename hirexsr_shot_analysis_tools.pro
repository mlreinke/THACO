;+
;NAME:
;	HIREXSR_IS_PLASMA
;PURPOSE:
;	This procedure checks if a shot is a plasma
;
;CALLING SEQUENCE
;	CHECK = HIREXSR_IS_PLASMA(shot, threshold = threshold)
;
;INPUTS:
;	shot        LONG  shot number
;       threshold   FLOAT mininum plasma current to exceed
;	
;OUPUTS:
;       Returns true(1) or false(0)
;
;
;MODIFICATION HISTORY:
;	Written by: 	Y. Podpaly  12/15/2010
;
;-


;note: basec on is_plasma by AIC
function HIREXSR_IS_PLASMA, shot, threshold = threshold

if keyword_set(threshold) then begin
    plasma_min = threshold 
endif else begin
    plasma_min = 0.2 ; MA
endelse

mdsopen, 'magnetics', shot, quiet=1
   ip   = mdsvalue('\MAGNETICS::TOP.PROCESSED.CURRENT_DATA:IP',status =status, quiet=1)
   ip_t = mdsvalue('dim_of(\MAGNETICS::TOP.PROCESSED.CURRENT_DATA:IP)',quiet=1)
mdsclose, quiet = 1

if n_elements(ip) EQ 0 then begin
    ; no ip data
    return, 0
endif

ip = abs(ip)*1e-6 ; convert to MA
if  max(ip) GE plasma_min then return, 1
return, 0

end



;+
;NAME:
;	HIREXSR_AR_EXIST
;PURPOSE:
;	This procedure checks if the argon signal is high enough to
;	justify running analysis codes. Checks how many ms the argon
;	is on B_side_lower or H_bottom when it is above a certain
;	number of standard devs about the mean
;
;CALLING SEQUENCE
;	CHECK = HIREXSR_AR_LOW(shot, time_thres = time_thresh,
;	lev_thresh = lev_thresh, volt_thresh=volt_thresh)
;
;INPUTS:
;	shot        LONG  shot number
;       time_thresh LONG  minimum number of ms for valve to be open (default - 35)
;	lev_thresh  LONG  minimum % to continue on (default - 30%)
;       volt_thresh LONG  minimum voltage to consider open (default -2)
;
;OUPUTS:
;       Returns true(1) or false(0)
;
;
;MODIFICATION HISTORY:
;	Written by: 	Y. Podpaly  12/15/2010]
;       Modified YP 12/21/2010 name changed to HIREXSR_AR_EXIST and
;       threshold check changed
;
;-

function HIREXSR_AR_EXIST, shot, time_thresh = time_thresh, lev_thresh = lev_thresh, volt_thresh=volt_thresh

if keyword_set(time_thresh) then min_time= time_thresh else min_time = 10 ; 10 ms min time
if keyword_set(lev_thresh) then min_perc= lev_thresh else min_perc = 30 ; 30%
if keyword_set(volt_thresh) then min_volt = volt_thresh else min_volt =2; 2 V


argon_in_bsdlo = 0
argon_in_hbot = 0

; checks if there is argon in B-side lower or H-bottom

mdsopen, 'engineering', LONG(shot)
label_bsdlo = mdsvalue('\ENGINEERING::TOP.TORVAC.GAS.PVALVE_3:LABEL')
label_hbot = mdsvalue('\ENGINEERING::TOP.TORVAC.GAS.PVALVE_5:LABEL')
mdsclose

; check b-side-lower
split = strsplit(label_bsdlo, /extract)
indarg = where(split EQ 'ARGON')
perc_argon = LONG(split(indarg+1))
if perc_argon GT min_perc then begin
    argon_in_bsdlo = 1
endif

; check h_bot
split = strsplit(label_hbot, /extract)
indarg = where(split EQ 'ARGON')
perc_argon = LONG(split(indarg+1))
if perc_argon GT min_perc then begin
    argon_in_hbot = 1
endif

; check for valve openings

if argon_in_bsdlo EQ 1 then begin
    mdsopen, 'spectroscopy', LONG(shot)
    waveform = mdsvalue('\SPECTROSCOPY::TOP.X_RAY_PHA:INCAA16:GAS_B_SD_LO')
    t_wave = mdsvalue('dim_of(\SPECTROSCOPY::TOP.X_RAY_PHA:INCAA16:GAS_B_SD_LO,0)')
    mdsclose
    
    sub_inds = where(t_wave GT 0. and t_wave LT 2.5) ; in shot
    wave = waveform(sub_inds)
    time = t_wave(sub_inds)
    
    on_indices = where(wave GT min_volt) ; threshold for on is above 3 standard devs
    time_on = n_elements(on_indices)*(time(1)-time(0))*1000. ; change to ms
    
    if time_on GT min_time then return, 1    
endif

if argon_in_hbot EQ 1 then begin
    mdsopen, 'spectroscopy', LONG(shot)
    waveform = mdsvalue('\SPECTROSCOPY::TOP.X_RAY_PHA:INCAA16:INPUT_05')
    t_wave = mdsvalue('dim_of(\SPECTROSCOPY::TOP.X_RAY_PHA:INCAA16:INPUT_05,0)')
    mdsclose
    
    sub_inds = where(t_wave GT 0 and t_wave LT 2.5) ; in shot
    wave = waveform(sub_inds)
    time = t_wave(sub_inds)
    
    on_indices = where(wave GT min_volt) ; threshold for on is above 3 standard devs
    time_on = n_elements(on_indices)*(time(1)-time(0))*1000. ; change to ms
    
    if time_on GT min_time then return, 1    
endif

return, 0
end

;+
;NAME:
;	HIREXSR_AUTO_TMAP
;PURPOSE:
;	This procedure returns a reasonable tmap and good channel
;	description. It is an automated first cut approach to this
;	process, so it does not do any data analysis. The initial good
;	array must be a real one, so that the code can grab a basic
;	'channels on' guess. 
;
;CALLING SEQUENCE
;	HIREXSR_AUTO_TMAP, shot, tmap, good, time = time, multi=multi
;
;INPUTS:
;	shot        LONG  shot number
;       tmap        INT[] array for the time maps
;
;
;OPTIONAL INPUTS:
;       time    FLOAT[] 2 element array for beginning and ending time
;                       of the shot defaults (0.5-1.7)
;       multi   INT     how many adjacent times to bin together
;       good    LONG[] array for the good channels MUST be a real
;                          good array from a previous shot (default = 1)
;       ipthresh FLOAT minimum acceptable current for disruption checking
;
;      
;KEYWORD PARAMETERS:
;	h	/h will use the default MAX_ALLOWED_CHANNELS=32 rather than 96
;
;OUPUTS:
;      
;       NONE
;
;LIMITATIONS:
;    - Currently requires an initial good channel guess
;    - Does NOT handle hydrogen-like module yet
;
;MODIFICATION HISTORY:
;	Written by: 	Y. Podpaly  01/10/2011
;	1/17/10:	M.L. Reinke - modified data check to use GETNCI to run quicker.  Added the /h keyword
;       1/20/11         Y. Podpaly - added disruption check to auto tmap
;-


pro HIREXSR_AUTO_TMAP, shot, tmap, good=good, time=time, multi = multi,h=h, ipthresh = ipthresh

; hard definitions
MINSHOTTIME = 0.0 ; start of discharge (s)
MAXSHOTTIME = 2.5 ; end of discharge  (s)
DEF_TIME_MIN = 0.5; default init time (s)
DEF_TIME_MAX = 1.7; default end time (s)
DEF_MULTI = 1 ; default time binning
IP_DEFAULT = 5.*10^4. ; 50 kA minimum plasma current
IF keyword_set(h) THEN MAX_ALLOWED_CHANNELS=32 ELSE MAX_ALLOWED_CHANNELS = 96 ; preset in Matt Reinke's codes

; check for disruptions
if keyword_set(ipthresh) then begin
    if ipthresh GT 0. then min_curr = ipthresh
endif else begin
    min_curr = IP_DEFAULT
endelse

mdsopen, 'magnetics', shot, quiet=1
   ip   = smooth(abs(mdsvalue('\MAGNETICS::TOP.PROCESSED.CURRENT_DATA:IP',status =status, quiet=1)), 5)
   ip_t = mdsvalue('dim_of(\MAGNETICS::TOP.PROCESSED.CURRENT_DATA:IP)',quiet=1)
mdsclose, quiet = 1

sub_ip = ip(where(ip_t GE DEF_TIME_MIN and ip_t LE DEF_TIME_MAX))
sub_ip_t = ip_t(where(ip_t GE DEF_TIME_MIN and ip_t LE DEF_TIME_MAX))

non_plasma = where(sub_ip LT min_curr)
if n_elements(non_plasma) GT 1 then begin
    time_disrupt = sub_ip_t(non_plasma(0))
    ;redefine DEF_TIME_MAX
    DEF_TIME_MAX = time_disrupt
endif


; get optional inputs first

; timing
if keyword_set(time) then begin
    if n_elements(time) EQ 2 then begin
        if time(0) GE MINSHOTTIME and time(0) LE MAXSHOTTIME then begin
            init_time = time(0)
        endif else begin
            print, 'INVALID TIMES SPECIFIED'
            return
        endelse
        
        if time(1) GE MINSHOTTIME and time(1) LE MAXSHOTTIME then begin
            end_time = time(1)
        endif else begin
            print, 'INVALID TIMES SPECIFIED'
            return
        endelse    
    endif else begin
        print, 'INVALID TIME ARRAY'
        return
    endelse
endif else begin
    init_time = DEF_TIME_MIN
    end_time = DEF_TIME_MAX
endelse

;time binning
if keyword_set(multi) then begin
    time_bin = FIX(multi)
endif else begin
    time_bin = FIX(DEF_MULTI)
endelse

; now get the shot information

mdsopen, 'spectroscopy', shot
;mod1 = mdsvalue('\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD1')
mod1=mdsvalue('getnci(\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD1, "length")')
;mod2 = mdsvalue('\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD2')
mod2=mdsvalue('getnci(\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD2, "length")')
;mod3 = mdsvalue('\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD3')
mod3=mdsvalue('getnci(\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD3, "length")')
;mod4 = mdsvalue('\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD4')
mod4=mdsvalue('getnci(\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD4, "length")')

mdsclose


;if not(array_equal(size(mod1), size(mod2)) and array_equal(size(mod2), size(mod3)) and array_equal(size(mod3), size(mod4))) then begin
print,'image sizes: '
print,mod1,mod2,mod3,mod4
IF mod1-mod2 NE 0 OR mod1-mod3 NE 0 OR mod1-mod4 NE 0 THEN BEGIN
    print, 'WARNING: Likely module malfunction on this shot'
ENDIF


srtimes = hirexsr_get_time(shot, /quiet)


if n_elements(srtimes) EQ  1 then begin
    print, 'Problem getting times'
    return
endif

number_times = n_elements(srtimes)

tmap = lonarr(number_times)-1 ; initialize tmap
valid_indices = where(srtimes GE init_time and srtimes LE end_time)
start_index = valid_indices(0)
end_index = valid_indices(n_elements(valid_indices)-1)

; deal with good values
if keyword_set(good) then begin
    newgood = lonarr(MAX_ALLOWED_CHANNELS, number_times) ;  preset new good values
    channelson = good(*, 0) ; grab the first time point from good which is *almost* always on
endif


; cycle through and create the tmap
count = 0 ; counter index
for index = start_index, end_index do begin
    tmap(index) = FIX(FIX(index-start_index)/time_bin)
    if keyword_set(good) then newgood(*, index-start_index) = channelson
end

if keyword_set(good) then good = newgood

end


;+
;NAME:
; HIREXSR_CHECK_SIZE
;-
function HIREXSR_CHECK_SIZE, shot
mdsopen, 'spectroscopy', shot
mod1=mdsvalue('getnci(\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD1, "length")')
mod2=mdsvalue('getnci(\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD2, "length")')
mod3=mdsvalue('getnci(\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD3, "length")')
mod4=mdsvalue('getnci(\SPECTROSCOPY::TOP.HIREX_SR.RAW_DATA:MOD4, "length")')

print,'image sizes: '
print,mod1,mod2,mod3,mod4
IF mod1-mod2 NE 0 OR mod1-mod3 NE 0 OR mod1-mod4 NE 0 THEN BEGIN
    ;print, 'Likely module malfunction on this shot'
    return, 0
ENDIF 
    
return, 1

end
;+
;NAME:
; HIREXSR_FIND_DEC
;-
pro HIREXSR_FIND_DEC, shot, timerange, dec_arr, beta = beta

test = hirexsr_check_size(shot)
if test EQ 0 then begin
    print, 'Likely module malfunction: returning'
    return
endif

if timerange(1) LE timerange(0) then begin
    print, 'Invalid Time Range'
    return
endif

mdsopen, 'spectroscopy', shot
original_str = mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.INFO.MOD2.MIRROR:ROT')
time = mdsvalue('dim_of(\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.PROFILES.Z:PRO, 1)')
mdsclose

original_beta = original_str(1)

beta_error_arr = fltarr(n_elements(dec_arr)) 

for index = 0, n_elements(dec_arr) -1 do begin
    if keyword_set(beta) then begin
        hirexsr_change_crystal_beta, shot, dec_arr(index)
    endif else begin
        hirexsr_change_crystal_decl, shot, dec_arr(index)
    endelse
    
    ;hirexsr_invert_data, shot, 0 ; w
    ;hirexsr_invert_data, shot, 1 ; x
    hirexsr_invert2tree, shot, 2 ; z
    
    ; now get the fit and moment
    hirexsr_load_momentptr, shot, 2, momptr,tau,pos,tpos,lam_o,z
    hirexsr_load_binning, shot, chmao, tch, tmap, good, chmax

    lowindex = where(time GE timerange(0))
    lowindex = lowindex(0)
    highindex =  where (time LE timerange(1))
    highindex = highindex(n_elements(highindex)-1)
    
    ; check for inversion reality
    if tmap(lowindex) EQ -1 then begin
        lowindex = where(tmap EQ 0)
        lowindex = lowindex(0)
    endif
    if tmap(highindex) EQ -1 then begin
        highindex = where(tmap EQ max(tmap))
        highindex = highindex(n_elements(highindex)-1)
    endif
    
    lowtime = tmap(lowindex)
    hightime = tmap(highindex)

    error = 0
    badframe = 0
    for interindex = lowtime, hightime do begin
        moments2 = *(momptr(interindex))
        if n_elements(moments2) GT 1 then begin
            temperr = total(abs(smooth(moments2(*, 9), 1, /nan)-smooth(moments2(*, 2), 1, /nan)))
            if not finite(temperr) then begin
                badframe++
            endif else begin
                error += temperr
            endelse
        endif else begin
            badframe++
            print, badframe
        endelse
    endfor

    beta_error_arr(index) = error/(hightime-lowtime-badframe+1) ; averaged error 

    ptr_free, momptr ; free the pointers to avoid memory leaks

endfor

hirexsr_change_crystal_beta, shot, original_beta
hirexsr_invert2tree, shot, 2

plot, dec_arr, beta_error_arr, psym = 2

; now fit a polynomial to this

fit = poly_fit(dec_arr, beta_error_arr, 2) ;  second degree
minima = -1.*fit(1)/(2.*fit(2))
print, 'Approximate minima: ' + string(minima)



;hirexsr_change_crystal_beta, shot, original_beta
;hirexsr_invert_data, shot, 2

;stop
end


;+
;NAME:
; HIREXSR_READ_BETA
;-

pro hirexsr_read_beta, shot, dec = dec

mdsopen, 'spectroscopy', LONG(shot)

if keyword_set(dec) then begin
case dec OF

    1: begin
        rot = mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.INFO.MOD1.MIRROR:ROT')
    end
    2: begin
        rot = mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.INFO.MOD2.MIRROR:ROT')
    end
        3: begin
        rot = mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.INFO.MOD3.MIRROR:ROT')
    end
        4: begin
        rot = mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.INFO.MOD4.MIRROR:ROT')
    end
    ELSE: begin
        rot = mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.INFO.MOD2.MIRROR:ROT')
    end
endcase
endif else begin
     rot = mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.INFO.MOD2.MIRROR:ROT')
 endelse


mdsclose

print, rot(1)

end

;end
