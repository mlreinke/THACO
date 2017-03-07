; code to take a hirex_sr inversion and normalize it with Jr data
;ypodpaly@mit.edu for edit requests
; updated  07/2011 to implement new data from HIREXSR

; this is the under construction normalizer
;@/usr/local/cmod/codes/spectroscopy/hirex_sr/hirex_sr_inversion.pro 
;procedure that prints string only if quiet is not set
pro quietprint, string, quiet = quiet
  if  not keyword_set(quiet) then begin;quiet NE 0 AND quiet NE 0./0. then begin
    print, string
   endif
end


function jrnormalizer,filename = filename, struct = struct, shot = shot, onlyw=onlyw, quiet=quiet, trange = trange, smooth = smooth, noinv = noinv
; filename - hirex inversion file
; struct - hirex inverted struct
; shot - provide a shot and it will invert it for you
; onlyw - only calibrate w line for debugging, not currently operational
; quiet, do not print warnings
; trange - range of calibration
; smooth - box car smooth on data value, only used to calculate the offset
; noinv - do not include the original inversion struct in the return

SIZE_STRUCT_VALUE = 8    ; size returns 8 if it is given a struct
jrtoreal = 10^(-5.)      ;  John stores things in cm/s; this converts it to to km/s
ERROR_CODE = -1          ; return on failure

;redefine quiet for easier if then statements
if keyword_set(quiet) then begin
        quiet = 1
    endif else begin
        quiet = 0
    endelse

if (keyword_set(filename) + keyword_set(struct) + keyword_set(shot)) NE 1 then begin
    quietprint, 'EITHER FILENAME OR STRUCT OPR SHOT MUST BE SET', quiet = quiet
    return, ERROR_CODE
endif

;create the smoothing var
if keyword_set(smooth) then begin
    smoothing = smooth
endif else begin
    smoothing = 0 ; do not smooth
end

if keyword_set(filename) then begin
    if quiet EQ 1 then begin
       restore, filename
    endif else begin
       restore, filename, /verb
    endelse
    inversionstr = inversion ; inversion is the general name that we use
endif 

if keyword_set(struct) then begin ; get the struct
    inversionstr = struct
    quietprint, 'STRUCT FOUND', quiet = quiet
endif 


if keyword_set(shot) then begin
    if quiet EQ 0 then begin
        inversionstr = hirexsr_get_profile(shot)
    endif else begin
        inversionstr = hirexsr_get_profile(shot, /quiet)
    endelse
endif



;restored file must have struct called inversion
if size(inversionstr, /type) EQ SIZE_STRUCT_VALUE then begin
	shot = LONG(inversionstr.shot) ; get shot number

	; get junior data
	mdsopen, 'spectroscopy', shot 
        jrrotation = mdsvalue('\SPECTROSCOPY::TOP.HIREXJR.ANALYSIS:ROT_VELOCITY', quiet = quiet)*jrtoreal
        jrtime = mdsvalue('\SPECTROSCOPY::TOP.HIREXJR.ANALYSIS:TIME', quiet = quiet)
	mdsclose

	;until further evaluation, assume that rotation = on axis rotation
	;this can be repaired by using the efit reconstruction

	if keyword_set(trange) then begin ; if trange specified do the renormalization only there
		time_range = trange
	endif else begin
		quietprint, 'UNSPECIFIED TIME RANGE: Taking [0.5, 1.5] s range', quiet = quiet
		time_range = [.5, 1.5]	; generally good time range
        endelse

        
       ; now cycle through w, z, lya1, lya2,
       ; m32 and renormalize to jr on axis

        
        ;w line
        ;test_w = size(inversionstr.w, /type) ; check if struck
        ;if test_w EQ SIZE_STRUCT_VALUE then begin 
        ;    quietprint, 'FOUND W LINE', quiet = quiet
        ;    
        ;    times = inversionstr.w.times
        ;    ;find low and high time
        ;    lowind = where(times GT time_range(0))
        ;    lowind = lowind(0)
        ;    highind = where(times GT time_range(1))
        ;    highind = highind(0)

        ;    ; take the whole range
        ;    if highind EQ -1 then begin
        ;        highind = n_elements(times)-1 
        ;    endif

            
        ;    wfix = fltarr(highind-lowind+1)
        ;    omega = smooth(inversionstr.w.omega(0, *), smoothing)
        ;    count = 0;
        ;    for index = lowind, highind do begin
        ;        jrind = where(jrtime GT times(index))
        ;        jrind = jrind(0)

         ;       ; select the two indices around 
         ;       jrindlow = jrind-1
         ;       jrindhigh = jrind
         ;       deltv = jrrotation(jrindhigh)-jrrotation(jrindlow)
                ; project the velocity
         ;       velocity_project = jrrotation(jrindlow) + deltv*(jrtime(jrindhigh)-times(index))/(jrtime(jrindhigh)-jrtime(jrindlow)) 
                
                ; now find how to offset the w line
         ;       omegaaxis = omega(index);inversionstr.w.omega(0, index)           
         ;       vaxis = omegaaxis*2*!Pi*inversionstr.w.r_major(0, index)
         ;       wfix(count) = velocity_project-vaxis
         ;       count++
         ;   endfor
         ;   avewoffset = mean(wfix)

         ;   ; now create a velocity array for the corrected data
         ;   wvel = inversionstr.w.omega
         ;   for index = 0, n_elements(wvel(*, 0)) -1 do begin
         ;       wvel(index, *) = 2*!Pi*inversionstr.w.omega(index,*)*inversionstr.w.r_major(index, *)+avewoffset
         ;   endfor

         ;   quietprint, 'CONTINUING', quiet = quiet
        ;endif else begin
        ;    quietprint, 'NO W LINE', quiet = quiet
        ;    quietprint, 'CONTINUING', quiet = quiet
        ;    avewoffset=0./0.
        ;    wvel=0./0.
        ;endelse	


     ;z line
        ;test_z = size(inversionstr, /type) ; check if struck
        ;if test_z EQ SIZE_STRUCT_VALUE then begin 
            ;quietprint, 'FOUND Z LINE', quiet = quiet
            
            times = inversionstr.time
            ;find low and high time
            lowind = where(times GT time_range(0))
            lowind = lowind(0)
            highind = where(times GT time_range(1))
            highind = highind(0)

                        ; take the whole range
            if highind EQ -1 then begin
                highind = n_elements(times)-1 
            endif
            zfix = fltarr(highind-lowind+1)
            omega = smooth(inversionstr.omg(0, *), smoothing)
            count = 0;
            for index = lowind, highind do begin
                jrind = where(jrtime GT times(index))
                jrind = jrind(0)

                ; select the two indices around 
                jrindlow = jrind-1
                jrindhigh = jrind
                deltv = jrrotation(jrindhigh)-jrrotation(jrindlow)
                ; project the velocity
                velocity_project = jrrotation(jrindlow) + deltv*(jrtime(jrindhigh)-times(index))/(jrtime(jrindhigh)-jrtime(jrindlow)) 
                
                ; now find how to offset the w line
                omegaaxis = omega(index);inversionstr.z.omega(0, index)
                vaxis = omegaaxis*2*!Pi*inversionstr.r_maj(0, index)
                zfix(count) = velocity_project-vaxis
                count++
            endfor
            aveoffset = mean(zfix)
            aveomegaoffset = aveoffset/(2. *!Pi*mean(inversionstr.r_maj(0, lowind:highind)))

            ; now create a velocity array for the corrected data
            corrvel = inversionstr.omg
            for index = 0, n_elements(corrvel(*, 0)) -1 do begin
                corrvel(index, *) = 2*!Pi*inversionstr.omg(index,*)*inversionstr.r_maj(index, *)+aveoffset
            endfor

            quietprint, 'CONTINUING',quiet = quiet
        ;endif else begin
        ;    quietprint, 'NO Z LINE', quiet = quiet
        ;    quietprint, 'CONTINUING', quiet = quiet
        ;    avezoffset=0./0.
        ;    zvel=0./0.
        ;endelse	

       ;test_lya1 = size(inversionstr.lya1, /type) ; check if struck
       ; if test_lya1 EQ SIZE_STRUCT_VALUE then begin 
       ;     quietprint, 'FOUND LYA1 LINE', quiet = quiet
       ;     
       ;     times = inversionstr.lya1.times
       ;     ;find low and high time
       ;     lowind = where(times GT time_range(0))
       ;     lowind = lowind(0)
       ;     highind = where(times GT time_range(1))
       ;     highind = highind(0)

       ;     ; take the whole range
       ;     if highind EQ -1 then begin
       ;         highind = n_elements(times)-1 
       ;     endif

       ;     lya1fix = fltarr(highind-lowind+1)
       ;     omega = smooth(inversionstr.lya1.omega(0, *), smoothing)
       ;     count = 0;
       ;     for index = lowind, highind do begin
       ;         jrind = where(jrtime GT times(index))
       ;         jrind = jrind(0)

       ;         ; select the two indices around 
       ;         jrindlow = jrind-1
       ;         jrindhigh = jrind
       ;         deltv = jrrotation(jrindhigh)-jrrotation(jrindlow)
                ; project the velocity
       ;         velocity_project = jrrotation(jrindlow) + deltv*(jrtime(jrindhigh)-times(index))/(jrtime(jrindhigh)-jrtime(jrindlow)) 
                
       ;         ; now find how to offset the w line
       ;         omegaaxis =omega(index);inversionstr.lya1.omega(0, index)
       ;         vaxis = omegaaxis*2*!Pi*inversionstr.w.r_major(0, index)
       ;         lya1fix(count) = velocity_project-vaxis
       ;         count++
       ;     endfor
       ;     avelya1offset = mean(lya1fix)

       ;    ; now create a velocity array for the corrected data
       ;    lya1vel = inversionstr.lya1.omega
       ;     for index = 0, n_elements(lya1vel(*, 0)) -1 do begin
       ;         lya1vel(index, *) = 2*!Pi*inversionstr.lya1.omega(index,*)*inversionstr.lya1.r_major(index, *)+avelya1offset
       ;     endfor
       ;     quietprint, 'CONTINUING', quiet = quiet
       ; endif else begin
       ;     quietprint, 'NO LYA1 LINE', quiet = quiet
       ;     quietprint, 'CONTINUING', quiet = quiet
       ;     avelya1offset=0./0.
       ;     lya1vel=0./0.
       ; endelse	

        ; lyman alpha 2
       ;test_lya2 = size(inversionstr.lya2, /type) ; check if struck
       ; if test_lya2 EQ SIZE_STRUCT_VALUE then begin 
        ;    quietprint, 'FOUND LYA2 LINE', quiet = quiet
        ;    
        ;    times = inversionstr.lya2.times
        ;    ;find low and high time
        ;    lowind = where(times GT time_range(0))
        ;    lowind = lowind(0)
        ;    highind = where(times GT time_range(1))
        ;    highind = highind(0)
        ;    
        ;    ; take the whole range
        ;    if highind EQ -1 then begin
        ;        highind = n_elements(times)-1 
        ;    endif

        ;    lya2fix = fltarr(highind-lowind+1)
        ;    omega = smooth(inversionstr.lya2.omega(0, *), smoothing)
        ;    count = 0;
        ;    for index = lowind, highind do begin
        ;        jrind = where(jrtime GT times(index))
        ;        jrind = jrind(0);

        ;        ; select the two indices around 
        ;        jrindlow = jrind-1
        ;        jrindhigh = jrind
        ;        deltv = jrrotation(jrindhigh)-jrrotation(jrindlow)
        ;        ; project the velocity
        ;        velocity_project = jrrotation(jrindlow) + deltv*(jrtime(jrindhigh)-times(index))/(jrtime(jrindhigh)-jrtime(jrindlow)) 
                
        ;        ; now find how to offset the w line
        ;        omegaaxis = omega(index);inversionstr.lya2.omega(0, index)
        ;        vaxis = omegaaxis*2*!Pi*inversionstr.lya2.r_major(0, index)
        ;        lya2fix(count) = velocity_project-vaxis
        ;        count++
        ;    endfor
        ;    avelya2offset = mean(lya2fix)

        ;    ; now create a velocity array for the corrected data
        ;    lya2vel = inversionstr.lya2.omega
        ;    for index = 0, n_elements(lya2vel(*, 0)) -1 do begin
        ;        lya2vel(index, *) = 2*!Pi*inversionstr.lya2.omega(index,*)*inversionstr.lya2.r_major(index, *)+avelya2offset
        ;    endfor

        ;    quietprint, 'CONTINUING', quiet = quiet
        ;endif else begin

        ;    quietprint, 'NO LYA2 LINE', quiet = quiet
        ;    quietprint, 'CONTINUING', quiet = quiet
        ;    avelya2offset=0./0.
        ;    lya2vel=0./0.
        ;endelse	

        ;moly calib
       ;test_m32 = size(inversionstr.m32, /type) ; check if struck
       ; if test_m32 EQ SIZE_STRUCT_VALUE then begin 
       ;     quietprint, 'FOUND M32 LINE', quiet = quiet
            
       ;     times = inversionstr.m32.times
       ;     ;find low and high time
       ;     lowind = where(times GT time_range(0))
       ;     lowind = lowind(0)
       ;     highind = where(times GT time_range(1))
       ;     highind = highind(0)
            
       ;                         ; take the whole range
       ;     if highind EQ -1 then begin
       ;         highind = n_elements(times)-1 
       ;     endif


      ;      m32fix = fltarr(highind-lowind+1)
      ;      omega = smooth(inversionstr.m32.omega(0, *), smoothing)
      ;      count = 0;
      ;      for index = lowind, highind do begin
      ;          jrind = where(jrtime GT times(index))
      ;          jrind = jrind(0)

      ;          ; select the two indices around 
      ;          jrindlow = jrind-1
      ;          jrindhigh = jrind
      ;          deltv = jrrotation(jrindhigh)-jrrotation(jrindlow)
                ; project the velocity
      ;          velocity_project = jrrotation(jrindlow) + deltv*(jrtime(jrindhigh)-times(index))/(jrtime(jrindhigh)-jrtime(jrindlow)) 
                
      ;          ; now find how to offset the w line
      ;          omegaaxis = omega(index);inversionstr.m32.omega(0, index)
      ;          vaxis = omegaaxis*2*!Pi*inversionstr.m32.r_major(0, index)
      ;          m32fix(count) = velocity_project-vaxis
      ;          count++
      ;      endfor
      ;      avem32offset = mean(m32fix)

      ;      ; now create a velocity array for the corrected data
      ;      m32vel = inversionstr.m32.omega
      ;      for index = 0, n_elements(m32vel(*, 0)) -1 do begin
      ;          m32vel(index, *) = 2*!Pi*inversionstr.m32.omega(index,*)*inversionstr.m32.r_major(index, *)+avem32offset
      ;      endfor

      ;      quietprint, 'CONTINUING', quiet = quiet
      ;  endif else begin
      ;      quietprint, 'NO M32 LINE', quiet = quiet
      ;      quietprint, 'CONTINUING', quiet = quiet
      ;      avem32offset=0./0.
       ;     m32vel=0./0.
      ;  endelse	

    endif else begin
	quietprint, 'DATA OF WRONG FORM', quiet = quiet
	quietprint, 'ABORTING', quiet = quiet
        return, ERROR_CODE
    endelse

    if keyword_set(noinv) then begin
        srvelstruct = {aveoffset:aveoffset, corrvel:corrvel, srtime:inversionstr.time, srrmaj:inversionstr.r_maj, jrrotation:jrrotation, jrtime:jrtime, shot:shot}
    endif else begin
        srvelstruct = {aveoffset:aveoffset, corrvel:corrvel, srtime:inversionstr.time, srrmaj:inversionstr.r_maj, jrrotation:jrrotation, jrtime:jrtime, shot:shot, inversionstr:inversionstr}
     endelse 
    return, srvelstruct
       
end

