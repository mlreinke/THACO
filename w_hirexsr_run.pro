; program to simplify running hirexsr

; event handler auxiliary widget
; written by YP 10/2010



; takes two help outputs and loads the differences into
; v.ptrs
pro created_diff_array, loadedvarsoriginal, loadedvarsnew, v

; get pointers
ptrvarnames = v.ptrvarnames
ptrvarvalues = v.ptrvarvalues

if n_elements(loadedvarsoriginal) EQ  n_elements(loadedvarsnew) then begin
    if array_equal(loadedvarsoriginal, loadedvarsnew) EQ 1 then begin
        return ; don't waste time if the arrays are the same
    endif
endif

; cycle through new vars
for index =0, n_elements(loadedvarsnew)-1 do begin
    curvar = loadedvarsnew(index)
    where_equal = where(loadedvarsoriginal EQ curvar)
    if n_elements(where_equals) GT 1 OR where_equal(0) EQ -1 then begin
        ; not in original array and should be added
        split_arr = strsplit(curvar, /extract)
        ;0 - name
        ;1 - type
        ;2 - =
        ;3 - value
        
        
        
    endif
endfor


end

pro auxiliary_create_event, event

tag = tag_names(event,/st)
WIDGET_CONTROL, event.top, get_uvalue = v
case tag OF
    'WIDGET_BUTTON':begin
         WIDGET_CONTROL, event.id, get_value=value
         case value OF
             'CLOSE': begin
                 WIDGET_CONTROL, event.top, /destroy
                 ptr_free, v.ptrvarnames
                 ptr_free, v.ptrvarvalues
             end
             'LOAD BETA': begin
                 WIDGET_CONTROL, /hourglass
                 u = v.u
                 activemodules = *u.ptr_activemodules
                 ; get shot
                 WIDGET_CONTROL, u.shottext, get_value = shot_string
                 cshot = LONG(shot_string(0))
                 ; get morder
                WIDGET_CONTROL, u.hemordertext, get_value = strhemorder
                WIDGET_CONTROL, u.hmordertext, get_value = strhmorder
                strhemorder = strhemorder(0); simplify
                strhmorder = strhmorder(0)
                heresult = execute('hemorder = ' + strcompress(strhemorder, /remove_all))
                hresult  = execute('hmorder = ' + strcompress(strhmorder, /remove_all))
                ;get beta
                WIDGET_CONTROL, v.betatext, get_value = beta_string
                beta = FLOAT(beta_string(0))
                
                hirexsr_change_crystal_beta, cshot, beta
                if activemodules(3) EQ 1 then  hirexsr_change_crystal_beta, cshot, beta, /h 
            end

            'CREATE TMAP': begin
                print, 'NOT HOOKED UP YET'
                WIDGET_CONTROL, v.tmintext, get_value = tmin
                WIDGET_CONTROL, v.tmaxtext, get_value = tmax
                WIDGET_CONTROL, v.tmulttext,get_value = tmult
                tmin = FLOAT(tmin(0))
                tmax = FLOAT(tmax(0))
                tmult = FLOAT(tmult(0))
            end

             ELSE: begin
                 ;do nothing
             end
         endcase
             
     end

     'WIDGET_TEXT_CH': begin

         case event.id OF
             v.commandtext: begin
                 WIDGET_CONTROL, v.commandtext, get_value = commandtext
                 WIDGET_CONTROL, v.outputtext, set_value = commandtext, /append
                 commandtext = strcompress(commandtext(0))
                 
                 ; get loaded vars
                 help, names='*', output = loadedvarsoriginal
                 
                 ;execute command
                 result = execute(commandtext)
                 
                 ;get loaded var
                 help, names='*', output=loadedvarsnew
                 
                 ;save/delete changed
                 created_diff_array, loadedvarsoriginal, loadedvarsnew, v
                 
                 if result EQ 1 then begin
                     WIDGET_CONTROL, v.outputtext, set_value = 'COMMAND EXECUTED SUCCESSFULLY', /append
                 endif else begin
                     WIDGET_CONTROL, v.outputtext, set_value = 'COMMAND FAILED', /append
                 endelse
             end
             ELSE: begin
                 ;print, event.id
             End
         endcase


     end

     ELSE: begin
         print, tag
     end
endcase

end

;create auxiliary widget
pro auxiliary_create, u

auxbase = WIDGET_BASE(group_leader = u.base, xsize = 400,/column, title= 'AUXILIARY HIREX COMMANDS')

; create two pointers for arrays of variables loaded by the text base
; interface

varnames = strarr(1)
varvalues = fltarr(1)
ptrvarnames = ptr_new(varnames)
ptrvarvalues = ptr_new(varvalues)

;tmap setup
tmaplabel = WIDGET_LABEL(auxbase, value ='TMAP')
tmapmasterbase = WIDGET_BASE(auxbase, /column, frame = 1)
tminbase = WIDGET_BASE(tmapmasterbase, /row)
tminlabel = WIDGET_LABEL(tminbase, value = 'TMIN: ')
tmintext = WIDGET_TEXT(tminbase, value = '', /editable)
tmaxbase = WIDGET_BASE(tmapmasterbase, /row)
tmaxlabel = WIDGET_LABEL(tmaxbase, value = 'TMAX: ')
tmaxtext = WIDGET_TEXT(tmaxbase, value = '', /editable)
tmultbase = WIDGET_BASE(tmapmasterbase, /row)
tmultlabel = WIDGET_LABEL(tmaxbase, value = 'TMULT: ')
tmulttext = WIDGET_TEXT(tmaxbase, value = '', /editable)
tmapbutton = WIDGET_BUTTON(tmapmasterbase, value = 'CREATE TMAP')

;change crystal beta
betamaxlabel = WIDGET_LABEL(auxbase, value ='CHANGE BETA')
betabase = WIDGET_BASE(auxbase, /row, frame = 1)
betalabel = WIDGET_LABEL(betabase, value = 'BETA: ')
betatext = WIDGET_TEXT(betabase, value = '', /editable)
betabutton = WIDGET_BUTTON(betabase, value='LOAD BETA')

; text interface
textlabel = WIDGET_LABEL(auxbase, value='TEXT INTERFACE')
commandtext = WIDGET_TEXT(auxbase, value = '', /editable)
outputtext = WIDGET_TEXT(auxbase, value = '', ysize = 6, /scroll)
quitbutton = WIDGET_BUTTON(auxbase, value = 'CLOSE')

v = {auxbase:auxbase, commandtext:commandtext, outputtext:outputtext, betabutton:betabutton, betatext:betatext, ptrvarnames:ptrvarnames, ptrvarvalues:ptrvarvalues, tmintext:tmintext, tmaxtext:tmaxtext, tmulttext:tmulttext,  u:u}
WIDGET_CONTROL, auxbase, set_uvalue = v
WIDGET_CONTROL, auxbase, /realize
xmanager, 'auxiliary_create', auxbase

end

; primary event handler
pro w_hirexsr_run_event, event

tag = tag_names(event,/st)
WIDGET_CONTROL, event.top, get_uvalue = u
WIDGET_CONTROL, u.shottext, get_value = shot_string
cshot = LONG(shot_string)
activemodules = *u.ptr_activemodules
activelines = *u.ptr_activelines
;print, 'Active lines', activelines
;print, 'Active Modules', activemodules

case tag OF
    'WIDGET_BUTTON':begin

        WIDGET_CONTROL, event.id, get_value=value
        case value OF
            'QUIT': begin
                ptr_free, u.ptr_activemodules
                WIDGET_CONTROL, event.top, /destroy
            end
            'STOP':begin
                stop ; for debugging only
            end

            'CURRENT':begin
                update_shot = mdsvalue('current_shot("cmod")')
                WIDGET_CONTROL, u.shottext, set_value = STRCOMPRESS(STRING(update_shot), /remove_all)
                cshot = update_shot
            end

            'RAW':begin
                hirex_look, shot = cshot
            end

            'CALIB':begin
                w_hirexsr_calib, shot = cshot
            end
            'DET_ALIGN':begin
                w_hirexsr_det_align, shot = cshot
            end
            'HE_MOM':begin
                w_hirexsr_he_moments, shot = cshot
            end

            'H_MOM': begin
                 w_hirexsr_h_moments, shot = cshot
            end
            
            'PROFILES': begin
                w_hirexsr_profiles, shot = cshot
            end

            'COMPARE': begin
                w_hirexsr_compare, shot = cshot
            end

            ; module buttons
            '1': begin
                buttonval = WIDGET_INFO(u.mod1button, /button_set)
                activemodules(0) = buttonval
                *(u.ptr_activemodules) = activemodules
            end
            '2': begin 
                buttonval = WIDGET_INFO(u.mod2button, /button_set)
                activemodules(1) = buttonval 
                *(u.ptr_activemodules) = activemodules
            end
            '3':begin
                buttonval = WIDGET_INFO(u.mod3button, /button_set)
                activemodules(2) = buttonval 
                *(u.ptr_activemodules) = activemodules
            end
            '4':begin
                buttonval = WIDGET_INFO(u.mod4button, /button_set)
                activemodules(3) = buttonval  
                *(u.ptr_activemodules) = activemodules
            end
            
            ; line buttons
            'w': begin
                buttonval = WIDGET_INFO(u.line1button, /button_set)
                activelines(0) = buttonval
                *(u.ptr_activelines) = activelines
            end
            'y': begin
                buttonval = WIDGET_INFO(u.line2button, /button_set)
                activelines(1) = buttonval   
                *(u.ptr_activelines) = activelines
            end
            'z': begin
                buttonval = WIDGET_INFO(u.line3button, /button_set)
                activelines(2) = buttonval  
                *(u.ptr_activelines) = activelines
            end
            'lya1': begin
                buttonval = WIDGET_INFO(u.line4button, /button_set)
                activelines(3) = buttonval                
                *(u.ptr_activelines) = activelines
            end
            
            'UPDATE TREE':begin
                hirexsr_update_tree, cshot
            end

            'COPY INFO': begin
                WIDGET_CONTROL, u.bshottext, get_value = bshot_string
                if bshot_string NE '' then begin
                    bshot = LONG(bshot_string)
                    if activemodules(0) EQ 1 then hirexsr_copy_info, bshot, cshot, module = 1 
                    if activemodules(1) EQ 1 then hirexsr_copy_info, bshot, cshot, module = 2
                    if activemodules(2) EQ 1 then hirexsr_copy_info, bshot, cshot, module = 3
                    if activemodules(3) EQ 1 then hirexsr_copy_info, bshot, cshot, module = 4
                endif else begin
                    ; write warning message here
                endelse
            end
            
            'WRITE WAVELENGTH': begin
                hirexsr_write_wavelengths, cshot
            end

            'COPY CALIB': begin
                WIDGET_CONTROL, u.bshottext, get_value = bshot_string
                if bshot_string NE '' then begin
                    bshot = LONG(bshot_string)
                    if activemodules(0) EQ 1 then hirexsr_copy_calib, bshot, cshot, module = 1 
                    if activemodules(1) EQ 1 then hirexsr_copy_calib, bshot, cshot, module = 2
                    if activemodules(2) EQ 1 then hirexsr_copy_calib, bshot, cshot, module = 3
                    if activemodules(3) EQ 1 then hirexsr_copy_calib, bshot, cshot, module = 4
                endif else begin
                    ; write warning message here
                endelse            
            end


            'WRITE WHITE': begin
                if activemodules(0) EQ 1 then hirexsr_write_white, cshot, 1
                if activemodules(1) EQ 1 then hirexsr_write_white, cshot, 2
                if activemodules(2) EQ 1 then hirexsr_write_white, cshot, 3
                if activemodules(3) EQ 1 then hirexsr_write_white, cshot, 4
            end

             'WRITE TRANS': begin
                if activemodules(0) EQ 1 then hirexsr_write_trans, cshot, 1
                if activemodules(1) EQ 1 then hirexsr_write_trans, cshot, 2
                if activemodules(2) EQ 1 then hirexsr_write_trans, cshot, 3
                if activemodules(3) EQ 1 then hirexsr_write_trans, cshot, 4
            end  

            'WRITE MORDER': begin
                WIDGET_CONTROL, u.hemordertext, get_value = strhemorder
                WIDGET_CONTROL, u.hmordertext, get_value = strhmorder
                strhemorder = strhemorder(0); simplify
                strhmorder = strhmorder(0)
                ; nice way of getting the results
                heresult = execute('hemorder = ' + strcompress(strhemorder, /remove_all))
                hresult  = execute('hmorder = ' + strcompress(strhmorder, /remove_all))
                ; load data
                if heresult EQ 1 then hirexsr_write_morder, cshot, morder = hemorder else print, 'ONE DOES NOT SIMPLY WALK INTO MORDER'
                if hresult EQ 1 then hirexsr_write_morder, cshot, morder = hmorder, /h else print, 'ONE DOES NOT SIMPLY WALK INTO MORDER'
            end
            
            'COPY BINNING':begin
                WIDGET_CONTROL, u.bshottext, get_value = bshot_string
                if bshot_string NE '' then begin
                    chmapval = WIDGET_INFO(u.chmapbutton, /button_set)
                    tmapval = WIDGET_INFO(u.tmapbutton, /button_set)
                    goodval = WIDGET_INFO(u.goodbutton, /button_set)
                    hbinval = WIDGET_INFO(u.hbinbutton, /button_set)
                    
                                ; set up string to execute
                    init_string = 'hirexsr_copy_binning,'+STRCOMPRESS(bshot_string, /remove_all)+', '+ STRCOMPRESS(STRING(cshot), /remove_all)
                    if chmapval EQ 1 then init_string  = init_string + ', /chmap'
                    if tmapval EQ 1 then init_string = init_string +', /tmap'
                    if goodval EQ 1 then init_string  = init_string +', /good'
                    if hbinval EQ 1 then init_string = init_string +', /h'
                    print, init_string
                    init_string = init_string(0)
                    result = execute(init_string)
                endif
                
            end
            

            'AVESPEC': begin
                WIDGET_CONTROL, /hourglass
                hirexsr_avespec2tree, cshot
            end

            'FITSPEC': begin
               nback = 0
               buttonval = WIDGET_INFO(u.nbackbutton, /button_set)
               if buttonval EQ 1 then nback = 1
                WIDGET_CONTROL, /hourglass
                if activelines(0) EQ 1 then hirexsr_fitspec2tree, cshot, 0, nback = nback
                if activelines(1) EQ 1 then hirexsr_fitspec2tree, cshot, 1, nback = nback
                if activelines(2) EQ 1 then hirexsr_fitspec2tree, cshot, 2, nback = nback
                if activelines(3) EQ 1 then hirexsr_fitspec2tree, cshot, 3, nback = nback
            end

            'INVERT':begin
                
                sine = 0
                buttonval = WIDGET_INFO(u.sinebutton, /button_set)
                if buttonval EQ 1 then sine = 1
                WIDGET_CONTROL, /hourglass
                if activelines(0) EQ 1 then hirexsr_invert2tree, cshot, 0, sine = sine
                if activelines(1) EQ 1 then hirexsr_invert2tree, cshot, 1, sine = sine
                if activelines(2) EQ 1 then hirexsr_invert2tree, cshot, 2, sine = sine
                if activelines(3) EQ 1 then hirexsr_invert2tree, cshot, 3, sine = sine
            end


            'AUXILIARY': begin
                auxiliary_create, u                
            end

            'CHMAP':begin
                ;nothing
            end
            
            'TMAP': begin
                ; nothing
            end
            
            'GOOD': begin
                ;nothing
            end

            'H BIN': begin
                ; nothing
            end

            'SIN':begin
                ; do nothing
            end

     
            ELSE: begin
                print, 'UNHANDLED BUTTON'
            end

        endcase
    end

    'WIDGET_TEXT_CH': begin
        ; nothing here yet
    end

    ELSE: begin
        print, 'UNHANDLED EXCEPTION: ' + tag
    END
endcase

end


pro w_hirexsr_run, debug = debug

user = REFORM(logname()) ; get user name. Why not? 
user = user(0)           ; recast

; user check; really if you are not a spectroscopy person you should
; not be using this
if user NE 'ypodpaly' and user NE 'mlreinke' and user NE 'rice' and user NE 'ldelgado' and user NE 'cgao' and user NE 'npablant' then begin
    print, 'NOT A MEMBER OF THE SPECTROSCOPY GROUP'
    return
endif


if keyword_set(debug) then begin
    shotnumber = 1100325022
endif else begin
    cshot = mdsvalue('current_shot("cmod")')
    shotnumber = REFORM(cshot(0))
endelse

; initializing
activemodules = [1, 1, 1, 1] ; all modules are on
ptr_activemodules = ptr_new(activemodules)
activelines = [1, 0, 1, 1]
ptr_activelines = ptr_new(activelines)

; ****************widget building***************

base = WIDGET_BASE(/column, xsize = 250, title = 'HiReX Sr Run'); top level base

; identify user
userbase  = WIDGET_BASE(base, /row, frame = 1)
userlabel = WIDGET_LABEL(userbase, value = 'USER: '+ strcompress(string(user), /remove_all)) 

; identify shot
shotbase  = WIDGET_BASE(base, /row, frame = 1)
shotlabel = WIDGET_LABEL(shotbase, value = 'SHOT: ')
shottext  = WIDGET_TEXT(shotbase, value = strcompress(string(shotnumber), /remove_all), /editable, xsize = 18)
currbutton = WIDGET_BUTTON(shotbase, value = 'CURRENT')

; identify locked mode shot shot
bshotbase  = WIDGET_BASE(base, /row, frame = 1)
bshotlabel = WIDGET_LABEL(bshotbase, value = 'BASE SHOT: ')
bshottext  = WIDGET_TEXT(bshotbase, value = '', /editable, xsize = 13)

; identify active modules
modlabel = WIDGET_LABEL(base, value = 'ACTIVE MODULES: ')
modbase = WIDGET_BASE(base, /row, /nonexclusive, frame = 1)
mod1button = WIDGET_BUTTON(modbase, value = '1')
mod2button = WIDGET_BUTTON(modbase, value = '2')
mod3button = WIDGET_BUTTON(modbase, value = '3')
mod4button = WIDGET_BUTTON(modbase, value = '4')

; identify active lines
linelabel = WIDGET_LABEL(base, value = 'ACTIVE LINES: ')
linebase = WIDGET_BASE(base, /row, /nonexclusive, frame = 1)
line1button = WIDGET_BUTTON(linebase, value = 'w')
line2button = WIDGET_BUTTON(linebase, value = 'y')
line3button = WIDGET_BUTTON(linebase, value = 'z')
line4button = WIDGET_BUTTON(linebase, value = 'lya1')

; individual steps button
updatebutton = WIDGET_BUTTON(base, value='UPDATE TREE')
infobutton = WIDGET_BUTTON(base, value='COPY INFO')
wavelengthsbutton = WIDGET_BUTTON(base, value='WRITE WAVELENGTH')
calibbutton = WIDGET_BUTTON(base, value = 'COPY CALIB')
whitebutton = WIDGET_BUTTON(base, value = 'WRITE WHITE')
transbutton = WIDGET_BUTTON(base, value = 'WRITE TRANS')
morderbutton = WIDGET_BUTTON(base, value= 'WRITE MORDER')

; morder alignment
hemorderbase = WIDGET_BASE(base, /row)
hmorderbase = WIDGET_BASE(base, /row)
hemorderlabel = WIDGET_LABEL(hemorderbase, value = 'He MORDER')
hemordertext = WIDGET_TEXT(hemorderbase, value = '[1,2,3]', /editable)
hmorderlabel = WIDGET_LABEL(hmorderbase, value = 'H  MORDER')
hmordertext = WIDGET_TEXT(hmorderbase, value = '[4]', /editable)

; binning step
binbutton = WIDGET_BUTTON(base, value = 'COPY BINNING')

; binning set-up
binlabel = WIDGET_LABEL(base, value = 'BINNING COPY: ')
binbase = WIDGET_BASE(base, /row, /nonexclusive, frame = 1)
chmapbutton = WIDGET_BUTTON(binbase, value = 'CHMAP')
tmapbutton = WIDGET_BUTTON(binbase, value = 'TMAP')
goodbutton = WIDGET_BUTTON(binbase, value = 'GOOD')
hbinbutton = WIDGET_BUTTON(binbase, value = 'H BIN')

; create data steps

avespecbutton = WIDGET_BUTTON(base, value = 'AVESPEC')
fitbase = WIDGET_BASE(base, /row)
fitspecbutton = WIDGET_BUTTON(fitbase, value = 'FITSPEC')
nbackbase = WIDGET_BASE(fitbase, frame = 1, /nonexclusive, xsize = 50)
nbackbutton=WIDGET_BUTTON(nbackbase, value = 'NBACK')
invertbase = WIDGET_BASE(base, /row)
invertbutton = WIDGET_BUTTON(invertbase, value = 'INVERT', xsize = 150)
sinebase = WIDGET_BASE(invertbase, frame = 1, /nonexclusive, xsize = 50)
sinebutton = WIDGET_BUTTON(sinebase, value = 'SIN')

; four widget base
widgetlabel = WIDGET_LABEL(base, value = 'WIDGETS:') 
widgbase = WIDGET_BASE(base, xsize = 250,/column, frame = 1, /align_center)
rawbutton = WIDGET_BUTTON(widgbase, value = 'RAW', /align_center, xsize = 200)
calibbutton = WIDGET_BUTTON(widgbase, value = 'CALIB', /align_center, xsize = 200)
detalignbutton = WIDGET_BUTTON(widgbase, value='DET_ALIGN', /align_center, xsize = 200)
hemombutton = WIDGET_BUTTON(widgbase, value='HE_MOM', /align_center, xsize = 200)
hmombutton = WIDGET_BUTTON(widgbase, value='H_MOM', /align_center, xsize = 200)
profbutton = WIDGET_BUTTON(widgbase, value='PROFILES', /align_center, xsize = 200)
compbutton = WIDGET_BUTTON(widgbase, value='COMPARE', /align_center, xsize =200)

;key buttons
controllabel = WIDGET_LABEL(base, value ='CONTROLS: ')
auxbutton = WIDGET_BUTTON(base, value = 'AUXILIARY')
if keyword_set(debug) then stopbutton = WIDGET_BUTTON(base, value = 'STOP')
quitbutton = WIDGET_BUTTON(base, value = 'QUIT')

; ******************** widget building end ***********************

u = {shottext:shottext, bshottext:bshottext, calibbutton:calibbutton, detalignbutton:detalignbutton, hemombutton:hemombutton, profbutton:profbutton, quitbutton:quitbutton, mod1button:mod1button, mod2button:mod2button, mod3button:mod3button, mod4button:mod4button, line1button:line1button, line2button:line2button, line3button:line3button, line4button:line4button, ptr_activemodules:ptr_activemodules, ptr_activelines:ptr_activelines, hemordertext:hemordertext, hmordertext:hmordertext, auxbutton:auxbutton, base:base, chmapbutton:chmapbutton, tmapbutton:tmapbutton, goodbutton:goodbutton, hbinbutton:hbinbutton, sinebutton:sinebutton, nbackbutton:nbackbutton}



WIDGET_CONTROL, base, set_uvalue = u
WIDGET_CONTROL, mod1button, /set_button
WIDGET_CONTROL, mod2button, /set_button
WIDGET_CONTROL, mod3button, /set_button
WIDGET_CONTROL, mod4button, /set_button
WIDGET_CONTROL, line1button, /set_button
;WIDGET_CONTROL, line2button, /set_button
WIDGET_CONTROL, line3button, /set_button
WIDGET_CONTROL, line4button, /set_button
WIDGET_CONTROL, chmapbutton, /set_button
WIDGET_CONTROL, tmapbutton, /set_button
WIDGET_CONTROL, goodbutton, /set_button

WIDGET_CONTROL, base, /realize
xmanager, 'w_hirexsr_run', base


end
