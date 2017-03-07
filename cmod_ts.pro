;ML Reinke - added some conditional statements to prevent loading
;            errors from Thomson reorganization

pro cmod_ts,shot,ts_data,status=status,core=core

tree='electrons'

mdsopen,tree,shot,status=sts

if not sts then begin
  print,'No TS data on shot '+string(shot)
  status=sts
  return
endif

signal='\'+tree+'::top.yag.results.details:timebase' 
time=mdsvalue(signal,stat=stat,/quiet)

;read edge data
signal='\'+tree+'::top.yag_edgets.results:rmid' 
r_edge=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag_edgets.results:rho' 
rho_edge=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag_edgets.results:te' 
Te_edge=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag_edgets.results.te:error' 
dTe_edge=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag_edgets.results:ne' 
ne_edge=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag_edgets.results.ne:error' 
dne_edge=mdsvalue(signal,stat=stat,/quiet)
z_edge=mdsvalue('\top.yag_edgets.data:fiber_z',/quiet)

;read 'old system' core data
signal='\'+tree+'::top.yag.results.details.chord_b:r_mid_t' 
r_b=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag.results.details.chord_b:rho_t' 
rho_b=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag.results.details.chord_b:te_t' 
Te_b=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag.results.details.chord_b:te_err' 
dTe_b=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag.results.details.chord_b:ne_t' 
ne_b=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag.results.details.chord_b:ne_err' 
dne_b=mdsvalue(signal,stat=stat,/quiet)

signal='\'+tree+'::top.yag.results.details.chord_f:r_mid_t' 
r_f=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag.results.details.chord_f:rho_t' 
rho_f=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag.results.details.chord_f:te_t' 
Te_f=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag.results.details.chord_f:te_err' 
dTe_f=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag.results.details.chord_f:ne_t' 
ne_f=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag.results.details.chord_f:ne_err' 
dne_f=mdsvalue(signal,stat=stat,/quiet)

signal='\'+tree+'::top.yag.results.details.chord_g:r_mid_t' 
r_g=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag.results.details.chord_g:rho_t' 
rho_g=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag.results.details.chord_g:te_t' 
Te_g=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag.results.details.chord_g:te_err' 
dTe_g=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag.results.details.chord_g:ne_t' 
ne_g=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag.results.details.chord_g:ne_err' 
dne_g=mdsvalue(signal,stat=stat,/quiet)

;read 'new system' core data
signal='\'+tree+'::top.yag_new.results.profiles:r_mid_t' 
r_core=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag_new.results.profiles:rho_t' 
rho_core=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag_new.results.profiles:te_rz' 
Te_core=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag_new.results.profiles:te_err' 
dTe_core=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag_new.results.profiles:ne_rz' 
ne_core=mdsvalue(signal,stat=stat,/quiet)
signal='\'+tree+'::top.yag_new.results.profiles:ne_err' 
dne_core=mdsvalue(signal,stat=stat,/quiet)
z_new=mdsvalue('\ELECTRONS::TOP.YAG_NEW.RESULTS.PROFILES:Z_SORTED',/quiet)

;combine the data
nt=n_elements(time)
IF NOT keyword_set(core) THEN n_edge=n_elements(rho_edge[0,*]) ELSE n_edge=0
n_core_new=n_elements(rho_core[0,*])
IF shot LT 1070806000 THEN n_core_old=3. ELSE n_core_old = 0.0

r=fltarr(nt,n_edge+n_core_new+n_core_old)
rho=fltarr(nt,n_edge+n_core_new+n_core_old)
Te=fltarr(nt,n_edge+n_core_new+n_core_old)
dTe=fltarr(nt,n_edge+n_core_new+n_core_old)
n_e=fltarr(nt,n_edge+n_core_new+n_core_old)
dn_e=fltarr(nt,n_edge+n_core_new+n_core_old)

IF NOT keyword_set(core) THEN BEGIN
	r[*,0:n_edge-1]=r_edge
        rho[*,0:n_edge-1]=rho_edge
        Te[*,0:n_edge-1]=Te_edge/1000. ;in keV
        dTe[*,0:n_edge-1]=dTe_edge/1000. ;in keV
        n_e[*,0:n_edge-1]=ne_edge
        dn_e[*,0:n_edge-1]=dne_edge
ENDIF
r[*,n_edge:n_edge+n_core_new-1]=r_core
rho[*,n_edge:n_edge+n_core_new-1]=rho_core
Te[*,n_edge:n_edge+n_core_new-1]=Te_core
dTe[*,n_edge:n_edge+n_core_new-1]=dTe_core
n_e[*,n_edge:n_edge+n_core_new-1]=ne_core
dn_e[*,n_edge:n_edge+n_core_new-1]=dne_core

IF shot LT 1070806000 THEN BEGIN
	r[*,n_edge+n_core_new]=r_b
	r[*,n_edge+n_core_new]=r_f
	r[*,n_edge+n_core_new]=r_g

        rho[*,n_edge+n_core_new]=rho_b
	rho[*,n_edge+n_core_new]=rho_f
	rho[*,n_edge+n_core_new]=rho_g

	Te[*,n_edge+n_core_new]=Te_b
	Te[*,n_edge+n_core_new]=Te_f
	Te[*,n_edge+n_core_new]=Te_g

	dTe[*,n_edge+n_core_new]=dTe_b
	dTe[*,n_edge+n_core_new]=dTe_f
	dTe[*,n_edge+n_core_new]=dTe_g

        n_e[*,n_edge+n_core_new]=ne_b
	n_e[*,n_edge+n_core_new]=ne_f
	n_e[*,n_edge+n_core_new]=ne_g

        dn_e[*,n_edge+n_core_new]=dne_b
        dn_e[*,n_edge+n_core_new]=dne_f
        dn_e[*,n_edge+n_core_new]=dne_g
ENDIF

for i=0,nt-1 do begin
  w=sort(rho[i,*])
  if i ne 0 then begin
    if total(rho[i,*])/n_elements(rho[i,*]) lt 0 then w=w_old
  endif
  r_temp=r[i,w]
  rho_temp=rho[i,w]
  Te_temp=Te[i,w]
  dTe_temp=dTe[i,w]
  ne_temp=n_e[i,w]
  dne_temp=dn_e[i,w]
  r[i,*]=r_temp
  rho[i,*]=rho_temp
  Te[i,*]=Te_temp
  dTe[i,*]=dTe_temp
  n_e[i,*]=ne_temp
  dn_e[i,*]=dne_temp
  w_old=w
endfor

;get some geometry data
z_ts=[z_new,z_edge]	;9/23/10 modified so that edge and new were being loaded and compile properly
r_ts=fltarr(n_elements(z_ts))+mdsvalue('.yag.results.param:r')


ts_data={time:time,r:r,rho:rho,Te:Te,dTe:dTe,n_e:n_e,dn_e:dn_e,$
	  r_ts:r_ts,z_ts:z_ts}

status=1
mdsclose
end
