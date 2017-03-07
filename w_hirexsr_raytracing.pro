;     this version adds writting-to-tree button comparing with version 5.

;     .compile /home/mlreinke/idl/general/mlr_functions
;     .compile /home/mlreinke/idl/genie/genie_line
;     .compile /home/mlreinke/idl/genie/genpos.pro
;     .compile /home/mlreinke/idl/genie/genpos_sphere
;     .compile /home/mlreinke/idl/hirexsr/hirexsr_load_data
;     .compile /home/mlreinke/idl/genie/genpos_sphere_scripts

; modified on Sep 2nd, 2010 by Chi Gao

FUNCTION tran_position,X,X0,angle=angle

;   x is [3,3] in tokamak fram
;   x0 is the origin of crystal in tokamak frame
;   angle = phi - beta
;   return [3,3] in crystal frame
    y=fltarr(3,3)
   
    y(0)=(x(0)-x0(0))*(-cos(angle))+(x(1)-x0(1))*(sin(angle))
    y(1)=(x(0)-x0(0))*(-sin(angle))+(x(1)-x0(1))*(-cos(angle))
    y(2)=x(2)-x0(2)
    
    y(3)=(x(3)-x0(0))*(-cos(angle))+(x(4)-x0(1))*(sin(angle))
    y(4)=(x(3)-x0(0))*(-sin(angle))+(x(4)-x0(1))*(-cos(angle))
    y(5)=x(5)-x0(2)

    y(6)=(x(6)-x0(0))*(-cos(angle))+(x(7)-x0(1))*(sin(angle))
    y(7)=(x(6)-x0(0))*(-sin(angle))+(x(7)-x0(1))*(-cos(angle))
    y(8)=x(8)-x0(2)

    return,Y
END

PRO w_helike_ar_vignetting,u
  
    shot=u.info.shot
    module=u.info.module
    help,shot,module
 
    t1=u.info.t1
    t2=u.info.t2
    widget_control,u.id.xi_slider,get_value=xi
    widget_control,u.id.zeta_slider,get_value=zeta   
    det=module

;setup detector points
    info=u.load_info
    x0=info.det.x0
    x1=info.det.x1
    x2=info.det.x2
    det_pos=[xi[0],zeta[0]]*172.0e-6
    d_pt=genpos_det2xyz(x0,x1,x2,det_pos)
	
;setup mirror points
    ny=u.info.he_ny
    nz=u.info.he_nz
    ;ny=50
    ;nz=50
    m_pts=genpos_grid_sphere(info,ny=ny,nz=nz)
    n_mir=ny*nz

    trans_info=fltarr(12,n_mir)
    trans_info(11,*)=1
	
;define reflected points
    r_pts=fltarr(3,n_mir)
    FOR i=0,n_mir-1 DO BEGIN
        xs=m_pts[*,i]
	x1=d_pt
	n_vec=[info.m.rad,0,0]-xs
	n_mag=sqrt(total(n_vec*n_vec))
	n=n_vec/n_mag
	r_pts[*,i]=genpos_spherical_refl(x1,xs,n)
    ENDFOR
;*****************************************************************

;*****************************************************************
;b-port points
    bpts=u.info.position.he.bport 
;	bpts=[[1648.44,932.16,6.7],$
;	      [1664.9,890.4,6.74],$
;	      [1648.36,932.3,190.85]]
    bpts/=1.0e3
    a=bpts
    dx=sqrt((a[0,0]-a[0,1])^2+(a[1,0]-a[1,1])^2+(a[2,0]-a[2,1])^2)
    dy=sqrt((a[0,2]-a[0,1])^2+(a[1,2]-a[1,1])^2+(a[2,2]-a[2,1])^2)

    wid=89.77/1.0e3
    ht=279.40/1.0e3
		
;make puncture plot on the b-port plane
    widget_control,u.id.draw1,get_value=draw_win
    wset,draw_win
    plot, [0],[0],xr=[-0.25,0.25],yr=[-0.25,0.25],/xsty,/ysty,position=[0.1,0.1,0.9,0.9];,tit='B-Port Flange'
    uvloc=fltarr(2,n_mir)
    area_through=0.
    area_total=0.
    FOR i=0,n_mir-1 DO BEGIN
        l=[[m_pts[*,i]],$
           [r_pts[*,i]]]
        out=line_plane_int(l,bpts)
        
        trans_info(0,i)=a[0,0]+(a[0,1]-a[0,0])*out[1]+(a[0,2]-a[0,0])*out[2]
        trans_info(1,i)=a[1,0]+(a[1,1]-a[1,0])*out[1]+(a[1,2]-a[1,0])*out[2]
        trans_info(2,i)=a[2,0]+(a[2,1]-a[2,0])*out[1]+(a[2,2]-a[2,0])*out[2]

        uvloc[*,i]=out[1:2]
	oplot, [out[1]]*dx,[out[2]]*dy,psym=3,color=200
        xx=out[1]*dx
        yy=out[2]*dy

        if i eq 0 then help,xx,l,yy

        tmp_through=area_through
        if xx gt -wid/2. && xx lt wid/2. then begin
            if yy gt -ht/2. && yy lt ht/2. then area_through+=1
            if yy lt -ht/2. && (xx^2+(yy+ht/2.)^2) lt (wid/2.)^2 then area_through+=1
            if yy gt  ht/2. && (xx^2+(yy-ht/2.)^2) lt (wid/2.)^2 then area_through+=1
        endif        
        if area_through eq tmp_through then trans_info(11,i)=0
        area_total+=1
    ENDFOR
    area_shade=area_total-area_through
    attenuation=area_shade/area_total
    *u.atnu1=attenuation
    widget_control,u.id.atnu1,set_value=num2str(*u.atnu1)
    
    oplot,wid/2*[1.0,1.0],ht/2.0*[-1.0,1.0]
    oplot,wid/2*[-1.0,-1.0],ht/2.0*[-1.0,1.0]
    th=make(0,!pi,50)
    oplot,wid/2.0*cos(th),ht/2.0+wid/2.0*sin(th)
    oplot,wid/2.0*cos(th),-1.0*ht/2.0-wid/2.0*sin(th)
;*****************************************************************

;*****************************************************************
;vessel points
    vpts=u.info.position.he.vessel
    vpts/=1.0e3
    a=vpts
    dx=sqrt((a[0,0]-a[0,1])^2+(a[1,0]-a[1,1])^2+(a[2,0]-a[2,1])^2)
    dy=sqrt((a[0,2]-a[0,1])^2+(a[1,2]-a[1,1])^2+(a[2,2]-a[2,1])^2)
    wid=203.20/1.0e3
    ht=425.45/1.0e3
    area_through=0.
    area_total=0.		
;make puncture plot on the vacuum plane
    widget_control,u.id.draw2,get_value=draw_win
    wset,draw_win
    plot, [0],[0],xr=[-0.4,0.4],yr=[-0.4,0.4],/xsty,/ysty,position=[0.1,0.1,0.9,0.9];,tit='Vaccum Vessel'
    uvloc=fltarr(2,n_mir)
    FOR i=0,n_mir-1 DO BEGIN
	l=[[m_pts[*,i]],$
           [r_pts[*,i]]]
        out=line_plane_int(l,bpts)
        uvloc[*,i]=out[1:2]
        oplot, [out[1]]*dx,[out[2]]*dy,psym=3,color=200
        xx=out[1]*dx
        yy=out[2]*dy
        tmp_through=area_through
        if xx gt -wid/2. && xx lt wid/2. then begin
            if yy gt -ht/2. && yy lt ht/2. then area_through+=1
            if yy lt -ht/2. && (xx^2+(yy+ht/2.)^2) lt (wid/2.)^2 then area_through+=1
            if yy gt  ht/2. && (xx^2+(yy-ht/2.)^2) lt (wid/2.)^2 then area_through+=1
        endif
        if area_through eq tmp_through then trans_info(11,i)=0
        area_total+=1
    ENDFOR

    area_shade=area_total-area_through
    attenuation=area_shade/area_total
    *u.atnu2=attenuation
    widget_control,u.id.atnu2,set_value=num2str(*u.atnu2)

    oplot,wid/2*[1.0,1.0],ht/2.0*[-1.0,1.0]
    oplot,wid/2*[-1.0,-1.0],ht/2.0*[-1.0,1.0]
    th=make(0,!pi,50)
    oplot,wid/2.0*cos(th),ht/2.0+wid/2.0*sin(th)
    oplot,wid/2.0*cos(th),-1.0*ht/2.0-wid/2.0*sin(th)
;*****************************************************************

;*****************************************************************
;make blank plot on the hemount plane
    widget_control,u.id.draw3,get_value=draw_win
    wset,draw_win   
    plot, [0],[0],xr=[-0.1,0.05],yr=[-0.2,0.2],/xsty,/ysty,position=[0.1,0.1,.9,.9],nodata=1

;*****************************************************************

;*****************************************************************
;reducer points
    redpts=u.info.position.he.reducer
    redpts/=1.0e3
    a=redpts
    dx=sqrt((a[0,0]-a[0,1])^2+(a[1,0]-a[1,1])^2+(a[2,0]-a[2,1])^2)
    dy=sqrt((a[0,2]-a[0,1])^2+(a[1,2]-a[1,1])^2+(a[2,2]-a[2,1])^2)
    area_shade=0.
    area_total=0.		

;reducer plane equation A*x+B*y+C*z+D=0
;see http://en.wikipedia.org/wiki/Plane_(geometry) and http://local.wasp.uwa.edu.au/~pbourke/geometry/planeeq/
    plane_x=[0,a[0,0],a[0,1],a[0,2]]
    plane_y=[0,a[1,0],a[1,1],a[1,2]]
    plane_z=[0,a[2,0],a[2,1],a[2,2]]

    Plane_A=plane_y(1)*(plane_z(2)-plane_z(3))+plane_y(2)*(plane_z(3)-plane_z(1))+plane_y(3)*(plane_z(1)-plane_z(2))
    Plane_B=plane_z(1)*(plane_x(2)-plane_x(3))+plane_z(2)*(plane_x(3)-plane_x(1))+plane_z(3)*(plane_x(1)-plane_x(2))
    Plane_C=plane_x(1)*(plane_y(2)-plane_y(3))+plane_x(2)*(plane_y(3)-plane_y(1))+plane_x(3)*(plane_y(1)-plane_y(2))
    Plane_D=plane_x(1)*(plane_y(2)*plane_z(3)-plane_y(3)*plane_z(2))+plane_x(2)*(plane_y(3)*plane_z(1)-plane_y(1)*plane_z(3))+plane_x(3)*(plane_y(1)*plane_z(2)-plane_y(2)*plane_z(1))
    Plane_D=-Plane_D
  
    
;make puncture plot on the reducer flange plane
    widget_control,u.id.draw4,get_value=draw_win
    wset,draw_win
    plot, [0],[0],xr=[-1.5,1.5]*dx,yr=[-1.5,1.5]*dy,/xsty,/ysty,position=[0.1,0.1,0.9,0.9]  ;,tit='Reducer Flange'
    uvloc=fltarr(2,n_mir)

    widget_control,u.id.xi_slider,get_value=xi
    widget_control,u.id.zeta_slider,get_value=zeta

    lambda=u.lambda(xi,zeta)*1.e3 ;in mAngstrom
    yy=(lambda-u.info.lamd)/u.info.delta_lamd  

    ;u.info,lamd=3000
    ;u.info.delta_lamd=2
    
    FOR i=0,n_mir-1 DO BEGIN
        l=[[m_pts[*,i]],$
           [r_pts[*,i]]]
        out=line_plane_int(l,redpts)


        trans_info(3,i)=a[0,0]+(a[0,1]-a[0,0])*out[1]+(a[0,2]-a[0,0])*out[2]
        trans_info(4,i)=a[1,0]+(a[1,1]-a[1,0])*out[1]+(a[1,2]-a[1,0])*out[2]
        trans_info(5,i)=a[2,0]+(a[2,1]-a[2,0])*out[1]+(a[2,2]-a[2,0])*out[2]

        zbport=trans_info(2,i)
        zflang=trans_info(5,i)
        ybport=trans_info(1,i)
        yflang=trans_info(4,i)        
        xbport=trans_info(0,i)
        xflang=trans_info(3,i)
        
        
;to get the incident angle; we have 3 points on the reducer plane, and
;we have the 2 points of the ray(one is in the bport plane, the other
;is in the reducer plane. see http://en.wikipedia.org/wiki/Plane_(geometry)

        dist_pts2pts=sqrt((xbport-xflang)^2+(ybport-yflang)^2+(zbport-zflang)^2)
        dist_pts2plane=abs(Plane_A*xbport+Plane_B*ybport+Plane_C*zbport+Plane_D)/sqrt(Plane_A^2+Plane_B^2+Plane_C^2)
        trans_info(6,i)=dist_pts2plane/dist_pts2pts   ;cos(theta)


        trans_info(7,i)=u.info.h1/trans_info(6,i)
        xx=(trans_info(7,i)-u.info.h1)/(u.info.h1*u.info.delta_h1)
        trans_info(8,i)=interpolate(u.info.transmission_1,xx,yy)       

        trans_info(9,i)=u.info.h2*sqrt(1+(trans_info(6,i)^2))               
        xx=(trans_info(9,i)-u.info.h2)/(u.info.h2*u.info.delta_h2)
        trans_info(10,i)=interpolate(u.info.transmission_2,xx,yy)

        uvloc[*,i]=out[1:2]
        oplot, [out[1]]*dx,[out[2]]*dy,psym=3,color=200
        tmp_shade=area_shade
        if ([out[1]])^2+([out[2]])^2 gt 1 then area_shade+=1
        if area_shade eq tmp_shade+1 then trans_info(11,i)=0
        area_total+=1
    ENDFOR

    attenuation=area_shade/area_total
    *u.atnu4=attenuation
    widget_control,u.id.atnu4,set_value=num2str(*u.atnu4)    

;here in calculating the transmission I have omitted the blocked/vignetting rays 
 
    subscp=where(trans_info(11,*) ne 0)
    
    if subscp(0) eq -1 then begin 
        total_not_vignet=0. 
        mean_trans=0
    endif else begin
        total_not_vignet=(n(subscp)+1.)/(n_mir)
        mean_trans=mean(trans_info(8,subscp))*u.info.p_h1+mean(trans_info(10,subscp))*u.info.p_h2
    endelse
    
    *u.trans_info_h=trans_info
    tmp_trans=*u.transmission
    total_trans=total_not_vignet*mean_trans
    tmp_trans(module-1,xi,zeta)=total_trans
    *u.transmission=tmp_trans
   
    min_cos=min(trans_info(6,*))
    max_angle=acos(min_cos)*180/!pi
    ave_angle=mean(acos(trans_info(6,*)))*180/!pi
    
    widget_control,u.id.trans,set_value=num2str(mean_trans)
    widget_control,u.id.max_angle,set_value=num2str(max_angle)
    widget_control,u.id.ave_angle,set_value=num2str(ave_angle)
    widget_control,u.id.total_atnu,set_value=num2str(1-total_not_vignet)
    widget_control,u.id.total_trans,set_value=num2str(total_trans)


    th=make(0,2.0*!pi,100)
    oplot, cos(th)*dx,sin(th)*dy
 
END

;*****************************************************************
   
PRO w_hlike_ar_vignetting,u

    shot=u.info.shot
    module=u.info.module
    help,shot,module	

    t1=u.info.t1
    t2=u.info.t2
    widget_control,u.id.xi_slider,get_value=xi
    widget_control,u.id.zeta_slider,get_value=zeta   

    det=module

    ;setup detector points
    info=u.load_info
    x0=info.det.x0
    x1=info.det.x1
    x2=info.det.x2
    det_pos=[xi[0],zeta[0]]*172.0e-6
    d_pt=genpos_det2xyz(x0,x1,x2,det_pos)
        
    ;setup mirror points
    ny=u.info.h_ny
    nz=u.info.h_nz
    ;ny=80
    ;nz=80
    m_pts=genpos_grid_sphere(info,ny=ny,nz=nz,/circ)
    n_mir=n(m_pts[0,*])

    trans_info=fltarr(12,n_mir)
    trans_info(11,*)=1

    ;define reflected points
    r_pts=fltarr(3,n_mir)
    FOR i=0L,n_mir-1 DO BEGIN
        xs=m_pts[*,i]
        x1=d_pt
        n_vec=[info.m.rad,0,0]-xs
        n_mag=sqrt(total(n_vec*n_vec))
        n=n_vec/n_mag
        r_pts[*,i]=genpos_spherical_refl(x1,xs,n)
    ENDFOR

;*****************************************************************
    ;b-port points
    bpts=u.info.position.h.bport	
    bpts/=1.0e3
    a=bpts
    dx=sqrt((a[0,0]-a[0,1])^2+(a[1,0]-a[1,1])^2+(a[2,0]-a[2,1])^2)
    dy=sqrt((a[0,2]-a[0,1])^2+(a[1,2]-a[1,1])^2+(a[2,2]-a[2,1])^2)

    wid=89.77/1.0e3
    ht=279.40/1.0e3
		
    ;make puncture plot on the b-port plane
    widget_control,u.id.draw1,get_value=draw_win
    wset,draw_win
    plot, [0],[0],xr=[-0.25,0.25],yr=[-0.25,0.25],/xsty,/ysty,position=[0.1,0.1,.9,.9] ;,tit='B-Port Flange'
    uvloc=fltarr(2,n_mir)
    area_through=0.
    area_total=0.
    FOR i=0L,n_mir-1 DO BEGIN
        l=[[m_pts[*,i]],$
           [r_pts[*,i]]]
        out=line_plane_int(l,bpts)

        trans_info(0,i)=a[0,0]+(a[0,1]-a[0,0])*out[1]+(a[0,2]-a[0,0])*out[2]
        trans_info(1,i)=a[1,0]+(a[1,1]-a[1,0])*out[1]+(a[1,2]-a[1,0])*out[2]
        trans_info(2,i)=a[2,0]+(a[2,1]-a[2,0])*out[1]+(a[2,2]-a[2,0])*out[2]

        uvloc[*,i]=out[1:2]
        oplot, [out[1]]*dx,[out[2]]*dy,psym=3,color=200
        xx=out[1]*dx
        yy=out[2]*dy
        tmp_through=area_through
        if xx gt -wid/2. && xx lt wid/2. then begin
            if yy gt -ht/2. && yy lt ht/2. then area_through+=1
            if yy lt -ht/2. && (xx^2+(yy+ht/2.)^2) lt (wid/2.)^2 then area_through+=1
            if yy gt  ht/2. && (xx^2+(yy-ht/2.)^2) lt (wid/2.)^2 then area_through+=1
        endif
        if area_through eq tmp_through then trans_info(11,i)=0
        area_total+=1
    ENDFOR
    
    area_shade=area_total-area_through
    attenuation=area_shade/area_total
    *u.atnu1=attenuation
    widget_control,u.id.atnu1,set_value=num2str(*u.atnu1)
    
    oplot,wid/2*[1.0,1.0],ht/2.0*[-1.0,1.0]
    oplot,wid/2*[-1.0,-1.0],ht/2.0*[-1.0,1.0]
    th=make(0,!pi,50)
    oplot,wid/2.0*cos(th),ht/2.0+wid/2.0*sin(th)
    oplot,wid/2.0*cos(th),-1.0*ht/2.0-wid/2.0*sin(th)
;******************************************************************

;******************************************************************
    ;vessel points
    vpts=u.info.position.h.vessel		
    vpts/=1.0e3
    a=vpts
    dx=sqrt((a[0,0]-a[0,1])^2+(a[1,0]-a[1,1])^2+(a[2,0]-a[2,1])^2)
    dy=sqrt((a[0,2]-a[0,1])^2+(a[1,2]-a[1,1])^2+(a[2,2]-a[2,1])^2)
    wid=203.20/1.0e3
    ht=425.45/1.0e3
    area_through=0.
    area_total=0.
    ;make puncture plot on the vacuum plane
    widget_control,u.id.draw2,get_value=draw_win
    wset,draw_win
    plot, [0],[0],xr=[-0.4,0.4],yr=[-0.4,0.4],/xsty,/ysty,position=[0.1,0.1,.9,.9] ;,tit='Vaccum Vessel'
    uvloc=fltarr(2,n_mir)
    FOR i=0L,n_mir-1 DO BEGIN
        l=[[m_pts[*,i]],$
           [r_pts[*,i]]]
        out=line_plane_int(l,bpts)
        
        uvloc[*,i]=out[1:2]
        oplot, [out[1]]*dx,[out[2]]*dy,psym=3,color=200
        xx=out[1]*dx
        yy=out[2]*dy
        tmp_through=area_through
        if xx gt -wid/2. && xx lt wid/2. then begin
            if yy gt -ht/2. && yy lt ht/2. then area_through+=1
            if yy lt -ht/2. && (xx^2+(yy+ht/2.)^2) lt (wid/2.)^2 then area_through+=1
            if yy gt  ht/2. && (xx^2+(yy-ht/2.)^2) lt (wid/2.)^2 then area_through+=1
        endif
        if area_through eq tmp_through then trans_info(11,i)=0
        area_total+=1                
    ENDFOR

    area_shade=area_total-area_through
    attenuation=area_shade/area_total
    *u.atnu2=attenuation
    widget_control,u.id.atnu2,set_value=num2str(*u.atnu2)
 
   
    oplot,wid/2*[1.0,1.0],ht/2.0*[-1.0,1.0]
    oplot,wid/2*[-1.0,-1.0],ht/2.0*[-1.0,1.0]
    th=make(0,!pi,50)
    oplot,wid/2.0*cos(th),ht/2.0+wid/2.0*sin(th)
    oplot,wid/2.0*cos(th),-1.0*ht/2.0-wid/2.0*sin(th)

;******************************************************************

;******************************************************************

    ;he-like points
    hemount=u.info.position.h.he_crystal
    hemount/=1.0e3
    a=hemount
    dx=-1.0*sqrt((a[0,0]-a[0,1])^2+(a[1,0]-a[1,1])^2+(a[2,0]-a[2,1])^2)
    dy=-1.0*sqrt((a[0,2]-a[0,1])^2+(a[1,2]-a[1,1])^2+(a[2,2]-a[2,1])^2)
    
    ;make puncture plot on the hemount plane
    widget_control,u.id.draw3,get_value=draw_win
    wset,draw_win
    plot, [0],[0],xr=[2*dx,-dx],yr=[2*dy,-2*dy],/xsty,/ysty,position=[0.1,0.1,.9,.9] ;,tit='He-like Crystal Moun
    uvloc=fltarr(2,n_mir)
    
    area_shade=0.
    area_total=0.
    FOR i=0L,n_mir-1 DO BEGIN
        l=[[m_pts[*,i]],$
           [r_pts[*,i]]]
        out=line_plane_int(l,hemount)
        uvloc[*,i]=out[1:2]
        oplot, [out[1]]*dx,[out[2]]*dy,psym=3,color=200
        tmp_shade=area_shade
        if [out[2]] gt 0 && [out[2]] lt 1 && [out[1]] gt 0 && [out[1]] lt 1 then area_shade+=1
        if area_shade eq tmp_shade+1 then trans_info(11,i)=0
        area_total+=1
    ENDFOR
    
    attenuation=area_shade/area_total
    *u.atnu3=attenuation
    widget_control,u.id.atnu3,set_value=num2str(*u.atnu3)
    
    tmp=where(uvloc[1,*] GE 0)
    IF tmp[0] EQ -1 THEN frac=0.0 ELSE frac=(n(tmp)+1.0)/n_mir

    
    oplot, [0,1,1,0,0]*dx,[0,0,1,1,0]*dy

;******************************************************************

;******************************************************************

;reducer points
    redpts=u.info.position.h.reducer
    redpts/=1.0e3
    a=redpts
    dx=sqrt((a[0,0]-a[0,1])^2+(a[1,0]-a[1,1])^2+(a[2,0]-a[2,1])^2)
    dy=sqrt((a[0,2]-a[0,1])^2+(a[1,2]-a[1,1])^2+(a[2,2]-a[2,1])^2)
    area_shade=0.
    area_total=0.

;reducer plane equation A*x+B*y+C*z+D=0
;see http://en.wikipedia.org/wiki/Plane_(geometry) and http://local.wasp.uwa.edu.au/~pbourke/geometry/planeeq/
    plane_x=[0,a[0,0],a[0,1],a[0,2]]
    plane_y=[0,a[1,0],a[1,1],a[1,2]]
    plane_z=[0,a[2,0],a[2,1],a[2,2]]

    Plane_A=plane_y(1)*(plane_z(2)-plane_z(3))+plane_y(2)*(plane_z(3)-plane_z(1))+plane_y(3)*(plane_z(1)-plane_z(2))
    Plane_B=plane_z(1)*(plane_x(2)-plane_x(3))+plane_z(2)*(plane_x(3)-plane_x(1))+plane_z(3)*(plane_x(1)-plane_x(2))
    Plane_C=plane_x(1)*(plane_y(2)-plane_y(3))+plane_x(2)*(plane_y(3)-plane_y(1))+plane_x(3)*(plane_y(1)-plane_y(2))
    Plane_D=plane_x(1)*(plane_y(2)*plane_z(3)-plane_y(3)*plane_z(2))+plane_x(2)*(plane_y(3)*plane_z(1)-plane_y(1)*plane_z(3))+plane_x(3)*(plane_y(1)*plane_z(2)-plane_y(2)*plane_z(1))
    Plane_D=-Plane_D
  
     
;make puncture plot on the reducer flange plane
    widget_control,u.id.draw4,get_value=draw_win
    wset,draw_win
    plot, [0],[0],xr=[-1.5,1.5]*dx,yr=[-1.5,1.5]*dy,/xsty,/ysty,position=[0.1,0.1,.9,.9] ;,tit='Reducer Flange'
    uvloc=fltarr(2,n_mir)
    

    widget_control,u.id.xi_slider,get_value=xi
    widget_control,u.id.zeta_slider,get_value=zeta 

    lambda=u.lambda(xi,zeta)*1.e3 ;in mAngstrom
    yy=(lambda-u.info.lamd)/u.info.delta_lamd



    FOR i=0L,n_mir-1 DO BEGIN
        l=[[m_pts[*,i]],$
           [r_pts[*,i]]]
        out=line_plane_int(l,redpts)

        trans_info(3,i)=a[0,0]+(a[0,1]-a[0,0])*out[1]+(a[0,2]-a[0,0])*out[2]
        trans_info(4,i)=a[1,0]+(a[1,1]-a[1,0])*out[1]+(a[1,2]-a[1,0])*out[2]
        trans_info(5,i)=a[2,0]+(a[2,1]-a[2,0])*out[1]+(a[2,2]-a[2,0])*out[2]

        zbport=trans_info(2,i)
        zflang=trans_info(5,i)
        ybport=trans_info(1,i)
        yflang=trans_info(4,i)        
        xbport=trans_info(0,i)
        xflang=trans_info(3,i)
        
;to get the incident angle; we have 3 points on the reducer plane, and
;we have the 2 points of the ray(one is in the bport plane, the other
;is in the reducer plane. see http://en.wikipedia.org/wiki/Plane_(geometry))

        dist_pts2pts=sqrt((xbport-xflang)^2+(ybport-yflang)^2+(zbport-zflang)^2)
        dist_pts2plane=abs(Plane_A*xbport+Plane_B*ybport+Plane_C*zbport+Plane_D)/sqrt(Plane_A^2+Plane_B^2+Plane_C^2)
        trans_info(6,i)=dist_pts2plane/dist_pts2pts   ;cos(theta)

        trans_info(7,i)=u.info.h1/trans_info(6,i)
        xx=(trans_info(7,i)-u.info.h1)/(u.info.h1*u.info.delta_h1)
        trans_info(8,i)=interpolate(u.info.transmission_1,xx,yy)       

        trans_info(9,i)=u.info.h2*sqrt(1+(trans_info(6,i)^2))               
        xx=(trans_info(9,i)-u.info.h2)/(u.info.h2*u.info.delta_h2)
        trans_info(10,i)=interpolate(u.info.transmission_2,xx,yy)

       
        uvloc[*,i]=out[1:2]
        oplot, [out[1]]*dx,[out[2]]*dy,psym=3,color=200
        tmp_shade=area_shade
        if ([out[1]])^2+([out[2]])^2 gt 1 then area_shade+=1
        if area_shade eq tmp_shade+1 then trans_info(11,i)=0
        area_total+=1
    ENDFOR
    
    attenuation=area_shade/area_total
    *u.atnu4=attenuation
    widget_control,u.id.atnu4,set_value=num2str(*u.atnu4)

;here in calculating the transmission I have omitted the blocked/vignetting rays 

    subscp=where(trans_info(11,*) ne 0)
    
    if subscp(0) eq -1 then begin 
        total_not_vignet=0. 
        mean_trans=0
    endif else begin
        total_not_vignet=(n(subscp)+1.)/(n_mir)
        mean_trans=mean(trans_info(8,subscp))*u.info.p_h1+mean(trans_info(10,subscp))*u.info.p_h2
    endelse
    
    *u.trans_info_h=trans_info
    tmp_trans=*u.transmission
    total_trans=total_not_vignet*mean_trans
    tmp_trans(module-1,xi,zeta)=total_trans
    *u.transmission=tmp_trans

    min_cos=min(trans_info(6,*))
    max_angle=acos(min_cos)*180/!pi
    ave_angle=mean(acos(trans_info(6,*)))*180/!pi
    
    widget_control,u.id.trans,set_value=num2str(mean_trans)
    widget_control,u.id.max_angle,set_value=num2str(max_angle)
    widget_control,u.id.ave_angle,set_value=num2str(ave_angle)
    widget_control,u.id.total_atnu,set_value=num2str(1-total_not_vignet)
    widget_control,u.id.total_trans,set_value=num2str(total_trans)


    th=make(0,2.0*!pi,100)
    oplot, cos(th)*dx,sin(th)*dy
    reducer=frac	

END




PRO plot_vignetting,u

    shot=u.info.shot
    module=u.info.module
    t1=u.info.t1
    t2=u.info.t2
    widget_control,u.id.xi_slider,get_value=xi
    widget_control,u.id.zeta_slider,get_value=zeta 
    if module eq 4 then w_hlike_ar_vignetting,u else w_helike_ar_vignetting,u
END

    


PRO reset1_u,u
    widget_control,u.id.he_ny,set_value=strtrim(u.reset.he_ny,2)
    widget_control,u.id.he_nz,set_value=strtrim(u.reset.he_nz,2)
    widget_control,u.id.h_ny,set_value=strtrim(u.reset.h_ny,2)
    widget_control,u.id.h_nz,set_value=strtrim(u.reset.h_nz,2)
   
    u.info.he_ny=u.reset.he_ny
    u.info.he_nz=u.reset.he_nz
    u.info.h_ny=u.reset.h_ny
    u.info.h_nz=u.reset.h_nz
    
END


PRO reset2_u,u
    widget_control,u.id.shotid,set_value=strtrim(u.reset.shot,2)
    widget_control,u.id.module,set_value=strtrim(u.reset.module,2)
    widget_control,u.id.t1,set_value=num2str(u.reset.t1,dp=2)
    widget_control,u.id.t2,set_value=num2str(u.reset.t2,dp=2)
    widget_control,u.id.xi_slider,set_value=u.reset.xi
    widget_control,u.id.zeta_slider,set_value=u.reset.zeta

    load_info=hirexsr_load_info(u.reset.module,shot=u.reset.shot)

    u.info.shot=u.reset.shot
    u.info.module=u.reset.module
    u.info.t1=u.reset.t1
    u.info.t2=u.reset.t2
    u.load_info=load_info
    *u.image_index=1    
END

PRO reset3_u,u

    widget_control,u.id.he_re_x1,set_value=num2str(u.reset.position.he.reducer(0),dp=2)
    widget_control,u.id.he_re_y1,set_value=num2str(u.reset.position.he.reducer(1),dp=2)
    widget_control,u.id.he_re_z1,set_value=num2str(u.reset.position.he.reducer(2),dp=2)
    widget_control,u.id.he_re_x2,set_value=num2str(u.reset.position.he.reducer(3),dp=2)
    widget_control,u.id.he_re_y2,set_value=num2str(u.reset.position.he.reducer(4),dp=2)
    widget_control,u.id.he_re_z2,set_value=num2str(u.reset.position.he.reducer(5),dp=2)
    widget_control,u.id.he_re_x3,set_value=num2str(u.reset.position.he.reducer(6),dp=2)
    widget_control,u.id.he_re_y3,set_value=num2str(u.reset.position.he.reducer(7),dp=2)
    widget_control,u.id.he_re_z3,set_value=num2str(u.reset.position.he.reducer(8),dp=2)

    widget_control,u.id.he_bp_x1,set_value=num2str(u.reset.position.he.bport(0),dp=2)
    widget_control,u.id.he_bp_y1,set_value=num2str(u.reset.position.he.bport(1),dp=2)
    widget_control,u.id.he_bp_z1,set_value=num2str(u.reset.position.he.bport(2),dp=2)
    widget_control,u.id.he_bp_x2,set_value=num2str(u.reset.position.he.bport(3),dp=2)
    widget_control,u.id.he_bp_y2,set_value=num2str(u.reset.position.he.bport(4),dp=2)
    widget_control,u.id.he_bp_z2,set_value=num2str(u.reset.position.he.bport(5),dp=2)
    widget_control,u.id.he_bp_x3,set_value=num2str(u.reset.position.he.bport(6),dp=2)
    widget_control,u.id.he_bp_y3,set_value=num2str(u.reset.position.he.bport(7),dp=2)
    widget_control,u.id.he_bp_z3,set_value=num2str(u.reset.position.he.bport(8),dp=2)
    
    widget_control,u.id.he_ve_x1,set_value=num2str(u.reset.position.he.vessel(0),dp=2)
    widget_control,u.id.he_ve_y1,set_value=num2str(u.reset.position.he.vessel(1),dp=2)
    widget_control,u.id.he_ve_z1,set_value=num2str(u.reset.position.he.vessel(2),dp=2)
    widget_control,u.id.he_ve_x2,set_value=num2str(u.reset.position.he.vessel(3),dp=2)
    widget_control,u.id.he_ve_y2,set_value=num2str(u.reset.position.he.vessel(4),dp=2)
    widget_control,u.id.he_ve_z2,set_value=num2str(u.reset.position.he.vessel(5),dp=2)
    widget_control,u.id.he_ve_x3,set_value=num2str(u.reset.position.he.vessel(6),dp=2)
    widget_control,u.id.he_ve_y3,set_value=num2str(u.reset.position.he.vessel(7),dp=2)
    widget_control,u.id.he_ve_z3,set_value=num2str(u.reset.position.he.vessel(8),dp=2)

    widget_control,u.id.h_re_x1,set_value=num2str(u.reset.position.h.reducer(0),dp=2)
    widget_control,u.id.h_re_y1,set_value=num2str(u.reset.position.h.reducer(1),dp=2)
    widget_control,u.id.h_re_z1,set_value=num2str(u.reset.position.h.reducer(2),dp=2)
    widget_control,u.id.h_re_x2,set_value=num2str(u.reset.position.h.reducer(3),dp=2)
    widget_control,u.id.h_re_y2,set_value=num2str(u.reset.position.h.reducer(4),dp=2)
    widget_control,u.id.h_re_z2,set_value=num2str(u.reset.position.h.reducer(5),dp=2)
    widget_control,u.id.h_re_x3,set_value=num2str(u.reset.position.h.reducer(6),dp=2)
    widget_control,u.id.h_re_y3,set_value=num2str(u.reset.position.h.reducer(7),dp=2)
    widget_control,u.id.h_re_z3,set_value=num2str(u.reset.position.h.reducer(8),dp=2)
    
    widget_control,u.id.h_bp_x1,set_value=num2str(u.reset.position.h.bport(0),dp=2)
    widget_control,u.id.h_bp_y1,set_value=num2str(u.reset.position.h.bport(1),dp=2)
    widget_control,u.id.h_bp_z1,set_value=num2str(u.reset.position.h.bport(2),dp=2)
    widget_control,u.id.h_bp_x2,set_value=num2str(u.reset.position.h.bport(3),dp=2)
    widget_control,u.id.h_bp_y2,set_value=num2str(u.reset.position.h.bport(4),dp=2)
    widget_control,u.id.h_bp_z2,set_value=num2str(u.reset.position.h.bport(5),dp=2)
    widget_control,u.id.h_bp_x3,set_value=num2str(u.reset.position.h.bport(6),dp=2)
    widget_control,u.id.h_bp_y3,set_value=num2str(u.reset.position.h.bport(7),dp=2)
    widget_control,u.id.h_bp_z3,set_value=num2str(u.reset.position.h.bport(8),dp=2)
    
    widget_control,u.id.h_ve_x1,set_value=num2str(u.reset.position.h.vessel(0),dp=2)
    widget_control,u.id.h_ve_y1,set_value=num2str(u.reset.position.h.vessel(1),dp=2)
    widget_control,u.id.h_ve_z1,set_value=num2str(u.reset.position.h.vessel(2),dp=2)
    widget_control,u.id.h_ve_x2,set_value=num2str(u.reset.position.h.vessel(3),dp=2)
    widget_control,u.id.h_ve_y2,set_value=num2str(u.reset.position.h.vessel(4),dp=2)
    widget_control,u.id.h_ve_z2,set_value=num2str(u.reset.position.h.vessel(5),dp=2)
    widget_control,u.id.h_ve_x3,set_value=num2str(u.reset.position.h.vessel(6),dp=2)
    widget_control,u.id.h_ve_y3,set_value=num2str(u.reset.position.h.vessel(7),dp=2)
    widget_control,u.id.h_ve_z3,set_value=num2str(u.reset.position.h.vessel(8),dp=2)

    widget_control,u.id.h_cr_x1,set_value=num2str(u.reset.position.h.he_crystal(0),dp=2)
    widget_control,u.id.h_cr_y1,set_value=num2str(u.reset.position.h.he_crystal(1),dp=2)
    widget_control,u.id.h_cr_z1,set_value=num2str(u.reset.position.h.he_crystal(2),dp=2)
    widget_control,u.id.h_cr_x2,set_value=num2str(u.reset.position.h.he_crystal(3),dp=2)
    widget_control,u.id.h_cr_y2,set_value=num2str(u.reset.position.h.he_crystal(4),dp=2)
    widget_control,u.id.h_cr_z2,set_value=num2str(u.reset.position.h.he_crystal(5),dp=2)
    widget_control,u.id.h_cr_x3,set_value=num2str(u.reset.position.h.he_crystal(6),dp=2)
    widget_control,u.id.h_cr_y3,set_value=num2str(u.reset.position.h.he_crystal(7),dp=2)
    widget_control,u.id.h_cr_z3,set_value=num2str(u.reset.position.h.he_crystal(8),dp=2)   

    u.info.position=u.reset.position
END


PRO reset5_u,u
    widget_control,u.id.he_l_crystal,set_value=strtrim(u.reset.he_l_0,2)
    widget_control,u.id.he_phi_crystal,set_value=strtrim(u.reset.he_phi_0,2)
    widget_control,u.id.he_z_crystal,set_value=strtrim(u.reset.he_z_0,2)

    widget_control,u.id.he_l_slider,set_value=0
    widget_control,u.id.he_phi_slider,set_value=0
    widget_control,u.id.he_z_slider,set_value=0

    u.info.he_l_c=u.reset.he_l_0
    u.info.he_phi_c=u.reset.he_phi_0
    u.info.he_z_c=u.reset.he_z_0
    
    ;help,u.reset.he_z_0,u.info.he_z_c

    u.info.position=u.reset.position
END

PRO reset6_u,u
    widget_control,u.id.h_l_crystal,set_value=strtrim(u.reset.h_l_0,2)
    widget_control,u.id.h_phi_crystal,set_value=strtrim(u.reset.h_phi_0,2)
    widget_control,u.id.h_z_crystal,set_value=strtrim(u.reset.h_z_0,2)

    widget_control,u.id.h_l_slider,set_value=0
    widget_control,u.id.h_phi_slider,set_value=0
    widget_control,u.id.h_z_slider,set_value=0

    u.info.h_l_c=u.reset.h_l_0
    u.info.h_phi_c=u.reset.h_phi_0
    u.info.h_z_c=u.reset.h_z_0
    
    ;help,u.reset.h_z_0,u.info.h_z_c

    u.info.position=u.reset.position
END


PRO get_image,u,image=image

    shot=u.info.shot
    module=u.info.module
    t1=u.info.t1
    t2=u.info.t2

    hirexsr_load_image,shot,module,image,t
    t_start=min(where(t gt t1))
    t_end=min(where(t gt t2))
    image=total(image(*,*,t_start:t_end),3)
    
    temp=image
    temp2 = temp[sort(temp)]
    i95 = temp2[0.98*n_elements(temp)]*1.0
    image(where(image GT i95)) = i95
    scale = i95/255
    if scale EQ 0 then scale = 1.
    image=image/scale
    image=reverse(reverse(transpose(image),2)) 
    *u.image=image
    *u.image_index=0
END


PRO plot_tv,u

    widget_control,u.id.xi_slider,get_value=xi
    widget_control,u.id.zeta_slider,get_value=zeta

    loadct,0,/silent
    widget_control,u.id.draw5,get_value=draw_win
    wset,draw_win
    tvimage,*u.image,position=[0,0,1,1]
    loadct,12,/silent
    xx=[1,195] 
    yy=[1,487]
    plot,xx,yy,xstyle=1,ystyle=1,nodata=1,noerase=1,position=[0,0,1,1]
    plots,[xi,xi],[0,486],color=200,thick=1    
    plots,[0,194],[zeta,zeta],color=200,thick=1

END




PRO plot_detector,u
    if *u.image_index eq 1 then get_image,u
    plot_tv,u

END

PRO PLOT_CRYSTAL_HE,u


    widget_control,u.id.he_l_slider,get_value=he_l_slider
    widget_control,u.id.he_phi_slider,get_value=he_phi_slider
    widget_control,u.id.he_z_slider,get_value=he_z_slider    
   
    delta_l=he_l_slider*0.05
    delta_phi=he_phi_slider*0.01*!pi/180
    delta_z=he_z_slider*0.05

    he_l=u.reset.he_l_0+ delta_l

    he_phi=u.reset.he_phi_0+delta_phi
    he_z=u.reset.he_z_0+delta_z

    widget_control,u.id.he_l_crystal,set_value=num2str(he_l,dp=2)
    widget_control,u.id.he_phi_crystal,set_value=num2str(he_phi,dp=2)
    widget_control,u.id.he_z_crystal,set_value=num2str(he_z,dp=2)

    widget_control,u.id.draw6,get_value=draw_win
    wset,draw_win

    l1=u.info.he_l1
    l2=he_l

    phi=he_phi
    z=he_z 
    beta=u.info.beta
    theta=atan(l2*sin(beta)/(l1+l2*cos(beta)))
    alpha=beta-theta
    R=sqrt(l1^2+l2^2+2*l1*l2*cos(alpha+theta))
    
    X=l2*cos(beta)+l1
    Y=l2*sin(beta)

;=============
    angle=phi-beta
    he_ori=[x,y,z]
    
    vpts=u.info.position.tok.vessel
    bpts=u.info.position.tok.bport
    redpts=u.info.position.tok.reducer
     
    vpts_he=tran_position(vpts,he_ori,angle=angle)
    bpts_he=tran_position(bpts,he_ori,angle=angle)
    redpts_he=tran_position(redpts,he_ori,angle=angle)
    position_He={reducer:redpts_he,bport:bpts_he,vessel:vpts_he}
    u.info.position.he=position_he
    widget_control,u.id.he_re_x1,set_value=num2str(u.info.position.he.reducer(0),dp=2)
    widget_control,u.id.he_re_y1,set_value=num2str(u.info.position.he.reducer(1),dp=2)
    widget_control,u.id.he_re_z1,set_value=num2str(u.info.position.he.reducer(2),dp=2)
    widget_control,u.id.he_re_x2,set_value=num2str(u.info.position.he.reducer(3),dp=2)
    widget_control,u.id.he_re_y2,set_value=num2str(u.info.position.he.reducer(4),dp=2)
    widget_control,u.id.he_re_z2,set_value=num2str(u.info.position.he.reducer(5),dp=2)
    widget_control,u.id.he_re_x3,set_value=num2str(u.info.position.he.reducer(6),dp=2)
    widget_control,u.id.he_re_y3,set_value=num2str(u.info.position.he.reducer(7),dp=2)
    widget_control,u.id.he_re_z3,set_value=num2str(u.info.position.he.reducer(8),dp=2)
    
    widget_control,u.id.he_bp_x1,set_value=num2str(u.info.position.he.bport(0),dp=2)
    widget_control,u.id.he_bp_y1,set_value=num2str(u.info.position.he.bport(1),dp=2)
    widget_control,u.id.he_bp_z1,set_value=num2str(u.info.position.he.bport(2),dp=2)
    widget_control,u.id.he_bp_x2,set_value=num2str(u.info.position.he.bport(3),dp=2)
    widget_control,u.id.he_bp_y2,set_value=num2str(u.info.position.he.bport(4),dp=2)
    widget_control,u.id.he_bp_z2,set_value=num2str(u.info.position.he.bport(5),dp=2)
    widget_control,u.id.he_bp_x3,set_value=num2str(u.info.position.he.bport(6),dp=2)
    widget_control,u.id.he_bp_y3,set_value=num2str(u.info.position.he.bport(7),dp=2)
    widget_control,u.id.he_bp_z3,set_value=num2str(u.info.position.he.bport(8),dp=2)
    
    widget_control,u.id.he_ve_x1,set_value=num2str(u.info.position.he.vessel(0),dp=2)
    widget_control,u.id.he_ve_y1,set_value=num2str(u.info.position.he.vessel(1),dp=2)
    widget_control,u.id.he_ve_z1,set_value=num2str(u.info.position.he.vessel(2),dp=2)
    widget_control,u.id.he_ve_x2,set_value=num2str(u.info.position.he.vessel(3),dp=2)
    widget_control,u.id.he_ve_y2,set_value=num2str(u.info.position.he.vessel(4),dp=2)
    widget_control,u.id.he_ve_z2,set_value=num2str(u.info.position.he.vessel(5),dp=2)
    widget_control,u.id.he_ve_x3,set_value=num2str(u.info.position.he.vessel(6),dp=2)
    widget_control,u.id.he_ve_y3,set_value=num2str(u.info.position.he.vessel(7),dp=2)
    widget_control,u.id.he_ve_z3,set_value=num2str(u.info.position.he.vessel(8),dp=2)

;=============

    u.info.position.he_ori=he_ori
    u.info.he_l_c=he_l
    u.info.he_phi_c=he_phi
    u.info.he_z_c=he_z

    
    widget_control,u.id.he_l_crystal,set_value=num2str(u.info.he_l_c,dp=2)
    widget_control,u.id.he_phi_crystal,set_value=num2str(u.info.he_phi_c*180/!pi,dp=2)
    widget_control,u.id.he_z_crystal,set_value=num2str(u.info.he_z_c,dp=2)


    l2=float(l2)

    radius_wall=2000            ;in mm

    theta_circ=findgen(101)*2*!pi/100
    x_circ=radius_wall*sin(theta_circ)
    y_circ=radius_wall*cos(theta_circ)

    plot,x_circ,y_circ,/iso,xrange=[0,2.5*radius_wall],yrange=[-1.25*radius_wall,1.25*radius_wall],xstyle=1,ystyle=1,title='He-like Cyrstal Layout'
    x_cr=l1+l2*cos(beta)
    y_cr=l2*sin(beta)

    plots,[0,x_cr],[0,y_cr]
    plots,[l1,x_cr],[0,y_cr]
    plots,[0,l1],[0,0]

    length=1000
    plots,[x_cr,x_cr-length*cos(phi-beta)],[y_cr,y_cr+length*sin(phi-beta)],color=200
    plots,[x_cr,x_cr-length*sin(phi-beta)],[y_cr,y_cr-length*cos(phi-beta)],color=200

    delta_hemount=[delta_l*cos(beta),delta_l*sin(beta),delta_z]
    u.info.position.tok.he_crystal=u.reset.position.tok.he_crystal-[[delta_hemount],[delta_hemount],[delta_hemount]]
    
    he_ori=u.info.position.he_ori
    
    reducer=u.info.position.tok.reducer    
    bport=u.info.position.tok.bport
    vessel=u.info.position.tok.vessel
    hemount=u.info.position.tok.he_crystal

    u.info.position.he.reducer=tran_position(reducer,he_ori,angle=phi-beta)
    u.info.position.he.bport=tran_position(bport,he_ori,angle=phi-beta)
    u.info.position.he.vessel=tran_position(vessel,he_ori,angle=phi-beta)
    u.info.position.h.he_crystal=tran_position(hemount,he_ori,angle=phi-beta)


end



PRO PLOT_CRYSTAL_H,u


    widget_control,u.id.h_l_slider,get_value=h_l_slider
    widget_control,u.id.h_phi_slider,get_value=h_phi_slider
    widget_control,u.id.h_z_slider,get_value=h_z_slider    
   
    delta_l=h_l_slider*0.05
    delta_phi=h_phi_slider*0.01*!pi/180
    delta_z=h_z_slider*0.05

    h_l=u.reset.h_l_0+ delta_l

    h_phi=u.reset.h_phi_0+delta_phi
    h_z=u.reset.h_z_0+delta_z

    widget_control,u.id.h_l_crystal,set_value=num2str(h_l,dp=2)
    widget_control,u.id.h_phi_crystal,set_value=num2str(h_phi,dp=2)
    widget_control,u.id.h_z_crystal,set_value=num2str(h_z,dp=2)

    widget_control,u.id.draw7,get_value=draw_win
    wset,draw_win

    l1=u.info.h_l1
    l2=h_l

    phi=h_phi
    z=h_z 
    beta=u.info.beta
    theta=atan(l2*sin(beta)/(l1+l2*cos(beta)))
    alpha=beta-theta
    R=sqrt(l1^2+l2^2+2*l1*l2*cos(alpha+theta))
    
    X=l2*cos(beta)+l1
    Y=l2*sin(beta)

;=============
    angle=phi-beta
    h_ori=[x,y,z]
    
    vpts=u.info.position.tok.vessel
    bpts=u.info.position.tok.bport
    redpts=u.info.position.tok.reducer
    hemount=u.info.position.tok.he_crystal
 
    vpts_h=tran_position(vpts,h_ori,angle=angle)
    bpts_h=tran_position(bpts,h_ori,angle=angle)
    redpts_h=tran_position(redpts,h_ori,angle=angle)
    hemount_h=tran_position(hemount,h_ori,angle=angle)

    position_H={reducer:redpts_h,bport:bpts_h,vessel:vpts_h,he_crystal:hemount_h}
    u.info.position.h=position_h

    widget_control,u.id.h_re_x1,set_value=num2str(u.info.position.h.reducer(0),dp=2)
    widget_control,u.id.h_re_y1,set_value=num2str(u.info.position.h.reducer(1),dp=2)
    widget_control,u.id.h_re_z1,set_value=num2str(u.info.position.h.reducer(2),dp=2)
    widget_control,u.id.h_re_x2,set_value=num2str(u.info.position.h.reducer(3),dp=2)
    widget_control,u.id.h_re_y2,set_value=num2str(u.info.position.h.reducer(4),dp=2)
    widget_control,u.id.h_re_z2,set_value=num2str(u.info.position.h.reducer(5),dp=2)
    widget_control,u.id.h_re_x3,set_value=num2str(u.info.position.h.reducer(6),dp=2)
    widget_control,u.id.h_re_y3,set_value=num2str(u.info.position.h.reducer(7),dp=2)
    widget_control,u.id.h_re_z3,set_value=num2str(u.info.position.h.reducer(8),dp=2)
    
    widget_control,u.id.h_bp_x1,set_value=num2str(u.info.position.h.bport(0),dp=2)
    widget_control,u.id.h_bp_y1,set_value=num2str(u.info.position.h.bport(1),dp=2)
    widget_control,u.id.h_bp_z1,set_value=num2str(u.info.position.h.bport(2),dp=2)
    widget_control,u.id.h_bp_x2,set_value=num2str(u.info.position.h.bport(3),dp=2)
    widget_control,u.id.h_bp_y2,set_value=num2str(u.info.position.h.bport(4),dp=2)
    widget_control,u.id.h_bp_z2,set_value=num2str(u.info.position.h.bport(5),dp=2)
    widget_control,u.id.h_bp_x3,set_value=num2str(u.info.position.h.bport(6),dp=2)
    widget_control,u.id.h_bp_y3,set_value=num2str(u.info.position.h.bport(7),dp=2)
    widget_control,u.id.h_bp_z3,set_value=num2str(u.info.position.h.bport(8),dp=2)
    
    widget_control,u.id.h_ve_x1,set_value=num2str(u.info.position.h.vessel(0),dp=2)
    widget_control,u.id.h_ve_y1,set_value=num2str(u.info.position.h.vessel(1),dp=2)
    widget_control,u.id.h_ve_z1,set_value=num2str(u.info.position.h.vessel(2),dp=2)
    widget_control,u.id.h_ve_x2,set_value=num2str(u.info.position.h.vessel(3),dp=2)
    widget_control,u.id.h_ve_y2,set_value=num2str(u.info.position.h.vessel(4),dp=2)
    widget_control,u.id.h_ve_z2,set_value=num2str(u.info.position.h.vessel(5),dp=2)
    widget_control,u.id.h_ve_x3,set_value=num2str(u.info.position.h.vessel(6),dp=2)
    widget_control,u.id.h_ve_y3,set_value=num2str(u.info.position.h.vessel(7),dp=2)
    widget_control,u.id.h_ve_z3,set_value=num2str(u.info.position.h.vessel(8),dp=2)

    widget_control,u.id.h_cr_x1,set_value=num2str(u.info.position.h.he_crystal(0),dp=2)
    widget_control,u.id.h_cr_y1,set_value=num2str(u.info.position.h.he_crystal(1),dp=2)
    widget_control,u.id.h_cr_z1,set_value=num2str(u.info.position.h.he_crystal(2),dp=2)
    widget_control,u.id.h_cr_x2,set_value=num2str(u.info.position.h.he_crystal(3),dp=2)
    widget_control,u.id.h_cr_y2,set_value=num2str(u.info.position.h.he_crystal(4),dp=2)
    widget_control,u.id.h_cr_z2,set_value=num2str(u.info.position.h.he_crystal(5),dp=2)
    widget_control,u.id.h_cr_x3,set_value=num2str(u.info.position.h.he_crystal(6),dp=2)
    widget_control,u.id.h_cr_y3,set_value=num2str(u.info.position.h.he_crystal(7),dp=2)
    widget_control,u.id.h_cr_z3,set_value=num2str(u.info.position.h.he_crystal(8),dp=2) 
;=============

    u.info.position.h_ori=h_ori
    u.info.h_l_c=h_l
    u.info.h_phi_c=h_phi
    u.info.h_z_c=h_z
    
    widget_control,u.id.h_l_crystal,set_value=num2str(u.info.h_l_c,dp=2)
    widget_control,u.id.h_phi_crystal,set_value=num2str(u.info.h_phi_c*180/!pi,dp=2)
    widget_control,u.id.h_z_crystal,set_value=num2str(u.info.h_z_c,dp=2)


    l2=float(l2)

    radius_wall=2000            ;in mm

    theta_circ=findgen(101)*2*!pi/100
    x_circ=radius_wall*sin(theta_circ)
    y_circ=radius_wall*cos(theta_circ)

    plot,x_circ,y_circ,/iso,xrange=[0,2.5*radius_wall],yrange=[-1.25*radius_wall,1.25*radius_wall],xstyle=1,ystyle=1,title='H-like Cyrstal Layout'

    x_cr=l1+l2*cos(beta)
    y_cr=l2*sin(beta)

    plots,[0,x_cr],[0,y_cr]
    plots,[l1,x_cr],[0,y_cr]
    plots,[0,l1],[0,0]

    length=1000
    plots,[x_cr,x_cr-length*sin(phi+beta)],[y_cr,y_cr+length*cos(phi+beta)],color=200
    plots,[x_cr,x_cr-length*cos(phi+beta)],[y_cr,y_cr-length*sin(phi+beta)],color=200

    delta_hemount=[delta_l*cos(beta),delta_l*sin(beta),delta_z]
    u.info.position.tok.he_crystal=u.reset.position.tok.he_crystal-[[delta_hemount],[delta_hemount],[delta_hemount]]
    
    h_ori=u.info.position.h_ori
    
    reducer=u.info.position.tok.reducer    
    bport=u.info.position.tok.bport
    vessel=u.info.position.tok.vessel
    hemount=u.info.position.tok.he_crystal
    
    
    u.info.position.h.reducer=tran_position(reducer,h_ori,angle=phi-beta)
    u.info.position.h.bport=tran_position(bport,h_ori,angle=phi-beta)
    u.info.position.h.vessel=tran_position(vessel,h_ori,angle=phi-beta)
    u.info.position.h.he_crystal=tran_position(hemount,h_ori,angle=phi-beta)

end


PRO plot_raytracing,u

    if u.info.module eq 4 then plot_crystal_h,u else  plot_crystal_he,u
    plot_detector,u
    plot_vignetting,u

END



PRO HE2TREE,U

    shot=u.info.shot
    ;module=u.info.module

    l1=1414.15
    beta=8.00004*!pi/180.
    l=u.info.he_l_c
    phi=u.info.he_phi_c

    z=u.info.he_z_c/1000.
    R=sqrt(l1^2+l^2+2*l1*l*cos(beta))/1000.
    gama=!pi/2.+phi-beta

    RAD=[R,0,Z]
    ROT=[0,0,GAMA]

    BPORT=u.info.position.He.bport
    REDUCE=u.info.position.He.reducer
    VACVES=u.info.position.He.vessel

    path='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.RAYTRACING'
    mdsopen,'spectroscopy',shot
    mdsput,path+'.BPORT','build_with_units($,"mm")',BPORT
    mdsput,path+'.REDUCE','build_with_units($,"mm")',VACVES
    mdsput,path+'.VACVES','build_with_units($,"mm")',BPORT

    for module =1,3 do begin
        path='\SPECTROSCOPY::TOP.HIREXSR.INFO.MOD'+num2str(module,1)
        mdsput,path+'MIRROR.RAD','build_with_units($,"m")',RAD
        mdsput,path+'MIRROR.ROT','build_with_units($,"radians")',ROT
    endfor

    mdsclose,'spectroscopy',shot
END



PRO H2TREE,U
    shot=u.info.shot
    module=4

    l1=1414.15
    beta=8.00004*!pi/180.
    l=u.info.he_l_c
    phi=u.info.he_phi_c

    z=u.info.he_z_c/1000.
    R=sqrt(l1^2+l^2+2*l1*l*cos(beta))/1000.
    gama=!pi/2.+phi-beta

    RAD=[R,0,Z]
    ROT=[0,0,GAMA]

    BPORT=u.info.position.H.bport
    REDUCE=u.info.position.H.reducer
    VACVES=u.info.position.H.vessel
    HECRYS=u.info.position.H.he_crystal

    path='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.RAYTRACING'
    mdsopen,'spectroscopy',shot
    mdsput,path+'.BPORT','build_with_units($,"mm")',BPORT
    mdsput,path+'.REDUCE','build_with_units($,"mm")',REDUCE
    mdsput,path+'.VACVES','build_with_units($,"mm")',VACVES
    mdsput,path+'.HECRYS','build_with_units($,"mm")',HECRYS
        
    path='\SPECTROSCOPY::TOP.HIREXSR.INFO.MOD'+num2str(module,1)
    mdsput,path+'MIRROR.RAD','build_with_units($,"m")',RAD
    mdsput,path+'MIRROR.ROT','build_with_units($,"radians")',ROT    
    
    mdsclose,'spectroscopy',shot
END


PRO w_hirexsr_raytracing_event,event

        widget_control,event.top,get_uvalue=u
        id = u.id
        tag = tag_names(event,/st)        
	button=' '
	idtags=tag_names(id)
        FOR i=0,n(idtags) DO IF id.(i) EQ event.id THEN ename=idtags[i]
	CASE tag OF
		"WIDGET_BUTTON": BEGIN
			widget_control,event.id,get_value=button,get_uvalue=uvalue
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF
				
	 			"QUIT":  widget_control,event.top,/destroy
                                "APPLY1": plot_raytracing,u                                
                                "APPLY2": plot_raytracing,u
                                "APPLY3": plot_raytracing,u
                                "APPLY4": plot_raytracing,u

                                "RESET1": BEGIN
                                            reset1_u,u
                                            plot_raytracing,u
                                        END
                                "RESET2": BEGIN
                                            reset2_u,u
                                            plot_raytracing,u
                                        END
                                "RESET3": BEGIN
                                            reset3_u,u
                                            plot_raytracing,u
                                        END  
                                "RESET4": BEGIN
                                            reset4_u,u
                                            plot_raytracing,u
                                        END      
                                "RESET5": BEGIN
                                            reset5_u,u
                                            plot_raytracing,u
                                        END 
                                "RESET6": BEGIN
                                            reset6_u,u
                                            plot_raytracing,u
                                        END    
                                "WRI_HE": BEGIN
                                            he2tree,u
                                        END   
                                "WRI_H": BEGIN
                                            h2tree,u
                                        END   
                         
				ELSE:
                            ENDCASE
                    END

  		"WIDGET_SLIDER": BEGIN
                    widget_control,event.id,get_value=slider,get_uvalue=uvalue
                    IF NOT keyword_set(uvalue) THEN uvalue='none'
                    *u.image_index=0
                    CASE ename OF                                       
                        'XI_SLIDER'   : BEGIN
                            plot_raytracing,u
                            widget_control,u.id.xi_slider,get_value=xi
                            widget_control,u.id.zeta_slider,get_value=zeta
                            tmp=*u.image
                            widget_control,u.id.counts,set_value=num2str(tmp(xi,zeta))
                        END
                        'ZETA_SLIDER' : BEGIN
                            plot_raytracing,u         
                            widget_control,u.id.xi_slider,get_value=xi
                            widget_control,u.id.zeta_slider,get_value=zeta
                            tmp=*u.image
                            widget_control,u.id.counts,set_value=num2str(tmp(xi,zeta))              
                        END   
                        'HE_L_SLIDER' : plot_raytracing,u
                        'HE_PHI_SLIDER': plot_raytracing,u
                        'HE_Z_SLIDER' : plot_raytracing,u
                        'H_L_SLIDER' : plot_raytracing,u
                        'H_PHI_SLIDER': plot_raytracing,u
                        'H_Z_SLIDER' : plot_raytracing,u
                        ELSE:
                    ENDCASE
                END
                "WIDGET_TEXT_CH": BEGIN
                    CASE event.id OF
                        u.id.he_ny : BEGIN
                            widget_control,u.id.he_ny,get_value=he_ny
                            u.info.he_ny=he_ny
                        END
                        u.id.he_nz : BEGIN
                                widget_control,u.id.he_nz,get_value=he_nz
                                u.info.he_nz=he_nz
                            END
                        u.id.h_ny : BEGIN
                            widget_control,u.id.h_ny,get_value=h_ny
                            u.info.h_ny=h_ny
                        END
                        u.id.h_nz : BEGIN
                            widget_control,u.id.h_nz,get_value=h_nz
                            u.info.h_nz=h_nz
                        END
                        
                        u.id.shotid : BEGIN
                            widget_control,u.id.shotid,get_value=shot
                            shot=long(shot(0))
                            print,'shot:'+strtrim(shot,2)
                            u.info.shot=shot
                            u.load_info=hirexsr_load_info(u.info.module,shot=shot)
                            hirexsr_load_lambda,shot,u.info.module,lambda ;lambda=[487,195]   stored format in tree
                            u.lambda=reverse(reverse(transpose(lambda),2)) ;u.lambda=[195,487] in detector coordinate
                            *u.image_index=1
                        END
                        u.id.module : BEGIN
                            widget_control,u.id.module,get_value=module
                            module=int(module(0))
                            print,'module:'+strtrim(module,2)
                            u.info.module=module
                            u.load_info=hirexsr_load_info(module,shot=u.info.shot)
                            hirexsr_load_lambda,u.info.shot,module,lambda
                            u.lambda=reverse(reverse(transpose(lambda),2)) 
                            *u.image_index=1
                        END
                            u.id.t1 : BEGIN
                                widget_control,u.id.t1,get_value=t1
                                t1=float(t1(0))                                        
                                print,'t1:'+strtrim(t1,2)
                                u.info.t1=t1
                                *u.image_index=1
                            END
                            u.id.t2 : BEGIN
                                widget_control,u.id.t2,get_value=t2
                                t2=float(t2(0))
                                print,'t2:'+strtrim(t2,2)
                                u.info.t2=t2
                                *u.image_index=1
                            END
;==================================================
                            
                            u.id.he_re_x1 : BEGIN
                                widget_control,u.id.he_re_x1,get_value=x1
                                u.info.position.he.reducer(0)=x1
                            END
                            
                            u.id.he_re_y1 : BEGIN
                                widget_control,u.id.he_re_y1,get_value=y1
                                u.info.position.he.reducer(1)=y1
                            END
                            
                            u.id.he_re_z1 : BEGIN
                                widget_control,u.id.he_re_z1,get_value=z1
                                u.info.position.he.reducer(2)=z1
                            END
                            
                            u.id.he_re_x2 : BEGIN
                                widget_control,u.id.he_re_x2,get_value=x2
                                u.info.position.he.reducer(3)=x2
                            END
                            
                            u.id.he_re_y2 : BEGIN
                                widget_control,u.id.he_re_y2,get_value=y2
                                u.info.position.he.reducer(4)=y2
                            END
                            
                            u.id.he_re_z2 : BEGIN
                                widget_control,u.id.he_re_z2,get_value=z2
                                u.info.position.he.reducer(5)=z2
                            END
                            
                            u.id.he_re_x3 : BEGIN
                                widget_control,u.id.he_re_x3,get_value=x3
                                u.info.position.he.reducer(6)=x3
                            END
                            
                            
                            u.id.he_re_y3 : BEGIN
                                widget_control,u.id.he_re_y3,get_value=y3
                                u.info.position.he.reducer(7)=y3
                            END
                            

                            u.id.he_re_z3 : BEGIN
                                widget_control,u.id.he_re_z3,get_value=z3
                                u.info.position.he.reducer(8)=z3
                            END
                            
                            

                            u.id.he_bp_x1 : BEGIN
                                widget_control,u.id.he_bp_x1,get_value=x1
                                u.info.position.he.bport(0)=x1
                            END
                            

                            u.id.he_bp_y1 : BEGIN
                                widget_control,u.id.he_bp_y1,get_value=y1
                                u.info.position.he.bport(1)=y1
                                print,y1
                            END


                            u.id.he_bp_z1 : BEGIN
                                widget_control,u.id.he_bp_z1,get_value=z1
                                u.info.position.he.bport(2)=z1
                            END

                            
                            u.id.he_bp_x2 : BEGIN
                                widget_control,u.id.he_bp_x2,get_value=x2
                                u.info.position.he.bport(3)=x2
                            END

                            
                            u.id.he_bp_y2 : BEGIN
                                widget_control,u.id.he_bp_y2,get_value=y2
                                u.info.position.he.bport(4)=y2
                            END
                                    

                            u.id.he_bp_z2 : BEGIN
                                widget_control,u.id.he_bp_z2,get_value=z2
                                u.info.position.he.bport(5)=z2
                            END
                            

                                u.id.he_bp_x3 : BEGIN	
                                    	widget_control,u.id.he_bp_x3,get_value=x3
                                        u.info.position.he.bport(6)=x3
                                    END


                                u.id.he_bp_y3 : BEGIN
                                    	widget_control,u.id.he_bp_y3,get_value=y3
                                        u.info.position.he.bport(7)=y3
                                    END


                                u.id.he_bp_z3 : BEGIN
                                    	widget_control,u.id.he_bp_z3,get_value=z3
                                        u.info.position.he.bport(8)=z3
                                    END



                                u.id.he_ve_x1 : BEGIN
                                    	widget_control,u.id.he_ve_x1,get_value=x1
                                        u.info.position.he.vessel(0)=x1
                                    END


                                u.id.he_ve_y1 : BEGIN
                                    	widget_control,u.id.he_ve_y1,get_value=y1
                                        u.info.position.he.vessel(1)=y1
                                    END


                                u.id.he_ve_z1 : BEGIN
                                    	widget_control,u.id.he_ve_z1,get_value=z1
                                        u.info.position.he.vessel(2)=z1
                                    END


                                u.id.he_ve_x2 : BEGIN
                                    	widget_control,u.id.he_ve_x2,get_value=x2
                                        u.info.position.he.vessel(3)=x2
                                    END


                                u.id.he_ve_y2 : BEGIN
                                    	widget_control,u.id.he_ve_y2,get_value=y2
                                        u.info.position.he.vessel(4)=y2
                                    END


                                u.id.he_ve_z2 : BEGIN
                                    	widget_control,u.id.he_ve_z2,get_value=z2
                                        u.info.position.he.vessel(5)=z2
                                    END
   
                                u.id.he_ve_x3 : BEGIN
                                    	widget_control,u.id.he_ve_x3,get_value=x3
                                        u.info.position.he.vessel(6)=x3
                                    END


                                u.id.he_ve_y3 : BEGIN
                                    	widget_control,u.id.he_ve_y3,get_value=y3
                                        u.info.position.he.vessel(7)=y3
                                    END


                                u.id.he_ve_z3 : BEGIN
                                    	widget_control,u.id.he_ve_z3,get_value=z3
                                        u.info.position.he.vessel(8)=z3
                                    END

;=======================================================

                                u.id.h_re_x1 : BEGIN
                                	widget_control,u.id.h_re_x1,get_value=x1
                                        u.info.position.h.reducer(0)=x1
                                    END


                                u.id.h_re_y1 : BEGIN
                                	widget_control,u.id.h_re_y1,get_value=y1
                                        u.info.position.h.reducer(1)=y1
                                    END


                                u.id.h_re_z1 : BEGIN
                                	widget_control,u.id.h_re_z1,get_value=z1
                                        u.info.position.h.reducer(2)=z1
                                    END


                                u.id.h_re_x2 : BEGIN
                                	widget_control,u.id.h_re_x2,get_value=x2
                                        u.info.position.h.reducer(3)=x2
                                    END


                                u.id.h_re_y2 : BEGIN
                                	widget_control,u.id.h_re_y2,get_value=y2
                                        u.info.position.h.reducer(4)=y2
                                    END


                                u.id.h_re_z2 : BEGIN
                                	widget_control,u.id.h_re_z2,get_value=z2
                                        u.info.position.h.reducer(5)=z2
                                    END


                                u.id.h_re_x3 : BEGIN
                                	widget_control,u.id.h_re_x3,get_value=x3
                                        u.info.position.h.reducer(6)=x3
                                    END


                                u.id.h_re_y3 : BEGIN
                                	widget_control,u.id.h_re_y3,get_value=y3
                                        u.info.position.h.reducer(7)=y3
                                    END


                                u.id.h_re_z3 : BEGIN
                                	widget_control,u.id.h_re_z3,get_value=z3
                                        u.info.position.h.reducer(8)=z3
                                    END

                                                                                                                     
                                u.id.h_bp_x1 : BEGIN
                                	widget_control,u.id.h_bp_x1,get_value=x1
                                        u.info.position.h.bport(0)=x1
                                    END


                                u.id.h_bp_y1 : BEGIN
                                	widget_control,u.id.h_bp_y1,get_value=y1
                                        u.info.position.h.bport(1)=y1
                                    END


                                u.id.h_bp_z1 : BEGIN
                                	widget_control,u.id.h_bp_z1,get_value=z1
                                        u.info.position.h.bport(2)=z1
                                    END


                                u.id.h_bp_x2 : BEGIN
                                	widget_control,u.id.h_bp_x2,get_value=x2
                                        u.info.position.h.bport(3)=x2
                                    END


                                u.id.h_bp_y2 : BEGIN
                                	widget_control,u.id.h_bp_y2,get_value=y2
                                        u.info.position.h.bport(4)=y2
                                    END


                                u.id.h_bp_z2 : BEGIN
                                	widget_control,u.id.h_bp_z2,get_value=z2
                                        u.info.position.h.bport(5)=z2
                                    END


                                u.id.h_bp_x3 : BEGIN
                                	widget_control,u.id.h_bp_x3,get_value=x3
                                        u.info.position.h.bport(6)=x3
                                    END

                                u.id.h_bp_y3 : BEGIN
                                	widget_control,u.id.h_bp_y3,get_value=y3
                                        u.info.position.h.bport(7)=y3
                                    END


                                u.id.h_bp_z3 : BEGIN
                                	widget_control,u.id.h_bp_z3,get_value=z3
                                        u.info.position.h.bport(8)=z3
                                    END
                                        
                                u.id.h_ve_x1 : BEGIN
                                	widget_control,u.id.h_ve_x1,get_value=x1
                                        u.info.position.h.vessel(0)=x1
                                    END

                                        
                                u.id.h_ve_y1 : BEGIN
                                	widget_control,u.id.h_ve_y1,get_value=y1
                                        u.info.position.h.vessel(1)=y1
                                    END

                                        
                                u.id.h_ve_z1 : BEGIN
                                	widget_control,u.id.h_ve_z1,get_value=z1
                                        u.info.position.h.vessel(2)=z1
                                    END

                                        
                                u.id.h_ve_x2 : BEGIN
                                	widget_control,u.id.h_ve_x2,get_value=x2
                                        u.info.position.h.vessel(3)=x2
                                    END

                                        
                                u.id.h_ve_y2 : BEGIN
                                	widget_control,u.id.h_ve_y2,get_value=y2
                                        u.info.position.h.vessel(4)=y2
                                    END

                                        
                                u.id.h_ve_z2 : BEGIN
                                	widget_control,u.id.h_ve_z2,get_value=z2
                                        u.info.position.h.vessel(5)=z2
                                    END

                                        
                                u.id.h_ve_x3 : BEGIN
                                	widget_control,u.id.h_ve_x3,get_value=x3
                                        u.info.position.h.vessel(6)=x3
                                    END


                                u.id.h_ve_y3 : BEGIN
                                	widget_control,u.id.h_ve_y3,get_value=y3
                                        u.info.position.h.vessel(7)=y3
                                    END

                                        
                                u.id.h_ve_z3 : BEGIN
                                	widget_control,u.id.h_ve_z3,get_value=z3
                                        u.info.position.h.vessel(8)=z3
                                    END

                                        
                                u.id.h_cr_x1 : BEGIN
                                	widget_control,u.id.h_cr_x1,get_value=x1
                                        u.info.position.h.he_crystal(0)=x1
                                    END

                                        
                                u.id.h_cr_y1 : BEGIN
                                	widget_control,u.id.h_cr_y1,get_value=y1
                                        u.info.position.h.he_crystal(1)=y1
                                    END

                                        
                                u.id.h_cr_z1 : BEGIN
                                	widget_control,u.id.h_cr_z1,get_value=z1
                                        u.info.position.h.he_crystal(2)=z1
                                    END

                                        
                                u.id.h_cr_x2 : BEGIN
                                	widget_control,u.id.h_cr_x2,get_value=x2
                                        u.info.position.h.he_crystal(3)=x2
                                    END

                                        
                                u.id.h_cr_y2 : BEGIN
                                	widget_control,u.id.h_cr_y2,get_value=y2
                                        u.info.position.h.he_crystal(4)=y2
                                    END

                                        
                                u.id.h_cr_z2 : BEGIN
                                	widget_control,u.id.h_cr_z2,get_value=z2
                                        u.info.position.h.he_crystal(5)=z2
                                    END

                                        
                                u.id.h_cr_x3 : BEGIN
                                	widget_control,u.id.h_cr_x3,get_value=x3
                                        u.info.position.h.he_crystal(6)=x3
                                    END

                                        
                                u.id.h_cr_y3 : BEGIN
                                	widget_control,u.id.h_cr_y3,get_value=y3
                                        u.info.position.h.he_crystal(7)=y3
                                    END

                                        
                                u.id.h_cr_z3 : BEGIN
                                	widget_control,u.id.h_cr_z3,get_value=z3
                                        u.info.position.h.he_crystal(8)=z3
                                    END


				ELSE: 
                                ENDCASE
                            END
		ELSE:
                        ENDCASE
           
	IF button NE 'Quit' THEN widget_control,event.top,set_uvalue=u

END






PRO w_hirexsr_raytracing
	user=logname()
	
	loadct,12,/silent
	base=widget_base(title='HIREXSR RAY TRACING',/row)
	left=widget_base(base,/column)
	center=widget_base(base,/column)
	right=widget_base(base,/column)
        rright=widget_base(base,/column)
        rrright=widget_base(base,/column)

;==============
	;setup left windows
        scale=1.5
	l1=widget_base(left,frame=5,/row)
	l1p1=widget_base(l1,frame=5,/column)
 	dum = widget_label(l1p1,value='B Port Flange')       
        draw1=widget_draw(l1p1,xsize=150*scale,ysize=200*scale)
        dum  =widget_label(l1p1,value='Attentuation')
        atnu1=widget_text(l1p1,xsize=7)
	l1p2=widget_base(l1,frame=5,/column)
        dum = widget_label(l1p2,value='Vaccum Vessel')
	draw2=widget_draw(l1p2,xsize=150*scale,ysize=200*scale)	
        dum  =widget_label(l1p2,value='Attentuation')
        atnu2=widget_text(l1p2,xsize=7)
        l2=widget_base(left,frame=5,/column)
        dum = widget_label(l2,value='He-like Crystal Mount(Module 4 only)')        
        draw3=widget_draw(l2,xsize=487,ysize=195)
        l2p1 = widget_base(l2,frame=5,/row)
        dum  =widget_label(l2p1,value='Attentuation')
        atnu3=widget_text(l2p1,xsize=7)
        l3=widget_base(left,frame=5,/row)
        l3p1=widget_base(l3,/column)
        dum = widget_label(l3p1,value='Reducer Flange')        
	draw4=widget_draw(l3p1,xsize=487*0.7,ysize=487*0.7)
        l3p1p1 = widget_base(l3p1,frame=5,/row)
        dum  =widget_label(l3p1p1,value='Attentuation')
        atnu4=widget_text(l3p1p1,xsize=7)
        l3p2=widget_base(l3,/column)
        l3p2p1=widget_base(l3p2)
	dum = widget_label(l3p2p1,value='He-like Crystal Points: ')        
        l3p2p2=widget_base(l3p2,/row)
        dum = widget_label(l3p2p2,value='nY=')
        he_ny = widget_text(l3p2p2,xsize=4,/edit)
	dum = widget_label(l3p2p2,value=' nZ=')
  	he_nz = widget_text(l3p2p2,xsize=4,/edit)      
        l3p2p3=widget_base(l3p2)
	dum = widget_label(l3p2p3,value='H-like Crystal Points: ')        
        l3p2p4=widget_base(l3p2,/row)
        dum = widget_label(l3p2p4,value='nY=')
	h_ny = widget_text(l3p2p4,xsize=4,/edit)
	dum = widget_label(l3p2p4,value=' nZ=')
	h_nz = widget_text(l3p2p4,xsize=4,/edit)  
        l3p2p5=widget_base(l3p2,/row)
        Apply1 = widget_button(l3p2p5,value='Apply')
        Reset1= widget_button(l3p2p5,value='Reset')   

;=========================     
	;setup center windows
        
	zeta_slider=widget_slider(center,ysize=487*1.,min=0,max=486,value=487/2.,/vert,/drag)

;=========================
	;setup right windows
        r1=widget_base(right,frame=5,/column)
        draw5=widget_draw(r1,xsize=195*1.5,ysize=487*1.)
        dum = widget_label(r1,value='Detector')
	xi_slider=widget_slider(r1,xsize=195*1.5,min=0,max=194,value=195/2.,/drag)
;=========================
	;setup rright base

        w_tab=widget_tab(right,location=0)

        rr1=widget_base(w_tab,title='Shot Info',/column)
        rr1p1=widget_base(rr1,/row)
        dum = widget_label(rr1p1,value='shot:')
	shotid = widget_text(rr1p1,xsize=12,/edit)        
        dum = widget_label(rr1p1,value='module')
	module = widget_text(rr1p1,xsize=7,/edit) 
        rr1p2=widget_base(rr1,/row)        
        dum = widget_label(rr1p2,value='t1(s)')
	t1 = widget_text(rr1p2,xsize=7,/edit)        
        dum = widget_label(rr1p2,value='t2(s)')
	t2 = widget_text(rr1p2,xsize=7,/edit) 
        rr1p3=widget_base(rr1,/row)        
        Apply2 = widget_button(rr1p3,value='Apply')
        Reset2= widget_button(rr1p3,value='Reset')
        Quit = widget_button(rr1p3,value='Quit')
                
        rr1p9=widget_base(rr1,/row)
        dum = widget_label(rr1p9,value='counts=')
	counts = widget_text(rr1p9,xsize=7)  
        rr1p4=widget_base(rr1,/row)        
        dum = widget_label(rr1p4,value='transmission=')
	trans = widget_text(rr1p4,xsize=7)      
        rr1p5=widget_base(rr1,/row)       
        dum = widget_label(rr1p5,value='average angle=')
	ave_angle = widget_text(rr1p5,xsize=7) 
        rr1p6=widget_base(rr1,/row)       
        dum = widget_label(rr1p6,value='maximum angle=')
	max_angle = widget_text(rr1p6,xsize=7) 
        rr1p7=widget_base(rr1,/row)       
        dum = widget_label(rr1p7,value='total attenuation=')
	total_atnu = widget_text(rr1p7,xsize=7)
        rr1p8=widget_base(rr1,/row)       
        dum = widget_label(rr1p8,value='total trans=')
	total_trans = widget_text(rr1p8,xsize=7)

;==========================================

        rr2=widget_base(w_tab,title='He Vector',/column)       
	dum = widget_label(rr2,value='He like Positions')        
        rr2p1  =widget_base(rr2,/row)
        dum      = widget_label(rr2p1,value='Reducer: ')
        dum      = widget_label(rr2p1,value='x1=') 
	he_re_x1 = widget_text(rr2p1,xsize=7,/edit)      
        dum      = widget_label(rr2p1,value='y1=')
	he_re_y1 = widget_text(rr2p1,xsize=7,/edit)       
        dum      = widget_label(rr2p1,value='z1=')
	he_re_z1 = widget_text(rr2p1,xsize=7,/edit)

	rr2p2 =widget_base(rr2,/row)
        dum     = widget_label(rr2p2,value='         ')
        dum     = widget_label(rr2p2,value='x2=') 
	he_re_x2= widget_text(rr2p2,xsize=7,/edit)      
        dum     = widget_label(rr2p2,value='y2=')
	he_re_y2= widget_text(rr2p2,xsize=7,/edit)       
        dum     = widget_label(rr2p2,value='z2=')
	he_re_z2= widget_text(rr2p2,xsize=7,/edit)

	rr2p3 =widget_base(rr2,/row)
        dum     = widget_label(rr2p3,value='         ')
        dum     = widget_label(rr2p3,value='x3=') 
	he_re_x3= widget_text(rr2p3,xsize=7,/edit)      
        dum     = widget_label(rr2p3,value='y3=')
	he_re_y3= widget_text(rr2p3,xsize=7,/edit)       
        dum     = widget_label(rr2p3,value='z3=')
	he_re_z3= widget_text(rr2p3,xsize=7,/edit)       
       
        rr2p4 =widget_base(rr2,/row)
        dum     = widget_label(rr2p4,value='B Port : ')
        dum     = widget_label(rr2p4,value='x1=') 
	he_bp_x1= widget_text(rr2p4,xsize=7,/edit)      
        dum     = widget_label(rr2p4,value='y1=')
	he_bp_y1= widget_text(rr2p4,xsize=7,/edit)       
        dum     = widget_label(rr2p4,value='z1=')
	he_bp_z1= widget_text(rr2p4,xsize=7,/edit)
	
        rr2p5 =widget_base(rr2,/row)
        dum     = widget_label(rr2p5,value='         ')
        dum     = widget_label(rr2p5,value='x2=') 
	he_bp_x2= widget_text(rr2p5,xsize=7,/edit)      
        dum     = widget_label(rr2p5,value='y2=')
	he_bp_y2= widget_text(rr2p5,xsize=7,/edit)       
        dum     = widget_label(rr2p5,value='z2=')
	he_bp_z2= widget_text(rr2p5,xsize=7,/edit)
	
        rr2p6 =widget_base(rr2,/row)
        dum     = widget_label(rr2p6,value='         ')
        dum     = widget_label(rr2p6,value='x3=') 
	he_bp_x3= widget_text(rr2p6,xsize=7,/edit)      
        dum     = widget_label(rr2p6,value='y3=')
	he_bp_y3= widget_text(rr2p6,xsize=7,/edit)       
        dum     = widget_label(rr2p6,value='z3=')
	he_bp_z3= widget_text(rr2p6,xsize=7,/edit)               

        rr2p7 =widget_base(rr2,/row)
        dum     = widget_label(rr2p7,value='Vessel : ')
        dum     = widget_label(rr2p7,value='x1=') 
	he_ve_x1= widget_text(rr2p7,xsize=7,/edit)      
        dum     = widget_label(rr2p7,value='y1=')
	he_ve_y1= widget_text(rr2p7,xsize=7,/edit)       
        dum     = widget_label(rr2p7,value='z1=')
	he_ve_z1= widget_text(rr2p7,xsize=7,/edit)

	rr2p8 =widget_base(rr2,/row)
        dum     = widget_label(rr2p8,value='         ')
        dum     = widget_label(rr2p8,value='x2=') 
	he_ve_x2= widget_text(rr2p8,xsize=7,/edit)      
        dum     = widget_label(rr2p8,value='y2=')
	he_ve_y2= widget_text(rr2p8,xsize=7,/edit)       
        dum     = widget_label(rr2p8,value='z2=')
	he_ve_z2= widget_text(rr2p8,xsize=7,/edit)

	rr2p9 =widget_base(rr2,/row)
        dum     = widget_label(rr2p9,value='         ')
        dum     = widget_label(rr2p9,value='x3=') 
	he_ve_x3= widget_text(rr2p9,xsize=7,/edit)      
        dum     = widget_label(rr2p9,value='y3=')
	he_ve_y3= widget_text(rr2p9,xsize=7,/edit)       
        dum     = widget_label(rr2p9,value='z3=')
	he_ve_z3= widget_text(rr2p9,xsize=7,/edit)  

	rr2p10 =widget_base(rr2,/row)
        Apply3  = widget_button(rr2p10,value='Apply')
        Reset3  = widget_button(rr2p10,value='Reset')       
        wri_he= widget_button(rr2p10,value='Write to He Tree')

;===============================
        rr3=widget_base(w_tab,title='H Vector',/column)
	dum = widget_label(rr3,value='H like Positions')        
	
        rr3p1 = widget_base(rr3,/row)
        dum     = widget_label(rr3p1,value='Reducer: ')
        dum     = widget_label(rr3p1,value='x1=') 
	h_re_x1 = widget_text(rr3p1,xsize=7,/edit)      
        dum     = widget_label(rr3p1,value='y1=')
	h_re_y1 = widget_text(rr3p1,xsize=7,/edit)       
        dum     = widget_label(rr3p1,value='z1=')
	h_re_z1 = widget_text(rr3p1,xsize=7,/edit)

	rr3p2 = widget_base(rr3,/row)
        dum     = widget_label(rr3p2,value='         ')
        dum     = widget_label(rr3p2,value='x2=') 
	h_re_x2 = widget_text(rr3p2,xsize=7,/edit)      
        dum     = widget_label(rr3p2,value='y2=')
	h_re_y2 = widget_text(rr3p2,xsize=7,/edit)       
        dum     = widget_label(rr3p2,value='z2=')
	h_re_z2 = widget_text(rr3p2,xsize=7,/edit)

	rr3p3 =widget_base(rr3,/row)
        dum     = widget_label(rr3p3,value='         ')
        dum     = widget_label(rr3p3,value='x3=') 
	h_re_x3 = widget_text(rr3p3,xsize=7,/edit)      
        dum     = widget_label(rr3p3,value='y3=')
	h_re_y3 = widget_text(rr3p3,xsize=7,/edit)       
        dum     = widget_label(rr3p3,value='z3=')
	h_re_z3 = widget_text(rr3p3,xsize=7,/edit)       
       
        rr3p4 =widget_base(rr3,/row)
        dum     = widget_label(rr3p4,value='B Port : ')
        dum     = widget_label(rr3p4,value='x1=') 
	h_bp_x1 = widget_text(rr3p4,xsize=7,/edit)      
        dum     = widget_label(rr3p4,value='y1=')
	h_bp_y1 = widget_text(rr3p4,xsize=7,/edit)       
        dum     = widget_label(rr3p4,value='z1=')
	h_bp_z1 = widget_text(rr3p4,xsize=7,/edit)
	
        rr3p5 =widget_base(rr3,/row)
        dum     = widget_label(rr3p5,value='         ')
        dum     = widget_label(rr3p5,value='x2=') 
	h_bp_x2 = widget_text(rr3p5,xsize=7,/edit)      
        dum     = widget_label(rr3p5,value='y2=')
	h_bp_y2 = widget_text(rr3p5,xsize=7,/edit)       
        dum     = widget_label(rr3p5,value='z2=')
	h_bp_z2 = widget_text(rr3p5,xsize=7,/edit)
	
        rr3p6 =widget_base(rr3,/row)
        dum     = widget_label(rr3p6,value='         ')
        dum     = widget_label(rr3p6,value='x3=') 
	h_bp_x3 = widget_text(rr3p6,xsize=7,/edit)      
        dum     = widget_label(rr3p6,value='y3=')
	h_bp_y3 = widget_text(rr3p6,xsize=7,/edit)       
        dum     = widget_label(rr3p6,value='z3=')
	h_bp_z3 = widget_text(rr3p6,xsize=7,/edit)               

        rr3p7 =widget_base(rr3,/row)
        dum     = widget_label(rr3p7,value='Vessel : ')
        dum     = widget_label(rr3p7,value='x1=') 
	h_ve_x1 = widget_text(rr3p7,xsize=7,/edit)      
        dum     = widget_label(rr3p7,value='y1=')
	h_ve_y1 = widget_text(rr3p7,xsize=7,/edit)       
        dum     = widget_label(rr3p7,value='z1=')
	h_ve_z1 = widget_text(rr3p7,xsize=7,/edit)

	rr3p8 =widget_base(rr3,/row)
        dum     = widget_label(rr3p8,value='         ')
        dum     = widget_label(rr3p8,value='x2=') 
	h_ve_x2 = widget_text(rr3p8,xsize=7,/edit)      
        dum     = widget_label(rr3p8,value='y2=')
	h_ve_y2 = widget_text(rr3p8,xsize=7,/edit)       
        dum     = widget_label(rr3p8,value='z2=')
	h_ve_z2 = widget_text(rr3p8,xsize=7,/edit)

	rr3p9 =widget_base(rr3,/row)
        dum     = widget_label(rr3p9,value='         ')
        dum     = widget_label(rr3p9,value='x3=') 
	h_ve_x3 = widget_text(rr3p9,xsize=7,/edit)      
        dum     = widget_label(rr3p9,value='y3=')
	h_ve_y3 = widget_text(rr3p9,xsize=7,/edit)       
        dum     = widget_label(rr3p9,value='z3=')
	h_ve_z3 = widget_text(rr3p9,xsize=7,/edit)  

        rr3p10=widget_base(rr3,/row)
        dum     = widget_label(rr3p10,value='He Crst: ')
        dum     = widget_label(rr3p10,value='x1=') 
	h_cr_x1 = widget_text(rr3p10,xsize=7,/edit)      
        dum     = widget_label(rr3p10,value='y1=')
	h_cr_y1 = widget_text(rr3p10,xsize=7,/edit)       
        dum     = widget_label(rr3p10,value='z1=')
	h_cr_z1 = widget_text(rr3p10,xsize=7,/edit)

	rr3p11=widget_base(rr3,/row)
        dum     = widget_label(rr3p11,value='         ')
        dum     = widget_label(rr3p11,value='x2=') 
	h_cr_x2 = widget_text(rr3p11,xsize=7,/edit)      
        dum     = widget_label(rr3p11,value='y2=')
	h_cr_y2 = widget_text(rr3p11,xsize=7,/edit)       
        dum     = widget_label(rr3p11,value='z2=')
	h_cr_z2 = widget_text(rr3p11,xsize=7,/edit)

	rr3p12=widget_base(rr3,/row)
        dum     = widget_label(rr3p12,value='         ')
        dum     = widget_label(rr3p12,value='x3=') 
	h_cr_x3 = widget_text(rr3p12,xsize=7,/edit)      
        dum     = widget_label(rr3p12,value='y3=')
	h_cr_y3 = widget_text(rr3p12,xsize=7,/edit)       
        dum     = widget_label(rr3p12,value='z3=')
	h_cr_z3 = widget_text(rr3p12,xsize=7,/edit)

        rr3p13=widget_base(rr3,/row)       
        Apply4  = widget_button(rr3p13,value='Apply')
        Reset4  = widget_button(rr3p13,value='Reset')        
        wri_h = widget_button(rr3p13,value='Write to H Tree')        


;=========================
	;setup rrright base
        w_tab=widget_tab(rrright,location=0)
        
        rrr1=widget_base(w_tab,title='He Crystal Layout',/column)

        
        draw6=widget_draw(rrr1,xsize=500,ysize=500)
        rrr1p1=widget_base(rrr1,frame=5,/column)

        rrr1p1p1=widget_base(rrr1p1,frame=5,/row)
        dum =widget_label(rrr1p1p1,value=' L :')        
        he_l_crystal=widget_text(rrr1p1p1,xsize=6)
        he_l_slider=widget_slider(rrr1p1p1,xsize=400,min=-2000,max=2000,value=0,/drag,/suppress_value)

        rrr1p1p2=widget_base(rrr1p1,frame=5,/row)
        dum =widget_label(rrr1p1p2,value='phi:')        
        he_phi_crystal=widget_text(rrr1p1p2,xsize=6)
        he_phi_slider=widget_slider(rrr1p1p2,xsize=400,min=-2000,max=2000,value=0,/drag,/suppress)

        rrr1p1p3=widget_base(rrr1p1,frame=5,/row)
        dum =widget_label(rrr1p1p3,value=' Z :')        
        he_z_crystal=widget_text(rrr1p1p3,xsize=6)
        he_z_slider=widget_slider(rrr1p1p3,xsize=400,min=-2000,max=2000,value=0,/drag,/suppress)

        rrr1p1p4=widget_base(rrr1p1,/row)
        Reset5= widget_button(rrr1p1p4,value='Reset')


;=============
        rrr2=widget_base(w_tab,title='H Crystal Layout',/column)

        
        draw7=widget_draw(rrr2,xsize=500,ysize=500)
        rrr2p1=widget_base(rrr2,frame=5,/column)

        rrr2p1p1=widget_base(rrr2p1,frame=5,/row)
        dum =widget_label(rrr2p1p1,value=' L :')        
        h_l_crystal=widget_text(rrr2p1p1,xsize=6)
        h_l_slider=widget_slider(rrr2p1p1,xsize=400,min=-2000,max=2000,value=0,/drag,/suppress)

        rrr2p1p2=widget_base(rrr2p1,frame=5,/row)
        dum =widget_label(rrr2p1p2,value='phi:')        
        h_phi_crystal=widget_text(rrr2p1p2,xsize=6)
        h_phi_slider=widget_slider(rrr2p1p2,xsize=400,min=-2000,max=2000,value=0,/drag,/suppress)

        rrr2p1p3=widget_base(rrr2p1,frame=5,/row)
        dum =widget_label(rrr2p1p3,value=' Z :')        
        h_z_crystal=widget_text(rrr2p1p3,xsize=6)
        h_z_slider=widget_slider(rrr2p1p3,xsize=400,min=-2000,max=2000,value=0,/drag,/suppress)

        rrr2p1p4=widget_base(rrr2p1,/row)
        Reset6= widget_button(rrr2p1p4,value='Reset')


	;build u structure
	id={base:base,draw1:draw1,draw2:draw2,draw3:draw3,draw4:draw4,draw5:draw5,draw6:draw6,draw7:draw7,$  ;,draw6:draw6
            xi_slider:xi_slider,zeta_slider:zeta_slider, $  ;,y_slider:y_slider,z_slider:z_slider
            he_re_x1:he_re_x1,he_re_x2:he_re_x2,he_re_x3:he_re_x3,he_re_y1:he_re_y1,he_re_y2:he_re_y2,he_re_y3:he_re_y3,he_re_z1:he_re_z1,he_re_z2:he_re_z2,he_re_z3:he_re_z3,$
            he_bp_x1:he_bp_x1,he_bp_x2:he_bp_x2,he_bp_x3:he_bp_x3,he_bp_y1:he_bp_y1,he_bp_y2:he_bp_y2,he_bp_y3:he_bp_y3,he_bp_z1:he_bp_z1,he_bp_z2:he_bp_z2,he_bp_z3:he_bp_z3,$
            he_ve_x1:he_ve_x1,he_ve_x2:he_ve_x2,he_ve_x3:he_ve_x3,he_ve_y1:he_ve_y1,he_ve_y2:he_ve_y2,he_ve_y3:he_ve_y3,he_ve_z1:he_ve_z1,he_ve_z2:he_ve_z2,he_ve_z3:he_ve_z3,$
            h_re_x1:h_re_x1,h_re_x2:h_re_x2,h_re_x3:h_re_x3,h_re_y1:h_re_y1,h_re_y2:h_re_y2,h_re_y3:h_re_y3,h_re_z1:h_re_z1,h_re_z2:h_re_z2,h_re_z3:h_re_z3,$
            h_bp_x1:h_bp_x1,h_bp_x2:h_bp_x2,h_bp_x3:h_bp_x3,h_bp_y1:h_bp_y1,h_bp_y2:h_bp_y2,h_bp_y3:h_bp_y3,h_bp_z1:h_bp_z1,h_bp_z2:h_bp_z2,h_bp_z3:h_bp_z3,$
            h_ve_x1:h_ve_x1,h_ve_x2:h_ve_x2,h_ve_x3:h_ve_x3,h_ve_y1:h_ve_y1,h_ve_y2:h_ve_y2,h_ve_y3:h_ve_y3,h_ve_z1:h_ve_z1,h_ve_z2:h_ve_z2,h_ve_z3:h_ve_z3,$
            h_cr_x1:h_cr_x1,h_cr_x2:h_cr_x2,h_cr_x3:h_cr_x3,h_cr_y1:h_cr_y1,h_cr_y2:h_cr_y2,h_cr_y3:h_cr_y3,h_cr_z1:h_cr_z1,h_cr_z2:h_cr_z2,h_cr_z3:h_cr_z3,$
            Apply1:apply1, reset1:reset1, Apply2:apply2, reset2:reset2, Apply3:apply3, reset3:reset3,Apply4:apply4, reset4:reset4,Quit:Quit, reset5:reset5, reset6:reset6, $
            shotid:shotid, module:module, t1:t1,t2:t2,$
            he_ny:he_ny,he_nz:he_nz,h_ny:h_ny,h_nz:h_nz,$
            trans:trans,ave_angle:ave_angle,max_angle:max_angle,total_trans:total_trans,$
            atnu1:atnu1,atnu2:atnu2,atnu3:atnu3,atnu4:atnu4,total_atnu:total_atnu,$
            he_l_slider:he_l_slider,he_phi_slider:he_phi_slider,he_z_slider:he_z_slider,$
            he_l_crystal:he_l_crystal,he_phi_crystal:he_phi_crystal,he_z_crystal:he_z_crystal, $
            h_l_slider:h_l_slider,h_phi_slider:h_phi_slider,h_z_slider:h_z_slider,$
            h_l_crystal:h_l_crystal,h_phi_crystal:h_phi_crystal,h_z_crystal:h_z_crystal, $
            counts:counts, wri_he:wri_he, wri_h:wri_h $
        }
        
	;set defaults	
	shot=1070830020
	module=4
        load_info=hirexsr_load_info(module,shot=shot)
        hirexsr_load_lambda,shot,module,lambda  ;fltarr(487,195)
        lambda=reverse(reverse(transpose(lambda),2)) 
	
        t1=1.05000
	t2=1.45000
        image=fltarr(195,487)
        trans=fltarr(4,195,487)
        image_index=0
        atnu1=0.
        atnu2=0.
        atnu3=0.
        atnu4=0.       
        ;for Helike cyrstal, mirror points
        he_ny=50
        he_nz=50
        he_trans_info=fltarr(12,he_ny*he_nz)
        ;he_trans_info[0,*]: b port point x coord        
        ;he_trans_info[1,*]: b port point y coord
        ;he_trans_info[2,*]: b port point z coord
        ;he_trans_info[3,*]: flange point x coord        
        ;he_trans_info[4,*]: flange point y coord
        ;he_trans_info[5,*]: flange point z coord
        ;he_trans_info[6,*]: tan(angle)
        ;he_trans_info[7,*]: thick(h1)
        ;he_trans_info[8,*]: transmission(h1)
        ;he_trans_info[9,*]: thick(h2)
        ;he_trans_info[10,*]: transmission(h2)        
        ;he_trans_info[11,*]: 0/1.initial value is 1. 0 means this ray is vignetting/blocked       


        ;for H-like cyrstal, mirror points
        h_ny=80
        h_nz=80         
        h_trans_info=fltarr(11,h_ny*h_nz)

        trans_info={he:he_trans_info,h:h_trans_info}

        widget_control,id.xi_slider,get_value=xi
        widget_control,id.zeta_slider,get_value=zeta	
        
;=================

; get the position info 
;       <l,phi,z>-><X,Y,Z>
; position_tok is in Tokamak frame
; position_he is in He-Crystal frame
; position_h is in H-Crystal frame

; position_Tok
        he_ori=[3672.98,317.46,-6.43]
        h_ori=[3767.36,334.25,19.93]

        hemount=[[3691.84,360.86,8.54],[3656.91,272.24,8.54],[3691.84,360.86,-23.21]]
 	vpts=[[1061.72,0.,0.],[1061.72,101.6,0.],[1061.72,0.,314.32]]  
	bpts=[[1797.55,54.73,0.05],[1797.55,99.61,0.08],[1797.55,54.57,184.2]]	
	redpts=[[3455.92,287.80,0.24],[3448.85,338.10,0.24],[3455.92,287.79,51.04]]
        position_Tok={reducer:redpts,bport:bpts,vessel:vpts,he_crystal:hemount}


        xr=3455.92
        yr=287.8
        xb=1797.55
        yb=54.73
        beta=atan((yr-yb)/(xr-xb))  ;beta=8.00004 for both he and h-like

        he_xc=3672.98
        he_yc=317.46
        he_l1=he_xc-he_yc/tan(beta)   ;he_l1=1414.15
        he_l_0=he_yc/sin(beta)  ;he_l_0=2281.03
        he_phi_0=29.51*!pi/180.

        h_xc=3767.36
        h_yc=334.25
        h_l1=h_xc-h_yc/tan(beta)   ;he_l1=1389.06
        h_l_0=h_yc/sin(beta)  ;he_l_0=2401.67
        h_phi_0=34.89*!pi/180.

;above paras are all editable

    
; position_He
        angle=he_phi_0-beta
        vpts_he=tran_position(vpts,he_ori,angle=angle)
        bpts_he=tran_position(bpts,he_ori,angle=angle)
        redpts_he=tran_position(redpts,he_ori,angle=angle)
        position_He={reducer:redpts_he,bport:bpts_he,vessel:vpts_he}
        
; position_H
        angle=h_phi_0-beta
        vpts_h=tran_position(vpts,h_ori,angle=angle)
        bpts_h=tran_position(bpts,h_ori,angle=angle)
        redpts_h=tran_position(redpts,h_ori,angle=angle)
        hemount_h=tran_position(hemount,h_ori,angle=angle)
        position_H={reducer:redpts_h,bport:bpts_h,vessel:vpts_h,he_crystal:hemount_h}

        position={Tok:position_Tok,He:position_He,H:position_H,he_ori:he_ori,h_ori:h_ori}
        

;=================

        ;waverange=[3000:4000](mA)
        ;thickrange=[1:1.09]*50.8 um
        filter_thick_1=fltarr(10)
        filter_thick_2=fltarr(10)
        wavelength_1=fltarr(500)
        wavelength_2=fltarr(500)
        transmission_1=fltarr(10,500)
        transmission_2=fltarr(10,500)

        for i=0,9 do begin 
            path='/usr/local/cmod/idl/HIREXSR/be_data/filter_h1x1.0'+strtrim(i,2)+'.dat'
            tmp=read_xray_data(path)
            filter_thick_1(i)=tmp.thick
            transmission_1(i,*)=tmp.TR
        endfor
            wavelength_1=tmp.E*1.e4

        for i=0,9 do begin 
            path='/usr/local/cmod/idl/HIREXSR/be_data/filter_h2x1.0'+strtrim(i,2)+'.dat'
            tmp=read_xray_data(path)
            filter_thick_2(i)=tmp.thick
            transmission_2(i,*)=tmp.TR
        endfor
            wavelength_2=tmp.E*1.e4
            

        u_info={shot:shot,module:module,t1:t1,t2:t2,he_ny:he_ny,he_nz:he_nz,h_ny:h_ny,h_nz:h_nz,position:position,$
                filter_thick_1:filter_thick_1,wavelength_1:wavelength_1,transmission_1:transmission_1 ,delta_h1:0.01,h1:50.8, $
                filter_thick_2:filter_thick_2,wavelength_2:wavelength_2,transmission_2:transmission_2 ,delta_h2:0.01,h2:1066.8, $
                delta_lamd:2,lamd:3000,p_h1:0.659973,p_h2:0.340027,$
                he_l_c:he_l_0,he_phi_c:he_phi_0,he_z_c:he_ori(2),$
                h_l_c:h_l_0,h_phi_c:h_phi_0,h_z_c:h_ori(2), $
                he_l1:he_l1,h_l1:h_l1,beta:beta $
}
        u_reset={shot:shot,module:module,t1:t1,t2:t2,xi:xi,zeta:zeta,he_ny:he_ny,he_nz:he_nz,h_ny:h_ny,h_nz:h_nz,position:position,$
                 he_l_0:he_l_0,he_phi_0:he_phi_0,he_z_0:he_ori(2),$
                 h_l_0:h_l_0,h_phi_0:h_phi_0,h_z_0:h_ori(2) $
                }
        u={id:id,info:u_info,reset:u_reset,$
           image:ptr_new(image),image_index:ptr_new(image_index),transmission:ptr_new(trans),$
           load_info:load_info,lambda:lambda,atnu1:ptr_new(atnu1),atnu2:ptr_new(atnu2),atnu3:ptr_new(atnu3),atnu4:ptr_new(atnu4),$
           trans_info_h:ptr_new(h_trans_info),trans_info_he:ptr_new(he_trans_info)}


	widget_control,base,set_uvalue=u
	widget_control,u.id.shotid,set_value=strtrim(u.info.shot,2)
	widget_control,u.id.module,set_value=strtrim(u.info.module,2)
	widget_control,u.id.t1,set_value=num2str(u.info.t1,dp=2)
	widget_control,u.id.t2,set_value=num2str(u.info.t2,dp=2)
        
        widget_control,u.id.he_ny,set_value=num2str(u.info.he_ny)
        widget_control,u.id.he_nz,set_value=num2str(u.info.he_nz)
        widget_control,u.id.h_ny,set_value=num2str(u.info.h_ny)
        widget_control,u.id.h_nz,set_value=num2str(u.info.h_nz)

        widget_control,u.id.he_re_x1,set_value=num2str(u.info.position.he.reducer(0),dp=2)
        widget_control,u.id.he_re_y1,set_value=num2str(u.info.position.he.reducer(1),dp=2)
        widget_control,u.id.he_re_z1,set_value=num2str(u.info.position.he.reducer(2),dp=2)
        widget_control,u.id.he_re_x2,set_value=num2str(u.info.position.he.reducer(3),dp=2)
        widget_control,u.id.he_re_y2,set_value=num2str(u.info.position.he.reducer(4),dp=2)
        widget_control,u.id.he_re_z2,set_value=num2str(u.info.position.he.reducer(5),dp=2)
        widget_control,u.id.he_re_x3,set_value=num2str(u.info.position.he.reducer(6),dp=2)
        widget_control,u.id.he_re_y3,set_value=num2str(u.info.position.he.reducer(7),dp=2)
        widget_control,u.id.he_re_z3,set_value=num2str(u.info.position.he.reducer(8),dp=2)

        widget_control,u.id.he_bp_x1,set_value=num2str(u.info.position.he.bport(0),dp=2)
        widget_control,u.id.he_bp_y1,set_value=num2str(u.info.position.he.bport(1),dp=2)
        widget_control,u.id.he_bp_z1,set_value=num2str(u.info.position.he.bport(2),dp=2)
        widget_control,u.id.he_bp_x2,set_value=num2str(u.info.position.he.bport(3),dp=2)
        widget_control,u.id.he_bp_y2,set_value=num2str(u.info.position.he.bport(4),dp=2)
        widget_control,u.id.he_bp_z2,set_value=num2str(u.info.position.he.bport(5),dp=2)
        widget_control,u.id.he_bp_x3,set_value=num2str(u.info.position.he.bport(6),dp=2)
        widget_control,u.id.he_bp_y3,set_value=num2str(u.info.position.he.bport(7),dp=2)
        widget_control,u.id.he_bp_z3,set_value=num2str(u.info.position.he.bport(8),dp=2)

        widget_control,u.id.he_ve_x1,set_value=num2str(u.info.position.he.vessel(0),dp=2)
        widget_control,u.id.he_ve_y1,set_value=num2str(u.info.position.he.vessel(1),dp=2)
        widget_control,u.id.he_ve_z1,set_value=num2str(u.info.position.he.vessel(2),dp=2)
        widget_control,u.id.he_ve_x2,set_value=num2str(u.info.position.he.vessel(3),dp=2)
        widget_control,u.id.he_ve_y2,set_value=num2str(u.info.position.he.vessel(4),dp=2)
        widget_control,u.id.he_ve_z2,set_value=num2str(u.info.position.he.vessel(5),dp=2)
        widget_control,u.id.he_ve_x3,set_value=num2str(u.info.position.he.vessel(6),dp=2)
        widget_control,u.id.he_ve_y3,set_value=num2str(u.info.position.he.vessel(7),dp=2)
        widget_control,u.id.he_ve_z3,set_value=num2str(u.info.position.he.vessel(8),dp=2)

;=============

        widget_control,u.id.h_re_x1,set_value=num2str(u.info.position.h.reducer(0),dp=2)
        widget_control,u.id.h_re_y1,set_value=num2str(u.info.position.h.reducer(1),dp=2)
        widget_control,u.id.h_re_z1,set_value=num2str(u.info.position.h.reducer(2),dp=2)
        widget_control,u.id.h_re_x2,set_value=num2str(u.info.position.h.reducer(3),dp=2)
        widget_control,u.id.h_re_y2,set_value=num2str(u.info.position.h.reducer(4),dp=2)
        widget_control,u.id.h_re_z2,set_value=num2str(u.info.position.h.reducer(5),dp=2)
        widget_control,u.id.h_re_x3,set_value=num2str(u.info.position.h.reducer(6),dp=2)
        widget_control,u.id.h_re_y3,set_value=num2str(u.info.position.h.reducer(7),dp=2)
        widget_control,u.id.h_re_z3,set_value=num2str(u.info.position.h.reducer(8),dp=2)

        widget_control,u.id.h_bp_x1,set_value=num2str(u.info.position.h.bport(0),dp=2)
        widget_control,u.id.h_bp_y1,set_value=num2str(u.info.position.h.bport(1),dp=2)
        widget_control,u.id.h_bp_z1,set_value=num2str(u.info.position.h.bport(2),dp=2)
        widget_control,u.id.h_bp_x2,set_value=num2str(u.info.position.h.bport(3),dp=2)
        widget_control,u.id.h_bp_y2,set_value=num2str(u.info.position.h.bport(4),dp=2)
        widget_control,u.id.h_bp_z2,set_value=num2str(u.info.position.h.bport(5),dp=2)
        widget_control,u.id.h_bp_x3,set_value=num2str(u.info.position.h.bport(6),dp=2)
        widget_control,u.id.h_bp_y3,set_value=num2str(u.info.position.h.bport(7),dp=2)
        widget_control,u.id.h_bp_z3,set_value=num2str(u.info.position.h.bport(8),dp=2)

        widget_control,u.id.h_ve_x1,set_value=num2str(u.info.position.h.vessel(0),dp=2)
        widget_control,u.id.h_ve_y1,set_value=num2str(u.info.position.h.vessel(1),dp=2)
        widget_control,u.id.h_ve_z1,set_value=num2str(u.info.position.h.vessel(2),dp=2)
        widget_control,u.id.h_ve_x2,set_value=num2str(u.info.position.h.vessel(3),dp=2)
        widget_control,u.id.h_ve_y2,set_value=num2str(u.info.position.h.vessel(4),dp=2)
        widget_control,u.id.h_ve_z2,set_value=num2str(u.info.position.h.vessel(5),dp=2)
        widget_control,u.id.h_ve_x3,set_value=num2str(u.info.position.h.vessel(6),dp=2)
        widget_control,u.id.h_ve_y3,set_value=num2str(u.info.position.h.vessel(7),dp=2)
        widget_control,u.id.h_ve_z3,set_value=num2str(u.info.position.h.vessel(8),dp=2)

        widget_control,u.id.h_cr_x1,set_value=num2str(u.info.position.h.he_crystal(0),dp=2)
        widget_control,u.id.h_cr_y1,set_value=num2str(u.info.position.h.he_crystal(1),dp=2)
        widget_control,u.id.h_cr_z1,set_value=num2str(u.info.position.h.he_crystal(2),dp=2)
        widget_control,u.id.h_cr_x2,set_value=num2str(u.info.position.h.he_crystal(3),dp=2)
        widget_control,u.id.h_cr_y2,set_value=num2str(u.info.position.h.he_crystal(4),dp=2)
        widget_control,u.id.h_cr_z2,set_value=num2str(u.info.position.h.he_crystal(5),dp=2)
        widget_control,u.id.h_cr_x3,set_value=num2str(u.info.position.h.he_crystal(6),dp=2)
        widget_control,u.id.h_cr_y3,set_value=num2str(u.info.position.h.he_crystal(7),dp=2)
        widget_control,u.id.h_cr_z3,set_value=num2str(u.info.position.h.he_crystal(8),dp=2)              
        
        widget_control,base,/realize

        get_image,u

        plot_raytracing,u
        plot_crystal_he,u

        widget_control,u.id.atnu1,set_value=num2str(*u.atnu1)
        widget_control,u.id.atnu2,set_value=num2str(*u.atnu2)
        widget_control,u.id.atnu3,set_value=num2str(*u.atnu3)
        widget_control,u.id.atnu4,set_value=num2str(*u.atnu4)

        tmp=*u.image
        widget_control,u.id.counts,set_value=num2str(tmp(xi,zeta))

	xmanager,'w_hirexsr_raytracing',base

END
