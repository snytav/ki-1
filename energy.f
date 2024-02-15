!17.08.06 - введено u0 стр.36 - вычитается при вычислении энергии электронов фона
!Программа вычисления энергии, мат. ожидания, дисперсии
!Вычисление энергии
      subroutine energy
      implicit real*8(a-h,o-z)
      include 'part.pf'
      
      parameter(fm=1000) 
      real*8 ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp),f(30),
     *hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp)
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      real*8 vec(1000),hh,elex,eley,elez,pi,epl,vecb(1000)
      real*8 efx,efy,efz,ebmx,ebmy,ebmz,cix,ciy,ciz
      real*8 imp_b,ibx,iby,ibz,imp_i,iix,iiy,iiz,imp_f,ifx,ify,ifz,eex,eey,eez
      common/a/tx,nt,ml,sst,ns,nt1,nt2		
      common/d/h1,h2,h3
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
      common/f/tau
      common/g/hx,hy,hz    
      common/h/ex,ey,ez
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/iff/xf,yf,zf,uf,vf,wf,amf
      common/c1/alpha,t0,tD
      common/et/ei0,ef0,eh0,ee0,ew0
      common/p/p_tau
      common/cc/ncc
       
      if(enout.ne.1) return
!----------------------------------------

      imp_b=0.d0
      imp_i=0.d0    
      imp_f=0.d0    
      ibx=0.d0
      iby=0.d0
      ibz=0.d0
      iix=0.d0
      iiy=0.d0
      iiz=0.d0
      ifx=0.d0
      ify=0.d0
      ifz=0.d0
     
         
      do j=1,jmi
      ibx=ibx+ub(j)
      iby=iby+vb(j)
      ibz=ibz+wb(j)
      imp_b=imp_b + sqrt(ub(j)**2+vb(j)**2+wb(j)**2)
      iix=iix+ui(j)
      iiy=iiy+vi(j)
      iiz=iiz+wi(j)
      imp_i=imp_i + sqrt(ui(j)**2+vi(j)**2+wi(j)**2)
      ifx=ifx+uf(2*j)+uf(2*j-1)
      ify=ify+vf(2*j)+vf(2*j-1)
      ifz=ifz+wf(2*j)+wf(2*j-1)
    
      imp_f=imp_f + sqrt(uf(2*j)**2+vf(2*j)**2+wf(2*j)**2)+
     + sqrt(uf(2*j-1)**2+vf(2*j-1)**2+wf(2*j-1)**2)
      enddo

      imp_b=imp_b/jmb
      imp_i=imp_i/jmi
      imp_f=imp_f/jmf
!Energy f
      ef=0.d0
      efg=0d0
      
      efx=0d0
      efy=0d0
      efz=0d0
c       ONLY PLASMA ELECTRONS      
      do j=1,jmf
c       ef=ef+(uf(j)**2+vf(j)**2+wf(j)**2)
         efg=efg+dsqrt(1.d0+beta0**2*
     *   (uf(j)**2+vf(j)**2+wf(j)**2))-1.d0
        
         efx=efx+dsqrt(1.d0+beta0**2*(uf(j)**2))-1.d0
         efy=efy+dsqrt(1.d0+beta0**2*(vf(j)**2))-1.d0
         efz=efz+dsqrt(1.d0+beta0**2*(wf(j)**2))-1.d0
      
 490     format('ef ',i10,6e15.5)
      enddo
      ef=-ef*amf/((imp-2)*(lmp-2)*(lmp-2))
      efg=-efg*amf/((imp-2)*(lmp-2)*(lmp-2)) 
      
      efx=-efx*amf/((imp-2)*(lmp-2)*(lmp-2))
      efy=-efy*amf/((imp-2)*(lmp-2)*(lmp-2))
      efz=-efz*amf/((imp-2)*(lmp-2)*(lmp-2))

      ew = 0d0
      ewg=0d0
c     BEAM ELECTRONS      
      
      ebmx=0d0
      ebmy=0d0
      ebmz=0d0
      
      do j=1,jmb
         ew =ew+(ub(j)**2+vb(j)**2+wb(j)**2)
         ewg=ewg+dsqrt(1.d0+beta0**2*
     *   (ub(j)**2+vb(j)**2+wb(j)**2))-1.d0

         ebmx=ebmx+dsqrt(1.d0+beta0**2*(ub(j)**2))-1.d0
         ebmy=ebmy+dsqrt(1.d0+beta0**2*(vb(j)**2))-1.d0
         ebmz=ebmz+dsqrt(1.d0+beta0**2*(wb(j)**2))-1.d0

      enddo
      ew= -ew*amb/((imp-2)*(lmp-2)*(lmp-2))
      ewg=-ewg*amb/((imp-2)*(lmp-2)*(lmp-2))
      ebmx=-ebmx*amb/((imp-2)*(lmp-2)*(lmp-2))
      ebmy=-ebmy*amb/((imp-2)*(lmp-2)*(lmp-2))
      ebmz=-ebmz*amb/((imp-2)*(lmp-2)*(lmp-2))
c     ############# Energy of ions #################
      eion=0d0
      do j=1,jmi
         eion=eion+dsqrt(1.d0+beta0**2*
     *   (ui(j)**2+vi(j)**2+wi(j)**2))-1.d0
      enddo
      eion= eion*ami/((imp-2)*(lmp-2)*(lmp-2))

c     ##############################################
c      epl=-epl*amf*h1*h2*h3
	ns1=idint(sst)
!----------------------------------------
      ee=0.d0
      eex=0.d0
      eey=0.d0
      eez=0.d0
      
      hh=0d0
      do k=2,kmp-1
         do l=2,lmp-1
            do i=2,imp-1
               ee=ee+ex(i,l,k)**2+ey(i,l,k)**2+ez(i,l,k)**2
               hh=hh+hx(i,l,k)**2+hy(i,l,k)**2+hz(i,l,k)**2
               eex=eex+ex(i,l,k)
               eey=eey+ey(i,l,k)
               eez=eez+ez(i,l,k)
            end do
         end do
      end do
      pi = 2d0*dasin(1d0)
      ee =ee /(4d0*(imp-2)*(lmp-2)*(lmp-2)*beta0**2)
      eex=eex/(4d0*(imp-2)*(lmp-2)*(lmp-2)*beta0**2)
      eey=eey/(4d0*(imp-2)*(lmp-2)*(lmp-2)*beta0**2)
      eez=eez/(4d0*(imp-2)*(lmp-2)*(lmp-2)*beta0**2)
      hh=hh/(4d0*(imp-2)*(lmp-2)*(lmp-2)) !*h1*h2*h3*beta0**2/2.d0
!----------------------------------------
            if(deb.ge.4) then
               write(37,*) 'ewfull ',ew 
            endif
        
       	if(mod((nt-1),ns1).eq.0) then
           vec(1) = ef
      	   vec(2) = eion
      	   vec(3) = imp_b
      	   vec(4) = imp_i
      	   vec(5) = imp_f
      	   vec(6) = ee
      	   vec(7) = hh
      	   vec(8) = elez
      	   vec(11) = efg
      	   vec(12) = ewg
      	   
      	   vec(21) = efx
      	   vec(22) = efy
      	   vec(23) = efz
      	   vec(24) = ebmx
      	   vec(25) = ebmy
      	   vec(26) = ebmz
        
           do i = 1,1000
              vecb(i) = vec(i)
           enddo
      call PARAreduceLMP2(vec)
           
           imp_b=vec(3)/jmb
           imp_i=vec(4)/jmi
           imp_f=vec(5)/jmf
           
       
       	   ef = vec(1)
      	   eion = vec(2)
      
      	   efx = vec(21) 
      	   efy = vec(22) 
      	   efz = vec(23) 
      
      	   ebmx = vec(24) 
      	   ebmy = vec(25) 
      	   ebmz = vec(26) 
      
      	   efg=vec(11)
      	   ewg=vec(12) 
           ew=ei+ef+ee
      	   if(nt.eq.1) then
      	      ef0 = ef
      	      ei0 = ei
      	      ee0 = ee
      	      ew0 = ew
      	   endif
            epl = ee+hh+efg+ewg+eion

      if(me.eq.0) then
      write(22,302) nt*tau,ibx,iby,ibz,iix,iiy,iiz,ifx,ify,ifz,
     + sqrt((ibx+iix+ifx)**2+(iby+iiy+ify)**2+(ibz+iiz+ifz)**2)
      write(21,300) nt*tau,ee,hh,efg,ewg,eion,epl,
     +imp_f,imp_b,imp_i,imp_f+imp_b+imp_i
	   endif      
        end if
  300 format(11e20.10)           
  302 format(11e20.10)           
        call curbalance
      
      end subroutine energy

!*******************************************************************
!Мат. ожидание и дисперсия на каждом шаге по времени
      subroutine M_D(uf,vf,wf,iflag,jm)
      implicit real*8(a-h,o-z)
      include 'part.pf'
      parameter(nx=1000,ny=20,nz=20)
      real*8 uxmin,uxmax,upmax,hxg,hyg,hzg,si1,si2,si3
      real*8 uf(jmp),vf(jmp),wf(jmp)	
      real*8 Npx(nx+1),Npy(ny+1),Npz(nz+1)
      real*8 Mx(nx+1),My(ny+1),Mz(nz+1),Mx0(nx+1),My0(ny+1),Mz0(nz+1),
     *	Mpx,Mpy,Mpz,Mpx0(2),Mpy0(2),Mpz0(2),Mp,Mpx1,Mpy1,Mpz1	! Mpy0(iflag)
      real*8 Dx(nx+1),Dy(ny+1),Dz(nz+1),Dx0(nx+1),Dy0(ny+1),Dz0(nz+1),
     *	Dpx,Dpy,Dpz,Dpx0,Dpy0,Dpz0,Dp
      real*8 ksi,vec(1000)
      integer jmall,nt1,nt2
      common/Mpp/Mpx0,Mpy0,Mpz0	
      common/a/tx,nt,ml,sst,ns,nt1,nt2
      common/p/p_tau
      common/f/tau
      common/c1/alpha,t0,tD
      common/start_p/Te0,ksi,ak0
!------------------------------------------------------------------	
      	Mpx=0.d0
      	Mpy=0.d0
      	Mpz=0.d0
      do 106 j=1,jm
	   Mpx=Mpx+uf(j)	
	   Mpy=Mpy+vf(j)		
	   Mpz=Mpz+wf(j)	
  106 continue
	   vec(1) = Mpx
	   vec(2) = Mpy
	   vec(3) = Mpz
	   vec(4) = jm
	   call PARAreduceLMP2(vec)
	   Mpx   = vec(1)
	   Mpy   = vec(2)
	   Mpz   = vec(3)
	   jmall = vec(4)
	Mpx=Mpx/jmall
	Mpy=Mpy/jmall
	Mpz=Mpz/jmall
c******************************	
	if(nt.eq.1)then
   	   Mpx0(iflag) = Mpx
	   Mpy0(iflag) = Mpy
	   Mpz0(iflag) = Mpz
	end if
	Mpx1=Mpx-Mpx0(iflag) 
	Mpy1=Mpy-Mpy0(iflag) 
	Mpz1=Mpz-Mpz0(iflag) 
c******************************	
      	Mp=(Mpx+Mpy+Mpz)/3.d0
!------------------------------------------------------------------	
       	Dpx=0.d0
      	Dpy=0.d0
      	Dpz=0.d0
      do 107 j=1,jm
      	   Dpx=Dpx+(uf(j)-Mpx)**2	
       	   Dpy=Dpy+(vf(j)-Mpy)**2	
      	   Dpz=Dpz+(wf(j)-Mpz)**2	
  107 continue
	   vec(1) = Dpx
	   vec(2) = Dpy
	   vec(3) = Dpz
	   call PARAreduceLMP2(vec)
	   Dpx   = vec(1)
	   Dpy   = vec(2)
	   Dpz   = vec(3)

	Dpx=Dpx/jmall
	Dpy=Dpy/jmall
	Dpz=Dpz/jmall
	Dp=(Dpx+Dpy+Dpz)/3.d0

	ss1=(nt-1.d0)/sst
	ns1=ss1
	ss2=ss1-ns1
	
	if((ss2.le.1e-8).and.(me.eq.0).and.(outterm.eq.1)) then
   	   if(iflag.EQ.1)then  
	      write(49,308) nt*tau,Mp,Dp	!/ak0      
	   endif
	   if(iflag.EQ.2)then 
	      write(50,308) nt*tau,Mp,Dp/alpha	!/ak0/alpha       
	   endif
	endif
       
  308 format(3e25.10)           
      
      end subroutine M_D

