      subroutine four3D(vl,vh,imx,imy,imz,m7,dnm,nm,dflg,dfhf,dfst)
      implicit none
      include 'part.pf'
      integer imx,imy,imz,i,j,k,m7,l
      real*8 vl(imp,lmp,kmp),dflg,dfhf,dfst,df
      real*8 v(imp,lmpf,kmp),vh(imp,lmpf,kmp),kr,kr2,dsqrt,Km
      real*8 hx(imp),hy(lmpf),hz(kmp),vx(imp),vy(lmpf),vz(kmp) 
      character*20 tr,tr1
      character*4  nm,dnm

      do l = 1,lmp
         !write(37,*) l,vl(2,l,2)
      enddo 

      call PARAgather(vl,v)
      if(deb.ge.2) then
         write(37,*) 'X'    
      endif


         if(me.eq.0) then       
            tr = dnm(1:4)    
            write(tr1,'(i6.6)') m7 
            tr=tr(1:4)//tr1 
            tr=tr(1:10)//'.dat' 
      
            open(91,file=tr,form='formatted')

            do i = 2,imx-1
               do j = 2,lmpf-1
                  do k = 2,imz-1
  161                format(3i10,e30.20)	    
 	             write(91,161) i,j,k,v(i,j,k)
       	          enddo	    
               enddo
            enddo 
            close(91)
         endif	    
         return

      	 
      do i = 1,imx
         do j = 1,lmpf
            do k = 1,imz
               vz(k) = v(i,j,k)	    
       	    enddo
c            call fourier(vz,hz,kmp) 
            do k = 1,imz
               vh(i,j,k) = hz(k) 	    
       	    enddo	    
         enddo
      enddo 
      if(deb.ge.2) then
         write(37,*) 'Y'   
      endif	 
      
      do i = 1,imx
         do k = 1,imz
            do j = 1,lmpf
               vy(j) = vh(i,j,k)	    
       	    enddo
c            call fourier(vy,hy,lmpf) 
            do j = 1,lmpf
               vh(i,j,k) = hy(j) 	    
       	    enddo	    
         enddo
      enddo 

      if(deb.ge.2) then
         write(37,*) 'Z'   
      endif	 

      do k = 1,imz
         do j = 1,lmpf
            do i = 1,imx
               vx(i) = vh(i,j,k)	    
       	    enddo
c            call fourier(vx,hx,imp) 
            do i = 1,imx
               vh(i,j,k) = hx(i) 	    
       	    enddo	    
         enddo
      enddo 

      if(me.eq.0) then
      
         df   = 0d0
         dflg = 0d0
      	 dfhf = 0d0
      	 dfst = 0d0 
         
	 kr = (imx-2)**2+(lmpf-2)**2+(imz-2)**2  	 
	 Km = dsqrt(kr)
         do i = 2,imx-1
            do j = 2,lmpf-1
               do k = 2,imz-1
                  kr2 = (i-1)**2+(j-1)**2+(k-1)**2
                  kr  = dsqrt(kr2)
        
       	          df = df + vh(i,j,k)**2
        
                  if(kr.lt.(Km/4d0)) then
      		     dflg = dflg + vh(i,j,k)**2
      		  else
      		     if(kr.gt.(Km/2d0)) then
      		        dfhf = dfhf + vh(i,j,k)**2
      		     endif
      		     if(kr.gt.(Km*0.9d0/dsqrt(3d0))) then
      		           dfst = dfst + vh(i,j,k)**2
      		     endif
      		  endif         		  
       	       enddo	    
            enddo
         enddo 

         if((outf.eq.1).and.(me.eq.0)) then       
            tr = nm(1:4)    
            write(tr1,'(i6.6)') m7 
            tr=tr(1:4)//tr1 
            tr=tr(1:10)//'.dat' 
      
            open(91,file=tr,form='formatted')

            do i = 2,imx-1
               do j = 2,lmpf-1
                  do k = 2,imz-1
  160                format(3i10,2e30.20)	    
                     kr2 = (i-1)**2+(j-1)**2+(k-1)**2
                     kr  = dsqrt(kr2)
c      		     if(dabs(vh(i,j,k)).gt.1.0d-8) then
       	                write(91,160) i,j,k,kr,dabs(vh(i,j,k))
c       		     endif	
       	          enddo	    
               enddo
            enddo 
            close(91)
         endif	    
           
      endif
      
      
      if(df.gt.1.0d-8) then
         dflg = dflg/df
         dfhf = dfhf/df
         dfst = dfst/df
      endif	 
      
      end

      subroutine fourier(en,enh,num1)
      implicit none
      include 'part.pf'
      integer im
      parameter(im=10000) 
      real*8 en(im),pi,tau,h(im+1),dlog,enh(im),dn
      complex*16 harm(im)
      integer num,i,k1,k,nlog,num1
      
      num = num1-2
      do i = 2,num+1
         harm(i-1) = en(i)
c         print*,i,harm(i)
      enddo
      dn = num
      nlog = idint(dlog(dn)/dlog(2d0))
      call fftc(harm,nlog,-1,num)
c      print*,'++++++++++++++++++'


      do i = num,1,-1
c         print*,i,nlog,en(i),harm(i)
      enddo 
c      stop
      k1 = num/2
      
      do k=1,k1+1 
         enh(k+1)=dreal(harm(k))/dn
      enddo 
c      stop
      do k=2,k1 
         enh(k1+k+1)=dimag(harm(k))/dn
      enddo 
c      stop
      end
         
      subroutine fftc(a,n,isi,np) 
      implicit real*8(a-h,o-z)
      parameter(im=10000) 
      integer np,isi,i,nn,j,m,mm,ii,n,isign
      real*8  pi2,th
      complex*16 a(im),t,w,w1 
      
      do 7 i=1,np 
      if(isi.gt.0) a(i)=a(i)/np 
    7 continue 
      pi2=8.d0*datan(1.d0) 
      nn=np 
      j=1 
      do 3 i=1,nn 
      if(i.ge.j) go to 1 
      t=a(j) 
      a(j)=a(i) 
      a(i)=t 
    1 m=nn/2 
    2 if(j.le.m) go to 3 
      j=j-m 
      m=m/2 
      if(m.ge.1) go to 2 
    3 j=j+m 
      mm=1 
    4 if(mm.ge.nn) return 
      ii=2*mm 
c      th= pi2/dsign(ii,n*isi) 
      th= pi2/isign(ii,n*isi) 
      w1=dcmplx(-2.0d0*dsin(th/2)**2,dsin(th)) 
      w=1 
      do 6 m=1,mm 
      do 5 i=m,nn,ii 
      t=w*a(i+mm) 
      a(i+mm)=a(i)-t 
      a(i)=a(i)+t 
    5 continue 
      w=w1*w+w 
    6 continue 
      mm=ii 
      go to 4 
      end 

c     xf - wf - coordinates and impulses
c     y0 - lower subdomain boundary
c     jm - number of particles (local)
c     Dx  - temperature distribution (out)
      subroutine tempSort(xf,yf,zf,uf,vf,wf,y0,jm,Dx,Dy,Dz,Mx,My,Mz)
      implicit none
      include 'part.pf'
      real*8 uf(jmp),vf(jmp),wf(jmp),h1,h2,h3,y0,minT,maxT	
      real*8 Mx(imp,lmp,kmp),Dx(imp,lmp,kmp),s2,s4,s6
      real*8 My(imp,lmp,kmp),Mz(imp,lmp,kmp),
     +       Dy(imp,lmp,kmp),Dz(imp,lmp,kmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),pi,tt,tfl2,tfl3,tfl0,tfl
      real*8 vec(1000),u,v,w,ps
      integer i,l,k,j,Num(imp,lmp,kmp),jm
      integer lp,jmf,jmi,n1,n2,m5,m7
      common/d/h1,h2,h3
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7
      
      if(deb.ge.4) write(37,*) 'in tempSort'
      do i = 1,imp
         do l = 1,lmp
            do k = 1,kmp
               Num(i,l,k) = 0	    
      	       Mx(i,l,k) = 0d0
      	       My(i,l,k) = 0d0
      	       Mz(i,l,k) = 0d0
      	       Dx(i,l,k) = 0d0
      	       Dy(i,l,k) = 0d0
      	       Dz(i,l,k) = 0d0
      	    enddo
      	 enddo
      enddo     
      if(deb.ge.4) write(37,*) 'in tempSort2'
      
      do j = 1,jm
         s2 = xf(j)/h1
         i  = idint(s2+1.5d0)
         s4 = (yf(j)-y0)/h2
         l  = idint(s4+1.5d0)
         s6 = zf(j)/h3
         k  = idint(s6+1.5d0)
	 if(deb.ge.7) write(37,*) i,l,k,j,xf(j),yf(j),y0
   
         ps=1.d0/dsqrt(1.d0+beta0**2*(uf(j)**2+vf(j)**2+wf(j)**2))
	 u = ps*uf(j)
	 v = ps*vf(j)
	 w = ps*wf(j)
	 
         Num(i,l,k) = Num(i,l,k) + 1
         Mx(i,l,k)  = Mx(i,l,k)  + u  	  
         My(i,l,k)  = My(i,l,k)  + v  	  
         Mz(i,l,k)  = Mz(i,l,k)  + w  	  
      enddo       
      if(deb.ge.4) write(37,*) 'after Mxyz'
       
      do i = 2,imp-1
         do l = 2,lmp-1
            do k = 2,kmp-1
               if(Num(i,l,k).gt.0) then
                  Mx(i,l,k) = Mx(i,l,k)/Num(i,l,k)
                  My(i,l,k) = My(i,l,k)/Num(i,l,k)
                  Mz(i,l,k) = Mz(i,l,k)/Num(i,l,k)
      	       endif	       
      	    enddo
      	 enddo
      enddo     
        
      do j = 1,jm
         s2 = xf(j)/h1
         i  = idint(s2+1.5d0)
         s4 = (yf(j)-y0)/h2
         l  = idint(s4+1.5d0)
         s6 = zf(j)/h3
         k  = idint(s6+1.5d0)

         ps=1.d0/dsqrt(1.d0+beta0**2*(uf(j)**2+vf(j)**2+wf(j)**2))
	 u = ps*uf(j)
	 v = ps*vf(j)
	 w = ps*wf(j)
         
         Dx(i,l,k) = Dx(i,l,k) + (u - Mx(i,l,k))**2	 	 
         Dy(i,l,k) = Dy(i,l,k) + (v - My(i,l,k))**2	 	 
         Dz(i,l,k) = Dz(i,l,k) + (w - Mz(i,l,k))**2	 	 
      enddo       
      if(deb.ge.4) write(37,*) 'after Dxyz'
c      call PARAFinal
c      stop    
      do i = 2,imp-1
         do l = 2,lmp-1
            do k = 2,kmp-1
               if(Num(i,l,k).gt.0) then
                  Dx(i,l,k) = Dx(i,l,k)/Num(i,l,k)
                  Dy(i,l,k) = Dy(i,l,k)/Num(i,l,k)
                  Dz(i,l,k) = Dz(i,l,k)/Num(i,l,k)
      	       endif	       
      	    enddo
      	 enddo
      enddo     

      end

      subroutine densSort(xf,yf,zf,y0,jm,a,md,px,label)
      implicit none
      include 'part.pf'
      real*8 uf(jmp),vf(jmp),wf(jmp),h1,h2,h3,y0	
      real*8 px(imp,lmp,kmp),s2,s4,s6,a,s21,s41,s61,s,a1
      real*8 xf(jmp),yf(jmp),zf(jmp),sst,tx
      integer i,l,k,j,Num(imp,lmp,kmp),jm,md,ml,nt,ns,nt1,nt2,label
      common/d/h1,h2,h3
      common/a/tx,nt,ml,sst,ns,nt1,nt2
      

      if(deb.ge.4) write(37,*) 'jm_dens ',jm,a,md      

      if(md.eq.0) then
         do i = 1,imp
            do l = 1,lmp
               do k = 1,kmp
      	          px(i,l,k) = 0d0
      	       enddo	  
      	    enddo
      	 enddo
      endif     
      if(deb.ge.4) write(37,*) 'jm_dens ',jm,a,md      
      do j = 1,jm
	 if(deb.ge.2) then
c	    write(37,367) i,l,k,j,xf(j),yf(j),zf(j),s2,s4,s6,h1,h2
	 endif
         s2 = xf(j)/h1
         i=idint(s2+1.5d0)
         s2=s2+1.5d0-i
         s4=(yf(j)-y0)/h2
         l=idint(s4+1.5d0)
         s4=s4+1.5d0-l
         s6=zf(j)/h3
         k=idint(s6+1.5d0)
         s6=s6+1.5d0-k
         
      	 s21=1.d0-s2
         s41=1.d0-s4
         s61=1.d0-s6

 367     format('dens367 ',i5,'lab ',i3,e15.5,3i5,i10,8e30.20)         
	 if(deb.ge.2) then
	    write(37,367) nt,label,a,i,l,k,j,xf(j),
     +                    yf(j),zf(j),s2,s4,s6,h1,h2
	 endif
	 a1 = a    
	 
      	 s=s21*s41*a1
         px(i,l,k)=px(i,l,k)+s*s61
         px(i,l,k+1)=px(i,l,k+1)+s*s6
         s=s21*s4*a1
         px(i,l+1,k)=px(i,l+1,k)+s*s61
         px(i,l+1,k+1)=px(i,l+1,k+1)+s*s6
         s=s2*s41*a1
         px(i+1,l,k)=px(i+1,l,k)+s*s61
         px(i+1,l,k+1)=px(i+1,l,k+1)+s*s6
         s=s2*s4*a1
         px(i+1,l+1,k)=px(i+1,l+1,k)+s*s61
         px(i+1,l+1,k+1)=px(i+1,l+1,k+1)+s*s6
c	 if(deb.ge.2) write(37,367) i,l,k,j,xf(j),yf(j),zf(j),y0,
c     +	              px(i,l,k),px(i,l+1,k),px(i+1,l,k),px(i+1,l+1,k)

      enddo       
      
      do i=1,imp-1
         do l=1,lmp-1
c            px(i,l,2)     = px(i,l,2)+px(i,l,kmp-1)
c            px(i,l,kmp-1) = px(i,l,2)
         enddo
      enddo

      do i=1,imp-1
         do l=1,kmp-1
c            px(i,2,l)     = px(i,2,l)+px(i,lmp-1,l)
c            px(i,lmp-1,2) = px(i,2,l)
         enddo
      enddo
      do i=1,kmp-1
         do l=1,lmp-1
c            px(2,l,i)     = px(2,l,i)+px(imp-1,l,i)
c            px(imp-1,l,i) = px(2,l,i)
         enddo
      enddo
      
       
      do i = 1,imp
         do l = 1,lmp
         
  223       format(2f10.3,e15.5)
        	    
            if(deb.ge.7) then
c	       write(37,223) i*h1,l*h2, px(i,l,5)
	    endif   
      	 enddo
      enddo     
           
      end

      subroutine getMaxHarm(ex,ki,kl,kk,kr,exm)
      implicit none
      include 'part.pf'
      real*8 ex(imp,lmpf,kmp),kr,exm,kr2
      integer ki,kl,kk,i,j,k
      
      exm = 0d0
      ki = 0
      kl = 0
      kk = 0
      kr = 0
      
         do i = 2,imp-1
            do j = 2,lmpf-1
               do k = 2,kmp-1

                  if(((i.gt.2).or.(k.gt.2).or.(j.gt.2))
     +		     .and.(exm.lt.ex(i,j,k))) then
                     exm = ex(i,j,k)
		     ki = i
		     kl = j
		     kk = k 
                  endif
       	       enddo	    
            enddo
         enddo 
         kr2 = (ki-1)**2+(kl-1)**2+(kk-1)**2
         kr  = dsqrt(kr2)
      end

      subroutine four3Dgrow(vl,vh,vh1,imx,imy,imz,
     +                      m7,dnm,nm,dflg,dfhf,dfst)
      implicit none
      include 'part.pf'
      integer imx,imy,imz,i,j,k,m7,l
      real*8 vl(imp,lmp,kmp),vh1(imp,lmp,kmp),dflg,dfhf,dfst,df
      real*8 v(imp,lmpf,kmp),vh(imp,lmpf,kmp),kr,kr2,dsqrt,Km
      real*8 hx(imp),hy(lmpf),hz(kmp),vx(imp),vy(lmpf),vz(kmp) 
      character*20 tr,tr1,nm,dnm

      call PARAgatherGrow(vl,v)
      
      if(deb.ge.2) then
         write(37,*) 'X ',v(3,3,3)    
      endif
      	 
      do i = 1,imx
         do j = 1,lmpf
            do k = 1,imz
               vz(k) = v(i,j,k)	    
       	    enddo
c            call fourier(vz,hz,kmp) 
            do k = 1,imz
               vh(i,j,k) = hz(k) 	    
       	    enddo	    
         enddo
      enddo 
      if(deb.ge.2) then
         write(37,*) 'Y ',vh(3,3,3)   
      endif	 
      
      do i = 1,imx
         do k = 1,imz
            do j = 1,lmpf
               vy(j) = vh(i,j,k)	    
       	    enddo
c            call fourier(vy,hy,lmpf) 
            do j = 1,lmpf
               vh(i,j,k) = hy(j) 	    
       	    enddo	    
         enddo
      enddo 

      if(deb.ge.2) then
         write(37,*) 'Z ',vh(3,3,3)   
      endif	 

      do k = 1,imz
         do j = 1,lmpf
            do i = 1,imx
               vx(i) = vh(i,j,k)	    
       	    enddo
c            call fourier(vx,hx,imp) 
            do i = 1,imx
               vh(i,j,k) = hx(i) 	    
       	    enddo	    
         enddo
      enddo 
c      write(37,*) vh1(3,3,3),vh(3,3,3)

      do k = 1,imz
         do j = 2,lmp-1
            do i = 1,imx
               vh1(i,j,k) = vh(i,meh*(lmp-2)+j,k)
c	       write(37,*) i,j,k,vh1(i,j,k),vh(i,me*(lmp-2)+j,k)
       	    enddo
         enddo
      enddo 
c      write(37,*) vh1(3,3,3),vh(3,3,3)

      if(me.eq.0) then
      
         df   = 0d0
         dflg = 0d0
      	 dfhf = 0d0
      	 dfst = 0d0 
         
	 kr = (imx-2)**2+(lmpf-2)**2+(imz-2)**2  	 
	 Km = dsqrt(kr)
         do i = 2,imx-1
            do j = 2,lmpf-1
               do k = 2,imz-1
                  kr2 = (i-1)**2+(j-1)**2+(k-1)**2
                  kr  = dsqrt(kr2)
        
       	          df = df + vh(i,j,k)**2
        
                  if(kr.lt.(Km/4d0)) then
      		     dflg = dflg + vh(i,j,k)**2
      		  else
      		     if(kr.gt.(Km/2d0)) then
      		        dfhf = dfhf + vh(i,j,k)**2
      		     endif
      		     if(kr.gt.(Km*0.9d0/dsqrt(3d0))) then
      		           dfst = dfst + vh(i,j,k)**2
      		     endif
      		  endif         		  
       	       enddo	    
            enddo
         enddo 

         if((outf.eq.1).and.(me.eq.0)) then       
            tr = nm(1:4)    
            write(tr1,'(i3.3)') m7 
            tr=tr(1:4)//tr1 
            tr=tr(1:7)//'.dat' 
      
            open(91,file=tr,form='formatted')

            do i = 2,imx-1
               do j = 2,lmpf-1
                  do k = 2,imz-1
  160                format(3i10,2e30.20)	    
                     kr2 = (i-1)**2+(j-1)**2+(k-1)**2
                     kr  = dsqrt(kr2)
c      		     if(dabs(vh(i,j,k)).gt.1.0d-8) then
       	                write(91,160) i,j,k,kr,dabs(vh(i,j,k))
c       		     endif	
       	          enddo	    
               enddo
            enddo 
            close(91)
         endif	    
           
      endif
      
      if(deb.ge.4) write(37,*) 'end four3Dgrow ',vh1(3,3,3)
      
      if(df.gt.1.0d-8) then
         dflg = dflg/df
         dfhf = dfhf/df
         dfst = dfst/df
      endif	 
      
      end

      subroutine selectGrow(vh1,imx,imy,imz,m7)
      implicit none
      include 'part.pf'
      integer imx,imy,imz,i,j,k,m7,l,knum,topn(grtot),nmax,grnum
      real*8 vh2(imp,lmp,kmp),vh1(imp,lmp,kmp),t,hmax,topres(2*grtot)
      real*8 kr,kr2,dsqrt,Km,tmax,toph(grtot),topg(grtot)
      common/vvhh/vh2,grnum
      common/top/toph,topg,topn
      
      if(deb.ge.4) write(37,*) 'in selectGrow ',vh1(3,3,3)
      if(m7.eq.1) then
         grnum = 0
         do i = 1,grtot
	    topn(i) = 0 
	    topg(i) = 0d0
	    toph(i) = 0d0
      	 enddo  	 
      	 if(me.eq.0) then
            open(87,file='grow.dat',form='formatted')	 
      	 endif    
      endif
      
      if(m7.gt.1) then
         tmax = 0d0
      if(deb.ge.4) write(37,*) 'i,topn ',i,topn(1)
       
      if(deb.ge.4) write(37,*) 'in if m7'
         do k = 2,imz-1
            do j = 2,lmp-1
               do i = 2,imx-1
                  knum = 1000000*i + 1000*(meh*(lmp-2) + j)+k	       
         	  t  = 0d0
                  if(dabs(vh2(i,j,k)).gt.1d-12) then
      		    t = dabs(vh1(i,j,k)/vh2(i,j,k))	    
      		  endif   
 620              format(3i5,3e15.5)         
                  if(deb.ge.4) then
        	     write(37,620) i,j,k,t,vh1(i,j,k),vh2(i,j,k)	    
		  endif     
       		  
      		  if(tmax.lt.t) then
      		     tmax = t
      		     nmax = knum
		     if(deb.ge.4) write(37,*) 'max ',nmax,tmax
      		  endif 
       	       enddo
            enddo
         enddo 
      if(deb.ge.4) write(37,*) 'i,topn ',i,topn(1)
      if(deb.ge.4) write(37,*) 'after if m7'
      	 
         do i = 1,grtot
       	    if(deb.ge.4) write(37,*) 'i,topn ',i,topn(i)
      
            if((topn(i).eq.nmax).or.(topn(i).eq.0)) then
      	       if(deb.ge.4) write(37,*) 'topn,nmax',topn(i),nmax
      	       if(topn(i).eq.0) then   
	          if((nmax.gt.2002002).and.
     +   	     (nmax.lt.(1000000*(imp-1) + 
     +                         1000*(lmpf-1)+kmp-1))) then
       	             topn(i) = nmax
		  else
		     goto 600 
		  endif     
       	       endif	  
      	       grnum = grnum + 1
       	       goto 600
      	    endif 	    
      	 enddo  	 
  600    continue
         if(deb.ge.4) write(37,*) 'nmax,grnum ',nmax,grnum
      
         do j = 1,grnum
            i = topn(j)/1000000
            l = mod(topn(j),1000000)/1000 - meh*(lmp-2)
            k = mod(topn(j),1000)
            if(deb.ge.4) write(37,*) 'j,i,l,k ',j,i,l,k,topn(j)

            toph(j) = vh1(i,l,k)        
      	    if(dabs(vh2(i,l,k)).gt.1d-12) then
               topg(j) = dabs(vh1(i,l,k)/vh2(i,l,k))
      	    endif   
	     if(deb.ge.4) write(37,*) 'topg ',topn(j),topg(j)

      	    topn(j) = 1000000*i + 1000*(l+ meh*(lmp-2))+k	       
	     if(deb.ge.4) write(37,*) 'topn(j) ',topn(j)

      	 enddo  	 
         if(deb.ge.4) write(37,*) 'b gath3'
         
         call PARAgather3(topn,topg,toph,m7,topres)
         if(deb.ge.4) write(37,*) 'a gath3'
          
         if(me.eq.0) then
 659        format(60e20.8)
            write(87,659) (topres(i),i=1,2*grtot) 	    
         endif    
        
      endif
          
      do k = 1,imz
         do j = 2,lmp-1
            do i = 1,imx
               vh2(i,j,k) = vh1(i,j,k)	    
       	    enddo
         enddo
      enddo 
           
      end

      subroutine selectHarmList(vh1,filename,fid)
      implicit none
      include 'part.pf'
      integer imx,imy,imz,i,j,k,m7,l,n,knum,delta,fid
      real*8 vh1(imp,lmpf,kmp)
      real*8 lst(90),pi,s27,mxh
      character*20 filename
      integer lp,jmf,jmi,n1,n2,m5,m27,jmb
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb

      if(deb.ge.2) write(37,*) 'in selectList ',m7,vh1(3,3,3),filename
      
      if(me.ne.0) return
      
      delta = (imp-2)*(lmpf-2)*(kmp-2)/grtot

      if(deb.ge.4) write(37,*) 'delta ',delta,(lmpf-2)*(kmp-2),m7

      
      if(m7.eq.0) then
            open(fid,file=filename,form='formatted')	 
            if(deb.ge.4) write(37,*) 'delta ',delta,(lmpf-2)*(kmp-2)
      endif
      if(deb.ge.4) write(37,*) 'delta ',delta,(lmpf-2)*(kmp-2)

      do i = 1,10 
         do l = 1,3 
            do k = 1,3 
       	       n = (i-1)*9 + (l-1)*3 + k
      	       lst(n) = dabs(vh1(i+1,l+1,k+1))
            enddo
         enddo	                  
      enddo 
          
 727  format(90e20.8)
      if(me.eq.0) then
         write(fid,727) (lst(i),i=1,90) 	    
      endif     
           
      end      
      
      subroutine getMaxLongHarm(vh1,s27,m27)
      implicit none
      include 'part.pf'
      real*8 vh1(imp,lmpf,kmp),s27,mxh,lst(27)
      integer m27,i,l,k,n
      
      s27 = 0d0
      m27 = -1
      mxh = 0d0
      do i = 1,3 
         do l = 1,3 
            do k = 1,3 
       	       n = (i-1)*9 + (l-1)*3 + k
      	       lst(n) = dabs(vh1(i+1,l+1,k+1))
	       s27 = s27 + dabs(vh1(i+1,l+1,k+1))
	       if((lst(n).gt.mxh).and.(n.gt.1)) then
	          mxh = lst(n)
		  m27 = n
	       endif
            enddo
         enddo	                  
      enddo 
      
      end 
      
      subroutine fourierLong(vne,vni,velf)
      implicit none
      include 'part.pf'
      real*8 vne(imp,lmpf,kmp),vni(imp,lmpf,kmp),velf(imp,lmpf,kmp)
      real*8 ave,avi,avf,pi,maxe,maxi,maxf
      integer ixe,iye,ixi,iyi,ixf,iyf,i,l,k
      character*20 tr,tr1
      integer lp,jmf,jmi,n1,n2,m5,m27,m7
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7

      tr = 'furx'    
      write(tr1,'(i3.3)') m7 
      tr=tr(1:4)//tr1 
      tr=tr(1:7)//'.dat' 
      
      if(me.eq.0) then 
         open(92,file=tr,form='formatted')
      endif     
      
      do i = 2,imp-1
         ave = 0d0
         avi = 0d0
         avf = 0d0     
         
         maxe = 0d0 
         maxi = 0d0 
         maxf = 0d0 
         do l = 2,lmpf-1
            do k = 2,kmp-1
               ave = ave + dabs(vne(i,l,k))
               avi = avi + dabs(vni(i,l,k))
               avf = avf + dabs(velf(i,l,k))
               
               if(maxe.lt.dabs(vne(i,l,k))) then
                  maxe = dabs(vne(i,l,k))
                  ixe = l
                  iye = k
               endif
               if(maxi.lt.dabs(vni(i,l,k))) then
                  maxi = dabs(vni(i,l,k))
                  ixi = l
                  iyi = k
               endif
               if(maxf.lt.dabs(velf(i,l,k))) then
                  maxe = dabs(velf(i,l,k))
                  ixf = l
                  iyf = k
               endif
            enddo
         enddo
         ave = ave /(lmpf-2)/(kmp-2)
         avi = avi /(lmpf-2)/(kmp-2)
         avf = avf /(lmpf-2)/(kmp-2)
         
 782     format(i10,2i5,2e15.5,2i5,2e15.5,2i5,2e15.5)          
         if(me.eq.0) then 
            write(92,782) i,ixe,iye,ave,maxe,ixe,iye,ave,maxe,
     +                      ixe,iye,ave,maxe
         endif         
      enddo
         if(me.eq.0) then 
            close(92)
         endif         
      end

      real*8 function getMidHarmonic(v)
      implicit none
      include 'part.pf'
      integer imx,imy,imz,i,j,k,m7,l
      real*8 v(imp,lmpf,kmp),kr,kr2,dsqrt,Km,df,dmid

      df   = 0d0
      dmid = 0d0
      
         do i = 2,imp-1
            do j = 2,lmpf-1
               do k = 2,kmp-1
                  kr2 = (i-1)**2+(j-1)**2+(k-1)**2
                  kr  = dsqrt(kr2)
        
       	          df = df + v(i,j,k)**2
        
                  if((i.ge.8).and.(i.le.(imp/8)).and.
     +               (l.le.(lmpf/8)).and.(k.le.(kmp/8))) then
      		     dmid = dmid + v(i,j,k)**2
      		  endif         		  
       	       enddo	    
            enddo
         enddo 
          
      getMidHarmonic  = dmid/df
      end       
      
      
