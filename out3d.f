      subroutine tempFlow3d(xf,yf,zf,uf,vf,wf,y0,jm)
      implicit none
      include 'part.pf'
      real*8 uf(jmp),vf(jmp),wf(jmp),h1,h2,h3,y0,minT,maxT	
      real*8 Mx(imp,lmp,kmp),Dx(imp,lmp,kmp),s2,s4,s6
      real*8 My(imp,lmp,kmp),Mz(imp,lmp,kmp),
     +       Dy(imp,lmp,kmp),Dz(imp,lmp,kmp)
      real*8 Fx(imp,lmp,kmp),globFx(imp,lmpf,kmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),pi,tt,tfl2,tfl3,tfl0,tfl
      real*8 vec(lmp2)
      character*20 tr,tr1
      integer i,l,k,j,Num(imp,lmp,kmp),jm
      integer lp,jmf,jmi,n1,n2,m5,m7
      common/d/h1,h2,h3
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7
      
      if(deb.ge.4) write(37,*) 'in temp3d'
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
      
      do j = 1,jm
         s2 = xf(j)/h1
         i  = idint(s2+1.5d0)
         s4 = (yf(j)-y0)/h2
         l  = idint(s4+1.5d0)
         s6 = zf(j)/h3
         k  = idint(s6+1.5d0)
	 if(deb.ge.7) write(37,*) i,l,k,j,xf(j),yf(j),y0
	 
         Num(i,l,k) = Num(i,l,k) + 1
         Mx(i,l,k)  = Mx(i,l,k)  + uf(j)  	  
         My(i,l,k)  = My(i,l,k)  + vf(j)  	  
         Mz(i,l,k)  = Mz(i,l,k)  + wf(j)  	  
      enddo       
      if(deb.ge.4) write(37,*) 'after Mxyz'
       
      do i = 1,imp
         do l = 1,lmp
            do k = 1,kmp
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
         
         Dx(i,l,k) = Dx(i,l,k) + (uf(j) - Mx(i,l,k))**2	 	 
         Dy(i,l,k) = Dy(i,l,k) + (vf(j) - My(i,l,k))**2	 	 
         Dz(i,l,k) = Dz(i,l,k) + (wf(j) - Mz(i,l,k))**2	 	 
      enddo       
      if(deb.ge.4) write(37,*) 'after Dxyz'
c      call PARAFinal
c      stop    
      do i = 1,imp
         do l = 1,lmp
            do k = 1,kmp
               if(Num(i,l,k).gt.0) then
                  Dx(i,l,k) = Dx(i,l,k)/Num(i,l,k)
                  Dy(i,l,k) = Dy(i,l,k)/Num(i,l,k)
                  Dz(i,l,k) = Dz(i,l,k)/Num(i,l,k)
                  Dx(i,l,k) = (Dx(i,l,k)+Dy(i,l,k)+Dz(i,l,k))/3.0d0		  
      	       endif	       
      	    enddo
      	 enddo
      enddo     

      do i=1,imp-1
         do l=1,lmp-1
c            DX(i,l,2)     = DX(i,l,2)+DX(i,l,kmp-1)
c            DX(i,l,kmp-1) = DX(i,l,2)
         enddo
      enddo

      do i=1,imp-1
         do l=1,kmp-1
c            DX(i,2,l)     = DX(i,2,l)+DX(i,lmp-1,l)
c            DX(i,lmp-1,2) = DX(i,2,l)
         enddo
      enddo
      do i=1,kmp-1
         do l=1,lmp-1
c            DX(2,l,i)     = DX(2,l,i)+DX(imp-1,l,i)
c            DX(imp-1,l,i) = DX(2,l,i)
         enddo
      enddo

      	 tfl2 = 0d0
      	 tfl3 = 0d0
         
         do i = 2,imp-1
            do l = 2,lmp-1
               do k = 2,kmp-1
                  tt = dabs(Dx(i,l,k))*
     +                dsqrt(Mx(i,l,k)**2+My(i,l,k)**2+Mz(i,l,k)**2)
                      Fx(i,l,k) = tt
                      if(deb.ge.4) write(37,*) 'curfl ',i,l,k,tt,tfl0
               enddo    
      	    enddo
      	 enddo
      
      call PARAgather(Fx,globFx)

      if(me.eq.0) then
         
         if((out3d.eq.1).and.(me.eq.0)) then       
            tr = 'flow'    
            write(tr1,'(i3.3)') m7 
            tr=tr(1:4)//tr1 
            tr=tr(1:7)//'.dat' 
      
            open(91,file=tr,form='formatted')

            do i = 2,imp-1
               do j = 2,lmpf-1
                  do k = 2,kmp-1
  174                format(4e15.5)	    
c                     if(globFx(i,l,k).gt.0) then
                        write(91,174) h1*(i-1.5d0),h2*(j-1.5d0),
     +                        h3*(k-1.5d0),globFx(i,l,k)
!                     ,Dx(i,l,k),
!     +               dsqrt(Mx(i,l,k)**2+My(i,l,k)**2+Mz(i,l,k)**2)
c                     endif 
       	          enddo	    
               enddo
            enddo 
            close(91)
         endif	    
           
      endif
      
      
           
      end
