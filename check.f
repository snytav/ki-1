
	!17.08.06 - ������� u0 ���.36 - ���������� ��� ���������� ������� ���������� ����
	!��������� ���������� �������, ���. ��������, ���������

	!���������� �������
      subroutine curbalance
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

      if(curbal.ne.1) return
c       CURRENT BALANCE ######################
      cfx = 0d0
      cfy = 0d0
      cfz = 0d0
      
c     ONLY PLASMA ELECTRONS      
      do j=1,jmf
         ps=1.d0/dsqrt(1.d0+
     +        beta0**2*(uf(j)*uf(j)+vf(j)*vf(j)+wf(j)*wf(j)))
         cfx = cfx + uf(j)*ps     
         cfy = cfy + vf(j)*ps     
         cfz = cfz + wf(j)*ps     
      end do
      cfx=cfx*amf
      cfy=cfy*amf
      cfz=cfz*amf

      cbx = 0d0
      cby = 0d0
      cbz = 0d0
      
c     BEAM ELECTRONS      
      do j=1,jmb
         ps=1.d0/dsqrt(1.d0+
     +        beta0**2*(ub(j)*ub(j)+vb(j)*vb(j)+wb(j)*wb(j)))
         cbx = cbx + ub(j)*ps     
         cby = cby + vb(j)*ps     
         cbz = cbz + wb(j)*ps     
      end do
      cbx=cbx*amb
      cby=cby*amb
      cbz=cbz*amb


      cix = 0d0
      ciy = 0d0
      ciz = 0d0
      
c     IONS      
      do j=1,jmi
         ps=1.d0/dsqrt(1.d0+
     +        beta0**2*(ui(j)*ui(j)+vi(j)*vi(j)+wi(j)*wi(j)))
         cix = cix + ui(j)*ps     
         ciy = ciy + vi(j)*ps     
         ciz = ciz + wi(j)*ps     
      end do
      cix=cix*ami
      ciy=ciy*ami
      ciz=ciz*ami

      
c#################################################
      vec(1) = cfx
      vec(2) = cfy
      vec(3) = cfz

      vec(4) = cbx
      vec(5) = cby
      vec(6) = cbz

      vec(7) = cix
      vec(8) = ciy
      vec(9) = ciz
      
      vec(10) = jmf

      call PARAreduceLMP2(vec)
      
      cfx = vec(1)
      cfy = vec(2)
      cfz = vec(3)
      
      cbx = vec(4)
      cby = vec(5)
      cbz = vec(6)

      cix = vec(7)
      ciy = vec(8)
      ciz = vec(9)

c######## ELECTRON TEMPERATURES #######################
      avfx=cfx/amf/vec(10)
      avfy=cfy/amf/vec(10)
      avfz=cfz/amf/vec(10)
      tempx=0d0
      tempy=0d0
      tempz=0d0
      do j=1,jmf
         ps=1.d0/dsqrt(1.d0+
     +        beta0**2*(uf(j)*uf(j)+vf(j)*vf(j)+wf(j)*wf(j)))
         tempx = tempx + (uf(j)*ps-avfx)**2
         tempy = tempy + (vf(j)*ps-avfy)**2 
         tempz = tempz + (wf(j)*ps-avfz)**2     

 492     format('tmp ',i10,9e15.5)
         if(deb.ge.4) then
            write(37,492) j,tempx,tempy,tempz,avfx,avfy,avfz,
     +                     (uf(j)*ps-avfx),(vf(j)*ps-avfy),
     +                      (wf(j)*ps-avfz) 
         endif


      enddo
      vec(1) = tempx
      vec(2) = tempy
      vec(3) = tempz

      vec(10) = jmf
ccccccccccccccccccccccccccccccccc
      do i = 1,1000
         vecb(i) = vec(i)
      enddo
      call PARAreduceLMP2(vec)
      do i  = 1,10
  56    format('reduceEnergy ',i10,2e30.20)
c        write(37,56) i,vecb(i),vec(i)
      enddo


c      call PARAreduceLMP2(vec)
      tempx = vec(1)/vec(10)
      tempy = vec(2)/vec(10)
      tempz = vec(3)/vec(10)
      


3012  format(14e20.10)
      if(me.eq.0) then
          if(nt.eq.1) then
             open(66,file='cbalance.txt',form='formatted')
          endif
          write(66,3012) nt*tau,cfx+cbx+cix,cfx,cfy,cfz,cbx,cby,cbz,
     +                          cix,ciy,ciz,tempx,tempy,tempz     
      endif      


c#############################################          
	end subroutine curbalance

!*******************************************************************

      
      subroutine ChainOut(nt)
      implicit none
      include 'part.pf'
      integer i,l,k,nt
      real*8 jxa
      real*8 jxb(imp,lmp,kmp),jxf(imp,lmp,kmp),ex(imp,lmp,kmp)
      real*8 ey(imp,lmp,kmp),ez(imp,lmp,kmp)
      real*8 hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp)
      real*8 jx(imp,lmp,kmp),jy(imp,lmp,kmp),jz(imp,lmp,kmp)
      character*40 tr,tr1
      common/h/ex,ey,ez
      common/g/hx,hy,hz    
      common/jtotal/jxf,jxb
      common/j/jx,jy,jz

      if(chnout.eq.1) then
         tr='jx' 
         write(tr1,'(i5.5)') nt 
         tr=tr(1:2)//tr1 
         tr=tr(1:7)//'.dat' 
        
         if(me.eq.0) then
            open(91,file=tr,form='formatted')
 320        format(3i5,7e25.10)        	 
            do i = 1,imp
               do l = 1,lmp
                  do k = 1,kmp
                     write(91,320) i,l,k,jxb(i,l,k),jxf(i,l,k),
     +                             jxf(i,l,k)-jxb(i,l,k),
     +                             jx(i,l,k),ex(i,l,k),
     +                             hy(i,l,k),hz(i,l,k)
                  enddo
               enddo
            enddo
            close(91)
	 endif
      endif

      if(curout.eq.1) then
         tr='jxal' 
         write(tr1,'(i3.3)') nt 
         tr=tr(1:4)//tr1 
         tr=tr(1:7)//'.dat' 
        
         if(me.eq.0) then
            open(91,file=tr,form='formatted')
 321        format(3i5,3e25.10)        	 
            do i = 2,imp-1
               do l = 2,lmp-1
                  do k = 2,kmp-1
                     write(91,321) i,l,k,jx(i,l,k),jy(i,l,k),jz(i,l,k)
                  enddo
               enddo
            enddo
            close(91)
	 endif
      endif

      end            
      
      subroutine poisson(y0p)
      implicit none
      include 'part.pf'
      real*8 ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp)
      real*8 jx(imp,lmp,kmp),jy(imp,lmp,kmp),jz(imp,lmp,kmp)
      real*8 hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp)
      real*8 dne(imp,lmp,kmp),dnb(imp,lmp,kmp),dni(imp,lmp,kmp)
      real*8 dn(imp,lmp,kmp),dive(imp,lmp,kmp),dn1(imp,lmp,kmp)
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      real*8 vec(1000),y0p,t
      real*8 tx,sst,h1,h2,h3,pi,amb,amf,ami,tau
      integer nt,ml,nt1,nt2,ns,lp,jmf,jmi,n1,n2,m5,m7,jmb,i,l,k
      common/a/tx,nt,ml,sst,ns,nt1,nt2		
      common/d/h1,h2,h3
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
      common/f/tau
      common/h/ex,ey,ez
      common/g/hx,hy,hz    
      common/j/jx,jy,jz
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/iff/xf,yf,zf,uf,vf,wf,amf
      common/den1/dn1


c      call densSort(xf,yf,zf,y0p,jmf,amf,0,dne,0)
c      call densSort(xb,yb,zb,y0p,jmb,amb,0,dnb,0)
c      call densSort(xi,yi,zi,y0p,jmi,ami,0,dni,0)
      
      do i = 1,imp
         do l = 1,lmp
            do k = 1,kmp
               dn(i,l,k)=dne(i,l,k)+dni(i,l,k)+dnb(i,l,k)
            enddo
         enddo
      enddo
      if(deb.ge.1) then
         write(37,*) 'jmi,jmf,jmb ',jmi,jmf,jmb
         write(37,*) 'ami,amf,amb ',ami,amf,amb
      endif
c      write(37,*) 'jXYZ 11_3 ',jx(11,11,11),jy(11,11,11),jz(11,11,11)
      
      do i = 2,imp-1
         do l = 2,lmp-1
            do k = 2,kmp-1
               dive(i,l,k)=(ex(i+1,l,k)-ex(i,l,k))/h1+
     +                     (ey(i,l+1,k)-ey(i,l,k))/h2+
     +                     (ez(i+1,l,k)-ez(i,l,k))/h3 
     
 237           format('div',3i3,3e8.1,' d ',3e8.1,' e ',
     +                           3e8.1,' h ',3e8.1,' j ',3e8.1)
               if((dabs(dive(i,l,k)-dn(i,l,k)).gt.1e-3)
     +            .and.(deb.ge.2)) then
                  write(37,*) i,l,k,'==========================='
                  write(37,237) i,l,k,dive(i,l,k),
     +                       dabs(dive(i,l,k)-dn(i,l,k)),
     +                       dn(i,l,k),
     +                       dne(i,l,k),dni(i,l,k),dnb(i,l,k),
     +                       ex(i,l,k),ey(i,l,k),ez(i,l,k),
     +                       hx(i,l,k),hy(i,l,k),hz(i,l,k),
     +                       jx(i,l,k),jy(i,l,k),jz(i,l,k)
                  write(37,237) i+1,l,k,dive(i+1,l,k),
     +                       dabs(dive(i+1,l,k)-dn(i+1,l,k)),
     +                       dn(i+1,l,k),
     +                       dne(i+1,l,k),dni(i+1,l,k),dnb(i+1,l,k),
     +                       ex(i+1,l,k),ey(i+1,l,k),ez(i+1,l,k),
     +                       hx(i+1,l,k),hy(i+1,l,k),hz(i+1,l,k),
     +                       jx(i+1,l,k),jy(i+1,l,k),jz(i+1,l,k)
                  write(37,237) i,l+1,k,dive(i,l+1,k),
     +                       dabs(dive(i,l+1,k)-dn(i,l+1,k)),
     +                       dn(i,l+1,k),
     +                       dne(i,l+1,k),dni(i,l+1,k),dnb(i,l+1,k),
     +                       ex(i,l+1,k),ey(i,l+1,k),ez(i,l+1,k),
     +                       hx(i,l+1,k),hy(i,l+1,k),hz(i,l+1,k),
     +                       jx(i,l+1,k),jy(i,l+1,k),jz(i,l+1,k)
                  write(37,237) i,l,k+1,dive(i,l,k+1),
     +                       dabs(dive(i,l,k+1)-dn(i,l,k+1)),
     +                       dn(i,l,k+1),
     +                       dne(i,l,k+1),dni(i,l,k+1),dnb(i,l,k+1),
     +                       ex(i,l,k+1),ey(i,l,k+1),ez(i,l,k+1),
     +                       hx(i,l,k+1),hy(i,l,k+1),hz(i,l,k+1),
     +                       jx(i,l,k+1),jy(i,l,k+1),jz(i,l,k+1)
                  write(37,*) '==========================='
               endif
            enddo
         enddo
      enddo

      if(nt.gt.1) then
         do i = 1,imp-1
            do l = 1,lmp-1
               do k = 1,kmp-1
                  t = (dn(i,l,k) - dn1(i,l,k))/tau+
     +                (jx(i+1,l,k)-jx(i,k,l))/h1+
     +                (jy(i,l+1,k)-jy(i,k,l))/h2+
     +                (jz(i,l,k+1)-jz(i,k,l))/h3
  295 format('cont ',3i5,4e9.1,' divj ',e8.1,' dn/dt ',e9.1,2e15.5)
                  if(deb.ge.2) then
                     write(37,295) i,l,k,t,
     +                          jx(i,l,k),jy(i,l,k),jz(i,l,k),
     +                          (jx(i+1,l,k)-jx(i,k,l))/h1+
     +                          (jy(i,l+1,k)-jy(i,k,l))/h2+
     +                          (jz(i,l,k+1)-jz(i,k,l))/h3,
     +                          (dn(i,l,k) - dn1(i,l,k))/tau,
     +                          dn1(i,l,k),dn(i,l,k)
                  endif
               enddo
            enddo
         enddo
      endif
      
      do i = 1,imp
         do l = 1,lmp
            do k = 1,kmp
               dn1(i,l,k) = dn(i,l,k)
            enddo
         enddo
      enddo
c      write(37,*) 'end poisson'
      end subroutine poisson

      subroutine traceparticle(y0p,sx,sy,sz)
      implicit none
      include 'part.pf'
      integer nt
      real*8 sx,sy,sz
      real*8 ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp)
      real*8 dne(imp,lmp,kmp),dnb(imp,lmp,kmp),dni(imp,lmp,kmp)
      real*8 dn(imp,lmp,kmp),dive(imp,lmp,kmp)
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      real*8 vec(1000),y0p
      real*8 tx,sst,h1,h2,h3,pi,amb,amf,ami,tau
      integer ml,nt1,nt2,ns,lp,jmf
      integer jmi,n1,n2,m5,m7,jmb,i,l,k
      common/a/tx,nt,ml,sst,ns,nt1,nt2		
      common/d/h1,h2,h3
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
      common/f/tau
      common/h/ex,ey,ez
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/iff/xf,yf,zf,uf,vf,wf,amf

      if((deb.lt.1).or.(jdeb.lt.1)) then
         return
      endif

      call poisson(y0p)

      if((nt.eq.1).and.(me.eq.0)) then
         open(34,file='tracepartE.dat',
     +         form='formatted',action='write',position='append')
c         uf(jdeb) = 10d0
      endif
c      if(deb.ge.2) then
 283     format(2i10,9e15.5)

      if(me.eq.0) then
         write(34,283) nt,jdeb,xf(jdeb),yf(jdeb),zf(jdeb),
     +                         uf(jdeb),vf(jdeb),wf(jdeb),
     +                         sx,sy,sz
      endif
      
      end subroutine traceparticle

