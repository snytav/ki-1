c      include 'Energy.f'
c      include 'Vvod.f'
!      include 'Dens.for'

      implicit real*8(a-h,o-z)
      include 'part.pf'
      include 'mass.par'
!      include 'mpif.h'
      real*8 ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp),
     *hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp),
     *jx(imp,lmp,kmp),jy(imp,lmp,kmp),jz(imp,lmp,kmp),
     *jp(imp,lmp,kmp),f(30),
     *qx(imp,lmp,kmp),qy(imp,lmp,kmp),qz(imp,lmp,kmp)
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      real*8 ufp(jmp),up(jmp),vp(jmp),wp(jmp)
      real*8 ei,eb,ef,eh,ee,ew,ei0,eb0,ef0,eh0,ee0
      real*8 mee,nee
      real*8 ni,ne,Te,hx0,hx00,ni0
      real*8 pp(imp,lmp,kmp)
      real*8 disp_old,gamma_t
      real*8 Mw(imp)	
      real*8 Dw(imp) 
      real*8 q, ksi
      real*8 xp,yp,zp,px,py,pz,qp,ambp
      integer jpp,rank,m7,m33,m57,nt3
      character*2 st,ci
      character*10 str, str1 
      character*40 st1,tr,tr1
      integer vnum,i1,i2,i33,i4	
      common/a/tx,nt,ml,sst,ns,nt1,nt2
      common/b/xm,ym,zm
      common/c/im,lm,km
      common/d/h1,h2,h3
      common/g/hx,hy,hz
      common/f/tau,c1,c2,c3
      common/h/ex,ey,ez
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/iff/xf,yf,zf,uf,vf,wf,amf
      common/o/qx,qy,qz
      common/j/jx,jy,jz
      common/p/p_tau
      common/c1/alpha,t0,tD
      common/c2/ni,ne,Te,hx0
      common/cc/ncc
      common/et/ei0,ef0,eh0,ee0,ew0
      common/den0/pp
      common/dispold/disp_old,gamma_t	    
      common/dwmw/Mw,Dw
      common/start_p/Te0,ksi,ak0
      common/vvv/vnum,i1,i2,i33,i4	
       	
       	
c      real*8, pointer :: atrib    ! (3*jmp*attr)
c      common/ctrl/atrib


      call PARAinit
      if(deb.ge.1) write(37,*) 'ONE'

      call allocControlAttribute

c      allocate(atrib(3*jmp*attr))
      
      if(deb.ge.2) then
         write(37,*) '0 ',nt,jdeb,uf(jdeb)
      endif
      
      if(outf.ne.0) then	
         open(25,file='laser018.lst',form='formatted')
         open(9,file='lasem18.dat',form='formatted')
         open(18,file='nn18.dat',form='formatted')
         open(66,file='ewfour.dat',form='formatted')
      endif    	 
      open(12,file='lstart18.dat',form='formatted')

      read(12,302) nt1	!���������� ���������� ������
      read(12,302) nt2	!���������� ������� ������
      read(12,302) nt3	!���������� ������� ������
      read(12,302) ml		!����� ������ �� ����
      read(12,301) tau	!��������� ���

      if(outf.ne.0) then
	if(ml.le.1)then
		ncc=0
	else
	  read(9,306)ncc
	  st=ci(ncc)
        st1='laser'//st(1:2)//'.lst'		!���� ��� ����� ����������
        open(10,file=st1,form='unformatted')
        end if
      endif	

  301 format(e10.3,40x)
  302 format(i10,40x) 
 9001 format(i10,i4,f14.6,i4)
      
      pi=3.14159265358979d0
      
      
      
      call start	!��������� ��������
      print*,'after start ',me
      call writeParticleList(xf,yf,zf,uf,vf,wf,nt,jmf,'strt','el')
      call writeParticleList(xi,yi,zi,ui,vi,wi,nt,jmi,'strt','in')
      call writeParticleList(xb,yb,zb,ub,vb,wb,nt,jmb,'strt','bm')
! -------------------------------- my
      if(beamf.eq.2) then
         if(deb.ge.1) write(37,*) 'ONE'
         call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
         if(deb.ge.1) write(37,*) 'ONE'
         if (rank.eq.0) then
             xp=0.d0
             yp=5.d0
             zp=5.d0
             px=0.699d0/(1-0.699d0**2)
             py=0.d0
             pz=0.d0
             qp=-1.896d0*100000/melb
             ambp=1.d0
             jpp=1
         open(12321,file='particle.dat',form='formatted')
         nt=0
         write(12321,*) 'begin'
         write(12321,987) nt*tau,xp,yp,zp,px,py,pz
  987    format(7e20.10)
          
         else
             qp=0.d0
             amb=0.d0
             jpp=0.d0
             xp=0.d0
             yp=0.d0
             zp=0.d0
             px=0.d0
             py=0.d0
             pz=0.d0
         endif
      endif
!---------------------------------

      if(deb.ge.2) then
         write(37,*) 'c ',nt,jdeb,uf(jdeb)
      endif
      
      
      
      y0p = ym*meh
      
      if(deb.ge.2) write(37,*) 'a st '
c      stop
c      if(m5.gt.2) stop
      
	!------------------------------------------------
	call open_f		!�������� ������
	!------------------------------------------------
      if(deb.ge.2) write(37,*) 'a open '
      im1=imp-1
      lm1=lmp-1
      km1=kmp-1
      im2=imp
      lm2=lmp
      km2=kmp
      itot=is+60*(imm+60*ih)
	j1=0
	j2=jm

      m7  = 0
      m5  = 0
      m57 = 0
      m33 = nt3
      call output(y0p,nt1,m33,nt3,0)
     
      if(deb.ge.2) write(37,*) 'b loop  '
c      stop
      
      call printTime(0)
      m7 = 1
      m33= 1
      do 1 while(m7.le.nt2)
      m5 = 1
      m57 = 1
      do 2 while(m5.le.nt1)
      if(deb.ge.2) then
         write(37,*) 'D ',nt,jdeb,uf(jdeb)
      endif
            
      call timeBegin(3)
      nt=nt+1
      tx=tx+tau
	ss1=(nt-1.d0)/sst
	ns1=ss1
	ss2=ss1-ns1
 9002 format(/8x,i10,1x,f7.4)
         call timeBegin(2)
         call writeallbinaryoutput(0,nt)
 1853 format('writeallbinaryoutput ',i5,3e30.20)
         if(deb.ge.1) write(37,1853) 0,xf(1),yf(1),zf(1)
         call emh1(nt)
         call writeFieldArrays(ex,ey,ez,hx,hy,hz,
     +                             qx,qy,qz,jx,jy,jz,
     +                             nt,'emh1')
         call timeEnd(2)
         call writeallbinaryoutput(50,nt)
         
         do 101 k=1,kmp
         do 101 l=1,lmp
         do 101 i=1,imp
            jx(i,l,k)=0.d0
            jy(i,l,k)=0.d0
            jz(i,l,k)=0.d0
  101    continue
      if(deb.ge.2) then
         write(37,*) 'E ',nt,jdeb,uf(jdeb)
      endif
         if(deb.ge.2) write(37,*) 'b M_D ',m5
c	 stop
c        if(m5.gt.2) stop
	
	call M_D(uf,vf,wf,1,jmf)		!���������� ����� ���������
							! ��� ���������� ����
	call M_D(ui,vi,wi,2,jmi)		!���������� ����� ���������

        if(deb.ge.2) write(37,*) 'b mo ',m5
c       if(m5.gt.2) stop
      if(deb.ge.2) then
         write(37,*) 'F ',nt,jdeb,uf(jdeb)
      endif

      call timeBegin(1)

      if(deb.ge.1) write(37,313) -1,xm,ym,zm,jmi,jmf,jmb

      if(beamf.eq.1) then
c         q=0.d0
    	 q=-1.d0/melb
         if(deb.ge.2) write(37,*) 'nt,amb,q,jmb ',nt,amb,q,jmb
c         stop
         call writeParticleList(xb,yb,zb,ub,vb,wb,nt,jmb,'afm0','bm')
         call sort(xb,yb,zb,ub,vb,wb,jmb)
         call writeParticleList(xb,yb,zb,ub,vb,wb,nt,jmb,'sort','bm')

         call writeFieldArrays(ex,ey,ez,hx,hy,hz,
     +                             qx,qy,qz,jx,jy,jz,
     +                             nt,'bmob')

         call move3(nt,xb,yb,zb,ub,vb,wb,amb,q,jmb,'bm')	!�������� ���������� ����
         
         call writeFieldArrays(ex,ey,ez,hx,hy,hz,
     +                             qx,qy,qz,jx,jy,jz,
     +                             nt,'amob')
         



         call writeParticleList(xb,yb,zb,ub,vb,wb,nt,jmb,'afm1','bm')
      endif
      call timeEnd(1)	
      if(deb.ge.1) write(37,313) 0,xm,ym,zm,jmi,jmf,jmb
      
! -------------------------  my - one particle -------------------------------
c      if(beamf.eq.2) then
c         call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
c         if (rank.eq.0) then
c           qp=-1.896d0*100000/melb
c           ambp=1.d0
c           call move3(nt,xp,yp,zp,px,py,pz,ambp,qp,jpp)      !�������� ���������� ����
c  160      format(7e20.10)
c                  write(12321,160) nt*tau,xp,yp,zp,px,py,pz
c         else
c           qp=0.d0
c           amb=0.d0
c           jpp=0.d0
c         endif
c      endif
c       q=0.d0
       q=alpha
      if(deb.ge.6) write(37,*) 'bmo ami,amf ',ami,amf
      if(deb.ge.2) then
         write(37,*) 'G ',nt,jdeb,uf(jdeb)
      endif
        	
      call timeBegin(1)	
c         write(37,564) ex(2,2,2),ey(2,2,2),ez(2,2,2)

         call writeallbinaryoutput(100,nt)
         if(deb.ge.1) write(37,1853) 100,xf(1),yf(1),zf(1)

      call writeFieldArrays(ex,ey,ez,hx,hy,hz,
     +                             qx,qy,qz,jx,jy,jz,
     +                             nt,'beme')

      call eme(nt)	!���������� ������������� ����
!      stop

      call writeFieldArrays(ex,ey,ez,hx,hy,hz,
     +                             qx,qy,qz,jx,jy,jz,
     +                             nt,'aeme')



 2233 format('ex ',e11.3)
      if(deb.ge.1) write(37,2233) ex(1,2,2)
         call writeallbinaryoutput(150,nt)
      
      q = -1.0
      if(deb.ge.1) write(37,*) 'b-pri yi(2) ',yi(2)
      
c      call printallpart(nt)

          call output(y0p,nt,m33,nt3,1)
      call writeallbinaryoutput(200,nt)
         if(deb.ge.1) write(37,1853) 200,xf(1),yf(1),zf(1)
         call writeParticleList(xb,yb,zb,ub,vb,wb,nt,jmb,'afm2','bm')

      q=alpha
      call writeParticleList(xi,yi,zi,ui,vi,wi,nt,jmi,'bmoi','in')
      call move3(nt,xi,yi,zi,ui,vi,wi,ami,q,jmi,'in')	!�������� �����
     
      call writeParticleList(xi,yi,zi,ui,vi,wi,nt,jmi,'amoi','in')
      call writeFieldArrays(ex,ey,ez,hx,hy,hz,
     +                             qx,qy,qz,jx,jy,jz,
     +                             nt,'amoi')


         call writeallbinaryoutput(250,nt)
          call output(y0p,nt,m33,nt3,2)
      if(ctrl_attr.eq.1) then
         write(37,*) 'ATTRIBUTES ',jmi,jmf,jmb
         call saveControlAttributesToFile(nt)
      endif
         call writeParticleList(xb,yb,zb,ub,vb,wb,nt,jmb,'afm3','bm')

          
      if(jx_four.eq.1) then
         call four3D(jx,dni2,imp,lmp,kmp,nt,
     +	 'jxmo','frjx',0d0,0d0,0d0)
      endif

      if(deb.ge.1) write(37,*) 'amo yi(2) ',yi(2)
c      call nonpercur(y0p,nt1,m33,nt3)
      if(deb.ge.2) then
         write(37,*) 'H ',nt,jdeb,uf(jdeb)
      endif
c        q=0.d0           
      	q=-1.d0/mel

      if(deb.ge.2) then
         write(37,*) 'mel ',mel
      endif
         call writeParticleList(xb,yb,zb,ub,vb,wb,nt,jmb,'afm4','bm')

      call writeParticleList(xf,yf,zf,uf,vf,wf,nt,jmf,'afm4','el')
      !stop

      write(37,*) 'afm5 jmf before ',jmf
      call move3(nt,xf,yf,zf,uf,vf,wf,amf,q,jmf,'el')	!�������� ���������� ����
      call writeParticleList(xf,yf,zf,uf,vf,wf,nt,jmf,'af45','el')
!      call PARAfinal
!      stop

      write(37,*) 'afm5 jmf after ',jmf
      if(deb.ge.1) write(37,1853) 270,xf(1),yf(1),zf(1)
         call writeallbinaryoutput(270,nt)
         call writeParticleList(xb,yb,zb,ub,vb,wb,nt,jmb,'afm5','bm')

      !call addbeam
      call writeParticleList(xb,yb,zb,ub,vb,wb,nt,jmb,'aadd','bm')

         call writeallbinaryoutput(275,nt)
      !endif

      call writeParticleList(xf,yf,zf,uf,vf,wf,nt,jmf,'afm5','el')

 313  format('clean lx,ly,lz ',i3,3e15.5,3i10)
      if(deb.ge.1) write(37,313) 1,xm,ym,zm,jmi,jmf,jmb
      call cleanparticles(xi,yi,zi,ui,vi,wi,
     +                    0d0,xm,0d0,ym,0d0,zm,jmi,nt,'in','move')
      if(deb.ge.1) write(37,313) 2,xm,ym,zm,jmi,jmf,jmb
      call cleanparticles(xf,yf,zf,uf,vf,wf,
     +                    0d0,xm,0d0,ym,0d0,zm,jmf,nt,'el','move')
      if(deb.ge.1) write(37,313) 3,xm,ym,zm,jmi,jmf,jmb

      call writeParticleList(xb,yb,zb,ub,vb,wb,nt,jmb,'becl','bm')

      call cleanparticles(xb,yb,zb,ub,vb,wb,
     +                    0d0,xm,0d0,ym,0d0,zm,jmb,nt,'bm','move')
      if(deb.ge.1) write(37,313) 4,xm,ym,zm,jmi,jmf,jmb

      call writeParticleList(xf,yf,zf,uf,vf,wf,nt,jmf,'afmo','el')
      call writeParticleList(xi,yi,zi,ui,vi,wi,nt,jmi,'afmo','in')
      call writeParticleList(xb,yb,zb,ub,vb,wb,nt,jmb,'afmo','bm')






      call timeEnd(1)	
      if(deb.ge.2) then
         write(37,*) 'I...jm,m5,nt1-1 ',nt,jdeb,uf(jdeb),jm,m5,nt1-1
      endif

      if(deb.ge.2) then
         write(37,*) 'J1 ',nt,jdeb,uf(jdeb)
         write(37,*) 'px(2,2,2) after move all',jx(2,2,2)
      endif
      

      if(deb.ge.2) then
         write(37,*) 'J2 ',nt,jdeb,uf(jdeb)
      endif


      if(m5.eq.nt1-1)then
      if(deb.ge.2) then
         write(37,*) 'J3 ',nt,jdeb,uf(jdeb)
      endif
      
c         do 103 i=1,jm
c            up(i)=uf(i)
c            vp(i)=vf(i)
c            wp(i)=wf(i)
c  103    continue
      end if
      if(deb.ge.2) write(37,*) 'b cur ',m5
c      if(m5.gt.2) stop
      if(deb.ge.2) then
         write(37,*) 'J ',nt,jdeb,uf(jdeb)
      endif

!      do 8 k=1,km2
!      do 8 l=1,lm2
!      jx(1,l,k)=jx(1,l,k)+jx(im1,l,k)
!      jx(im1,l,k)=jx(1,l,k)
!      jy(2,l,k)=jy(2,l,k)+jy(im2,l,k)
!      jy(im1,l,k)=jy(im1,l,k)+jy(1,l,k)
!      jz(2,l,k)=jz(2,l,k)+jz(im2,l,k)
!      jz(im1,l,k)=jz(im1,l,k)+jz(1,l,k)
!      jx(im2,l,k)=jx(2,l,k)
!      jy(1,l,k)=jy(im1,l,k)
!      jy(im2,l,k)=jy(2,l,k)
!      jz(1,l,k)=jz(im1,l,k)
!      jz(im2,l,k)=jz(2,l,k)
!    8 continue
    
      if(deb.ge.2) then
         write(37,*) 'px(2,2,2) bnd1 ',jx(2,2,2),jx(2,2,km2)
         write(37,*) 'K ',nt,jdeb,uf(jdeb)
      endif

      do 10 l=1,lm2
      do 10 i=1,im2
      jx(i,l,2)=jx(i,l,2)+jx(i,l,km2)
      jx(i,l,km1)=jx(i,l,km1)+jx(i,l,1)
      jy(i,l,2)=jy(i,l,2)+jy(i,l,km2)
      jy(i,l,km1)=jy(i,l,km1)+jy(i,l,1)
      jz(i,l,1)=jz(i,l,1)+jz(i,l,km1)
      jz(i,l,km1)=jz(i,l,1)
      jx(i,l,1)=jx(i,l,km1)
      jx(i,l,km2)=jx(i,l,2)
      jy(i,l,1)=jy(i,l,km1)
      jy(i,l,km2)=jy(i,l,2)
      jz(i,l,km2)=jz(i,l,2)
   10 continue
       if(deb.ge.2) write(37,*) 'b send ',m5
c      if(m5.gt.2) stop
      if(deb.ge.2) then
         write(37,*) 'px(2,2,2) bnd2 ',jx(2,2,2)
         write(37,*) 'L ',nt,jdeb,uf(jdeb)
      endif

      
c************************************
      call PARAreduceP(jx,jy,jz)
      
      
      call PARAsendP1(jx,jy,jz)
         
      call PARAsendP2(jx,jy,jz)
      if(deb.ge.2) then
         write(37,*) 'px(2,2,2) send ',jx(2,2,2)
      endif
c      write(37,*) 'p jXYZ 11_3 ',jx(11,11,11),jy(11,11,11),jz(11,11,11)

c*********************************
      
      if(deb.ge.2) write(37,*) 'b field ',m5
!      stop
c      if(m5.gt.2) stop

      call timeBegin(2)
      if(deb.ge.1) write(37,*) 'b fields ',yi(2)
         call writeallbinaryoutput(400,nt)

      call emh2(nt)	!���������� ���������� ����
      call writeFieldArrays(ex,ey,ez,hx,hy,hz,
     +                             qx,qy,qz,jx,jy,jz,
     +                             nt,'emh2')

         call writeallbinaryoutput(500,nt)
      !call PARAfinal
      !top
      
      if(deb.ge.2) then
         write(37,*) 'px(2,2,2) b eme  ',jx(2,2,2)
      endif
      
      !call eme(nt)	!���������� ������������� ����
      call writeFieldArrays(ex,ey,ez,hx,hy,hz,
     +                             qx,qy,qz,jx,jy,jz,
     +                             nt,'emee')

      call writeParticleList(xb,yb,zb,ub,vb,wb,nt,jmb,'bbma','bm')
      call addbeam
      call writeParticleList(xb,yb,zb,ub,vb,wb,nt,jmb,'abma','bm')
      call writeParticleList(xb,yb,zb,ub,vb,wb,nt,jmb,'abma','el')
      
      call timeEnd(2)
      call TimeEnd(3) 
      call printTime(nt)
      call writeallbinaryoutput(600,nt)

      !stop

cc my energy
      if(m57.eq.1000)then
      call energy	!���������� �������
      m57=0
      endif
      m57=m57+1
cc     
      if(nt.eq.ntd) then
         nt2=nt2*nt1/ntd1
         nt1=ntd1
         m7=ntd/ntd1+ntd1
         m5 = 1
      endif
      m5=m5+1   


    2 continue
      if(m7.lt.1) goto 1		!1->max
 9005 format(/8x,'nt=',i10,10x,'�६�=',f7.4)
          call output(y0p,nt1,m33,nt3,3)
    
      m7=m7+1      
      m33=m33+1
    1 continue
    
!------------------------------------- my -----------------
      if(beamf.eq.2) then
         call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
         if (rank.eq.0) then
          close(12321)
         endif
      endif
!----------------------------------------------------------

      
      if((me.eq.0).and.(ons.eq.0)) then
         close(94)
      endif	 
	!----------------------------------------
      if(ml.ne.1.and.ml.ne.3) goto 5
	call vivod			!����� � ����
  	!----------------------------------------
    5 continue


  366 format(i8,'	',e16.9)
  367 format(i8,'	',i8,'	',e16.9)

      if(outf.ne.0) then
         close(10)
  	 close(11)
         close(12)
         close(18)
 3600 format(2x,'calculation time=',i8)     
      
         close(25)
         close(9)
      endif
  300 format(i4,'   ',e12.4)           
  304 format(e14.6)
  306 format(i4)
      call PARAfinal
      stop
      close(22)
      close(21)
      
      end

!------------------------------------------------------------------
	!������� ��������� ��������
      subroutine start
!	USE DFLIB 
      implicit real*8(a-h,o-z)
      include 'part.pf'
	include 'mass.par'
c      parameter(imp=42,lmp=22,kmp=3,jmp=224000)
      real*8 qx(imp,lmp,kmp),qy(imp,lmp,kmp),qz(imp,lmp,kmp),
     *hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp),
     *jx(imp,lmp,kmp),jy(imp,lmp,kmp),jz(imp,lmp,kmp),
     *ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp),f(30)
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
	real*8 ak(jmp),trg,g05,vb0,vf01,vf02,pinv1,pinv2,gb0
      real*8 mee,nee,mii,ne0,kl,kl_m,tau_e
      real*8 ni,ne,Te,hx0,lx,ly,lz,hx00,ni0
      real*8 dT,dn,dE,du,kk,const1,gam
      real*8 dln,eps,dx1,dx2,eps1,ak0
      real*8 ksi,termx
      real*8 N_D, N_D1
      real*8 lambda,lambda_m,gamma_t
      real*8 xpmin,xpmax,ypmin,ypmax,zpmin,zpmax
      real*8 xbmin,xbmax,ybmin,ybmax,zbmin,zbmax
      integer vnum,i1,i2,i33,i4,irg,nt1,nt2
      real*8 tf0,Tbb,ux(jmp),uy(jmp),uz(jmp),vy,vz
      common/a/tx,nt,ml,sst,ns,nt1,nt2
      common/b/lx,ly,lz
      common/c/im,lm,km
      common/d/h1,h2,h3
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
      common/g/hx,hy,hz
	common/f/tau,c1,c2,c3
      common/j/jx,jy,jz
      common/h/ex,ey,ez
      common/o/qx,qy,qz
	common/p/p_tau
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/iff/xf,yf,zf,uf,vf,wf,amf
      common/c1/alpha,t0,tD
      common/c2/ni,ne,Te,hx0
	common/cc/ncc
	common/et/ei0,ef0,eh0,ee0,ew0
	common/dd/dT,dn,dE,du,kk,gam,tau_e
	common/dispold/disp_old,gamma_t
	common/start_p/Te0,ksi,ak0
      common/vvv/vnum,i1,i2,i33,i4	
      
      !!write(37,*) 'start in'
      !!print *,'START '
 9006 format(2x,'pa�otaet �o��po�pamma  start')
	eps=1.d-10
      pi=3.14159265358979d0
	!------------------------------------
      if(ml.le.1) goto 1
        if(outf.ne.0) then
	   call vvod			!���� ��������� �������� �� �����
	endif   
      goto 2
    1 continue

!---------------------------------------------------------------
	ns=0
	!���������
!      k=1.3807d-16       !erg/K
      zk=1.6d-12          !1�� � ���
      ee=4.8032d-10       !����� ��������� ��. ����
      mee=mel*9.1094d-28      !����� ��������� � �������
      mii=mi*9.1094d-28      !����� ����� � ������� ��� ����� (������������!)
!      mii=1.6726d-24      !�����  ����� � �������
      cc=2.9979d+10       !�������� ����� � ��/�
	
      !!write(37,*) 'start 11'
	
      im = imp-2
      km = kmp-2
      lm = lmp-2
      read(12,302) lp		!���������� ������ � ������ ����� lp*lp*lp

!----------------------------------------------------------------
      read(12,301)ni0     !��������� ����� (1/cm^3)
!      read(12,301)nb0    !��������� �����  (1/cm^3)
!      read(12,301)ne0    !��������� ����������
      ne0=ni0	!-nb0		!��������� ����������
	!!print*,'ne0',ne0
      read(12,301)Te0     !����������� ����������   (eV)
      read(12,301)hx00    !��������� ���������� ��������� ����  (gauss)
	!!print*,'hx00',hx00
!      read(12,301)gamma  !�������� �����  
      read(12,301)lx      !������ ������� �� x (� �������� �����)
      read(12,301)ly      !������ ������� �� y (� �������� �����)
      read(12,301)lz      !������ ������� �� z (� �������� �����)
      read(12,302)nak		!�������� ��� ak
	read(12,301)sst		
	read(12,301)u0		!��������� �������� �������� ���������� ����
!	!!print*,'u0',u0
	read(12,301)dT		!��������� �����������
	!!print*,'dT',dT
!	read(12,301)kk		!��������� k
	kk=2*pi/lx
	!!print*,'kk',kk
!	read*
c*******************************************************************
      lx = lx0
 1889 format('lx,lx0 ')
      ly = ly0
c*******************************************************************      
c       Changing parameters depending on the variant number
  391 format(i3)
  
      open(92,file='var.dat',form='formatted')
      read(92,391) vnum
      close(92)    
      
      i1  = vnum/125
      i2  = (vnum-125*i1)/25
      i33 = (vnum-125*i1-25*i2)/5
      i4  = vnum - 125*i1 - 25*i2 - 5*i33
      
      lx = lx/2**i1
      lz = ly
      
      ni = ni0*2**i2
c       ni = 1d0 
      
      tf0 = Te0*2**i33
c       tf0 = 1d0
      
c       lp = 50
      lp = lp*(i4+1)
c      print*,vnum,i1,i2,i33,i4
c      stop
c*******************************************************************      
      ne = ni
      lz=ly               !������ ������� �� z
      ly = ly/nproc
      nee=ni0         
      !write(37,*) 'start 1.5 ','lx ',lx,'ly ',ly
      if(deb.ge.1) print*,'start 1.5 ','lx ',lx,'ly ',ly,'lz ',lz
!----------------------------------------------------------------
      
      jmf=(imp-2)*(lmp-2)*nproc*(kmp-2)*lp/fproc
      jmi=(imp-2)*(lmp-2)*nproc*(kmp-2)*lp/fproc/2
      if(deb.ge.4) write(37,*) 'lp,jmf,jmi  ',lp,jmf,jmi,imp,lmp,kmp
      if(deb.ge.4) write(37,*) 'nproc,fproc ',
     +               nproc,fproc,(lmp-2)*nproc*(kmp-2)/fproc*lp
      
      s=jm/(lx0*ly0*lz0)
      s1=ni0/s
	ksi=s1 !!!
	!!print*,'ksi',ksi
      alpha=mee/mii      !me/mi
	!!print*,'alpha= ',alpha

	!��������� ����������� ��������
	kl_m=dlog(3.d0*2.d0*pi*   
     *	dsqrt(2.d0*T**3/(nee*dsqrt(ksi**5)*(4.d0*pi*ee**2)**3)))	

!----------------------------------------------------------------
      if(outf.ne.0) then     
!   	 write(25,*)'--------------------------------'
!c	 write(25,*)'������� ���������:'
!	 write(25,532)im,lm,km
!c  532    format('������� ����� - ',i10,i10,i10)
!	 write(25,533)lp
!  533    format('Number of particles - ',i10)
!	 write(25,534)jm
!c  534    format('������ ����� ������ - ',i10)
!	 write(25,*)'--------------------------------'
!c	 write(25,*)'��������� ������:'
!	 write(25,535)lx,ly,lz
!  535    format('Domain size - ',f10.3,f10.3,f10.3)
!	 write(25,536)lx0,ly0,lz0
!c  536    format('������� ������� � �� - ',f12.6,f12.6,f12.6)
!	 write(25,537)alpha
!  537    format('alpha - ',f12.6)
!  	 write(25,538)ksi
!  538    format('ksi - ',f12.6)
!	 rD=rD1/dsqrt(ksi)
!	 !!print*,'rD',rD
!         t0=dsqrt(mee/(ppi*nee*ee**2))      !����� � ����
!         vv0=rD1/t0                          !��/���
	 vv0_m=rD/t0
         h0=dsqrt(T*ppi*nee)	
         !!print*,'h0',h0
         hx0=hx00
         !!print*,'hx0',hx0
         j0=ee*nee*dsqrt(T/mee)
         ro0=ee*nee
 	 N_D1=ne0*rD1**3
	 N_D=ne0*rD**3
	 !!print*,'N_D1',N_D1
	 !!print*,'N_D',N_D
	 !electron collision rate
	 vee=2.91d0*(1.0d-6)*ne0*kl*dsqrt(Te0**3)
	     !��������� electron collision rate
	 vee_m=2.91d0*(1.0d-6)*s*kl_m*dsqrt(Te0**3)
	 !ion collision rate
	 vii=4.8d0*(1.0d-8)*ni0*kl*dsqrt(Te0**3)
	 wce=ee*hx00/(mee*cc)
	 wci=ee*hx00/(mii*cc)
	 wpe=dsqrt(4.d0*pi*nee*ee**2/mee)
	 wpi=dsqrt(4.d0*pi*ni0*ee**2/mii)
	 fce=wce/2.d0/pi
	 fci=wci/2.d0/pi
	 fpe=wpe/2.d0/pi
  	 fpi=wpi/2.d0/pi
         !!write(37,*) 'start 33'
	 !!print*,'wpe',wpe
	 p_tau=3.d0*dsqrt(2.d0)*pi**(3.d0/2.d0)/kl/alpha*N_D1/wpe
	 !!print*,'p_tau',p_tau
	 te=1.d0/vee
  	 lambda=vv0*te
	 !!print*,'lambda',lambda
	 lambda_m=vv0_m/vee_m	!��������� ����� ���������� �������
	 !+++++++++++++++++++++++++++++
	 const1=1.71d0
	 dn=const1*kk**2*dT
	 dE=const1*kk*dT
	 gam=2.d0*3.16d0/3.d0*kk**2*tau_e
	 du=const1*gam*kk*dT
	 !+++++++++++++++++++++++++++++
	 !!print*,'du',du
	 !!print*,'const1',const1
	 !!print*,'dn',dn
	 !!print*,'gam',gam
	 !!print*,'dE',dE
         !!write(37,*) 'start 4444'
  301    format(f10.4,40x)   
  302    format(i10,40x)   
	 !����� ��������� ��������
         !write(18,401) ni0,ni
  401	 format('ni',2x,e10.3,2x,f10.4)
!        !write(18,402) nb0,nb
!  402	 format("nb",2x,e10.4,2x,e10.4)
         !write(18,403) ne0,ne
  403	 format('ne',2x,e10.3,2x,e10.3)
         !write(18,404) Te0,Te
  404	 format('Te',2x,f10.4,2x,f10.4)
         !write(18,405) hx00,hx0
  405	 format('Hx0',2x,f10.4,2x,e10.3)
!        !write(18,406) gamma,vb0
!  406	 format('gamma',2x,f10.4,2x,"vb",2x,e14.8)
         !write(18,407) lx0,lx
  407	 format('lx',2x,f10.4,2x,e14.7)
         !write(18,408) rD
  408	 format('rD',2x,f14.8)
         !write(18,409) ly
  409	 format('ly',2x,f14.8)
         !!write(37,*) 'start 5555'
!	 write(25,501) ni0
!  501    format('  ��������� �����  - ',e12.4,'  1/��^3')
!	 write(25,502) Te0
!  502    format('  ����������� ���������� - ',f6.1,' ��')
!	 write(25,*) '------------------------------'
!c         write(25,*) '������� ���������:'
!	 write(25,503) rD1
!  503    format('����� - ���������� ������ - ',f12.4,' ��')
!         ss1=1.d12*t0
!	 write(25,504) t0,ss1
!  504    format('t0 - ',e12.4,' ss1 =',e12.4,' picosec')
!	 write(25,505) vv0,beta0
!  505    format('vv0 - ',e12.4,' cm/s ',10x,'beta0=',f10.3)
!	 write(25,506) h0
!  506    format('h0 - ',e12.4,' Gs')
!         write(25,*) '-------------------------------'
!	 write(25,*) '��������� ���������: '
!	 write(25,539) rD
  539    format('Debye radius - ',f12.4,' cm')
   	 write(25,540) kl_m
  540    format('Coulomb logarythm - ',e12.4)
         write(25,*) '-------------------------------'
!------------------------------------------------	
      endif
      !write(37,*) 'jmi ',jmi
 
      tx=0.d0
      nt=0
  746 format('h1,hx,lm ',e30.20,i10,e30.20)
      h1=lx/im		!��� �� ���������� x
      if(deb.ge.1) write(37,746) h1,lm,lx
      h2=ly/lm		!��� �� ���������� y
      h3=lz/km		!��� �� ���������� z
      !!print*,'h1',h1
      !!print*,'h2',h2
      !!print*,'h3',h3
      im1=im+1
      lm1=lm+1
      km1=km+1
      im2=im+2
      lm2=lm+2
      km2=km+2
!------------------------------------------------
      !!write(37,*) 'start 666'

      i3=lp*im
      l3=lp*lm
      k3=lp*km
!      jm=i3*l3*k3
      jm=im*lm*km*lp/fproc
	!!print*,'jmp',jmp
	!!print*,'jm',jm
      if(jm.gt.jmp) then
         print*,'jmp < jm.  jm=',jm,' jmp=',jmp
	 if(outf.ne.0) then
            write(25,*) 'jmp < jm.  jm=',jm,' jmp=',jmp
	 endif    
         stop
      endif
      s=jm/(lx0*ly0*lz0)
!  507 format('  ��������� ��������� ������  - ',e12.4,'  1/��^3')
!      s1=ni0/s
!	ak0=1.d0/s1 !!!
!      beta1=dsqrt(s1)
!	!!print*,'beta1',beta1
!  508 format('s1,beta1 -',2f10.4)
!      ss=jm/(lx*ly*lz)
!
!      !���������� ������������� �������
!	s1=vv0*mee*cc/(ee*hx00)
!      s2=s1/rD	!��������� ������������� ������� � ������� �����
!c      write(25,510) s1,s2
!  510 format(' s1,s2 - ',e12.4,'      ?/cm - ',f10.4)
!	!����������� �������
!c	write(25,*)
!c      write(25,511)
!  511 format('Frequencies')
!
!	!electron gyrofrequency
!c	write(25,5121)
! 5121 format('electron gyrofrequency')
!c      write(25,512) wce,fce
!  512 format('wce - ',e12.4,' ���/�       fce - ',e12.4,'  ��')
!
!	!ion gyrofrequency
!c	write (25,5131)
! 5131 format('ion gyrofrequency')
!c      write(25,513) wci,fci
!  513 format('wci - ',e12.4,' ���/�       fci - ',e12.4,'  ��')
!
!	!electron plasma frequency
!c	write(25,5141)
! 5141 format('electron plasma frequency')
!c      write(25,514) wpe,fpe
!  514 format('wpe - ',e16.8,' ���/�       fpe - ',e16.8,'  ��')
!
!	!ion plasma frequency
!c	write(25,5151)
! 5151 format('ion plasma frequency')
!c      write(25,515) wpi,fpi
!  515 format('wpi - ',e12.4,' ���/�       fpi - ',e12.4,'  ��')
!
!	!electron collision rate
!c	write(25,5161)
! 5161 format('electron collision rate')
!c      write(25,516) vee,vii
!  516 format('ve - ',e12.4,' 1/c       vi - ',e12.4,'  1/c')
!
!
!
!	!����������� �����
!c	write(25,*)
!c      write(25,517)
!  517 format('����������� �����')
!	!������������ ������ ������
!	s3=s1*mii
!	!����������� ������������ �����
!	s4=cc/wpe
!	!������ ������������ �����
!	s5=cc/wpi

c      write(25,518) s1
  518 format(' s1 - ',e12.4,' ')

c      write(25,519) s3
  519 format('s3 - ',e12.4,'  ')

c      write(25,520) s4
  520 format('s4 - ',e12.4,'  ')

c      write(25,521) s5
  521 format('s5 - ',e12.4,' ' )

	!Velocities
	!electron thermal velocity
	vte=dsqrt(T/mee)		!T - ?
	!ion thermal velocity
	vti=dsqrt(T/mii)		!T - ?
	!Alfven velocity
	vA=hx00/dsqrt(4*pi*ni0*mii)

c      write(25,522) vte
  522 format('electron thermal velocity - ',e12.4,' cm/c ')

c      write(25,523) vti
  523 format('ion thermal velocity - ',e12.4,' cm/c ')

c      write(25,524) vA
  524 format('Alfven velocity - ',e12.4,' cm/c ')

	!����� ���������� �������


	!����� ������� ������ �������
c	write(25,*)

	t1=lx0/cc

c      if(deb.ge.4) then
c         write(37,*) 'RNG ',jmi,jmf,beamf,me*jmi*(3+12+beamf*12)
c      endif
c     MAKING THE RANDOM GENERATOR WORK THE SAME
      call PARAsetRandomIon  
c     END RANDOM GENERATOR
!---------------------------------------------------------------------
      if(deb.ge.1) print*,'start 1.7 ','lx ',lx,'ly ',ly,'lz ',lz
      !!write(37,*) 'start 777'
	!������������� �����
      ami=2.d0*ni/lp*(1d0+rbd)	!/(lp*lp*lp)
      !!print*,'ami',ami
      
      j=0
      
c      write(37,*) 'b ions ',jmi,jmf       
      
      do k=1,jmi
  	 z= lz*wrapg05cae (lz)
   	 y=meh*ly+ly* wrapg05cae (ly)
         x=lx* wrapg05cae (lx)
 2242       format('ini-c0 ',i10,3e30.20)
            if(deb.ge.4) then
               write(37,2242) k,x,y,z
            endif
            j=j+1      
            xi(j)=x
            yi(j)=y
            zi(j)=z
            ui(j)=0.d0	!g05dde(0.d0,dsqrt(1.d0*alpha)) 
            vi(j)=0.d0	!g05dde(0.d0,dsqrt(1.d0*alpha)) 
            wi(j)=0.d0	!g05dde(0.d0,dsqrt(1.d0*alpha)) 
 2241       format('ini-c ',i10,3e30.20)
            if(deb.ge.4) then
               write(37,2241) j,x,y,z
            endif
      enddo
!------------------------------------------------
      !!write(37,*) '8888888'
      
	!������������� ���������� ����
      amf=-ne/lp	!/(lp*lp*lp)
      if(deb.ge.4) write(37,*) 'ami,amf ',ami,amf
!      ak=nak*1.d0/beta1 !�����������
!      ak0=1.0d0 !�����������

	!!print*,'ak0',ak0
      j=0
	!-----------
	x=0.d0
!	dx1=0.d0
	dln=2.d0*pi*ne/kk/i3
	!!print*,'dln',dln
!      ak0=1.0d0/alpha !����������� ����������

	!------------------------------
c     MAKING THE RANDOM GENERATOR WORK THE SAME
      call PARAsetRandomBeam  
c     END RANDOM GENERATOR
	
c****************** BEAM ****************************************	
      if(beamf.eq.1) then
        j =0
        jmb = jmf/2
c         amb=0.d0
        amb = amf*rbd*2d0
          
        Tbb=(Tb*rimp)**2/(1+rimp**2)
        do i=1,jmb
           j=j+1
           xb(j)=xi(j)
           yb(j)=yi(j)
           zb(j)=zi(j)
           vb0  =wrapg05dde(0.d0,Tb*rimp)
           ux(j)=rimp+vb0
           uy(j)=wrapg05dde(0.d0,Tb*rimp) 
           uz(j)=wrapg05dde(0.d0,Tb*rimp) 
 992       format('beam-imp ',i10,5e30.20)
           write(37,992) j,xi(j),vb0,ux(j),uy(j),uz(j)
        enddo
        
        do j=1,jmb
           vb0=sqrt(1.d0-ux(j)**2-uy(j)**2-uz(j)**2)
           ub(j)=ux(j)/vb0
           vb(j)=uy(j)/vb0
           wb(j)=uz(j)/vb0
 993       format('ini-beam ',i10,7e30.20)
           if(deb.ge.4) then
              write(37,993) j+me*jmb,ux(j),uy(j),uz(j),
     +                           vb0,ub(j),vb(j),wb(j)
           endif
        enddo
      
      else
         jmb = 0
         amb = 0d0
      endif	

c     MAKING THE RANDOM GENERATOR WORK THE SAME
      call PARAsetRandomElectrons  
c     END RANDOM GENERATOR
c           open(91,file='Tex0.dat',form='formatted')       
c           write(91,*),'tex0=',tex0,'tey0=',tey0,'tez0=',tez0    
cc           print*,'tex0=',tex0,'tey0=',tey0,'tez0=',tez0    
c           close(91)

      j=0
      do i=1,jmb
           j=j+1
    	   
           xf(2*j-1)=xi(j)
           yf(2*j-1)=yi(j)
           zf(2*j-1)=zi(j)

           xf(2*j)=xi(j)
           yf(2*j)=yi(j)
           zf(2*j)=zi(j)


c          FIRST SETTING TRANVERSE
c razbros v skorostyax
           
           vy=wrapg05dde(0.d0,tey0)    
           vz=wrapg05dde(0.d0,tez0)        
 137       format('el-vel ',3e30.20)
           write(37,137) xi(j),vy,vz

c          INVERSE CURRENT
           

           termx = wrapg05dde(0.d0,tex0)
           gb0=1d0/dsqrt(1d0+ub(j)**2+vb(j)**2+wb(j)**2)
           vb0=ub(j)*gb0
 145       format('el-temp ',4e30.20)
           write(37,145) xi(j),termx,gb0,vb0
      	   if (beamf.eq.1) then
        	vf01=-rbd*vb0+termx  
                vf02=-rbd*vb0-termx
       	   else
      	        vf01=+termx
      	        vf02=-termx
      	   endif
!           format('el-vf ',)       	   
          pinv1= vf01/dsqrt((1d0-vf01**2-vy**2-vz**2)) 
          pinv2= vf02/dsqrt((1d0-vf02**2-vy**2-vz**2))
 158      format('el-vf ',5e30.20) 
          write(37,158) xi(j),vf01,vf02,pinv1,pinv2
            
          vf(2*j-1)= vy/dsqrt((1d0-vf01**2-vy**2-vz**2))
          vf(2*j)  =-vy/dsqrt((1d0-vf02**2-vy**2-vz**2))
          wf(2*j-1)= vz/dsqrt((1d0-vf01**2-vy**2-vz**2))
          wf(2*j)  =-vz/dsqrt((1d0-vf02**2-vy**2-vz**2))
          
c my correct end
         
          uf(2*j-1)= pinv1+
     +               0.01d0*dsin(mfrq*2d0*pi*xf(2*j-1)/lx)          
          uf(2*j)  = pinv2+
     +               0.01d0*dsin(mfrq*2d0*pi*xf(2*j)/lx)            
        
        
 897       format('ini-electron ',i10,7e30.20) 
           if(deb.ge.4) then
               write(37,897) j,xf(2*j),yf(2*j),zf(2*j),uf(2*j),
     +                         uf(2*j-1),vf(2*j),wf(2*j)

 344          format('ini-electron1 ',9e30.20)           
              write(37,344) xf(2*j),ub(j),uf(2*j),gb0,vb0,
     +                      vf01,vf02,pinv1,pinv2
           endif
        enddo

      xpmin = 0d0
      xpmax = xmp
      ypmin = (ly - ymp)/2
      ypmax  = ly - ypmin
      zpmin = (lz - zmp)/2
      zpmax  = lz - zpmin

      xbmin = 0d0
      xbmax = xmb
      ybmin = (ly - ymb)/2
      ybmax  = ly - ybmin
      zbmin = (lz - zmb)/2
      zbmax  = lz - zbmin
      
      if(deb.ge.2) then
         write(37,*) 'xpmin,xpmax ',xpmin,xpmax
         write(37,*) 'ypmin,ypmax ',ypmin,ypmax
         write(37,*) 'zpmin,zpmax ',zpmin,zpmax         
         write(37,*) 'xbmin,xbmax,ybmin,ybmax,zbmin,zbmax ',
     +                xbmin,xbmax,ybmin,ybmax,zbmin,zbmax         
      endif

      call writeParticleList(xf,yf,zf,uf,vf,wf,nt,jmf,'astr','el')
      call writeParticleList(xi,yi,zi,ui,vi,wi,nt,jmi,'astr','in')
      call writeParticleList(xb,yb,zb,ub,vb,wb,nt,jmb,'astr','bm')


      call cleanparticles(xi,yi,zi,ui,vi,wi,
     +     xpmin,xpmax,ypmin,ypmax,zpmin,zpmax,jmi,nt,'in','ends')
      call cleanparticles(xf,yf,zf,uf,vf,wf,
     +     xpmin,xpmax,ypmin,ypmax,zpmin,zpmax,jmf,nt,'el','ends')
      call cleanparticles(xb,yb,zb,ub,vb,wb,
     +     xbmin,xbmax,ybmin,ybmax,zbmin,zbmax,jmb,nt,'bm','ends')
!-------------------------------------------------

      do 7 k=1,kmp
      do 7 l=1,lmp
      do 7 i=1,imp
!	   x=(i-1)*h1	
         ex(i,l,k)=0.d0	!		-dE*dcos(kk*x)
         ey(i,l,k)=0.d0
         ez(i,l,k)=0.d0
!	   jx(i,l,k)=ne-dn*dsin(kk*x)	
    7 continue

      do 8 k=1,kmp
      do 8 l=1,lmp
      do 8 i=1,imp
         hx(i,l,k)=hx0
         hy(i,l,k)=0.d0
         hz(i,l,k)=0.d0
    8 continue
    !------------------------------------------------
    2 c1=tau/h1
      c2=tau/h2
      c3=tau/h3
c      write(37,855) tau,c1,c2,c3,h1,h2,h3,jmf
 855  format('tauCh,jmf ',7e10.3,i10)
      
      !print*,h1,h2,h3,tf0,ami,amf,jmi,jmf
      !!write(37,*) '999'
      if(deb.ge.1) print*,'start 1.9 ','lx ',lx,'ly ',ly,'lz ',lz
      if(deb.ge.1) print*,'ni,ni0,Te0 ',ni,ni0,Te0
      
      !!print 9007
 9007 format(2x,'kohe�  pa�ot�  start')
      ! SETTING TEST DISTRIBUTION
      call init(0)
      return
      end subroutine start

!------------------------------------------------


	!�������� ������ 

      subroutine move3(nt,xi,yi,zi,ui,vi,wi,ami,q,jm,sort)	!�������� �������� q
      implicit real*8(a-h,o-z)
      include 'part.pf'
      include 'mass.par'
c      parameter(imp=42,lmp=22,kmp=3,jmp=224000)
      real*8 ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp),
     *hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp),
     *px(imp,lmp,kmp),py(imp,lmp,kmp),pz(imp,lmp,kmp),
     *pp(imp,lmp,kmp)
      real*8 pxf(imp,lmp,kmp),pxb(imp,lmp,kmp)
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
	real*8 q,wk,y0,dbu0,dbv0,dbw0
	real*8 kappa,dwp,pi
      character*2 sort
!      integer sort
      common/b/xm,ym,zm
      common/c/im,lm,km
      common/g/hx,hy,hz
      common/d/h1,h2,h3
      common/f/tau,c1,c2,c3
      common/h/ex,ey,ez
      common/j/px,py,pz
      common/jtotal/pxf,pxb
      common/c1/alpha,t0,tD
	common/kap/kappa
 1011 format(2x,' ࠡ�⠥� ��楤��   move3')
      tau1=q*tau*0.5d0
      if(deb.ge.2) then
         write(37,*) '0 ',nt,jdeb,ui(jdeb)
      endif
      dwp = 0      
      if(magf.eq.0) then
         do i = 1,imp
            do l = 1,lmp
               do k = 1,kmp
    	          hx(i,l,k) = 0d0
	          hy(i,l,k) = 0d0
	          hz(i,l,k) = 0d0
	       enddo
            enddo    	    
         enddo
      endif
                
      y0 = ym*meh
c------------------------------------------      
      im1=im+1
      lm1=lm+1
      km1=km+1
      im2=im+2
      lm2=lm+2
      km2=km+2

      if(deb.ge.3) then
         do i = 1,imp
            do l = 1,lmp
               do k = 1,kmp
 1145             format('pxpypzmove3 ',i10,3i5,e15.5,3e30.20)
    	          write(37,1145) nt,i,l,k,ami,
     +                  px(i,l,k),py(i,l,k),pz(i,l,k)
	       enddo
            enddo    	    
         enddo
      endif

      if(deb.ge.2) then
 1306    format('sizeparticles ami,jm ',e12.3,i15)
         write(37,1306) ami,jm
      endif

      do 100 j=1,jm
      x=xi(j)
      y=yi(j)
      z=zi(j)
      pu=ui(j)
      pv=vi(j)
      pw=wi(j)
      if(deb.ge.4) then
 1521    format('crd-start a',e12.3,' x ',e30.20,' y ',e30.20,
     +          ' z ',e30.20,' pu ',e30.20,
     +          ' pv ',e30.20,' pw ',e30.20,' nt ',i5)
         write(37,1521) ami,x,y,z,pu,pv,pw,nt
      endif

      call writeParticleCurrents(px,py,pz,j,nt,'d100',
     +                                 x,y,z,pu,pv,pw)

      call writeControlAttribute(j,ami,1,x,nt)
      call writeControlAttribute(j,ami,2,y,nt)
      call writeControlAttribute(j,ami,3,z,nt)
      call writeControlAttribute(j,ami,4,pu,nt)
      call writeControlAttribute(j,ami,5,pv,nt)
      call writeControlAttribute(j,ami,6,pw,nt)


      if(deb.ge.4) then
         write(37,*) '3 ',nt,jdeb,ui(jdeb)
      endif
	!---------------------------------------------------------------------------------  
	!------------------- ������ ����� ����� ------------------------------------------
	!---------------------------------------------------------------------------------
      s2=x/h1
      call writeControlAttribute(j,ami,32,z,nt)
      call writeControlAttribute(j,ami,31,y,nt)
      call writeControlAttribute(j,ami,30,x,nt)

      call writeControlAttribute(j,ami,251,x,nt)


      i=idint(s2+1.d0)
      i1=idint(s2+1.5d0)
      s1=i-s2

       call writeControlAttribute(j,ami,14,s2,nt)
      call writeControlAttribute(j,ami,15,s2+1.0,nt)
      s_i = i
      call writeControlAttribute(j,ami,16,s_i,nt)
      call writeControlAttribute(j,ami,17,i-s2,nt)


      s2=i1-0.5d0-s2
      s_i = i
      s_i1 = i1
      call writeControlAttribute(j,ami,29,s_i1,nt)
      call writeControlAttribute(j,ami,33,s_i,nt)
      call writeControlAttribute(j,ami,34,s1,nt)
      call writeControlAttribute(j,ami,35,s2,nt)
           
       s4=(y-y0)/h2
       l=idint(s4+1.d0)
       l1=idint(s4+1.5d0)
       s_l1 = l1
       s3=l-s4
       s4=l1-0.5d0-s4
       s6=z/h3
       k=idint(s6+1.d0)
       k1=idint(s6+1.5d0)
       s5=k-s6
       s6=k1-0.5d0-s6
      
      s11=1.d0-s1
      s21=1.d0-s2
      s31=1.d0-s3
      s41=1.d0-s4
      s51=1.d0-s5
      s61=1.d0-s6
      if((deb.ge.4).and.(j.eq.jdeb).and.(ami.lt.0d0)) then
  376    format('s1,s2,s3,s4,s5,s6 ',6e15.5)
         write(37,376) s1,s2,s3,s4,s5,s6
      endif

      
      if((deb.ge.4)) then
 165     format('sx-field y,ami ',e30.20,e12.3,' tau1 ',e30.20,' s1 ',
     +           e30.20,' s4 ',e30.20,' s6 ',e30.20, ' s61 ',e30.20,
     +           ' s41 ',e30.20,' s11 ',e30.20)
 166     format('electric y,ami ',e30.20,e12.3,
     +          ' ex(i,l1,k1) ',e30.20,' ex(i,l1,k1+1) ',e30.20,
     +          ' ex(i,l1+1,k1) ',e30.20,' ex(i,l1+1,k1+1) ',e30.20,
     +          ' ex(i+1,l1,k1) ',e30.20,' ex(i+1,l1,k1+1) ',e30.20,
     +          ' ex(i+1,l1+1,k1) ',e30.20,' ex(i+1,l1+1,k1+1) ',e30.20
     +          ' l ',i3,' l1 ',i3,' nt ',i5)

         !write(37,*) 'nt,j,zi(j),ui(j) ',nt,prev_jm+j,zi(prev_jm+j),ui(prev_jm+j)
         write(37,165) y,ami,tau1,s1,s4,s6,s61,s41,s11
         write(37,166) y,ami,
     +          ex(i,l1,k1), ex(i,l1,k1+1),
     +          ex(i,l1+1,k1),ex(i,l1+1,k1+1),
     +          ex(i+1,l1,k1),ex(i+1,l1,k1+1),
     +          ex(i+1,l1+1,k1),ex(i+1,l1+1,k1+1),
     +          l+(lmp-2)*meh,l1+(lmp-2)*meh,nt

      endif
      
      sx=tau1*(s1*(s4*(s6*ex(i,l1,k1)+s61*ex(i,l1,k1+1))+
     +s41*(s6*ex(i,l1+1,k1)+s61*ex(i,l1+1,k1+1)))+
     +s11*(s4*(s6*ex(i+1,l1,k1)+s61*ex(i+1,l1,k1+1))+
     +s41*(s6*ex(i+1,l1+1,k1)+s61*ex(i+1,l1+1,k1+1))))

      call writeControlAttribute(j,ami,110,(s6*ex(i,l1,k1)+
     +                          s61*ex(i,l1,k1+1)),nt)
      call writeControlAttribute(j,ami,111,(s6*ex(i,l1+1,k1)+
     +                          s61*ex(i,l1+1,k1+1)),nt)
      call writeControlAttribute(j,ami,112,(s6*ex(i+1,l1,k1)+
     + s61*ex(i+1,l1,k1+1)),nt)
      call writeControlAttribute(j,ami,113,(s6*ex(i+1,l1+1,k1)+
     +s61*ex(i+1,l1+1,k1+1)),nt)

      call writeControlAttribute(j,ami,114,s1,nt)
      call writeControlAttribute(j,ami,115,s11,nt)
      call writeControlAttribute(j,ami,116,s4,nt)
      call writeControlAttribute(j,ami,117,s41,nt)

      call writeControlAttribute(j,ami,128,
     +     s1*(s4*(s6*ex(i,l1,k1)+s61*ex(i,l1,k1+1))+
     +s41*(s6*ex(i,l1+1,k1)+s61*ex(i,l1+1,k1+1))),nt)
 
      call writeControlAttribute(j,ami,129,
     +s11*(s4*(s6*ex(i+1,l1,k1)+s61*ex(i+1,l1,k1+1))+
     +s41*(s6*ex(i+1,l1+1,k1)+s61*ex(i+1,l1+1,k1+1))),nt)

      call writeControlAttribute(j,ami,120,ex(i,l1,k1),nt)
      call writeControlAttribute(j,ami,121,ex(i,l1+1,k1),nt)
      call writeControlAttribute(j,ami,122,ex(i+1,l1,k1),nt)
      call writeControlAttribute(j,ami,123,ex(i+1,l1+1,k1),nt)
      call writeControlAttribute(j,ami,124,ex(i,l1,k1+1),nt)
      call writeControlAttribute(j,ami,125,ex(i,l1+1,k1+1),nt)
      call writeControlAttribute(j,ami,126,ex(i+1,l1,k1+1),nt)
      call writeControlAttribute(j,ami,127,ex(i+1,l1+1,k1+1),nt)

      call writeControlAttribute(j,ami,252,x,nt)

      s_i  = i
      s_l1 = l1
      s_k1 = k1

      call writeControlAttribute(j,ami,118,s_i,nt)
      call writeControlAttribute(j,ami,119,s_l1,nt)
      call writeControlAttribute(j,ami,128,s_k1,nt)

      call writeControlAttribute(j,ami,18,s1,nt)
      call writeControlAttribute(j,ami,19,s4,nt)
      call writeControlAttribute(j,ami,20,s6,nt)
      call writeControlAttribute(j,ami,21,s61,nt)
      call writeControlAttribute(j,ami,22,s41,nt)


      sy=tau1*(s2*(s3*(s6*ey(i1,l,k1)+s61*ey(i1,l,k1+1))+
     +s31*(s6*ey(i1,l+1,k1)+s61*ey(i1,l+1,k1+1)))+
     +s21*(s3*(s6*ey(i1+1,l,k1)+s61*ey(i1+1,l,k1+1))+
     +s31*(s6*ey(i1+1,l+1,k1)+s61*ey(i1+1,l+1,k1+1))))

      call writeControlAttribute(j,ami,23,s2,nt)
      call writeControlAttribute(j,ami,24,s3,nt)
      call writeControlAttribute(j,ami,25,s21,nt)
      call writeControlAttribute(j,ami,26,s31,nt)


      sz=tau1*(s2*(s4*(s5*ez(i1,l1,k)+s51*ez(i1,l1,k+1))+
     +s41*(s5*ez(i1,l1+1,k)+s51*ez(i1,l1+1,k+1)))+
     +s21*(s4*(s5*ez(i1+1,l1,k)+s51*ez(i1+1,l1,k+1))+
     +s41*(s5*ez(i1+1,l1+1,k)+s51*ez(i1+1,l1+1,k+1))))

      call writeControlAttribute(j,ami,27,s5,nt)
      call writeControlAttribute(j,ami,28,s51,nt)


      call writeControlAttribute(j,ami,33,sx/tau1,nt)
      call writeControlAttribute(j,ami,34,sy/tau1,nt)
      call writeControlAttribute(j,ami,35,sz/tau1,nt)


      s_i  = i
      s_l1 = l1
      s_k1 = k1
      call writeControlAttribute(j,ami,7,s_i,nt)
      call writeControlAttribute(j,ami,8,s_l1,nt)
      call writeControlAttribute(j,ami,9,s_k1,nt)

      call writeControlAttribute(j,ami,10,ex(i,l1,k1),nt)
      call writeControlAttribute(j,ami,11,h1,nt)
      call writeControlAttribute(j,ami,12,h2,nt)

      call writeControlAttribute(j,ami,13,h3,nt)



      if((j.eq.jdeb).and.(deb.ge.2).and.(ami.lt.0d0)) then

      endif    
      call writeControlAttribute(j,ami,253,x,nt)

      if(deb.ge.4) then
 1695    format('pu-sx x,a',e30.20,e12.3,' pu ',e30.20,' pv ',e30.20,
     +          ' pw ',e30.20,' sx ',e30.20,
     +          ' sy ',e30.20,' sz ',e30.20,' nt ',i5)
         write(37,1695) x,ami,pu,pv,pw,sx,sy,sz,nt         
      endif
      
      call traceparticle(y0p,sx,sy,sz)


      pu = pu + sx
      pv = pv + sy
      pw = pw + sz


	!----------------------------------------------------------------------------------
	!------------------������ ����� ����� ---------------------------------------------
	!----------------------------------------------------------------------------------
      ps=tau1*beta0/dsqrt(1.d0+beta0**2*(pu*pu+pv*pv+pw*pw))

      call writeControlAttribute(j,ami,48,tau1,nt)

      if((j.eq.jdeb).and.(deb.ge.4)) then
	write(37,*) 'ps',ps,s2,s3,s5,s51,s31
      endif	
      bx=ps*(s2*(s3*(s5*hx(i1,l,k)+s51*hx(i1,l,k+1))+
     +s31*(s5*hx(i1,l+1,k)+s51*hx(i1,l+1,k+1)))+
     +s21*(s3*(s5*hx(i1+1,l,k)+s51*hx(i1+1,l,k+1))+
     +s31*(s5*hx(i1+1,l+1,k)+s51*hx(i1+1,l+1,k+1))))
      by=ps*(s1*(s4*(s5*hy(i,l1,k)+s51*hy(i,l1,k+1))+
     +s41*(s5*hy(i,l1+1,k)+s51*hy(i,l1+1,k+1)))+
     +s11*(s4*(s5*hy(i+1,l1,k)+s51*hy(i+1,l1,k+1))+
     +s41*(s5*hy(i+1,l1+1,k)+s51*hy(i+1,l1+1,k+1))))
      bz=ps*(s1*(s3*(s6*hz(i,l,k1)+s61*hz(i,l,k1+1))+
     +s31*(s6*hz(i,l+1,k1)+s61*hz(i,l+1,k1+1)))+
     +s11*(s3*(s6*hz(i+1,l,k1)+s61*hz(i+1,l,k1+1))+
     +s31*(s6*hz(i+1,l+1,k1)+s61*hz(i+1,l+1,k1+1))))

      call writeControlAttribute(j,ami,36,bx/ps,nt)
      call writeControlAttribute(j,ami,37,by/ps,nt)
      call writeControlAttribute(j,ami,38,bz/ps,nt)


      su=pu+pv*bz-pw*by
      sv=pv+pw*bx-pu*bz
      sw=pw+pu*by-pv*bx



      s1=bx*bx
      s2=by*by
      s3=bz*bz
      s4=bx*by
      s5=by*bz
      s6=bz*bx
      s=1.d0+s1+s2+s3



      if((deb.ge.4)) then
         write(37,*) '2 nt,j,zi(j),ui(j) ',nt,j,zi(j),ui(j)
      endif

	!----------------------------------------------------------------------------------
	!--------- ������ ����� ����� -----------------------------------------------------
	!----------------------------------------------------------------------------------
	!��� �� ����������� �������� � �������� ��� ���������� �������
      pu1=((1.d0+s1)*su+(s4+bz)*sv+(s6-by)*sw)/s	
      pv1=((s4-bz)*su+(1.d0+s2)*sv+(s5+bx)*sw)/s
      pw1=((s6+by)*su+(s5-bx)*sv+(1.d0+s3)*sw)/s

      if((j.eq.jdeb).and.(deb.ge.2).and.(ami.lt.0d0)) then
         write(37,*) 'pw1,sz,s ',nt,j,pw1,sz,s
	 write(37,*) 'bx,by,bz ',bx,by,bz
	 write(37,*) 'sx,sy,sz ',sx,sy,sz
      endif

      pu=pu1 + sx	
      pv=pv1 + sy
      pw=pw1 + sz

	  ps=(pu*pu+pv*pv+pw*pw)
      ps=1.d0/dsqrt(1.d0+beta0**2*(pu*pu+pv*pv+pw*pw))
      if((j.eq.jdeb).and.(deb.ge.4)) then
         write(37,*) 'beta0 ',beta0
         write(37,*) me,nt,j,zi(j),wi(j),z
      endif
      if((deb.ge.4)) then
         write(37,*) '3 nt,j,zi(j),ui(j) ',nt,j,zi(j),ui(j)
         write(37,*) '3 ps,pu,pv,pw ',ps,pu,pv,pw
      endif
      call writeControlAttribute(j,ami,254,x,nt)
      
      u=ps*pu
      v=ps*pv
      w=ps*pw
      x1=x+tau*u
      y1=y+tau*v
      z1=z+tau*w

      if(deb.ge.4) then
 1410    format('x1y1z1 j,x1,x,u,ps,pu,y1,y,v,pv,z1,z,w,pw,nt ',
     +          e12.3,' ',i5,' ',13e30.20,' nt ',i5)
         write(37,1410) ami,j,x1,x,u,ps,pu,y1,
     +   y,v,pv,z1,z,w,pw,nt
      endif

      call writeControlAttribute(j,ami,255,x1,nt)

      call writeControlAttribute(j,ami,39,x1,nt)
      call writeControlAttribute(j,ami,40,y1,nt)
      call writeControlAttribute(j,ami,41,z1,nt)
      call writeControlAttribute(j,ami,42,pu,nt)
      call writeControlAttribute(j,ami,43,pv,nt)
      call writeControlAttribute(j,ami,44,pw,nt)


c*************************************************      
 160  format(2i10,4e10.3)
  
      if((ami.gt.-1e-3).and.(ami.lt.0).and.(deb.ge.4)) then
         write(37,160) nt,j,x1,ui(j),u,sx
      endif
!---------------------------------------------------------------------------------
      s2=x1/h1
      i2=idint(s2+1.5d0)
      s4=(y1-y0)/h2
      l2=idint(s4+1.5d0)
      if((deb.ge.4)) then
         write(37,*) '4 nt,j,zi(j),wi(j) ',nt,j,zi(j),ui(j)
      endif

      s6=z1/h3
      k2=idint(s6+1.5d0)
      i=abs(i2-i1)
      l=abs(l2-l1)
      k=abs(k2-k1)
      s_i = i
      s_l = l
      s_k = k
      call writeControlAttribute(j,ami,45,s_i,nt)
      call writeControlAttribute(j,ami,46,s_l,nt)
      call writeControlAttribute(j,ami,47,s_k,nt)
      m=4*i+2*l+k
      if((deb.ge.4)) then
 1459    format('j,x,x1,y,y1,y0,z,z1,i,l,k,i1,l1,k1,i2,l2,k2,m ',
     +           i15,e12.3,7e30.20,10i4,' nt ',i5)
         write(37,1459) j,ami,x,x1,y,y1,y0,z,z1,
     +                  i,l,k,i1,l1,k1,i2,l2,k2,m,nt
      endif


      if((nt.eq.5).and.(deb.ge.1).and.(ami.gt.0d0)) then
         write(37,*) 'x,x1 ',x,x1
         write(37,*) 'y,y1,y0 ',y,y1,y0
         write(37,*) 'z,z1 ',z,z1
         write(37,*) 'nt,j,i,l,k,i1,l1,k1,i2,l2,k2,m ',nt,j,i,l
	 write(37,*)  k,i1,l1,k1,i2,l2,k2,m
      endif
      goto(1,2,3,4,5,6,7),m
      call pqr(i1,l1,k1,x,y,z,y0,x1,y1,z1,ami,j,jm,0,nt)
      goto 18
    1 z2=h3*(0.5d0*(k1+k2)-1.d0)
      s=(z2-z)/(z1-z)
      x2=x+(x1-x)*s  
      y2=y+(y1-y)*s  
      goto 11
    2 y2=h2*(0.5d0*(l1+l2)-1.d0)+y0
      s=(y2-y)/(y1-y)
      x2=x+(x1-x)*s  
      z2=z+(z1-z)*s  
      goto 11
    3 y2=h2*(0.5d0*(l1+l2)-1.d0)+y0
      z2=h3*(0.5d0*(k1+k2)-1.d0)
      s=((z1-z)*(z2-z)+(y1-y)*(y2-y))/((z1-z)**2+(y1-y)**2)
      x2=x+(x1-x)*s  
      goto 11
    4 x2=h1*(0.5d0*(i1+i2)-1.d0)
      s=(x2-x)/(x1-x)
      y2=y+(y1-y)*s  
      z2=z+(z1-z)*s  
      goto 11
    5 x2=h1*(0.5d0*(i1+i2)-1.d0)
      z2=h3*(0.5d0*(k1+k2)-1.d0)
      s=((z1-z)*(z2-z)+(x1-x)*(x2-x))/((z1-z)**2+(x1-x)**2)
      y2=y+(y1-y)*s  
      goto 11
    6 x2=h1*(0.5d0*(i1+i2)-1.d0)
      y2=h2*(0.5d0*(l1+l2)-1.d0)+y0
      s=((y1-y)*(y2-y)+(x1-x)*(x2-x))/((y1-y)**2+(x1-x)**2)
      z2=z+(z1-z)*s  
      goto 11
    7 x2=h1*(0.5d0*(i1+i2)-1.d0)
      y2=h2*(0.5d0*(l1+l2)-1.d0)+y0
      z2=h3*(0.5d0*(k1+k2)-1.d0)
  11  if((j.eq.jdeb).and.(deb.ge.4)) then
         write(37,*) me,nt,j,z,z2,z1,zi(j),wi(j)
      endif
      call pqr(i1,l1,k1,x,y,z,y0,x2,y2,z2,ami,j,jm,0,nt)
      call pqr(i2,l2,k2,x2,y2,z2,y0,x1,y1,z1,ami,j,jm,1,nt)
   18 continue   
      call writeControlAttribute(j,ami,256,x1,nt)

      s1=dsign(1.d0,x1*(xm-x1))
      s1=s1+dsign(1.d0,y1*(ym-y1))
      s1=s1+dsign(1.d0,z1*(zm-z1))
      if(s1.gt.2.5d0) goto 15


  112 if(x1.gt.0.d0) goto 12
      !x1=x1+xm
      call writeControlAttribute(j,ami,257,x1,nt)
      goto 13
   12 if(x1.le.xm) goto 13
      !x1=x1-xm
   13 if(z1.gt.0.d0) goto 14
      z1=z1+zm
      goto 15
   14 if(z1.le.zm) goto 15
      z1=z1-zm
   15 continue

      if((z1.lt.0).or.(z1.gt.zm)) goto 112

      call beamboundcheck(xi(j),x1,sort,j,nt)


      xi(j)=x1
      yi(j)=y1
      zi(j)=z1
      ui(j)=pu
      vi(j)=pv
      wi(j)=pw
      call writeControlAttribute(j,ami,258,xi(j),nt)

      if((deb.ge.4)) then
         write(37,*) '5 nt,j,zi(j),ui(j) ',nt,j,zi(j),ui(j)
      endif
        
         
  100 continue

      if(deb.ge.4) then
         do i = 1,imp
            do l = 1,lmp
               do k = 1,kmp
 1553             format('p-in-result ',i3,3i5,e30.20,e12.3,' nt ',i5)
                  write(37,1553) 0,
     +             i,l+(lmp-2)*meh,k,px(i,l,k),ami,nt
	       enddo
            enddo    	    
         enddo
      endif
        
       
      

      if(deb.ge.2)  write(37,*) 'nt,ym,jm1 ',nt,ym,jm,sort,ui(1)

      if(deb.ge.4) then
         do i = 1,imp
            do l = 1,lmp
               do k = 1,kmp
 1554             format('p-result ',i3,3i5,e30.20,e12.3,' nt ',i5)
                  write(37,1554) 0,
     +              i,l+(lmp-2)*meh,k,px(i,l,k),ami,nt
	       enddo
            enddo    	    
         enddo
      endif


      if((deb.ge.4)) then
         write(37,*) '6 nt,j,zi(j),ui(j) ',nt,j,zi(j),ui(jdeb)
      endif
     
 
      if(deb.ge.2)  write(37,*) 'nt,ym,jm ',nt,ym,jm,sort
      !stop

!      if(sort.eq.2) then
!         chsort='bm'
!      endif
!      if(sort.eq.1) then
!         chsort='in'
!      endif
!      if(sort.eq.0) then
!         chsort='el'
!      endif
!      print *,chsort

!      call PARAfinal
!      stop

      
      call writeParticleList(xi,yi,zi,ui,vi,wi,nt,jm,'btra',sort)

!      call PARAfinal
!      stop

      call PARAtransferPart(nt,ym,jm,xi,yi,zi,ui,vi,wi) 
 
!      if(sort.eq.2) then
!         chsort='bm'
!      endif
!      if(sort.eq.1) then
!         chsort='in'
!      endif
!      if(sort.eq.0) then
!         chsort='el'
!      endif


      call writeParticleList(xi,yi,zi,ui,vi,wi,nt,jm,'atra',sort)

!      call PARAfinal
!      stop
 
c      stop
      
      if((deb.ge.4)) then
         write(37,*) 'AT nt,j,zi(j),ui(j) ',nt,j,zi(j),ui(jdeb)
      endif


      if((ami.lt.0).and.(deb.ge.2)) then
         write(37,*) 'dwp ',dwp,ami
      endif
        
      if(deb.ge.4) then
        write(37,*) 'b CURRENTS =================================== ',nt
         do i = 1,imp
            do l = 1,lmp
               do k = 1,kmp
  453             format('curr-end-move-X ',3i5,e12.3,e30.20,' nt ',i5,
     +                                   ' meh ',i4 )
                  write(37,453) i,l+(lmp-2)*meh,k,ami,px(i,l,k),nt,meh
               enddo
            enddo
         enddo 
         write(37,*) 'px(2,2,2) ',px(2,2,2)
         write(37,*) 'END CURRENTS ================================ ',nt
       
      endif

      if(deb.ge.4) then
        write(37,*) 'b CURRENTS =================================== ',nt
         do i = 1,imp
            do l = 1,lmp
               do k = 1,kmp
 4531             format('curr-end-move-Y ',3i5,e12.3,e30.20,' nt ',i5,
     +                                   ' meh ',i4 )
                  write(37,4531) i,l+(lmp-2)*meh,k,ami,py(i,l,k),nt,meh
               enddo
            enddo
         enddo
         write(37,*) 'px(2,2,2) ',px(2,2,2)
         write(37,*) 'END CURRENTS ================================ ',nt

      endif

      if(deb.ge.4) then
        write(37,*) 'b CURRENTS =================================== ',nt
         do i = 1,imp
            do l = 1,lmp
               do k = 1,kmp
 4532             format('curr-end-move-Z ',3i5,e12.3,e30.20,' nt ',i5,
     +                                   ' meh ',i4 )
                  write(37,4532) i,l+(lmp-2)*meh,k,ami,pz(i,l,k),nt,meh
               enddo
            enddo
         enddo
         write(37,*) 'px(2,2,2) ',px(2,2,2)
         write(37,*) 'END CURRENTS ================================ ',nt

      endif



        
       
c  kohe� c�mm�pobah��.
       
      if(deb.ge.4) write(37,*) 'end mo ',xi(jdeb)
      return
      end subroutine move3

c------------------------------------------------------------
      subroutine pqr(i,l,k,x,y,z,y0,x1,y1,z1,a1,j,jmf,num,nt)
      implicit none
      include 'part.pf'
c      parameter(imp=42,lmp=22,kmp=3)
      real*8 p(imp,lmp,kmp),q(imp,lmp,kmp),r(imp,lmp,kmp),y0
      real*8 dx1,dy1,dz1,s1,s2,s3,su,sv,sw,h1,h2,h3,c1,c2,c3
      real*8 dx,dy,dz,a,tau,x,x1,y,y1,z,z1,a1,t_ilk
      integer j,i,l,k,jmf,num,nt
      real*8  s_i,s_l,sk
      common/d/h1,h2,h3
      common/f/tau,c1,c2,c3
      common/j/p,q,r
      dx=0.5d0*(x+x1)-h1*(i-1.5d0)
      dy=0.5d0*(y+y1)-y0-h2*(l-1.5d0)
      dz=0.5d0*(z+z1)-h3*(k-1.5d0)

      if(deb.ge.4) then
 1641    format('dxyz j,dx ',i10,e30.20,' x ',e30.20,' x1 ',e30.20,
     +          ' h1 ',e30.20,' i ',i5)
 1642    format('dxyz j,dy,a ',i10,e30.20,e12.3,' y ',e30.20,' y1 '
     +          ,e30.20,' h2 ',e30.20,' l ',i5,' y0 ',e30.20)
 1643    format('dxyz j,dz ',i10,e30.20,' z ',e30.20,' z1 ',e30.20,
     +          ' h3 ',e30.20,' k ',i5)
         write(37,1641) j,dx,x,x1,h1,i
         write(37,1642) j,dy,a1,y,y1,h2,l,y0
         write(37,1643) j,dz,z,z1,h3,k
      endif

      call writeControlAttribute(j,a1,48+num,dx,nt)
      call writeControlAttribute(j,a1,50+num,dy,nt)
      call writeControlAttribute(j,a1,52+num,dz,nt)

c      dx=0.5d0*(x+x1)-h1*(i-1.0d0)
c      dy=0.5d0*(y+y1)-y0-h2*(l-1.0d0)
c      dz=0.5d0*(z+z1)-h3*(k-1.0d0)
      
      a = a1
      
      dx1=h1-dx
      dy1=h2-dy
      dz1=h3-dz
      call writeControlAttribute(j,a1,90+num,dx1,nt)
      call writeControlAttribute(j,a1,92+num,dy1,nt)
      call writeControlAttribute(j,a1,94+num,dz1,nt)
      su=x1-x
      sv=y1-y
      sw=z1-z
      if(deb.ge.4) then
 1629   format('uvw j,a1,su,sv,sw,x,x1,y,y1,z,z1 ',i5,e12.3,9e30.20,i2)
         write(37,1629) j,a1,su,sv,sw,x,x1,y,y1,z,z1,num
      endif
      s1=sv*sw/12.d0
      s2=su*sw/12.d0
      s3=su*sv/12.d0

      call writeControlAttribute(j,a1,54+num,su,nt)
      call writeControlAttribute(j,a1,56+num,sv,nt)
      call writeControlAttribute(j,a1,58+num,sw,nt)


      if(deb.ge.4) then
         write(37,*) 'j,i,l,k ',j,i,l,k
         write(37,*) 'y1,y,z1,z ',y1,y,z1,z
         write(37,*) 'sv,sw ',sv,sw
         write(37,369) j,x1,x,su,a,tau,h1,h2
         write(37,3701) j,s1,s2,s3, dx1,dy1,dz1, su,x1,h2,h3
!         write(37,371) z,a1,p(i,l,k),p(i,l,k+1),p(i,l+1,k),p(i,l+1,k+1)
      endif

      call writeControlAttribute(j,a1,106,su*a,nt)
      call writeControlAttribute(j,a1,107,tau*h2*h3,nt)
      call writeControlAttribute(j,a1,108,su*a/(tau*h2*h3),nt)
      call writeControlAttribute(j,a1,109,1.0d0/(tau*h2*h3),nt)
      call writeControlAttribute(j,a1,110,su*a/(tau*h2*h3),nt)


      su=su*a/(tau*h2*h3)
      sv=sv*a/(tau*h1*h3)
      sw=sw*a/(tau*h1*h2)

      call writeControlAttribute(j,a1,96+num,su,nt)
      call writeControlAttribute(j,a1,98+num,sv,nt)
      call writeControlAttribute(j,a1,100+num,sw,nt)
      call writeControlAttribute(j,a1,102,a,nt)
      call writeControlAttribute(j,a1,103,tau,nt)
      call writeControlAttribute(j,a1,104,h2,nt)
      call writeControlAttribute(j,a1,105,h3,nt)

      if(deb.ge.4) then
  369    format('pqqr th,j,x1,x,su,a,tau,h1,h2 ',i10,7e12.3)     
 1369    format('PQR th,j,x1,x,su,a,tau,h1,h2 ',i10,7e12.3)     
         write(37,1369) j,x1,x,su,a,tau,h1,h2
         write(37,3702) j,s1,s2,s3, dx1,dy1,dz1, su,x1,x
!         write(37,371) j,a1,p(i,l,k),p(i,l,k+1),
!     +                 p(i,l+1,k),p(i,l+1,k+1)
      endif
      s_i = i
      s_l = l
      sk = k
      call writeControlAttribute(j,a1,60+num,s_i,nt)
      call writeControlAttribute(j,a1,62+num,s_l,nt)
      call writeControlAttribute(j,a1,64+num,sk,nt)

      call writeControlAttribute(j,a1,250+num,p(i,l,k),nt)
      call writeControlAttribute(j,a1,252+num,p(i,l,k+1),nt)
      call writeControlAttribute(j,a1,254+num,p(i,l+1,k),nt)
      call writeControlAttribute(j,a1,256+num,p(i,l+1,k+1),nt)
      
      if(deb.ge.4) then
 1695    format('p-triple ',e30.20,e12.3,
     +          ' (',3i4,') ',
     +          ' (',3i4,') ',
     +          ' (',3i4,') ',
     +          ' (',3i4,') ')
         write(37,1695) x,a1,i,l,k,i,l,k+1,i,l+1,k,i,l+1,k+1

 1696    format('ppp j,su ',e15.5,e30.20,' sv ',e30.20,
     +          ' sw ',e30.20,' dy1 ',e30.20,' dz1 ',e30.20,' s1 ',
     +          e30.20,
     +          ' dz ',e30.20,' dy ',e30.20,' dx ',e30.20,' dx1 ',
     +          e30.20,
     +          ' s2 ',e30.20,' s3 ',e30.20,
     +          ' x  ',e30.20,' x1 ',e30.20,
     +          ' y  ',e30.20,' y1 ',e30.20,
     +          ' z  ',e30.20,' z1 ',e30.20,
     +          ' p  ',e30.20,' ilk ',3i5,' nt ',i5
     +)
         write(37,1696) a1,su,sv,sw,dy1,dz1,s1,dz,dy,dx,dx1,s2,s3,
     +                  x,x1,y,y1,z,z1,q(i,l,k),i,l,k,nt
  361    format('befor-p y ',e30.20,' j ',i10,' a ',e12.3,
     +         ' p ',e30.20,
     +         ' pk ',e30.20,
     +         ' pl ',e30.20,
     +         ' plk ',e30.20,
     +         ' ilk ',3i5,
     +         ' nt ',i5
     +         )
         write(37,361) y,j,a1,p(i,l,k),p(i,l,k+1),p(i,l+1,k),
     +         p(i,l+1,k+1), 
     +         i,l+(lmp-2)*meh,k,nt

      endif
      t_ilk = p(i,l,k) 

!$omp critical
      p(i,l,k)=p(i,l,k)+su*(dy1*dz1+s1)
      p(i,l,k+1)=p(i,l,k+1)+su*(dy1*dz-s1)
      p(i,l+1,k)=p(i,l+1,k)+su*(dy*dz1-s1)
      p(i,l+1,k+1)=p(i,l+1,k+1)+su*(dy*dz+s1)
!$omp end critical

      call writeControlAttribute(j,a1,260+num,p(i,l,k),nt)
      call writeControlAttribute(j,a1,262+num,p(i,l,k+1),nt)
      call writeControlAttribute(j,a1,264+num,p(i,l+1,k),nt)
      call writeControlAttribute(j,a1,266+num,p(i,l+1,k+1),nt)



      call writeControlAttribute(j,a1,66+num,su*(dy1*dz1+s1),nt)
      call writeControlAttribute(j,a1,68+num,su*(dy1*dz-s1),nt)
      call writeControlAttribute(j,a1,70+num,su*(dy*dz1-s1),nt)
      call writeControlAttribute(j,a1,72+num,su*(dy*dz+s1),nt)


      if(deb.ge.4) then
 3701   format('p1qr th,j,s1,s2,s3,dx1,dy1,dz1,su,sv,sw,x1,x,s3,h2,h3 ',
     +          i10,11e30.20)
 370     format('pqr ',3e30.20,i10,12e30.20)
!                        3  4  5  6   7   8   9 10 11 12 13 14 15
 3702    format('pq1r th,j,s1,s2,s3,dx1,dy1,dz1,su,sv,sw,x1,x,s3 ',
     +          i10,9e30.20)
         write(37,370) a1,z,z1,j,s1,s2,s3,dx1,
     +                 dy1,dz1,su,sv,sw,x1,x,s3

 3611    format('after-pX y ',e30.20,' j ',i10,' a ',e12.3,
     +         ' p ',e30.20,
     +         ' pk ',e30.20,
     +         ' pl ',e30.20,
     +         ' plk ',e30.20,
     +         ' ilk ',3i5,
     +         ' nt ',i5
     +         )
         write(37,3611) y,j,a1,p(i,l,k),p(i,l,k+1),p(i,l+1,k),
     +         p(i,l+1,k+1),
     +         i,l+(lmp-2)*meh,k,nt

  372    format('q q(k+1) q(l+1) q(l+1,k+1)  ',i10,4e12.3)
  373    format('r r(k+1) r(l+1) r(l+1,k+1)  ',i10,4e12.3)
         write(37,372) j,q(i,l,k),q(i,l,k+1),q(i,l+1,k),
     +q(i,l+1,k+1)
      endif

      
      q(i,l,k)=q(i,l,k)+sv*(dx1*dz1+s2)
      q(i,l,k+1)=q(i,l,k+1)+sv*(dx1*dz-s2)
      q(i+1,l,k)=q(i+1,l,k)+sv*(dx*dz1-s2)
      q(i+1,l,k+1)=q(i+1,l,k+1)+sv*(dx*dz+s2)
!$omp end critical

      if(deb.ge.4) then
3612    format('after-pY y ',e30.20,' j ',i10,' a ',e12.3,
     +         ' p ',e30.20,
     +         ' pk ',e30.20,
     +         ' pl ',e30.20,
     +         ' plk ',e30.20,
     +         ' ilk ',3i5,
     +         ' nt ',i5
     +         )
         write(37,3612) y,j,a1,q(i,l,k),q(i,l,k+1),q(i,l+1,k),
     +         q(i,l+1,k+1),
     +         i,l+(lmp-2)*meh,k,nt

      endif

      call writeControlAttribute(j,a1,74+num,sv*(dx1*dz1+s2),nt)
      call writeControlAttribute(j,a1,76+num,sv*(dx1*dz-s2),nt)
      call writeControlAttribute(j,a1,78+num,sv*(dx*dz1-s2),nt)
      call writeControlAttribute(j,a1,80+num,sv*(dx*dz+s2),nt)

      
      r(i,l,k)=r(i,l,k)+sw*(dx1*dy1+s3)
      r(i,l+1,k)=r(i,l+1,k)+sw*(dx1*dy-s3)
      r(i+1,l,k)=r(i+1,l,k)+sw*(dx*dy1-s3)
      r(i+1,l+1,k)=r(i+1,l+1,k)+sw*(dx*dy+s3)

      call writeControlAttribute(j,a1,82+num,sw*(dx1*dy1+s3),nt)
      call writeControlAttribute(j,a1,84+num,sw*(dx1*dy-s3),nt)
      call writeControlAttribute(j,a1,86+num,sw*(dx*dy1-s3),nt)
      call writeControlAttribute(j,a1,88+num,sw*(dx*dy+s3),nt)

      if(deb.ge.4) then
 3613    format('after-pZ y ',e30.20,' j ',i10,' a ',e12.3,
     +         ' p ',e30.20,
     +         ' pk ',e30.20,
     +         ' pl ',e30.20,
     +         ' plk ',e30.20,
     +         ' ilk ',3i5,
     +         ' nt ',i5
     +         )
         write(37,3613) y,j,a1,r(i,l,k),r(i,l,k+1),r(i,l+1,k),
     +         r(i,l+1,k+1),
     +         i,l+(lmp-2)*meh,k,nt

      endif

      if((deb.ge.4).and.(j.eq.jdeb).and.(a.lt.0d0)) then
         write(37,*) 'AFTER==========================='
         
!        write(37,370) a1,z,z1,j,s1,s2,s3,dx1,dy1,dz1,su,sv,sw,x1,x,s3
!         write(37,371) j,a1,p(i,l,k),p(i,l,k+1),p(i,l+1,k),
!     +                 p(i,l+1,k+1)
         write(37,372) j,q(i,l,k),q(i,l,k+1),q(i,l+1,k),
     +                 q(i,l+1,k+1)
         write(37,373) j,r(i,l,k),r(i,l,k+1),r(i,l+1,k),
     +                 r(i,l+1,k+1)
      endif

!    ! stop

      return
      end subroutine pqr
c------------------------------------------------------------
	!1-�� ���� ���������� ���������� ����
      subroutine emh1(nt)
      implicit real*8(a-h,o-z)
       include 'part.pf'
c      include 'laser1.par'
c      include 'laser2.par'
c      include 'laser3.par'
c      parameter(imp=42,lmp=22,kmp=3)
      real*8 hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp),
     =ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp),
     =qx(imp,lmp,kmp),qy(imp,lmp,kmp),qz(imp,lmp,kmp)

      common/c/im,lm,km
      common/d/h1,h2,h3
      common/f/tau,c1,c2,c3
      common/g/hx,hy,hz
      common/h/ex,ey,ez
      common/o/qx,qy,qz
c********************************************************
      integer nt
      real*8 hx1(imp,lmp,kmp),hy1(imp,lmp,kmp),hz1(imp,lmp,kmp)
      common/dbgh/hx1,hy1,hz1
      do i = 1,imp
         do l = 1,lmp
            do k = 1,kmp
               hx1(i,l,k) = hx(i,l,k)
               hy1(i,l,k) = hy(i,l,k)
               hz1(i,l,k) = hz(i,l,k)
            enddo
         enddo
      enddo

c********************************************************

!      !!print 1010
 1010 format(2x,' ࠡ�⠥� emh1')
      im1=im+1
      lm1=lm+1
      km1=km+1
      im2=im+2
      lm2=lm+2
      km2=km+2

      c11=c1/beta0
      c21=c2/beta0
      c31=c3/beta0
      do 1 k=1,km1
      do 1 l=1,lm1
      do 1 i=1,im2
      qx(i,l,k)=0.5d0*(c31*(ey(i,l,k+1)-ey(i,l,k))-
     -c21*(ez(i,l+1,k)-ez(i,l,k)))
      hx(i,l,k)=hx(i,l,k)+qx(i,l,k)
      if(deb.ge.4) then
 218     format('QX ',4i4,7e30.20)
         write(37,218) nt,i,l,k,hx(i,l,k),hx1(i,l,k),
     +       qx(i,l,k),ey(i,l,k+1),ey(i,l,k),
     +       ez(i,l+1,k),ez(i,l,k)
      endif
    1 continue
      do 2 k=1,km1
      do 2 l=1,lm2
      do 2 i=1,im1
      qy(i,l,k)=0.5d0*(c11*(ez(i+1,l,k)-ez(i,l,k))-
     -c31*(ex(i,l,k+1)-ex(i,l,k)))
      hy(i,l,k)=hy(i,l,k)+qy(i,l,k)
      if(deb.ge.4) then
 219     format('QY ',4i4,7e30.20)
         write(37,219) nt,i,l,k,hy(i,l,k),hy1(i,l,k),
     +       qy(i,l,k),ez(i+1,l,k),ez(i,l,k),
     +       ex(i,l,k+1),ex(i,l,k)
      endif
    2 continue
    
      do 3 k=1,km2
      do 3 l=1,lm1
      do 3 i=1,im1
      qz(i,l,k)=0.5d0*(c21*(ex(i,l+1,k)-ex(i,l,k))-
     -c11*(ey(i+1,l,k)-ey(i,l,k)))
      hz(i,l,k)=hz(i,l,k)+qz(i,l,k)
      if(deb.ge.4) then
 220     format('QZ ',4i4,7e30.20)
         write(37,220) nt,i,l,k,hz(i,l,k),hz1(i,l,k),
     +       qz(i,l,k),ex(i,l+1,k),ex(i,l,k),
     +       ey(i+1,l,k),ey(i,l,k)
      endif
    3 continue
      return
      end subroutine emh1
c------------------------------------------------------------
	!2-�� ���� ���������� ���������� ����
      subroutine emh2(nt)
      implicit real*8(a-h,o-z)
      include 'part.pf'
c      include 'laser1.par'
c      include 'laser2.par'
c      include 'laser3.par'
c      parameter(imp=42,lmp=22,kmp=3)
      real*8 hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp),
     =qx(imp,lmp,kmp),qy(imp,lmp,kmp),qz(imp,lmp,kmp)
      real*8 hx1(imp,lmp,kmp),hy1(imp,lmp,kmp),hz1(imp,lmp,kmp)
      common/c/im,lm,km
      common/g/hx,hy,hz
      common/o/qx,qy,qz
      common/dbgh2/hx1,hy1,hz1
      integer nt
!      !!print 1010
 1010 format(2x,' ࠡ�⠥� emh2')
      write(37,*) 'in emh2'
      im1=im+1
      lm1=lm+1
      km1=km+1
      im2=im+2
      lm2=lm+2
      km2=km+2
      do 1 k=1,km1
      do 1 l=1,lm1
      do 1 i=1,im2
      hx(i,l,k)=hx(i,l,k)+qx(i,l,k)
      if(deb.ge.4) then
 229     format('Q2X ',4i4,3e30.20)
         write(37,229) nt,i,l,k,hx(i,l,k),hx1(i,l,k),
     +       qx(i,l,k)
      endif
    1 continue
      do 2 k=1,km1
      do 2 l=1,lm2
      do 2 i=1,im1
      hy(i,l,k)=hy(i,l,k)+qy(i,l,k)
      if(deb.ge.4) then
 218     format('Q2Y ',4i4,3e30.20)
         write(37,218) nt,i,l,k,hy(i,l,k),hy1(i,l,k),
     +       qy(i,l,k)
      endif

    2 continue
      do 3 k=1,km2
      do 3 l=1,lm1
      do 3 i=1,im1
      hz(i,l,k)=hz(i,l,k)+qz(i,l,k)
      if(deb.ge.4) then
 220     format('Q2Z ',4i4,7e30.20)
         write(37,220) nt,i,l,k,hz(i,l,k),hz1(i,l,k),
     +       qz(i,l,k)
      endif

    3 continue
      return
      end subroutine emh2
c------------------------------------------------------------
	!���������� �������������� ����
      subroutine eme(nt)
      implicit real*8(a-h,o-z)
      include 'part.pf'
c      include 'laser1.par'
c      include 'laser2.par'
c      include 'laser3.par'
c      parameter(imp=42,lmp=22,kmp=3)
      real*8 hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp),
     =ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp),
     =jx(imp,lmp,kmp),jy(imp,lmp,kmp),jz(imp,lmp,kmp)
      integer nt
      real*8 ex1(imp,lmp,kmp),de,ey1(imp,lmp,kmp),ez1(imp,lmp,kmp)
      common/c/im,lm,km
      common/d/h1,h2,h3
      common/f/tau,c1,c2,c3
      common/g/hx,hy,hz
      common/h/ex,ey,ez
      common/j/jx,jy,jz
      common/dbge/ex1,ey1,ez1
      integer lms,lmf
      character*40 fn
      character*20 tr,tr1,num
c b��c�eh�e �ektp��eck�x �o�e�
!      !!print 1010
      fn='ex'
      write(num,'(i5.5)') me
      fn=fn(1:2)//num
      fn=fn(1:7)//'nt'
      write(tr1,'(i3.3)') nt
      fn=fn(1:9)//tr1
      fn=fn(1:12)//'.dat'

 1010 format(2x,' ࠡ�⠥� eme')
      im1=im+1
      lm1=lm+1
      km1=km+1
      im2=im+2
      lm2=lm+2
      km2=km+2
c******************************************************
      if(magf.eq.0) then
         do i = 1,imp
            do l = 1,lmp
               do k = 1,kmp
    	          hx(i,l,k) = 0d0
	          hy(i,l,k) = 0d0
	          hz(i,l,k) = 0d0
	       enddo
            enddo    	    
         enddo
      endif

c******************************************************
      if(deb.ge.2) then
         write(37,855) tau,c1,c2,c3,h1,h2,h3,beta0
      endif
 855  format('eme ',8e10.3)
      do i = 1,imp
         do l = 1,lmp
            do k = 1,kmp
               ex1(i,l,k) = ex(i,l,k)
               ey1(i,l,k) = ey(i,l,k)
               ez1(i,l,k) = ez(i,l,k)
            enddo
         enddo
      enddo
c******************************************************
      c11=c1/beta0
      c21=c2/beta0
      c31=c3/beta0

      open(38,file=fn,form='formatted')

      do 1 k=2,km1
      do 1 l=2,lm1
      do 1 i=1,im1
      ex(i,l,k)=ex(i,l,k)+c21*(hz(i,l,k)-hz(i,l-1,k))-
     *c31*(hy(i,l,k)-hy(i,l,k-1))-tau*jx(i,l,k)
      if(deb.ge.4) then
 147     format('EX ',4i4,7e30.20)
         write(37,147) nt,i,l,k,ex(i,l,k),ex1(i,l,k),hz(i,l,k),
     +                    hz(i,l-1,k),
     +                    hy(i,l,k),hy(i,l,k-1),jx(i,l,k)
      endif
         write(38,147) nt,i,l,k,ex(i,l,k),
     +                    ex1(i,l,k),hz(i,l,k),
     +                    hz(i,l-1,k),
     +                    hy(i,l,k),hy(i,l,k-1),jx(i,l,k)
!     +                    l,meh,lmp-2,l+(lmp-2)*meh

      
    1 continue
      close(38)

      fn='ey'
      write(num,'(i5.5)') me
      fn=fn(1:2)//num
      fn=fn(1:7)//'nt'
      write(tr1,'(i3.3)') nt
      fn=fn(1:9)//tr1
      fn=fn(1:12)//'.dat'
!      fn=fn(1:16)//'.dat'
      open(38,file=fn,form='formatted')

      do 2 k=2,km1
      do 2 l=1,lm1
      do 2 i=2,im1
      ey(i,l,k)=ey(i,l,k)+c31*(hx(i,l,k)-hx(i,l,k-1))-
     *c11*(hz(i,l,k)-hz(i-1,l,k))-tau*jy(i,l,k)
      if(deb.ge.4) then
 148     format('EY ',4i4,7e30.20)
         write(37,148) nt,i,l,k,ey(i,l,k),ey1(i,l,k),hx(i,l,k),
     +                    hx(i,l,k-1),
     +                    hz(i,l,k),hz(i-1,l,k),jy(i,l,k)
      endif

         write(38,148) nt,i,l,k,ey(i,l,k),ey1(i,l,k),hx(i,l,k),
     +                    hx(i,l,k-1),
     +                    hz(i,l,k),hz(i-1,l,k),jy(i,l,k)
    2 continue   

      close(38)

      fn='ez'
      write(num,'(i5.5)') me
      fn=fn(1:2)//num
      fn=fn(1:7)//'nt'
      write(tr1,'(i3.3)') nt
      fn=fn(1:9)//tr1
      fn=fn(1:12)//'.dat'
      open(38,file=fn,form='formatted')
 

!    2 continue
      do 3 k=1,km1
      do 3 l=2,lm1
      do 3 i=2,im1
      ez(i,l,k)=ez(i,l,k)+c11*(hy(i,l,k)-hy(i-1,l,k))-
     *c21*(hx(i,l,k)-hx(i,l-1,k))-tau*jz(i,l,k)
      if(deb.ge.4) then
 149     format('EZ ',4i4,7e30.20)
         write(37,149) nt,i,l,k,ez(i,l,k),ez1(i,l,k),hy(i,l,k),
     +                    hy(i-1,l,k),
     +                    hx(i,l,k),hx(i,l-1,k),jz(i,l,k)
      endif
      write(38,149) nt,i,l,k,ez(i,l,k),ez1(i,l,k),hy(i,l,k),
     +                    hy(i-1,l,k),
     +                    hx(i,l,k),hx(i,l-1,k),jz(i,l,k)
    3 continue
      close(38)

!    3 continue
      do 4 k=2,km1
      do 4 i=1,im1
      ex(i,1,k)=ex(i,lm1,k)
      ex(i,lm2,k)=ex(i,2,k)
    4 continue
      do 5 l=1,lm2
      do 5 i=1,im1
      ex(i,l,1)=ex(i,l,km1)
      ex(i,l,km2)=ex(i,l,2)
    5 continue
      do 6 k=2,km1
      do 6 l=1,lm1
      ey(1,l,k)=0   !ey(im1,l,k)
      ey(im2,l,k)=0 !ey(2,l,k)
    6 continue
      do 7 l=1,lm1
      do 7 i=1,im2
      ey(i,l,1)=ey(i,l,km1)
      ey(i,l,km2)=ey(i,l,2)
    7 continue
      do 8 k=1,km1
      do 8 l=2,lm1
      ez(1,l,k)=0 !   ez(im1,l,k)
      ez(im2,l,k)=0 ! ez(2,l,k) 
    8 continue
      do 9 k=1,km1
      do 9 i=1,im2
      ez(i,1,k)=ez(i,lm1,k)
      ez(i,lm2,k)=ez(i,2,k)
    9 continue

      call PARAsendE1(ex,ey,ez)
      call PARAsendE2(ex,ez)
    
      return
      end subroutine eme
c------------------------------------------------------------
!      subroutine pr3(b,l1,m1,l2,m2,l3,m3)
!      include 'laser1.par'
!      include 'laser2.par'
!      include 'laser3.par'
!c      parameter(imp=42,lmp=22,kmp=3)
!      real*8 b(imp,lmp,kmp),v(8)
!      integer mv(8)
!      do 1 k=l3,m3
!      k1=k
!      !write(25,901) k1
!!  901 format(9x,5hc�o� ,i2)
!!      n1=l1
!!      n2=l1+7
!!    5 if(n2.gt.m1) n2=m1
!!      n3=n2-n1+1
!!      do 2 i=n1,n2
!!    2 mv(i+1-n1)=i
!!      !write(25,902) (mv(i),i=1,n3)
!!  902 format(8x,8(i2,7x))
!!      do 4 l0=l2,m2
!!      l=m2-l0+1
!!      l02=l
!!      do 3 i=n1,n2
!!    3 v(i+1-n1)=b(i,l,k)
!!      !write(25,903) l02,(v(i),i=1,n3)
!!  903 format(1x,i2,1x,8e9.3)
!    4 continue
!      if(n2.eq.m1) goto 1
!      n1=n1+8
!      n2=n2+8
!      goto 5
!    1 continue
!      return
!      end subroutine pr3
c------------------------------------------------------
	!������� �.�., �������������� �� ����������� ������  
      real*8 function g05dde(a,b)
      implicit real*8 (a-h,o-z)
      real*8 a, b
      real*8 store1, store2
      real*8 half, one, t, u, v, w, x
      integer n
      real*8 d(41)
      real*8 wrapg05cae
      common /cag05b/ store1, store2
      data one /1.0d0/, half /0.5d0/
      data d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), d(9),
     * d(10), d(11), d(12), d(13), d(14) /0.0d0,0.674489750196082d0,
     * 1.150349380376008d0,1.534120544352546d0,1.862731867421652d0,
     * 2.153874694061456d0,2.417559016236505d0,2.660067468617460d0,
     * 2.885634912426757d0,3.097269078198785d0,3.297193345691964d0,
     * 3.487104104114431d0,3.668329285121323d0,3.841930685501911d0/
      data d(15), d(16), d(17), d(18), d(19), d(20), d(21), d(22),
     * d(23), d(24), d(25), d(26), d(27) /4.008772594168585d0,
     * 4.169569323349106d0,4.324919040826046d0,4.475328424654204d0,
     * 4.621231001499247d0,4.763001034267814d0,4.900964207963193d0,
     * 5.035405969463927d0,5.166578119728753d0,5.294704084854598d0,
     * 5.419983174916868d0,5.542594057802940d0,5.662697617459439d0/
      data d(28), d(29), d(30), d(31), d(32), d(33), d(34), d(35),
     * d(36), d(37), d(38), d(39), d(40) /5.780439324478935d0,
     * 5.895951216739571d0,6.009353565530745d0,6.120756285971941d0,
     * 6.230260137989044d0,6.337957754553790d0,6.443934526538564d0,
     * 6.548269367831731d0,6.651035379893011d0,6.752300431407015d0,
     * 6.852127665896068d0,6.950575947916750d0,7.047700256664409d0/
      data d(41) /7.143552034352190d0/
      u = store1
      do 20 n=1,39
         if (u.gt.half) go to 40
         u = u + u
   20 continue
      n = 40
   40 t = d(n)
      u = wrapg05cae(x)
   60 w = (d(n+1)-t)*u
      v = w*(w*half+t)
   80 u = wrapg05cae(x)
      if (v.le.u) go to 100
      v = wrapg05cae(x)
      if (u.gt.v) go to 80
      u = (v-u)/(one-u)
      go to 60
  100 u = (u-v)/(one-v)
      if (u.gt.half) go to 120
      store1 = u + u
      g05dde = a + b*(w+t)
      return
  120 store1 = u + u - one
      g05dde = a - b*(w+t)
      return
      end

!------------------------------------------------------------------
      real*8 function g05cae(x)
      implicit real*8 (a-h,o-z)
      real*8 x
      integer ix,iy,iz
      common /cag05a/ ix,iy,iz
      ix = 171*(ix-(ix/177)*177) - 2*(ix/177)
      iy = 172*(iy-(iy/176)*176) - 2*(iy/176)
      iz = 170*(iz-(iz/178)*178) - 2*(iz/178)
      if (ix.lt.0) ix = ix+30269
      if (iy.lt.0) iy = iy+30307
      if (iz.lt.0) iz = iz+30323
      ax=ix
      ay=iy
      az=iz
      ai = ax/30269.0d0+ay/30307.0d0+az/30323.0d0
      ii=ai
      g05cae = ai-ii
      return
      end
      block data
      integer ix,iy,iz
      real*8 normal,gamma
      common /cag05a/ ix,iy,iz
      common /cag05b/ normal,gamma
      data ix/1/,iy/255/,iz/25555/
      data normal/1.0d0/,gamma/-1.0d0/
      end
!------------------------------------------------------------------
	!������� �������� ������ � �����
      integer function ich(st)
      integer i,k
      character*2 st
      character*1 a
      a=st(1:1)
      i=ichar(a)-48
      a=st(2:2)
      k=ichar(a)-48
      ich=10*i+k
      return
      end
!------------------------------------------------------------------
	!������� �������� ����� � ������
      character*2 function ci(i)
      integer i,j,k
      character*1 a,b
      k=i/10
      j=i-10*k+48
      k=k+48
      a=char(k)
      b=char(j)
      ci=a//b
      return
      end
      
      subroutine sort(xi,yi,zi,ui,vi,wi,jm)
      implicit none
      include 'part.pf'
       real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 x,y,z,px,py,pz
      integer jm,i,j

      do i = 1,jm
         do j = i+1,jm
            if(yi(i).gt.yi(j)) then
               x = xi(i)
               y = yi(i)
               z = zi(i)
               px = ui(i)
               py = vi(i)
               pz = wi(i)

               xi(i) = xi(j)
               yi(i) = yi(j)
               zi(i) = zi(j)
               ui(i) = ui(j)
               vi(i) = vi(j)
               wi(i) = wi(j)

               xi(j) = x
               yi(j) = y
               zi(j) = z
               ui(j) = px
               vi(j) = py
               wi(j) = pz

            endif
         enddo
      enddo


      end subroutine sort
      
