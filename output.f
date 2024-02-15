c     sort 0 - electrons, 1 - ions, 2 - beam,22 - extra beam particles added at each timestep
      subroutine writeallparticles(sort,x,y,z,pu,pv,pw,jm,m7)
      implicit none
      include 'part.pf'
      real*8 getMidHarmonic
c      integer m33,nt3
      character*40 tr,tr1
c      common/flow/etfl0,itfl0
c      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
c      common/vvv/vnum,i1,i2,i33,i4	
c      common/j/jx,jy,jz
c      real*8 vel(3)
      real *8 x(jmp),y(jmp),z(jmp),pu(jmp),pv(jmp),pw(jmp)
      integer i,j,sort,jm,m7
      
      if(sort.eq.0) then
         tr='elpp'
      endif
      if(sort.eq.1) then
         tr='inpp'
      endif
      if(sort.eq.2) then
         tr='bmpp'
      endif

      if(sort.eq.22) then
         tr='badd'
      endif

      if(sort.eq.23) then
         tr='bfly'
      endif
      
      write(tr1,'(i3.3)') meh
      tr=tr(1:4)//tr1
c      print*,,tr
      write(tr1,'(i5.5)') m7
c      print*,stage,tr,tr1
      tr=tr(1:7)//tr1 
c      print*,stage,tr,tr1
      tr=tr(1:12)//'.dat' 
c      print*,stage,tr,tr1
      open(94,file=tr,form='formatted')
      
  127 format(i10,6e30.20)      
      do j=1,jm
         write(94,127) j,x(j),y(j),z(j),pu(j),pv(j),pw(j)
      enddo
      
      close(94)
      
      end

      subroutine writeGridarray(fd,a)
      implicit none
      include 'part.pf'
      real*8 a(imp,lmp,kmp)
      integer fd,i,l,k
      
      write(fd) (((a(i,l,k),k=1,kmp),l=1,lmp),i=1,imp)

      

     
      end

      subroutine writeParticlearray(fd,a,jm)
      implicit none
      include 'part.pf'
      real*8 a(jmp)
      integer fd,i,jm
      write(37,*) 'particles ',jm,a(1)
      write(fd) (a(i),i=1,jm)
     
      end
      
      subroutine cleanparticles(x,y,z,pu,pv,pw,
     +                xmin,xmax,ymin,ymax,zmin,zmax,jm,nt,sort,wh)
      implicit none
      include 'part.pf'
      real*8 x(jmp),y(jmp),z(jmp),pu(jmp),pv(jmp),pw(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),puf(jmp),pvf(jmp),pwf(jmp)
      real*8 x1(jmp),y1(jmp),z1(jmp),pu1(jmp),pv1(jmp),pw1(jmp)
      real*8 xmin,xmax,ymin,ymax,zmin,zmax
      integer j,jm,jm1,i,nt,k
      character*2 sort
      character*4 wh

 660  format('begin_clean ',6e10.3,i10) 
      write(37,660) xmin,xmax,ymin,ymax,zmin,zmax,jm

      
      i=0
      k=0
      do j = 1,jm
 661     format('in_clean ',i10,9e30.20)
         if(deb.ge.2) then
            write(37,661) j,xmin,x(j),xmax,
     +                      ymin,y(j),ymax,zmin,z(j),zmax
         endif
         
         if(((x(j).gt.xmin).and.(x(j).lt.xmax)).and.
     +      ((y(j).gt.ymin).and.(y(j).lt.ymax)).and.
     +      ((z(j).gt.zmin).and.(z(j).lt.zmax))) then
        
            if(deb.ge.2) write(37,*) j

            i = i + 1
            x1(i) = x(j)
            y1(i) = y(j)
            z1(i) = z(j)
            pu1(i) = pu(j)
            pv1(i) = pv(j)
            pw1(i) = pw(j)
         else
            k = k + 1
            xf(k)  = x(j)
            yf(k)  = y(j)
            zf(k)  = z(j)
            puf(k) = pu(j)
            pvf(k) = pv(j)
            pwf(k) = pw(j)
         endif
      enddo
      jm=i
      print*,'GOT fly LISTS ',me
      print*,'sort ',sort
      call writeParticleList(xf,yf,zf,puf,pvf,pwf,nt,k,wh,sort)
      print*,'WRITE fly LISTS ',me
      
      do j = 1,jm
            x(j)  = x1(j)
            y(j)  = y1(j)
            z(j)  = z1(j)
            pu(j) = pu1(j)
            pv(j) = pv1(j)
            pw(j) = pw1(j)
 662     format('clean_assign ',i10,9e30.20)
         if(deb.ge.2) then
            write(37,662) j,xmin,x(j),xmax,
     +                      ymin,y(j),ymax,zmin,z(j),zmax
         endif
      enddo
      print*,'END fly LISTS ',me
      
      
      
      end      


      subroutine writeallbinaryoutput(stage,nt)
      implicit none
      include 'part.pf'
      include 'mass.par'
      integer lp,jmf,jmi,n1,n2,m5,m7,m27,jmb,nt1,nt
      integer vnum,i1,i2,i33,i4,ki,kl,kk,kid,kld,kkd
      integer kih,klh,kkh,mee27,mhh27,mne27,mni27
      real*8 fd(18),tp(4),tp0(4),y0p,pi,ami,amf,krd,exmd,krh,exmh,nueff
      real*8 enb,dnb,see27,shh27,sne27,sni27,nu
      real*8 dne(imp,lmp,kmp),dni(imp,lmp,kmp),dTe(imp,lmp,kmp),
     +       dTi(imp,lmp,kmp),kr,exm,dee(imp,lmp,kmp)
      real*8 dne2(imp,lmpf,kmp),dni2(imp,lmpf,kmp),dTe2(imp,lmpf,kmp),
     +       dTi2(imp,lmpf,kmp),dee2(imp,lmpf,kmp)	   	   
      real*8 dnbm(imp,lmp,kmp),dnb2(imp,lmpf,kmp),amb
      real*8 dhh(imp,lmp,kmp),dhh2(imp,lmpf,kmp)
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp),s27
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      real*8 etfl2,etfl3,etfl0,etfl,itfl2,itfl3,itfl0,itfl,mhd
      real*8 dTeX(imp,lmpf,kmp),DTeY(imp,lmpf,kmp),DTeZ(imp,lmpf,kmp)
      real*8 dVX(imp,lmpf,kmp),DVY(imp,lmpf,kmp),DVZ(imp,lmpf,kmp)
      real*8 jx(imp,lmp,kmp),jy(imp,lmp,kmp),jz(imp,lmp,kmp)
      integer stage
c================================================================
      real*8 ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp),
     *hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp)
      common/g/hx,hy,hz
      common/h/ex,ey,ez
c================================================================      
      real*8 qx(imp,lmp,kmp),qy(imp,lmp,kmp),qz(imp,lmp,kmp)
      common/o/qx,qy,qz
      
      real*8 getMidHarmonic
      integer m33,nt3
      character*40 tr,tr1
      common/flow/etfl0,itfl0
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/iff/xf,yf,zf,uf,vf,wf,amf
      common/vvv/vnum,i1,i2,i33,i4	
      common/j/jx,jy,jz
      real*8 vel(3)
      integer i
      
      if(write_all.eq.0) then
         return
      endif
      
      tr='mumu'
      write(tr1,'(i3.3)') stage
      tr=tr(1:4)//tr1
      print*,stage,tr
      write(tr1,'(i5.5)') nt
      print*,stage,tr,tr1
      tr=tr(1:7)//tr1 
      print*,stage,tr,tr1
      tr=tr(1:12)//'.dat' 
      print*,stage,tr,tr1
      open(94,file=tr,form='unformatted')

      call writegridarray(94,ex)
      call writegridarray(94,ey)
      call writegridarray(94,ez)

      call writegridarray(94,hx)
      call writegridarray(94,hy)
      call writegridarray(94,hz)

      call writegridarray(94,jx)
      call writegridarray(94,jy)
      call writegridarray(94,jz)

      vel(1) = jmi
      vel(2) = mel/mi
      vel(3) = ami
      write(94) (vel(i),i=1,3)
      call writeparticlearray(94,xi,jmi)
      call writeparticlearray(94,yi,jmi)
      call writeparticlearray(94,zi,jmi)
      call writeparticlearray(94,ui,jmi)
      call writeparticlearray(94,vi,jmi)
      call writeparticlearray(94,wi,jmi)

      vel(1) = jmf
      vel(2) = -1.0d0
      vel(3) = amf
      write(94) (vel(i),i=1,3)

 111  format('x,y,z 111 ',3e30.20)
      write(37,111) xf(1),yf(1),zf(1)
      call writeparticlearray(94,xf,jmf)
      call writeparticlearray(94,yf,jmf)
      call writeparticlearray(94,zf,jmf)
      call writeparticlearray(94,uf,jmf)
      call writeparticlearray(94,vf,jmf)
      call writeparticlearray(94,wf,jmf)

      vel(1) = jmb
      vel(2) = -1.0d0
      vel(3) = amb
      write(94) (vel(i),i=1,3)
      call writeparticlearray(94,xb,jmb)
      call writeparticlearray(94,yb,jmb)
      call writeparticlearray(94,zb,jmb)
      call writeparticlearray(94,ub,jmb)
      call writeparticlearray(94,vb,jmb)
      call writeparticlearray(94,wb,jmb)



      close(94)

         call four3D(ex,dni2,imp,lmp,kmp,nt,
     +	 'muex','frex',fd(10),fd(11),fd(12))


      end subroutine writeallbinaryoutput


      subroutine writecurrents(j1,nt)
      implicit none
      include 'part.pf'
      include 'mass.par'
      integer lp,jmf,jmi,n1,n2,m5,m7,m27,jmb,nt1,nt
      integer vnum,i1,i2,i33,i4,ki,kl,kk,kid,kld,kkd
      integer kih,klh,kkh,mee27,mhh27,mne27,mni27
      real*8 fd(18),tp(4),tp0(4),y0p,pi,ami,amf,krd,exmd,krh,exmh,nueff
      real*8 enb,dnb,see27,shh27,sne27,sni27,nu
      real*8 dne(imp,lmp,kmp),dni(imp,lmp,kmp),dTe(imp,lmp,kmp),
     +       dTi(imp,lmp,kmp),kr,exm,dee(imp,lmp,kmp)
      real*8 dne2(imp,lmpf,kmp),dni2(imp,lmpf,kmp),dTe2(imp,lmpf,kmp),
     +       dTi2(imp,lmpf,kmp),dee2(imp,lmpf,kmp)	   	   
      real*8 dnbm(imp,lmp,kmp),dnb2(imp,lmpf,kmp),amb
      real*8 dhh(imp,lmp,kmp),dhh2(imp,lmpf,kmp)
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp),s27
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      real*8 etfl2,etfl3,etfl0,etfl,itfl2,itfl3,itfl0,itfl,mhd
      real*8 dTeX(imp,lmpf,kmp),DTeY(imp,lmpf,kmp),DTeZ(imp,lmpf,kmp)
      real*8 dVX(imp,lmpf,kmp),DVY(imp,lmpf,kmp),DVZ(imp,lmpf,kmp)
      real*8 jx(imp,lmp,kmp),jy(imp,lmp,kmp),jz(imp,lmp,kmp)
      integer stage,j1
c================================================================
      real*8 ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp),
     *hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp)
      common/g/hx,hy,hz
      common/h/ex,ey,ez
c================================================================      
      real*8 qx(imp,lmp,kmp),qy(imp,lmp,kmp),qz(imp,lmp,kmp)
      common/o/qx,qy,qz
      
      real*8 getMidHarmonic
      integer m33,nt3
      character*40 tr,tr1
      common/flow/etfl0,itfl0
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/iff/xf,yf,zf,uf,vf,wf,amf
      common/vvv/vnum,i1,i2,i33,i4	
      common/j/jx,jy,jz
      real*8 vel(3)
      integer i
      
      tr='curr'
      write(tr1,'(i3.3)') nt
      tr=tr(1:4)//tr1
      print*,stage,tr
      write(tr1,'(i5.5)') j1
      print*,stage,tr,tr1
      tr=tr(1:7)//tr1 
      print*,stage,tr,tr1
      tr=tr(1:12)//'.dat' 
      print*,stage,tr,tr1
      open(94,file=tr,form='unformatted')

      call writegridarray(94,jx)
      call writegridarray(94,jy)
      call writegridarray(94,jz)

      end subroutine writecurrents

      subroutine output(y0p,nt1,m33,nt3,label)
      implicit none
      include 'part.pf'
      integer lp,jmf,jmi,n1,n2,m5,m7,m27,jmb,nt1
      integer vnum,i1,i2,i33,i4,ki,kl,kk,kid,kld,kkd
      integer kih,klh,kkh,mee27,mhh27,mne27,mni27,label
      real*8 fd(18),tp(4),tp0(4),y0p,pi,ami,amf,krd,exmd,krh,exmh,nueff
      real*8 enb,dnb,see27,shh27,sne27,sni27,nu
      real*8 dne(imp,lmp,kmp),dni(imp,lmp,kmp),dTe(imp,lmp,kmp),
     +       dTi(imp,lmp,kmp),kr,exm,dee(imp,lmp,kmp)
      real*8 dne2(imp,lmpf,kmp),dni2(imp,lmpf,kmp),dTe2(imp,lmpf,kmp),
     +       dTi2(imp,lmpf,kmp),dee2(imp,lmpf,kmp)	   	   
      real*8 dnbm(imp,lmp,kmp),dnb2(imp,lmpf,kmp),amb
      real*8 dhh(imp,lmp,kmp),dhh2(imp,lmpf,kmp)
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp),s27
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      real*8 etfl2,etfl3,etfl0,etfl,itfl2,itfl3,itfl0,itfl,mhd
      real*8 dTeX(imp,lmpf,kmp),DTeY(imp,lmpf,kmp),DTeZ(imp,lmpf,kmp)
      real*8 dVX(imp,lmpf,kmp),DVY(imp,lmpf,kmp),DVZ(imp,lmpf,kmp)
      real*8 jx(imp,lmp,kmp),jy(imp,lmp,kmp),jz(imp,lmp,kmp)
c================================================================
      real*8 ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp),
     *hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp)
      common/g/hx,hy,hz
      common/h/ex,ey,ez
c================================================================      
      real*8 qx(imp,lmp,kmp),qy(imp,lmp,kmp),qz(imp,lmp,kmp)
      common/o/qx,qy,qz


      
      real*8 getMidHarmonic
      integer m33,nt3
      character*40 tr,tr1
      common/flow/etfl0,itfl0
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/iff/xf,yf,zf,uf,vf,wf,amf
      common/vvv/vnum,i1,i2,i33,i4	
      common/j/jx,jy,jz
      

      call ChainOut(m7)
      
      if(deb.ge.2) write(37,*) 'in output ',y0p,m7
c      if(m5.gt.2) stop
      if(deb.ge.4) write(37,*) 'begin out ',xi(jdeb)
          
      if(m7.eq.0) then   
         if(deb.ge.2) write(37,*) 'b temp ',ons,me,m7
c        if(m5.gt.2) stop
            if(out3d.eq.1) then
               call tempFlow3d(xf,yf,zf,uf,vf,wf,y0p,jmf)
            endif   
            
         if(deb.ge.2) write(37,*) 'b temp,fl0,fl0i ',ons,me,m7,y0p,
     +     etfl0,itfl0

         if((ons.eq.0).and.(me.eq.0)) then
            tr = 'type'    
            write(tr1,'(i3.3)') vnum 
            tr=tr(1:4)//tr1 
            tr=tr(1:7)//'.dat' 
            open(94,file=tr,form='formatted')
         endif	 
      endif      

            if(out3d.eq.1) then
               call tempFlow3d(xf,yf,zf,uf,vf,wf,y0p,jmf)
            endif   


      if(deb.ge.2) write(37,*) 'a type ',m5
      if(deb.ge.2) write(37,*) 'temp1 ',m5

      if(wpbeam.eq.1) then
         call writeallparticles(2,xb,yb,zb,ub,vb,wb,jmb,m7)
      endif
      if(wpelec.eq.1) then
         call writeallparticles(0,xf,yf,zf,uf,vf,wf,jmf,m7)
      endif
      if(wpion.eq.1) then
         call writeallparticles(1,xi,yi,zi,ui,vi,wi,jmi,m7)
      endif

c      stop
      if(fte.eq.1) then
         call tempSort(xf,yf,zf,uf,vf,wf,y0p,jmf,
     +	               dTeX,DTeY,DTeZ,dVX,DVY,DVZ)

         call four3D(dTeX,dTe2,imp,lmp,kmp,m7,
     +	 'dnTX','frTe',fd(1),fd(2),fd(3))
         call four3D(dTeY,dTe2,imp,lmp,kmp,m7,
     +	 'dnTY','frTe',fd(1),fd(2),fd(3))
         call four3D(dTeZ,dTe2,imp,lmp,kmp,m7,
     +	 'dnTZ','frTe',fd(1),fd(2),fd(3))
         call four3D(dVX,dTe2,imp,lmp,kmp,m7,
     +	 'dnVX','frTe',fd(1),fd(2),fd(3))
         call four3D(dVY,dTe2,imp,lmp,kmp,m7,
     +	 'dnVY','frTe',fd(1),fd(2),fd(3))
         call four3D(dVZ,dTe2,imp,lmp,kmp,m7,
     +	 'dnVZ','frTe',fd(1),fd(2),fd(3))
      endif
      if(deb.ge.2) write(37,*) 'dnTe ',m5,fd(1),fd(2),fd(3)
c      stop

      call densSort(xf,yf,zf,y0p,jmf,amf,0,dne,label)

      if(deb.ge.2) write(37,*) 'dnne ',m5
c      stop
      if(fne.eq.1) then
         call four3D(dne,dne2,imp,lmp,kmp,m7,
     +	 'dnne','frne',fd(4),fd(5),fd(6))
c         call getMaxHarm(dne2,kid,kld,kkd,krd,exmd)
c         call selectHarmList(dne2,'harmdne.dat',87)
      endif 	 

      if(fnb.eq.1) then
         call densSort(xb,yb,zb,y0p,jmb,amb,0,dnbm,label)

         call four3D(dnbm,dnb2,imp,lmp,kmp,m7,
     +	 'dnbm','frbm',fd(4),fd(5),fd(6))
c         call getMaxHarm(dne2,kid,kld,kkd,krd,exmd)
c         call selectHarmList(dnb2,'harmdnb.dat',88)
      endif 	 
c      call tempSort(xi,yi,zi,ui,vi,wi,y0p,jmi,dTi,tp(3),tp(4),
c     +      itfl2,itfl3,itfl0,itfl)
      if(fti.eq.1) then	  	    
         call four3D(dTi,dTi2,imp,lmp,kmp,m7,'dnTi','frTi',
     +          	 fd(7),fd(8),fd(9))
      endif	 

      call densSort(xi,yi,zi,y0p,jmi,ami,0,dni,label)
      if(fni.eq.1) then
         call four3D(dni,dni2,imp,lmp,kmp,m7,
     +	 'dnni','frni',fd(10),fd(11),fd(12))
c         call selectHarmList(dni2,'harmdni.dat',89)
      endif 	 
        
      if(deb.ge.2) write(37,*) 'dnTi ',m5
       
      call densField(dee,0)
      if(fee.eq.1) then
         call four3D(dee,dee2,imp,lmp,kmp,m7,
     +	 'dnee','free',fd(13),fd(14),fd(15))
c         mhd = getMidHarmonic(dee2) 
c         call getMaxHarm(dee2,ki,kl,kk,kr,exm)
c         call selectHarmList(dee2,'harmelf.dat',90)
      endif 	 
         
         
      call densField(dhh,1)
      if(fhh.eq.1) then
         call four3D(dhh,dhh2,imp,lmp,kmp,m7,
     +	 'dnhh','frhh',fd(16),fd(17),fd(18))
c         call getMaxHarm(dhh2,kih,klh,kkh,krh,exmh)
      endif 	 
          
      if(radf.eq.1) then
         call avRadius(dne,dni,dTe,dTi,y0p)
      endif	 
     
      if(map.eq.1) then
         call outDens2D(y0p)
      endif    
         call outDens1D(y0p)
c         call phasePlane
      
      if(disf.eq.1) then
         call distrFunc1D(0,1)
c         call distrFunc2D(0,1)
      endif
       
      call getBeamEnergy(enb,dnb)
      
c      call getAnomFreq(nueff)
      nueff = 0d0
      call getFreq(nu,m7)
                  
      if(tp0(1).gt.1d-8) then
         tp(1) = tp(1)/tp0(1)
      endif

      if(tp0(2).gt.1d-8) then
         tp(2) = tp(2)/tp0(2)
      endif
        
      if(tp0(3).gt.1d-8) then
         tp(3) = tp(3)/tp0(3)
      endif
        
      if(tp0(4).gt.1d-8) then
         tp(4) = tp(4)/tp0(4)
      endif
      
      if (m33.eq.nt3)then 
       call printv
       m33=0
      endif    
      end subroutine output
!-------------------------------------------------------------------------------------
      subroutine avRadius(dne,dni,dTe,dTi,y0)
      implicit none
      include 'part.pf'
      integer lp,jmf,jmi,n1,n2,m5,m7,i,l,k
      real*8 dne(imp,lmp,kmp),dni(imp,lmp,kmp),dTe(imp,lmp,kmp),
     +       dTi(imp,lmp,kmp),pi,ami,amf,t,h1,h2,h3   
      real*8 rde(imp),rle(imp),rdi(imp),rli(imp),dte1(imp),dti1(imp),
     +       dne1(imp),dni1(imp),y0 
      character*40 tr,tr1
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/iff/xf,yf,zf,uf,vf,wf,amf
      common/d/h1,h2,h3
      
      do i  = 1,imp
         rde(i)   = 0d0 
         rdi(i)   = 0d0 
         dte1(i)  = 0d0 
         dti1(i)  = 0d0 
         dni1(i)  = 0d0 
         dne1(i)  = 0d0 
      enddo
      
      do i = 1,imp
         do l = 1,lmp 
            do k = 1,kmp
               if(dabs(dne(i,l,k)).gt.1.0d-8) then
	          rde(i) = rde(i) + dsqrt(dTe(i,l,k)/dabs(dne(i,l,k)))
	       else
	          rde(i) = 0d0    	  
	       endif                        	    

               if(dabs(dni(i,l,k)).gt.1.0d-8) then
	          rdi(i) = rdi(i) + dsqrt(dTi(i,l,k)/dabs(dni(i,l,k)))
	       else
	          rdi(i) = 0d0    	  
	       endif                        	    
               dte1(i) = dte1(i) + dTe(i,l,k)       
               dti1(i) = dti1(i) + dTi(i,l,k)       
               dne1(i) = dne1(i) + dne(i,l,k)       
               dni1(i) = dni1(i) + dni(i,l,k)       
            enddo   	    	     
      	 enddo
	 rde(i)  = rde(i)/lmp/kmp
	 rdi(i)  = rdi(i)/lmp/kmp
	 dne1(i) = dne1(i)/lmp/kmp
	 dni1(i) = dni1(i)/lmp/kmp
	 dte1(i) = dte1(i)/lmp/kmp
	 dti1(i) = dti1(i)/lmp/kmp
      enddo 
      
      call getLarmor(xf,yf,zf,uf,vf,wf,y0,jmf,rle)
      call getLarmor(xi,yi,zi,ui,vi,wi,y0,jmi,rli)
       
      call PARAreduce8(rde,rdi,dne1,dni1,dte1,dti1,rli,rle)
       
      if((outf.eq.1).and.(me.eq.0)) then
         tr = 'rdrl'    
         write(tr1,'(i3.3)') m7 
         tr=tr(1:4)//tr1 
         tr=tr(1:7)//'.dat' 
         open(95,file=tr,form='formatted')
 156     format(e10.3,8e15.5)        	 
         do i = 1,imp
	    write(95,156) i*h1,dte1(i),dne1(i),rde(i),dti1(i),
     +     	          dni1(i),rdi(i),rle(i),rli(i)
	 enddo	 
      endif
      
      end subroutine avRadius
!------------------------------------------------------------------------------      
      real*8 function minreal(a,b)
      implicit none
      real*8 a,b
      
      if(a.lt.b) then
         minreal = a
      else
         minreal = b
      endif	  	  
       	  
      end
!-----------------------------------------------------------------------------
      subroutine getLarmor(xf,yf,zf,uf,vf,wf,y0,jm,rle)     
      include 'part.pf'
      integer lp,jmf,jmi,n1,n2,m5,m7,i,l,k
      real*8 vx(imp,lmp,kmp),vy(imp,lmp,kmp),vz(imp,lmp,kmp)
      real*8 hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      real*8 rle(imp),rx,ry,rz,rl,minreal,y0
      common/g/hx,hy,hz

      if(deb.ge.2) write(37,*) 'in getLarmor'
      
      call getVelDens(xf,yf,zf,uf,vf,wf,y0,jm,uf,vx)
      call getVelDens(xf,yf,zf,uf,vf,wf,y0,jm,vf,vy)
      call getVelDens(xf,yf,zf,uf,vf,wf,y0,jm,wf,vz)
      
      do i = 1,imp
         rl = 0d0
         do l = 1,lmp
            do k = 1,kmp
	       if(dabs(hx(i,l,k)).gt.1.0d-8) then
	          rx = dsqrt(vy(i,l,k)**2+vz(i,l,k)**2)/dabs(hx(i,l,k))
	       else
	          rx = 0d0 
	       endif	  
	       if(dabs(hy(i,l,k)).gt.1.0d-8) then
	          ry = dsqrt(vx(i,l,k)**2+vz(i,l,k)**2)/dabs(hy(i,l,k))
	       else
	          ry = 0d0
	       endif    
	       if(dabs(hz(i,l,k)).gt.1.0d-8) then	  
	          rz = dsqrt(vy(i,l,k)**2+vx(i,l,k)**2)/dabs(hz(i,l,k))
     	       else    	  
	          rz = 0d0
	       endif        	  
	       rl = rl + minreal(rx,minreal(ry,rz))
            enddo	  
         enddo
	 rle(i) = rl/lmp/kmp
	 if(deb.ge.4) write(37,*) 'i,rle(i) ',i,rle(i)
      enddo
      
      if(deb.ge.2) write(37,*) 'end getLarmor'

      end subroutine getLarmor
!-------------------------------------------------------------------------------
      subroutine getVelDens(xf,yf,zf,uf,vf,wf,y0,jm,vv,px)
      implicit none
      include 'part.pf'
      real*8 px(imp,lmp,kmp),s2,s4,s6,a,s21,s41,s61,s
      real*8 xf(jmp),yf(jmp),zf(jmp),ps,pu,pv,pw,y0
      real*8 uf(jmp),vf(jmp),wf(jmp),vv(jmp),h1,h2,h3
      integer jm,i,k,l,j
      common/d/h1,h2,h3
      
      if(deb.ge.2) write(37,*) 'in getVelDens'
      
      do i = 1,imp
         do l = 1,lmp
            do k = 1,kmp
               px(i,l,k) = 0d0
            enddo	  
         enddo
      enddo
      
      do j = 1,jm
         pu = uf(j)
         pv = vf(j)
         pw = wf(j)
         ps=1.d0/dsqrt(1.d0+beta0**2*(pu*pu+pv*pv+pw*pw))
         a=ps*vv(j)
      
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
         
      	 s=s21*s41*a
         px(i,l,k)=px(i,l,k)+s*s61
         px(i,l,k+1)=px(i,l,k+1)+s*s6
         s=s21*s4*a
         px(i,l+1,k)=px(i,l+1,k)+s*s61
         px(i,l+1,k+1)=px(i,l+1,k+1)+s*s6
         s=s2*s41*a
         px(i+1,l,k)=px(i+1,l,k)+s*s61
         px(i+1,l,k+1)=px(i+1,l,k+1)+s*s6
         s=s2*s4*a
         px(i+1,l+1,k)=px(i+1,l+1,k)+s*s61
         px(i+1,l+1,k+1)=px(i+1,l+1,k+1)+s*s6
      enddo       

      if(deb.ge.2) write(37,*) 'end getVelDens ',px(2,2,2)

      end
!----------------------------------------------------------------------------------      
      subroutine partj(num,srt,j,a)
      implicit none
      include 'part.pf'
      real*8 s2,s4,s6,s21,s41,s61,s,a,ms,ami,amf,h1,h2,h3
      integer jm,i,k,l,j,m7,i2,l2,num,srt,lp,jmf,jmi,n1,n2,m5
      real*8 pi,pu,pv,pw,ps,q0
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp),amb
      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/iff/xf,yf,zf,uf,vf,wf,amf
      common/d/h1,h2,h3
      
      if(srt.eq.0) then
         pu = uf(j)
         pv = vf(j)
         pw = wf(j)
         ps=1.d0/dsqrt(1.d0+beta0**2*(pu*pu+pv*pv+pw*pw))
	 ms = amf
      else
         if(srt.eq.1) then
            pu = ui(j)
            pv = vi(j)
            pw = wi(j)
            ps=1.d0/dsqrt(1.d0+beta0**2*(pu*pu+pv*pv+pw*pw))
   	    ms = ami
   	 endif
      endif
      
      if(srt.eq.2) then
         pu = ub(j)
         pv = vb(j)
         pw = wb(j)
         ps=1.d0/dsqrt(1.d0+beta0**2*(pu*pu+pv*pv+pw*pw))
	 ms = amb
      endif	    
      
      
      if(num.eq.0) then
         a = ms
      endif
      if(num.eq.1) then
         a = ps*pu
      endif
      if(num.eq.2) then
         a = ps*pv 
      endif
      if(num.eq.3) then
         a = ps*pw 
      endif
      if(num.eq.4) then
         a = ps**2*(pu*pu+pv*pv+pw*pw)/2d0*ps*pu 
      endif
      if(num.eq.5) then
         a = ps**2*(pu*pu+pv*pv+pw*pw)/2d0*ps*pv 
      endif
      if(num.eq.6) then
         a = ps**2*(pu*pu+pv*pv+pw*pw)/2d0*ps*pw 
      endif

      if(num.eq.7) then
         if((pu*pu+pv*pv+pw*pw).gt.1.0d-8) then
         
            a = dasin(dsqrt((pv*pv+pw*pw)/(pu*pu+pv*pv+pw*pw)))
            if(deb.ge.7) then
               write(37,*) 'tu ',j,a,pv*pv+pw*pw,pu*pu+pv*pv+pw*pw,pu*pu
            endif
         else
            a = 0d0
         endif              
      endif
      
      if((abs(j-jmf/2).lt.2).and.(deb.ge.4)) then
         write(37,*) 'partj j,jmf,a ',j,jmf,a
      endif	 
      
      end
!------------------------------------------------------------------------------------------------
c     px  - resulting matrix
c     num - code of quantity (1,2,3 for velocity, 4,5,6 for energy flow)
c     srt - 0 : electrons, other - ions 
c     y0  - lower y boundary for processor     
      subroutine calcDens2D(px,num,srt,y0)
      implicit none
      include 'part.pf'
      real*8 s2,s4,s6,a,s21,s41,s61,s,pi,y0
      real*8 hx,hy,xm,ym,zm,h1,h2,h3,x1,y1,z1
      integer jm,i,k,l,j,m7,i2,l2,num,srt,lp,jmf,jmi,n1,n2,m5,m33
      real*8 px(omp,omp),ami,amf
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7
      common/b/xm,ym,zm
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/iff/xf,yf,zf,uf,vf,wf,amf
         
      if(deb.ge.2) write(37,*) 'qq3 ',hx,hy,y0
      
      do i = 1,omp
         do l = 1,omp
	    px(i,l) = 0d0
	 enddo
      enddo
        
      hy = ym*nproc/(omp-2)
      hx = xm/(omp-2)
      if(deb.ge.2) write(37,*) 'qq4 ',hx,hy,y0
      if(srt.eq.0) then
         jm = jmf
      else
         jm = jmi
      endif
      
      do j = 1,jm
         call partj(num,srt,j,a)
      
         if(srt.eq.0) then
            x1 = xf(j)
      	    y1 = yf(j)
      	    z1 = zf(j)
         else
            x1 = xi(j)
	    y1 = yi(j)
	    z1 = zi(j)
         endif
      
         s2 = x1/hx
         i=idint(s2+1.5d0)
         s2=s2+1.5d0-i
         s4=y1/hy
         l=idint(s4+1.5d0)
         s4=s4+1.5d0-l
         
      	 s21=1.d0-s2
         s41=1.d0-s4
         if((deb.ge.7).and.(num.eq.7)) write(37,*) 'qq4.5 ',j,y1,l,a,
     +       px(2,2)
         
      	 s=s21*s41*a
         px(i,l)=px(i,l)+s
         s=s21*s4*a
         px(i,l+1)=px(i,l+1)+s
         s=s2*s41*a
         px(i+1,l)=px(i+1,l)+s
         s=s2*s4*a
         px(i+1,l+1)=px(i+1,l+1)+s
      enddo       
      if(deb.ge.2) write(37,*) 'qq5 ',px(1,8)
c      stop
c      call PARAreduce(px)
      
      end subroutine calcDens2D
!------------------------------------------------------------------------------------------          
      subroutine outDens2D(y0)
      implicit none
      include 'part.pf'
      real*8 px(imp,lmpf,kmp),s2,s4,s6,a,s21,s41,s61,s,pi
      real*8 p2(omp,omp),hx,hy,xm,ym,zm,h1,h2,h3,ami,amf,y0
      integer jm,i,k,l,j,i2,l2,lp,jmf,jmi,n1,n2,m5,m7
      real*8 vx(omp,omp),vy(omp,omp),vz(omp,omp),ne(omp,omp)
      real*8 vxi(omp,omp),vyi(omp,omp),vzi(omp,omp)
      real*8 tfx(omp,omp),tfy(omp,omp),tfz(omp,omp)
      real*8 tfxi(omp,omp),tfyi(omp,omp),tfzi(omp,omp)
      real*8 ang(omp,omp)
      character*40 tr,tr1
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/iff/xf,yf,zf,uf,vf,wf,amf
      common/d/h1,h2,h3
      common/b/xm,ym,zm
      
      hy = ym*nproc/(omp-2)
      hx = xm/(omp-2)
      
      if(deb.ge.2) write(37,*) 'qq1,x,hy ',hx,hy  

      call calcDens2D(ne,0,0,y0)

      call calcDens2D(vx,1,0,y0)
      call calcDens2D(vy,2,0,y0)
      call calcDens2D(vz,3,0,y0)
      call calcDens2D(vxi,1,1,y0)
      call calcDens2D(vyi,2,1,y0)
      call calcDens2D(vzi,3,1,y0)

      call calcDens2D(tfx,4,0,y0)
      call calcDens2D(tfy,5,0,y0)
      call calcDens2D(tfz,6,0,y0)
      call calcDens2D(tfxi,4,1,y0)
      call calcDens2D(tfyi,5,1,y0)
      call calcDens2D(tfzi,6,1,y0)
      
      call calcDens2D(ang,7,0,y0)

      call PARAreduce12(ne,vx,vy,vz,vxi,vyi,vzi,
     +                  tfx,tfy,tfz,tfxi,tfyi,tfzi,ang)
      if(me.eq.0) then
      
         tr='dens' 
         write(tr1,'(i3.3)') m7 
         tr=tr(1:4)//tr1 
         tr=tr(1:7)//'.dat' 
         open(96,file=tr,form='formatted')
 320     format(2e12.3,14e15.5)        	 

  
         do i = 1,omp
            do l = 1,omp
	       write(96,320) i*hx,l*hy,ne(i,l),vx(i,l),vy(i,l),vz(i,l),
     +	                         vxi(i,l),vyi(i,l),vzi(i,l),  
     +	                         tfx(i,l),tfy(i,l),tfz(i,l),  
     +	                         tfxi(i,l),tfyi(i,l),tfzi(i,l),ang(i,l)
            enddo
         enddo
         close(96)
      endif
      if(deb.ge.2) write(37,*) 'qq2'  
c      call PARAFinal
c      stop
      
      end
!---------------------------------------------------------------------------
      subroutine densField(dee,num)
      implicit none
      include 'part.pf'
      real*8 ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp),dsqrt,
     *hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp),dee(imp,lmp,kmp)
      integer i,l,k,num
      common/g/hx,hy,hz
      common/h/ex,ey,ez
      
      do i = 1,imp
         do l = 1,lmp
      	    do k = 1,kmp
	       if(num.eq.0) then
                  dee(i,l,k) = ex(i,l,k)
c                  dee(i,l,k) = dsqrt(ex(i,l,k)**2+
c     +		                     ey(i,l,k)**2+ez(i,l,k)**2)
	       else    	  
                  dee(i,l,k) = dsqrt(hx(i,l,k)**2+
     +		                     hy(i,l,k)**2+hz(i,l,k)**2)
	       endif    	  
            enddo
      	 enddo
      enddo	             
      end
!--------------------------------------------------------------------------
c     px  - resulting matrix
c     num - code of quantity (1,2,3 for velocity, 4,5,6 for energy flow)
c     srt - 0 : electrons, other - ions 
c     y0  - lower y boundary for processor     
      subroutine distrFunc1D(srt,num)
      implicit none
      include 'part.pf'
      integer dfsize
      parameter(dfsize=10000)
      real*8 s2,s4,s6,a,s21,s41,s61,s,pi,y0,ms,hp
      real*8 hx,hy,xm,ym,zm,h1,h2,h3,x1,y1,z1
      integer jm,i,k,l,j,m7,i2,l2,num,srt,lp,jmf,jmi,n1,n2,m5,jmb
      real*8 px(dfsize),pimp(dfsize),ami,amf
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp),amb
      character*20 tr,tr1
      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
      common/b/xm,ym,zm
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/iff/xf,yf,zf,uf,vf,wf,amf
      
      print*,'in disf'

      do i = 1,dfsize
	 px(i) = 0d0
	 pimp(i) = 0d0
      enddo
      print*,'in disf 1'
        
      hy = 4d0/(dfsize-2)
      hp = 20.0d0/(dfsize-2)
      print*,'in disf2 ',jm

      do j = 1,jmf
         call partj(1,0,j,a)
         s4=(a+2d0)/hy
         l=idint(s4+1.5d0)
         s4=s4+1.5d0-l
         s41=1.d0-s4
         
      	 s=s41
         px(l)=px(l)+s
         s=s4
         px(l+1)=px(l+1)+s
         
         s4 = (uf(j)+10.0d0)/hp
         l=idint(s4+1.5d0)
         s4=s4+1.5d0-l
         
         s41=1.d0-s4
         
      	 s=s41
         pimp(l)=pimp(l)+s
         s=s4
         pimp(l+1)=pimp(l+1)+s

      enddo       
      
      do j = 1,jmb
         call partj(1,2,j,a)
         s4=(a+2d0)/hy
         l=idint(s4+1.5d0)
         s4=s4+1.5d0-l
         
         s41=1.d0-s4
      	 s=s41
         px(l)=px(l)+s
         s=s4
         px(l+1)=px(l+1)+s
         
         s4 = (ub(j)+10.0d0)/hp
         l=idint(s4+1.5d0)
         s4=s4+1.5d0-l
         
         s41=1.d0-s4
      	 s=s41
         pimp(l)=pimp(l)+s
         s=s4
         pimp(l+1)=pimp(l+1)+s
  344    format('bDF ',i10,3e15.5,i5,2e15.5)
      enddo       
      
      
      
      call PARAreduce(px)
c      stop
c      call PARAreduce(px)
      if((me.eq.0).and.(disf.eq.1)) then
             tr='disf' 
             write(tr1,'(i3.3)') m7 
             tr=tr(1:4)//tr1 
             tr=tr(1:7)//'.dat' 
             open(86,file=tr,form='formatted')
 320         format(i10,4e30.20)        	 
  
             do l = 2,dfsize-1
	        write(86,320) l,(l-dfsize/2)*hy,px(l), 
     +	        (l-dfsize/2)*hp,pimp(l)
             enddo
             print*,'out'
             close(86)
      endif
             print*,'out1'

      end
!--------------------------------------------------------------    
      subroutine distrFunc2D(srt,num)
      implicit none
      include 'part.pf'
      real*8 s2,s4,s6,a,s21,s41,s61,s,pi,y0,ms
      real*8 hx,hy,xm,ym,zm,h1,h2,h3,x1,y1,z1
      integer jm,i,k,l,j,m7,i2,l2,num,srt,lp,jmf,jmi,n1,n2,m5
      real*8 px(omp,omp),ami,amf
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      character*20 tr,tr1
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7
      common/b/xm,ym,zm
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/iff/xf,yf,zf,uf,vf,wf,amf
         
      
      do i = 1,omp
         do l = 1,omp
	    px(i,l) = 0d0
	 enddo
      enddo
        
      hy = 4d0/(omp-2)
      hx = xm/(omp-2)
      if(deb.ge.2) write(37,*) 'qq4 ',hx,hy,y0
      if(srt.eq.0) then
         jm = jmf
	 ms = -amf
      else
         jm = jmi
	 ms = ami
      endif
      
      do j = 1,jm
         call partj(num,srt,j,a)
      
         if(srt.eq.0) then
            x1 = xf(j)
      	    y1 = yf(j)
      	    z1 = zf(j)
         else
            x1 = xi(j)
	    y1 = yi(j)
	    z1 = zi(j)
         endif
      
         s2 = x1/hx
         i=idint(s2+1.5d0)
         s2=s2+1.5d0-i
         s4=(a+2d0)/hy
         l=idint(s4+1.5d0)
         s4=s4+1.5d0-l
         
      	 s21=1.d0-s2
         s41=1.d0-s4
         
      	 s=s21*s41*ms
         px(i,l)=px(i,l)+s
         s=s21*s4*ms
         px(i,l+1)=px(i,l+1)+s
         s=s2*s41*ms
         px(i+1,l)=px(i+1,l)+s
         s=s2*s4*ms
         px(i+1,l+1)=px(i+1,l+1)+s
      enddo       
      call PARAreduce(px)
      if(deb.ge.2) write(37,*) 'qq5 ',px(1,8)
c      call PARAreduce(px)
      if((me.eq.0).and.(disf.eq.1)) then
             tr='disf' 
             write(tr1,'(i3.3)') m7 
             tr=tr(1:4)//tr1 
             tr=tr(1:7)//'.dat' 
             open(96,file=tr,form='formatted')
 320         format(2e12.3,e15.5)        	 
  
             do i = 2,omp-1
                do l = 2,omp-1
	           write(96,320) (i-1)*hx,(l-omp/2)*hy,px(i,l)
                enddo
             enddo
      endif

      end
!-------------------------------------------------------------------
      subroutine getBeamEnergy(enb,dnb)
      implicit none
      include 'part.pf'
      integer bmst,bmvli,lp,jmf,jmi,n1,n2,m5,m7,j
      real*8 enb,dnb,pi,t,h1,h2,h3,vec(1000)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp),amf
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7
      common/iff/xf,yf,zf,uf,vf,wf,amf
      common/d/h1,h2,h3
      
      enb = 0d0
      dnb = 0d0
      if(beamf.eq.0) return
      
      bmst=(imp-2)*(lmp-2)*(kmp-2)*lp
      
      do j=jmf/2+1,jmf
         enb=enb+((uf(j))**2+(vf(j))**2+(wf(j))**2)
         t = dsqrt(uf(j)**2+vf(j)**2+wf(j)**2)
         if(((0.8d0*rimp).lt.t).and.((1.2d0*rimp).gt.t)) then
	    dnb = dnb+1d0
	 endif
      enddo
      enb=amf
      vec(1) = enb
      vec(2) = dnb
      call PARAreduceLMP2(vec)
      enb = vec(1)
      dnb = vec(2)
      
      end
!---------------------------------------------------------------------
      subroutine timeBegin(num)
      implicit none
      external amicro 
      real*8 amicro
      integer num
      real*8 tm(3),pTime,fTime,fullTime,cpTime,cfTime,cfullTime
      common/timeval/tm,pTime,fTime,fullTime,cpTime,cfTime,cfullTime
      
      tm(num) = amicro()
      end
!---------------------------------------------------------------------
      subroutine timeEnd(num)
      implicit none
      external amicro 
      real*8 amicro
      integer num
      real*8 tm(3),pTime,fTime,fullTime,cpTime,cfTime,cfullTime
      common/timeval/tm,pTime,fTime,fullTime,cpTime,cfTime,cfullTime
      
      tm(num) = - (tm(num) - amicro())
      
      if(num.eq.1) then
         cpTime    = cpTime   + tm(1)
      endif

      if(num.eq.2) then
         cfTime    = cfTime   + tm(2)
      endif

      if(num.eq.3) then
         cfullTime    = cfullTime   + tm(3)
      endif
      end
      subroutine cumultime(tm1,time)
      implicit none
      external amicro
      real*8 amicro,tm1,tm2,time
      
      tm2 = amicro()
      time = time+tm2-tm1
      
      end subroutine cumultime
!--------------------------------------------------------------------
      subroutine printTime(nt)
      implicit none
      include 'part.pf'
      integer num
      integer bmst,bmvli,lp,jmf,jmi,n1,n2,m5,m7,nt
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7
      real*8 tm(3),pTime,fTime,fullTime,cpTime,cfTime,cfullTime,pi
      ! particle send-receive time, field send-receive time, 
      ! collective communications time FOR CURRENT TIMESTEP
      real*8 psrt,fsrt,colt,gatht
      integer lc,rc
      real*8  rlc,rrc
      common/ptranstat/lc,rc
      common/mpitime/psrt,fsrt,colt
      common/timeval/tm,pTime,fTime,fullTime,cpTime,cfTime,cfullTime

      character*40 st1,st,fn
      external  writecpuinfo
    
      call  writecpuinfo(me)
      
      rlc = lc    ! tranforming number of flown-away particles into real
      rrc = rc

      if(nt.eq.0) then
         pTime    = 0d0
	 fTime    = 0d0
	 fullTime = 0d0
	 
         cpTime    = 0d0
	 cfTime    = 0d0
	 cfullTime = 0d0
            
         psrt  = 0d0
         fsrt  = 0d0
         colt  = 0d0
	 gatht = 0d0
         
         write(st,'(i5.5)') me

         st1='time'//st
         fn=st1(1:9)//'.dat'

	 open(93,file=fn,form='formatted')
      else	 

         pTime    = pTime    + cpTime
         fTime    = fTime    + cfTime
         fullTime = fullTime + cfullTime
      
  950    format(i10,9e15.6)
         if(me.eq.0) then
            write(93,950) nt,cpTime,cfTime,cfullTime,fsrt,colt,psrt,
     +                    gatht,rlc,rrc
         endif
      
         if((nt.eq.totSteps).and.(me.eq.0)) then
            write(93,950) nt,pTime,fTime,fullTime 
         endif
         cpTime    = 0d0
	 cfTime    = 0d0
	 cfullTime = 0d0

      endif
         
      end
!-----------------------------------------------------------------------------
      subroutine getFreq(nu,nt)
      implicit none
      include 'part.pf'
      real*8 nueff,nucl
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp),amf
      real*8 gam(omp,omp),numb(omp,omp),taue,vec(lmp2)
      real*8 hix,hiy,imaxx,imaxy,iminx,iminy,s,s2,s21,s41
      real*8 uf1(jmp),vf1(jmp),wf1(jmp),pi,avx,avy,a,s4,nu,temp(1000)
      real*8 datan,dacos,alpha,pi2
      integer j,lp,jmf,jmi,n1,n2,m5,m7,i,k,l,npi2,npi2c,nt
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7
      common/iff/xf,yf,zf,uf,vf,wf,amf
      common/imp1/uf1,vf1,wf1,hix,hiy,iminx,iminy
      common/turnprev/npi2

      pi2   = 2d0*datan(1.0d0)
      npi2c = 0
      do j = 1,jmf
         alpha = dacos(uf(j)/dsqrt(uf(j)**2+vf(j)**2+wf(j)**2))
         if(alpha.gt.pi2) npi2c = npi2c+1
      enddo
      if(deb.ge.4) then
         write(37,*) 'npi2c,nt,npi2 ',npi2c,nt,npi2
      endif      
      
      temp(1) = npi2c
      call PARAreduceLMP2(temp)
      npi2c = temp(1)
      if(deb.ge.4) then
         write(37,*) 'PAara npi2c ',npi2c
      endif      
      
      if(nt.gt.0) then
         nu = (npi2c - npi2)
         nu = nu/jmf/fproc
         npi2 = npi2c
      else
         nu = npi2c
         nu = nu/jmf/fproc
      endif

      if(deb.ge.4) then
         write(37,*) 'npi2,nu,jmf ',npi2,nu,jmf
      endif      
      
      npi2 = npi2c

      end
!-------------------------------------------------------------------------
      subroutine interpolate(y0,sec)
      implicit none
      include 'part.pf'
      real*8 sec(omp,omp,6),x,z,sx,sy,sz,bx,by,bz,hx0,hz0,y,y0
      real*8 s1,s2,s3,s4,s5,s6,s11,s21,s31,s41,s51,s61
      real*8 ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp),
     *hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp),h1,h2,h3
      common/g/hx,hy,hz
      common/h/ex,ey,ez
      common/d/h1,h2,h3

      integer i,k,i1,i2,l1,l2,k1,k2,p,q,l
      
      hx0=lx0/omp
      hz0=lz0/omp
      
      y  = ly0/2d0
      y0 = meh*(ly0/nproc)
      
      do p = 1,omp
         do q = 1,omp
            x = (p-0.5d0)*hx0
            z = (q-0.5d0)*hz0
        
c***************************************
            s2=x/h1
            i=idint(s2+1.d0)
            i1=idint(s2+1.5d0)
            s1=i-s2
            s2=i1-0.5d0-s2
            
            s4=(y-y0)/h2
            l=idint(s4+1.d0)
            l1=idint(s4+1.5d0)
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
      
            sx=(s1*(s4*(s6*ex(i,l1,k1)+s61*ex(i,l1,k1+1))+
     +         s41*(s6*ex(i,l1+1,k1)+s61*ex(i,l1+1,k1+1)))+
     +         s11*(s4*(s6*ex(i+1,l1,k1)+s61*ex(i+1,l1,k1+1))+
     +         s41*(s6*ex(i+1,l1+1,k1)+s61*ex(i+1,l1+1,k1+1))))

            sy=(s2*(s3*(s6*ey(i1,l,k1)+s61*ey(i1,l,k1+1))+
     +         s31*(s6*ey(i1,l+1,k1)+s61*ey(i1,l+1,k1+1)))+
     +         s21*(s3*(s6*ey(i1+1,l,k1)+s61*ey(i1+1,l,k1+1))+
     +         s31*(s6*ey(i1+1,l+1,k1)+s61*ey(i1+1,l+1,k1+1))))
     
            sz=(s2*(s4*(s5*ez(i1,l1,k)+s51*ez(i1,l1,k+1))+
     +         s41*(s5*ez(i1,l1+1,k)+s51*ez(i1,l1+1,k+1)))+
     +         s21*(s4*(s5*ez(i1+1,l1,k)+s51*ez(i1+1,l1,k+1))+
     +         s41*(s5*ez(i1+1,l1+1,k)+s51*ez(i1+1,l1+1,k+1))))
      
            bx=(s2*(s3*(s5*hx(i1,l,k)+s51*hx(i1,l,k+1))+
     +         s31*(s5*hx(i1,l+1,k)+s51*hx(i1,l+1,k+1)))+
     +         s21*(s3*(s5*hx(i1+1,l,k)+s51*hx(i1+1,l,k+1))+
     +         s31*(s5*hx(i1+1,l+1,k)+s51*hx(i1+1,l+1,k+1))))
      
            by=(s1*(s4*(s5*hy(i,l1,k)+s51*hy(i,l1,k+1))+
     +         s41*(s5*hy(i,l1+1,k)+s51*hy(i,l1+1,k+1)))+
     +         s11*(s4*(s5*hy(i+1,l1,k)+s51*hy(i+1,l1,k+1))+
     +         s41*(s5*hy(i+1,l1+1,k)+s51*hy(i+1,l1+1,k+1))))
      
            bz=(s1*(s3*(s6*hz(i,l,k1)+s61*hz(i,l,k1+1))+
     +         s31*(s6*hz(i,l+1,k1)+s61*hz(i,l+1,k1+1)))+
     +         s11*(s3*(s6*hz(i+1,l,k1)+s61*hz(i+1,l,k1+1))+
     +         s31*(s6*hz(i+1,l+1,k1)+s61*hz(i+1,l+1,k1+1))))
c***************************************            
            sec(p,q,1) = sx            
            sec(p,q,2) = sy            
            sec(p,q,3) = sz
            sec(p,q,4) = bx            
            sec(p,q,5) = by            
            sec(p,q,6) = bz            
         enddo
      enddo
      
      call PARAreduce6(sec)

      end
!---------------------------------------------------------------------------------
      subroutine outdens1d(y0)
      implicit none
      include 'part.pf'
      real*8 s2,s4,s6,a,s21,s41,s61,s,pi,y0
      real*8 hx,hy,xm,ym,zm,h1,h2,h3,x1,y1,z1,amf,e,amb
      real*8 sec(omp,omp,6)
      integer jm,i,k,l,j,m7,i2,l2,num,srt,lp,jmf,jmi,n1,n2,m5,jmb
      character*20 tr,tr1
      real*8 be(1000),ee(1000),fe(1000),bn(1000),fn(1000)
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/b/xm,ym,zm
      common/iff/xf,yf,zf,uf,vf,wf,amf
        
      if(out1dd.ne.1) return
      
      if(deb.ge.2) write(37,*) 'qq3 ',hx,hy,y0
      
      do i = 1,omp
         be(i) = 0d0
         fe(i) = 0d0
         ee(i) = 0d0
         bn(i) = 0d0
         fn(i) = 0d0
      enddo
        
      hx = xm/(omp-2)
      if(deb.ge.2) write(37,*) '1d hx ',hx
      
      do j = 1,jmf
            x1 = xf(j)
      	    y1 = yf(j)
      	    z1 = zf(j)
      
         s2 = x1/hx
         i=idint(s2+1.5d0)
         s2=s2+1.5d0-i
         
      	 s21=1.d0-s2

      	 s=s21*amf
         fn(i)=fn(i)+s
         s=s2*amf
         fn(i+1)=fn(i+1)+s

         e=amf*(dsqrt(1.d0+beta0**2*(uf(j)**2+vf(j)**2+wf(j)**2))-1.d0)
      	 s = s21*e
         fe(i)=fe(i)+s
         s = s2*e
         fe(i+1)=fe(i+1)+s


      enddo       
      if(deb.ge.2) write(37,*) 'fn ',j,s2,s21,fn(i),fn(i+1)
      call PARAreduceLMP2(fn)
      call PARAreduceLMP2(fe)

      do j = 1,jmb
            x1 = xb(j)
      	    y1 = yb(j)
      	    z1 = zb(j)
      
         s2 = x1/hx
         i=idint(s2+1.5d0)
         s2=s2+1.5d0-i
         
      	 s21=1.d0-s2

      	 s=s21*amb
         bn(i)=bn(i)+s
         s=s2*amb
         bn(i+1)=bn(i+1)+s

         e=amb*(dsqrt(1.d0+beta0**2*(ub(j)**2+vb(j)**2+wb(j)**2))-1.d0)
      	 s = s21*e
         be(i)=be(i)+s
         s = s2*e
         be(i+1)=be(i+1)+s


      enddo       
      if(deb.ge.2) write(37,*) 'fn ',j,s2,s21,fn(i),fn(i+1)
      call PARAreduceLMP2(bn)
      call PARAreduceLMP2(be)

      call interpolate(y0,sec)
      
         tr='pn1d' 
         write(tr1,'(i3.3)') m7 
         tr=tr(1:4)//tr1 
         tr=tr(1:7)//'.dat' 
         open(96,file=tr,form='formatted')
 320     format(e12.5,5e15.5)        	 

      do i = 2,omp-1
         e = 0
         do k = 2,omp-1
            e = e + sec(i,k,1)**2+sec(i,k,2)**2+sec(i,k,3)**2
         enddo
         e = e/(omp-2)/8d0
         write(96,320) hx*(i-1.5d0),fn(i),fe(i),bn(i),be(i),e
      enddo
      
      end
!------------------------------------------------------------------------------------  
c     px  - resulting matrix
c     num - code of quantity (1,2,3 for velocity, 4,5,6 for energy flow)
c     srt - 0 : electrons, other - ions 
c     y0  - lower y boundary for processor     
      subroutine phasePlane
      implicit none
      include 'part.pf'
      real*8 s2,s4,s6,a,s21,s41,s61,s,pi,y0
      real*8 hx,hy,xm,ym,zm,h1,h2,h3,x1,y1,z1
      integer jm,i,k,l,j,m7,i2,l2,num,srt,lp,jmf,jmi,n1,n2,m5,nt,jmb
      real*8 px(omp,omp),ami,amf,tau,c1,c2,c3,amb
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      common/f/tau,c1,c2,c3
      common/b/xm,ym,zm
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
      common/ib/xb,yb,zb,ub,vb,wb,amb
      real*8 ps,v,v0       
      
      character*20 tr,tr1
         
      hy = 2d0/(omp-2)
      hx = lx0/(omp-2) 
      if(deb.ge.2) write(37,*) 'qq3 phase ',hx,hy,y0
c      stop
      
      if(m7.eq.1) then
         tr='phas' 
         write(tr1,'(i3.3)') me
         tr=tr(1:4)//tr1 
         tr=tr(1:7)//'.dat' 
         open(38,file=tr,form='formatted')
 320     format(2i10,4e15.5)        	 

      do j = 1,jmb
         write(38,320) m7,j,xb(j),ub(j),vb(j),wb(j)
      enddo
       close(38)
      endif
      
      end
c_____________________________________________________________________
      subroutine printv
      implicit none
      include 'part.pf'
      integer j,lp,jmf,jmi,n1,n2,m5,m7,jmb
      real*8 pi,ux,uy,uz,gb0,amb,ami,amf
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      character*40 tr,tr1
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/iff/xf,yf,zf,uf,vf,wf,amf
      
      if(beam_write.eq.0) then
         return
      endif
      
      tr = 'beam'    
            write(tr1,'(i6.6,i3.3)') m7,me 
            tr=tr(1:4)//tr1 
            tr=tr(1:13)//'.dat' 
            open(38,file=tr,form='formatted')
        
           do j=1,jmb
              gb0=1d0/dsqrt(1d0+ub(j)**2+vb(j)**2+wb(j)**2)
              ux=ub(j)*gb0
              uy=vb(j)*gb0
              uz=wb(j)*gb0
             
              write(38,321) xb(j),yb(j),zb(j),ux,uy,uz
c              write(38,322) xb(j),ux
c              write(38,321)ub(j),vb(j),wb(j)
           enddo 
 321     format(6e30.20)
         close(38)


c      tr = 'efon'    
c            write(tr1,'(i6.6,i3.3)') m7,me 
c            tr=tr(1:4)//tr1 
c            tr=tr(1:13)//'.dat' 
c            open(38,file=tr,form='formatted')
c        
c           do j=1,jmf
c              write(38,321)uf(j),vf(j),wf(j)
c           enddo 
c         close(38)
c         
c      tr = 'ions'    
c            write(tr1,'(i6.6,i3.3)') m7,me 
c            tr=tr(1:4)//tr1 
c            tr=tr(1:13)//'.dat' 
c            open(38,file=tr,form='formatted')
c        
c           do j=1,jmi
c              write(38,321)ui(j),vi(2*j),wi(2*j)
c           enddo 
c         close(38)
               

c         tr = 'efon'    
c            write(tr1,'(i6.6,i3.3)') m7,me 
c            tr=tr(1:4)//tr1 
c            tr=tr(1:13)//'.dat' 
c            open(38,file=tr,form='formatted')
c        
c           do j=1,jmi
c              gb0=1d0/dsqrt(1d0+ui(j)**2+vi(j)**2+wi(j)**2)
c              ux=ui(j)*gb0
c              uy=vi(j)*gb0
c              uz=wi(j)*gb0
c              write(38,321) xf(j),yf(j),zf(j),ux,uy,uz
c           enddo 
c         close(38) 

      end  subroutine printv
