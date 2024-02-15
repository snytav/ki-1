





      subroutine allocControlAttribute
      implicit none
      include 'part.pf'
      include 'mass.par'

c***********************************************************

      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      real*8 ami,amb,amf
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/iff/xf,yf,zf,uf,vf,wf,amf

c**********************************************************
      real*8 a,t
      integer j,num,pos,nt
      real*8, pointer :: atrib(:)    ! (3*jmp*attr)
      common/ctrl/atrib

      if(ctrl_attr.eq.0) then
         return
      endif

      allocate(atrib(3*jmp*attr))

      end subroutine allocControlAttribute

      subroutine writeControlAttribute(j,a,num,t,nt)
      implicit none
      include 'part.pf'
      include 'mass.par'

c***********************************************************

      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      real*8 ami,amb,amf
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/iff/xf,yf,zf,uf,vf,wf,amf

c**********************************************************
      real*8 a,t
      integer j,num,pos,nt
      real*8, pointer :: atrib(:)    ! (3*jmp*attr)
      common/ctrl/atrib

      if(ctrl_attr.eq.0) then
         return
      endif


      if(dabs(ami-a).lt.1.0d-15) then
         pos = num + (j-1)*attr
      endif

      if(dabs(amf-a).lt.1.0d-15) then
         pos = num + (j-1)*attr + jmp*attr
      endif

      if(dabs(amb-a).lt.1.0d-15) then
         pos = num + (j-1)*attr + 2*jmp*attr
      endif

      atrib(pos) = t

      if((j.eq.18).and.(num.eq.44)) then
 445     format('writeControlAttribute44 ',i5,i10,e15.5,i10,i10,e25.17)

         write(37,445) nt,j,a,num,pos,t
      endif

      end subroutine writeControlAttribute
       

      subroutine saveControlAttributesToFile(nt)
      implicit none
      include 'part.pf'
      include 'mass.par'

c***********************************************************

      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      integer lp,jmf,jmi,n1,n2,m5,m7,jmb
      real*8 ami,amb,amf,pi
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/iff/xf,yf,zf,uf,vf,wf,amf
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb

c**********************************************************

      real*8 t,f(3)
      integer j,num,nt,l,i,fi(3)
      real*8, pointer ::  atrib(:) ! (3*jmp*attr)
      common/ctrl/atrib

      character*40 tr,tr1
      
      if(ctrl_attr.eq.0) then
         return
      endif

      tr='ctrl'
      write(tr1,'(i5.5)') nt
      tr=tr(1:4)//tr1

      f(1) = ami
      f(2) = amf
      f(3) = amb

      fi(1) = jmi
      fi(2) = jmf
      fi(3) = jmb

 137  format('jmp,attr jmi,jmf,jmb ',3i10,2i3,3i10)
      write(37,137) jmp,attr,3*jmp*attr,meh,me,jmi,jmf,jmb
      open(38,file=tr,form='unformatted')
      write(38) (f(i),i=1,3)
      write(38) (fi(i),i=1,3)
      write(38) (atrib(i),i=1,3*jmp*attr)
      close(38)

      end subroutine saveControlAttributesToFile
