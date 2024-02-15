	!Режимы записи на диск (параметр ml):
	! 0 - начальные данные не считываются из файла; конечные данные не выводятся в файл
	! 1 - начальные данные не считываются из файла; конечные данные    выводятся в файл
	! 2 - начальные данные    считываются из файла; конечные данные не выводятся в файл
	! 3 - начальные данные    считываются из файла; конечные данные    выводятся в файл
	!Открываются файлы для вывода
	subroutine open_f
      implicit real*8(a-h,o-z)
      include 'part.pf'
	if((outterm.eq.1).and.(me.eq.0)) then
		open(21,file='energy.dat',form='formatted')
		open(22,file='impuls.dat',form='formatted')
		open(49,file='thermf.dat',form='formatted')
		open(50,file='thermi.dat',form='formatted')
	endif
	!------------------------------------------------
  306 format(i4)
	end subroutine open_f

!-----------------------------------------------------------------------	
	!Ввод из файла информации для продолжения счета
	subroutine vvod
      implicit real*8(a-h,o-z)
      include 'laser1.par'
      include 'laser2.par'
      include 'laser3.par'
      include 'laser4.par'
c      parameter(imp=42,lmp=22,kmp=3,jmp=224000)
      real*8 ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp),f(30)
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
!      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp),
     *hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp)
      real*8 mee,nee,mii,ne0
      real*8 ni,ne,Te,hx0,lx,ly,lz,lx0,ly0,lz0,hx00,ni0
	real*8 dT,dn,dE,du,kk,gamma
      common/a/tx,nt,ml,sst,ns
      common/b/lx,ly,lz
      common/c/im,lm,km
      common/d/h1,h2,h3
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7
      common/f/tau,c1,c2,c3
      common/g/hx,hy,hz
      common/h/ex,ey,ez
      common/ii/xi,yi,zi,ui,vi,wi,ami
!      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/iff/xf,yf,zf,uf,vf,wf,amf
      common/c1/alpha,t0,tD
      common/c2/ni,ne,Te,hx0
	common/cc/ncc
	common/et/ei0,ef0,eh0,ee0,ew0
	common/dd/dT,dn,dE,du,kk,gam
	common/dispold/disp_old
	common/den1/j_pp,i_pp
      rewind 10

      read(10) (f(i),i=1,30)

      read(10) (((hx(i,l,k),i=1,imp),l=1,lmp),k=1,kmp)
      read(10) (((hy(i,l,k),i=1,imp),l=1,lmp),k=1,kmp)
      read(10) (((hz(i,l,k),i=1,imp),l=1,lmp),k=1,kmp)
      read(10) (((ex(i,l,k),i=1,imp),l=1,lmp),k=1,kmp)
      read(10) (((ey(i,l,k),i=1,imp),l=1,lmp),k=1,kmp)
      read(10) (((ez(i,l,k),i=1,imp),l=1,lmp),k=1,kmp)

      read(10) (xi(j),j=1,jmp)
      read(10) (yi(j),j=1,jmp)
      read(10) (zi(j),j=1,jmp)
      read(10) (ui(j),j=1,jmp)
      read(10) (vi(j),j=1,jmp)
      read(10) (wi(j),j=1,jmp)

      read(10) (xf(j),j=1,jmp)
      read(10) (yf(j),j=1,jmp)
      read(10) (zf(j),j=1,jmp)
      read(10) (uf(j),j=1,jmp)
      read(10) (vf(j),j=1,jmp)
      read(10) (wf(j),j=1,jmp)

      tx=f(1)
      nb=f(2)+0.1
      nt=f(3)+0.1
      im=f(4)+0.1
      lm=f(5)+0.1
      km=f(6)+0.1
      lp=f(7)+0.1
	ei0=f(8)
	ef0=f(9)
	eh0=f(10)
	ee0=f(11)
	ew0=f(12)
      lx=f(13)
      ly=f(14)
      lz=f(15)
      h1=f(16)
      h2=f(17)
      h3=f(18)
	ns=f(19)+0.1
	sst=f(20)
!	u0=f(21)
	jm=f(22)+0.1
	ami=f(23)
	amf=f(24)
	alpha=f(25)
	beta0=f(26)
	dT=f(27)
	kk=f(28)
	disp_old=f(29)	
	i_pp=f(30)+0.1	        
	end subroutine vvod

!-----------------------------------------------------------------------------
	!Вывод в файл при завершении программы для возможного продолжения счета
	subroutine vivod
      implicit real*8(a-h,o-z)
      include 'laser1.par'
      include 'laser2.par'
      include 'laser3.par'
      include 'laser4.par'
c      parameter(imp=42,lmp=22,kmp=3,jmp=224000)
      real*8 ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp),f(30)
      real*8 hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp)
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
!      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      real*8 mee,nee,mii,ne0
      real*8 ni,ne,Te,hx0,lx,ly,lz,lx0,ly0,lz0,hx00,ni0
	real*8 dT,dn,dE,du,kk,gamma
      common/a/tx,nt,ml,sst,ns
      common/b/lx,ly,lz
      common/c/im,lm,km
      common/d/h1,h2,h3
      common/g/hx,hy,hz
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7
      common/f/tau,c1,c2,c3
      common/h/ex,ey,ez
      common/ii/xi,yi,zi,ui,vi,wi,ami
!      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/iff/xf,yf,zf,uf,vf,wf,amf
      common/c1/alpha,t0,tD
      common/c2/ni,ne,Te,hx0
	common/cc/ncc
	common/et/ei0,ef0,eh0,ee0,ew0
	common/dd/dT,dn,dE,du,kk,gam
	common/dispold/disp_old
      common/den1/j_pp,i_pp
      rewind 10
      f(1)=tx
      f(3)=nt
      f(4)=im
      f(5)=lm
      f(6)=km
      f(7)=lp
	f(8)=ei0
	f(9)=ef0
	f(10)=eh0
	f(11)=ee0
	f(12)=ew0
      f(13)=lx
      f(14)=ly
      f(15)=lz
      f(16)=h1
      f(17)=h2
      f(18)=h3
	f(19)=ns
	f(20)=sst
!	f(21)=u0
	f(22)=jm
	f(23)=ami
	f(24)=amf
	f(25)=alpha
	f(26)=beta0
	f(27)=dT
	f(28)=kk
	f(29)=disp_old
	f(30)=i_pp
      write(11) (f(i),i=1,30)

 9014 format(2x,'maccЁb  f  §aЇЁcah')
      write(11) (((hx(i,l,k),i=1,imp),l=1,lmp),k=1,kmp)
      write(11) (((hy(i,l,k),i=1,imp),l=1,lmp),k=1,kmp)
      write(11) (((hz(i,l,k),i=1,imp),l=1,lmp),k=1,kmp)
      write(11) (((ex(i,l,k),i=1,imp),l=1,lmp),k=1,kmp)
      write(11) (((ey(i,l,k),i=1,imp),l=1,lmp),k=1,kmp)
      write(11) (((ez(i,l,k),i=1,imp),l=1,lmp),k=1,kmp)

      write(11) (xi(j),j=1,jmp)
      write(11) (yi(j),j=1,jmp)
      write(11) (zi(j),j=1,jmp)
      write(11) (ui(j),j=1,jmp)
      write(11) (vi(j),j=1,jmp)
      write(11) (wi(j),j=1,jmp)

      write(11) (xf(j),j=1,jmp)
      write(11) (yf(j),j=1,jmp)
      write(11) (zf(j),j=1,jmp)
      write(11) (uf(j),j=1,jmp)
      write(11) (vf(j),j=1,jmp)
      write(11) (wf(j),j=1,jmp)
	!-----------------------------------------
      end subroutine vivod
