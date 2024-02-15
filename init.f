      !num - number of initial set
      ! 0  - single electron with 0 impulse,
      !      zero magnetic field and constant electric field
      subroutine init(num)
!   USE DFLIB
      implicit real*8(a-h,o-z)
      include 'part.pf'
    include 'mass.par'
c      parameter(imp=42,lmp=22,kmp=3,jmp=224000)
       integer num ! number of initial set
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
      integer i,l,k

      if(num.eq.1) then
         do i = 1,imp
            do l = 1,lmp
               do k = 1,kmp
                  hx(i,l,k) = 0d0
                  hy(i,l,k) = 0d0
                  hz(i,l,k) = 0d0
                  jx(i,l,k) = 0d0
                  jy(i,l,k) = 0d0
                  jz(i,l,k) = 0d0
                  qx(i,l,k) = 0d0
                  qy(i,l,k) = 0d0
                  qz(i,l,k) = 0d0
                  ex(i,l,k) = 1d0
                  ey(i,l,k) = 0d0
                  ez(i,l,k) = 0d0
               enddo
            enddo
         enddo
      endif
      jmf = 1
      xf(1) = h1*0.5d0
      yf(1) = h2*0.5d0
      zf(1) = h3*0.5d0
      uf(1) = 0d0
      vf(1) = 0d0
      wf(1) = 0d0
      jmi = 0
      jmb = 0

      end subroutine init





