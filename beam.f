      subroutine beamboundcheck(x,x1,sort,j,nt)
      implicit none
      include 'part.pf'
      real*8 x,x1
      character*2 sort
      integer bcc,fst,j,nt
      common/bcc1/bcc

      if((j.eq.1).and.(sort.eq.'bm')) then
         bcc=0
      endif

      if((sort.eq.'bm').and.(deb.ge.2)) then
          write(37,*) 'beam doubt',j,x,x1
      endif

      if((sort.eq.'bm').and.(x.lt.xmb).and.(x1.gt.xmb)) then
         write(37,*) 'beam doubt solved',j,x,x1

         bcc = bcc + 1

         if((sort.eq.'bm').and.(deb.ge.2)) then
            write(37,*) 'beam runs ',j,bcc,nt
         endif
      endif

      
          
      end 

      subroutine addbeam
      implicit none
      include 'part.pf'
      real*8 x,x1
      integer bcc,j,sort,nt,jml
      real*8 ux(jmp/imp),uy(jmp/imp),uz(jmp/imp)
      real*8 xl(jmp/imp),yl(jmp/imp),zl(jmp/imp),ul(jmp/imp),
     +       vl(jmp/imp),wl(jmp/imp)
      real*8 xbl(jmp/imp),ybl(jmp/imp),zbl(jmp/imp),ubl(jmp/imp),
     +       vbl(jmp/imp),wbl(jmp/imp)
      real*8 xfl(jmp/imp),yfl(jmp/imp),zfl(jmp/imp),ufl(jmp/imp),
     +       vfl(jmp/imp),wfl(jmp/imp)
      real*8 xil(jmp/imp),yil(jmp/imp),zil(jmp/imp),uil(jmp/imp),
     +       vil(jmp/imp),wil(jmp/imp)

      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      real*8 xb(jmp),yb(jmp),zb(jmp),ub(jmp),vb(jmp),wb(jmp)
      real*8 xf(jmp),yf(jmp),zf(jmp),uf(jmp),vf(jmp),wf(jmp)
      real*8 ak(jmp),trg,g05,vb0,vf01,vf02,pinv1,pinv2,gb0
      real*8 mee,nee,mii,ne0,kl,kl_m,tau_e
      real*8 ni,ne,Te,hx0,lx,ly,lz,hx00,ni0
      real*8 dT,dn,dE,du,kk,const1,gam
      real*8 dln,eps,dx1,dx2,eps1,ak0
      real*8 ksi,termx,h1,h2,h3
      real*8 lambda,lambda_m
      real*8 tf0,Tbb,vy,vz,z,y,wrapg05cae,wrapg05dde,pi,vf0,pinv
      real*8 tx,sst,c1,c2,c3,tau,ami,amb,amf,p_tau
  
      integer lp,i,ml,ns,im,lm,km,jmf,jmi,jmb,n1,n2,m5,m7,jmil,jmfl
      integer ndde,ncae
      common/rng/ndde,ncae

      common/a/tx,nt,ml,sst,ns
      common/b/lx,ly,lz
      common/c/im,lm,km
      common/d/h1,h2,h3
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
      common/f/tau,c1,c2,c3
      common/p/p_tau
      common/ii/xi,yi,zi,ui,vi,wi,ami
      common/ib/xb,yb,zb,ub,vb,wb,amb
      common/iff/xf,yf,zf,uf,vf,wf,amf
      common/bcc1/bcc
      
      

      write(37,*) 'beamstatb ',bcc,ndde,ncae
      

      if(me.ne.0) return
 
      do j = 1,bcc
  	    z=(lz-zmb)*0.5d0 + zmb*wrapg05cae (lz)
   	    y=(ly-ymb)*0.5d0 + ymb* wrapg05cae (ly)
            x=xmb* wrapg05cae (lx)
 70         format('basiccorrds ',3e30.20)
            if(deb.ge.2) then
               write(37,70) x,y,z
            endif
            xil(j)=x
            yil(j)=y
            zil(j)=z
            uil(j)=0.d0	!g05dde(0.d0,dsqrt(1.d0*alpha)) 
            vil(j)=0.d0	!g05dde(0.d0,dsqrt(1.d0*alpha)) 
            wil(j)=0.d0	!g05dde(0.d0,dsqrt(1.d0*alpha)) 
      enddo

        do j=1,bcc
           xbl(j)=xil(j)
           ybl(j)=yil(j)
           zbl(j)=zil(j)
 71        format('beamcorrds ',3e30.20)
            if(deb.ge.2) then
               write(37,71) xbl(j),ybl(j),zbl(j)
            endif
           vb0  =wrapg05dde(0.d0,Tb*rimp)
           ux(j)=rimp+vb0
           uy(j)=wrapg05dde(0.d0,Tb*rimp) 
           uz(j)=wrapg05dde(0.d0,Tb*rimp) 
        enddo
        
        do j=1,bcc
           vb0=sqrt(1.d0-ux(j)**2-uy(j)**2-uz(j)**2)
           ubl(j)=ux(j)/vb0
           vbl(j)=uy(j)/vb0
           wbl(j)=uz(j)/vb0
        enddo
      
        j=0
        do i=1,bcc
           j=j+1
    	   
           xfl(2*j-1)=xil(j)
           yfl(2*j-1)=yil(j)
           zfl(2*j-1)=zil(j)

           xfl(2*j)=xil(j)
           yfl(2*j)=yil(j)
           zfl(2*j)=zil(j)

           vy=wrapg05dde(0.d0,tey0)    
           vz=wrapg05dde(0.d0,tez0)        

c          INVERSE CURRENT
           
           termx = wrapg05dde(0.d0,tex0)
           gb0=1d0/dsqrt(1d0+ubl(j)**2+vbl(j)**2+wbl(j)**2)
           vb0=ubl(j)*gb0
           
           vf01=-rbd*vb0+termx  
           vf02=-rbd*vb0-termx
           if(deb.ge.2) then
              write(37,*) 'vf01,rbd,vb0 ',vf01,rbd,vb0
           endif
       	   
           pinv1= vf01/dsqrt((1d0-vf01**2-vy**2-vz**2)) 
           pinv2= vf02/dsqrt((1d0-vf02**2-vy**2-vz**2))
            
           vfl(2*j-1)= vy/dsqrt((1d0-vf01**2-vy**2-vz**2))
           vfl(2*j)  =-vy/dsqrt((1d0-vf02**2-vy**2-vz**2))
           wfl(2*j-1)= vz/dsqrt((1d0-vf01**2-vy**2-vz**2))
           wfl(2*j)  =-vz/dsqrt((1d0-vf02**2-vy**2-vz**2))
          
c my correct end
         
           ufl(2*j-1)= pinv1
           ufl(2*j)  = pinv2
        
        
 897       format('electron ',i10,7e15.5) 
           if(deb.ge.2) then
              write(37,897) j,xfl(2*j),yfl(2*j),zfl(2*j),ufl(2*j),
     +                       ufl(2*j-1),vfl(2*j),wfl(2*j)
 344          format('electron ',i10,6e15.5)           
 345          format('inverse current: b,f,-f',i10,5e15.5)
              write(37,345) j,vb0,vf01,vf02,
     +                      (vf01-termx)/rbd,(vb0+(vf01-termx)/rbd)
           endif
         
      enddo
      jmil = j
      jmfl = 2*j

      do j=1,bcc
         xb(j+jmb) = xbl(j)
         yb(j+jmb) = ybl(j)
         zb(j+jmb) = zbl(j)
         ub(j+jmb) = ubl(j)
         vb(j+jmb) = vbl(j)
         wb(j+jmb) = wbl(j)
      enddo
      jmb=jmb+bcc
      if(deb.ge.2) then
        write(37,* ) 'badd ',m7,bcc,jmfl
      endif
      if(baddst.eq.1) then
         write(37,*) 'badd write ',nt,bcc
         call writeParticleList(xbl,ybl,zbl,ubl,vbl,wbl,nt,bcc,
     +                      'badd','bm')
      endif


      call writeParticleList(xfl,yfl,zfl,ufl,vfl,wfl,nt,jmfl,
     +                      'badd','el')

      do j=1,jmfl
         xf(j+jmf) = xfl(j)
         yf(j+jmf) = yfl(j)
         zf(j+jmf) = zfl(j)
         uf(j+jmf) = ufl(j)
         vf(j+jmf) = vfl(j)
         wf(j+jmf) = wfl(j)
      enddo
      jmf=jmf+jmfl

      end
