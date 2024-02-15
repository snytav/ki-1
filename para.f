      ! starting to make a sequential branch
      subroutine PARAinit
      implicit none
      !include 'mpif.h'
      include 'part.pf'
      integer ierr,sz
      character*20 st,st1
      common/sizz/sz

!      call MPI_Init(ierr)
!
!      call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
!
!      call MPI_Comm_size(MPI_COMM_WORLD, sz, ierr)

      write(st,'(2i2.2)') sz,me 
      
      st1='par'//st
      
      if(deb.gt.0) then
        open(37,file=st,form='formatted')
      endif    
      
      call PARAcreateTopology
      
      end

      subroutine PARAcreateTopology
      implicit none
!      include 'mpif.h'
      include 'part.pf'
      integer ierr,sz,grp,grp2,ranks(nproc),i,hsize,vsize,grnum
      character*20 st,st1
      
      hsize = nproc
      vsize = fproc/nproc
        
      meh = mod(me,nproc)
      mev = me/nproc
      grnum = meh
      if(deb.gt.0) write(37,*) grnum, mev
!      call MPI_Comm_split(MPI_COMM_WORLD, grnum,
!     +                     mev, vertComm,ierr)
!
!      call MPI_Comm_group(MPI_COMM_WORLD,grp,ierr)
       
      do i = 1,nproc
         ranks(i) = mev*hsize + (i-1)
      if(deb.gt.0) write(37,*) 'i,vsize,mev,ranks(i),vertComm ',
     +                          i,vsize,mev,ranks(i),vertComm	 
      enddo
!      call MPI_Group_incl(grp,nproc,ranks,grp2,ierr)
      
!      call MPI_Comm_create(MPI_COMM_WORLD,grp2,horComm,ierr)
      if(deb.gt.0) write(37,*) 'meh,horComm ',meh,horComm
      
c      call PARAfinal
c      stop
      end


      subroutine PARAfinal
      implicit none
!      include 'mpif.h'
      integer ierr

!      call MPI_Barrier(MPI_COMM_WORLD,ierr)
!      call MPI_Finalize(ierr)
       
      end 

      subroutine PARAsync()
      implicit none
!      include 'mpif.h'
      integer ierr
      
!      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      
      end

c      subroutine check3D1(me,ix1,iy1,iz1,
c     + imx,imy,imz,erix,eriy,eriz,flag)
c      implicit none
c      include 'part.pf'
c      integer ierr,me,lm1,nt,flag,endimy,imx,imy,imz,erix,eriy,eriz
c      integer i,l,k,ix1,iy1,iz1
c      real*8 ex(imp,lmp,kmp),foo(1000),exc(imp,lmp*nproc,kmp),t,tc,dabs
c
c      common/c3d/ex,exc
c
c      flag = -1
c      erix = -1
c      eriy = -1
c      eriz = -1
c      
c      if(me.eq.(nproc-1)) then
c         endimy = (lmp-2)*me + imy
c      else
c         endimy = (lmp-2)*(me+1)	       
c      endif
c      
c      write(47,*) imx,endimy,imz
c       
c         do i  = ix1,imx
c            do l = (lmp-2)*me+iy1,endimy
c               do k  = iz1,imz
c        	  t  = ex(i,l-(lmp-2)*me,k)       
c                  tc = exc(i,l,k)
c  588             format('i,l,k,t,tc ',3i6,2e25.8)		  
c		  !!!!!write(37,588) i,l,k,t,tc
c                  if(dabs(t-tc).gt.1.0d-15) then
c       		     erix = i
c      		     eriy = l
c      		     eriz = k
c      		     flag = 0
c      		     goto 398
c       		  endif		         	       
c      	       enddo
c            enddo	       
c      	 enddo   
c       	 
c 398 continue   	 
c      end 

          
c      subroutine check1D(me,jmpreal,ex,exc,erix,flag)
c      implicit none
c      include 'part.pf'
c      integer ierr,me,lm1,nt,flag,erix,j,jmpreal
c      real*8 ex(jmp),exc(jmp*nproc),t,tc,dabs
c          
c      flag = -1
c      
c      do j = jmpreal*me+1,(me+1)*jmpreal
c      
c         t  = ex(j-lmp*me)       
c         tc = exc(j)
c        	 
c         write(47,*)t,tc,j,me     
c        	  
c         if(dabs(t-tc).gt.1.0d-8) then
c            erix = j
c      	    flag = 0
c      	    goto 399
c         endif		         	       
c      enddo
c         	 
c  399 continue   	 
c      end 
c
      subroutine PARAsendE1(ex,ey,ez)
      implicit none
!      include 'mpif.h'
      include 'part.pf'
      integer ierr,lm1
!      integer stat(MPI_STATUS_SIZE)
      real*8 sbuf(4*imp*kmp),buf(4*imp*kmp),
     +ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp)
      integer sz,i,l,k,src,dst
      common/sizz/sz
      include 'timempi.par' 

      
      if(nproc.eq.1) then
         do l = 1,imp
            do k = 1,kmp
               ex(l,lmp,k) = 0 !ex(l,2,k)
               ey(l,lmp,k) = 0 !ey(l,2,k)
               ez(l,lmp,k) = 0 !ez(l,2,k)
            enddo
         enddo
         return
      endif
!      call MPI_Barrier(horComm,ierr)

         if(meh.gt.0) then
	    src = meh - 1
	 else
	    src = nproc - 1    
	 endif     
        
         if(meh.lt.(nproc-1)) then
      	    dst = meh + 1
      	 else
       	    dst = 0    
       	 endif     
         
            do l = 1,imp
               do k = 1,kmp
       	          sbuf((l - 1)*kmp + k)           = ex(l,2,k)
                  sbuf((l - 1)*kmp + k+imp*kmp)   = ey(l,2,k)
      	          sbuf((l - 1)*kmp + k+2*imp*kmp) = ez(l,2,k)
      	       enddo
            enddo
      	    
         tm1 = amicro()
!         call MPI_Sendrecv(sbuf,3*imp*kmp,MPI_DOUBLE_PRECISION,
!     +	                      src,888,
!     +                        buf, 3*imp*kmp,MPI_DOUBLE_PRECISION,
!     +                        dst,888,
!     +         	              horComm,stat,ierr)
         call cumultime(tm1,fsrt)         

            do l = 1,imp
               do k = 1,kmp
       	          ex(l,lmp,k) = buf((l - 1)*kmp + k)
      	          ey(l,lmp,k) = buf((l - 1)*kmp + k+imp*kmp)
       	          ez(l,lmp,k) = buf((l - 1)*kmp + k+2*imp*kmp)
      	       enddo
            enddo
        	    
      end

      subroutine PARAsendE2(ex,ez)
      implicit none
!      include 'mpif.h'
      include 'part.pf'
      integer ierr,lm1
!      integer stat(MPI_STATUS_SIZE)
      real*8 sbuf(4*imp*kmp),buf(4*imp*kmp),
     +ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp)
      integer sz,i,l,k,src,dst
      common/sizz/sz
      include 'timempi.par'       
 
      if(nproc.eq.1) then
         do l = 1,imp
            do k = 1,kmp
               ex(l,1,k) = 0 !ex(l,lmp-1,k)
               ez(l,1,k) = 0 !ez(l,lmp-1,k)
            enddo
         enddo
         return
      endif     
!      call MPI_Barrier(horComm,ierr)
          
         if(meh.gt.0) then
       	    src = meh - 1
	 else
	    src = nproc - 1    
	 endif     
        
         if(meh.lt.(nproc-1)) then
      	    dst = meh + 1
      	 else
       	    dst = 0    
       	 endif     
         
            do l = 1,imp
               do k = 1,kmp
       	          sbuf((l - 1)*kmp + k)           = ex(l,lmp-1,k)
                  sbuf((l - 1)*kmp + k+imp*kmp)   = ez(l,lmp-1,k)
      	       enddo
            enddo
      	   
         tm1 = amicro() 
!         call MPI_Sendrecv(sbuf,2*imp*kmp,MPI_DOUBLE_PRECISION,
!     +	                      dst,888,
!     +                        buf, 2*imp*kmp,MPI_DOUBLE_PRECISION,
!     +                        src,888,
!     +         	              horComm,stat,ierr)
         call cumultime(tm1,fsrt)

            do l = 1,imp
               do k = 1,kmp
       	          ex(l,1,k) = buf((l - 1)*kmp + k)
      	          ez(l,1,k) = buf((l - 1)*kmp + k+imp*kmp)
      	       enddo
            enddo
        	    
      end
      
      subroutine PARAreduceP(px,py,pz)
      implicit none
!      include 'mpif.h'
      include 'part.pf'
      integer ierr,lm1
!      integer stat(MPI_STATUS_SIZE)
      real*8 sbufx(lmp*imp*kmp),bufx(lmp*imp*kmp),s1,r1
      real*8 sbufy(lmp*imp*kmp),bufy(lmp*imp*kmp)
      real*8 sbufz(lmp*imp*kmp),bufz(lmp*imp*kmp)
      
      real*8 px(imp,lmp,kmp),py(imp,lmp,kmp),pz(imp,lmp,kmp)
      integer sz,i,l,k,src,dst
      common/sizz/sz
      include 'timempi.par'
c=========================================

c      write(37,*) 'beginning of PARAreduceP'
      
c      do i=1,imp 
c      do l=1,lmp 
c      do k=1,kmp 
c        write(37,*) i,l,k,px(i,l,k),py(i,l,k),pz(i,l,k)
c      enddo
c      enddo
c      enddo
              

c==========================================
      if(deb.ge.2) write(37,*) 'in PARA reduce'

      if(nproc.eq.1) return
      
         if(mev.gt.0) then
            s1 = 3.14e0   
       	    src = mev - 1
	 else
            s1 = 0   
	    src = fproc/nproc - 1    
	 endif     
        
         if(mev.lt.((fproc/nproc)-1)) then
      	    dst = mev + 1
      	 else
       	    dst = 0    
       	 endif     
c      write(37,*) 'b send,mev,src,dst ',s1,r1,mev,src,dst
c      call MPI_Sendrecv(s1,1,MPI_DOUBLE_PRECISION,
c     +	                     src,889,
c     +                        r1, 1,MPI_DOUBLE_PRECISION,
c     +                        dst,889,  
c     +         	              vertComm,stat,ierr)       
c      write(37,*) 'a send ',s1,r1
      
!      call MPI_Barrier(vertComm,ierr)
      
      do i = 1,imp
         do l = 1,lmp
            do k = 1,kmp
               sbufx((i - 1)*lmp*kmp+(l-1)*kmp+ k) = px(i,l,k)
               sbufy((i - 1)*lmp*kmp+(l-1)*kmp+ k) = py(i,l,k)
               sbufz((i - 1)*lmp*kmp+(l-1)*kmp+ k) = pz(i,l,k)
            enddo	       
         enddo
      enddo
      
      if(deb.ge.2) write(37,*) 'b reduce ',sbufx(10)
      tm1 = amicro()
!      call MPI_Allreduce(sbufx,bufx,imp*lmp*kmp,MPI_DOUBLE_PRECISION,
!     +                MPI_SUM,vertComm,ierr)
!      if(deb.ge.2) write(37,*) 'a reduce ',bufx(10)
!      call MPI_Allreduce(sbufy,bufy,imp*lmp*kmp,MPI_DOUBLE_PRECISION,
!     +                MPI_SUM,vertComm,ierr)
!      call MPI_Allreduce(sbufz,bufz,imp*lmp*kmp,MPI_DOUBLE_PRECISION,
!     +                MPI_SUM,vertComm,ierr)
      call cumultime(tm1,colt)
      if(deb.ge.2) write(37,*) 'a reduce ',bufx(10)
        
      do i = 1,imp
         do l = 1,lmp
            do k = 1,kmp
               px(i,l,k) = bufx((i - 1)*lmp*kmp+(l-1)*kmp+ k)
               py(i,l,k) = bufy((i - 1)*lmp*kmp+(l-1)*kmp+ k)
               pz(i,l,k) = bufz((i - 1)*lmp*kmp+(l-1)*kmp+ k)
            enddo	       
         enddo
      enddo
      
c=========================================

c      write(37,*) 'end of PARAreduceP'
      
c      do i=1,imp 
c      do l=1,lmp 
c      do k=1,kmp 
c        write(37,*) 'PARAreduceP',i,l,k,px(i,l,k),py(i,l,k),pz(i,l,k)
c      enddo
c      enddo
c      enddo
                     

c==========================================

       
      end      
      

      subroutine PARAsendP1(px,py,pz)
      implicit none
!      include 'mpif.h'
      include 'part.pf'
      integer ierr,lm1
!      integer stat(MPI_STATUS_SIZE)
      real*8 sbuf(5*imp*kmp),buf(5*imp*kmp),
     +px(imp,lmp,kmp),py(imp,lmp,kmp),pz(imp,lmp,kmp)
      integer sz,i,l,k,src,dst
      common/sizz/sz
      include 'timempi.par'
        
        
      if(deb.ge.3) then
         write(37,*) 'in send1 ',px(2,2,2)
      endif
      
      if(nproc.eq.1) then 
         do l = 1,imp
            do k = 1,kmp
               px(l,lmp-1,k) = px(l,lmp-1,k)+ px(l,1,k) 
               px(l,lmp,k)   = px(l,lmp,k)  + px(l,2,k)
               py(l,lmp-1,k) = py(l,lmp-1,k)+ py(l,1,k)
               pz(l,lmp-1,k) = pz(l,lmp-1,k)+ pz(l,1,k)
               pz(l,lmp,k)   = pz(l,lmp,k)  + pz(l,2,k)
            enddo
         enddo
         if(deb.ge.3) then
            write(37,*) 'leaving send1 ',px(2,2,2)
         endif
         
         return
      endif
      
c      call PARAreduceP(px,py,pz)
           
!      call MPI_Barrier(horComm,ierr)
          
         if(meh.gt.0) then
       	    src = meh - 1
	 else
	    src = nproc - 1    
	 endif     
        
         if(meh.lt.(nproc-1)) then
      	    dst = meh + 1
      	 else
       	    dst = 0    
       	 endif     
         
            do l = 1,imp
               do k = 1,kmp
       	          sbuf((l - 1)*kmp + k)             = px(l,1,k)
                  sbuf((l - 1)*kmp + k+  imp*kmp)   = px(l,2,k)
                  sbuf((l - 1)*kmp + k+2*imp*kmp)   = py(l,1,k)
                  sbuf((l - 1)*kmp + k+3*imp*kmp)   = pz(l,1,k)
                  sbuf((l - 1)*kmp + k+4*imp*kmp)   = pz(l,2,k)
      	       enddo
            enddo
      	    
         tm1 = amicro()
!         call MPI_Sendrecv(sbuf,5*imp*kmp,MPI_DOUBLE_PRECISION,
!     +	                      src,889,
!     +                        buf, 5*imp*kmp,MPI_DOUBLE_PRECISION,
!     +                        dst,889,
!     +         	              horComm,stat,ierr)
         call cumultime(tm1,fsrt)
      do l = 1,imp
         do k = 1,kmp
         px(l,lmp-1,k) = px(l,lmp-1,k)+ buf((l - 1)*kmp + k)
         px(l,lmp,k)   = px(l,lmp,k)  + buf((l - 1)*kmp + k +   imp*kmp)
         py(l,lmp-1,k) = py(l,lmp-1,k)+ buf((l - 1)*kmp + k + 2*imp*kmp)
         pz(l,lmp-1,k) = pz(l,lmp-1,k)+ buf((l - 1)*kmp + k + 3*imp*kmp)
         pz(l,lmp,k)   = pz(l,lmp,k)  + buf((l - 1)*kmp + k + 4*imp*kmp)
         enddo
      enddo
        	    
      end

      subroutine PARAsendP2(px,py,pz)
      implicit none
!      include 'mpif.h'
      include 'part.pf'
      integer ierr,lm1
!      integer stat(MPI_STATUS_SIZE)
      real*8 sbuf(5*imp*kmp),buf(5*imp*kmp),
     +px(imp,lmp,kmp),py(imp,lmp,kmp),pz(imp,lmp,kmp)
      integer sz,i,l,k,src,dst
      common/sizz/sz
      include 'timempi.par'

      if(deb.ge.3) then
         write(37,*) 'in send2 ',px(2,2,2),px(2,lmp,2)
      endif
      
      if(nproc.eq.1) then
         do l = 1,imp
            do k = 1,kmp
               px(l,1,k) = px(l,lmp-1,k)
               px(l,2,k) = px(l,lmp,k)
               py(l,1,k) = py(l,lmp-1,k)
               pz(l,1,k) = pz(l,lmp-1,k)
               pz(l,2,k) = pz(l,lmp,k)
            enddo
         enddo
         if(deb.ge.3) then
            write(37,*) 'leaving send2 ',px(2,2,2),px(2,lmp,2)
         endif
         return
      endif
      !!write(37,*) 'sp2 ',px(5,2,2) ,px(5,lmp,2)       
           
!      call MPI_Barrier(horComm,ierr)
          
         if(meh.gt.0) then
       	    src = meh - 1
	 else
	    src = nproc - 1    
	 endif     
        
         if(meh.lt.(nproc-1)) then
      	    dst = meh + 1
      	 else
       	    dst = 0    
       	 endif     
         
            do l = 1,imp
               do k = 1,kmp
       	          sbuf((l - 1)*kmp + k)             = px(l,lmp-1,k)
                  sbuf((l - 1)*kmp + k+  imp*kmp)   = px(l,lmp,k)
                  sbuf((l - 1)*kmp + k+2*imp*kmp)   = py(l,lmp-1,k)
                  sbuf((l - 1)*kmp + k+3*imp*kmp)   = pz(l,lmp-1,k)
                  sbuf((l - 1)*kmp + k+4*imp*kmp)   = pz(l,lmp,k)
      	       enddo
            enddo
         tm1 = amicro()
!         call MPI_Sendrecv(sbuf,5*imp*kmp,MPI_DOUBLE_PRECISION,
!     +	                      dst,889,
!     +                        buf, 5*imp*kmp,MPI_DOUBLE_PRECISION,
!     +                        src,889,
!     +         	              horComm,stat,ierr)
         call cumultime(tm1,fsrt)
      !!write(37,*) 'sp2 ',px(5,2,2),px(5,lmp,2)               
      
      do l = 1,imp
         do k = 1,kmp
            px(l,1,k) = buf((l - 1)*kmp + k)
            px(l,2,k) = buf((l - 1)*kmp + k +   imp*kmp)
            py(l,1,k) = buf((l - 1)*kmp + k + 2*imp*kmp)
            pz(l,1,k) = buf((l - 1)*kmp + k + 3*imp*kmp)
            pz(l,2,k) = buf((l - 1)*kmp + k + 4*imp*kmp)
         enddo
      enddo
        	    
      !!write(37,*) 'sp2 ',px(5,2,2),px(5,lmp,2)               
      end
       
      subroutine PARAreduce(pp)
      implicit none
!      include 'mpif.h'
      include 'part.pf'
      integer ierr,lm1
!      integer stat(MPI_STATUS_SIZE)
      real*8 pp(omp,omp),sbuf(omp*omp),buf(omp*omp)
      integer sz,i,l,k,src,dst
      common/sizz/sz
      include 'timempi.par'
      !!!!!write(37,*) 'in PARA reduce'

      if(fproc.eq.1) return
           
!      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      
      do i = 1,omp
         do l = 1,omp
           sbuf((i - 1)*omp + l) = pp(i,l)
c     	   !!!write(37,*) 'A i,l,pp(i,l) ',i,l,pp(i,l)
	 enddo
      enddo
         
      if(deb.ge.2) write(37,*) 'PARAreduce '
      tm1 = amicro()
!      call MPI_Allreduce(sbuf,buf,omp*omp,MPI_DOUBLE_PRECISION,
!     +                MPI_SUM,MPI_COMM_WORLD,ierr)
      call cumultime(tm1,colt)
      !!!!!write(37,*) 'a reduce'
        
      do i = 1,omp
         do l = 1,omp
            pp(i,l) = buf((i - 1)*omp + l)
c	    !!!write(37,*) 'B i,l,pp(i,l) ',i,l,pp(i,l)
         enddo
      enddo
      !!!!!write(37,*) 'leaving reduce'
       
      end      

      subroutine PARAreduce12(pp0,pp1,pp2,pp3,pp4,pp5,pp6,
     +                        pp7,pp8,pp9,pp10,pp11,pp12,pp13)
      implicit none
!      include 'mpif.h'
      include 'part.pf'
      integer ierr,lm1
!      integer stat(MPI_STATUS_SIZE)
      real*8 pp0(omp,omp)      
      real*8 pp1(omp,omp),pp2(omp,omp),pp3(omp,omp),pp4(omp,omp)
      real*8 pp5(omp,omp),pp6(omp,omp),pp7(omp,omp),pp8(omp,omp)
      real*8 pp9(omp,omp),pp10(omp,omp),pp11(omp,omp),pp12(omp,omp)
      real*8 pp13(omp,omp)
      real*8 sbuf(14*omp*omp),buf(14*omp*omp)
      integer sz,i,l,k,src,dst
      common/sizz/sz
      include 'timempi.par'
      
      if(deb.ge.2) write(37,*) 'in PARA reduce12 '
c	 call PARAFinal
c	 stop

      if(fproc.eq.1) return
           
!      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      
      do i = 1,omp
         do l = 1,omp
           sbuf((i - 1)*omp + l) = pp1(i,l)
           sbuf((i - 1)*omp + l + omp*omp) = pp2(i,l)
           sbuf((i - 1)*omp + l + 2*omp*omp) = pp3(i,l)
           sbuf((i - 1)*omp + l + 3*omp*omp) = pp4(i,l)
           sbuf((i - 1)*omp + l + 4*omp*omp) = pp5(i,l)
           sbuf((i - 1)*omp + l + 5*omp*omp) = pp6(i,l)
           sbuf((i - 1)*omp + l + 6*omp*omp) = pp7(i,l)
           sbuf((i - 1)*omp + l + 7*omp*omp) = pp8(i,l)
           sbuf((i - 1)*omp + l + 8*omp*omp) = pp9(i,l)
           sbuf((i - 1)*omp + l + 9*omp*omp) = pp10(i,l)
           sbuf((i - 1)*omp + l + 10*omp*omp) = pp11(i,l)
           sbuf((i - 1)*omp + l + 11*omp*omp) = pp12(i,l)
           sbuf((i - 1)*omp + l + 12*omp*omp) = pp0(i,l)
           sbuf((i - 1)*omp + l + 13*omp*omp) = pp13(i,l)
c     	   !!!write(37,*) 'A i,l,pp(i,l) ',i,l,pp(i,l)
	 enddo
      enddo
         
      if(deb.ge.3) then
         write(37,*) 'b reduce',buf(omp/4*3+1)
         do l = 1,omp
           write(37,*) l,buf((i - 1)*omp + l)
	 enddo  
      endif
      if(deb.ge.2) write(37,*) 'PARAreduce12 '
      tm1 = amicro()
!      call MPI_Allreduce(sbuf,buf,14*omp*omp,MPI_DOUBLE_PRECISION,
!     +                MPI_SUM,MPI_COMM_WORLD,ierr)
      call cumultime(tm1,colt)
      if(deb.ge.3) then
         write(37,*) 'a reduce',buf(omp/4*3+1)
         do l = 1,omp
           write(37,*) l,buf((i - 1)*omp + l)
	 enddo  
      endif
          
      do i = 1,omp
         do l = 1,omp
           pp1(i,l)  = buf((i - 1)*omp + l)
           pp2(i,l)  = buf((i - 1)*omp + l + omp*omp)
           pp3(i,l)  = buf((i - 1)*omp + l + 2*omp*omp)
           pp4(i,l)  = buf((i - 1)*omp + l + 3*omp*omp)
           pp5(i,l)  = buf((i - 1)*omp + l + 4*omp*omp)
           pp6(i,l)  = buf((i - 1)*omp + l + 5*omp*omp)
           pp7(i,l)  = buf((i - 1)*omp + l + 6*omp*omp)
           pp8(i,l)  = buf((i - 1)*omp + l + 7*omp*omp)
           pp9(i,l)  = buf((i - 1)*omp + l + 8*omp*omp)
           pp10(i,l) = buf((i - 1)*omp + l + 9*omp*omp)
           pp11(i,l) = buf((i - 1)*omp + l + 10*omp*omp)
           pp12(i,l) = buf((i - 1)*omp + l + 11*omp*omp)
           pp0(i,l)  = buf((i - 1)*omp + l + 12*omp*omp)
           pp13(i,l) = buf((i - 1)*omp + l + 13*omp*omp)
         enddo
      enddo
      if(deb.ge.2) write(37,*) 'leaving reduce 12'
       
      end      


      subroutine PARAreduce1(pp)
      implicit none
!      include 'mpif.h'
      include 'part.pf'
      integer ierr,lm1
!      integer stat(MPI_STATUS_SIZE)
      real*8 pp(imp-1,lmp2),sbuf((imp-1)*lmp2),buf((imp-1)*lmp2)
      integer sz,i,l,k,src,dst
      common/sizz/sz
      include 'timempi.par'
      !!!!!write(37,*) 'in PARA reduce'

      if(fproc.eq.1) return
           
!      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      
      do i = 1,imp-1
         do l = 1,lmp2
           sbuf((i - 1)*lmp2 + l) = pp(i,l)
c       	   !!!!write(37,*) 'AA i,l,pp(i,l) ',i,l,pp(i,l)
	 enddo
      enddo
         
      if(deb.ge.2) write(37,*) 'PARAreduce1 '
      tm1 = amicro()
!      call MPI_Allreduce(sbuf,buf,(imp-1)*lmp2,MPI_DOUBLE_PRECISION,
!     +                MPI_SUM,MPI_COMM_WORLD,ierr)
      call cumultime(tm1,colt)
      !!!!!write(37,*) 'a reduce'
        
      do i = 1,imp-1
         do l = 1,lmp2
            pp(i,l) = buf((i - 1)*lmp2 + l)
c	    !!!!write(37,*) 'BB i,l,pp(i,l) ',i,l,pp(i,l)
         enddo
      enddo
      !!!!!write(37,*) 'leaving reduce'
       
      end      

      subroutine delpart(j,n3,x,y,z,pu,pv,pw)
      implicit none
      include 'part.pf'
      real*8 x(jmp),y(jmp),z(jmp),pu(jmp),pv(jmp),pw(jmp)
      integer j,n3
      
      if(j.lt.n3) then
         x(j)   = x(n3)
         z(j)   = z(n3)
         y(j)   = y(n3)
         pu(j)   = pu(n3)
         pv(j)   = pv(n3)
         pw(j)   = pw(n3)
      endif
      n3 = n3-1
      end   

      
      subroutine PARAtransferPart(m4,ym,n3,xi,yi,zi,ui,vi,wi)
      implicit none
!      include 'mpif.h'
      include 'part.pf' 
      real*8 xi(jmp),yi(jmp),zi(jmp),ui(jmp),vi(jmp),wi(jmp)
      integer ierr,lm1,m4,flag,sp,n3
!      integer stat(MPI_STATUS_SIZE)
      real*8  lbuf(bjmp),rbuf(bjmp),lrbuf(bjmp),rrbuf(bjmp) 
      real*8  y00,ym,yg,yg1
      integer sz,i,l,k,src,dst,j,lc,rc,lguest,rguest,ff
      common/ptranstat/lc,rc
      include 'timempi.par'
      
      if(deb.ge.4) write(37,*) 'in PARAtransf ',n3,nproc
     
      if(nproc.eq.1) then
         do j = 1,n3 
            if(yi(j).lt.0) then
	       yi(j) = yi(j) + ym
 751           format('partDown x,y,y1',3e30.20)
               write(37,751) xi(j),yi(j),yi(j)+ym
	    endif
            if(yi(j).gt.ym) then
	       yi(j) = yi(j) - ym
 750           format('partUp x,y,y1',3e30.20)
               write(37,750) xi(j),yi(j),yi(j)-ym
	    endif
	 enddo    
         return
      endif
      y00 = ym*meh
      
      lc = 0
      rc = 0
      if(deb.ge.1) then
        !!write(37,*) m4, '=== in transfer======================= ', 
c     +	me,nproc,y00,ym     
      endif    
      
         if(meh.gt.0) then
       	    src = meh - 1
      	 else
       	    src = nproc - 1    
      	 endif     
        
         if(meh.lt.(nproc-1)) then
      	    dst = meh + 1
      	 else
       	    dst = 0    
       	 endif     
      
      if(deb.ge.1) then
        !!write(37,*) 'in transfer src,dst ',src,dst,me     
      endif    
      
      j = 1
 962     flag = 0
         if((yi(j).gt.(y00+ym)).and.((rc*8).lt.bjmp)) then
            if(meh.eq.(nproc-1)) then
	       yi(j) = yi(j) - ym*nproc
	    endif 

c	    if(m4.eq.60) then
   	       !!!!!write(37,*) 'R ',j,n3,x(j),y(j),pu(j),pv(j)
c	    endif   
      	    rbuf(rc*8+1) = xi(j)
      	    rbuf(rc*8+2) = yi(j)
      	    rbuf(rc*8+3) = zi(j)
      	    rbuf(rc*8+4) = ui(j)
      	    rbuf(rc*8+5) = vi(j)
      	    rbuf(rc*8+6) = wi(j)
            rc = rc + 1	    
            flag = 1
      	 else     
            if((yi(j).le.y00).and.((lc*8).lt.bjmp)) then
               if(meh.eq.0) then
	          yi(j) = yi(j) + ym*nproc
	       endif 
c	       if(m4.eq.60) then
   	          !!!!!write(37,*) 'L ',j,n3,x(j),y(j),pu(j),pv(j)
c	       endif   
      	       lbuf(lc*8+1) = xi(j)
      	       lbuf(lc*8+2) = yi(j)
      	       lbuf(lc*8+3) = zi(j)
      	       lbuf(lc*8+4) = ui(j)
      	       lbuf(lc*8+5) = vi(j)
      	       lbuf(lc*8+6) = wi(j)
       	       lc = lc + 1
       	       flag = 1
       	    endif   
      	 endif     
      	 if(flag.eq.1) then
            call delPart(j,n3,xi,yi,zi,ui,vi,wi)	    
      	 else
      	    j = j + 1
      	 endif
      
      if(j.le.n3) goto 962
      
      if(deb.ge.3) then
         write(37,*) 'in transfer lc,rc' ,lc,rc,src,dst,meh    
      endif     
      
!      call MPI_Sendrecv(lc,1,MPI_INTEGER,src,881,
!     +                  rguest,1,MPI_INTEGER,
!     +                  dst,881,horComm,stat,ierr)
!
!      call MPI_Sendrecv(rc,1,MPI_INTEGER,dst,881,
!     +                  lguest,1,MPI_INTEGER,
!     +                  src,881,horComm,stat,ierr)
        
      if(deb.ge.3) then
         write(37,*) 'lc,lguest,rc,rguest ',lc,lguest,rc,rguest
      endif     
      tm1 = amicro()
!      call MPI_Sendrecv(lbuf,(lc+1)*8,MPI_DOUBLE_PRECISION,src,883,
!     +                  rrbuf,(rguest+1)*8,MPI_DOUBLE_PRECISION,
!     +                  dst,883,horComm,stat,ierr)
!
!      call MPI_Sendrecv(rbuf,(rc+1)*8,MPI_DOUBLE_PRECISION,dst,883,
!     +                  lrbuf,(lguest+1)*8,MPI_DOUBLE_PRECISION,
!     +                  src,883,horComm,stat,ierr)
      call cumultime(tm1,psrt)
      if(deb.ge.1) then
         write(37,*) '2 send in transfer ',lc,lguest,rc,rguest,m4      
         write(37,*) 'present n3 ',n3,m4
      endif     
c      call PARAFinal 
c      stop
      do l = 1,lguest
         yg  = lrbuf((l-1)*8 + 2)
         sp  = lrbuf((l-1)*8 + 8)
         
c	 !!write(37,*) l,lguest,n3,jmp,yg
         if(n3.lt.(jmp-1)) then
	    n3 = n3 + 1
	    
	    xi(n3) = lrbuf((l-1)*8 + 1)
	    yi(n3) = yg
	    zi(n3) = lrbuf((l-1)*8 + 3)
	    ui(n3) = lrbuf((l-1)*8 + 4)
	    vi(n3) = lrbuf((l-1)*8 + 5)
	    wi(n3) = lrbuf((l-1)*8 + 6)
	 endif    
c	 !!write(37,*) 'yg ',yi(n3)
      enddo

      do l = 1,rguest
         yg  = rrbuf((l-1)*8 + 2)
         sp  = rrbuf((l-1)*8 + 8)

c	 !!write(37,*) l,lguest,n3,jmp,yg
         if(n3.lt.(jmp-1)) then
	    n3 = n3 + 1
	    
	    xi(n3) = rrbuf((l-1)*8 + 1)
	    yi(n3) = yg
	    zi(n3) = rrbuf((l-1)*8 + 3)
	    ui(n3) = rrbuf((l-1)*8 + 4)
	    vi(n3) = rrbuf((l-1)*8 + 5)
	    wi(n3) = rrbuf((l-1)*8 + 6)
	 endif    
c	 !!write(37,*) 'yg ',yi(n3)
      enddo
      
      if(deb.ge.1) then
         !!write(37,*) 'accepted ',lguest+rguest
         !!write(37,*) 'now ',n3,me
      endif     
      
      
      end
      
      subroutine PARAreduceLMP2(vec)
      implicit none
      include 'part.pf'
!      include 'mpif.h'
      real*8 vec1(1000),vec2(1000),vec(1000)
      integer ntot,nn3,ierr,pcmin,pcmax,l
      include 'timempi.par' 
      
      do l = 1,1000
         vec1(l) = vec(l)
      enddo
      call PARAsync
      if(deb.ge.3) then
         write(37,*) 'lmp2 1 ',vec1(1),vec1(4)
      endif
      
      if(deb.ge.2) write(37,*) ' reduceLMP2 '
      tm1 = amicro()
!      call MPI_Allreduce(vec1,vec2,1000,MPI_DOUBLE_PRECISION,
!     +                MPI_SUM,MPI_COMM_WORLD,ierr)
      call cumultime(tm1,colt)
      if(deb.ge.3) then
         write(37,*) 'lmp2 2 ',vec2(1),vec2(4)
      endif
      
      do l = 1,1000
         vec(l) = vec2(l)
      enddo

      end
      
      subroutine PARAgatherGrow(p,ppar)
      implicit none
!      include 'mpif.h'
      include 'part.pf'
      integer ierr,lm1
!      integer stat(MPI_STATUS_SIZE)
      real*8 p(imp,lmp,kmp),ppar(imp,lmpf,kmp),
     +       rbuf(imp*lmpf*kmp),buf(imp*lmpf*kmp)     
      integer sz,i,l,k,src,dst
      common/sizz/sz
      include 'timempi.par'

      if(deb.ge.2) then
         write(37,*) 'in PARA gatherGG'
      endif	 

      if(nproc.eq.1) then
         do i = 1,imp
            do l = 2,lmp-1 
               do k = 1,kmp
      	          ppar(i,l,k) = p(i,l,k)
	       enddo 
	    enddo
         enddo
      
         return
      endif 	 
           
!      call MPI_Barrier(horComm,ierr)
      
      do i = 1,imp
         do l = 2,lmp-1 
            do k = 1,kmp
      	       buf((l-2)*imp*kmp+(i-1)*kmp+k) = p(i,l,k)
      	    enddo 
      	 enddo
      enddo
      do l = 1,lmp-1
         if(deb.ge.2) then
c	    write(37,*) 
c     +	    l,buf((l+me*(lmp-2))*imp*kmp+(2-1)*kmp+2),p(2,l,2)
         endif
      enddo 
c      do l = lmp,lmpf-1
c         if(deb.ge.2) then
c	    write(37,*) l,
c     +	    buf(l*imp*kmp+(2-1)*kmp+2)
c         endif 
c      enddo
         
      if(deb.ge.2) then
         write(37,*) 'b gath ',meh,horComm
      endif	 
      tm1 = amicro()
!      call MPI_Allgather(buf, (lmp-2)*imp*kmp,MPI_DOUBLE_PRECISION,
!     +                rbuf,(lmp-2)*imp*kmp,MPI_DOUBLE_PRECISION,
!     +                horComm,ierr)
      call cumultime(tm1,gatht)
      if(deb.ge.2) then
         write(37,*) 'PARAGatherGrow a gath ',meh,horComm
      endif	 
          
      do i = 1,imp
         do k = 1,kmp
	    ppar(i,1,k)    = 0d0
	    ppar(i,lmpf,k) = 0d0
            do l = 2,lmpf-1 
	       ppar(i,l,k) = rbuf((l-2)*imp*kmp+(i-1)*kmp+k) 
 	    enddo 
	 enddo
      enddo
c      do l = lmp,lmpf-1
c         if(deb.ge.2) then 
c	    write(37,*) l,rbuf(l*imp*kmp+(2-1)*kmp+2)
c	 endif     
c      enddo
      do l = 1,lmpf-1
         if(deb.ge.2) then
	    write(37,*) l,ppar(2,l,2),rbuf(l*imp*kmp+(2-1)*kmp+2)
	 endif    
      enddo 

      if(deb.ge.2) then
         write(37,*) 'leaving gath'
      endif	 
       
      end      


      subroutine PARAgather(p,ppar)
      implicit none
!      include 'mpif.h'
      include 'part.pf'
      integer ierr,lm1
!      integer stat(MPI_STATUS_SIZE)
      real*8 p(imp,lmp,kmp),ppar(imp,lmpf,kmp),
     +       rbuf(imp*lmpf*kmp),buf(imp*lmpf*kmp)     
      integer sz,i,l,k,src,dst
      common/sizz/sz
      include 'timempi.par'
      
      if(deb.ge.2) then
         write(37,*) 'in PARA gather ',horComm
      endif	 

      if(nproc.eq.1) then
         do i = 1,imp
            do l = 2,lmp-1 
               do k = 1,kmp
      	          ppar(i,l,k) = p(i,l,k)
	       enddo 
	    enddo
         enddo
      
         return
      endif 	 

      if(deb.ge.2) then
         write(37,*) 'in PARA gather2 ',horComm
      endif	 
      
           
c      call MPI_Barrier(horComm,ierr)
c      call PARAfinal
c      stop
      
      do i = 1,imp
         do l = 2,lmp-1 
            do k = 1,kmp
      	       buf((l-2)*imp*kmp+(i-1)*kmp+k) = p(i,l,k)
      	    enddo 
      	 enddo
      enddo
         
      if(deb.ge.2) then
         write(37,*) 'b gath ',meh,horComm
      endif	 

      tm1 = amicro()
!      call MPI_Gather(buf, (lmp-2)*imp*kmp,MPI_DOUBLE_PRECISION,
!     +                rbuf,(lmp-2)*imp*kmp,MPI_DOUBLE_PRECISION,
!     +                0,horComm,ierr)
      call cumultime(tm1,gatht)
      if(deb.ge.2) then
         write(37,*) 'PARAgather a gath ',meh,horComm
      endif	 
          
      do i = 1,imp
         do k = 1,kmp
	    ppar(i,1,k)    = 0d0
	    ppar(i,lmpf,k) = 0d0
            do l = 2,lmpf-1 
	       ppar(i,l,k) = rbuf((l-2)*imp*kmp+(i-1)*kmp+k) 
	    enddo 
	 enddo
      enddo
c      do l = lmp,lmpf-1
c         if(deb.ge.2) then 
c	    write(37,*) l,rbuf(l*imp*kmp+(2-1)*kmp+2)
c	 endif     
c      enddo
      
      do l = 1,lmpf-1
         if(deb.ge.2) then
	    write(37,*) l,ppar(2,l,2),rbuf(l*imp*kmp+(2-1)*kmp+2)
	 endif    
      enddo 

      if(deb.ge.2) then
         write(37,*) 'leaving gath'
      endif	 
       
      end      

      subroutine PARAmm(pcmin,pcmax)
      implicit none
      include 'part.pf'
!      include 'mpif.h'
      real*8 ad,dl,pcav,tmass,buf(10),rbuf(10)
      integer ntot,nn3,ierr,pcmin,pcmax
      include 'timempi.par'

      if(fproc.eq.1) return
      
      buf(1)   = pcmin
      buf(2)   = -pcmax
      
      if(deb.ge.2) write(37,*) 'PARAmm '
      tm1 = amicro()
!      call MPI_Allreduce(buf,rbuf,2,MPI_DOUBLE_PRECISION,
!     +                MPI_MIN,MPI_COMM_WORLD,ierr)
      call cumultime()
      pcmin   =  rbuf(1)
      pcmax   = -rbuf(2)   

      end

      subroutine PARAreduce8(rde,rdi,dne1,dni1,dte1,dti1,rli,rle)
      implicit none
      include 'part.pf'
!      include 'mpif.h'
      real*8 rde(imp),rdi(imp),dne1(imp),dni1(imp),dte1(imp),dti1(imp),
     +       rli(imp),rle(imp),buf(8*imp),rbuf(8*imp)
      integer i,ierr 
      include 'timempi.par'
         
c      print*,'qq'
c      call PARAfinal
c      stop
      
      if(fproc.eq.1) return
      
      do i = 1,imp
         buf(i)       =  rde(i)
         buf(i+imp)   =  rdi(i)
         buf(i+2*imp) =  dne1(i)
         buf(i+3*imp) =  dni1(i)
         buf(i+4*imp) =  dte1(i)
         buf(i+5*imp) =  dti1(i)
         buf(i+6*imp) =  rli(i)
         buf(i+7*imp) =  rle(i)
      enddo

c      print*,'qq'
c      call PARAfinal
c      stop
      
      if(deb.ge.2) write(37,*) 'b reduce8 '
      tm1 = amicro()
!      call MPI_Allreduce(buf,rbuf,8*imp,MPI_DOUBLE_PRECISION,
!     +                   MPI_SUM,MPI_COMM_WORLD,ierr)
      call cumultime(tm1,colt)
      do i = 1,imp
         rde(i)  = rbuf(i)/nproc
         rdi(i)  = rbuf(i+imp)/nproc
         dne1(i) = rbuf(i+2*imp)/nproc
         dni1(i) = rbuf(i+3*imp)/nproc
         dte1(i) = rbuf(i+4*imp)/nproc
         dti1(i) = rbuf(i+5*imp)/nproc
         rli(i)  = rbuf(i+6*imp)/nproc
         rle(i)  = rbuf(i+7*imp)/nproc
      enddo


      end  
      
      subroutine PARAgather3(topn,topg,toph,m7,topres2)
      implicit none
      include 'part.pf'
!      include 'mpif.h'
      real*8 topg(grtot),toph(grtot),topres(2*grtot),topres2(2*grtot)
      real*8 sbuf(3*grtot),rbuf(3*grtot*nproc),gmax
      integer topn(grtot),m7,n1,n2,i,j,ierr,jmax
      common/tpr/topres
      include 'timempi.par'
      
      if(m7.eq.1) then
         do i  = 1,2*grtot
	    topres(i) = 0
	 enddo
      endif
      
      if(nproc.eq.1) then
         do j = 1,grtot
            topres(2*j-1) = topn(j)     
            topres(2*j)   = toph(j)     
         enddo
         
         return
      endif
      
      do j = 1,grtot
         sbuf((j-1)*3+1) = topn(j)     
         sbuf((j-1)*3+2) = topg(j)     
         sbuf((j-1)*3+3) = toph(j)     
      enddo

      tm1 = amicro()
!      call MPI_Gather(sbuf, 3*grtot,MPI_DOUBLE_PRECISION,
!     +                rbuf, 3*grtot,MPI_DOUBLE_PRECISION,
!     +                0,horComm,ierr)
      call cumultime(tm1,gatht)
      
      do i = 1,grtot
         if(topres(2*i-1).gt.0) then
	    n1 = topres(2*i-1)
	    do j = 1,grtot*nproc
	       n2 = rbuf((j-1)*3+1)
	       if(n1.eq.n2) then
	          topres(2*i) = rbuf((j-1)*3+3)
	       endif
	    enddo
 119        continue	    
	 else
            gmax = 0d0
	    do j = 1,grtot*nproc
	       n2 = rbuf((j-1)*3+1)
	       if((n2.gt.0).and.(gmax.lt.rbuf((j-1)*3+2))) then
		  gmax = rbuf((j-1)*3+2)
		  jmax = j
	       endif
	    enddo    
	    if(gmax.gt.0) then
	       topres(2*i-1) = rbuf((jmax-1)*3+1)
	       topres(2*i)   = rbuf((jmax-1)*3+3)
	       gmax = 0d0
	       rbuf((jmax-1)*3+1) = 0
            endif
	    
  120       continue	    
	 endif
      enddo      
      
      do i  = 1,2*grtot
         topres2(i) = topres(i)
      enddo

      end 


      subroutine PARAreduce6(pp)
      implicit none
!      include 'mpif.h'
      include 'part.pf'
      integer ierr,lm1
!      integer stat(MPI_STATUS_SIZE)
      real*8 pp(omp,omp,6)      
      real*8 sbuf(6*omp*omp),buf(6*omp*omp)
      integer sz,i,l,k,src,dst
      common/sizz/sz
      include 'timempi.par'
      
      if(deb.ge.2) write(37,*) 'in PARA reduce12 '
c	 call PARAFinal
c	 stop

      if(fproc.eq.1) return
           
!      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      
      do i = 1,omp
         do l = 1,omp
            do k=1,6
               sbuf((i - 1)*omp + l + (k-1)*omp*omp) = pp(i,l,k)
            enddo
	 enddo
      enddo
         
      if(deb.ge.2) write(37,*) 'b reduce6 '
      tm1 = amicro()
!      call MPI_Allreduce(sbuf,buf,6*omp*omp,MPI_DOUBLE_PRECISION,
!     +                MPI_SUM,MPI_COMM_WORLD,ierr)
      call cumultime(tm1,colt)
      do i = 1,omp
         do l = 1,omp
            do k=1,6
               pp(i,l,k) = buf((i - 1)*omp + l + (k-1)*omp*omp)
            enddo
	 enddo
      enddo

       
      end      
    
      subroutine PARAsetRandomIon
      implicit none
!      include 'mpif.h'
      include 'part.pf'
      integer times,i,k,j,lp,jmf,jmi,n1,n2,m5,m7,jmb
      real*8 t,pi,lx,ly,lz,x,wrapg05cae,wrapg05dde
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
      common/b/lx,ly,lz

      if(fproc.eq.1) return


      call setrngdebug
c      write(37,*) 'jmf,jmb,jmi ',jmf,jmb,jmi
         do k=1,jmi*me
  	    t = lz*wrapg05cae (lz)
   	    t = wrapg05cae (ly)
            x = wrapg05cae (lx)
         enddo
c      if(deb.ge.2) write(37,*) 'end setrng ion'         

      end

      subroutine PARAsetRandomBeam
      implicit none
!      include 'mpif.h'
      include 'part.pf'
      integer times,i,k,j,lp,jmf,jmi,n1,n2,m5,m7,jmb
      real*8 t,pi,lx,ly,lz,x,wrapg05cae,wrapg05dde
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
      common/b/lx,ly,lz

      if(fproc.eq.1) return
      jmb = jmf/2
      
c      do times = 1,me
c ostatok ot ionov
        
         do k=1,jmi*(fproc - (me+1))
       	    t = lz*wrapg05cae (lz)
      	    t = wrapg05cae (ly)
            x = wrapg05cae (lx)
         enddo
c
         if(beamf.eq.1) then
            do i=1,jmb*me
               j=j+1
               t=wrapg05dde(0.d0,Tb*rimp)        
               t=wrapg05dde(0.d0,Tb*rimp)        
               t=wrapg05dde(0.d0,Tb*rimp)        
            enddo
         endif	
c       if(deb.ge.2) write(37,*) 'end setrng beam'

      end
        
      subroutine PARAsetRandomElectrons
      implicit none
!      include 'mpif.h'
      include 'part.pf'
      integer times,i,k,j,lp,jmf,jmi,n1,n2,m5,m7,jmb
      real*8 t,pi,lx,ly,lz,x,wrapg05cae,wrapg05dde
      common/e/pi,lp,jmf,jmi,n1,n2,m5,m7,jmb
      common/b/lx,ly,lz
        
      if(fproc.eq.1) return
      jmb = jmf/2
      
         if(beamf.eq.1) then
            do i=1,jmb*(fproc - (me+1))
               t = wrapg05dde(0.d0,Tb*rimp)        
               t = wrapg05dde(0.d0,Tb*rimp)        
               t = wrapg05dde(0.d0,Tb*rimp)        
            enddo
            do i=1,jmb*me
               j=j+1
               t = wrapg05dde(0.d0,tey0)        
               t = wrapg05dde(0.d0,tez0)        
               t = wrapg05dde(0.d0,tex0)        
            enddo
         endif	
       if(deb.ge.2) write(37,*) 'end setrng electrons'
      end
        
        
      subroutine setrngdebug
      implicit none
      real*8 t
      integer ndde,ncae
      common/rng/ndde,ncae
      
      ndde = 0
      ncae = 0  
      end

      real*8 function wrapg05dde(a,b) 
      implicit none
      include 'part.pf'
      real*8 t,a,b,g05dde
      integer ndde,ncae
      common/rng/ndde,ncae

      t = g05dde(a,b)
      ndde = ndde + 1
      if(deb.ge.2) write(37,787) ndde,t
 787  format('g05dde ',i10,e30.20)         
 
      wrapg05dde = t
      end

      real*8 function wrapg05cae(a) 
      implicit none
      include 'part.pf'
      real*8 t,a,g05cae
      integer ndde,ncae
      common/rng/ndde,ncae

      t = g05cae(a)
      ncae = ncae + 1
      if(deb.ge.3) write(37,788) ncae,t
 788  format('g05cae ',i10,e30.20)         
 
      wrapg05cae = t
      end

      subroutine getParticlesOnMinorPEs(amf,jmf,jmfm)
      implicit none
      include 'part.pf'
!      include 'mpif.h'
      real*8 amf
      integer all_jmf(fproc),ierr,i,jmf,jmfm
      integer all_jmf_res(fproc)
      include 'timempi.par'

      write(37,*) 'begin minorPE ',me,amf,jmf

      do i = 1,fproc
         all_jmf(i) = 0
      enddo

      all_jmf(me+1) = jmf

      do i = 1,fproc
 1399    format('minorPEb ',e12.3,2i5)
         write(37,1399) amf,me,all_jmf(i)
      enddo

      tm1 = amicro()
!      call MPI_Allreduce(all_jmf,all_jmf_res,fproc,MPI_INT,
!     +                   MPI_SUM,MPI_COMM_WORLD,ierr)
      call cumultime(tm1,colt)

      do i = 1,fproc
 1398    format('minorPEa ',e12.3,2i5)
         write(37,1398) amf,me,all_jmf_res(i)
      enddo


      jmfm = 0
      do i = 1,me
         jmfm=jmfm+all_jmf_res(i)
         write(37,*) 'minorPE ',amf,i,me,jmfm,all_jmf_res(i)
      enddo
      write(37,*) 'end minorPE ',me,amf,jmfm
      

      end subroutine getParticlesOnMinorPEs   

      subroutine writeParticleList(x,y,z,pu,pv,pw,nt,jm,wh,sort)
      implicit none
      include 'part.pf'
      real*8 x(jmp),y(jmp),z(jmp),pu(jmp),pv(jmp),pw(jmp)
      character*2  sort
      character*4  wh
      character*5  num
      character*40 fn,tr1
      integer      nt,jm,j

      if(jm.eq.0) return

      if(deb.ge.4) then
         write(37,*) 'writeParticleList ',me,sort,jm,wh
      else
         return
      endif


      !print*,'sort ',sort


      fn=sort
      write(num,'(i5.5)') me
      fn=fn(1:2)//num
      fn=fn(1:7)//'nt'
      write(tr1,'(i3.3)') nt
      fn=fn(1:9)//tr1
      fn=fn(1:12)//wh
      fn=fn(1:16)//'.dat'

      open(38,file=fn,form='formatted')
      
 1479 format(6e30.20)

      do j = 1,jm
         write(38,1479) x(j),y(j),z(j),pu(j),pv(j),pw(j)
      enddo

      close(38)

      end subroutine writeParticleList

      subroutine writeFieldArrays(ex,ey,ez,hx,hy,hz,
     +                             qx,qy,qz,jx,jy,jz,
     +                             nt,wh)
      implicit none
      include 'part.pf'
      real*8 ex(imp,lmp,kmp),ey(imp,lmp,kmp),ez(imp,lmp,kmp)
      real*8 hx(imp,lmp,kmp),hy(imp,lmp,kmp),hz(imp,lmp,kmp)
      real*8 qx(imp,lmp,kmp),qy(imp,lmp,kmp),qz(imp,lmp,kmp)
      real*8 jx(imp,lmp,kmp),jy(imp,lmp,kmp),jz(imp,lmp,kmp)
      integer nt,i,l,k
      character*4  wh
      character*5  num
      character*40 fn,tr1


      if(deb.lt.4) return

      if(deb.ge.2) write(37,*) 'writeParticleList ',me

      fn='fd'
      write(num,'(i5.5)') me
      fn=fn(1:2)//num
      fn=fn(1:7)//'nt'
      write(tr1,'(i3.3)') nt
      fn=fn(1:9)//tr1
      fn=fn(1:12)//wh
      fn=fn(1:16)//'.dat'



      open(38,file=fn,form='formatted')

      do i = 1,imp
         do l = 1,lmp
            do k = 1,kmp
               write(38,1476) i,l,k,ex(i,l,k),ey(i,l,k),ez(i,l,k),
     +                              hx(i,l,k),hy(i,l,k),hz(i,l,k),
     +                              qx(i,l,k),qy(i,l,k),qz(i,l,k),
     +                              jx(i,l,k),jy(i,l,k),jz(i,l,k)
            enddo
         enddo
      enddo
      close(38)

 1476 format(3i5,12e30.20)

      end subroutine writeFieldArrays

      subroutine writeParticleCurrents(jx,jy,jz,j,nt,wh,
     +                                 x,y,z,px,py,pz)
      implicit none
      include 'part.pf'
      real*8 jx(imp,lmp,kmp),jy(imp,lmp,kmp),jz(imp,lmp,kmp)
      integer nt,i,l,k,j
      character*4  wh
      character*5  num
      character*40 fn,tr1
      real*8       x,y,z,px,py,pz


      if(deb.lt.5) return

      !if(deb.ge.2) write(37,*) 'writeParticleList ',me

      fn='pj'
      write(num,'(i5.5)') me
      fn=fn(1:2)//num
      fn=fn(1:7)//'nt'
      write(tr1,'(i3.3)') nt
      fn=fn(1:9)//tr1
      fn=fn(1:12)//wh
      write(tr1,'(i8.8)') j
      fn=fn(1:16)//tr1
      fn=fn(1:24)//'.dat'



      open(38,file=fn,form='formatted')
1527  format(6e30.20)
      write(38,1527) x,y,z,px,py,pz

      do i = 1,imp
         do l = 1,lmp
            do k = 1,kmp
               write(38,1476) i,l,k,jx(i,l,k),jy(i,l,k),jz(i,l,k)
            enddo
         enddo
      enddo
      close(38)

 1476 format(3i5,12e30.20)

      end subroutine writeParticleCurrents



      
      
