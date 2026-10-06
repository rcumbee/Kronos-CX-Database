cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c radcash                                                          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    calculates cascade of radiative transitions for one-electron  
c    excited state of the hydrogenic atom               
c           
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Copyright 2017 P. D. Mullen, R. S. Cumbee, P. C. Stancil
c
c                        Kronos version 3.1   (released May 16, 2017)
c
c    This file is part of Kronos. Kronos is a series of python and fortran codes 
c    utilizing a comprehensive charge exchange database at its backbone to generate 
c    theoretical line ratios and spectra for many single electron capture systems. 
c
c    You are morally obligated to cite the following paper in his/her
c    scientific literature that results from use of any part of Kronos:
c
c    [1] Mullen et al. (2016), ApJS 224, 31
c    [2] Cumbee et al. (2017), ApJ, in prep
c    [3] Lyons et al. (2017), ApJS, in prep
c    [4] Mullen et al. (2017), ApJ, in prep
c        
c
c    Kronos is free software; you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation; either version 3 of the License, or
c    (at your option) any later version.
c
c    Kronos is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Milun Rakovic  July 2006
c    pcs 5-28-12 add Hardness ratio calculation and comments 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    input - init.dat  - relative cross sections for each nl level
c    output - many files, but most important are:
c             phot* and ratios.dat
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Initializations
c    nnn is the maximum possible n value in this code.
c    energy is the collision energy 
c    potential is the ionization energy of species
c    rmax is the maximum radius of radial integrations
c    q is the charge 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit real*8(a-h,o-z)
      real sump, sump2, energy
      parameter (nnn=40)
      dimension c(nnn,nnn), cnew(nnn,nnn)
      dimension phot(nnn,nnn,nnn,nnn)
      common rmax,q,potential

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
c     initialize three arrays, c, cnew, and phot 
c     to be 0 to max possible n. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

      sump=0.0d0
      sump2=0.0d0

      do i=1,nnn
      do j=1,nnn
      c(i,j)=0.d0
      cnew(i,j)=0.d0
      do i1=1,nnn
      do j1=1,nnn
      phot(i,j,i1,j1)=0.d0
      enddo
      enddo
      enddo
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
c     Now, read in the initial state population 
c     through the cross section file.
c     N, L, c(n,l)=cross sections, determine nmax
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    

      open(12,file='./init.dat', status='old')
          nmax=1
 2    read(12,*)n,l,c(n,l+1)
       if (n.lt.40) then
       if (n.gt.nmax) nmax=n
       goto 2
       endif
       
       read(12,*)q,rmax,energy,potential

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc          
c     file cevol.dat follows evolution of 
c     the atomic state contained in array c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c First all n,l,c is written to file, and the normalization. 
c cnew is the array of the population through each transition. 
       cnew(2,1)= c(2,1)
       cnew(1,1)= c(1,1)
       open(55,file='./cevol.dat')
c file cevolh.dat follows evolution of the atomic state contained in array c

       nmax0= nmax
       do ntest=1,nmax
       do ltest=0,ntest-1
       write(55,*)ntest,ltest, c(ntest,ltest+1)
       tnrm=tnrm+ c(ntest,ltest+1)
       enddo
       enddo
       write(55,*)'norm=', tnrm

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    
c       ---------start cascade-------------
c
c      nmax starts at the highest n and is subtracted by 1 for each cascade
c	        to lower energy levels.
c      c is modified after each cascade to the newest population.
c	        those values higher than nmax are left alone.
c      cnew is set to zero for recalculation up to nmax for each
c	        cascade
c      phot is the summation of the populations for each n and l level.
c         icasc is the counter for the cascade, the final value should 
c          be nmax-1.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

        icasc=0
c stable is the sum of the n=1 and 2 levels through each cascade.
        stable= c(2,1)+c(1,1)
		

 3      if (nmax.eq.1) goto 100
        if(nmax.eq.2.and.c(2,2).eq.0.d0) goto 100        
c nmax is reduced after each cascade since the transitions are moving
c to the lowest energy level.The cascade continues until only n=2-1 is left.
        icasc= icasc+1

       do n=1,nmax
       do l=0,n-1
       do np= 1,n-1
       do lp= 0, np-1
       cnew(np,lp+1)= cnew(np,lp+1) + prob1(np,lp,n,l)*c(n,l+1)
       phot(np,lp+1,n,l+1)= phot(np,lp+1,n,l+1)
     &  + prob1(np,lp,n,l)*c(n,l+1)
       enddo
       enddo
       enddo
       enddo

      
       nmax=nmax-1
c c is set equal to cnew and cnew is set to 0 after each cascade. 
c This only counts for n = 0 to nmax. 
       do n=1,nmax
       do l= 0,n-1
       c(n,l+1)= cnew(n,l+1) 
       cnew(n,l+1)=0.d0
       enddo
       enddo

       cnew(2,1)= c(2,1)
       cnew(1,1)= c(1,1)

       stable= stable+ c(2,1)+c(1,1)

       write(55,*) 'icasc=', icasc
       write(55,*) 'nmax=', nmax
c While nmax>1, write n,l,c (which is the current relative intensity) to file. 
c Write from n=0 to nmax. 
       if (nmax.gt.1) then
       tnrm=0.d0      
       do ntest=1,nmax
       do ltest=0,ntest-1
       write(55,*)ntest,ltest, c(ntest,ltest+1)
       tnrm=tnrm+ c(ntest,ltest+1)
       enddo
       enddo
       write(55,*)'norm=', tnrm
       else
       write(55,*)1,0, c(1,1)
       write(55,*)2,0, c(2,1)
       write(55,*)'norm=', c(1,1)+c(2,1)      
       endif

       goto 3
         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c            
c       ---------end cascade-------------
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Calculate Ratios
c    ephot is the energy difference of the initial and final state
c    ephot0 is the energy difference of the 2-1 state.
c    rad is equivalent to the population of a state. The relative 
c    intensity is equal to rad/phot
c    rad1 is summation of the l states for each n. 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 100   open(2,file='./phot.dat')
       open(3,file='./phot1.dat')
       open(7,file='./photout.dat')
       open(20,file='./ratios.dat')
       write(7,*) energy

       tot=0.d0

       do n=nmax0,2,-1
       do np=1,nmax0-1
       ephot= enr(n)-enr(np)
       ephot0= enr(2)-enr(1)
c Want to calculate Einstein A Coefficients here. 
       rad1=0.d0
       do lp= 0, np-1
       do l= 0,n-1
       rad= phot(np,lp+1,n,l+1)
	   
       if (rad.gt.0.d0) then 
       write(2,*)n,l,  '-->'  ,np,lp,'  phot(eV)=', ephot,
     &   '   ',rad,' relat= ',rad/phot(1,1,2,2)

       rad1= rad1+rad	 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	 

c pcs added to compute line ratios
       if((np.eq.1).and.(l.eq.1)) then
c           write(20,*)n,l, '-->'  ,np,lp,'  phot(eV)=', ephot,
c     .      '   ',rad,' relat= ',rad/phot(1,1,2,2)
           if(n.gt.2) sump=sump + rad/phot(1,1,2,2)
           if(n.eq.2) sump2=sump2 + rad/phot(1,1,2,2)
       end if
c           write(20,999) sump


       endif
       enddo
       enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Write resulting photons energies and Relative Intensities to file. 
       if (rad1.gt.0.d0)  then
          write(3,*)n, '-->'  ,np,'  phot(eV)=', ephot,
     &   '   ',rad1,' relat= ',rad1/phot(1,1,2,2)

          write(7,*)ephot, rad1/phot(1,1,2,2),rad1


       tot= tot+rad1
       endif
       enddo
	   enddo
       

   	   write(20,997) 'Energy(eV)', energy
   	   write(20,997) 'H ', sump/sump2
   	   write(20,*) '# H-like Hardness ratio as described in P. ',
     .   'Beiersdorfer, et al., Astrophys. J. 549, L147 (2001).'
       write(3,*)'totrad=', tot
       write(3,*)'totrad+stable=', tot+stable

       write(7,*)-1.d0, 0.d0
997    FORMAT(A14,F10.4)


       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c            Functions and Subroutines Below   
c
c            calculates energy levels in eV
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function enr(n)
c        calculates energy levels in eV
      implicit real*8(a-h,o-z)
      common rmax,q,potential
      enr = potential / n**2
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    calculates transition probabilities in terms of ordinary quant. num.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    

      function prob1(np,lp,n,l)
c    calculates transition probabilities in terms of ordinary quant. num.
      implicit real*8(a-h,o-z)
     
      prob1= 0.d0
      if  (np.ge.n.or.lp.eq.l) return
      if (lp-l.gt.1.or.l-lp.gt.1) return
      if (trans(np,lp,n,l).eq.0.d0) return
      anorm=0.d0
      do nn=1,n-1
      do ll= 0,nn-1 
      anorm= anorm + trans(nn,ll,n,l)
      enddo
      enddo

      prob1= trans(np,lp,n,l)/anorm
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   radiative probability  transition between hydrogen eigenstates in 1/sec
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function  trans(np,lp,n,l)
c   radiative probability  transition between hydrogen eigenstates in 1/sec
      implicit real*8 (a-h,o-z)
      dimension Y(100000),Z(100000)
      common rmax,q,potential
      trans=0.d0
      
      if(n.le.np.or.lp.eq.l) return
      if (lp.gt.l+1) return
      if (l.gt.lp+1) return

      ndim= 100000
      rmin=0.d0                

      do i=1,100000
         y(i)=0.d0      
         z(i)=0.d0
      enddo

      h= (rmax-rmin)/ndim

      do i=0,ndim-1
         r= rmin+i*h
         y(i+1)= hrad(n,l,r)*hrad(np,lp,r)*r**3
      enddo

      call  QSF(H,Y,Z,NDIM)
	  
      trans= z(ndim)/q
 

      lmax= l
      if (lmax.lt.lp) lmax=lp
         c= 1000.d0/7.2973d0
         trans= trans**2/(2*l+1)*lmax/6
         trans=trans*(1.d0/np**2-1.d0/n**2)**3/c**3
c         conversion from a.u. to 1/sec
         trans= trans*2*13.6058d16/6.5822d0
c       write(90,*) np,lp,n,l,q,trans
       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    radial eigenfunction of hydrogen atom in atomic units (Messiah & Sobelman)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       double precision function hrad(n,l,r)
c    radial eigenfunction of hydrogen atom in atomic units (Messiah & Sobelman)
       implicit real*8 (a-h,o-z)
       hrad=0.d0
       if (l.gt.n-1) return
       c= 1.d0
       do i=n-l,n+l
          c=c*i
       enddo

       c= dsqrt(.5d0/c/n)*(2.d0/n)**1.5d0
       x= r/n
       c= c*(2*x)**l*dexp(-x)
       hrad= c*funlag(n-l-1,2*l+1,2*x)
       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    lebedev's laguerre polynomial
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       double precision function funlag(n,m,x)
c    lebedev's laguerre polynomial
       implicit real*8 (a-h,o-z)
       funlag= 0.d0
       if (n.lt.0) return
       if (n.eq.0) then 
          funlag= 1.d0
          return
       endif
       if (n.eq.1) then
           funlag= 1.d0*(1+m)-x
           return
       endif    
       flag0=1.d0
       flag1=  1.d0*(1+m)-x
       k=1

  1    k=k+1
       funlag= (x-m*1.d0-2.d0*(k-1)-1.d0)*flag1
       funlag= funlag+ (k-1+m)*flag0
       funlag= -funlag/k
       flag0=flag1
       flag1=funlag             
       if(k.lt.n) goto 1
          return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE QSF(H,Y,Z,NDIM)
      implicit real*8 (a-h,o-z)
      dimension Y(*),Z(*)
      HT=H/3
      IF(NDIM-5)7,8,1
    1 SUM1=Y(2)+Y(2) 			     
      SUM1=SUM1+SUM1
      SUM1=HT*(Y(1)+SUM1+Y(3))
      AUX1=Y(4)+Y(4)
      AUX1=AUX1+AUX1
      AUX1=SUM1+HT*(Y(3)+AUX1+Y(5))
      aux2=ht*(y(1)+3.875d0*(y(2)+y(5))+2.625d0*(y(3)+y(4))+y(6))
      SUM2=Y(5)+Y(5)
      SUM2=SUM2+SUM2
      SUM2=AUX2-HT*(Y(4)+SUM2+Y(6))
      Z(1)=0.d0
      AUX=Y(3)+Y(3)
      AUX=AUX+AUX
      Z(2)=SUM2-HT*(Y(2)+AUX+Y(4))
      Z(3)=SUM1
      Z(4)=SUM2
      IF(NDIM-6)5,5,2
    2 DO 4 I=7,NDIM,2
      SUM1=AUX1
      SUM2=AUX2
      AUX1=Y(I-1)+Y(I-1)
      AUX1=AUX1+AUX1
      AUX1=SUM1+HT*(Y(I-2)+AUX1+Y(I))
      Z(I-2)=SUM1
      IF(I-NDIM)3,6,6
    3 AUX2=Y(I)+Y(I)
      AUX2=AUX2+AUX2
      AUX2=SUM2+HT*(Y(I-1)+AUX2+Y(I+1))
    4 Z(I-1)=SUM2                
    5 Z(NDIM-1)=AUX1
      Z(NDIM)=AUX2
      RETURN
    6 Z(NDIM-1)=SUM2
      Z(NDIM)=AUX1
      RETURN
    7 IF(NDIM-3)12,11,8
    8 SUM2=1.125d0*HT*(Y(1)+Y(2)+Y(2)+Y(2)+Y(3)+Y(3)+Y(3)+Y(4))
      SUM1=Y(2)+Y(2)
      SUM1=SUM1+SUM1
      SUM1=HT*(Y(1)+SUM1+Y(3))
      Z(1)=0.d0				 			
      AUX1=Y(3)+Y(3)
      AUX1=AUX1+AUX1
      Z(2)=SUM2-HT*(Y(2)+AUX1+Y(4))
      IF(NDIM-5)10,9,9
    9 AUX1=Y(4)+Y(4)
      AUX1=AUX1+AUX1
      Z(5)=SUM1+HT*(Y(3)+AUX1+Y(5))
   10 Z(3)=SUM1
      Z(4)=SUM2
      RETURN
   11 SUM1=HT*(1.25d0*Y(1)+Y(2)+Y(2)-.25d0*Y(3))
      SUM2=Y(2)+Y(2)
      SUM2=SUM2+SUM2
      Z(3)=HT*(Y(1)+SUM2+Y(3))
      Z(1)=0.d0
      Z(2)=SUM1
   12 RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


