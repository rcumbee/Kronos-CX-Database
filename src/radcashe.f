cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c radcashe                                                         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    calculates cascade of radiative transitions for one-electron  
c    excited state of helium-like ions
c                         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Copyright 2017 P. D. Mullen, R. S. Cumbee, P. C. Stancil
c
c                        Kronos version 3.0   (released Feb 7, 2017)
c
c    This file is part of Kronos. Kronos is a series of python and fortran codes 
c    utilizing a comprehensive charge exchange database at its backbone to generate 
c    theoretical line ratios and spectra for many single electron capture systems. 
c
c    You are morally obligated to cite the following paper in his/her
c    scientific literature that results from use of any part of HAM:
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
c    pcs 6-28-12 add Hardness, G, R ratio calculations and comments
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  input - init.dat  - relative cross sections for each nlS level
c  output - many files, but most important are:
c           phot* and ratios.dat
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit real*8(a-h,o-z)
      real sump, hard, interc, forb, res, sump2, energy
      parameter (nl=15)
      dimension trans(nl,nl,2,nl,nl,2),alam(nl,nl,2,nl,nl,2)
      dimension c(nl,nl,2), cnew(nl,nl,2)
      dimension prob(nl,nl,2,nl,nl,2)
      dimension phot(nl,nl,2,nl,nl,2)
      
      sump=0.0d0
      sump2=0.0d0
      hard=0.0d0

      do i=1,nl
      do j=1,nl
      do k=1,2
      do i1=1,nl
      do j1=1,nl
      do k1=1,2
      trans(i,j,k,i1,j1,k1)=0.d0
      prob(i,j,k,i1,j1,k1)=0.d0
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo


      do i=1,nl
      do j=1,nl
      do k=1,2
      c(i,j,k)=0.d0
      cnew(i,j,k)=0.d0
      do i1=1,nl
      do j1=1,nl
      do k1=1,2
      alam(i,j,k,i1,j1,k1)=0.d0
      phot(i,j,k,i1,j1,k1)=0.d0
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  
c  read transition rates for n,l,i --> np,lp,ip (singl i,ip=0, or tripl i,ip=1) 
c  and the corresp. wavelengths 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

      open(11,file='./radrathe.dat',status='old')
 1    read(11,*)n,l,i,np,lp,ip,
     &  alam(np,lp+1,ip+1,n,l+1,i+1),trans(np,lp+1,ip+1,n,l+1,i+1)
      if (n.lt.nl+1) goto 1
c    end read

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  
c    Now, read in the initial state population through the cross section file.
c    N, L, c(n,l)=cross sections, determine nmax
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

      open(12,file='./init.dat', status='old')
      nmax=1
 2    read(12,*)n,l,i,c(n,l+1,i+1)
      if (n.lt.nl+1) then
      if (n.gt.nmax) nmax=n
      goto 2
      endif

      read(12,*)energy

c     end read

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c             create probability matrix
c	  
c 	Here we calculate the branching ratios.	  
c     Prob matrix is the matrix of branching ratios.
c     anorm is the normalized Einstein A coefficient for
c        each upper nls level
c     the prob matrix is the the Einstein A coefficient
c        for each transition divided by the anorm for each
c        upper level--- AKA the branching ratio.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	

       do n=1,nmax
       do l=0,n-1
       do i= 0,1

       anorm= 0.d0
       do np=1,n
       do lp=0,np-1
       do ip= 0,1
       anorm=anorm+trans(np,lp+1,ip+1,n,l+1,i+1)
       enddo
       enddo
       enddo

       if (anorm.gt.0.d0) then
       do ip= 0,1
       do np=1,n
       do lp=0,np-1
       prob(np,lp+1,ip+1,n,l+1,i+1)
     &                =trans(np,lp+1,ip+1,n,l+1,i+1)/anorm
c234567
       enddo
       enddo
       enddo
       else
       prob(n,l+1,i+1,n,l+1,i+1)=1.d0
       endif

       enddo
       enddo
       enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	  
c 	End calculating brancing ratios
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       do n=1,nmax
       do l=0,n-1
       do i= 0,1

       aaa=0.d0
       do np=1,n
       do lp=0,np-1
       do ip= 0,1
       test=prob(np,lp+1,ip+1,n,l+1,i+1)
       aaa=aaa+test
c       if (test.ne.0.d0) 
c     &   write(33,*)np,lp,ip,n,l,i,prob(np,lp+1,ip+1,n,l+1,i+1)
       enddo
       enddo
       enddo
c       write(33,*)'norm=',aaa
c       write(33,*)'------------------------'

       enddo
       enddo
       enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   end create matrix
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c   
c   file cevol.dat follows evolution of the atomic state contained in array c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

c First all n,l,c is written to file, and the normalization. 
c cnew is the array of the population through each transition. 

       open(55,file='./cevol.dat')
       nmax0= nmax
       tnrm=0.d0      
       do i=0,1
       do ntest=1,nmax
       do ltest=0,ntest-1
       test= c(ntest,ltest+1,i+1)
       if (test.ne.0.d0) write(55,*)ntest,ltest,i,test
       tnrm=tnrm+ c(ntest,ltest+1,i+1)
       enddo
       enddo
       enddo
       write(55,*)'norm=', tnrm

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c
c          --------start cascade----------
c
c    nmax starts at the highest n and is subtracted by 1 for each cascade
c	     to lower energy levels.
c    c is modified after each cascade to the newest population.
c	     those values higher than nmax are left alone.
c    cnew is set to zero for recalculation up to nmax for each
c	     cascade
c    phot is the summation of the populations for each n and l level.
c    icasc is the counter for the cascade, the final value should be nmax-1.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

       stable=0.d0

       icasc=0

 3     if(nmax.eq.2.and.c(2,2,1).eq.0.d0.and.c(2,2,2).eq.0.d0) goto 100        
       
       icasc= icasc+1

       if(icasc.gt.50) goto 100

       do n=1,nmax
       do l=0,n-1
       do i=0,1
       if (prob(n,l+1,i+1,n,l+1,i+1).lt.1.d0) then
       do np= 1,n
       do lp= 0, np-1
       do ip= 0,1
       cnew(np,lp+1,ip+1)= cnew(np,lp+1,ip+1)
     &  + prob(np,lp+1,ip+1,n,l+1,i+1)*c(n,l+1,i+1)
       if((n-np)**2+(l-lp)**2+(i-ip)**2.gt.0) then
       phot(np,lp+1,ip+1,n,l+1,i+1)= phot(np,lp+1,ip+1,n,l+1,i+1)
     &  + prob(np,lp+1,ip+1,n,l+1,i+1)*c(n,l+1,i+1)
       endif
       enddo
       enddo
       enddo
       else
       cnew(n,l+1,i+1)= c(n,l+1,i+1)
       stable=stable+c(n,l+1,i+1)
       endif
       enddo
       enddo
       enddo       

       do n=1,nmax
       do l= 0,n-1
       do i=0,1
       c(n,l+1,i+1)= cnew(n,l+1,i+1) 
       cnew(n,l+1,i+1)=0.d0
       enddo
       enddo
       enddo

       nmax1=nmax
       nmax=1
      
       do n=1,nmax1 
       do l=0,n-1
       do i=0,1
       if(c(n,l+1,i+1).gt.0.d0.and.n.gt.nmax) nmax=n
       enddo
       enddo
       enddo

       write(55,*) 'icasc=', icasc
       write(55,*) 'nmax=', nmax

       tnrm=0.d0      
       do ntest=1,nmax
       do ltest=0,ntest-1
       do itest=0,1
       test= c(ntest,ltest+1,itest+1)
       if (test.ne.0.d0) write(55,*)ntest,ltest,itest,test
       tnrm=tnrm+ c(ntest,ltest+1,itest+1)
       enddo
       enddo
       enddo
       write(55,*)'norm=', tnrm

       goto 3
         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c          
c          --------end cascade----------
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Calculate Ratios
c ephot is the energy difference of the initial and final state
c ephot0 is the energy difference of the 2-1 state.
c rad is equivalent to the population of a state. The relative intensity
c is equal to rad/phot
c rad1 is summation of the l states for each n. 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 100   open(2,file='./phot.dat')
       open(7,file='./photout.dat')
       open(27,file='./ratios.dat')
       write(7,*) energy
	   

       tot=0.d0
       phot0= phot(1,1,1,2,2,1)

       do n=nmax0,2,-1
       do np=1,n
       do lp= 0, np-1
       do l= 0,n-1
       do i=0,1
       do ip=0,1
       rad= phot(np,lp+1,ip+1,n,l+1,i+1)
       ephot= 1.d8/alam(np,lp+1,ip+1,n,l+1,i+1)/8065.5d0

       if (rad.gt.0.d0) then 
       write(2,*)n,l,i , '-->'  ,np,lp,ip,'  phot(eV)=', ephot,
     &   '   ',rad,' relat= ',rad/phot0
c
c pcs added to compute line ratios
c        if((np.eq.1).and.(i.eq.0).and.(l.eq.1)) then
       if((np.eq.1)) then
c           write(20,*)n,l,i , '-->'  ,np,lp,ip,'  phot(eV)=', ephot,
c     .      '   ',rad,' relat= ',rad/phot0
       if(n.gt.2) sump=sump + rad/phot0
       if(n.eq.2) sump2=sump2 + rad/phot0
c           write(20,*) sump,sump2
       end if 
       if((np.eq.1).and.(n.eq.2)) then
c           write(20,*)n,l,i , '-->'  ,np,lp,ip,'  phot(eV)=', ephot,
c     .      '   ',rad,' relat= ',rad/phot0
       if((l.eq.1).and.(ip.eq.0).and.(i.eq.0)) res=rad/phot0 
       if((l.eq.0).and.(ip.eq.0).and.(i.eq.1)) forb=rad/phot0 
       if((l.eq.1).and.(ip.eq.0).and.(i.eq.1)) interc=rad/phot0 
       endif
       write(7,*)ephot, rad/phot0, rad

       tot= tot+rad
       endif
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo

c       write(20,999)
c       write(20,998)
c	   write(20,*) energy, res, forb, interc, forb/interc, res/interc, 
c     .              (forb+interc)/res,sump/sump2
   	   write(27,997) 'Energy(eV)', energy
   	   write(27,997) 'H ', sump/(res + forb + interc)
   	   write(27,997) "H'", sump/res
   	   write(27,997) "R ", forb/interc
   	   write(27,997) 'G ', (forb + interc)/res
   	   write(27,997) "G'", (forb + interc)/(res + sump)
   	   write(27,*) '# He-like Line ratios are provided as described',
     .   ' in Leutenegger, M. A., et al. PRL 105.6 (2010): 063201.'
   	   write(27,*) '# H = K beta+ / SUM( K alpha F + I + R)'
   	   write(27,*) "# H' = K beta+ / K alpha R"
   	   write(27,*) "# R = K alpha F / K alpha I"
   	   write(27,*) "# G = SUM(K alpha F + I) / K alpha R"
   	   write(27,*) "# G' = SUM(K alpha F + I)/SUM( K alpha R + ",
     .              "K beta+)"
	   
	   
       write(7,*)-1.d0, 0.d0
999    format('# Collision        He-like line intensities and ratios')
998    format('# energy (eV/u)    Res    Forb   Interc   R   L   G  H')
997    FORMAT(A14,F10.4)

       end


