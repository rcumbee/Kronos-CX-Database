cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c radcash                                                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    calculates cascade of radiative transitions for one-electron  c
c    excited state of the hydrogenic atom                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Milun Rakovic  July 2006
c pcs 5-28-12 add Hardness ratio calculation and comments 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  input - inith.dat  - relative cross sections for each nl level
c          obtained from, for example, casmainh
c  output - many files, but most important are:
c           phot* and ratios.dat
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c   Initializations
c   nnn is the maximum possible n value in this code.
c   rmax - maximum radius of radial integrations
c   z is the charge 
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      implicit real*8(a-h,o-z)
      real sump, sump2, energy
      parameter (nnn=40)
      dimension c(nnn,nnn), cnew(nnn,nnn)
      dimension phot(nnn,nnn,nnn,nnn)
      common rmax,z,potential
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
c     initialize three arrays, c, cnew, and phot to be 0 to max possible n. 
c-----------------------------------------  
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
c-----------------------------------------  
c Now, read in the initial state population through the cross section file.
c N, L, c(n,l)=cross sections, determine nmax
c-----------------------------------------  
      open(12,file='SourceCode/Bare/inith.dat', status='old')
          nmax=1
 2    read(12,*)n,l,c(n,l+1)
       if (n.lt.40) then
       if (n.gt.nmax) nmax=n
       goto 2
       endif
       
       read(12,*)z,rmax,energy,potential
c-----------------------------------------        
c file cevolh.dat follows evolution of the atomic state contained in array c
c-----------------------------------------  
c First all n,l,c is written to file, and the normalization. 
c cnew is the array of the population through each transition. 
       cnew(2,1)= c(2,1)
       cnew(1,1)= c(1,1)

      open(55,file='cevolh.dat')

c file cevolh.dat follows evolution of the atomic state contained in array c

       nmax0= nmax
       tnrm=0.d0      
       do ntest=1,nmax
       do ltest=0,ntest-1
       write(55,*)ntest,ltest, c(ntest,ltest+1)
       tnrm=tnrm+ c(ntest,ltest+1)
       enddo
       enddo
       write(55,*)'norm=', tnrm

c-----------------------------------------  
c~~~~~~start cascade~~~~~~~~~~~~~~~~~~~~~~
c      nmax starts at the highest n and is subtracted by 1 for each cascade
c	 to lower energy levels.
c      c is modified after each cascade to the newest population.
c	  those values higher than nmax are left alone.
c      cnew is set to zero for recalculation up to nmax for each
c	  cascade
c      phot is the summation of the populations for each n and l level.
c  icasc is the counter for the cascade, the final value should be nmax-1.
c-----------------------------------------  
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
         
c-----------------------------------------           
c       end cascade
c-----------------------------------------

c-----------------------------------------
C Calculate Ratios
c ephot is the energy difference of the initial and final state
c ephot0 is the energy difference of the 2-1 state.
c rad is equivalent to the population of a state. The relative intensity
c is equal to rad/phot
c rad1 is summation of the l states for each n. 

c-----------------------------------------
 100   open(2,file='photh.dat')
       open(3,file='phot1h.dat')
       open(7,file='photouth.dat')
       open(20,file='ratios.dat')
       write(7,*) energy

       tot=0.d0

       do n=nmax0,2,-1
       do np=1,nmax0-1
       ephot= enr(n,z)-enr(np,z)
       ephot0= enr(2,10.d0)-enr(1,10.d0)
c Want to calculate Einstein A Coefficients here. 
       rad1=0.d0
       do lp= 0, np-1
       do l= 0,n-1
       rad= phot(np,lp+1,n,l+1)
       if (rad.gt.0.d0) then 
       write(2,*)n,l,  '-->'  ,np,lp,'  phot(eV)=', ephot,
     &   '   ',rad,' relat= ',rad/phot(1,1,2,2)
c23456
       write(25,*)n,l,  '-->'  ,np,lp,'  phot(eV)=', ephot,
     &   '   ',rad,' relat= ',rad/phot(1,1,2,2)*ephot/ephot0

       rad1= rad1+rad
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	   
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
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c Write resulting photons energies and Relative Intensities to file. 
       if (rad1.gt.0.d0)  then
          write(3,*)n, '-->'  ,np,'  phot(eV)=', ephot,
     &   '   ',rad1,' relat= ',rad1/phot(1,1,2,2)
c23456
          write(35,*)n, '-->'  ,np,'  phot(eV)=', ephot,
     &   '   ',rad1,' relat= ',rad1/phot(1,1,2,2)*ephot/ephot0
          write(7,*)ephot, rad1/phot(1,1,2,2),rad1


       tot= tot+rad1
       endif
       enddo
	   enddo
       

       write(20,999)
       write(20,*) energy, sump/sump2
       write(3,*)'totrad=', tot
       write(3,*)'totrad+stable=', tot+stable

       write(7,*)-1.d0, 0.d0
999    format('# Collision energy (eV/u), Hardness ratio')

       end
c###################################################################
c########## Functions and Subroutines Below   ######################
c###################################################################
cc**************************************************************************************
c        calculates energy levels in eV
c**************************************************************************************
      function enr(n,z)
c        calculates energy levels in eV
      implicit real*8(a-h,o-z)
      common rmax,q,potential
      enr = potential / n**2
      return
      end
cc************************************************************************************
c    calculates transition probabilities in terms of ordinary quant. num.
c************************************************************************************    
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
c######################################################
c************************************************************************************
c   radiative probability  transition between hydrogen eigenstates in 1/sec
c************************************************************************************
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
       write(90,*) np,lp,n,l,q,trans
       return
       end
c######################################################
c************************************************************************************
c    radial eigenfunction of hydrogen atom in atomic units (Messiah & Sobelman)
c************************************************************************************
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
c######################################################
c************************************************************************************
c    lebedev's laguerre polynomial
c************************************************************************************
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
c######################################################
c**************************************************************************************
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
c######################################################
c######################################################


