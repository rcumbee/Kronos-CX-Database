cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c radcashe                                                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    calculates cascade of radiative transitions for one-electron  c
c    excited state of helium-like ions                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Milun Rakovic  July 2006
c pcs 6-28-12 add Hardness, G, R ratio calculations and comments
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  input - inithe.dat  - relative cross sections for each nl level
c          obtained from, for example, casmainh
c  output - many files, but most important are:
c           phot* and ratios.dat
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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

c  read transition rates for n,l,i --> np,lp,ip (singl i,ip=0, or tripl i,ip=1) and the corresp. wavelengths 
      open(11,file='SourceCode/NonBare/radrathe.dat',status='old')
 1    read(11,*)n,l,i,np,lp,ip,
     &  alam(np,lp+1,ip+1,n,l+1,i+1),trans(np,lp+1,ip+1,n,l+1,i+1)
      if (n.lt.nl+1) goto 1
c    end read



c    read initial state population

      open(12,file='SourceCode/NonBare/inithe.dat', status='old')
      nmax=1
 2    read(12,*)n,l,i,c(n,l+1,i+1)
       if (n.lt.nl+1) then
       if (n.gt.nmax) nmax=n
       goto 2
       endif

       read(12,*)z,energy

c     end read

c   create probability matrix

c         write(*,*)'nmax=',nmax
     

  

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

       do n=1,nmax
       do l=0,n-1
       do i= 0,1

       aaa=0.d0
       do np=1,n
       do lp=0,np-1
       do ip= 0,1
       test=prob(np,lp+1,ip+1,n,l+1,i+1)
       aaa=aaa+test
       if (test.ne.0.d0) 
     &   write(33,*)np,lp,ip,n,l,i,prob(np,lp+1,ip+1,n,l+1,i+1)
       enddo
       enddo
       enddo
       write(33,*)'norm=',aaa
       write(33,*)'------------------------'

       enddo
       enddo
       enddo



c   end create matrix


      open(55,file='cevolhe.dat')

c file cevol.dat follows evolution of the atomic state contained in array c

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

c      start cascade

         stable=0.d0

        icasc=0

 3     if(nmax.eq.2.and.c(2,2,1).eq.0.d0.and.c(2,2,2).eq.0.d0) goto 100        
       
        icasc= icasc+1

        if(icasc.gt.20) goto 100

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
     &   + prob(np,lp+1,ip+1,n,l+1,i+1)*c(n,l+1,i+1)
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
         
c       end cascade

 100   open(2,file='phothe.dat')
       open(7,file='photout.dat')
       open(20,file='ratios.dat')
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
cpcs added to compute line ratios
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
         write(7,*)ephot, rad/phot0

        tot= tot+rad
        endif
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo

        write(20,999)
        write(20,998)
        write(20,*) energy, res, forb, interc, forb/interc, res/interc, 
     .              (forb+interc)/res,sump/sump2
        write(2,*)'totrad=', tot
        write(2,*)'totrad+stable=', tot+stable

         write(7,*)-1.d0, 0.d0
999    format('# Collision        He-like line intensities and ratios')
998    format('# energy (eV/u)    Res    Forb   Interc   R   L   G  H')

       end


