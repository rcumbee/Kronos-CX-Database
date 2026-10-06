ccccccccccccccccccccccccccccccccccccccccccccccc
c gauss                                       c
ccccccccccccccccccccccccccccccccccccccccccccccc
c code to compute convolved spectrum assuming c
c gaussian line shape                         c
ccccccccccccccccccccccccccccccccccccccccccccccc
c Original version by Milun Rakovic July 2006 c
c 06-01-2012 PCS added some features and      c
c            documentation                    c
ccccccccccccccccccccccccccccccccccccccccccccccc
c input -  photouth.dat - output from         c
c       radcash.f                             c 
c input  - gauss.in  - control parameters     c
c output - exper.dat - convolved spectrum     c
c output - photline.dat - normalized line     c
c          spectrum                           c
ccccccccccccccccccccccccccccccccccccccccccccccc   
c
      implicit real*8(a-h,o-z)
          dimension ephot(6000), rad(6000) 
        double precision energy
        open(7,file='photout.dat',status='old')
        open(11,file='SourceCode/NonBare/gauss.in',status='old')
        open(12,file='photline.dat')

        read(7,*) energy
        write(12,996) 
        write(12,995) energy
        write(12,999)
        write(12,997)

        do i=1,6000
         ephot(i)=0.d0
         rad(i)=0.d0
        enddo

            read(11,*)delta,nn,emin,emax

           dpi= delta*dsqrt(dacos(-1.d0))


           i=0
  1      i=i+1
           read(7,*)ephot(i), rad(i)
           if (ephot(i).gt.0.d0) then
              write(12,*)ephot(i),rad(i)/dpi
              goto 1
           end if
    
           n=i-1


        open(8,file='spectrum.dat')
        write(8,998) delta
        write(8,995) energy
        write(8,999)
        write(8,997)

        do ii=0,nn
        e= emin+ ii*(emax-emin)/nn
         emis=0.d0
         do i=1,n
         emis= emis + rad(i)*dexp(-((e-ephot(i))/delta)**2)
         enddo
         emis= emis/dpi
         if (emis.gt.1.0d-99) write(8,*)e,emis
        enddo
998     format('# Convolved spectrum, gaussian FWHM = ',f8.2,' eV')
999     format('# Photon energy                  Intensity')
997     format('#     (eV)                       (arbitrary)') 
996     format('# Line spectrum')
995     format('# Collision energy = ',f8.3,' eV/u')
      end


