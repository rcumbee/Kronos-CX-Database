cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c gauss                                       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      code to compute convolved spectrum 
c      assuming gaussian line shape
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
c    Orig. version by Milun Rakovic July 2006 
c    06-01-2012 PCS added some features and      
c       documentation  
c                  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   input -  photout.dat - output from          
c   input  - gauss.in  - control parameters     
c   output - spectrum.dat - convolved spectrum  
c   output - photline.dat - normalized line     
c          spectrum                           
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   

      implicit real*8(a-h,o-z)
          dimension ephot(6000), rad(6000) 
        double precision energy
        open(7,file='./photout.dat',status='old')
        open(11,file='./gauss.in',status='old')
        open(12,file='./photline.dat')

        read(7,*) energy
        write(12,996) 
        write(12,995) energy
        write(12,999)
        write(12,997)

        do i=1,1000
         ephot(i)=0.d0
         rad(i)=0.d0
        enddo

            read(11,*)delta,nn,emin,emax

           dpi= delta*dsqrt(dacos(-1.d0))


           i=0
  1      i=i+1
           read(7,*)ephot(i), rad(i)
           if (ephot(i).gt.0.d0) then
              write(12,*)ephot(i),rad(i)
              goto 1
           end if
    
           n=i-1


        open(8,file='./spectrum.dat')
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
         if (emis.gt.1.0d-99) write(8,*)e,emis
        enddo
998     format('# Convolved spectrum, gaussian FWHM = ',f8.2,' eV')
999     format('# Photon energy                  Intensity')
997     format('#     (eV)                       (arbitrary)') 
996     format('# Line spectrum')
995     format('# Collision energy = ',f8.3,' eV/u')
      end


