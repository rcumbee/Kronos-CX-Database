# Developed by: Patrick Mullen

#####################################################################################

What is Kronos?


Kronos is a series of python and fortran codes utilizing a comprehensive 
charge exchange database at its backbone to generate theoretical line ratios 
and spectra for many single electron capture systems.  Kronos’ power can be 
found in that it is 1) very easy to use, 2) it allows for the user to specify 
a variety of input parameters of interest 3) it is quick, and 4) above all,  
it contains a fully accessible and comprehensive database of multi-channel 
Landau-Zener single electron capture nlS-resolved cross sections for many 
projectile ion/ neutral target pairs at a wide variety of collision energies.  
These cross sections can be used in other Cascade models outside of Kronos or 
studied directly — in other words, for some users, Kronos may never be run 
and it serves the sole purpose of acting as a charge exchange cross section database.  

#####################################################################################


What is contained within the Kronos distribution?

I. CXDatabase:
	1. excitation energy files (*.nist.in)
	2. n-resolved cross section files (*_nres.cs)
	3. nlS-resolved cross section files
	4. Bare Ion Projectiles
		Low Energy Distribution (*_lowe.cs)
		Statistical Distribution (*_stat.cs)
	5. Non-Bare Ion Projectiles
		No Distribution Function (*_mclz.cs)
	6. Cascade Files
		Transitions 
			nlS initial—> nlS final
		Photon Energies
	7.Transition Probabilities (Einstein A Coefficients)
II. kronos.py 
	-the main script to be run
III. README.txt
	-Contains this information.
IV. Results
	-Empty until you run kronos.py
	-After Running kronos.py:
		1. Line Ratios ASCII File (LineRatios.dat)
		2. Line Ratios Tex File (LineRatios.tex)
		3. Line Ratios PDF File (LineRatios.pdf)
		4. Spectrum File (spectrum.dat)
		5. Miscellaneous
V. Source Code
	-Contains other files/codes called upon by kronos.py that are mandatory for proper running.  

#####################################################################################

How do I run Kronos?

-It is worth noting that many people will not need to run Kronos; to some, 
it may only be a source of charge exchange data (in CXDatabase directory).  
Otherwise, see below:

-You MUST have Pylab installed on your machine.  I recommend obtaining Pylab 
from http://continuum.io/downloads (Python 2.7).  You must also have gfortran compilers.  

1. To run, in working directory simply type:

./kronos.py -- help

2. This brings up a help menu.  It details how to properly run the program, i.e. we must 
add tags/command-line-options (-I, -Z, -A, -E, etc.).  

3. An example of running kronos.py for C^{5+} + H2O CX collisions at 10 eV/u mimicking a detector resolution of 10 eV:

	./kronos.py -I C -Z 5 -A H2O -E 10 -R 10

4. Some tags/command-line-options are required (example -I, -Z, -A, -E, etc.) while others 
are not (-G, -R, -X, -N, etc.)
	-This is also explained in ./kronos.py --help

5. After running program, view results in the Results directory.

#####################################################################################

