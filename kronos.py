#! /usr/bin/env python

##################################################################################
#
#    Copyright 2017 P. D. Mullen, R. S. Cumbee, P. C. Stancil
#
#                        Kronos version 3.1   (released May 16, 2017)
#
#    This file is part of Kronos. Kronos is a series of python and fortran codes 
#    utilizing a comprehensive charge exchange database at its backbone to generate 
#    theoretical line ratios and spectra for many single electron capture systems. 
#
#    You are morally obligated to cite the following paper in his/her
#    scientific literature that results from use of any part of Kronos:
#
#    [1] Mullen et al. (2016), ApJS 224, 31
#    [2] Cumbee et al. (2017), ApJ, in prep
#    [3] Lyons et al. (2017), ApJS, in prep
#    [4] Mullen et al. (2017), ApJ, in prep
#        
#
#    Kronos is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.
#
#    Kronos is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##################################################################################

#Add source directory of cx_functions and import
import sys
sys.path.append('./src/')
from cx_functions import *

#Main Function
def main():
	getFiles() #Gets the CSFILE relevant to query
	cascade() #Runs the Cascade Model for X-ray Emission
	getResults() #Gets results by interacting with output Cascade Files.
	cleanUp() #Move files to results directory and empty working directory
	end() #Close Program

#Call main function
if __name__ == "__main__":
	main()
	