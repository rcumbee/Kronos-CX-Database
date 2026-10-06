#! /usr/bin/env python
# Developed by: Patrick Mullen

##################################################################################
#
#    Copyright 2017 P. D. Mullen, R. S. Cumbee, P. C. Stancil
#
#                        Kronos version 3.0   (released Feb 7, 2017)
#
#    This file is part of Kronos. Kronos is a series of python and fortran codes 
#    utilizing a comprehensive charge exchange database at its backbone to generate 
#    theoretical line ratios and spectra for many single electron capture systems. 
#
#    You are morally obligated to cite the following paper in his/her
#    scientific literature that results from use of any part of HAM:
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

from numpy import loadtxt
import argparse,os,sys,math,subprocess,shutil,time, matplotlib.pyplot

#Gets the CSFILE relevant to query
def getFiles():
	parseCommandLineOptions()
	problemGenerator()
	availabilityCheck()

#Runs the Cascade Model for X-ray Emission
def cascade():
	prepCascade()
	doCascade()
	makeSpectrum()

#Gets results by interacting with output Cascade Files.  
#Plots if set by -G parser
def getResults():
	lineRatios()
	plotSpectrum()

#Move files to results directory and empty working directory
def cleanUp():
	moveFiles()
	deleteBin()

#Inputs information specified in command line options
def parseCommandLineOptions():
	#Declare Globals
	global projectile, charge, target, collEnergy, distro
	global resolution, rmax, methodCheck, graphCheck, numberOfPoints

	#Parser Information
	parser=argparse.ArgumentParser()
	parser.add_argument("-I", "--projectile", help="Specify Projectile Ion -- Required", required=True)
	parser.add_argument("-Z", "--charge", help="Specify Projectile Charge -- Required", required=True)
	parser.add_argument("-A", "--target", help="Specify Target Neutral -- Required", required=True)
	parser.add_argument("-E", "--collisionEnergy", help="Specify Collision Energy in eV/u -- Required", \
		required=True)
	parser.add_argument("-D", "--distribution", help="Specify distribution function -- \
		Only Required for n-resolved data (Options: stat, lowe, lowemod, separable, SL1)", required=False)
	parser.add_argument("-R", "--resolution", help="Specify the FWHM Resolution in eV -- Optional", \
		required=False)
	parser.add_argument("-P", "--numberOfPoints", help="Specify number of data points -- Optional", \
		required=False)
	parser.add_argument("-X", "--rmax", help= "Specify rmax -- Optional (not recommended)", required=False)
	parser.add_argument("-M", "--method", help="Specify method of Charge Exchange Data Calculation \
		-- Optional (Options: MCLZ, CTMC, AOCC, MOCC, QMOCC, RCMD...check Data Availability)", required=False)
	parser.add_argument("-G", "--graph", help="Specify \"yes\" or \"no\" to graph -- Optional", required=False)
	args=parser.parse_args()

	#Required Command Line Options
	projectile = str(args.projectile)
	charge = int(args.charge)
	target = str(args.target)		
	collEnergy = float(args.collisionEnergy)
	distro = args.distribution #only required if bare ion collision

	#Optional Command Line Options
	resolution = 10.0 if args.resolution == None else float(args.resolution)
	numberOfPoints = 4200 if args.numberOfPoints == None else int(args.numberOfPoints)
	rmax = 200.0 if args.rmax == None else float(args.rmax)
	methodCheck = args.method #set to None if not specified
	graphCheck = args.graph   #set to None if not specified

def problemGenerator():
	#Macros
	global hashLine, newLine
	hashLine = "###################################################\n" #line of hashes for organized output
	newLine  = "\n" #newLine

	#Print Problem To User
	print newLine, hashLine, newLine, "Collision System:", newLine
	print projectile + "^{" + str(charge) + "+}" + " + " + target + " --> " + projectile + "^{" + \
		str(charge-1) + "+}" + " + " + target + "^{+}" + newLine
	print "Projectile Ion: ", projectile
	print "Projectile Ion Charge: ", str(charge)
	print "Target Neutral: ", target
	print "Collision Energy Requested: ", str(collEnergy), "eV/u"
	print "FWHM Resolution: ", str(resolution), "eV"
	print newLine

	#Read in elements.dat
	global atomicZ
	AtomicNumbers, Elements = loadtxt("CXDatabase/elements.dat", dtype = str, unpack = True, usecols=[0,1])
	Elements = Elements.tolist() #indexing requires list format
	elementInFile = projectile in Elements #Make sure specified element exists in elements.dat

	#Determine if bare or non-bare collision
	global collisionType
	if elementInFile == True: 
		projectileIndex = Elements.index(projectile)  #find the index of projectile in elements.dat
		atomicZ = int(AtomicNumbers[projectileIndex]) #find the atomic number associated with projectile

		if atomicZ == charge:
			print "...This is a bare ion collision.\n"
			print hashLine
			collisionType = "bare" #set collision type to bare
		else:
			print "...This is a non-bare ion collision.\n"
			print hashLine
			collisionType = "non-bare" #set collision type to non-bare

	else: 
		print hashLine, newLine, "Error: Cannot Find Projectile Type in ./CXDatabase/elements.dat\n", \
			newLine, hashLine
		sys.exit("Exiting Unsuccessfully.\n" + newLine + hashLine)

def availabilityCheck():
	#Is data file for collision in CXDatabase?
	global directoryTree
	directoryTree=[os.path.join(dp, f) \
		for dp, dn, fn in os.walk(os.path.expanduser("./CXDatabase")) for f in fn]

	if collisionType == "bare":
		findFile = projectile.lower()+str(charge)+"+"+target.lower()+"_"
		FoundFiles = [i for i in directoryTree if findFile in i and ".cs" in i]
		if len(FoundFiles)==0:
			print "Error: Cannot Find Data for this Collision System.\n", newLine, hashLine
			sys.exit("Exiting Unsuccessfully.\n" + newLine + hashLine)

	if collisionType == "non-bare":
		findFile = projectile.lower()+str(charge)+"+"+target.lower()+"_"
		FoundFiles = [i for i in directoryTree if findFile in i and ".cs" in i and "nres" not in i]
		if len(FoundFiles)==0:
			print "Error: Cannot Find Data for this Collision System.\n", newLine, hashLine
			sys.exit("Exiting Unsuccessfully.\n" + newLine + hashLine)

	#Which methods are avaible? 
	global method
	method=None
	checkMCLZ = "Yes" if any("mclz" in i for i in FoundFiles) else "No"
	checkCTMC = "Yes" if any("ctmc" in i for i in FoundFiles) else "No"
	checkAOCC = "Yes" if any("aocc" in i for i in FoundFiles) else "No"
	checkQMOCC= "Yes" if any("qmocc" in i for i in FoundFiles) else "No"
	checkRCMD = "Yes" if any("rcmd" in i for i in FoundFiles) else "No"
	for i in range(len(FoundFiles)): #Special attention here to discern QMOCC and MOCC
		if "mocc" in FoundFiles[i]:
			if "qmocc" not in FoundFiles[i]: checkMOCC="Yes"
		else: checkMOCC="No" 

	print "Data Availability:"
	print "MCLZ:", checkMCLZ
	print "CTMC:", checkCTMC
	print "AOCC:", checkAOCC
	print "MOCC:", checkMOCC
	print "QMOCC:", checkQMOCC
	print "RCMD:", checkRCMD
	print newLine

	#Set hierarchy for method usage if not manually set. 
	#Hierarchy (starting with most preferred): RCMD -> QMOCC -> MOCC -> AOCC -> CTMC -> MCLZ
	if checkMCLZ == "Yes": method = "mclz"
	if checkCTMC == "Yes": method = "ctmc"
	if checkAOCC == "Yes": method = "aocc"
	if checkMOCC == "Yes": method = "mocc" 
	if checkQMOCC == "Yes": method = "qmocc"
	if checkRCMD == "Yes": method = "rcmd"


	#Check if method was manually set
	if methodCheck != None:
		method=methodCheck.lower()
		if not any(method in i for i in FoundFiles) or (method=="mocc" and checkMOCC=="No"): 
			print "Error: No", method.upper(), "Data for this System."
			sys.exit("Exiting Unsuccessfully.\n" + newLine + hashLine)
	print "...Proceeding with", method.upper(), "Cross Sections."
	print newLine, hashLine

	#Check to see if distribution function assigned for bare and not assigned for non-bare
	if collisionType == "bare" and distro==None and method=="mclz": 
		print "Must specify distribution function for bare collisions."
		sys.exit("Exiting Unsuccessfully.\n" + newLine + hashLine)
	if collisionType == "bare" and distro!=None and method!="mclz": 
		print "\n...Do not specifiy distribution function for bare collisions using non-MCLZ data."
		print "...Ignoring -D " + str(distro) + newLine
	if collisionType == "non-bare" and distro!=None:
		print "\n...Do not specifiy distribution function for non-bare collisions."
		print "...Ignoring -D " + str(distro) + newLine

	#CS File to be used for remainder of code
	global CSFILE
	if method!="mocc":
		CSFileCheck = [i for i in FoundFiles if method in i]
	else:
		CSFileCheck = [i for i in FoundFiles if "mocc" in i and "qmocc" not in i]
	if len(CSFileCheck) > 1:
		print "Error: Conflicting Files in CXDatabase for this collision system and method:"
		for i in range(len(CSFileCheck)): print CSFileCheck[i]
		sys.exit("\nExiting Unsuccessfully.\n" + newLine + hashLine)
	CSFILE = CSFileCheck[0]

	#NIST File to be used for remainder of code
	global NISTFILE 
	NistFileCheck = [i for i in os.listdir(os.path.dirname(CSFILE)) if "nist.in" in i]
	if len(NistFileCheck) > 1:
		print "Error: Conflicting Files in CXDatabase for this collision system and method:"
		for i in range(len(NistFileCheck)): print NistFileCheck[i]
		sys.exit("\nExiting Unsuccessfully.\n" + newLine + hashLine)
	if len(NistFileCheck) == 0:
		print "Error: No NIST File for this Collision.  Add NIST file to specific collision directory."
		sys.exit("\nExiting Unsuccessfully.\n" + newLine + hashLine)
	NISTFILE = os.path.join(os.path.dirname(CSFILE), NistFileCheck[0])

	#Radrathe File to be used for remainder of code if non-bare ion collision
	global RADRATHE
	RADRATHE=None
	if collisionType == "non-bare":
		RadratheFileCheck = [i for i in os.listdir(os.path.dirname(CSFILE)) if "radrathe.dat" in i]
		if len(RadratheFileCheck) == 0:
			print str("Error: No radrathe file for this Collision. Add radrathe file to specific" +\
				" collision directory.")
			sys.exit("\nExiting Unsuccessfully.\n" + newLine + hashLine)
		RADRATHE = os.path.join(os.path.dirname(CSFILE), RadratheFileCheck[0])

def prepCascade():
	#Find closest collision energy in CSFILE
	global closestCollEnergy
	FileCollEnergies = loadtxt(CSFILE, dtype=float, unpack = True, usecols=[0])
	FileExcitEnergies = loadtxt(NISTFILE, dtype=float, unpack=True, usecols=[-1])
	closestCollEnergy = min(enumerate(FileCollEnergies), key=lambda x: abs(x[1]-collEnergy))
	ionizationEnergy = FileExcitEnergies[-1] #ionization energy from NIST File
	print "...Proceeding with Collision Energy", str(closestCollEnergy[1]), "eV/u.\n"
	if abs(closestCollEnergy[1]-collEnergy) > 500.0:
		print "...Caution: This differs by", str(abs(closestCollEnergy[1]-collEnergy)), \
			"eV/u than collision energy requested. Perhaps specifiy different method.\n"

	#Find CS data associated with closest collision energy
	f=open(CSFILE, "r")
	for i in range(7): f.readline() #skip first 7 lines of File (Assuming MCLZ header structure)
	StateStrings = f.readline()
	StateStrings = StateStrings.replace("#", "\t"); StateStrings = StateStrings.replace("(eV/u)","\t")
	StateStrings = StateStrings.split() #strings corresponding to the capturing states
	del StateStrings[-1] # remove Total label
	for i in range(closestCollEnergy[0]): f.readline() #skip lines until reaching closest energy
	CrossSections = (f.readline()).split(); del CrossSections[0] #remove first entry which is coll energy
	del CrossSections[-1] #remove final entry which is total cross section
	f.close()

	#If collision is bare, apply distribution functions
	if collisionType == "bare" and method == "mclz":
		#Distribution Functions
		def statistical(n, l): #See Krasnopolsky et al. (2004)
			return ((2.0*float(l)+1.0)/(float(n)**2.0))

		def lowEnergy(n, l): #See Krasnopolsky et al. (2004)
			return ((2.0*float(l)+1.0)*((math.factorial(float(n)-1.0))**2.0)/((math.factorial(float(n)\
				+float(l)))*(math.factorial(float(n)-1.0-float(l)))))

		def lowEnergyModifiedACX(n, l): #See Smith et al. (2014)
			if n==1: return 0 #no population of n=1 states allowed
			else: return ((float(l))*(float(l) + 1.0)*(2.0*float(l)+1.0) * ((math.factorial(float(n) - 1.0)) \
				* (math.factorial(float(n) - 2.0)))/((math.factorial(float(n) + float(l))) * \
				(math.factorial(float(n)-float(l) - 1.0))))
				
		def separableACX(l, q): #See Smith et al. (2014)...only distribution function with charge.
			return (((2.0*float(l) +1.0)/q)*math.exp((-float(l)*(float(l)+1.0))/q))

		def SL1(n, lPrime): #See Mullen et al. (2016)
			l = lPrime + 1
			if lPrime == n-1: 
				return 0
			elif lPrime == 1: 
				return (lowEnergy(n, l) + lowEnergy(n,0)) #maintain normalization of distr. function
			else: 
				return lowEnergy(n, l)

		#Set the Distribution Function that is used
		if str(distro).lower()=="stat": 
			def lDistro(n,l): return statistical(n,l)
		elif str(distro).lower()=="lowe":
			def lDistro(n,l): return lowEnergy(n,l)
		elif str(distro).lower()=="lowemod":
			def lDistro(n,l): return lowEnergyModifiedACX(n,l)
		elif str(distro).lower()=="separable":
			def lDistro(l,q): return separableACX(l,q)
		elif str(distro).lower()=="sl1":
			def lDistro(n,l): return SL1(n,l)
		else: 
			print "Error: Do not recognize distribution function."
			sys.exit("Exiting Unsuccessfully.\n" + newLine + hashLine)

		nData = [(int(StateStrings[i][2:]), float(CrossSections[i])) for i in range(len(StateStrings))]
		csData =[]
		#csData Structure (n,l,S, crosssection)
		for i in range(len(nData)):
			for j in range(nData[i][0]):
				if str(distro)!="separable": # not equal to separable because separable takes l & charge as arg.
					csData.append((nData[i][0], j, 2, lDistro(nData[i][0], j)*nData[i][1]))
				else:
					csData.append((nData[i][0], j, 2, lDistro(j, charge)*nData[i][1])) # arguments l & charge
		csData = sorted(csData, key=lambda x: (x[0], x[1])) # sort data by n and l


	#Sort data by ascending nlS.  Only for non-bare collisions.
	#Use dictionary to take l state and go to corresponding number (note that there is no j, p, s)
	#Assumes that all data will have capture into n-state <= 15.  Used for non-bare collisions.  
	if collisionType == "non-bare" or method != "mclz":

		dict = {"s":0, "p":1, "d":2, "f":3, "g":4, "h":5, "i":6, "k":7, "l":8, "m":9, \
		"n":10, "o":11, "q":12, "r":13, "t":14, "1":0, "3":1, "2":2}

		#csData Structure (n,l,S, crosssection)
		csData = sorted([(int(StateStrings[i][0:-5]), dict[StateStrings[i][-5]], dict[StateStrings[i][-3]], \
			float(CrossSections[i])) for i in range(len(StateStrings))], key=lambda x: (x[0], x[1], x[2]))

	#Output Data to File Required by Cascade Model
	if not os.path.exists("./bin"):
		os.makedirs("./bin")
	cascadeFile = open("./bin/init.dat", "w")

	if collisionType == "bare":
		for i in range(0, len(csData)):
			cascadeFile.write( str(csData[i][0]) + "\t" + str(csData[i][1]) + "\t" + \
				str(csData[i][3]) + "\n" )
		cascadeFile.write(str(100) + "\t" + str(0) + "\t" + str(0.0) + "\n")
		cascadeFile.write(str(charge) + "\t" + str(rmax) + "\t" + str(closestCollEnergy[1]) + "\t" + \
			str(ionizationEnergy*(-1.23981E-4)) + "\n") #Final entry is -E_0 for enr in radcash 
		cascadeFile.close()

	if collisionType == "non-bare":
		for i in range(0, len(csData)):
			cascadeFile.write( str(csData[i][0]) + "\t" + str(csData[i][1]) + "\t" + \
				str(csData[i][2]) + "\t" + str(csData[i][3]) + "\n" )
		cascadeFile.write(str(100) + "\t" + str(0) + "\t" + str(0) + "\t" + str(0.0) + "\n")
		cascadeFile.write(str(closestCollEnergy[1]) + "\n") #Final entry is coll Energy for use by radcash
		cascadeFile.close()

def doCascade():
	#Compile and Run radcash
	os.chdir("./bin/") # Change working directory to ./bin/
	if collisionType == "bare":
		subprocess.call(["gfortran", "-o", "cascade", "../src/radcash.f"])
		print "...Now Running Cascade Model for X-ray Emission."
		print "...This May Take a While.\n"
		subprocess.call(["./cascade"])

	#Copy Radrathe File; Compile and Run radcashe
	if collisionType == "non-bare":
		shutil.copyfile("."+RADRATHE, "./radrathe.dat")
		subprocess.call(["gfortran", "-o", "cascade", "../src/radcashe.f"])
		print "...Now Running Cascade Model for X-ray Emission.\n"
		subprocess.call(["./cascade"])

def makeSpectrum():
	#Create gauss.in file for gauss.f
	gaussFile = open("./gauss.in", "w")
	gaussFile.write(str(resolution) + "\t" + str(numberOfPoints) + "\t" + "0.d0" + "\t" + "1000.d0" + "\n")
	gaussFile.close()

	#Compile and Run gauss
	subprocess.call(["gfortran", "-o", "gauss", "../src/gauss.f"])
	subprocess.call(["./gauss"])

def lineRatios():
	#Turns integer into Roman numeral for results output
	def int_to_roman(integer):
	    returnstring=""
	    table=[["X",10],["IX",9],["V",5],["IV",4],["I",1]]
	    for pair in table:
	        while integer-pair[1]>=0:
	            integer-=pair[1]
	            returnstring+=pair[0]
	    return returnstring
	
	#Converts Collision Energy in (eV/u) to (km/s)
	def convertSpeed(E):
		return (((E/1000.0)/25.0)**0.5 * 2.1877E3)

	#Dictionary that converts transition to line Name (i.e. n=2, l=1 -> n=1, l=0 goes to Ly alpha)
	dict = {"2110":"Ly $\\alpha$", "3110":"Ly $\\beta$", "4110":"Ly $\\gamma$", "5110":"Ly $\\delta$",\
	"6110":"Ly $\\epsilon$", "7110":"Ly $\\zeta$", "8110":"Ly $\\eta$", "9110":"Ly $\\theta$",\
	"10110":"Ly $\\iota$", "11110":"Ly $\\kappa$", "12110":"Ly $\\lambda$", "13110":"Ly $\\mu$",\
	"14110":"Ly $\\nu$", "15110":"Ly $\\xi$", "16110":"Ly $\\omicron$", "17110":"Ly $\\pi$", \
	"18110":"Ly $\\rho$", "19110":"K $\\sigma$","210100":"K $\\alpha$ r", "211100":"K $\\alpha$ i" ,\
	"201100":"K $\\alpha$ f", "310100":"K $\\beta$", "410100":"K $\\gamma$", "510100":"K $\\delta$",\
	"610100":"K $\\epsilon$", "710100":"K $\\zeta$", "810100":"K $\\eta$", "910100":"K $\\theta$",\
	"1010100":"K $\\iota$", "1110100":"K $\\kappa$", "1210100":"K $\\lambda$", "1310100":"K $\\mu$",\
	"1410100":"K $\\nu$", "1510100":"K $\\xi$", "1610100":"K $\\omicron$", "1710100":"K $\\pi$",\
	"1810100":"K $\\rho$", "1910100":"K $\\sigma$", "2010100":"K $\\tau$"}

	#Generate Cascade Data from radcash(e) Output Files
	if collisionType == "bare":
		ni,li,nf,lf,energy,line,linenorm = loadtxt("./phot.dat",dtype=str, unpack=True,\
			usecols=[0,1,3,4,6,7,9])
		transition=[(i+j+k+l) for i,j,k,l in zip(ni,li,nf,lf)]
	if collisionType == "non-bare":
		ni,li,si,nf,lf,sf,energy,line,linenorm = loadtxt("./phot.dat",dtype=str, unpack=True,\
			usecols=[0,1,2,4,5,6,8,9,11])
		transition=[(i+j+k+l+m+n) for i,j,k,l,m,n in zip(ni,li,si,nf,lf,sf)]
	CascadeData = sorted(zip(transition, energy, line, linenorm), key=lambda x: x[0])
	CascadeData = [(dict[i[0]],i[1],i[2],i[3]) if i[0] in dict else (i[0],i[1],i[2],i[3]) for i in CascadeData]

	#Write Line Ratios File
	lineRatiosFile = open("./line_ratios.dat","w")
	lineRatiosFile.write("#Collision Energy = " + str(closestCollEnergy[1]) + " eV/u" + "\t" + \
		"Collision Velocity = " + str(convertSpeed(closestCollEnergy[1])) + "\n")
	lineRatiosFile.write("#Method = " + method.upper() + "\n")
	if collisionType == "bare" and method == "mclz": lineRatiosFile.write("#Distrib = " + distro + "\n")
	lineRatiosFile.write("#EmissionLine" + "\t" + "Photon Energy (eV)" + "\t" + "LineRatio" + "\t"\
		+ "NormLineRatio" + "\n")
	for i in range(len(CascadeData)):
		if CascadeData[i][0][0:1]=="L" or CascadeData[i][0][0:1]=="K":
			lineRatiosFile.write(CascadeData[i][0] + "\t" + CascadeData[i][1] + "\t" + \
				CascadeData[i][2] + "\t" + CascadeData[i][3] + "\n")
	lineRatiosFile.close()

	#Write Latex Table
	shutil.copyfile("../src/aastex.cls", "./aastex.cls")
	texFile = open("./line_ratios.tex","w")

	#Preamble and Begin Table
	texFile.write("\\documentclass{aastex}"+"\n"+"\\begin{document}"+"\n"+"\\begin{deluxetable}{cccc}"+"\n"+\
		"\\tablecolumns{4}"+"\n"+"\\tablewidth{0pc}"+"\n"+"\\tablecaption{"+projectile.upper()+" " +\
		int_to_roman(charge)+" "+method.upper()+" Line Ratios for "+target+" Targets}" +\
		"\\tablehead{\multicolumn{4}{c}{Collision Energy = "+"%.1f" % float(closestCollEnergy[1])+\
		" eV/u $\\vert$ Collsion Velocity = "+"%.1f" % float(convertSpeed(closestCollEnergy[1]))+" km/s}"+\
		"\\\\"+"\n"+"\\colhead{Emission Line} & \colhead{Photon Energy (eV)} & \\colhead{Line Ratio}"+\
		"& \\colhead{Norm. Line Ratio}}"+"\n"+"\\startdata"+"\n")

	#Cascade Data
	for i in range(len(CascadeData)):
		if CascadeData[i][0][0:1]=="L" or CascadeData[i][0][0:1]=="K":
			texFile.write(CascadeData[i][0] + "\t&\t" + ("%.2f" % float(CascadeData[i][1])) + "\t&\t" + \
				("%.4f" % float(CascadeData[i][2]))+"\t&\t"+("%.4f" % float(CascadeData[i][3]))+"\\\\"+"\n")

	#End Table and Close Document
	texFile.write("\\\\"+"\n"+"\\hline"+"\n"+"\\enddata"+"\n")
	texFile.write("\\tablecomments{These line ratios were generated using " + method.upper() + " Data. ")
	if collisionType=="bare" and method == "mclz":
		texFile.write("The " + distro + " distribution function was applied to the MCLZ $n$-resolved \
			cross sections.}\n")
	else: 
		texFile.write("}\n")
	texFile.write("\\end{deluxetable}"+"\n"+"\\end{document}")
	texFile.close()

	#Build, Compile and Run TeX Document
	FNULL = open(os.devnull, "w") #pdflatex output to dev/null
	subprocess.call(["pdflatex", "line_ratios.tex"], stdout=FNULL)
	FNULL.close()

def plotSpectrum():
	#Plots Normalized Spectrum if -G is set to yes
	if graphCheck != None and 'y' in graphCheck.lower():
		photEnergiesSpec, IntensitySpec = loadtxt("spectrum.dat",dtype=float, unpack=True, usecols=[0,1])
		photEnergiesLines, IntensityLines = loadtxt("photline.dat", dtype=float, unpack=True, usecols=[0,1])
		matplotlib.pyplot.plot(photEnergiesSpec, IntensitySpec)
		matplotlib.pyplot.vlines(photEnergiesLines,[0],IntensityLines)
		matplotlib.pyplot.legend(['spectrum.dat', 'photline.dat'])
		matplotlib.pyplot.suptitle("Projectile Ion: " + projectile + "   " + "Charge: " + str(charge) + "   " \
			+ "Target: " + target + "   " + "Coll. Energy: " + str(closestCollEnergy[1]) + " eV/u" + "   " + \
			"FWHM: " + str(resolution) + " eV")
		matplotlib.pyplot.xlabel ("Photon Energy (eV)")
		matplotlib.pyplot.ylabel("Relative Intensity (arb. units)")
		matplotlib.pyplot.draw()

def moveFiles():
	#Make Directory Tree for Results
	resultsPath1 = "../results/" + projectile.lower()  + str(charge) + "+" + target.lower()
	if not os.path.exists(resultsPath1):
		os.makedirs(resultsPath1)
	resultsTree = resultsPath1 + "/" + "%.1f" % float(closestCollEnergy[1]) + "eV-u_" + \
		str(time.strftime("%d.%m.%y@%H.%M.%S") + "/TeXTable")
	os.makedirs(resultsTree)

	#Move Files to this Directory
	for filename in os.listdir("."):
		if filename.endswith(".dat"):
			shutil.move(filename, resultsTree + "/../")
		if filename.endswith(".pdf") or filename.endswith(".tex") or filename.endswith(".cls"):
			shutil.move(filename, resultsTree)

def deleteBin():
	#Delete remaining contents of ./bin/
	for filename in os.listdir("."):
		os.remove(filename)

def end():
	#Keep Program Running Until Plot Closed (if -G yes)
	if graphCheck != None and 'y' in graphCheck.lower():
		matplotlib.pyplot.show()

	#Close Program Successfully
	print hashLine
	print "...Exiting Successfully.\n"
	print hashLine
