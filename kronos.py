#! /usr/bin/env python
# Developed by: Patrick Mullen

###################################################################

import matplotlib
import argparse 

###################################################################

parser=argparse.ArgumentParser()
parser.add_argument("-I", "--projectile", help="Specify Projectile Ion -- Required", required=True)
parser.add_argument("-Z", "--charge", help="Specify Projectile Charge -- Required", required=True)
parser.add_argument("-A", "--target", help="Specify Target Neutral -- Required", required=True)
parser.add_argument("-E", "--collisionenergy", help="Specify Collision Energy in eV/u -- Required", required=True)
parser.add_argument("-D", "--distribution", help="Choose lowe or stat distribution -- Only Required for MCLZ bare ion collisions", required=False)
parser.add_argument("-R", "--resolution", help="Specify the FWHM Resolution in eV -- Optional", required=False)
parser.add_argument("-P", "--numberofpoints", help="Specify number of data points -- Optional", required=False)
parser.add_argument("-N", "--name", help="Specify Your Name -- Optional", required=False)
parser.add_argument("-L", "--initials", help="Specify Your Initials -- Optional", required=False)
parser.add_argument("-X", "--rmax", help= "Specify rmax -- Optional", required=False)
parser.add_argument("-M", "--method", help="Specify method of Charge Exchange Data Calculation -- Optional", required=False)
parser.add_argument("-G", "--graph", help="Specify yes or no to graph -- Optional", required=False)
args=parser.parse_args()

from pylab import *
import os
import sys

projectile = args.projectile
charge = args.charge
target = args.target
distro = str(args.distribution)
collenergy = args.collisionenergy
methodcheck = str(args.method)
graphcheck = str(args.graph)
if args.resolution == None:
	resolution = 150.0
else:
	resolution = float(args.resolution)
if args.numberofpoints == None:
	numberofpoints = int(4200.0)
else:
	numberofpoints = int(args.numberofpoints)

if args.rmax == None:
	rmax = 200
else:
	rmax = float(args.rmax)

print '\n'

print '###################################################'
print '\n'
if args.name != None:
	print '...Welcome,', args.name
print '...This is Kronos.'
print '\n'
print '###################################################'
print '\n'
print '...You have selected the CX SEC collision:', '\n'
print str(projectile) + '^{' + str(charge) + '+} + ' + str(target) +  ' --> ' + str(projectile) + '^{' + str(int(charge) -1) + \
'+} + ' + str(target) + '^{+}'
print 'Projectile Ion: ', str(projectile)
print 'Projectile Ion Charge: ', str(charge)
print 'Target Neutral: ', str(target)
print 'Collision Energy: ', str(collenergy), 'eV/u'
print 'FWHM Resolution: ', str(resolution), 'eV'
print '\n'

atomicnumber,Ion = loadtxt('CXDatabase/elements.dat', dtype = str, unpack = True, usecols=[0,1])

arrayION = np.array(Ion)
ions = arrayION.tolist()
foundit = str(projectile) in Ion

if foundit == True: 
	founditindex = ions.index(str(projectile))
	founditcharge = atomicnumber[int(founditindex)]
	if int(founditcharge) == int(charge):
		print '...I have detected that this is a bare ion collision.'
		barisity = 'bare'
	else:
		print '...I have detected that this is a non-bare ion collision.'
		barisity = 'nonbare'
		if str(distro) != 'None':
			print '\n'
			print '...There is no need to specify distribution function! This is a non-bare ion collision!'
if foundit == False:
	print '...This projectile ion does not exist in my universe.'
	sys.exit("...Exiting Kronos Unsuccessfully" + '\n' + '###################################################')

print '\n'
print '###################################################'

os.system('ls CXDatabase/Projectile_Ions >> availability.txt')
databaseprojectiles = loadtxt('availability.txt', unpack = True, dtype=str)
os.remove('availability.txt')
check1 = str(projectile) in databaseprojectiles
if check1 != True:
	print '\n'
	print 'No Data for this System.  Please provide charge exchange cross sections.'
	print '\n'
	sys.exit("...Exiting Kronos Unsuccessfully"  + '\n' + '###################################################')

os.system('ls CXDatabase/Projectile_Ions/' + str(projectile) +  '/Charge/ >> availability.txt')
databasecharges = loadtxt('availability.txt', unpack = True, dtype=str)
os.remove('availability.txt')
check2 = str(charge) in databasecharges
if check2 != True:
	print '\n'
	print '...No Data for this System.  Please provide charge exchange cross sections.'
	print '\n'
	sys.exit("...Exiting Kronos Unsuccessfully"  + '\n' + '###################################################')

os.system('ls CXDatabase/Projectile_Ions/' + str(projectile) +  '/Charge/' + str(charge) + '/Targets/ >> availability.txt')
databasetargets = loadtxt('availability.txt', unpack = True, dtype=str)
os.remove('availability.txt')
check3 = str(target) in databasetargets
if check3 != True:
	print '\n'
	print '...No Data for this System.  Please provide charge exchange cross sections.'
	print '\n'
	sys.exit("...Exiting Kronos Unsuccessfully"  + '\n' + '###################################################')

os.system('find CXDatabase/Projectile_Ions/' + str(projectile) +  '/Charge/' + str(charge) + '/Targets/ -name "*.cs" >> availability.txt')
check4 = os.stat('availability.txt').st_size == 0
os.remove('availability.txt')
if check4 == True:
	print '\n'
	print '...No Data for this System.  Please provide charge exchange cross sections.'
	print '\n'
	sys.exit("...Exiting Kronos Unsuccessfully"  + '\n' + '###################################################')

methods = list()
os.system('find CXDatabase/Projectile_Ions/' + str(projectile) +  '/Charge/' + str(charge) + '/Targets/' + str(target) + \
	' -name "*.cs" >> availability.txt')
databasefiles = loadtxt('availability.txt', unpack = True, dtype=str)
os.remove('availability.txt')
lengthcheck = len('CXDatabase/Projectile_Ions/' + str(projectile) +  '/Charge/' + str(charge) + '/Targets/' + str(target) + '/' \
	+ str(projectile.lower()) + str(charge) + '+' + str(target.lower()) + '_sec_')

if int(databasefiles.ndim) == 0:
	methods.append(str(databasefiles)[(int(lengthcheck)):(int(lengthcheck)+4)])	
else:
	for i in range(0, len(databasefiles)):
		methods.append(databasefiles[i][(int(lengthcheck)):(int(lengthcheck)+4)])

check5 = 'mclz' in methods
check6 = 'ctmc' in methods
check7 = 'aocc' in methods
check8 = 'mocc' in methods
check9 = 'qmoc' in methods
check10 = 'rcmd' in methods

print '\n'
print '...Checking availability of SRG Recommended:', 
if check10 == True:
	print 'Yes'
else:
	print 'No'

print '...Checking availability of QMOCC:', 
if check9 == True:
	print 'Yes'
else:
	print 'No'

print '...Checking availability of MOCC:', 
if check8 == True:
	print 'Yes'
else:
	print 'No'

print '...Checking availability of AOCC:', 
if check7 == True:
	print 'Yes'
else:
	print 'No'

print '...Checking availability of CTMC:', 
if check6 == True:
	print 'Yes'
else:
	print 'No'

print '...Checking availability of MCLZ:', 
if check5 == True:
	print 'Yes'
else:
	print 'No'

if str(args.method) == 'RCMD' and check10 != True:
	print '\n'
	print '...No SRG Recommended Data for this System. Specify a different method, for which there is data.'
	print '\n'
	sys.exit("...Exiting Kronos Unsuccessfully"  + '\n' + '###################################################')

if str(args.method) == 'QMOCC' and check9 != True:
	print '\n'
	print '...No QMOCC Data for this System. Specify a different method, for which there is data.'
	print '\n'
	sys.exit("...Exiting Kronos Unsuccessfully"  + '\n' + '###################################################')

if str(args.method) == 'MOCC'  and check8 != True:
	print '\n'
	print '...No MOCC Data for this System. Specify a different method, for which there is data.'
	print '\n'
	sys.exit("...Exiting Kronos Unsuccessfully"  + '\n' + '###################################################')

if str(args.method) == 'AOCC'  and check7 != True:
	print '\n'
	print '...No AOCC Data for this System. Specify a different method, for which there is data.'
	print '\n'
	sys.exit("...Exiting Kronos Unsuccessfully"  + '\n' + '###################################################')

if str(args.method) == 'CTMC' and check6 != True:
	print '\n'
	print '...No CTMC Data for this System. Specify a different method, for which there is data.'
	print '\n'
	sys.exit("...Exiting Kronos Unsuccessfully"  + '\n' + '###################################################')

if str(args.method) == 'MCLZ' and check5 != True:
	print '\n'
	print '...No MCLZ Data for this System. Specify a different method, for which there is data.'
	print '\n'
	sys.exit("...Exiting Kronos Unsuccessfully"  + '\n' + '###################################################')

print '\n'

if check5 == True:
	method = 'mclz'
	distro = str(args.distribution) 
if check6 == True:
	method = 'ctmc'
if check7 == True:
	method = 'aocc'
if check8 == True:
	method = 'mocc'
if check9 == True:
	method = 'qmocc'
if check10 == True:
	method = 'rcmd'
if methodcheck == 'None':
	pass
else:
	method=str(args.method).lower()

if method == 'rcmd':
	if str(args.distribution) != 'None':
		print '...There is no need to specify a distribution function for'
		print '...RCMD/QMOCC/MOCC/AOCC/CTMC data. Ignoring "-D ' + str(args.distribution) + '" option!'
		print '\n'
	
if method == 'qmocc':
	if str(args.distribution) != 'None':
		print '...There is no need to specify a distribution function for'
		print '...QMOCC/MOCC/AOCC/CTMC data. Ignoring "-D ' + str(args.distribution) + '" option!'
		print '\n'
if method == 'mocc':
	if str(args.distribution) != 'None':
		print '...There is no need to specify a distribution function for'
		print '...QMOCC/MOCC/AOCC/CTMC data. Ignoring "-D ' + str(args.distribution) + '" option!'
		print '\n'
if method == 'aocc':
	if str(args.distribution) != 'None':
		print '...There is no need to specify a distribution function for'
		print '...QMOCC/MOCC/AOCC/CTMC data. Ignoring "-D ' + str(args.distribution) + '" option!'
		print '\n'
if method == 'ctmc':
	if str(args.distribution) != 'None':
		print '...There is no need to specify a distribution function for'
		print '...QMOCC/MOCC/AOCC/CTMC data. Ignoring "-D ' + str(args.distribution) + '" option!'
		print '\n'
if method == 'mclz':
	distro = str(args.distribution)

if str(distro) == 'None' and str(barisity) == 'bare' and method == 'mclz':
	print '...Please Specify Distribution Function!'
	print '...If MCLZ data, you must specify distribution as lowe or stat.'
	print '\n'
	sys.exit("...Exiting Kronos Unsuccessfully" + '\n' + '###################################################')

if barisity == 'bare':
	print '...Proceeding with ' + str(method).upper() + ' Cross Sections.'
	if method == 'ctmc' or method == 'aocc' or method == 'mocc' or method == 'qmocc' or method == 'rcmd':
		line = 'python SourceCode/Bare/initmaker_h.py -q ' + str(collenergy) + ' -charge ' + str(charge) + ' -i ' + \
		'CXDatabase/Projectile_Ions/' + str(projectile) + '/Charge/' + str(charge) + '/Targets/'  \
		+  str(target) + '/' + str(projectile.lower()) + str(charge) + '+' + str(target.lower()) + '_sec_' + str(method) + '.cs -n '\
		+ 'CXDatabase/Projectile_Ions/' + str(projectile) + '/Charge/' + str(charge) + '/Targets/'  +  str(target)\
		+ '/' + str(projectile.lower()) + str(charge) + '.nist.in -rmax ' + str(rmax)
	else:
		line = 'python SourceCode/Bare/initmaker_h.py -q ' + str(collenergy) + ' -charge ' + str(charge) + ' -i ' + \
		'CXDatabase/Projectile_Ions/' + str(projectile) + '/Charge/' + str(charge) + '/Targets/'  \
		+  str(target) + '/' + str(projectile.lower()) + str(charge) + '+' + str(target.lower()) + '_sec_' + str(method) + '_' + \
		str(distro).lower() + '.cs -n '\
		+ 'CXDatabase/Projectile_Ions/' + str(projectile) + '/Charge/' + str(charge) + '/Targets/'  +  str(target)\
		+ '/' + str(projectile.lower()) + str(charge) + '.nist.in -rmax ' + str(rmax)
	os.system(str(line))
	res, pts, start,end = loadtxt('SourceCode/Bare/gauss.in', dtype = str, unpack = True, usecols=[0,1,2,3])
	os.remove('SourceCode/Bare/gauss.in')
	outputfile = open('SourceCode/Bare/gauss.in', 'w')
	outputfile.write(str(float(resolution)) + '\t' + str(int(numberofpoints)) + '\t' + str(start) + '\t' + str(end))
	outputfile.close()
	print '\n'
	print '###################################################'
	print '\n'
	print '...Now running the Cascade Model for X-ray Emission.'
	os.system('mv inith.dat SourceCode/Bare/inith.dat')
	os.system('gfortran -o SourceCode/Bare/radcash_h SourceCode/Bare/radcash_h.f')
	os.system('./SourceCode/Bare/radcash_h')
	os.system('gfortran -o SourceCode/Bare/gauss_h SourceCode/Bare/gauss_h.f')
	os.system('./SourceCode/Bare/gauss_h')
	os.remove('fort.25')
	os.remove('fort.35')
	os.remove('fort.90')
	closestenergy = loadtxt('closestenergy.dat', unpack = True, dtype=float)
	os.remove('closestenergy.dat')
	collenergy=closestenergy
	os.system('mkdir -p Results/' + str(projectile) + '_' + str(charge))
	os.system('mkdir -p Results/' + str(projectile) + '_' + str(charge) + '/' + str(target))
	os.system('mkdir -p Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy)) + 'eV-u@' + str(resolution) + 'eV_FWHM')
	os.system('mkdir -p Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
		'eV_FWHM/TeX-Table/')
	if method == 'mclz':
		os.system('mv cevolh.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) \
			+ 'eV_FWHM/cevolh_' + str(method) + '_' + str(distro) +'.dat')
		os.system('mv photh.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/photh_' +str(method) + '_' + str(distro) +'.dat')
		os.system('mv phot1h.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/phot1h_' + str(method) + '_' + str(distro)+ '.dat')
		os.system('mv photline.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) +\
		 'eV_FWHM/photline_' +str(method) + '_' + str(distro) + '.dat')
		os.system('mv photouth.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/photouth_' + str(method) + '_' + str(distro)+ '.dat')
		os.system('mv ratios.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/ratios_' + str(method) + '_' + str(distro)+ '.dat')
		os.system('mv spectrum.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/spectrum_' + str(method) + '_' + str(distro) + '.dat')
		os.system('mv SourceCode/Bare/inith.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + \
			str(resolution) + 'eV_FWHM/inith_' + str(method) + '_' + str(distro) +  '.dat')
		x,y = loadtxt('Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + 'eV_FWHM/spectrum_' \
			+ str(method) + '_' + str(distro) + '.dat', dtype = float, unpack = True, usecols=[0,1])
		z,xi = loadtxt('Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + 'eV_FWHM/photline_' + \
			str(method) + '_' + str(distro) + '.dat', dtype = float, unpack = True, usecols=[0,1])
	else:
		os.system('mv cevolh.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) \
			+ 'eV_FWHM/cevolh_' + str(method) +'.dat')
		os.system('mv photh.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/photh_' +str(method) +'.dat')
		os.system('mv phot1h.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/phot1h_' + str(method)+ '.dat')
		os.system('mv photline.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) +\
		 'eV_FWHM/photline_' +str(method) + '.dat')
		os.system('mv photouth.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/photouth_' + str(method)+ '.dat')
		os.system('mv ratios.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/ratios_' + str(method)+ '.dat')
		os.system('mv spectrum.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/spectrum_' + str(method) + '.dat')
		os.system('mv SourceCode/Bare/inith.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + \
			str(resolution) + 'eV_FWHM/inith_' + str(method) +  '.dat')
		x,y = loadtxt('Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + 'eV_FWHM/spectrum_' \
			+ str(method) + '.dat', dtype = float, unpack = True, usecols=[0,1])
		z,xi = loadtxt('Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + 'eV_FWHM/photline_' + \
			str(method) + '.dat', dtype = float, unpack = True, usecols=[0,1])		
	plot(x, y)
	vlines(z, [0], xi)
	legend(['spectrum.dat', 'photline.dat'])
	suptitle('Projectile Ion: ' + str(projectile) + '   '  'Charge: ' + str(charge) + '   ' + 'Target: ' + str(target) + '   ' + 'Collision Energy: ' + \
		'%s' % float('%.5g' % float(collenergy))  + ' eV/u' + '   ' + 'FWHM Resolution: ' + str(resolution) + ' eV')
	xlabel('Photon Energy (eV)')
	ylabel('Relative Intensty (Arbitrary)')
	print '\n'
	os.system('cp SourceCode/Bare/aastex.cls aastex.cls')
	if method == 'mclz':
		os.system('python SourceCode/Bare/zeus_h.py -I' + str(projectile) + ' -Z ' + str(charge) + ' -A ' + str(target) + ' -E ' + str(collenergy) + ' -F Results/' \
			+ str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  \
			+ 'eV-u@' + str(resolution) + 'eV_FWHM/phot1h_' + str(method) + '_' + str(distro)+ '.dat')
		os.system('mv LineRatios.aux Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/TeX-Table/LineRatios_' + str(method) + '_' + str(distro) + '.aux')
		os.system('mv LineRatios.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/LineRatios_' + str(method) + '_' + str(distro) + '.dat')
		os.system('mv LineRatios.log Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/TeX-Table/LineRatios_' + str(method) + '_' + str(distro) + '.log')
		os.system('mv LineRatios.pdf Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/TeX-Table/LineRatios_' + str(method) + '_' + str(distro) + '.pdf')
		os.system('mv LineRatios.tex Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/TeX-Table/LineRatios_' + str(method) + '_' + str(distro) + '.tex')
	else:
		os.system('./SourceCode/Bare/zeus_h.py -I' + str(projectile) + ' -Z ' + str(charge) + ' -A ' + str(target) + ' -E ' + '%s' % float('%.5g' % float(collenergy))  + ' -F Results/' \
			+ str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  \
			+ 'eV-u@' + str(resolution) + 'eV_FWHM/phot1h_' + str(method)+ '.dat')
		os.system('mv LineRatios.aux Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/TeX-Table/LineRatios_' + str(method) + '.aux')
		os.system('mv LineRatios.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/LineRatios_' + str(method) + '.dat')
		os.system('mv LineRatios.log Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/TeX-Table/LineRatios_' + str(method) + '.log')
		os.system('mv LineRatios.pdf Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/TeX-Table/LineRatios_' + str(method) + '.pdf')
		os.system('mv LineRatios.tex Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy))  + 'eV-u@' + str(resolution) + \
			'eV_FWHM/TeX-Table/LineRatios_' + str(method) + '.tex')
	os.remove('aastex.cls')
	print '\n'
	print "...Successfully Exiting Kronos."
	print '\n'
	print '###################################################'


if barisity == 'nonbare':
	print '...Proceeding with ' + str(method).upper() + ' Cross Sections.'
	line = 'python SourceCode/NonBare/initmaker_he.py -q ' + str(collenergy) + ' -charge ' + str(charge) + ' -i ' + 'CXDatabase/Projectile_Ions/' + \
	str(projectile) + '/Charge/' + str(charge) + '/Targets/'  +  \
	str(target) + '/' + str(projectile.lower()) + str(charge) + '+' + str(target.lower()) + '_sec_' + str(method) + '.cs' 
	os.system(str(line))
	res, pts, start,end = loadtxt('SourceCode/NonBare/gauss.in', dtype = str, unpack = True, usecols=[0,1,2,3])
	os.remove('SourceCode/NonBare/gauss.in')
	outputfile = open('SourceCode/NonBare/gauss.in', 'w')
	outputfile.write(str(float(resolution)) + '\t' + str(int(numberofpoints)) + '\t' + str(start) + '\t' + str(end))
	outputfile.close()
	print '\n'
	print '###################################################'
	print '\n'
	print '...Now running the Cascade Model for X-ray Emission.'
	os.system('cp CXDatabase/Projectile_Ions/' + str(projectile) + '/Charge/' + str(charge) + '/Targets/'  +  str(target) + '/radrathe.dat SourceCode/NonBare/radrathe.dat')
	os.system('mv inithe.dat SourceCode/NonBare/inithe.dat')
	os.system('gfortran -o SourceCode/NonBare/radcashe_he SourceCode/NonBare/radcashe_he.f')
	os.system('./SourceCode/NonBare/radcashe_he')
	os.system('gfortran -o SourceCode/NonBare/gauss_he SourceCode/NonBare/gauss_he.f')
	os.system('./SourceCode/NonBare/gauss_he')
	os.remove('fort.33')
	os.remove('SourceCode/NonBare/radrathe.dat')
	closestenergy = loadtxt('closestenergy.dat', unpack = True, dtype=float)
	os.remove('closestenergy.dat')
	collenergy=closestenergy
	os.system('mkdir -p Results/' + str(projectile) + '_' + str(charge))
	os.system('mkdir -p Results/' + str(projectile) + '_' + str(charge) + '/' + str(target))
	os.system('mkdir -p Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy)) + 'eV-u@' + str(resolution) + 'eV_FWHM')
	os.system('mkdir -p Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy)) + 'eV-u@' + str(resolution) + 'eV_FWHM/TeX-Table/')
	os.system('mv cevolhe.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy)) + 'eV-u@' + str(resolution) + 'eV_FWHM/cevolhe.dat')
	os.system('mv phothe.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy)) + 'eV-u@' + str(resolution) + 'eV_FWHM/phothe.dat')
	os.system('mv photline.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy)) + 'eV-u@' + str(resolution) + 'eV_FWHM/photline.dat')
	os.system('mv photout.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy)) + 'eV-u@' + str(resolution) + 'eV_FWHM/photout.dat')
	os.system('mv ratios.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy)) + 'eV-u@' + str(resolution) + 'eV_FWHM/ratios.dat')
	os.system('mv spectrum.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy)) + 'eV-u@' + str(resolution) + 'eV_FWHM/spectrum.dat')
	os.system('mv SourceCode/NonBare/inithe.dat Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy)) + 'eV-u@' + str(resolution) + \
		'eV_FWHM/inithe.dat')
	x,y = loadtxt('Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy)) + 'eV-u@' + str(resolution) + 'eV_FWHM/spectrum.dat', dtype = float, \
		unpack = True, usecols=[0,1])
	z,xi = loadtxt('Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy)) + 'eV-u@' + str(resolution) + 'eV_FWHM/photline.dat', dtype = float, \
		unpack = True, usecols=[0,1])
	plot(x, y)
	vlines(z, [0], xi)
	legend(['spectrum.dat', 'photline.dat'])
	suptitle('Projectile Ion: ' + str(projectile) + '   '  'Charge: ' + str(charge) + '   ' + 'Target: ' + str(target) + '   ' + 'Collision Energy: ' + '%s' % float('%.5g' % float(collenergy)) + ' eV/u' + '   ' \
		+ 'FWHM Resolution: ' + str(resolution) + ' eV')
	xlabel('Photon Energy (eV)')
	ylabel('Relative Intensity (Arbitrary)')
	print '\n'
	os.system('cp SourceCode/NonBare/aastex.cls aastex.cls')
	os.system('python SourceCode/NonBare/zeus_he.py -I' + str(projectile) + ' -Z ' + str(charge) + ' -A ' + str(target) + ' -E ' + '%s' % float('%.5g' % float(collenergy)) + ' -F Results/' + str(projectile) + '_' + \
		str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy)) \
		+ 'eV-u@' + str(resolution) + 'eV_FWHM/phothe.dat')
	os.system('mv Line* Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy)) + 'eV-u@' + str(resolution) + 'eV_FWHM/TeX-Table/')
	os.system('mv Results/' + str(projectile) + '_' + str(charge) + '/' + str(target) + '/' + '%s' % float('%.5g' % float(collenergy)) + 'eV-u@' + str(resolution) + 'eV_FWHM/TeX-Table/LineRatios.dat Results/' + \
		str(projectile) + '_' + str(charge) + '/' + str(target) \
		+ '/' + '%s' % float('%.5g' % float(collenergy)) + 'eV-u@' + str(resolution) + 'eV_FWHM/LineRatios.dat')
	os.remove('aastex.cls')
	print '\n'
	print "...Successfully Exiting Kronos."
	print '\n'
	print '###################################################'
if str(graphcheck).lower() == 'yes':
	show()
