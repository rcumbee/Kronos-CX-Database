#! /usr/bin/env python
# Developed by: Patrick Mullen

###################################################################

# Necessary imports are made here.  

from pylab import *
import argparse 
import os

###################################################################

# This provides feedback to the user that the program has initialized.  

print '...Creating Line Ratio Table and TeX File/PDF'

parser=argparse.ArgumentParser()
parser.add_argument("-I", "--projectile", help="Specify Projectile Ion")
parser.add_argument("-Z", '--charge' ,help="Specify Projectile Charge")
parser.add_argument("-E", "--collisionenergy", help = "Specify the Collision Energy")
parser.add_argument("-A", "--target", help = "Specify the Target")
parser.add_argument("-F", "--filename", help = "Specify the Filename")
args=parser.parse_args()

CollisionEnergy = float(args.collisionenergy)
Projectile = str(args.projectile)
Charge = int(args.charge)
Target = str(args.target)
filename = str(args.filename)

###################################################################

def int_to_roman (integer):

    returnstring=''
    table=[['X',10],['IX',9],['V',5],['IV',4],['I',1]]

    for pair in table:

        while integer-pair[1]>=0:

            integer-=pair[1]
            returnstring+=pair[0]

    return returnstring

def convertspeed(E):
	z = ((E/1000.0)/25.0)**0.5 * 2.1877E3
	return z

lineslist = list()

f = open(filename, 'r')
for line in f:
	lineslist.append(line)

del lineslist[-1]
del lineslist[-1]

for i in range(0, len(lineslist)): 
	lineslist[i] = lineslist[i].replace('phot(eV)=', '\t', 10)
	lineslist[i] = lineslist[i].replace('-->', '\t', 10)
	lineslist[i] = lineslist[i].replace('relat=', '\t', 10)
	lineslist[i] = lineslist[i].replace('\n', '\t', 10)

f = open('hephaistus.dat', 'w')
for i in range(0, len(lineslist)): 
	f.write(lineslist[i] + '\n')
f.close()

initial, final, photenergy, relative = loadtxt('hephaistus.dat', dtype=str, unpack = True, usecols=[0,1,2,-1])

os.remove('hephaistus.dat')

initial = initial.tolist()
final = final.tolist()
photenergy = photenergy.tolist()
relative = relative.tolist()

if len(initial) > 1:
	initial.reverse()
	final.reverse()
	photenergy.reverse()
	relative.reverse()

dict = {'21':'Lyman alpha', '31':'Lyman beta ', '41':'Lyman gamma', '51':'Lyman delta', '61':'Lyman epsilon', '71':'Lyman zeta', '81':'Lyman eta', '91':\
'Lyman theta', '101':'Lyman iota', '111': 'Lyman kappa', '121': 'Lyman lambda', '131': 'Lyman mu', '141': 'Lyman nu', '151':'Lyman xi', '161': 'Lyman omicron',\
'171': 'Lyman pi', '181': 'Lyman rho', '191': 'Lyman sigma', '201': 'Lyman tau'}

indices = [i for i, x in enumerate(final) if x == "1"]

tags = list()

for i in range(0, len(indices)):
	tags.append(dict[initial[indices[i]] + final[indices[i]]])

if len(initial) > 1:
	outputfile = open('LineRatios.dat', 'w')
	outputfile.write('#Collision Energy = ' + str(CollisionEnergy) + ' eV/u   Collision Velocity = ' + '%.1f' % float(convertspeed(float(CollisionEnergy))) + ' km/s '+'\n')
	outputfile.write('#Photon Energy (eV) ' + '\t' + 'Emission Line' + '\t' + 'Line Ratio' + '\n')
	for i in range(0, len(tags)):
		outputfile.write(str(photenergy[indices[i]]) + '\t' + str(tags[i]) + '\t' + relative[indices[i]] + '\n')
	outputfile.close()
else: 
	outputfile = open('LineRatios.dat', 'w')
	outputfile.write('#Collision Energy = ' + str(CollisionEnergy) + ' eV/u   Collision Velocity = ' + '%.1f' % float(convertspeed(float(CollisionEnergy))) + ' km/s '+'\n')
	outputfile.write('#Photon Energy (eV) ' + '\t' + 'Emission Line' + '\t' + 'Line Ratio' + '\n')
	for i in range(0, len(tags)):
		outputfile.write(str(photenergy) + '\t' + str(tags[i]) + '\t' + str(relative) + '\n')
	outputfile.close()

for i in range(0, len(tags)):
	tags[i] = tags[i].replace('Lyman ', 'Ly $\\')

if len(initial) > 1:
	outputfile = open('LineRatios.tex', 'w')
	outputfile.write('\documentclass{aastex}' + '\n')
	outputfile.write('\\begin{document}' + '\n')
	outputfile.write('\\begin{deluxetable}{ccc}' + '\n')
	outputfile.write('\\tablecolumns{3}' + '\n')
	outputfile.write('\\tablewidth{0pc}' + '\n')
	outputfile.write('\\tablecaption{' + str(Projectile) + ' ' + str(int_to_roman(int(Charge))) + ' Charge Exchange X-ray Emission with ' + str(Target) + ' Targets}' + '\n')
	outputfile.write('\\tablehead{\multicolumn{3}{c}{Collision Energy = ' + str(CollisionEnergy)+ ' eV/u $\\vert$ Collsion Velocity = ' + '%.1f' % float(convertspeed(float(CollisionEnergy))) + ' km/s}' +  ' \\' + '\\' + ' \colhead{Photon Energy (eV)} & \colhead{Emission Line} & \colhead{Line Ratio}}')
	outputfile.write('\startdata' + '\n')
	for i in range(0, len(tags)):
		outputfile.write(str(tags[i]) + '$   &   ' +'%.2f' % float(photenergy[indices[i]]) + '   &   ' +   '%.4f' % float(relative[indices[i]]) + '\t' + '\\' + '\\' + '\n')
	outputfile.write('\\' + '\\' + '\n')
	outputfile.write('\hline' + '\n')
	outputfile.write('\enddata' + '\n')
	outputfile.write('\end{deluxetable}' + '\n') 
	outputfile.write('\end{document}')
	outputfile.close()
else:
	outputfile = open('LineRatios.tex', 'w')
	outputfile.write('\documentclass{aastex}' + '\n')
	outputfile.write('\\begin{document}' + '\n')
	outputfile.write('\\begin{deluxetable}{ccc}' + '\n')
	outputfile.write('\\tablecolumns{3}' + '\n')
	outputfile.write('\\tablewidth{0pc}' + '\n')
	outputfile.write('\\tablecaption{' + str(Projectile) + ' ' + str(int_to_roman(int(Charge))) + ' Charge Exchange X-ray Emission with ' + str(Target) + ' Targets}' + '\n')
	outputfile.write('\\tablehead{\multicolumn{3}{c}{Collision Energy = ' + str(CollisionEnergy)+ ' eV/u $\\vert$ Collsion Velocity = ' + '%.1f' % float(convertspeed(float(CollisionEnergy))) + ' km/s}' +  ' \\' + '\\' + ' \colhead{Emission Line} & \colhead{Photon Energy (eV)} &  \colhead{Line Ratio}}')
	outputfile.write('\startdata' + '\n')
	for i in range(0, len(tags)):
		outputfile.write('%.2f' % float(photenergy) + '   &   ' + str(tags[i]) + '$   &   ' +  '%.4f' % float(relative) + '\t' + '\\' + '\\' + '\n')
	outputfile.write('\\' + '\\' + '\n')
	outputfile.write('\hline' + '\n')
	outputfile.write('\enddata' + '\n')
	outputfile.write('\end{deluxetable}' + '\n') 
	outputfile.write('\end{document}')
	outputfile.close()

os.system('pdflatex LineRatios.tex >/dev/null')
