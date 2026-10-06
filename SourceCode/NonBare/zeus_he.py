#! /usr/bin/env python
# Developed by: Patrick Mullen

###################################################################

# Necessary imports are made here.

from pylab import *
import argparse
import os

###################################################################

# This provides feedback to the user that the program has initialized.

print '...Creating Latex Figure and PDF'

parser=argparse.ArgumentParser()
parser.add_argument("-I", "--projectile", help="Specify Projectile Ion")
parser.add_argument("-Z", '--charge' ,help="Specify Projectile Charge")
parser.add_argument("-E", "--collisionenergy", help = "Specify the Collision Energy")
parser.add_argument("-A", "--target", help = "Specify the Target")
parser.add_argument("-F", "--filename", help = "Specify the Filename")
parser.add_argument("-M", "--method", help="Specify method ")

args=parser.parse_args()

CollisionEnergy = float(args.collisionenergy)
Projectile = str(args.projectile)
Charge = int(args.charge)
Target = str(args.target)
filename = str(args.filename)
method = str(args.method)

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
newlineslist = list()

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

initialn, initiall, initials, finaln, finall, finals, photenergy, relative = loadtxt('hephaistus.dat', dtype=str, unpack = True, usecols=[0,1,2,3,4,5,6,-1])

os.remove('hephaistus.dat')

initialn = initialn.tolist()
initiall = initiall.tolist()
initials = initials.tolist()
finaln = finaln.tolist()
finall = finall.tolist()
finals = finals.tolist()
photenergy = photenergy.tolist()
relative = relative.tolist()

initialn.reverse()
initiall.reverse()
initials.reverse()
finaln.reverse()
finall.reverse()
finals.reverse()
photenergy.reverse()
relative.reverse()

dict = {'210100':'K alpha r', '211100':'K alpha i' , '201100':'K alpha f','310100':'K beta  ', '410100':'K gamma ', '510100':'K delta ', '610100':'K epsilon ', '710100':'K zeta  ', '810100':'K eta   ', '910100':\
'K theta ', '1010100':'K iota  ', '1110100': 'K kappa ', '1210100': 'K lambda ', '1310100': 'K mu    ', '1410100': 'K nu    ', '1510100':'K xi    ', '1610100': 'K omicron',\
'1710100': 'K pi    ', '1810100': 'K rho   ', '1910100': 'K sigma', '2010100': 'K tau   '}

checklist = ['201100','211100','210100','310100', '410100', '510100', '610100', '710100', '810100', '910100'\
, '1010100', '1110100', '1210100', '1310100', '1410100', '1510100', '1610100',\
'1710100', '1810100', '1910100', '2010100']

newlist = list()

for i in range(0, len(finals)):
	newlist.append(initialn[i] + initiall[i] + initials[i] + finaln[i] + finall[i] + finals[i])

index = list()

for i in range(0, len(checklist)):
	index.append(newlist.index(checklist[i]) if checklist[i] in newlist else 'nothere')

finalindex = list()

for i in index:
	if i not in finalindex:
		finalindex.append(i)

finalindex.remove('nothere')

finaltags = list()

for i in range(0, len(finalindex)):
	finaltags.append(str(dict[newlist[finalindex[i]]]))

outputfile = open('LineRatios.dat', 'w')
outputfile.write('#Collision Energy = ' + str(CollisionEnergy) + ' eV/u   Collision Velocity = ' + '%.1f' % float(convertspeed(float(CollisionEnergy))) +'\n')
outputfile.write('#Photon Energy (eV) ' + '\t' + 'Emission Line' + '\t' + 'Line Ratio' + '\n')
for i in range(0, len(finalindex)):
	outputfile.write(str(photenergy[finalindex[i]]) + '\t' + str(finaltags[i]) + '\t' + str(relative[finalindex[i]]) + '\n')
outputfile.close()

for i in range(0, len(finaltags)):
	finaltags[i] = finaltags[i].replace('K ', 'K $\\')

outputfile = open('LineRatios.tex', 'w')
outputfile.write('\documentclass{aastex}' + '\n')
outputfile.write('\\begin{document}' + '\n')
outputfile.write('\\begin{deluxetable}{ccc}' + '\n')
outputfile.write('\\tablecolumns{3}' + '\n')
outputfile.write('\\tablewidth{0pc}' + '\n')
outputfile.write('\\tablecaption{'  + str(Projectile) + ' ' + str(int_to_roman(int(Charge))) + ' Charge Exchange X-ray Emission with ' + str(Target) + ' Targets using the '  + str(method) + ' ' + 'method }' + '\n')
outputfile.write('\\tablehead{\multicolumn{3}{c}{Collision Energy = ' + str(CollisionEnergy)+ ' eV/u $\\vert$ Collsion Velocity = ' + '%.1f' % float(convertspeed(float(CollisionEnergy))) + ' km/s}' +  ' \\' + '\\' + ' \colhead{Photon Energy (eV)} & \colhead{Emission Line} & \colhead{Line Ratio}}')
outputfile.write('\startdata' + '\n')
for i in range(0, len(finaltags)):
	outputfile.write('%.2f' % float(photenergy[finalindex[i]]) + '   &   ' + str(finaltags[i][0:-1]) + '$ ' + str(finaltags[i][-1]) + ' &   ' +  '%.4f' % float(relative[finalindex[i]]) + '\t' + '\\' + '\\' + '\n')
outputfile.write('\\' + '\\' + '\n')
outputfile.write('\hline' + '\n')
outputfile.write('\enddata' + '\n')
outputfile.write('\end{deluxetable}' + '\n')
outputfile.write('\end{document}')
outputfile.close()

os.system('pdflatex LineRatios.tex >/dev/null')
