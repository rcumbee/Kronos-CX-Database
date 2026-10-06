#! /usr/bin/env python
# Developed by: Patrick Mullen
#######################################################################

from pylab import *
import argparse

#######################################################################

parser=argparse.ArgumentParser()
parser.add_argument("-q", "--collisionenergy", help="Specify Projectile Ion")
parser.add_argument("-charge", '--projcharge' ,help="Specify Projectile Charge")
parser.add_argument("-i", "--CXFILE", help = "Specify the CXFILE")
parser.add_argument("-n", "--NISTFILE", help = "Specify the NISTFILE")
parser.add_argument("-rmax", "--rmax", help= "Specify rmax -- Optional", required=False)

args=parser.parse_args()

CXFILE = args.CXFILE
NISTFILE = args.NISTFILE
rmax = args.rmax
q = float(args.collisionenergy)
charge = int(args.projcharge)

energies = loadtxt(CXFILE, dtype = float, unpack = True, usecols=[0])
closest = min(enumerate(energies), key=lambda x: abs(x[1]-q))
excites = loadtxt(NISTFILE, dtype = float, unpack = True, usecols=[-1])

#######################################################################

print '\n', '...Proceeding with collision energy,', str(closest[1]), 'eV/u.', '\n'

#######################################################################

dict = {'0':'s', '1':'p', '2':'d', '3':'f', '4':'g', '5':'h', '6':'i', '7':'k', '8':'l', '9':'m', '10':'n', '11':'o', '12':'q', '13':'r', '14':'t'}

#######################################################################

f = open(CXFILE, "r")
for i in range(0, 7):
	f.readline()
x = str(f.readline())
for i in range(0, int(closest[0])):
	f.readline()
y = str(f.readline())
f.close()

#######################################################################

CXSTRINGS = x.replace(",", "      ")
CXSTRINGS = CXSTRINGS.split()

CROSSSECS = y.split()
CROSSSECSnorm = list(CROSSSECS)

CXSTRINGS.remove('#')
CXSTRINGS.remove('(eV/u)')
CXSTRINGSUPDATE = list()

for i in range(0, len(CXSTRINGS)):
	if i == 0 or i % 2 == 0:
		CXSTRINGSUPDATE.append(CXSTRINGS[i])
	else:
		CXSTRINGSUPDATE.append(dict[CXSTRINGS[i]])

del CXSTRINGS[-1]

#######################################################################

almostlist = list()

counter = 0
while counter < len(CXSTRINGSUPDATE)-1:
	almostlist.append(str(CXSTRINGSUPDATE[counter]) + str(CXSTRINGSUPDATE[counter + 1]))
	counter += 2

finallist = list()

for i in range(0, len(almostlist)):
	finallist.append(str(almostlist[i]) + '(2' + str(almostlist[i][-1]).capitalize() + ')')

for i in range(0, len(CROSSSECS)):
	CROSSSECSnorm[i] = float(CROSSSECS[i])/float(CROSSSECS[-1])
	CROSSSECS[i] = float(CROSSSECS[i])/float(1.0)

del CROSSSECS[0]
del CROSSSECS[-1]
del CROSSSECSnorm[0]
del CROSSSECSnorm[-1]

dict = {'s':'0', 'p':'1', 'd':'2', 'f':'3', 'g':'4', 'h':'5', 'i':'6', 'k':'7', 'l':'8', 'm':'9', 'n':'10', 'o':'11', 'q':'12', 'r':'13', 't':'14', '1':'0', '3':'1'}

nlist = list()
llist = list()
slist = list()
tobe = list()
ascent = list()
ascentcross = list()
ascentcrossnorm = list()
indexing = list()

for i in range(0, len(finallist)):
	nlist.append(str(finallist[i][0:-5]))
	llist.append(str(dict[finallist[i][-5]]))
	slist.append(str(finallist[i][-3]))

for i in range(0, len(nlist)):
	tobe.append(int(nlist[i] + llist[i] + slist[i]))
	ascent.append(int(nlist[i] + llist[i] + slist[i]))

ascent.sort()
ascentcross = [x for (y,x) in sorted(zip(tobe,CROSSSECS))]
ascentcrossnorm = [x for (y,x) in sorted(zip(tobe,CROSSSECSnorm))]

for i in range(0, len(tobe)):
	indexing.append(tobe.index(ascent[i]))

print '...Sum of normalized cross sections: ', sum(ascentcross)

#######################################################################

outputfile = open('closestenergy.dat', 'w')
outputfile.write(str(closest[1]))
outputfile.close()

#######################################################################

outputfile = open('inith.dat', 'w')
for i in range(0, len(indexing)):
	outputfile.write(str(nlist[indexing[i]]) + '\t' + str(llist[indexing[i]]) + '\t' + str(ascentcross[i]) +  '\t' + str(ascentcrossnorm[i]) + '\n')
outputfile.write(str(int(100.0)) + '\t' + '0' + '\t'  + '0' + '\t' + '0.0d0' + '\t' + '0.0d0'  + '\n')
outputfile.write(str(charge) + '\t' + str(rmax) + '\t' + str(closest[1]) + '\t' + str(excites[-1] * float(-1.23984E-4)))
outputfile.close()


#######################################################################
