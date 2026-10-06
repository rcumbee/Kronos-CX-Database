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
args=parser.parse_args()

NISTFILE = args.CXFILE
q = float(args.collisionenergy)
charge = int(args.projcharge)

#######################################################################

dict = {'s':'0', 'p':'1', 'd':'2', 'f':'3', 'g':'4', 'h':'5', 'i':'6', 'k':'7', 'l':'8', 'm':'9', 'n':'10', 'o':'11', 'q':'12', 'r':'13', 't':'14', '1':'0', '3':'1'}

#######################################################################

energies = loadtxt(NISTFILE, dtype = float, unpack = True, usecols=[0])
closest = min(enumerate(energies), key=lambda x: abs(x[1]-q))

#######################################################################

print '\n', '...Proceeding with collision energy,', str(closest[1]), 'eV/u.', '\n'

#######################################################################

f = open(NISTFILE, "r")
for i in range(0, 7):
	f.readline()
x = str(f.readline())
for i in range(0, int(closest[0])):
	f.readline()
y = str(f.readline())
f.close()

#######################################################################

CXSTRINGS = x.split()
CROSSSECS = y.split()
CROSSSECSnorm = list(CROSSSECS)

CXSTRINGS.remove('#')
CXSTRINGS.remove('(eV/u)')

for i in range(0, len(CROSSSECS)):
	CROSSSECSnorm[i] = float(CROSSSECS[i])/float(CROSSSECS[-1])
	CROSSSECS[i] = float(CROSSSECS[i])/float(1.0)


del CROSSSECS[0]
del CROSSSECS[-1]
del CROSSSECSnorm[0]
del CROSSSECSnorm[-1]

nlist = list()
llist = list()
slist = list()
tobe = list()
ascent = list()
ascentcross = list()
ascentcrossnorm = list()
indexing = list()

for i in range(0, len(CXSTRINGS)-1):
	nlist.append(str(CXSTRINGS[i][0:-5]))
	llist.append(str(dict[CXSTRINGS[i][-5]]))
	slist.append(str(CXSTRINGS[i][-3]))

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

outputfile = open('inithe.dat', 'w')
for i in range(0, len(indexing)):
	outputfile.write(str(nlist[indexing[i]]) + '\t' + str(llist[indexing[i]]) + '\t' + str(dict[slist[indexing[i]]]) + '\t' + str(ascentcross[i])  + '\n')
outputfile.write(str(int(100)) + '\t' + str(0) + '\t' + str(0) + '\t' + '0.0d0' + '\t' + '0.0d0' + '\n')
outputfile.write(str(charge) + '\t' + str(closest[1]))
outputfile.close()

#######################################################################
