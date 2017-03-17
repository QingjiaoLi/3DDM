##This incremental blocks from .prob file takes care one block each run
#choose 2/4
#update probability peffect by excluding pairs that automatically in contact without restraints
#iterative update of p-effective
#take care of Nuc & Het extra
#for pwish < 1, since pwish=1 has been imposed in fmax_restraint
import datetime
import sys
import re
import string
import math
from StringIO import StringIO
import math
import time
import numpy
from collections import defaultdict

getnum = re.compile(r'[^\d.]+')
#---------freq info of optimization step---------
#blocksize = 1000             #the number of lines read in probability file
targetfb = sys.argv[1]       #number in p% for  the next optimization
lastfName = sys.argv[2]       #last pymfile prefix
probfile = sys.argv[3]       #probability filename
strucDir = sys.argv[4]       #pymfile directory
nstruct  = int(sys.argv[5])  #number of structures examined
lstart   = int(float(sys.argv[6]))  #readline start of probability file
beadfile = sys.argv[7]       #bead index file
blocksize = int(sys.argv[8]) #number of lines read
outDir = sys.argv[9]         #output directory
runName = sys.argv[10]       #model name
lad_p = sys.argv[11]         #current lad probability for filename fout

pymname = lastfName+runName+'.pym'
fout=open(outDir+'/'+targetfb+"."+lad_p+runName+'.peffect_ActDist.'+str(lstart),'w')
targetfreq = float(getnum.sub('',targetfb))/100.
p100 = re.search("p100",lastfName)
if(p100):              #only after p100
    haveActDist = False
else:
    haveActDist = True
    plastfile = strucDir+'/'+lastfName+runName+".peffect_ActDist" 
#------------------------------------------------------------
tol = 0.999
t1=time.time()
#---------------------------------------------------------
def compute_distance(b1,b2,coor):
    '''Calculate distance between bead1 and bead2 from the coor[]'''
    vd = coor[b1]-coor[b2]
    dist = numpy.sqrt(numpy.dot(vd,vd))
    return dist

def readpym_diplo(pymfile,nhap=1169):
    '''Given pymfile returns xyzr'''
    fcoor = open(pymfile,'r')
    n = 2*nhap #first n beads are diploid domains
    lastcoor = []
    radii = []
    for line in fcoor.readlines():
        if line[:7]=="SPHERE,":
            xyz=string.split(line,',')
            x = float(xyz[1])
            y = float(xyz[2])
            z = float(xyz[3])
            radii.append(float(xyz[4]))
            lastcoor.append(numpy.array([x,y,z]))
    fcoor.close()
    return lastcoor[:n],radii[:n]

def existingPortion(v, rsum):
    '''Given a vector of distances, return a fraction that is in contact'''
    #rsum is sum of 2 radii, and v is a list of surface-surface distance 
    incontact = [s for s in v if s < rsum]
    n = len(v)*1.0
    return len(incontact)/n

def cleanProbability(pij,pexist):
    '''Equation for effective P'''
    try:
        pclean = (pij-pexist)/(1.0-pexist) #can be 0 division
    except:
        pclean = pij #reset to factory setting
    #return pclean
    return max(0,pclean) #non-negative correcting matp

def readP0(fname):
    '''Given a filename, read the previous P-effective'''
    lastp = {}
    for line in open(fname).readlines():
        tmp = line.split()
        b1 = int(tmp[0])
        b2 = int(tmp[1])
        #pwish = float(tmp[2])
        lastp[(b1,b2)] = float(tmp[4])
    return lastp

#-------------------------------------------------
#READ the coarse bead info to find out intrachr
fsize = open(beadfile,'r')
bead_to_chr = []
for line in fsize.readlines()[1:]:
    cols = string.split(line)
    ch = cols[1]
    bead_to_chr.append(ch)
fsize.close()


#Get the previous p-effective
if(haveActDist):              #only after p100
    plast = readP0(plastfile)
else:
    plast = {}

#registering the pairs and P-wish
targetpairs = []
for line in open(probfile).readlines()[lstart:lstart+blocksize]:
    tmp = line.split()
    b1 = int(tmp[0])
    b2 = int(tmp[1])
    p = float(tmp[2])
    b = (b1,b2)
    try:
        p0 = plast[b]
    except KeyError:
        p0 = 0.0
    if targetfreq <= p < 1:
        bp = (b1,b2,p,p0)
        if bead_to_chr[b1] != bead_to_chr[b2]:
            targetpairs.append(bp)
        elif abs(b2-b1) > 1:
            targetpairs.append(bp)


distances = defaultdict(list)
for i in xrange(nstruct):
    fname = "%s/ensemble/copy%d/%s"%(strucDir,i,pymname)
    coor,r = readpym_diplo(fname)  #coordinates and radii
    nbead = len(r)/2
    for b1b2c in targetpairs:
        b1 = b1b2c[0]
        b2 = b1b2c[1]
        b3 = b1 + nbead
        b4 = b2 + nbead
        vd = coor[b1]-coor[b2]
        d1 = numpy.sqrt(numpy.dot(vd,vd))-r[b1]-r[b2]
        vd = coor[b3]-coor[b4]
        d2 = numpy.sqrt(numpy.dot(vd,vd))-r[b1]-r[b2]
        vd = coor[b1]-coor[b4]
        d3 = numpy.sqrt(numpy.dot(vd,vd))-r[b1]-r[b2]
        vd = coor[b3]-coor[b2]
        d4 = numpy.sqrt(numpy.dot(vd,vd))-r[b1]-r[b2]
        alld = sorted([d1,d2,d3,d4])
        if bead_to_chr[b1] != bead_to_chr[b2]:
            distances[b1b2c].append(alld[0])
            distances[b1b2c].append(alld[1])
        else:
            distances[b1b2c].append(d1)
            distances[b1b2c].append(d2)

radii = r

for k,v in distances.items():
    s = sorted(v)
    b1 = k[0]
    b2 = k[1]
    pwish = k[2]
    plast = k[3]
    pnow = existingPortion(v,radii[b1]+radii[b2])
    t = cleanProbability(pnow,plast) #effective background now
    #if (t != pwish):
    #    p = cleanProbability(pwish,t) #iterative while doing gradual steps
    #else:
    #    p = pwish
    p = cleanProbability(pwish,t) #iterative while doing gradual steps
    if p>0:
        o=min((2*nstruct)-1,int(round(2*p*nstruct))) #intrachr diploid assumption 2v instead of v
        fout.write('%4d %4d %5.3f %7.1f %5.3f %5.3f\n'%(b1,b2,pwish,s[o],p,pnow))

#-------------------------
t2=time.time()
#print 'time spend is ', t2-t1, ' s'
