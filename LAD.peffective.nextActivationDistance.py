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

NEregion = 50.0 #in nm
rNE = 2000.0 #in nm nuclear radius
getnum = re.compile(r'[^\d.]+')
#---------freq info of optimization step---------
#blocksize = 1000             #the number of lines read in probability file
targetfb = sys.argv[1]       #number in p% for  the next optimization
lastfName = sys.argv[2]       #last pymfile prefix
strucDir = sys.argv[3]       #pymfile directory
nstruct  = int(sys.argv[4])  #number of structures examined
beadfile = sys.argv[5]       #bead index file
outDir = sys.argv[6]         #out
runName = sys.argv[7]         #out
lad_p = sys.argv[8]         #current lad probability
#runName='.smallNucleol.DamID' #with dot first is better
pymname = lastfName+runName+'.pym'
ladprobfile = "/panfs/cmb-panasas2/qjl_001/Drosophila/code/modeling/fly_uniform/withDamID/lamDamID.prob"

targetfreq = float(getnum.sub('',lad_p))/100. #for DamID lamin freq
p100 = re.search("p100",lastfName)
if(p100):              #only after p100
    haveActDist = False
else:
    haveActDist = True
    plastfile = strucDir+'/'+lastfName+runName+".peffect_lad_ActDist"
fout=open(outDir+'/'+targetfb+"."+lad_p+runName+'.peffect_lad_ActDist','w')
#------------------------------------------------------------
tol = 0.999
t1=time.time()
#---------------------------------------------------------
def getBead(query,beadInfo):
    '''given a (chr, start, end) in query, returns the order within beadInfo'''
    c1 = query[0]
    p1 = 0.5*(query[1]+query[2])
    i = 0
    for c,s,e in beadInfo:
        if c == c1 and s<p1<e:
            break
        i += 1
    return i

def compute_distance(b1,b2,coor):
    '''Calculate distance between bead1 and bead2 from the coor[]'''
    vd = coor[b1]-coor[b2]
    dist = numpy.sqrt(numpy.dot(vd,vd))
    return dist

def readpym_diploid(pymfile,nhap=1169):
    '''Given pymfile returns xyzr'''
    fcoor = open(pymfile,'r')
    n = 2*nhap
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

def existingPortionNE(v, r=250.0):
    '''Given a vector of radial positions v, return a fraction that is within r'''
    inside = [s for s in v if s < r]
    n = len(v)*1.0
    return len(inside)/n

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
        #pwish = float(tmp[2])
        lastp[b1] = float(tmp[3])
    return lastp

#-------------------------------------------------
#READ the coarse bead info to find out intrachr
fsize = open(beadfile,'r')
bead_to_chr = []
ordered_beadInfo = []
for line in fsize.readlines()[1:]:
    cols = string.split(line)
    ch = cols[1]
    bead_to_chr.append(ch)
    sta = float(cols[2])
    end = float(cols[3])
    ordered_beadInfo.append((ch,sta,end))
fsize.close()


#Get the previous p-effective
if(haveActDist):              #only after p100
    plast = readP0(plastfile)
else:
    plast = {}

#registering the LAD indexes and P-wish
flad =open(ladprobfile, 'r') #this is 1-based indexing
targetPloci = []
for line in flad.readlines():
    cols = string.split(line)
    b = int(cols[0])-1 #translate to 0-based id
    p = float(cols[1]) #fraction in population
    try:
        p0 = plast[b]
    except KeyError:
        p0 = 0.0
    if targetfreq <= p < 1: #just include above p
        bp = (b,p,p0)
        targetPloci.append(bp)


distNE = defaultdict(list)
if(len(targetPloci) > 0):
    for i in xrange(nstruct):
        fname = "%s/ensemble/copy%d/%s"%(strucDir,i,pymname)
        #print "struc",i
        coor,r = readpym_diploid(fname)  #coordinates and radii
        nbead = len(r)/2
        for bpP0 in targetPloci:
            b1 = bpP0[0]
            b3 = b1 + nbead #the diploid copy
            vd = coor[b1]
            d1 = numpy.sqrt(numpy.dot(vd,vd))
            vd = coor[b3]
            d2 = numpy.sqrt(numpy.dot(vd,vd))
            r1 = r[b1]  #the radius of bead
            #collect position of 2 homolog copies
            distNE[bpP0].append(rNE-d1-r1) #distance from surface
            distNE[bpP0].append(rNE-d2-r1) #from surf

for k,v in distNE.items():
    s = sorted(v)
    b1 = k[0]
    pwish = k[1]
    plast = k[2]
    pnow = existingPortionNE(v,NEregion)
    t = cleanProbability(pnow,plast) #effective background now
    p = cleanProbability(pwish,t) #peffective iterative while doing gradual steps
    o=min((2*nstruct)-1,int(round(2*p*nstruct))) #total 2M positions
    fout.write('%4d %5.3f %7.1f %5.3f %5.3f\n'%(b1,pwish,s[o],p,pnow)) #beadID,experiment_p,dist,current_peffective,existed_p

#-------------------------
t2=time.time()
print 'time spend is ', t2-t1, ' s'
sys.exit()
