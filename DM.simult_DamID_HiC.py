# Lamin-TAD + HiC restraints;
# Lamin usually don't have 100% value, so lad_p=l100 will generate population with only Hi-C
# p=1 means 2 out of 4 possible pairwise are considered as restraints; for intra p=1 is intra_both
# step by step contact restraints application
# HET is attached to each pericentromeric arm, then CEN is added for each chr; 2*6 + 4 beads
# homolog pairing applies to centromeres, no need to enforce in HET
import IMP
import IMP.core
import IMP.atom
import IMP.display
import IMP.algebra
import datetime
import time
import os
import string
import random
import sys
import math
import itertools
import numpy
import re
from operator import itemgetter
from collections import defaultdict

getnum = re.compile(r'[^\d.]+')
eps = 1.e-4 #small fraction
offset=0.999  #to get rid of =
########### modify these block every stepwise ##############
#currentfb = 'p010' #output name
RunName='.DamID_HiC' #with dot first is better
#lastfb = 'p050'
############################################################
pymname=RunName+'.pym'

'''##########INPUT ARGS############'''
workdir = sys.argv[1]+"/"
inputDir = sys.argv[2]+"/" #need full path
modelName = sys.argv[3]
currentfb = sys.argv[4] #e.g. 'p100' as output name & stage; the last 3 characters must be percentage digits!!!
targetfreq = float(getnum.sub('',currentfb))/100.
p100 = re.search("p100",currentfb)
if(p100):
    useLastPopulation = False
    lad_p='l100' # no lad restraints added
else:
    useLastPopulation = True

if(useLastPopulation):
    lastfb = sys.argv[5]  #must input the last coordinate filename prefix, e.g. p100
    pymdir = sys.argv[6]+"/" #need full path
    lad_p = sys.argv[7] #probability level of lamin-TAD, just a string to read filename
    lad_lastp = sys.argv[8] #last probability level of lamin-TAD, just a string to read filename
    #lastpym = pymdir+lastfb+pymname #input coordinate file name
    lastpym = "%s%s.%s%s"%(pymdir,lastfb,lad_lastp,pymname) #input coordinate file name
    hlistfile = workdir+currentfb+'.'+lad_p+RunName+'.peffect_ActDist'
    ladActivationDist = workdir+currentfb+"."+lad_p+RunName+'.peffect_lad_ActDist'
    #print "Last structure: ",lastpym
    print "The following distance threshold files are used:\n %s\n %s\n"%(hlistfile, ladActivationDist)
'''###########################################'''

ktouch = 1.e-4  # touching spring constant 1e-5 is good
xradsum = 1.0  #(1+xradsum)*radii_sum center-center
# Reading model setup
size_file = inputDir+modelName+'.ind'
probability_file = inputDir+modelName+'.from40k.prob'
damID_probfile = inputDir+"lamDamID.prob"
print "Optimizing at Hi-C prob %.2f with input:\n %s \n %s"%(targetfreq,size_file,probability_file)
scalestep = 1/math.sqrt(targetfreq)
#Drosophila chromosomes to be modeled
chrname = ['2L','2R','3L','3R','4','X']

sys.stdout.flush()
#inactiveX = []   # list of bead_id for chrX
onethird = 1./3.
# WARNING: bead_id starts from 0 here, but 1 in external input files
bead_bp = []  # bp in each bead
Rb = []    # Bead radius
interact = {}
#touchRestraints = [] # touching spheres restraints

begin_chromosome = []
nbead_chr = []  # #of_beads in chr
bead_to_chr = []
bead_to_hires = [] #(i138_start,i138_end) for each bead, important for a_ij
ordered_beadInfo = []

fsize = open(size_file,'r')
chr0=''
#arm0='L'
countbead = 0
subcentro = []
#actina = []
for line in fsize.readlines()[1:]:
    cols = string.split(line)
    ch = cols[1]
    arm = list(ch)[-1]
    #arm0 = arm
    bead_to_chr.append(ch)
    sta = float(cols[2])
    end = float(cols[3])
    ordered_beadInfo.append((ch,sta,end))
    bp = (end-sta)
    bead_bp.append(bp)
    if ch != chr0:
        begin_chromosome.append(countbead)
        chr0 = ch
    countbead = countbead+1
fsize.close()
#counting beads for each chromosome
nbead_chr = [bead_to_chr.count(x) for x in chrname]
ndomain = sum(nbead_chr) #nhaploid 
nhaploid = ndomain
end_chromosome = [i+j-1 for (i,j) in zip(begin_chromosome,nbead_chr)]
at_end = [0,2,5] #select chromosomes which centro is at end
at_first = [1,3,4] #select chromosomes which centro is at first idx
subcentro = sorted([end_chromosome[i] for i in at_end]+\
  [begin_chromosome[i] for i in at_first])
subcentro = subcentro+[ndomain + i for i in subcentro]
cenX = [end_chromosome[5],ndomain + end_chromosome[5]] #bead of cenX

# Constants
rscale = 1.35 #occupancy 15% including Heterochromatin
rscale = 1.38 #occupancy 15% including Heterochromatin
cenhet_scale = 0.05 #let centromere is 1% of het volume
#rscale = 1.4 #
rad_nucleus = 2000.0 # nm
laminR = rad_nucleus-50.0 #lowerbound distance from center to place LAD near NE
rad_hetroClus = rad_nucleus/3.0 #nm; #about 1/3 rad_nucleus
rad_nucleolus = 0.5*rad_hetroClus #nm; observation from paper figures
hetsizes = [5.4, 11.0, 8.2, 8.2, 3.1, 20.0] #relative het size following chrname order
hetsum = sum(hetsizes)
hetrads = [rad_hetroClus*(0.5*h/hetsum)**onethird for h in hetsizes] #0.5 because of diploid share
hetvol = [h**3 for h in hetrads]
hetvol_chr = [hetvol[0]+hetvol[1],hetvol[2]+hetvol[3],hetvol[4],hetvol[5]]
#print 2*sum(hetvol),rad_hetroClus**3
cenrads = [(cenhet_scale*h)**onethird for h in hetvol_chr]
cenvol = [h**3 for h in cenrads]


#hetroClus_center = IMP.algebra.Vector3D((rad_nucleus-rad_hetroClus),0,0) #on x-axis like yeast
#print 'Position of HET',hetroClus_center
cdensity = 107.45 # bp/nm assuming 197 bp/nucleosomes and 6nucleosomes/11nm
r_n2n = {} # sqrt<r^2> end-to-end distance in nm
kscale = (0.75*15**2)**onethird # 3/4*r**2 where r=15nm
rad_randomCap = 1*rad_nucleus # nm
r_n2n = {} # sqrt<r^2> end-to-end distance in nm
#rscale = 1.105928 #occupancy 10% hard volume nuc occupancy

totvol = 0.
for i in xrange(ndomain):
    Lc = bead_bp[i]/cdensity # in nm
    rad_bead = rscale*kscale*Lc**onethird
    Rb.append(rad_bead) # nm
    totvol = totvol + 4.*3.1415/3.*rad_bead**3
    #print i,bead_bp[i]/1.e6,rad_bead
dnavol = totvol*2
dnahet = dnavol+2*4*math.pi/3*(sum(hetvol)+sum(cenvol))
totvol = dnahet + 4*math.pi/3*(rad_nucleolus**3)
nucvol = (4*math.pi/3)*rad_nucleus**3
dnaocc = dnavol/nucvol
dnahet = dnahet/nucvol
print 'occupancy: %.2f euchromatin, %.2f with het, and %.2f with Nucleolus in Rnuc %d'%(dnaocc, dnahet, totvol/nucvol,rad_nucleus)
sys.stdout.flush()


#diploid Rb; 2xtotal haploid beads 
Rb=Rb+Rb

# Read the hires probability
fpread = open(probability_file,'r')
for line in fpread.readlines():
    cols = string.split(line)
    b1 = int(cols[0]) #prob file is 0-based
    b2 = int(cols[1])
    pread = float(cols[2])
    #if pread > 0:#remove this is safe as prob file contains nonzero column3
    #    interact[(b1,b2)] = pread
    interact[(b1,b2)] = pread
fpread.close()
t1=time.time()
#___________________________ IMP starts ___________________
m = IMP.Model()
chain = IMP.container.ListSingletonContainer()
IMP.set_check_level(IMP.NONE)
IMP.set_log_level(IMP.SILENT)
#----------------------------------------------------------  

#### FORMAT order of beads : 1169+1169+6+6+4+4+1 (2*domains+2*HET+2*CEN+nucleolus)
nhet = len(hetrads) #each haploid
ncen = len(cenrads) #each haploid
hetstart = 2*ndomain
censtart = hetstart + 2*nhet
hetbeads = [2*nhaploid + i for i in range(2*nhet)] #0-based index pointer to HETs
cenbeads = [2*(nhaploid + nhet) + i for i in range(2*ncen)] #0-based index pointer to CEN
nucbead = 2*(nhaploid + nhet + ncen) #0-based index of nucleolus
for h in hetrads:
    Rb.append(h)
for h in hetrads:
    Rb.append(h)
for h in cenrads:
    Rb.append(h)
for h in cenrads:
    Rb.append(h)
Rb.append(rad_nucleolus)

totalBead = len(Rb)
#NOW Rb is complete explicit radii for [2*chromosomes+2*6het+2*4cen+nucleolus]

#-----------------------------------------------------------------------------------
def timespend(t0):
    t1 = time.time()
    return int(round(t1-t0,0))

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

def mdstep(t,step):
    xyzr = chain.get_contained_particles()
    o = IMP.atom.MolecularDynamics()
    o.set_model(m)
    md = IMP.atom.VelocityScalingOptimizerState(xyzr,t,10)  # replace 300 K with 500 K
    o.add_optimizer_state(md)
    #print 'optimizing with temperature',t,'and',step,'steps'
    s=o.optimize(step)
    o.remove_optimizer_state(md)
    #print 'MD',step,'steps done @',datetime.datetime.now()
    return s

def cgstep(step):
    o = IMP.core.ConjugateGradients()
    o.set_model(m)
    s=o.optimize(step)
    #print 'CG',step,'steps done @',datetime.datetime.now()
    return s

def simulated_anneal(hot,cold,nc=10,nstep=500):
    """perform a cycle of simulated annealing from hot to cold"""
    dt = (hot-cold)/nc
    for i in range(nc):
        t = hot-dt*i
        mdstep(t,nstep)
    mdstep(cold,300)
    cgstep(100)

def simulated_anneal_scored(hot,cold,nc=10,nstep=500, lowscore=10):
    """perform a cycle of simulated annealing but stop if reach low score"""
    dt = (hot-cold)/nc
    for i in range(nc):
        t = hot-dt*i
        mdstep(t,nstep)
        score = cgstep(100)
        if score < lowscore:
            return t,score
            break
    mdstep(cold,300)
    score = cgstep(100)
    return t,score

def compute_distance(b1,b2):
    '''Caculate surface-to-surface distance between bead1 and bead2'''
    p1 = IMP.core.XYZR(chain.get_particle(b1))
    p2 = IMP.core.XYZR(chain.get_particle(b2))
    checkdist = IMP.core.get_distance(p1,p2)
    return checkdist

def euclideanDist(p1,p2):
    '''Just point to point distance'''
    #p1 and p2 must be [x,y,z]
    d = numpy.sqrt(sum([(x1-x2)**2 for x1,x2 in zip(p1,p2)]))
    return d

def is_inCap(b,rcap=rad_nucleus,rbuf=20):
    '''Check if bead center is within a radius from nucleus center'''
    p = IMP.core.XYZR(chain.get_particle(b))
    r = euclideanDist([p.get_x(),p.get_y(),p.get_z()],[0,0,0])
    ####if rcap > (r + p.get_radius()+rbuf):
    if rcap > (r):
        return True
    else:
        return False

def is_outCap(b,rcap=rad_nucleus):
    '''Check if bead center is within a radius from nucleus center'''
    p = IMP.core.XYZR(chain.get_particle(b))
    r = euclideanDist([p.get_x(),p.get_y(),p.get_z()],[0,0,0])
    if rcap < r:
        return True
    else:
        return False

def outCapScore(b,rcap=rad_nucleus, k=1.0):
    '''Check if bead center is within a radius from nucleus center'''
    p = IMP.core.XYZR(chain.get_particle(b))
    r = euclideanDist([p.get_x(),p.get_y(),p.get_z()],[0,0,0])
    r = r+(p.get_radius()) #reach the surface of this bead
    if rcap < r:
        d = r-rcap
    else:
        d = 0.0
    return(k*d*d/2.0)

def inCapScore(b,rcap=rad_nucleus, k=1.0, rbuf=5.0):
    '''Check if bead center is within a radius from nucleus center'''
    p = IMP.core.XYZR(chain.get_particle(b))
    r = euclideanDist([p.get_x(),p.get_y(),p.get_z()],[0,0,0])
    r = r+(p.get_radius()*(1+rbuf/100)) #reach the surface of this bead
    if rcap < r:
        d = 0.0
    else:
        d = r-rcap
    return(k*d*d/2.0)

def lastcoorDistance(beadpair):
    '''Caculate surface-to-surface distance between a beadpair in the lastcoor'''
    b1 = beadpair[0]
    b2 = beadpair[1]
    vd = lastcoor[b1]-lastcoor[b2]
    dist = numpy.sqrt(numpy.dot(vd,vd))-Rb[b1]-Rb[b2]
    return dist #in nm

def lastcoorRpos(b):
    '''Returns radial position of bead b from the lastcoor'''
    vd = lastcoor[b]
    dist = numpy.sqrt(numpy.dot(vd,vd))
    return dist #in nm

def pym(filename):
    '''get the pym file'''
    pym2 = IMP.display.PymolWriter(filename)
    g3 = IMP.display.XYZRsGeometry(chain,IMP.core.XYZR.get_default_radius_key())
    g3.set_name("beads")
    g3.set_color(IMP.display.Color(1,1,1))
    pym2.add_geometry(g3)

def inDiameterRestraint(beadset,diameter, kspring=1):
    lc = IMP.container.ListSingletonContainer()
    for b in beadset:
        lc.add_particle(chain.get_particle(b))
    hu = IMP.core.HarmonicUpperBound(0,kspring)
    dr = IMP.core.DiameterRestraint(hu,lc,diameter)
    m.add_restraint(dr)

def intra_both(b1,b2,dscale, kspring=1,nbead=ndomain):
    '''give f100 restraint to both copies'''
    upperdist = dscale*(Rb[b1]+Rb[b2])
    ds = IMP.core.SphereDistancePairScore\
       (IMP.core.HarmonicUpperBound(upperdist,kspring))
    dst = IMP.core.SphereDistancePairScore\
       (IMP.core.HarmonicUpperBound(.1,ktouch))
    p1 = chain.get_particle(b1)
    p2 = chain.get_particle(b2)
    apair = IMP.ParticlePair(p1,p2)
    pr = IMP.core.PairRestraint(ds,apair)
    prt = IMP.core.PairRestraint(dst,apair)
    m.add_restraint(pr)
    #m.add_restraint(prt)
    #touchRestraints.append(prt)
    c1 = b1 + nbead
    c2 = b2 + nbead
    p1 = chain.get_particle(c1)
    p2 = chain.get_particle(c2)
    apair = IMP.ParticlePair(p1,p2)
    pr = IMP.core.PairRestraint(ds,apair)
    prt = IMP.core.PairRestraint(dst,apair)
    m.add_restraint(pr)
    #m.add_restraint(prt)
    #touchRestraints.append(prt)
    #print "Add fmax",b1,"-",b2
    inbothRestraints.append([b1,b2])

def restrain_one(b1,b2,dscale, kspring=1.0):
    '''give f100 restraint to one pair'''
    upperdist = dscale*(Rb[b1]+Rb[b2])
    ds = IMP.core.SphereDistancePairScore\
       (IMP.core.HarmonicUpperBound(upperdist,kspring))
    dst = IMP.core.SphereDistancePairScore\
       (IMP.core.HarmonicUpperBound(.1,ktouch))
    p1 = chain.get_particle(b1)
    p2 = chain.get_particle(b2)
    apair = IMP.ParticlePair(p1,p2)
    pr = IMP.core.PairRestraint(ds,apair)
    prt = IMP.core.PairRestraint(dst,apair)
    m.add_restraint(pr)
    #m.add_restraint(prt)
    #touchRestraints.append(prt)
    #print "Add fmax",b1,"-",b2

def minpair(lpair,to_n,dscale, kspring=1.0):
    '''give restraint to to_n=1,2,3,or 4 of lpair possible pairs'''
    ambi = IMP.container.ListPairContainer()
    for p in lpair:
        pi = chain.get_particle(p[0])
        pj = chain.get_particle(p[1])
        pair = IMP.ParticlePair(pi,pj)
        ambi.add_particle_pair(pair)
    upperdist = dscale*(Rb[p[0]]+Rb[p[1]])
    ds = IMP.core.SphereDistancePairScore\
    (IMP.core.HarmonicUpperBound(upperdist,kspring))
    dst = IMP.core.SphereDistancePairScore\
    (IMP.core.HarmonicUpperBound(.1,ktouch))
    minpr = IMP.container.MinimumPairRestraint(ds,ambi,to_n)
    minprt = IMP.container.MinimumPairRestraint(dst,ambi,to_n)
    m.add_restraint(minpr)
    #m.add_restraint(minprt)
    #touchRestraints.append(minprt)
    minpairRestraints.append([lpair,to_n])


def fmax_restraints4out4(dscale,rnum,duptcc, intercount=0, nbead=ndomain):
    '''assign restraints for prob=1.0'''
    linear_dict = {} #will be sorted by distance
    for k,v in interact.items():
        c1 = k[0]
        c2 = k[1]
        c3 = c1+nbead
        c4 = c2+nbead
        if v > 0.999:
            linear_dict[(c1,c2,v)] = abs(c2-c1)
    listsort = sorted(linear_dict.items(),key=itemgetter(1),reverse=False)
    for i,j in listsort:
        c1 = i[0]
        c2 = i[1]
        c3 = c1+nbead
        c4 = c2+nbead
        bpair = [(c1,c2),(c1,c4),(c2,c3),(c3,c4)]
        if (bead_to_chr[c1] == bead_to_chr[c2]): #INTRA
            if abs(c1-c2)>1: # only none consecutive
                minpair(bpair,4,dscale)
                rnum=rnum+1
                duptcc=duptcc+3
                pair_enforced[(c1,c2)]=['all4',dscale*(Rb[c1]+Rb[c2])]
        else:                                    #INTER
            minpair(bpair,4,dscale)
            pair_enforced[(c1,c2)]=['all4_inter',dscale*(Rb[c1]+Rb[c2])]
            rnum=rnum+1
            duptcc=duptcc+4
            intercount += 1

    return rnum,duptcc,intercount

def fmax_restraints(dscale, rnum, duptcc, intercount=0, nbead=ndomain):
    '''assign restraints for prob=1.0'''
#i.e.: restrain_one, intra_both, minpair
    linear_dict = {} #will be sorted by distance
    for k,v in interact.items():
        c1 = k[0]
        c2 = k[1]
        c3 = c1+nbead
        c4 = c2+nbead
        if v > 0.999:
            linear_dict[(c1,c2,v)] = abs(c2-c1)
    listsort = sorted(linear_dict.items(),key=itemgetter(1),reverse=False)
    for i,j in listsort:
        c1 = i[0]
        c2 = i[1]
        c3 = c1+nbead
        c4 = c2+nbead
        bpair = [(c1,c2),(c1,c4),(c2,c3),(c3,c4)]
        if (bead_to_chr[c1] == bead_to_chr[c2]): #INTRA
            if abs(c1-c2)>1: # only none consecutive
                intra_both(c1,c2,dscale)
                rnum=rnum+1
                duptcc=duptcc+1
                pair_enforced[(c1,c2)]=['min2',dscale*(Rb[c1]+Rb[c2])]
        else:                                    #INTER
            minpair(bpair,2,dscale)
            pair_enforced[(c1,c2)]=['min2_inter',dscale*(Rb[c1]+Rb[c2])]
            rnum=rnum+1
            duptcc=duptcc+1
            intercount += 1
    return rnum,duptcc, intercount


def apply_contact(contact_data,dscale,rnum,duptcc,intercount=0, nbead=ndomain):
    '''Based on the last conformations, 1-2 possible pairs can be selected'''
    fhlist = open(contact_data,'r')
    intraTargetDistances = {} #will be applied contacts
    interTarget = {} #will be applied contacts
    for line in fhlist.readlines():
        tmp=string.split(line)
        dcutoff = float(tmp[3])
        bprob = float(tmp[2])
        c1 = int(tmp[0]); c2 = int(tmp[1])
        c3 = c1+nbead; c4 = c2+nbead
        bpair = [(c1,c2),(c1,c4),(c2,c3),(c3,c4)]
        lastdist1 = lastcoorDistance(bpair[0])
        lastdist2 = lastcoorDistance(bpair[3])
        lastdist3 = lastcoorDistance(bpair[1])
        lastdist4 = lastcoorDistance(bpair[2])
        if (bead_to_chr[c1] == bead_to_chr[c2]): #INTRACHR
            min2 = min(lastdist1,lastdist2)
            max2 = max(lastdist1,lastdist2)
            dcross = min(lastdist3, lastdist4)
            if max2 <= dcutoff:
                b = (c1,c2,2,bprob) #both pairs will be in contact
                intraTargetDistances[b] = abs(c2-c1)
                pair_enforced[(c1,c2)]=['min2',dscale*(Rb[c1]+Rb[c2])]
            elif min2 <= dcutoff:
                b = (c1,c2,1,bprob) #only 1 pair will be in contact
                intraTargetDistances[b] = abs(c2-c1)
                pair_enforced[(c1,c2)]=['min1',dscale*(Rb[c1]+Rb[c2])]
            #else none will be in contact
        else: #INTERCHR
            b = (c1,c2,bprob)
            alld = sorted([lastdist1,lastdist2,lastdist3,lastdist4])
            min2 = alld[1] #second minimum
            min4 = alld[0] #minimum all
            if min2 <= dcutoff:
                interTarget[b] = 2 #
                pair_enforced[(c1,c2)]=['min2_inter',dscale*(Rb[c1]+Rb[c2])]
            elif min4 <= dcutoff:
                interTarget[b] = 1 #
                pair_enforced[(c1,c2)]=['min1_inter',dscale*(Rb[c1]+Rb[c2])]
    fhlist.close()
    listsort = sorted(intraTargetDistances.items(),key=itemgetter(1),reverse=False)
    #sort the intra based on distance first
    for i,j in listsort:
        c1 = i[0]
        c2 = i[1]
        c3=c1+nbead
        c4=c2+nbead
        nchoose = i[2]
        sqp = 10*math.sqrt(i[3])
        if nchoose > 1:
            intra_both(c1,c2,dscale, kspring=sqp)
            duptcc=duptcc+1
        else:
            ambiPairList=[(c1,c2),(c3,c4)]
            minpair(ambiPairList,1,dscale,kspring=sqp)
        rnum=rnum+1
    for bpair,nchoose in interTarget.items():
        c1 = bpair[0]; c2 = bpair[1]
        c3 = c1+nbead; c4 = c2+nbead
        ambiPairList = [(c1,c2),(c1,c4),(c2,c3),(c3,c4)]
        minpair(ambiPairList,nchoose,dscale,kspring=sqp)
        duptcc=duptcc+nchoose-1
        rnum=rnum+1
        intercount += 1

    return rnum,duptcc,intercount

def apply_contact4out4(contact_data,dscale,rnum,duptcc,intercount=0, nbead=ndomain):
    '''Based on the last conformations, 0-4 possible pairs can be selected'''
    fhlist = open(contact_data,'r')
    intraTargetDistances = {} #will be applied contacts
    interTarget = {} #will be applied contacts
    for line in fhlist.readlines():
        tmp=string.split(line)
        dcutoff = float(tmp[3])
        bprob = float(tmp[2])
        c1 = int(tmp[0]); c2 = int(tmp[1])
        c3 = c1+nbead; c4 = c2+nbead
        bpair = [(c1,c2),(c1,c4),(c2,c3),(c3,c4)]
        lastdist1 = lastcoorDistance(bpair[0])
        lastdist2 = lastcoorDistance(bpair[3])
        lastdist3 = lastcoorDistance(bpair[1])
        lastdist4 = lastcoorDistance(bpair[2])
        alld = sorted([lastdist1,lastdist2,lastdist3,lastdist4])#low to high
        if (bead_to_chr[c1] == bead_to_chr[c2]): #INTRACHR
            if alld[3] < dcutoff:
                b = (c1,c2,4,bprob) #all pairs will be in contact
                intraTargetDistances[b] = abs(c2-c1)
                pair_enforced[(c1,c2)]=['all4',dscale*(Rb[c1]+Rb[c2])]
            elif alld[2] < dcutoff:
                b = (c1,c2,3,bprob) #3pairs will be in contact
                intraTargetDistances[b] = abs(c2-c1)
                pair_enforced[(c1,c2)]=['min3',dscale*(Rb[c1]+Rb[c2])]
            elif alld[1] < dcutoff:
                b = (c1,c2,2,bprob) #both pairs will be in contact
                intraTargetDistances[b] = abs(c2-c1)
                pair_enforced[(c1,c2)]=['min2',dscale*(Rb[c1]+Rb[c2])]
            elif alld[0] < dcutoff:
                b = (c1,c2,1,bprob) #only 1 pair will be in contact
                intraTargetDistances[b] = abs(c2-c1)
                pair_enforced[(c1,c2)]=['min1',dscale*(Rb[c1]+Rb[c2])]
            #else none will be in contact
        else: #INTERCHR
            b = (c1,c2,bprob)
            if alld[3] < dcutoff:
                interTarget[b] = 4 #
                pair_enforced[(c1,c2)]=['all4_inter',dscale*(Rb[c1]+Rb[c2])]
            elif alld[2] < dcutoff:
                interTarget[b] = 3 #
                pair_enforced[(c1,c2)]=['min3_inter',dscale*(Rb[c1]+Rb[c2])]
            elif alld[1] < dcutoff:
                interTarget[b] = 2 #
                pair_enforced[(c1,c2)]=['min2_inter',dscale*(Rb[c1]+Rb[c2])]
            elif alld[0] < dcutoff:
                interTarget[b] = 1 #
                pair_enforced[(c1,c2)]=['min1_inter',dscale*(Rb[c1]+Rb[c2])]
    fhlist.close()
    listsort = sorted(intraTargetDistances.items(),key=itemgetter(1),reverse=False)
    #sort the intra based on distance first
    for i,j in listsort:
        c1 = i[0]
        c2 = i[1]
        c3=c1+nbead
        c4=c2+nbead
        nchoose = i[2]
        sqp = 10*math.sqrt(i[3]) #spring constant proportional to cube_root of prob
        if nchoose == 2:
            intra_both(c1,c2,dscale, kspring=sqp)
            duptcc=duptcc+1
        elif nchoose == 1:
            ambiPairList = [(c1,c2),(c3,c4)]
            minpair(ambiPairList,1,dscale,kspring=sqp)
        else:
            ambiPairList = [(c1,c2),(c1,c4),(c2,c3),(c3,c4)]
            minpair(ambiPairList,nchoose,dscale,kspring=sqp)
            duptcc=duptcc+nchoose-1
        rnum=rnum+1
    for bpair,nchoose in interTarget.items():
        c1 = bpair[0]; c2 = bpair[1]
        c3 = c1+nbead; c4 = c2+nbead
        ambiPairList = [(c1,c2),(c1,c4),(c2,c3),(c3,c4)]
        minpair(ambiPairList,nchoose,dscale,kspring=sqp)
        duptcc=duptcc+nchoose-1
        rnum=rnum+1
        intercount += 1

    return rnum,duptcc,intercount

def condenseChromosome(nbead=ndomain):
    '''collapsing chains into the cenbeads'''
    nchr = len(begin_chromosome)
    for i in range(nchr):
        b1=begin_chromosome[i]
        b2=end_chromosome[i]
        p0A=IMP.core.XYZ(chain.get_particle(subcentro[i])) #subcen
        p0B=IMP.core.XYZ(chain.get_particle(subcentro[i+nchr]))
        coorA = p0A.get_coordinates()
        coorB = p0B.get_coordinates()
        for j in range(b1,b2+1): #from b1 to (and include) b2
            p1A=IMP.core.XYZ(chain.get_particle(j))
            p1B=IMP.core.XYZ(chain.get_particle(j+nbead))
            #p1A=chain.get_particle(j)
            #p1B=chain.get_particle(j+nbead)
            p1A.set_coordinates(coorA) #same as subcen coordinates
            p1B.set_coordinates(coorB) #same as subcen coordinates

def homologuePairing(bstart=0,n=1169,dscale=4):
    '''Restraint all chr pairs'''
    #n is the number of elements in haploid group
    #bstart is the starting bead index, 0-based.
    for j in range(n):
        i = j+bstart
        restrain_one(i,i+n,dscale)
        pair_enforced[(i,i+n)]=['homopair',dscale*(2*Rb[i])]

def addFISH_nearNE(locus, nbead=ndomain):
    '''given locus info, put harmonicLowerBound to NE'''
    rfish = rad_nucleus*locus[0]
    lbfish = IMP.core.HarmonicLowerBound(rfish,1.0)
    bfish = IMP.container.ListSingletonContainer()
    bfish.add_particle(chain.get_particle(getBead(locus[1],ordered_beadInfo)))
    bfish.add_particle(chain.get_particle(nbead+getBead(locus[1],ordered_beadInfo)))
    ssfish = IMP.core.DistanceToSingletonScore(lbfish,center)
    resfish = IMP.container.SingletonsRestraint(ssfish,bfish)
    m.add_restraint(resfish)

def addFISH_farNE(locus, nbead=ndomain):
    '''given locus info, put harmonicUpperBound to NE'''
    rfish = rad_nucleus*locus[0]
    lbfish = IMP.core.HarmonicUpperBound(rfish,1.0)
    bfish = IMP.container.ListSingletonContainer()
    bfish.add_particle(chain.get_particle(getBead(locus[1],ordered_beadInfo)))
    bfish.add_particle(chain.get_particle(nbead+getBead(locus[1],ordered_beadInfo)))
    ssfish = IMP.core.DistanceToSingletonScore(lbfish,center)
    resfish = IMP.container.SingletonsRestraint(ssfish,bfish)
    m.add_restraint(resfish)

def addFISH_actDist(fish_actDist_data, nbead=ndomain):
    '''given FISH_ActDist file, put harmonicLowerBound to NE'''
    fhlist = open(fish_actDist_data,'r')
    b = []
    bnot = []
    for line in fhlist.readlines():
        tmp=string.split(line)
        dcutoff = float(tmp[2]) #rtoNE position threshold
        c1 = int(tmp[0])
        c3 = c1+nbead
        lastdist1 = rad_nucleus-lastcoorRpos(c1)
        lastdist3 = rad_nucleus-lastcoorRpos(c3)
        if(lastdist1 < dcutoff):
            b.append(c1)
        else:
            bnot.append(c1)
        if(lastdist3 < dcutoff):
            b.append(c3)
        else:
            bnot.append(c3)
    lowerbound_fromCenter(b)
    upperbound_fromCenter(bnot)
    return b,bnot

def stayNE_actDist(actDist_file, nbead=ndomain, NEshell=1950.0):
    '''given lad_ActDist file, put harmonicLowerBound to NE'''
    #also returns the list of all the diploid domains
    b = [] #bead IDs range(2*nbead) that will be imposed as LADs (lamin associated domain)
    bnot = [] #beads that will not be imposed to be near NE
    fhlist = open(actDist_file,'r') #must be 0-based indexing!
    for line in fhlist.readlines():
        tmp=string.split(line)
        c1 = int(tmp[0])
        c3 = c1+nbead
        dcutoff = float(tmp[2]) #rtoNE position threshold
        lastdist1 = rad_nucleus-lastcoorRpos(c1)-Rb[c1] #surface measure
        lastdist3 = rad_nucleus-lastcoorRpos(c3)-Rb[c3] #surface measure
        if(lastdist1 < dcutoff):
            b.append(c1)
        else:
            bnot.append(c1)
        if(lastdist3 < dcutoff):
            b.append(c3)
        else:
            bnot.append(c3)
    if(len(b) > 0):
        lowerbound_fromCenter(b, r=NEshell) #close to NE restraint
        print "NE restraints assigned for LADs:",b #2*nbead IDs
    else:
        print 'No NE restraints applied'
    return b,bnot

def stayNE_random(ladprobfile, nbead=ndomain, NEshell=1950.0):
    '''Place random beads to the NE with probability according to ladprobfile'''
    flad =open(ladprobfile, 'r') #this is 1-based indexing
    targetPloci = []
    for line in flad.readlines():
        cols = string.split(line)
        b = int(cols[0])-1 #translate to 0-based id
        p = float(cols[1]) #fraction in population
        pran = random.random()
        bb = [b,b+nbead]
        if pran < p:
            lowerbound_fromCenter(bb, r=NEshell) #close to NE restraint
            print "DamID NE restraints at (random starting config) assigned for LADs:%d with p=%.3f and chance %.3f"%(b,p,pran) #2*nbead IDs

def lowerbound_fromCenter(bset,r=1750.0):
    '''apply lowerBound restraint for beads in bset r away from nuclear center'''
    for b in bset:
        contmp = IMP.container.ListSingletonContainer()
        contmp.add_particle(chain.get_particle(b))
        hlb = IMP.core.HarmonicLowerBound(r-Rb[b],1.0)
        ss = IMP.core.DistanceToSingletonScore(hlb,center)
        res = IMP.container.SingletonsRestraint(ss,contmp)
        m.add_restraint(res)

def upperbound_fromCenter(bset,r=1750.0):
    '''apply upperBound restraint for beads in bset r from nuclear center'''
    bcontainer = IMP.container.ListSingletonContainer()
    hlb = IMP.core.HarmonicUpperBound(r,1.0)
    for b in bset:
        bcontainer.add_particle(chain.get_particle(b))
    ss = IMP.core.DistanceToSingletonScore(hlb,center)
    res = IMP.container.SingletonsRestraint(ss,bcontainer)
    m.add_restraint(res)

def returnPrint(d):
    if d < 1:
        return 0
    else:
        return d

def consec_dist_from_prob(r1,r2,p,xcontact=2):
    '''xcontact is the scaling (r1+r2) where a contact is defined; default=2'''
    try:
        d = (r1+r2)*(1.+(xcontact**3-1)/p)**(1./3)
    except:
        d = 100*(r1+r2) #just a big number
    return d-r1-r2 #will be used as surface-to-surface upperbound


def restrain_d(b1,b2,d,kspring=1):
    '''give d restraint to both copies'''
    upperdist = d
    ds = IMP.core.SphereDistancePairScore(IMP.core.HarmonicUpperBound(upperdist,kspring))
    p1 = chain.get_particle(b1)
    p2 = chain.get_particle(b2)
    apair = IMP.ParticlePair(p1,p2)
    pr = IMP.core.PairRestraint(ds,apair)
    m.add_restraint(pr)

def consec_restraints(nbead=ndomain):
    """ Assign consecutive restraints with distance expected random model"""
    consec_interact={}
    for i in xrange(nbead-1):
        if bead_to_chr[i]==bead_to_chr[i+1]:
            consec_interact[(i,i+1)]=interact[(i,i+1)] #extract the probability input
    consec_num=0
    for k,v in consec_interact.items():
        c1 = k[0]
        c2 = k[1]
        c3 = c1+nbead
        c4 = c2+nbead
        r1=Rb[c1]
        r2=Rb[c2]
        consec_d=consec_dist_from_prob(r1,r2,v)
        restrain_d(c1,c2,consec_d)
        restrain_d(c3,c4,consec_d)
        consec_num+=1
        consecRestraints.append([c1,c2,consec_d])
        pair_enforced[(c1,c2)]=['consec',consec_d] #[type,surface-to-surface]
    return consec_num

def checkNEbeads(beads,rNE=1750):
    unsatisFISH = 0
    s = 0.0 #score for NE restraints
    print 'Number of TADs with NE restraints:',len(beads)
    for b in beads:
        scap = inCapScore(b, rcap=rNE)
        if(scap > 0 ):
            s += scap
            unsatisFISH += 1
            print "Surface of bead %d is under %.1f from center, score= %.4f"%(b,rNE,s)
    print "Total beads violate LAD-NE restraints:",unsatisFISH
    return unsatisFISH


def pEnforced_Distance_wHets(c1,c2,tag, nbead=ndomain):
    '''given bead pair and tag, returns the appropriate distance'''
    #for model with HET on each arm + cen per chr
    tagNoDiplo=['hetnuc','hetcen','nuccen','subcen-het','cen-nuc','rDNA','cen2','cen3','cen4','cenX','homopair']
    d12 = compute_distance(c1,c2)
    if tag in tagNoDiplo:
        d = d12
    else:
        c3 = c1+nbead; c4 = c2+nbead
        d14 = compute_distance(c1,c4)
        d34 = compute_distance(c3,c4)
        d32 = compute_distance(c3,c2)
        alld = sorted([d12,d34,d14,d32])
        if tag == 'min1_inter':
            d = alld[0]
        elif tag == 'min2_inter':
            d = alld[1]
        elif tag in ['min3_inter','min3']:
            d = alld[2]
        elif tag in ['all4_inter','all4']:
            d = alld[3]
        elif tag == 'min1':  #min1 intrachr
            d = min(d12,d34)
        elif tag in ['consec','min2']: #intrachr
            d = max(d12,d34)
    return d

def checkTCC_restraints(dtol=20):
    unsatisfied = defaultdict(list)
    contactScore = 0
    dviolated = []
    for k,v in pair_enforced.items():
        tag = v[0]
        maxdis = v[1]
        c1 = k[0]
        c2 = k[1]
        d = pEnforced_Distance_wHets(c1,c2,tag) - maxdis #surface-to-surface
        #drelative = (d+maxdis)/(maxdis+maxdis) #relative distance
        drelative = d/maxdis #relative distance
        if d > dtol:
            contactScore += 1
            unsatisfied[tag].append([c1,c2])
            dviolated.append(drelative)
            print "Reldist violation pair: %d-%d %s %.4f"%(c1,c2,tag,drelative)
    #dscore = sum(x*x*0.5 for x in dviolated) #use spring equation with k=1
    dscore = sum(x for x in dviolated) #use relative to r1+r2
    print "------- Detected unsatisfied pair enforced %d -------"%(contactScore)
    for k,v in unsatisfied.items():
        print k,len(v)
    print "Contact distance violation score: %.4f"%(dscore)
    return unsatisfied

def checkRestraints(dtol=20, n=totalBead, nucprint=True, beads_NE=[]):
    #excluded volume is not included here
    checkTCC_restraints(dtol=dtol)
    unsatisNuc = 0
    s = 0. #score for nucleus restraint
    for i in range(n):
        scap = outCapScore(i)
        if(scap > eps):
            s += scap
            unsatisNuc += 1
            if(nucprint):
                print "bead %d outside nucleus with score %.4f"%(i, scap)
    print "Total beads outside nucleus: %d with score %.4f"%(unsatisNuc,s)
    if(beads_NE):
        checkNEbeads(beads_NE, rNE=laminR)
    print "==================================="


def getInitialRandom(x0=2000,y0=0,z0=0,magx=500,magy=500,magz=500):
    dx=int((2*random.random()-1)*magx*random.random())
    dy=int((2*random.random()-1)*magy*random.random())
    dz=int((2*random.random()-1)*magz*random.random())
    initialr = IMP.algebra.Vector3D(x0+dx,y0+dy,z0+dz) #initial random 
    return initialr

#===================functions done=====================================

##### Get the last conformation ######################
if(useLastPopulation):
    print "Reading last coordinates",lastpym
    fcoor = open(lastpym,'r')
    lastcoor = []
    for line in fcoor.readlines():
        if line[:7]=="SPHERE,":
            xyz=string.split(line,',')
            x = float(xyz[1])
            y = float(xyz[2])
            z = float(xyz[3])
            lastcoor.append(numpy.array([x,y,z]))
    fcoor.close()

center = IMP.algebra.Vector3D(0,0,0)
randomCap = IMP.algebra.Sphere3D(center,rad_randomCap)
containerNucleolus = IMP.container.ListSingletonContainer()
containerHET = IMP.container.ListSingletonContainer()
containerCEN = IMP.container.ListSingletonContainer()

random.seed()
rdummy=int(random.random()*1000*random.random())
for i in xrange(rdummy):
    ranvec = IMP.algebra.get_random_vector_in(randomCap)
    random.seed()
rdummy=int(random.random()*10000+rdummy)
for i in xrange(rdummy):
    ranvec = IMP.algebra.get_random_vector_in(randomCap)

for i in xrange(totalBead):
    if(useLastPopulation):
        #error will occur if totalBead is not the same in lastcoor
        coor = IMP.algebra.Vector3D(lastcoor[i][0],lastcoor[i][1],lastcoor[i][2])
    else:
        coor = IMP.algebra.get_random_vector_in(randomCap)
    sph = IMP.algebra.Sphere3D(coor, Rb[i])
    p0 = IMP.Particle(m)
    chain.add_particle(p0)
    IMP.atom.Mass.setup_particle(p0,1)
    sp = IMP.core.XYZR.setup_particle(p0,sph)
    sp.set_coordinates_are_optimized(True)
    if i in hetbeads:
        containerHET.add_particle(p0)
    elif i in cenbeads:
        containerCEN.add_particle(p0)
    elif i == nucbead:
        containerNucleolus.add_particle(p0)


#!!!!!!!!!!!!!!!!!!!!!!! RESTRAINTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

consecRestraints=[] #list of consecutive restraints [b1,b2,d]
minpairRestraints=[] #list of minimum pairs restraints [[(b1,b2),(b3,b4),(etc)],to_n]
inbothRestraints=[] #list of minimum pairs restraints [b1,b2]
pair_enforced = {}#register all pair in restrains
homologuePairing(bstart=0, n=nhaploid, dscale=3) #domains
#homologuePairing(bstart=0, n=nhaploid, dscale=4) #domains
homologuePairing(bstart=hetstart, n=nhet, dscale=2) #hetro cluster
homologuePairing(bstart=censtart, n=ncen, dscale=2) #CENs

#connecting diploid domain chains to their het
rsubhet = 0.5
for i in range(2*len(chrname)):
    restrain_one(subcentro[i],hetbeads[i],rsubhet) # regular surface to surface 1
    pair_enforced[(subcentro[i],hetbeads[i])]=['subcen-het',rsubhet*(Rb[subcentro[i]]+Rb[hetbeads[i]])] #optional

#tethering CEN to nucleolus
rcennuc = 0.1
for b in cenbeads:
    restrain_one(b,nucbead,rcennuc) #rscale modify!
    pair_enforced[(b,nucbead)]=['cen-nuc',rcennuc*(rad_nucleolus+Rb[b])] #optional

#rDNA
rnuchet = 0.1
restrain_one(hetbeads[5],nucbead,rnuchet) #Nucleolus touches rDNA X as het!
restrain_one(hetbeads[5+nhet],nucbead,rnuchet) #Nucleolus touches rDNA X as het!

#CEN connecting arms via het
rcenhet = 0.1
restrain_one(hetbeads[0],cenbeads[0],rcenhet) # 2L
restrain_one(hetbeads[0+nhet],cenbeads[0+ncen],rcenhet) # 2L
restrain_one(hetbeads[1],cenbeads[0],rcenhet) # 2R
restrain_one(hetbeads[1+nhet],cenbeads[0+ncen],rcenhet) # 2R
restrain_one(hetbeads[2],cenbeads[1],rcenhet) # 3L
restrain_one(hetbeads[2+nhet],cenbeads[1+ncen],rcenhet) # 3L
restrain_one(hetbeads[3],cenbeads[1],rcenhet) # 3R
restrain_one(hetbeads[3+nhet],cenbeads[1+ncen],rcenhet) # 3R
restrain_one(hetbeads[4],cenbeads[2],rcenhet) # 4
restrain_one(hetbeads[4+nhet],cenbeads[2+ncen],rcenhet) # 4
restrain_one(hetbeads[5],cenbeads[3],rcenhet) # X
restrain_one(hetbeads[5+nhet],cenbeads[3+ncen],rcenhet) # X

#OPTIONAL: pair related to het, cen, nucleolus for scoring check
pair_enforced[(hetbeads[5],nucbead)]=['rDNA',rnuchet*(rad_nucleolus+hetrads[5])]
pair_enforced[(hetbeads[5+nhet],nucbead)]=['rDNA',rnuchet*(rad_nucleolus+hetrads[5])]
pair_enforced[(hetbeads[0],cenbeads[0])]=['cen2',rcenhet*(cenrads[0]+hetrads[0])]
pair_enforced[(hetbeads[0+nhet],cenbeads[0+ncen])]=['cen2',rcenhet*(cenrads[0]+hetrads[0])]
pair_enforced[(hetbeads[1],cenbeads[0])]=['cen2',rcenhet*(cenrads[0]+hetrads[1])]
pair_enforced[(hetbeads[1+nhet],cenbeads[0+ncen])]=['cen2',rcenhet*(cenrads[0]+hetrads[1])]
pair_enforced[(hetbeads[2],cenbeads[1])]=['cen3',rcenhet*(cenrads[1]+hetrads[2])]
pair_enforced[(hetbeads[2+nhet],cenbeads[1+ncen])]=['cen3',rcenhet*(cenrads[1]+hetrads[2])]
pair_enforced[(hetbeads[3],cenbeads[1])]=['cen3',rcenhet*(cenrads[1]+hetrads[3])]
pair_enforced[(hetbeads[3+nhet],cenbeads[1+ncen])]=['cen3',rcenhet*(cenrads[1]+hetrads[3])]
pair_enforced[(hetbeads[4],cenbeads[2])]=['cen4',rcenhet*(cenrads[2]+hetrads[4])]
pair_enforced[(hetbeads[4+nhet],cenbeads[2+ncen])]=['cen4',rcenhet*(cenrads[2]+hetrads[4])]
pair_enforced[(hetbeads[5],cenbeads[3])]=['cenX',rcenhet*(cenrads[3]+hetrads[5])]
pair_enforced[(hetbeads[5+nhet],cenbeads[3+ncen])]=['cenX',rcenhet*(cenrads[3]+hetrads[5])]

# Set up hetClus; all bead is outside hetClus
#lbn = IMP.core.HarmonicLowerBound(rad_hetroClus,1.0)
#dtssn = IMP.core.SphereDistanceToSingletonScore(lbn,hetroClus_center)
#res_hetroClus = IMP.container.SingletonsRestraint(dtssn,chain)
#m.add_restraint(res_hetroClus) #beads in chain do not enter hetroClus

#condenseChromosome()
#---------------------------------------------------------------------------------

t2=time.time()
#print 'score', score
#print "total time spent is %.1fs"%(t2-t1)

# Set up excluded volume
evr = IMP.core.ExcludedVolumeRestraint(chain)
m.add_restraint(evr) #1

# Set up Nucleus cap  #
ubnuc = IMP.core.HarmonicUpperBound(rad_nucleus,1.0)
ssnuc = IMP.core.DistanceToSingletonScore(ubnuc,center) #center-to-center distance
#ssnuc = IMP.core.SphereDistanceToSingletonScore(ubnuc,center) #surface-to-surface distance
rnuc = IMP.container.SingletonsRestraint(ssnuc,chain)
m.add_restraint(rnuc) #2

#smaller cap for nucleolus so it totally will be inside NE
ubnuch = IMP.core.HarmonicUpperBound(rad_nucleus-rad_nucleolus, 1)
ssnuch = IMP.core.DistanceToSingletonScore(ubnuch,center) #center-to-center distance
rnuch = IMP.container.SingletonsRestraint(ssnuch,containerNucleolus)
m.add_restraint(rnuch) #3

#smaller cap for HETs so they totally will be inside NE
for b in hetbeads:
    contmp = IMP.container.ListSingletonContainer()
    contmp.add_particle(chain.get_particle(b))
    ubnuch = IMP.core.HarmonicUpperBound(rad_nucleus-Rb[b], 1)
    ssnuch = IMP.core.DistanceToSingletonScore(ubnuch,center) #center-to-center distance
    m.add_restraint(IMP.container.SingletonsRestraint(ssnuch,contmp))

#Need to cluster centromeres ???
#inDiameterRestraint(cenbeads,diameter=0.5*rad_nucleolus, kspring=1)

#--------!!!!!!!!!!!!!!!!!!!!!!!-------------- TCC ---------!!!!!!!!!!!!!!!!!!!!!!!-----------
minscore = 10*totalBead
added=0
intertcc = 0
doubletcc = 0 #1 is from the top1 higher_freq > fmax 
added = consec_restraints()
nr_consec = added #number of consecutive restraints in haploid
doubletcc = added
print "Total consec restraints",added, doubletcc
added,doubletcc, intertcc = fmax_restraints(xradsum,added,doubletcc)
print "Total fmax restraints %d, duplicates %d, inter %d"%(added-nr_consec,doubletcc-nr_consec, intertcc)

if(useLastPopulation):
    added,doubletcc,intertcc = apply_contact(hlistfile,xradsum,added,doubletcc, intercount=intertcc)
    print '#of TCC restraints %d with duplicates %d and interchr %d'%(added,doubletcc,intertcc)  #number of TCC translated as restraints in a model
#else:
#    stayNE_random(damID_probfile) #DamID is applied at the begining with exact probability, after p100 it will use activation distance

sys.stdout.flush()
#print len(consecRestraints)
#print len(minpairRestraints)
#print len(inbothRestraints)
print "Total pairs enforced including homolog-pair, het/nuc related:", len(pair_enforced)
#print "===========  Optimization Starts =========="
shrinkscore=20
if(useLastPopulation):
    dr = 0.1 #not more than 1!!!
    radStartShrink = int((1+dr)*rad_nucleus)
    radEndShrink = int((1-dr)*rad_nucleus)
    interScale=0.10; shrinkscore=20; highscore=minscore;nucExpand=dr+1.3
    radExpand = int(nucExpand*rad_nucleus)
    incr = int(interScale*rad_nucleus)
    nucrads = range(radEndShrink,radStartShrink,incr)
    nucrads.append(radStartShrink)
    nucrads = sorted(nucrads, reverse=True)
    print "Optimization with decreasing NE %s"%(nucrads)
    expanded=False
    print "\t--- Start shrinking ---"
    for r_nuc in nucrads:
        #print "\t Temporary NE radius %.1f"%(r_nuc)
        m.remove_restraint(rnuc) #will be replaced by temporary
        ubnuc = IMP.core.HarmonicUpperBound(r_nuc,1.0)
        ssnuc = IMP.core.DistanceToSingletonScore(ubnuc,center)
        #ssnuc = IMP.core.SphereDistanceToSingletonScore(ubnuc,center)
        rnuc = IMP.container.SingletonsRestraint(ssnuc,chain)
        m.add_restraint(rnuc) # temporary NE placed
        temp, score = simulated_anneal_scored(2000,300, nc=5, lowscore=shrinkscore)
        print "\t  Shrinking time elapsed %s s, score= %.1f at T=%.1f"%(timespend(t1),score,temp)
        score = cgstep(500)
        print "\t Score optimizing temporary NE %.1f done: %.1f"%(r_nuc,score)
        if score > highscore:
            if expanded == False: #haven't done expansion steps before
                print "\t    --- Start expanding ---"
                m.remove_restraint(rnuc)
                ubexpend = IMP.core.HarmonicUpperBound(radExpand,1.0)
                ssexpend = IMP.core.DistanceToSingletonScore(ubexpend,center)
                rexpend = IMP.container.SingletonsRestraint(ssexpend,chain)
                m.add_restraint(rexpend)
                expanded = True
                temp, score = simulated_anneal_scored(75000,500,nstep=1000, nc=5, lowscore=shrinkscore)
                print "\t     Expansion time elapsed %ss... score= %.1f at T=%.1f"%(timespend(t1),score,temp)
                m.remove_restraint(rexpend)
                m.add_restraint(rnuc)
                print "\t     ...back to shrinking..."
            else:
                print '\t     Score is still high after expansion:',score
                break
    #checkRestraints(dtol=20)
else:
    temp, score = simulated_anneal_scored(50000,10000, nc=2, lowscore=minscore)
    if score > minscore:
        #checkRestraints(dtol=20)
        temp, score = simulated_anneal_scored(500000,10000, nc=5, lowscore=minscore)
    print "p100 T=%.1f K, score now: %.1f"%(temp,score)

# Putting back original radius cap  #
m.remove_restraint(rnuc)
#Adjust cap for surface of domains inside nuclear space
for b in range(ndomain):
    contmp = IMP.container.ListSingletonContainer()
    contmp.add_particle(chain.get_particle(b))
    contmp.add_particle(chain.get_particle(b+nhaploid))
    ubnuc = IMP.core.HarmonicUpperBound(rad_nucleus-Rb[b], 1)
    ssnuc = IMP.core.DistanceToSingletonScore(ubnuc,center) #center-to-center distance
    m.add_restraint(IMP.container.SingletonsRestraint(ssnuc,contmp))


#restraints domains near NE
if(useLastPopulation):
    ladbeads,tmplads = stayNE_actDist(ladActivationDist, NEshell=laminR)

temp, score = simulated_anneal_scored(2000,300, nc=4, nstep=int(2000*scalestep), lowscore=20)
print "Recover nucleus %.1f nm at T=%.1f K, score: %.1f"%(rad_nucleus,temp,score)

print "^_^  relaxing structure  *^_^*"
sys.stdout.flush()

simulated_anneal(10000,5000,nc=5, nstep=int(1000*scalestep))
simulated_anneal(5000,500,nc=5, nstep=int(1000*scalestep))
mdstep(500,500)
score=cgstep(int(1000*scalestep))
print "score T=500K %.1f"%(score)
mdstep(300,500)# from 1000 to 5000
score1=cgstep(int(1000*scalestep)) #from 1000 to 5000
print "score T=300K %.1f"%(score1)
pym("%s.%s%s"%(currentfb,lad_p,pymname))

mdstep(273,1000)
score=cgstep(2000)

#checkRestraints(dtol=20, nucprint=True)
if(useLastPopulation):
    checkRestraints(dtol=1, beads_NE=ladbeads)
else:
    checkRestraints(dtol=1)

print "Final score %.1f" %(score)
if score < score1:
    #pym(currentfb+pymname)
    pym("%s.%s%s"%(currentfb,lad_p,pymname))
    print 'pym file written'

t2=time.time()
print "total time spent is %.1fs"%(t2-t1)
#print the pairs


#--------------------------------------------------------
#m.evaluate(False)
#m.show()
