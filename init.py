import os, sys
sys.path.append('/n/ghernquist/kchua/Orbit/01-code')
sys.path.append('/n/ghernquist/kchua/Shape/02-Runs/code')
import readsubfHDF5
from common import *
import loadhalo
import random as rand
import axial_ratios as ar
import conversions

## Writes coefficient of HO expansion into file
def getK(pos,mass,rvir,filename='./var',nlim=16,llim=10):
    postemp=pos.copy()
    a=finda( np.sqrt((postemp**2).sum(axis=1)),mass)
    postemp,mass=postosph(postemp, mass)
    postemp[:,0]/=a ###Normalize to a=1
    #if (self.pos[:,1]==1).sum():
    #raise Exception( "Error: cos(theta)=1 will cause procedure to produce NANs")
        #self.mass=np.delete(mass,self.pos[:,0]==0)
    #self.pos=np.delete(self.pos,self.pos[:,0]==0,axis=0)
        #self.phibar=Phibarlist(self.pos[:,0],self.nlim,self.llim)
    K=knlm(postemp,mass,nlim,llim)
    np.savez(filename,knlm=K,a=a,rvir=rvir)

## Writes initial positions (km) and velocities (km/s) into file:
def getinit(pos,vel,rvir,numpart,savedir='./',nbins=15):
    edges=np.logspace(np.log10(0.02),np.log10(0.5),nbins+1)
    partnum=np.zeros((nbins,numpart),dtype='int')
    dist= np.sqrt((pos**2).sum(axis=1))/rvir
    for i in np.arange(nbins):
        partnum[i]=rand.sample(np.nonzero( (dist>edges[i]) & (dist<edges[i+1]))[0],numpart)
    #temp=partnum[:5].ravel()
    #xa=np.hstack((pos[temp],vel[temp]))
    #temp=partnum[5:10].ravel()
    #xb=np.hstack((pos[temp],vel[temp]))
    #temp=partnum[10:].ravel()
    #xc=np.hstack((pos[temp],vel[temp]))
    temp=partnum.ravel()
    x=np.hstack((pos[temp],vel[temp])).reshape(15,numpart,6)
    #np.savez(savedir+'init',x=(xa,xb,xc),pos=pos[partnum.ravel()],vel=vel[partnum.ravel()])
    np.savez(savedir+'init',x=x,pos=pos[partnum.ravel()],vel=vel[partnum.ravel()])

def getrotmat(cat,grp,dir,savedir='./'):
    print dir,grp
    q,s,n,rotmat=ar.axial(cat,dir,135,grp,15,rmin=2.24e-2,rmax=0.452,binwidth=0.15)
    np.savez(savedir+'qs',q=q,s=s,n=n,rotmat=rotmat,mvir=cat.Group_M_Crit200[grp],rvir=cat.Group_R_Crit200[grp])

kpc=3.08568e+16 #km
h=0.704 #Hubble parameter
fields=["Coordinates","Velocities"]

def loadandwritedm(cat,groupid,dir,savedir='./'):
    grouppos=cat.GroupPos[groupid]
    groupvel=cat.SubhaloVel[cat.GroupFirstSub[groupid]]
    rvir=cat.Group_R_Crit200[groupid]/h*kpc

    part=loadhalo.loadSnapSubset(dir,135,cat.GroupFirstSub[groupid],0,1,fields)
    pos=(image(grouppos,part.get("Coordinates"),75000)-grouppos)/h*kpc
    vel=part.get("Velocities")-groupvel
    #mass=np.ones(len(pos))*0.00052946428432085776/h
    mass=np.ones(len(pos))*0.00044089652436109581/h

    getK(pos,mass,cat.Group_R_Crit200[groupid],savedir+'var.npz')
    #getinit(pos,vel,rvir,100,savedir)

def loadandwritehot(cat,groupid,dir,savedir='./'):
    grouppos=cat.GroupPos[groupid]
    groupvel=cat.SubhaloVel[cat.GroupFirstSub[groupid]]
    rvir=cat.Group_R_Crit200[groupid]/h*kpc
    subid=cat.GroupFirstSub[groupid]

    partdm=loadhalo.loadSnapSubset(dir,135,cat.GroupFirstSub[groupid],0,1,fields)
    posdm=(image(grouppos,partdm.get("Coordinates"),75000)-grouppos)/h*kpc
    veldm=partdm.get("Velocities")-groupvel
    massdm=np.ones(len(posdm))*0.00044089652436109581/h

    partgas=loadhalo.loadSnapSubset(dir,135,subid,0,0,["Coordinates","InternalEnergy","Masses","ElectronAbundance"])
    gastemp=conversions.GetTemp(partgas.get("InternalEnergy"),partgas.get("ElectronAbundance"),5./3.)
    posh=image(np.zeros(3),partgas.get("Coordinates")[gastemp >= 1e5]-grouppos,75000)/h*kpc
    massh=partgas.get("Masses")[gastemp >= 1e5]/h

    getinit(posdm,veldm,rvir,100,savedir)
    getK(np.vstack((posdm,posh)),np.hstack((massdm,massh)),cat.Group_R_Crit200[groupid],savedir+'varh.npz')

def loadandwritecold(cat,groupid,dir,savedir='./'):
    grouppos=cat.GroupPos[groupid]
    groupvel=cat.SubhaloVel[cat.GroupFirstSub[groupid]]
    rvir=cat.Group_R_Crit200[groupid]/h*kpc
    subid=cat.GroupFirstSub[groupid]

    partstar=loadhalo.loadSnapSubset(dir,135,subid,0,4,["Coordinates","Velocities","Masses"])
    poss=image(np.zeros(3),partstar.get("Coordinates")-grouppos,75000)/h*kpc
    vels=partstar.get("Velocities")-groupvel
    masss=partstar.get("Masses")

    partgas=loadhalo.loadSnapSubset(dir,135,subid,0,0,["Coordinates","InternalEnergy","Masses","ElectronAbundance"])
    gastemp=conversions.GetTemp(partgas.get("InternalEnergy"),partgas.get("ElectronAbundance"),5./3.)
    posc=image(np.zeros(3),partgas.get("Coordinates")[gastemp <= 1e5]-grouppos,75000)/h*kpc
    massc=partgas.get("Masses")[gastemp <= 1e5]/h

    getK(np.vstack((poss,posc)),np.hstack((masss,massc)),cat.Group_R_Crit200[groupid],savedir+'varc.npz')

catdir='/n/hernquistfs1/Illustris/Runs/L75n1820FP/'
cat=readsubfHDF5.subfind_catalog(catdir+'/output/',135,keysel=["Group_M_Crit200","Group_R_Crit200",  "GroupFirstSub",  "GroupPos", "SubhaloVel","SubhaloPos"])
#all=np.load('/n/ghernquist/kchua/Orbit/10-match/matchedgrps.npz')
#List=np.unique(np.sort(all['DM'][all['m200fp'] > 50]))
List=[237,305,317,362,372,377,378,408,416,417,418,419,438,445,448,457,458,472,481,484,488,533,536,548,551,569,579,590,633,646,650,656,661,664,670,699,701,726,738,742,748,753,774,791,807,820,854,858,879,956]

for i in range(len(List)):
    savedir = "GROUP_"+str(List[i]).zfill(3)+"/"
    if not os.path.isdir(savedir):os.makedirs(savedir)
    loadandwritedm(cat,List[i],catdir,savedir)
    loadandwritehot(cat,List[i],catdir,savedir)
    loadandwritecold(cat,List[i],catdir,savedir)
    getrotmat(cat,List[i],catdir,savedir)
