import os, sys
#sys.path.append('/n/ghernquist/kchua/Orbit/01-code/')
sys.path.append('/n/ghernquist/kchua/Orbit/201-code-C-Mar2016/Python/')
sys.path.append('/n/ghernquist/kchua/Shape/02-Runs/code')
#sys.path.append('/n/home04/kchua/PYTHON/Illustris/illustris_python')
import numpy as np
import readsubfHDF5
import common
import tables
import snapshot
import random as rand
import axial_ratios as ar
import conversions

kpc=3.08568e+16 #km
h=0.704 #Hubble parameter
fields={1:["Coordinates","Velocities"],4:["Coordinates","Velocities","Masses"],\
        0:["Coordinates","InternalEnergy","Masses","ElectronAbundance"]}
keysel=["Group_M_Crit200","Group_R_Crit200","GroupFirstSub","GroupPos","SubhaloVel",           "SubhaloPos"]

class Init:
    def __init__(self,catdir,savedir='./',snapnum=135,savepartids=False,RotateParticles=False):
        self.catdir=catdir
        self.snapnum=snapnum
        self.savedir=savedir

        self.cat=readsubfHDF5.subfind_catalog(catdir,snapnum,keysel=keysel)
        #self.rvir=self.cat.Group_R_Crit200[groupid]/h*kpc
        self.savepartids=savepartids
        self.RotateParticles=RotateParticles

        self.fields={1:["Coordinates","Velocities"],4:["Coordinates","Velocities","Masses"],\
        0:["Coordinates","InternalEnergy","Masses","ElectronAbundance"]}

        if catdir.find('1820DM') >0:
            self.whichdm='1820DM'
        elif catdir.filebase.find('1820FP') >0:
            self.whichdm='1820FP'
        else:
            print "Cant find DM or FP!"
            sys.exit(0)

        print self.catdir,self.snapnum

        #self.catdir=cat.filebase.split(self.whichdm)[0]+self.whichdm+'/output/'
        N=self.cat.filebase.find('/groups_')
        if int(self.cat.filebase[N+8:N+11]) != self.snapnum:
            print "Mismatching specification of SnapNum with Subfind catalogue"
            sys.exit(0)

    ## Writes coefficient of HO expansion into file
    def getK(self,partdata,filename,nlim=16,llim=10):
        pos=partdata["Coordinates"]
        mass=partdata["Masses"]
        postemp=pos.copy()
        a=common.finda( np.sqrt((postemp**2).sum(axis=1)),mass)
        postemp,mass=common.postosph(postemp, mass)
        postemp[:,0]/=a ###Normalize to a=1
        K=common.knlm(postemp,mass,nlim,llim)
        with tables.open_file(self.savedir+filename,'w') as f:
            f.create_carray("/","a",tables.Float64Atom(),(1,))
            f.create_carray("/","Rvir",tables.Float64Atom(),(1,))
            f.create_carray("/","Knlm",tables.Float64Col(),(nlim,llim,llim,2))
            f.root.a[:]=a
            f.root.Rvir[:]=self.cat.Group_R_Crit200[partdata["groupid"]]
            f.root.Knlm[:]=K
        return K

    ## Writes initial positions (km) and velocities (km/s) into file:
    def getinit(self,partdata,numpart,filename='init.hdf5',nbins=15):
        edges=np.logspace(np.log10(0.02),np.log10(0.5),nbins+1)
        partnum=np.zeros((nbins,numpart),dtype='int')
        pos=partdata["Coordinates"]
        vel=partdata["Velocities"]
        dist= np.sqrt((pos**2).sum(axis=1))/partdata["rvir"]
        for i in np.arange(nbins):
            partnum[i]=rand.sample(np.nonzero( (dist>edges[i]) & (dist<edges[i+1]))[0],numpart)
        temp=partnum.ravel()
        x=np.hstack((pos[temp],vel[temp]))#.reshape(15*numpart,6)
        #np.savez(savedir+'init',x=x,pos=pos[partnum.ravel()],vel=vel[partnum.ravel()])
        with tables.open_file(self.savedir+'init.hdf5','w') as f:
            f.create_carray("/","x",tables.Float64Col(),(15*numpart,6))
            f.root.x[:]=x
            if self.savepartids == True:
                f.create_carray("/","IDs",tables.UInt64Col(),(15*numpart,))
                f.root.IDs[:]=partdata["ParticleIDs"][temp]

    def getrotmat(self,groupid,saving=True):
        print self.catdir,groupid
        q,s,n,rotmat=ar.axial(self.cat,self.catdir,self.snapnum,groupid,2,rmin=0.02,rmax=0.5)
        if saving==True:
            with tables.open_file(self.savedir+'shape.hdf5','w') as f:
                f.create_carray("/","qs",tables.Float64Col(),(len(q),2))
                f.create_carray("/","rotmat",tables.Float64Col(),(len(q),3,3))
                f.create_carray("/","Rvir",tables.Float64Atom(),(1,))
                f.create_carray("/","Mvir",tables.Float64Atom(),(1,))
                f.root.qs[:,0]=q
                f.root.qs[:,1]=s
                f.root.rotmat[:]=rotmat
                f.root.Rvir[:]=self.cat.Group_R_Crit200[groupid]
                f.root.Mvir[:]=self.cat.Group_M_Crit200[groupid]
        return rotmat

    #loadHalo(basePath,snapNum,id,partType,fields=None)
    def loadandwriteDM(self,groupid):
        cat=self.cat
        grouppos=cat.GroupPos[groupid]
        groupvel=cat.SubhaloVel[cat.GroupFirstSub[groupid]]
        rotmat=self.getrotmat(groupid,saving=True)[1]

        if self.savepartids==True:
            self.fields[1]=self.fields[1]+["ParticleIDs"]
        part=snapshot.loadSubhalo(self.catdir,self.snapnum,cat.GroupFirstSub[groupid],1,self.fields[1])
        part["Coordinates"]=(common.image(grouppos,part["Coordinates"],75000)-grouppos)/h*kpc
        part["Velocities"]=part["Velocities"]-groupvel
        #pos=(image(grouppos,part.get("Coordinates"),75000)-grouppos)/h*kpc
        #vel=part.get("Velocities")-groupvel
        if cat.filebase.find('DM') >0:
            part["Masses"]=np.ones(part["count"])*0.00052946428432085776/h
        elif cat.filebase.find('FP') >0:
            part["Masses"]=np.ones(part["count"])*0.00044089652436109581/h
        else:
            print "Cant find DM or FP!"
        part["rvir"]=self.cat.Group_R_Crit200[groupid]/h*kpc
        part["groupid"]=groupid

        #Rotate particles if needed:
        if self.RotateParticles==True:
            part["Coordinates"]=(np.dot(part["Coordinates"],rotmat))#.reshape(nn,N,3)
            part["Velocities"]=(np.dot(part["Velocities"],rotmat))#.reshape(nn,N,3)

        #Get test particles
        self.getinit(part,100,'init.hdf5')

        #Get Knlm
        return self.getK(part,self.savedir+'var.hdf5')

    def loadandwriteFP(self,component,groupid):
        cat=self.cat
        grouppos=cat.GroupPos[groupid]
        groupvel=cat.SubhaloVel[cat.GroupFirstSub[groupid]]
        subid=cat.GroupFirstSub[groupid]
        rotmat=self.getrotmat(groupid,saving=True)[1]

        if component=='hot':
            savefile=self.savedir+'varh.hdf5'
            parttype=1
        elif component=='cold':
            savefile=self.savedir+'varc.hdf5'
            parttype=4
        part1=snapshot.loadSubhalo(self.catdir,self.snapnum,subid,parttype,self.fields[parttype])
        part1["Coordinates"]=(common.image(grouppos,part1["Coordinates"],75000)-grouppos)/h*kpc
        part1["Velocities"]=part1["Velocities"]-groupvel

        partgas=snapshot.loadSubhalo(self.catdir,self.snapnum,subid,0,self.field[0])
        gastemp=conversions.GetTemp(partgas["InternalEnergy"],partgas["ElectronAbundance"],5./3.)
        if component=='hot':
            sel=gastemp>=1e5
        elif component=='cold':
            self=gastemp<=1e5
        partgas["Coordinates"]=(common.image(grouppos,partgas["Coordinates"][sel],75000)-grouppos)/h*kpc
        partgas["Masses"]=partgas["Masses"][sel]/h

        #Rotate particles if needed:
        if self.RotateParticles==True:
            part1["Coordinates"]=(np.dot(part1["Coordinates"],rotmat))#.reshape(nn,N,3)
            part1["Velocities"]=(np.dot(part1["Velocities"],rotmat))#.reshape(nn,N,3)
            partgas["Coordinates"]=(np.dot(partgas["Coordinates"],rotmat))#.reshape(nn,N,3)
            partgas["Velocities"]=(np.dot(partgas["Velocities"],rotmat))#.reshape(nn,N,3)

        #Save test particles:
        if component=='hot':
            part1["Masses"]=np.ones(part1["count"])*0.00044089652436109581/h
            self.getinit(part1,100,'initdm.hdf5')
        elif component=='cold':
            self.getinit(part1,100,'initstar.hdf5')

        #Save Knlm:
        part={}
        part["Coordinates"]=np.vstack((part1["Coordinates"],partgas["Coordinates"]))
        part["Masses"]=np.hstack((part1["Masses"],partgas["Masses"]))
        part["rvir"]=self.cat.Group_R_Crit200[groupid]/h*kpc
        part["groupid"]=groupid
        self.getK(part,savefile)

#catdir='/n/hernquistfs2/Illustris/Runs/L75n1820DM/output/'
#List=[879]

#for i in range(len(List)):
    #savedir = "GROUP_"+str(List[i]).zfill(3)+"/"
    #savedir = "./"
    #if not os.path.isdir(savedir):os.makedirs(savedir)
    #init=Init(cat,catdir,savedir,snapnum=103)
    #init.loadandwriteDM(List[i])
    #init.loadandwriteFP('hot',List[i])
    #init.loadandwriteFP('cold',List[i])
    #getrotmat(cat,List[i],catdir,savedir,saving=True)
