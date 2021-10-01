import os, sys
#sys.path.append('/n/ghernquist/kchua/Orbit/201-code-C-Mar2016/Python/')
sys.path.append('/u/eddc/Codes/HaloShape')
import numpy as np
import readsubfHDF5, readsnapHDF5
import common, utils
import tables
import snapshot
import random as rand
import numexpr as ne
#import axial_ratios as ar
import ellipsoid_numba as ell
import conversions
import argparse

kpc = 3.08568e+16 #km
fields = {1:["Coordinates","Velocities"],4:["Coordinates","Velocities","Masses"],\
          0:["Coordinates","InternalEnergy","Masses","ElectronAbundance"]}
keysel = ["Group_M_Crit200", "Group_R_Crit200", "GroupFirstSub", "GroupPos", "SubhaloVel", "SubhaloPos"]

class Init:

    def __init__(self, catdir, snapnum, savepartids=True,\
            GetRotMat=True, RotateParticles=True, GetCoef=True,
            GetInit=True, GetInitStar=False, relaxedonly=None):

        if catdir.find('/output') < 0:
            self.snapdir = os.path.join(catdir,'output/')
        else:
            self.snapdir = catdir
        print(self.snapdir)
        self.snapnum = snapnum
        
        self.savepartids=savepartids
        self.RotateParticles=RotateParticles
        self.GetRotMat=GetRotMat
        self.GetCoef=GetCoef
        self.GetInit=GetInit
        self.GetInitStar=GetInitStar

        # Read snapshot header for boxsize
        snapstr = str(snapnum).zfill(3)
        self.header = readsnapHDF5.snapshot_header(os.path.join(self.snapdir,'snapdir_'+snapstr,'snap_'+snapstr))
        self.hubble = self.header.hubble
        self.boxsize = self.header.boxsize
        print('\tBoxsize =', self.boxsize)

        self.parttypes = [1]
        if self.header.cooling == 1:
            self.parttypes.append(0)
        if self.header.sfr == 1:
            self.parttypes.append(4)
        print('\tAvailable parttypes =',self.parttypes)

        print('You have chosen the following:')
        if self.RotateParticles: print('Rotate particles using ', self.RotateParticles)
        if self.GetRotMat: print('Save rotation matrix')
        if self.GetCoef: print('Save Knlm coefs')
        if self.GetInit: print('Save DM init file')
        if self.GetInitStar: print('Save stellar init file')

        self.dmpartmass = self.header.massarr[1] #0.00044089652436109581

        if self.GetInitStar and (self.header.sfr != 1):
            self.GetInitStar=False
            print('Warning: no stellar init for DMO runs!')

        # Read SUBFIND catalogue
        self.cat = readsubfHDF5.subfind_catalog(self.snapdir, snapnum, keysel=keysel)
        N = self.cat.filebase.find('/groups_')
        if int(self.cat.filebase[N+8:N+11]) != self.snapnum:
            raise ValueError('Mismatching specification of SnapNum with Subfind catalogue')

    @staticmethod
    def getK(halo, nlim=20, llim=12):
        """
        Get coefficients of HO expansion
        """

        pos = halo["Coordinates"].copy()
        mass = halo["Masses"]

        # Get half-mass radius
        a = common.finda(np.linalg.norm(pos,axis=1),mass)
        
        # Convert to spherical coordinates
        pos = common.postosph(pos)

        # Normalize r to a=1
        pos[:,0] /= a

        # Get coefficients
        K = common.knlm(pos, mass, nlim, llim)

        return K,a

    ## Get initial positions (km) and velocities (km/s):
    @staticmethod
    def getinit(halo, numpart=150, range=(0.05,1), nbins=10):

        out = {}

        edges = np.logspace(np.log10(range[0]), np.log10(range[1]), nbins+1)
        count = np.zeros(nbins, dtype='int')
        pos = halo["Coordinates"]
        vel = halo["Velocities"]
        R200 = halo["R200"]
        dist = np.linalg.norm(pos,axis=1)/R200

        for i in np.arange(nbins):
            bin = ((dist>edges[i]) & (dist<edges[i+1]))
            if bin.sum() < numpart:
                count[i] = bin.sum()
            else:
                count[i] = numpart
        out['count'] = count

        partnum = np.zeros(count.sum(), dtype='int64')
        cumcount = np.zeros(nbins+1, dtype='int64')
        cumcount[1:] = count.cumsum()
        #print(cumcount)

        for i in np.arange(nbins):
            bin = ((dist>edges[i]) & (dist<edges[i+1]))
            partnum[cumcount[i]:cumcount[i+1]] = rand.sample(np.nonzero(bin)[0],count[i])
        temp = partnum

        out['edges'] = edges
        out['Coordinates'] = np.hstack((pos[temp],vel[temp]))#.reshape(15*numpart,6)
        out['ParticleIDs'] = halo["ParticleIDs"][temp]

        return out

    @staticmethod
    def getrotmat(halo):

        R200 = halo["R200"]
        pos = halo["Coordinates"]/R200
        mass = halo["Masses"]
        if mass.std() < 1e-6: mass=None

        rotmat = np.zeros((3,3,3))
        qs = np.zeros((3,2))

        qs[0,0],qs[0,1],n,rotmat[0],_ = ell.ellipsoidfit(pos, R200, 0, 0.2, mass, weighted=True)
        qs[1,0],qs[1,1],n,rotmat[1],_ = ell.ellipsoidfit(pos, R200, 0, 0.5, mass, weighted=True)
        qs[2,0],qs[2,1],n,rotmat[2],_ = ell.ellipsoidfit(pos, R200, 0, 1.0, mass, weighted=True)

        return qs, rotmat


    def readhalo(self, groupid, parttype=1):
        assert parttype in self.parttypes, 'Specified parttype is not available in this run'

        snapdir = self.snapdir
        snapnum = self.snapnum
        cat = self.cat
        centre = cat.GroupPos[groupid]
        subnum = cat.GroupFirstSub[groupid]
        groupvel = cat.SubhaloVel[subnum]
        

	readfields = fields[parttype]
        if self.savepartids is True:
            readfields = readfields + ["ParticleIDs"]
        #print(readfields)	
        halo = snapshot.loadSubhalo(snapdir, snapnum, subnum, parttype, readfields)

        try:
            halo["Coordinates"] = utils.image(halo["Coordinates"]-centre, None, self.boxsize)/self.hubble #in kpc
            halo["Velocities"] = halo["Velocities"] - groupvel #in km/s
        except:
            print('readhalo failed:', groupid, centre)
            return -1
  
        # Set mass for DM particles
        if parttype == 1:
            halo["Masses"] = np.ones(halo["count"]) * self.dmpartmass/self.hubble
        
        # Calculate gastemp for gas elements  
        elif parttype == 0:
            halo["temp"] = conversions.GetTemp(halo["InternalEnergy"], halo["ElectronAbundance"], 5./3.)

        halo["R200"] = cat.Group_R_Crit200[groupid]/self.hubble #in kpc
        halo["groupid"] = groupid

        return halo


    def loadandwriteDM(self, groups, savedir):

        cat = self.cat

        for grp in groups:
            # Read halo
            halo = self.readhalo(grp, 1)

            # Get rotation matrices
            if self.GetRotMat:
                qs, rotmat = self.getrotmat(halo)

            # Rotate particles:
            if self.RotateParticles:
                halo["Coordinates"] = np.dot(halo["Coordinates"], rotmat[2])
                halo["Velocities"] = np.dot(halo["Velocities"], rotmat[2])

            # Get test particles
            if self.GetInit:
                initout = self.getinit(halo)

            # Get Knlm
            if self.GetCoef:
                K,a = self.getK(halo)

    def loadandwriteFP(self, groups, savedir):

        assert 4 in self.parttypes
        cat=self.cat

        for grp in groups:
            halodm = self.loadhalo(grp, 1)
            halostars = self.loadload(grp, 4)
            halogas = self.loadhalo(grp, 0)
            selh=gastemp>=1e5
            selc=gastemp<=1e5
    
            #Get rotation matrices
            halo = {}
            if self.GetRotMat:
                halo["Coordinates"] = np.vstack((halodm["Coordinates"], halostars["Coordinates"], halogas["Coordinates"]))
                halo["Masses"] = np.hstack((halodm["Masses"], halostars["Masses"], halogas["Masses"]))
                halo["R200"] = halodm["R200"]
                halo["groupid"] = groupid
                rotmat = self.getrotmat()
    
            #Rotate particles if needed:
            if self.RotateParticles:
                halodm["Coordinates"] = (np.dot(halodm["Coordinates"],rotmat))
                halodm["Velocities"] = (np.dot(halodm["Velocities"],rotmat))
                halostars["Coordinates"] = (np.dot(halostars["Coordinates"],rotmat))
                halostars["Velocities"] = (np.dot(halostars["Velocities"],rotmat))
                halogas["Coordinates"] = (np.dot(halogas["Coordinates"],rotmat))
   
            #Save test particles:
            if self.GetInit:
                xinit, count, edges, partids = self.getinit(halodm, 100)
    
            if self.GetInitStar:
                xinit, count, edges, partids = self.getinit(halostars, 100, range=[0.05,0.5])
    
            #Save Knlm:
            if self.GetCoef:
                halo["Coordinates"] = np.vstack((halodm["Coordinates"],halogas["Coordinates"][selh]))
                halo["Masses"] = np.hstack((halodm["Masses"],halogas["Masses"][selh]))
                Kh, ah = self.getK(halo)
    
                halo["Coordinates"] = np.vstack((halostars["Coordinates"],halogas["Coordinates"][selc]))
                halo["Masses"] = np.hstack((halostars["Masses"],halogas["Masses"][selc]))
                Kc, ac = self.getK(halo)

if __name__ == '__main__':

    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("snap", type=int, help="snapshot number")
    parser.add_argument("t", help="Sim type: FP, DM or NR")
    parser.add_argument("-c","--chunksize", type=int, default=100, help="Chunksize for parallel calculation")
    parser.add_argument("--minmass", type=float, default=11, help="Min halo mass to be considered")
    parser.add_argument("--maxmass", type=float, default=100, help="Min halo mass to be considered")
    parser.add_argument("--nthreads", type=int, default=4, help="Number of threads for numexpr")
    parser.add_argument("--test", action="store_true", help="Turn on testing mode")
    parser.add_argument("-b","--base", type=str, default='/virgo/simulations/IllustrisTNG/', 
			help="Base directory. Default = /virgo/simulations/IllustrisTNG/")
    args = parser.parse_args()
    
    print('######################################################')
    print('Calculating coefficients!\n')
    fbase = args.base
    fdir = os.path.join(fbase,args.t)
    assert os.path.isdir(fdir), fdir+" does not exist!"
    
    snap = args.snap
    snapstr=str(snap).zfill(3)
    fout = args.t+'_HOcoef_'+snapstr+'.hdf5'
    
    header = readsnapHDF5.snapshot_header(fdir+'/output/snapdir_'+snapstr+'/snap_'+snapstr)
    boxsize = header.boxsize
    hubble = header.hubble
    
    chunksize = args.chunksize
    ell.numba.config.NUMBA_NUM_THREADS = args.nthreads
    ne.set_num_threads(args.nthreads)
    mass=[args.minmass,args.maxmass]
    if args.test:
        mass=[12,12.01]
        chunksize=5
    
    print('Directory: ',fdir)
    print('Snapshot: ', snap)
    print('Boxsize: ', boxsize)
    print('Hubble parameter: ', hubble)
    print('Halo mass range for calculation: 10^'+str(mass),'M_sun')
    #print('Number of cores for sharedmem: ', sharedmem.cpu_count())
    print('Number of threads for numexpr: ', args.nthreads)
    print(' ')


    grp = 1000
    coef = Init(fdir, args.snap)
    print(coef.cat.Group_M_Crit200[grp])

    halo = coef.readhalo(grp, 1)
    print(halo["count"])

    qs, rotmat = coef.getrotmat(halo)
    halo["Coordinates"] = np.dot(halo["Coordinates"], rotmat[2])
    halo["Velocities"] = np.dot(halo["Velocities"], rotmat[2])

    initout = coef.getinit(halo)
    Knlm,a = coef.getK(halo)

    nbodypot, nbodyforce = common.nbodypot(halo['Coordinates'], halo['Masses'], initout['Coordinates'][:,:3])
    
    sph = common.cart2sph(initout["Coordinates"][:,:3])
    sph[:,0] /= a

    pot1 = common.reconstruct_potden(Knlm, a, sph, 20, 12)
    force, pot = common.getforce(Knlm, a, sph, 20, 12)
    force = common.xyzforce(force, sph)

    print(np.allclose(pot1[:,0], pot))
    print('Potential recontruction', np.abs((pot-nbodypot)/nbodypot).mean() * 100)
    print('Force reconstruction', np.linalg.norm(force-nbodyforce)/np.linalg.norm(nbodyforce))

    np.savez('test_{}.npz'.format(grp), 
             grpid=grp, subid=coef.cat.GroupFirstSub[grp], Knlm=Knlm, a=a, 
             rotmat=rotmat, qs=qs, init=initout, nbodypot=nbodypot, nbodyforce=nbodyforce,
             pot=pot, force=force)

    #SSE = np.zeros((8,11))
    #ferr = np.zeros_like(SSE)
    #for i,nlim in enumerate(range(4,20,2)):
    #    for j,llim in enumerate(range(1,12,1)):
#	    pot = common.reconstruct_potden(Knlm, a, sph, nlim, llim)
#            SSE[i,j] = ((pot[:,0]-nbodypot)**2).sum()
#            ferr[i,j] = np.abs((pot[:,0]-nbodypot)/nbodypot).mean() * 100


