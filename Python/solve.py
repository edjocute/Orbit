### Following Hernquist-Ostriker 1991# Functions for DMONLY

## Version of force calculation from (16 Oct 2014)
from common import *
import cfunc
import tables
#import gc
#gc.disable()

class Problem():
    def __init__(self, filename='var',nlim=12,llim=6):
        self.G=6.6738e-20 * 1.989e40 #km^3/M_sun/s^2
        self.nlim=nlim
        self.llim=llim
        #temp=np.load(filename)
        with tables.open_file('./'+filename+'.hdf5','r') as f:
            self.K = f.root.Knlm[:]
            self.a = f.root.a[:]
        #postemp=pos.copy()
        #self.a=finda(np.sqrt( (postemp**2).sum(axis=1)),mass)
        #postemp,mass=postosph(postemp, mass)
        #postemp[:,0]/=self.a ###Normalize to a=1
        #if (self.pos[:,1]==1).sum():
    #       raise Exception( "Error: cos(theta)=1 will cause procedure to produce NANs")
                #self.mass=np.delete(mass,self.pos[:,0]==0)
        #self.pos=np.delete(self.pos,self.pos[:,0]==0,axis=0)
                #self.phibar=Phibarlist(self.pos[:,0],self.nlim,self.llim)
        #self.K=knlm(postemp,mass,self.nlim,self.llim)
        #del postemp
        #del mass

    def potential(self,xyz):
        sph= xyz.copy().reshape(len(xyz),3)
        sph= cart2sph(sph/self.a)
        #temp[:,0]=temp[:,0]/self.a ##Normalize input coordinates to a=1##
        return vspotden(self, sph)

    def f(self, xyz, t):
        """xyz is the particle input array [x,y,z,vx,vy,vz] centered around potential minimum"""
        vec= xyz.copy().reshape(len(xyz)/6,6)
        #vec=np.array([xyz])
        diff=np.empty_like(vec)
        ##diff= derivative array [x',y',z',vx',vy',vz']##
        diff[:,:3]=vec[:,3:]
        sph=cart2sph(vec[:,:3]/self.a)
        ##sph= Particle position in spherical coordinates normalized to a=1 ##
        diff[:,3:]=xyzforce(vsforce(self, sph),sph)
        #return diff[0]
        return diff.ravel()
        ##returns 1d array

    def acc(self, pos):
        sph=cart2sph(pos/self.a)
        #aa=vsforce(self,sph)
        aaa=cfunc.vsforce(sph,self.K,self.llim,self.nlim,self.G,self.a)
        return xyzforce(aaa,sph)
        #return xyzforce(vsforce(self,sph),sph)

def vABCD(var,l,m,r,gegen):
    #pos should be converted to r,cos(theta),phi in units of a
    abcd=np.zeros((len(r),4))
    xi=(r-1.)/(r+1.)
    for n in xrange(0,var.nlim):
        #gegentest=special.eval_gegenbauer(n,2.*l+1.5,xi)
        #if not np.allclose(gegentest,gegen[n]):
        #    print 'gegen error!'
        #    print n,l,gegentest[0],gegen[n][0],xi[0]
        rhobartemp=vrhobar_nl(n,l,r,gegen[n])
        Phibartemp=vPhibar_nl(l,r,gegen[n])
        abcd[:,0]+= rhobartemp*var.K[n,l,m,0]
        abcd[:,1]+= rhobartemp*var.K[n,l,m,1]
        abcd[:,2]+= Phibartemp*var.K[n,l,m,0]
        abcd[:,3]+= Phibartemp*var.K[n,l,m,1]
    return abcd[:,0],abcd[:,1],abcd[:,2],abcd[:,3]

def vsforce(var,sph):
    ###Requires positions to be normalized to a=1
    temp=np.zeros((3,len(sph)))
    sintheta=np.sqrt(1-sph[:,1]**2)
    sinthetainv=1./sintheta
    r=sph[:,0]
    rinv=1./r
    rp1inv=1./(1.+r)
    xi=(r-1.)*rp1inv
    sph2=sph[:,2]
    cosmall=[np.cos(m*sph2) for m in range(var.llim)]
    sinmall=[np.sin(m*sph2) for m in range(var.llim)]
    Phi=np.empty((2,var.nlim,len(r)))
    for ll in xrange(0,var.llim):
        gegen1=cfunc.gegen_n(var.nlim,2*ll+1.5,xi)
        #gegen2=cfunc.gegen_nm1(var.nlim,2*ll+2.5,xi)
        Phifac = -sqrt4pi *r**ll * rp1inv**(2*ll+1)
        #Phidifffac1= (8*ll+6)*rp1inv**2 * Phifac
        #Phidifffac2= (ll*rinv-(2*ll+1)*rp1inv) * Phifac
        Phi[0]= Phifac * gegen1
        Phi[1]= Phifac* ((8*ll+6)*rp1inv**2*cfunc.gegen_nm1(var.nlim,2*ll+2.5,xi) + (ll*rinv-(2*ll+1)*rp1inv) *     gegen1)
        for mm in xrange(0,ll+1):
            #### New summation over n using numpy.einsum
            Plmtheta,Plmthetap=cfunc.genlegendrediff(ll,mm,sph[:,1])
            CDEF=np.einsum('ijk,jm->imk',Phi,var.K[:var.nlim,ll,mm])
            #EFCD=np.tensordot(var.K[:var.nlim,ll,mm],Phi,axes=([0],[1])) #tensordot is slower than einsum

            print ll,mm,Plmtheta[0],Plmthetap[0]
            cosmphi,sinmphi=cosmall[mm],sinmall[mm]
            temp[0]-= Plmtheta *(CDEF[1,0]* cosmphi+ CDEF[1,1]*sinmphi )
            temp[1]+= rinv* sintheta * Plmthetap * (CDEF[0,0]*cosmphi +CDEF[0,1]*sinmphi)
            if mm:
                temp[2]-= rinv* mm * Plmtheta*sinthetainv * (CDEF[0,1]*cosmphi -CDEF[0,0]*sinmphi)
    return temp.T*var.G/(var.a**2) ## Normalization to account for G and to scale our results back to a!=1

def vspotden(var,sph):
    """Requires positions to be normalized to a=1"""
    temp=np.zeros((len(sph),2))
    xi=(sph[:,0]-1.)/(1.+sph[:,0])
    for ll in xrange(var.llim):
        gegen=cfunc.gegen_n(var.nlim,2.*ll+1.5,xi)
        for mm in xrange(ll+1):
            A,B,C,D=vABCD(var,ll,mm,sph[:,0],gegen)
            Plm=genlegendre(ll,mm,sph[:,1])
            cosmphi,sinmphi=np.cos(mm*sph[:,2]),np.sin(mm*sph[:,2])
            temp[:,1]+=  Plm* (A*cosmphi + B*sinmphi)  #density
            temp[:,0]+=  Plm* (C*cosmphi + D*sinmphi)  #potential
    return temp*var.G/var.a
