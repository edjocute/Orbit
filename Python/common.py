### Following Hernquist-Ostriker 1991
import numpy as np
from scipy import special
import numexpr as ne
from numba import jit
import numba

#ne.set_num_threads(1)

##################
### Conversion ###
##################

# Constants
sqrt4pi = np.sqrt(4*np.pi)
G = 6.6738e-20 * 1.989e40 #km^3/M_sun/s^2

def sph2cart(sph):
    """
    Convert from spherical to cartesian coordinates
    """ 
    ptsnew = np.zeros(sph.shape)
    sintheta=np.sqrt(1-sph[:,1]**2)
    ptsnew[:,0] = sph[:,0]*sintheta*np.cos(sph[:,2])
    ptsnew[:,1] = sph[:,0]*sintheta*np.sin(sph[:,2])
    ptsnew[:,2] = sph[:,0]*sph[:,1]

    return ptsnew

def cart2sph(xyz):
    """
    Convert from cartesian to spherical coordinates
    """ 
    if len(xyz.shape)==1:
        xy = xyz[0]**2 + xyz[1]**2
        r = np.sqrt(xy + xyz[2]**2)
        if r==0:
            return np.array([0.,1e-5,1e-5])
        else:
            return np.array([r,xyz[2]/r, np.arctan2(xyz[1], xyz[0]) ])
    else:
        ptsnew = np.zeros(xyz.shape)
        xy = xyz[:,0]**2 + xyz[:,1]**2
        ptsnew[:,0] = np.sqrt(xy + xyz[:,2]**2) #r
    #ptsnew[:,1] = np.arctan2(np.sqrt(xy), xyz[:,2]) #theta defined from Z-axis down
        ptsnew[:,1] = xyz[:,2]/ptsnew[:,0] # returns cos(theta) instead
        ptsnew[:,2] = np.arctan2(xyz[:,1], xyz[:,0]) #phi
        return ptsnew

def finda(dist,mass):
    #sph = particle spherical coordinates (r, cos(theta), phi)
    #mass = array of particle masses
    #dist=sph[:,0]
    argsort = dist.argsort()
    halfmass = dist[argsort][(mass[argsort].cumsum() >= mass.sum()/2.)][0:2]
    halfmass = halfmass.mean()

    return halfmass/(1.+np.sqrt(2))

##Returns positions normalized to a=1 in spherical coordinates
def postosph(pos):
    #pos=pos-pos[np.where(pot==pot.min())[0]] #Shift potential minimum to center
    pos[np.where( (pos**2).sum(axis=1)==0    )[0]]+=[1e-15,0,0] # Give potential minimum nonzero angle
    #mask=np.where((pos**2).sum(axis=1)==0)[0]
    #mass=np.delete(mass,mask)
    #pos=np.delete(pos,mask,axis=0)
    pos=cart2sph(pos)

    return pos

def xyzforce(fvec,sph):
    fxyz=np.empty_like(fvec)
    costheta=sph[:,1]
    sintheta=np.sqrt(1-costheta**2)
    cosphi=np.cos(sph[:,2])
    sinphi=np.sin(sph[:,2])
    fxyz[:,0] = sintheta*cosphi*fvec[:,0] + costheta*cosphi*fvec[:,1] - sinphi*fvec[:,2]
    fxyz[:,1] = sintheta*sinphi*fvec[:,0] + costheta*sinphi*fvec[:,1] + cosphi*fvec[:,2]
    fxyz[:,2] = costheta*fvec[:,0] - sintheta*fvec[:,1]

    return fxyz

#########################
### Special Functions ###
#########################

def factorial2(m):
    # computes (m)!! only for odd m
    res = 1
    for i in xrange(1,m+1,2):
        res *= i
    return res

def factorial(m):
    # computes (m)!
    res = 1
    for i in xrange(1,m+1):
        res *= i
    return res

def genlegendre(l,m,x):
    ### Returns the associated Legendre polynomial P_lm(x)
    assert 0<=m<=l and all(abs(x)<=1.)

    pmm= (-1.)**m * factorial2(2*m-1) * (1.-x**2)**(m/2.)
    if l==m:
        return pmm
    pmm1= x*(2.*m+1)*pmm
    if l == (m+1):
        return pmm1

    for ll in xrange(m+2,l+1):
        temp=(x* (2.*ll-1.) * pmm1 - (ll+m-1.)*pmm)/(ll-m)
        pmm=pmm1
        pmm1=temp
    return temp

def genlegendrelist(l,m,x):
    Plist=np.zeros(llim,llim,len(x))

    for ll in xrange(llim):
        for mm in xrange(ll+1):
            Plist[ll,mm]=genlegendre(ll,mm,x)
    return Plist

#@profile
def genlegendrediff(l,m,x):
    """ 
    Returns:
    (1) associated Legendre polynomial P_lm(x)
    (2) First derivative of the associated Legendre polynomial dP_lm/dx(x)
    """
    #assert 0<=m<=l and abs(x).max() <= 1

    xx = (1.-x**2)
    fac = factorial2(2*m-1)* xx**(0.5*m-1)
    pmm =  (-1)**m * fac * xx
    #pmm= (-1)**m * fac * (1.-x**2)**(m/2.)
    #ne.evaluate("(-1)**m * fac * (1.-x**2)**(m/2.)")
    if m==0:
        pmmd=np.zeros_like(x)
    else:
		#pmmd=(-1)**(m+1) * factorial2(2*m-1) * m*x *(1.-x**2)**(m/2.-1)
        #pmmd=(-1)**(m+1) * fac * m*x *(1.-x**2)**(m/2.-1)a
        pmmd=(-1)**(m+1) *m*x * fac
    if l==m:
        return pmm,pmmd

    pmm1 = x*(2*m+1)*pmm
    pmm1d = (2*m+1)*(pmm+x*pmmd)
    if l==(m+1):
        return pmm1,pmm1d

    # Recursion
    for ll in xrange(m+2,l+1):
        tempd = ((2*ll-1)*pmm1 - (ll+m-1)*pmmd + x*(2*ll-1)*pmm1d)/(ll-m)
        pmmd = pmm1d
        pmm1d = tempd
        temp = (x*(2*ll-1)*pmm1 - (ll+m-1)*pmm)/(ll-m)
        pmm = pmm1
        pmm1 = temp

    return temp,tempd

def gegenbauerlist2(xi,nlim,llim): #Tabulates gegenbauer values for all halo particles (memory intensive)
    Clist=np.zeros((nlim,llim,len(xi))    )

    alpha=np.array([2*i+1.5*np.ones(len(xi)) for i in xrange(llim)])
    Clist[0]= 1.
    if nlim==1: return Clist

    Clist[1]= 2*alpha*xi
    if nlim==2: return Clist

    for n in xrange(2,nlim):
        Clist[n]= (2*(n+alpha-1)*xi*Clist[n-1] - (n+2*alpha-2)*Clist[n-2])/n
    return Clist

# Tabulates ultraspherical harmonics C^a_n for given a (or l)
#@profile
def gegen_n(nlim,a,x):
    assert nlim > 0 and a > 0
    arrayall=np.zeros((nlim,len(x)))
    arrayall[0] = 1.
    arrayall[1] = 2.*a*x
    for nn in xrange(2,nlim):
        arrayall[nn] = (2*(nn+a-1)*x*arrayall[nn-1] - (nn+2*a-2)*arrayall[nn-2])/nn
    return arrayall

# Tabulates ultraspherical harmonics C^a_(n-1) for given a (or l)
#@profile
def gegen_nm1(nlim,a,x):
    assert nlim > 0 and a > 0
    arrayall=np.zeros((nlim,len(x)))
    arrayall[1] = 1.
    arrayall[2] = 2.*a*x
    for nn in xrange(3,nlim):
        arrayall[nn] = (2*(nn+a-2)*x*arrayall[nn-1] - (nn+2*a-3)*arrayall[nn-2])/(nn-1)
    return arrayall

##########################################################################################################
###                 Basis Functions and Expansion coefficient calculation                              ###
##########################################################################################################

def K_nl(n,l):
    return 0.5 * n *(n+4.*l+3) + (l+1.)*(2.*l+1)

def I_nl(n,l):
    return -K_nl(n,l) *4*np.pi/2**(8.*l+6) * factorial(n+4*l+2)/(factorial(n)* (n+2.*l+1.5)* special.gamma(2.*l + 1.5)**2)

def Abar_nl(n,l):
    return 1./I_nl(n,l)

def N_lm(l,m):
    temp=(2.*l+1.)/(4*np.pi) * np.double(factorial(l-m))/factorial(l+m)
    if m==0:
        return temp
    elif m!=0:
        return temp*2.

def vPhibar_nl(l,r,gegen):

    #return -r**l/(1+r)**(2.*l+1)  * gegen * sqrt4pi
    return -np.power(r,l)/np.power(1+r, 2*l+1)  * gegen * sqrt4pi

def vrhobar_nl(n,l,r,gegen):

    #return -K_nl(n,l)/(2*np.pi)*r**(l-1.)/(1+r)**(2.*l+3)  * gegen * sqrt4pi
    return -K_nl(n,l)/(2*np.pi) * np.power(r,l-1)/np.power(1+r, 2*l+3) * gegen * sqrt4pi

def Phibarlist(r,nlim,llim): #Tabulates Phibar basis functions for all halo particles (Decaprecated because memory intensive)
    gegen=gegenbauerlist2((r-1.)/(r+1.),nlim,llim)
    temp=np.zeros_like(gegen)
    for n in xrange(nlim):
        for l in xrange(llim):
            temp[n,l]= -r**l/(1.+r)**(2*l+1.) * gegen[n,l] *sqrt4pi
    return temp

#@profile
def calc_knlm(pos,mass,nlim,llim):

    """
    Calculate Hernquist-Ostriker coefficients
    """ 
    
    # Create empty arrays
    K = np.zeros((nlim,llim,llim,2))
    r = pos[:,0]
    phi = pos[:,2]

    # New radial coordinate
    xi = (r-1.)/(r+1.)

    for ll in  xrange(0,llim):
        phibarfac = - ne.evaluate("sqrt4pi*(r**ll)/((1.+r)**(2*ll+1.))")
        #phibarfac=-sqrt4pi*(pos[:,0]**ll)/((1.+pos[:,0])**(2*ll+1.))
        gegen = gegen_n(nlim,2*ll+1.5,xi)

        for mm in xrange(0,ll+1):
            Plm=genlegendre(ll,mm,pos[:,1])
            cosm = ne.evaluate("cos(mm*phi)")
            sinm = ne.evaluate("sin(mm*phi)")

            for nn in xrange(nlim):
                #base= mass*  Plm * phibarfac*special.eval_gegenbauer(nn,2*ll+1.5,xi)
                base = mass*Plm*phibarfac*gegen[nn]
                abartemp = Abar_nl(nn,ll)*N_lm(ll,mm)
                K[nn,ll,mm,0] = abartemp*(base * cosm).sum()
                K[nn,ll,mm,1] = abartemp*(base * sinm).sum()
    return K


###########################################################################################################


def vABCD(Knlm, l, m, r, gegen, nlim=None):
    """ 
    Calulcate quantities A,B,C,D
    
    positions should be in (r,cos(theta),phi), with r normalized to a=1
    """
    if nlim is None:
        nlim = Knlm.shape[0]

    ABCD = np.zeros((4, len(r)))
    #xi = (r-1.)/(r+1.)

    for n in xrange(0, nlim):
        rhobartemp = vrhobar_nl(n,l,r,gegen[n])
        Phibartemp = vPhibar_nl(l,r,gegen[n])

        ABCD[0] += rhobartemp * Knlm[n,l,m,0]
        ABCD[1] += rhobartemp * Knlm[n,l,m,1]
        ABCD[2] += Phibartemp * Knlm[n,l,m,0]
        ABCD[3] += Phibartemp * Knlm[n,l,m,1]

    return ABCD


def reconstruct_potden(Knlm, a, sph, nlim, llim):

    """
    Reconstruct potential from HO coefficients

    - Requires positions (sph) to be normalized to a=1
    - knlm: Coefficients
    - sph: spherical coordinates for locations
    """

    temp = np.zeros((len(sph),2))
    xi = (sph[:,0]-1.)/(1.+sph[:,0])

    for ll in xrange(llim):
        gegen = gegen_n(nlim,2.*ll+1.5,xi)

        for mm in xrange(ll+1):
            A,B,C,D = vABCD(Knlm,ll,mm,sph[:,0],gegen, nlim=nlim)
            Plm = genlegendre(ll,mm,sph[:,1])
            cosmphi,sinmphi = np.cos(mm*sph[:,2]),np.sin(mm*sph[:,2])
            temp[:,1] +=  Plm* (A*cosmphi + B*sinmphi)  #density
            temp[:,0] +=  Plm* (C*cosmphi + D*sinmphi)  #potential

    return temp/a #*G

def getforce(Knlm, a, sph, nlim, llim):
    """
    Get SCF forces from HO coefficients

    - Requires positions (sph) to be normalized to a=1
    - knlm: Coefficients
    - sph: spherical coordinates for locations
    """

    force = np.zeros((3,len(sph)))
    pot = np.zeros(len(sph))

    sintheta = np.sqrt(1-sph[:,1]**2)
    sinthetainv = 1./sintheta
    r = sph[:,0]
    rinv = 1./r
    rp1inv = 1./(1.+r)
    xi = (r-1.)*rp1inv
    sph2 = sph[:,2]
    cosmall = [np.cos(m*sph2) for m in range(llim)]
    sinmall = [np.sin(m*sph2) for m in range(llim)]
    Phi = np.zeros((2,nlim,len(r)))

    for ll in xrange(llim):
        gegen1 = gegen_n(nlim, 2*ll+1.5, xi)
        Phifac = -sqrt4pi * np.power(r,ll) * np.power(rp1inv,2*ll+1)
        Phi[0] = Phifac * gegen1
        Phi[1] = Phifac* ((8*ll+6)* np.power(rp1inv,2) * gegen_nm1(nlim, 2*ll+2.5, xi) \
                 + (ll*rinv-(2*ll+1)*rp1inv) * gegen1)

        for mm in xrange(0,ll+1):
            #### New summation over n using numpy.einsum
            Plmtheta,Plmthetap = genlegendrediff(ll,mm,sph[:,1])
            CDEF = np.einsum('ijk,jm->imk', Phi, Knlm[:nlim,ll,mm])

            #print ll,mm,Plmtheta[0],Plmthetap[0]
            cosmphi,sinmphi = cosmall[mm],sinmall[mm]
            pot +=  Plmtheta * (CDEF[0,0]*cosmphi+ CDEF[0,1]*sinmphi )
            force[0] -= Plmtheta *(CDEF[1,0]* cosmphi+ CDEF[1,1]*sinmphi )
            force[1] += rinv* sintheta * Plmthetap * (CDEF[0,0]*cosmphi +CDEF[0,1]*sinmphi)
            if mm:
                force[2] -= rinv* mm * Plmtheta*sinthetainv * (CDEF[0,1]*cosmphi - CDEF[0,0]*sinmphi)

    return force.T/(a**2), pot/a #*G ## Normalization to account for G and to scale our results back to a!=1



def nbodypot(xyz, mass, outpos):

    """
    Calculate potential and force using N-body interation
    """
    assert xyz.shape[1] == 3
    assert outpos.shape[1] == 3

    pot = np.zeros(len(outpos))
    force = np.zeros((len(outpos),3))

    for i in xrange(len(outpos)):
        if i%100==0: print(i)

        dr = xyz - outpos[i]
        d2 = ne.evaluate("sum(dr**2, axis=1)")
        sel = d2 > 0
        invr = ne.evaluate("1./sqrt(d2)")[sel, np.newaxis]

        m = mass[sel, np.newaxis]
        #print(m.shape, invr.shape, dr.shape)
        pot[i] = -ne.evaluate("sum(m*invr)")
        dr = dr[sel]
        force[i] = ne.evaluate("sum(m*invr**3 * dr, axis=0)")

    print
    return pot, force

def nbody(xyz):
    pot=np.zeros(len(xyz))
    force =np.zeros((len(xyz),3))
    for i in xrange(len(xyz)):
        if i%200==0:
            print(i)
        for j in xrange(i+1,len(xyz)):
            dr=xyz[j]-xyz[i]
            r=np.sqrt((dr**2).sum())
            rinv=1./r
            pot[i]-=rinv
            pot[j]-=rinv
            force[i]+=rinv**3 * dr
            force[j]-=rinv**3 * dr

    return force,pot

## For Subhalos
def image(a,b,boxsize):#a and b are 3-d vectors
    diff=b-a
    b[(b-a)>(boxsize/2.)]-=boxsize
    b[(b-a)<(-boxsize/2.)]+=boxsize
    return b

def shortestdist(a,b,boxsize):#a and b are 3-d vectors
    if len(b.shape)==1:
        return np.sqrt(((image(a,b,boxsize)-a)**2).sum())
    else:
        return np.sqrt(((image(a,b,boxsize)-a)**2).sum(axis=1))

def getsubhalos(cat,group):
    grouppos = cat.GroupPos[group]
    count = np.array([],dtype='uint32')
    for j in xrange(1,cat.GroupNsubs[group]):
        sub = j+cat.GroupFirstSub[group]
        if cat.SubhaloLen[sub] >= 100:
            dist=shortestdist(grouppos,cat.SubhaloPos[sub],75000)
            if (dist <= cat.Group_R_Crit200[group]):
                count=np.append(count,sub)
    return count
