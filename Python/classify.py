import numpy as np
import tables
#from multiprocessing import Pool, Array
import ctypes
from functools import partial
import os,sys
sys.path.append('/n/ghernquist/kchua/Orbit/201-code-C-Mar2016/taxon')
#import taxon
#import taxon32768
#import taxon8192
#import solve
import matplotlib.pyplot as plt
import sharedmem

Myr = 3.15569e16 #s
kpc=3.08568e+16 #km
h=0.704

def classall(infile,components,rotation,outfile='classout.hdf5',npoints=16384,nn=16384,lim=[12,6],end=False):
    sys.path.append('/n/ghernquist/kchua/Orbit/201-code-C-Mar2016/Python')
    import taxon
    import solve
    def calcclass(i,t,x,v,nn):
        #taxon(t,x,v,n,jsub,jdim,jcla,jcl,jpan,jlin,jcom,arch)
        tt=np.asfortranarray(t[i])
        xx=np.asfortranarray(x[i].transpose())
        vv=np.asfortranarray(v[i].transpose())*Myr
        if nn==32768:
            out=taxon32768.taxon(tt,xx,vv,nn,1,3,0,0,0,0,0,'test')
        elif nn==16384:
            #out=taxon.taxon(tt,xx,vv,nn,1,3,0,0,0,0,0,'test')
            out=taxon.taxon(tt,xx,vv,nn,1,3,0,0,0,0,0,'test')
        elif nn==8192:
            out=taxon8192.taxon(tt,xx,vv,nn,1,3,0,0,0,0,0,'test')
        return out
    ## N = no. of particles in each bin
    ## nn = no. of points to be used in classification
    with tables.open_file('shape.hdf5','r') as rotfile:
        rotmat=rotfile.root.rotmat[rotation]
    with tables.open_file(infile,'r') as u:
        #t=u.root.t[:]/Myr
        npart=u.root.x.shape[0]/(npoints+1)/6
        print 'analyzing orbits for',infile
        classout=np.zeros((npart,2))
        x=u.root.x[:].reshape(npart,npoints+1,6)
        t=u.root.t[:].reshape(npart,npoints+1)/Myr
    if end:
        print 'using only ',end,' points'
        t=t[:end]
        x=x[:,:end,:]
    else:
        print 'using all points'

    interval=npoints/nn
    print 'interval = ', interval
    with tables.open_file(outfile,'w') as file:
        file.create_carray("/","classorb",tables.Int32Col(),(npart,2))
        file.create_carray("/","classification",tables.Int32Col(),(npart,))
        distout=file.create_carray("/","avgdist",tables.Float64Col(),(npart,2))
        file.create_carray("/","totE",tables.Float64Col(),(npart,2))
        file.create_carray("/","Ncomponents",tables.Int32Atom(),(1,))
        file.root.Ncomponents[0]=components
        v=sharedmem.copy(x[:,:,3:])
        x=sharedmem.copy(x[:,:,:3])

        distout[:,0]=np.sqrt((x[:,0]**2).sum(axis=1))
        with sharedmem.MapReduce() as pool:
            def avgr(i,xx,vv):
                d=np.sqrt((xx[i]**2).sum(axis=1))
                vr=abs(np.einsum('ij,ij->i',xx[i],vv[i])/d)
                return np.sum(d/vr)/np.sum(1./vr)
            partialfunc=partial(avgr,xx=x,vv=v)
            out=pool.map(partialfunc,xrange(npart))
            distout[:,1]=out

        if components==2:
            A=solve.Problem('varh',nlim=lim[0],llim=lim[1])
            B=solve.Problem('varc',nlim=lim[0],llim=lim[1])
            file.root.totE[:,0]=A.potential(x[:, 0])[:,0]+ B.potential(x[:, 0])[:,0] + \
                    (v[:,0]**2).sum(axis=1)/2.
            file.root.totE[:,1]=A.potential(x[:,-1])[:,0]+ B.potential(x[:,-1])[:,0] + \
                    (v[:,-1]**2).sum(axis=1)/2.
            del A; del B
        elif components==1:
            A=solve.Problem('var',nlim=lim[0],llim=lim[1]) #calculate initial and final energy
            file.root.totE[:,0]=A.potential(x[:,0])[:,0]+(v[:,0]**2).sum(axis=1)/2.
            file.root.totE[:,1]=A.potential(x[:,-1])[:,0]+(v[:,-1]**2).sum(axis=1)/2.
            del A
        elif components==3:
            A=solve.Problem('varc',nlim=lim[0],llim=lim[1]) #calculate initial and final energy
            file.root.totE[:,0]=A.potential(x[:,0])[:,0]+(v[:,0]**2).sum(axis=1)/2.
            file.root.totE[:,1]=A.potential(x[:,-1])[:,0]+(v[:,-1]**2).sum(axis=1)/2.
            del A
        else:
            print "Number of components not specified"

        #t_base = Array(ctypes.c_double, npart*nn)
        #tt = np.ctypeslib.as_array(t_base.get_obj()).reshape(npart,nn)
        #x_base = Array(ctypes.c_double, npart*nn*3)
        #xx = np.ctypeslib.as_array(x_base.get_obj()).reshape(npart,nn,3)
        #v_base = Array(ctypes.c_double, npart*nn*3)
        #vv = np.ctypeslib.as_array(v_base.get_obj()).reshape(npart,nn,3)
        #tt[:]=t[:,::interval][:,:nn]
        #xx[:]=x[:,::interval,:3][:,:nn]
        #vv[:]=x[:,::interval,3:][:,:nn]
        #del t
        #del x
        t=sharedmem.copy(t[:,::interval][:,:nn])
        x=x[:,::interval][:,:nn]
        v=v[:,::interval][:,:nn]
        x=np.dot(x,rotmat)
        v=np.dot(v,rotmat)
        # Use 8 threads
        #pool=Pool(8)
        with sharedmem.MapReduce() as pool:
            partialfunc=partial(calcclass,t=t,x=x,v=v,nn=nn)
            out=pool.map(partialfunc,xrange(npart))
        #pool.close()
        #pool.join()
            file.root.classorb[:]=np.array(out)

        with sharedmem.MapReduce() as pool:
            out=pool.map(getclass,file.root.classorb[:,0])
            file.root.classification[:]=out

#def calcclass(i,t,x,v,nn):
#    #taxon(t,x,v,n,jsub,jdim,jcla,jcl,jpan,jlin,jcom,arch)
#    tt=np.asfortranarray(t[i])
#    xx=np.asfortranarray(x[i].transpose())
#    vv=np.asfortranarray(v[i].transpose())*Myr
#    if nn==32768:
#        out=taxon32768.taxon(tt,xx,vv,nn,1,3,0,0,0,0,0,'test')
#    elif nn==16384:
#        #out=taxon.taxon(tt,xx,vv,nn,1,3,0,0,0,0,0,'test')
#        out=taxon.taxon(tt,xx,vv,nn,1,3,0,0,0,0,0,'test')
#    elif nn==8192:
#        out=taxon8192.taxon(tt,xx,vv,nn,1,3,0,0,0,0,0,'test')
#    return out

def getclass(classorb):
    xtube=[121,211,221,311,321]
    ytube=[122,212,222,312,322]
    ztube=[123,213,223,313,323]
    box=[300]
    resbox=[310,320,210,120,220]
    All=[[],xtube,ytube,ztube,box,resbox]
    for n in xrange(len(All)):
        if classorb in All[n]:
            return n
    return 0

def avgdist(filename):
    with tables.open_file(filename,'r') as u:
        t=u.root.t[:]
        npart=u.root.x.shape[0]/len(t)/6
        distout=np.zeros((npart,2))
        x=u.root.x[:].reshape(npart,len(t),6)
        distout[:,0]=np.sqrt((x[:,0,:3]**2).sum(axis=1))
        for i in xrange(npart):
            xx=x[i,:,:3]
            vv=x[i,:,3:]
            d=np.sqrt((xx**2).sum(axis=1))
            vr=abs(np.einsum('ij,ij->i',xx,vv)/d)
            distout[i,1]=np.sum(d/vr)/np.sum(1./vr)
        return distout

def getperiod(x):
    A=np.zeros((x.shape[0],5),dtype='int32')
    print A.shape
    for n in np.arange(x.shape[0]):
        for i in np.arange(3):
            temp=np.fft.fft(x[n,:,i]+1j*x[n,:,i+3]*Myr)
            temp=np.abs(temp)
            temp=np.argsort(temp)[::-1][0]
            if temp>x.shape[1]/2:
                temp=x.shape[1]-temp
            A[n,i+1]=temp
    A[:,0]=np.median(A[:,1:4],axis=1)
    return A
    #return np.argsort(A)[::-1][0]

def extractlines(FX):
    N=FX.shape[0]
    nlim=10

    def h(mm):
        m=mm.copy()
        if m.__class__==np.int64:
            m=np.array([m])
        m[m>N/2]=m[m>N/2]-N
        return mm

    def extractline(Fx):
        Ax=np.abs(Fx)
        m1=h(np.argsort(Ax)[::-1][0])
        if Ax[h(m1-1)]>Ax[h(m1+1)]:
            m2=m1-1
        else:
            m2=m1+1
        k=np.pi/N*(m2-m1)
        #a=np.arctan2(np.sin(k),np.cos(k)+Ax[m1]/Ax[m2])
        a=np.arctan2(np.sin(k),(np.cos(k)+Ax[h(m1)]/Ax[h(m2)]))
        s=a*N/np.pi+m1

        rho = lambda m:1-np.cos(2*np.pi*(s-m)/N) - np.cos(2*np.pi*(s-m)) + np.cos(2*np.pi*(s-m)*(N-1)/N)
        sigma = lambda m:np.sin(2*np.pi*(s-m)/N) - np.sin(2*np.pi*(s-m)) + np.sin(2*np.pi*(s-m)*(N-1)/N)
        cosN=1.-np.cos(2*np.pi*(s-m1)/N)
        A=N*Ax[h(m1)]*np.sqrt(cosN/(1.-np.cos(2*np.pi*(s-m1))  ))

        p=rho(m1)
        sg=sigma(m1)
        g=np.angle(Fx[h(m1)])
        phi=np.arctan2( (cosN* (p*np.sin(g) - sg*np.cos(g)) ),
                        (cosN* (p*np.cos(g) + sg*np.sin(g)) ))

        line= lambda m: A*np.exp(1j*phi)*(rho(m)+1j*sigma(m))/(2.*N*(1.-np.cos(2*np.pi*(s-m)/N)))

        l=line(np.arange(N))
        if m1>N/2:
            l[(m1+1):]=-l[(m1+1):]
        else:
            l[:m1]=-l[:m1]

        #l=-l
        #l[m1]=-l[m1]
        print Fx[h(np.arange(m1-1,m1+2))]
        print l[np.arange(m1-1,m1+2)]

        Fx=Fx-l

        return Fx,[s,A,phi]

    #Fx=np.fft.fft(x)
    Fx=FX.copy()
    P=np.zeros((nlim,3))
    for i in np.arange(nlim):
        print i
        Fx,P[i]=extractline(Fx)

    return P






#def plotclass(dist,C,nbins=10):
def plotclass(dir,dir2='',dir3='',nbins=10):
    import string
    if dir[-5:]=='.hdf5':
        dirt=string.join(dir.split('/')[:-1],'/')
    else:
        dirt=dir
    if os.path.isfile(dirt+'/var.hdf5'):
        varfile=dirt+'/var.hdf5'
    elif os.path.isfile(dirt+'/varh.hdf5'):
        varfile=dirt+'/varh.hdf5'
    else:
        print 'Variable file not found'
        sys.exit(2)
    with tables.open_file(varfile,'r') as var:
        rvir=var.root.Rvir[:]/h*kpc
    edges=np.logspace(-1.7,0,nbins+1)
    x=(edges[:-1]+edges[1:])/2
    edges*=rvir

    lines=['-','--',':']
    dirs=[dir]
    if dir2: dirs=[dir,dir2]
    if dir3: dirs=[dir,dir2,dir3]
    print dirs
    plt.figure()
    for i in xrange(len(dirs)):
        dirt=dirs[i]
        if dirt[-5:]!='.hdf5':
            dirt=str(dirt)+'/classout.hdf5.'
        with tables.open_file(dirt,'r') as f:
            dist=f.root.avgdist[:,1]
            C=f.root.classification[:]
        wbox=np.where( (C==4) | (C==5))[0]
        wtube=np.where( (C==1) | (C==2) | (C==3))[0]
        wirr=np.where(C==0)[0]
        box=np.array([bin.sum() for bin in ((edges[i]<dist[wbox]) & (dist[wbox]<edges[i+1])  for i in np.arange(len(edges)-1))],dtype='float')
        tube=np.array([bin.sum() for bin in ((edges[i]<dist[wtube]) & (dist[wtube]<edges[i+1])  for i in np.arange(len(edges)-1))],dtype='float')
        irr=np.array([bin.sum() for bin in ((edges[i]<dist[wirr]) & (dist[wirr]<edges[i+1])  for i in np.arange(len(edges)-1))],dtype='float')
        plt.semilogx(x[box>3],(box/(box+tube+irr))[box>3],c='k',ls=lines[i])
        plt.semilogx(x[tube>3],(tube/(box+tube+irr))[tube>3],c='b',ls=lines[i])
        plt.semilogx(x,1.-(box+tube)/(box+tube+irr),c='r',ls=lines[i])
