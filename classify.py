import numpy as np
import tables
from multiprocessing import Pool, Array
import ctypes
from functools import partial
import os,sys
sys.path.append('/n/ghernquist/kchua/Orbit/201-code-C-Mar2016/compare')
sys.path.append('/n/ghernquist/kchua/Orbit/201-code-C-Mar2016/taxon')
import taxon
import taxonopt
import taxon32768
import taxon8192
import solve
import matplotlib.pyplot as plt

Myr = 3.15569e16 #s

def classall(infile,outfile='classout.hdf5',rotmat=False,npoints=16384,nn=16384,end=False):
    ## N = no. of particles in each bin
    ## nn = no. of points to be used in classification
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
        distout[:,0]=np.sqrt((x[:,0,:3]**2).sum(axis=1))
        for i in xrange(npart): #calculate average distance
            xx=x[i,:,:3]
            vv=x[i,:,3:]
            d=np.sqrt((xx**2).sum(axis=1))
            vr=abs(np.einsum('ij,ij->i',xx,vv)/d)
            distout[i,1]=np.sum(d/vr)/np.sum(1./vr)

        A=solve.Problem('var') #calculate initial and final energy
        file.root.totE[:,0]=A.potential(x[:,0,:3])[:,0]+(x[:,0,3:]**2).sum(axis=1)/2.
        file.root.totE[:,1]=A.potential(x[:,-1,:3])[:,0]+(x[:,-1,3:]**2).sum(axis=1)/2.
        del A

        t_base = Array(ctypes.c_double, npart*nn)
        tt = np.ctypeslib.as_array(t_base.get_obj()).reshape(npart,nn)
        x_base = Array(ctypes.c_double, npart*nn*3)
        xx = np.ctypeslib.as_array(x_base.get_obj()).reshape(npart,nn,3)
        v_base = Array(ctypes.c_double, npart*nn*3)
        vv = np.ctypeslib.as_array(v_base.get_obj()).reshape(npart,nn,3)
        tt[:]=t[:,::interval][:,:nn]
        xx[:]=x[:,::interval,:3][:,:nn]
        vv[:]=x[:,::interval,3:][:,:nn]
        del t
        del x

        # Use 8 threads
        pool=Pool(8)
        partialfunc=partial(calcclass,t=tt,x=xx,v=vv,nn=nn)
        out=pool.map(partialfunc,xrange(npart))
        pool.close()
        pool.join()
        file.root.classorb[:]=np.array(out)

        pool2=Pool(8)
        out=pool2.map(getclass,file.root.classorb[:,0])
        pool2.close()
        pool2.join()
        file.root.classification[:]=out

def calcclass(i,t,x,v,nn):
    tt=np.asfortranarray(t[i])
    xx=np.asfortranarray(x[i].transpose())
    vv=np.asfortranarray(v[i].transpose())*Myr
    if nn==32768:
        out=taxon32768.taxon(tt,xx,vv,nn,1,3,0,0,1,0,0,'test')
    elif nn==16384:
        out=taxon.taxon(tt,xx,vv,nn,1,3,0,0,1,0,0,'test')
    elif nn==8192:
        out=taxon8192.taxon(tt,xx,vv,nn,1,3,0,0,1,0,0,'test')
    return out

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

def plotclass(dist,C,nbins=10):
    edges=np.logspace(-2,0,nbins+1)
    x=(edges[:-1]+edges[1:])/2
    wbox=np.where( (C==4) | (C==5))[0]
    wtube=np.where( (C==1) | (C==2) | (C==3))[0]
    wirr=np.where(C==0)[0]
    box=np.array([bin.sum() for bin in ((edges[i]<dist[wbox]) & (dist[wbox]<edges[i+1])  for i in np.arange(len(edges)-1))],dtype='float')
    tube=np.array([bin.sum() for bin in ((edges[i]<dist[wtube]) & (dist[wtube]<edges[i+1])  for i in np.arange(len(edges)-1))],dtype='float')
    irr=np.array([bin.sum() for bin in ((edges[i]<dist[wirr]) & (dist[wirr]<edges[i+1])  for i in np.arange(len(edges)-1))],dtype='float')
    plt.figure()
    plt.semilogx(x,box/(box+tube+irr),'k-')
    plt.semilogx(x,tube/(box+tube+irr),'b-')
    plt.semilogx(x,irr/(box+tube+irr),'r-')
