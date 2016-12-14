import numpy as np
import matplotlib.pyplot as plt
import tables
from scipy import stats
import os,sys, getopt
sys.path.append('/n/ghernquist/kchua/Orbit/201-code-C-Mar2016/python')
import classify

#plt.interactive(False)
def intclass(filedir,components,doclass=True,makeplot=False):
    if doclass:
        classify.classall(filedir+'/save_temp.hdf5',filedir+'/classout_temp.hdf5',npoints=16384,\
                components=components)

    with tables.open_file(filedir+'/classout_temp.hdf5','r') as f:
        E=f.root.totE[:,0]
        C=f.root.classorb[:,1]
        def fitfunc(x,a,b,c):
            return a*(abs(x)/1.e4)**2 + b*abs(x)/1.e4 + c
        #plt.plot(abs(E),10**fitfunc(abs(E[C!=0]),p[0],p[1],p[2]),'-')
        lin=stats.linregress(abs(E[C!=0])/1.e4,np.log10(C[C!=0]))
        with open(filedir+'params.txt','w') as paramf:
            paramf.write('%.10f \t %.10f' % (lin[0],lin[1]))

        if makeplot:
            fig,ax=plt.subplots(2,2,figsize=(12,10))

            ax[0,0].semilogy(abs(f.root.totE[:,0]),f.root.classorb[:,1],'.',ms=1)
            ax[0,0].plot(abs(E),10**(lin[0]*abs(E)/1e4+lin[1]),'-')
            ax[0,0].set_xlabel('abs(E)')
            ax[0,0].set_ylabel('Test Count')

            ax[0,1].plot(abs(E),C/10**(lin[0]*abs(E)/1e4+lin[1])*70,'.')
            #The below statement is no longer true.
            #This line is actually not accurate since I'm saving only particular intervals
            #The actual prediction requires accounting for a floor in multiples of savenpoints
            #ax[0,1].plot(abs(E),C*np.floor(70/10**(lin[0]*abs(E)/1e4+lin[1])*5e18/3e13/16384)*16384*3e13/5e18,'.')
            ax[0,1].set_xlabel('abs(E)')
            ax[0,1].set_ylabel('Predicted Count')

            ax[1,0].semilogy(abs(f.root.totE[:,0]),f.root.avgdist[:,1],'.',ms=1)
            ax[1,0].set_xlabel('abs(E)')
            ax[1,0].set_ylabel('R_avg')

            ax[1,1].plot(f.root.classorb[:,1],1-f.root.totE[:,1]/f.root.totE[:,0],'.')
            ax[1,1].set_xlabel('Test Count')
            ax[1,1].set_ylabel(r'$1-\Delta E/E_0$')

            fig.savefig(filedir+'/intermediate_results.pdf',bbox_true=True)

def finclass(filedir,components,makeplot=False):
    classify.classall(filedir+'/save.hdf5',filedir+'/classout.hdf5',npoints=16384,\
            components=components)

def plotsingle(filedir):
    if filedir[-5:]=='.hdf5':
        filename=filedir
    else:
        filename=filedir+'/classout.hdf5'
    with tables.open_file(filename,'r') as f:
        E=f.root.totE[:,0]
        C=f.root.classorb[:,1]

        fig,ax=plt.subplots(2,2,figsize=(12,10))

        ax[0,0].plot(abs(f.root.totE[:,0]),f.root.classorb[:,1],'.',ms=1)
        ax[0,0].set_xlabel('abs(E)')
        ax[0,0].set_ylabel('Count')

        ax[0,1].loglog(f.root.avgdist[:,1],f.root.avgdist[:,0],'.',ms=1)
        ax[0,1].set_xlabel('R_avg')
        ax[0,1].set_ylabel('R_0')

        ax[1,0].semilogy(abs(f.root.totE[:,0]),f.root.avgdist[:,1],'.',ms=1)
        ax[1,0].set_xlabel('abs(E)')
        ax[1,0].set_ylabel('R_avg')

        ax[1,1].plot(f.root.classorb[:,1],1-f.root.totE[:,1]/f.root.totE[:,0],'.')
        ax[1,1].set_xlabel('Test Count')
        ax[1,1].set_ylabel(r'$1-\Delta E/E_0$')

        fig.tight_layout()
        #fig.savefig(filedir+'/check_results.pdf',bbox_true=True)

def plotall(filedir,whichdist,nbins=10,components=1):
    import glob
    per=np.percentile
    allgroups=glob.glob(filedir+'/GROUP_*/')
    allbox=np.zeros((len(allgroups),nbins))
    alltube=np.zeros((len(allgroups),nbins))
    allirr=np.zeros((len(allgroups),nbins))
    edges=np.logspace(-1.7,0,nbins+1)
    x=(edges[:-1]+edges[1:])/2

    for i in xrange(len(allgroups)):
        dirt=allgroups[i]
        if components==1:
            varfile='var.hdf5'
        elif components==2:
            varfile='varh.hdf5'
        with tables.open_file(dirt+varfile,'r') as var:
            rvir=var.root.Rvir[:]/classify.h*classify.kpc
        with tables.open_file(dirt+'/classout.hdf5','r') as f:
            dist=f.root.avgdist[:,whichdist]/rvir
            C=f.root.classification[:]

        #print dist.min()
        wbox=np.where( (C==4) | (C==5))[0]
        wtube=np.where( (C==1) | (C==2) | (C==3))[0]
        wirr=np.where(C==0)[0]
        box=np.array([bin.sum() for bin in ((edges[i]<dist[wbox]) & (dist[wbox]<edges[i+1])  for i in np.arange(len(edges)-1))],dtype='float')
        tube=np.array([bin.sum() for bin in ((edges[i]<dist[wtube]) & (dist[wtube]<edges[i+1])  for i in np.arange(len(edges)-1))],dtype='float')
        irr=np.array([bin.sum() for bin in ((edges[i]<dist[wirr]) & (dist[wirr]<edges[i+1])  for i in np.arange(len(edges)-1))],dtype='float')
        allbox[i]=box/(box+tube+irr)
        alltube[i]=tube/(box+tube+irr)
        allirr[i]=1.-(box+tube)/(box+tube+irr)

    #print allbox
    plt.figure()
    plt.semilogx(x,np.median(allbox,axis=0),'k-')
    plt.semilogx(x,np.median(alltube,axis=0),'b-')
    plt.semilogx(x,np.median(allirr,axis=0),'r-')

    plt.legend(['Box','Tube','Irr'],loc=2)

    plt.fill_between(x,per(allbox,25,axis=0),per(allbox,75,axis=0),color='k',alpha=0.1)
    plt.fill_between(x,per(alltube,25,axis=0),per(alltube,75,axis=0),color='b',alpha=0.1)
    plt.fill_between(x,per(allirr,25,axis=0),per(allirr,75,axis=0),color='r',alpha=0.1)

    plt.axvline(0.5,c='k',ls='--')
    plt.xlabel(r'$r/R_{200}$')
    plt.ylabel('Orbit Fraction')


if __name__ == "__main__":
    indir=''
    step=''
    c=1
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hs:c:",["step=","components="])
    except getopt.GetoptError:
        print 'intermediate.py -s first/second -c 1/2 <i/o dir>'
        sys.exit(2)
    for opt,arg in opts:
        if opt == '-h':
            print 'intermediate.py -s first/second -c 1/2 <i/o dir>'
            sys.exit()
        elif opt in ("-s", "--step"):
            step = arg
        elif opt in ("-c", "--components"):
            c=int(arg)
            print "no. of components = ",c

    dir = args[0]

    if step == 'first':
        intclass(dir,c,doclass=True,makeplot=True)
    elif step == 'second':
        finclass(dir,c)
    elif step == 'plotall':
        plotall(dir,whichdist=1,nbins=15,components=c)
    elif step =='plotsingle':
        plotsingle(dir)
    else:
        print 'Invalid specification of step'
        sys.exit(2)
