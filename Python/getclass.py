import numpy as np
import tables
from scipy import stats
import os,sys, getopt
sys.path.append('/n/ghernquist/kchua/Orbit/201-code-C-Mar2016/python')
import classify

varfile=['var.hdf5','varh.hdf5']

#plt.interactive(False)
def intclass(filedir,components,rotation,doclass=True,makeplot=False):
    if doclass:
        classify.classall(filedir+'/save_temp.hdf5',components,rotation, \
                outfile=filedir+'/classout_temp.hdf5',npoints=16384)

    with tables.open_file(filedir+'/classout_temp.hdf5','r') as f:
        E=f.root.totE[:,0]
        C=f.root.classorb[:,1]
        def fitfunc(x,a,b,c):
            return a*(abs(x)/1.e4)**2 + b*abs(x)/1.e4 + c
        #plt.plot(abs(E),10**fitfunc(abs(E[C!=0]),p[0],p[1],p[2]),'-')
        lin=stats.linregress(abs(E[C!=0])/1.e4,np.log10(C[C!=0]))
        with open(filedir+'params.txt','w') as paramf:
            paramf.write('%.10f \t %.10f' % (lin[0],lin[1]))
    print lin[0],lin[1]

    if makeplot:
        fig=plotsingle(filedir+'/classout_temp.hdf5',lin=lin)
        fig.savefig(filedir+'/intermediate_results.pdf',bbox_tight=True)

def finclass(filedir,components,rotation,outfile,makeplot=False):
    classify.classall(filedir+'/save.hdf5',components,rotation,\
            outfile=filedir+outfile,npoints=16384)

def plotsingle(filedir,lin=False):
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    if filedir[-5:]=='.hdf5':
        filename=filedir
    else:
        filename=filedir+'/classout_temp.hdf5'
    with tables.open_file(filename,'r') as f:
        E=f.root.totE[:,0]
        C=f.root.classorb[:,1]

        fig,ax=plt.subplots(2,2,figsize=(12,10))

        ax[0,0].semilogy(abs(f.root.totE[:,0]),f.root.classorb[:,1],'.',ms=1)
        if lin:
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

        fig.tight_layout()
        return fig

def plotfin(filedir,components,nbins=10,whichdist=1,ls='-',\
        range=[0.05,1.]):
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

    if filedir[-5:]=='.hdf5':
        filename=filedir
        dirt=filedir.split('classout')[0]
    else:
        filename=filedir+'/classout.hdf5'

    with tables.open_file(filename,'r') as f:
        try:
            Ncom=f.root.Ncomponents[0]
        except tables.NoSuchNodeError:
            Ncom=components
        assert Ncom in [0,1]

        with tables.open_file(dirt+varfile[Ncom],'r') as var:
            rvir=var.root.Rvir[:]/classify.h*classify.kpc
        dist=f.root.avgdist[:,whichdist]/rvir
        C=f.root.classification[:]

        edges=np.logspace(np.log10(range[0]),np.log10(range[1]),nbins+1)
        x=(edges[:-1]+edges[1:])/2
        allbox,alltube,allirr=calcfrac(C,dist,edges)

        fig,ax=plt.subplots(2,1,figsize=(6,10))
        for w,c in [[allbox,'k'],[alltube,'b'],[allirr,'r']]:
            ax[0].semilogx(x,w,c=c,ls=ls)
        ax[0].legend(['Box','Tube','Irr'],loc=2)
        ax[0].set_xlabel(r'$r/R_{200}$')
        ax[0].set_ylabel('Orbit Fraction')
        ax[1].loglog(dist,1-f.root.totE[:,1]/f.root.totE[:,0],'.')
        ax[1].set_xlabel('Avg_Dist')
        ax[1].set_ylabel(r'$1-\Delta E/E_0$')
        fig.tight_layout()
        return fig

def plotall(filedir,components,nbins=10,whichdist=1,\
        filename='classout.hdf5',ls='-',ploterr=False,range=[0.05,1.]):
    import matplotlib.pyplot as plt
    import glob
    per=np.percentile
    allgroups=glob.glob(filedir+'/GROUP_*/')
    edges=np.logspace(np.log10(range[0]),np.log10(range[1]),nbins+1)
    x=(edges[:-1]+edges[1:])/2

    if len(allgroups)>0:
        N=len(allgroups)
    else:
        N=1
        ploterr=False
    allbox=np.zeros((len(allgroups),nbins))
    alltube=np.zeros((len(allgroups),nbins))
    allirr=np.zeros((len(allgroups),nbins))

    for i in xrange(N):
        if len(allgroups)>0:
            dirt=allgroups[i]
        else:
            dirt='./'
        with tables.open_file(dirt+'/'+filename,'r') as f:
            try:
                Ncom=f.root.Ncomponents[0]
            except tables.NoSuchNodeError:
                Ncom=components
            assert Ncom in [0,1]

            with tables.open_file(dirt+varfile[Ncom],'r') as var:
                rvir=var.root.Rvir[:]/classify.h*classify.kpc
                dist=f.root.avgdist[:,whichdist]/rvir
                C=f.root.classification[:]
        allbox[i],alltube[i],allirr[i]=calcfrac(C,dist,edges)

    fig=plt.figure()
    for w,c in [[allbox,'k'],[alltube,'b'],[allirr,'r']]:
        plt.semilogx(x,np.median(w,axis=0),c=c,ls=ls)
    plt.legend(['Box','Tube','Irr'],loc=2)

    if ploterr==True:
        for w,c in [[allbox,'k'],[alltube,'b'],[allirr,'r']]:
            plt.fill_between(x,per(w,25,axis=0),per(w,75,axis=0),color=c,alpha=0.1)
    #plt.axvline(0.5,c='k',ls='--')
    plt.xlabel(r'$r/R_{200}$')
    plt.ylabel('Orbit Fraction')
    return fig

def pdfmulti(dir,components,filename='classout3.hdf5'):
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    import glob

    alldir=glob.glob(dir+'/GROUP_*')

    with PdfPages(dir+'allplots.pdf') as pdf:
        fig=plotall(dir,components,filename='classout3.hdf5',ploterr=True)
        pdf.savefig(fig)
        plt.close()
        for i in alldir:
            print i
            fig=plotfin(i+'/classout3.hdf5',components)
            pdf.savefig(fig)
            plt.close()

def calcfrac(C,dist,edges):
    wbox=np.where( (C==4) | (C==5))[0]
    wtube=np.where( (C==1) | (C==2) | (C==3))[0]
    wirr=np.where(C==0)[0]
    nbox=np.array([bin.sum() for bin in ((edges[i]<dist[wbox]) & (dist[wbox]<edges[i+1])  for i in np.arange(len(edges)-1))],dtype='float')
    ntube=np.array([bin.sum() for bin in ((edges[i]<dist[wtube]) & (dist[wtube]<edges[i+1])  for i in np.arange(len(edges)-1))],dtype='float')
    nirr=np.array([bin.sum() for bin in ((edges[i]<dist[wirr]) & (dist[wirr]<edges[i+1])  for i in np.arange(len(edges)-1))],dtype='float')
    fbox=nbox/(nbox+ntube+nirr)
    ftube=ntube/(nbox+ntube+nirr)
    firr=1.-(fbox+ftube)
    return fbox,ftube,firr

if __name__ == "__main__":
    step=''
    c=1
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hs:c:r:f:",["step=","components=","rotation=","filename="])
    except getopt.GetoptError:
        print 'intermediate.py -s first/second -c 1/2 <i/o dir> -r # -f filename'
        sys.exit(2)
    for opt,arg in opts:
        if opt == '-h':
            print 'intermediate.py -s first/second -c 1/2 <i/o dir> -r # -f filename'
            sys.exit()
        elif opt in ("-s", "--step"):
            step = arg
            print "step=",step
        elif opt in ("-c", "--components"):
            c=int(arg)
            print "no. of components = ",c
        elif opt in ("-r", "--rotation"):
            r=int(arg)
        elif opt in ("-f", "--filename"):
            f=arg


    dir = args[0]

    if step == 'first':
        intclass(dir,c,r,doclass=True,makeplot=True)
    elif step == 'second':
        finclass(dir,c,r,f)
    elif step == 'plotall':
        plotall(dir,whichdist=1,nbins=15,components=c)
    elif step =='plotsingle':
        plotsingle(dir)
    else:
        print 'Invalid specification of step'
        sys.exit(2)
