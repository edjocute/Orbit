import numpy as np
import tables
from scipy import stats
import os,sys, getopt
sys.path.append('/n/ghernquist/kchua/Orbit/201-code-C-Mar2016/python')
import classify

varfile={1:'var.hdf5',2:'varh.hdf5',3:'varc.hdf5'}

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
        assert Ncom in [1,2,3]

        with tables.open_file(dirt+varfile[Ncom],'r') as var:
            rvir=var.root.Rvir[:]/classify.h*classify.kpc
        dist=f.root.avgdist[:,whichdist]/rvir
        C=f.root.classification[:]

        edges=np.logspace(np.log10(range[0]),np.log10(range[1]),nbins+1)
        x=(edges[:-1]+edges[1:])/2
        allfrac=calcfrac(C,dist,edges)

        col=['k','b','r','b','b','b']
        ls=['-','-','-','--','-.',':']

        fig,ax=plt.subplots(2,1,figsize=(6,10))
        for i in np.arange(6):
            ax[0].semilogx(x,allfrac[i],c=col[i],ls=ls[i])
        ax[0].legend(['Box','Tube','Irr'],loc=2)
        ax[0].set_xlabel(r'$r/R_{200}$')
        ax[0].set_ylabel('Orbit Fraction')
        ax[1].loglog(dist,1-f.root.totE[:,1]/f.root.totE[:,0],'.')
        ax[1].set_xlabel('Avg_Dist')
        ax[1].set_ylabel(r'$1-\Delta E/E_0$')
        fig.tight_layout()
        return fig

def plotall(filedir,components,nbins=10,whichdist=1,\
        filename='classout.hdf5',ls='-',ploterr=True,plottubes=False,range=[0.05,1.]):
    import matplotlib.pyplot as plt
    f,ax=plt.subplots(1,1)
    plotdir(ax,filedir,components,nbins,whichdist,filename,ls,ploterr,plottubes,range)
    return f

def plotdir(ax,filedir,components,nbins=10,whichdist=1,\
        filename='classout.hdf5',ls='-',ploterr=True,plottubes=False,range=[0.05,1.]):
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
    allfrac=np.zeros((len(allgroups),6,nbins))

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
            assert Ncom in [1,2,3]

            with tables.open_file(dirt+varfile[Ncom],'r') as var:
                rvir=var.root.Rvir[:]/classify.h*classify.kpc
                dist=f.root.avgdist[:,whichdist]/rvir
                C=f.root.classification[:]
        allfrac[i]=calcfrac(C,dist,edges)

    col=['k','b','r','b','b','b']
    lines=['-','-','-','--','-.',':']
    if plottubes==True:
        for i in np.arange(6):
            ax.semilogx(x,np.median(allfrac[:,i],axis=0),c=col[i],ls=lines[i],lw=2)
        leg=ax.legend(['Box','Tube','Irr','Xtube','Ytube','Ztube'],loc=3)
    elif plottubes==False:
        for i in np.arange(3):
            ax.semilogx(x,np.median(allfrac[:,i],axis=0),c=col[i],ls=ls,lw=2)
        leg=ax.legend(['Box','Tube','Irr'],loc=3)
    ax.add_artist(leg)

    if ploterr==True:
        for i in np.arange(3):
            w=allfrac[:,i]
            ax.fill_between(x,per(w,25,axis=0),per(w,75,axis=0),\
                    color=col[i],alpha=0.1)
    ax.set_xlabel(r'$r/R_{200}$')
    ax.set_ylabel('Orbit Fraction')

def pdfmulti(dir,components,filename='classout3.hdf5',range=[0.05,1.]):
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    import glob

    alldir=glob.glob(dir+'/GROUP_*')

    with PdfPages(dir+'allplots.pdf') as pdf:
        fig=plotall(dir,components,filename='classout3.hdf5',ploterr=True,plottubes=True,range=range)
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
    wx=np.where(C==1)[0]
    wy=np.where(C==2)[0]
    wz=np.where(C==3)[0]

    nbox=np.array([bin.sum() for bin in ((edges[i]<dist[wbox]) & (dist[wbox]<edges[i+1])  for i in np.arange(len(edges)-1))])
    ntube=np.array([bin.sum() for bin in ((edges[i]<dist[wtube]) & (dist[wtube]<edges[i+1])  for i in np.arange(len(edges)-1))])
    nirr=np.array([bin.sum() for bin in ((edges[i]<dist[wirr]) & (dist[wirr]<edges[i+1])  for i in np.arange(len(edges)-1))])
    ntot=(nbox+ntube+nirr).astype('float')

    nx=np.array([bin.sum() for bin in ((edges[i]<dist[wx]) & (dist[wx]<edges[i+1])  for i in np.arange(len(edges)-1))])
    ny=np.array([bin.sum() for bin in ((edges[i]<dist[wy]) & (dist[wy]<edges[i+1])  for i in np.arange(len(edges)-1))])
    nz=np.array([bin.sum() for bin in ((edges[i]<dist[wz]) & (dist[wz]<edges[i+1])  for i in np.arange(len(edges)-1))])

    fbox=nbox/ntot
    ftube=ntube/ntot
    firr=1.-(fbox+ftube)
    fx=nx/ntot
    fy=ny/ntot
    fz=nz/ntot
    try:
        assert np.allclose(fx+fy+fz,ftube)
    except AssertionError:
        print fx,fy,fz,ftube

    return np.vstack((fbox,ftube,firr,fx,fy,fz))

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
        pdfmulti(dir,c,f,range=[0.05,1.])
    elif step == 'plotallstars':
        pdfmulti(dir,c,f,range=[0.05,0.5])
    elif step =='plotsingle':
        plotsingle(dir)
    else:
        print 'Invalid specification of step'
        sys.exit(2)
