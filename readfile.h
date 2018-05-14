#ifndef READFILE_H_
#define READFILE_H_

class Readparams{
    public: 
        double endTime, dt, tol;
        int nPoints, Npart, saveint, NumComponents=1;
        int NumSavepoints=16384; //no. of time points to be saved in output file per particle
        bool verbose, twocomp=false, firstpass, test;
        char *infilename1, *infilename2, *outfilename, *initfilename;
        int read(int argc, char *argv[]);
        void operator() (int argc, char *argv[]){
            read(argc,argv);
        }
};  

#endif
