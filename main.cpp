#include <math.h>
#include <stdio.h>
#include "main.h"
#include <omp.h>
#include <string>
#include <fstream>

Readparams allparams;

int main( int argc, char* argv[]){
    struct Indata infile,infile2; //to store coefficients
    std::vector<state_type> xinit,xinit_run; //to store array of initial positions and velocities

    allparams.read(argc,argv);

    /* Call either the 1-component or 2-component routine 
     * based on whether 2nd component input is present
     */
    loadHdf5Input(allparams.infilename1, &infile);
    std::cout << "First component read\n";
    if (allparams.twocomp){
        loadHdf5Input(allparams.infilename2, &infile2);
        std::cout << "Second component read\n";
    } 

    loadHdf5Init(allparams.initfilename, xinit);
    //Npart = (vm.count("N")) ? vm["N"].as<int>() : xinit.size();
    /* Check whether 1st or 2nd pass and adjust the integration points accordingly. */


    /* FOR TESTING ONLY:
     * Save test potentials to file and exit */
    if (allparams.test){
        calcOrb Orbit = calcOrb(infile,infile2);
        //calcOrb Orbit = calcOrb(&allparams,infile,infile2);
        Orbit.test(xinit);
        return 0;
    }


    if (allparams.firstpass){
        allparams.Npart=100;
        xinit_run = std::vector<state_type>();
        std::cout << "FIRST PASS. Points to integrate = " << allparams.Npart << "\n";
        std::cout << "T_end=" << allparams.endTime <<'\n';
        int Interval= xinit.size()/allparams.Npart;
        //std::cout << "interval = " << Interval << "\n";
        for (int n=0; n<allparams.Npart; n++){
            xinit_run.push_back(xinit[n*Interval]);
        }
    }
    else{
        std::cout << "SECOND PASS. Integrating all points." << "\n";
        allparams.Npart=xinit.size();
        xinit_run = xinit;
    }

    /* Create output arrays based on number of integration points */
    state_type XX(allparams.Npart*(allparams.NumSavepoints)*6); //output array
    state_type T(allparams.Npart*(allparams.NumSavepoints));
    std::cout << "Total points in input/Points to integrate = " << xinit.size() << "/" << allparams.Npart << '\n';

    //calcAcc ACC = (allparams.twocomp) ? calcAcc(&infile,&infile2) : calcAcc(&infile);
    //calcOrb Orbit = calcOrb(&allparams,infile,infile2);
    calcOrb Orbit = calcOrb(infile,infile2);
    double time1=omp_get_wtime();
    Orbit.integrate(XX,T,xinit_run);
    double time2=omp_get_wtime();

    std::cout << "... done!" << "\n";

    /* Print some info to screen for checking */
    if (allparams.verbose){
    for (int n=0; n< ((allparams.Npart<10) ? allparams.Npart : 10); n++){
        std::cout << n << '\n';
        std::cout << T[n*(allparams.NumSavepoints)] << '\t' ;
        for ( size_t j=0; j<=5; j++){
            std::cout <<  XX[n*6*(allparams.NumSavepoints)+j] << '\t';
        }
        std::cout << '\n';

        std::cout <<  T[(n+1)*(allparams.NumSavepoints)-1] << '\t' ;
        for ( size_t j=0; j<=5; j++){
            std::cout <<  XX[n*6*(allparams.NumSavepoints)+6*(allparams.NumSavepoints)+j] << '\t';
        }
        std::cout << '\n';
    }
    }

    /* Save data to output file */
    saveHdf5(XX,T);
    fprintf(stdout,"All done! This took %.4g mins\n", (time2-time1)/60 );
return 0;
}
