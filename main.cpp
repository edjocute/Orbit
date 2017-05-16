#include <math.h>
#include <stdio.h>
#include "main.h"
#include <boost/numeric/odeint.hpp>
#include <boost/program_options.hpp>
#include <omp.h>
#include <string>
#include <fstream>

using namespace boost::numeric::odeint;

//typedef symplectic_rkn_sb3a_mclachlan< state_type > stepper_type;
//typedef runge_kutta4< state_type > stepper_type;
typedef runge_kutta_cash_karp54< state_type > ck54_type;
typedef runge_kutta_dopri5< state_type > dopri5_type;
typedef controlled_runge_kutta< dopri5_type > controlled_dopri5_type;
typedef dense_output_runge_kutta< controlled_dopri5_type > dense_output_dopri5_type;
typedef runge_kutta_fehlberg78< state_type > rkf78_type;
//typedef velocity_verlet< state_type > stepper_type;

int main( int argc, char* argv[]){
    struct Indata infile,infile2; //to store coefficients
    std::vector<state_type> xinit,xinit_run; //to store array of initial positions and velocities
    int Npart;

    Readparams allparams;
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
    calcAcc ACC = (allparams.twocomp) ? calcAcc(&infile,&infile2) : calcAcc(&infile);

    loadHdf5Init(allparams.initfilename, xinit);
    //Npart = (vm.count("N")) ? vm["N"].as<int>() : xinit.size();
    /* Check whether 1st or 2nd pass and adjust the integration points accordingly. */
    if (allparams.firstpass){
        Npart=100;
        xinit_run = std::vector<state_type>();
        std::cout << "FIRST PASS. Points to integrate = " << Npart << "\n";
        std::cout << "T_end=" << allparams.endTime <<'\n';
        int Interval= xinit.size()/Npart;
        //std::cout << "interval = " << Interval << "\n";
        for (int n=0; n<Npart; n++){
            xinit_run.push_back(xinit[n*Interval]);
        }
    }
    else{
        std::cout << "SECOND PASS. Integrating all points." << "\n";
        Npart=xinit.size();
        xinit_run = xinit;
    }

    /* Create output arrays based on number of integration points */
    state_type XX(Npart*(allparams.NumSavepoints+1)*6); //output array
    state_type T(Npart*(allparams.NumSavepoints+1));
    std::cout << "Total points in input/Points to integrate = " << xinit.size() << "/" << Npart << '\n';

    /* Calculate end time array for particles  */
    state_type endTimeAll(Npart), dtAll(Npart);
    if (allparams.firstpass){
        for (int n=0; n<Npart; n++){
            endTimeAll[n]=allparams.endTime;
            dtAll[n]=allparams.endTime/allparams.nPoints;
        }
    }
    else{
        /* Read fitting params from linear regression from file */
        std::ifstream paramfile("params.txt");
        double m,c,E,endtemp;
        //int ntemp;
        if (paramfile.is_open()){
            paramfile >> m >> c;
            paramfile.close();
        }
        else std::cout << "Unable to open params.txt";

        std::cout << "Fitting params m,c=" << m << "," << c << "\n";
        for (int n=0; n<Npart; n++){
            ACC.getEnergy(xinit_run[n],E);
            endTimeAll[n]=70./pow(10,fabs(E)*m/10000+c) * allparams.endTime;
            dtAll[n]=endTimeAll[n]/allparams.nPoints;
            /*ntemp=endTimeAll[n]/dt;
            if ( double(ntemp % NumSavepoints)/ntemp > 0.1){
                ntemp--;                ntemp |= ntemp >> 1;
                ntemp |= ntemp >> 2;    ntemp |= ntemp >> 4; 
                ntemp |= ntemp >> 8;    ntemp |= ntemp >> 16;
                ntemp++;
                dtAll[n]=endTimeAll[n]/ntemp;
            }
            else{
                ntemp = (ntemp/NumSavepoints)*NumSavepoints;
                dtAll[n]=dt;
                endTimeAll[n]=ntemp*dt;
            }
            std::cout << "ntemp/dt=" << ntemp << "," << dtAll[n] << "\n";*/
        }
    }

    double time1=omp_get_wtime();
    #pragma omp parallel
    {
        if (omp_get_thread_num()==0){
            std::cout << "Num_threads = " << omp_get_num_threads() << "\n";
        }
        #pragma omp for schedule(dynamic,2)
        for (int n=0; n<Npart; n++){
            std::vector<state_type> Xpriv;
            state_type Tpriv;
            integrate_const( make_controlled(1.0e-10,allparams.tol,rkf78_type())
                                     //dopri5
                    ,ACC , xinit_run[n] , 0.0 , endTimeAll[n], dtAll[n] , saveStates(Xpriv,Tpriv) );
            if (allparams.verbose){
                std::cout << "n,npoints,size(T) = " << n <<","<< allparams.nPoints <<","<< Tpriv.size() <<"\n";
            }
            for (int x=0;x<=allparams.NumSavepoints; x++){
                T[n*(allparams.NumSavepoints+1)+x] = Tpriv[allparams.saveint*x];
                for (int y=0; y<6; y++){
                    XX[n*6*(allparams.NumSavepoints+1)+6*x+y]=Xpriv[allparams.saveint*x][y];
                }
            }
        }//end omp for
    }//end omp parallel
    double time2=omp_get_wtime();


    /* Print some info to screen for checking */
    if (allparams.verbose){
    for (int n=0; n< ((Npart<10) ? Npart : 10); n++){
        std::cout << n << '\n';
        std::cout << T[n*(allparams.NumSavepoints+1)] << '\t' ;
        for ( size_t j=0; j<=5; j++){
            std::cout <<  XX[n*6*(allparams.NumSavepoints+1)+j] << '\t';
        }
        std::cout << '\n';

        std::cout <<  T[(n+1)*(allparams.NumSavepoints+1)-1] << '\t' ;
        for ( size_t j=0; j<=5; j++){
            std::cout <<  XX[n*6*(allparams.NumSavepoints+1)+6*(allparams.NumSavepoints)+j] << '\t';
        }
        std::cout << '\n';
    }
    }

    /* Save data to output file */
    saveHdf5(allparams,XX,T);
    fprintf(stdout,"Done! This took %.4g mins\n", (time2-time1)/60 );

    //free(filename);
return 0;
}
