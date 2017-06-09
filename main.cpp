#include <math.h>
#include <stdio.h>
#include "main.h"
#include <boost/numeric/odeint.hpp>
#include <boost/program_options.hpp>
#include <omp.h>
#include <string>
#include <fstream>
#include <complex>
#include <fftw3.h>

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
    //calcAcc ACC = (allparams.twocomp) ? calcAcc(&infile,&infile2) : calcAcc(&infile);

    loadHdf5Init(allparams.initfilename, xinit);
    //Npart = (vm.count("N")) ? vm["N"].as<int>() : xinit.size();
    /* Check whether 1st or 2nd pass and adjust the integration points accordingly. */


    /* FOR TESTING ONLY:
     * Save test potentials to file and exit */
    if (allparams.test){
        calcAcc ACC = (allparams.twocomp) ? calcAcc(&infile,&infile2) : calcAcc(&infile);
        double pot=0;
        std::cout << "Saving potentials to file:" << "\n";
        std::ofstream potfile("Potentials.txt");
        if (potfile.is_open()){
            for (int n=0; n<xinit.size(); n++){
                ACC.getPotential(xinit[n],pot);
                potfile << pot << std::endl;
            }
        }
        potfile.close();

        state_type testacc(6);
        std::cout << "Saving accelerations to file:" << std::endl;
        std::ofstream accfile("Acc.txt");
        if (accfile.is_open()){
            for (int n=0; n<xinit.size(); n++){
                ACC.getCartAcc(xinit[n],testacc);
                for (int m=0; m<6; m++){ 
                    accfile << testacc[m] << "\t";
                }
                accfile << "\n";
            }
        }
        std::cout << "done!" << std::endl;
        return 0;
    }


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
    state_type XX(Npart*(allparams.NumSavepoints)*6); //output array
    state_type T(Npart*(allparams.NumSavepoints));
    std::cout << "Total points in input/Points to integrate = " << xinit.size() << "/" << Npart << '\n';

    /* Calculate end time array for particles  */
    state_type endTimeAll(Npart), dtAll(Npart);
    if (allparams.firstpass){
        #pragma omp parallel for
        for (int n=0; n<Npart; n++){
            endTimeAll[n]=allparams.endTime;
            dtAll[n]=allparams.endTime/(allparams.nPoints-1);
        }
    }
    //else{
        /* Read fitting params from linear regression from file */
        //std::ifstream paramfile("params.txt");
        //double m,c,E,endtemp;
        /*if (paramfile.is_open()){
            paramfile >> m >> c;
            paramfile.close();
        }
        else std::cout << "Unable to open params.txt";
        std::cout << "Fitting params m,c=" << m << "," << c << "\n";*/
/*
        for (int n=0; n<Npart; n++){
            ACC.getEnergy(xinit_run[n],E);
            endTimeAll[n]=70./pow(10,fabs(E)*m/10000+c) * allparams.endTime;
            dtAll[n]=endTimeAll[n]/(allparams.nPoints-1);*/
            /*if ( double(ntemp % NumSavepoints)/ntemp > 0.1){
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
        //}
    //}

    double time1=omp_get_wtime();
    #pragma omp parallel shared(XX,T,allparams)
    {
        if (omp_get_thread_num()==0){
            std::cout << "Starting loops..." << "\n";
        }

        std::vector<state_type> Xpriv;
        state_type Tpriv, v(allparams.nPoints);
        std::vector<std::complex<double>> fv(allparams.nPoints/2+1);
        int maxi;
        double endtimepriv,dtpriv,normv,maxamp;
        fftw_plan plan;

        calcAcc ACC = (allparams.twocomp) ? calcAcc(&infile,&infile2) : calcAcc(&infile);

        /*std::vector<std::complex<double>> v(allparams.nPoints);*/

        /*FFTW plan is NOT thread safe. Do fftw_plan one thread at a time.*/
        #pragma omp critical(createplan)
        plan = fftw_plan_dft_r2c_1d(allparams.nPoints,v.data(),reinterpret_cast<fftw_complex*>(fv.data()),FFTW_ESTIMATE);
        
            /*fftw_plan plan = fftw_plan_dft_1d(allparams.nPoints,reinterpret_cast<fftw_complex*>(&v[0]),
                                reinterpret_cast<fftw_complex*>(&v[0]),FFTW_FORWARD,FFTW_ESTIMATE);*/
        if (omp_get_thread_num()==0){
            std::cout << "Num_threads = " << omp_get_num_threads() << "\n";
        }

        #pragma omp for schedule(dynamic,1)
        for (int n=0; n<Npart; n++){

            /*Coarse integration*/
            Xpriv.clear(); Tpriv.clear();
            dtpriv = allparams.endTime/(allparams.nPoints-1);
            integrate_const( make_controlled(1.0e-10,1e-7,rkf78_type()) //dopri5
                    ,ACC , xinit_run[n] , 0.0 , allparams.endTime, dtpriv, saveStates(Xpriv,Tpriv) );

            maxamp=0; 
            for (int k=0; k<3; k++){
                #pragma omp simd
                for (unsigned int i=0; i<allparams.nPoints; i++)
                    v[i]=Xpriv[allparams.saveint*i][k];
                    //v[i]=std::complex<double> (Xpriv[allparams.saveint*i][k],Xpriv[allparams.saveint*i][k+3]);
                fftw_execute(plan);

                for (int i=1; i<allparams.nPoints/2+1; i++){
                    normv=std::norm(fv[i]);
                    //normv=std::norm(v[i])+std::norm(v[allparams.nPoints-i]);
                    if (normv > maxamp){
                        maxamp=normv;
                        maxi=i;
                    }
                }
            }
            //std::cout << n << "," << maxi << "\n";

            /*Fine integration*/
            if (!allparams.firstpass){
                Xpriv.clear(); Tpriv.clear();
                endtimepriv = 70./maxi * allparams.endTime;
                dtpriv = endtimepriv/(allparams.nPoints-1);

                integrate_const( make_controlled(1.0e-10,allparams.tol,rkf78_type())
                            ,ACC , xinit_run[n] , 0.0 , endtimepriv, dtpriv , saveStates(Xpriv,Tpriv) );

                /*Estimate periods again and redo if periods<40 or periods>100*/
                maxamp=0;
                for (int k=0; k<3; k++){
                    #pragma omp simd
                    for (unsigned int i=0; i<allparams.nPoints; i++)
                        v[i]=Xpriv[allparams.saveint*i][k];
                    fftw_execute(plan);

                    for (int i=1; i<allparams.nPoints/2+1; i++){
                        normv=std::norm(fv[i]);
                        if (normv > maxamp){
                            maxamp=normv;
                            maxi=i;
                        }
                    }
                }

                /*Final integration if needed*/
                if ( (maxi<40) || (maxi>100) ){
                    Xpriv.clear(); Tpriv.clear();
                    endtimepriv = 70./maxi * allparams.endTime;
                    dtpriv = endtimepriv/(allparams.nPoints-1);
                    std::cout << n << "," << maxi << "," << endtimepriv << "," << dtpriv << "\n";
                    integrate_const( make_controlled(1.0e-10,allparams.tol,rkf78_type())
                            ,ACC , xinit_run[n] , 0.0 , endtimepriv, dtpriv , saveStates(Xpriv,Tpriv) );
                }
            }

            /*Save private variables to global*/
            if (allparams.verbose){
                std::cout << "n,npoints,size(T) = " << n <<","<< allparams.nPoints <<","<< Tpriv.size() <<"\n";
            }
            for (int x=0;x<allparams.NumSavepoints; x++){
                T[n*(allparams.NumSavepoints)+x] = Tpriv[allparams.saveint*x];
                for (int y=0; y<6; y++){
                    XX[n*6*(allparams.NumSavepoints)+6*x+y]=Xpriv[allparams.saveint*x][y];
                }
            }
        }//end omp for
        fftw_destroy_plan(plan);
    }//end omp parallel
    double time2=omp_get_wtime();
    std::cout << "... done!" << "\n";

    /* Print some info to screen for checking */
    if (allparams.verbose){
    for (int n=0; n< ((Npart<10) ? Npart : 10); n++){
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
    saveHdf5(allparams,XX,T);
    fprintf(stdout,"All done! This took %.4g mins\n", (time2-time1)/60 );
return 0;
}
