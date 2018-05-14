#include "main.h"
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <fstream>
#include <complex>
#include <fftw3.h>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;
//typedef symplectic_rkn_sb3a_mclachlan< state_type > stepper_type;
//typedef runge_kutta4< state_type > stepper_type;
//typedef runge_kutta_cash_karp54< state_type > ck54_type;
//typedef runge_kutta_dopri5< state_type > dopri5_type;
//typedef controlled_runge_kutta< dopri5_type > controlled_dopri5_type;
//typedef dense_output_runge_kutta< controlled_dopri5_type > dense_output_dopri5_type;
typedef runge_kutta_fehlberg78<array6> rkf78_type;
//typedef velocity_verlet< state_type > stepper_type;
//Readparams allparams;

void calcOrb::test(std::vector<array6> const &xinit){
        double pot=0;
        //calcAcc ACC=(allparams.twocomp) ? calcAcc(&invar,&invarfp) : calcAcc(&invar);

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
                /*ACC.getCartAcc(xinit[n],testacc);
                for (int m=0; m<6; m++){ 
                    accfile << testacc[m] << "\t";
                }*/
                accfile << "\n";
            }
        }
        accfile.close();
        std::cout << "done!" << std::endl;
}

void calcOrb::integrate(state_type &XX, state_type &T, std::vector<array6> &xinit_run){
    /* Calculate end time array for particles  */
    std::cout << allparams.Npart << "\n";
    state_type endTimeAll(allparams.Npart), dtAll(allparams.Npart);

    //calcAcc ACC=(allparams.twocomp) ? calcAcc(&invar,&invarfp) : calcAcc(&invar);
    if (allparams.firstpass){
        #pragma omp parallel for
        for (int n=0; n<allparams.Npart; n++){
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
            std::cout << "Starting loops... for " << allparams.Npart << "particles\n";
        }

        std::vector<array6> Xpriv;
        state_type Tpriv, v(allparams.nPoints);
        Xpriv.reserve(16385);
        Tpriv.reserve(16385);

        complex_type fv(allparams.nPoints/2+1);
        int maxi;
        double endtimepriv,dtpriv,normv,maxamp;
        fftw_plan plan;

        //calcAcc ACC = (allparams.twocomp) ? calcAcc(&infile,&infile2) : calcAcc(&infile);
        /*std::vector<std::complex<double>> v(allparams.nPoints);*/

        /*FFTW plan is NOT thread safe. Do fftw_plan one thread at a time.*/
        #pragma omp critical(createplan)
        plan = fftw_plan_dft_r2c_1d(allparams.nPoints,v.data(),
                reinterpret_cast<fftw_complex*>(fv.data()),FFTW_ESTIMATE);
        /*fftw_plan plan = fftw_plan_dft_1d(allparams.nPoints,reinterpret_cast<fftw_complex*>(&v[0]),
            reinterpret_cast<fftw_complex*>(&v[0]),FFTW_FORWARD,FFTW_ESTIMATE);*/

        if (omp_get_thread_num()==0){
            std::cout << "Num_threads = " << omp_get_num_threads() << "\n";
        }

        #pragma omp for schedule(dynamic,2)
        for (int n=0; n<allparams.Npart; n++){

            /*Coarse integration*/
            Xpriv.clear(); Tpriv.clear(); //IMPORTANT: DO NOT DELETE THIS LINE
            endtimepriv = allparams.endTime;
            dtpriv = endtimepriv/(allparams.nPoints-1);
            integrate_const( make_controlled(1.0e-10,1e-7,rkf78_type())
                    ,ACC , xinit_run[n] , 0.0 , endtimepriv, dtpriv, saveStates(Xpriv,Tpriv) );
            getperiod(plan, Xpriv, v, fv, maxamp, maxi);
            //std::cout << n << "," << maxi << "\n";

            /*Fine integration*/
            if (!allparams.firstpass){
                Xpriv.clear(); Tpriv.clear();
                endtimepriv = 70./maxi * endtimepriv;
                dtpriv = endtimepriv/(allparams.nPoints-1);
                integrate_const( make_controlled(1.0e-10,allparams.tol,rkf78_type())
                            ,ACC , xinit_run[n] , 0.0 , endtimepriv, dtpriv , saveStates(Xpriv,Tpriv) );
                getperiod(plan, Xpriv, v, fv, maxamp, maxi);

                /*Final integration if needed*/
                if ( (maxi<40) || (maxi>100) ){
                    Xpriv.clear(); Tpriv.clear();
                    std::cout << n << "," << maxi << "," << endtimepriv << "," << dtpriv << "\n";
                    endtimepriv *= 70./maxi;
                    dtpriv = endtimepriv/(allparams.nPoints-1);
                    integrate_const( make_controlled(1.0e-10,0.1*allparams.tol,rkf78_type())
                        , ACC , xinit_run[n] , 0.0 , endtimepriv, dtpriv , saveStates(Xpriv,Tpriv) );
                    getperiod(plan,Xpriv,v,fv,maxamp,maxi);
                    std::cout << n << "," << maxi << "," << endtimepriv << "," << dtpriv << "\n";
                }
            }

            /*Save private variables to global*/
            if (allparams.verbose){
                std::cout << "n,npoints,size(T) = " << n <<","<< allparams.nPoints <<","<< Tpriv.size()  <<"\n";
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
}

void calcOrb::getperiod(fftw_plan plan, const std::vector<array6>& Xpriv, state_type &v, 
        complex_type &fv, double &maxamp, int &maxi){
    maxamp=0;
    double normv;
    for (int k=0; k<3; k++){
        #pragma omp simd
        for (unsigned int i=0; i<allparams.nPoints; i++)
            v[i]=Xpriv[allparams.saveint*i][k];
            //v[i]=std::complex<double> (Xpriv[allparams.saveint*i][k],Xpriv[allparams-         >saveint*  i][k+3]);
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
}

void calcOrb::integrateDecay(state_type &XX, state_type &T, std::vector<array6> &xinit_run){
    //calcAcc ACC=(allparams.twocomp) ? calcAcc(&invar,&invarfp) : calcAcc(&invar);
    std::cout << allparams.Npart << "\n";
    double time1=omp_get_wtime();
    #pragma omp parallel shared(XX,T,allparams)
    {
        if (omp_get_thread_num()==0){
            std::cout << "Starting loops... for " << allparams.Npart << "particles\n";
        }

        std::vector<array6> Xpriv;
        state_type Tpriv;
        Xpriv.reserve(16385);
        Tpriv.reserve(16385);

        double endtimepriv,dtpriv,normv,maxamp;

        if (omp_get_thread_num()==0){
            std::cout << "Num_threads = " << omp_get_num_threads() << "\n";
        }

        #pragma omp for schedule(dynamic,2)
        for (int n=0; n<allparams.Npart; n++){
            Xpriv.clear(); Tpriv.clear(); //IMPORTANT: DO NOT DELETE THIS LINE
            endtimepriv = allparams.decaytime;
            dtpriv = endtimepriv/(allparams.nPoints-1);
            integrate_const( make_controlled(1.0e-10,1e-12,rkf78_type())
                    , ACC , xinit_run[n] , 0.0 , allparams.endTime, dtpriv, saveStates(Xpriv,Tpriv) );

            for (int i = 0; i<6; i++){
                xinit_run[n][i] = Xpriv.back()[i];
            }
        }//end omp for
    }//end omp parallel
}
