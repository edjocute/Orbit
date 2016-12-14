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
    char *filename;
    struct Indata infile,infile2; //to store coefficients
    double endTime, dt, tol;
    int nPoints, Npart, saveint;
    int savenpoints=16384; //no. of time points to be saved in output file per particle
    std::vector<state_type> xinit,xinit_run; //to store array of initial positions and velocities
    bool verbose;

    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()
        ("infile1", po::value<std::string>()->required(),"input file 1")
        ("infile2", po::value<std::string>(),"input file 2")
        ("o,o", po::value<std::string>()->default_value("save.hdf5"),"output file")
        ("n,n", po::value<int>(& nPoints)->default_value(16384),"no. of integration points")
        ("i,i", po::value<std::string>()->default_value("init.hdf5"),"File with initial positions")
        ("N,N", po::value<int>(), "no. of particles to integrate")
        ("e,e", po::value<double>( &endTime)->default_value(5e18), "End time of integration in s")
        ("Firstpass,F", po::value<bool> () -> default_value(true), "First or second pass")
        ("dt,dt", po::value<double>(& dt) -> default_value(3e13), "Time Step")
        ("verbose,v", po::value<bool> (& verbose) -> default_value(false), "Verbosity")
        ("tol,", po::value<double>(& tol) -> default_value(1e-12), "Relative tolerance for integration");
    po::positional_options_description positionalOptions;
    positionalOptions.add("infile1",1);
    positionalOptions.add("infile2",1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc,argv).options(desc).positional(positionalOptions).run(),vm);
    po::notify(vm);

    std::cout << "Output file is: " << vm["o"].as<std::string>() << "\n";
    std::cout << "No. of time points/ save points = " << nPoints <<","<< savenpoints <<"\n";
    std::cout << "Rtol = " << tol << "\n";
    saveint = nPoints / savenpoints;
    std::cout << "Saveint = " << saveint << "\n" << "\n";


    filename=strdup(vm["infile1"].as<std::string>().c_str());
    loadHdf5Input(filename, &infile);
    
    /* Call either the 1-component or 2-component routine 
     * based on whether 2nd component input is present
     */
    if (vm.count("infile2")){
        std::cout << "2 coefficient files given, 2 component model assumed\n";
        filename=strdup(vm["infile2"].as<std::string>().c_str());
        loadHdf5Input(strdup(vm["infile2"].as<std::string>().c_str()), &infile2);
    }
    calcAcc ACC = (vm.count("infile2")) ? calcAcc(&infile,&infile2) : calcAcc(&infile);

    filename=strdup(vm["i"].as<std::string>().c_str());
    loadHdf5Init(filename, xinit);
    Npart = (vm.count("N")) ? vm["N"].as<int>() : xinit.size();

    /* Print some info */
    std::cout << "G="<< GRAVITY << "\n" << "\n";
    //fprintf(stdout,"a=%e\n",infile.scalerad);
    //fprintf(stdout,"X0=%e,%e,%e,%e,%e,%e\n",xinit[0][0],xinit[0][1],xinit[0][2],xinit[0][3],xinit[0][4],xinit[0][5]);

    /* Check whether 1st or 2nd pass and adjust the integration points accordingly. */
    if (vm["Firstpass"].as<bool>()){
        Npart=100;
        xinit_run = std::vector<state_type>();
        std::cout << "FIRST PASS. Points to integrate = " << Npart << "\n";
        std::cout << "T_end=" << endTime <<'\n';
        int Interval= xinit.size()/Npart;
        //std::cout << "interval = " << Interval << "\n";
        for (int n=0; n<Npart; n++){
            xinit_run.push_back(xinit[n*Interval]);
        }
    }
    else{
        std::cout << "SECOND PASS. Integrating all points." << "\n";
        xinit_run = xinit;
    }

    /* Create output arrays based on number of integration points */
    state_type XX(Npart*(savenpoints+1)*6); //output array
    state_type T(Npart*(savenpoints+1));
    std::cout << "Total points in input/Points to integrate = " << xinit.size() << "/" << Npart << '\n';

    /* Calculate end time array for particles  */
    state_type endTimeAll(Npart), dtAll(Npart);
    if (vm["Firstpass"].as<bool>()){
        for (int n=0; n<Npart; n++){
            endTimeAll[n]=endTime;
            dtAll[n]=endTime/nPoints;
        }
    }
    else{
        /* Read fitting params from linear regression from file */
        std::ifstream paramfile("params.txt");
        double m,c,E,endtemp;
        //int ntemp;
        paramfile >> m >> c;
        std::cout << "Fitting params m,c=" << m << "," << c << "\n";
        for (int n=0; n<Npart; n++){
            ACC.getEnergy(xinit_run[n],E);
            endTimeAll[n]=70/pow(10,fabs(E)*m/10000+c) * endTime;
            dtAll[n]=endTimeAll[n]/nPoints;
            /*ntemp=endTimeAll[n]/dt;
            if ( double(ntemp % savenpoints)/ntemp > 0.1){
                ntemp--;                ntemp |= ntemp >> 1;
                ntemp |= ntemp >> 2;    ntemp |= ntemp >> 4; 
                ntemp |= ntemp >> 8;    ntemp |= ntemp >> 16;
                ntemp++;
                dtAll[n]=endTimeAll[n]/ntemp;
            }
            else{
                ntemp = (ntemp/savenpoints)*savenpoints;
                dtAll[n]=dt;
                endTimeAll[n]=ntemp*dt;
            }
            std::cout << "ntemp/dt=" << ntemp << "," << dtAll[n] << "\n";*/
        }
    }


    //dense_output_dopri5_type dopri5 = make_dense_output( 1E-10 , tol , dopri5_type() );
    double time1=omp_get_wtime();
    #pragma omp parallel
    {
        if (omp_get_thread_num()==0){
            std::cout << "Num_threads = " << omp_get_num_threads() << "\n";
        }
        #pragma omp for schedule(dynamic,2)
        for (int n=0; n<Npart; n++)
        {
            //std::cout << "endtime=" << n << ","  << endTimeAll[n] << "\n";
            std::vector<state_type> Xpriv;
            state_type Tpriv;
            integrate_const( make_controlled(1.0e-10,tol,rkf78_type())
                                     //dopri5
                    ,ACC , xinit_run[n] , 0.0 , endTimeAll[n], dtAll[n] , saveStates(Xpriv,Tpriv) );
            if (verbose){
                std::cout << "n,npoints,size(T) = " << n <<","<< nPoints <<","<< Tpriv.size() <<"\n";
            }

            //saveint=nPoints/savenpoints;
            //std::cout << "nPoints=" << nPoints << "\n";
            for (int x=0;x<=savenpoints; x++)
            {
                T[n*(savenpoints+1)+x] = Tpriv[saveint*x];
                for (int y=0; y<6; y++)
                {
                    XX[n*6*(savenpoints+1)+6*x+y]=Xpriv[saveint*x][y];
                }
            }
            //if (n==0) T=Tpriv;
        }//end omp for
    }//end omp parallel
    double time2=omp_get_wtime();


    /* Print some info to screen for checking */
    if (verbose){
    for (int n=0; n< ((Npart<10) ? Npart : 10); n++)
    {
        std::cout << n << '\n';
        std::cout << T[n*(savenpoints+1)] << '\t' ;
        for ( size_t j=0; j<=5; j++)
        {
            std::cout <<  XX[n*6*(savenpoints+1)+j] << '\t';
        }
        std::cout << '\n';

        std::cout <<  T[(n+1)*(savenpoints+1)-1] << '\t' ;
        for ( size_t j=0; j<=5; j++)
        {
            std::cout <<  XX[n*6*(savenpoints+1)+6*(savenpoints)+j] << '\t';
        }
        std::cout << '\n';
    }
    }

    /* Save data to output file */
    filename=strdup(vm["o"].as<std::string>().c_str());
    saveHdf5(filename,XX,T);
    fprintf(stdout,"Done! This took %.4g mins\n", (time2-time1)/60 );

    free(filename);
return 0;
}
