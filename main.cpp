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
//typedef runge_kutta_cash_karp54< state_type > stepper_type;
//typedef runge_kutta_dopri5< state_type > stepper_type;
typedef runge_kutta_fehlberg78< state_type > stepper_type;
//typedef velocity_verlet< state_type > stepper_type;

int main( int argc, char* argv[]){
    char *filename;
    struct Indata infile,infile2; //to store coefficients
    double endTime, dt;
    int nPoints, Npart, saveint, savenpoints=16384;
    std::vector<state_type> xinit,xinit_run; //to store array of initial positions and velocities

    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()
        ("infile1", po::value<std::string>()->required(),"input file 1")
        ("infile2", po::value<std::string>(),"input file 2")
        ("o,o", po::value<std::string>()->default_value("save.hdf5"),"output file")
        ("n,n", po::value<int>(& nPoints)->default_value(32768),"no. of integration points")
        ("i,i", po::value<std::string>()->default_value("init.hdf5"),"File with initial positions")
        ("N,N", po::value<int>(), "no. of particles to integrate")
        ("e,e", po::value<double>( &endTime)->default_value(5e18), "End time of integration in s")
        ("Firstpass,F", po::value<bool> () -> default_value(true), "First or second pass")
        ("dt,dt", po::value<double>(& dt) -> default_value(3e13), "Time Step");
    po::positional_options_description positionalOptions;
    positionalOptions.add("infile1",1);
    positionalOptions.add("infile2",1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc,argv).options(desc).positional(positionalOptions).run(),vm);
    po::notify(vm);

    std::cout << "Output file is: " << vm["o"].as<std::string>() << "\n";
    std::cout << "Npoints=" << nPoints << "\t";
    
    //dt=endTime/nPoints;
    std::cout << "dt=" << dt << '\t' << "endTime=" <<endTime <<'\n';

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
    //int saveint=nPoints/16384;
    calcAcc ACC = (vm.count("infile2")) ? calcAcc(&infile,&infile2) : calcAcc(&infile);

    filename=strdup(vm["i"].as<std::string>().c_str());
    loadHdf5Init(filename, xinit);
    Npart = (vm.count("N")) ? vm["N"].as<int>() : xinit.size();

    /* Print some info */
    fprintf(stdout,"G=%e\n",GRAVITY);
    fprintf(stdout,"a=%e\n",infile.scalerad);
    fprintf(stdout,"X0=%e,%e,%e,%e,%e,%e\n",xinit[0][0],xinit[0][1],xinit[0][2],xinit[0][3],xinit[0][4],xinit[0][5]);

    /* Check whether 1st or 2nd pass and adjust the integration points accordingly. */
    if (vm["Firstpass"].as<bool>()){
        Npart=100;
        xinit_run = std::vector<state_type>();
        std::cout << "Doing FIRST PASS. Num points to integrate = " << Npart << "\n";
        int Interval= xinit.size()/Npart;
        std::cout << "interval = " << Interval << "\n";
        for (int n=0; n<Npart; n++){
            xinit_run.push_back(xinit[n*Interval]);
        }
    }
    else{
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
        }
    }
    else{
        /* Read fitting params from linear regression from file */
        std::ifstream paramfile("params.txt");
        double m,c,E,endtemp;
        int ntemp;
        paramfile >> m >> c;
        std::cout << "m,c=" << m << "," << c << "\n";
        for (int n=0; n<Npart; n++){
            ACC.getEnergy(xinit_run[n],E);
            endTimeAll[n]=70/pow(10,fabs(E)*m/10000+c) * endTime;
            ntemp=endTimeAll[n]/dt;
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
            std::cout << "ntemp/dt=" << ntemp << "," << dtAll[n] << "\n";
        }
    }


    double time1=omp_get_wtime();
    #pragma omp parallel
    {
        if (omp_get_thread_num()==0){
            std::cout << "Num_threads = " << omp_get_num_threads() << "\n";
        }
        //std::vector<state_type> Xpriv;
        //state_type Tpriv;
        #pragma omp for schedule(static,2) private(saveint,nPoints)
        for (int n=0; n<Npart; n++)
        {
            //std::cout << "endtime=" << n << ","  << endTimeAll[n] << "\n";
            std::vector<state_type> Xpriv;
            state_type Tpriv;
            nPoints=integrate_const( stepper_type() , ACC , xinit_run[n] ,
                0.0 , endTimeAll[n], dtAll[n] , saveStates(Xpriv,Tpriv) );
            

            saveint=nPoints/savenpoints;
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
    for (int n=0; n<10; n++)
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

    /* Save data to output file */
    filename=strdup(vm["o"].as<std::string>().c_str());
    saveHdf5(filename,XX,T);
    fprintf(stdout,"Done! This took %.4g mins\n", (time2-time1)/60 );

    free(filename);
return 0;
}
