#include <math.h>
#include <stdio.h>
#include "main.h"
#include <boost/numeric/odeint.hpp>
#include <boost/program_options.hpp>
#include <omp.h>
#include <string>

using namespace boost::numeric::odeint;

//typedef symplectic_rkn_sb3a_mclachlan< state_type > stepper_type;
//typedef runge_kutta4< state_type > stepper_type;
//typedef runge_kutta_cash_karp54< state_type > stepper_type;
//typedef runge_kutta_dopri5< state_type > stepper_type;
typedef runge_kutta_fehlberg78< state_type > stepper_type;
//typedef velocity_verlet< state_type > stepper_type;

int main( int argc, char* argv[]){
    char *filename;
    struct Indata infile; //to store coefficients
    //state_type acc(6), sph(6);
    const double endTime=1e18;
    //double dt=std::stod(argv[1]);
    int nPoints,Npart;//=std::stoi(argv[1]);
    std::vector<state_type> xinit;
    //std::vector<state_type> X;
    state_type T;

    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()
        ("infile1", po::value<std::string>()->required(),"input file 1")
        ("infile2", po::value<std::string>(),"input file 2")
        ("o,o", po::value<std::string>()->default_value("save.hdf5"),"output file")
        ("n,n", po::value<int>(& nPoints)->default_value(32768),"no. of integration points")
        ("i,i", po::value<std::string>()->default_value("init.hdf5"),"File with initial positions")
        ("N,N", po::value<int>(), "no. of particles to integrate");
    po::positional_options_description positionalOptions;
    positionalOptions.add("infile1",1);
    positionalOptions.add("infile2",1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc,argv).options(desc).positional(positionalOptions).run(),vm);
    po::notify(vm);
    if ( vm.count("infile2") ){
        std::cout << "2 coefficient files given, 2 component model assumed\n";
    }
    std::cout << "Output file is: " << vm["o"].as<std::string>() << "\n";
    std::cout << "Npoints=" << nPoints << "\t";

    double dt=endTime/nPoints;
    std::cout << "dt=" << dt << '\t' << "endTime=" <<endTime <<'\n';

    //snprintf(filename,150,vm["infile"].as<std::string>());
    filename=strdup(vm["infile1"].as<std::string>().c_str());
    loadHdf5Input(filename, &infile);
    calcAcc ACC(&infile);
    /*
    if (vm.count("infile2")){
        //snprintf(filename,150,vm["infile2"].as<std::string>());
        struct Indata infile2;
        loadHdf5Input(strdup(vm["infile2"].as<std::string>().c_str()), &infile2);
        calcAcc ACC(&infile,&infile2);
    }
    else{
        calcAcc ACC(&infile);
    }*/

    filename=strdup(vm["i"].as<std::string>().c_str());
    loadHdf5Init(filename, xinit);
    if (vm.count("N")){
        Npart=vm["N"].as<int>();
    }
    else {
        Npart=xinit.size();
    }
    state_type XX(Npart*(nPoints+1)*6);

    std::cout << "Num points to integrate = " << Npart << '\n';
    fprintf(stdout,"G=%e\n",GRAVITY);
    fprintf(stdout,"a=%e\n",infile.scalerad);
    fprintf(stdout,"X0=%e,%e,%e,%e,%e,%e\n",xinit[0][0],xinit[0][1],xinit[0][2],xinit[0][3],xinit[0][4],xinit[0][5]);


    #pragma omp parallel
    {
        if (omp_get_thread_num()==0){
            std::cout << "Num_threads = " << omp_get_num_threads() << "\n";
        }
        //std::vector<state_type> Xpriv;
        //state_type Tpriv;
        #pragma omp for
        for (int n=0; n<Npart; n++)
        {
            std::vector<state_type> Xpriv;
            state_type Tpriv;
            nPoints=integrate_const( stepper_type() , ACC , xinit[n] ,
                0.0 , endTime, dt , saveStates(Xpriv,Tpriv) );
            
            #pragma omp critical
            {
            for (int x=0;x<=nPoints; x++)
            {
                for (int y=0; y<6; y++)
                {
                    XX[n*6*(nPoints+1)+6*x+y]=Xpriv[x][y];
                }

            }
            }//end omp critical
            //
            if (n==0) T=Tpriv;

        }//end omp parallel

        /*
        #pragma omp critical
        X.insert(X.end(), Xpriv.begin(), Xpriv.end());
        T.insert(T.end(),Tpriv.begin(),Tpriv.end());
        */
    }

    /* Print info to screen for checking */
    for (int n=0; n<Npart; n++)
    {
        std::cout << n << '\n';
        std::cout << T[0] << '\t' ;
        for ( size_t j=0; j<=5; j++)
        {
            std::cout <<  XX[n*6*(nPoints+1)+j] << '\t';
        }
        std::cout << '\n';

        std::cout <<  T[nPoints] << '\t' ;
        for ( size_t j=0; j<=5; j++)
        {
            std::cout <<  XX[n*6*(nPoints+1)+6*(nPoints)+j] << '\t';
        }
        std::cout << '\n';
    }

    //snprintf(filename,150,"./save.hdf5");
    filename=strdup(vm["o"].as<std::string>().c_str());
    saveHdf5(filename,XX,T);
    fprintf(stdout,"Done!!\n");

    free(filename);
return 0;
}

