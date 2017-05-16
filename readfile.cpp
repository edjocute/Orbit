#include "main.h"
#include <math.h>
#include <hdf5.h>
#include <stdio.h>
#include <boost/program_options.hpp>
#include <string>
#include <fstream>


/* Read and load coefficients Knlm and scale radius a */
int loadHdf5Input(char *filename, struct Indata *var){
    hid_t   hdf_file,hdf_group,hdf_data;
    herr_t  status;
    double  temp[NMAX][LMAX][LMAX][2];

    fprintf(stdout,"Reading file %s ...",filename);
    hdf_file = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
    if (hdf_file < 0){
        return -1;
    }

    
    if ( (hdf_group=H5Gopen2(hdf_file,"/",H5P_DEFAULT)) < 0){
        H5Gclose(hdf_file);
        return -1;
    }

    if ( (hdf_data=H5Dopen2(hdf_file,"/Knlm",H5P_DEFAULT)) < 0){
        H5Dclose(hdf_data);
        return -1;
    }
    //status=H5Dread(hdf_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, var->Knlm);
    status=H5Dread(hdf_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
    for (int n=0; n<NMAX; n++){
        for (int l=0; l<LMAX; l++){
            for (int m=0; m<LMAX; m++){
                var->Knlm[l][m][n][0]=temp[n][l][m][0]; // I reorder the matrix since the summation
                var->Knlm[l][m][n][1]=temp[n][l][m][1]; // over n is done first
                                                        // [0]/[1] -> cosine/sine terms
            }
        }
    }
    
    /* Read virial radius (not really needed) */
    if ( (hdf_data=H5Dopen2(hdf_file,"/Rvir",H5P_DEFAULT)) < 0){
        H5Dclose(hdf_data);
        return -1;
    }
    status=H5Dread(hdf_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var->virialrad);

    /* Read Hernquist scale radius needed to normalize positions */
    if ( (hdf_data=H5Dopen2(hdf_file,"/a",H5P_DEFAULT)) < 0){
        H5Dclose(hdf_data);
        return -1;
    }
    status=H5Dread(hdf_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var->scalerad);

    //H5Gclose(hdf_group);
    H5Fclose(hdf_file);
    H5Dclose(hdf_data);
    fprintf(stdout," file read successfully!\n");
    return 0;
}

/* Load and read file with initial positions and velocites */
int loadHdf5Init(char *filename, std::vector<state_type> &init){
    hid_t   hdf_file,hdf_group,hdf_data,hdf_dspace;
    herr_t  status;
    state_type x(6);

    fprintf(stdout,"Reading file %s ...\n",filename);
    hdf_file = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
    if (hdf_file < 0){
        return -1;
    }

    if ( (hdf_group=H5Gopen2(hdf_file,"/",H5P_DEFAULT)) < 0){
        H5Fclose(hdf_file);
        return -1;
    }

    if ( (hdf_data=H5Dopen2(hdf_file,"/x",H5P_DEFAULT)) < 0){
        H5Dclose(hdf_data);
        return -1;
    }
    hdf_dspace = H5Dget_space(hdf_data);
    hsize_t dims[2];
    H5Sget_simple_extent_dims(hdf_dspace,dims,NULL);
    fprintf(stdout,"No. of test particles read = %d\n",dims[0]);
    double *temp = new double[dims[0]*dims[1]];
    status=H5Dread(hdf_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
    for (int i=0; i<dims[0]; i++){
        for (int j=0; j<dims[1]; j++){
            x[j]=temp[i*dims[1]+j];
        }
        init.push_back(x);
    }

    delete [] temp;
    H5Gclose(hdf_group);
    H5Fclose(hdf_file);
    H5Dclose(hdf_data);
    fprintf(stdout," file read successfully!\n");
    return 0;
}


int saveHdf5(Readparams &allparams, state_type &OUT, state_type &TIME){
    hid_t   hdf_file,hdf_group,hdf_data,dataspace_id;
    herr_t  status;

    fprintf(stdout,"Writing file %s ...",allparams.outfilename);
    hdf_file = H5Fcreate(allparams.outfilename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    if (hdf_file < 0){
        return -1;
    }

    /*
    if ( (hdf_group=H5Gopen2(hdf_file,"/",H5P_DEFAULT)) < 0){
        H5Fclose(hdf_file);
        return -1;
    }*/
    
    /* Write particle positions and velocities.
     * Ordered in chunk where first Nstep+1 lines correspond to particle 1,
     * step Nstep+1 chunk correspond to particle 2 etc. */
    hsize_t dims[1]={OUT.size()};
    dataspace_id=H5Screate_simple(1,dims,NULL);
    if ( (hdf_data=H5Dcreate2(hdf_file,"x",H5T_NATIVE_DOUBLE,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT)) < 0){
        H5Dclose(hdf_data);
        return -1;
    }
    status=H5Dwrite(hdf_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &OUT[0]);

    /*Write times*/
    hsize_t dims2[1]={TIME.size()};
    dataspace_id=H5Screate_simple(1,dims2,NULL);
    if ( (hdf_data=H5Dcreate2(hdf_file,"t",H5T_NATIVE_DOUBLE,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT)) < 0){
        H5Dclose(hdf_data);
        return -1;
    }
    status=H5Dwrite(hdf_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &TIME[0]);

    /*Write no. of components*/
    dims2[0]={1};
    dataspace_id=H5Screate_simple(1,dims2,NULL);
    if ( (hdf_data=H5Dcreate2(hdf_file,"NumComponents",H5T_NATIVE_INT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT)) < 0){
        H5Dclose(hdf_data);
        return -1;
    }
    status=H5Dwrite(hdf_data, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &allparams.NumComponents);


    H5Fclose(hdf_file);
    H5Dclose(hdf_data);
    H5Sclose(dataspace_id);
    fprintf(stdout," file written successfully!\n");
    return 0;
}

int Readparams::read(int argc, char *argv[]){
    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()
        ("infile1", po::value<std::string>()->required(),"input file 1")
        ("infile2", po::value<std::string>(),"input file 2")
        ("o,o", po::value<std::string>()->default_value("save.hdf5"),"output file")
        ("n,n", po::value<int>(&nPoints)->default_value(16384),"no. of integration points")
        ("i,i", po::value<std::string>()->default_value("init.hdf5"),"File with initial positions")
        ("N,N", po::value<int>(&Npart), "no. of particles to integrate")
        ("e,e", po::value<double>( &endTime)->default_value(5e18), "End time of integration in s")
        ("Firstpass,F", po::value<bool> () -> default_value(true), "First or second pass")
        //("dt,dt", po::value<double>(& dt) -> default_value(3e13), "Time Step")
        ("verbose,v", po::value<bool> (& verbose) -> default_value(false), "Verbosity")
        ("tol,", po::value<double>(& tol) -> default_value(1e-12), "Relative tolerance for              integration");
    po::positional_options_description positionalOptions;
    positionalOptions.add("infile1",1);
    positionalOptions.add("infile2",1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc,argv).options(desc).positional(positionalOptions).run(),vm);
    po::notify(vm);

    infilename1=strdup(vm["infile1"].as<std::string>().c_str());
    initfilename=strdup(vm["i"].as<std::string>().c_str());
    outfilename=strdup(vm["o"].as<std::string>().c_str());
    saveint = nPoints / NumSavepoints;

    /* Print some info */
    std::cout << "Output file is: " << vm["o"].as<std::string>() << "\n";
    std::cout << "No. of time points/ save points = " << nPoints <<","<< NumSavepoints <<"\n";
    std::cout << "Save Interval = " << saveint << "\n";
    std::cout << "Rtol = " << tol << "\n";
    std::cout << "NLIM,LLIM,NDIM=" << NLIM <<","<< LLIM << "," << NDIM << "\n";
    std::cout << "G="<< GRAVITY << "\n" << "\n";


    if (vm.count("infile2")){
        std::cout << "2 coefficient files given, 2 component model assumed\n";
        NumComponents=2;
        twocomp=true;
        infilename2=strdup(vm["infile2"].as<std::string>().c_str());
    }

    firstpass=vm["Firstpass"].as<bool>();

    std::ofstream compfile("components.param");
    if (compfile.is_open()){
        if (vm.count("infile2"))
            compfile << "2";
        else
            compfile << "1";
        compfile.close();
    }

    return 0;
}
