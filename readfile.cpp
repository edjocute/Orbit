#include "main.h"
#include <math.h>
#include <hdf5.h>

int loadHdf5Input(char *filename, struct Indata *var){

    hid_t   hdf_file,hdf_group,hdf_data;
    herr_t  status;
    char    name[150];
    double  temp[NMAX][LMAX][LMAX][2];

    //snprintf(name,150,"%s.hdf5",filename);
    fprintf(stdout,"Reading file %s ...",filename);
    hdf_file = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
    if (hdf_file < 0){
        return -1;
    }

    //snprintf(name,150,"/root");
    if ( (hdf_group=H5Gopen2(hdf_file,"/",H5P_DEFAULT)) < 0){
        H5Fclose(hdf_file);
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
                var->Knlm[l][m][n][0]=temp[n][l][m][0];
                var->Knlm[l][m][n][1]=temp[n][l][m][1];
            }
        }
    }


    
    if ( (hdf_data=H5Dopen2(hdf_file,"/rvir",H5P_DEFAULT)) < 0){
        H5Dclose(hdf_data);
        return -1;
    }
    status=H5Dread(hdf_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var->virialrad);

    if ( (hdf_data=H5Dopen2(hdf_file,"/a",H5P_DEFAULT)) < 0){
        H5Dclose(hdf_data);
        return -1;
    }
    status=H5Dread(hdf_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &var->scalerad);

    H5Gclose(hdf_group);
    H5Fclose(hdf_file);
    H5Dclose(hdf_data);
    fprintf(stdout," file read successfully!\n");
    return 0;
}


int loadHdf5Init(char *filename, std::vector<state_type> &init){

    hid_t   hdf_file,hdf_group,hdf_data,hdf_dspace;
    herr_t  status;
    char    name[150];
    state_type x(6);

    //snprintf(name,150,"%s.hdf5",filename);
    fprintf(stdout,"Reading file %s ...\n",filename);
    hdf_file = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
    if (hdf_file < 0){
        return -1;
    }

    //snprintf(name,150,"/root");
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
    //status=H5Dread(hdf_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, var->Knlm);
    status=H5Dread(hdf_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
    for (int i=0; i<dims[0]; i++){
        for (int j=0; j<dims[1]; j++){
            x[j]=temp[i*dims[1]+j];
        }
        //if (i==0) fprintf(stdout,"%e,%e,%e,%e,%e,%e\n",x[0],x[1],x[2],x[3],x[4],x[5]);
        init.push_back(x);
    }

    delete [] temp;
    H5Gclose(hdf_group);
    H5Fclose(hdf_file);
    H5Dclose(hdf_data);
    fprintf(stdout," file read successfully!\n");
    return 0;
}


int saveHdf5(char *filename, state_type &OUT, state_type &TIME){

    hid_t   hdf_file,hdf_group,hdf_data,dataspace_id;
    herr_t  status;
    char    name[150];

    //snprintf(name,150,"%s.hdf5",filename);
    fprintf(stdout,"Writing file %s ...",filename);
    hdf_file = H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    if (hdf_file < 0){
        return -1;
    }

    //snprintf(name,150,"/root");
    /*
    if ( (hdf_group=H5Gopen2(hdf_file,"/",H5P_DEFAULT)) < 0){
        H5Fclose(hdf_file);
        return -1;
    }*/
    
    /* Write particle positions and velocities.
     * Ordered in chunk where first Nstep+1 lines correspond to particle 1,
     * step Nstep+1 chunk correspond to particle 2 etc.
     */
    
    hsize_t dims[1]={OUT.size()};
    //fprintf(stdout,"Dim of particle data = %d,%d\n",dims[0],dims[1]);
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

    H5Fclose(hdf_file);
    H5Dclose(hdf_data);
    H5Sclose(dataspace_id);
    fprintf(stdout," file written successfully!\n");
    return 0;
}
