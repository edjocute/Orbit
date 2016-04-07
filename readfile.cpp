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



