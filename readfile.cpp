#include "main.h"
#include <math.h>
#include <hdf5.h>

int load_hdf5_input(char *filename, struct Indata *var){

    hid_t   hdf_file,hdf_group,hdf_data;
    herr_t  status;
    char    name[150];

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
    status=H5Dread(hdf_data, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, var->Knlm);
    
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



