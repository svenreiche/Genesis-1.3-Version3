#include "HDF5base.h"

HDF5Base::HDF5Base(){}
HDF5Base::~HDF5Base(){}

//----------------------
// generating expandable dataset

void HDF5Base::createExpDataset(hid_t fid, char *name, int nz, int ns)
{
  hsize_t dims[2] = {nz, ns};
  hsize_t maxdims[2]={H5S_UNLIMITED, H5S_UNLIMITED};
  hid_t dataspace = H5Screate_simple(2,dims,maxdims);
  hid_t cparms=H5Pcreate(H5P_DATASET_CREATE);
  hsize_t chunk_dims[2]={10,10};
  herr_t status=H5Pset_chunk( cparms, 2, chunk_dims);
  hid_t dataset=H5Dcreate(fid,name,H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT,cparms,H5P_DEFAULT);

  H5Dclose(dataset);
  H5Sclose(dataspace);
  
  return;

}

void HDF5Base::expandDataset(hid_t fid, vector<double> *rec, int pos, int recsize, int slice, char *name)
{

  double *data= new double[recsize];
  
  int stride=rec->size()/recsize;
  for(int i=0;i<recsize;i++){
    data[i]=rec->at(i*stride+pos);
    //    cout << data[i] << endl;
  }
  //  return;

  hid_t did=H5Dopen(fid,name,H5P_DEFAULT);
  hsize_t size[2] = {recsize, slice};  
  herr_t status=H5Dset_extent(did,size); 

  hid_t filespace=H5Dget_space(did);  
  hsize_t offset[2] = {0,slice-1};
  hsize_t dims[2]={recsize,1};
  hid_t dataspace=H5Screate_simple(2,dims,dims);  

  status=H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,dims,NULL);
  status=H5Dwrite(did,H5T_NATIVE_DOUBLE,dataspace,filespace,H5P_DEFAULT,data);
  
  status=H5Sclose(dataspace);
  status=H5Sclose(filespace);
  status=H5Dclose(did);

  delete [] data;
  return;
}


//------------------------
// writing procedures

void HDF5Base::writeDataDouble(hid_t fid, char *name, double *data, int size)
{
  hsize_t dims[1];
  dims[0]=size;
  hid_t dataspace_id=H5Screate_simple(1,dims,NULL);
  hid_t dataset_id=H5Dcreate(fid,name,H5T_NATIVE_DOUBLE,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  H5Dwrite(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
}

void HDF5Base::writeDataInt(hid_t fid, char *name, int *data, int size)
{
  hsize_t dims[1];
  dims[0]=size;
  hid_t dataspace_id=H5Screate_simple(1,dims,NULL);
  hid_t dataset_id=H5Dcreate(fid,name,H5T_NATIVE_INT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  H5Dwrite(dataset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
}

//--------------------- 
// reading procedures


void HDF5Base::readDataDouble(hid_t fid, char *name, double *data, int size)
{

  hsize_t dims[1];
  dims[0]=size;
  hid_t dataspace_id=H5Screate_simple(1,dims,NULL);
  hid_t dataset_id=H5Dopen(fid,name,H5P_DEFAULT);
  H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset_id);     
  H5Sclose(dataspace_id);
  return;
}

void HDF5Base::readDataInt(hid_t fid, char *name, int *data, int size)
{
  hsize_t dims[1];
  dims[0]=size;
  hid_t dataspace_id=H5Screate_simple(1,dims,NULL);
  hid_t dataset_id=H5Dopen(fid,name,H5P_DEFAULT);
  H5Dread(dataset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset_id);     
  H5Sclose(dataspace_id);
  return;
}


int HDF5Base::getDatasetSize(hid_t fid, char *name)
{

  hsize_t dims[1],maxdims[1];
  hid_t  dsid=H5Dopen(fid,name,H5P_DEFAULT);
  hid_t spaceid=H5Dget_space(dsid);
  H5Sget_simple_extent_dims(spaceid,dims,maxdims);
  
  return dims[0];

}

 
