
#ifndef __GEN_HDF5BASE__
#define __GEN_HDF5BASE__


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

#include "hdf5.h"


using namespace std;

class HDF5Base{
 public:
  HDF5Base();
  virtual ~HDF5Base();
 protected:
  void createExpDataset(hid_t fid, char *name, int nz, int ns);
  void expandDataset(hid_t fid, vector<double> *rec, int pos, int recsize, int slice, char *name);
  void writeDataDouble(hid_t fid, char *name, double *data, int size);
  void writeDataInt(hid_t fid, char *name, int *data, int size);
  void readDataDouble(hid_t fid, char *name, double *data, int size);
  void readDataInt(hid_t fid, char *name, int *data, int size);
  int getDatasetSize(hid_t fid, char *name);
};



#endif

