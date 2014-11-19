

#ifndef __GEN_OUTPUTHDF5__
#define __GEN_OUTPUTHDF5__


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

#include "mpi.h"
#include "HDF5base.h"
#include "hdf5.h"

#include "outputbase.h"

#include "genesis_fortran_calls.h"

using namespace std;

class OutputHDF5: public OutputBase, public HDF5Base{
 public:
  OutputHDF5();
  virtual ~OutputHDF5();
  void open(string );
  void close();
  void global();
  void input(); 
  void init(string, int);
  void finalize(int rank, bool isparallel);
  void writeRecord(int, bool, int, int, int);
 private:
  string filename;
  vector<double> profile;
  hid_t fid;
};



#endif
