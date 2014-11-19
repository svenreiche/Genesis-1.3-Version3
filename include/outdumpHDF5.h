
#ifndef __GEN_OUTDUMPHDF5__
#define __GEN_OUTDUMPHDF5__


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

#include "outdumpbase.h"
#include "genesis_fortran_calls.h"

#include "HDF5base.h"
#include "hdf5.h"
#include "mpi.h"

using namespace std;

class OutdumpHDF5 : public OutdumpBase, public HDF5Base {
 public:
  OutdumpHDF5();
  virtual ~OutdumpHDF5();
  void init(string, int );
  void finalize(int, bool);
  void writeRecord(int, int, bool,int,int,int);

 private:
  hid_t createFieldFile(string, const char *,int);
  hid_t createPartFile(string, const char *);
  void writeBeam(hid_t, int, int, double *, double);
  void writeField(hid_t, int, int, double *);
  void updateWorkArray(int);
  void dumpSlippageField();

  int npart,nfield, nslp,nwork;  
  double *work;
  hid_t fiddpa, fidpar;
  vector<hid_t> fiddfl,fidfld;

};


#endif
