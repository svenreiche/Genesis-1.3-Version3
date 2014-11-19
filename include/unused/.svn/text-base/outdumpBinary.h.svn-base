
#ifndef __GEN_OUTDUMPBIN__
#define __GEN_OUTDUMPBIN__


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "outdumpbase.h"

#include "genesis_fortran_calls.h"


using namespace std;

class OutdumpBinary : public OutdumpBase {
 public:
  OutdumpBinary();
  virtual ~OutdumpBinary();
  void init(string,int );
  void finalize(int, bool);
  void writeRecord(int, bool);
  void writeRecord(int, bool, int, int, int);
  void dumpSlippageField();
 private:
  string filename;  
  double *workp,*workf;
  int npart,nfield, nslp;
  ofstream parout, fldout;
};



#endif
