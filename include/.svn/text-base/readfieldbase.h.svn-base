#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>

#ifndef __GEN_READFIELDBASE__
#define __GEN_READFIELDBASE__

#include "genesis_fortran_common.h"
#include "genesis_fortran_calls.h"

using namespace std;

class ReadFieldBase{
 public:
  ReadFieldBase();
  virtual ~ReadFieldBase();
  virtual void init(string, int)=0;
  virtual void readfield(int)=0;
  virtual void readslippage()=0;
  //  virtual void finalize(int, bool)=0;
  //  virtual void writeRecord(int, bool)=0;
 protected:
  bool check();
  bool isOpen;
  double slicespacing,wavelength,refposition,gridsize;
  int slicecount, offset;
  double *work;
  int nwork;
  int harm;
};



#endif
