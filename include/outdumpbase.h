#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>


#ifndef __GEN_OUTDUMPBASE__
#define __GEN_OUTDUMPBASE__

#include "genesis_fortran_common.h"


using namespace std;

class OutdumpBase{
 public:
  OutdumpBase();
  virtual ~OutdumpBase();
  virtual void init(string, int)=0;
  virtual void finalize(int, bool)=0;
  virtual void writeRecord(int, int, bool,int,int,int)=0;
 protected:
  int firstout;
  int hcount,hmax;
  bool dumppar,dumpfld, writepar, writefld;
  int ipradi,ippart,isradi,ispart;
};



#endif
