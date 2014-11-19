

#ifndef __GEN_OUTPUTASCII__
#define __GEN_OUTPUTASCII__


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "mpi.h"

#include "outputbase.h"

#include "genesis_fortran_calls.h"

using namespace std;

class OutputASCII: public OutputBase{
 public:
  OutputASCII();
  virtual ~OutputASCII();
  void open(string );
  void close();
  void global();
  void input(); 
  void init(string, int);
  void finalize(int rank, bool isparallel);
  void writeRecord(int, bool, int, int, int);
  void writeRecordHeader(int, double); 
 private:
  ofstream fout; 
  string filename;  

};



#endif
