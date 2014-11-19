#include <stdlib.h>
#include <string>

#ifndef __GEN_OUTPUTBASE__
#define __GEN_OUTPUTBASE__

#include "genesis_fortran_common.h"

using namespace std;

typedef struct Flag{
  bool power[HARMMAX]; // power
  bool signal[HARMMAX];   // amp for spec
  bool phase[HARMMAX];
  bool bunch[HARMMAX];
  bool farfield[HARMMAX];
  bool diver[HARMMAX];
  bool size[HARMMAX];  
  bool inc[HARMMAX];
  bool bunchphase[HARMMAX];
  bool energy;
  bool xrms;
  bool yrms;
  bool xpos;
  bool ypos;
  bool espread;
  bool error;
} Flag;

class OutputBase{
 public:
  OutputBase();
  virtual ~OutputBase();
  virtual void init(string, int)=0;
  virtual void finalize(int rank, bool isparallel)=0;
  virtual void writeRecord(int, bool, int, int, int)=0;

 protected:
  Flag flag;
  bool isscan;
  int nstepz,firstout;
  // output flags
};


#endif


