
#include "outputbase.h"

// ----- constructor/destructor

OutputBase::OutputBase()
{
  flag.energy=false;
  flag.xrms=false;
  flag.yrms=false;
  flag.xpos=false;
  flag.ypos=false;
  flag.espread=false;
  flag.error=false;
  for (int i=0; i< HARMMAX; i++){
    flag.power[i]=false;
    flag.signal[i]=false;
    flag.phase[i]=false;
    flag.bunch[i]=false;
    flag.bunchphase[i]=false;
    flag.farfield[i]=false;
    flag.diver[i]=false;
    flag.size[i]=false;
    flag.inc[i]=false;
  }
  return;
}


OutputBase::~OutputBase(){}
