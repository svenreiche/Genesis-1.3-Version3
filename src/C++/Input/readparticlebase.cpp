
#include "readparticlebase.h"

ReadParticleBase::ReadParticleBase(){}
ReadParticleBase::~ReadParticleBase(){}

bool ReadParticleBase::check()
{
  if (inputcom_.nbins!=beamletsize) { // mismatch between bins requested and in file
    cout << "NBINS does not match beamlet size in partfile" << endl;
    return false; 
   }
  double ds=inputcom_.xlamds*inputcom_.zsep;
  if (inputcom_.itdp==0){
    ds=slicespacing;    // for no time-dependent simulation the slice length doesn-t matter
  }

  if ((fabs(ds-slicespacing)/ds)>1e-3) {  // check whether first required slice is within the dum 
    cout << "Selected slice length does not match length of input file" << endl;
    return false; 
  } 

  if ((inputcom_.nslice+inputcom_.ntail)>slicecount) { // check whether last slice is inside the dump
    cout << "Time window larger than partfile time-window" << endl;
    return false; 
   }
  if ((inputcom_.ntail)<0) { // check whether last slice is inside the dump
    cout << "First slice before first entry of partfile time-window" << endl;
    return false; 
   }

  double tmp=(slicelength-inputcom_.xlamds*inputcom_.convharm)/slicelength;
  if (fabs(tmp)>0.01) { // check for wavelength match
    cout << "Slice length not consistent with slice length in partfile" << endl;
    return false; 
  }
  return true;
}
