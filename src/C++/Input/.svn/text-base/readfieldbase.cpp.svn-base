
#include "readfieldbase.h"

ReadFieldBase::ReadFieldBase(){}
ReadFieldBase::~ReadFieldBase(){}

bool ReadFieldBase::check()
{

  double ds=inputcom_.xlamds*inputcom_.zsep;
  double dd=ds/slicespacing-round(ds/slicespacing);
  if (fabs(dd/ds)>1e-3) {  // check whether first required slice is within the dump 
    cout << "Selected slice length does not match length of input file" << endl;
    return false; 
  } 

  ds=inputcom_.xlamds/inputcom_.nharm;
  if ((fabs(ds-wavelength)/ds)>1e-3) {  // check whether first required slice is within the dump 
    cout << "Selected wave length does not match wave length of input file" << endl;
    return false; 
  } 

  // this is not quite needed
  //if (inputcom_.nslice>slicecount) { // check whether last slice is inside the dump
  //  cout << "Time window larger than partfile time-window" << endl;
  //  return false; 
  // }


  double tmp=f2cgetdxy_()/f2cgetxkw0_();  
  if ((fabs(tmp-gridsize)/gridsize)>0.001) { // check for wavelength match
    cout << "Slice length not consistent with slice length in partfile" << endl;
    return false; 
  }
  return true;
}
