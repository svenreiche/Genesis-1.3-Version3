
#include "outdumpBinary.h"
#include <stdlib.h>

// constructor destructor
OutdumpBinary::OutdumpBinary(){}
OutdumpBinary::~OutdumpBinary(){}

void OutdumpBinary::init(string outputfile,int rank)
{
  filename=outputfile;
  npart=inputcom_.npart;
  nfield=2*inputcom_.ncar*inputcom_.ncar;
  nslp=f2cgetnslp_()*inputcom_.itdp;
  dumppar=(inputcom_.idmppar != 0);
  dumpfld=(inputcom_.idmpfld != 0);
  firstout=f2cgetnslp_()*(1-inputcom_.iotail)*inputcom_.itdp;
  if (dumppar) { workp=new double [npart]; }
  if (dumpfld) { workf=new double [nfield];}
  return;
}


void OutdumpBinary::finalize(int rank, bool isparallel)
{
  
  char cfile[100];
  ifstream fin; 
  int count;
  
  // particle ----------------------------
  if (dumppar){
    if (parout.is_open()) { parout.close();} 
    if ((isparallel)&&(rank==0)) { 
      sprintf(cfile,"%s.dpa",filename.c_str());
      parout.open(cfile,ios::binary); 
      count=npart*sizeof(double);
      for (int islice=1+firstout; islice <=inputcom_.nslice ; islice++){
         sprintf(cfile,"%s.dpa.slice%6.6d",filename.c_str(),islice);
         fin.open(cfile,ios::binary);
         for (int i=0;i<6;i++){
            fin.read((char *)workp,count);
            parout.write((char *)workp,count);
         }
         fin.close();
         sprintf(cfile,"rm %s.dpa.slice%6.6d",filename.c_str(),islice);
	 system(cfile);
      }
      parout.close();
    }
    delete [] workp;
  }

  // field ----------------------------------
  if (dumpfld){
    if (fldout.is_open()) { fldout.close();} 
    if ((isparallel)&&(rank==0)) { 
      sprintf(cfile,"%s.dfl",filename.c_str());
      fldout.open(cfile,ios::binary); 
      count=nfield*sizeof(double);
      for (int islice=1+firstout; islice <=inputcom_.nslice ; islice++){
         sprintf(cfile,"%s.dfl.slice%6.6d",filename.c_str(),islice);
         fin.open(cfile,ios::binary);
         fin.read((char *)workf,count);
         fin.close();
         fldout.write((char *)workf,count);
         sprintf(cfile,"rm %s.dfl.slice%6.6d",filename.c_str(),islice);
	 system(cfile);
      }
      this->dumpSlippageField();
      fldout.close();

    }
    delete [] workf;
  }
  return;
}

void OutdumpBinary::dumpSlippageField()
{
  if (inputcom_.nslice <= firstout) { return ; }
  int count=nfield*sizeof(double);
  for (int i=nslp-1;i>0;i--){
    f2cgetslippagefield_(workf,&i,&i);
    fldout.write((char *)workf,count);
  }
}


void OutdumpBinary::writeRecord(int islice, bool isparallel, int a, int b , int c)
{
  this->writeRecord(islice, isparallel);
  return;
}


void OutdumpBinary::writeRecord(int islice, bool isparallel)
{

  char cfile[100];
  int count;

  // particle ------------------------
  if ((dumppar)&&(islice > firstout)){
    if (isparallel){
      if (parout.is_open()) {parout.close();}
      sprintf(cfile,"%s.dpa.slice%6.6d",filename.c_str(),islice);
      parout.open(cfile,ios::binary);    
    } else {
      if (!parout.is_open()){
        sprintf(cfile,"%s.dpa",filename.c_str());
        parout.open(cfile,ios::binary);    
      }
    }

    for (int i=0; i<npart;i++){
      workp[i]=-1.;
    }
    count=npart*sizeof(double);  
    for (int i=1;i<=6;i++){
        f2cgetbeam_(workp,&i);
        parout.write((char *)workp,count);
     }
  }

  // field ---------------------------

  if ((dumpfld)&&(islice > firstout)){
    if (isparallel){
      if (fldout.is_open()) {fldout.close();}
      sprintf(cfile,"%s.dfl.slice%6.6d",filename.c_str(),islice);
      fldout.open(cfile,ios::binary);    
    } else {
      if (!fldout.is_open()){
        sprintf(cfile,"%s.dfl",filename.c_str());
        fldout.open(cfile,ios::binary);    
      }
    }
    int i=0;
    count=nfield*sizeof(double);  
    f2cgetfield_(workf,&i);
    fldout.write((char *)workf,count);
  }


  return; 

  

}
