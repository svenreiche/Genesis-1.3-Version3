#include "readfieldHDF5.h"

ReadFieldHDF5::ReadFieldHDF5(){}

ReadFieldHDF5::~ReadFieldHDF5()
{
  if(isOpen){
    H5Fclose(fid);
  }
}

void ReadFieldHDF5::init(string in_filename, int in_harm)
{
  isOpen=false;
  filename="";
  if (in_filename.size() < 2) { return; }

  filename=in_filename;
  isOpen=true;
  // allocate some work arrays
  harm=in_harm;
  nwork=inputcom_.ncar;
  nwork=nwork*nwork*2;
  work=new double [nwork];
  // open file and read some parameter description
  fid=H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  readDataDouble(fid,(char *)"slicespacing",&slicespacing,1);
  readDataDouble(fid,(char *)"wavelength",&wavelength,1);
  readDataDouble(fid,(char *)"refposition",&refposition,1);
  readDataDouble(fid,(char *)"gridsize",&gridsize,1);
  readDataInt(fid,(char *)"slicecount",&slicecount,1);
  
  // check for consistent input parameters
  if(this->check()==false){
    isOpen=false;   
    H5Fclose(fid);
  };
  
  // calculate the pointed in the dump field for the first slice
  if (inputcom_.alignradf==0){
    offset=f2cgetnslp_();   // add slippage field
  }else{
    offset=inputcom_.offsetradf; // enforce an offset
  }
  return; 
}


void ReadFieldHDF5::readfield(int islice){
  
  if(!isOpen){ return; } // skip if partfile option is not selected


  int irec=islice+offset-1;
  
  double spos=inputcom_.zsep*inputcom_.xlamds*irec;
  
  if ((spos < 0) || (spos >= slicespacing*slicecount)) {
    cout << " Slice outside of field file - using internal generation of wavefront" << endl;
    return;
  }

  irec=round(spos/slicespacing)+1;
  char name[20];
  sprintf(name,"slice%6.6d/field",irec);
  int nsize=getDatasetSize(fid, name);
  cout << "Reading slice: " << irec << " from field dump file" << endl;
  if (nsize!=nwork){ // check of mismatch in grid sizes
    cout << "Grid size NCAR does not match with size of input field" << endl;
    H5Fclose(fid);
    isOpen=false;
  }

  readDataDouble(fid,name,work,nsize);
  int i=1; // number of harmonics
  f2cputfullfield_(work,&i,&harm);

  return;
}

void ReadFieldHDF5::readslippage(){
  
  if(!isOpen){ return; } // skip if partfile option is not selected

  int nslp=f2cgetnslp_();
  for (int islp=1; islp <= nslp; islp++){
    int err=swapfield_(&islp);    // copy slippage field into current field
    int islice=1-islp;            
    this->readfield(islice);      // update current field
    err=swapfield_(&islp);    // copy back the field into the slippage array
  }

  return;

 
}

