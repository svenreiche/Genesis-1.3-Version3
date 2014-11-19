#include "readparticleHDF5.h"

ReadParticleHDF5::ReadParticleHDF5(){}

ReadParticleHDF5::~ReadParticleHDF5()
{
  if (isOpen){
    H5Fclose(fid);
  }
}

void ReadParticleHDF5::init(string in_filename)
{
  isOpen=false;
  filename="";
  if (in_filename.size() < 2) { return; }

  filename=in_filename;
  isOpen=true;
  // allocate some work arrays
  nwork=inputcom_.npart;
  work=new double [nwork];
  // open file and read some parameter description
  fid=H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  readDataDouble(fid,(char *)"slicespacing",&slicespacing,1);
  readDataDouble(fid,(char *)"slicelength",&slicelength,1);
  readDataDouble(fid,(char *)"refposition",&refposition,1);
  readDataInt(fid,(char *)"slicecount",&slicecount,1);
  readDataInt(fid,(char *)"beamletsize",&beamletsize,1);
  
  // check for consistent input parameters
  if(this->check()==false){
    isOpen=false;   
    H5Fclose(fid);
  };
  
  return; 
}

void ReadParticleHDF5::readpart(int islice){
  
  if(!isOpen){ return; } // skip if partfile option is not selected

  double current;
  //  double spos=(islice-1+inputcom_.ntail)*inputcom_.xlamds*inputcom_.zsep-refposition;
  //  int index=static_cast<int> (rint(spos/slicespacing))+1;
  char name[20];
  sprintf(name,"slice%6.6d/gamma",islice);
  int nsize=getDatasetSize(fid, name);
  cout << "Reading slice: " << islice << " from particle dump file" << endl;
  if (nsize>nwork){ // allocate extra work array to hold field
    delete [] work;
    nwork=nsize;
    work=new double [nwork];
  }

  sprintf(name,"slice%6.6d/current",islice);
  readDataDouble(fid,name,&current,1);
  cout << "Current: " << current << endl;
  f2cputcurrent_(&current);


  for (int i=1;i<=6;i++){
     switch (i){
   	 case 1:
            sprintf(name,"slice%6.6d/gamma",islice);
	    break;
	 case 2:
            sprintf(name,"slice%6.6d/theta",islice);
	    break;
	 case 3:
 	    sprintf(name,"slice%6.6d/x",islice);
	    break;
	 case 4:
	    sprintf(name,"slice%6.6d/y",islice);
	    break;
	 case 5:
	    sprintf(name,"slice%6.6d/px",islice);
 	    break;
	 case 6:
	    sprintf(name,"slice%6.6d/py",islice);
	    break;
     }
     readDataDouble(fid,name,work,nsize);
     f2cputbeam_(work,&i,&nsize);
  }
  readpart_(&nsize);
  return;
}

