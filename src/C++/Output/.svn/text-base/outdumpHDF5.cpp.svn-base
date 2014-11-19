
#include "outdumpHDF5.h"
#include <stdlib.h>

// constructor destructor
OutdumpHDF5::OutdumpHDF5()
{
  work=NULL;
}

OutdumpHDF5::~OutdumpHDF5()
{
  if (work!=NULL){
    delete [] work;
  }
}



void OutdumpHDF5::init(string outputfile, int rank)
{
  nwork=0;
  npart=6*inputcom_.npart;      // initial guess for size of macro particle record
  nfield=2*inputcom_.ncar*inputcom_.ncar;   // size for field record
  nslp=f2cgetnslp_()*inputcom_.itdp;        // guess for the slippage field
  dumppar=(inputcom_.idmppar != 0);         // flag whether particle dump is requested
  dumpfld=(inputcom_.idmpfld != 0);         // flag whether field dump is requested
  writepar=(inputcom_.ippart != 0);         // flag whether particle output is req.  
  writefld=(inputcom_.ipradi != 0);         // flag whether field output is requested
  ipradi=inputcom_.ipradi;
  isradi=inputcom_.isradi;
  ippart=inputcom_.ippart;
  ispart=inputcom_.ispart;

  firstout=f2cgetnslp_()*(1-inputcom_.iotail)*inputcom_.itdp; // record of first output (not dumps)
  int nharm=inputcom_.nharm;  // harmonics
  int iallharm=inputcom_.iallharm;

  hcount=1;
  hmax=1;
  if (nharm>1){
    hmax=nharm;
    if (iallharm==0){
      hcount=2;
    }else{
      hcount=nharm;
    }
  }

  // allocate working memory
  if ((dumppar) || (writepar)) { 
    this->updateWorkArray(npart);              // header for dumping particle
  }
  if ((dumpfld) || (writefld)){
    this->updateWorkArray(nfield*hcount);
  }
  if (rank>0) { return; }



  //---------------------------------
  // master node is opening dump files

  hid_t fid;


  if (dumppar){
    fiddpa=this->createPartFile(outputfile,".dpa.h5");
  }
  if (writepar){
    fidpar=this->createPartFile(outputfile,".par.h5");
  }
 
  if (dumpfld){
    fid=this->createFieldFile(outputfile,".dfl.h5",1);
    fiddfl.push_back(fid); 

    for (int i=2;i<=hcount;i++){
      int harm;
      if (i==hcount){
	harm=hmax;
      }else{
	harm=i;
      }
      char ext[15];
      sprintf(ext,".dfl.harm%d.h5",harm);
      fid=this->createFieldFile(outputfile,ext,harm);
      fiddfl.push_back(fid);    
    }
  }

  if (writefld){
    fid=this->createFieldFile(outputfile,".fld.h5",1);
    fidfld.push_back(fid); 

    for (int i=2;i<=hcount;i++){
      int harm;
      if (i==hcount){
	harm=hmax;
      }else{
	harm=i;
      }
      char ext[15];
      sprintf(ext,".fld.harm%d.h5",harm);
      fid=this->createFieldFile(outputfile,ext,harm);
      fidfld.push_back(fid);    
    }
  }
  return;
}

//------------------------------------
// creating the dumps files

  hid_t OutdumpHDF5::createPartFile(string filename,const char *ext){

    filename.append(ext);

    hid_t fid=H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT); // open HDF5 file
    double tmp=inputcom_.xlamds*inputcom_.zsep;  // slice separation XLAMDS*ZSEP
    writeDataDouble(fid,(char *)"slicespacing",&tmp,1);
    tmp=inputcom_.xlamds;                        // slice length = XLAMDS 
    writeDataDouble(fid,(char *)"slicelength",&tmp,1);
    tmp=(inputcom_.ntail+firstout)*inputcom_.xlamds*inputcom_.zsep;  // positin of bunch for first slice
    writeDataDouble(fid,(char *)"refposition",&tmp,1);
    int inttmp=inputcom_.nslice-firstout;       // total number of records
    writeDataInt(fid,(char *)"slicecount",&inttmp,1);
    inttmp=inputcom_.nbins;                     // size of beamlet
    writeDataInt(fid,(char *)"beamletsize",&inttmp,1);
    return fid;
}

hid_t OutdumpHDF5::createFieldFile(string filename, const char *ext,  int harm){

      filename.append(ext);

      hid_t  fid=H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      double tmp=inputcom_.xlamds*inputcom_.zsep;  // slice separation XLAMDS*ZSEP
      writeDataDouble(fid,(char *)"slicespacing",&tmp,1);
      tmp=inputcom_.xlamds/static_cast<double>(harm);                 // reference wavelength 
      writeDataDouble(fid,(char *)"wavelength",&tmp,1);
      tmp=(inputcom_.ntail+firstout)*inputcom_.xlamds*inputcom_.zsep;  // positin of bunch for first slice
      writeDataDouble(fid,(char *)"refposition",&tmp,1);
      int inttmp=inputcom_.nslice-firstout+nslp;       // total number of records
      writeDataInt(fid,(char *)"slicecount",&inttmp,1);
      tmp=f2cgetdxy_()/f2cgetxkw0_();             // grid spacing
      writeDataDouble(fid,(char *)"gridsize",&tmp,1); 
      return fid;
}

//--------------------------
// cleaning up at the end of the fun

void OutdumpHDF5::finalize(int rank, bool isparallel)
{
  if (rank==0){
    if (dumppar){ H5Fclose(fiddpa); }
    if (writepar){ H5Fclose(fidpar); }
    if (dumpfld){
      this->dumpSlippageField();
      for (int i=0;i<fiddfl.size();i++){
	H5Fclose(fiddfl.at(i));
      }
    }
    if (writefld){
      for (int i=0;i<fidfld.size();i++){
	H5Fclose(fidfld.at(i));
      }
    }
  }
  
  delete [] work;
  work=NULL;

  return;  
}



//------------------------------
// slippage field

void OutdumpHDF5::dumpSlippageField()
{

  for (int j=nslp-1;j>0;j--){
    f2cgetslippagefield_(work,&j,&hcount);
    for (int i=0; i<fiddfl.size();i++){
        this->writeField(fiddfl.at(i), j+inputcom_.nslice, nfield, &work[i*nfield]);
    }
  }
  return;
}

void OutdumpHDF5::updateWorkArray(int ntmp){
   
  if (ntmp>nwork){
    delete [] work;
    work = new double [ntmp];
    nwork=ntmp;
  }
  return;

}

//----------------------------------------------------------------------------------
// write single record

void OutdumpHDF5::writeRecord(int istepz, int islice, bool isparallel, int rank, int size, int nslice)
{

  MPI::Status status;
  int tag=1;
  int nstop=islice+size;
  if (nstop>nslice) {
    nstop=nslice+1;  //  because loop goes till "< nstop"
  }

  //------------------------
  // dumping at end of undulator

  if (istepz<0){

   //-----------
    // dump particle
   if (dumppar){
  
    npart=6*inputcom_.npart;   // number of valid macro particles
    this->updateWorkArray(npart);
    double current=f2cgetcurrent_();
    f2cgetbeam_(work);
    
    if (rank >0) {
       MPI::COMM_WORLD.Send(&current, 1, MPI::DOUBLE, 0, tag);
       MPI::COMM_WORLD.Send(&npart, 1, MPI::INT, 0, tag);
       MPI::COMM_WORLD.Send(work, npart, MPI::DOUBLE, 0, tag);
    } else {
       this->writeBeam(fiddpa, islice, npart, work, current);
       for (int j=islice+1; j < nstop; j++){
          cout << "MPI: Gathering particle dumps for slices " << islice +1 << " to " << nstop-1 << endl;
          MPI::COMM_WORLD.Recv(&current, 1, MPI::DOUBLE, j-islice, tag,status);
          MPI::COMM_WORLD.Recv(&npart, 1, MPI::INT, j-islice, tag,status);
          this->updateWorkArray(npart); 
	  MPI::COMM_WORLD.Recv(work, npart, MPI::DOUBLE, j-islice, tag, status);
          this->writeBeam(fiddpa, j, npart, work, current);
       }
    }      
   }
 
  //--------------
  // dump field
   if (dumpfld){

    int first=1;
    f2cgetfullfield_(work,&hcount,&first);

    if (rank>0) {
       MPI::COMM_WORLD.Send(work, nfield*hcount, MPI::DOUBLE, 0, tag);
    } else {
      for (int i=0; i<fiddfl.size();i++){
        this->writeField(fiddfl.at(i), islice, nfield, &work[i*nfield]);
      }
      for (int j=islice+1; j < nstop; j++){
         cout << "MPI: Gathering field dumps for slices " << islice +1 << " to " << nstop-1 << endl;
	 MPI::COMM_WORLD.Recv(work, nfield*hcount, MPI::DOUBLE, j-islice, tag, status);
         for (int i=0; i<fiddfl.size();i++){
           this->writeField(fiddfl.at(i), j, nfield, &work[i*nfield]);
         }
      }
    }
   }

   return;  
  }

  //--------------------------
  // output during the integration   
  
  //-----------------
  // writing particles

  if ((ippart>0) && ((istepz % ippart) == 0)){
    
    npart=6*inputcom_.npart;   // number of valid macro particles
    this->updateWorkArray(npart);
    double current=f2cgetcurrent_();
    f2cgetbeam_(work);
    if (rank >0) {
      if ((islice>firstout) && ((islice % ispart)==0)){
      MPI::COMM_WORLD.Send(&current, 1, MPI::DOUBLE, 0, tag);
      MPI::COMM_WORLD.Send(&npart, 1, MPI::INT, 0, tag);
      MPI::COMM_WORLD.Send(work, npart, MPI::DOUBLE, 0, tag);
      }
    } else {
      hid_t gid;
      char name[20];
      sprintf(name,"step%6.6d",istepz);  
      if (islice==1){
        gid=H5Gcreate(fidpar,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      } else {
        gid=H5Gopen(fidpar,name,H5P_DEFAULT);
      }
      if ((islice>firstout) && ((islice % ispart)==0)){
	this->writeBeam(gid, islice, npart, work, current);
      }
      for (int j=islice+1; j < nstop; j++){
        if ((j>firstout) && ((j % ispart)==0)){
          cout << "MPI: Gathering particle dumps for slices " << islice +1 << " to " << nstop-1 << endl;
          MPI::COMM_WORLD.Recv(&current, 1, MPI::DOUBLE, j-islice, tag,status);
          MPI::COMM_WORLD.Recv(&npart, 1, MPI::INT, j-islice, tag,status);
          this->updateWorkArray(npart); 
	  MPI::COMM_WORLD.Recv(work, npart, MPI::DOUBLE, j-islice, tag, status);
	  this->writeBeam(gid, j, npart, work, current);
        }
      }
      H5Gclose(gid);
    }      
  }

  //--------------
  // writing field
  if ((ipradi>0)&&((istepz % ipradi) == 0)){
    int first=1;
    f2cgetfullfield_(work,&hcount,&first);

    if (rank>0) {
      if ((islice>firstout) && ((islice % isradi)==0)){
        MPI::COMM_WORLD.Send(work, nfield*hcount, MPI::DOUBLE, 0, tag);
      }
    } else {
      vector<hid_t> gidfld;
      hid_t gid;
      char name[20];
      sprintf(name,"step%6.6d",istepz);  
      if (islice==1){
        for (int i=0;i<fidfld.size();i++){
          gid=H5Gcreate(fidfld.at(i),name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          gidfld.push_back(gid);
	}
      } else {
        for (int i=0;i<fidfld.size();i++){
          gid=H5Gopen(fidfld.at(i),name,H5P_DEFAULT);
          gidfld.push_back(gid);
	}
      }

      if ((islice>firstout) && ((islice % isradi)==0)){
        for (int i=0; i<fiddfl.size();i++){
          this->writeField(gidfld.at(i), islice, nfield, &work[i*nfield]);
	}
      }
      for (int j=islice+1; j < nstop; j++){
        if ((j>firstout) && ((j % isradi)==0)){
          cout << "MPI: Gathering field dumps for slices " << islice +1 << " to " << nstop-1 << endl;
	  MPI::COMM_WORLD.Recv(work, nfield*hcount, MPI::DOUBLE, j-islice, tag, status);
          for (int i=0; i<fiddfl.size();i++){
             this->writeField(gidfld.at(i), j, nfield, &work[i*nfield]);
          }
        }
      }
      for (int i=0;i<fidfld.size();i++){
	H5Gclose(gidfld.at(i));
      }
      gidfld.clear();
    }

  }
  return;
}


void OutdumpHDF5::writeBeam(hid_t fid, int islice, int npart, double *beam, double current)
{

  int nbeam=npart/6;
   char name[20];
   sprintf(name,"slice%6.6d",islice);
   hid_t gid=H5Gcreate(fid,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
   writeDataDouble(gid,(char *)"current",&current,1);
   writeDataDouble(gid,(char *)"gamma",&beam[0],nbeam);
   writeDataDouble(gid,(char *)"theta",&beam[nbeam],nbeam);
   writeDataDouble(gid,(char *)"x",&beam[2*nbeam],nbeam);
   writeDataDouble(gid,(char *)"y",&beam[3*nbeam],nbeam);
   writeDataDouble(gid,(char *)"px",&beam[4*nbeam],nbeam);
   writeDataDouble(gid,(char *)"py",&beam[5*nbeam],nbeam);
   H5Gclose(gid);
   return;
    
}



void OutdumpHDF5::writeField(hid_t fid, int islice, int nsize, double *field)
{

   char name[20];
   sprintf(name,"slice%6.6d",islice);
   hid_t gid=H5Gcreate(fid,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
   writeDataDouble(gid,(char *)"field",field,nsize);
   H5Gclose(gid);
   return;
    
}
