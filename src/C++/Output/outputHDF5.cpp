// class outputASCII
// controls main output in ASCII format

#include "outputHDF5.h"

// constructor/destructor

OutputHDF5::OutputHDF5()
{
  fid=-1;
  profile.clear();
}

OutputHDF5::~OutputHDF5()
{
  //  if (fout.is_open()) { fout.close(); }
  // return;
}


//----------------------------------------------------------------------------
// open/close

void OutputHDF5::open(string file)
{
  fid=H5Fcreate(file.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT); // open HDF5 file
  return;
}

void OutputHDF5::close(){
  if (fid>0) { 
    H5Fclose(fid);
    fid=-1;
  }
  return;
}

//---------------------------------------------------------------------------------------
// some final processing/output of the output file - merging in parallel operation

void OutputHDF5::finalize(int rank, bool isparallel)
{
  if (rank==0) {
    int nprof=profile.size();
    double *prof=new double [nprof];
    for (int i=0;i<nprof;i++){
      prof[i]=profile.at(i);
    }
    if (nprof>1){
      if (isscan){
        writeDataDouble(fid,(char *)"scan",prof,nprof);
      } else {
        writeDataDouble(fid,(char *)"current",prof,nprof);
      }
      this->close();   // close all open files
    }
  }
  return;
}



//----------------------------------------------------------------------------
// write record

void OutputHDF5::writeRecord(int islice, bool isparallel, int rank, int size, int nslice)
{

  double *value=new double[10*HARMMAX+10];
  int valsize=10*HARMMAX+10;
  vector<double> record;
  int tag=1;
  int icount;
  double tmp;
  MPI::Status status;

  if (rank >0) {
      if (isscan){
       tmp=f2cgetscanvalue_();
      } else{
       tmp=f2cgetcurrent_();
      }
      MPI::COMM_WORLD.Send(&tmp, 1, MPI::DOUBLE, 0, tag);

      for (int i=0; i<nstepz+1;i+=inputcom_.iphsty){
         int ir=i/inputcom_.iphsty+1;
         f2cgetdiagnorecord_(&ir,value); // get current data record
	 MPI::COMM_WORLD.Send(value, valsize, MPI::DOUBLE, 0, tag);
      }
  } else {
    int nstop=islice+size;
    if (nstop>nslice) {
      nstop=nslice+1;
    }
    for (int j=islice; j < nstop; j++){
      if (j> firstout){
	 if ( j==islice){
            if (isscan){
               tmp=f2cgetscanvalue_();
            } else{
               tmp=f2cgetcurrent_();
            }
         }else{
	   MPI::COMM_WORLD.Recv(&tmp, 1, MPI::DOUBLE, j-islice, tag,status);
         }
         profile.push_back(tmp);
	 //         this->writeRecordHeader(j,tmp);  // write header : labels, current/scan parameter
         record.clear();
         icount=0;
         for (int i=0; i<nstepz+1;i+=inputcom_.iphsty){
           icount++;
           if (j==islice){
             int ir=i/inputcom_.iphsty+1;
             f2cgetdiagnorecord_(&ir,value); // get current data record
	   } else {

	     MPI::COMM_WORLD.Recv(value, valsize, MPI::DOUBLE, j-islice, tag, status);
           }
           for (int i=0;i<valsize;i++){
             record.push_back(value[i]);
           }
	 }
         // write new record
         if (flag.power[0])  {this->expandDataset(fid,&record,0,icount,j,(char *)"power");}
 	 if (flag.inc[0])    {this->expandDataset(fid,&record,1,icount,j,(char *)"increment");}
         if (flag.signal[0]) {this->expandDataset(fid,&record,2,icount,j,(char *)"signalamp");}
         if (flag.phase[0])  {this->expandDataset(fid,&record,3,icount,j,(char *)"signalphase");}
         if (flag.size[0])   {this->expandDataset(fid,&record,4,icount,j,(char *)"radsize");}
         if (flag.diver[0])  {this->expandDataset(fid,&record,5,icount,j,(char *)"divergence");}
         if (flag.energy)    {this->expandDataset(fid,&record,6,icount,j,(char *)"energy");}
         if (flag.bunch[0])  {this->expandDataset(fid,&record,7,icount,j,(char *)"bunching");}
         if (flag.xrms)      {this->expandDataset(fid,&record,8,icount,j,(char *)"xbeamsize");}
         if (flag.yrms)      {this->expandDataset(fid,&record,9,icount,j,(char *)"ybeamsize");}
         if (flag.error)     {this->expandDataset(fid,&record,10,icount,j,(char *)"error");}
         if (flag.xpos)      {this->expandDataset(fid,&record,11,icount,j,(char *)"xposition");}
         if (flag.ypos)      {this->expandDataset(fid,&record,12,icount,j,(char *)"yposition");}
         if (flag.espread)   {this->expandDataset(fid,&record,13,icount,j,(char *)"energyspread");}
         if (flag.farfield[0]) {this->expandDataset(fid,&record,14,icount,j,(char *)"farfield");}
         if (flag.bunchphase[0]) {this->expandDataset(fid,&record,15,icount,j,(char *)"bunchphase");}
         for (int i = 1; i<HARMMAX;i++){
	   char name[20];
           sprintf(name,"bunching%d",i+1);
           if (flag.bunch[i])  {this->expandDataset(fid,&record,16+(i-1)*5,icount,j,name);}
           sprintf(name,"power%d",i+1);
           if (flag.power[i])  {this->expandDataset(fid,&record,17+(i-1)*5,icount,j,name);}
           sprintf(name,"bunchphase%d",i+1);
           if (flag.bunchphase[i])  {this->expandDataset(fid,&record,18+(i-1)*5,icount,j,name);}
           sprintf(name,"signalamp%d",i+1);
           if (flag.signal[i])  {this->expandDataset(fid,&record,19+(i-1)*5,icount,j,name);}
           sprintf(name,"signalphase%d",i+1);
           if (flag.phase[i])  {this->expandDataset(fid,&record,20+(i-1)*5,icount,j,name);}
	 }        
      }
    }
  }
  

  delete [] value;
  return;

}



//--------------------------------------------------------------------------
// global + input are writing the header for the main output file. input can also be used to create a template file.


void OutputHDF5::init(string outputfile, int rank)
{
  // get some local information out of input parameters

  filename=outputfile+".h5";   // extension for hdf5 output file

  // extract some information from common block
  flag.power[0]=inputcom_.lout[0];   //1
  flag.inc[0]=inputcom_.lout[1];     //2
  flag.signal[0]=inputcom_.lout[2];  //3
  flag.phase[0]=inputcom_.lout[3];   //4
  flag.size[0]=inputcom_.lout[4];    //5
  flag.diver[0]=inputcom_.lout[5];   //6
  flag.energy=inputcom_.lout[6];     //7
  flag.bunch[0]=inputcom_.lout[7];   //8
  flag.xrms=inputcom_.lout[8];       //9
  flag.yrms=inputcom_.lout[9];       //10
  flag.error=inputcom_.lout[10];      //11
  flag.xpos=inputcom_.lout[11];       //12
  flag.ypos=inputcom_.lout[12];       //13
  flag.espread=inputcom_.lout[13];    //14
  flag.farfield[0]=inputcom_.lout[14];//15  
  flag.bunchphase[0]=inputcom_.lout[15];//16  
  for (int i=1;i<HARMMAX;i++){  
    flag.power[i]=(inputcom_.lout[i+15]!=0);
    flag.signal[i]=(inputcom_.lout[i+15]!=0);
    flag.phase[i]=(inputcom_.lout[i+15]!=0);
    flag.bunch[i]=(inputcom_.lout[i+15]!=0);
    flag.bunchphase[i]=(inputcom_.lout[i+15]!=0);
  }

  isscan=false;
  if ((inputcom_.iscan > 0)  && ( inputcom_.iscan < 23)){ isscan=true; } // check for scan feature
  nstepz=f2cgetnstepz_();                                       // get integration steps
  firstout=f2cgetnslp_()*(1-inputcom_.iotail)*inputcom_.itdp;   // calculate the first output

  // write global variables
  if (rank==0) {
      this->global(); 
  }


  return;
}


//-------------------------------------------------------------------
// write header of genesis main output file

void OutputHDF5::global()
{

  this->open(filename);
 
  int nrec,nslc;

  // global parameters
  hid_t gid=H5Gcreate(fid,(char *)"global",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
 
  for (int i=15; i< 15+HARMMAX; i++){ // version 1.0 format: Higher harmonics are indicated with a 4, indicating 4 outputs per harmonics
    inputcom_.lout[i]*=4;
  }
  writeDataInt(gid,(char *)"lout",inputcom_.lout,15+HARMMAX);
 
  int itmp=0;
  if (inputcom_.iphsty > 0){ itmp=nstepz/inputcom_.iphsty +1; }// note integer division  
  writeDataInt(gid,(char *)"z-steps",&itmp,1);
  nrec=itmp; 

  itmp=inputcom_.nslice-firstout;
  itmp=itmp/inputcom_.ishsty;
  writeDataInt(gid,(char *)"slices",&itmp,1);
  nslc=itmp;  

  writeDataDouble(gid,(char *)"wavelength",&inputcom_.xlamds,1);
  double ztmp=inputcom_.xlamds*inputcom_.zsep*inputcom_.ishsty;

  writeDataDouble(gid,(char *)"slicelength",&ztmp,1);

  ztmp=f2cgetdxy_()/f2cgetxkw0_();
  writeDataDouble(gid,(char *)"meshsize",&ztmp,1);

  writeDataInt(gid,(char *)"npart",&inputcom_.npart,1);

  int nz,nt;
  if (inputcom_.ippart > 0){
    nz= nstepz/inputcom_.ippart+1;
    nt= itmp/inputcom_.ispart;
  }
  else{
    nz=0;
    nt=0;
  }
  writeDataInt(gid,(char *)"part-nz",&nz,1);
  writeDataInt(gid,(char *)"part-ns",&nz,1);
  
  if (inputcom_.ipradi > 0){
    nz= nstepz/inputcom_.ipradi+1;
    nt= itmp/inputcom_.isradi;
  }
  else{
    nz=0;
    nt=0;
  }
  writeDataInt(gid,(char *)"field-nz",&nz,1);
  writeDataInt(gid,(char *)"field-ns",&nz,1);
  H5Gclose(gid); 

  //
  // input parameters  
  this->input();

  //
  // magnetic lattice 
  double *aw = new double [NZMAX +1];
  double *qf = new double [NZMAX +1];
  double *zpos = new double [NZMAX+1];
  double dz=f2cgetmaglattice_(aw,qf);  
  int icount=0;
  for (int i=0;i<nstepz+1;i+=inputcom_.iphsty){
    zpos[icount]=dz*i;
    icount++;
  }
  gid=H5Gcreate(fid,(char *)"lattice",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  writeDataDouble(gid,(char *)"z",zpos,icount);
  writeDataDouble(gid,(char *)"aw",aw,icount);
  writeDataDouble(gid,(char *)"qfld",qf,icount);
  H5Gclose(gid); 
  delete [] aw;
  delete [] qf;
  delete [] zpos; 

  // create expendable datasets

  if (flag.power[0])  {this->createExpDataset(fid,(char *)"power",nrec,nslc);}
  if (flag.inc[0])    {this->createExpDataset(fid,(char *)"increment",nrec,nslc);}
  if (flag.signal[0]) {this->createExpDataset(fid,(char *)"signalamp",nrec,nslc);}
  if (flag.phase[0])  {this->createExpDataset(fid,(char *)"signalphase",nrec,nslc);}
  if (flag.size[0])   {this->createExpDataset(fid,(char *)"radsize",nrec,nslc);}
  if (flag.diver[0])  {this->createExpDataset(fid,(char *)"divergence",nrec,nslc);}
  if (flag.energy)    {this->createExpDataset(fid,(char *)"energy",nrec,nslc);}
  if (flag.bunch[0])  {this->createExpDataset(fid,(char *)"bunching",nrec,nslc);}
  if (flag.xrms)      {this->createExpDataset(fid,(char *)"xbeamsize",nrec,nslc);}
  if (flag.yrms)      {this->createExpDataset(fid,(char *)"ybeamsize",nrec,nslc);}
  if (flag.error)     {this->createExpDataset(fid,(char *)"error",nrec,nslc);}
  if (flag.xpos)      {this->createExpDataset(fid,(char *)"xposition",nrec,nslc);}
  if (flag.ypos)      {this->createExpDataset(fid,(char *)"yposition",nrec,nslc);}
  if (flag.espread)   {this->createExpDataset(fid,(char *)"energyspread",nrec,nslc);}
  if (flag.farfield[0]) {this->createExpDataset(fid,(char *)"farfield",nrec,nslc);}
  if (flag.bunchphase[0]) {this->createExpDataset(fid,(char *)"bunchphase",nrec,nslc);}
  for (int i = 1; i<HARMMAX;i++){
	  char name[20];
          sprintf(name,"bunching%d",i+1);
          if (flag.bunch[i])  {this->createExpDataset(fid,name,nrec,nslc);}
          sprintf(name,"power%d",i+1);
          if (flag.power[i])  {this->createExpDataset(fid,name,nrec,nslc);}
          sprintf(name,"signalphase%d",i+1);
          if (flag.phase[i])  {this->createExpDataset(fid,name,nrec,nslc);}
          sprintf(name,"signalamp%d",i+1);
          if (flag.signal[i])  {this->createExpDataset(fid,name,nrec,nslc);}
          sprintf(name,"bunchphase%d",i+1);
          if (flag.bunchphase[i])  {this->createExpDataset(fid,name,nrec,nslc);}
   }   
  return;
}

//----------------------------------------------------------------------------------
// repead of main input file

void OutputHDF5::input()
{

  hid_t gid=H5Gcreate(fid,(char *)"input",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  writeDataDouble(gid,(char *)"aw0",&inputcom_.aw0,1);
  writeDataDouble(gid,(char *)"xkx",&inputcom_.xkx,1);
  writeDataDouble(gid,(char *)"xky",&inputcom_.xky,1);
  writeDataDouble(gid,(char *)"wcoefz",inputcom_.wcoefz,3);
  writeDataDouble(gid,(char *)"xlamd",&inputcom_.xlamd,1);
  writeDataDouble(gid,(char *)"fbess0",&inputcom_.fbess0,1);
  writeDataDouble(gid,(char *)"delaw",&inputcom_.delaw,1);
  writeDataInt(gid,(char *)"iertyp",&inputcom_.iertyp,1);
  writeDataInt(gid,(char *)"iwityp",&inputcom_.iwityp,1);
  writeDataDouble(gid,(char *)"awd",&inputcom_.awd,1);
  writeDataDouble(gid,(char *)"awx",&inputcom_.awx,1);
  writeDataDouble(gid,(char *)"awy",&inputcom_.awy,1);
  writeDataInt(gid,(char *)"iseed",&inputcom_.iseed,1);
  writeDataInt(gid,(char *)"npart",&inputcom_.npart,1);
  writeDataDouble(gid,(char *)"gamma0",&inputcom_.gamma0,1);
  writeDataDouble(gid,(char *)"delgam",&inputcom_.delgam,1);
  writeDataDouble(gid,(char *)"rxbeam",&inputcom_.rxbeam,1);
  writeDataDouble(gid,(char *)"rybeam",&inputcom_.rybeam,1);
  writeDataDouble(gid,(char *)"alphax",&inputcom_.alphax,1);
  writeDataDouble(gid,(char *)"alphay",&inputcom_.alphay,1);
  writeDataDouble(gid,(char *)"emitx",&inputcom_.emitx,1);
  writeDataDouble(gid,(char *)"emity",&inputcom_.emity,1);
  writeDataDouble(gid,(char *)"xbeam",&inputcom_.xbeam,1);
  writeDataDouble(gid,(char *)"ybeam",&inputcom_.ybeam,1);
  writeDataDouble(gid,(char *)"pxbeam",&inputcom_.pxbeam,1);
  writeDataDouble(gid,(char *)"pybeam",&inputcom_.pybeam,1);
  writeDataDouble(gid,(char *)"conditx",&inputcom_.conditx,1);
  writeDataDouble(gid,(char *)"condity",&inputcom_.condity,1);
  writeDataDouble(gid,(char *)"bunch",&inputcom_.bunch,1);
  writeDataDouble(gid,(char *)"bunchphase",&inputcom_.bunchphase,1);
  writeDataDouble(gid,(char *)"emod",&inputcom_.emod,1);
  writeDataDouble(gid,(char *)"emodphase",&inputcom_.emodphase,1);
  writeDataDouble(gid,(char *)"xlamds",&inputcom_.xlamds,1);
  writeDataDouble(gid,(char *)"prad0",&inputcom_.prad0,1);
  writeDataDouble(gid,(char *)"pradh0",&inputcom_.pradh0,1);
  writeDataDouble(gid,(char *)"zrayl",&inputcom_.zrayl,1);
  writeDataDouble(gid,(char *)"zwaist",&inputcom_.zwaist,1);
  writeDataInt(gid,(char *)"ncar",&inputcom_.ncar,1);
  writeDataInt(gid,(char *)"lbc",&inputcom_.lbc,1);
  writeDataDouble(gid,(char *)"rmax0",&inputcom_.rmax0,1);
  writeDataDouble(gid,(char *)"dgrid",&inputcom_.dgrid,1);
  writeDataInt(gid,(char *)"nscr",&inputcom_.nscr,1);
  writeDataInt(gid,(char *)"nscz",&inputcom_.nscz,1);
  writeDataInt(gid,(char *)"nptr",&inputcom_.nptr,1);
  writeDataInt(gid,(char *)"nwig",&inputcom_.nwig,1);
  writeDataDouble(gid,(char *)"zsep",&inputcom_.zsep,1);
  writeDataDouble(gid,(char *)"delz",&inputcom_.delz,1);
  writeDataInt(gid,(char *)"nsec",&inputcom_.nsec,1);
  writeDataInt(gid,(char *)"iorb",&inputcom_.iorb,1);
  writeDataDouble(gid,(char *)"zstop",&inputcom_.zstop,1);
  writeDataInt(gid,(char *)"magin",&inputcom_.magin,1);
  writeDataInt(gid,(char *)"magout",&inputcom_.magout,1);
  writeDataDouble(gid,(char *)"quadf",&inputcom_.quadf,1);
  writeDataDouble(gid,(char *)"quadd",&inputcom_.quadd,1);
  writeDataDouble(gid,(char *)"fl",&inputcom_.fl,1);
  writeDataDouble(gid,(char *)"dl",&inputcom_.dl,1);
  writeDataDouble(gid,(char *)"drl",&inputcom_.drl,1);
  writeDataDouble(gid,(char *)"f1st",&inputcom_.f1st,1);
  writeDataDouble(gid,(char *)"qfdx",&inputcom_.qfdx,1);
  writeDataDouble(gid,(char *)"qfdy",&inputcom_.qfdy,1);
  writeDataDouble(gid,(char *)"solen",&inputcom_.solen,1);
  writeDataDouble(gid,(char *)"sl",&inputcom_.sl,1);
  writeDataInt(gid,(char *)"ildgam",&inputcom_.ildgam,1);
  writeDataInt(gid,(char *)"ildpsi",&inputcom_.ildpsi,1);
  writeDataInt(gid,(char *)"ildx",&inputcom_.ildx,1);
  writeDataInt(gid,(char *)"ildy",&inputcom_.ildy,1);
  writeDataInt(gid,(char *)"ildpx",&inputcom_.ildpx,1);
  writeDataInt(gid,(char *)"ildpy",&inputcom_.ildpy,1);
  writeDataInt(gid,(char *)"itgaus",&inputcom_.itgaus,1);
  writeDataInt(gid,(char *)"nbins",&inputcom_.nbins,1);
  writeDataInt(gid,(char *)"igamgaus",&inputcom_.igamgaus,1);
  writeDataInt(gid,(char *)"inverfc",&inputcom_.inverfc,1);
  writeDataInt(gid,(char *)"ione4one",&inputcom_.ione4one,1);
  writeDataInt(gid,(char *)"lout",inputcom_.lout,15+HARMMAX) ;
  writeDataInt(gid,(char *)"iphsty",&inputcom_.iphsty,1);
  writeDataInt(gid,(char *)"ishsty",&inputcom_.ishsty,1);
  writeDataInt(gid,(char *)"ippart",&inputcom_.ippart,1);
  writeDataInt(gid,(char *)"ispart",&inputcom_.ispart,1);
  writeDataInt(gid,(char *)"ipradi",&inputcom_.ipradi,1);
  writeDataInt(gid,(char *)"isradi",&inputcom_.isradi,1);
  writeDataInt(gid,(char *)"idump",&inputcom_.idump,1);
  writeDataInt(gid,(char *)"iotail",&inputcom_.iotail,1);
  writeDataInt(gid,(char *)"nharm",&inputcom_.nharm,1);
  writeDataInt(gid,(char *)"iallharm",&inputcom_.iallharm,1);
  writeDataInt(gid,(char *)"iharmsc",&inputcom_.iharmsc,1);
  writeDataDouble(gid,(char *)"curpeak",&inputcom_.curpeak,1);
  writeDataDouble(gid,(char *)"curlen",&inputcom_.curlen,1);
  writeDataInt(gid,(char *)"ntail",&inputcom_.ntail,1);
  writeDataInt(gid,(char *)"nslice",&inputcom_.nslice,1);
  writeDataDouble(gid,(char *)"shotnoise",&inputcom_.shotnoise,1);
  writeDataInt(gid,(char *)"isntyp",&inputcom_.isntyp,1);
  writeDataInt(gid,(char *)"iall",&inputcom_.iall,1);
  writeDataInt(gid,(char *)"itdp",&inputcom_.itdp,1);
  writeDataInt(gid,(char *)"ipseed",&inputcom_.ipseed,1);
  writeDataInt(gid,(char *)"iscan",&inputcom_.iscan,1);
  writeDataInt(gid,(char *)"nscan",&inputcom_.nscan,1);
  writeDataDouble(gid,(char *)"svar",&inputcom_.svar,1);
  writeDataInt(gid,(char *)"isravg",&inputcom_.isravg,1);
  writeDataInt(gid,(char *)"isrsig",&inputcom_.isrsig,1);
  writeDataDouble(gid,(char *)"cuttail",&inputcom_.cuttail,1);
  writeDataDouble(gid,(char *)"eloss",&inputcom_.eloss,1);
  writeDataDouble(gid,(char *)"version",&inputcom_.version,1);
  writeDataInt(gid,(char *)"ndcut",&inputcom_.ndcut,1);
  writeDataInt(gid,(char *)"idmpfld",&inputcom_.idmpfld,1);
  writeDataInt(gid,(char *)"idmppar",&inputcom_.idmppar,1);
  writeDataInt(gid,(char *)"ilog",&inputcom_.ilog,1);
  writeDataInt(gid,(char *)"ffspec",&inputcom_.ffspec,1);
  writeDataDouble(gid,(char *)"ibfield",&inputcom_.ibfield,1);
  writeDataDouble(gid,(char *)"imagl",&inputcom_.imagl,1);
  writeDataDouble(gid,(char *)"idril",&inputcom_.idril,1);
  writeDataInt(gid,(char *)"alignradf",&inputcom_.alignradf,1);
  writeDataInt(gid,(char *)"offsetradf",&inputcom_.offsetradf,1);
  writeDataInt(gid,(char *)"multconv",&inputcom_.multconv,1);
  H5Gclose(gid); 
}


