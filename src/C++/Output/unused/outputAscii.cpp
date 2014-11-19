// class outputASCII
// controls main output in ASCII format

#include "outputAscii.h"

// constructor/destructor

OutputASCII::OutputASCII()
{
}

OutputASCII::~OutputASCII()
{
  if (fout.is_open()) { fout.close(); }
  return;
}

//----------------------------------------------------------------------------
// open/close

void OutputASCII::open(string file)
{
  fout.open(file.c_str()); 
  return;
}

void OutputASCII::close(){
  if (fout.is_open()) { fout.close(); }
  return;
}


//----------------------------------------------------------------------------
// write record

void OutputASCII::writeRecord(int islice, bool isparallel, int rank, int size, int nslice)
{
 
  double *value=new double[10*HARMMAX+10];
  int valsize=10*HARMMAX+10;
  int tag=1;
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
         this->writeRecordHeader(j,tmp);  // write header : labels, current/scan parameter

         for (int i=0; i<nstepz+1;i+=inputcom_.iphsty){
           if (j==islice){
             int ir=i/inputcom_.iphsty+1;
             f2cgetdiagnorecord_(&ir,value); // get current data record
	   } else {

	     MPI::COMM_WORLD.Recv(value, valsize, MPI::DOUBLE, j-islice, tag, status);
           }
           // write out by column
           if (flag.power[0])  { fout << scientific << setw(14) << setprecision(4) << value[0];}
           if (flag.inc[0])    { fout << scientific << setw(14) << setprecision(4) << value[1];}
           if (flag.signal[0]) { fout << scientific << setw(14) << setprecision(4) << value[2];}
           if (flag.phase[0])  { fout << scientific << setw(14) << setprecision(4) << value[3];}
           if (flag.size[0])   { fout << scientific << setw(14) << setprecision(4) << value[4];}
           if (flag.diver[0])  { fout << scientific << setw(14) << setprecision(4) << value[5];}
           if (flag.energy)    { fout << scientific << setw(14) << setprecision(4) << value[6];}
           if (flag.bunch[0])  { fout << scientific << setw(14) << setprecision(4) << value[7];}
           if (flag.xrms)      { fout << scientific << setw(14) << setprecision(4) << value[8];}
           if (flag.yrms)      { fout << scientific << setw(14) << setprecision(4) << value[9];}
           if (flag.error)     { fout << scientific << setw(14) << setprecision(4) << value[10];}
           if (flag.xpos)      { fout << scientific << setw(14) << setprecision(4) << value[11];}
           if (flag.ypos)      { fout << scientific << setw(14) << setprecision(4) << value[12];}
           if (flag.espread)   { fout << scientific << setw(14) << setprecision(4) << value[13];}
           if (flag.farfield[0]) { fout << scientific << setw(14) << setprecision(4) << value[14];}
           for (int i = 1; i<HARMMAX;i++){
              if (flag.bunch[i])  { fout << scientific << setw(14) << setprecision(4) << value[16+(i-1)*5];}
              if (flag.power[i])  { fout << scientific << setw(14) << setprecision(4) << value[17+(i-1)*5];}
              if (flag.phase[i])  { fout << scientific << setw(14) << setprecision(4) << value[20+(i-1)*5];}
              if (flag.signal[i]) { fout << scientific << setw(14) << setprecision(4) << value[19+(i-1)*5];}
           }
           fout << endl;
        }
      }
    }
  }
  

  delete [] value;
  return;
}


void OutputASCII::writeRecordHeader(int islice, double value)
{
  string valuetag;
  // check for slice parameter to be written
  if (isscan){
    valuetag="scan value";
  }
  else{
    valuetag="current";
  }
  fout << endl << "********** output: slice " << setw(5) << islice << endl;
  fout << "           =================" << endl;
  fout << " " << scientific << setw(14) << setprecision(4) << value << " " << valuetag << endl << endl << endl;

  if (flag.power[0])  { fout << "    " << setfill(' ') << setw(10) << right << "power";}
  if (flag.inc[0])    { fout << "    " << setfill(' ') << setw(10) << right << "increment";}
  if (flag.signal[0]) { fout << "    " << setfill(' ') << setw(10) << right << "p_mid";}
  if (flag.phase[0])  { fout << "    " << setfill(' ') << setw(10) << right << "phi_mid";}
  if (flag.size[0])   { fout << "    " << setfill(' ') << setw(10) << right << "r_size";}
  if (flag.diver[0])  { fout << "    " << setfill(' ') << setw(10) << right << "angle";}
  if (flag.energy)    { fout << "    " << setfill(' ') << setw(10) << right << "energy";}
  if (flag.bunch[0])  { fout << "    " << setfill(' ') << setw(10) << right << "bunching";}
  if (flag.xrms)      { fout << "    " << setfill(' ') << setw(10) << right << "xrms";}
  if (flag.yrms)      { fout << "    " << setfill(' ') << setw(10) << right << "yrms";}
  if (flag.error)     { fout << "    " << setfill(' ') << setw(10) << right << "error";}
  if (flag.xpos)      { fout << "    " << setfill(' ') << setw(10) << right << "<x>";}
  if (flag.ypos)      { fout << "    " << setfill(' ') << setw(10) << right << "<y>";}
  if (flag.espread)   { fout << "    " << setfill(' ') << setw(10) << right << "e-spread";}
  if (flag.farfield[0]){ fout << "    " << setfill(' ') << setw(10) << right << "far_field";}
  if (flag.bunch[1])  { fout << "    " << setfill(' ') << setw(10) << right << "2nd_bunching";}
  if (flag.power[1])  { fout << "    " << setfill(' ') << setw(10) << right << "2nd_power";}
  if (flag.phase[1])  { fout << "    " << setfill(' ') << setw(10) << right << "2nd_phase";}
  if (flag.signal[1]) { fout << "    " << setfill(' ') << setw(10) << right << "2nd_p-mid";}
  if (flag.bunch[2])  { fout << "    " << setfill(' ') << setw(10) << right << "3rd_bunching";}
  if (flag.power[2])  { fout << "    " << setfill(' ') << setw(10) << right << "3rd_power";}
  if (flag.phase[2])  { fout << "    " << setfill(' ') << setw(10) << right << "3rd_phase";}
  if (flag.signal[2]) { fout << "    " << setfill(' ') << setw(10) << right << "3rd_p-mid";}
  for (int i=3;i<HARMMAX;i++){
    if (flag.bunch[i])  { fout << "    " << setfill(' ') << setw(10) << right << i+1 << "th_bunching";}
    if (flag.power[i])  { fout << "    " << setfill(' ') << setw(10) << right << i+1 << "th_power";}
    if (flag.phase[i])  { fout << "    " << setfill(' ') << setw(10) << right << i+1 << "th_phase";}
    if (flag.signal[i]) { fout << "    " << setfill(' ') << setw(10) << right << i+1 << "th_p-mid";}
  }
  fout << endl;
  return;
}

//---------------------------------------------------------------------------------------
// some final processing/output of the output file - merging in parallel operation

void OutputASCII::finalize(int rank, bool isparallel)
{
  if (rank==0) {
    this->close();   // close all open files
  }
  return;
}


//--------------------------------------------------------------------------
// global + input are writing the header for the main output file. input can also be used to create a template file.


void OutputASCII::init(string outputfile, int rank)
{
  // get some local information out of input parameters

  filename=outputfile;
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
  for (int i=1;i<HARMMAX;i++){  
    flag.power[i]=(inputcom_.lout[i+16]!=0);
    flag.signal[i]=(inputcom_.lout[i+16]!=0);
    flag.phase[i]=(inputcom_.lout[i+16]!=0);
    flag.bunch[i]=(inputcom_.lout[i+16]!=0);
  }
  isscan=false;
  if ((inputcom_.iscan > 0)  && ( inputcom_.iscan < 23)){ isscan=true; } // check for scan feature
  nstepz=f2cgetnstepz_();                                       // get integration steps
  firstout=f2cgetnslp_()*(1-inputcom_.iotail)*inputcom_.itdp;   // calculate the first output

  if (rank==0) { this->global(); } // write header of output file
  return;
}

//-------------------------------------------------------------------
// write header of genesis main output file

void OutputASCII::global()
{

  this->open(filename);

  // header
  fout << "--------------------------------------------------------" << endl;
  fout << "Genesis 1.3 output started..." << endl;
  fout << "Version: " << GENESIS_VERSION << endl;
  fout << "Platform: " <<endl;
  //  fout << "Input filename: " << inputcom_.inputfile << endl;
  fout << "Start: " << endl;

  // input parameter
  this->input();

  for (int i=15; i< 15+HARMMAX; i++){ // version 1.0 format: Higher harmonics are indicated with a 4, indicating 4 outputs per harmonics
    inputcom_.lout[i]*=4;
  }

  // global parameter
  fout << "flags for output parameter" << endl;
  for (int i=0;i<15+HARMMAX;i++){
    fout << " " << inputcom_.lout[i];
  }
  fout << endl;

  // output of record information
  int itmp=0;
  if (inputcom_.iphsty > 0){ itmp=nstepz/inputcom_.iphsty +1; }// note integer division  
  fout << setw(5) << itmp << " entries per record" << endl;  // output in z
  
  itmp=inputcom_.nslice-firstout;
  itmp=itmp/inputcom_.ishsty;
  fout << setw(5) << itmp << " history records " << endl;  // output in t (slices)

  // reference wavelength
  fout << scientific << setw(14) << setprecision(4) << inputcom_.xlamds << " wavelength" << endl;
  double ztmp=inputcom_.xlamds*inputcom_.zsep*inputcom_.ishsty;
  fout << scientific << setw(14) << setprecision(4) << ztmp; 
  fout << " separation of output slices" << endl;   // separation of slices

  fout << setw(5) << inputcom_.ncar << " number of gridpoints" << endl;   // field grid size
  ztmp=f2cgetdxy_()/f2cgetxkw0_();
  fout << scientific << setw(14) << setprecision(4) << ztmp; 
  fout << " meshsize " << endl;  // grid spacing
  fout << setw(7) << inputcom_.npart << " number of particles" << endl; // number particles

  if (inputcom_.ippart > 0){
     fout << setw(5) << nstepz/inputcom_.ippart+1 << " particle: records in z" << endl; 
     fout << setw(5) << itmp/inputcom_.ispart << " particle: records in t" << endl;     
  }
  else{
     fout << "    0 particle: records in z" << endl; 
     fout << "    0 particle: records in t" << endl;     
  }
  if (inputcom_.ipradi > 0){
     fout << setw(5) << nstepz/inputcom_.ipradi+1 << " field: records in z" << endl; 
     fout << setw(5) << itmp/inputcom_.isradi << " field: records in t" << endl;     
  }
  else{
     fout << "    0 field: records in z" << endl; 
     fout << "    0 field: records in t" << endl;     
  }

  //
  // magnetic lattice 
  double *aw = new double [NZMAX +1];
  double *qf = new double [NZMAX +1];
  double dz=f2cgetmaglattice_(aw,qf);  
  fout << "    z[m]      ";
  fout << "    aw        ";
  fout << "    qfld [T/m]"<< endl;
  for (int i=0;i<nstepz+1;i+=inputcom_.iphsty){
    fout << scientific << setw(14) << setprecision(4) << dz*i; 
    fout << scientific << setw(14) << setprecision(4) << aw[i];
    fout << scientific << setw(14) << setprecision(4) << qf[i] << endl;
  }
  delete [] aw;
  delete [] qf;
 
  return;
}

//----------------------------------------------------------------------------------
// repead of main input file

void OutputASCII::input()
{
  fout << " $newrun" << endl;
  fout << " aw0      = " << inputcom_.aw0 << endl;
  fout << " xkx      = " << inputcom_.xkx << endl;
  fout << " xky      = " << inputcom_.xky << endl;
  fout << " wcoefz   = " << inputcom_.wcoefz[0] << " " << inputcom_.wcoefz[1] ;
  fout << " " << inputcom_.wcoefz[2] << endl;
  fout << " xlamd    = " << inputcom_.xlamd << endl;
  fout << " fbess0   = " << inputcom_.fbess0 << endl;
  fout << " delaw    = " << inputcom_.delaw << endl;
  fout << " iertyp   = " << inputcom_.iertyp << endl;
  fout << " iwityp   = " << inputcom_.iwityp << endl;
  fout << " awd      = " << inputcom_.awd << endl;
  fout << " awx      = " << inputcom_.awx << endl;
  fout << " awy      = " << inputcom_.awy << endl;
  fout << " iseed    = " << inputcom_.iseed << endl;
  fout << " npart    = " << inputcom_.npart << endl;
  fout << " gamma0   = " << inputcom_.gamma0 << endl;
  fout << " delgam   = " << inputcom_.delgam << endl;
  fout << " rxbeam   = " << inputcom_.rxbeam << endl;
  fout << " rybeam   = " << inputcom_.rybeam << endl;
  fout << " alphax   = " << inputcom_.alphax << endl;
  fout << " alphay   = " << inputcom_.alphay << endl;
  fout << " emitx    = " << inputcom_.emitx << endl;
  fout << " emity    = " << inputcom_.emity << endl;
  fout << " xbeam    = " << inputcom_.xbeam << endl;
  fout << " ybeam    = " << inputcom_.ybeam << endl;
  fout << " pxbeam   = " << inputcom_.pxbeam << endl;
  fout << " pybeam   = " << inputcom_.pybeam << endl;
  fout << " conditx  = " << inputcom_.conditx << endl;
  fout << " condity  = " << inputcom_.condity << endl;
  fout << " bunch    = " << inputcom_.bunch << endl;
  fout << " bunchphase = " << inputcom_.bunchphase << endl;
  fout << " emod     = " << inputcom_.emod << endl;
  fout << " emodphase= " << inputcom_.emodphase << endl;
  fout << " xlamds   = " << inputcom_.xlamds << endl;
  fout << " prad0    = " << inputcom_.prad0 << endl;
  fout << " pradh0   = " << inputcom_.pradh0 << endl;
  fout << " zrayl    = " << inputcom_.zrayl << endl;
  fout << " zwaist   = " << inputcom_.zwaist << endl;
  fout << " ncar     = " << inputcom_.ncar << endl;
  fout << " lbc      = " << inputcom_.lbc << endl;
  fout << " rmax0    = " << inputcom_.rmax0 << endl;
  fout << " dgrid    = " << inputcom_.dgrid << endl;
  fout << " nscr     = " << inputcom_.nscr << endl;
  fout << " nscz     = " << inputcom_.nscz << endl;
  fout << " nptr     = " << inputcom_.nptr << endl;
  fout << " nwig     = " << inputcom_.nwig << endl;
  fout << " zsep     = " << inputcom_.zsep << endl;
  fout << " delz     = " << inputcom_.delz << endl;
  fout << " nsec     = " << inputcom_.nsec << endl;
  fout << " iorb     = " << inputcom_.iorb << endl;
  fout << " zstop    = " << inputcom_.zstop << endl;
  fout << " magin    = " << inputcom_.magin << endl;
  fout << " magout   = " << inputcom_.magout << endl;
  fout << " quadf    = " << inputcom_.quadf << endl;
  fout << " quadd    = " << inputcom_.quadd << endl;
  fout << " fl       = " << inputcom_.fl << endl;
  fout << " dl       = " << inputcom_.dl << endl;
  fout << " drl      = " << inputcom_.drl << endl;
  fout << " f1st     = " << inputcom_.f1st << endl;
  fout << " qfdx     = " << inputcom_.qfdx << endl;
  fout << " qfdy     = " << inputcom_.qfdy << endl;
  fout << " solen    = " << inputcom_.solen << endl;
  fout << " sl       = " << inputcom_.sl << endl;
  fout << " ildgam   = " << inputcom_.ildgam << endl;
  fout << " ildpsi   = " << inputcom_.ildpsi << endl;
  fout << " ildx     = " << inputcom_.ildx << endl;
  fout << " ildy     = " << inputcom_.ildy << endl;
  fout << " ildpx    = " << inputcom_.ildpx << endl;
  fout << " ildpy    = " << inputcom_.ildpy << endl;
  fout << " itgaus   = " << inputcom_.itgaus << endl;
  fout << " nbins    = " << inputcom_.nbins << endl;
  fout << " igamgaus = " << inputcom_.igamgaus << endl;
  fout << " inverfc  = " << inputcom_.inverfc << endl;
  fout << " ione4one = " << inputcom_.ione4one << endl;
  fout << " lout     = " ;
  for (int i =0; i<25;i++) { fout << inputcom_.lout[i] << " " ; }
  fout << endl;
  fout << " iphsty   = " << inputcom_.iphsty << endl;
  fout << " ishsty   = " << inputcom_.ishsty << endl;
  fout << " ippart   = " << inputcom_.ippart << endl;
  fout << " ispart   = " << inputcom_.ispart << endl;
  fout << " ipradi   = " << inputcom_.ipradi << endl;
  fout << " isradi   = " << inputcom_.isradi << endl;
  fout << " idump    = " << inputcom_.idump << endl;
  fout << " iotail   = " << inputcom_.iotail << endl;
  fout << " nharm    = " << inputcom_.nharm << endl;
  fout << " iallharm = " << inputcom_.iallharm << endl;
  fout << " iharmsc  = " << inputcom_.iharmsc << endl;
  fout << " curpeak  = " << inputcom_.curpeak << endl;
  fout << " curlen   = " << inputcom_.curlen << endl;
  fout << " ntail    = " << inputcom_.ntail << endl;
  fout << " nslice   = " << inputcom_.nslice << endl;
  fout << " shotnoise= " << inputcom_.shotnoise << endl;
  fout << " isntyp   = " << inputcom_.isntyp << endl;
  fout << " iall     = " << inputcom_.iall << endl;
  fout << " itdp     = " << inputcom_.itdp << endl;
  fout << " ipseed   = " << inputcom_.ipseed << endl;
  fout << " iscan    = " << inputcom_.iscan << endl;
  fout << " nscan    = " << inputcom_.nscan << endl;
  fout << " svar     = " << inputcom_.svar << endl;
  fout << " isravg   = " << inputcom_.isravg << endl;
  fout << " isrsig   = " << inputcom_.isrsig << endl;
  fout << " cuttail  = " << inputcom_.cuttail << endl;
  fout << " eloss    = " << inputcom_.eloss << endl;
  fout << " version  = " << inputcom_.version << endl;
  fout << " ndcut    = " << inputcom_.ndcut << endl;
  fout << " idmpfld  = " << inputcom_.idmpfld << endl;
  fout << " idmppar  = " << inputcom_.idmppar << endl;
  fout << " ilog     = " << inputcom_.ilog << endl;
  fout << " ffspec   = " << inputcom_.ffspec << endl;
  fout << " ibfield  = " << inputcom_.ibfield << endl;
  fout << " imagl    = " << inputcom_.imagl << endl;
  fout << " idril    = " << inputcom_.idril << endl;
  fout << " alignradf= " << inputcom_.alignradf << endl;
  fout << " offsetradf = " << inputcom_.offsetradf << endl;
  fout << " multconv = " << inputcom_.multconv << endl;
  fout << " $end" << endl;
}

