#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

#include <stdio.h>
#include <cstring>
#include <ctime>

#include "mpi.h"

#ifdef MPE
#include "mpe.h"
#endif

// genesis headerfiles & classes

#include "readparticleHDF5.h"
#include "readfieldHDF5.h"
#include "outdumpHDF5.h"
#include "outputHDF5.h"
#include "RandomU.h"

// interface with fortran

#include "genesis_fortran_common.h"
#include "genesis_fortran_calls.h"

using namespace std;

int main (int argc, char *argv[]) {


        //-------------------------------------------------------
        // init MPI and get size etc.
        //

        MPI::Status status; //MPI
        MPI::Init(argc, argv); //MPI

        int size=1;
        int rank=0;

        size=MPI::COMM_WORLD.Get_size(); // get size of cluster
        rank=MPI::COMM_WORLD.Get_rank(); // assign rank to node

        bool isparallel=(size > 1) ; // flag for parallel support
        mpisetstatus_(&size,&rank); // distribute MPI to fortran code - should be obsolete

	//-------------------------------------------------------------
        // some profiling options, includen when the macro MPE is defined in the input deck

#ifdef MPE
        MPE_Init_log();

        int event1a = MPE_Log_get_event_number(); 
        int event1b = MPE_Log_get_event_number(); 
        int event2a = MPE_Log_get_event_number(); 
        int event2b = MPE_Log_get_event_number(); 
        int event3a = MPE_Log_get_event_number(); 
        int event3b = MPE_Log_get_event_number(); 
        int event4a = MPE_Log_get_event_number(); 
        int event4b = MPE_Log_get_event_number(); 

        if (rank == 0) {
	   MPE_Describe_state(event1a, event1b, "Computation", "red");
	   MPE_Describe_state(event2a, event2b, "Loading Field",   "blue");
	   MPE_Describe_state(event3a, event3b, "Swapping Field",    "green");
	   MPE_Describe_state(event4a, event4b, "Output",   "orange");
        }

#endif
        

        //------------------------------------------------------------
        // Parsing Input Argument and Input file


        if (argc==1) {
          return 0; 
        } // create template

        int nlen=strlen(argv[argc-1]);

        int err=readin_(argv[argc-1],&nlen);
	if (err < 0){
          return err;
        }
        err=initio_();

        time_t timer;
	if (rank==0) {
          time(&timer);
	  cout << "Starting Time: " << ctime(&timer) << endl;
        }

        

        //--------------------------------------------------------------
        // Parse magnetic file + create magnetic lattice

        err=initrun_();
       
 	//----------------------------------------------------------------
        // some fortran-C++ crap


        int nslice=inputcom_.nslice;
        int ntail=inputcom_.ntail;
        int itdp=inputcom_.itdp;
        double ku=f2cgetxkw0_();
        int nstepz=f2cgetnstepz_();
	int nsep=f2cgetnsep_();
	int nslp=f2cgetnslp_();
        int nharm=1;
        if (inputcom_.nharm>1){
           if (inputcom_.iallharm==0){
	      nharm=2;
            }
            else {
	      nharm=(inputcom_.nharm);
            }
	 }

        //--------------------------------------------------------------
        // Open outputfiles etc
        //
	//        err=outglob_(); // old fortran routine
	//        string outputfilename=inputcom_.outputfile;
        char file[30];
	int len=f2cgetoutputfilename_(file);
        file[len]='\0';
        string outputfilename=file;

	len=f2cgetpartfilename_(file);
        file[len]='\0';
        string partfilename=file;

        len=f2cgetfieldfilename_(file);
        file[len]='\0';
        string fieldfilename=file;  

        len=f2cgetfieldharmfilename_(file);
        file[len]='\0';
        string fieldharmfilename=file;  


        //       OutputBase  *outputasc= new OutputASCII;
	//       OutdumpBase *outdump = new OutdumpBinary;
        OutputBase  *output= new OutputHDF5;
        OutdumpBase *outdump = new OutdumpHDF5;
        ReadParticleBase *inpart = new ReadParticleHDF5;
        ReadFieldBase *infield  = new ReadFieldHDF5;
        ReadFieldBase *infieldh = new ReadFieldHDF5;


	output->init(outputfilename,rank);
        outdump->init(outputfilename, rank);

        inpart->init(partfilename);
	infield->init(fieldfilename,1);
	infieldh->init(fieldharmfilename,nharm);
        if (rank==0){
	  infield->readslippage();    
	  infieldh->readslippage();    
        }

        double *value=new double[10*HARMMAX+10];
        int valsize=10*HARMMAX+10;


        //-----------------------------------
        // generate unique seeds

        unsigned int ipseed=abs(inputcom_.ipseed);
        RandomU ran2;
        ran2.set(ipseed); 
        vector<int> seeds;
        int seedloc;
        if (inputcom_.iall==0){
          for (int i=0;i<nslice;i++){
	    ran2.getElement();
            seedloc=ran2.getSeed();
            seeds.push_back(seedloc);
          }
        }else{
          for (int i=0;i<nslice;i++){
            seeds.push_back(ipseed);
          }
	}
  
        //------------------------------------------
        // allocate some working space
 
        int worksize=0;
        double *workf, *workg;
        if (itdp!=0){      // allocate memory for MPI communication (radiation field)
          worksize=16*inputcom_.ncar*inputcom_.ncar*nharm;
	  workf=new double[worksize];
          workg=new double[worksize];
          worksize=worksize/8;   // number of field elements in record
	}


	//----------------------------------------------------------------
        // main loop

        for(int islice=1+rank; islice <= nslice; islice+=size){

#ifdef MPE 
	  MPE_Log_event(event2a, 0, "start loading");
#endif

	  int mpi_loop=size; 
          if ((islice-rank) >= (nslice-size+1)) { mpi_loop= (nslice % size) ; }
          if (mpi_loop == 0) { mpi_loop=size; }
 
          //  beginning of undulatorstart new slice, reset position in undulator    
          int istepz=0;

          //---------------------------------------
          // load field and beam for given slice
   
          seedloc=-abs(seeds.at(islice-1));
    
	  err=doscan_(&islice);
          err=dotime_(&islice);

          err=loadrad_(&islice);
          err=loadbeam_(&islice,&ku,&seedloc);
          // overwrite internal generation
          inpart->readpart(islice+ntail);  
	  infield->readfield(islice+ntail);  // has to be after loadrad
	  infieldh->readfield(islice+ntail);  // has to be after loadrad
           
	  //          cout << "Current Slice: " << islice << endl;
          err=output_(&istepz,&islice,&ku); // diagnostic & status

#ifdef MPE 
	  MPE_Log_event(event2b, 0, "end loading");
#endif
#ifdef MPE 
	  MPE_Log_event(event4a, 0, "start output");
#endif
          outdump->writeRecord(istepz, islice,isparallel,rank, size, nslice);  
#ifdef MPE 
	  MPE_Log_event(event4b, 0, "end output");
#endif
          
	  //-----------------------------
          // start of integration
	  //          if (rank==0){ cout << nslp << " " << nsep << endl; }
          for (int islp=1; islp <= nslp; islp++){
            int lstepz=nsep;
            if (islp == nslp){ lstepz=nstepz-(islp-1)*nsep;}  // correct length at undulator end
            for (int isep=1;isep<=lstepz;isep++){
#ifdef MPE 
	  MPE_Log_event(event1a, 0, "start calculation");
#endif

              istepz++;
	      err=stepz_(&istepz,&ku);
	      err=output_(&istepz,&islice,&ku);

#ifdef MPE 
	  MPE_Log_event(event1b, 0, "end calculation");
#endif
#ifdef MPE 
	  MPE_Log_event(event4a, 0, "start output");
#endif
	      outdump->writeRecord(istepz, islice,isparallel,rank, size, nslice);  
#ifdef MPE 
	  MPE_Log_event(event4b, 0, "end output");
#endif


            } 
      
            //-----------------------------------------------
            // shift radiation field in time dependent simulations
            //

#ifdef MPE 
	  MPE_Log_event(event3a, 0, "start field swapping");
#endif

            if ((itdp!=0) && (mpi_loop > 1) && (islp < nslp)){
              int first=1;
	        err=f2cgetfullfield_(workf,&nharm,&first); 
		int rank_next=rank+1;
                int rank_prev=rank-1;
                int tag=1;
                if (rank_next >= mpi_loop ) { rank_next=0; }
                if (rank_prev < 0 ) { rank_prev = mpi_loop-1; }	
	        if ( (rank % 2)==0 ){
		  MPI::COMM_WORLD.Send(workf, worksize, MPI::DOUBLE, rank_next, tag);
		  MPI::COMM_WORLD.Recv(workg, worksize, MPI::DOUBLE, rank_prev, tag,status);
		}
		else {
		  MPI::COMM_WORLD.Recv(workg, worksize, MPI::DOUBLE, rank_prev, tag,status);
		  MPI::COMM_WORLD.Send(workf, worksize, MPI::DOUBLE, rank_next, tag);
                }
                err=f2cputfullfield_(workg,&nharm,&first);
	    }
            if ((itdp!=0) && (rank==0) && (islp<nslp)){
  	          err=swapfield_(&islp);   // swapping out current field with time record
	    }
#ifdef MPE 
	  MPE_Log_event(event3b, 0, "end field swapping");
#endif

          }


	  //-------------------------------------------------------------
	  // write out the results for each slice 
          //
#ifdef MPE 
	  MPE_Log_event(event4a, 0, "start output");
#endif
          output->writeRecord(islice,isparallel, rank, size, nslice);  
	  outdump->writeRecord(-1, islice,isparallel,rank, size, nslice);  
#ifdef MPE 
	  MPE_Log_event(event4b, 0, "end output");
#endif
                 
          cout << "Node "<<rank<< ": Calculation done for slice "<< islice  << endl; 

	  if (rank==0) {
	     time(&timer);
	     cout << endl << "Inner Loop " << (int((islice-1)/size)+1) << " of " << (int((nslice-1)/size)+1) << " done: " << ctime(&timer) << endl;
          }


        } // end of loop over slices

#ifdef MPE 
	  MPE_Log_event(event4a, 0, "start output");
#endif

	outdump->finalize(rank,isparallel); 
        output->finalize(rank,isparallel);
#ifdef MPE 
	  MPE_Log_event(event4b, 0, "end output");
#endif

        last_();

        
        //---------------------------------------------------------------
        // clean up allocated memory
        seeds.clear();
        if (itdp!=0){ 
            delete [] workf ;
            delete [] workg ;
        }    // release memory for swaping fields between nodes
   

        //---------------------------------------------------------------
        // MPI Termination

        cout << "Node "<<rank<< " is terminating the program"  << endl; 

 	if (rank==0) {
          time(&timer);
	  cout << "Ending Time: " << ctime(&timer) << endl;
        }

#ifdef MPE
        MPE_Finish_log(outputfilename.c_str());
#endif
        MPI::Finalize(); // node turned off

        return 0;
}
