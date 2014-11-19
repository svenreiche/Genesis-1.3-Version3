extern "C"{
  int readin_(char *name,int *len);
  int initio_();
  int initrun_();
  //  int outglob_();
  void mpisetstatus_(int *size, int *id);
  //int mpimerge_();
  int doscan_(int *slice);
  int dotime_(int *slice);
  int loadrad_(int *slice);
  int loadbeam_(int *slice,double *ku, int *seed);
  int readpart_(int *npart);
  int output_(int *step,int *slice, double *ku);
  //  int openoutputbinmpi_(int *slice);
  //  int closeoutputbinmpi_();
  //  int outdumpslippage_();
  //  int outhist_(int *slice);
  //  int outdump_(int *slice);
  int stepz_(int *step,double *ku);
  int swapfield_(int *slp);
  int f2cgetnslp_();
  int f2cgetnsep_();
  int f2cgetnstepz_();
  double f2cgetxkw0_();
  double f2cgetdxy_();
  double f2cgetscanvalue_();
  double f2cgetcurrent_();
  int f2cputcurrent_(double *);
  int f2cputnpart0_(int *);
  double f2cgetmaglattice_(double *,double *);
  int f2cgetdiagnorecord_(int *, double *);
  int f2cgetbeam_(double *);
  int f2cputbeam_(double *, int *,int *);
  int f2cgetfullfield_(double *, int *, int *);
  int f2cputfullfield_(double *, int *, int *);
  int f2cgetslippagefield_(double *, int *, int *);
  void last_();
  int f2cgetoutputfilename_(char *);
  int f2cgetpartfilename_(char *);
  int f2cgetfieldfilename_(char *);
  int f2cgetfieldharmfilename_(char *);
}

