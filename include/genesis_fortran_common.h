#ifndef __GENESIS_FORTRAN_COMMON__
#define __GENESIS_FORTRAN_COMMON__

extern "C"{
  extern struct{
    double aw0,xkx,xky,wcoefz[3],xlamd,fbess0,delaw,awd,awx,awy;
    double gamma0,delgam,rxbeam,rybeam,alphax,alphay,emitx;
    double emity,xbeam,ybeam,pxbeam,pybeam,cuttail,curpeak;
    double conditx,condity,bunch,bunchphase,emod,emodphase;
    double xlamds,prad0,zrayl,zwaist,rmax0,zsep,delz,zstop;
    double quadf,quadd,fl,dl,drl,f1st,qfdx,qfdy,sl,solen;
    double curlen,shotnoise,svar,dgrid,eloss,version;
    double ibfield,imagl,idril,igamref,pradh0;
    double itram11,itram12,itram13;
    double itram14,itram15,itram16,itram21,itram22,itram23;
    double itram24,itram25,itram26,itram31,itram32,itram33;
    double itram34,itram35,itram36,itram41,itram42,itram43;
    double itram44,itram45,itram46,itram51,itram52,itram53;
    double itram54,itram55,itram56,itram61,itram62,itram63;
    double itram64,itram65,itram66,rmax0sc, phaseshift;
    int iseed,nwig,nsec,npart,ncar,lbc,nscr,nscz,nptr;
    int ildgam,ildpsi,ildx,ildy,ildpx,ildpy,itgaus,nbins;
    int iphsty,ishsty,ippart,ispart,ipradi,isradi,iertyp,iwityp;
    int idump,iotail,nharm,iallharm,iharmsc,magin,magout,lout[40];
    int ffspec,ntail,nslice,iall,itdp,ipseed,iscan,nscan,isntyp;
    int isravg,isrsig,iorb,ndcut,idmppar,idmpfld,ilog, igamgaus;
    int convharm,alignradf,offsetradf,multconv,trama;
    int iscrkup,inverfc,ione4one,iautophase;
  } inputcom_;
}

#endif
