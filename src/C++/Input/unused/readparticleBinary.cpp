
#include "readparticleBinary.h"

ReadParticleBinary::ReadParticleBinary(){}
ReadParticleBinary::~ReadParticleBinary()
{
  cout << "Destructor is called" << endl;
  if (isOpen){
    cout << "closing partfile: " << filename << endl;
    fin.close();
  }
}

void ReadParticleBinary::init(string in_filename)
{
  isOpen = false;
  filename="";
  if (in_filename.size()> 1){
    cout << "opening: " << in_filename << endl;
    filename=in_filename;
    fin.open(filename.c_str(),ios::binary);
    isOpen=true;
  }
  return; 
  

}
