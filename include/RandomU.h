/*
 *  RandomU.h
 *  Genesis
 *
 *  Created by Sven Reiche on 5/26/10.
 *  Copyright 2010 Paul Scherrer Institut. All rights reserved.
 *
 */

#include <math.h>
#include <stdlib.h>


#ifndef __GENESIS_RANDOMU__
#define __GENESIS_RANDOMU__



class RandomU{
public:
	RandomU(unsigned int = 0 );
	~RandomU();
	void set(unsigned int);
	double getElement();
        int getSeed();
private:
	int iv[32],iy;
	int iseed,iseed2;
};

inline int RandomU::getSeed(){
  return iy;
}

#endif


