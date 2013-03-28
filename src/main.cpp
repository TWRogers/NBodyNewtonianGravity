//============================================================================
// Name        : gravity.cpp
// Author      : Thomas Rogers
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================


#include "Gravity.h"
#include "Gravity_2D.h"
#include "cpbitmap.h"
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>

int main() {
	/* Please look at Gravity.cpp and Gravity_2D.cpp for source code*/
    /* 2-D CODE  */

	unsigned long int Ng = 512;
	int Np = 51;
	double dt = 0.01;		   //set time step
	srand ( time(NULL) );      //Set random seed
	Gravity_2D b(Np,Ng,dt);
//	b.addParticle(0,0.50,0.50,0.1,0.,0.);
//	b.addParticle(1,0.49,0.49,0.000001,-0.0892,0.0892);
//	b.addStep();
//	b.calcAcceleration();
//	std::cout << b.getAccelerationx()[0] << "\t" << b.getAccelerationx()[1] << std::endl;
//	std::cout << b.getAccelerationy()[0] << "\t" << b.getAccelerationy()[1] << std::endl;
//	b.printDoubleArray(b.getPotential(),(char*)"2dPot.txt",Ng);
//	b.addParticle(1,0.6,0.5,0.01,0.,0.);
//	b.addParticle(2,0.44,0.57,0.03,0.,0.);
//	b.addParticle(0,0.5,0.5,0.1,0.,0.);

	for(int i = 0; i < Np-1; i++)
	{
		double x,y,m;//,vx,vy;
		x = .5*(double)rand()/RAND_MAX+0.25;
		y = .5*(double)rand()/RAND_MAX+0.25;
		m = 0.0001*(double)rand()/RAND_MAX;
		//vx = 0.1*(x-0.5);
		//vy = 0.1*(y-0.5);
		b.addParticle(i+1,x,y,m,0.,0.);
	}

//	std::ofstream position2D;
//	position2D.open("position2D1024.txt");
	int iMax =2000;
	for(int i = 0; i < iMax; i++)     //loop over time
	{
		b.addStep();
//      position2D << sqrt(pow((b.getx()[0]-0.5),2)+pow((b.gety()[0]-0.5),2)) << "\t" << sqrt(pow((b.getx()[1]-0.5),2)+pow((b.gety()[1]-0.5),2)) << std::endl;
//		position2D << i*dt << "\t" << b.getx()[0] << "\t" << b.getx()[1] << "\t" << b.getvx()[0] << "\t" << b.getvx()[1] << "\t" << b.getAccelerationx()[0] << "\t" << b.getAccelerationx()[1]
//		<< "\t" << b.gety()[0] << "\t" << b.gety()[1] << "\t" << b.getvy()[0] << "\t" << b.getvy()[1] << "\t" << b.getAccelerationy()[0] << "\t" << b.getAccelerationy()[1] << std::endl;

		if (i%20 == 0)				   //print out bitmaps
		{
			std::ostringstream fileName;
			fileName << "./t_" << i*dt << ".bmp";
			writebmp(b.getPotential(), Ng, Ng,(char*)fileName.str().c_str());
		}

	}
//	writebmp(b.getTrace(), Ng, Ng,(char*)"trace.bmp"); //Print trace .bmp

//std::cout << b.getAccelerationx()[0] << std::endl;
//std::cout << b.getAccelerationy()[0] << std::endl;
//std::cout << b.getAccelerationx()[1] << std::endl;
//std::cout << b.getAccelerationy()[1] << std::endl;

//writebmp(b.getPotential(), Ng, Ng, (char*)"2dPot");
//writebmp(b.getDensityField(), Ng, Ng, (char*)"2dDen");

//b.printDoubleArray(b.getDensityField(),(char*)"2d.txt",Ng);
//b.printDoubleArray(b.getPotential(),(char*)"2dPot.txt",Ng);




    /* 1-D CODE  */

//	int Ng = 16;
//	int Np = 2;
//	double dt = 0.01;
//	Gravity a(Np,Ng,dt);
//	a.addParticle(0,0.4,0.01,0.);
//	a.addParticle(1,0.6,0.01,0.);
//	std::ofstream dataOut;
//	dataOut.open("16n3.txt");
//
//for(int i = 0; i < 10000; i++)
//{
//	a.addStep();
//	dataOut << i*dt << "\t";
//	//dataOut << a.calcKineticEnergy() << "\t";
//	dataOut << a.getPosition()[0] << "\t" << a.getVelocity()[0] << "\t" << a.getAcceleration()[0] << "\t";
//	dataOut << a.getPosition()[1] << "\t" << a.getVelocity()[1] << "\t" << a.getAcceleration()[1] << std::endl;
//}


//	a.addStep();

//  std::cout << a.getAcceleration()[0] << "\t" << a.getAcceleration()[1] << std::endl;
//  a.printArray(a.getPotential(),(char*)"potential1d2p.txt",Ng);

//  a.printDoubleArray(a.getDensityField(),(char*)"rho.txt",Ng);
//  a.printArray(a.getPotential(),(char*)"./equilibrium.txt",Ng);

	return 0;
}
