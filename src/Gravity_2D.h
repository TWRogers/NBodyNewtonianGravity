/*
 * Gravity_2D.h
 *
 *  Created on: 8 Dec 2010
 *      Author: trogers
 */

#ifndef GRAVITY_2D_H_
#define GRAVITY_2D_H_

#include "Gravity.h"

class Gravity_2D
{
public:
	Gravity_2D(const unsigned int Np, const unsigned long int Ng, double dt);
	void addParticle(const unsigned int i, const double x, const double y, const double m, double vx, double vy);
	void calcPotential();		//Calculates the potential from density field
	void calcAcceleration();	//Calculates acceleration components
	void addStep();			    //Adds time step- calculates new velocities and
	                            //positions and updates them

	//Access methods
	double * getx();
	double * gety();
	double * getvx();
	double * getvy();
	double * getPotential();
	double * getDensityField();
	double * getAccelerationy();
	double * getAccelerationx();
	double * getTrace();

	//Prints double array in text file in format required for gnuplot splot function
	void printDoubleArray(double * const array, char* filename, const int size) const;
	~Gravity_2D();
private:

	fftw_plan _backward_Phi_k;  //FFTW Plan for transforming potential into spatial domain
	fftw_plan _forward_rho; 	//FFTW plan for transforming density field into Fourier domain

	void calc_Phi_k();
	void initDensityField();
	void normalise(double * array, unsigned long int size);

	const unsigned int _Np;     //Number of particles
	const unsigned int _Ng;     //Number of grid units
	const double _dt;
	unsigned int _i;

	double*    _trace;

	double*	_M;                 //Masses of the particles
	double*	_x;                 //x position of the particles
	double*	_y;                 //y position of the particles
	double*	_rho;               //Density field, row major format
	double*	_Phi;				//Potential, row major format
	double*	_Ax;				//Acceleration in x of the particles
	double*	_Ay;				//Acceleration in y of the particles
	double*	_vx;				//Velocity in x of the particles
	double*	_vy;				//Velocity in y of the particles
	double*	_vxPrev;			//Previous velocity in x for Leapfrog in x
	double*	_vyPrev;			//Previous velocity in y for Leapfrog in y
	double*	_xPrev;				//Previous position in x for Leapfrog in x
	double*	_yPrev;				//Previous position in y for Leapfrog in y

	unsigned int* _i_in_grid;	//Index in x of particles
	unsigned int* _j_in_grid;	//Index in y of particles

	fftw_complex* _Phi_k;		//Potential in the Fourier domain, row major format
	fftw_complex* _rho_k;		//Density field in the Fourier domain, row major format
};

#endif /* GRAVITY_2D_H_ */
