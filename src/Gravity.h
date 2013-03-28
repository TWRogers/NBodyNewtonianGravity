/*
 * Gravity.h
 *
 *  Created on: 1 Dec 2010
 *      Author: trogers
 */

#ifndef GRAVITY_H_
#define GRAVITY_H_

#include <complex.h>
#include <fftw3.h>


class Gravity {
public:
	Gravity(const unsigned int Np, const unsigned int Ng, const double dt);

	void addParticle(const unsigned int i,
					 const double x,
					 const double m,
					 double p);             //Defines particle i

	void calcPotential();                   //Calculates the new potential
	void calcAcceleration();				//Calculates new accelerations
	double calcKineticEnergy();				//Calculates the total kinetic energy
	void addStep();							//Adds time step- calculates new velocities and
											//positions and updates them

	/* Access methods */
	double * getVelocity();					//Returns particle velocities
	double * getMass();						//Returns the particle masses
	double * getAcceleration();				//Returns the particle accelerations
	double * getPosition();					//Returns the particles positions

	double * getDensityField();				//Returns the density field
	double * getPotential();				//Returns the potential grid

	~Gravity();

private:

	void initDensityField();
	void calc_Phi_k();
	void normalise(double * array, unsigned int size); //Function to normalise result of FFT

	const unsigned int _Np;     //Number of particles
	const unsigned int _Ng;     //Number of grid units

	fftw_plan _backward_Phi_k;  //FFTW Plan for transforming potential into spatial domain
	fftw_plan _forward_rho; 	//FFTW plan for transforming density field into Fourier domain

	const double _dt;			//Constant time-step for Leapfrog
	unsigned int _i;			//Counter for number of steps

	double * _vPrev;			//Previous velocities of particles- for Leapfrog
	double * _xPrev;			//Previous positions of particle for Leapfrog

	unsigned int * _i_in_grid;  //Array for storing a particles index in the grid

	double * _A;  				//Accelerations of the particles
	double * _M;  				//Masses of the particles
	double * _v;  				//Momentums of the particles
	double * _x;  				//Positions of the particles
	double * _rho;				//Density field
	double * _Phi;				//Potential
	fftw_complex * _Phi_k;		//Potential in the Fourier domain
	fftw_complex * _rho_k;		//Density field in the Fourier domain

public:
	//Prints the array of doubles to a text file
	void printArray(double * const array, char * filename, const int size) const;
};

#endif /* GRAVITY_H_ */
