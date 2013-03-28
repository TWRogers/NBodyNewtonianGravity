/*
 * Gravity.cpp
 *
 *  Created on: 1 Dec 2010
 *      Author: trogers
 */

#include <new>
#include <math.h>
#include "Gravity.h"
#include <iostream>
//#include "matrix.h"

using namespace std;
const double pi = acos(-1.);

Gravity::Gravity(const unsigned int Np, const unsigned int Ng, const double dt):_Np(Np),_Ng(Ng),_dt(dt),_i(0)
{
	//Assign memory for arrays as required
	_i_in_grid = new unsigned int[_Np];

	_M         = new double[_Np];
	_v         = new double[_Np];
	_vPrev     = new double[_Np];
	_xPrev     = new double[_Np];
	_x         = new double[_Np];
	_rho       = new double[_Ng];
	_Phi       = new double[_Ng];
	_A         = new double[_Np];

	_Phi_k = (fftw_complex*) fftw_malloc(sizeof(*_Phi_k)*(_Ng/2+1));
	_rho_k = (fftw_complex*) fftw_malloc(sizeof(*_rho_k)*(_Ng/2+1));

	_forward_rho = fftw_plan_dft_r2c_1d(_Ng, _rho, _rho_k, FFTW_MEASURE);//Create plan
	_backward_Phi_k = fftw_plan_dft_c2r_1d(_Ng, _Phi_k, _Phi, FFTW_MEASURE);
}

void Gravity::addParticle(const unsigned int i, const double x, const double m, const double v)
{
	if(i < _Np) //Checks that requested particle is within range
	{
		//Checks whether mass needs replacing or adding
		if(!_M[i]) cout << "Adding new mass (" << i << ") with: M = " << m << ", v = " << v << " x = " << x << endl;
		else  cout << "Replacing mass (" << i << ") with: M = " << m << ", v = " << v << " x = " << x << endl;
		_M[i] = m;                  //Store mass in array
		_v[i] = v;					//Store velocity in array
		_x[i] = x;					//Store position in array
		cout << "DONE!" << endl;
	}
	else
	{	//Gives error if particle not allowed
		cout << "Error:::> Particle not in range 0 < i < " << _Np-1;
		return;
	}

}
void Gravity::calcPotential()
{
	initDensityField();								//Update density field
	fftw_execute(_forward_rho);						//Execute FFT
	calc_Phi_k();									//Calculate Phi in Fourier domian

	fftw_execute(_backward_Phi_k);					//Transform back into spatial domain
	normalise(_Phi,_Ng);							//Normalise the result
}
void Gravity::calcAcceleration()
{
	for(unsigned int i = 0; i < _Np; i++)//Loops over particles
	{
		if(_i_in_grid[i] == 0)          _A[i] = -0.5*_Ng*(_Phi[1]-_Phi[_Ng-1]);//Takes care of periodic BCs
		else if(_i_in_grid[i] == _Ng-1) _A[i] = -0.5*_Ng*(_Phi[0]-_Phi[_Ng-2]);
		else                            _A[i] = -0.5*_Ng*(_Phi[_i_in_grid[i]+1]-_Phi[_i_in_grid[i]-1]);
	}									//Calculates acceleration by central approximation
}
double Gravity::calcKineticEnergy()
{
	double kineticEnergy = 0.;
	for(unsigned int i = 0; i < 1; i++)      //Loops over the particles
	{
		kineticEnergy += 0.5*_v[i]*_v[i]*_M[i];//Adds the kinetic energy of ith mass
	//	kineticEnergy += _M[i]*_Phi[i];
	}
	return kineticEnergy;
}
void Gravity::addStep()
{

	for(unsigned int i = 0; i < _Np; i++)
	{
		if(_i == 0)				      //Checks for first step
		{
			_xPrev[i] = _x[i];		  //Stores the first "previous x"
			_vPrev[i] = _v[i];		  //Stores the first "previous v"
			_x[i] += _v[i]*_dt;		  //Performs initial Euler step for x
			_v[i] += _A[i]*_dt;       //Performs initial Euler step for v
		}
		else
		{
			double x_temp = _x[i];    //Temporary store of x before update
			double v_temp = _v[i];	  //Temporary store of v before update

			_x[i] = _x[i] + _v[i]*_dt;//Calculates next x via Leapfrog
			_v[i] = _v[i] + _A[i]*_dt;//Calculates next v via Leapfrog


			_xPrev[i] = x_temp;       //Sets x n-1 value for next step
			_vPrev[i] = v_temp;       //Sets x n-1 value for next step

		}
	}
	calcPotential();
	calcAcceleration();
	_i++;
}

Gravity::~Gravity()
{
	//Looks after dynamically assigned memory
	delete [] _A;		  _A         = NULL;
	delete [] _M;		  _M         = NULL;
	delete [] _v;		  _v         = NULL;
	delete [] _x;		  _x         = NULL;
	delete [] _i_in_grid; _i_in_grid = NULL;
	delete [] _rho;		  _rho       = NULL;
	delete [] _Phi;       _Phi       = NULL;
	delete [] _vPrev;	  _vPrev	 = NULL;
	delete [] _xPrev;	  _xPrev 	 = NULL;

	fftw_free(_rho_k);	  _rho_k     = NULL;
	fftw_free(_Phi_k);	  _Phi_k     = NULL;

	fftw_destroy_plan(_backward_Phi_k);				//Destroy plans
	fftw_destroy_plan(_forward_rho);
}

double * Gravity::getVelocity()
{
	return _v;
}
double * Gravity::getDensityField()
{
	return _rho;
}
double * Gravity::getAcceleration()
{
	return _A;
}
double * Gravity::getPosition()
{
	return _x;
}
double * Gravity::getPotential()
{
	return _Phi;
}

void Gravity::normalise(double * array, unsigned int size)
{
	for(unsigned int i = 0; i < size; i++)
	{
		array[i] = array[i]/size; //Divides each element by size to normalise
	}
}


void Gravity::calc_Phi_k()
{
	//_Phi_k[0] = 0;
	for(unsigned int k = 1; k < _Ng/2+1; k++) //Loops over k avoiding infinity at k=0
	{
		_Phi_k[k] = _rho_k[k]/(2.*_Ng*_Ng*(cos(2*pi*k/_Ng)-1.)); //Calculate k'th value of Phi
	}
}
void Gravity::initDensityField()
{
	double h = 1./_Ng;                        //Defines cell separation
	unsigned int grid_index;				  //Index in density field to be calculated
	for(unsigned int i = 0; i < _Ng; i++)
	{
		_rho[i] = 0.;                         //Initialise density grid to 0
	}
	for(unsigned int i = 0; i < _Np; i++)
	{
		double diff = _x[i]/h-static_cast<int>(_x[i]/h); //Find difference between actual value and rounded down
		if (diff < 0.5)                             //Check whether the position should be rounded down
		{
			grid_index = static_cast<int>(_x[i]/h);
		}
		else
		{
			grid_index = static_cast<int>(_x[i]/h+1);
		}
		if(grid_index > _Ng-1)
		{
			grid_index = 0;//If the new grid_index is out of grid, it loops round
			_x[i] = 0.;                   //Resets x so that it stays in range 0 < x < L
		}
		else if (grid_index <= 0.)
		{
			grid_index = _Ng -1;
			_x[i] = 1.;
		}
		_rho[grid_index] += _M[i]/h;
		_i_in_grid[i] = grid_index;       //Stores grid index
	}
}

void Gravity::printArray(double * const array, char* filename, const int size) const
{
	FILE *outFile;
	outFile = fopen(filename,"wt");
	for(int i = 0; i < size; i++)
	{
		fprintf(outFile, "%f\t%f\n", (double)i/size, array[i]);
	}
	fclose(outFile);
	outFile = NULL;
}


