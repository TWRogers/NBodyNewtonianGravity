/*
 * Gravity_2D.cpp
 *
 *  Created on: 8 Dec 2010
 *      Author: trogers
 */

#include "Gravity_2D.h"
#include <new>
#include <math.h>
#include <iostream>
//#include "matrix.h"

using namespace std;
const double pi = acos(-1.);

Gravity_2D::Gravity_2D(const unsigned int Np, const unsigned long int Ng, double dt):_Np(Np),_Ng(Ng),_dt(dt),_i(0)
{
	//Assign memory for arrays as required
	_M         = new double[_Np];
	_x         = new double[_Np];
	_y         = new double[_Np];
	_rho       = new double[_Ng*_Ng];
	_Phi       = new double[_Ng*_Ng];
	_Ax		   = new double[_Np];
	_Ay		   = new double[_Np];
	_vx        = new double[_Np];
	_vy        = new double[_Np];
	_vxPrev    = new double[_Np];
	_vyPrev    = new double[_Np];
	_xPrev     = new double[_Np];
	_yPrev     = new double[_Np];

	_trace     = new double[_Ng*_Ng];

	_i_in_grid = new unsigned int[_Np];
	_j_in_grid = new unsigned int[_Np];

	_Phi_k = (fftw_complex*) fftw_malloc(sizeof(*_Phi_k)*(_Ng/2+1)*_Ng);
	_rho_k = (fftw_complex*) fftw_malloc(sizeof(*_rho_k)*(_Ng/2+1)*_Ng);

	_forward_rho = fftw_plan_dft_r2c_2d(_Ng, _Ng, _rho, _rho_k, FFTW_ESTIMATE);
	_backward_Phi_k = fftw_plan_dft_c2r_2d(_Ng, _Ng, _Phi_k, _Phi, FFTW_ESTIMATE);

}
double * Gravity_2D::getx()
{
	return _x;
}
double * Gravity_2D::gety()
{
	return _y;
}
double * Gravity_2D::getvx()
{
	return _vx;
}
double * Gravity_2D::getvy()
{
	return _vy;
}
double * Gravity_2D::getDensityField()
{
	return _rho;
}
void Gravity_2D::initDensityField()
{
	double h = 1./_Ng; //Calculate the size of each grid cell
	unsigned int grid_index_x, grid_index_y;
	for(unsigned int i = 0; i < _Ng; i++)
	{
		for(unsigned int j = 0; j < _Ng; j++)
		{
			_rho[i*_Ng + j] = 0.;  //Initialise density grid to 0
		}
	}

	for(unsigned int i = 0; i < _Np; i++)
	{
		double diff_x = _x[i]/h-static_cast<int>(_x[i]/h); //Find difference between actual value and rounded down
		double diff_y = _y[i]/h-static_cast<int>(_y[i]/h); //Find difference between actual value and rounded down
		//determines grid index i
		if (diff_x < 0.5) //Check whether the position should be rounded down
		{
			grid_index_x = static_cast<int>(_x[i]/h);
		}
		else
		{
			grid_index_x = static_cast<int>(_x[i]/h+1);
		}
		//determines grid index j
		if (diff_y < 0.5) //Check whether the position should be rounded down
		{
			grid_index_y = static_cast<int>(_y[i]/h);
		}
		else
		{
			grid_index_y = static_cast<int>(_y[i]/h+1);
		}

		if(grid_index_x > _Ng-1)
		{
			grid_index_x = 0;				   //If the new grid_index is out of grid, it loops round
			_x[i] = 0.;                        //Resets x so that it stays in range 0 < x < L
		}
		else if(grid_index_x < 0)
		{
			grid_index_x = _Ng-1;				   //If the new grid_index is out of grid, it loops round
			_x[i] = 1.;                        //Resets x so that it stays in range 0 < x < L
		}
		if(grid_index_y > _Ng-1)
		{
			grid_index_y = 0;				   //If the new grid_index is out of grid, it loops round
			_y[i] = 0.;                        //Resets x so that it stays in range 0 < x < L
		}
		else if(grid_index_y < 0)
		{
			grid_index_y = _Ng-1;				   //If the new grid_index is out of grid, it loops round
			_y[i] = 1.;                        //Resets x so that it stays in range 0 < x < L
		}

		_i_in_grid[i] = grid_index_x;
		_j_in_grid[i] = grid_index_y;
		if(i == 1)_trace[grid_index_x*_Ng + grid_index_y] ++;
		//cout << _i_in_grid[i] << "\t" << _j_in_grid[i] << endl;
		_rho[grid_index_x*_Ng + grid_index_y] += _M[i]/(h*h);
	}
}
double * Gravity_2D::getTrace()
{
	return _trace;
}
void Gravity_2D::calcAcceleration()
{
	//Central approximation for acceleration in x direction
	for(unsigned int i = 0; i < _Np; i++)
	{
		if(_i_in_grid[i] == 0)          _Ax[i] = -0.5*_Ng*(_Phi[_Ng+_j_in_grid[i]]-_Phi[_Ng*(_Ng-1)+_j_in_grid[i]]);
		else if(_i_in_grid[i] == _Ng-1) _Ax[i] = -0.5*_Ng*(_Phi[_j_in_grid[i]]-_Phi[_Ng*(_Ng-2)+_j_in_grid[i]]);
		else                            _Ax[i] = -0.5*_Ng*(_Phi[_Ng*(_i_in_grid[i]+1)+_j_in_grid[i]]-_Phi[_Ng*(_i_in_grid[i]-1)+_j_in_grid[i]]);
	}
	//Central approximation for acceleration in y direction
	for(unsigned int i = 0; i < _Np; i++)
	{
		if(_j_in_grid[i] == 0)          _Ay[i] = -0.5*_Ng*(_Phi[_Ng*_i_in_grid[i]+1]-_Phi[_Ng*_i_in_grid[i]+_Ng-1]);
		else if(_j_in_grid[i] == _Ng-1) _Ay[i] = -0.5*_Ng*(_Phi[_Ng*_i_in_grid[i]]-_Phi[_Ng*_i_in_grid[i]+_Ng-2]);
		else                            _Ay[i] = -0.5*_Ng*(_Phi[_Ng*_i_in_grid[i]+_j_in_grid[i]+1]-_Phi[_Ng*_i_in_grid[i]+_j_in_grid[i]-1]);
	}
}
void Gravity_2D::calcPotential()
{
	initDensityField();
	fftw_execute(_forward_rho);
	calc_Phi_k();                        //Calulate phi in Fourier domain

	fftw_execute(_backward_Phi_k);		 //Transform back into spatial domain
	normalise(_Phi,_Ng);
}
void Gravity_2D::addStep()
{

	for(unsigned int i = 0; i < _Np; i++)
	{
		if(_i == 0) //Does Euler step for first time
		{
			_xPrev[i] = _x[i];
			_yPrev[i] = _y[i];
			_vxPrev[i] = _vx[i];
			_vyPrev[i] = _vy[i];
			_x[i] += _vx[i]*_dt;
			_vx[i] += _Ax[i]*_dt;
			_y[i] += _vy[i]*_dt;
			_vy[i] += _Ay[i]*_dt;
		}
		else
		{
			double x_temp = _x[i];             //Temporarily store the N-1 values
			double vx_temp = _vx[i];
			double y_temp = _y[i];             //Temporarily store the N-1 values
			double vy_temp = _vy[i];

			_x[i] = _xPrev[i] + 2*_vx[i]*_dt;
			_vx[i] = _vxPrev[i] + 2*_Ax[i]*_dt;         //Calculate N+1 values
			_y[i] = _yPrev[i] + 2*_vy[i]*_dt;
			_vy[i] = _vyPrev[i] + 2*_Ay[i]*_dt;         //Calculate N+1 values


			_xPrev[i] = x_temp; 					  //Sets N-1 values
			_vxPrev[i] = vx_temp;
			_yPrev[i] = y_temp; 					  //Sets N-1 values
			_vyPrev[i] = vy_temp;

		}
	}
	calcPotential();
	calcAcceleration();
	_i++;
}
void Gravity_2D::calc_Phi_k()
{
	int j = 0; //Use as index to simplify
	for(unsigned int w = 0; w < _Ng; w++)
	{
		for(unsigned int k = 0; k < _Ng/2+1; k++)
		{
			if(k+w != 0) //Checks for DC component to avoid infinity
			_Phi_k[j] = _rho_k[j]/(2.*_Ng*_Ng*(cos(2.*pi*w/(_Ng))+cos(2.*pi*k/_Ng)-2.));
			j++;		 //Increments index counter
		}
	}

}
void Gravity_2D::normalise(double * array, unsigned long int size)
{
	int j = 0; //Use as index to simplify
	for(unsigned int w = 0; w < size; w++)
	{
		for(unsigned int k = 0; k < size; k++)
		{
			array[j] = array[j]/(size*size);
			j++;		 //Increments index counter
		}
	}
}
double * Gravity_2D::getPotential()
{
	return _Phi;
}
double * Gravity_2D::getAccelerationy()
{
	return _Ay;
}
double * Gravity_2D::getAccelerationx()
{
	return _Ax;
}
void Gravity_2D::addParticle(const unsigned int i, const double x, const double y, const double m, double vx, double vy)
{
	if(i < _Np) //Checks that requested particle is within range
	{
		//Checks whether mass needs replacing or adding
		if(!_M[i]) cout << "Adding new mass (" << i << ") with: M = " << m << ", x = " << x << ", y = " << y << endl;
		else  cout << "Replacing mass (" << i << ") with: M = " << m <<  ", x = " << x << ", y = " << y  << endl;
		_M[i] = m;                  //Store mass in array
		_x[i] = x;					//Store x position in array
		_y[i] = y;					//Store y position in array
		_vx[i] = vx;				//Store x velocity in array
		_vy[i] = vy;				//Store y velocity in array
		cout << "DONE!" << endl;
	}
	else
	{
		//Gives error if particle not allowed
		cout << "Error:::> Particle not in range 0 < i < " << _Np-1;
		return;
	}
}
void Gravity_2D::printDoubleArray(double * const array, char* filename, const int size) const
{
	FILE *outFile;
	outFile = fopen(filename,"wt");
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			fprintf(outFile, "%f\n", array[i*size + j]);
		}
	fprintf(outFile, "\n");
	}
	fclose(outFile);
	outFile = NULL;
}
Gravity_2D::~Gravity_2D()
{
	//Looks after dynamically assigned memory
	delete [] _M;		  _M         = NULL;
	delete [] _x;		  _x         = NULL;
	delete [] _y;		  _y         = NULL;
	delete [] _rho;		  _rho       = NULL;
	delete [] _Phi;       _Phi       = NULL;
	delete [] _Ay;		  _Ay		 = NULL;
	delete [] _Ax;        _Ax        = NULL;
	delete [] _i_in_grid; _i_in_grid = NULL;
	delete [] _j_in_grid; _j_in_grid = NULL;
	delete [] _vx;     	  _vx        = NULL;
	delete [] _vy;        _vy        = NULL;
	delete [] _vxPrev;    _vxPrev    = NULL;
	delete [] _vyPrev;    _vyPrev    = NULL;
	delete [] _xPrev;     _xPrev     = NULL;
	delete [] _yPrev;     _yPrev     = NULL;

	fftw_free(_rho_k);	  _rho_k     = NULL;
	fftw_free(_Phi_k);	  _Phi_k     = NULL;

	fftw_destroy_plan(_backward_Phi_k);
	fftw_destroy_plan(_forward_rho);
}
