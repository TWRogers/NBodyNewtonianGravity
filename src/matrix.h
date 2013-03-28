#include <new>
#include <vector>
#include<iostream>
#include <cassert>

/*
 * Matrix class, defines Matrix and vector operations required for this analysis
 */
class Matrix
{
public:
Matrix();
Matrix(const int rows, const int cols):m_cols(cols),m_rows(rows)
{
  m_matrix=new double * [m_rows];                                 //creates 1st dimension of 2 dimensional dynamic array
  for(int i=0; i<m_rows; i++)                                     //loops over first dimension creating dynamic arrays inside each
  *(m_matrix+i)=new double [m_cols];
  set_zero();                                                     //always intialise the matrix to zero matrix
}
Matrix(const Matrix& a):m_cols(a.m_cols),m_rows(a.m_rows)         //constructor for copying a matrix
{
  m_matrix=new double * [m_rows];                                 //takes on same dimensions
  for(int i=0; i<m_rows; i++)
  *(m_matrix+i)=new double [m_cols];
  set_elements(a);                                                //copies values
}
~Matrix()                                                         //destructor- important to correctly delete array to prevent segmentation faults
{
  for(int i=0; i<m_rows; i++)
  {
    delete [] *(m_matrix+i);
    *(m_matrix+i) = NULL;
  }
	delete [] m_matrix;
  m_matrix = NULL;
}
void set_zero()                                                    //loops through and sets all elements to zero
{
    for(int i = 0; i < m_rows; i++)
    {
      for(int j = 0; j < m_cols; j++)
      {
        m_matrix[i][j] = 0;
      }
    }
}

void operator+= (Matrix &other)                                     //overloads += so that one can increment a matrix by other
{ 
 if(other.m_rows != m_rows || other.m_cols != m_cols)               //sanity check
  {
	  std::cout << "Cannot increment matrix, check dimensions" << std::endl;
    assert(0);
  }
  else
  {
    for(int i = 0; i < m_rows; i++)
    {
      for(int j = 0; j < m_cols; j++)
      {
        m_matrix[i][j] += other.get_element(i,j);
      }
    }
  }
}


Matrix operator* (double x)                                           //overloads the * operator for multiplying matrix by doubles
{ 
    Matrix temp(m_rows,m_cols);
    for(int i = 0; i < m_rows; i++)
    {
      for(int j = 0; j < m_cols; j++)
      {
        temp.set_element(i,j,m_matrix[i][j]*x);
      }
    }
  return temp;
}
void scale(double x)                                                   //allows the return of a scaled matrix
{ 
    for(int i = 0; i < m_rows; i++)
    {
      for(int j = 0; j < m_cols; j++)
      {
        m_matrix[i][j] = m_matrix[i][j]*x;
      }
    }
}
std::vector<double> operator* (std::vector<double> vec)                 //overload * operator for multiplying a std::vector<double> by matrix
{ 
  std::vector<double> temp_vec(vec.size());
  if(vec.size() != m_cols)
  {
	 std::cout << "Cannot multiply vector by matrix, check dimensions" << std::endl; //sanity check to see whether multiplication is compatible
   assert(0);
  }
  else
  {
      std::vector<double> temp_vec(vec.size());
      for(int k = 0; k < temp_vec.size(); k++)
      {
        double temp_element = 0.;
        for(int i = 0; i < m_cols; i++)
        {
	  temp_element += m_matrix[k][i]*vec[i];
        }
        temp_vec[k] = temp_element;
      }
  return temp_vec;
  }

}

void print_row(int i)                                                       //prints matrix row i to screen
{
  for(int j = 0; j < m_cols; j++)
  {
    std::cout << m_matrix[i][j] << "\t";
  }
  std::cout << std::endl;
}
void print_matrix()                                                          //loops over rows printing each to screen
{
  for(int i = 0; i < m_rows; i++) print_row(i);  
  std::cout << std::endl;
}

void set_element(int i, int j, double x)                                     //set element in row i, col j equal to x
{            
  m_matrix[i][j] = x;
}
void set_elements(Matrix a)                                                  //copy elements from Matrix a to this matrix
{
    for(int i = 0; i < m_rows; i++)
    {
      for(int j = 0; j < m_cols; j++)
      {
        m_matrix[i][j] = a.get_element(i,j);
      }
    }
}
void set_unit()                                                                //set matrix to the unit matrix                                                       
{
   if(m_rows != m_cols)
  {
	 std::cout << "Cannot create a " << m_rows << "x" << m_cols << " matrix!" << std::endl; //sanity check to see whether multiplication is compatible
   assert(0);
  }
  else
  {
    for(int i = 0; i < m_rows; i++)
    {
      for(int j = 0; j < m_cols; j++)
      {
        if( i == j ) m_matrix[i][j] = 1;
        else  m_matrix[i][j] = 0;
      }
    }
  }
}
double get_element(int i, int j)                                                 //return value in row i collumn j
{
  return m_matrix[i][j];
}
public:
const int m_rows;                                                                //const number of rows to prevent outside tampering
const int m_cols;
private:
double ** m_matrix;
};





