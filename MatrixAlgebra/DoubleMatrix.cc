/************************************************************************
 *                                                                      *
 * This subclass of object implements a double precision matrix.        *
 * The matrix is represented as a double array, with the elements       *
 * increasing by rows.                                                  *
 *                                                                      *
 * File:DoubleMatrix.cc                                                 *
 *                                                                      *
 ************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <C_Libraries/ftoc.h>
#include <AlgebraRoutines/algebra_routines.h>
#include "DoubleMatrix.h"
#include "matrix_routines.h"


#define MAX(a,b)  ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))
#ifndef YES
#define YES     1
#define NO      0
#endif

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the DoubleMatrix class.
// Input:       numberRows:             Rows in the matrix
//              numberCols:             Columns in the matrix
//              input_elements:         Input matrix elements
//                      
// Output:                              None
//
// Notes:
// ############################# Class Constructor ###############################
DoubleMatrix::DoubleMatrix(int numberRows, int numberCols, const double *input_elements)
{
//
// Initialize arrays:
//
  _matrixString                         = NULL;
  _matrix                               = NULL;
  _singularValues                       = NULL;
//
// Use assign function to initialize matrix:
//
  assign(numberRows, numberCols, input_elements);
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Copy constructor for the DoubleMatrix class.
// Input:       matrix:                 A previously allocated DoubleMatrix
//                      
// Output:                              None
//
// Notes:
// ############################# Class Constructor ###############################
DoubleMatrix::DoubleMatrix(const DoubleMatrix& inputMatrix)
{
//
// Initialize arrays, if arrays are added, check subclass, DoubleSquareMatrix
// to make sure they are initialized.
//
  _matrixString                         = NULL;
  _matrix                               = NULL;
  _singularValues                       = NULL;
//
// Use assign function to initialize matrix:
//
  assign(inputMatrix.getRows(), inputMatrix.getColumns(), inputMatrix.getMatrix());
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the DoubleMatrix class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
DoubleMatrix::~DoubleMatrix()
{
  delete [] _matrix;
  delete [] _matrixString;
  delete [] _singularValues;
  return;
}

// ############################# Public Method ###############################
// assign -- Sets the DoubleMatrix's element to the input element array
// Input:       inputMatrix:            array of elements for the DoubleMatrix
//              numberRows:             new number of rows (if 0 use old number)
//              numberCols:             new number of columns (if 0 use old number)
//                      
// Output:                              None
//
// NOTES:
// 1. The element array needs to be ordered by rows: [0,0], [0,1], [0,2] ...
// ############################# Public Method ###############################
void DoubleMatrix::assign(int numberRows, int numberCols, const double *inputMatrix)
{
  int   i;

  if(numberRows > 0)
  {
    _rows               = numberRows;
    if(numberCols > 0)
      _columns          = numberCols;
    delete [] _matrix;
    _matrix             = new double [_rows*_columns];
  }
//
// Copy elements
//
  if(inputMatrix != NULL)
  {
    for(i=0; i<_rows*_columns; i++)
      _matrix[i]        = inputMatrix[i];
  }
  return;
}

#define NO_U_MATRIX     0
#define NO_V_MATRIX     0
// ############################# Public Method ###############################
// singularValues  Returns the singular values from largest to smallest
// Input:               none
//                      
// Output:              singular values
//
// ############################# Public Method ###############################
const double *DoubleMatrix::singularValues()
{
  int   n;
  double        *u, *v;
 
  delete [] _singularValues;

  n                     = MAX(_rows, _columns);
  _singularValues       = new double[n];
  u                     = new double[n*n];
  v                     = new double[n*n];
  svd(_matrix, _rows, _columns, _singularValues, NO_U_MATRIX, u, NO_V_MATRIX, v);

  delete [] u;
  delete [] v;
  return _singularValues;
}

// ############################# Public Method ###############################
// cond  Returns the condition number of the matrix
// Input:               none
//                      
// Output:              condition number
//
// ############################# Public Method ###############################
double DoubleMatrix::cond()
{
  int           n;
  double        cond;
  const double *singular_values;

  singular_values       = singularValues();

  n                     = MIN(_rows, _columns);
  if(singular_values[n-1] > 0.)
  {
    cond                = singular_values[0]/singular_values[n-1];
  }
  else
    cond                = 0.;

  return cond;
}


// ############################# Public Method ###############################
// transpose  Performs the transpose on the matrix
// Input:               none
//                      
// Output:              the transposed matrix
//
// ############################# Public Method ###############################
DoubleMatrix& DoubleMatrix::transpose()
{
  int           i, j, old_rows, old_columns, number_columns;
  double        *temp;


//
// Allocate temp and copy matrix
//
  temp                  = new double[_rows*_columns];
  for(i=0; i<_rows*_columns; i++)
    temp[i]             = _matrix[i];

  old_columns           = _columns;
  old_rows              = _rows;
  _columns              = _rows;
  _rows                 = old_columns;
  number_columns        = _columns;
  for(i=0; i<_rows; i++)
  {
    for(j=0; j<_columns; j++)
    {
      _matrix INDEX(i,j) = temp[i+j*old_columns];
    }
  }
  delete [] temp;
  return *this;
}

#define MAX_COEFF_SIZE  64
// ############################# Public Method ###############################
// getMatrixString()  - Get a string for printing out the matrix
// Input:                       None
//          
// Output:                      The output string
//
// ############################# Public Method ###############################
const char *DoubleMatrix::getMatrixString()
{
  int           i,j,cols_lin, n;
  int           number_rows, number_columns;
  const         double  *elements;
  char          temp[MAX_COEFF_SIZE];

  if(_matrixString == NULL)
    _matrixString       = new char[MAX_MATRIX_STRING_SIZE];
  _matrixString[0]      = '\0';
  n                     = MAX_MATRIX_STRING_SIZE - 1;
//
// Get the matrix parameters:
//
  number_columns        = _columns;
  number_rows           = _rows;
  elements              = _matrix;
  cols_lin = 5;                 // 5 per line for double types

  for(i = 0 ; i < number_rows ; i++)
  {
    for(j = 0 ; j<number_columns ; j++)
    {
      if(j%cols_lin == 0)
      {                         // newline every cols_lin
        if(j == 0)              // start of row
        {
          sprintf(temp, "\nRow %d:", i);
          strncat(_matrixString, temp, n);
          n             -= strlen(temp);
        }
        else
        {
          sprintf(temp, "\n       ");
          strncat(_matrixString, temp, n);
          n             -= strlen(temp);
        }
       }
       sprintf(temp, "%g ", elements INDEX(i,j));
       strncat(_matrixString, temp, n);
       n                -= strlen(temp);
    }
  }
  sprintf(temp,"\n");
  strncat(_matrixString, temp, n);
  return _matrixString;
}

// ############################# Public Method ###############################
// Operator =  Sets the DoubleMatrix equal to the input DoubleMatrix
// Input:       y:              input DoubleMatrix
//                      
// Output:                      the result DoubleMatrix
//
// ############################# Public Method ###############################
DoubleMatrix& DoubleMatrix::operator = (const DoubleMatrix& y)
{
  
  if( this == &y )                                      // Check for x=x
    return *this;
//
// Perform the copy using the assign function.
//
  assign(y.getRows(), y.getColumns(), y.getMatrix());

  return *this;
}

// ############################# Public Method ###############################
// Operator +=  Add the input DoubleMatrix to this.
// Input:       y:      input DoubleMatrix
//                      
// Output:                              the result DoubleMatrix
//
// NOTES:
// ############################# Public Method ###############################
DoubleMatrix& DoubleMatrix::operator += (const DoubleMatrix& y)
{
  int   i, min_size;
  const double *input_elements;
//
// Find minimum matrix size, so we don't overrun
// It really doesn't make sense to add when the matrices
// aren't the same size, but this works.
//
  min_size                      = MIN( (_rows*_columns), ( y.getRows()*y.getColumns() ) );
  input_elements        = y.getMatrix();
  for(i=0; i<min_size; i++)
    _matrix[i]          += input_elements[i];

  return *this;
}

// ############################# Public Method ###############################
// Operator -=  Sub the input DoubleMatrix to this.
// Input:       y:      input DoubleMatrix
//                      
// Output:                              the result DoubleMatrix
//
// NOTES:
// ############################# Public Method ###############################
DoubleMatrix& DoubleMatrix::operator -= (const DoubleMatrix& y)
{
  int   i, min_size;
  const double *input_elements;
//
// Find minimum matrix size, so we don't overrun
// It really doesn't make sense to add when the matrices
// aren't the same size, but this works.
//
  min_size                      = MIN( (_rows*_columns), ( y.getRows()*y.getColumns() ) );
  input_elements        = y.getMatrix();
  for(i=0; i<min_size; i++)
    _matrix[i]          -= input_elements[i];

  return *this;
}

// ############################# Public Method ###############################
// Operator *=  Postmultiply this by the input DoubleMatrix.
// Input:       y:                      input DoubleMatrix
//                      
// Output:                              the result DoubleMatrix
//
// NOTES:
// 1. This translates into *this = *this*y
// ############################# Public Method ###############################
DoubleMatrix& DoubleMatrix::operator *= (const DoubleMatrix& y)
{
  int           i, j, k, number_columns;
  double        *result, *a_row_ptr;
  const         double *y_elements;
//
// Check dimensions
//
  if(y.getRows() != _columns)
  {
    fprintf(stderr, "Bad multiplicand in multMatrix:\n");
    return *this;
  }
//
// Perform the multiplication, store in result:
//
  result                = new double[_rows*y.getColumns()];
  number_columns        = y.getColumns();
  y_elements            = y.getMatrix();
  for(i=0; i<_rows; i++)
  {
    a_row_ptr = &_matrix[i*_columns];
    for(j=0; j<y.getColumns(); j++)
    {
      result INDEX(i,j) = 0.;
      for(k=0; k<_columns; k++)
        result INDEX(i,j) += a_row_ptr[k] * y_elements INDEX(k,j);
    }
  }
//
// Re-size matrix, copy result:
//
  if(_columns != y.getColumns())
  {
    _columns    = y.getColumns();
    delete [] _matrix;
    _matrix     = new double[_rows*_columns];
  }
  for(i=0; i<_rows*_columns; i++)
    _matrix[i]  = result[i];

  delete        [] result;
  return        *this;
}

// ############################# Public Method ###############################
// Operator =  Sets the DoubleMatrix equal to the input double
// Input:       y:              input double
//                      
// Output:                      the result DoubleMatrix
//
// ############################# Public Method ###############################
DoubleMatrix& DoubleMatrix::operator = (double y)
{
  _rows         = _columns      = 1;
  _matrix[0]    = y;
  return *this;
}

// ############################# Public Method ###############################
// Operator +=  Adds a double to the DoubleMatrix.
// Input:       y:      Double to add
//                      
// Output:              the result DoubleMatrix
//
// NOTES:
//
// ############################# Public Method ###############################
DoubleMatrix& DoubleMatrix::operator += (double y)
{
  for(int i=0; i<_rows*_columns; i++)
    _matrix[i]  += y;
  return *this;
}

// ############################# Public Method ###############################
// Operator -=  Subtracts a double from the DoubleMatrix.
// Input:       y:      Double to subtract
//                      
// Output:              the result DoubleMatrix
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
DoubleMatrix& DoubleMatrix::operator -= (double y)
{
  for(int i=0; i<_rows*_columns; i++)
    _matrix[i]  -= y;
  return *this;
}

// ############################# Public Method ###############################
// Operator *=  Multiplies a double by the DoubleMatrix.
// Input:       y:      Double to mult
//                      
// Output:              the result DoubleMatrix
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
DoubleMatrix& DoubleMatrix::operator *= (double y)
{
  for(int i=0; i<_rows*_columns; i++)
    _matrix[i]  *= y;
  return *this;
}

// ############################# Public Method ###############################
// Operator /=  Divides the DoubleMatrix by a double
// Input:       y:      Double to divide into
//                      
// Output:              the result DoubleMatrix
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
DoubleMatrix& DoubleMatrix::operator /= (double y)
{
  if(y != 0.)
  {
    for(int i=0; i<_rows*_columns; i++)
      _matrix[i]        /= y;
  }
  return *this;
}

/************************************* Won't compile on NT ************************
// ############################ Public Function ###################################
// << - print out the matrix
//
// Input:       s:      Open stream
//              x:      Matrix to print
//
// Output:              None
//
// NOTES:
// ############################ Public Function ###################################
ostream&  operator << (ostream& s, const DoubleMatrix& x)
{
  int           i,j,cols_lin;
  int           number_rows, number_columns;
  const         double  *elements;
//
// Get the matrix parameters:
//
  number_columns        = x.getColumns();
  number_rows           = x.getRows();
  elements                      = x.getMatrix();
  cols_lin = 5;                 // 5 per line for double types

  for(i = 0 ; i < number_rows ; i++)
  {
    for(j = 0 ; j<number_columns ; j++)
    {
      if(j%cols_lin == 0)
      {                         // newline every cols_lin
        if(j == 0)              // start of row
          s << "\nRow " << i << ":";
        else
          s << "\n       ";
       }
       s << elements INDEX(i,j) << " ";
    }
  }
  s << "\n";
  return s;
}
********************************************************************************************/

// ############################# Public Method ###############################
// Operator -  Returns the negative of the input DoubleMatrix
// Input:       x:      input DoubleMatrix
//                      
// Output:              the output DoubleMatrix
//
// NOTES
// ############################# Public Method ###############################
DoubleMatrix operator - (const DoubleMatrix& x)
{
  return (-1.*x);
}

// ############################# Public Method ###############################
// Operator +  Returns the sum of the two input DoubleMatrixs
// Input:       x:              input DoubleMatrix
//              y:              second DoubleMatrix
//                      
// Output:              the sum DoubleMatrix
//
// ############################# Public Method ###############################
DoubleMatrix operator + (const DoubleMatrix& x, const DoubleMatrix& y)
{
  DoubleMatrix  sum = x;
  sum   += y;                   //Stroustrup p. 302
  return sum;
}

// ############################# Public Method ###############################
// Operator -  Returns the difference of the two input DoubleMatrixs
// Input:       x:              input DoubleMatrix
//              y:              second DoubleMatrix
//                      
// Output:              the difference DoubleMatrix
//
// NOTES:
// 1. For modulo-2 addition, + equals -.
// ############################# Public Method ###############################
DoubleMatrix operator - (const DoubleMatrix& x, const DoubleMatrix& y)
{
  DoubleMatrix  difference = x;
  difference    -= y;                                   //Stroustrup p. 302
  return difference;
}

// ############################# Public Method ###############################
// Operator *  Returns the product of the two input DoubleMatrixs
// Input:       x:              input DoubleMatrix
//              y:              second DoubleMatrix
//                      
// Output:                      the product DoubleMatrix
//
// ############################# Public Method ###############################
DoubleMatrix operator * (const DoubleMatrix& x, const DoubleMatrix& y)
{
  DoubleMatrix  prod = x;
  prod  *= y;                   //Stroustrup p. 302
  return prod;
}

