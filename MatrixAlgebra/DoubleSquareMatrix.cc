/************************************************************************
 *                                                                      *
 * This subclass of object implements a double precision square matrix. *
 * The matrix is represented as a double array, with the elements       *
 * increasing by rows.                                                  *
 * The usual overloaded operators are included: +, -, *, /, etc.        *
 * This class is a subclass of the DoubleMatrix class.                  *
 *                                                                      *
 * File:DoubleSquareMatrix.cc                                           *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 07/02/04  - Extracted from DoubleMatrix.                         *
 *  2. 05/19/05  - Don't allocate luMatrix until needed.                *
 ************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <C_Libraries/ftoc.h>
#include <AlgebraRoutines/algebra_routines.h>
#include "DoubleSquareMatrix.h"
#include "matrix_routines.h"


#define MAX(a,b)  ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))
#ifndef YES
#define YES     1
#define NO      0
#endif

// ############################ Private Function ###################################
// luDecomp - Calculate the LU decomposition.
//            Taken from: C Tools for Engineers and Scientists.
//
// Input:                               None
// Output:                  YES on succes
//
// NOTES:
// 1. This routine must be called before findDeterminant(),  etc.
// 2. _luMatrix, and _pivot are modified.
// ############################ Private Function ###################################
LOGICAL DoubleSquareMatrix::luDecomp()
{

  if(_luDecompPerformed)                                                // Check if we already factored
    return YES;
  if(_oldSize < _rows*_columns)
  {
    _oldSize            = _rows*_columns;
    delete [] _luMatrix;
    delete [] _pivot;
    _luMatrix           = new double [_rows*_columns];
    _pivot              = new int[_rows];
  }
  _luDecompPerformed    = lu_decomp(_matrix, _columns, _luMatrix, _pivot);
  return _luDecompPerformed;
}

// ############################ Private Function ###################################
// resetInstanceVariables - Resets instance variables.
//
// Input:                               None
// Output:                              None
//
// NOTES:
// ############################ Private Function ###################################
void DoubleSquareMatrix::resetInstanceVariables()
{
  _matrixString                         = NULL;
  _matrix                               = NULL;
  _singularValues                       = NULL;
  _luMatrix                             = NULL;
  _pivot                                = NULL;
  _oldSize                              = 0;
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the DoubleSquareMatrix class.
// Input:       size:                   Number of rows or columns
//              input_elements:         Input matrix elements
//                      
// Output:                              None
//
// Notes:
// ############################# Class Constructor ###############################
DoubleSquareMatrix::DoubleSquareMatrix(int size, const double *input_elements)
{
//
// Initialize all arrays, base class constructor not called
//
  resetInstanceVariables();

  assign(size, size, input_elements);

  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Copy constructor for the DoubleSquareMatrix class.
// Input:       matrix:                 A previously allocated DoubleSquareMatrix
//                      
// Output:                              None
//
// Notes:
// ############################# Class Constructor ###############################
DoubleSquareMatrix::DoubleSquareMatrix(const DoubleSquareMatrix& inputMatrix)
{
//
// Initialize all arrays, base class constructor not called
//
  resetInstanceVariables();

  assign(inputMatrix.getRows(), inputMatrix.getRows(), inputMatrix.getMatrix());

  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Copy constructor for the DoubleSquareMatrix class.
// Input:       matrix:                 A previously allocated DoubleMatrix
//                      
// Output:                              None
//
// Notes:
// ############################# Class Constructor ###############################
DoubleSquareMatrix::DoubleSquareMatrix(const DoubleMatrix& inputMatrix)
{
//
// Force square matrix:
//
  _columns                              = _rows = MIN((inputMatrix.getRows()), (inputMatrix.getColumns()));
//
// Initialize all arrays, base class constructor not called
//
  resetInstanceVariables();

  assign(_rows, _rows, inputMatrix.getMatrix());

  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Sets up a diagonal square matrix
// Input:       x:              double input
//              size:           number of rows and columns, default = 1
//                      
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
DoubleSquareMatrix::DoubleSquareMatrix(double x, int size)
{
  int           i, j, k;
  double        *temp;
//
//
// Initialize arrays:
//
  resetInstanceVariables();
//
//
  temp          = new double[size*size];
  k             = 0;
  for(i=0; i<size; i++)
  {
    for(j=0; j<size; j++)
    {
      if(i==j)
        temp[k] = x;
      else
        temp[k] = 0.;
      k++;
    }
  }
  assign(size, size, temp);

  return;
}


// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the DoubleSquareMatrix class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
DoubleSquareMatrix::~DoubleSquareMatrix()
{
  delete [] _luMatrix;
  delete [] _pivot;
  return;
}

// ############################# Public Method ###############################
// assign -- Sets the DoubleSquareMatrix's element to the input element array
// Input:       inputMatrix:            array of elements for the DoubleSquareMatrix
//              rows:                   new number of rows (if 0 use old number)
//              cols:                   new number of columns (if 0 use old number)
//                                      not really used, just included to override
//                                      base class
//                      
// Output:                              None
//
// NOTES:
// 1. The element array needs to be ordered by rows: [0,0], [0,1], [0,2] ...
// ############################# Public Method ###############################
void DoubleSquareMatrix::assign(int rows, int cols, const double *inputMatrix)
{
  int   i;

  if(rows > 0)
  {
    _rows               = _columns      = rows;
    delete [] _matrix;
    _matrix             = new double [_rows*_columns];
  }
  _luDecompPerformed    = NO;
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


// ############################ Public Function ###################################
// findDeterminant - Calculate the determinant.  Taken from: C Tools for Engineers
//                   and Scientists.
//
// Input:               None
// Output               Matrix determinant
//
// ############################ Public Function ###################################
double DoubleSquareMatrix::findDeterminant()
{
  int number_columns, i, sign;
  double d;

  if(!luDecomp())
      return 0.;

  number_columns = _columns;
  sign = 0;
  d    = 1.;
  for(i=0; i<number_columns; i++)
  {
    if(_pivot[i] != i)
      sign++;
    d *= _luMatrix INDEX(i,i);
  }
  sign = sign - ((sign>>1)<<1);
  if(sign)
    d = -d;
  return d;
}

// ############################ Public Function ###################################
// invert - Inverts the matrix, using the LU decomposition.
//          Taken from: C Tools for Engineers and Scientists.
//
// Input:               None
// Output:              flag = YES on success
//                      Return pointer to inverted matrix
//
// NOTES:
// 1. This routine calls luDecomp().
// ############################ Public Function ###################################
DoubleSquareMatrix &DoubleSquareMatrix::invert(LOGICAL *flag)
{
  int           i;
  double        *inverse_matrix;


  if(!luDecomp())
  {
    flag        = NO;
    return      *this;
  }

  inverse_matrix                = new double [_rows*_columns];
  *flag                         = invert_using_LU(_luMatrix, _pivot, _columns, inverse_matrix);
//
// Copy inverse matrix to *this
//
  for(i=0; i<_rows*_columns; i++)
    _matrix[i]  = inverse_matrix[i];
//
// delete the temp arrays:
//
  delete [] inverse_matrix;

  return *this;
}

// ############################ Public Function ###################################
// linearSolve - Finds the solution to Ax=b, using the LU decomposition.
//               Taken from: C Tools for Engineers and Scientists.
//
// Input:       bVector         The b vector
// Output:      xResult         Solution vector
//
// Returns                      YES on success
//
// NOTES:
// 1. This routine calls luDecomp().
// ############################ Public Function ###################################
LOGICAL DoubleSquareMatrix::linearSolve(const DoubleMatrix &bVector, DoubleMatrix &xResult)
{
  int           i, result;
  double        *b;
  const         double *input_elements;
//
// Check for illegal size
//
  if(_rows != bVector.getRows())
  {
    return NO;
  }
//
// Perform LU decomposition
//
  if(!luDecomp())                                                       // Find the LU decomposition
  {
    return NO;
  }
//
// Allocate temp vector
//
  b                     = new double[_columns];
  input_elements        = bVector.getMatrix();
  for(i=0; i<_rows; i++)
    b[i]                = input_elements[i];
//
// Solve, using LU
//
  result = linear_solve(_luMatrix, _pivot, b, _columns);
//
// Copy solution into x:
//
  xResult.assign(_rows, 1, b);

  delete [] b;

  return        result;
}

// ############################# Public Method ###############################
// transpose  Performs the transpose on the matrix
// Input:               none
//                      
// Output:              the transposed matrix
//
// ############################# Public Method ###############################
DoubleSquareMatrix& DoubleSquareMatrix::transpose()
{
  int           i, j, number_columns;
  double        temp;

  number_columns        = _columns;

//
// In place transpose
//
  for(i=0; i<_columns; i++)
  {
    for(j=i; j<_rows; j++)
    {
      temp                      = _matrix INDEX(i,j);
      _matrix INDEX(i,j)        = _matrix INDEX(j, i);
      _matrix INDEX(j,i)        = temp;
    }
  }
  return *this;
}

// ############################# Public Method ###############################
// Operator /=  Divide this by the input DoubleSquareMatrix.
// Input:       y:      input DoubleSquareMatrix
//                      
// Output:              the result DoubleSquareMatrix
//
// NOTES:
// 1. This translates into *this = *this/y = y_inverse*(*this)
// ############################# Public Method ###############################
DoubleSquareMatrix& DoubleSquareMatrix::operator /= (const DoubleSquareMatrix& y)
{
  LOGICAL                       flag;
  DoubleSquareMatrix            temp_mat, y_copy;
//

//
  y_copy        = y;
//
// Do invert:
//
  DoubleSquareMatrix inverse_y  = y_copy.invert(&flag);
  if(flag == NO)
  {
      fprintf(stderr, "Couldn't find inverse in /=\n");
      return *this;
  }
  *this         = inverse_y * (*this);

  return *this;
}

// ############################# Public Method ###############################
// Operator -  Returns the negative of the input DoubleSquareMatrix
// Input:       x:      input DoubleSquareMatrix
//                      
// Output:              the output DoubleSquareMatrix
//
// NOTES
// ############################# Public Method ###############################
DoubleSquareMatrix operator - (const DoubleSquareMatrix& x)
{
  return (-1.*x);
}

// ############################# Public Method ###############################
// Operator +  Returns the sum of the two input DoubleSquareMatrixs
// Input:       x:              input DoubleSquareMatrix
//              y:              second DoubleSquareMatrix
//                      
// Output:              the sum DoubleSquareMatrix
//
// ############################# Public Method ###############################
DoubleSquareMatrix operator + (const DoubleSquareMatrix& x, const DoubleSquareMatrix& y)
{
  DoubleSquareMatrix    sum = x;
  sum   += y;                   //Stroustrup p. 302
  return sum;
}

// ############################# Public Method ###############################
// Operator -  Returns the difference of the two input DoubleSquareMatrixs
// Input:       x:              input DoubleSquareMatrix
//              y:              second DoubleSquareMatrix
//                      
// Output:              the difference DoubleSquareMatrix
//
// NOTES:
// 1. For modulo-2 addition, + equals -.
// ############################# Public Method ###############################
DoubleSquareMatrix operator - (const DoubleSquareMatrix& x, const DoubleSquareMatrix& y)
{
  DoubleSquareMatrix    difference = x;
  difference    -= y;                                   //Stroustrup p. 302
  return difference;
}

// ############################# Public Method ###############################
// Operator *  Returns the product of the two input DoubleSquareMatrixs
// Input:       x:              input DoubleSquareMatrix
//              y:              second DoubleSquareMatrix
//                      
// Output:                      the product DoubleSquareMatrix
//
// ############################# Public Method ###############################
DoubleSquareMatrix operator * (const DoubleSquareMatrix& x, const DoubleSquareMatrix& y)
{
  DoubleSquareMatrix    prod = x;
  prod  *= y;                   //Stroustrup p. 302
  return prod;
}

// ############################# Public Method ###############################
// Operator /  Returns the division of the two input DoubleSquareMatrixs
// Input:       x:              input DoubleSquareMatrix
//              y:              second DoubleSquareMatrix
//                      
// Output:              the result DoubleSquareMatrix
//
// ############################# Public Method ###############################
DoubleSquareMatrix operator / (const DoubleSquareMatrix& x, const DoubleSquareMatrix& y)
{
  DoubleSquareMatrix    result = x;
  result                /= y;                   //Stroustrup p. 302
  return                result;
}

