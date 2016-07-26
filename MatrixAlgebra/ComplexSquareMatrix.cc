/************************************************************************
 *                                                                      *
 * This subclass of ComplexSquareMatrix implements a square Complex     *
 * matrix.                                                              *
 * The matrix is represented as a Complex array, with the elements      *
 * increasing by rows.                                                  *
 * The usual overloaded operators are included: +, -, *, /, etc.        *
 *                                                                      *
 * File:ComplexSquareMatrix.cc                                          *
 *                                                                      *
 ************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#if defined(WIN32)
#include <GNU/Complex.h>
#include <AlgebraRoutines/algebra_routines.h>
#include <C_Libraries/ftoc.h>
#else
#include "Complex.h"
#include "algebra_routines.h"
#include "ftoc.h"
#endif
#include "ComplexSquareMatrix.h"
#include "matrix_routines.h"


#define MAX(a,b)  ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))
#ifndef YES
#define YES     1
#define NO      0
#endif

// ############################ Public Function ###################################
// luDecomp - Calculate the LU decomposition.
//            Taken from: C Tools for Engineers and Scientists.
//
// Input:                               None
// Output:                  YES on succes
//
// NOTES:
// 1. This routine must be called before findDeterminant(),  etc.
// ############################ Public Function ###################################
LOGICAL ComplexSquareMatrix::luDecomp()
{
  if(_luDecompPerformed)                                                // Check if we already factored
    return YES;
  if(_oldSize < _rows*_columns)
  {
    _oldSize            = _rows*_columns;
    delete [] _luMatrix;
    delete [] _pivot;
    _luMatrix           = new Complex [_rows*_columns];
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
void ComplexSquareMatrix::resetInstanceVariables()
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
// Class Constructor -- Constructor for the ComplexSquareMatrix class.
// Input:       size:                   # of rows or columns
//              input_elements:         Input matrix elements
//                      
// Output:                                                      None
//
// Notes:
// ############################# Class Constructor ###############################
ComplexSquareMatrix::ComplexSquareMatrix(int size, const Complex *input_elements)
{
//
// Initialize all arrays, base class constructor not called
//
  resetInstanceVariables();

  assign(size, size, input_elements);
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Copy constructor for the ComplexSquareMatrix class.
// Input:       matrix:                 A previously allocated ComplexSquareMatrix
//                      
// Output:                                                      None
//
// Notes:
// ############################# Class Constructor ###############################
ComplexSquareMatrix::ComplexSquareMatrix(const ComplexSquareMatrix&  inputMatrix)
{
//
// Initialize all arrays, base class constructor not called
//
  resetInstanceVariables();

  assign(inputMatrix.getRows(), inputMatrix.getRows(), inputMatrix.getMatrix());
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Copy constructor for the ComplexSquareMatrix class.
// Input:       matrix:                 A previously allocated ComplexSquare
//                      
// Output:                              None
//
// Notes:
// ############################# Class Constructor ###############################
ComplexSquareMatrix::ComplexSquareMatrix(const ComplexMatrix & inputMatrix)
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
// Class Constructor -- Sets a matrix equal to a diagonal matrix
// Input:       x:              Diagonal element
//              size:           # of rows or columns, default = 1
//                      
// Output:                                                      None
//
// Notes:
// ############################# Class Constructor ###############################
ComplexSquareMatrix::ComplexSquareMatrix(Complex x, int size)
{
  int           i, j, k;
  Complex       *temp;
//
// Initialize all arrays, base class constructor not called
//
  resetInstanceVariables();
//
//
//
  temp          = new Complex[size*size];
  k             = 0;
  for(i=0; i<size; i++)
  {
    for(j=0; j<size; j++)
    {
      if(i==j)
        temp[k] = x;
      else
        temp[k] = Complex(0., 0.);
      k++;
    }
  }
  assign(size, size, temp);

  delete [] temp;
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Sets a matrix equal to a diagonal matrix
// Input:       x:              array of diagonal elements
//              size:           # of rows or columns, default = 1
//                      
// Output:                                                      None
//
// Notes:
// ############################# Class Constructor ###############################
ComplexSquareMatrix::ComplexSquareMatrix(const Complex *elements, int size)
{
  int           i, j, k;
  Complex       *temp;
//
// Initialize all arrays, base class constructor not called
//
  resetInstanceVariables();
//
  temp          = new Complex[size*size];
  k             = 0;
  for(i=0; i<size; i++)
  {
    for(j=0; j<size; j++)
    {
      if(i==j)
        temp[k] = elements[j];
      else
        temp[k] = Complex(0., 0.);
      k++;
    }
  }
  assign(size, size, temp);

  delete [] temp;
  return;
}

// ############################# Class Destructor #################################
// ComplexSquareMatrix -- Destructor for the ComplexSquareMatrix class
// Input:               none
// Output:              none
//
// ############################# Class Destructor #################################
ComplexSquareMatrix::~ComplexSquareMatrix()
{
  delete [] _luMatrix;
  delete [] _pivot;
  return;
}


// ############################# Public Method ###############################
// assign -- Sets the ComplexSquareMatrix's element to the input element array
//
// Input:       inputMatrix:            array of elements for the ComplexSquareMatrix
//              numberRows:             new number of rows (if 0 use old number)
//              numberCols:             new number of columns (if 0 use old number)
//                      
// Output:                              None
//
// NOTES:
// 1. The element array needs to be ordered by rows: [0,0], [0,1], [0,2] ...
// ############################# Public Method ###############################
void ComplexSquareMatrix::assign(int numberRows, int numberCols, const Complex *inputMatrix)
{
  int   i;

  numberRows            = MIN(numberRows, numberCols);
  if(numberRows > 0)
  {
    _rows               = _columns      = numberRows;
    delete [] _matrix;
    _matrix             = new Complex [_rows*_columns];
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
// findDeterminant - Calculate the determinant.  Assumes LU factorization has been
//               performed. Taken from: C Tools for Engineers and Scientists.
//
// Input:                               None
// Output:                  Matrix determinant
//
// ############################ Public Function ###################################
Complex ComplexSquareMatrix::findDeterminant()
{
  int number_columns, i, sign;
  Complex d;

  if(_rows != _columns)
  {
    fprintf(stderr, "\ndeterminant ERROR: non-square, size = %d x %d\n", _rows, _columns);
    return Complex(0., 0.);
  }

  if(!luDecomp())
    return Complex(0., 0.);

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
// Input:                               None
//
// Output:              flag:           YES on success
//                                      Return pointer to inverted matrix
//
// NOTES:
// 1. This routine calls luDecomp().
// ############################ Public Function ###################################
ComplexSquareMatrix &ComplexSquareMatrix::invert(int *flag)
{
  int           i;
  Complex       *inverse_matrix;


  if(!luDecomp())
  {
    flag        = NO;
    return      *this;
  }

  inverse_matrix                = new Complex [_rows*_columns];
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
// Input:       input           The b vector
//                              
// Output:      flag:           = YES on success
//                              return pointer to result
//
// NOTES:
// 1. This routine calls luDecomp().
// ############################ Public Function ###################################
LOGICAL ComplexSquareMatrix::linearSolve(const ComplexMatrix &bVector, ComplexMatrix &xResult)
{
  int           i, result;
  Complex       *b;
  const         Complex *input_elements;

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
  if(!luDecomp())                                       // Find the LU decomposition
  {
    return NO;
  }
//
// Allocate temp vector
//
  b                     = new Complex[_columns];
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
// Operator transpose  Performs the transpose on the matrix
// Input:                                               none
//                      
// Output:              the transposed matrix
//
// ############################# Public Method ###############################
ComplexSquareMatrix& ComplexSquareMatrix::transpose()
{
  int           i, j, number_columns;
  Complex       temp;

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
// Operator hermitian  Performs the conjugate transpose on the matrix
// Input:               none
//                      
// Output:              the conjugate transposed matrix
//
// ############################# Public Method ###############################
ComplexSquareMatrix& ComplexSquareMatrix::hermitian()
{
  int           i, j, number_columns;
  Complex       temp;

  number_columns        = _columns;

//
// In place transpose
//
  for(i=0; i<_columns; i++)
  {
    for(j=i; j<_rows; j++)
    {
      temp                      = _matrix INDEX(i,j);
      _matrix INDEX(i,j)        = conj(_matrix INDEX(j, i));
      _matrix INDEX(j,i)        = conj(temp);
    }
  }
  return *this;
}


// ############################# Public Method ###############################
// Operator /=  Divide this by the input ComplexSquareMatrix.
// Input:       y:      input ComplexSquareMatrix
//                      
// Output:              the result ComplexSquareMatrix
//
// NOTES:
// 1. This translates into *this = *this/y = y_inverse*(*this)
// ############################# Public Method ###############################
ComplexSquareMatrix& ComplexSquareMatrix::operator /= (const ComplexSquareMatrix&  y)
{
  int                   flag;
  ComplexSquareMatrix   y_copy;
//
// If *this has a single column, than we can use linear solve instead
// of calculating the inverse:
//
  y_copy                = y;
//
// Do invert:
//
  ComplexSquareMatrix inverse_y = y_copy.invert(&flag);
  if(flag == NO)
  {
      fprintf(stderr, "Couldn't find inverse in /=\n");
      return *this;
  }
  *this         = inverse_y * (*this);

  return *this;
}


// ############################# Public Method ###############################
// Operator <<=  Shifts the matrix element locations by an integer
// Input:       shift:          Shift value
//                      
// Output:                      the result ComplexMatrix
//
// NOTES:
//   1. Zeros are shifted in at the left.
//
// ############################# Public Method ###############################
ComplexSquareMatrix& ComplexSquareMatrix::operator <<= (int shift)
{
  int   i, j, k, number_zeros, number_elements, number_shifts;

  number_elements       = _rows*_columns;
//
// Handle boundary cases
//
  if((shift >= number_elements) || (shift <= -number_elements))
  {
    for(i=0; i<number_elements; i++)
      _matrix[i]                = Complex(0., 0.);
    return *this;
  }
  if(shift == 0)
    return *this;
//
// Handle positive shifts
//
  if(shift > 0)
  {
    number_shifts       = number_elements - shift;
    j                   = shift;
    for(i=0; i<number_shifts; i++)
    {
      _matrix[i]        = _matrix[j++];
    }
    number_zeros        = shift;
    j                   = number_shifts;
    for(i=0; i<number_zeros; i++)
    {
      _matrix[j++]      = Complex(0., 0.);
    }
  }
  else
  {
    number_shifts       = number_elements + shift;
    j                   = number_elements-1;
    k                   = number_elements-1+shift;
    for(i=0; i<number_shifts; i++)
    {
      _matrix[j--]      = _matrix[k--];
    }
    number_zeros        = -shift;
    for(i=0; i<number_zeros; i++)
    {
      _matrix[i]                = Complex(0., 0.);
    }
  }
  return *this;
}

// ############################# Public Method ###############################
// Operator -  Returns the negative of the input ComplexSquareMatrix
// Input:       x:                      input ComplexSquareMatrix
//                      
// Output:              the output ComplexSquareMatrix
//
// NOTES
// ############################# Public Method ###############################
ComplexSquareMatrix operator - (const ComplexSquareMatrix&  x)
{
  return (Complex(-1.,0.)*x);
}

// ############################# Public Method ###############################
// Operator +  Returns the sum of the two input ComplexMatrices
// Input:       x:                      input ComplexSquareMatrix
//              y:                      second ComplexSquareMatrix
//                      
// Output:              the sum ComplexSquareMatrix
//
// ############################# Public Method ###############################
ComplexSquareMatrix operator + (const ComplexSquareMatrix&  x, const ComplexSquareMatrix&  y)
{
  ComplexSquareMatrix   sum = x;
  sum   += y;                                   //Stroustrup p. 3022
  return sum;
}

// ############################# Public Method ###############################
// Operator -  Returns the difference of the two input ComplexMatrices
// Input:       x:                      input ComplexSquareMatrix
//              y:                      second ComplexSquareMatrix
//                      
// Output:              the difference ComplexSquareMatrix
//
// NOTES:
// 1. For modulo-2 addition, + equals -.
// ############################# Public Method ###############################
ComplexSquareMatrix operator - (const ComplexSquareMatrix&  x, const ComplexSquareMatrix&  y)
{
  ComplexSquareMatrix   difference = x;
  difference    -= y;                                   //Stroustrup p. 302
  return difference;
}

// ############################# Public Method ###############################
// Operator *  Returns the product of the two input ComplexMatrices
// Input:       x:                      input ComplexSquareMatrix
//              y:                      second ComplexSquareMatrix
//                      
// Output:              the product ComplexSquareMatrix
//
// ############################# Public Method ###############################
ComplexSquareMatrix operator * (const ComplexSquareMatrix&  x, const ComplexSquareMatrix&  y)
{
  ComplexSquareMatrix   prod = x;
  prod  *= y;                                   //Stroustrup p. 302
  return prod;
}

// ############################# Public Method ###############################
// Operator /  Returns the division of the two input ComplexMatrices
// Input:       x:                      input ComplexSquareMatrix
//              y:                      second ComplexSquareMatrix
//                      
// Output:              the result ComplexSquareMatrix
//
// ############################# Public Method ###############################
ComplexSquareMatrix operator / (const ComplexSquareMatrix&  x, const ComplexSquareMatrix&  y)
{
  ComplexSquareMatrix   result = x;
  result        /= y;                                   //Stroustrup p. 3022
  return        result;
}

