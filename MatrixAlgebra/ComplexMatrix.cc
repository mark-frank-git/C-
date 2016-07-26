/************************************************************************
 *                                                                      *
 * This subclass of object implements a complex matrix.                 *
 * The matrix is represented as a Complex array, with the elements      *
 * increasing by rows.                                                  *
 * The usual overloaded operators are included: +, -, *, /, etc.        *
 *                                                                      *
 * File:ComplexMatrix.cc                                                *
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
#include "ComplexMatrix.h"
#include "matrix_routines.h"


#define MAX(a,b)  ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))
#ifndef YES
#define YES     1
#define NO      0
#endif

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the ComplexMatrix class.
// Input:       numberRows:             Rows in the matrix
//              numberCols:             Columns in the matrix
//              input_elements:         Input matrix elements
//                      
// Output:                                                      None
//
// Notes:
// ############################# Class Constructor ###############################
ComplexMatrix::ComplexMatrix(int numberRows, int numberCols, const Complex *input_elements)
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
// Class Constructor -- Sets all elements equal to input complex element
// Input:       numberRows:     # of rows
//              numberCols:     # of cols
//              x:              Complex element
//                      
// Output:                                                      None
//
// Notes:
// ############################# Class Constructor ###############################
ComplexMatrix::ComplexMatrix(int numberRows, int numberCols, Complex x)
{
  int           i;
  Complex       *elements;
//
// Initialize arrays:
// NOTE: if we add any more arrays, we need to modify constructors in ComplexSquareMatrix
//
  _matrixString                 = NULL;
  _matrix                       = NULL;
  _singularValues               = NULL;
  _rows                         = MAX(1, numberRows);
  _columns                      = MAX(1, numberCols);
//
// Set all elements equal to input element:
//
  elements                      = new Complex[_rows*_columns];
  for(i=0; i<_rows*_columns; i++)
  {
    elements[i]                 = x;
  }
//
// Use assign operator to initialize matrix:
//
  assign(_rows, _columns, elements);
  delete [] elements;
  return;
}


// ############################# Class Constructor ###############################
// Class Constructor -- Copy constructor for the ComplexMatrix class.
// Input:       matrix:                 A previously allocated ComplexMatrix
//                      
// Output:                              None
//
// Notes:
// ############################# Class Constructor ###############################
ComplexMatrix::ComplexMatrix(const ComplexMatrix&  inputMatrix)
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
  assign(inputMatrix.getRows(), inputMatrix.getColumns(), inputMatrix.getMatrix());
  return;
}

// ############################# Class Destructor #################################
// ComplexMatrix -- Destructor for the ComplexMatrix class
// Input:               none
// Output:              none
//
// ############################# Class Destructor #################################
ComplexMatrix::~ComplexMatrix()
{
  delete [] _matrix;
  delete [] _singularValues;
  delete [] _matrixString;
  return;
}


// ############################# Public Method ###############################
// assign -- Sets the ComplexMatrix's element to the input element array
//
// Input:       inputMatrix:            array of elements for the ComplexMatrix
//              numberRows:             new number of rows (if 0 use old number)
//              numberCols:             new number of columns (if 0 use old number)
//                      
// Output:                              None
//
// NOTES:
// 1. The element array needs to be ordered by rows: [0,0], [0,1], [0,2] ...
// ############################# Public Method ###############################
void ComplexMatrix::assign(int numberRows, int numberCols, const Complex *inputMatrix)
{
  int   i;

  if(numberRows > 0)
  {
    _rows               = numberRows;
    if(numberCols > 0)
      _columns  = numberCols;
    delete [] _matrix;
    _matrix             = new Complex [_rows*_columns];
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

// ############################# Public Method ###############################
// singularValues  Returns the singular values from largest to smallest
// Input:               none
//                      
// Output:              singular values
//
// NOTES:
// 1. Not yet working for Complex matrices.
// ############################# Public Method ###############################
const double *ComplexMatrix::singularValues()
{
  int           n;
  Complex       *u, *v;
 
  delete [] _singularValues;

  n                     = MAX(_rows, _columns);
  _singularValues       = new double[n];
  u                     = new Complex[n*n];
  v                     = new Complex[n*n];
  svd(_matrix, _rows, _columns, _singularValues, u, v);

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
double ComplexMatrix::cond()
{
  int           n;
  double        cond;
  const double *singular_values;

  singular_values       = singularValues();

  n                     = MIN(_rows, _columns);
  if(fabs(singular_values[n-1]) > 0.)
  {
    cond                = fabs(singular_values[0])/fabs(singular_values[n-1]);
  }
  else
    cond                = 0.;

  return cond;
}

// ############################ Public Function ###################################
// toeplitzSolve - Finds the solution to Ax=b, using the Toeplitz matrix algorithm
//                 as given in "Digital Signal Processing," by Roberts and Mullis.
//
// Input:       r:                      The first row of the Toeplitz matrix
//              b:                      The result of Ax=b
//              n:                      The dimension of the Toeplitz matrix
//
// Output:      conda:                  Estimate of the 2-norm of the matrix
//              flag:                   0 = normal, 1 = illegal n, 2 = singular matrix found
//
// NOTES:
// 1. This function taken from the Netlib function dsytcl()
// ############################ Public Function ###################################
ComplexMatrix & ComplexMatrix::toeplitzSolve(Complex *r, Complex *b, int n, double *conda, int *flag)
{
  int           i,  k, ip;
  Complex       *y, *x, *w;
  Complex       alpha, beta, c, eta, gamma, t;
  double        s, smin, yn2;
  LOGICAL       compute_error;
//
// Initialization:
//
  k                     = 0;
  *flag                 = 0;
  if(n <= 0)
  {
    *flag               = 1;
    return *this;
  }
//
// Allocate arrays:
//

  x                     = new Complex[n];
  y                     = new Complex[n];
  w                     = new Complex[n];
  compute_error         = NO;
  while (YES)                           // break out of this loop on error
  {
    if(abs(r[0]) == 0.)
    {
      compute_error     = YES;
      break;
    }
//
// Compute the first elements of x and y:
//
    x[0]                = b[0]/r[0];
    y[0]                = -r[1]/r[0];
    gamma               = r[0] + r[1]*y[0];
    smin                = abs(r[0]);
    s                   = smin;
    k                   = 0;            // # of elements to process
    while((n>1) && (k<(n-1)))
    {
//
// Estimate the smallest singular value of leading submatrix
//
      ip                = izamax(k+1, y, 1);
      yn2               = norm(y[ip]);
      s                 = abs(gamma)/MAX(1., yn2);
      if(s < smin)
        smin            = s;
//
// Perform a classical Levinson step and update the algorithm:
//
      beta              = b[k+1] - zdot(k+1, &r[1], 1, &x[0], -1);
      c                 = -r[(k+2)%n] - zdot(k+1, &r[1], 1, &y[0], -1);
      if(abs(gamma) == 0.)
        break;
      alpha             = beta/gamma;
      eta               = c/gamma;
      for(i=0; i<=k; i++)
      {
        x[i]            += alpha*y[k-i];
        w[i]            = y[i] + eta*y[k-i];
      }
      for(i=0; i<=k; i++)
        y[i]            = w[i];
      k++;
      x[k]              = alpha;
      y[k]              = eta;
      gamma             -= c*eta;
    }
    t                   = zasum(n, &r[0], 1);
    *conda              = abs(t)/smin;
    break;
  }
  if(compute_error)
    *flag               = 2;
//
// Re-size *this, and copy solution:
//
  _columns      = 1;
  _rows         = n;
  delete [] _matrix;
  _matrix       = new Complex[_rows];
  for(i=0; i<_rows; i++)
    _matrix[i]  = x[i];

  delete [] x;
  delete [] y;
  delete [] w;
  return        *this;
}


// ############################# Public Method ###############################
// Operator transpose  Performs the transpose on the matrix
// Input:                                               none
//                      
// Output:              the transposed matrix
//
// ############################# Public Method ###############################
ComplexMatrix& ComplexMatrix::transpose()
{
  int           i, j, old_rows, old_columns, number_columns;
  Complex       *temp;

  number_columns        = _columns;

//
// Allocate temp and copy matrix
//
  temp                  = new Complex[_rows*_columns];
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

// ############################# Public Method ###############################
// Operator hermitian  Performs the conjugate transpose on the matrix
// Input:               none
//                      
// Output:              the conjugate transposed matrix
//
// ############################# Public Method ###############################
ComplexMatrix& ComplexMatrix::hermitian()
{
  int           i, j, old_rows, old_columns, number_columns;
  Complex       *temp;

  number_columns        = _columns;

//
// Allocate temp and copy matrix
//
  temp                  = new Complex[_rows*_columns];
  for(i=0; i<_rows*_columns; i++)
    temp[i]             = conj(_matrix[i]);

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


// ############################# Public Method ###############################
// Operator reverseElements  Reverses the elements (indices) of a matrix
// Input:               none
//                      
// Output:              the reversed matrix
//
// Notes:
// 1. This may only make sense for a vector.
//
// ############################# Public Method ###############################
ComplexMatrix& ComplexMatrix::reverseElements()
{
  int           i, j;
  Complex       *temp;

//
// Copy to temp, and then reverse.  There is probably a more
// efficient way to implement this.
//
  temp                  = new Complex[_rows*_columns];
  for(i=0; i<_rows*_columns; i++)
    temp[i]             = _matrix[i];
  j                     = _rows*_columns-1;
  for(i=0; i<_rows*_columns;i++)
    _matrix[i]          = temp[j--];

  delete [] temp;
  return *this;
}


// ############################# Public Method ###############################
// Operator replaceRow  Replaces a row in the matrix
// Input:       rowNumber       : rowNumber to replace
//              newRow:         : new row of matrix
//                      
// Output:                      the row replaced matrix
//
// Notes:
// 1. The new row matrix should be a 1xN matrix, otherwise, we just take the
//    first row of the newRow matrix.
//
// ############################# Public Method ###############################
ComplexMatrix&  ComplexMatrix::replaceRow(int rowNumber, const ComplexMatrix& newRow)
{
  int           j, k, new_cols;
  const Complex *new_elements;

  new_cols      = newRow.getColumns();
  if(new_cols != _columns)
  {
    fprintf(stderr,"Error: new cols = %d, old cols = %d in replaceRow\n", new_cols, _columns);
    return *this;
  }
  if(rowNumber > (_rows-1))
  {
    fprintf(stderr,"Error: replace row = %d, number rows = %d in replaceRow\n", rowNumber, _rows);
    return *this;
  }
//
// Replace the row with the new row:
//
  new_elements          = newRow.getMatrix();
  k                     = rowNumber*_columns;
  for(j=0; j<_columns; j++)
  {
      _matrix[k++]      = new_elements[j];
  }

  return        *this;
}


// ############################# Public Method ###############################
// Operator postTimesDiagonal: post multiplies the matrix times a diagonal matrix.
// Input:       y               : post multiplying matrix
//                      
// Output:                      the multiplied matrix
//
// Notes:
// 1. It would probably be cleaner to develop a ComplexDiagonalMatrix class and
//    use it instead.
// ############################# Public Method ###############################
ComplexMatrix&  ComplexMatrix::postTimesDiagonal(ComplexMatrix &y)
{
  int           i, j, number_columns, result_index, y_columns;
  Complex       *result;
  Complex       dot_product;
  const         Complex *y_elements;
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
  y_columns             = y.getColumns();
  result                = new Complex[_rows*y_columns];
  number_columns        = y.getColumns();
  y_elements            = y.getMatrix();
  for(i=0; i<_rows; i++)
  {
    for(j=0; j<y_columns; j++)
    {
      result_index              = i*y_columns + j;
      result[result_index]      = _matrix[i*_columns + j]*y_elements[j*y_columns + j];
    }
  }
//
// Re-size matrix, copy result:
//
  if(_columns != y.getColumns())
  {
    _columns    = y.getColumns();
    delete [] _matrix;
    _matrix     = new Complex[_rows*_columns];
  }
  for(i=0; i<_rows*_columns; i++)
    _matrix[i]  = result[i];

  delete        [] result;
  return        *this;
}


#define MAX_COEFF_SIZE  64
// ############################# Public Method ###############################
// getMatrixString()  - Get a string for printing out the matrix
// Input:                       None
//          
// Output:                      The output string
//
// ############################# Public Method ###############################
const char *ComplexMatrix::getMatrixString()
{
  int           i,j,cols_lin, n;
  int           number_rows, number_columns;
  const         Complex *elements;
  char          temp[MAX_COEFF_SIZE];

  if(_matrixString == NULL)
    _matrixString       = new char[MAX_CMPLX_STRING_SIZE];
  _matrixString[0]      = '\0';
  n                     = MAX_CMPLX_STRING_SIZE - 1;
//
// Get the matrix parameters:
//
  number_columns        = _columns;
  number_rows           = _rows;
  elements              = _matrix;
  cols_lin = 5;                 // 5 per line for Complex types

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
       sprintf(temp, "%5.2e + %5.2e i  ", real(elements INDEX(i,j)), imag(elements INDEX(i,j)));
       strncat(_matrixString, temp, n);
       n                -= strlen(temp);
    }
  }
  sprintf(temp,"\n");
  strncat(_matrixString, temp, n);
  return _matrixString;
}


// ############################# Public Method ###############################
// Operator =  Sets the ComplexMatrix equal to the input ComplexMatrix
// Input:       y:      input ComplexMatrix
//                      
// Output:              the result ComplexMatrix
//
// ############################# Public Method ###############################
ComplexMatrix& ComplexMatrix::operator = (const ComplexMatrix&  y)
{  
  if( this == &y )                                                      // Check for x=x
    return *this;
//
// Copy using assign function.
//
  assign(y.getRows(), y.getColumns(), y.getMatrix());

  return *this;
}

// ############################# Public Method ###############################
// Operator +=  Add the input ComplexMatrix to this.
// Input:       y:      input ComplexMatrix
//                      
// Output:              the result ComplexMatrix
//
// NOTES:
// ############################# Public Method ###############################
ComplexMatrix& ComplexMatrix::operator += (const ComplexMatrix&  y)
{
  int   i, min_size;
  const Complex *input_elements;
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
// Operator -=  Sub the input ComplexMatrix to this.
// Input:       y:      input ComplexMatrix
//                      
// Output:              the result ComplexMatrix
//
// NOTES:
// ############################# Public Method ###############################
ComplexMatrix& ComplexMatrix::operator -= (const ComplexMatrix&  y)
{
  int   i, min_size;
  const Complex *input_elements;
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
// Operator *=  Postmultiply this by the input ComplexMatrix.
// Input:       y:      input ComplexMatrix
//                      
// Output:              the result ComplexMatrix
//
// NOTES:
// 1. This translates into *this = *this*y
// ############################# Public Method ###############################
ComplexMatrix& ComplexMatrix::operator *= (const ComplexMatrix&  y)
{
  int           i, j, k, number_columns, result_index, y_columns;
  int           y_index;
  Complex       *result, *a_row_ptr;
  Complex       dot_product;
  const         Complex *y_elements, *y_ptr;
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
  y_columns             = y.getColumns();
  result                = new Complex[_rows*y_columns];
  number_columns        = y.getColumns();
  y_elements            = y.getMatrix();
  for(i=0; i<_rows; i++)
  {
    for(j=0; j<y_columns; j++)
    {
      a_row_ptr                 = &_matrix[i*_columns];
      dot_product               = Complex(0., 0.);
      y_index                   = j;
      y_ptr                     = &y_elements[y_index];
      for(k=0; k<_columns; k++)
      {
        dot_product             += (*a_row_ptr++) * (*y_ptr);
        y_ptr                   += y_columns;
      }
      result_index              = i*y_columns + j;
      result[result_index]      = dot_product;
    }
  }
//
// Re-size matrix, copy result:
//
  if(_columns != y.getColumns())
  {
    _columns    = y.getColumns();
    delete [] _matrix;
    _matrix     = new Complex[_rows*_columns];
  }
  for(i=0; i<_rows*_columns; i++)
    _matrix[i]  = result[i];

  delete        [] result;
  return        *this;
}

// ############################# Public Method ###############################
// Operator =  Sets the ComplexMatrix equal to the input Complex
// Input:       y:      input Complex
//                      
// Output:              the result ComplexMatrix
//
// ############################# Public Method ###############################
ComplexMatrix& ComplexMatrix::operator = (Complex y)
{
  _rows         = _columns      = 1;
  _matrix[0]    = y;
  return *this;
}

// ############################# Public Method ###############################
// Operator +=  Adds a Complex to the ComplexMatrix.
// Input:               y:              Complex to add
//                      
// Output:                              the result ComplexMatrixx
//
// NOTES:
//
// ############################# Public Method ###############################
ComplexMatrix& ComplexMatrix::operator += (Complex y)
{
  for(int i=0; i<_rows*_columns; i++)
    _matrix[i]  += y;
  return *this;
}

// ############################# Public Method ###############################
// Operator -=  Subtracts a Complex from the ComplexMatrix.
// Input:               y:              Complex to subtract
//                      
// Output:                              the result ComplexMatrixx
//
// NOTES:
//   1. 
//
// ############################# Public Method ###############################
ComplexMatrix& ComplexMatrix::operator -= (Complex y)
{
  for(int i=0; i<_rows*_columns; i++)
    _matrix[i]  -= y;
  return *this;
}

// ############################# Public Method ###############################
// Operator *=  Multiplies a Complex by the ComplexMatrix.
// Input:               y:              Complex to mult
//                      
// Output:                              the result ComplexMatrixx
//
// NOTES:
//   1. 
//
// ############################# Public Method ###############################
ComplexMatrix& ComplexMatrix::operator *= (Complex y)
{
  for(int i=0; i<_rows*_columns; i++)
    _matrix[i]  *= y;
  return *this;
}

// ############################# Public Method ###############################
// Operator /=  Divides the ComplexMatrix by a Complex
// Input:                       y:      Complex to divide into
//                      
// Output:                              the result ComplexMatrix
//
// NOTES:
//   1. 
//
// ############################# Public Method ###############################
ComplexMatrix& ComplexMatrix::operator /= (Complex y)
{
  if(y != 0.)
  {
    for(int i=0; i<_rows*_columns; i++)
      _matrix[i]        /= y;
  }
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
ComplexMatrix& ComplexMatrix::operator <<= (int shift)
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


/********************** Won't compile on NT ***************************************
// ############################ Public Function ###################################
// << - print out the matrix
//
// Input:       s:                      Open stream
//              x:                      Matrix to print
//
// Output:                  None
//
// NOTES:
// ############################ Public Function ###################################
ostream&  operator << (ostream& s, const ComplexMatrix&  x)
{
  int           i,j,cols_lin;
  int           number_rows, number_columns;
  const         Complex *elements;
//
// Get the matrix parameters:
//
  number_columns        = x.getColumns();
  number_rows           = x.getRows();
  elements              = x.getMatrix();
  cols_lin = 5;                 // per line for Complex types

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
*********************************************************************************/

// ############################# Public Method ###############################
// Operator -  Returns the negative of the input ComplexMatrix
// Input:       x:                      input ComplexMatrix
//                      
// Output:              the output ComplexMatrix
//
// NOTES
// ############################# Public Method ###############################
ComplexMatrix operator - (const ComplexMatrix&  x)
{
  return (Complex(-1.,0.)*x);
}

// ############################# Public Method ###############################
// Operator +  Returns the sum of the two input ComplexMatrices
// Input:       x:                      input ComplexMatrix
//              y:                      second ComplexMatrix
//                      
// Output:              the sum ComplexMatrix
//
// ############################# Public Method ###############################
ComplexMatrix operator + (const ComplexMatrix&  x, const ComplexMatrix&  y)
{
  ComplexMatrix sum = x;
  sum   += y;                                   //Stroustrup p. 3022
  return sum;
}

// ############################# Public Method ###############################
// Operator -  Returns the difference of the two input ComplexMatrices
// Input:       x:                      input ComplexMatrix
//              y:                      second ComplexMatrix
//                      
// Output:              the difference ComplexMatrix
//
// NOTES:
// 1. For modulo-2 addition, + equals -.
// ############################# Public Method ###############################
ComplexMatrix operator - (const ComplexMatrix&  x, const ComplexMatrix&  y)
{
  ComplexMatrix difference = x;
  difference    -= y;                                   //Stroustrup p. 302
  return difference;
}

// ############################# Public Method ###############################
// Operator *  Returns the product of the two input ComplexMatrices
// Input:       x:                      input ComplexMatrix
//              y:                      second ComplexMatrix
//                      
// Output:              the product ComplexMatrix
//
// ############################# Public Method ###############################
ComplexMatrix operator * (const ComplexMatrix&  x, const ComplexMatrix&  y)
{
  ComplexMatrix prod = x;
  prod  *= y;                                   //Stroustrup p. 302
  return prod;
}


// ############################# Public Method ###############################
// appendColumns  Appends columns to the first input Matrix
//
// Input:       x:                      input ComplexMatrix
//              y:                      second ComplexMatrix
//                      
// Output:              the result ComplexMatrix
//
// ############################# Public Method ###############################
ComplexMatrix appendColumns(const ComplexMatrix&  x, const ComplexMatrix&  y)
{
  int           i, j, k, x_index, y_index;
  int           x_m, x_n, y_m, y_n, result_m, result_n;
  ComplexMatrix result = x;
  Complex       *temp;
  const Complex *x_elements, *y_elements;

  x_m           = x.getRows();
  x_n           = x.getColumns();
  y_m           = y.getRows();
  y_n           = y.getColumns();
  result_m      = x_m;
  result_n      = x_n + y_n;
  if(x_m != y_m)
  {
    fprintf(stderr,"Error: x rows = %d, y rows = %d in appendColumns\n", x_m, y_m);
    return result;
  }
//
// Get elements of x and y, and copy into temp:
//
  temp          = new Complex[result_m*result_n];
  x_elements    = x.getMatrix();
  y_elements    = y.getMatrix();
  k             = 0;
  x_index       = y_index       = 0;
  for(i=0; i<result_m; i++)
  {
    for(j=0; j<x_n; j++)
    {
      temp[k++] = x_elements[x_index++];
    }
    for(j=0; j<y_n; j++)
    {
      temp[k++] = y_elements[y_index++];
    }
  }
//
// Assign the elements to the result matrix
//
  result.assign(result_m, result_n, temp);

  delete [] temp;
  return        result;
}


// ############################# Public Method ###############################
// appendRows  Appends rows to the first input Matrix
//
// Input:       x:                      input ComplexMatrix
//              y:                      second ComplexMatrix
//                      
// Output:              the result ComplexMatrix
//
// ############################# Public Method ###############################
ComplexMatrix appendRows(const ComplexMatrix&  x, const ComplexMatrix&  y)
{
  int           j, k;
  int           x_m, x_n, y_m, y_n, result_m, result_n;
  ComplexMatrix result = x;
  Complex       *temp;
  const Complex *x_elements, *y_elements;

  x_m           = x.getRows();
  x_n           = x.getColumns();
  y_m           = y.getRows();
  y_n           = y.getColumns();
  result_m      = x_m + y_m;
  result_n      = x_n;
  if(x_n != y_n)
  {
    fprintf(stderr,"Error: x cols = %d, y cols = %d in appendRows\n", x_n, y_n);
    return result;
  }
//
// Get elements of x and y, and copy into temp:
//
  temp          = new Complex[result_m*result_n];
  x_elements    = x.getMatrix();
  y_elements    = y.getMatrix();
  k             = 0;
   for(j=0; j<x_n*x_m; j++)
  {
      temp[k++] = x_elements[j];
  }
  for(j=0; j<x_n*y_n; j++)
  {
      temp[k++] = y_elements[j];
  }
//
// Assign the elements to the result matrix
//
  result.assign(result_m, result_n, temp);

  delete        [] temp;
  return        result;
}

