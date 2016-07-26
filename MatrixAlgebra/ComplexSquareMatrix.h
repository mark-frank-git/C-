#ifndef _COMPLEXSQUAREMATRIX_H
#define _COMPLEXSQUAREMATRIX_H
/************************************************************************
 *                                                                      *
 * This subclass of ComplexMatrix implements a square Complex matrix.   *
 * The matrix is represented as a Complex array, with the elements      *
 * increasing by rows.                                                  *
 * The usual overloaded operators are included: +, -, *, /, etc.        *
 *                                                                      *
 * File:ComplexSquareMatrix.h                                           *
 *                                                                      *
 ************************************************************************/
#include        "ComplexMatrix.h"

class ComplexSquareMatrix: public ComplexMatrix
{
private:
  int           *_pivot;                                // pivot for LU decomp
  int           _oldSize;                               // old size for LU decomp alloc
  Complex       *_luMatrix;                             // LU decomposition
  LOGICAL       _luDecompPerformed;                     // YES = LU decomp was performed

//
// Private methods:
//
  LOGICAL       luDecomp();                             // Calculate LU decomp, rtn YES on success
  void          resetInstanceVariables();


public:
/*****************************
 * Constructor/Destructor    *
 *****************************/
  ComplexSquareMatrix(int size=1, const Complex *elements=NULL);
  ComplexSquareMatrix(const ComplexSquareMatrix &matrix);               // Copy constructor
  ComplexSquareMatrix(const ComplexMatrix &matrix);                     // Copy constructor
  ComplexSquareMatrix(Complex x, int size=1);                           // diagonal matrix, all elements the same
  ComplexSquareMatrix(const Complex *elements, int size=1);             // diagonal matrix, with array input
  ~ComplexSquareMatrix();

/**********************
 * Set parameters:    *
 **********************/
  void assign(int numberRows=0, int numberCols=0, const Complex *inputMatrix=NULL);
                                                                        // Sets a new ComplexMatrix

/**********************
 * Get parameters:    *
 **********************/
  LOGICAL       getLUDecompState()              const;                  // Get status of LU decomposition
  const Complex *getLuMatrix()                  const;                  // Return a pointer to the LU matrix

/**********************
 * Calculations:      *
 **********************/
  Complex       findDeterminant();                                      // Find and return matrix determinant
  ComplexSquareMatrix   & invert(int *flag);                            // Invert the matrix, rtn the inverse
  LOGICAL       linearSolve(const ComplexMatrix &b, ComplexMatrix &x);  // Solve Ax = b for x, rtn YES on success


/****************************************
 * The following operators              *
 * modify the ComplexSquareMatrix       *
 ****************************************/
  ComplexSquareMatrix&  transpose();                                    // Transpose operator
  ComplexSquareMatrix&  hermitian();                                    // Conjugate transpose

/****************************************
 * The following operators              *
 * modify the ComplexSquareMatrix       *
 * by a second ComplexSquareMatrix      *
 ****************************************/
  ComplexSquareMatrix&  operator /= (const ComplexSquareMatrix& y);

/****************************
 * The following operators      *
 * modify the ComplexMatrix     *
 * by a integer.                *
 ****************************/
  ComplexSquareMatrix&  operator <<=  (int shift);

};

/****************************
 * The following operators      *
 * take 1 or 2 ComplexMatrices  *
 * and return a third           *
 ****************************/
  ComplexSquareMatrix   operator - (const ComplexSquareMatrix& x);
  ComplexSquareMatrix   operator + (const ComplexSquareMatrix& x, const ComplexSquareMatrix& y);
  ComplexSquareMatrix   operator - (const ComplexSquareMatrix& x, const ComplexSquareMatrix& y);
  ComplexSquareMatrix   operator * (const ComplexSquareMatrix& x, const ComplexSquareMatrix& y);
  ComplexSquareMatrix   operator / (const ComplexSquareMatrix& x, const ComplexSquareMatrix& y);


/****************************************
 * The following operators              *
 * modify a ComplexSquareMatrix by a    *
 * Complex, and return a                *
 * ComplexSquareMatrix.                 *
 ****************************************/
  inline ComplexSquareMatrix operator + (Complex x, const ComplexSquareMatrix&  y);
  inline ComplexSquareMatrix operator - (Complex x, const ComplexSquareMatrix&  y);
  inline ComplexSquareMatrix operator * (Complex x, const ComplexSquareMatrix&  y);

  inline ComplexSquareMatrix operator + (const ComplexSquareMatrix&  x, Complex y);
  inline ComplexSquareMatrix operator - (const ComplexSquareMatrix&  x, Complex y);
  inline ComplexSquareMatrix operator * (const ComplexSquareMatrix&  x, Complex y);
  inline ComplexSquareMatrix operator / (const ComplexSquareMatrix&  x, Complex y);
  inline ComplexSquareMatrix operator / (Complex x, const ComplexSquareMatrix& y);


/****************************************
 * Inline operator functions defined    *
 * above.                               *
 ****************************************/
inline ComplexSquareMatrix operator + (Complex x, const ComplexSquareMatrix &y)
{
  ComplexSquareMatrix   sum = y;
  sum   += x;                                   //Stroustrup p. 302
  return sum;
}
inline ComplexSquareMatrix operator - (Complex x, const ComplexSquareMatrix &y)
{
  ComplexSquareMatrix   result = -y;
  result        += x;                           //Stroustrup p. 302
  return result;
}
inline ComplexSquareMatrix operator * (Complex x, const ComplexSquareMatrix &y)
{
  ComplexSquareMatrix   product = y;
  product       *= x;                           //Stroustrup p. 302
  return product;
}
inline ComplexSquareMatrix operator / (Complex x, const ComplexSquareMatrix &y)
{
  ComplexSquareMatrix   dividend = x;
  dividend      /= y;                           //Stroustrup p. 302
  return dividend;
}
inline ComplexSquareMatrix operator + (const ComplexSquareMatrix &x, Complex y)
{
  ComplexSquareMatrix   sum = x;
  sum   += y;                                   //Stroustrup p. 302
  return sum;
}
inline ComplexSquareMatrix operator - (const ComplexSquareMatrix &x, Complex y)
{
  ComplexSquareMatrix   difference = x;
  difference    -= y;                           //Stroustrup p. 302
  return difference;
}
inline ComplexSquareMatrix operator * (const ComplexSquareMatrix &x, Complex y)
{
  ComplexSquareMatrix   product = x;
  product       *= y;                           //Stroustrup p. 302
  return product;
}
inline ComplexSquareMatrix operator / (const ComplexSquareMatrix &x, Complex y)
{
  ComplexSquareMatrix   dividend = x;
  dividend      /= y;                           //Stroustrup p. 302
  return dividend;
}
inline ComplexSquareMatrix operator << (const ComplexSquareMatrix &x, int shift)
{
  ComplexSquareMatrix   dividend = x;
  dividend      <<= shift;                      //Stroustrup p. 302
  return dividend;
}

// Return the LU decomposition state
  inline LOGICAL        ComplexSquareMatrix::getLUDecompState()  const  {return _luDecompPerformed;}
// Return a pointer to the LU matrix
  inline const Complex  *ComplexSquareMatrix::getLuMatrix()     const   {return _luMatrix;} 


#endif
