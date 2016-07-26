#ifndef _DOUBLESQUAREMATRIX_H
#define _DOUBLESQUAREMATRIX_H
/************************************************************************
 *                                                                      *
 * This subclass of DoubleMatrix implements a double precision square   *
 * matrix.                                                              *
 * The matrix is represented as a double array, with the elements       *
 * increasing by rows.                                                  *
 * The usual overloaded operators are included: +, -, *, /, etc.        *
 * This class is a subclass of the DoubleMatrix class.                  *
 *                                                                      *
 * File:DoubleSquareMatrix.h                                            *
 *                                                                      *
************************************************************************/

#include        "DoubleMatrix.h"

class DoubleSquareMatrix: public DoubleMatrix
{
private:
  int           *_pivot;                                        // pivot for LU decomp
  int           _oldSize;                                       // old size for LU allocate
  double        *_luMatrix;                                     // LU decomposition
  LOGICAL       _luDecompPerformed;                             // YES = LU decomp was performed

//
// Private methods:
//
  LOGICAL       luDecomp();                                     // Calculate LU decomp, rtn YES on success
  void          resetInstanceVariables();

public:
//
// Public methods:
//
/*********************************
 * Constructors, destructors:    *
 *********************************/
  DoubleSquareMatrix(int size=1, const double *elements=NULL);  // Initializes a new matrix
  DoubleSquareMatrix(const DoubleSquareMatrix &matrix);         // Copy constructor
  DoubleSquareMatrix(const DoubleMatrix &matrix);               // Copy constructor from general matrix
  DoubleSquareMatrix(double x, int size=1);                     // diagonal matrix
  ~DoubleSquareMatrix();

/**********************
 * Set parameters:    *
 **********************/
  void assign(int rows=0, int cols=0, const double *inputMatrix=NULL);  // Sets a new DoubleSquareMatrix

/**********************
 * Get parameters:    *
 **********************/
  LOGICAL       getLUDecompState() const;                       // Get status of LU decomposition
  const double *getLuMatrix() const;                            // Return a pointer to the LU matrix
                
/**********************
 * Calculations:      *
 **********************/
  double                findDeterminant();                      // Find and return matrix determinant
  DoubleSquareMatrix    &invert(LOGICAL *flag);                 // Invert the matrix, rtn the inverse
  LOGICAL               linearSolve(const DoubleMatrix &b, DoubleMatrix &x);    // Solve Ax = b for x, rtn YES on success

/****************************************
 * The following operators              *
 * modify the DoubleSquareMatrix        *
 ****************************************/
  DoubleSquareMatrix&   transpose();                            // Transpose operator

/****************************
 * The following operators      *
 * modify the DoubleMatrix      *
 * by a second DoubleMatrix     *
 ****************************/
  DoubleSquareMatrix&   operator /= (const DoubleSquareMatrix& y);

};
//
// The following are nonmember functions taking two arguments
//

/****************************
 * The following operators      *
 * take 1 or 2 DoubleMatrixs    *
 * and return a third           *
 ****************************/
  DoubleSquareMatrix    operator - (const  DoubleSquareMatrix& x);
  DoubleSquareMatrix    operator + (const  DoubleSquareMatrix& x, const  DoubleSquareMatrix& y);
  DoubleSquareMatrix    operator - (const  DoubleSquareMatrix& x, const  DoubleSquareMatrix& y);
  DoubleSquareMatrix    operator * (const  DoubleSquareMatrix& x, const  DoubleSquareMatrix& y);
  DoubleSquareMatrix    operator / (const  DoubleSquareMatrix& x, const  DoubleSquareMatrix& y);


/********************************
 * The following operators      *
 * modify a double by a         *
 * DoubleMatrix, and return a   *
 * DoubleMatrix.                *
 ********************************/

/********************************
 * The following operators      *
 * modify a double by a         *
 * DoubleMatrix, and return a   *
 * DoubleMatrix.                *
 ********************************/
  inline DoubleSquareMatrix     operator + (double x, const  DoubleSquareMatrix &y);
  inline DoubleSquareMatrix     operator - (double x, const  DoubleSquareMatrix &y);
  inline DoubleSquareMatrix     operator * (double x, const  DoubleSquareMatrix &y);

  inline DoubleSquareMatrix     operator + (const  DoubleSquareMatrix &x, double y);
  inline DoubleSquareMatrix     operator - (const  DoubleSquareMatrix &x, double y);
  inline DoubleSquareMatrix     operator * (const  DoubleSquareMatrix &x, double y);
  inline DoubleSquareMatrix     operator / (const  DoubleSquareMatrix &x, double y);
  inline DoubleSquareMatrix     operator / (double x, const DoubleSquareMatrix &y);

/************************
 * Inline functions     *
 * that are defined     *
 * above.               *
 ************************/
inline DoubleSquareMatrix  operator + (double x, const DoubleSquareMatrix &y)
{
  DoubleSquareMatrix    sum = y;
  sum   += x;                           //Stroustrup p. 302
  return sum;
}
inline DoubleSquareMatrix  operator - (double x, const DoubleSquareMatrix &y)
{
  DoubleSquareMatrix    result = -y;
  result        += x;                   //Stroustrup p. 302
  return result;
}
inline DoubleSquareMatrix  operator * (double x, const DoubleSquareMatrix &y)
{
  DoubleSquareMatrix    product = y;
  product       *= x;                   //Stroustrup p. 302
  return product;
}
inline DoubleSquareMatrix  operator + (const DoubleSquareMatrix &x, double y)
{
  DoubleSquareMatrix    sum = x;
  sum   += y;                           //Stroustrup p. 302
  return sum;
}
inline DoubleSquareMatrix  operator - (const DoubleSquareMatrix &x, double y)
{
  DoubleSquareMatrix    difference = x;
  difference    -= y;                   //Stroustrup p. 302
  return difference;
}
inline DoubleSquareMatrix  operator * (const DoubleSquareMatrix &x, double y)
{
  DoubleSquareMatrix    product = x;
  product       *= y;                   //Stroustrup p. 302
  return product;
}
inline DoubleSquareMatrix  operator / (const DoubleSquareMatrix &x, double y)
{
  DoubleSquareMatrix    dividend = x;
  dividend      /= y;                   //Stroustrup p. 302
  return dividend;
}
inline DoubleSquareMatrix       operator / (double x, const DoubleSquareMatrix &y)
{
  DoubleSquareMatrix    dividend = x;
  dividend              /= y;           //Stroustrup p. 302
  return dividend;
}

// Return the LU decomposition state
  inline LOGICAL        DoubleSquareMatrix::getLUDecompState()  const   {return _luDecompPerformed;}
// Return a pointer to the LU matrix
  inline const double   *DoubleSquareMatrix::getLuMatrix()      const   {return _luMatrix;} 

#endif
