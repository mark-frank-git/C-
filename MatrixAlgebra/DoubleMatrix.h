#ifndef _DOUBLEMATRIX_H
#define _DOUBLEMATRIX_H
/************************************************************************
 *                                                                      *
 * This subclass of object implements a double precision matrix.        *
 * The matrix is represented as a double array, with the elements       *
 * increasing by rows.                                                  *
 * The usual overloaded operators are included: +, -, *, /, etc.        *
 *                                                                      *
 * File:DoubleMatrix.h                                                  *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/05/93  - Started                                              *
 *  2. 08/31/95  - Objective C -> C++.                                  *
 *  3. 02/27/97  - Added operator overloading.                          *
 *  4. 03/10/97  - print() -> <<.                                       *
 *  5. 05/19/05  - add singular values.                                 *
 ************************************************************************/
/////////////// Won't compile on NT:       #include <fstream.h>

#ifndef LOGICAL
#define LOGICAL                 char
#endif
#ifndef YES
#define YES                     1
#define NO                      0
#endif
#define MAX_MATRIX_STRING_SIZE  8192            // Matrix string size

class DoubleMatrix
{
protected:
  int           _rows, _columns;                // # of rows, cols in matrix
  double        *_matrix;                       // matrix elements
  double        *_singularValues;               // singular values
  char          *_matrixString;                 // For printing out the matrix

//
// Private methods:
//

public:
//
// Public methods:
//
/*********************************
 * Constructors, destructors:    *
 *********************************/
  DoubleMatrix(int numberRows=1, int numberCols=1, const double *elements=NULL); // Initializes a new matrix
  DoubleMatrix(const DoubleMatrix &matrix);                             // Copy constructor
  virtual ~DoubleMatrix();

/**********************
 * Set parameters:    *
 **********************/
  virtual void assign(int numberRows=0, int numberCols=0, const double *inputMatrix=NULL);      // Sets a new DoubleMatrix

/**********************
 * Get parameters:    *
 **********************/
  int   getRows() const;                                                // Get the number of rows in the matrix
  int   getColumns() const;                                             // Get the number of columns
  const double *getMatrix() const;                                      // Return a pointer to the matrix
                
/**********************
 * Calculations:      *
 **********************/
  const double *singularValues();
  double        cond();

/****************************
 * The following operators      *
 * modify the DoubleMatrix      *
 ****************************/
  virtual DoubleMatrix& transpose();                                            // Transpose operator

/********************************
 * The following functions      *
 * aid in printing the          *
 * matrix                       *
 ********************************/
  const char *getMatrixString();

/****************************
 * The following operators      *
 * modify the DoubleMatrix      *
 * by a second DoubleMatrix     *
 ****************************/
  DoubleMatrix& operator =  (const DoubleMatrix& y);
  DoubleMatrix& operator += (const DoubleMatrix& y);
  DoubleMatrix& operator -= (const DoubleMatrix& y);
  DoubleMatrix& operator *= (const DoubleMatrix& y);

/****************************
 * The following operators      *
 * modify the DoubleMatrix      *
 * by a double.                 *
 ****************************/
  DoubleMatrix& operator =  (double y);
  DoubleMatrix& operator += (double y);
  DoubleMatrix& operator -= (double y);
  DoubleMatrix& operator *= (double y);
  DoubleMatrix& operator /= (double y);

};
//
// The following are nonmember functions taking two arguments
//
/************** Won't compile on NT *********************
  ostream&  operator << (ostream& s, DoubleMatrix& x);                  // Outputs a Matrix
 *******************************************************/

/****************************
 * The following operators      *
 * take 1 or 2 DoubleMatrixs    *
 * and return a third           *
 ****************************/
  DoubleMatrix  operator - (const DoubleMatrix& x);                             // Negation operator
  DoubleMatrix  operator + (const DoubleMatrix& x, const DoubleMatrix& y);      // Adds two DoubleMatrixs
  DoubleMatrix  operator - (const DoubleMatrix& x, const DoubleMatrix& y);      // Subtracts two DoubleMatrixs
  DoubleMatrix  operator * (const DoubleMatrix& x, const DoubleMatrix& y);      // Multiplies two DoubleMatrixs


/********************************
 * The following operators      *
 * modify a double by a         *
 * DoubleMatrix, and return a   *
 * DoubleMatrix.                *
 ********************************/
  inline DoubleMatrix   operator + (double x, const DoubleMatrix &y);
  inline DoubleMatrix   operator - (double x, const DoubleMatrix &y);
  inline DoubleMatrix   operator * (double x, const DoubleMatrix &y);

  inline DoubleMatrix   operator + (const DoubleMatrix &x, double y);
  inline DoubleMatrix   operator - (const DoubleMatrix &x, double y);
  inline DoubleMatrix   operator * (const DoubleMatrix &x, double y);
  inline DoubleMatrix   operator / (const DoubleMatrix &x, double y);


/************************
 * Inline functions     *
 * that are defined     *
 * above.               *
 ************************/
  inline int    DoubleMatrix::getRows()                 const   {return _rows;}                 // Get the number of rows in the matrix
  inline int    DoubleMatrix::getColumns()              const   {return _columns;}              // Get the number of columns
  inline const double   *DoubleMatrix::getMatrix()      const   {return _matrix;}               // Return a pointer to the matrix

/****************************************
 * Inline operator functions defined    *
 * above.                               *
 ****************************************/
inline DoubleMatrix  operator + (double x, const DoubleMatrix &y)
{
  DoubleMatrix  sum = y;
  sum   += x;                   //Stroustrup p. 302
  return sum;
}
inline DoubleMatrix  operator - (double x, const DoubleMatrix &y)
{
  DoubleMatrix  result = -y;
  result        += x;                   //Stroustrup p. 302
  return result;
}
inline DoubleMatrix  operator * (double x, const DoubleMatrix &y)
{
  DoubleMatrix  product = y;
  product       *= x;                   //Stroustrup p. 302
  return product;
}
inline DoubleMatrix  operator + (const DoubleMatrix &x, double y)
{
  DoubleMatrix  sum = x;
  sum   += y;                   //Stroustrup p. 302
  return sum;
}
inline DoubleMatrix  operator - (const DoubleMatrix &x, double y)
{
  DoubleMatrix  difference = x;
  difference    -= y;                   //Stroustrup p. 302
  return difference;
}
inline DoubleMatrix  operator * (const DoubleMatrix &x, double y)
{
  DoubleMatrix  product = x;
  product       *= y;                   //Stroustrup p. 302
  return product;
}
inline DoubleMatrix  operator / (const DoubleMatrix &x, double y)
{
  DoubleMatrix  dividend = x;
  dividend      /= y;                   //Stroustrup p. 302
  return dividend;
}


#endif
