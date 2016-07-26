#ifndef _COMPLEXMATRIX_H
#define _COMPLEXMATRIX_H
/************************************************************************
 *                                                                      *
 * This subclass of object implements a complex matrix.                 *
 * The matrix is represented as a Complex array, with the elements      *
 * increasing by rows.                                                  *
 * The usual overloaded operators are included: +, -, *, /, etc.        *
 *                                                                      *
 * File:ComplexMatrix.h                                                 *
 *                                                                      *
 ************************************************************************/
/////// Won't compile on NT #include <fstream.h>

#ifndef LOGICAL
#define LOGICAL char
#endif
#ifndef YES
#define YES             1
#define NO              0
#endif

#define MAX_CMPLX_STRING_SIZE   32768                   // Matrix string size

class           Complex;                                // Class description

class ComplexMatrix
{
protected:
  int           _rows, _columns;                        // # of rows, cols in matrix
  Complex       *_matrix;                               // matrix elements
  double        *_singularValues;                       // Singulare values
  char          *_matrixString;                         // For printing out the matrix

//
// Private methods:
//

public:
/*****************************
 * Constructor/Destructor    *
 *****************************/
  ComplexMatrix(int numberRows=1, int numberCols=1, const Complex *elements=NULL);
  ComplexMatrix(int numberRows, int numberCols, Complex element);       // Set all elements equal to input element
  ComplexMatrix(const ComplexMatrix &matrix);                           // Copy constructor
  virtual ~ComplexMatrix();

/**********************
 * Set parameters:    *
 **********************/
  virtual void assign(int numberRows=0, int numberCols=0, const Complex *inputMatrix=NULL);
                                                                                // Sets a new ComplexMatrix

/**********************
 * Get parameters:    *
 **********************/
  int getRows()                                 const;                  // Get # of rows in the matrix
  int getColumns()                              const;                  // Get the number of columns
  const Complex *getMatrix()                    const;                  // Return a pointer to the matrix

/**********************
 * Calculations:      *
 **********************/
  const double  *singularValues();
  double        cond();
  ComplexMatrix & ComplexMatrix::toeplitzSolve(Complex *r, Complex *b, int n, double *conda, int *flag);

/****************************
 * The following operators      *
 * modify the ComplexMatrix     *
 ****************************/
  virtual ComplexMatrix&        transpose();                            // Transpose operator
  virtual ComplexMatrix&        hermitian();                            // Conjugate transpose
  ComplexMatrix&                reverseElements();                      // Reverse elements in Matrix.  This probably
                                                                        // only makes sense for a vector
  ComplexMatrix&        replaceRow(int rowNumber, const ComplexMatrix& newRow);
  ComplexMatrix&        postTimesDiagonal(ComplexMatrix &diagonalMatrix); // post mult by a diagonal matrix

/********************************
 * The following functions      *
 * aid in printing the          *
 * matrix                       *
 ********************************/
  const char *getMatrixString();

/****************************
 * The following operators      *
 * modify the ComplexMatrix     *
 * by a second ComplexMatrix    *
 ****************************/
  ComplexMatrix&        operator =  (const ComplexMatrix& y);
  ComplexMatrix&        operator += (const ComplexMatrix& y);
  ComplexMatrix&        operator -= (const ComplexMatrix& y);
  ComplexMatrix&        operator *= (const ComplexMatrix& y);

/****************************
 * The following operators      *
 * modify the ComplexMatrix     *
 * by a Complex.                *
 ****************************/
  ComplexMatrix&        operator =  (Complex y);
  ComplexMatrix&        operator += (Complex y);
  ComplexMatrix&        operator -= (Complex y);
  ComplexMatrix&        operator *= (Complex y);
  ComplexMatrix&        operator /= (Complex y);

/****************************
 * The following operators      *
 * modify the ComplexMatrix     *
 * by a integer.                *
 ****************************/
  ComplexMatrix&        operator <<=  (int shift);

};
//
// The following are nonmember functions taking two arguments
//

/****************************
 * The following operators      *
 * are for printing the         *
 * matrix.                      *
 ****************************/
/***************** Won't compile on NT *********************************
  ostream&  operator << (ostream& s, ComplexMatrix& x);                         // Outputs a Matrix
 ***********************************************************************/

/****************************
 * The following operators      *
 * take 1 or 2 ComplexMatrices  *
 * and return a third           *
 ****************************/
  ComplexMatrix operator - (const ComplexMatrix& x);                            // Negation operator
  ComplexMatrix operator + (const ComplexMatrix& x, const ComplexMatrix& y);    // Adds two ComplexMatrices
  ComplexMatrix operator - (const ComplexMatrix& x, const ComplexMatrix& y);    // Subtracts two ComplexMatrices
  ComplexMatrix operator * (const ComplexMatrix& x, const ComplexMatrix& y);    // Multiplies two ComplexMatrices

  ComplexMatrix appendColumns(const ComplexMatrix& x, const ComplexMatrix& y);  // Appends a matrix to a matrix
  ComplexMatrix appendRows(const ComplexMatrix& x, const ComplexMatrix& y);     // Appends a matrix to a matrix, rows


/********************************
 * The following operators      *
 * modify a ComplexMatrix by a  *
 * Complex, and return a        *
 * ComplexMatrix.               *
 ********************************/
  inline ComplexMatrix operator + (Complex x, const ComplexMatrix&  y);
  inline ComplexMatrix operator - (Complex x, const ComplexMatrix&  y);
  inline ComplexMatrix operator * (Complex x, const ComplexMatrix&  y);

  inline ComplexMatrix operator + (const ComplexMatrix&  x, Complex y);
  inline ComplexMatrix operator - (const ComplexMatrix&  x, Complex y);
  inline ComplexMatrix operator * (const ComplexMatrix&  x, Complex y);
  inline ComplexMatrix operator / (const ComplexMatrix&  x, Complex y);


/********************************
 * The following operators      *
 * modify a ComplexMatrix by a  *
 * integer, and return a        *
 * ComplexMatrix.               *
 ********************************/
  inline ComplexMatrix operator << (const ComplexMatrix&  y, int shift);

/************************
 * Inline functions     *
 * that are defined     *
 * above.               *
 ************************/
  inline int    ComplexMatrix::getRows()                const   {return _rows;}         // Get the number of rows in the matrix
  inline int    ComplexMatrix::getColumns()             const   {return _columns;}      // Get the number of columns
  inline const  Complex *ComplexMatrix::getMatrix()     const   {return _matrix;}       // Return a pointer to the matrix

/****************************************
 * Inline operator functions defined    *
 * above.                               *
 ****************************************/
inline ComplexMatrix operator + (Complex x, const ComplexMatrix &y)
{
  ComplexMatrix sum = y;
  sum   += x;                                   //Stroustrup p. 302
  return sum;
}
inline ComplexMatrix operator - (Complex x, const ComplexMatrix &y)
{
  ComplexMatrix result = -y;
  result        += x;                           //Stroustrup p. 302
  return result;
}
inline ComplexMatrix operator * (Complex x, const ComplexMatrix &y)
{
  ComplexMatrix product = y;
  product       *= x;                           //Stroustrup p. 302
  return product;
}
inline ComplexMatrix operator + (const ComplexMatrix &x, Complex y)
{
  ComplexMatrix sum = x;
  sum   += y;                                   //Stroustrup p. 302
  return sum;
}
inline ComplexMatrix operator - (const ComplexMatrix &x, Complex y)
{
  ComplexMatrix difference = x;
  difference    -= y;                           //Stroustrup p. 302
  return difference;
}
inline ComplexMatrix operator * (const ComplexMatrix &x, Complex y)
{
  ComplexMatrix product = x;
  product       *= y;                           //Stroustrup p. 302
  return product;
}
inline ComplexMatrix  operator / (const ComplexMatrix &x, Complex y)
{
  ComplexMatrix dividend = x;
  dividend      /= y;                   //Stroustrup p. 302
  return dividend;
}

inline ComplexMatrix operator << (const ComplexMatrix &x, int shift)
{
  ComplexMatrix dividend = x;
  dividend      <<= shift;                      //Stroustrup p. 302
  return dividend;
}


#endif
