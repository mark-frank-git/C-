#ifndef _FLOATVECTOR_H
#define _FLOATVECTOR_H  1
/************************************************************************
 *                                                                      *
 * This class implements a floating point vector object.                *
 * coefficients of the FloatVector are single precision variables.      *
 * The usual operations are defined for this class, +, x, -, /, +=, etc.*
 * This class should be useful for digital signal processing.           *
 *                                                                      *
 * File:FloatVector.h                                                   *
 *                                                                      *
 *                                                                      *
 * NOTES:                                                               *
 *   1. vector = a[0] a[1]  ....  a[n]                                  *
 *                                                                      *
 ************************************************************************/

#ifdef  USE_IO_STREAM                           // Won't compile on NT *******************************
#include <fstream.h>
#endif                                          // ***************************************************

#ifndef NULL
#define NULL    0
#endif
#ifndef LOGICAL
#define LOGICAL    char
#endif


#define MAX_SHIFT           1024                    // Maximum shift right or left
#define MAX_STRING_SIZE     2048                    // FloatVector string size

class FloatVector
{
private:
  int       _size;                                  // Reported length of the vector
  int       _allocatedSize;                         // actual length of the _coefficients array
  int       _tempSize;                              // Size of the temp coefficients array
  float     *_coefficients;                         // array containing the FloatVector coefficients
  float     *_tempCoefficients;                     // temporary array for storing coefficients
//
// Private methods:
//
  void updateMemory(int newLength);

public:
//
// Public methods:
//
/*********************************
 * Constructors, destructors:    *
 *********************************/
  FloatVector(int size=0, float *coeff=NULL);               // Initializes a FloatVector
  FloatVector(const FloatVector &vector);
  FloatVector(float x);
  ~FloatVector();

/**********************
 * Set parameters:    *
 **********************/
  void  assign(int size, const float *coeff);               // Sets a new FloatVector
  void  fill(int size, float value);                        // Sets all elements of vector = value

/**********************
 * Get parameters:    *
 **********************/
  const float *coefficients() const  {return _coefficients;} // Return the FloatVector coefficients
  int   size() const                 {return _size;}         // Return the length of the FloatVector

/****************************
 * Modifying the vector     *
 ****************************/
  void  reverse();                                      // Reverse the order of coefficients of the vector
  void  normalize();                                    // Make the vector have unity norm
  void  chop(int n);                                    // Chop off n points from end of vector

/****************************
 * Calculating properties   *
 * of the vector:           *
 ****************************/
  float  norm() const;                             // Returns sqrt(sum(coeffs^2))
  float  max() const;                              // Returns the maximum element
  float  min() const;                              // Returns the minimum element
  float  mean() const;                             // Returns the average value
  float  get(int index) const;                     // Returns the element at the input index
  int    maxIndex() const;                         // Returns the index of the maximum element
  int    minIndex() const;                         // Returns the index of the minimum element


/****************************
 * The following operators  *
 * modify the FloatVector   *
 * by a second FloatVector  *
 ****************************/
  FloatVector&  operator =  (const FloatVector& y);
  FloatVector&  operator += (const FloatVector& y);
  FloatVector&  operator -= (const FloatVector& y);
  FloatVector&  operator *= (const FloatVector& y);
  FloatVector&  operator /= (const FloatVector& y);
  FloatVector&  append(const FloatVector& y);
  FloatVector&  prepend(const FloatVector& y);

/****************************
 * The following operators  *
 * modify the FloatVector   *
 * by a float.              *
 ****************************/
  FloatVector&  operator =  (float y);
  FloatVector&  operator += (float y);
  FloatVector&  operator -= (float y);
  FloatVector&  operator *= (float y);
  FloatVector&  operator /= (float y);
  FloatVector&  operator >= (float y);
  FloatVector&  operator <= (float y);

/****************************
 * The following operators  *
 * modify the FloatVector   *
 * by an int.               *
 ****************************/
  FloatVector&  operator >>= (int numberShifts);
  FloatVector&  operator <<= (int numberShifts);

};
//
// The following are nonmember functions
//

/****************************
 * The following operators  *
 * are for printing the     *
 * FloatVector string.      *
 ****************************/

#ifdef  USE_IO_STREAM                           // Won't compile on NT *******************************
ostream&  operator << (ostream& s, FloatVector& x);       // Outputs a FloatVector
#endif                                          // ***************************************************

/****************************
 * The following operators  *
 * make comparisons between *
 * two FloatVectors.            *
 ****************************/
  LOGICAL  operator == (const FloatVector& x, const FloatVector& y);       // Checks if two FloatVectors are equal
  LOGICAL  operator != (const FloatVector& x, const FloatVector& y);       // Checks if two FloatVectors are not equal

/****************************
 * The following operators  *
 * take 1 or 2 FloatVectors *
 * and return a third       *
 ****************************/
  FloatVector   operator - (const FloatVector& x);                          // Negation operator
  FloatVector   operator >> (const FloatVector& x, int numberShift);        // Shift right, 0 fill
  FloatVector   operator << (const FloatVector& x, int numberShift);        // Shift left, 0 fill
  FloatVector   operator + (const FloatVector& x, const FloatVector& y);    // Adds two FloatVectors
  FloatVector   operator - (const FloatVector& x, const FloatVector& y);    // Subtracts two FloatVectors
  FloatVector   operator * (const FloatVector& x, const FloatVector& y);    // Multiply two FloatVectors element wise
  FloatVector   operator / (const FloatVector& x, const FloatVector& y);    // Divides two FloatVectors
  FloatVector   convolve(const FloatVector& x, const FloatVector& y);       // Convolves the two vectors
  FloatVector   decimateVector(const FloatVector& x, int m);                // Decimate a vector by m
  FloatVector   interpolateVector(const FloatVector& x, int m);             // Interpolate a vector by m
  FloatVector   diff(const FloatVector& x);                                 // Returns first order difference of elements
  FloatVector   find(const FloatVector& x);                                 // Returns a vector of indices of non-zero
                                                                            // elements


/********************************
 * The following operators      *
 * modify a float by a          *
 * FloatVector, and return a    *
 * FloatVector.                 *
 ********************************/
  FloatVector   operator + (float x, const FloatVector &y);
  FloatVector   operator - (float x, const FloatVector &y);
  FloatVector   operator * (float x, const FloatVector &y);
  FloatVector   operator / (float x, const FloatVector &y);

/********************************
 * The following operators      *
 * modify a FloatVector by a    *
 * float, and return a          *
 * FloatVector.                 *
 ********************************/
  FloatVector   operator + (const FloatVector &x, float y);
  FloatVector   operator - (const FloatVector &x, float y);
  FloatVector   operator * (const FloatVector &x, float y);
  FloatVector   operator / (const FloatVector &x, float y);

/********************************
 * The following operators      *
 * return a FloatVector with    *
 * elements of 0 or 1 depending *
 * on the result of a           *
 * comparison with a float.     *
 ********************************/
  FloatVector   operator > (const FloatVector &x, float y);
  FloatVector   operator < (const FloatVector &x, float y);




#endif


