#ifndef _INTVECTOR_H
#define _INTVECTOR_H  1
/************************************************************************
 *                                                                      *
 * This class implements a crude version of the STL <int>vector         *
 *                                                                      *
 * File:IntVector.h                                                     *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. April 4, 2015 - Started                                          *
 *                                                                      *
 * NOTES:                                                               *
 *   1. vector = a[0] a[1]  ....  a[n]                                  *
 *                                                                      *
 ************************************************************************/

class IntVector
{
private:
  int       _size;                                  // Reported length of the vector
  int       _allocatedSize;                         // actual length of the _coefficients array
  int       _tempSize;                              // Size of the temp coefficients array
  int      *_coefficients;                          // array containing the IntVector coefficients
  int      *_tempCoefficients;                      // temporary array for storing coefficients
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
  IntVector(int size=0, int *coeff=NULL);               // Initializes a IntVector
  IntVector(const IntVector &vector);
  IntVector(int x);
  ~IntVector();

/**********************
 * Set parameters:    *
 **********************/
  void  assign(int size, const int *coeff);             // Sets a new IntVector
  void  fill(int size, int value);                      // Sets all elements of vector = value
  void  push_back(const int value);                     // Puts a new element in vector
  void  clear() {_size=0;}                              // Clear the vector

/**********************
 * Get parameters:    *
 **********************/
  const int *coefficients() const  {return _coefficients;} // Return the IntVector coefficients
  const int find(const int value) const;                // Finds value in the vector
  const int size()                const {return _size;} // Return the length of the IntVector

  
/****************************
 * The following operators  *
 * modify the IntVector   *
 * by a second IntVector  *
 ****************************/
  IntVector&  operator =  (const IntVector& y);

};

#endif

