/************************************************************************
 *                                                                      *
 * This class implements a crude version of the STL <int>vector         *
 *                                                                      *
 * File:IntVector.h                                                     *
 *                                                                      *
 *                                                                      *
 * NOTES:                                                               *
 *   1. vector = a[0] a[1]  ....  a[n]                                  *
 *                                                                      *
 ************************************************************************/


#include <stdio.h>
#include <string.h>
#include "IntVector.h"                                                // Class prototypes


#define MAX(a,b)  ( ((a)<(b)) ? (b) : (a) )

#define GROW_SIZE   10
// ############################# Private Method ###############################
// updateMemory -- Method to perform common memory check/adjustment.
// Input:           newSize         New length of vector
//          
// Output:                          None
//
// Notes:
// ############################# Private Method ###############################
void IntVector::updateMemory (int newSize)
{
  int       old_size;
  int      *temp;                          // Pointer to old memory

  _size         = newSize;
  old_size      = _allocatedSize;           // Save old size
  if(_allocatedSize >= _size)
    return;                                 // No need to reallocate
  temp              = _coefficients;
  _allocatedSize    = newSize + GROW_SIZE;
  _coefficients     = new int[_allocatedSize];
  for(int i=0; i<old_size; i++)             // Copy old coefficients
    _coefficients[i] = temp[i];
  delete [] temp;
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the IntVector class.
// Input:           size:           length of the IntVector
//                  coeff:          IntVector coefficients
//          
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
IntVector::IntVector(int size, int *coeff)
{
// init instance variables:
  _allocatedSize        = _size     = MAX(1, size);
  _coefficients         = new int[_size];
  _tempCoefficients     = NULL;
  _tempSize             = 0;
  for(int i=0; i<_size; i++)
  {
    if(coeff == NULL)
      _coefficients[i]  = 1;
    else
      _coefficients[i]  = coeff[i];
  }
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Copy constructor for the IntVector class.
// Input:           vector:         A previously allocated IntVector
//          
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
IntVector::IntVector(const IntVector& vector)
{
// init instance variables:
  _size                   = vector.size();
  _size                   = MAX(1, _size);
  _allocatedSize          = _size;
  _coefficients           = new int[_size];
  _tempCoefficients       = NULL;
  _tempSize               = 0;
  const int *input_coeff  = vector.coefficients();
  for(int i=0; i<_size; i++)
  {
    _coefficients[i]      = input_coeff[i];
  }
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Copy constructor for the IntVector class, for input float
// Input:           x:              float input
//          
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
IntVector::IntVector(int x)
{
// init instance variables:
  _size                 = _allocatedSize    = 1;
  _coefficients         = new int[_size];
  _coefficients[0]      = x;
  _tempCoefficients     = NULL;
  _tempSize             = 0;
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the IntVector class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
IntVector::~IntVector()
{
  delete [] _coefficients;
  delete [] _tempCoefficients;
  return;
}

// ############################# Public Method ###############################
// assign -- Sets the IntVector's _coefficients to the input coefficient array
// Input:           size:           new size of the vector
//                  coeff:          coefficients to set vector equal to
//          
// Output:                      None
// ############################# Public Method ###############################
void IntVector::assign(int size, const int *coeff)
{
  updateMemory(size);
  
  for(int i=0; i<_size; i++)
  {
    if(coeff == NULL)
      _coefficients[i]  = 1;
    else
      _coefficients[i]  = coeff[i];
  }
  return;
}

// ############################# Public Method ###############################
// fill -- Sets the IntVector's coefficients to the input float value
// Input:           size:           new size of the vector
//                  value:          float value to set all elements
//          
// Output:                      None
// ############################# Public Method ###############################
void IntVector::fill(int size, int value)
{
  updateMemory(size);
  
  for(int i=0; i<_size; i++)
  {
    _coefficients[i]    = value;
  }
  return;
}

// ############################# Public Method ###############################
// push_back -- Adds an element to the vector.
// Input:       value:          Element to add
//          
// Output:                      None
//
// Notes:
//  1. The instance variables, _coefficients[] are modified.
// ############################# Public Method ###############################
void IntVector::push_back(const int value)
{
  int oldSize   = _size;
  updateMemory(_size+1);
  _coefficients[oldSize]        = value;
  return;
}

// ############################# Public Method ###############################
// normalize -- Normalizes the vector to have unity norm.
// Input:                       None
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
const int IntVector::find(const int value) const
{
  for(int i=0; i<_size; i++)
  {
    if(_coefficients[i] == value)
      return i;
  }
  return -1;
}

// ############################# Public Method ###############################
// Operator =  Sets the IntVector equal to the input IntVector
// Input:           y:          input IntVector
//          
// Output:                      the result IntVector
//
// ############################# Public Method ###############################
IntVector& IntVector::operator = (const IntVector& y)
{
  const int *input_coeff;

  if( this == &y )                          // Check for x=x
    return *this;
  input_coeff   = y.coefficients();
  updateMemory(y.size());
  for(int i=0; i<_size; i++)
    _coefficients[i] = input_coeff[i];
  return *this;
}

