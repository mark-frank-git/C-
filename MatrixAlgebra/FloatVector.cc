/************************************************************************
 *                                                                      *
 * This class implements a floating point vector object.                *
 * coefficients of the FloatVector are single precision variables.      *
 * The usual operations are defined for this class, +, x, -, /, +=, etc.*
 * This class should be useful for digital signal processing.           *
 *                                                                      *
 * File:FloatVector.cc                                                  *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 03/27/97  - Started                                              *
 *                                                                      *
 * NOTES:                                                               *
 *   1. vector = a[0] a[1]  ....  a[n]                                  *
 ************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "FloatVector.h"                                                // Class prototypes


#define MAX(a,b)  ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))
#ifndef YES
#define YES 1
#define NO  0
#endif

#define GROW_SIZE   10
// ############################# Private Method ###############################
// updateMemory -- Method to perform common memory check/adjustment.
// Input:           newSize         New length of vector
//          
// Output:                          None
//
// Notes:
// ############################# Private Method ###############################
void FloatVector::updateMemory (int newSize)
{
  int       old_size;
  float     *temp;                          // Pointer to old memory

  _size         = newSize;
  old_size      = _allocatedSize;           // Save old size
  if(_allocatedSize >= _size)
    return;                                 // No need to reallocate
  temp              = _coefficients;
  _allocatedSize    = newSize + GROW_SIZE;
  _coefficients     = new float[_allocatedSize];
  for(int i=0; i<old_size; i++)             // Copy old coefficients
    _coefficients[i] = temp[i];
  delete [] temp;
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the FloatVector class.
// Input:           size:           length of the FloatVector
//                  coeff:          FloatVector coefficients
//          
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
FloatVector::FloatVector(int size, float *coeff)
{
  int i;
// init instance variables:
  _allocatedSize        = _size     = MAX(1, size);
  _coefficients         = new float[_size];
  _tempCoefficients     = NULL;
  _tempSize             = 0;
  for(i=0; i<_size; i++)
  {
    if(coeff == NULL)
      _coefficients[i]  = 1;
    else
      _coefficients[i]  = coeff[i];
  }
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Copy constructor for the FloatVector class.
// Input:           vector:         A previously allocated FloatVector
//          
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
FloatVector::FloatVector(const FloatVector& vector)
{
  int   i;
  const float *input_coeff;

// init instance variables:
  _size                 = vector.size();
  _size                 = MAX(1, _size);
  _allocatedSize        = _size;
  _coefficients         = new float[_size];
  _tempCoefficients     = NULL;
  _tempSize             = 0;
  input_coeff           = vector.coefficients();
  for(i=0; i<_size; i++)
    _coefficients[i]    = input_coeff[i];
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Copy constructor for the FloatVector class, for input float
// Input:           x:              float input
//          
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
FloatVector::FloatVector(float x)
{

// init instance variables:
  _size                 = _allocatedSize    = 1;
  _coefficients         = new float[_size];
  _coefficients[0]      = x;
  _tempCoefficients     = NULL;
  _tempSize             = 0;
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the FloatVector class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
FloatVector::~FloatVector()
{
  delete [] _coefficients;
  delete [] _tempCoefficients;
  return;
}

// ############################# Public Method ###############################
// assign -- Sets the FloatVector's _coefficients to the input coefficient array
// Input:           size:           new size of the vector
//                  coeff:          coefficients to set vector equal to
//          
// Output:                      None
// ############################# Public Method ###############################
void FloatVector::assign(int size, const float *coeff)
{
  int i;

  updateMemory(size);
  
  for(i=0; i<_size; i++)
  {
    if(coeff == NULL)
      _coefficients[i]  = 1;
    else
      _coefficients[i]  = coeff[i];
  }
  return;
}

// ############################# Public Method ###############################
// fill -- Sets the FloatVector's coefficients to the input float value
// Input:           size:           new size of the vector
//                  value:          float value to set all elements
//          
// Output:                      None
// ############################# Public Method ###############################
void FloatVector::fill(int size, float value)
{
  int i;

  updateMemory(size);
  
  for(i=0; i<_size; i++)
  {
    _coefficients[i]    = value;
  }
  return;
}

// ############################# Public Method ###############################
// reverse -- Reverse the FloatVector.
// Input:                       None
//          
// Output:                      None
//
// Notes:
//  1. The instance variables, _coefficients[] are modified.
// ############################# Public Method ###############################
void FloatVector::reverse()
{
  int   i, n;
  float temp;
  
  n = _size - 1;
  for(i=0; i<_size/2; i++)
  {
    temp                = _coefficients[i];
    _coefficients[i]    = _coefficients[n];
    _coefficients[n--]  = temp;
  }
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
void FloatVector::normalize()
{
  int       i;
  float     vec_norm;
  
  vec_norm  = norm();
  if(vec_norm > 0.)
  {
    for(i=0; i<_size; i++)
    {
      _coefficients[i]  /= vec_norm;
    }
  }
  return;
}

// ############################# Public Method ###############################
// chop -- Chops off n points from end of vector.
// Input:       n:              Number of points to chop off
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
void FloatVector::chop(int n)
{
  _size     -= n;
  _size     = MAX(1, _size);
  return;
}

// ############################# Public Method ###############################
// norm -- Returns the L2 norm of the vector.
// Input:                       None
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
float FloatVector::norm() const
{
  int       i;
  double    sum;
  
  sum       = 0.;
  for(i=0; i<_size; i++)
  {
    sum     += _coefficients[i]*_coefficients[i];
  }
  return sqrt(sum);
}

// ############################# Public Method ###############################
// max -- Returns the maximum element of the vector.
// Input:                       None
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
float FloatVector::max() const
{
  int       i;
  float     max_element;

  max_element   = _coefficients[0];
  for(i=1; i<_size; i++)
  {
    max_element = MAX(_coefficients[i], max_element);
  }
  return max_element;
}

// ############################# Public Method ###############################
// min -- Returns the minimum element of the vector.
// Input:                       None
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
float FloatVector::min() const
{
  int       i;
  float     min_element;

  min_element   = _coefficients[0];
  for(i=1; i<_size; i++)
  {
    min_element = MIN(_coefficients[i], min_element);
  }
  return min_element;
}

// ############################# Public Method ###############################
// mean -- Returns the average value of all the elements.
// Input:                       None
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
float FloatVector::mean() const
{
  int       i;
  double    sum;

  sum       = 0.;
  for(i=0; i<_size; i++)
  {
    sum     += _coefficients[i];
  }
  sum       /= _size;
  return (float) sum;
}

// ############################# Public Method ###############################
// get -- Returns the element at the given index
// Input:                       None
//          
// Output:                      None
//
// Notes:
// ############################# Public Method ###############################
float FloatVector::get(int index) const
{
//
// Error checking:
//
  index = MAX(0, index);
  index = MIN(index, _size-1);

  return _coefficients[index];
}

// ############################# Public Method ###############################
// maxIndex -- Returns the index of the maximum element of the vector.
// Input:                       None
//          
// Output:                      None
//
// Notes:
// 1. Indexing starts at 0.
// ############################# Public Method ###############################
int FloatVector::maxIndex() const
{
  int       i, max_index;
  float     max_element;

  max_element   = _coefficients[0];
  max_index     = 0;
  for(i=1; i<_size; i++)
  {
    if(_coefficients[i] > max_element)
    {
      max_element   =_coefficients[i];
      max_index     = i;
    }
  }
  return max_index;
}

// ############################# Public Method ###############################
// minIndex -- Returns the index of the minimum element of the vector.
// Input:                       None
//          
// Output:                      None
//
// Notes:
// 1. Indexing starts at 0.
// ############################# Public Method ###############################
int FloatVector::minIndex() const
{
  int       i, min_index;
  float     min_element;

  min_element   = _coefficients[0];
  min_index     = 0;
  for(i=1; i<_size; i++)
  {
    if(_coefficients[i] < min_element)
    {
      min_element   =_coefficients[i];
      min_index     = i;
    }
  }
  return min_index;
}

// ############################# Public Method ###############################
// Operator =  Sets the FloatVector equal to the input FloatVector
// Input:           y:          input FloatVector
//          
// Output:                      the result FloatVector
//
// ############################# Public Method ###############################
FloatVector& FloatVector::operator = (const FloatVector& y)
{
  const float *input_coeff;

  if( this == &y )                          // Check for x=x
    return *this;
  input_coeff   = y.coefficients();
  updateMemory(y.size());
  for(int i=0; i<_size; i++)
    _coefficients[i] = input_coeff[i];
  return *this;
}

// ############################# Public Method ###############################
// Operator +=  Add the input FloatVector to this.
// Input:           y:          input FloatVector
//          
// Output:                      the result FloatVector
//
// NOTES:
// ############################# Public Method ###############################
FloatVector& FloatVector::operator += (const FloatVector& y)
{
  int   i, y_size, old_size, min_size;
  const float *input_coeff;

  y_size        = y.size();
  old_size      = _size;
  min_size      = MIN(y_size, old_size);
  _size         = MAX(_size, y_size);
  updateMemory(_size);
  input_coeff   = y.coefficients();
  for(i=0; i<min_size; i++)
    _coefficients[i]    += input_coeff[i];
  if(y_size > old_size)
  {
    for(i=old_size; i<y_size; i++)
    {
      _coefficients[i]  = input_coeff[i];
    }
  }
  return *this;
}

// ############################# Public Method ###############################
// Operator -=  Sub the input FloatVector from this.
// Input:           y:          input FloatVector
//          
// Output:                      the result FloatVector
//
// NOTES:
// ############################# Public Method ###############################
FloatVector& FloatVector::operator -= (const FloatVector& y)
{
  int   i, y_size, old_size, min_size;
  const float *input_coeff;

  y_size        = y.size();
  old_size      = _size;
  min_size      = MIN(y_size, old_size);
  _size         = MAX(_size, y_size);
  updateMemory(_size);
  input_coeff   = y.coefficients();
  for(i=0; i<min_size; i++)
    _coefficients[i]    -= input_coeff[i];
  if(y_size > old_size)
  {
    for(i=old_size; i<y_size; i++)
    {
      _coefficients[i]  = -input_coeff[i];
    }
  }
  return *this;
}

// ############################# Public Method ###############################
// Operator *=  Mult this by the input FloatVector.
// Input:           y:          input FloatVector
//          
// Output:                      the result FloatVector
//
// NOTES:
// 1. Multiplication is performed element-wise
// ############################# Public Method ###############################
FloatVector& FloatVector::operator *= (const FloatVector& y)
{
  int   i, min_size;
  const float *input_coeff;

  min_size      = MIN(_size, y.size());
  input_coeff   = y.coefficients();
  updateMemory(min_size);
  for(i=0; i<min_size; i++)
  {
    _coefficients[i]    *= input_coeff[i];
  }
  return *this;
}

// ############################# Public Method ###############################
// Operator /=  Divide this by the input FloatVector.
// Input:           y:          input FloatVector
//          
// Output:                      the result FloatVector
//
// NOTES:
// 1. Division is performed element-wise.
// ############################# Public Method ###############################
FloatVector& FloatVector::operator /= (const FloatVector& y)
{
  int   i, min_size;
  const float *input_coeff;

  min_size      = MIN(_size, y.size());
  input_coeff   = y.coefficients();
  updateMemory(min_size);
  for(i=0; i<min_size; i++)
  {
    if(input_coeff[i] != 0.)
      _coefficients[i]  /= input_coeff[i];
  }
  return *this;
}

// ############################# Public Method ###############################
// append:  Appends the input FloatVector onto the end of this
// Input:           y:          input FloatVector
//          
// Output:                      the result FloatVector
//
// NOTES:
// ############################# Public Method ###############################
FloatVector& FloatVector::append(const FloatVector& y)
{
  int   i, k, new_size;
  const float *y_coeff;
//
// Find the new size of the vector, and allocate memory for it:
//
  new_size      = _size + y.size();
  if(new_size > _tempSize)
  {
    delete [] _tempCoefficients;
    _tempSize           = new_size;
    _tempCoefficients   = new float[new_size];
  }
//
// Copy the old coefficients:
//
  for(i=0; i<_size; i++)
    _tempCoefficients[i]    = _coefficients[i];
//
// Copy the new:
//
  y_coeff   = y.coefficients();
  k         = _size;
  for(i=0; i<y.size(); i++)
    _tempCoefficients[k++]  = y_coeff[i];
//
// Store temp coefficients:
//
  assign(new_size, _tempCoefficients);

  return *this;
}

// ############################# Public Method ###############################
// prepend:  prepends the input FloatVector onto the beginning of this
// Input:           y:          input FloatVector
//          
// Output:                      the result FloatVector
//
// NOTES:
// ############################# Public Method ###############################
FloatVector& FloatVector::prepend(const FloatVector& y)
{
  int   i, k, new_size;
  const float *y_coeff;
//
// Find the new size of the vector, and allocate memory for it:
//
  new_size      = _size + y.size();
  if(new_size > _tempSize)
  {
    delete [] _tempCoefficients;
    _tempSize           = new_size;
    _tempCoefficients   = new float[new_size];
  }
//
// Copy the new:
//
  y_coeff   = y.coefficients();
  for(i=0; i<y.size(); i++)
    _tempCoefficients[i]    = y_coeff[i];
//
// Copy the old coefficients:
//
  k         = y.size();
  for(i=0; i<_size; i++)
    _tempCoefficients[k++]  = _coefficients[i];
//
// Store temp coefficients:
//
  assign(new_size, _tempCoefficients);

  return *this;
}


// ############################# Public Method ###############################
// Operator =  Sets the FloatVector equal to the input float
// Input:           y:          input float
//          
// Output:                      the result FloatVector
//
// ############################# Public Method ###############################
FloatVector& FloatVector::operator = (float y)
{

  updateMemory(1);
  _coefficients[0]  = y;
  return *this;
}

// ############################# Public Method ###############################
// Operator +=  Adds a float to the FloatVector.
// Input:           y:              float to add
//          
// Output:                          the result FloatVector
//
// NOTES:
//   1. vector = a[n] a[n-1]  ....  a[0]
//
// ############################# Public Method ###############################
FloatVector& FloatVector::operator += (float y)
{
  for(int i=0; i<_size; i++)
    _coefficients[i]    += y;
  return *this;
}

// ############################# Public Method ###############################
// Operator -=  Subtracts a float from the FloatVector.
// Input:           y:              float to subtract
//          
// Output:                          the result FloatVector
//
// NOTES:
//   1. vector = a[n] a[n-1]  ....  a[0]
//
// ############################# Public Method ###############################
FloatVector& FloatVector::operator -= (float y)
{
  for(int i=0; i<_size; i++)
    _coefficients[i]    -= y;
  return *this;
}

// ############################# Public Method ###############################
// Operator *=  Multiplies a float by the FloatVector.
// Input:           y:              float to mult
//          
// Output:                          the result FloatVector
//
// NOTES:
//   1. vector = a[n] a[n-1]  ....  a[0]
//
// ############################# Public Method ###############################
FloatVector& FloatVector::operator *= (float y)
{
  for(int i=0; i<_size; i++)
    _coefficients[i]    *= y;
  return *this;
}

// ############################# Public Method ###############################
// Operator /=  Divides the FloatVector by a float
// Input:           y:              float to divide into
//          
// Output:                          the result FloatVector
//
// NOTES:
//   1. vector = a[n] a[n-1]  ....  a[0]
//
// ############################# Public Method ###############################
FloatVector& FloatVector::operator /= (float y)
{
  if(y != 0.)
  {
    for(int i=0; i<_size; i++)
      _coefficients[i]  /= y;
  }
  return *this;
}

// ############################# Public Method ###############################
// Operator >=  Checks whether each element of the vector is greater than the
//              input float.  If greater, set the element to 1, else set to 0.
// Input:           y:              float to compare to
//          
// Output:                          the result FloatVector
//
// NOTES:
//   1. vector = a[n] a[n-1]  ....  a[0]
//
// ############################# Public Method ###############################
FloatVector& FloatVector::operator >= (float y)
{
  for(int i=0; i<_size; i++)
    _coefficients[i]    = (float)(_coefficients[i] > y);
  return *this;
}

// ############################# Public Method ###############################
// Operator <=  Checks whether each element of the vector is less than the
//              input float.  If greater, set the element to 1, else set to 0.
// Input:           y:              float to compare to
//          
// Output:                          the result FloatVector
//
// NOTES:
//   1. vector = a[n] a[n-1]  ....  a[0]
//
// ############################# Public Method ###############################
FloatVector& FloatVector::operator <= (float y)
{
  for(int i=0; i<_size; i++)
    _coefficients[i]    = (float)(_coefficients[i] < y);
  return *this;
}

// ############################# Public Method ###############################
// Operator >>=  Shift the FloatVector right by n shifts.
// Input:           numberShifts:   Number of places to shift right
//          
// Output:                          the result FloatVector
//
// NOTES:
//   1. vector = a[0] a[1]  ....  a[n]
//
// ############################# Public Method ###############################
FloatVector& FloatVector::operator >>= (int numberShifts)
{
  int   i, j, old_size;
//
// Bound the # of shifts, and update the order
//
  old_size      = _size;
  numberShifts  = MAX(0, numberShifts);
  numberShifts  = MIN(numberShifts, MAX_SHIFT);
  _size         += numberShifts;
  updateMemory(_size);
//
// Copy to the right:
//
  j = _size - 1;
  for(i=0; i<old_size; i++)
  {
    _coefficients[j]    = _coefficients[j-numberShifts];
    j--;
  }
//
// Now, zero fill:
//
  for(i=0; i<numberShifts; i++)
  {
    _coefficients[i]    = 0.;
  }
  return *this;
}

// ############################# Public Method ###############################
// Operator <<=  Shift the FloatVector left by n shifts.
// Input:           numberShifts:   Number of places to shift left
//          
// Output:                          the result FloatVector
//
// NOTES:
//   1. vector = a[0] a[1]  ....  a[n]
//
// ############################# Public Method ###############################
FloatVector& FloatVector::operator <<= (int numberShifts)
{
  int   i;

//
// Bound the # of shifts, and update the order
//
  numberShifts  = MAX(0, numberShifts);
  numberShifts  = MIN(numberShifts, MAX_SHIFT);
  _size         -= numberShifts;
  _size         = MAX(1, _size);
//
// Do the shifting by copying:
//
  for(i=0; i<_size; i++)
    _coefficients[i] = _coefficients[i+numberShifts];
  return *this;
}

#ifdef  USE_IO_STREAM                           // Won't compile on NT *******************************
// ############################# Public Method ###############################
// <<   - Print the FloatVector
// Input:           s:              ofstream for printing
//                  x:              the FloatVector to print
//          
// Output:                          None
//
// ############################# Public Method ###############################
ostream& operator << (ostream& s, const FloatVector& x)
{
  const float *coeff;

  coeff = x.coefficients();
  s << "[ ";
  for(int i=0; i<x.size(); i++)
    s << coeff[i] << " ";
  return s <<  "]\n";
}
#endif                                          // ***************************************************

// ############################# Public Method ###############################
// Operator ==  Checks whether the two FloatVectors are equal
// Input:           x:          first FloatVector
//                  y:          second FloatVector
//          
// Output:                      YES if two FloatVectors are equal, else NO
//
// NOTES:
// 1. This is probably not very useful for floating point _coefficients
// ############################# Public Method ###############################
LOGICAL operator == (const FloatVector& x, const FloatVector& y)
{
  int   i, n;
  const float *x_coeff, *y_coeff;
  LOGICAL  equal;
  
  if( (n = x.size()) != y.size())
    return NO;

  equal     = YES;
  x_coeff   = x.coefficients();
  y_coeff   = y.coefficients();
  for(i=0; i<n; i++)
  {
    if(x_coeff[i] != y_coeff[i])
    {
      equal = NO;
      break;
    }
  }
  return equal;
}

// ############################# Public Method ###############################
// Operator !=  Checks whether the two FloatVectors are not equal
// Input:           x:          first FloatVector
//                  y:          second FloatVector
//          
// Output:                      YES if two FloatVectors are not equal, else NO
//
// ############################# Public Method ###############################
LOGICAL operator != (const FloatVector& x, const FloatVector& y)
{
  return (!(y==x));
}

// ############################# Public Method ###############################
// Operator -  Returns the negative of the input FloatVector
// Input:           x:          input FloatVector
//          
// Output:                      the output FloatVector
//
// NOTES
// ############################# Public Method ###############################
FloatVector operator - (const FloatVector& x)
{
  return (-1.*x);
}

// ############################# Public Method ###############################
// Operator >>  Returns the right shift of the input FloatVector (0 fill)
// Input:           x:          input FloatVector
//                  numberShift # of places to shift left
//          
// Output:                      the output FloatVector
//
// NOTES:
// 1. The maximum shift is equal to MAX_SHIFT
// ############################# Public Method ###############################
FloatVector operator >> (const FloatVector& x, int numberShift)
{
  FloatVector   shift = x;
  shift >>= numberShift;                            //Stroustrup p. 302
  return shift;
}

// ############################# Public Method ###############################
// Operator <<  Returns the left shift of the input FloatVector (0 fill)
// Input:           x:          input FloatVector
//                  numberShift # of places to shift left
//          
// Output:                      the output FloatVector
//
// NOTES:
// 1. The maximum shift is equal to MAX_SHIFT
// ############################# Public Method ###############################
FloatVector operator << (const FloatVector& x, int numberShift)
{
  FloatVector   shift = x;
  shift <<= numberShift;                            //Stroustrup p. 302
  return shift;
}

// ############################# Public Method ###############################
// Operator +  Returns the sum of the two input FloatVectors
// Input:           x:          input FloatVector
//                  y:          second FloatVector
//          
// Output:                      the sum FloatVector
//
// ############################# Public Method ###############################
FloatVector operator + (const FloatVector& x, const FloatVector& y)
{
  FloatVector   sum = x;
  sum   += y;                           //Stroustrup p. 302
  return sum;
}

// ############################# Public Method ###############################
// Operator -  Returns the difference of the two input FloatVectors
// Input:           x:          input FloatVector
//                  y:          second FloatVector
//          
// Output:                      the difference FloatVector
//
// NOTES:
// 1. For modulo-2 addition, + equals -.
// ############################# Public Method ###############################
FloatVector operator - (const FloatVector& x, const FloatVector& y)
{
  FloatVector   difference = x;
  difference    -= y;                   //Stroustrup p. 302
  return difference;
}

// ############################# Public Method ###############################
// Operator *  Returns the product of the two input FloatVectors
// Input:           x:          input FloatVector
//                  y:          second FloatVector
//          
// Output:                      the product FloatVector
//
// ############################# Public Method ###############################
FloatVector operator * (const FloatVector& x, const FloatVector& y)
{
  FloatVector   prod = x;
  prod  *= y;                           //Stroustrup p. 302
  return prod;
}

// ############################# Public Method ###############################
// Operator /  Returns the division of the two input FloatVectors
// Input:           x:          input FloatVector
//                  y:          second FloatVector
//          
// Output:                      the result FloatVector
//
// ############################# Public Method ###############################
FloatVector operator / (const FloatVector& x, const FloatVector& y)
{
  FloatVector   result = x;
  result    /= y;                           //Stroustrup p. 302
  return    result;
}

// ############################# Public Method ###############################
// convolve  Convolves the two input vectors, and returns the result
// Input:           x:          input FloatVector
//                  y:          second FloatVector
//          
// Output:                      the result FloatVector
//
// NOTES:
// 1. The y vector is slid across the x vector.  Note that this convolution is
//    not commutative.
// ############################# Public Method ###############################
FloatVector convolve(const FloatVector& x, const FloatVector& y)
{
  int           k, n, max_index;
  int           x_size, y_size, output_size;
  static        int temp_size   = 0;
  static        float *temp_coeff = NULL;
  const         float *x_coeff, *y_coeff;
  double        sum;
  FloatVector   result;

  x_size        = x.size();
  y_size        = y.size();
  output_size   = x_size - y_size;
  output_size   = ABS( output_size ) + 1;
  x_coeff       = x.coefficients();
  y_coeff       = y.coefficients();
//
// Allocate output array
//
  if(output_size > temp_size)
  {
    delete [] temp_coeff;
    temp_size   = output_size;
    temp_coeff  = new float[output_size];
  }

  for(k=0; k<output_size; k++)
  {
    max_index       = MIN(y_size, x_size-k);
    sum             = 0.;
    for(n=0; n<max_index; n++)
      sum           += x_coeff[n+k]*y_coeff[n];
    temp_coeff[k]   = (float) sum;
  }

  result.assign(output_size, temp_coeff);
  return    result;
}

// ############################# Public Method ###############################
// decimateVector: Decimates the input vector by selecting every mth coefficient.
// Input:           x:          input FloatVector
//                  m:          the decimation factor
//          
// Output:                      the result FloatVector
//
// NOTES:
// ############################# Public Method ###############################
FloatVector decimateVector(const FloatVector& x, int m)
{
  int           i, k;
  int           x_size, output_size;
  const         float *x_coeff;
  static        int temp_size   = 0;
  static        float *temp_coeff = NULL;
  FloatVector   result;

  m             = MAX(1, m);
  x_size        = x.size();
  output_size   = x_size/m;
  x_coeff       = x.coefficients();
//
// Allocate output array
//
  if(output_size > temp_size)
  {
    delete [] temp_coeff;
    temp_size   = output_size;
    temp_coeff  = new float[output_size];
  }

  k             = 0;
  for(i=0; i<output_size; i++)
  {
    temp_coeff[i]   = x_coeff[k];
    k                       += m;
  }

  result.assign(output_size, temp_coeff);
  return    result;
}

// ############################# Public Method ###############################
// interpolateVector: Interpolates the input vector by inserting zeros.
// Input:           x:          input FloatVector
//                  m:          the interpolation factor
//          
// Output:                      the result FloatVector
//
// NOTES:
// ############################# Public Method ###############################
FloatVector interpolateVector(const FloatVector& x, int m)
{
  int           i, j, k;
  int           x_size, output_size;
  const         float *x_coeff;
  static        int temp_size   = 0;
  static        float *temp_coeff = NULL;
  FloatVector   result;

  m             = MAX(1, m);
  x_size        = x.size();
  output_size   = x_size*m;
  x_coeff       = x.coefficients();
//
// Allocate output array
//
  if(output_size > temp_size)
  {
    delete [] temp_coeff;
    temp_size   = output_size;
    temp_coeff  = new float[output_size];
  }

  k             = 0;
  for(i=0; i<x_size; i++)
  {
    temp_coeff[k++]     = x_coeff[i];
    for(j=1; j<m; j++)
      temp_coeff[k++]   = 0.;
  }

  result.assign(output_size, temp_coeff);
  return    result;
}

// ############################# Public Method ###############################
// diff: Returns the first order difference of the input vector.
// Input:           x:          input FloatVector
//          
// Output:                      the result FloatVector
//
// NOTES:
// 1. result = [ (x[1]-x[0])  (x[2]-x[1]) ... (x[n-1]-x[n-2]) ]
// ############################# Public Method ###############################
FloatVector diff(const FloatVector& x)
{
  int           k;
  int           x_size, output_size;
  const         float *x_coeff;
  static        int temp_size   = 0;
  static        float *temp_coeff = NULL;
  FloatVector   result;

  x_size        = x.size();
  output_size   = x_size - 1;
  output_size   = MAX(1, output_size);
  x_coeff       = x.coefficients();
//
// Allocate output array
//
  if(output_size > temp_size)
  {
    delete [] temp_coeff;
    temp_size   = output_size;
    temp_coeff  = new float[output_size];
  }

  for(k=1; k<x_size; k++)
  {
    temp_coeff[k-1] = x_coeff[k] - x_coeff[k-1];
  }

  result.assign(output_size, temp_coeff);
  return    result;
}

// ############################# Public Method ###############################
// find: Returns a vector of indices of the non-zero elements of input vector.
// Input:           x:          input FloatVector
//          
// Output:                      the result FloatVector
//
// NOTES:
// 1. Indices start at zero
// ############################# Public Method ###############################
FloatVector find(const FloatVector& x)
{
  int           k, number_non_zero;
  int           x_size, output_size;
  const         float *x_coeff;
  static        int temp_size   = 0;
  static        float *temp_coeff = NULL;
  FloatVector   result;

  x_size        = x.size();
  output_size   = x_size;
  output_size   = MAX(1, output_size);
  x_coeff       = x.coefficients();

//
// Allocate output array
//
  if(output_size > temp_size)
  {
    delete [] temp_coeff;
    temp_size   = output_size;
    temp_coeff  = new float[output_size];
  }
  number_non_zero   = 0;
  for(k=0; k<x_size; k++)
  {
    if(x_coeff[k] != 0.)
    {
      temp_coeff[number_non_zero++] = k;
    }
  }

  result.assign(number_non_zero, temp_coeff);
  return    result;
}

// ############################# Public Method ###############################
// Operator +  Returns the sum of a float and a FloatVector
// Input:           x:          input float
//                  y:          input FloatVector
//          
// Output:                      the sum FloatVector
//
// ############################# Public Method ###############################
FloatVector operator + (float x, const FloatVector &y)
{
  FloatVector   sum = y;
  sum   += x;                           //Stroustrup p. 302
  return sum;
}

// ############################# Public Method ###############################
// Operator -  Returns the difference of a float and a FloatVector
// Input:           x:          input float
//                  y:          input FloatVector
//          
// Output:                      the difference FloatVector
//
// ############################# Public Method ###############################
FloatVector operator - (float x, const FloatVector &y)
{
  FloatVector   result = -y;
  result    += x;                           //Stroustrup p. 302
  return result;
}

// ############################# Public Method ###############################
// Operator *  Returns the product of a float and a FloatVector
// Input:           x:          input float
//                  y:          input FloatVector
//          
// Output:                      the product FloatVector
//
// ############################# Public Method ###############################
FloatVector operator * (float x, const FloatVector &y)
{
  FloatVector   product = y;
  product   *= x;                           //Stroustrup p. 302
  return product;
}

// ############################# Public Method ###############################
// Operator /  Returns the division of a float by a FloatVector
// Input:           x:          input float
//                  y:          input FloatVector
//          
// Output:                      the sum FloatVector
//
// ############################# Public Method ###############################
FloatVector operator / (float x, const FloatVector &y)
{
  FloatVector   dividend = x;
  dividend  /= y;                           //Stroustrup p. 302
  return dividend;
}

// ############################# Public Method ###############################
// Operator +  Returns the sum of a FloatVector and a float
// Input:           x:          input FloatVector
//                  y:          input float
//          
// Output:                      the sum FloatVector
//
// ############################# Public Method ###############################
FloatVector operator + (const FloatVector &x, float y)
{
  FloatVector   sum = x;
  sum   += y;                           //Stroustrup p. 302
  return sum;
}


// ############################# Public Method ###############################
// Operator -  Returns the difference of a FloatVector and a float
// Input:           x:          input FloatVector
//                  y:          input float
//          
// Output:                      the difference FloatVector
//
// ############################# Public Method ###############################
FloatVector operator - (const FloatVector &x, float y)
{
  FloatVector   difference = x;
  difference    -= y;                           //Stroustrup p. 302
  return difference;
}


// ############################# Public Method ###############################
// Operator *  Returns the product of a FloatVector and a float
// Input:           x:          input FloatVector
//                  y:          input float
//          
// Output:                      the product FloatVector
//
// ############################# Public Method ###############################
FloatVector operator * (const FloatVector &x, float y)
{
  FloatVector   product = x;
  product   *= y;                           //Stroustrup p. 302
  return product;
}


// ############################# Public Method ###############################
// Operator /  Returns the dividend of a FloatVector by a float
// Input:           x:          input FloatVector
//                  y:          input float
//          
// Output:                      the dividend FloatVector
//
// ############################# Public Method ###############################
FloatVector operator / (const FloatVector &x, float y)
{
  FloatVector   dividend = x;
  dividend  /= y;                           //Stroustrup p. 302
  return dividend;
}

// ############################# Public Method ###############################
// Operator >  Returns a FloatVector whose elements are 1 if x[i] > y, or
//             0 otherwise
// Input:           x:          input FloatVector
//                  y:          input float
//          
// Output:                      the 0, 1 vector
//
// ############################# Public Method ###############################
FloatVector operator > (const FloatVector &x, float y)
{
  FloatVector   result = x;
  result    >= y;                           //Stroustrup p. 302
  return result;
}

// ############################# Public Method ###############################
// Operator <  Returns a FloatVector whose elements are 1 if x[i] < y, or
//             0 otherwise
// Input:           x:          input FloatVector
//                  y:          input float
//          
// Output:                      the 0, 1 vector
//
// ############################# Public Method ###############################
FloatVector operator< (const FloatVector &x, float y)
{
  FloatVector   result = x;
  result    <= y;                           //Stroustrup p. 302
  return result;
}

