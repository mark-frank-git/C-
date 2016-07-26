/************************************************************************
 *                                                                      *
 * This subclass of object implements a polynomial object.  Operations  *
 * on the polynomial are taken over the binary number field, GF(2).     *
 * The usual operations are defined for this class, +, x, -, /, +=, etc.*
 * This class should be useful for analyzing linear block codes.        *
 *                                                                      *
 * File:Mod2Poly.cc                                                     *
 *                                                                      *
 * Note: for p(x) = a[0]*x^n + a[1]*x^(n-1) + . . . + a[n]              *
 *       p(x) has order n, but the size of the array is n+1.            *
 *                                                                      *
 ************************************************************************/

#include <stdio.h>
#include <string.h>
#include "Mod2Poly.h"                                           // Class prototypes


#define MAX(a,b)  ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))

#define GROW_SIZE       10
// ############################# Private Method ###############################
// updateMemory -- Method to perform common memory check/adjustment.
// Input:                       newOrder                New order of polynomial
//                      
// Output:                                      None
//
// Notes:
// ############################# Private Method ###############################
void Mod2Poly::updateMemory (int newOrder)
{
  int   old_size;
  BIT   *temp;                                          // Pointer to old memory

  _order                = newOrder;
  old_size      = _size;                                        // Save old size
  if(_size > newOrder)
    return;                                             // No need to reallocate
  temp          = _coefficients;
  _size         = newOrder + GROW_SIZE;
  _size         = MAX(_size, 1);
  _coefficients = new BIT[_size];
  for(int i=0; i<old_size; i++)                         // Copy old coefficients
    _coefficients[i] = temp[i];
  delete [] temp;
  return;
}

// ############################# Private Method ###############################
// updateOrder -- Check that x[order+1] != 0, if 0, decrement order accordingly
// Input:                                       None
//                      
// Output:                                      None
//
// Notes:
// ############################# Private Method ###############################
void Mod2Poly::updateOrder ()
{
  while(_coefficients[0] == 0)
  {
    if(_order < 1)
      break;
    _order--;                           // reduce order
    for(int i=0; i<=_order; i++)                // shift coefficients
      _coefficients[i]  = _coefficients[i+1];
  }
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the Polynomial class.
// Input:               order:                  order of the polynomial
//                      coeff:                  polynomial coefficients
//                      
// Output:                                      None
//
// Notes:
// ############################# Class Constructor ###############################
Mod2Poly::Mod2Poly(int pOrder, const BIT *coeff)
{
  int i;
// init instance variables:
  _order        = pOrder;
  _size         = _order+1;
  _size         = MAX(1, _size);
  _polyString   = NULL;
  _coefficients = new BIT[_size];
  for(i=0; i<_size; i++)
  {
    if(coeff == NULL)
      _coefficients[i]  = 1;
    else
      _coefficients[i]  = coeff[i];
  }
  updateOrder();                                        // Remove leading zeros
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the Polynomial class.
// Input:               input:                  An input bit
//                      
// Output:                                      None
//
// Notes:
// ############################# Class Constructor ###############################
Mod2Poly::Mod2Poly(BIT input)
{
//
// init instance variables:
//
  _order                = 0;
  _size                 = _order+1;
  _size                 = MAX(1, _size);
  _coefficients         = new BIT[_size];
  _polyString           = NULL;
  _coefficients[0]      = input;
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the Polynomial class.
// Input:               poly:                   A previously allocated polynomial
//                      
// Output:                                      None
//
// Notes:
// ############################# Class Constructor ###############################
Mod2Poly::Mod2Poly(const Mod2Poly& poly)
{
  int   i;
  const BIT *input_coeff;

// init instance variables:
  _order        = poly.getOrder();
  _size         = _order+1;
  _size         = MAX(1, _size);
  _coefficients = new BIT[_size];
  _polyString   = NULL;
  input_coeff   = poly.getCoefficients();
  for(i=0; i<_size; i++)
    _coefficients[i]    = input_coeff[i];
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the Polynomial class.
//
// Input:                                               None
//
// Output:                                              None
//
// Notes:
// ############################# Class Destructor ###############################
Mod2Poly::~Mod2Poly()
{
  delete [] _coefficients;
  delete [] _polyString;
  return;
}

// ############################# Public Method ###############################
// assign -- Sets the polynomial's coefficients to the input coefficient array
// Input:                       aPoly:                  array of coefficients for the A polynomial
//                              order:                  order of the A polynomial
//                      
// Output:                                              None
// ############################# Public Method ###############################
void Mod2Poly::assign(int pOrder, const BIT *coeff)
{
  int i;

  updateMemory(pOrder);
  for(i=0; i<(pOrder+1); i++)
  {
    if(coeff == NULL)
      _coefficients[i]  = 1;
    else
      _coefficients[i]  = coeff[i];
  }
  updateOrder();                                                        // Remove leading zeros
  return;
}

// ############################# Public Method ###############################
// hammingWeight -- Return the number of non-zero coefficients
// Input:                               None
//                      
// Output:                              The hamming weight
// ############################# Public Method ###############################
int Mod2Poly::hammingWeight()
{
  int i, weight;

  weight        = 0;
  for(i=0; i<=_order; i++)
    if(_coefficients[i] == 1)
      weight++;
  return weight;
}

// ############################# Public Method ###############################
// reverse -- Reverse the polynomial.
// Input:                       None
//
// Output:                      None
//
// Notes:
//  1. The instance variables, coefficients[] are modified.
// ############################# Public Method ###############################
void Mod2Poly::reverse()
{
  int i, n;
  BIT temp;

  n = _order;
  for(i=0; i<(_order+1)/2; i++)
  {
    temp                = _coefficients[i];
    _coefficients[i]    = _coefficients[n];
    _coefficients[n--]  = temp;
  }
  return;
}

// ############################# Public Method ###############################
// getPolyString()      - Returns a string for printing the polynomial
// Input:                                       None
//                      
// Output:                                      Pointer to character string
//
// ############################# Public Method ###############################
const char *Mod2Poly::getPolyString()
{
  int i;
  char temp[20];

  if(_polyString == NULL)
    _polyString = new char[MAX_STRING_SIZE];
//
// Check for boundary case:
//
  if(_order == 0)
  {
    sprintf(_polyString,"%d", _coefficients[0]);
    return _polyString;
  }

  if(_coefficients[0] == 1)
    sprintf(_polyString, "x^%d", _order);
  else
    strcpy(_polyString, "0");
  for(i=1; i<_order; i++)
  {
    if(_coefficients[i] != 0)
      sprintf(temp, " + x^%d", _order-i);
    else
      sprintf(temp, " +    ");
    strcat(_polyString, temp);
  }
  if(_coefficients[_order] == 1)
    sprintf(temp," + 1");
  else
    sprintf(temp," + 0");
  strcat(_polyString, temp);
  return _polyString;
}

// ############################# Public Method ###############################
// Operator =  Sets the polynomial equal to the input polynomial
// Input:                       y:                      input polynomial
//                      
// Output:                                              the result polynomial
//
// ############################# Public Method ###############################
Mod2Poly& Mod2Poly::operator = (const Mod2Poly& y)
{
  const BIT *input_coeff;

  if( this == &y )                                                      // Check for x=x
    return *this;
  input_coeff   = y.getCoefficients();
  updateMemory(y.getOrder());
  for(int i=0; i<_order+1; i++)
    _coefficients[i] = input_coeff[i];
  return *this;
}

// ############################# Public Method ###############################
// Operator =  Sets the polynomial equal to the input BIT
// Input:           y:          input BIT
//
// Output:                      the result polynomial
//
// ############################# Public Method ###############################
Mod2Poly& Mod2Poly::operator = (BIT y)
{

  updateMemory(0);
  _coefficients[0]   = y;
  return *this;
}

// ############################# Public Method ###############################
// Operator +=  Add the input polynomial to this.
// Input:                       y:                      input polynomial
//                      
// Output:                                              the result polynomial
//
// NOTES:
// 1. Mod 2 addition uses exclusive or operator (^)
// ############################# Public Method ###############################
Mod2Poly& Mod2Poly::operator += (const Mod2Poly& y)
{
  int   i, j, y_order, old_order;
  const BIT *input_coeff;

  old_order     = _order;
  y_order       = y.getOrder();
  _order        = MAX(_order, y_order);
  updateMemory(_order);
  input_coeff   = y.getCoefficients();
  if(_order > old_order)
  {
    j                   = _order;
    for(i=old_order; i>=0; i--)
    {
      _coefficients[j]  = _coefficients[i] ^ input_coeff[j];
      j--;
    }
    for(i=0; i<(_order-old_order); i++)
      _coefficients[i]  = input_coeff[i];
  }
  else
  {
    j                   = _order;
    for(i=y_order; i>=0; i--)
    {
      _coefficients[j]  ^= input_coeff[i];
      j--;
    }
  }
  updateOrder();                                // Check that x[0] != 0
  return *this;
}

// ############################# Public Method ###############################
// Operator -=  Sub the input polynomial from this.
// Input:                       y:                      input polynomial
//                      
// Output:                                              the result polynomial
//
// NOTES:
// 1. For modulo-2 addition, + equals -.
// ############################# Public Method ###############################
Mod2Poly& Mod2Poly::operator -= (const Mod2Poly& y)
{
  *this += y;
  return *this;
}

// ############################# Public Method ###############################
// Operator *=  Mult this by the input polynomial.
// Input:                       y:                      input polynomial
//                      
// Output:                                              the result polynomial
//
// NOTES:
// 1. Mod 2 multiplication uses and operator (&)
// 2. Mod 2 addition uses exclusive or operator (^)
// ############################# Public Method ###############################
Mod2Poly& Mod2Poly::operator *= (const Mod2Poly& y)
{
  int   i, j, k, m, n, old_order, y_order;
  int   output_pwr, temp_pwr, y_pwr;
  const BIT *y_coeff;
  BIT   *temp;

  y_order               = y.getOrder();
  old_order             = _order;
  _order                += y_order;
  updateMemory(_order);
  y_coeff               = y.getCoefficients();
  temp                  = new BIT[old_order+1];
  for(i=0; i<old_order+1; i++)
    temp[i]             = _coefficients[i];                                     // Save old coefficients
  for(i=0; i<=_order; i++)
  {
    output_pwr          = _order-i;
    _coefficients[i] = 0;
    m                   = MIN(i, old_order);
    for(j=0; j<=m; j++)
    {
      temp_pwr          = old_order-j;
      n                 = MIN(i, y_order);
      for(k=0; k<=n; k++)
      {
        y_pwr           = y_order-k;
        if(output_pwr == (y_pwr+temp_pwr) )
          _coefficients[i] ^= temp[j]&y_coeff[k];
      }
    }
  }
  delete [] temp;
 
  return *this;
}

// ############################# Public Method ###############################
// Operator /=  Divide this by the input polynomial.
// Input:                       y:                      input polynomial
//                      
// Output:                                              the result polynomial
//
// NOTES:
// 1. Mod 2 subtraction = addition.
// ############################# Public Method ###############################
Mod2Poly& Mod2Poly::operator /= (const Mod2Poly& y)
{
  int           i, j, y_order, number_shifts;
  static        int old_shifts = -1;
  static        BIT *dividend;

// First check for boundary case
  y_order               = y.getOrder();
  if(y_order > _order)
  {
    _order              = 0;
    _coefficients[0]    = 0;
    return *this;
  }

//
// Allocate the new poly coefficients
//
  number_shifts = _order - y_order;
  if(old_shifts < number_shifts)                                        // Reduce calls to new
  {
    delete [] dividend;
    old_shifts  = number_shifts;
    dividend    = new BIT[number_shifts+1];
  }
  for(i=0; i<number_shifts+1; i++)
    dividend[i] = 0;
//
// Now, do long division
//
  Mod2Poly      shifted_poly = y << number_shifts;
  for(i=0; i<number_shifts+1; i++)
  {
    if(_order < y_order)                                        // Time to exit?
      break;
    if(_order >= shifted_poly.getOrder())                       // Is add possible
    {
      *this             += shifted_poly;                        // Yes update this by subtracting
      dividend[i]       = 1;
    }
    else
      dividend[i]       = 0;
    shifted_poly >>= 1;
  }
//
// Modify this, and return it
//
  _order        = number_shifts;
  for(j=0; j<_order+1; j++)
    _coefficients[j] = dividend[j];
    
  return *this;
}

// ############################# Public Method ###############################
// Operator %=  Return remainder of this divided by the input polynomial.
// Input:                       y:                      input polynomial
//                      
// Output:                                              the result polynomial
//
// NOTES:
// 1. Mod 2 subtraction = addition.
// ############################# Public Method ###############################
Mod2Poly& Mod2Poly::operator %= (const Mod2Poly& y)
{
  int           y_order;
  Mod2Poly      dividend;

// First check for boundary case
  y_order       = y.getOrder();
  if(_order < y_order)
  {
    return *this;
  }

//
// Remainder = *this + (y * (*this/y));
//
  dividend      = *this / y;
  *this         += y * dividend;
//
  updateOrder();        // Check that we don't have initial coefficients equal to zero
    
  return *this;
}

// ############################# Public Method ###############################
// Operator ^= Bit wise exclusive or the input polynomial to this.
// Input:                       y:                      input polynomial
//                      
// Output:                                              the result polynomial
//
// NOTES:
// 1. Mod 2 addition uses exclusive or operator (^), this operator is a dupli-
//    cation of += for convenience.
// ############################# Public Method ###############################
Mod2Poly& Mod2Poly::operator ^= (const Mod2Poly& y)
{
  *this += y;
  return *this;
}

// ############################# Public Method ###############################
// Operator &= Bit wise and the input polynomial to this.
// Input:                       y:                      input polynomial
//                      
// Output:                                              the result polynomial
//
// NOTES:
// ############################# Public Method ###############################
Mod2Poly& Mod2Poly::operator &= (const Mod2Poly& y)
{
  int   i, j, k, y_order, old_order;
  const BIT *input_coeff;

  old_order             = _order;
  y_order               = y.getOrder();
  _order                = MIN(_order, y_order);
  updateMemory(_order);
  input_coeff   = y.getCoefficients();
//
// & the input with this
//
  j                     = y_order;
  k                     = old_order;
  for(i=_order; i>=0; i--)
  {
      _coefficients[i]  = _coefficients[k] & input_coeff[j];
      j--;
      k--;
  }
  updateOrder();                                                // Check that x[0] != 0
  return *this;
}

// ############################# Public Method ###############################
// Operator >>=  Shift the polynomial right by n shifts.
// Input:                       numberShifts:   Number of places to shift right
//                      
// Output:                                      the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + . . . + a[n]
//
// ############################# Public Method ###############################
Mod2Poly& Mod2Poly::operator >>= (int numberShifts)
{
//
// Bound the # of shifts, and update the order
//
  numberShifts  = MAX(0, numberShifts);
  numberShifts  = MIN(numberShifts, MAX_SHIFT);
  _order        -= numberShifts;
  _order        = MAX(0, _order);
  return        *this;
}

// ############################# Public Method ###############################
// Operator <<=  Shift the polynomial left by n shifts.
// Input:                       numberShifts:   Number of places to shift left
//                      
// Output:                                      the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + . . . + a[n]
//
// ############################# Public Method ###############################
Mod2Poly& Mod2Poly::operator <<= (int numberShifts)
{
  int   i, old_order;

//
// Bound the # of shifts, and update the order
//
  numberShifts  = MAX(0, numberShifts);
  numberShifts  = MIN(numberShifts, MAX_SHIFT);
  old_order     = _order;
  _order        += numberShifts;
  updateMemory(_order);
//
// Do the shifting by filling zeros
//
  for(i=old_order+1; i<(_order+1); i++)
    _coefficients[i] = 0;
  updateOrder();                                                                        // Remove leading zeros
  return *this;
}

#ifdef USE_IO_STREAM
// ############################# Public Method ###############################
// <<   - Print the polynomial
// Input:       s:                              ofstream for printing
//              x:                              the polynomial to print
//                      
// Output:                                      None
//
// ############################# Public Method ###############################
ostream& operator << (ostream& s, const Mod2Poly& x)
{
  return s << x.getPolyString() << "\n";
}
#endif

// ############################# Public Method ###############################
// Operator ==  Checks whether the two polynomials are equal
// Input:                       x:                      first polynomial
//                              y:                      second polynomial
//                      
// Output:                      YES if two polynomials are equal, else NO
//
// ############################# Public Method ###############################
LOGICAL operator == (const Mod2Poly& x, const Mod2Poly& y)
{
  int   i, n;
  const BIT *x_coeff, *y_coeff;
  LOGICAL       equal;
  
  if( (n = x.getOrder()) != y.getOrder())
    return NO;

  equal         = YES;
  x_coeff       = x.getCoefficients();
  y_coeff       = y.getCoefficients();
  for(i=0; i<(n+1); i++)
  {
    if(x_coeff[i] != y_coeff[i])
    {
      equal     = NO;
      break;
    }
  }
  return equal;
}

// ############################# Public Method ###############################
// Operator !=  Checks whether the two polynomials are not equal
// Input:                       x:                      first polynomial
//                              y:                      second polynomial
//                      
// Output:                      YES if two polynomials are not equal, else NO
//
// ############################# Public Method ###############################
LOGICAL operator != (const Mod2Poly& x, const Mod2Poly& y)
{
  return (!(y==x));
}

// ############################# Public Method ###############################
// Operator >  Checks whether the first polynomial is greater than the second
// Input:                       a:                      first polynomial
//                              b:                      second polynomial
//                      
// Output:                      YES if first is greater than second
//
// NOTES:
// 1. For example, if a(x) = x^2 + 1, b(x) = x + 1, return YES, similarly,
//    if a(x) = x^2 + x, b(x) = x^2 + 1, return YES, etc.
//
// ############################# Public Method ###############################
LOGICAL operator > (const Mod2Poly& a, const Mod2Poly& b)
{
  int   i, a_order, b_order;
  const BIT *a_coeff, *b_coeff;
  LOGICAL       a_greater;
  
  if( (a_order = a.getOrder()) > (b_order = b.getOrder()) )
    return YES;
  else if( a_order < b_order )
    return NO;

  a_greater     = NO;
  a_coeff       = a.getCoefficients();
  b_coeff       = b.getCoefficients();
  for(i=0; i<(a_order+1); i++)
  {
    if(a_coeff[i] > b_coeff[i])
    {
      a_greater = YES;
      break;
    }
  }
  return a_greater;
}

// ############################# Public Method ###############################
// Operator >=  Checks whether the first polynomial is greater or equal than the second
// Input:                       a:                      first polynomial
//                              b:                      second polynomial
//                      
// Output:                      YES if first is greater than second
//
// NOTES:
// 1. For example, if a(x) = x^2 + 1, b(x) = x + 1, return YES, similarly,
//    if a(x) = x^2 + x, b(x) = x^2 + 1, return YES, etc.
//
// ############################# Public Method ###############################
LOGICAL operator >= (const Mod2Poly& a, const Mod2Poly& b)
{
  return !(a<b);
}

// ############################# Public Method ###############################
// Operator <  Checks whether the first polynomial is less than the second
// Input:                       a:                      first polynomial
//                              b:                      second polynomial
//                      
// Output:                                              YES if first is less than second
//
// NOTES:
// 1. For example, if a(x) = x^2 + 1, b(x) = x + 1, return NO, similarly,
//    if a(x) = x^2 + x, b(x) = x^2 + 1, return NO, etc.
//
// ############################# Public Method ###############################
LOGICAL operator < (const Mod2Poly& a, const Mod2Poly& b)
{
  int   i, a_order, b_order;
  const BIT *a_coeff, *b_coeff;
  LOGICAL       a_less;
  
  if( (a_order = a.getOrder()) < (b_order = b.getOrder()) )
    return YES;
  else if( a_order > b_order )
    return NO;

  a_less        = NO;
  a_coeff       = a.getCoefficients();
  b_coeff       = b.getCoefficients();
  for(i=0; i<(a_order+1); i++)
  {
    if(a_coeff[i] < b_coeff[i])
    {
      a_less    = YES;
      break;
    }
  }
  return a_less;
}

// ############################# Public Method ###############################
// Operator <=  Checks whether the first polynomial is less or equal than the second
// Input:                       a:                      first polynomial
//                              b:                      second polynomial
//                      
// Output:                                              YES if first is less or equal than second
//
// NOTES:
// 1. For example, if a(x) = x^2 + 1, b(x) = x + 1, return NO, similarly,
//    if a(x) = x^2 + x, b(x) = x^2 + 1, return NO, etc.
//
// ############################# Public Method ###############################
LOGICAL operator <= (const Mod2Poly& a, const Mod2Poly& b)
{
  return !(a>b);
}

// ############################# Public Method ###############################
// hammingDistance: Returns the hamming distance between two input polynomials
// Input:                       a:                      first polynomial
//                              b:                      second polynomial
//                      
// Output:                                              Hamming distance between two code words
//
//
// ############################# Public Method ###############################
int hammingDistance(const Mod2Poly& a, const Mod2Poly& b)
{
  Mod2Poly      c = a + b;
  return c.hammingWeight();
}

// ############################# Public Method ###############################
// Operator -  Returns the negative of the input polynomial
// Input:                       x:                      input polynomial
//                      
// Output:                                              the output polynomial
//
// NOTES
// 1. The additive inverse for modulo-2 addition is the element itself, so
//    we just call the constructor.
//
// ############################# Public Method ###############################
Mod2Poly operator - (const Mod2Poly& x)
{
  return Mod2Poly(x);
}

// ############################# Public Method ###############################
// Operator <<  Returns the left shift of the input polynomial (0 fill)
// Input:                       x:                      input polynomial
//                                      numberShift     # of places to shift left
//                      
// Output:                                              the output polynomial
//
// NOTES:
// 1. The maximum shift is equal to MAX_SHIFT
// ############################# Public Method ###############################
Mod2Poly operator << (const Mod2Poly& x, int numberShift)
{
  Mod2Poly      shift = x;
  shift <<= numberShift;                                //Stroustrup p. 302
  return shift;
}

// ############################# Public Method ###############################
// Operator >>  Returns the right shift of the input polynomial (0 fill)
// Input:                       x:                      input polynomial
//                              numberShift     # of places to shift left
//                      
// Output:                                              the output polynomial
//
// NOTES:
// 1. The maximum shift is equal to MAX_SHIFT
// ############################# Public Method ###############################
Mod2Poly operator >> (const Mod2Poly& x, int numberShift)
{
  Mod2Poly      shift = x;
  shift >>= numberShift;                                        //Stroustrup p. 302
  return shift;
}

// ############################# Public Method ###############################
// Operator +  Returns the sum of the two input polynomials
// Input:                       x:                      input polynomial
//                                      y:                      second polynomial
//                      
// Output:                                              the sum polynomial
//
// ############################# Public Method ###############################
Mod2Poly operator + (const Mod2Poly& x, const Mod2Poly& y)
{
  Mod2Poly      sum = x;
  sum   += y;                                                   //Stroustrup p. 302
  return sum;
}

// ############################# Public Method ###############################
// Operator -  Returns the difference of the two input polynomials
// Input:                       x:                      input polynomial
//                              y:                      second polynomial
//                      
// Output:                                              the difference polynomial
//
// NOTES:
// 1. For modulo-2 addition, + equals -.
// ############################# Public Method ###############################
Mod2Poly operator - (const Mod2Poly& x, const Mod2Poly& y)
{
  return (x + y);
}

// ############################# Public Method ###############################
// Operator *  Returns the product of the two input polynomials
// Input:                       x:                      input polynomial
//                              y:                      second polynomial
//                      
// Output:                                              the product polynomial
//
// ############################# Public Method ###############################
Mod2Poly operator * (const Mod2Poly& x, const Mod2Poly& y)
{
  Mod2Poly      prod = x;
  prod  *= y;                                           //Stroustrup p. 302
  return prod;
}

// ############################# Public Method ###############################
// Operator /  Returns the division of the two input polynomials
// Input:                       x:                      input polynomial
//                              y:                      second polynomial
//                      
// Output:                                              the result polynomial
//
// ############################# Public Method ###############################
Mod2Poly operator / (const Mod2Poly& x, const Mod2Poly& y)
{
  Mod2Poly      result = x;
  result        /= y;                                           //Stroustrup p. 302
  return        result;
}

// ############################# Public Method ###############################
// Operator %  Returns the remainder of the division of the two input polynomials
// Input:                       x:                      input polynomial
//                              y:                      second polynomial
//                      
// Output:                                              the result polynomial
//
// ############################# Public Method ###############################
Mod2Poly operator % (const Mod2Poly& x, const Mod2Poly& y)
{
  Mod2Poly      result = x;
  result        %= y;                                           //Stroustrup p. 302
  return        result;
}

// ############################# Public Method ###############################
// Operator ^  Returns the exclusive or of the two input polynomials
// Input:                       x:                      input polynomial
//                              y:                      second polynomial
//                      
// Output:                                              the result polynomial
//
// NOTES:
// 1. Exclusive or is the same as +, this is added as a convenience.
// ############################# Public Method ###############################
Mod2Poly operator ^ (const Mod2Poly& x, const Mod2Poly& y)
{
  Mod2Poly      sum = x;
  sum   ^= y;                                                   //Stroustrup p. 302
  return sum;
}

// ############################# Public Method ###############################
// Operator &  Returns the bitwise 'and' of the two input polynomials
// Input:                       x:                      input polynomial
//                              y:                      second polynomial
//                      
// Output:                                              the result polynomial
//
// NOTES:
// ############################# Public Method ###############################
Mod2Poly operator & (const Mod2Poly& x, const Mod2Poly& y)
{
  Mod2Poly      sum = x;
  sum   &= y;                                                   //Stroustrup p. 302
  return sum;
}
