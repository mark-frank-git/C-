/************************************************************************
 *                                                                      *
 * This subclass of object implements a polynomial object.  The         *
 * coefficients of the polynomial are Complex variables.                *
 * The usual operations are defined for this class, +, x, -, /, +=, etc.*
 * This class should be useful for digital signal processing.           *
 *                                                                      *
 * File:ComplexPoly.h                                                   *
 *                                                                      *
 * Note: for p(x) = a[0]*x^n + a[1]*x^(n-1) + . . . + a[n]              *
 *       p(x) has order n, but the length of the array is n+1.          *
 *                                                                      *
 ************************************************************************/

#include <stdio.h>
#include <string.h>
#include "ComplexPoly.h"                                                // Class prototypes
#include "DoublePoly.h"
#if defined(WIN32)
#include <GNU/Complex.h>
#include <PolyRoots/header.h>
#include <AlgebraRoutines/algebra_routines.h>
#else
#include "Complex.h"
#endif

#define ABS_EPS         1.e-20                          // For check if leading coeff == 0.
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
// Input:           newOrder        New order of polynomial
//          
// Output:                          None
//
// Notes:
// ############################# Private Method ###############################
void ComplexPoly::updateMemory (int newOrder)
{
  int       old_size;
  Complex   *temp;                          // Pointer to old memory

  _order     = newOrder;
  old_size  = _size;                     // Save old size
  if(_size > newOrder)
    return;                             // No need to reallocate
  temp      = _coefficients;
  _size      = newOrder + GROW_SIZE;
  _coefficients  = new Complex[_size];
  for(int i=0; i<old_size; i++)         // Copy old coefficients
    _coefficients[i] = temp[i];
  delete [] temp;
  return;
}

// ############################# Private Method ###############################
// maxCoeffAbs -- Finds and returns the abs of the maximum coefficient
// Input:                       None
//          
// Output:                      Norm of maximum coefficient
//
// Notes:
// ############################# Private Method ###############################
double ComplexPoly:: maxCoeffAbs ()
{
  int           i;
  double        max_coeff, abs_coeff;

  max_coeff     = abs(_coefficients[0]);
  for(i=1; i<=_order; i++)
  {
    abs_coeff   = abs(_coefficients[i]);
    if(abs_coeff > max_coeff)
      max_coeff = abs_coeff;
  }
  return        max_coeff;
}

// ############################# Private Method ###############################
// updateOrder -- Check that x[order+1] != 0, if 0, decrement order accordingly
// Input:                           None
//          
// Output:                          None
//
// Notes:
// ############################# Private Method ###############################
void ComplexPoly::updateOrder ()
{
  double        max_coeff_abs;

  max_coeff_abs         = maxCoeffAbs();
  if(max_coeff_abs == 0.)
  {
    _order      = 0;
    return;
  }
  while( abs(_coefficients[0]) < max_coeff_abs * ABS_EPS)
  {
    _order--;
    for(int i=0; i<=_order; i++)
      _coefficients[i]  = _coefficients[i+1];
    if(_order == -1)
      break;
  }
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the ComplexPoly class.
// Input:           order:          order of the polynomial
//                  coeff:          polynomial coefficients
//          
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
ComplexPoly::ComplexPoly(int pOrder, const Complex *coeff)
{
  int i;
// init instance variables:
  _order         = pOrder;
  _size          = _order+1;
  _size          = MAX(1, _size);
  _polyString    = NULL;
  _polyRoots     = NULL;
  _coefficients  = new Complex[_size];
  for(i=0; i<_size; i++)
  {
    if(coeff == NULL)
      _coefficients[i]   = 1;
    else
      _coefficients[i]   = coeff[i];
  }
  updateOrder();                            // Remove leading zeros
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Copy constructor for the ComplexPoly class.
// Input:           poly:           A previously allocated polynomial
//          
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
ComplexPoly::ComplexPoly(const ComplexPoly& poly)
{
  int   i;
  const Complex *input_coeff;

// init instance variables:
  _order         = poly.getOrder();
  _size          = _order+1;
  _size          = MAX(1, _size);
  _coefficients  = new Complex[_size];
  _polyString    = NULL;
  _polyRoots     = NULL;
  input_coeff   = poly.getCoefficients();
  for(i=0; i<_size; i++)
    _coefficients[i] = input_coeff[i];
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Copy constructor for the ComplexPoly class.
// Input:           x:          A Complex number
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
ComplexPoly::ComplexPoly(Complex x)
{
// init instance variables:
 _order                 = 0;
 _size                  = _order+1;
 _coefficients          = new Complex[_size];
 _polyString            = NULL;
 _polyRoots             = NULL;
 _coefficients[0]       = x;
 return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Copy constructor for the ComplexPoly class.
// Input:           x:          A double number
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
ComplexPoly::ComplexPoly(double x)
{
// init instance variables:
  _order                        = 0;
  _size                 = _order+1;
  _coefficients         = new Complex[_size];
  _polyString           = NULL;
  _polyRoots            = NULL;
  _coefficients[0]      = Complex(x, 0.);
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the ComplexPoly class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
ComplexPoly::~ComplexPoly()
{
  delete [] _coefficients;
  delete [] _polyString;
  delete [] _polyRoots;
  return;
}

// ############################# Public Method ###############################
// assign -- Sets the polynomial's coefficients to the input coefficient array
// Input:           aPoly:          array of new coefficients for the polynomial
//                  order:          order of the A polynomial
//          
// Output:                      None
// ############################# Public Method ###############################
void ComplexPoly::assign(int pOrder, const Complex *coeff)
{
  int i;

  updateMemory(pOrder);
  for(i=0; i<(pOrder+1); i++)
  {
    if(coeff == NULL)
      _coefficients[i]   = 1;
    else
      _coefficients[i]   = coeff[i];
  }
  updateOrder();                            // Remove leading zeros
  return;
}


// ############################# Public Method ###############################
// setRoots -- Sets the polynomial's coefficients according to a set of roots
// Input:           roots:      new set of roots.
//                  order:      order of the polynomial
//          
// Output:                      None
//
// Notes:
// 1. It may be more efficient/computationally accurate to try to match up
//    complex roots before multiplying out the polynomial, but this was not
//    done for simplicity.
// ############################# Public Method ###############################
void ComplexPoly::setRoots(int pOrder, const Complex *roots)
{
  int i;
  ComplexPoly   new_root;
  Complex       _coefficients[2];

//
// Set polynomial to first root:
//
  _coefficients[0]      = Complex(1., 0.);
  _coefficients[1]      = Complex(-roots[0].real(), -roots[0].imag());
  assign(1, _coefficients);
//
// Multiply by the rest of the roots
  for(i=1; i<pOrder; i++)
  {
    _coefficients[0]    = Complex(1., 0.);
    _coefficients[1]    = Complex(-roots[i].real(), -roots[i].imag());
    new_root.assign(1, _coefficients);
    *this               *= new_root;
  }
  return;
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
void ComplexPoly::reverse()
{
  int i, n;
  Complex temp;
  
  n = _order;
  for(i=0; i<(_order+1)/2; i++)
  {
    temp                = _coefficients[i];
    _coefficients[i]     = _coefficients[n];
    _coefficients[n--]   = temp;
  }
  updateOrder();
  return;
}

// ############################# Public Method ###############################
// evaluateAtComplexPoint -- Evaluate the polynomial in a at a complex point.
// Input:       point:          The complex point
//          
// Output:                      The complex polynomial value
//
// Notes:
//  1. The instance variables, a[] are modified.
// ############################# Public Method ###############################
Complex ComplexPoly::evaluateAtComplexPoint(Complex &point) const
{
  int     i;
  Complex result, prod, sum; 

/****************************
 * Evaluate polynomial using*
 * Horner's method:         *
 ****************************/
  if(_order>0)
  {
   prod         = _coefficients[0];
   prod         *= point;
   for(i=1; i<=_order; i++)
   {
     sum        = _coefficients[i] + prod;
     prod       = sum*point; 
   } 
   result   = sum;
  }
  else
    result = _coefficients[0];  
  return result;
}

#define REAL_SCALE      1.e5
#define IMAG_SCALE      1.e6
// ############################# Public Method ###############################
// roots -- Find the roots of the polynomial. Return roots in the form:
//             c[0], c[1] .. r[0], r[1] ... c*[0], c*[1] ... where the c's are
//             the complex roots, and the r's are the real roots
// Input:                       None
//          
// Output:                      The complex roots of A
//
// Notes:
//  1. Use Jenkins Traub root finder.
// ############################# Public Method ###############################
Complex *ComplexPoly::roots()
{
#if defined(WIN32)
  int           i, n;
  double        newerr;         // error of actual root
  double        maxerr;         // max error of all the roots
  dcomplex      ns;             // root determined by Muller's method
  dcomplex      *complex_a, *rice_roots, *pred, *rice_ptr;
  int           nred;           // highest exponent of deflated polynomial
  unsigned      char error;     // indicates an error in poly_check
  unsigned      char flag = 1;  // Complex coefficients
  int           red,
                diff;           // number of roots at 0

  if(_polyRoots != NULL)
    delete [] _polyRoots;
  _polyRoots = new Complex[MAX(1, _order)];
//
// Check boundary case:
//
  if(_order < 1)
  {
    _polyRoots[0] = Complex(0., 0.);
    return _polyRoots;
  }
//
// Allocate arrays:
//
  complex_a     = new dcomplex[(_order+1)];
  rice_roots    = new dcomplex[(_order+1)];
  rice_ptr      = rice_roots;
  pred          = new dcomplex[(_order+1)];
  for(i=0; i<=_order; i++)
  {
    complex_a[i].r      = _coefficients[_order-i].real();       // Save the coefficients backwards for Rice routines
    complex_a[i].i      = _coefficients[_order-i].imag();
  }

/*******************************
 * Use Rice' routines to        *
 * find roots.                  *
 *******************************/
  maxerr        = 0.;           //initialize max. error of determined roots
  nred          = _order;       // At the beginning: degree defl. polyn. =
                                // degree of original polyn.
  n             = _order;
                                // check input of the polynomial
  error         = poly_check(complex_a,&nred,&n,rice_roots);
  diff          = (n-nred);     // reduce polynomial, if roots at 0
  rice_roots    += diff;
  n             = nred;         // update n

  if (error)
  {
    printf("error in roots\n");
    delete [] complex_a;
    delete [] rice_roots;
    delete [] pred;
    return _polyRoots;          // error in poly_check(); return error 
  }

                                // polynomial is linear or quadratic
  if (lin_or_quad(complex_a,nred,rice_roots)==0)
  {
    n           += diff;        // remember roots at 0
    maxerr      = DBL_EPSILON; 
  }
  else
  {
    monic(complex_a,&n);        // get monic polynom

    for (i=0;i<=n;i++)
      pred[i]=complex_a[i];     // original polynomial
                                // = deflated polynomial at beginning
                                // of Muller
     do
     { 
                                                // Muller method
          ns                    = muller(pred,nred); 
                                                // Newton method
          rice_roots[nred-1]    = newton(complex_a,n,ns,&newerr,flag);
                                                // stores max. error of all roots
          if (newerr>maxerr)  
               maxerr   = newerr;
                                                // deflate polynomial
          red = poldef(pred,nred,rice_roots,flag);
          pred += red;                          // forget lowest coefficients
          nred -= red;                          // reduce degree of polynomial
     } while (nred>2); 
                                                // last one or two roots
     (void) lin_or_quad(pred,nred,rice_roots);
     if (nred==2)
     {
          rice_roots[1] = newton(complex_a,n,rice_roots[1],&newerr,flag);
          if (newerr>maxerr)
               maxerr   =newerr;
     }
     rice_roots[0]      = newton(complex_a,n,rice_roots[0],&newerr,flag);
     if (newerr>maxerr)
          maxerr=newerr;

     n  += diff;                        // remember roots at 0
  }  
/*******************************
 * Can't re-order roots, since  *
 * they might not appear in     *
 * complex pairs.               *
 *******************************/
  for(i=0; i<_order; i++)
  {
    _polyRoots[i]       = Complex(rice_ptr[i].r, rice_ptr[i].i);
  }
  delete [] complex_a;
  delete [] rice_roots;
  delete [] pred;

#endif
  return _polyRoots;
}

// ############################# Public Method ###############################
// getPolyString()  - Print the polynomial
// Input:                       None
//          
// Output:                      Polynomial string
//
// ############################# Public Method ###############################
const char *ComplexPoly::getPolyString()
{
  int i;
  char temp[64];

  if(_polyString == NULL)
    _polyString  = new char[MAX_STRING_SIZE];
  _polyString[0] = '\0';
  for(i=0; i<_order; i++)
  {
    sprintf(temp, "(%g+%gi)*x^%d + ", real(_coefficients[i]), imag(_coefficients[i]), _order-i);
    strcat(_polyString, temp);
  }
  sprintf(temp, "(%g+%gi)",  real(_coefficients[_order]), imag(_coefficients[_order]));
  strcat(_polyString, temp);
  return _polyString;
}

#ifdef USE_IO_STREAM
// ############################# Public Method ###############################
// printPoly()  - Print the polynomial
// Input:           fout:           ofstream for printing
//          
// Output:                          None
//
// ############################# Public Method ###############################
void ComplexPoly::printPoly(ofstream &fout) const
{
  const char    *output;

  output    = getPolyString();
  fout << output << "\n";
  return;
}
#endif

// ############################# Public Method ###############################
// Operator =  Sets the polynomial equal to the input polynomial
// Input:           y:          input polynomial
//          
// Output:                      the result polynomial
//
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator = (const ComplexPoly& y)
{
  const Complex *input_coeff;

  if( this == &y )                          // Check for x=x
    return *this;
  input_coeff   = y.getCoefficients();
  updateMemory(y.getOrder());
  for(int i=0; i<_order+1; i++)
    _coefficients[i] = input_coeff[i];
  return *this;
}

// ############################# Public Method ###############################
// Operator =  Sets the polynomial equal to the input real polynomial
// Input:           y:          input polynomial
//          
// Output:                      the result polynomial
//
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator = (const DoublePoly& y)
{
  const double *input_coeff;

  input_coeff   = y.getCoefficients();
  updateMemory(y.getOrder());
  for(int i=0; i<_order+1; i++)
    _coefficients[i] = Complex(input_coeff[i], 0.);
  return *this;
}

// ############################# Public Method ###############################
// Operator =  Sets the polynomial equal to the input Complex
// Input:           y:          input Complex
//
// Output:                      the result polynomial
//
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator = (Complex y)
{
  updateMemory(0);
  _coefficients[0]      = y;
  return *this;
}

// ############################# Public Method ###############################
// Operator =  Sets the polynomial equal to the input double
// Input:           y:          input double
//
// Output:                      the result polynomial
//
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator = (double y)
{
  updateMemory(0);
  _coefficients[0]      = Complex(y, 0.);
  return *this;
}

// ############################# Public Method ###############################
// Operator +=  Add the input polynomial to this.
// Input:           y:          input polynomial
//          
// Output:                      the result polynomial
//
// NOTES:
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator += (const ComplexPoly& y)
{
  int           i, j, old_order, y_order;
  const         Complex *input_coeff;
  Complex       *save_coeff;
//
// Get and save y coefficients in case this == y
//
  y_order       = y.getOrder();
  input_coeff   = y.getCoefficients();
  save_coeff    = new Complex[y_order+1];
  for(i=0; i<y_order+1; i++)
    save_coeff[i]       = input_coeff[i];
//
// update order:
//
  old_order     = _order;
  _order         = MAX(_order, y_order);
  updateMemory(_order);
  if(_order > old_order)
  {
    j                   = _order;
    for(i=old_order; i>=0; i--)
    {
      _coefficients[j]  = _coefficients[i] + save_coeff[j];
      j--;
    }
    for(i=0; i<(_order-old_order); i++)
      _coefficients[i]  = save_coeff[i];
  }
  else
  {
    j                   = _order;
    for(i=y_order; i>=0; i--)
    {
      _coefficients[j]  += save_coeff[i];
      j--;
    }
  }
  updateOrder();                                // Check that x[0] != 0
  delete [] save_coeff;
  return *this;
}

// ############################# Public Method ###############################
// Operator -=  Add the input polynomial to this.
// Input:           y:          input polynomial
//          
// Output:                      the result polynomial
//
// NOTES:
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator -= (const ComplexPoly& y)
{
  int           i, j, old_order, y_order;
  const         Complex *input_coeff;
  Complex       *save_coeff;

//
// Get and save y coefficients in case this == y
//
  y_order       = y.getOrder();
  input_coeff   = y.getCoefficients();
  save_coeff    = new Complex[y_order+1];
  for(i=0; i<y_order+1; i++)
    save_coeff[i]       = input_coeff[i];
//
// update order:
//
  old_order     = _order;
  _order        = MAX(_order, y_order);
  updateMemory(_order);
  if(_order > old_order)
  {
    j                   = _order;
    for(i=old_order; i>=0; i--)
    {
      _coefficients[j]  = _coefficients[i]- save_coeff[j];
      j--;
    }
    for(i=0; i<(_order-old_order); i++)
      _coefficients[i]  = -save_coeff[i];
  }
  else
  {
    j                   = _order;
    for(i=y_order; i>=0; i--)
    {
      _coefficients[j]  -= save_coeff[i];
      j--;
    }
  }
  updateOrder();                                // Check that x[0] != 0
  delete [] save_coeff;
  return *this;
}

// ############################# Public Method ###############################
// Operator *=  Mult this by the input polynomial.
// Input:           y:          input polynomial
//          
// Output:                      the result polynomial
//
// NOTES:
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator *= (const ComplexPoly& y)
{
  int       i, j, k, m, n, old_order, y_order;
  const     Complex *input_coeff;
  Complex   *temp, *save_coeff;

//
// Get and save y coefficients in case this == y
//
  y_order       = y.getOrder();
  input_coeff   = y.getCoefficients();
  save_coeff    = new Complex[y_order+1];
  for(i=0; i<y_order+1; i++)
    save_coeff[i]       = input_coeff[i];
//
// Update this order
//
  old_order     = _order;
  _order         += y_order;
  updateMemory(_order);
  temp          = new Complex[old_order+1];
  for(i=0; i<old_order+1; i++)
    temp[i]     = _coefficients[i];                  // Save old coefficients
  for(i=0; i<=_order; i++)
  {
    _coefficients[i]            = 0;
    m                           = MIN(i, old_order);
    for(j=0; j<=m; j++)
    {
      n                         = MIN(i, y_order);
      for(k=0; k<=n; k++)
      {
        if(j+k == i)
          _coefficients[i]      += temp[j]*save_coeff[k];
      }
    }
  }
  delete [] temp;
  delete [] save_coeff;
 
  return *this;
}

// ############################# Public Method ###############################
// Operator /=  Divide this by the input polynomial.
// Input:           y:          input polynomial
//          
// Output:                      the result polynomial
//
// NOTES:
// 1. I am not sure if this works.
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator /= (const ComplexPoly& y)
{
  int           i, y_order, number_shifts;
  static        int old_shifts = -1;
  static        Complex *dividend;
  const         Complex *shifted_coeff;
  Complex       shifted_xn_coeff;
  ComplexPoly   sub_poly;

// First check for boundary case
  y_order       = y.getOrder();
  if(y_order > _order)
  {
    _order           = 0;
    _coefficients[0] = 0;
    return *this;
  }
//
// Allocate the new poly coefficients
//
  number_shifts     = _order - y_order;
  if(old_shifts < number_shifts)                            // Reduce calls to new
  {
    delete [] dividend;
    old_shifts  = number_shifts;
    dividend    = new Complex[number_shifts+1];
  }
  for(i=0; i<number_shifts+1; i++)
  {
    dividend[i] = Complex(0.,0.);
  }
//
// Now, do long division
//
  ComplexPoly   shifted_poly = y << number_shifts;
  shifted_coeff         = shifted_poly.getCoefficients();
  shifted_xn_coeff      = shifted_coeff[0];
  if(norm(shifted_xn_coeff)== 0.)
    return *this;                                       // Error, shouldn't happen
  for(i=0; i<number_shifts+1; i++)
  {
    if(_order < y_order)                                 // Time to exit?
      break;
    if(_order >= shifted_poly.getOrder())                // Is subtract possible
    {
      dividend[i]       = _coefficients[0]/shifted_xn_coeff;
      sub_poly          = shifted_poly * dividend[i];
      *this             -= sub_poly;                    // Yes update this
    }
    else
      dividend[i]       = Complex(0.,0.);
    shifted_poly >>= 1;
  }
//
// Modify this, and return it
//
  _order = number_shifts;
  for(i=0; i<_order+1; i++)
    _coefficients[i] = dividend[i];
    
  return *this;
}

// ############################# Public Method ###############################
// Operator %=  Return remainder of this divided by the input polynomial.
// Input:           y:          input polynomial
//          
// Output:                      the result polynomial
//
// NOTES:
// 1. Mod 2 subtraction = addition.
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator %= (const ComplexPoly& y)
{
  ComplexPoly   div;

// First check for boundary case
  if(_order < y.getOrder())
  {
    return *this;
  }

//
// Modulo equals remainder after division:
//
  div   = *this / y;
  div   *= y;
  *this -= div;
  
//
  updateOrder();        // Check that we don't have initial coefficients equal to zero
    
  return *this;
}

// ############################# Public Method ###############################
// Operator >>=  Shift the polynomial right by n shifts.
// Input:           numberShifts:   Number of places to shift right
//          
// Output:                          the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + . . . + a[n]
//
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator >>= (int numberShifts)
{
//
// Bound the # of shifts, and update the order
//
  numberShifts  = MAX(0, numberShifts);
  numberShifts  = MIN(numberShifts, MAX_SHIFT);
  _order         -= numberShifts;
  _order         = MAX(0, _order);
  return *this;
}

// ############################# Public Method ###############################
// Operator <<=  Shift the polynomial left by n shifts.
// Input:           numberShifts:   Number of places to shift left
//          
// Output:                          the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + . . . + a[n]
//
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator <<= (int numberShifts)
{
  int   i, old_order;

//
// Bound the # of shifts, and update the order
//
  numberShifts  = MAX(0, numberShifts);
  numberShifts  = MIN(numberShifts, MAX_SHIFT);
  old_order     = _order;
  _order         += numberShifts;
  updateMemory(_order);
//
// Do the shifting by filling in zeros:
//
  for(i=old_order+1; i<(_order+1); i++)
    _coefficients[i] = Complex(0.,0.);
  updateOrder();                                    // Remove leading zeros
  return *this;
}

// ############################# Public Method ###############################
// Operator +=  Adds a Complex to the polynomial.
// Input:           y:              Complex to add
//
// Output:                          the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator += (Complex y)
{
  _coefficients[_order]   += y;
  return *this;
}

// ############################# Public Method ###############################
// Operator -=  Subtracts a Complex from the polynomial.
// Input:           y:              Complex to sub
//
// Output:                          the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator -= (Complex y)
{
 _coefficients[_order]   -= y;
 return *this;
}

// ############################# Public Method ###############################
// Operator *=  Mults a Complex by the polynomial.
// Input:           y:              Complex to mult
//
// Output:                          the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator *= (Complex y)
{
  for(int i=0; i<=_order; i++)
  {
    _coefficients[i]    *= y;
  }
  return *this;
}

// ############################# Public Method ###############################
// Operator *=  Divides the polynomial by a Complex.
// Input:           y:              Complex to divide
//
// Output:                          the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator /= (Complex y)
{
  if(norm(y) > 0)
  {
    for(int i=0; i<=_order; i++)
    {
      _coefficients[i]  /= y;
    }
  }
  return *this;
}

// ############################# Public Method ###############################
// Operator +=  Adds a double to the polynomial.
// Input:           y:              double to add
//
// Output:                          the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator += (double y)
{
  _coefficients[_order]   += Complex(y, 0.);
  return *this;
}

// ############################# Public Method ###############################
// Operator -=  Subtracts a double from the polynomial.
// Input:           y:              double to sub
//
// Output:                          the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator -= (double y)
{
  _coefficients[_order]   -= Complex(y, 0.);
  return *this;
}

// ############################# Public Method ###############################
// Operator *=  Mults a double by the polynomial.
// Input:           y:              double to mult
//
// Output:                          the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator *= (double y)
{
 for(int i=0; i<=_order; i++)
 {
   _coefficients[i]     *= Complex(y, 0.);
 }
 return *this;
}

// ############################# Public Method ###############################
// Operator *=  Divides the polynomial by a double.
// Input:           y:              double to divide
//
// Output:                          the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
ComplexPoly& ComplexPoly::operator /= (double y)
{
 if(norm(y) > 0)
 {
   for(int i=0; i<=_order; i++)
   {
     _coefficients[i]   /= Complex(y, 0.);
   }
 }
 return *this;
}


// ############################# Public Method ###############################
// Operator ==  Checks whether the two polynomials are equal
// Input:           x:          first polynomial
//                  y:          second polynomial
//          
// Output:                      YES if two polynomials are equal, else NO
//
// NOTES:
// 1. This is probably not useful for floating point coefficients
// ############################# Public Method ###############################
LOGICAL operator == (const ComplexPoly& x, const ComplexPoly& y)
{
  int   i, n;
  const Complex *x_coeff, *y_coeff;
  LOGICAL  equal;
  
  if( (n = x.getOrder()) != y.getOrder())
    return NO;

  equal     = YES;
  x_coeff   = x.getCoefficients();
  y_coeff   = y.getCoefficients();
  for(i=0; i<(n+1); i++)
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
// Operator !=  Checks whether the two polynomials are not equal
// Input:           x:          first polynomial
//                  y:          second polynomial
//          
// Output:                      YES if two polynomials are not equal, else NO
//
// ############################# Public Method ###############################
LOGICAL operator != (const ComplexPoly& x, const ComplexPoly& y)
{
  return (!(y==x));
}

// ############################# Public Method ###############################
// Operator >  Checks whether the first polynomial is greater than the second
// Input:           a:          first polynomial
//                  b:          second polynomial
//          
// Output:                      YES if first is greater than second
//
// NOTES:
// 1. For example, if a(x) = x^2 + 1, b(x) = x + 1, return YES, similarly,
//    if a(x) = x^2 + x, b(x) = x^2 + 1, return YES, etc.
//
// ############################# Public Method ###############################
LOGICAL operator > (const ComplexPoly& a, const ComplexPoly& b)
{
  int   i, a_order, b_order;
  const Complex *a_coeff, *b_coeff;
  LOGICAL  a_greater;
  
  if( (a_order = a.getOrder()) > (b_order = b.getOrder()) )
    return YES;
  else if( a_order < b_order )
    return NO;

  a_greater = NO;
  a_coeff   = a.getCoefficients();
  b_coeff   = b.getCoefficients();
  for(i=0; i<(a_order+1); i++)
  {
    if( abs(a_coeff[i]) > abs(b_coeff[i]))
    {
      a_greater = YES;
      break;
    }
  }
  return a_greater;
}

// ############################# Public Method ###############################
// Operator >=  Checks whether the first polynomial is greater or equal than the second
// Input:           a:          first polynomial
//                  b:          second polynomial
//          
// Output:                      YES if first is greater than second
//
// NOTES:
// 1. For example, if a(x) = x^2 + 1, b(x) = x + 1, return YES, similarly,
//    if a(x) = x^2 + x, b(x) = x^2 + 1, return YES, etc.
//
// ############################# Public Method ###############################
LOGICAL operator >= (const ComplexPoly& a, const ComplexPoly& b)
{
  return !(a<b);
}

// ############################# Public Method ###############################
// Operator <  Checks whether the first polynomial is less than the second
// Input:           a:          first polynomial
//                  b:          second polynomial
//          
// Output:                      YES if first is less than second
//
// NOTES:
// 1. For example, if a(x) = x^2 + 1, b(x) = x + 1, return NO, similarly,
//    if a(x) = x^2 + x, b(x) = x^2 + 1, return NO, etc.
//
// ############################# Public Method ###############################
LOGICAL operator < (const ComplexPoly& a, const ComplexPoly& b)
{
  int   i, a_order, b_order;
  const Complex *a_coeff, *b_coeff;
  LOGICAL  a_less;
  
  if( (a_order = a.getOrder()) < (b_order = b.getOrder()) )
    return YES;
  else if( a_order > b_order )
    return NO;

  a_less    = NO;
  a_coeff   = a.getCoefficients();
  b_coeff   = b.getCoefficients();
  for(i=0; i<(a_order+1); i++)
  {
    if( abs(a_coeff[i]) < abs(b_coeff[i]) )
    {
      a_less    = YES;
      break;
    }
  }
  return a_less;
}

// ############################# Public Method ###############################
// Operator <=  Checks whether the first polynomial is less or equal than the second
// Input:           a:          first polynomial
//                  b:          second polynomial
//          
// Output:                      YES if first is less or equal than second
//
// NOTES:
// 1. For example, if a(x) = x^2 + 1, b(x) = x + 1, return NO, similarly,
//    if a(x) = x^2 + x, b(x) = x^2 + 1, return NO, etc.
//
// ############################# Public Method ###############################
LOGICAL operator <= (const ComplexPoly& a, const ComplexPoly& b)
{
  return !(a>b);
}

// ############################# Public Method ###############################
// Operator -  Returns the negative of the input polynomial
// Input:           x:          input polynomial
//          
// Output:                      the output polynomial
//
// NOTES
// 1. The additive inverse for modulo-2 addition is the element itself, so
//    we just call the constructor.
//
// ############################# Public Method ###############################
ComplexPoly operator - (const ComplexPoly& x)
{
  return ComplexPoly(x);
}

// ############################# Public Method ###############################
// Operator >>  Returns the right shift of the input polynomial (0 fill)
// Input:           x:          input polynomial
//                  numberShift # of places to shift left
//          
// Output:                      the output polynomial
//
// NOTES:
// 1. The maximum shift is equal to MAX_SHIFT
// ############################# Public Method ###############################
ComplexPoly operator >> (const ComplexPoly& x, int numberShift)
{
  int       i, new_order;
  static    int old_order = 0;
  static    Complex *shift  = NULL;
  const     Complex *input_coeff;
//
// Bound the # of shifts
//
  numberShift   = MAX(0, numberShift);
  numberShift   = MIN(numberShift, MAX_SHIFT);
  new_order     = x.getOrder() - numberShift;
  new_order     = MAX(0, new_order);
//
// Allocate the new poly coefficients and load them
//
  if(old_order < new_order)                         // Reduce calls to new
  {
    delete [] shift;
    old_order   = new_order;
    shift       = new Complex[new_order+1];
  }
//
// Check boundary condition:
//
  if(new_order==0)
    shift[0] = 0.;
  else
  {
    input_coeff = x.getCoefficients();
    for(i=0; i<new_order+1; i++)
      shift[i]  = input_coeff[i];
  }
    
  return ComplexPoly(new_order, shift);
}

// ############################# Public Method ###############################
// Operator <<  Returns the left shift of the input polynomial (0 fill)
// Input:           x:          input polynomial
//                  numberShift # of places to shift left
//          
// Output:                      the output polynomial
//
// NOTES:
// 1. The maximum shift is equal to MAX_SHIFT
// ############################# Public Method ###############################
ComplexPoly operator << (const ComplexPoly& x, int numberShift)
{
  int       i, new_order;
  static    int old_order = 0;
  static    Complex *shift  = NULL;
  const     Complex *input_coeff;
//
// Bound the # of shifts
//
  numberShift   = MAX(0, numberShift);
  numberShift   = MIN(numberShift, MAX_SHIFT);
  new_order     = x.getOrder() + numberShift;
//
// Allocate the new poly coefficients and load them
//
  if(old_order < new_order)                         // Reduce calls to new
  {
    delete [] shift;
    old_order   = new_order;
    shift       = new Complex[new_order+1];
  }
  input_coeff   = x.getCoefficients();
  for(i=0; i<=x.getOrder(); i++)
    shift[i]    = input_coeff[i];
  for(i=(x.getOrder()+1); i<=new_order; i++)
    shift[i]    = Complex(0.,0.);
    
  return ComplexPoly(new_order, shift);
}

// ############################# Public Method ###############################
// Operator +  Returns the sum of the two input polynomials
// Input:           x:          input polynomial
//                  y:          second polynomial
//          
// Output:                      the sum polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator + (const ComplexPoly& x, const ComplexPoly& y)
{
  ComplexPoly   sum = x;
  sum   += y;                           //Stroustrup p. 302
  return sum;
}

// ############################# Public Method ###############################
// Operator -  Returns the difference of the two input polynomials
// Input:           x:          input polynomial
//                  y:          second polynomial
//          
// Output:                      the difference polynomial
//
// NOTES:
// 1. For modulo-2 addition, + equals -.
// ############################# Public Method ###############################
ComplexPoly operator - (const ComplexPoly& x, const ComplexPoly& y)
{
  return (x + y);
}

// ############################# Public Method ###############################
// Operator *  Returns the product of the two input polynomials
// Input:           x:          input polynomial
//                  y:          second polynomial
//          
// Output:                      the product polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator * (const ComplexPoly& x, const ComplexPoly& y)
{
  ComplexPoly   prod = x;
  prod  *= y;                           //Stroustrup p. 302
  return prod;
}

// ############################# Public Method ###############################
// Operator /  Returns the division of the two input polynomials
// Input:           x:          input polynomial
//                  y:          second polynomial
//          
// Output:                      the result polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator / (const ComplexPoly& x, const ComplexPoly& y)
{
  ComplexPoly   result = x;
  result    /= y;                           //Stroustrup p. 302
  return    result;
}

// ############################# Public Method ###############################
// Operator %  Returns the remainder of the division of the two input polynomials
// Input:           x:          input polynomial
//                  y:          second polynomial
//          
// Output:                      the result polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator % (const ComplexPoly& x, const ComplexPoly& y)
{
  ComplexPoly   result = x;
  result    %= y;                           //Stroustrup p. 302
  return    result;
}



// ############################# Public Method ###############################
// Operator +  Returns the sum of a Complex and a polynomial
// Input:           x:          input Complex
//                  y:          input polynomial
//
// Output:                      the sum polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator + (Complex x, const ComplexPoly &y)
{
 ComplexPoly    sum = y;
 sum   += x;                           //Stroustrup p. 302
 return sum;
}

// ############################# Public Method ###############################
// Operator -  Returns the difference of a Complex and a polynomial
// Input:           x:          input Complex
//                  y:          input polynomial
//
// Output:                      the difference polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator - (Complex x, const ComplexPoly &y)
{
 ComplexPoly    result = -y;
 result    += x;                           //Stroustrup p. 302
 return result;
}

// ############################# Public Method ###############################
// Operator *  Returns the product of a Complex and a polynomial
// Input:           x:          input Complex
//                  y:          input polynomial
//
// Output:                      the product polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator * (Complex x, const ComplexPoly &y)
{
 ComplexPoly    product = y;
 product   *= x;                           //Stroustrup p. 302
 return product;
}

// ############################# Public Method ###############################
// Operator /  Returns the division of a Complex by a polynomial
// Input:           x:          input Complex
//                  y:          input polynomial
//
// Output:                      the sum polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator / (Complex x, const ComplexPoly &y)
{
 ComplexPoly    dividend = x;
 dividend  /= y;                           //Stroustrup p. 302
 return dividend;
}

// ############################# Public Method ###############################
// Operator +  Returns the sum of a polynomial and a Complex
// Input:           x:          input polynomial
//                  y:          input Complex
//
// Output:                      the sum polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator + (const ComplexPoly &x, Complex y)
{
 ComplexPoly    sum = x;
 sum   += y;                           //Stroustrup p. 302
 return sum;
}


// ############################# Public Method ###############################
// Operator -  Returns the difference of a polynomial and a Complex
// Input:           x:          input polynomial
//                  y:          input Complex
//
// Output:                      the difference polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator - (const ComplexPoly &x, Complex y)
{
 ComplexPoly    difference = x;
 difference    -= y;                           //Stroustrup p. 302
 return difference;
}

// ############################# Public Method ###############################
// Operator *  Returns the product of a polynomial and a Complex
// Input:           x:          input polynomial
//                  y:          input Complex
//
// Output:                      the product polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator * (const ComplexPoly &x, Complex y)
{
 ComplexPoly    product = x;
 product   *= y;                           //Stroustrup p. 302
 return product;
}

// ############################# Public Method ###############################
// Operator /  Returns the dividend of a polynomial by a Complex
// Input:           x:          input polynomial
//                  y:          input Complex
//
// Output:                      the dividend polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator / (const ComplexPoly &x, Complex y)
{
 ComplexPoly    dividend = x;
 dividend  /= y;                           //Stroustrup p. 302
 return dividend;
}

// ############################# Public Method ###############################
// Operator +  Returns the sum of a double and a polynomial
// Input:           x:          input double
//                  y:          input polynomial
//
// Output:                      the sum polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator + (double x, const ComplexPoly &y)
{
  ComplexPoly    sum = y;
  sum   += x;                           //Stroustrup p. 302
  return sum;
}

// ############################# Public Method ###############################
// Operator -  Returns the difference of a double and a polynomial
// Input:           x:          input double
//                  y:          input polynomial
//
// Output:                      the difference polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator - (double x, const ComplexPoly &y)
{
  ComplexPoly    result = -y;
  result    += x;                           //Stroustrup p. 302
  return result;
}

// ############################# Public Method ###############################
// Operator *  Returns the product of a double and a polynomial
// Input:           x:          input double
//                  y:          input polynomial
//
// Output:                      the product polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator * (double x, const ComplexPoly &y)
{
  ComplexPoly    product = y;
  product   *= x;                           //Stroustrup p. 302
  return product;
}

// ############################# Public Method ###############################
// Operator /  Returns the division of a double by a polynomial
// Input:           x:          input double
//                  y:          input polynomial
//
// Output:                      the sum polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator / (double x, const ComplexPoly &y)
{
  ComplexPoly    dividend = x;
  dividend  /= y;                           //Stroustrup p. 302
  return dividend;
}

// ############################# Public Method ###############################
// Operator +  Returns the sum of a polynomial and a double
// Input:           x:          input polynomial
//                  y:          input double
//
// Output:                      the sum polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator + (const ComplexPoly &x, double y)
{
  ComplexPoly    sum = x;
  sum   += y;                           //Stroustrup p. 302
  return sum;
}

// ############################# Public Method ###############################
// Operator -  Returns the difference of a polynomial and a double
// Input:           x:          input polynomial
//                  y:          input double
//
// Output:                      the difference polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator - (const ComplexPoly &x, double y)
{
  ComplexPoly    difference = x;
  difference    -= y;                           //Stroustrup p. 302
  return difference;
}

// ############################# Public Method ###############################
// Operator *  Returns the product of a polynomial and a double
// Input:           x:          input polynomial
//                  y:          input double
//
// Output:                      the product polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator * (const ComplexPoly &x, double y)
{
  ComplexPoly    product = x;
  product   *= y;                           //Stroustrup p. 302
  return product;
}

// ############################# Public Method ###############################
// Operator /  Returns the dividend of a polynomial by a double
// Input:           x:          input polynomial
//                  y:          input double
//
// Output:                      the dividend polynomial
//
// ############################# Public Method ###############################
ComplexPoly operator / (const ComplexPoly &x, double y)
{
  ComplexPoly    dividend = x;
  dividend  /= y;                           //Stroustrup p. 302
  return dividend;
}

// ############################# Public Method ###############################
// pow()  Returns a polynomial raised to a power
// Input:           x:          input polynomial
//                  n:          input power
//
// Output:                      the polynomial raised to a power
//
// ############################# Public Method ###############################
ComplexPoly pow(const ComplexPoly &x, int n)
{
  ComplexPoly    pow_result = x;
  ComplexPoly   temp_pol   = x;
  while( --n > 0 )
    pow_result  *= temp_pol;
  return pow_result;
}
