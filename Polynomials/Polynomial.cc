/************************************************************************
 *                                                                      *
 * This subclass of object implements a polynomial object.  Three       *
 * polynomials are contained in this object: a, b, c.  Most of the      *
 * operations take place on the a polynomial, with the b and c poly-    *
 * nomials used for multiplying, etc.                                   *
 *                                                                      *
 * File:Polynomial.h                                                    *
 *                                                                      *
 * Note: for p(x) = a[0]*x^n + a[1]*x^(n-1) + . . . + a[n]              *
 *       p(x) has order n, but the length of the array is n+1.          *
 *                                                                      *
 ************************************************************************/

#include <stdio.h>
#include <string.h>
#include "Polynomial.h"                                             // Class prototypes
#if defined(WIN32)
#include <GNU/Complex.h>
#include <AlgebraRoutines/algebra_routines.h>
#include <PolyRoots/header.h>
#else
#include "Complex.h"
#endif

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
void Polynomial::updateMemory (int newOrder)
{
  int       old_size;
  double    *temp;                          // Pointer to old memory

  order     = newOrder;
  old_size  = size;                     // Save old size
  if(size > newOrder)
    return;                             // No need to reallocate
  temp      = coefficients;
  size      = newOrder + GROW_SIZE;
  coefficients  = new double[size];
  for(int i=0; i<old_size; i++)         // Copy old coefficients
    coefficients[i] = temp[i];
  delete [] temp;
  return;
}

// ############################# Private Method ###############################
// updateOrder -- Check that x[0] != 0, if 0, decrement order accordingly
// Input:                           None
//          
// Output:                          None
//
// Notes:
// ############################# Private Method ###############################
void Polynomial::updateOrder ()
{
  while(coefficients[0] == 0.)
  {
    order--;
    for(int i=0; i<=order; i++)
       coefficients[i]  = coefficients[i+1];
    if(order == -1)
      break;
  }
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the Polynomial class.
// Input:           order:          order of the polynomial
//                  coeff:          polynomial coefficients
//          
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
Polynomial::Polynomial(int pOrder, const double *coeff)
{
  int i;
// init instance variables:
  order         = pOrder;
  size          = order+1;
  size          = MAX(1, size);
  polyString    = NULL;
  polyRoots     = NULL;
  coefficients  = new double[size];
  for(i=0; i<size; i++)
  {
    if(coeff == NULL)
      coefficients[i]   = 1;
    else
      coefficients[i]   = coeff[i];
  }
  updateOrder();                            // Remove leading zeros
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Copy constructor for the Polynomial class.
// Input:           poly:           A previously allocated polynomial
//          
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
Polynomial::Polynomial(const Polynomial& poly)
{
  int   i;
  const double *input_coeff;

// init instance variables:
  order         = poly.getOrder();
  size          = order+1;
  size          = MAX(1, size);
  coefficients  = new double[size];
  polyString    = NULL;
  polyRoots     = NULL;
  input_coeff   = poly.getCoefficients();
  for(i=0; i<size; i++)
    coefficients[i] = input_coeff[i];
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Copy constructor for the Polynomial class, for input double
// Input:           x:              double input
//          
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
Polynomial::Polynomial(double x)
{

// init instance variables:
  order             = 0;
  size              = order+1;
  size              = MAX(1, size);
  coefficients      = new double[size];
  polyString        = NULL;
  polyRoots         = NULL;
  coefficients[0]   = x;
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the Polynomial class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
Polynomial::~Polynomial()
{
  delete [] coefficients;
  delete [] polyString;
  delete [] polyRoots;
  return;
}

// ############################# Public Method ###############################
// assign -- Sets the polynomial's coefficients to the input coefficient array
// Input:           pOrder:         order of the input set of coefficients
//                  coeff:          array of coefficients for the A polynomial
//          
// Output:                      None
// ############################# Public Method ###############################
void Polynomial::assign(int pOrder, const double *coeff)
{
  int i;

  updateMemory(pOrder);
  for(i=0; i<(pOrder+1); i++)
  {
    if(coeff == NULL)
      coefficients[i]   = 1;
    else
      coefficients[i]   = coeff[i];
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
void Polynomial::reverse()
{
  int i, n;
  double temp;
  
  n = order;
  for(i=0; i<(order+1)/2; i++)
  {
    temp                = coefficients[i];
    coefficients[i]     = coefficients[n];
    coefficients[n--]   = temp;
  }
  return;
}

// ############################# Public Method ###############################
// evaluateAtRealPoint -- Evaluate the polynomial in a at a real point.
// Input:       point:          The real point
//          
// Output:                      The real polynomial value
//
// Notes:
// ############################# Public Method ###############################
double Polynomial::evaluateAtRealPoint(double point)
{
  int           i;
  double        result, prod, sum; 

/****************************
 * Evaluate polynomial using*
 * Horner's method:         *
 ****************************/
  if(order>0)
  {
   prod         = coefficients[0];
   prod         *= point;
   sum          = 0.;
   for(i=1; i<=order; i++)
   {
     sum        = coefficients[i];
     sum        += prod;
     prod       = sum*point; 
   } 
   result   = sum;
  }
  else
    result = coefficients[0];  
  return result;
}

// ############################# Public Method ###############################
// evaluateAtComplexPoint -- Evaluate the polynomial in a at a complex point.
// Input:       point:          The complex point
//          
// Output:                      The complex polynomial value
//
// Notes:
// ############################# Public Method ###############################
Complex Polynomial::evaluateAtComplexPoint(Complex &point) const
{
  int     i;
  Complex result, cprod, ctemp, csum; 

/****************************
 * Evaluate polynomial using*
 * Horner's method:         *
 ****************************/
  if(order>0)
  {
   ctemp    = Complex(coefficients[0],0.);
   cprod    = ctemp*point;  
   for(i=1; i<=order; i++)
   {
     ctemp  = Complex(coefficients[i],0.);
     csum   = ctemp+cprod;
     cprod  = csum*point; 
   } 
   result   = csum;
  }
  else
    result = Complex(coefficients[0], 0.);  
  return result;
}

#define REAL_SCALE      1.e5
#define IMAG_SCALE      1.e6
// ############################# Public Method ###############################
// roots -- Find the roots of the polynomial in A. Return roots in the form:
//             c[0], c[1] .. r[0], r[1] ... c*[0], c*[1] ... where the c's are
//             the complex roots, and the r's are the real roots
// Input:                       None
//          
// Output:                      The complex roots of A
//
// Notes:
//  1. Use Rice University's root finder.
//  2. Rice's code defines its own dcomplex structure.
// ############################# Public Method ###############################
Complex *Polynomial::roots()
{
#if defined(WIN32)
  int           i, j, n;
  int           number_complex, number_real;
  LOGICAL       new_complex_root;
  Complex       *complex_roots, *real_roots;
  double        newerr;         // error of actual root
  double        maxerr;         // max error of all the roots
  dcomplex      ns;             // root determined by Muller's method
  dcomplex      *complex_a, *rice_roots, *pred, *rice_ptr;
  int           nred;           // highest exponent of deflated polynomial
  unsigned      char error;     // indicates an error in poly_check
  unsigned      char flag = 0;  // Real coefficients
  int           red,
                diff;           // number of roots at 0

  if(polyRoots != NULL)
    delete [] polyRoots;
  polyRoots = new Complex[MAX(1, order)];
//
// Check boundary case:
//
  if(order < 1)
  {
    polyRoots[0] = Complex(0., 0.);
    return polyRoots;
  }
//
// Allocate arrays:
//
  complex_a     = new dcomplex[(order+1)];
  rice_roots    = new dcomplex[(order+1)];
  rice_ptr      = rice_roots;
  pred          = new dcomplex[(order+1)];
  for(i=0; i<=order; i++)
  {
    complex_a[i].r      = coefficients[order-i];        // Save the coefficients backwards for Rice routines
    complex_a[i].i      = 0.;
  }

/*******************************
 * Use Rice' routines to        *
 * find roots.                  *
 *******************************/
  maxerr        = 0.;           //initialize max. error of determined roots
  nred          = order;        // At the beginning: degree defl. polyn. =
                                // degree of original polyn.
  n             = order;
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
    return polyRoots;           // error in poly_check(); return error 
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
 * Now, re-order roots, to      *
 * match up complex pairs.      *
 *******************************/
  for(i=0; i<order; i++)
  {
    polyRoots[i]        = Complex(rice_ptr[i].r, rice_ptr[i].i);
  }
  real_roots            = new Complex[order];
  complex_roots         = new Complex[order];
  number_complex        = number_real = 0;
  for(i=0; i<order; i++)
  {
    if( ABS(real(polyRoots[i])) > REAL_SCALE*ABS(imag(polyRoots[i])) )
    {
      real_roots[number_real] = polyRoots[i];
      number_real++;
    }
    else
    {
      new_complex_root = YES;
      for(j=0; j<number_complex; j++)
      {
        if( ABS((polyRoots[i].imag() + complex_roots[j].imag())) <=
                            (ABS((polyRoots[i].imag())))/IMAG_SCALE)
        {
          new_complex_root = NO;
          break;
        }
      }
      if( new_complex_root )
      {
        complex_roots[number_complex] = polyRoots[i];
        number_complex++;
      }
    }
  }
  number_complex = order - number_real;                 // error check
  number_complex /= 2;
  for(i=0; i<number_complex; i++)
  {
    polyRoots[i] = complex_roots[i];
  }
  for(j=0; j<number_real; j++)
  {
    polyRoots[i] = real_roots[j];
    i++;
  }
  for(j=0; j<number_complex; j++)                           // conjugate roots
  {
    polyRoots[i] = conj(complex_roots[j]);
    i++;
  }
  delete [] complex_a;
  delete [] rice_roots;
  delete [] pred;
  delete [] real_roots;
  delete [] complex_roots;
#endif
  return polyRoots;
}

#define MAX_COEFF_SIZE  64
// ############################# Public Method ###############################
// getPolyString()  - Get a string for printing out the polynomial
// Input:                       None
//          
// Output:                      Polynomial string
//
// ############################# Public Method ###############################
const char *Polynomial::getPolyString()
{
  int   i, n;
  char  temp[MAX_COEFF_SIZE];

  if(polyString == NULL)
    polyString  = new char[MAX_STRING_SIZE];
  polyString[0] = '\0';
  n             = MAX_STRING_SIZE - 1;
  for(i=0; i<order; i++)
  {
    sprintf(temp, "%g*x^%d + ", coefficients[i], order-i);
    strncat(polyString, temp, n);
    n           -= strlen(temp);
  }
  sprintf(temp, "%g", coefficients[order]);
  strncat(polyString, temp, n);
  return polyString;
}


// ############################# Public Method ###############################
// Operator =  Sets the polynomial equal to the input polynomial
// Input:           y:          input polynomial
//          
// Output:                      the result polynomial
//
// ############################# Public Method ###############################
Polynomial& Polynomial::operator = (const Polynomial& y)
{
  const double *input_coeff;

  if( this == &y )                          // Check for x=x
    return *this;
  input_coeff   = y.getCoefficients();
  updateMemory(y.getOrder());
  for(int i=0; i<order+1; i++)
    coefficients[i] = input_coeff[i];
  return *this;
}

// ############################# Public Method ###############################
// Operator =  Sets the polynomial equal to the input double
// Input:           y:          input double
//          
// Output:                      the result polynomial
//
// ############################# Public Method ###############################
Polynomial& Polynomial::operator = (double y)
{

  updateMemory(0);
  coefficients[0]   = y;
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
Polynomial& Polynomial::operator += (const Polynomial& y)
{
  int           i, j, y_order, old_order;
  const         double *input_coeff;
  double    *save_coeff;
//
// Get and save y coefficients in case this == y
//
  y_order       = y.getOrder();
  input_coeff   = y.getCoefficients();
  save_coeff    = new double[y_order+1];
  for(i=0; i<y_order+1; i++)
    save_coeff[i]       = input_coeff[i];
//
// update order:
//
  old_order     = order;
  order         = MAX(order, y_order);
  updateMemory(order);
  if(order > old_order)
  {
    j           = order;
    for(i=old_order; i>=0; i--)
    {
      coefficients[j]   = coefficients[i] + save_coeff[j];
      j--;
    }
    for(i=0; i<(order-old_order); i++)
      coefficients[i]   = save_coeff[i];
  }
  else
  {
    j           = order;
    for(i=y_order; i>=0; i--)
    {
      coefficients[j]   += save_coeff[i];
      j--;
    }
  }
  updateOrder();                                // Check that x[0] != 0
  delete [] save_coeff;
  return *this;
}

// ############################# Public Method ###############################
// Operator -=  Sub the input polynomial from this.
// Input:           y:          input polynomial
//          
// Output:                      the result polynomial
//
// NOTES:
// ############################# Public Method ###############################
Polynomial& Polynomial::operator -= (const Polynomial& y)
{
  int   i, j, y_order, old_order;
  const double *input_coeff;
  double    *save_coeff;
//
// Get and save y coefficients in case this == y
//
  y_order       = y.getOrder();
  input_coeff   = y.getCoefficients();
  save_coeff    = new double[y_order+1];
  for(i=0; i<y_order+1; i++)
    save_coeff[i]       = input_coeff[i];
//
// update order:
//
  old_order     = order;
  order         = MAX(order, y_order);
  updateMemory(order);
  if(order > old_order)
  {
    j           = order;
    for(i=old_order; i>=0; i--)
    {
      coefficients[j]   = coefficients[i]- save_coeff[j];
      j--;
    }
    for(i=0; i<(order-old_order); i++)
      coefficients[i]   = -save_coeff[i];
  }
  else
  {
    j           = order;
    for(i=y_order; i>=0; i--)
    {
      coefficients[j]   -= save_coeff[i];
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
Polynomial& Polynomial::operator *= (const Polynomial& y)
{
  int       i, j, k, m, n, old_order, y_order;
  const     double *input_coeff;
  double    *temp, *save_coeff;
//
// Get and save y coefficients in case this == y
//
  y_order       = y.getOrder();
  input_coeff   = y.getCoefficients();
  save_coeff    = new double[y_order+1];
  for(i=0; i<y_order+1; i++)
    save_coeff[i]       = input_coeff[i];
//
// Update this order
//
  old_order     = order;
  order         += y_order;
  updateMemory(order);
  temp          = new double[old_order+1];
  for(i=0; i<old_order+1; i++)
    temp[i]     = coefficients[i];                  // Save old coefficients
  for(i=0; i<=order; i++)
  {
    coefficients[i] = 0;
    m               = MIN(i, old_order);
    for(j=0; j<=m; j++)
    {
      n             = MIN(i, y_order);
      for(k=0; k<=n; k++)
      {
        if(j+k == i)
          coefficients[i] += temp[j]*save_coeff[k];
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
Polynomial& Polynomial::operator /= (const Polynomial& y)
{
  int           i, y_order, number_shifts;
  static        int old_shifts = -1;
  static        double *dividend;
  const         double *shifted_coeff;
  double        shifted_xn_coeff;
  Polynomial    sub_poly;

// First check for boundary case
  y_order       = y.getOrder();
  if(y_order > order)
  {
    order           = 0;
    coefficients[0] = 0;
    return *this;
  }

//
// Allocate the new poly coefficients
//
  number_shifts     = order - y_order;
  if(old_shifts < number_shifts)                            // Reduce calls to new
  {
    delete [] dividend;
    old_shifts  = number_shifts;
    dividend    = new double[number_shifts+1];
  }
//
// Now, do long division
//
  Polynomial    shifted_poly = y << number_shifts;
  shifted_coeff     = shifted_poly.getCoefficients();
  shifted_xn_coeff  = shifted_coeff[0];
  if(shifted_xn_coeff == 0.)
    return *this;                                       // Error, shouldn't happen
  for(i=0; i<number_shifts+1; i++)
  {
    if(order < y_order)                                 // Time to exit?
      break;
    if(order >= shifted_poly.getOrder())                // Is subtract possible
    {
      dividend[i]   = coefficients[0]/shifted_xn_coeff;
      sub_poly      = shifted_poly * dividend[i];
      *this         -= sub_poly;                        // Yes update this by subtracting
    }
    else
      dividend[i]   = 0;
    shifted_poly >>= 1;
  }
//
// Modify this, and return it
//
  order = number_shifts;
  for(i=0; i<order+1; i++)
    coefficients[i] = dividend[i];
    
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
Polynomial& Polynomial::operator %= (const Polynomial& y)
{
  Polynomial    div;

// First check for boundary case
  if(order < y.getOrder())
  {
    return *this;
  }
//
// Modulo equals remainder after division:
//
  div   = *this / y;
  div   *= y;
  *this -= div;
    
  return *this;
}

// ############################# Public Method ###############################
// Operator >>=  Shift the polynomial right by n shifts.
// Input:           numberShifts:   Number of places to shift right
//          
// Output:                          the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
Polynomial& Polynomial::operator >>= (int numberShifts)
{
//
// Bound the # of shifts, and update the order
//
  numberShifts  = MAX(0, numberShifts);
  numberShifts  = MIN(numberShifts, MAX_SHIFT);
  order         -= numberShifts;
  order         = MAX(0, order);
  return *this;
}

// ############################# Public Method ###############################
// Operator <<=  Shift the polynomial left by n shifts.
// Input:           numberShifts:   Number of places to shift left
//          
// Output:                          the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
Polynomial& Polynomial::operator <<= (int numberShifts)
{
  int   i, old_order;

//
// Bound the # of shifts, and update the order
//
  numberShifts  = MAX(0, numberShifts);
  numberShifts  = MIN(numberShifts, MAX_SHIFT);
  old_order     = order;
  order         += numberShifts;
  updateMemory(order);
//
// Do the shifting by filling in zeros:
//
  for(i=old_order+1; i<(order+1); i++)
    coefficients[i] = 0;
  updateOrder();                                    // Remove leading zeros
  return *this;
}

// ############################# Public Method ###############################
// Operator +=  Adds a double to the polynomial.
// Input:           y:              Double to add
//          
// Output:                          the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
Polynomial& Polynomial::operator += (double y)
{
  coefficients[order]   += y;
  return *this;
}

// ############################# Public Method ###############################
// Operator -=  Subtracts a double from the polynomial.
// Input:           y:              Double to subtract
//          
// Output:                          the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
Polynomial& Polynomial::operator -= (double y)
{
  coefficients[order]   -= y;
  return *this;
}

// ############################# Public Method ###############################
// Operator *=  Multiplies a double by the polynomial.
// Input:           y:              Double to mult
//          
// Output:                          the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
Polynomial& Polynomial::operator *= (double y)
{
  for(int i=0; i<=order; i++)
    coefficients[i] *= y;
  return *this;
}

// ############################# Public Method ###############################
// Operator /=  Divides the polynomial by a double
// Input:           y:              Double to divide into
//          
// Output:                          the result polynomial
//
// NOTES:
//   1. p(x) = a[0]*x^n + a[1]*x^(n-1) + .... + a[n]
//
// ############################# Public Method ###############################
Polynomial& Polynomial::operator /= (double y)
{
  if(y != 0.)
  {
    for(int i=0; i<=order; i++)
      coefficients[i]   /= y;
  }
  return *this;
}


#ifdef  USE_IO_STREAM                           // Won't compile on NT *******************************
// ############################# Public Method ###############################
// <<   - Print the polynomial
// Input:           s:              ofstream for printing
//                  x:              the polynomial to print
//          
// Output:                          None
//
// ############################# Public Method ###############################
ostream& operator << (ostream& s, const Polynomial& x)
{
  return s << x.getPolyString() << "\n";
}
#endif                                          // ***************************************************

// ############################# Public Method ###############################
// Operator ==  Checks whether the two polynomials are equal
// Input:           x:          first polynomial
//                  y:          second polynomial
//          
// Output:                      YES if two polynomials are equal, else NO
//
// NOTES:
// 1. This is probably not very useful for floating point coefficients
// ############################# Public Method ###############################
LOGICAL operator == (const Polynomial& x, const Polynomial& y)
{
  int   i, n;
  const double *x_coeff, *y_coeff;
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
LOGICAL operator != (const Polynomial& x, const Polynomial& y)
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
LOGICAL operator > (const Polynomial& a, const Polynomial& b)
{
  int   i, a_order, b_order;
  const double *a_coeff, *b_coeff;
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
LOGICAL operator >= (const Polynomial& a, const Polynomial& b)
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
LOGICAL operator < (const Polynomial& a, const Polynomial& b)
{
  int   i, a_order, b_order;
  const double *a_coeff, *b_coeff;
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
LOGICAL operator <= (const Polynomial& a, const Polynomial& b)
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
// ############################# Public Method ###############################
Polynomial operator - (const Polynomial& x)
{
  return (-1.*x);
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
Polynomial operator >> (const Polynomial& x, int numberShift)
{
  int       i, new_order;
  static    int old_order = 0;
  static    double *shift   = NULL;
  const     double *input_coeff;
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
    shift       = new double[new_order+1];
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
    
  return Polynomial(new_order, shift);
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
Polynomial operator << (const Polynomial& x, int numberShift)
{
  int       i, new_order;
  static    int old_order = 0;
  static    double *shift   = NULL;
  const     double *input_coeff;
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
    shift       = new double[new_order+1];
  }
  input_coeff   = x.getCoefficients();
  for(i=0; i<=x.getOrder(); i++)
    shift[i]    = input_coeff[i];
  for(i=(x.getOrder()+1); i<=new_order; i++)
    shift[i]    = 0.;
    
  return Polynomial(new_order, shift);
}

// ############################# Public Method ###############################
// Operator +  Returns the sum of the two input polynomials
// Input:           x:          input polynomial
//                  y:          second polynomial
//          
// Output:                      the sum polynomial
//
// ############################# Public Method ###############################
Polynomial operator + (const Polynomial& x, const Polynomial& y)
{
  Polynomial    sum = x;
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
Polynomial operator - (const Polynomial& x, const Polynomial& y)
{
  Polynomial    difference = x;
  difference    -= y;                   //Stroustrup p. 302
  return difference;
}

// ############################# Public Method ###############################
// Operator *  Returns the product of the two input polynomials
// Input:           x:          input polynomial
//                  y:          second polynomial
//          
// Output:                      the product polynomial
//
// ############################# Public Method ###############################
Polynomial operator * (const Polynomial& x, const Polynomial& y)
{
  Polynomial    prod = x;
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
Polynomial operator / (const Polynomial& x, const Polynomial& y)
{
  Polynomial    result = x;
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
Polynomial operator % (const Polynomial& x, const Polynomial& y)
{
  Polynomial    result = x;
  result    %= y;                           //Stroustrup p. 302
  return    result;
}

// ############################# Public Method ###############################
// Operator +  Returns the sum of a double and a polynomial
// Input:           x:          input double
//                  y:          input polynomial
//          
// Output:                      the sum polynomial
//
// ############################# Public Method ###############################
Polynomial operator + (double x, const Polynomial &y)
{
  Polynomial    sum = y;
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
Polynomial operator - (double x, const Polynomial &y)
{
  Polynomial    result = -y;
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
Polynomial operator * (double x, const Polynomial &y)
{
  Polynomial    product = y;
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
Polynomial operator / (double x, const Polynomial &y)
{
  Polynomial    dividend = x;
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
Polynomial operator + (const Polynomial &x, double y)
{
  Polynomial    sum = x;
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
Polynomial operator - (const Polynomial &x, double y)
{
  Polynomial    difference = x;
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
Polynomial operator * (const Polynomial &x, double y)
{
  Polynomial    product = x;
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
Polynomial operator / (const Polynomial &x, double y)
{
  Polynomial    dividend = x;
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
Polynomial pow(const Polynomial &x, int n)
{
  Polynomial    pow_result = x;
  Polynomial    temp_pol   = x;
  while( --n > 0 )
    pow_result  *= temp_pol;
  return pow_result;
}
