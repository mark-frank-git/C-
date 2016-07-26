#ifndef _NUMERICAL_H
#define _NUMERICAL_H    1
/************************************************************************
s *                                                                     *
 * This subclass of object implements certain numerical math methods.   *
 * The methods implemented include a function zero finder, and a        *
 * Simpson's integrator.                                                *
 *                                                                      *
 * File:Numerical.h                                                     *
 *                                                                      *
 ************************************************************************/

class Numerical
{
private:
  double  machineEps;           /* Machine epsilon */

public:
/**********************
 * Constructors:    *
 **********************/
  Numerical();
  ~Numerical();
        
/**********************
 * Numerical functions*
 **********************/
  double zeroIn(double ax, double bx, double (*function)(double), double tol=1.e-5);
                                                            // Find function zero in an interval
  double functionMin(double a, double b, double (*f)(double x), double tol=1.e-5);
  double integrateFrom(double  ax, double bx, double (*integrand)(double));
                                                            // Integrate a function
/**********************
 * Set parameters:    *
 **********************/

/**********************
 * Get parameters:    *
 **********************/


};
#endif
