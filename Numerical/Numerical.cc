/************************************************************************
 *                                                                      *
 * This subclass of object implements certain numerical math methods.   *
 * The methods implemented include a function zero finder, a one        *
 * dimensional integrator, and a two dimensional integrator.            *
 *                                                                      *
 * File:Numerical.cc                                                    *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 01/05/95  - Objective-C -> C++                                   *
 *  2. 08/30/05  - Added functionMin().                                 *
 ************************************************************************/

#include "Numerical.h"
#include <math.h>

#define ABS(a)      ((a) >= 0 ? (a) : (-a))
#define SIGN(a,b)   ( ((b)>=0 ) ? ABS((a)) : -ABS((a)) )
#define MAX(a,b)    ((a)>(b)?(a):(b))

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the Numerical class.
// Input:                           None
//          
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
Numerical::Numerical()
{
  double tol1;

/*******************
 * Find machine eps*
 *******************/
  machineEps = 1.;
  tol1 = 1.1;
  while(tol1 > 1.)
  {
    machineEps /= 2.;
    tol1 = 1. + machineEps;
  }
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the Numerical class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
Numerical::~Numerical()
{
  return;
}

// ############################# Public Method ###############################
// zeroIn -- Find a zero in [ax, bx] of a function supplied by the sender
// Input:           tol:        tolerance on zero finder
//                  ax:         lower limit of search region
//                  bx:         upper limit of search region
//                  function:   pointer to function to find zero
//          
// Output:                      Value of x such that f(x) = 0.
// ############################# Public Method ###############################
double Numerical::zeroIn(double ax, double bx, double (*function)(double), double tol)
{
  int    begin_step, converged, bisection;
  double a, b, c, d, e, fa, fb, fc, tol1;
  double xm, p, q, r, s;

/********************
 * Initialization   *
 ********************/
  c = d = e = fc = 0.;
  a  = ax;
  b  = bx;
  fa = function(a);
  fb = function(b);
  
/*******************
 * Begin step:     *
 *******************/
  converged  = 0;
  begin_step = 1;
  while(!converged)
  {
    if(begin_step)
    {
      c  = a;
      fc = fa;
      d  = b - a;
      e  = d;
    }
    if( ABS(fc) < ABS(fb) )
    {
      a  = b;
      b  = c;
      c  = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
  
/*******************
 * Convergence test*
 *******************/
  tol1 = 2.*machineEps*ABS(b) + 0.5*tol;
  xm   = 0.5*(c-b);
  if( (ABS(xm)<=tol1) || (fb == 0.) )
  {
    converged = 1;
    break;
  }

/********************
 * Bisection        *
 * necessary?       *
 ********************/
  bisection = 1;
  if( (ABS(e)>=tol1) && (ABS(fa)>ABS(fb)) )
  {
  
/********************
 * quadratic interp *
 * possible?        *
 ********************/
     if( a == c)
     {
/********************
 * linear interp.   *
 ********************/
        s = fb/fa;
        p = 2.*xm*s;
        q = 1. - s;
     }
     else
     {
/*******************
 * Inverse quadrat.*
 * interp.         *
 *******************/
       q = fa/fc;
       r = fb/fc;
       s = fb/fa;
       p = s*(2.*xm*q*(q-r) - (b-a)*(r-1.));
       q = (q-1.)*(r-1.)*(s-1.);
     }
/*******************
 * Adjust signs:   *
 *******************/
     q = (p > 0.) ? -q : q;
     p = ABS(p);
/*******************
 * Is interpolation*
 * acceptable?     *
 *******************/
     if(  ( 2.*p < (2.*xm*q - ABS(tol1*q)) )  &&
          ( p < ABS(0.5*e*q) )                 )
     {
       e = d;
       d = p/q;
       bisection = 0;
     }
     else
       bisection = 1;
   }
/*******************
 * Bisection       *
 *******************/
   if(bisection)
   {
     d = xm;
     e = d;
   }
/*******************
 * Complete step   *
 *******************/
     a  = b;
     fa = fb;
     if( ABS(d) > tol1)
       b += d;
     if( ABS(d) <= tol1)
       b += SIGN(tol1, xm);
     fb = function(b);
     if(fb*fc/ABS(fc) > 0.)
       begin_step = 1;
     else
       begin_step = 0;
   }
   return b;
}

#define GSL_SUCCESS     1
#define GSL_FAILURE     0
// ############################# Public Method ###############################
// functionMin -- Find a minimum in [a b] of a function supplied by the sender
// Input:           a:          left border of the range
//                  b:          right border of the range
//                  f:          pointer to function to find min
//                  tol:        tolerance on minimum
//          
// Output:                      Value of x such that f(x) = min.
//      functionMin returns an estimate for the minimum location with accuracy
//      3//SQRT_EPSILON//abs(x) + tol.
//      The function always obtains a local minimum which coincides with
//      the global one only if a function under investigation being
//      unimodular.
//      If a function being examined possesses no local minimum within
//      the given range, Fminbr returns 'a' (if f(a) < f(b)), otherwise
//      it returns the right range boundary value b.
//
// Algorithm
//      G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
//      computations. M., Mir, 1980, p.202 of the Russian edition
//
//      The function makes use of the "gold section" procedure combined with
//      the parabolic interpolation.
//      At every step program operates three abscissae - x,v, and w.
//      x - the last and the best approximation to the minimum location,
//          i.e. f(x) <= f(a) or/and f(x) <= f(b)
//          (if the function f has a local minimum in (a,b), then the both
//          conditions are fulfiled after one or two steps).
//      v,w are previous approximations to the minimum location. They may
//      coincide with a, b, or x (although the algorithm tries to make all
//      u, v, and w distinct). Points x, v, and w are used to construct
//      interpolating parabola whose minimum will be treated as a new
//      approximation to the minimum location if the former falls within
//      [a,b] and reduces the range enveloping minimum more efficient than
//      the gold section procedure. 
//      When f(x) has a second derivative positive at the minimum location
//      (not coinciding with a or b) the procedure converges superlinearly
//      at a rate order about 1.324
// ############################# Public Method ###############################
double Numerical::functionMin(double a, double b, double (*f)(double x), double tol)
{
  double x,v,w;                         /* Abscissae, descr. see above  */
  double fx;                            /* f(x)                         */
  double fv;                            /* f(v)                         */
  double fw;                            /* f(w)                         */
  const double r = (3.-sqrt(5.0))/2;    /* Gold section ratio           */

//  assert( tol > 0 && b > a );

  v = a + r*(b-a);  fv = (*f)(v);       /* First step - always gold section*/
  x = v;  w = v;
  fx=fv;  fw=fv;

  for(;;)               /* Main iteration loop  */
  {
    double range = b-a;                 /* Range over which the minimum */
                                        /* is seeked for                */
    double middle_range = (a+b)/2;
    double tol_act =                    /* Actual tolerance             */
                machineEps*fabs(x) + tol/3;
    double new_step;                    /* Step at this iteration       */

       

    if( fabs(x-middle_range) + range/2 <= 2*tol_act )
      return x;                         /* Acceptable approx. is found  */

                                        /* Obtain the gold section step */
    new_step = r * ( x<middle_range ? b-x : a-x );


                        /* Decide if the interpolation can be tried     */
    if( fabs(x-w) >= tol_act  )         /* If x and w are distinct      */
    {                                   /* interpolatiom may be tried   */
        register double p;              /* Interpolation step is calcula-*/
        register double q;              /* ted as p/q; division operation*/
                                        /* is delayed until last moment */
        register double t;

        t = (x-w) * (fx-fv);
        q = (x-v) * (fx-fw);
        p = (x-v)*q - (x-w)*t;
        q = 2*(q-t);

        if( q>(double)0 )               /* q was calculated with the op-*/
          p = -p;                       /* posite sign; make q positive */
        else                            /* and assign possible minus to */
          q = -q;                       /* p                            */

        if( fabs(p) < fabs(new_step*q) &&       /* If x+p/q falls in [a,b]*/
            p > q*(a-x+2*tol_act) &&            /* not too close to a and */
            p < q*(b-x-2*tol_act)  )            /* b, and isn't too large */
          new_step = p/q;                       /* it is accepted         */
                                        /* If p/q is too large then the */
                                        /* gold section procedure can   */
                                        /* reduce [a,b] range to more   */
                                        /* extent                       */
    }

    if( fabs(new_step) < tol_act )      /* Adjust the step to be not less*/
      if( new_step > (double)0 )        /* than tolerance               */
        new_step = tol_act;
      else
        new_step = -tol_act;

                                /* Obtain the next approximation to min */
    {                           /* and reduce the enveloping range      */
      register double t = x + new_step; /* Tentative point for the min  */
      register double ft = (*f)(t);
      if( ft <= fx )
      {                                 /* t is a better approximation  */
        if( t < x )                     /* Reduce the range so that     */
          b = x;                        /* t would fall within it       */
        else
          a = x;
      
        v = w;  w = x;  x = t;          /* Assign the best approx to x  */
        fv=fw;  fw=fx;  fx=ft;
      }
      else                              /* x remains the better approx  */
      {                              
        if( t < x )                     /* Reduce the range enclosing x */
          a = t;                   
        else
          b = t;
      
        if( ft <= fw || w==x )
        {
           v = w;  w = t;
           fv=fw;  fw=ft;
        }
        else if( ft<=fv || v==x || v==w )
        {
           v = t;
           fv=ft;
        }
      }
      
    }                   /* ----- end-of-block ----- */
  }             /* ===== End of loop ===== */

}


#define NUMSEG  100
// ############################# Public Method ###############################
// integrateFrom -- This a Simpson's integrator from Quinn Curtis
// Input:           ax:         lower limit of integration
//                  bx:         upper limit of integration
//                  function:   pointer to function to integrate
//          
// Output:                      Integration value
// ############################# Public Method ###############################
double Numerical::integrateFrom(double  ax, double bx, double (*integrand)(double))
{
  int    numseg;
  double xl, xh;
  double strtpnt, endpnt, area1;

  double area2, area, segwidth;


   xl = ax;
   xh = bx;
   numseg = NUMSEG;
   segwidth = (xh - xl) / numseg;

   if(segwidth<=machineEps) 
    return 0.;
   endpnt = xh;

   area1 = 0.0;

   area2 = 0.0;

   area = 0.0;

   if ( (numseg%2) != 0 ) {

      area1 = 3.0 / 8.0 * segwidth * 
             (      integrand(endpnt - 3.0 * segwidth) +

              3.0 * integrand(endpnt - 2.0 * segwidth) +

              3.0 * integrand(endpnt - segwidth) + 
                    integrand(endpnt)
             );

      endpnt = endpnt - 3.0 * segwidth;

   }

   else

      area1 = 0.0;



   if ( numseg != 3 ) {

      strtpnt = xl;

      while ( (strtpnt < endpnt - segwidth) )
      {

         area2 += 1.0 / 3.0 * (segwidth * 
                 (      integrand(strtpnt) + 
                  4.0 * integrand(strtpnt + segwidth) + 
                        integrand(strtpnt + 2.0 * segwidth)
                 ));

         strtpnt += 2.0 * segwidth;
     }

   }

   area = area1 + area2;

   return area;

}



