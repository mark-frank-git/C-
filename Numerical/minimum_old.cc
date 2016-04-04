// ############################# Public Method ###############################
// functionMin -- Find a minimum in [ax, bx] of a function supplied by the sender
// Input:           tol:        tolerance on zero finder
//                  ax:         lower limit of search region
//                  bx:         upper limit of search region
//                  function:   pointer to function to find zero
//          
// Output:                      Value of x such that f(x) = 0.
// ############################# Public Method ###############################
int Numerical::brent_iterate (double (*function)(double), double *x_minimum, double * f_minimum, double x_lower, double x_upper)
{
  brent_state_t *state = (brent_state_t *) vstate;

  const double x_left   = x_lower;
  const double x_right  = x_upper;
  const double golden   = 0.3819660;      /* golden = (3 - sqrt(5))/2 */

  double v = x_lower + golden * (x_upper - x_lower);
  double w = v;

  const double z        = *x_minimum    = x_lower;
  double d = 0.;
  double e = 0.;
  double u, f_u;
  const double f_v      = function(v);
  const double f_w      = f_v;
  const double f_z      = *f_minimum;

  double f_lower        = function(x_lower);
  double f_upper        = function(x_upper);

  const double w_lower  = (z - x_left);
  const double w_upper  = (x_right - z);

  const double tolerance =  machineEps * fabs (z);

  double p = 0, q = 0, r = 0;

  const double midpoint = 0.5 * (x_left + x_right);

  if (fabs (e) > tolerance)
    {
      /* fit parabola */

      r = (z - w) * (f_z - f_v);
      q = (z - v) * (f_z - f_w);
      p = (z - v) * q - (z - w) * r;
      q = 2 * (q - r);

      if (q > 0)
        {
          p = -p;
        }
      else
        {
          q = -q;
        }

      r = e;
      e = d;
    }

  if (fabs (p) < fabs (0.5 * q * r) && p < q * w_lower && p < q * w_upper)
    {
      double t2 = 2 * tolerance ;

      d = p / q;
      u = z + d;

      if ((u - x_left) < t2 || (x_right - u) < t2)
        {
          d = (z < midpoint) ? tolerance : -tolerance ;
        }
    }
  else
    {
      e = (z < midpoint) ? x_right - z : -(z - x_left) ;
      d = golden * e;
    }


  if (fabs (d) >= tolerance)
    {
      u = z + d;
    }
  else
    {
      u = z + ((d > 0) ? tolerance : -tolerance) ;
    }

  state->e = e;
  state->d = d;

  f_u   = function(u);

  if (f_u > f_z)
    {
      if (u < z)
        {
          *x_lower = u;
          *f_lower = f_u;
          return GSL_SUCCESS;
        }
      else
        {
          *x_upper = u;
          *f_upper = f_u;
          return GSL_SUCCESS;
        }
    }
  else if (f_u < f_z)
    {
      if (u < z)
        {
          *x_upper = z;
          *f_upper = f_z;
        }
      else
        {
          *x_lower = z;
          *f_lower = f_z;
        }

      *x_minimum = u;
      *f_minimum = f_u;
      return GSL_SUCCESS;
    }
  else if (f_u <= f_w || w == z)
    {
      state->v = w;
      state->f_v = f_w;
      state->w = u;
      state->f_w = f_u;
      return GSL_SUCCESS;
    }
  else if (f_u <= f_v || v == z || v == w)
    {
      state->v = u;
      state->f_v = f_u;
      return GSL_SUCCESS;
    }
  else
    {
      return GSL_FAILURE;
    }
}
