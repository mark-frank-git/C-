/************************************************************************
 *                                                                      *
 *  This file contains subroutines for performing Fourier transfomrs.   *
 *                                                                      *
 *  Author: Mark Frank                                                  *
 *                                                                      *
 *  File: fft_routines.c                                                *
 *                                                                      *
 *************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <GNU/Complex.h>
#include "fft_routines.h"

#ifndef TWOPI
#define TWOPI                   6.283185307             /* 2*PI         */
#define PI                      3.141592654
#endif

#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))

#ifndef EPS
#define EPS     0.0001
#endif

/************************************************************************
 *                                                                      *
 * void complex_fft(Complex *x, int m, int operation)                   *
 *                                                                      *
 *                                                                      *
 *    integer input variables                                           *
 *    -----------------------                                           *
 *    m = log2(number of points in array)                               *
 *    operation = 1 for FFT, -1 for IFFT                                *
 *                                                                      *
 *    float input variables                                             *
 *    ---------------------                                             *
 *    x = floating point complex array                                  *
 *                                                                      *
 *                                                                      *
 * Requires pointer to complex array, x and power of 2 size of FFT, m   *
 * (size of FFT = 2**m).  Places FFT output on top of input array.      *
 ** Modified from "C Language Algorithms for Digital Signal Processing,"*
 * by Embree and Kimble                                                 *
 ************************************************************************/
void complex_fft(Complex *x, int m, int operation)
{
  static Complex *w     = NULL;         // Used to store the w complex array
  static int mstore     = 0;            // stores m for future reference
  static int opstore    = 0;
  static int n          = 0;

  Complex   u,temp,tm;
  Complex   *xi,*xip,*xj,*wptr;

  int       i,j,k,l,le,windex;

  double    scale;
  double    arg,w_real,w_imag,wrecur_real,wrecur_imag,wtemp_real;
    
//
// The following gets run if the size of the fft has
// changed or we want a new operation.
//
  if((m != mstore) || (operation != opstore))
  {
    delete [] w;
    opstore = operation;
    mstore  = m;
    if(m == 0) 
      return;                               /* if m=0 then done */

//
// n = 2**m = fft length
//
    n   = 1 << m;
    le  = n/2;

//
// allocate the storage for w
//
    w   = new Complex[le];

/* calculate the w values recursively */

    arg             = 4.0*atan(1.0)/le;         /* PI/le calculation */
    wrecur_real     = w_real = cos(arg);
    if(operation == INVERSE_FFT)
      wrecur_imag   = w_imag = sin(arg);
    else
      wrecur_imag   = w_imag = -sin(arg);
    xj = w;
    for (j = 1 ; j < le ; j++)
    {
      *xj           = Complex(wrecur_real, wrecur_imag);
      xj++;
      wtemp_real    = wrecur_real*w_real - wrecur_imag*w_imag;
      wrecur_imag   = wrecur_real*w_imag + wrecur_imag*w_real;
      wrecur_real   = wtemp_real;
    }
  }

/* start fft */

  le        = n;
  windex    = 1;
  for (l = 0 ; l < m ; l++)
  {
    le >>= 1;

/* first iteration with no multiplies */

    for(i = 0 ; i < n ; i += 2*le)
    {
      xi        = x + i;
      xip       = xi + le;
      temp      = *xi + *xip;
      *xip      = *xi - *xip;
      *xi       = temp;
    }

/* remaining iterations use stored w */

    wptr = w + windex - 1;
    for (j = 1 ; j < le ; j++)
    {
      u = *wptr;
      for (i = j ; i < n ; i += 2*le)
      {
        xi      = x + i;
        xip     = xi + le;
        temp    = *xi + *xip;
        tm      = *xi - *xip;
        *xip    = tm * u;
        *xi     = temp;
      }
      wptr += windex;
    }
    windex *= 2;
  }      

/* rearrange data by bit reversing */

  j = 0;
  for (i = 1 ; i < (n-1) ; i++)
  {
    k   = n/2;
    while(k <= j)
    {
      j -= k;
      k >>= 1;                  // Divide by 2
    }
    j += k;
    if (i < j)
    {
      xi    = x + i;
      xj    = x + j;
      temp  = *xj;
      *xj   = *xi;
      *xi   = temp;
    }
  }
/************************
 * If IFFT scale by     *
 * 1/n                  *
 ************************/
  if( (operation == INVERSE_FFT) && (n>0) )
  {
     scale  = (double) 1.0/n;
     for(i=0; i<n; i++)
     {
       *x   *= scale;
       x++;
     }
  }
  return;
       
}

/************************************************************************
 *                                                                      *
 * void fft_real(Complex *x, int m)                                     *
 *                                                                      *
 *                                                                      *
 *    integer input variables                                           *
 *    -----------------------                                           *
 *    m = log2(number of points in array)                               *
 *                                                                      *
 *    Complex input/output variables                                    *
 *    ---------------------                                             *
 *    x = complex array = {real[0] + j real[1]}, {real[2] + j real[3]}, *
 *                      = {real[0] + j imag[0]}, etc. on output.        *
 *                                                                      *
 **** Modified from "C Language Algorithms for Digital Signal           *
 * Processing," by Embree and Kimble                                    *
 *                                                                      *
 * NOTE: the original routine did not return value for k = N/2??        *
 *                                                                      *
 **************************************************************************/
void fft_real(Complex *x, int m)
{
  static    Complex  *cf = NULL;
  static    Complex  *y  = NULL;
  static    int      mstore = 0;
  int       p,num,k,num_output;
  int       i, j, n;
  float     Realsum, Realdif, Imagsum, Imagdif;
  double    factor, arg;
  Complex   *ck, *xk, *xnk;

// First call the fft routine using the x array but with
// half the size of the real fft.
  n     = 1 << m;
  p     = m - 1;
  complex_fft(x, p, FORWARD_FFT);

//
// Next create the coefficients for recombination, if required
//
  num           = 1 << p;    /* num is half the real sequence length.  */
  num_output    = num+1;

  if (m!=mstore)
  {
    mstore = m;
    delete [] cf;
    delete [] y;
    cf      = new Complex[num];
    y       = new Complex[num_output];

    factor = 4.0*atan(1.0)/num;
    for (k = 1; k <= num; k++)
    {
      arg       = factor*k;
      cf[k-1]   = Complex(cos(arg), sin(arg));
    }
  }  

//
// DC component, no multiplies
//
  y[0]      = Complex(real(x[0]) + imag(x[0]), 0.0);
//
// N/2 component
//
  y[num]    = Complex(real(x[0]) - imag(x[0]), 0.);

//
// other frequencies by trig recombination
//
  ck        = cf;
  xk        = x + 1;
  xnk       = x + num - 1;
  for (k = 1; k < num; k++)
  {
    Realsum = ( real(*xk) + real(*xnk) );
    Imagsum = ( imag(*xk) + imag(*xnk) );
    Realdif = ( real(*xk) - real(*xnk) );
    Imagdif = ( imag(*xk) - imag(*xnk) );

    y[k]    = Complex(Realsum + real(*ck) * Imagsum - imag(*ck) * Realdif, 
                      Imagdif - imag(*ck) * Imagsum - real(*ck) * Realdif)/2.;
    ck++;
    xk++;
    xnk--;
  }
//
// Store results in output and reflect across num_output
//
  for(i=0; i<num_output; i++)
    x[i]    = y[i];
  j         = num;
  for(i=num; i<n; i++)              // spectrum is reflected about num_output
    x[i]    = x[j--];
  return;
}

/****************************************************************************
 * void complex_dft(COMPLEX *data_in, COMPLEX *data_out, int number_pts,    *
 *          int operation)                                                  *
 *                                                                          *
 *    integer input variables                                               *
 *    -----------------------                                               *
 *    number_pts = number input points                                      *
 *    operation = 1 for DFT, -1 for IDFT                                    *
 *                                                                          *
 *    COMPLEX input variables                                               *
 *    -----------------------------                                         *
 *    data_in[number_pts] = double Complex input array                      *
 *                                                                          *
 *    COMPLEX output variables                                              *
 *    -----------------------------                                         *
 *    data_out[number_pts] = double Complex output array                    *
 *                                                                          *
 * This function performs a straight DFT of N points on an array of         *
 * complex numbers whose first member is pointed to by data_in.  The        *
 * output is placed in an array pointed to by data_out.                     *
 **** Modified from "C Language Algorithms for Digital Signal Processing,"  *
 * by Embree and Kimble                                                     *
 ****************************************************************************/
void complex_dft(Complex *data_in, Complex *data_out, int number_pts,
          int operation)
{
  int       i,k,n,p;
  static    int nstore = 0;    /* store N for future use */
  static    int optstore = 1;
  static    Complex *cf;     /* coefficient storage */
  Complex   *cfptr,*Dinptr;
  double    arg, arg_i;

/* Create the coefficients if N has changed */

  if( (number_pts != nstore) || (optstore != operation) )
  {
    if(nstore != 0)
    delete [] cf;                                       /* free previous */
    nstore      = number_pts;
    optstore    = operation;

    cf          = new Complex[number_pts];

    arg         = 8.0*atan(1.0)/(double)number_pts;
    arg_i       = 0.;
    for (i=0 ; i<number_pts ; i++)
    {
      if(operation == FORWARD_FFT)
        cf[i]   = Complex(cos(arg_i), -sin(arg_i));
      else
        cf[i]   = Complex(cos(arg_i), -sin(arg_i))/(double)number_pts;
      arg_i     += arg;
    }
  }

/* Perform the DFT calculation */

  for (k=0 ; k<number_pts ; k++) 
  {
    Dinptr      = data_in;
    if(operation == FORWARD_FFT)
    {
      *data_out = *Dinptr;
    }
    else
    {
      *data_out = *Dinptr * cf[0].real();
    }
    Dinptr++;
    for (n=1; n<number_pts; n++) 
    {
       p            = (n*k) % number_pts;
       cfptr        = cf + p;                       /* pointer to cf modulo N */
       *data_out    += *Dinptr * (*cfptr);
       Dinptr++;
    }
    data_out++;                                     /* next output */
  }
}

#define FORWARD_DCT 0
#define INVERSE_DCT 1
/**************************************************************************
 *                                                                       *
 * void dct(float *data_in, float *data_out, int number_pts,             *
 *           int number_output)                                           *
 *                                                                        *
 *    integer input variables                                             *
 *    -----------------------                                             *
 *    number_pts = number input points                                    *
 *    number_output = number output points                                *
 *                                                                        *
 *    float input variables                                               *
 *    -----------------------------                                       *
 *    data_in[number_pts] = Real input array                              *
 *                                                                        *
 *    float output variables                                              *
 *    -----------------------------                                       *
 *    data_out[number_pts] = float  output array                          *
 *                                                                        * 
 * This function performs a straight DCT of N points on an array of      *
 * float   numbers whose first member is pointed to by data_in.  The     *
 * output is placed in an array pointed to by data_out.                   *
 **** Modified from "C Language Algorithms for Digital Signal Processing,"*
 * by Embree and Kimble                                                   *
 **************************************************************************/
void dct(float *data_in, float *data_out, int number_pts,
          int number_output, int dct_type)
{
    int i,k,n,p;
    static int nstore = 0;      /* store N for future use */
    static float *cf;         /* coefficient storage */
    float *cfptr,*Dinptr;
    double arg;

/* Create the coefficients if N has changed */

    if( (number_pts != nstore) )
    {
      if(nstore != 0)
        free((char *) cf);    /* free previous */
      nstore = number_pts;

      cf = (float  *) calloc(number_pts, sizeof(float));

      arg = TWOPI/(double)number_pts;
      if(dct_type == FORWARD_DCT)
      {
       for (i=0 ; i<number_pts ; i++)
        cf[i] = 2.*cos(arg*i);
      }
      else
      {
       for (i=0 ; i<number_pts ; i++)
        cf[i] = 2.*(float)cos(arg*i)/number_pts;
      }
    }

/* Perform the DFT calculation */
    for (k=0 ; k<number_output ; k++)
    {
      Dinptr = data_in;
      *data_out = *Dinptr;
      Dinptr++;
      for (n=1; n<number_pts; n++)
      {
        p = (n*k) % number_pts;
        cfptr = cf + p;         /* pointer to cf modulo N */

        *data_out += (*Dinptr) * (*cfptr);
        Dinptr++;
      }
      data_out++;          /* next output */
    }
}

/************************************************************************
 *                                                                      *
 * void real_czt(float *input_data, int number_input_pts,               *
 *                float omega_start, float omega_end,                   *
 *                Complex *output_data, int number_output_pts)          *
 *                                                                      *
 *                                                                      *
 *    float input variables                                             *
 *    -----------------------                                           *
 *    input_data:       array of input data samples                     *
 *    omega_start:      starting point of evaluation [0, 2*PI]          *
 *    omega_end:        ending point of evaluation [omega_start, 2*PI]  *
 *                                                                      *
 *    int input variables                                               *
 *    ---------------------                                             *
 *    number_input_pts:     Number of points in input array             *
 *    number_output_pts:    Number of points to evaluate chirp Z        *
 *                                                                      *
 *    Complex output variables                                          *
 *    ---------------------                                             *
 *    output_data:          Complex DFT over [omega_start, omega_end].  *
 *                          This array must be allocated in calling     *
 *                          function.
 *                                                                      *
 *                                                                      *
 *  Based on the algorithm given in "Digital Signal Analysis," by       *
 *  Stearns and Hush.                                                   *
 ************************************************************************/
void real_czt(const float *input_data, int number_input_pts, float omega_start, float omega_end,
                    Complex *output_data, int number_output_pts)
{
//
// We use static variables for efficient computation when multiple
// calls to this function are made with the same input paramters
//
  static int n_store            = 0;
  static int m_store            = 0;
  static int l_store            = 0;

  static float theta_0_store    = -1000.;
  static float theta_1_store    = -1000.;

  static Complex *AnBn          = NULL;
  static Complex *Gm            = NULL;
  static Complex *Bm            = NULL;

  int       n, m;
  float     phi_0;
  Complex   A, B, A_minus_n, sqrt_B, B_n_over_2;
  Complex   B_l_sq_over_2, B_to_l, B_to_nl;


//
// phi_0 is the output frequency spacing:
// Probably should error check omega_start and omega_end
//
  number_output_pts             = MAX(2, number_output_pts);
  omega_start                   = MAX((-PI), omega_start);
  omega_end                     = MIN(omega_end, TWOPI);
  phi_0                         = (omega_start - omega_end)/(number_output_pts-1);
  A                             = exp( Complex(0., theta_0_store));
  B                             = exp( Complex(0., -phi_0) );
//
// First check if we need to re-calculate stored values:
//
  if( (n_store != number_input_pts) || (number_output_pts != m_store) || 
      ( ABS( (theta_0_store-omega_start) ) > EPS ) || ( ABS( (theta_1_store-omega_end) ) > EPS )   )
  {
    delete [] AnBn;
    delete [] Gm;
    delete [] Bm;
    n_store         = number_input_pts;
    m_store         = number_output_pts;
    theta_0_store   = omega_start;
    theta_1_store   = omega_end;
//
// Calculate l, a power of 2
//
    
    l_store         = log2(m_store + n_store - 1);
    l_store         = 1 << l_store;
//
// Allocate the arrays:
//
    AnBn            = new Complex[n_store];
    Gm              = new Complex[l_store];
    Bm              = new Complex[l_store];
//
// Initialize the AnBn array:
//
    A               = exp( Complex(0., theta_0_store));
    A_minus_n       = Complex(1., 0.);
    sqrt_B          = sqrt(B);
    B_n_over_2      = Complex(1., 0.);
    for(n=0; n<n_store; n++)
    {
      AnBn[n]       = A_minus_n*pow(B_n_over_2, n);
      A_minus_n     /= A;
      B_n_over_2    *= sqrt_B;
    }
//
// Initialize the Bm array:
//
    B_n_over_2      = Complex(1., 0.);          // This will actually be B^-n_sq/2
    for(n=0; n< m_store; n++)
    {
      Bm[n]         = pow(B_n_over_2, n);
      B_n_over_2    /= sqrt_B;
    }
    for(n=m_store; n<=(l_store-n_store); n++)
      Bm[n]         = Complex(0., 0.);
    B_n_over_2      = pow(sqrt_B, -n);          // This will actually be B^-n_sq/2
    B_l_sq_over_2   = pow(B, -l_store*l_store/2);
    B_to_l          = pow(B, l_store);
    B_to_nl         = pow(B_to_l, n);
    for(; n<l_store; n++)
    {
      Bm[n]         = B_l_sq_over_2 * B_to_nl * pow(B_n_over_2, n);
      B_n_over_2    /= sqrt_B;
      B_to_nl       *= B_to_l;
    }
    m               = log2(l_store);
    complex_fft(Bm, m, FORWARD_FFT);
  }
//
// Multiply the input by the AnBn array to form the gn array
// Then, FFT to compute Gm
//
  for(n=0; n<n_store; n++)
  {
    Gm[n]   = AnBn[n]*input_data[n];
  }
  for(; n<l_store; n++)
    Gm[n]   = Complex(0., 0.);
  m         = log2(l_store);
  complex_fft(Gm, m, FORWARD_FFT);
//
// Multiply the FFTs together and invert to perform the convolution:
//
  for(n=0; n<l_store; n++)
  {
    Gm[n]   *= Bm[n];
  }
  complex_fft(Gm, m, INVERSE_FFT);
//
// Finally, multiply by B^(n/2) to get result
//
  for(n=0; n<m_store; n++)
  {
    output_data[n]  = Gm[n]*pow(B, (double)n*n/2);
  }

  return;
}

/************************************************************************
 *                                                                      *
 * void complex_czt(Complex *input_data, int number_input_pts,          *
 *                  float omega_start, float omega_end,                 *
 *                  Complex *output_data, int number_output_pts)        *
 *                                                                      *
 *    float input variables                                             *
 *    -----------------------                                           *
 *    input_real:       real input array                                *
 *    input_imag:       imag input array                                *
 *    omega_start:      starting point of evaluation [0, 2*PI]          *
 *    omega_end:        ending point of evaluation [omega_start, 2*PI]  *
 *                                                                      *
 *    int input variables                                               *
 *    ---------------------                                             *
 *    number_input_pts:     Number of points in input array             *
 *    number_output_pts:    Number of points to evaluate chirp Z        *
 *                                                                      *
 *    Complex output variables                                          *
 *    ---------------------                                             *
 *    output_data:          Complex DFT over [omega_start, omega_end].  *
 *                          This array must be allocated in calling     *
 *                          function.
 *                                                                      *
 *                                                                      *
 *  Based on the algorithm given in "Digital Signal Analysis," by       *
 *  Stearns and Hush.                                                   *
 ************************************************************************/
void complex_czt(const float *input_real, const float *input_imag, int number_input_pts, float omega_start,
                 float omega_end, Complex *output_data, int number_output_pts)
{
//
// We use static variables for efficient computation when multiple
// calls to this function are made with the same input paramters
//
  static int n_store            = 0;
  static int m_store            = 0;
  static int l_store            = 0;

  static float theta_0_store    = -1000.;
  static float theta_1_store    = -1000.;

  static Complex *AnBn          = NULL;
  static Complex *Gm            = NULL;
  static Complex *Bm            = NULL;

  int           n, m;
  float         phi_0;
  double        theta, phi;
  Complex       A, B, A_minus_n, sqrt_B, B_n_over_2;
  Complex       B_l_sq_over_2, B_to_l, B_to_nl;


//
// phi_0 is the output frequency spacing:
// Probably should error check omega_start and omega_end
//
  number_output_pts             = MAX(2, number_output_pts);
  omega_start                   = MAX((-PI), omega_start);
  omega_end                     = MIN(omega_end, TWOPI);
  phi_0                         = (omega_end - omega_start)/(number_output_pts-1);
  A                             = exp( Complex(0., theta_0_store));
  B                             = exp( Complex(0., -phi_0) );
  sqrt_B                        = exp( Complex(0., -phi_0/2.));
//
// First check if we need to re-calculate stored values:
//
  if( (n_store != number_input_pts) || (number_output_pts != m_store) || 
      ( ABS( (theta_0_store-omega_start) ) > EPS ) || ( ABS( (theta_1_store-omega_end) ) > EPS )   )
  {
    delete [] AnBn;
    delete [] Gm;
    delete [] Bm;
    n_store         = number_input_pts;
    m_store         = number_output_pts;
    theta_0_store   = omega_start;
    theta_1_store   = omega_end;
//
// Calculate l, a power of 2
//
    
    l_store         = log2(m_store + n_store - 1);
    l_store         = 1 << l_store;
//
// Allocate the arrays:
//
    AnBn            = new Complex[n_store];
    Gm              = new Complex[l_store];
    Bm              = new Complex[l_store];
//
// Initialize the AnBn array:
//
    theta               = 0.;
    phi                 = 0.;
    for(n=0; n<n_store; n++)
    {
      A_minus_n         = exp(Complex(0., theta));
      B_n_over_2        = exp(Complex(0., phi));
      AnBn[n]           = A_minus_n*B_n_over_2;
      theta             -= theta_0_store;
      if(theta < -TWOPI)
        theta           += TWOPI;
      phi               = -(n+1)*(n+1)*phi_0/2.;
      while(phi < -TWOPI)
        phi             += TWOPI;
    }
//
// Initialize the Bm array:
//
    phi                 = 0.;
    for(n=0; n< m_store; n++)
    {
      Bm[n]             = exp(Complex(0., phi));
      phi               = (n+1)*(n+1)*phi_0/2.;
      while(phi > TWOPI)
        phi             -= TWOPI;
    }
    for(n=m_store; n<=(l_store-n_store); n++)
      Bm[n]         = Complex(0., 0.);
    for(; n<l_store; n++)
    {
      Bm[n]         = exp(Complex(0., (l_store-n)*(l_store-n)*phi_0/2.));
    }
    m               = log2(l_store);
    complex_fft(Bm, m, FORWARD_FFT);
  }
//
// Multiply the input by the AnBn array to form the gn array
// Then, FFT to compute Gm
//
  for(n=0; n<n_store; n++)
  {
    Gm[n]   = AnBn[n]*Complex(input_real[n], input_imag[n]);
  }
  for(; n<l_store; n++)
    Gm[n]   = Complex(0., 0.);
  m         = log2(l_store);
  complex_fft(Gm, m, FORWARD_FFT);
//
// Multiply the FFTs together and invert to perform the convolution:
//
  for(n=0; n<l_store; n++)
  {
    Gm[n]   *= Bm[n];
  }
  complex_fft(Gm, m, INVERSE_FFT);
//
// Finally, multiply by B^(n/2) to get result
//
  for(n=0; n<m_store; n++)
  {
    output_data[n]  = Gm[n]*exp(Complex(0., -n*n*phi_0/2.));
  }

  return;
}

/**********************************************************************
 * int log2(int n)                                                      *
 *                                                                      *
 * log2 - base 2 logarithm                                              *
 *                                                                      *
 * Returns base 2 log such that n = 2**ans where ans = log2(n).         *
 * if log2(i) is between two values, the larger is returned.            *
 *                                                                      *
 **********************************************************************/
int log2(int n)
{
  unsigned int mask,i;

  if(n == 0) 
    return -1;                /* zero is an error, return -1 */

  n--;                        /* get the max index, x-1     */

  for(mask = 1 , i = 0 ; ; mask *= 2 , i++)
  {
    if(n == 0) 
      return i;             /* return log2 if all zero   */
    n &= (~mask);           /* AND off a bit             */
  }
  return n;
}

// ############################ Public Function ###################################
// powerOfTwo - Returns whether or not the input integer is a power of 2.
//
// Input:           n:      integer to be tested
// Output:                  YES if n is a power of 2
//
// ############################ Public Function ###################################
LOGICAL powerOfTwo(int n)
{
  int   old_n, two_to_n;
  int   _powerOfTwo = -1;

  old_n = n;
  while(1)
  {
    if(n==0)
      break;
    n >>= 1;            // Divide by 2
    _powerOfTwo++;
  }
  if(_powerOfTwo < 1)
    return NO;
  two_to_n  = 1<<_powerOfTwo;
  if(old_n == two_to_n)
    return YES;
  return NO;
}
