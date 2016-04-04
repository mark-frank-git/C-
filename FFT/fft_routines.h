#ifndef _FFT_ROUTINE_H
#define _FFT_ROUTINE_H  1
/************************************************************************
 *                                                                      *
 *  This file contains subroutines for performing Fourier transfomrs.   *
 *                                                                      *
 *  Author: Mark Frank                                                  *
 *                                                                      *
 *  File: fft_routines.h                                                *
 *                                                                      *
 *  Revisions History:                                                  * 
 *   08/27/95  - started.                                               *
 *   09/15/95  - c.real() -> real(c).                                   *
 *   08/25/99  - Added complex input CZT.                               *
 ************************************************************************/

#define INVERSE_FFT             -1              // See _operationType below
#define FORWARD_FFT             1


#ifndef LOGICAL
#define LOGICAL                 char
#endif
#ifndef YES
#define YES                     1
#define NO                      0
#endif


void    complex_fft(Complex *x, int m, int operation);
void    fft_real(Complex *x, int m);
void    complex_dft(Complex *data_in, Complex *data_out, int number_pts,
                    int operation);
void    dct(float *data_in, float *data_out, int number_pts, int number_output,
            int dct_type);
void    real_czt(const float *input_data, int number_input_pts,
                 float omega_start, float omega_end,
                    Complex *output_data, int number_output_pts);
void    complex_czt(const float *input_real, const float *input_imag,
                 int number_input_pts, float omega_start,
                 float omega_end, Complex *output_data,
                 int number_output_pts);
int     log2(int x);
LOGICAL powerOfTwo(int n);

#endif