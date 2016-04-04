/************************************************************************
 *                                                                      *
 * This file contains some general purpose matrix routines.             *
 *                                                                      *
 * File:matrix_routines.cc                                              *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 07/06/04  - Extracted from DoubleMatrix.                         *
 ************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#if defined(WIN32)
#include <C_Libraries/ftoc.h>
#include <AlgebraRoutines/algebra_routines.h>
#include <AlgebraRoutines/misc.h>
#else
#include "ftoc.h"
#include "algebra_routines.h"
#include "misc.h"
#endif
#include "matrix_routines.h"

#define MAX(a,b)  ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))
#ifndef YES
#define YES     1
#define NO      0
#endif

// ############################ Private Function ###################################
// luDecomp - Calculate the LU decomposition.
//            Taken from: C Tools for Engineers and Scientists.
//
// Input:               inputMatrix             matrix to perform lu matrix
//                      size                    number of rows or columns
//
// Output:              luMatrix                output LU decomposition
//                      pivot                   pivot vector
//
// Returns:                                     Returns YES on success, NO, otherwise
//
// NOTES:
// 1. This routine must be called before findDeterminant(),  etc.
// ############################ Private Function ###################################
LOGICAL lu_decomp(const double *inputMatrix, int size, double *luMatrix, int *pivot)
{
  int           number_columns, number_rows;
  int           i,j,n,k,l,kp1,nm1;
  double        *a, t;

  number_rows                           = number_columns        = n     = size;
  a                                     = luMatrix;
  for(i=0; i<size*size; i++)
    luMatrix[i]                         = inputMatrix[i];
  nm1                                   = n-1;

  if (nm1>=1)                                                   // Non-trivial problem
  {
    DOFOR(k,nm1)
    {
      kp1=k+1;
//
// Partial pivoting ROW exchanges-search over column
// in FORTRAN, the increment would be 1 not n in izamax call
//
      pivot [k]                         = l = isamax((n-k),&(a INDEX(k,k)),number_columns)+k;
      if(a INDEX(l,k) != 0.)
      {                                                         // nonsingular pivot found
        if(l!=k)
        {                                                       // interchange needed
          t                             = a INDEX(l,k);
          a INDEX(l,k)                  = a INDEX(k,k);
          a INDEX(k,k)                  = t;
        }
        t                               = -1./a INDEX(k,k);     // scale row
        sscal(nm1-k,t,&(a INDEX(k+1,k)),number_columns);
        DOBYYY(j,kp1,n)
        {
          t                             = a INDEX(l,j);
          if(l!=k)
          {
            a INDEX(l,j)                = a INDEX(k,j);
            a INDEX(k,j)                = t;
          }
          saxpy(nm1-k,t,&(a INDEX(k+1,k)), n, &(a INDEX(k+1,j)), n);
        }
      }
      else /*pivot singular*/
      { 
        fprintf(stderr, "Singular matrix found in LU decomp\n");
        return NO;
      }
    }                                                           // main loop over k
  }
  pivot [nm1]=nm1;
  if (a INDEX(nm1,nm1) ==0.0)
  {
    fprintf(stderr, "Singular matrix found in LU decomp\n");
    return NO;
  }
  return YES;
}

// ############################ Private Function ###################################
// linear_solve - Finds the solution to Ax=b, using the LU decomposition.
//               Taken from: C Tools for Engineers and Scientists.
//
// Input:               luMatrix        LU decomposition of A
//                      pivot           Pivot vector from LU decomposition
//                      size            Size of A matrix
//
// Input/Output:        bVector         b vector on input, x result on output
//
// Returns                      YES on success
//
// NOTES:
// 1. This routine calls luDecomp().
// ############################ Private Function ###################################
LOGICAL linear_solve(const double *luMatrix, const int *pivot, double *bVector, int size)
{

  int           number_columns;
  int           n,k,l,nm1;
  const double  *a;
  double        t;
//
// Change names
//
  a                     = luMatrix;
  n                     = number_columns = size;
  nm1                   = n-1;
//
// solve ly=b first
//
  DOFOR(k,nm1)
  {
    l                   = pivot[k];
    t                   = bVector[l];
    if(l!=k)
    {
      bVector[l]        = bVector[k];
      bVector[k]        = t;
    }
    saxpy( nm1-k,t, &(a INDEX(k+1,k)),n,&(bVector[k+1]),1);
  }

/* solve Ux=y*/
  DOFOR(l,n)
  {
    k                   = nm1-l;
    bVector[k]          = bVector[k]/ a INDEX(k,k);
    t                   = -bVector[k];
    saxpy(k,t,&(a INDEX(0,k)),n,bVector,1);
  }

  return        YES;
}

// ############################ Public Function ###################################
// invert_using_LU - Inverts the matrix, using the LU decomposition.
//          Taken from: C Tools for Engineers and Scientists.
//
// Input:       luMatrix        Previously calculated LU decomposition
//              pivot           Pivot vector from LU decomposition
//              size            number of rows or columns
//
// Output:      inverseMatrix   inverse
//
// Returns                      YES on success
//
// NOTES:
// ############################ Public Function ###################################
LOGICAL invert_using_LU(const double *luMatrix, const int *pivot, int size, double *inverseMatrix)
{
  int           number_columns, n;
  int           i,j,k,l,kb,kp1,nm1;
  double        *work, *a, t;

  j                             = size*size;
  for(i=0; i<j; i++)
    inverseMatrix[i]            = luMatrix[i];
  a                             = inverseMatrix;
  number_columns                = n = size;
  nm1                           = n-1;
  work                          = new double[number_columns];
  DOFOR(k,n)
  {
        a INDEX(k,k)=t=1./ a INDEX(k,k);
        t = -t;
        sscal(k,t,&(a INDEX(0,k)),number_columns);
        kp1 =k+1;
        if (nm1>=kp1)
        {
            DOBYYY(j,kp1,n)
            {
                t =a INDEX(k,j);
                a INDEX(k,j) = 0.0;
                saxpy(k+1,t,&(a INDEX(0,k)),number_columns,&(a INDEX(0,j)),
                      number_columns);
            }
        }

  }
/*inv(u)*inv(l)*/
  if (nm1>=1)
  {
        DOFOR(kb,nm1)
        {
            k   =nm1-kb-1;
            kp1 =k+1;
            DOBYYY(i,kp1,n)
            {
                work[i]         =a INDEX(i,k);
                a INDEX(i,k)    =0.0;
            }
            DOBYYY(j,kp1,n)
            {
                t=work [j];
                saxpy(n,t,&(a INDEX(0,j)),number_columns,&(a INDEX(0,k))
                           ,number_columns);
            }
            l   = pivot [k];
            if(l!=k)
            {
              sswap(n,&(a INDEX(0,k)),number_columns,&(a INDEX(0,l))
                          ,number_columns);
            }
        }
  }
//
// delete the temp arrays:
//
  delete [] work;

  return YES;
}

// ############################ Public Function ###################################
// toeplitz_solve - Finds the solution to Ax=b, using the Toeplitz matrix algorithm
//                 as given in "Digital Signal Processing," by Roberts and Mullis.
//
// Input:       r:                      The first row of the Toeplitz matrix (may be modified!!)
//              b:                      The result of Ax=b
//              size:                   The dimension of the Toeplitz matrix
//
// Output:      x:                      The resulting vector
//              conda:                  Estimate of the 2-norm of the matrix
//              flag:                   0 = normal, 1 = illegal n, 2 = singular matrix found
//
// NOTES:
// 1. This function taken from the Netlib function dsytcl()
// ############################ Public Function ###################################
void toeplitz_solve(double *r, double *b, int size, double *x, double *conda, int *flag)
{
  int           i,  k, ip, n;
  double        *y, *w;
  double        alpha, beta, c, eta, gamma, t, yn2;
  double        s, smin;
  LOGICAL       compute_error;
//
// Initialization:
//
  k                     = 0;
  n                     = size;
  *flag                 = 0;
  if(n <= 0)
  {
    *flag               = 1;
    return;
  }
//
// Allocate arrays:
//
  x                     = new double[n];
  y                     = new double[n];
  w                     = new double[n];
  compute_error         = NO;
  while (YES)                           // break out of this loop on error
  {
    if(fabs(r[0]) == 0.)
    {
      compute_error     = YES;
      break;
    }
//
// Compute the first elements of x and y:
//
    x[0]                = b[0]/r[0];
    y[0]                = -r[1]/r[0];
    gamma               = r[0] + r[1]*y[0];
    smin                = fabs(r[0]);
    s                   = smin;
    k                   = 0;            // # of elements to process
    while((n>1) && (k<(n-1)))
    {
//
// Estimate the smallest singular value of leading submatrix
//
      ip                = isamax(k+1, y, 1);
      yn2               = y[ip]*y[ip];
      s                 = fabs(gamma)/MAX(1., (fabs(yn2)));
      if(s < smin)
        smin            = s;
//
// Perform a classical Levinson step and update the algorithm:
//
      beta              = b[k+1] - sdot(k+1, &r[1], 1, &x[0], -1);      // fixed
      c                 = -r[(k+2)%n] - sdot(k+1, &r[1], 1, &y[0], -1);
      if(fabs(gamma) == 0.)
        break;
      alpha             = beta/gamma;
      eta               = c/gamma;
      for(i=0; i<=k; i++)
      {
        x[i]            += alpha*y[k-i];                // Should it be [k-i+1]?
        w[i]            = y[i] + eta*y[k-i];            // same comment
      }
      for(i=0; i<=k; i++)
        y[i]            = w[i];
      k++;
      x[k]              = alpha;
      y[k]              = eta;
      gamma             -= c*eta;
    }
    t                   = sasum(n, &r[0], 1);
    *conda              = fabs(t)/smin;
    break;
  }
  if(compute_error)
    *flag               = 2;

  delete [] y;
  delete [] w;
  return;
}

#define SORT_DESCENDING 1
// ############################ Public Function ###################################
// svd - Finds the singular value decomposition of an input matrix.
//
// Input:       a:                      The matrix to decompose
//              aRows:                  number of rows in a
//              aCols:                  number of cols in a
//              uFlag                   = (0: do not compute u); (1: compute u[m,n]); (-1: compute u[m,m])
//              vFlag                   = (0: do not compute v); (1: compute v)
//
// Outputs:     singularValues:         Singular values from largest to smallest
//              u:                      U matrix (see uFlag)
//              v:                      V matrix (see vFlag)
//
// Returns:                             error flag
//
// NOTES:
// 1. This function taken from C Tools for Engineers and Scientists.
// ############################ Public Function ###################################
int     svd(const double *a, int aRows, int aCols, double *singularValues, int uFlag, double *u, int vFlag, double *v)
{
  int           m, n, matu, matv, ierr, coln, number_columns;
  int           i,j,k,mn,l,m1,n1,top,colu;
  double        dum,f,g,h,s,y,scale,anorm;
  double        *rv1;
//
// Initialize:
//
  m     = aRows;
  n     = coln  = number_columns        = aCols;
  matu  = uFlag;
  matv  = vFlag;
  ierr=0;
  colu=n;
  rv1   = new double[n];
  if (matu<0)
    colu=m;

  m1=m-1;
  n1=n-1;
  DOFOR(i,m)
  {
    DOFOR(j,n)
      u IU(i,j)= a INDEX(i,j);
  }
/* reduction to bidiagonal form via Householder xfms*/
g       = 0.;
scale   = 0.;
anorm   = 0.;
l       = 0;
DOFOR(i,n)
        {
        l=i+1;
        rv1[i]=g*scale;
        g=0.;
        s=0.;
        scale=0.;
        if(i<=m1)/* in fortran, <= but c i runs 0-m1*/
                {
                DFFOR(k,i,m)scale+= fabs(u IU(k,i));    
                if(scale!=0.)
                        {
                        DFFOR(k,i,m)
                                {
                                dum=u IU(k,i)/scale;
                                s+=dum*dum;     
                                u IU(k,i)=dum;
                                }
                        f= u IU(i,i);
                        g= -dsign(sqrt(s),f);
                        h=1./(f*g-s);/* 1/h of others- * faster than / */
                        u IU(i,i)=f-g;
                        if (i!=n1)
                                {
                                DFFOR(j,l,n)
                                        {
                                        s=0.;
                                        DFFOR(k,i,m)s+=
                                         u IU(k,i)*u IU(k,j);
                                        f=s*h;
                                        DFFOR(k,i,m)
                                                u IU(k,j)+=f*u IU(k,i);
                                        }
                                }/*i!=n*/
                        DFFOR(k,i,m)u IU(k,i)*=scale;
                        }/*scale!=0*/
                }/* i<=m*/
singularValues[i]=scale*g;
g=0.;
s=0.;
scale=0.;
if(i<=m1 && i!=n1)
        {
        DFFOR(k,l,n)scale+= fabs(u IU(i,k));
        if(scale!=0.)
                {
                DFFOR(k,l,n)
                        {
                        dum= u IU(i,k)/scale;
                        s+=dum*dum;
                        u IU(i,k)=dum;
                        }
                f=u IU(i,l);
                g=-dsign(sqrt(s),f);
                h=1./(f*g-s);
                u IU(i,l)=f-g;
                DFFOR(k,l,n)rv1[k]=u IU(i,k)*h;
                if(i!=m1)
                        {
                        DFFOR(j,l,m)
                                {
                                s=0.;
                                DFFOR(k,l,n)
                                        s+=u IU(j,k)*u IU(i,k);
                                DFFOR(k,l,n)u IU(j,k)+=s*rv1[k];

                                }

                        }
                DFFOR(k,l,n)u IU(i,k)*=scale;
                }
        }
y=fabs(singularValues[i])+fabs(rv1[i]) ;
anorm= MAX( anorm , y);
        }
/*accum. rhs*/
/*printf(" acc rhs\n");*/
if(matv)
        {
        for(i=n-1;i>=0;i--)
                {
                if(i!=n1)
                        {
                        if (g!=0.)
                                {/*here FMM do division by g,u indiv*/
                                /*h=1./(g*u IU(i,l));*/
                                DFFOR(j,l,n)v INDEX(j,i)
                                        =(u IU(i,j)/u IU(i,l) )/g;
                                DFFOR(j,l,n)
                                        {s=0.;
                                        DFFOR(k,l,n)s+=
                                                u IU(i,k)*v INDEX(k,j);
                                        DFFOR(k,l,n)v INDEX(k,j)+=
                                                s* v INDEX(k,i);
                                        }
                                }
                        DFFOR(j,l,n)
                                {
                                v INDEX(i,j)=0.;
                                v INDEX(j,i)=0.;
                                }
                        }
                v INDEX(i,i)=1.;
                g=rv1[i];
                l=i;
                }
        }
/*accum. lhs*/
/*printf(" accum lhs\n");*/
if(matu)
        {
        if(matu==-1)
                {
                DFFOR(i,n,m)
                        {
                        DFFOR(j,n,m) u IU(i,j)=0.;
                        u IU(i,i)=1.;
                        }
                }       
        mn=n;
        if(m<n) mn=m;
        top=n;
        if(matu==-1)top=m;
        for(i=mn-1;i>=0;i--)
                {
                l=i+1;
                g=singularValues[i];
                if (i!=n1)
                        {
                        DFFOR(j,l,top)u IU(i,j)=0.;
                        }
                if(g!=0.)
                        {
                        /*h= 1./(g* u IU(i,i));*/
                        if(i!=(mn-1))
                        {
                        DFFOR(j,l,top)
                                {
                                s=0.;
                                DFFOR(k,l,m)
                                        s += u IU(k,i) * u IU(k,j);
                                f=(s/u IU(i,i))/g;
                                DFFOR(k,i,m)u IU(k,j)+=f*u IU(k,i);
                                }/* DFFOR j */
                        }/* i!=mn-1*/
                                DFFOR(j,i,m)u IU(j,i)/=g;
                        
                        }
                else/*g=0*/
                        {
                        DFFOR(j,i,m)u IU(j,i)=0.;
                        }               
                u IU(i,i)+=1.;
                }
        }
  ierr  = diagon(u,v,rv1,singularValues,anorm,coln,colu,n,m);
//
// Sort the singular values:
//
  sort(n, singularValues, SORT_DESCENDING);  delete [] rv1;
  return ierr;
}

#define ITMAX   30
// ############################ Public Function ###################################
// diagon - diagonalize in binary form.
//
// Input:       u:                      U matrix
//              v:                      V matrix
//              rv1:                    work vector
//              w:                      singular values
//              anorm:                  Norm of A
//
// Returns:                             error flag
//
// NOTES:
// 1. This function taken from C Tools for Engineers and Scientists.
// ############################ Public Function ###################################
int diagon(double *u, double *v, double *rv1, double *w, double anorm, int coln, int colu, int n, int m)
{
  int           i,j,k,l,l1,k1,i1,matu,matv,its, ierr;
  int           number_columns;
  double        c,s,h,f,g,x,z,y;

  number_columns        = coln;
  matu                  = 1;
  matv                  = 1;
  ierr                  = 0;
  l1                    = 0;
/* diagonalize binary form*/
  for(k=n-1;k>=0;k--)
        {
        k1=k-1;
        its=0;
testsplit:
        for(l=k;l>=0;l--)
                {
                l1=l-1;
                if ( fabs(rv1[l])+anorm==anorm)goto conv;
                if(fabs(w[l1])+anorm==anorm)break;
                }
        c=0.;
        s=1.;
        DFFOR(i,l,k+1) /* want i to range from l to k inclusive */
                {
                f=s*rv1[i];
                rv1[i]*=c;      
        if(fabs(f)+anorm!= anorm )
                        {
                        g=w[i];
                        h=sqrt(f*f+g*g);
                        w[i]=h;
                        c=g/h;
                        s=-f/h;
                        if(matu)
                                {
                                DOFOR(j,m)
                                        {
                                        y=u IU(j,l1) ;
                                        z=u IU(j,i) ;
                                        u IU(j,l1)=y*c+z*s;
                                        u IU(j,i)=-y*s+z*c;
                                        }
                                }/*matu*/
                        }/* if */
                }/* DFFOR i,l,k*/
                conv:;
                z=w[k];
                if(l!=k)
                        {
                        if (its>=ITMAX)
                                {
                                ierr=k;
                                return ierr;
                                }
                        its++;
                        x=w[l];
                        y=w[k1];
                        g=rv1[k1];
                        h=rv1[k];
                        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.*h*y);
                        g=sqrt(f*f+1.);
                        f=((x-z)*(x+z) +h*(y/(f+dsign(g,f))-h))/x;
                        /* do another QR xfm*/
                        c=1.;
                        s=1.;
                        DFFOR(i1,l,k1+1)/* i1 ranges l to k1 inclusive*/
                                {
                                i=i1+1;
                                g=rv1[i];
                                y=w[i];
                                h=s*g;
                                g=c*g;
                                z=sqrt(f*f+h*h);
                                rv1[i1]=z;
                                c=f/z;
                                s=h/z;
                                f=x*c+g*s;
                                g=-x*s+g*c;
                                h=y*s;
                                y=y*c;
                                if(matv)
                                        {
                                        DOFOR(j,n)
                                                {
                                                x= v INDEX(j,i1);
                                                z= v INDEX(j,i);
                                                v INDEX(j,i1)=
                                                   x*c+z*s;
                                                v INDEX(j,i)=
                                                   -x*s+z*c;
                                                }
                                        }
                                z=sqrt(h*h+f*f);
                                w[i1]=z;
                                if(z!=0.)
                                        {
                                        c=f/z;
                                        s=h/z;
                                        }
                                f=c*g+s*y;
                                x=-s*g+c*y;
                                if(matu)
                                        {
                                        DOFOR(j,m)
                                                {
                                                y= u IU(j,i1);
                                                z= u IU(j,i);
                                                u IU(j,i1)=
                                                   y*c+z*s;
                                                u IU(j,i)=
                                                   -y*s+z*c;
                                                }
                                        }       
                                }
                        rv1[l]=0.;
                        rv1[k]=f;
                        w[k]=x;
                        goto testsplit; 
                        }/*l!=k*/               
        /*convergence*/
        if(z<0.)
                {/* make w[k] non-neg*/
                w[k]=-z;
                if(matv)
                {
                DOFOR(j,n) v INDEX(j,k)=-v INDEX(j,k);
                }
                }       
        }
return ierr;
}


// ############################ Public Function ###################################
// dsign - Sign function.
//
//
// NOTES:
// 1. This function taken from C Tools for Engineers and Scientists.
// ############################ Public Function ###################################
double dsign(double to, double from)
{
  double x;
  x=0.;
  if(from>0.)x=1.;
  if(from<0.)x=-1.;
  return( x*fabs(to));
}

















// ############################ Private Function ###################################
// luDecomp - Calculate the LU decomposition of a Complex matrix
//            Taken from: C Tools for Engineers and Scientists.
//
// Input:               inputMatrix             matrix to perform lu matrix
//                      size                    number of rows or columns
//
// Output:              luMatrix                output LU decomposition
//                      pivot                   pivot vector
//
// Returns:                                     Returns YES on success, NO, otherwise
//
// NOTES:
// 1. This routine must be called before findDeterminant(),  etc.
// ############################ Private Function ###################################
LOGICAL lu_decomp(const Complex *inputMatrix, int size, Complex *luMatrix, int *pivot)
{
  int           number_columns, number_rows;
  int           i,j,n,k,l,kp1,nm1;
  Complex       *a, t;

  number_rows                           = number_columns        = n     = size;
  a                                     = luMatrix;
  for(i=0; i<size*size; i++)
    luMatrix[i]                         = inputMatrix[i];
  nm1                                   = n-1;

  if (nm1>=1)                                                   // Non-trivial problem
  {
    DOFOR(k,nm1)
    {
      kp1=k+1;
//
// Partial pivoting ROW exchanges-search over column
// in FORTRAN, the increment would be 1 not n in izamax call
//
      pivot [k]                         = l = izamax((n-k),&(a INDEX(k,k)),number_columns)+k;
      if(abs(a INDEX(l,k)) != 0.)
      {                                                         // nonsingular pivot found
        if(l!=k)
        {                                                       // interchange needed
          t                             = a INDEX(l,k);
          a INDEX(l,k)                  = a INDEX(k,k);
          a INDEX(k,k)                  = t;
        }
        t                               = -1./a INDEX(k,k);     // scale row
        zscal(nm1-k,t,&(a INDEX(k+1,k)),number_columns);
        DOBYYY(j,kp1,n)
        {
          t                             = a INDEX(l,j);
          if(l!=k)
          {
            a INDEX(l,j)                = a INDEX(k,j);
            a INDEX(k,j)                = t;
          }
          zaxpy(nm1-k,t,&(a INDEX(k+1,k)), n, &(a INDEX(k+1,j)), n);
        }
      }
      else /*pivot singular*/
      { 
        fprintf(stderr, "Singular matrix found in LU decomp\n");
        return NO;
      }
    }                                                           // main loop over k
  }
  pivot [nm1]=nm1;
  if (abs(a INDEX(nm1,nm1)) ==0.0)
  {
    fprintf(stderr, "Singular matrix found in LU decomp\n");
    return NO;
  }
  return YES;
}

// ############################ Private Function ###################################
// linear_solve - Finds the solution to Ax=b, using the LU decomposition for Complex matrices
//               Taken from: C Tools for Engineers and Scientists.
//
// Input:               luMatrix        LU decomposition of A
//                      pivot           Pivot vector from LU decomposition
//                      size            Size of A matrix
//
// Input/Output:        bVector         b vector on input, x result on output
//
// Returns                      YES on success
//
// NOTES:
// 1. This routine calls luDecomp().
// ############################ Private Function ###################################
LOGICAL linear_solve(const Complex *luMatrix, const int *pivot, Complex *bVector, int size)
{

  int           number_columns;
  int           n,k,l,nm1;
  const Complex *a;
  Complex       t;
//
// Change names
//
  a                     = luMatrix;
  n                     = number_columns = size;
  nm1                   = n-1;
//
// solve ly=b first
//
  DOFOR(k,nm1)
  {
    l                   = pivot[k];
    t                   = bVector[l];
    if(l!=k)
    {
      bVector[l]        = bVector[k];
      bVector[k]        = t;
    }
    zaxpy( nm1-k,t, &(a INDEX(k+1,k)),n,&(bVector[k+1]),1);
  }

/* solve Ux=y*/
  DOFOR(l,n)
  {
    k                   = nm1-l;
    bVector[k]          = bVector[k]/ a INDEX(k,k);
    t                   = -bVector[k];
    zaxpy(k,t,&(a INDEX(0,k)),n,bVector,1);
  }

  return        YES;
}

// ############################ Public Function ###################################
// invert_using_LU - Inverts the Complex matrix, using the LU decomposition.
//          Taken from: C Tools for Engineers and Scientists.
//
// Input:       luMatrix        Previously calculated LU decomposition
//              pivot           Pivot vector from LU decomposition
//              size            number of rows or columns
//
// Output:      inverseMatrix   inverse
//
// Returns                      YES on success
//
// NOTES:
// ############################ Public Function ###################################
LOGICAL invert_using_LU(const Complex *luMatrix, const int *pivot, int size, Complex *inverseMatrix)
{
  int           number_columns, n;
  int           i,j,k,l,kb,kp1,nm1;
  Complex       *work, *a, t;

  j                             = size*size;
  for(i=0; i<j; i++)
    inverseMatrix[i]            = luMatrix[i];
  a                             = inverseMatrix;
  number_columns                = n = size;
  nm1                           = n-1;
  work                          = new Complex[number_columns];
  DOFOR(k,n)
  {
        a INDEX(k,k)=t=1./ a INDEX(k,k);
        t = -t;
        zscal(k,t,&(a INDEX(0,k)),number_columns);
        kp1 =k+1;
        if (nm1>=kp1)
        {
            DOBYYY(j,kp1,n)
            {
                t =a INDEX(k,j);
                a INDEX(k,j) = 0.0;
                zaxpy(k+1,t,&(a INDEX(0,k)),number_columns,&(a INDEX(0,j)),
                      number_columns);
            }
        }

  }
/*inv(u)*inv(l)*/
  if (nm1>=1)
  {
        DOFOR(kb,nm1)
        {
            k   =nm1-kb-1;
            kp1 =k+1;
            DOBYYY(i,kp1,n)
            {
                work[i]         =a INDEX(i,k);
                a INDEX(i,k)    =0.0;
            }
            DOBYYY(j,kp1,n)
            {
                t=work [j];
                zaxpy(n,t,&(a INDEX(0,j)),number_columns,&(a INDEX(0,k))
                           ,number_columns);
            }
            l   = pivot [k];
            if(l!=k)
            {
              zswap(n,&(a INDEX(0,k)),number_columns,&(a INDEX(0,l))
                          ,number_columns);
            }
        }
  }
//
// delete the temp arrays:
//
  delete [] work;

  return YES;
}

// ############################ Public Function ###################################
// toeplitz_solve - Finds the solution to Ax=b, using the Toeplitz matrix algorithm
//                 as given in "Digital Signal Processing," by Roberts and Mullis.
//
// Input:       r:                      The first row of the Toeplitz matrix (may be modified!!)
//              b:                      The result of Ax=b
//              size:                   The dimension of the Toeplitz matrix
//
// Output:      x:                      The resulting vector
//              conda:                  Estimate of the 2-norm of the matrix
//              flag:                   0 = normal, 1 = illegal n, 2 = singular matrix found
//
// NOTES:
// 1. This function taken from the Netlib function dsytcl()
// ############################ Public Function ###################################
void toeplitz_solve(Complex *r, Complex *b, int size, Complex *x, Complex *conda, int *flag)
{
  int           i,  k, ip, n;
  Complex       *y, *w;
  Complex       alpha, beta, c, eta, gamma, t;
  double        s, smin, yn2;
  LOGICAL       compute_error;
//
// Initialization:
//
  k                     = 0;
  n                     = size;
  *flag                 = 0;
  if(n <= 0)
  {
    *flag               = 1;
    return;
  }
//
// Allocate arrays:
//
  x                     = new Complex[n];
  y                     = new Complex[n];
  w                     = new Complex[n];
  compute_error         = NO;
  while (YES)                           // break out of this loop on error
  {
    if(abs(r[0]) == 0.)
    {
      compute_error     = YES;
      break;
    }
//
// Compute the first elements of x and y:
//
    x[0]                = b[0]/r[0];
    y[0]                = -r[1]/r[0];
    gamma               = r[0] + r[1]*y[0];
    smin                = abs(r[0]);
    s                   = smin;
    k                   = 0;            // # of elements to process
    while((n>1) && (k<(n-1)))
    {
//
// Estimate the smallest singular value of leading submatrix
//
      ip                = izamax(k+1, y, 1);
      yn2               = norm(y[ip]);
      s                 = abs(gamma)/MAX(1., yn2);
      if(s < smin)
        smin            = s;
//
// Perform a classical Levinson step and update the algorithm:
//
      beta              = b[k+1] - zdot(k+1, &r[1], 1, &x[0], -1);      // fixed
      c                 = -r[(k+2)%n] - zdot(k+1, &r[1], 1, &y[0], -1);
      if(abs(gamma) == 0.)
        break;
      alpha             = beta/gamma;
      eta               = c/gamma;
      for(i=0; i<=k; i++)
      {
        x[i]            += alpha*y[k-i];                // Should it be [k-i+1]?
        w[i]            = y[i] + eta*y[k-i];            // same comment
      }
      for(i=0; i<=k; i++)
        y[i]            = w[i];
      k++;
      x[k]              = alpha;
      y[k]              = eta;
      gamma             -= c*eta;
    }
    t                   = zasum(n, &r[0], 1);
    *conda              = abs(t)/smin;
    break;
  }
  if(compute_error)
    *flag               = 2;

  delete [] y;
  delete [] w;
  return;
}

#define GSL_DBL_EPSILON        2.2204460492503131e-16
// ############################ Public Function ###################################
// svd - Finds the singular value decomposition of an input Complex matrix.
//
// Input:       aMatrix:                The matrix to decompose
//              aRows:                  number of rows in a
//              aCols:                  number of cols in a
//
// Outputs:     singularValues:         Singular values from largest to smallest
//              uMatrix:                U matrix (see uFlag)
//              qMatrix:                Q matrix (see vFlag)
//
// Returns:                             error flag
//
// NOTES:
// 1. This function taken from GSL.
// 2. Not yet working for complex matrices.
// where U is a an M x N matrix (U^T U = I), D is a diagonal N x N
// matrix and Q is an N x N orthogonal matrix (Q^T Q = I)
//
// U is stored in the original matrix A, which has the same size
//
// Q is stored as a separate matrix (not Q^T). You must take the
// transpose to form the product above.
//
// The diagonal matrix D is stored in the vector S,  D_ii = S_i
//
// The one-sided Jacobi algorithm is used (see Algorithm 4.1 in the
// following paper). 
//
// See James Demmel, Kresimir Veselic, "Jacobi's Method is more
// accurate than QR", Lapack Working Note 15 (LAWN15), October 1989.
// Available from netlib.  
// ############################ Public Function ###################################
int     svd(const Complex *aMatrix, int aRows, int aCols, double *singularValues, Complex *uMatrix, Complex *qMatrix)
{
  int           i, j, k, m, n;
  int           number_columns;
  int           rtn_value;
  Complex       a_ik, a_ij, q_ik, q_ij, new_ij, new_ik;

  m                     = aRows;
  n                     = aCols;
  number_columns        = n;
//
// Copy A to U
//
  for(i=0; i<aRows; i++)
  {
    for(j=0; j<aCols; j++)
    {
      uMatrix INDEX(i,j)        = aMatrix INDEX(i,j);
    }
  }

        /* Initialize the rotation counter and the sweep counter. */
      int count = 1;
      int sweep = 0;
      int sweepmax = n;

      double tolerance = 10 * GSL_DBL_EPSILON ;

      /* Always do at least 12 sweeps. */
      sweepmax = MAX (sweepmax, 12);

      /* Set Q to the identity matrix. */
//
//       gsl_matrix_set_identity (Q);
//
      k                 = 0;
      for(i=0; i<n; i++)
      {
        for(j=0; j<n; j++)
        {
          if (i==j)
            qMatrix[k++]        = Complex(1., 0.);
          else
            qMatrix[k++]        = Complex(0., 0.);
        }
      }

      /* Orthogonalize A by plane rotations. */

      while (count > 0 && sweep <= sweepmax)
        {
          /* Initialize rotation counter. */
          count = n * (n - 1) / 2;

          for (j = 0; j < n - 1; j++)
            {
              for (k = j + 1; k < n; k++)
                {
                  double a = 0.0;
                  double b = 0.0;
                  Complex p;
                  double q = 0.0;
                  double r = 0.0;
                  Complex cosine, sine;
                  double v;
//
//  gsl_vector cj = gsl_matrix_column (A, j);
//  gsl_vector ck = gsl_matrix_column (A, k);
//  gsl_blas_ddot (&cj, &ck, &p);
//
//  a = gsl_blas_dnrm2 (&cj);
//  b = gsl_blas_dnrm2 (&ck);
//  q = a*a;
//  r = b*b;
//
                  p     = complex_dot_columns(uMatrix, aRows, aCols, j, k);
                  a     = complex_norm_column(uMatrix, aRows, aCols, j);
                  b     = complex_norm_column(uMatrix, aRows, aCols, k);
                  q     = a ;
                  r     = b ;
                  a     = sqrt(a);
                  b     = sqrt(b);
                  /* NOTE: this could be handled better by scaling
                   * the calculation of the inner products above.
                   * But I'm too lazy. This will have to do. [GJ]
                   */

                  /* This is an adhoc method of testing for a "zero"
                     singular value. We consider it to be zero if it
                     is sufficiently small compared with the currently
                     leading column. Note that b <= a is guaranteed by
                     the sweeping algorithm. BJG */
                  
                  if (b <= tolerance * a) 
                    {
                      /* probably |b| = 0 */
                      count--;
                      continue;
                    }

                  if (abs(p) <= tolerance * a * b)
                    {
                      /* columns j,k orthogonal
                       * note that p*p/(q*r) is automatically <= 1.0
                       */
                      count--;
                      continue;
                    }

                  /* calculate rotation angles */
                  if (q < r)
                    {
                      cosine = 0.0;
                      sine = 1.0;
                    }
                  else
                    {
                      q -= r;
//
//  v = gsl_hypot (2.0 * p, q);
//
                      v      = gsl_hypot(2.*abs(p), q);
                      cosine = sqrt ((v + q) / (2.0 * v));
                      sine = p / (v * cosine);
                    }

                  /* apply rotation to A */
                  for (i = 0; i < m; i++)
                    {
//
//  const REAL Aik = gsl_matrix_get (A, i, k);
//  const REAL Aij = gsl_matrix_get (A, i, j);
//  gsl_matrix_set (A, i, j, Aij * cosine + Aik * sine);
//  gsl_matrix_set (A, i, k, -Aij * sine + Aik * cosine);
//
                    a_ik                = uMatrix INDEX(i, k);
                    a_ij                = uMatrix INDEX(i, j);
                    new_ij              = a_ij*cosine + a_ik*sine;
                    new_ik              = -a_ij*sine + a_ik*cosine;
                    uMatrix INDEX(i, k) = new_ik;
                    uMatrix INDEX(i, j) = new_ij;
                    }

                  /* apply rotation to Q */
                  for (i = 0; i < n; i++)
                    {
//
//                    const REAL Qij = gsl_matrix_get (Q, i, j);
//                    const REAL Qik = gsl_matrix_get (Q, i, k);
//                    gsl_matrix_set (Q, i, j, Qij * cosine + Qik * sine);
//                    gsl_matrix_set (Q, i, k, -Qij * sine + Qik * cosine);
//
                    q_ik                = qMatrix INDEX(i, k);
                    q_ij                = qMatrix INDEX(i, j);
                    new_ij              = q_ij*cosine + q_ik*sine;
                    new_ik              = -q_ij*sine + q_ik*cosine;
                    qMatrix INDEX(i, k) = new_ik;
                    qMatrix INDEX(i, j) = new_ij;
                    }
                }
            }

          /* Sweep completed. */
          sweep++;
        }

      /* 
       * Orthogonalization complete. Compute singular values.
       */

      for (j = 0; j < n; j++)
        {
          /* Calculate singular values. */
//  gsl_vector column = gsl_matrix_column (A, j);
//  double norm = gsl_blas_dnrm2 (&column);
//  gsl_vector_set (S, j, norm);
//
          singularValues[j]     = sqrt(complex_norm_column(uMatrix, aRows, aCols, j));

          /* Normalize vectors. */
          if(singularValues[j] > 0.)
          {
            for (i = 0; i < m; i++)
            {
//  const REAL Aij = gsl_matrix_get (A, i, j);
//  const REAL Sj = gsl_vector_get (S, j);
//  gsl_matrix_set (A, i, j, Aij / Sj);
              a_ij                      = uMatrix INDEX(i, j);
              new_ij                    = a_ij/singularValues[j];
              uMatrix INDEX(i, j)       = new_ij;
            }
          }
        }

      if (count > 0)
        {
          /* reached sweep limit */
          rtn_value     = -1;
        }
      else
      {
        rtn_value       = 1;
      }

      return rtn_value;
}


double gsl_hypot (const double x, const double y)
{
  double xabs = fabs(x) ;
  double yabs = fabs(y) ;
  double min, max;

  if (xabs < yabs) {
    min = xabs ;
    max = yabs ;
  } else {
    min = yabs ;
    max = xabs ;
  }

  if (min == 0) 
    {
      return max ;
    }

  {
    double u = min / max ;
    return max * sqrt (1 + u * u) ;
  }
}


// ############################ Public Function ###################################
// complex_dot_columns - Forms the dot product between two columns of matrix.
//
// Input:       aMatrix:                The matrix to decompose
//              aRows:                  number of rows in a
//              aCols:                  number of cols in a
//              column1:                1st column
//              column2:                2nd column
//
// Outputs:
//
// Returns:                             dot product
//
// NOTES:
// ############################ Public Function ###################################
Complex complex_dot_columns(const Complex *aMatrix, int aRows, int aCols, int column1, int column2)
{
  int           i, number_columns;
  Complex       dot_product;

  number_columns        = aCols;
  dot_product           = Complex(0., 0.);
  for(i=0; i<aRows; i++)
  {
    dot_product         += conj(aMatrix INDEX(i, column1)) * aMatrix INDEX(i, column2);
  }
  return dot_product;
}


// ############################ Public Function ###################################
// complex_norm_column - Finds the norm of a Complex matrix column.
//
// Input:       aMatrix:                The matrix to decompose
//              aRows:                  number of rows in a
//              aCols:                  number of cols in a
//              column:                 column to find norm
//
// Outputs:
//
// Returns:                             norm
//
// NOTES:
// ############################ Public Function ###################################
double complex_norm_column(const Complex *aMatrix, int aRows, int aCols, int column)
{
  int           i, number_columns;
  double        col_norm;

  number_columns        = aCols;
  col_norm              = 0.;
  for(i=0; i<aRows; i++)
  {
    col_norm            += norm(aMatrix INDEX(i, column) );
  }
  return col_norm;
}
