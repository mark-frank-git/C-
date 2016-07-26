/********************************************************************************
 *                                                                              *
 * This subclass of object implements a general probability distribution        *
 * function.                                                                    *
 *                                                                              *
 * File: /User/frank/Objc_Classes/Probability/Probability.h                     *
 *                                                                              *
 ********************************************************************************/

#include "Probability.h"

#if defined(WIN32)
#include <GNU/Random.h>
#include <GNU/MLCG.h>
#include <GNU/Uniform.h>
#else
#include "Random.h"
#include "MLCG.h"
#include "Uniform.h"
#endif

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ROUND(a) ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )

#define MIN_RHO         -0.99999999
#define MAX_RHO         0.99999999
#define MIN_VARIANCE    1.e-20
#define SQRT2           1.414213562             /* sqrt(2)      */
#define SQRT3           1.732050808             /* sqrt(3)      */
#define SQRT2PI         2.506628275             /* sqrt(2*PI)   */
#define YES             1
#define NO              0

#ifndef NO_EXTERN_C                 /* C libraries */
extern "C"  {
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifndef NO_EXTERN_C
}
#endif

// ############################# Private Function ###############################
// initMeanVariance --Initialize means and variances
//
// Input:                   None
//          
// Output:                  None
//
// Notes:
//  1. The instance variables, _digitalPoles, _digitalZeros are modified.
// ############################# Private Function ###############################
void Probability::initMeanVariance()
{
  int i;

  delete [] mean;
  delete [] variance;
  delete [] stdDev;
  
  mean          = new double[dimensions];
  variance      = new double[dimensions];
  stdDev        = new double[dimensions];
  for(i=0; i<dimensions; i++)
  {
    mean[i]     = 0.;
    variance[i] = stdDev[i] = 1.;
  }
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the Probability class.
//
// Input:           type:           window type
//
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################

Probability::Probability ()
{
// 
// Initialize instance variables:
//
  mean          = variance = stdDev = NULL;
  distribution  = NORMAL;
  probability   = 1.;
  rho           = 0.;
  doubleOutputs = NULL;
  setDimensions(2);
  initMeanVariance();
  
  linearCongruential            = new MLCG(SEED1, SEED2);
  uniformGenerator              = new Uniform(0.0, 1.0, linearCongruential);
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the Probability class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
Probability::~Probability ()
{
  delete linearCongruential;
  delete uniformGenerator;
  delete [] doubleOutputs;
  delete [] mean;
  delete [] variance;
  delete [] stdDev;
  return;
}

// ############################ Public Function ###################################
// initGenerator - Initializes the noise generator to a defined state.
// Input:               None
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void Probability::initGenerator()
{
  if(linearCongruential != NULL)
    linearCongruential->reset();
  return;
}

// ############################# Public Function ###############################
// setMean -- sets a new set of means
//
// Input:       newMean:        Array of means
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void Probability::setMean( double *newMean)
{
  int i;
  for(i=0; i<dimensions; i++)
  {
    mean[i] = newMean[i];
  }
}

// ############################# Public Function ###############################
// setVariance -- sets a new set of variances
//
// Input:       newVariance:    Array of variances
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void Probability::setVariance( double *newVariance)
{
  int i;
  for(i=0; i<dimensions; i++)
  {
    variance[i] = MAX(MIN_VARIANCE, newVariance[i]);
    stdDev[i]   = sqrt(variance[i]);
  }
}

// ############################# Public Function ###############################
// setCorrelationCoefficient -- sets a correlation coefficient for 2D distributions
//
// Input:       newRho:         New correlation coefficient
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void Probability::setCorrelationCoefficient(double newRho)
{
  rho = MAX(MIN_RHO,newRho);
  rho = MIN(MAX_RHO,rho);
}

// ############################# Public Function ###############################
// setProbDistribution -- sets a new type of probability distribution
//
// Input:       newDistribution:        New probability distribution
//          
// Output:                              None
//
// Notes:
// ############################# Public Function ###############################
void Probability::setProbDistribution(int newDistribution)
{
  distribution = newDistribution;
}


// ############################# Public Function ###############################
// setProbability -- sets a new probability of occurrence
//
// Input:       newProbability:         Probability of occurring
//          
// Output:                              None
//
// Notes:
// ############################# Public Function ###############################
void Probability::setProbability(double newProbability)
{
  probability = newProbability;
}

// ############################# Public Function ###############################
// setDimensions -- sets a new number of dimensions for sample outputs
//
// Input:       newDimensions:          New number of dimensions
//          
// Output:                              None
//
// Notes:
// ############################# Public Function ###############################
void Probability::setDimensions(int newDimensions)
{
  delete [] doubleOutputs;
  dimensions    = newDimensions;
  doubleOutputs = new double[dimensions];
  initMeanVariance();
}

// ############################# Public Function ###############################
// setNewSeeds -- sets a new set of random seeds for random object
//
// Input:       seed1->seed2:           New set of seeds
//          
// Output:                              None
//
// Notes:
// ############################# Public Function ###############################
void Probability::setNewSeeds(long seed1, long seed2)
{
  linearCongruential->seed1(seed1);
  linearCongruential->seed2(seed2);
  return;
}


/*##############################*
 * These methods get parameters *
 *##############################*/
double * Probability::getMean()
{
  return mean;
}

double * Probability::getVariance()
{
  return variance;
}

int Probability::getProbDistribution()
{
  return distribution;
}

double  Probability::getProbability()
{
  return probability;
}

int Probability::getDimensions()
{
  return dimensions;
}

// ############################# Public Function ###############################
// newSample -- Returns a set of samples from our distribution
//
// Input:                               None
//          
// Output:                              Array of samples, array size = dimensions
//
// Notes:
// ############################# Public Function ###############################
double * Probability::newSample()
{
  int    i, j, k;
  double *sample, range;
  double u1,s,ln_s;
  double v1,v1_sq,v2_sq, alpha, a, tr;
  double std_devx, std_devy, term, mean_y, sum;
  static double v2, sqrt_lns, temp;
  static int odd = 1;
  static double markov_old = 0.;
  
  sample        = doubleOutputs;
  switch(distribution)
  {
    case UNIFORM:
      for(i=0; i<dimensions; i++)
      {
        sample[i]  = (*uniformGenerator)();             // sample in [0,1]
        sample[i] -= 0.5;
        range      = SQRT3*stdDev[i];
        sample[i] *= range + range;
        sample[i] += mean[i];
      }
      break;
    case NORMAL:
    default:
      for(i=0; i<dimensions; i++)
      {
        if(odd)
        {
          s = 2.;
          while(s > 1.)
          {
            u1          = (*uniformGenerator)();        // returns 0<=x<=1
            v2          = u1+u1-1.;
            v2_sq       = v2*v2;
            u1          = (*uniformGenerator)();
            v1          = u1 + u1 - 1.;
            v1_sq       = v1*v1;
            s           = v1_sq + v2_sq;
          }
          ln_s          = log(s);
          sqrt_lns      = sqrt(-(ln_s+ln_s)/s);
          odd           = 0;
          sample[i]     = v1*sqrt_lns;
        }
        else
        {
          odd           = 1;
          sample[i]     = v2*sqrt_lns;
        }
        if(i==1)
        {
          sample[i]     = stdDev[i]*(rho*temp + sample[i]*sqrt(1.-rho*rho));
          sample[i]     += mean[i];
        }
        else
        {
          temp          = sample[i];
          sample[i]     *= stdDev[i];
          sample[i]     += mean[i];
        }
      }
      break;
    case EXPONENTIAL_DIST:
      for(i=0; i<dimensions; i++)
      {
        alpha   = mean[i] - stdDev[i];
        u1      = (*uniformGenerator)();
        if(u1>0.)
          sample[i] = -stdDev[i]*log(u1);
        else
          sample[i] = 0.;
          sample[i] += alpha;
      }
      break;
    case ERLANG_DIST:
      for(i=0; i<dimensions; i++)
      {
        a   = mean[i]/variance[i];
        k   = ROUND(mean[i]*a);
        tr  = 1.;
        for(j=0; j<k; j++)
          tr  *= (*uniformGenerator)();
        if(tr>0.)
          sample[i] = -log(tr)/a;
        else
          sample[i] = 0.;
      }
      break;
    case LOG_NORMAL:
      for(i=0; i<dimensions; i++)
      {
        std_devx        = stdDev[i];
        term            = log(variance[i]/mean[i]/mean[i] + 1.);
        std_devy        = sqrt(term);
        if(mean[i]>0.)
          mean_y        = log(mean[i]) - 0.5*term;
        else
          mean_y        = 0.;
        sum             = -6.0;
        for(j=0; j<12; j++)
          sum           += (*uniformGenerator)();
        sample[i]       = exp(mean_y + std_devy*sum);
      }
      break;
    case GAUSS_MARKOV:
      for(i=0; i<dimensions; i++)
      {
        if(odd)
        {
          s = 2.;
          while(s > 1.)
          {
            u1          = (*uniformGenerator)();
            v2          = u1+u1-1.;
            v2_sq       = v2*v2;
            u1          = (*uniformGenerator)();
            v1          = u1 + u1 - 1.;
            v1_sq       = v1*v1;
            s           = v1_sq + v2_sq;
          }
          ln_s          = log(s);
          sqrt_lns      = sqrt(-(ln_s+ln_s)/s);
          odd           = 0;
          sample[i]     = rho*markov_old + sqrt(1.-rho*rho)*v1*sqrt_lns;
          markov_old    = sample[i];
        }
        else
        {
          odd           = 1;
          sample[i]     = rho*markov_old + sqrt(1.-rho*rho)*v2*sqrt_lns;
          markov_old    = sample[i];
        }
      }
      break;
  }
  return sample;
}

/*###############################*
 * Return the value of our den-  *
 * sity function evaluated at x  *
 *###############################*/
double Probability::densityFunctionAtX( double * x)
{
  int    i, j, n;
  int    in_range;
  double p_of_x, range, volume, temp, mu, sigma_sq, sum;
  double alpha, c;
  
  switch(distribution)
  {
    case UNIFORM:
      in_range = YES;
      volume   = 1.;
      for(i=0; i<dimensions; i++)
      {
        range   = SQRT3*stdDev[i];
        volume *= range;
        if( (x[i]<(mean[i]+range)) && (x[i]>=(mean[i]-range)) )
          ;
        else
        {
          in_range = NO;
          break;
        }
      }
      if(in_range)
        p_of_x = 1./volume;
      else
        p_of_x = 0.;
      break;
    case NORMAL:
    default:
      temp   = (1-rho*rho);
      if(dimensions == 2)
        p_of_x = exp(rho*(x[0]-mean[0])*(x[1]-mean[1])/stdDev[0]/
                 stdDev[1]/temp);
      else
        p_of_x = 1.;
      temp *= 2.;
      for(i=0; i<dimensions; i++)
      {
        p_of_x *= exp(-(x[i]-mean[i])*(x[i]-mean[i])/temp/variance[i]);
        p_of_x /= SQRT2PI*stdDev[i];
      }
      break;
    case EXPONENTIAL_DIST:
      in_range = YES;
      for(i=0; i<dimensions; i++)
      {
        alpha   = mean[i] - stdDev[i];
        if( x[i]>=alpha )
          ;
        else
        {
          in_range = NO;
          break;
        }
      }
      if(in_range)
      {
        p_of_x = 1.;
        for(i=0; i<dimensions; i++)
        {
          p_of_x *= exp(-(x[i]-mean[i])/stdDev[i]);
          p_of_x /= stdDev[i];
        }
      }
      else
        p_of_x = 0.;
      break;
    case ERLANG_DIST:       /* see Papoulis p. 77 */
      in_range = YES;
      p_of_x = 1.;
      for(i=0; i<dimensions; i++)
      {
        if(x[i] <= 0.)
       {
         in_range = NO;
         break;
       }
       c   = mean[i]/variance[i];
       n   = ROUND(mean[i]*c);
       sum = 0.;
       for(j=1; j<n; j++)
         sum += log((double)j);
       temp = n*log(c) + (n-1)*log(x[i]) - c*x[i] - sum;
       p_of_x *= exp(temp);
      }
      if(!in_range)
        p_of_x = 0.;
      break;
    case LOG_NORMAL:        /* see Whalen p. 23 */
      in_range = YES;
      p_of_x = 1.;
      for(i=0; i<dimensions; i++)
      {
        if(x[i] <= 0.)
        {
          in_range = NO;
          break;
        }
        if(mean[i]>0.)
        {
          sigma_sq = log(variance[i] + mean[i]*mean[i]) - 2.*log(mean[i]);
          mu       = log(mean[i]) - sigma_sq/2;
        }
        else
          mu = sigma_sq = 0.;
        temp = log(x[i]) - mu;
        p_of_x *= exp(-temp*temp/2./sigma_sq)/SQRT2PI/sqrt(sigma_sq)/x[i];
      }
      if(!in_range)
        p_of_x = 0.;
      break;
  }
  return p_of_x;
}

