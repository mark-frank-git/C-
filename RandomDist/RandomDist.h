#ifndef _RANDOM_DIST_H
#define _RANDOM_DIST_H 1
/********************************************************************************
 *                                                                              *
 * This subclass of object implements a class for generating samples from       *
 * different types of distributions.                                            *
 *                                                                              *
 * File: /User/frank/Objc_Classes/RandomDist/RandomDist.h                       *
 *                                                                              *
 ********************************************************************************/

#define NORMAL                  0
#define UNIFORM                 1
#define EXPONENTIAL_DIST        2
#define ERLANG_DIST             3                       // special case of Gamma dist
#define LOG_NORMAL              4
#define GAUSS_MARKOV            5
#define RICIAN_DIST             6
#define RAYLEIGH_DIST           7
#define CHI_SQUARED_DIST        8
#define NON_CENTRAL_DIST        9
#define CONSTANT_DIST           10                      // Used for testing RadarScatterers

#ifndef LOGICAL
#define LOGICAL                 char
#endif

class   UniformNumber;

class RandomDist
{
private:
  UniformNumber *_uniformGenerator;             // Random number generator
  
  int           _distributionType;              // distribution type, see above
  int           _dimensions;                    // # of dimensions in distribution
  int           _oddSample;                     // 1 = odd sample
  double        *_mean;                         // mean of the distribution
  double        *_variance;                     // variance
  double        *_stdDev;                       // standard deviation
  double        *_ricianAlpha;                  // Calculated Rician alphas from means
  double        _rho;                           // correlation coefficient
  double        _probability;                   // Probability of occurrence
  double        _markovOld;                     // Old value of Markov output
  double        *_doubleOutputs;                // Used for array outputs
//
// Private functions:
//
  void          initMeanVariance();             // initialize for new means, vars
  double        findRiceAlphaFromMean(double mean);     // Local for Rician distribution

  double        *normalSample();                // Return samples from the various distributions
  double        *uniformSample();
  double        *exponentialSample();
  double        *erlangSample();
  double        *lognormalSample();
  double        *gaussMarkovSample();
  double        *ricianSample();
  double        *rayleighSample();
  double        *constantSample();

  double        normalDensityAt(double *x);     // Return density functions for the various distributions
  double        uniformDensityAt(double *x);
  double        exponentialDensityAt(double *x);
  double        erlangDensityAt(double *x);
  double        lognormalDensityAt(double *x);
  double        gaussMarkovDensityAt(double *x);
  double        ricianDensityAt(double *x);
  double        rayleighDensityAt(double *x);
  double        chiSquaredDensityAt(double *x);
  double        nonCentralDensityAt(double *x);
  double        oldnonCentralDensityAt(double *x);

public:
/*********************************
 * Constructors/destructors:     *
 *********************************/
  RandomDist();                            // Class Constructor
  ~RandomDist();                           // Class Destructor

/*******************************
 * These methods initialize:    *
 *******************************/
  void  initGenerator();

/*******************************
 * These methods set parameters:*
 *******************************/
  void setMean(double *newMean);
  void setVariance(double* newVariance);
  void setCorrelationCoefficient( double newRho);
  void setProbability( double probability);
  void setDistributionType( int newDistribution);
  void setDimensions( int newDimension);
  void setSeed(long seed);

/********************************
 * These methods get parameters *
 * from the object.             *
 ********************************/
  const double  *mean()                 {return _mean;}                         // Get parameters
  const double  *variance()             {return _variance;}
  double        correlationCoefficient(){return _rho;}
  double        probability()           {return _probability;}
  int           distributionType()      {return _distributionType;}
  int           dimensions()            {return _dimensions;}

//
// Getting outputs:
//
  double *newSample();                      // Return samples
  double  densityFunctionAtX(double *x);    // Return PDF
  
};

#endif
