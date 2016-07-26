/********************************************************************************
 *                                                                              *
 * This subclass of object implements a general probability distribution        *
 * function.                                                                    *
 *                                                                              *
 * File: /User/frank/Objc_Classes/Probability/Probability.m                     *
 *                                                                              *
 * Revision History:                                                            *
 ********************************************************************************/

#define NORMAL              0
#define UNIFORM             1
#define EXPONENTIAL_DIST    2
#define ERLANG_DIST         3                   // special case of Gamma dist
#define LOG_NORMAL          4
#define GAUSS_MARKOV        5

#define SEED1           2723                    // Default seeds
#define SEED2           177

class   MLCG;
class   Uniform;

class Probability
{
private:
  MLCG          *linearCongruential;            // Required for random number generators
  Uniform       *uniformGenerator;              // Uniform number generator
  double        *mean;                          // mean of the distribution
  double        *variance;                      // variance
  double        *stdDev;                        // standard deviation
  double        rho;                            // correlation coefficient
  double        probability;                    // Probability of occurrence
  int           distribution;                   // distribution type, see above
  int           dimensions;                     // # of dimensions in distribution
  double        *doubleOutputs;                 // Used for array outputs
//
// Private functions:
//
  void          initMeanVariance();             // initialize for new means, vars

public:
/*********************************
 * Constructors/destructors:     *
 *********************************/
  Probability();                            // Class Constructor
  ~Probability();                           // Class Destructor

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
  void setProbDistribution( int newDistribution);
  void setProbability( double newProbability);
  void setDimensions( int newDimension);
  void setNewSeeds(long seed1, long seed2);

/********************************
 * These methods generate       *
 * new PN data.                 *
 ********************************/
  double * getMean();                       // Get parameters
  double * getVariance();
  int      getProbDistribution();
  double   getProbability();
  int      getDimensions();

//
// Getting outputs:
//
  double *newSample();                      // Return samples
  double  densityFunctionAtX(double *x);    // Return PDF
  
};

