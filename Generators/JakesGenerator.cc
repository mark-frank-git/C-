/********************************************************************************
 * This class implements a fading simulator based on Jakes' model, resulting    *
 * in an envelope that is Rayleigh distributed. Provisions are made for the     *
 * user's _velocity.                                                            *
 * File: /User/frank/C++/Generators/JakesGenerator.cc                           *
 *                                                                              *
 ********************************************************************************/
#include <C_Libraries/constants.h>
#include "JakesGenerator.h"                                                     // Class definitions
#include <stdio.h>
#include <math.h>

#define RANDOM_MULT     11
#define RANDOM_OFFSET   7

#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define ROUND(a)        ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )


#define USE_SINE_TABLE          1
#define NUMBER_SINE_POINTS      72
double sine_table[NUMBER_SINE_POINTS+1] =
{
0, 0.0871557427590095, 0.173648177689373, 0.25881904513554, 0.342020143368499, 0.422618261792335, 0.500000000059208, 0.573576436416384, 0.64278760975637, 0.707106781259063, 0.766044443192222, 0.819152044360885, 0.866025403852807, 0.906307787099253, 0.93969262084047, 0.965925826333307, 0.984807753043868, 0.996194698108631, 1, 0.996194698072879, 0.984807752972637, 0.965925826227138, 0.939692620700171, 0.906307786925892, 0.866025403647704, 0.8191520441256, 0.766044442928546, 0.707106780969003, 0.642787609442133, 0.573576436080362, 0.499999999703959, 0.422618261420561, 0.34202014298303, 0.25881904473931, 0.173648177285398, 0.0871557423503633, -4.10206857022848e-10, -0.087155743167655, -0.173648178093348, -0.258819045531769, -0.342020143753967, -0.422618262164108, -0.500000000414458, -0.573576436752406, -0.642787610070607, -0.707106781549123, -0.766044443455898, -0.81915204459617, -0.866025404057911, -0.906307787272614, -0.939692620980769, -0.965925826439476, -0.9848077531151, -0.996194698144382, -1, -0.996194698037127, -0.984807752901405, -0.965925826120968, -0.939692620559872, -0.906307786752531, -0.8660254034426, -0.819152043890315, -0.76604444266487, -0.707106780678943, -0.642787609127897, -0.573576435744341, -0.499999999348709, -0.422618261048787, -0.342020142597562, -0.258819044343081, -0.173648176881423, -0.0871557419417174, 0.
};
double cosine_table[NUMBER_SINE_POINTS+1] =
{
1, 0.996194698090755, 0.984807753008253, 0.965925826280222, 0.939692620770321, 0.906307787012573, 0.866025403750255, 0.819152044243243, 0.766044443060384, 0.707106781114033, 0.642787609599252, 0.573576436248373, 0.499999999881584, 0.422618261606448, 0.342020143175764, 0.258819044937425, 0.173648177487386, 0.0871557425546864, -2.05103428511424e-10, -0.0871557429633323, -0.173648177891361, -0.258819045333654, -0.342020143561233, -0.422618261978222, -0.500000000236833, -0.573576436584394, -0.642787609913488, -0.707106781404093, -0.76604444332406, -0.819152044478527, -0.866025403955359, -0.906307787185934, -0.93969262091062, -0.965925826386392, -0.984807753079484, -0.996194698126507, -1, -0.996194698055003, -0.984807752937021, -0.965925826174053, -0.939692620630022, -0.906307786839212, -0.866025403545152, -0.819152044007958, -0.766044442796708, -0.707106780823973, -0.642787609285015, -0.573576435912351, -0.499999999526334, -0.422618261234674, -0.342020142790297, -0.258819044541196, -0.173648177083411, -0.0871557421460412, 6.1531028553512e-10, 0.0871557433719789, 0.173648178295336, 0.258819045729884, 0.342020143946701, 0.422618262349996, 0.500000000592083, 0.573576436920416, 0.642787610227725, 0.707106781694153, 0.766044443587735, 0.819152044713812, 0.866025404160462, 0.906307787359295, 0.939692621050919, 0.965925826492561, 0.984807753150716, 0.996194698162258, 1.
};

#define A    16807
#define M    2147483647
#define Q    127773                     /* M div A */
#define R    2836                       /* m mod A */
/***************************************************************
 *
 * double urand()
 *
 *   implicit long input variables
 *   -----------------------------
 *   seed = initialize seed the first time this routine is
 *           called.
 *
 *  This routine generates a random number in [0,1] based on
 *  the algorithm given in [PAR88] S.K. Park and K.W. Miller,
 *  "Random numver generators: good ones are hard to find,"
 *  Comm. ACM, vol. 32, no. 10, pp. 1192-1201, Oct. 1988.
 *
 **************************************************************/
double JakesGenerator::urand()
{
  long   lo,hi,test;
  double rand;

  hi            = _randomSeed/Q;
  lo            = _randomSeed%Q;
  test          = A*lo - R*hi;

  if(test > 0)
    _randomSeed = test;
  else
    _randomSeed = test + M;
  rand          = (double)_randomSeed/M;
  return rand;
}

// ############################# Private Method ###############################
// updateParameters -- Update generator parameters based on new _velocity
// Input:               None
//                      
// Output:              None
//
// Notes:
//  1. 
// ############################# Private Method ###############################
void JakesGenerator::updateParameters()
{
  int           n;
  double        doppler_omega, alpha_n_k;
//
// Allocate arrays:
//
  if(_oldNumberWaves < _numberSineWaves)
  {
    _oldNumberWaves             = _numberSineWaves;
    delete [] _cosThetas;
    delete [] _sinThetas;
    delete [] _cosPhiInits;
    delete [] _sinPhiInits;
    delete [] _omegaMCosAlphaDeltaT;
    delete [] _omegaMSinAlphaDeltaT;
    _cosThetas                  = new double[_numberSineWaves/NUMBER_QUADRANTS];
    _sinThetas                  = new double[_numberSineWaves/NUMBER_QUADRANTS];
    _cosPhiInits                = new double[_numberSineWaves/NUMBER_QUADRANTS];
    _sinPhiInits                = new double[_numberSineWaves/NUMBER_QUADRANTS];
    _omegaMCosAlphaDeltaT       = new double[_numberSineWaves/NUMBER_QUADRANTS];
    _omegaMSinAlphaDeltaT       = new double[_numberSineWaves/NUMBER_QUADRANTS];
  }
//
// initialize gain
//
  _normalizingGain              = 2./sqrt((double)_numberSineWaves);
  _tableScale                   = NUMBER_SINE_POINTS/TWOPI;
//
// initialize inital phases
//
  _randomSeed                   = RANDOM_MULT * _fadeNumber + RANDOM_OFFSET;
  for(n=0; n<_numberSineWaves/NUMBER_QUADRANTS;n++)
  {
    _cosPhiInits[n]             = TWOPI*urand();
    _sinPhiInits[n]             = TWOPI*urand();
  }
//
// Find Doppler frequency
//
  doppler_omega                 = TWOPI*_dopplerFreq;  

// 
// We need to check the new generator parameters, and update the jakes
// generator oscillators
//
  for(n=0; n<_numberSineWaves/NUMBER_QUADRANTS; n++)
  {
    alpha_n_k                   = TWOPI*n/_numberSineWaves + TWOPI*_fadeNumber/_numberSineWaves/_numberFaders
                                + PI/2./_numberSineWaves/_numberFaders;
    _omegaMCosAlphaDeltaT[n]    = doppler_omega*cos(alpha_n_k)*_samplingTime;
    _omegaMSinAlphaDeltaT[n]    = doppler_omega*sin(alpha_n_k)*_samplingTime;
     _cosThetas[n]              = _cosPhiInits[n];
     _sinThetas[n]              = _sinPhiInits[n];
  }

  return;
}

// ############################# Class Constructor ###############################
// JakesGenerator -- Constructor for the JakesGenerator class
// Input:       fadeNumber:             The # of the fader (for independent faders)
//              numberSineWaves:        The total # of sine waves for the generator
//              numberFaders:           The total # of faders (for independent faders)
//              doppler:                Doppler frequency in Hertz.
//              sampleTime:             Time between samples in seconds.
//                      
// Output:                              None
//
// Notes:
// 1. The walshNumber is used to insure orthogonal simulators.  E.g., use a different
//    walshNumber for different satellites.  See P. Dent, et. al., Electronics Letters,
//    June 24, 1993.
// ############################# Class Constructor ###############################
JakesGenerator::JakesGenerator(int fadeNumber, int numberSineWaves, int numberFaders, float doppler,
                 float sampleTime)
{
// 
// Initialize instance variables:
//
  _numberSineWaves              = 1;
  _oldNumberWaves               = 0;
  _cosThetas                    = NULL;
  _sinThetas                    = NULL;
  _cosPhiInits                  = NULL;
  _sinPhiInits                  = NULL;
  _omegaMCosAlphaDeltaT         = NULL;
  _omegaMSinAlphaDeltaT         = NULL;
  _iData                        = _qData                = _magnitudeData        = NULL;
  _numberSimulatePoints         = _oldNumberPoints      = 0;
 
  setFadeNumber(fadeNumber);
  setNumberSineWaves(numberSineWaves);
  setNumberFaders(numberFaders);
  setDopplerFreq(doppler);
  setSamplingTime(sampleTime);
  initGenerator();

  return;
}

// ############################# Class Destructor ###############################
// ~JakesGenerator -- Destructor for the JakesGenerator class
// Input:                       None
//                      
// Output:                      None
// ############################# Class Destructor ###############################
JakesGenerator::~JakesGenerator()
{
// 
// Delete allocated data:
//
  delete        [] _cosThetas;
  delete        [] _sinThetas;
  delete        [] _cosPhiInits;
  delete        [] _sinPhiInits;
  delete        [] _omegaMCosAlphaDeltaT;
  delete        [] _omegaMSinAlphaDeltaT;
  delete        [] _iData;
  delete        [] _qData;
  delete        [] _magnitudeData;

  return;
}

// ############################# Public Method ###############################
// initGenerator -- Initializes the generator
// Input:               None
//                      
// Output:              None
// ############################# Public Method ###############################
void JakesGenerator::initGenerator()
{
  int   i;
  for(i=0; i<_numberSineWaves/NUMBER_QUADRANTS; i++)
  {
    _cosThetas[i]       = _cosPhiInits[i];
    _sinThetas[i]       = _sinPhiInits[i];
  }

  return;
}

// ############################# Public Method ###############################
// setFadeNumber -- Sets a new fade number
// Input:       number:         New fade number for independent faders
//                      
// Output:                              None
// ############################# Public Method ###############################
void JakesGenerator::setFadeNumber(int number)
{
  _fadeNumber           = MAX(number, 0);
  updateParameters();

  return;
}

// ############################# Public Method ###############################
// setNumberSineWaves -- Sets a new number of sine waves
// Input:       number:         New number of sine waves must be a factor of 4
//                      
// Output:                              None
// ############################# Public Method ###############################
void JakesGenerator:: setNumberSineWaves(int number)
{
  _numberSineWaves      = MAX(number, 1);
  while(_numberSineWaves%4 != 0)
    _numberSineWaves++;
  updateParameters();

  return;
}

// ############################# Public Method ###############################
// setNumberFaders -- Sets a new number of independent faders
// Input:       number:         New number of independent faders
//                      
// Output:                              None
// ############################# Public Method ###############################
void JakesGenerator:: setNumberFaders(int number)
{
  _numberFaders         = MAX(number, 1);
  updateParameters();

  return;
}

// ############################# Public Method ###############################
// setSamplingTime -- Sets a new sampling time
// Input:       sampleTime:             Time between samples in seconds.
//                      
// Output:                              None
// ############################# Public Method ###############################
void JakesGenerator::setSamplingTime(float sampleTime)
{
  _samplingTime                 = sampleTime;
  updateParameters();

  return;
}

// ############################# Public Method ###############################
// setDopplerFreq -- Sets a new Doppler frequency in Hertz
// Input:       doppler:        Doppler frequency in Hertz.
//                      
// Output:                      None
// ############################# Public Method ###############################
void JakesGenerator::setDopplerFreq(float doppler)
{
  _dopplerFreq                  = MAX(doppler, MIN_DOPPLER);
  updateParameters();

  return;
}

// ############################# Public Method ###############################
// generateIAndQData -- Generates I and Q scattered data
// Input:       numberPoints:           Number of samples to generate
//                      
// Output:                              None
//
// Notes:
//  1. If you want I and Q data, call this function first, then request I and Q data.
//  2. This function is still non-functional
// ############################# Public Method ###############################
void JakesGenerator::generateIAndQData(long numberPoints)
{
  long i;
  int  n, index;
  double xc, xs;

  _numberSimulatePoints         = numberPoints;
  if(_oldNumberPoints < _numberSimulatePoints)
  {
    _oldNumberPoints            = _numberSimulatePoints;
    if(_iData != NULL)
    {
      delete [] _iData;
      delete [] _qData;
      delete [] _magnitudeData;
    }
    _iData                      = new float [_numberSimulatePoints];
    _qData                      = new float [_numberSimulatePoints];
    _magnitudeData              = new float [_numberSimulatePoints];
  }
//
// Here is where we generate the data:
//
  for(i=0; i<_numberSimulatePoints; i++)
  {
    xc = xs = 0.;
    for(n=0; n<_numberSineWaves/NUMBER_QUADRANTS; n++)
    {
#if USE_SINE_TABLE
      index                     = ROUND(_tableScale*_cosThetas[n]);
      xc                        += cosine_table[index];
      index                     = ROUND(_tableScale*_sinThetas[n]);
      xs                        += sine_table[index];
      _cosThetas[n]             += _omegaMCosAlphaDeltaT[n];
      _sinThetas[n]             += _omegaMSinAlphaDeltaT[n];
      while(_cosThetas[n] > TWOPI)
        _cosThetas[n]           -= TWOPI;
      while(_cosThetas[n] < 0.)
        _cosThetas[n]           += TWOPI;
      while(_sinThetas[n] > TWOPI)
        _sinThetas[n]           -= TWOPI;
      while(_sinThetas[n] < 0.)
        _sinThetas[n]           += TWOPI;
#else
      xs                        += sin(_sinThetas[n]);
      xc                        += cos(_cosThetas[n]);
      _cosThetas[n]             += _omegaMCosAlphaDeltaT[n];
      _sinThetas[n]             += _omegaMSinAlphaDeltaT[n];
      
#endif
    }
    _iData[i] = _normalizingGain*xc;                    /* Normalize power              */
    _qData[i] = _normalizingGain*xs;
  }

  return;
}

// ############################# Public Method ###############################
// complexFading -- Generates I and Q fading data
// Input:                               None
//                      
// Output:                              Complex I and Q fading data
//
// Notes:
// ############################# Public Method ###############################
Complex JakesGenerator::complexFading()
{
  double        xc, xs;
  Complex       fading;
  int           n, index;
//
// Sum the sinusoids:
//
  xc = xs = 0.;
  for(n=0; n<_numberSineWaves/NUMBER_QUADRANTS; n++)
  {
#if USE_SINE_TABLE
      index                     = ROUND(_tableScale*_cosThetas[n]);
      xc                        += cosine_table[index];
      index                     = ROUND(_tableScale*_sinThetas[n]);
      xs                        += sine_table[index];
      _cosThetas[n]             += _omegaMCosAlphaDeltaT[n];
      _sinThetas[n]             += _omegaMSinAlphaDeltaT[n];
      while(_cosThetas[n] > TWOPI)
        _cosThetas[n]           -= TWOPI;
      while(_cosThetas[n] < 0.)
        _cosThetas[n]           += TWOPI;
      while(_sinThetas[n] > TWOPI)
        _sinThetas[n]           -= TWOPI;
      while(_sinThetas[n] < 0.)
        _sinThetas[n]           += TWOPI;
#else
      xs                        += sin(_sinThetas[n]);
      xc                        += cos(_cosThetas[n]);
      _cosThetas[n]             += _omegaMCosAlphaDeltaT[n];
      _sinThetas[n]             += _omegaMSinAlphaDeltaT[n];
      
#endif
  }
  fading                        = _normalizingGain*Complex(xc, xs);                     /* Normalize power              */

  return fading;
}

// ############################# Public Method ###############################
// calculateFadeMagnitude -- Generates fade magnitude data from I and Q data
// Input:                       None
//                      
// Output:                      an array of magnitude data given in dB
//
// Notes:
//  1. You must first call generateIandQData() before calling this function
// ############################# Public Method ###############################
const float *JakesGenerator::calculateFadeMagnitude()
{
  double temp;

  for(long i=0; i<_numberSimulatePoints; i++)
  {
    temp                = _iData[i]*_iData[i] + _qData[i]*_qData[i];
    if(temp>0.)
      _magnitudeData[i] = 10.*log10(temp);
    else
      _magnitudeData[i] = 0.;
  }

  return _magnitudeData;
}
