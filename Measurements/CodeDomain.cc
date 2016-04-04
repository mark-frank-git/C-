/************************************************************************
 *                                                                      *
 * This class is used for making code domain power measurements on CDMA *
 * QPSK symbols.                                                        *
 *                                                                      *
 * File: /User/frank/C++/Measurements/CodeDomain.cc                     *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 03/20/99 - Started.                                              *
 *  2. 05/31/99 - Add calculateRho().                                   *
 ************************************************************************/

#include "CodeDomain.h"
#include <Generators/WalshCodes.h>
#include <C_Libraries/constants.h>
#include <stdio.h>
#include <math.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ROUND(a) ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )


// ############################# Class Constructor #################################
// CodeDomain -- Constructor for the CodeDomain class
//
// Input:       length:         Walsh length in chips
//              samples:        Samples/chip in input signal
//
// Output:              None
// ############################# Class Constructor #################################
CodeDomain::CodeDomain(int length, int samples)
{
// 
// Initialize instance variables:
//
  m_walshGen                            = NULL;
  m_codePowers                          = NULL;
  m_integrateIWalsh                     = NULL;
  m_integrateQWalsh                     = NULL;
  m_crossIWalsh                         = NULL;
  m_crossQWalsh                         = NULL;
  m_idealIWalsh                         = NULL;
  m_idealQWalsh                         = NULL;
  m_walshCodes                          = NULL;
  m_walshLength                         = 0;
  setMinCodePower(MIN_CODE_POWER);
  setSamplesPerChip(samples);
  setWalshLength(length);
  resetStatistics();
  return;
}


// ############################# Class Destructor ###############################
// CodeDomain -- Destructor for the CodeDomain class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
CodeDomain::~CodeDomain()
{
  delete [] m_codePowers;
  delete [] m_integrateIWalsh;
  delete [] m_integrateQWalsh;
  delete [] m_crossIWalsh;
  delete [] m_crossQWalsh;
  delete [] m_idealIWalsh;
  delete [] m_idealQWalsh;
  delete [] m_walshCodes;
  delete m_walshGen;
  
  return;
}

// ############################ Public Function ###################################
// resetStatistics - Reset the accumulated statistics counters
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void CodeDomain::resetStatistics()
{
  m_accumulatedWalshPeriods     = 0;

  m_meanErrorAngle              = 0.;
  m_varianceErrorAngle          = 0.;
  return;
}

// ############################ Public Function ###################################
// setMinCodePower - Sets a new minimum code power for regenerating Walsh signal
//
// Input:       power:          new minimum code power
//
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void CodeDomain::setMinCodePower(float power)
{
  m_minCodePower        = power;                        // No error check
  return;
}
// ############################ Public Function ###################################
// setSamplesPerChip - Sets a new number of samples/chip.
//
// Input:       samples:        new number of samples/chip in input signal
//
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void CodeDomain::setSamplesPerChip(int samples)
{
  m_samplesPerChip      = samples;
  setWalshLength(m_walshLength);                // allocate ideal arrays
  return;
}

// ############################ Public Function ###################################
// setWalshLength - Sets a new walsh code length.
//
// Input:       length:         new Walsh code length
//
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void CodeDomain::setWalshLength(int length)
{
  int   n;
//
// Delete old storage:
//
  delete m_walshGen;
  delete [] m_codePowers;
  delete [] m_integrateIWalsh;
  delete [] m_integrateQWalsh;
  delete [] m_crossIWalsh;
  delete [] m_crossQWalsh;
  delete [] m_idealIWalsh;
  delete [] m_idealQWalsh;
//
// Allocate new storage:
//
  m_walshLength         = length;
  m_walshGen            = new WalshCodes(length);
  m_codePowers          = new double[length];
  m_integrateIWalsh     = new double[length];
  m_integrateQWalsh     = new double[length];
  m_crossIWalsh         = new double[length];
  m_crossQWalsh         = new double[length];
  m_idealIWalsh         = new float[m_samplesPerChip*length];
  m_idealQWalsh         = new float[m_samplesPerChip*length];
  m_walshCodes          = new const char * [length];
  for(n=0; n<length; n++)
    m_walshCodes[n]     = m_walshGen->codeForIndex(n);
  return;
}

#define MIN_POWER       -100.

// ############################ Public Function ###################################
// processRealWalshPeriod - Process an array of real samples over a single
//                             Walsh period.
//
// Input:       inputI:         input real array after de-spreading
//
// Output:                      None
//
// Notes:
// 1. The size of the above arrays must be greater than or equal to the number of
//    chips in a Walsh period (i.e., 64 for IS-95) * number_samples_per_chip.
// 2. The outputs of this processing are stored in m_integrateIWalsh[].
// ############################ Public Function ###################################
void CodeDomain::processRealWalshPeriod(const float *inputI)
{
  int           n, i;
  int           walsh_index, chip_counter;
//
// Reset integrators:
//
  for(i=0; i<m_walshLength; i++)
    m_integrateIWalsh[i]        = m_integrateQWalsh[i]  = 0.;
  walsh_index   = 0;
  chip_counter  = 0;
  for(n=0; n<m_walshLength*m_samplesPerChip; n++)
  {
    for(i=0; i<m_walshLength; i++)                              // Loop over walsh codes
    {
      if(m_walshCodes[i][walsh_index] == 0)                     // Note: use (0,1)->(1,-1) mapping
      {
        m_integrateIWalsh[i]    += inputI[n];
      }
      else
      {
        m_integrateIWalsh[i]    -= inputI[n];
      }
    }
    if(++chip_counter == m_samplesPerChip)
    {
      walsh_index++;
      chip_counter      = 0;
    }
  }     // End loop over Walsh period
  return;
}

// ############################ Public Function ###################################
// processComplexWalshPeriod - Process an array of complex samples over a single
//                             Walsh period.
//
// Input:       inputI:         input real array after de-spreading
//              inputQ:         input imag array after de-spreading
//
// Output:                      None
//
// Notes:
// 1. The size of the above arrays must be greater than or equal to the number of
//    chips in a Walsh period (i.e., 64 for IS-95) * number_samples_per_chip.
// 2. The outputs of this processing are stored in m_integrateIWalsh[], and
//    m_integrateQWalsh[].
// ############################ Public Function ###################################
void CodeDomain::processComplexWalshPeriod(const float *inputI, const float *inputQ)
{
  int           n, i;
  int           walsh_index, chip_counter;
//
// Reset integrators:
//
  for(i=0; i<m_walshLength; i++)
    m_integrateIWalsh[i]        = m_integrateQWalsh[i]  = 0.;
  walsh_index   = 0;
  chip_counter  = 0;
  for(n=0; n<m_walshLength*m_samplesPerChip; n++)
  {
    for(i=0; i<m_walshLength; i++)                              // Loop over walsh codes
    {
      if(m_walshCodes[i][walsh_index] == 0)                     // Note: use (0,1)->(1,-1) mapping
      {
        m_integrateIWalsh[i]    += inputI[n];
        m_integrateQWalsh[i]    += inputQ[n];
      }
      else
      {
        m_integrateIWalsh[i]    -= inputI[n];
        m_integrateQWalsh[i]    -= inputQ[n];
      }
    }
    if(++chip_counter == m_samplesPerChip)
    {
      walsh_index++;
      chip_counter      = 0;
    }
  }     // End loop over Walsh period
  return;
}

// ############################ Public Function ###################################
// processCrossWalshPeriod - Process an array of complex samples over a single
//                           Walsh period.
//
// Input:       inputI:         input real array after cross-spreading
//              inputQ:         input imag array after cross-spreading
//
// Output:                      None
//
// Notes:
// 1. The size of the above arrays must be greater than or equal to the number of
//    chips in a Walsh period (i.e., 64 for IS-95) * number_samples_per_chip.
// 2. The outputs of this processing are stored in m_crossIWalsh[], and
//    m_crossQWalsh[].
// ############################ Public Function ###################################
void CodeDomain::processCrossWalshPeriod(const float *inputI, const float *inputQ)
{
  int           n, i;
  int           walsh_index, chip_counter;
//
// Reset integrators:
//
  for(i=0; i<m_walshLength; i++)
    m_crossIWalsh[i]    = m_crossQWalsh[i]      = 0.;
  walsh_index   = 0;
  chip_counter  = 0;
  for(n=0; n<m_walshLength*m_samplesPerChip; n++)
  {
    for(i=0; i<m_walshLength; i++)                              // Loop over walsh codes
    {
      if(m_walshCodes[i][walsh_index] == 0)                     // Note: use (0,1)->(1,-1) mapping
      {
        m_crossIWalsh[i]        += inputI[n];
        m_crossQWalsh[i]        += inputQ[n];
      }
      else
      {
        m_crossIWalsh[i]        -= inputI[n];
        m_crossQWalsh[i]        -= inputQ[n];
      }
    }
    if(++chip_counter == m_samplesPerChip)
    {
      walsh_index++;
      chip_counter      = 0;
    }
  }     // End loop over Walsh period
  return;
}

// ############################ Public Function ###################################
// findCodePowers - Process an array of input samples to find code domain powers.
//
// Input:       inputArray:     input samples array
//              size:           Size of the input array
//
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void CodeDomain::findCodePowers(const float *inputArray, int size)
{
  int           n, i, walsh_periods, index, increment;
  double        sum;
//
// Reset integrators:
//
  for(i=0; i<m_walshLength; i++)
    m_codePowers[i]     = 0.;
//
// calculate number of walsh periods:
//
  increment                     = m_walshLength*m_samplesPerChip;
  if(increment > 0)
    walsh_periods               = size/increment;
  else
    walsh_periods               = 0;
  index                         = 0;
  for(n=0; n<walsh_periods; n++)
  {
//
// Process the Walsh codes over a single period:
//
    processRealWalshPeriod(&inputArray[index]);
    index                       += increment;
//
// Sum the current Walsh period:
//
    for(i=0; i<m_walshLength; i++)                              // Average the code powers
    {
      sum               = m_integrateIWalsh[i];
      m_codePowers[i]   += sum*sum;
    }
  }     // End loop over number of Walsh periods
//
// Convert code domain to dB (referenced to pilot):
//
  if(m_codePowers[0] > 0.)
  {
    for(i=1; i<m_walshLength; i++)
    {
      if(m_codePowers[i] > 0.)
        m_codePowers[i] = 10.*log10(m_codePowers[i]/m_codePowers[0]);
      else
        m_codePowers[i] = MIN_POWER;
    }
    m_codePowers[0]     = 0.;
  }
  return;
}

// ############################ Public Function ###################################
// findCodePowers - Process an array of complex samples to find code domain powers.
//
// Input:       inputI:         input real array
//              inputQ:         input imag array
//              size:           Size of the input array
//
// Output:                      None
//
// Notes:
// 1. This is a rotation insensitive version of findCodePowers
// ############################ Public Function ###################################
void CodeDomain::findCodePowers(const float *inputI, const float *inputQ, int size)
{
  int           n, i, code, walsh_periods, index, increment;
  int           start, chip_counter, walsh_index;
  double        sum_i, sum_q;
//
// Reset integrators:
//
  for(i=0; i<m_walshLength; i++)
    m_codePowers[i]     = 0.;
//
// calculate number of walsh periods:
//
  increment                     = m_walshLength*m_samplesPerChip;
  if(increment > 0)
    walsh_periods               = size/increment;
  else
    walsh_periods               = 0;
  start                         = 0;
  for(n=0; n<walsh_periods; n++)
  {
//
// Loop over Walsh codes:
//
    for(code=0; code<m_walshLength; code++)
    {
      sum_i             = sum_q = 0.;
      index             = start;
      walsh_index       = 0;
      chip_counter      = 0;
//
// Sum the current Walsh period:
//
      for(i=0; i<increment; i++)
      {
        if(m_walshCodes[code][walsh_index] == 0)                // Note: use (0,1)->(1,-1) mapping
        {
          sum_i         += inputI[index];
          sum_q         += inputQ[index];
        }
        else
        {
          sum_i         -= inputI[index];
          sum_q         -= inputQ[index];
        }
        index++;
        if(++chip_counter == m_samplesPerChip)
        {
          walsh_index++;
          chip_counter  = 0;
        }
      }  // End loop over Walsh period
      m_codePowers[code]        += (sum_i+sum_q)*(sum_i+sum_q); // Find the power for the current code
    }  // End loop over Walsh codes
    start               += increment;                           // Move to next Walsh period
    
  }     // End loop over number of Walsh periods
//
// Convert code domain to dB (referenced to pilot):
//
  if(m_codePowers[0] > 0.)
  {
    for(i=1; i<m_walshLength; i++)
    {
      if(m_codePowers[i] > 0.)
        m_codePowers[i] = 10.*log10(m_codePowers[i]/m_codePowers[0]);
      else
        m_codePowers[i] = MIN_POWER;
    }
    m_codePowers[0]     = 0.;
  }
  return;
}

// ############################ Public Function ###################################
// findCodePowersEVM - Process an array of complex samples to find code domain powers.
//
// Input:       inputI:         input real array
//              inputQ:         input imag array
//              size:           Size of the input array
//
// Output:                      None
//
// Notes:
// 1. This is a rotation insensitive version of findCodePowers
// 2. This is a different version of findCodePowers() which may overlap work load
//    with EVM calculation.
// ############################ Public Function ###################################
void CodeDomain::findCodePowersEVM(const float *inputI, const float *inputQ, int size)
{
  int           n, i, walsh_periods, index, increment;
  double        sum;
//
// Reset integrators:
//
  for(i=0; i<m_walshLength; i++)
    m_codePowers[i]     = 0.;
//
// calculate number of walsh periods:
//
  increment                     = m_walshLength*m_samplesPerChip;
  if(increment > 0)
    walsh_periods               = size/increment;
  else
    walsh_periods               = 0;
  index                         = 0;
  for(n=0; n<walsh_periods; n++)
  {
//
// Process the Walsh codes over a single period:
//
    processComplexWalshPeriod(&inputI[index], &inputQ[index]);
    index                       += increment;
//
// Sum the current Walsh period:
//
    for(i=0; i<m_walshLength; i++)                              // Average the code powers
    {
      sum               = m_integrateIWalsh[i]+m_integrateQWalsh[i];
      m_codePowers[i]   += sum*sum;
    }
  }     // End loop over number of Walsh periods
//
// Convert code domain to dB (referenced to pilot):
//
  if(m_codePowers[0] > 0.)
  {
    for(i=1; i<m_walshLength; i++)
    {
      if(m_codePowers[i] > 0.)
        m_codePowers[i] = 10.*log10(m_codePowers[i]/m_codePowers[0]);
      else
        m_codePowers[i] = MIN_POWER;
    }
    m_codePowers[0]     = 0.;
  }
  return;
}

// ############################ Public Function ###################################
// generateWalshPeriod - Generate output I and Q data streams based on averaged
//                       input data.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// 1. You must first call processComplexWalshPeriod() to find the average Walsh
//    values.
// 2. The outputs of this processing are stored in m_idealIWalsh[], and
//    m_idealQWalsh[].
// 3. This function only works for IS-95 where I data = Q data
// ############################ Public Function ###################################
void CodeDomain::generateWalshPeriod()
{
  int           n, i;
  int           walsh_index, chip_counter;
  float         max_power;
  float         *powers;
  double        sum, sum_i, sum_q;
//
// Find maximum Walsh channel power, so that we can ignore low power
// codes. Power calculation is not rotation tolerant.
//
  powers        = new float[m_walshLength];
  powers[0]     = m_integrateIWalsh[0]*m_integrateIWalsh[0] + m_integrateQWalsh[0]*m_integrateQWalsh[0];
  max_power     = powers[0];
  for(i=1; i<m_walshLength; i++)
  {
    powers[i]   = m_integrateIWalsh[i]*m_integrateIWalsh[i] + m_integrateQWalsh[i]*m_integrateQWalsh[i];
    max_power   = MAX(powers[i], max_power);
  }
//
// Generate the output signal by summing all of higher power Walsh codes:
// Note: The following assumes the same bit sequence on I and Q
// channels (as for IS-95).
//
  walsh_index   = chip_counter  = 0;
  for(n=0; n<m_walshLength*m_samplesPerChip; n++)
  {
    sum_i       = sum_q = 0.;
    for(i=0; i<m_walshLength; i++)
    {
      if(powers[i] > max_power*m_minCodePower)
      {
        sum             = (m_integrateIWalsh[i]+m_integrateQWalsh[i])/2.;
        if(m_walshCodes[i][walsh_index] == 0)                   // Note: use (0,1)->(1,-1) mapping
        {
          sum_i         += sum;
          sum_q         += sum;
        }
        else
        {
          sum_i         -= sum;
          sum_q         -= sum;
        }
      }
    }
    m_idealIWalsh[n]    = sum_i/m_walshLength;                  // scale the sum, and save it
    m_idealQWalsh[n]    = sum_q/m_walshLength;
    if(++chip_counter == m_samplesPerChip)
    {
      walsh_index++;
      chip_counter      = 0;
    }
  }     // End loop over Walsh period
  delete [] powers;
  return;
}

// ############################ Public Function ###################################
// generateWalshPeriodGeneral - Generate output I and Q data streams based on averaged
//                       input data.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// 1. You must first call processComplexWalshPeriod() to find the average Walsh
//    values.
// 2. The outputs of this processing are stored in m_idealIWalsh[], and
//    m_idealQWalsh[].
// 3. This version of generateWalshPeriod does not compensate for phase rotation.
// ############################ Public Function ###################################
void CodeDomain::generateWalshPeriodGeneral()
{
  int           n, i;
  int           walsh_index, chip_counter;
  float         max_power;
  float         *powers;
  double        sum_i, sum_q;
//
// Find maximum Walsh channel power, so that we can ignore low power
// codes. The power calculation is rotation tolerant
//
  powers        = new float[m_walshLength];
  powers[0]     = (m_integrateIWalsh[0]+m_integrateQWalsh[0]) * (m_integrateIWalsh[0]+m_integrateQWalsh[0]);
  max_power     = powers[0];
  for(i=1; i<m_walshLength; i++)
  {
    powers[i]   = (m_integrateIWalsh[i]+m_integrateQWalsh[i]) * (m_integrateIWalsh[i]+m_integrateQWalsh[i]);
    max_power   = MAX(powers[i], max_power);
  }
//
// Generate the output signal by summing all of higher power Walsh codes:
// Note: The following does not assume the same bit sequence on I and Q
// channels (IS-95 does have the same bit sequence on both).
//
  walsh_index   = chip_counter  = 0;
  for(n=0; n<m_walshLength*m_samplesPerChip; n++)
  {
    sum_i       = sum_q = 0.;
    for(i=0; i<m_walshLength; i++)
    {
      if(powers[i] > max_power*m_minCodePower)
      {
        if(m_walshCodes[i][walsh_index] == 0)                   // Note: use (0,1)->(1,-1) mapping
        {
          sum_i         += m_integrateIWalsh[i];
          sum_q         += m_integrateQWalsh[i];
        }
        else
        {
          sum_i         -= m_integrateIWalsh[i];
          sum_q         -= m_integrateQWalsh[i];
        }
      }
    }
    m_idealIWalsh[n]    = sum_i/m_walshLength;                  // scale the sum, and save it
    m_idealQWalsh[n]    = sum_q/m_walshLength;
    if(++chip_counter == m_samplesPerChip)
    {
      walsh_index++;
      chip_counter      = 0;
    }
  }     // End loop over Walsh period
  delete [] powers;
  return;
}

#define PILOT_ONLY_ROTATION_ANGLE       0
// ############################ Public Function ###################################
// calculateRotationAngle - Calculate rotation angle of processed input signal for
//                          EVM calculation.
//
// Input:                       None
//
// Output:                      Rotation angle in radians
//
// Notes:
// 1. You must first call processComplexWalshPeriod() and processCrossWalshPeriod()
//    to find the average Walsh values.
// 2. This function only works for IS-95 where I data = Q data
// ############################ Public Function ###################################
double CodeDomain::calculateRotationAngle()
{
  int           i, count;
  float         max_power;
  float         *powers;
  double        cos_theta, sin_theta, theta, sum;
//
// Find maximum Walsh channel power, so that we can ignore low power
// codes.  These powers are already calculated in generateWalshPeriod,
// so we could use the results from there.
//
  powers        = new float[m_walshLength];
  powers[0]     = (m_integrateIWalsh[0]+m_integrateQWalsh[0]) * (m_integrateIWalsh[0]+m_integrateQWalsh[0]);
  max_power     = powers[0];
  for(i=1; i<m_walshLength; i++)
  {
    powers[i]   = (m_integrateIWalsh[i]+m_integrateQWalsh[i]) * (m_integrateIWalsh[i]+m_integrateQWalsh[i]);
    max_power   = MAX(powers[i], max_power);
  }
//
// Find the rotation angle by taking the arctangent of the sine over the cosine.
//
  sum   = 0.;
  count = 0;
  for(i=0; i<m_walshLength; i++)                        // Loop over Walsh codes
  {
#if PILOT_ONLY_ROTATION_ANGLE
    if(i==0)
#else
    if(powers[i] > max_power*m_minCodePower)
#endif
    {
      cos_theta         = (m_integrateIWalsh[i]+m_integrateQWalsh[i]);
      sin_theta         = (m_crossQWalsh[i] - m_crossIWalsh[i]);
      if(cos_theta > 0.)
      {
        sum             += sin_theta/cos_theta;
        count++;
      }
    }
  }     // End loop over Walsh period
  delete [] powers;
  if(count > 0)
    theta       = atan(sum/count);
  else
    theta       = 0.;
  return theta;
}

// ############################ Public Function ###################################
// accumulatePhaseJitter - Accumulate phase jitter statistics over an input Walsh
//                         period.
//
// Input:       inputI:         input real array
//              inputQ:         input imag array
//
// Output:                      None
//
// Notes:
// 1. The following routines must be called before calling this function:
//    * processComplexWalshPeriod()
//    * generateWalshPeriod()
// 2. The phase jitter statistics are accumulated in m_meanErrorAngle and
//    m_varianceErrorAngle.
// ############################ Public Function ###################################
void CodeDomain::accumulatePhaseJitter(const float *inputI, const float *inputQ)
{
  int           n, number_samples;
  double        phi1, phi2, error_angle;
//
// update counter:
//
  m_accumulatedWalshPeriods++;
//
// Calculate the angle to ideal and actual constellation points:
//
  number_samples        = m_walshLength*m_samplesPerChip;
  for(n=0; n<number_samples; n++)
  {
    phi1                        = atan2(m_idealQWalsh[n], m_idealIWalsh[n]);
    phi2                        = atan2(inputQ[n], inputI[n]);
    error_angle                 = phi1 - phi2;
    m_meanErrorAngle            += error_angle;
    m_varianceErrorAngle        += error_angle*error_angle;
  }
  return;
}

// ############################ Public Function ###################################
// calculateEVM - Calculate EVM based on data in input arrays and pre-calculated
//                data over the Walsh period.
//
// Input:       inputI:         input real array
//              inputQ:         input imag array
//
// Output:                      The EVM value over the Walsh period
//
// Notes:
// 1. The following routines must be called before calling this function:
//    * processComplexWalshPeriod()
//    * generateWalshPeriod()
// 2. This function is based on Dave Whipple's submission to TR45.5
// 3. This calculation does not account for carrier feedthrough
// ############################ Public Function ###################################
float CodeDomain::calculateEVM(const float *inputI, const float *inputQ)
{
  int           n, number_samples;
  float         float_samples;
  float         error_i, error_q, evm;
  double        error_power, reference_power;

  number_samples        = m_walshLength*m_samplesPerChip;
  float_samples         = (float)number_samples;
//
// Calculate the error power and the reference power:
//
  error_power           = 0.;
  reference_power       = 0.;
  for(n=0; n<number_samples; n++)
  {
    error_i             = inputI[n] - m_idealIWalsh[n];
    error_q             = inputQ[n] - m_idealQWalsh[n];
    error_power         += error_i*error_i + error_q*error_q;
    reference_power     += m_idealIWalsh[n]*m_idealIWalsh[n] + m_idealQWalsh[n]*m_idealQWalsh[n];
  }
//
// EVM is the square root of error power over reference * 100.
//
  if(reference_power > 0.)
    evm                 = 100.*sqrt(error_power/reference_power);
  else
    evm                 = 0.;
  return evm;
}

// ############################ Public Function ###################################
// calculateEVM - Calculate EVM based on data in input arrays and pre-calculated
//                data over the Walsh period.
//
// Input:       inputI:         input real array
//              inputQ:         input imag array
//              theta:          rotation angle in radians
//
// Output:                      The EVM value over the Walsh period
//
// Notes:
// 1. The following routines must be called before calling this function:
//    * processComplexWalshPeriod()
//    * generateWalshPeriod()
// 2. This function is based on Dave Whipple's submission to TR45.5
// 3. This calculation does not account for carrier feedthrough
// ############################ Public Function ###################################
float CodeDomain::calculateEVM(const float *inputI, const float *inputQ, double theta)
{
  int           n, number_samples;
  float         float_samples;
  float         error_i, error_q, evm;
  float         derotate_i, derotate_q;
  double        error_power, reference_power;
  float         cos_theta, sin_theta;

  number_samples        = m_walshLength*m_samplesPerChip;
  float_samples         = (float)number_samples;
//
// Calculate the error power and the reference power:
//
  error_power           = 0.;
  reference_power       = 0.;
  cos_theta             = cos(theta);
  sin_theta             = sin(theta);
  for(n=0; n<number_samples; n++)
  {
    derotate_i          = inputI[n]*cos_theta + inputQ[n]*sin_theta;    // mult by exp(-j*theta);
    derotate_q          = -inputI[n]*sin_theta + inputQ[n]*cos_theta;
    error_i             = derotate_i - m_idealIWalsh[n];
    error_q             = derotate_q - m_idealQWalsh[n];
    error_power         += error_i*error_i + error_q*error_q;
    reference_power     += m_idealIWalsh[n]*m_idealIWalsh[n] + m_idealQWalsh[n]*m_idealQWalsh[n];
  }
//
// EVM is the square root of error power over reference * 100.
//
  if(reference_power > 0.)
    evm                 = 100.*sqrt(error_power/reference_power);
  else
    evm                 = 0.;
  return evm;
}

// ############################ Public Function ###################################
// calculateRho - Calculate rho based on data in input arrays and pre-calculated
//                data over the Walsh period.
//
// Input:       inputI:         input real array
//              inputQ:         input imag array
//              iPN:            input I channel PN data
//              qPN:            input Q channel PN data
//
// Output:                      The rho value over the Walsh period
//
// Notes:
// 1. The input PN arrays are assumed to have values of +/- 1
// ############################ Public Function ###################################
float CodeDomain::calculateRho(const float *inputI, const float *inputQ, const char *iPN, const char *qPN)
{
  int           n, number_samples;
  int           index, pn_increment;
  float         rho;
  double        i_correlation_sum;
  double        q_correlation_sum;
  double        norm_sum, ideal_norm;

  number_samples        = m_walshLength*m_samplesPerChip;
//
// Initialize parameters:
//
  pn_increment          = index = 0;
  i_correlation_sum     = 0.;
  q_correlation_sum     = 0.;
  norm_sum              = 0.;
  for(n=0; n<number_samples; n++)
  {
    if(iPN[index] == 1)
    {
      i_correlation_sum += inputI[n];
      q_correlation_sum += inputQ[n];
    }
    else
    {
      i_correlation_sum -= inputI[n];
      q_correlation_sum -= inputQ[n];
    }
    if(qPN[index] == 1)
    {
      i_correlation_sum += inputQ[n];
      q_correlation_sum -= inputI[n];
    }
    else
    {
      i_correlation_sum -= inputQ[n];
      q_correlation_sum += inputI[n];
    }
    norm_sum            += inputI[n]*inputI[n] + inputQ[n]*inputQ[n];
    pn_increment++;
    if(pn_increment==(m_samplesPerChip))
    {
      index++;
      pn_increment      = 0;
    }
  }
//
// Rho is the correlation term over the norm term.
//
  ideal_norm            = m_walshLength*2.;                             // 64 * R^2
  if(norm_sum > 0.)
    rho                 = (i_correlation_sum*i_correlation_sum + q_correlation_sum*q_correlation_sum)/
                          norm_sum/ideal_norm;

  else
    rho                 = 0.;
  return rho;
}


// ############################ Public Function ###################################
// calculatePhaseJitter - Calculate phase jitter based on accumulated statistics.
//
// Input:                       None
//
// Output:                      The phase jitter in degrees RMS
//
// Notes:
// 1. The following routines must be called before calling this function:
//    * accumulatePhaseJitter()
// ############################ Public Function ###################################
float CodeDomain::calculatePhaseJitter()
{
  int           number_samples, samples_per_period;
  float         phase_jitter;
//
// Scale the accumulated statistics by the # of samples:
//
  samples_per_period    = m_walshLength*m_samplesPerChip;
  number_samples        = samples_per_period*m_accumulatedWalshPeriods;
  if(number_samples > 0)
  {
    m_meanErrorAngle            /= number_samples;
    m_varianceErrorAngle        /= number_samples;
  }
//
// calculate phase jitter:
//
  phase_jitter                  = sqrt(m_varianceErrorAngle - m_meanErrorAngle*m_meanErrorAngle);
  phase_jitter                  *= DEG_RAD;
  return phase_jitter;
}
