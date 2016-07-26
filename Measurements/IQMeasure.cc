/************************************************************************
 *                                                                      *
 * This is used for making measurements on QAM I and Q symbols.  The    *
 * measurements are geared for efficiency in that variables common to   *
 * multiple measurements should only be calculated once.                *
 *                                                                      *
 *                                                                      *
 * File: /User/frank/C++/QAM/IQMeasure.cc                               *
 *                                                                      *
 ************************************************************************/

#include "IQMeasure.h"
#include <stdio.h>
#include <math.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ROUND(a) ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )

//
// The offset for converting from index to I and Q values:
//
int     qam_offset[MAX_BITS+1]  = {0, 0, 1, 1, 3, 5, 7, 9, 15};



// ############################ Private Function ###################################
// updateTargetStatistics - Update the target statistics based on the accumulated
//                          counters.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// 1.  The instance variables, m_meanTargetI, m_meanTargetQ, m_targetVarianceI, etc. are updated.
// ############################ Private Function ###################################
void IQMeasure::updateTargetStatistics()
{
  int       i;
  double    dbl_points;
//
// Normalize the mean target values by number of points collected:
//
  for(i=0; i<m_qamPoints; i++)
  {
    if( m_numberTargetPoints[i] > 1)
    {
      dbl_points                = (double)m_numberTargetPoints[i];
      m_meanTargetI[i]          /= dbl_points;
      m_meanTargetQ[i]          /= dbl_points;
      m_meanTargetErrorI[i]     /= dbl_points;
      m_meanTargetErrorQ[i]     /= dbl_points;
      m_targetVarianceI[i]      = (m_targetVarianceI[i] - m_meanTargetErrorI[i]*m_meanTargetErrorI[i]*dbl_points)/
                                  (dbl_points-1.);
      m_targetVarianceQ[i]      = (m_targetVarianceQ[i] - m_meanTargetErrorQ[i]*m_meanTargetErrorQ[i]*dbl_points)/
                                  (dbl_points-1.);
    }
  }
  return;
}

#define INFINITE_SLOPE  1.e10
// ############################ Private Function ###################################
// leastSquaresFit - Do a linear least squares fit to the input data
//
// Input:       xData:              x data points array
//              yData:              y data points array
//              numberDataPoints:   size of the arrays
//
// Output:      lineSlope:          slope of the line
//              lineIntercept:      intercept of the line
//
// Notes:
// 1. From Forsythe, Malcolm, and Moler, "Computer methods for mathematical
//    computations."
// ############################ Private Function ###################################
void IQMeasure::leastSquaresFit(double *xData, double *yData, int numberDataPoints,
                                double *lineSlope, double *lineIntercept)
{
  int   i;
  double p[2][2];               // least squares matrix;
  double q[2], determinant;

//
// Error check:
//
  if(numberDataPoints < 3)
    return;
//
// Find matrix elements:
//
  p[0][0]   = numberDataPoints;
  p[0][1]   = p[1][0]   = p[1][1]   = 0.;
  q[0]      = q[1]      = 0.;
  for(i=0; i<numberDataPoints; i++)
  {
    p[0][1] += xData[i];
    p[1][1] += xData[i]*xData[i];
    q[0]    += yData[i];
    q[1]    += xData[i]*yData[i];
  }
  p[1][0]   = p[0][1];                  // symmetric matrix
//
// Find matrix inverse
//
  determinant   = p[0][0]*p[1][1] - p[0][1]*p[1][0];
  if(determinant == 0.)
  {
    *lineSlope      = INFINITE_SLOPE;
    *lineIntercept  = xData[0];
    return;
  }
//
// Solve for slope and intercept:
//
  *lineIntercept    = (q[0]*p[1][1] - q[1]*p[0][1])/determinant;
  *lineSlope        = (q[1]*p[0][0] - q[0]*p[1][0])/determinant;

  return;
}

// ############################ Private Function ###################################
// findImbalance - Find the current value of amplitude imbalance based on accumulated counters.
//
// Input:                       None
//
// Output:                      amplitude imbalance
//
// Notes:
// 1.  updateTargetStatistics must be called before this function.
// ############################ Private Function ###################################
double IQMeasure::findImbalance()
{
  int       i;
  double    vi, vq, v1, v2;
  double    amplitude_imbalance;
//
// Calculate amplitude imbalance:
//
  vi            = vq        = 0.;
  for(i=0; i<m_qamPoints; i++)
  {
    if( (m_idealTargetI[i]!=0.) && (m_idealTargetQ[i]!=0.) )
    {
      vi                    += (m_idealTargetI[i] + m_meanTargetErrorI[i])/m_idealTargetI[i];
      vq                    += (m_idealTargetQ[i] + m_meanTargetErrorQ[i])/m_idealTargetQ[i];
    }
  }
  vi                        /= m_qamPoints;
  vq                        /= m_qamPoints;
  v1                        = MIN(vi, vq);
  v2                        = MAX(vi, vq);

  if(v1 > 0.)
    amplitude_imbalance     = 100.*(v2/v1 - 1.);
  else
    amplitude_imbalance     = 0.;
  return amplitude_imbalance;
}


// ############################ Private Function ###################################
// findQuadError - Find the current value of quadrature error based on accumulated counters.
//
// Input:                       None
//
// Output:                      output quadrature error
//
// Notes:
// 1.  updateTargetStatistics must be called before this function.
// ############################ Private Function ###################################
double IQMeasure::findQuadError()
{
  int           i, j, k, l, number_data_pts;
  double        avg_i_slope, avg_q_slope, phi_1, phi_2, slope;
  double        quad_error;
//
// Quadrature Error:
//
#if LEAST_SQUARES_FIT
  double    *x, *y;
  x     = new double[m_qamLevels];
  y     = new double[m_qamLevels];
  for(i=0; i<m_qamLevels; i++)            // Calculate curve fit for lines parallel to I axis
  {
    k               = i*m_qamLevels;
    number_data_pts = 0;
    for(j=0; j<m_qamLevels; j++)
    {
      if( m_numberTargetPoints[k] > 0)
      {
        x[number_data_pts]  = m_meanTargetI[k];
        y[number_data_pts]  = m_meanTargetQ[k];
        number_data_pts++;
      }
      k++;
    }
    leastSquaresFit(x, y, number_data_pts, &m_slopeI[i], &m_interceptI[i]);
  }

  for(i=0; i<m_qamLevels; i++)            // Calculate curve fit lines parallel to Q axis
  {
    k               = i;
    number_data_pts = 0;
    for(j=0; j<m_qamLevels; j++)
    {
      if( m_numberTargetPoints[k] > 0)
      {
        x[number_data_pts]  = m_meanTargetI[k];
        y[number_data_pts]  = m_meanTargetQ[k];
        number_data_pts++;
      }
      k                     += m_qamLevels;
    }
    leastSquaresFit(x, y, number_data_pts, &m_slopeQ[i], &m_interceptQ[i]);
  }
  delete [] x;
  delete [] y;
#else
  for(i=0; i<m_qamLevels; i++)            // Calculate slopes for lines parallel to I axis
  {
    l               = k             = i*m_qamLevels;
    m_slopeI[i]     = 0.;
    number_data_pts = 0;
    for(j=1; j<m_qamLevels; j++)
    {
      k++;
      if( j == (m_qamLevels-1) )
      {
        slope   = m_meanTargetQ[k] - m_meanTargetQ[l];
        if( m_meanTargetI[k] != m_meanTargetI[l])
        {
          slope /= m_meanTargetI[k] - m_meanTargetI[l];
          number_data_pts++;
          m_slopeI[i]   += slope;
        }
      }
    }
    if(number_data_pts > 0)
      m_slopeI[i]       /= number_data_pts;
  }

  for(i=0; i<m_qamLevels; i++)            // Calculate slopes for lines parallel to Q axis
  {
    k               = l = i;
    m_slopeQ[i]     = 0.;
    number_data_pts = 0;
    for(j=1; j<m_qamLevels; j++)
    {
      k             += m_qamLevels;
      if( j == (m_qamLevels-1) )
      {
        slope       = m_meanTargetI[k] - m_meanTargetI[l];
        if( m_meanTargetQ[k] != m_meanTargetQ[l])
        {
          slope /= m_meanTargetQ[k] - m_meanTargetQ[l];
          number_data_pts++;
          m_slopeQ[i]   += slope;
        }
      }
    }
    if(number_data_pts > 0)
      m_slopeQ[i]       /= number_data_pts;
  }
#endif
//
// Average the slopes:
//
  avg_i_slope   = 0.;
  avg_q_slope   = 0.;
  for(i=0; i<m_qamLevels; i++)
  {
    avg_i_slope += m_slopeI[i];
    avg_q_slope += m_slopeQ[i];
  }
  avg_i_slope   /= m_qamLevels;
  avg_q_slope   /= m_qamLevels;
//
// Calculate angles from slopes, these angles may not be the same
// as in the DVB document, but they should work
//
#if LEAST_SQUARES_FIT
  if(avg_i_slope > 0.)
  {
    phi_1       = atan(avg_q_slope);
    phi_2       = atan(avg_i_slope);
    quad_error  = 90. + RAD_TO_DEG*(phi_2 - phi_1);
  }
  else if(avg_i_slope < 0.)
  {
    if(avg_q_slope != 0.)
      phi_1     = atan(1./avg_q_slope);
    else
      phi_1     = 0.;
    phi_2       = atan(avg_i_slope);
    quad_error  = - RAD_TO_DEG*(phi_1 + phi_2);
  }
#else
  phi_1         = atan(avg_i_slope);
  phi_2         = atan(avg_q_slope);
  quad_error    = RAD_TO_DEG*(phi_1 + phi_2);
#endif

  return quad_error;
}

#define EQUAL_CHECK     0.017                           // This corresponds to about 1 degree difference
// ############################ Private Function ###################################
// findNormalizedMERAndPJ - Uncouple the MER and phase jitter measurements.
//
// Input:                       None
//
// Output:              mer:    MER minus the phase jitter
//                      pj:     phase jitter minus the 'MER' component
//
// Notes:
// 1. updateTargetStatistics must be called before this function.
// 2. m_meanSigPower must be calculated before calling this function.
// 3. See the write up in "Reducing the impact of phase noise on Magnet measurements".
// ############################ Private Function ###################################
void IQMeasure::findNormalizedMERAndPJ(double *mer, double *pj)
{
  int           i, number_avgs;
  double        r_sq, i_sq, q_sq, cos_theta_sq, sin_theta_sq;
  double        sigma_i_sq, sigma_q_sq, sigma_phi_sq, two_sigma_n_sq;
  double        mean_pj, mean_noise;
//
// For each constellation point we assume that 
// sigma_I^2 = sigma_I_phi^2 + sigma_n^2, and
// sigma_Q^2 = sigma_Q_phi^2 + sigma_n^2
// where: sigma_I_phi = r sin(theta) sigma_phi; constellation point = (r, theta)
//
  mean_pj       = mean_noise    = 0.;
  number_avgs   = 0;
  for(i=0; i<m_qamPoints; i++)
  {
    if(m_targetVarianceI[i] > 0.)                       // have statistics been collected?
    {
      sigma_i_sq        = m_targetVarianceI[i];
      sigma_q_sq        = m_targetVarianceQ[i];         // previously calculated in updateTargetStatistics()
//
// Find polar form of constellation point:
//
      i_sq              = m_idealTargetI[i]*m_idealTargetI[i];
      q_sq              = m_idealTargetQ[i]*m_idealTargetQ[i];
      r_sq              = i_sq + q_sq;
      if(r_sq > 0.)                                     // This should always be true
      {
        cos_theta_sq    = q_sq/r_sq;
        sin_theta_sq    = i_sq/r_sq;
      }
      else
        cos_theta_sq    = sin_theta_sq  = 0.;
//
// Now, solve the two equations in two unknowns:
//
      if(fabs(cos_theta_sq-sin_theta_sq) > EQUAL_CHECK) // We can't use this method at the 45 degree points
      {
        sigma_phi_sq    = (sigma_i_sq - sigma_q_sq)/r_sq/(cos_theta_sq - sin_theta_sq);
        two_sigma_n_sq  = sigma_i_sq + sigma_q_sq - r_sq*sigma_phi_sq;
//
// Update the averages, weight by the number of measurements at each constellation point
//
        mean_pj         += sigma_phi_sq*m_numberTargetPoints[i];
        mean_noise      += two_sigma_n_sq*m_numberTargetPoints[i];
        number_avgs     += m_numberTargetPoints[i];
      }
    }
  }
  if(number_avgs > 0)
  {
    if(mean_pj > 0.)
      *pj       = sqrt(mean_pj/number_avgs);
    mean_noise  /= number_avgs;
    *mer        = m_meanSigPower/mean_noise;
  }
  return;
}

// ############################ Private Function ###################################
// processIQSymbolAtIndex - This method adds a new I/Q symbol to the statistics counters.
//
// Input:       iData:          I (in-phase) value of the input symbol
//              qData:          Q (quadrature) value of the input symbol
//              iIdeal:         I value of the ideal constellation point
//              qIdeal:         Q value of the ideal constellation point
//              qamIndex:       index to store measurement values
//
// Output:                      None
//
// Notes:
// ############################ Private Function ###################################
void IQMeasure::processIQSymbolAtIndex(double iData, double qData, double iIdeal, double qIdeal, int qamIndex)
{
  double        i_error, q_error;
  double        phi1, phi2, error_angle;
  double        var_i, var_q;
//
// Store the target points:
//
  m_idealTargetI[qamIndex]
                        = iIdeal;
  m_idealTargetQ[qamIndex]
                        = qIdeal;
  m_meanTargetI[qamIndex]
                        += iData;
  m_meanTargetQ[qamIndex]
                        += qData;

//
// Calculate modulation error
//
  i_error               = iIdeal - iData;
  q_error               = qIdeal - qData;
  var_i                 = i_error*i_error;
  var_q                 = q_error*q_error;
  m_meanI               += i_error;
  m_meanQ               += q_error;
  m_modulationError     += var_i + var_q;
//
// Store the target error for amp imbalance:
//
  m_meanTargetErrorI[qamIndex]
                        += -i_error;
  m_meanTargetErrorQ[qamIndex]
                        += -q_error;
//
// Store the target variance for normalized phase jitter measurement:
//
  m_targetVarianceI[qamIndex]
                        += var_i;
  m_targetVarianceQ[qamIndex]
                        += var_q;

//
// Calculate the error angles for phase jitter:
//
  phi1                  = atan2(qIdeal, iIdeal);
  phi2                  = atan2(qData, iData);
  error_angle           = phi1 - phi2;
  m_meanErrorAngle      += error_angle;
  m_varianceErrorAngle  += error_angle*error_angle;
//
// Update counters:
//
  m_numberTargetPoints[qamIndex]++;
  m_numberMeanPoints++;
  return;
}

// ############################ Private Function ###################################
// allocateArrays - This method should be called whenever the number of constellation
//                  points has changed.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################ Private Function ###################################
void IQMeasure::allocateArrays()
{
//
// delete arrays:
//
  delete [] m_numberTargetPoints;
  delete [] m_meanTargetI;
  delete [] m_meanTargetQ;
  delete [] m_meanTargetErrorI;
  delete [] m_meanTargetErrorQ;
  delete [] m_targetVarianceI;
  delete [] m_targetVarianceQ;
  delete [] m_idealTargetI;
  delete [] m_idealTargetQ;
  delete [] m_slopeI;
  delete [] m_interceptI;
  delete [] m_slopeQ;
  delete [] m_interceptQ;
  
//
// Allocate arrays:
//
  m_numberTargetPoints                  = new int[m_qamPoints];
  m_meanTargetI                         = new double[m_qamPoints];
  m_meanTargetQ                         = new double[m_qamPoints];
  m_meanTargetErrorI                    = new double[m_qamPoints];
  m_meanTargetErrorQ                    = new double[m_qamPoints];
  m_targetVarianceI                     = new double[m_qamPoints];
  m_targetVarianceQ                     = new double[m_qamPoints];
  m_idealTargetI                        = new double[m_qamPoints];
  m_idealTargetQ                        = new double[m_qamPoints];
  m_slopeI                              = new double[m_qamLevels];
  m_interceptI                          = new double[m_qamLevels];
  m_slopeQ                              = new double[m_qamLevels];
  m_interceptQ                          = new double[m_qamLevels];

  return;
}


// ############################# Class Constructor #################################
// IQMeasure -- Constructor for the IQMeasure class
//
// Input:       bits:   # of bits per constellation point
//              scale:  qam scale factor
//
// Output:              None
// ############################# Class Constructor #################################
IQMeasure::IQMeasure(int bits, float scale)
{

// 
// Initialize instance variables:
//
  m_numberTargetPoints                  = NULL;
  m_meanTargetI                         = NULL;
  m_meanTargetQ                         = NULL;
  m_meanTargetErrorI                    = NULL;
  m_meanTargetErrorQ                    = NULL;
  m_targetVarianceI                     = NULL;
  m_targetVarianceQ                     = NULL;
  m_idealTargetI                        = NULL;
  m_idealTargetQ                        = NULL;
  m_slopeI                              = NULL;
  m_interceptI                          = NULL;
  m_slopeQ                              = NULL;
  m_interceptQ                          = NULL;
  
  setQAMBits(bits);
  setQAMMinScale(scale);
  m_maxConstellationPt                  = -1.;
//
// Allocate statistics objects:
//

  m_merStatistics                       = new SampleStatistic();
  m_evmStatistics                       = new SampleStatistic();
  m_suppressionStatistics               = new SampleStatistic();
  m_imbalanceStatistics                 = new SampleStatistic();
  m_quadErrorStatistics                 = new SampleStatistic();
  m_phaseJitterStatistics               = new SampleStatistic();
  m_normalizedMERStatistics             = new SampleStatistic();
  m_normalizedPhaseJitterStatistics     = new SampleStatistic();

// 
// Reset all counters and statistics:
//
  resetCounters();
  resetStatistics();

  return;
}


// ############################# Class Destructor ###############################
// IQMeasure -- Destructor for the IQMeasure class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
IQMeasure::~IQMeasure()
{
//
// delete arrays:
//
  delete [] m_numberTargetPoints;
  delete [] m_meanTargetI;
  delete [] m_meanTargetQ;
  delete [] m_meanTargetErrorI;
  delete [] m_meanTargetErrorQ;
  delete [] m_targetVarianceI;
  delete [] m_targetVarianceQ;
  delete [] m_idealTargetI;
  delete [] m_idealTargetQ;
  delete []  m_slopeI;
  delete [] m_interceptI;
  delete [] m_slopeQ;
  delete [] m_interceptQ;

  delete [] m_merStatistics;
  delete [] m_evmStatistics;
  delete [] m_suppressionStatistics;
  delete [] m_imbalanceStatistics;
  delete [] m_quadErrorStatistics;
  delete [] m_phaseJitterStatistics;
  delete [] m_normalizedMERStatistics;
  delete [] m_normalizedPhaseJitterStatistics;
  
  return;
}

// ############################ Public Function ###################################
// resetCounters - Resets the counters between trials.
//
// Input:                       None
// Output:                      None
//
// Notes:
// 1. This function is also called from resetStatistics().
// ############################ Public Function ###################################
void IQMeasure::resetCounters()
{
  int   i;
  
  m_numberMeanPoints    = 0;
  m_numberPowerPoints   = 0;
  m_meanI               = m_meanQ   = 0.;
  m_modulationError     = 0.;
  m_currentMeanPower    = 0.;
  m_meanErrorAngle      = 0.;
  m_varianceErrorAngle  = 0.;

  for(i=0; i<m_qamPoints; i++)
  {
    m_numberTargetPoints[i]     = 0;
    m_meanTargetErrorI[i]       = m_meanTargetErrorQ[i]         = 0.;
    m_meanTargetI[i]            = m_meanTargetQ[i]              = 0.;
    m_idealTargetI[i]           = m_idealTargetQ[i]             = 0.;
    m_targetVarianceI[i]        = m_targetVarianceQ[i]          = 0.;
  }
  return;
}

// ############################ Public Function ###################################
// resetStatistics - Resets the statistics between channel changes.
//
// Input:                       None
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void IQMeasure::resetStatistics()
{
  m_merStatistics->reset();
  m_evmStatistics->reset();
  m_suppressionStatistics->reset();
  m_imbalanceStatistics->reset();
  m_quadErrorStatistics->reset();
  m_phaseJitterStatistics->reset();
  m_normalizedMERStatistics->reset();
  m_normalizedPhaseJitterStatistics->reset();
  resetCounters();
  return;
}

// ############################ Public Function ###################################
// setQAMMinScale - Sets a new constellation half spacing.
//
// Input:       scale:          half spacing between constellation points
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void IQMeasure::setQAMMinScale(float scale)
{
  m_qamMinScale = scale;
//
// Update mean signal power:
//
  m_meanSigPower    = idealSigPower();
  return;
}
  
// ############################ Public Function ###################################
// setQAMBits - Sets the # of bits per constellation point.
//
// Input:       bits:           # of bits
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void IQMeasure::setQAMBits(int bits)
{
  m_qamBits     = MAX(2, bits);
  m_qamBits     = MIN(m_qamBits, MAX_BITS);

  m_qamPoints   = 1 << m_qamBits;
  m_qamLevels   = 1 << (m_qamBits/2);                   // doesn't work for odd _qamBits
  allocateArrays();
//
// Update mean signal power:
//
  m_meanSigPower    = idealSigPower();
  return;
}

// ############################ Public Function ###################################
// processIQSymbol - This method adds a new I/Q symbol to the statistics counters.
//
// Input:       iData:          I (in-phase) value of the input symbol
//              qData:          Q (quadrature) value of the input symbol
//              iIdeal:         I value of the ideal constellation point
//              qIdeal:         Q value of the ideal constellation point
//
// Output:                      None
//
// Notes:
// 1. After calling this function for the desired number of input symbols, call
//    updateStatistics(), to update the statistics for the current trial, and then
//    call one of meanMER(), stdDevMER(), etc. to get the output statistics.
// ############################ Public Function ###################################
void IQMeasure::processIQSymbol(double iData, double qData, double iIdeal, double qIdeal)
{
  int           i_index, q_index, qam_point;
  float         offset;
  double        ideal_mag2;

//
// Update constellation point values:
//
  ideal_mag2            = iIdeal*iIdeal + qIdeal*qIdeal;
  m_currentMeanPower    += ideal_mag2;
  m_numberPowerPoints++;
  m_maxConstellationPt  = MAX(m_maxConstellationPt, ideal_mag2);

//
// Get constellation point index for
// amplitude imbalance, quadrature error:
//
  offset                = qam_offset[m_qamBits];
  i_index               = (ROUND( (iIdeal/m_qamMinScale + offset)) )/2;
  q_index               = (ROUND( (qIdeal/m_qamMinScale + offset)) )/2;
  qam_point             = i_index + q_index*m_qamLevels;
  qam_point             = MAX(0, qam_point);
  qam_point             = MIN(qam_point, (m_qamPoints-1));
//
// Store the target points, and the modulation error
//
  processIQSymbolAtIndex(iData, qData, iIdeal, qIdeal, qam_point);
  return;
}

// ############################ Public Function ###################################
// processIQSymbolInRing - This method adds a new I/Q symbol to the statistics counters,
//                         if the point lies in the input ring.
//
// Input:       iData:          I (in-phase) value of the input symbol
//              qData:          Q (quadrature) value of the input symbol
//              iIdeal:         I value of the ideal constellation point
//              qIdeal:         Q value of the ideal constellation point
//              ring:           1 = inner most ring, 2 = next outer ring, etc.
//
// Output:                      None
//
// Notes:
// 1. After calling this function for the desired number of input symbols, call
//    updateStatistics(), to update the statistics for the current trial, and then
//    call one of meanMER(), stdDevMER(), etc. to get the output statistics.
// ############################ Public Function ###################################
void IQMeasure::processIQSymbolInRing(double iData, double qData, double iIdeal, double qIdeal, int ring)
{
  int           i_index, q_index, qam_point;
  int           ring_low, ring_high;
  float         offset;
  double        ideal_mag2;

//
// Update constellation point values:
//
  ideal_mag2            = iIdeal*iIdeal + qIdeal*qIdeal;
  m_currentMeanPower    += ideal_mag2;
  m_numberPowerPoints++;
  m_maxConstellationPt  = MAX(m_maxConstellationPt, ideal_mag2);

//
// Get constellation point index for
// amplitude imbalance, quadrature error:
//
  offset                = qam_offset[m_qamBits];
  i_index               = (ROUND( (iIdeal/m_qamMinScale + offset)) )/2;
  q_index               = (ROUND( (qIdeal/m_qamMinScale + offset)) )/2;
  qam_point             = i_index + q_index*m_qamLevels;
  qam_point             = MAX(0, qam_point);
  qam_point             = MIN(qam_point, (m_qamPoints-1));
//
// Check if point lies in the ring:
//
  ring_low              = m_qamLevels/2 - 1 - ring;
  ring_high             = m_qamLevels/2 + ring;
  if( (i_index <= ring_low) || (i_index >= ring_high) )
    return;
  if( (q_index <= ring_low) || (q_index >= ring_high) )
    return;
//
// Point lies in ring, so store the target points, and the modulation error
//
  processIQSymbolAtIndex(iData, qData, iIdeal, qIdeal, qam_point);
  return;
}

// ############################ Public Function ###################################
// updateStatistics - Update the statistics based on current accumulated values.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void IQMeasure::updateStatistics()
{
  double    mer, pj, evm;
  double    carrier_suppression;
//
// Process statistics counters:
//
  if(m_numberMeanPoints > 0)
  {
    m_modulationError           /= m_numberMeanPoints;
    m_meanI                     /= m_numberMeanPoints;
    m_meanQ                     /= m_numberMeanPoints;
    m_meanErrorAngle            /= m_numberMeanPoints;
    m_varianceErrorAngle        /= m_numberMeanPoints;
  }
  if(m_numberPowerPoints > 0)
    m_currentMeanPower  /= m_numberPowerPoints;
//
// Process target statistics counters:
//
  updateTargetStatistics();

//
// Update MER:
//
  mer                   = evm   = 0.;                                           // keep compiler happy
  if(m_meanSigPower < 0.)
    m_meanSigPower      = m_currentMeanPower;                                   // First time through
  if(m_modulationError > 0.)
      mer               = m_meanSigPower/m_modulationError;
  if(mer > 0.)
     *m_merStatistics   += 10.*log10(mer);
 //
// EVM:
//
  if(m_maxConstellationPt > 0.)
    evm                 = sqrt(m_modulationError/m_maxConstellationPt);
  *m_evmStatistics      += 100.*evm;

//
// Carrier Suppression:
//
  carrier_suppression           = m_meanI*m_meanI + m_meanQ*m_meanQ;            // Referenced to 0 + j0
  if(m_meanSigPower > 0.)
    carrier_suppression         /= m_meanSigPower;
  if(carrier_suppression > 0.)
    *m_suppressionStatistics    += -10.*log10(carrier_suppression);

//
// Amplitude Imbalance and Quadrature error
//
  *m_imbalanceStatistics        += findImbalance();
  *m_quadErrorStatistics        += findQuadError();

//
// Phase jitter:
//
  *m_phaseJitterStatistics      += sqrt(m_varianceErrorAngle - m_meanErrorAngle*m_meanErrorAngle);

//
// Normalized MER and phase jitter:
//
  findNormalizedMERAndPJ(&mer, &pj);
  if(mer > 0.)
     *m_normalizedMERStatistics         += 10.*log10(mer);
  *m_normalizedPhaseJitterStatistics    += pj;

//
// Update mean signal power:
//
  m_meanSigPower    = FILTER_ALPHA*m_meanSigPower + (1. - FILTER_ALPHA)* m_currentMeanPower;

  return;
}

// ############################ Public Function ###################################
// idealSigPower - This method returns the ideal average constellation power.
//
// Input:                       None
//
// Output:                      average signal power in watts
//
// Notes:
// ############################ Public Function ###################################
double IQMeasure::idealSigPower()
{
  int           i_index, q_index;
  double        ideal_mag2, i_value, q_value, sig_power;

  sig_power     = 0.;
  q_value       = m_qamMinScale;
  for(q_index=0; q_index<m_qamLevels/2; q_index++)
  {
    i_value     = m_qamMinScale;
    for(i_index=0; i_index<m_qamLevels/2; i_index++)
    {
      ideal_mag2        = i_value*i_value + q_value*q_value;
      sig_power         += ideal_mag2;
      i_value           += 2.*m_qamMinScale;
    }
    q_value             += 2.*m_qamMinScale;
  }
  return (4.*sig_power/m_qamPoints);
}


