#ifndef _ADAPTIVE_NOTCH_H
#define _ADAPTIVE_NOTCH_H    1
/************************************************************************
 *                                                                      *
 * This class implements an adaptive notch filter.  It can be used,     *
 * for example, for detecting sine waves in noise.                      *
 *                                                                      *
 * File:AdaptiveNotch.h                                                 *
 *                                                                      *
 * The transfer function of the filter is:                              *
 *                                                                      *
 *           1 + a*z^(-1) + z^(-2)                                      *
 *    H(z) = --------------------                                       *
 *           1 + rho*a*z^(-1) + rho^2*z^(-2)                            *
 *                                                                      *
 * Where a = -2*cos(omega), with omega equal to the notch frequency,    *
 * and rho = contraction factor to ensure that the poles stay inside    *
 * the unit circle.                                                     *
 * See: Y. Xiao and K. Shida, "New gradient-based algorithms for        *
 * adaptive IIR notch filters, Proc. ICSP '98, pp. 453-456.             *
 *                                                                      *
 * Revision history:                                                    *
 *  1.  10/15/01  - Started.                                            *
 ************************************************************************/

//
// _adaptType
//
#define NOTCH_SG                0               // sign gradient
#define NOTCH_PG                1               // plain gradient
#define NOTCH_LMP               2               // least Mean p-power
#define NOTCH_MNG               3               // memoryless nonlinear gradient
#define NOTCH_SG_MNG            4               // Hybrid
#define NOTCH_LMP_MNG           5               // Hybrid

#define NOTCH_ORDER             2               // 2nd order notch

#define DEFAULT_MU              0.01            // default value of step size
#define DEFAULT_EPSILON         0.02            // default value of epsilon
#define DEFAULT_RHO             0.9             // default value of contaction factor
#define DEFAULT_A_COEFF         1.              // default value of a
#ifndef YES
#define YES             1
#define NO              0
#endif

class AdaptiveNotch
{
private:
  int           _adaptType;                     // Adaptation algorithm

  float         _stepSize;                      // adaptation parameter, mu
  float         _epsilon;                       // adaptation parameter, epsilon
  float         _contractionFactor;             // Stability factor, rho
  float         _aCoefficient;                  // a coefficient
  float         _storedAValue;                  // reset value of a coefficient

  float         _inputs[NOTCH_ORDER];           // past inputs
  float         _outputs[NOTCH_ORDER];          // past outputs
//
// Private functions:
//
  float         filterOutput(float newInput);   // Output of filter
  float         gradientOutput();               // Gradient signal calculation
  void          updateFilter(float output, float gradient, float inputData);
                                                // Update filter coefficients, and past inputs/outputs
// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  AdaptiveNotch(int type=NOTCH_SG, float stepSize=DEFAULT_MU, float epsilon = DEFAULT_EPSILON,
                float contraction=DEFAULT_RHO, float aCoeff=DEFAULT_A_COEFF);
  ~AdaptiveNotch();
  
/********************************
 * Constructors, destructors    *
 ********************************/
  void          resetFilter();

/**********************
 * Set parameters:    *
 **********************/
  void          setAdaptType(int type);
  void          setStepSize(float step);
  void          setEpsilon(float epsilon);
  void          setContractionFactor(float contraction);
  void          setACoefficient(float a);
  void          setNotchFrequency(float omega);

/**********************
 * Get parameters:    *
 **********************/
  int           adaptType()                     {return _adaptType;             }
  float         stepSize()                      {return _stepSize;              }
  float         contractionFactor()             {return _contractionFactor;     }
  float         aCoefficient()                  {return _aCoefficient;          }
  float         notchFrequency();                       // Calculate notch frequency

/************************
 * Calculating filters: *
 ************************/
  float         processInput(float inputData);          // process new sample through filter

};
#endif
