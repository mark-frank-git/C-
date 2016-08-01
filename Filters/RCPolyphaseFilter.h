#ifndef _RCPOLYPHASEFILTER_H
#define _RCPOLYPHASEFILTER_H    1
/************************************************************************
 *                                                                      *
 * This subclass of ComplexDigitalFilter adds functionality for         *
 * calculating the transfer function of RC polyphase filters.           *
 * The actual filtering functions are defined in the superclass.        *
 *                                                                      *
 * File:RCPolyphaseFilter.h                                             *
 *                                                                      *
 * The filter is stored in the following form:                          *
 *                                                                      *
 *           b[0] + b[1]z^(-1) + ... + b[n]z^-(n)                       *
 *  H(z)   = ------------------------------------                       *
 *           1    + a[1]z^(-1) + ... + a[n]z^-(n)                       *
 * Where the a and b polynomials are stored in the polynomial objects.  *
 * NOTE: The order of the filter is the # taps - 1.                     *
 *                                                                      *
 * SEE: "RC Sequence Asymmetric Polyphase Networks fro RF Integrated    *
 *       Transceivers," by Galal, Ragaie, Tawfik.                       *
 *                                                                      *
 ************************************************************************/


#ifndef YES
#define YES     1
#define NO      0
#endif

#define DEFAULT_STAGES          2
#define DEFAULT_R               16.e3
#define DEFAULT_C               2.43e-12
#define DEFAULT_SAMPLING_FREQ   16.368e6        // 16fo

class Complex;                                  // class prototype
class ComplexDigitalFilter;

class RCPolyphaseFilter: public ComplexDigitalFilter
{
private:
  int           _numberStages;                  // Number of RC stages
  float         _r1, _r2, _r3;                  // The resistor values in ohms
  float         _c1, _c2, _c3;                  // The capacity values in Farads

//
// Private functions:
//
  void          findSTransferFunction();        // Find the transfer function in the S domain

// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  RCPolyphaseFilter(int numberStages=DEFAULT_STAGES);
  ~RCPolyphaseFilter();

  /**********************
 * Set parameters:    *
 **********************/
  void          setNumberStages(int stages)     {_numberStages=stages;  return;}
  void          setR1(float r1)                 {_r1 = r1;              return;}
  void          setR2(float r2)                 {_r2 = r2;              return;}
  void          setR3(float r3)                 {_r3 = r3;              return;}
  void          setC1(float c1)                 {_c1 = c1;              return;}
  void          setC2(float c2)                 {_c2 = c2;              return;}
  void          setC3(float c3)                 {_c3 = c3;              return;}

/**********************
 * Get parameters:    *
 **********************/
  float         r1()                    {return _r1;}
  float         r2()                    {return _r2;}
  float         r3()                    {return _r3;}
  float         c1()                    {return _c1;}
  float         c2()                    {return _c2;}
  float         c3()                    {return _c3;}

/************************
 * Calculating filters: *
 ************************/
  void          findTransferFunction();                 // Find transfer function for RC polyphase

};
#endif
