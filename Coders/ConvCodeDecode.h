#ifndef _CONV_CODE_DECODE_H
#define _CONV_CODE_DECODE_H 1
/***************************************************************************//*!
**  \file         ConvCodeDecode.h
**  \brief        Abstract base class
**  \description  abstract base class for the convolutional
**                coder and decoder classes.
**  \project      EE 598 class
**  \author       Mark Frank
**  \Notes        The convolutional coders and decoders are defined by the output path
**                polynomials in octal.  So, for example, for the GSM Sync channel,
**                G0 = 1 + D^3 + D^4, so that the octal rep of the polynomial would
**                be, 23. Similarly for G1 = 1 + D + D^3 + D^4 -> 33 (octal)
********************************************************************************/
#define MAX_CODER_OUTPUTS       3                       //!<  Max value for k
#define MAX_CODER_INPUTS        1                       //!<  Max value for n
#define MAX_CODER_TRANSITIONS   4                       //!<  2**n
#define MAX_CODER_STAGES        32                      //!<  32 bit integers

class ConvCodeDecode
{
protected:
  int   m_numberOutputs;                         //!<  number of outputs, n
  int   m_numberInputs;                          //!<  number of inputs, k
  int   m_numberStages;                          //!<  memory size, m
  int   m_shiftPtr;                              //!<  Shift register pointer
  int   m_puncturePtr;                           //!<  puncture pointer
  int   m_punctureSize;                          //!<  puncture size

  int   m_g[MAX_CODER_OUTPUTS];                  //!<  generator polynomials in octal
  int   m_p[MAX_CODER_OUTPUTS];                  //!<  puncture polynomials in octal

  int   m_gBinary[MAX_CODER_OUTPUTS];            //!<  generator polynomials in binary
  int   m_pBinary[MAX_CODER_OUTPUTS];            //!<  puncture polynomials in binary


//
// Private functions:
//
  void  octalPolynomialsToBinary();             //!< Convert the polynomials in octal to binary representation

//
// Public functions:
//
public:
    
/********************************
 * Constructors, destructors    *
 ********************************/
  ConvCodeDecode(const int numberOutputs=2,            //!< Default constructor
                 const int numberInputs=1,
                 const int numberStages=4,
                 const int g1=25,
                 const int g2=37,
                 const int g3=0,
                 const int p1=1,
                 const int p2=17,
                 const int p3=0);
  virtual ~ConvCodeDecode();                     //!< Default destructor

/********************************
 * Initializing coder:          *
 ********************************/
  virtual       void    initCoder();             //!< Initialize coder before operation

/********************************
 * Setting coder parameters:    *
 ********************************/
  void  setNumberOutputs(const int n);                 //!< Set the number of output paths (polynomials)
  void  setNumberInputs(const int k);                  //!< Set the number of input channels
  void  setNumberStages(const int m);                  //!< Set the number of shift register stages
  void  setGPolynomial(const int gOctal,               //!< Set the polynomial for the number'th output
                       const int number);
  void  setPPolynomial(const int pOctal,               //!< Set the puncture polynomial for the number'th output
                       const int number);


/********************************
 * Getting coder parameters:    *
 ********************************/
  int   numberOutputs() const {return m_numberOutputs; }
  int   numberInputs()  const {return m_numberInputs;  }
  int   numberStages()  const {return m_numberStages;  }
  int   g1Octal()       const {return m_g[0];          }
  int   g2Octal()       const {return m_g[1];          }
  int   g3Octal()       const {return m_g[2];          }
  int   p1Octal()       const {return m_p[0];          }
  int   p2Octal()       const {return m_p[1];          }
  int   p3Octal()       const {return m_p[2];          }
  int   numberStates()  const;

  private:
  // These are not implemented, just listed here so they are not implicitly used:
  // lint !e1704
  ConvCodeDecode & operator = (const ConvCodeDecode & rhs);       //!< Disabled assignment operator
  // lint !e1704
  ConvCodeDecode (const ConvCodeDecode & rhs);                    //!< Disabled copy constructor

};
#endif



