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
#include "ConvCodeDecode.h"
#include <stdio.h>


#define MAX(a,b)        ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )

#ifndef YES
#define YES 1
#define NO  0
#endif

#define OCTAL_BITS  3


/********************** Private Function **********************************//**
**  \brief              Octal to Binary
**  \description        Converts the polynomials given in octal to binary representations
**  \param[in]          None
**  \return             None
**   \post         
**  <b>Modified:</b>    m_gBinary, m_pBinary, m_punctureSize
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
void ConvCodeDecode::octalPolynomialsToBinary()
{
  int   i, g_poly, p_poly, octal_power, octal_digit;
  int   puncture_bits;

  m_punctureSize = 0;
  for(i=0; i<m_numberOutputs; i++)
  {
//
// Convert octal representation to binary for g polynomials:
//
    g_poly      = m_g[i];
    p_poly      = m_p[i];
    m_gBinary[i] = 0;
    m_pBinary[i] = 0;
    octal_power = 1;
    while(g_poly/10)
    {
      octal_digit       = g_poly % 10;
      m_gBinary[i]       += octal_digit*octal_power;
      octal_power       *= 8;
      g_poly            /= 10;
    }
    m_gBinary[i]         += g_poly*octal_power;
//
// Repeat for puncture polynomial:
//
    puncture_bits       = 0;
    octal_power         = 1;
    while(p_poly/10)
    {
      octal_digit       = p_poly % 10;
      m_pBinary[i]       += octal_digit*octal_power;
      octal_power       *= 8;
      p_poly            /= 10;
      puncture_bits     += OCTAL_BITS;
    }
    m_pBinary[i]         += p_poly*octal_power;
    while(p_poly)
    {
      puncture_bits++;
      p_poly    /= 2;
    }
//
// Find the maximum number of bits in the puncture polynomials:
//
    m_punctureSize       = MAX(m_punctureSize, puncture_bits);
  }
  
  return;
}


/********************** Class Constructor **********************************//**
**  \brief              Constructor
**  \description        Constructor
**  \param[in]          numberOutputs  = number of output bits per input bit
**  \param[in]          numberInputs   = number of input shift registers
**  \param[in]          numberStages   = number of stages in coder
**  \param[in]          g1->g3         = generator polynomials in octal
**  \param[in]          p1->p3         = puncture polynomials in octal
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
ConvCodeDecode::ConvCodeDecode(const int numberOutputs,
                               const int numberInputs,
                               const int numberStages,
                               const int g1,
                               const int g2,
                               const int g3,
                               const int p1,
                               const int p2,
                               const int p3)
{
//
// Initialize instance variables:
//
  setNumberOutputs(numberOutputs);
  setNumberInputs(numberInputs);
  setNumberStages(numberStages);
  setGPolynomial(g1, 0);
  setGPolynomial(g2, 1);
  setGPolynomial(g3, 2);
  setPPolynomial(p1, 0);
  setPPolynomial(p2, 1);
  setPPolynomial(p3, 2);
//
// Initialize coder:
//
  initCoder();
  
  return;
}
  

/********************** Class Destructor **********************************//**
**  \brief              Class Destructor
**  \description        Class Destructor
**  \param[in]          None
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
ConvCodeDecode::~ConvCodeDecode()
{
  return;
}


/********************** Public Function **********************************//**
**  \brief              Initializes the coder to start accepting new input
**  \description        Initializes the coder to start accepting new input
**  \param[in]          None
**  \return             None
**   \post         
**  <b>Modified:</b>    m_shiftPtr, m_puncturePtr, m_gBinary, m_pBinary
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
void ConvCodeDecode::initCoder()
{
  m_shiftPtr     = 0;
  m_puncturePtr  = 0;
  octalPolynomialsToBinary();
  return;
}

/********************** Public Function **********************************//**
**  \brief              Sets a new number of outputs from coder
**  \description        Sets a new number of outputs from coder
**  \param[in]          n          = new number of outputs
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       m_numberOutputs
**
**  \TODO         None
*******************************************************************************/
void ConvCodeDecode::setNumberOutputs(const int n)
{
  m_numberOutputs        = MAX(n, 1);
  m_numberOutputs        = MIN(m_numberOutputs, MAX_CODER_OUTPUTS);
  return;
}

/********************** Public Function **********************************//**
**  \brief              Sets a new number of inputs to the coder
**  \description        Sets a new number of inputs to the coder
**  \param[in]          k          = new number of inputs
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       m_numberInputs
**
**  \TODO         None
*******************************************************************************/
void ConvCodeDecode::setNumberInputs(const int k)
{
  m_numberInputs         = MAX(k, 1);
  m_numberInputs         = MIN(m_numberInputs, MAX_CODER_INPUTS);
  return;
}

/********************** Public Function **********************************//**
**  \brief              Sets a new number of shift register stages for the coder
**  \description        Sets a new number of shift register stages for the coder
**  \param[in]          m          = new number of stages
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       m_numberStages 
**
**  \TODO         None
*******************************************************************************/
void ConvCodeDecode::setNumberStages(const int m)
{
  m_numberStages         = MAX(m, 1);
  m_numberStages         = MIN(m_numberStages, MAX_CODER_STAGES);
  return;
}

/********************** Public Function **********************************//**
**  \brief              Sets a new generator polynomial
**  \description        Sets a new generator polynomial, g
**  \param[in]          gOctal      = the g polynomial in octal format
**  \param[in]          number      = the number of the polynomial [0, MAX_CODER_OUTPUTS-1]
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       m_g 
**
**  \TODO         None
*******************************************************************************/
void ConvCodeDecode::setGPolynomial(const int gOctal,
                                    const int number)
{
  int n                 = MAX(0, number);
  n                     = MIN(n, (MAX_CODER_OUTPUTS-1));
  m_g[n]                = gOctal;
  return;
}


/********************** Public Function **********************************//**
**  \brief              Sets a new puncture polynomial
**  \description        Sets a new puncture polynomial, p
**  \param[in]          pOctal      = the p polynomial in octal format
**  \param[in]          number      = the number of the polynomial [0, MAX_CODER_OUTPUTS-1]
**  \return             None
**   \post         
**  <b>Modified:</b>    m_p
**  <b>Notes:</b>       None 
**
**  \TODO         None
*******************************************************************************/
void ConvCodeDecode::setPPolynomial(const int pOctal,
                                    const int number)
{
  int n                 = MAX(0, number);
  n                     = MIN(n, (MAX_CODER_OUTPUTS-1));
  m_p[n]                = pOctal;
  return;
}


/********************** Public Function **********************************//**
**  \brief              Returns the number of states in the (decoder) trellis
**  \description        Returns the number of states in the (decoder) trellis
**  \param[in]          None
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       m_p 
**
**  \TODO         None
*******************************************************************************/
int ConvCodeDecode::numberStates() const
{
  int   number_states;
  number_states         = 1 << (m_numberInputs * m_numberStages);
  return        number_states;
}

