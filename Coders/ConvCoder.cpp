/***************************************************************************//*!
**  \file         ConvCoder.cpp
**  \brief        Convolutional coder
**  \description  implements a (n, k, m) convolutional coder.
*                 Where: k = number inputs, n = number outputs, and m = memory size.
**  \project      EE 598 class
**  \author       Mark Frank
**  \Notes        The convolutional coders and decoders are defined by the output path
**                polynomials in octal.  So, for example, for the GSM Sync channel,
**                G0 = 1 + D^3 + D^4, so that the octal rep of the polynomial would
**                be, 23. Similarly for G1 = 1 + D + D^3 + D^4 -> 33 (octal)
********************************************************************************/
#include "ConvCoder.h"
#include <stdio.h>



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
ConvCoder::ConvCoder(int numberOutputs,
                     int numberInputs,
                     int numberStages,
                     int g1,
                     int g2,
                     int g3,
                     int p1,
                     int p2,
                     int p3)
          :ConvCodeDecode(numberOutputs, numberInputs, numberStages, g1, g2, g3, p1, p2, p3)
{
  m_shiftRegister        = NULL;
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
ConvCoder::~ConvCoder()
{
  delete [] m_shiftRegister;
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
void ConvCoder::initCoder()
{
//
// Call super's function:
//
  ConvCodeDecode::initCoder();
//
// Add this:
//
  delete [] m_shiftRegister;
  m_shiftRegister        = new int[m_numberStages];
  for(int i=0; i<m_numberStages; i++)
    m_shiftRegister[i]   = 0;
  return;
}

/********************** Public Function **********************************//**
**  \brief              Sets the shift register to the input state
**  \description        Sets the shift register to the input state
**  \param[in]          state     = new shift register state with bit 0 going to the
**                                  last stage in the shift register.
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       m_shiftPtr, m_shiftRegister
**
**  \TODO         None
*******************************************************************************/
void ConvCoder::setState(int state)
{
  int   i;
  
  m_shiftPtr                     = m_numberStages;
  for(i=0; i<m_numberStages; i++)
  {
    m_shiftPtr--;
    m_shiftRegister[m_shiftPtr]   = state & 1;
    state                       >>= 1;
  }
  return;
}

/********************** Public Function **********************************//**
**  \brief              Returns data out of the coder from an old data input
**  \description        Returns data out of the coder from an old data input
**  \param[in]          None
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
const int *ConvCoder::oldDataOutput() const
{
  return        m_coderOutput;
}


/********************** Public Function **********************************//**
**  \brief              Returns data out of the coder with a new data input
**  \description        Returns data out of the coder with a new data input
**  \param[in]          inputData    = new input data
**  \return             an array of output data, a -1 in the output data
**                      indicates punctured data
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
const int *ConvCoder::newDataOutput(const int inputData)
{
  int           i, puncture_bit;
//
// first, shift the puncture matrix
//
  m_puncturePtr--;
  if(m_puncturePtr < 0)
    m_puncturePtr        = m_punctureSize-1;
//
// Input new data to feed the shift register
//
  unpuncturedOutput(inputData);
//
// Puncture the data by inserting -1s:
//
  puncture_bit  = 1 << m_puncturePtr;                    // Selects puncture bit to check
  for(i=0; i<m_numberOutputs; i++)
  {
    if( (m_pBinary[i] & puncture_bit) == 0)              // Is the data punctured?
      m_coderOutput[i]           = PUNCTURED_OUTPUT;
  }
  return        m_coderOutput;
}

/********************** Public Function **********************************//**
**  \brief              Returns data out of the coder with a new data input
**  \description        Returns data out of the coder with a new data input
**  \param[in]          inputData    = new input data
**  \return             an array of output data, a -1 in the output data
**                      indicates punctured data
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       Similar to the above, but no puncturing
**
**  \TODO         None
*******************************************************************************/
const int *ConvCoder::unpuncturedOutput(const int inputData)
{
  int           i, j, shift_ptr, poly_ptr;
  int           shift_data;
//
// Input new data to feed the shift register
//
  poly_ptr              = 1<<m_numberStages;
  shift_ptr             = m_shiftPtr;
  for(j=0; j<m_numberOutputs; j++)
  {
      if((m_gBinary[j] & poly_ptr) != 0)
        m_coderOutput[j] = inputData;                    // Add the input data to the outputs
      else
        m_coderOutput[j] = 0;
  }
  for(i=0; i<m_numberStages; i++)                        // XOR the outputs of the shift register
  {
      poly_ptr                  >>= 1;
      shift_data                = m_shiftRegister[shift_ptr];
      shift_ptr++;
      if(shift_ptr == m_numberStages)
        shift_ptr               = 0;
      for(j=0; j<m_numberOutputs; j++)
      {
        if((m_gBinary[j] & poly_ptr) != 0)
          m_coderOutput[j]       ^= shift_data;
      }
  }
//
// Shift the new input data into the shift register:
//
  m_shiftPtr--;
  if(m_shiftPtr < 0)
    m_shiftPtr                   = m_numberStages-1;
  m_shiftRegister[m_shiftPtr]     = inputData;
    
  return        m_coderOutput;
}


/********************** Public Function **********************************//**
**  \brief              Returns the current state of the shift register
**  \description        Returns the current state of the shift register
**  \param[in]          None
**  \return             current shift register state with bit 0 equal to the
**                      last stage in the shift register
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
const int ConvCoder::state() const
{
  int   i, shift_ptr, power_of_2;
  int   state = 0;
  
  shift_ptr     = m_shiftPtr;
  power_of_2    = 1 << (m_numberStages-1);
  for(i=0; i<m_numberStages; i++)
  {
    if(m_shiftRegister[shift_ptr] > 0)
      state     |= power_of_2;
    power_of_2  >>= 1;
    shift_ptr++;
    if(shift_ptr == m_numberStages)
      shift_ptr = 0;
  }
  return state;
}

