/***************************************************************************//*!
**  \file         SymbolEstimate.cc
**  \brief        symbol estimator
**  \description  This class takes input I/Q data and outputs symbol estimates
**  \project      BlockDraw
**  \author       Mark Frank
** Revision history:
**  1. 10/17/09  - Started
*/ /*
*******************************************************************************/
// Include Files ///////////////////////////////////////////////////////////////
//
#include "SymbolEstimate.h"

// MACROS ///////////////////////////////////////////////////////////////////
//
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )


/********************** Constructor ****************************************//**
**  \brief        Constructor
**  \description  Class constructor for the Symbol estimator
**  \param[in]    modTYpe          = Type of modulation
**  \return       None
**   \post         
**  <b>Modified:</b>  m_symbolTypes, m_numberSymbols, m_estimatedSymbols, m_iData, m_qData 
**
**  \TODO         None
*******************************************************************************/
SymbolEstimate::SymbolEstimate (eMod_ModType modType)
               :m_symbolTypes(0),
                m_numberSymbols(0),
                m_estimatedSymbols(0),
                m_iData(0),
                m_qData(0)
{
  SetModType(modType);
  return;
}


/********************** Destructor ****************************************//**
**  \brief        Destructor
**  \description  Class destructor for the Symbol estimator
**  \return       None
**   \post         
**  <b>Modified:</b>  None
**
**  \TODO         None
*******************************************************************************/
SymbolEstimate::~SymbolEstimate()
{
  delete [] m_symbolTypes;
  delete [] m_estimatedSymbols;
  delete [] m_iData;
  delete [] m_qData;
  return;
}

/********************** Public Function **********************************//**
**  \brief              Resets the estimates
**  \description        Resets the accumulated statistics for the symbol estimates
**  \return             None
**   \post         
**  <b>Modified:</b>    m_estimatedSymbols
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
void SymbolEstimate::ResetEstimates(){
  m_numberSymbols                       = 0;
  for(int i=0; i<m_numberSymbolTypes; i++)
  {
    m_estimatedSymbols[i].iAverage      = 0.;
    m_estimatedSymbols[i].qAverage      = 0.;
    m_estimatedSymbols[i].numberSymbols = 0;
  }
  return;
}
  

/********************** Public Function **********************************//**
**  \brief              Accumulates statistics for an array of symbols
**  \description        Accumulates statistics for an array of symbols
**  \param[in]          *iData          = I values of input data
**  \param[in]          *qData          = Q values of input data
**  \param[in]          numberSymbols   = # of symbols to estimate
**  \return             None
**   \post         
**  <b>Modified:</b>    m_estimatedSymbols, m_numberSymbols
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
void SymbolEstimate::AccumulateSymbols(                 //!< Accumulate estimate statistics for an array of symbols
                  const double *iData,                  //!< Input I data
                  const double *qData,                  //!< Input Q data
                  const int    numberSymbols)           //!< Size of above arrays
{
  SetMaxNumberSymbols(numberSymbols);
  m_numberSymbols       = 0;
  for(int i=0; i<numberSymbols; i++)
  {
    AccumulateSymbol(iData[i], qData[i]);
  }
}


/********************** Public Function **********************************//**
**  \brief              Accumulates statistics for a single symbol
**  \description        Accumulates statistics for a single symbol
**  \param[in]          iData           = I value of input data
**  \param[in]          qData           = Q value of input data
**  \return             None
**   \post         
**  <b>Modified:</b>    m_estimatedSymbols, m_numberSymbols
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
void SymbolEstimate::AccumulateSymbol(                  //!< Accumulate estimate statistics for a single symbol
                  const double  iData,                  //!< Input I data
                  const double  qData)                  //!< Input Q data
{
  int symbolType                                = GetSymbolType(m_modulationType, iData, qData);
  int symbolNumber                              = MIN(m_maxNumberSymbols, m_numberSymbols);
  m_numberSymbols++;
  m_symbolTypes[symbolNumber]                   = symbolType;
  m_estimatedSymbols[symbolType].iAverage       += iData;
  m_estimatedSymbols[symbolType].qAverage       += qData;
  m_estimatedSymbols[symbolType].numberSymbols++;
}

/********************** Public Function **********************************//**
**  \brief              Estimate symbols
**  \description        Estimates I and Q data arrays for all of the symbols
**  \return             None
**   \post         
**  <b>Modified:</b>    m_estimatedSymbols, m_numberSymbols
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
void SymbolEstimate::EstimateSymbols()          //!< Estimate symbols
{
  delete [] m_iData;
  delete [] m_qData;
  int symbolsToAnalyze  = MIN(m_maxNumberSymbols, m_numberSymbols);
  m_iData               = new double[symbolsToAnalyze];
  m_qData               = new double[symbolsToAnalyze];
  //
  // Estimate the average symbol values:
  //
  for(int i=0; i<m_numberSymbolTypes; i++)
  {
    int numberAverages          = m_estimatedSymbols[i].numberSymbols;
    if(numberAverages > 0)
    {
      m_estimatedSymbols[i].iAverage    /= numberAverages;
      m_estimatedSymbols[i].qAverage    /= numberAverages;
    }
  }
  //
  // Now, fill in the symbols
  //
  for(int i=0; i<symbolsToAnalyze; i++)
  {
    int symbolType              = m_symbolTypes[i];
    m_iData[i]                  = m_estimatedSymbols[symbolType].iAverage;
    m_qData[i]                  = m_estimatedSymbols[symbolType].qAverage;
  }
  return;
}


/********************** Public Function **********************************//**
**  \brief              Set max number of symbols
**  \description        Sets the max number of symbols for m_symbolTypes array
**  \return             None
**   \post         
**  <b>Modified:</b>    m_symbolTypes, m_maxNumberSymbols
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
void SymbolEstimate::SetMaxNumberSymbols(               //!< Sets max number of symbols
                          const int maxNumber)          //!< New max number of symbols
{
  m_maxNumberSymbols    = maxNumber;
  delete [] m_symbolTypes;
  m_symbolTypes         = new int[m_maxNumberSymbols];
  return;
}



/********************** Protected Function **********************************//**
**  \brief              Sets modulation type
**  \description        Sets new modulation type
**  \param[in]          modType           = type of modulation to use in estimator
**  \return             None
**   \post         
**  <b>Modified:</b>    m_modulationType, m_numberSymbolTypes
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
void SymbolEstimate::SetModType(                        //!< Sets new mod type
                          const eMod_ModType modType)   //!< New modulation type
{
  m_modulationType      = modType;
  m_numberSymbolTypes   = 2;
  switch(modType)
  {
    case eMod_BPSK:
    default:
      m_numberSymbolTypes       = 2;
      break;
    case eMod_QPSK:
      m_numberSymbolTypes       = 4;
      break;
  }
  delete [] m_estimatedSymbols;
  m_estimatedSymbols            = new SymbolEstimateStruct[m_numberSymbolTypes];
}
      
  

/********************** Protected Function **********************************//**
**  \brief              Gets symbol type
**  \description        Gets the symbol type from the I and Q data
**  \param[in]          modType         = Modulation type to use in estimate
**  \param[in]          iData           = I value of input data
**  \param[in]          qData           = Q value of input data
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
const int SymbolEstimate::GetSymbolType(                //!< Gets the symbol type/number
                          const eMod_ModType modType,   //!< modulation type
                          double iData,                 //!< I data for symbol type to estimate
                          double qData) const           //!< Q data for symbol type to estimate
{
  int symbolType                = 0;
  switch(modType)
  {
    case eMod_BPSK:
    default:
      symbolType                = GetSymbolTypeForBPSK(iData, qData);
      break;
    case eMod_QPSK:
      symbolType                = GetSymbolTypeForQPSK(iData, qData);
      break;
  }
  return symbolType;
}


/********************** Protected Function **********************************//**
**  \brief              Gets symbol type
**  \description        Gets the symbol type from the I and Q data for BPSK modulation
**  \param[in]          iData           = I value of input data
**  \param[in]          qData           = Q value of input data
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       The symbol type is somewhat arbitrary
**
**  \TODO         None
*******************************************************************************/
const int SymbolEstimate::GetSymbolTypeForBPSK(         //!< Gets the symbol type/number
                          double iData,                 //!< I data for symbol type to estimate
                          double qData) const           //!< Q data for symbol type to estimate
{
  int symbolType                = 0;
  if(iData < 0.)
  {
    symbolType                  = 1;
  }
  return symbolType;
}


/********************** Protected Function **********************************//**
**  \brief              Gets symbol type
**  \description        Gets the symbol type from the I and Q data for QPSK modulation
**  \param[in]          iData           = I value of input data
**  \param[in]          qData           = Q value of input data
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       The symbol type is somewhat arbitrary
**
**  \TODO         None
*******************************************************************************/
const int SymbolEstimate::GetSymbolTypeForQPSK(         //!< Gets the symbol type/number
                          double iData,                 //!< I data for symbol type to estimate
                          double qData) const           //!< Q data for symbol type to estimate
{
  int symbolType                = 0;                    // 1st quadrant
  if(iData < 0.)                                        // 2nd or 3rd
  {
    if(qData > 0.)
    {
      symbolType                = 1;                    // 2nd quadrant
    }
    else
    {
      symbolType                = 2;                    // 3rd quadrant
    }
  }
  else if(qData < 0.)
  {
    symbolType                  = 3;                    // 4th quadrant
  }
  return symbolType;
}



  