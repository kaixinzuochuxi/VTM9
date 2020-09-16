/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.
 *
 * Copyright (c) 2010-2020, ITU/ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/** \file     RateCtrl.h
    \brief    Rate control manager class
*/

#ifndef __ENCRATECTRL__
#define __ENCRATECTRL__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "../CommonLib/CommonDef.h"

#include <vector>
#include <algorithm>

using namespace std;

//! \ingroup EncoderLib
//! \{

#include "../EncoderLib/EncCfg.h"
#include <list>

const int g_RCInvalidQPValue = -999;
const int g_RCSmoothWindowSize = 40;
const int g_RCMaxPicListSize = 32;
const double g_RCWeightPicTargetBitInGOP    = 0.9;
const double g_RCWeightPicRargetBitInBuffer = 1.0 - g_RCWeightPicTargetBitInGOP;
const int g_RCIterationNum = 20;
const double g_RCWeightHistoryLambda = 0.5;
const double g_RCWeightCurrentLambda = 1.0 - g_RCWeightHistoryLambda;
const int g_RCLCUSmoothWindowSize = 4;
const double g_RCAlphaMinValue = 0.05;
const double g_RCAlphaMaxValue = 500.0;
const double g_RCBetaMinValue  = -3.0;
const double g_RCBetaMaxValue  = -0.1;

#define ALPHA     6.7542;
#define BETA1     1.2517
#define BETA2     1.7860

struct TRCLCU
{
  int m_actualBits;
  int m_QP;     // QP of skip mode is set to g_RCInvalidQPValue
  int m_targetBits;
  double m_lambda;
  double m_bitWeight;
  int m_numberOfPixel;
  double m_costIntra;
  int m_targetBitsLeft;
  double m_actualSSE;
  double m_actualMSE;
};

struct TRCParameter
{
  double m_alpha;
  double m_beta;
  int    m_validPix;
  double m_skipRatio;
};

class EncRCSeq
{
public:
  EncRCSeq();
  ~EncRCSeq();

public:
  void create( int totalFrames, int targetBitrate, int frameRate, int GOPSize, int picWidth, int picHeight, int LCUWidth, int LCUHeight, int numberOfLevel, bool useLCUSeparateModel, int adaptiveBit );
  void destroy();
  void initBitsRatio( int bitsRatio[] );
  void initGOPID2Level( int GOPID2Level[] );
  void initPicPara( TRCParameter* picPara  = NULL );    // NULL to initial with default value
  void initLCUPara( TRCParameter** LCUPara = NULL );    // NULL to initial with default value
  void updateAfterPic ( int bits );
  void setAllBitRatio( double basicLambda, double* equaCoeffA, double* equaCoeffB );

public:
  int  getTotalFrames()                 { return m_totalFrames; }
  int  getTargetRate()                  { return m_targetRate; }
  int  getFrameRate()                   { return m_frameRate; }
  int  getGOPSize()                     { return m_GOPSize; }
  int  getPicWidth()                    { return m_picWidth; }
  int  getPicHeight()                   { return m_picHeight; }
  int  getLCUWidth()                    { return m_LCUWidth; }
  int  getLCUHeight()                   { return m_LCUHeight; }
  int  getNumberOfLevel()               { return m_numberOfLevel; }
  int  getAverageBits()                 { return m_averageBits; }
  int  getLeftAverageBits()             { CHECK(!( m_framesLeft > 0 ), "No frames left"); return (int)(m_bitsLeft / m_framesLeft); }
  bool getUseLCUSeparateModel()         { return m_useLCUSeparateModel; }

  int  getNumPixel()                    { return m_numberOfPixel; }
  int64_t  getTargetBits()                { return m_targetBits; }
  int  getNumberOfLCU()                 { return m_numberOfLCU; }
  int* getBitRatio()                    { return m_bitsRatio; }
  int  getBitRatio( int idx )           { CHECK(!( idx<m_GOPSize), "Idx exceeds GOP size"); return m_bitsRatio[idx]; }
  int* getGOPID2Level()                 { return m_GOPID2Level; }
  int  getGOPID2Level( int ID )         { CHECK(!( ID < m_GOPSize ), "Idx exceeds GOP size"); return m_GOPID2Level[ID]; }
  TRCParameter*  getPicPara()                                   { return m_picPara; }
  TRCParameter   getPicPara( int level )                        { CHECK(!( level < m_numberOfLevel ), "Level too big"); return m_picPara[level]; }
  void           setPicPara( int level, TRCParameter para )     { CHECK(!( level < m_numberOfLevel ), "Level too big"); m_picPara[level] = para; }
  TRCParameter** getLCUPara()                                   { return m_LCUPara; }
  TRCParameter*  getLCUPara( int level )                        { CHECK(!( level < m_numberOfLevel ), "Level too big"); return m_LCUPara[level]; }
  TRCParameter   getLCUPara( int level, int LCUIdx )            { CHECK(!( LCUIdx  < m_numberOfLCU ), "LCU id exceeds number of LCU"); return getLCUPara(level)[LCUIdx]; }
  void           setLCUPara( int level, int LCUIdx, TRCParameter para ) { CHECK(!( level < m_numberOfLevel ), "Level too big"); CHECK(!( LCUIdx  < m_numberOfLCU ), "LCU id exceeds number of LCU"); m_LCUPara[level][LCUIdx] = para; }

  int  getFramesLeft()                  { return m_framesLeft; }
  int64_t  getBitsLeft()                  { return m_bitsLeft; }

  double getSeqBpp()                    { return m_seqTargetBpp; }
  double getAlphaUpdate()               { return m_alphaUpdate; }
  double getBetaUpdate()                { return m_betaUpdate; }

  int    getAdaptiveBits()              { return m_adaptiveBit;  }
  double getLastLambda()                { return m_lastLambda;   }
  void   setLastLambda( double lamdba ) { m_lastLambda = lamdba; }
  void setBitDepth(int bitDepth) { m_bitDepth = bitDepth; }
  int getbitDepth() { return m_bitDepth; }

private:
  int m_totalFrames;
  int m_targetRate;
  int m_frameRate;
  int m_GOPSize;
  int m_picWidth;
  int m_picHeight;
  int m_LCUWidth;
  int m_LCUHeight;
  int m_numberOfLevel;
  int m_averageBits;

  int m_numberOfPixel;
  int64_t m_targetBits;
  int m_numberOfLCU;
  int* m_bitsRatio;
  int* m_GOPID2Level;
  TRCParameter*  m_picPara;
  TRCParameter** m_LCUPara;

  int m_framesLeft;
  int64_t m_bitsLeft;
  double m_seqTargetBpp;
  double m_alphaUpdate;
  double m_betaUpdate;
  bool m_useLCUSeparateModel;

  int m_adaptiveBit;
  double m_lastLambda;
  int m_bitDepth;
};

class EncRCGOP
{
public:
  EncRCGOP();
  ~EncRCGOP();

public:
  void create( EncRCSeq* encRCSeq, int numPic );
  void destroy();
  void updateAfterPicture( int bitsCost );

private:
  int  xEstGOPTargetBits( EncRCSeq* encRCSeq, int GOPSize );
  void   xCalEquaCoeff( EncRCSeq* encRCSeq, double* lambdaRatio, double* equaCoeffA, double* equaCoeffB, int GOPSize );
  double xSolveEqua(EncRCSeq* encRCSeq, double targetBpp, double* equaCoeffA, double* equaCoeffB, int GOPSize);

public:
  EncRCSeq* getEncRCSeq()        { return m_encRCSeq; }
  int  getNumPic()                { return m_numPic;}
  int  getTargetBits()            { return m_targetBits; }
  int  getPicLeft()               { return m_picLeft; }
  int  getBitsLeft()              { return m_bitsLeft; }
  int  getTargetBitInGOP( int i ) { return m_picTargetBitInGOP[i]; }
  double getMinEstLambda()        { return m_minEstLambda; }
  double getMaxEstLambda()        { return m_maxEstLambda; }

private:
  EncRCSeq* m_encRCSeq;
  int* m_picTargetBitInGOP;
  int m_numPic;
  int m_targetBits;
  int m_picLeft;
  int m_bitsLeft;
  double m_minEstLambda;
  double m_maxEstLambda;
};

class EncRCPic
{
public:
  EncRCPic();
  ~EncRCPic();

public:
  void create( EncRCSeq* encRCSeq, EncRCGOP* encRCGOP, int frameLevel, list<EncRCPic*>& listPreviousPictures );
  void destroy();

  int    estimatePicQP    ( double lambda, list<EncRCPic*>& listPreviousPictures );
  int    getRefineBitsForIntra(int orgBits);
  double calculateLambdaIntra(double alpha, double beta, double MADPerPixel, double bitsPerPixel);
  double estimatePicLambda( list<EncRCPic*>& listPreviousPictures, bool isIRAP);

  void   updateAlphaBetaIntra(double *alpha, double *beta);

  double getLCUTargetBpp(bool isIRAP);
  double getLCUEstLambdaAndQP(double bpp, int clipPicQP, int *estQP);
  double getLCUEstLambda( double bpp );
  int    getLCUEstQP( double lambda, int clipPicQP );
  void updateAfterCTU(int LCUIdx, int bits, int QP, double lambda, double skipRatio, bool updateLCUParameter = true);
  void updateAfterPicture( int actualHeaderBits, int actualTotalBits, double averageQP, double averageLambda, bool isIRAP);

  double clipRcAlpha(const int bitdepth, const double alpha);
  double clipRcBeta(const double beta);

  void addToPictureLsit( list<EncRCPic*>& listPreviousPictures );
  double calAverageQP();
  double calAverageLambda();

private:
  int xEstPicTargetBits( EncRCSeq* encRCSeq, EncRCGOP* encRCGOP );
  int xEstPicHeaderBits( list<EncRCPic*>& listPreviousPictures, int frameLevel );
#if V0078_ADAPTIVE_LOWER_BOUND
  int xEstPicLowerBound( EncRCSeq* encRCSeq, EncRCGOP* encRCGOP );
#endif

public:
  EncRCSeq*      getRCSequence()                         { return m_encRCSeq; }
  EncRCGOP*      getRCGOP()                              { return m_encRCGOP; }

  int  getFrameLevel()                                    { return m_frameLevel; }
  int  getNumberOfPixel()                                 { return m_numberOfPixel; }
  int  getNumberOfLCU()                                   { return m_numberOfLCU; }
  int  getTargetBits()                                    { return m_targetBits; }
  int  getEstHeaderBits()                                 { return m_estHeaderBits; }
  int  getLCULeft()                                       { return m_LCULeft; }
  int  getBitsLeft()                                      { return m_bitsLeft; }
  int  getPixelsLeft()                                    { return m_pixelsLeft; }
  int  getBitsCoded()                                     { return m_targetBits - m_estHeaderBits - m_bitsLeft; }
  int  getLCUCoded()                                      { return m_numberOfLCU - m_LCULeft; }
#if V0078_ADAPTIVE_LOWER_BOUND
  int  getLowerBound()                                    { return m_lowerBound; }
#endif
  TRCLCU* getLCU()                                        { return m_LCUs; }
  TRCLCU& getLCU( int LCUIdx )                            { return m_LCUs[LCUIdx]; }
  int  getPicActualHeaderBits()                           { return m_picActualHeaderBits; }
#if U0132_TARGET_BITS_SATURATION
  void setBitLeft(int bits)                               { m_bitsLeft = bits; }
#endif
  void setTargetBits( int bits )                          { m_targetBits = bits; m_bitsLeft = bits;}
  void setTotalIntraCost(double cost)                     { m_totalCostIntra = cost; }
  void getLCUInitTargetBits();

  int  getPicActualBits()                                 { return m_picActualBits; }
  int  getPicActualQP()                                   { return m_picQP; }
  double getPicActualLambda()                             { return m_picLambda; }
  int  getPicEstQP()                                      { return m_estPicQP; }
  void setPicEstQP( int QP )                              { m_estPicQP = QP; }
  double getPicEstLambda()                                { return m_estPicLambda; }
  void setPicEstLambda( double lambda )                   { m_picLambda = lambda; }
  double getPicMSE()                                      { return m_picMSE; }
  void  setPicMSE(double avgMSE)                           { m_picMSE = avgMSE; }

private:
  EncRCSeq* m_encRCSeq;
  EncRCGOP* m_encRCGOP;

  int m_frameLevel;
  int m_numberOfPixel;
  int m_numberOfLCU;
  int m_targetBits;
  int m_estHeaderBits;
  int m_estPicQP;
#if V0078_ADAPTIVE_LOWER_BOUND
  int m_lowerBound;
#endif
  double m_estPicLambda;

  int m_LCULeft;
  int m_bitsLeft;
  int m_pixelsLeft;

  TRCLCU* m_LCUs;
  int m_picActualHeaderBits;    // only SH and potential APS
  double m_totalCostIntra;
  double m_remainingCostIntra;
  int m_picActualBits;          // the whole picture, including header
  int m_picQP;                  // in integer form
  double m_picLambda;
  double m_picMSE;
  int m_validPixelsInPic;
};

class RateCtrl
{
public:
  RateCtrl();
  ~RateCtrl();

public:
  void init(int totalFrames, int targetBitrate, int frameRate, int GOPSize, int picWidth, int picHeight, int LCUWidth, int LCUHeight, int bitDepth, int keepHierBits, bool useLCUSeparateModel, GOPEntry GOPList[MAX_GOP]);
  void destroy();
  void initRCPic( int frameLevel );
  void initRCGOP( int numberOfPictures );
  void destroyRCGOP();

public:
  void       setRCQP ( int QP ) { m_RCQP = QP;   }
  int        getRCQP () const   { return m_RCQP; }
  EncRCSeq* getRCSeq()          { CHECK( m_encRCSeq == NULL, "Object does not exist" ); return m_encRCSeq; }
  EncRCGOP* getRCGOP()          { CHECK( m_encRCGOP == NULL, "Object does not exist" ); return m_encRCGOP; }
  EncRCPic* getRCPic()          { CHECK( m_encRCPic == NULL, "Object does not exist" ); return m_encRCPic; }
  list<EncRCPic*>& getPicList() { return m_listRCPictures; }
#if U0132_TARGET_BITS_SATURATION
  bool       getCpbSaturationEnabled()  { return m_CpbSaturationEnabled;  }
  uint32_t       getCpbState()              { return m_cpbState;       }
  uint32_t       getCpbSize()               { return m_cpbSize;        }
  uint32_t       getBufferingRate()         { return m_bufferingRate;  }
  int        updateCpbState(int actualBits);
  void       initHrdParam(const GeneralHrdParams* generalHrd, const OlsHrdParams* olsHrd, int iFrameRate, double fInitialCpbFullness);
#endif

private:
  EncRCSeq* m_encRCSeq;
  EncRCGOP* m_encRCGOP;
  EncRCPic* m_encRCPic;
  list<EncRCPic*> m_listRCPictures;
  int        m_RCQP;
#if U0132_TARGET_BITS_SATURATION
  bool       m_CpbSaturationEnabled;    // Enable target bits saturation to avoid CPB overflow and underflow
  int        m_cpbState;                // CPB State
  uint32_t       m_cpbSize;                 // CPB size
  uint32_t       m_bufferingRate;           // Buffering rate
#endif
};



//#if pre_ana
//
//
//class pre_analysis
//{
//public:
//  int m_size;
//  vector<PelBuf*> pre_ana_buf;
//  int curidx;
//  int CTUsize = 64;
//  int framew;
//  int frameh;
//
//  pre_analysis(int iframew, int iframeh) {
//    m_size = 0; curidx = 0; framew = iframew; frameh = iframeh;
//  }
//  //~pre_analysis();
//
//  void createbuf(int w,int h)
//  {
//    for (int i = 0; i < m_size; i++)
//    {
//      //pre_ana_buf.push_back( new PelBuf);
//      Pel *tbuf = new Pel [w*h];
//      
//      //PelBuf temp((Pel*)xMalloc(Pel, w*h), w, h);
//      PelBuf temp(tbuf, w, h);
//      pre_ana_buf.push_back(&temp);
//      
//    }
//  }
//
//  void addbuf(int n)
//  {
//    if (n > 0)
//    {
//      for (int i = 0; i < n; i++)
//      {
//        pre_ana_buf.push_back(new PelBuf);
//      }
//    }
//  }
//
//  //void update(int n)
//  //{
//  //  for (int i = 0; i < n; i++)
//  //  {
//  //    pre_ana_buf.push_back(pre_ana_buf.front()); pre_ana_buf.pop_front();
//  //  }
//  //}
//
//  void destroy()
//  {
//    while (pre_ana_buf.size()!=0)
//    {
//      delete pre_ana_buf.back()->buf;
//      pre_ana_buf.back()->buf = NULL;
//      delete pre_ana_buf.back();
//      pre_ana_buf.back() = NULL;
//      pre_ana_buf.pop_back();
//    }
//    
//  }
//
//  int xCalcHADs8x8_ISlice(const Pel *piOrg, const int iStrideOrg)
//  {
//    int k, i, j, jj;
//    int diff[64], m1[8][8], m2[8][8], m3[8][8], iSumHad = 0;
//
//    for (k = 0; k < 64; k += 8)
//    {
//      diff[k + 0] = piOrg[0];
//      diff[k + 1] = piOrg[1];
//      diff[k + 2] = piOrg[2];
//      diff[k + 3] = piOrg[3];
//      diff[k + 4] = piOrg[4];
//      diff[k + 5] = piOrg[5];
//      diff[k + 6] = piOrg[6];
//      diff[k + 7] = piOrg[7];
//
//      piOrg += iStrideOrg;
//    }
//
//    //horizontal
//    for (j = 0; j < 8; j++)
//    {
//      jj = j << 3;
//      m2[j][0] = diff[jj] + diff[jj + 4];
//      m2[j][1] = diff[jj + 1] + diff[jj + 5];
//      m2[j][2] = diff[jj + 2] + diff[jj + 6];
//      m2[j][3] = diff[jj + 3] + diff[jj + 7];
//      m2[j][4] = diff[jj] - diff[jj + 4];
//      m2[j][5] = diff[jj + 1] - diff[jj + 5];
//      m2[j][6] = diff[jj + 2] - diff[jj + 6];
//      m2[j][7] = diff[jj + 3] - diff[jj + 7];
//
//      m1[j][0] = m2[j][0] + m2[j][2];
//      m1[j][1] = m2[j][1] + m2[j][3];
//      m1[j][2] = m2[j][0] - m2[j][2];
//      m1[j][3] = m2[j][1] - m2[j][3];
//      m1[j][4] = m2[j][4] + m2[j][6];
//      m1[j][5] = m2[j][5] + m2[j][7];
//      m1[j][6] = m2[j][4] - m2[j][6];
//      m1[j][7] = m2[j][5] - m2[j][7];
//
//      m2[j][0] = m1[j][0] + m1[j][1];
//      m2[j][1] = m1[j][0] - m1[j][1];
//      m2[j][2] = m1[j][2] + m1[j][3];
//      m2[j][3] = m1[j][2] - m1[j][3];
//      m2[j][4] = m1[j][4] + m1[j][5];
//      m2[j][5] = m1[j][4] - m1[j][5];
//      m2[j][6] = m1[j][6] + m1[j][7];
//      m2[j][7] = m1[j][6] - m1[j][7];
//    }
//
//    //vertical
//    for (i = 0; i < 8; i++)
//    {
//      m3[0][i] = m2[0][i] + m2[4][i];
//      m3[1][i] = m2[1][i] + m2[5][i];
//      m3[2][i] = m2[2][i] + m2[6][i];
//      m3[3][i] = m2[3][i] + m2[7][i];
//      m3[4][i] = m2[0][i] - m2[4][i];
//      m3[5][i] = m2[1][i] - m2[5][i];
//      m3[6][i] = m2[2][i] - m2[6][i];
//      m3[7][i] = m2[3][i] - m2[7][i];
//
//      m1[0][i] = m3[0][i] + m3[2][i];
//      m1[1][i] = m3[1][i] + m3[3][i];
//      m1[2][i] = m3[0][i] - m3[2][i];
//      m1[3][i] = m3[1][i] - m3[3][i];
//      m1[4][i] = m3[4][i] + m3[6][i];
//      m1[5][i] = m3[5][i] + m3[7][i];
//      m1[6][i] = m3[4][i] - m3[6][i];
//      m1[7][i] = m3[5][i] - m3[7][i];
//
//      m2[0][i] = m1[0][i] + m1[1][i];
//      m2[1][i] = m1[0][i] - m1[1][i];
//      m2[2][i] = m1[2][i] + m1[3][i];
//      m2[3][i] = m1[2][i] - m1[3][i];
//      m2[4][i] = m1[4][i] + m1[5][i];
//      m2[5][i] = m1[4][i] - m1[5][i];
//      m2[6][i] = m1[6][i] + m1[7][i];
//      m2[7][i] = m1[6][i] - m1[7][i];
//    }
//
//    for (i = 0; i < 8; i++)
//    {
//      for (j = 0; j < 8; j++)
//      {
//        iSumHad += abs(m2[i][j]);
//      }
//    }
//    iSumHad -= abs(m2[0][0]);
//    iSumHad = (iSumHad + 2) >> 2;
//    return(iSumHad);
//  }
//
//  int CalIntraSATD(int idx)
//  {
//    /*for (int y = 0; y < frameh; y += CTUsize)
//    {
//      for (int x=0;x<framew;x+=CTUsize)
//      {
//        Area CurCU = Area(x, y, min(CTUsize, framew - x), min(CTUsize, frameh - y));
//
//      }
//        
//    }
//  */
//
//    /////default 
//    /*int  xBl, yBl;
//    const int iBlkSize = CTUsize;
//    const Pel* pOrgInit = pre_ana_buf[idx]->buf;
//    int  iStrideOrig = pre_ana_buf[idx]->stride;
//
//    int iSumHad = 0;
//    printf("\n");
//
//    for (yBl = 0; (yBl + iBlkSize) <= pre_ana_buf[idx]->height; yBl += iBlkSize)
//    {
//      for (xBl = 0; (xBl + iBlkSize) <= pre_ana_buf[idx]->width; xBl += iBlkSize)
//      {
//        const Pel* pOrg = pOrgInit + iStrideOrig * yBl + xBl;
//        printf("%d\t", xCalcHADs8x8_ISlice(pOrg, iStrideOrig));
//        iSumHad += xCalcHADs8x8_ISlice(pOrg, iStrideOrig);
//        
//      }
//    }
//    return(iSumHad);*/
//  
//    Pel* tbuf = (Pel*) xMalloc(Pel, CTUsize*CTUsize);
//    PelBuf tmp(tbuf, CTUsize, CTUsize);
//    for (int y = 0; y < frameh; y += CTUsize)
//    {
//      for (int x = 0; x < framew; x += CTUsize)
//      {
//        Area CurCU = Area(x, y, min(CTUsize, framew - x), min(CTUsize, frameh - y));
//        tmp.fill(pre_ana_buf[idx]->subBuf(x, y, CTUsize, CTUsize).computeAvg());
//        auto dist = m_pcRdCost->getDistPart(org, reco, cs.sps->getBitDepth(toChannelType(compID)), compID, DF_SSE);
//
//      }
//    }
//  }
//};
//
//
//
//#endif

#endif


