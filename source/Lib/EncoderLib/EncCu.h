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

/** \file     EncCu.h
    \brief    Coding Unit (CU) encoder class (header)
*/

#ifndef __ENCCU__
#define __ENCCU__

// Include files
#include "CommonLib/CommonDef.h"
#include "CommonLib/IntraPrediction.h"
#include "CommonLib/InterPrediction.h"
#include "CommonLib/TrQuant.h"
#include "CommonLib/Unit.h"
#include "CommonLib/UnitPartitioner.h"
#include "CommonLib/IbcHashMap.h"
#include "CommonLib/LoopFilter.h"

#include "DecoderLib/DecCu.h"

#include "CABACWriter.h"
#include "IntraSearch.h"
#include "InterSearch.h"
#include "RateCtrl.h"
#include "EncModeCtrl.h"

#if yang2019content
//#include "CommonLib/TrQuant.cpp"
#include "CommonLib/TrQuant_EMT.h"
#include <numeric>
#endif
#if UsePipe
#if iswindows
#include "CommonLib/TypeDef.h"
#endif
#endif
//! \ingroup EncoderLib
//! \{

class EncLib;
class HLSWriter;
class EncSlice;

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// CU encoder class
struct GeoMergeCombo
{
  int splitDir;
  int mergeIdx0;
  int mergeIdx1;
  double cost;
  GeoMergeCombo() : splitDir(), mergeIdx0(-1), mergeIdx1(-1), cost(0.0) {};
  GeoMergeCombo(int _splitDir, int _mergeIdx0, int _mergeIdx1, double _cost) : splitDir(_splitDir), mergeIdx0(_mergeIdx0), mergeIdx1(_mergeIdx1), cost(_cost) {};
};
struct GeoMotionInfo
{
  uint8_t   m_candIdx0;
  uint8_t   m_candIdx1;

  GeoMotionInfo(uint8_t candIdx0, uint8_t candIdx1) : m_candIdx0(candIdx0), m_candIdx1(candIdx1) { }
  GeoMotionInfo() { m_candIdx0 = m_candIdx1 = 0; }
};
struct SmallerThanComboCost
{
  inline bool operator() (const GeoMergeCombo& first, const GeoMergeCombo& second)
  {
      return (first.cost < second.cost);
  }
};
class GeoComboCostList
{
public:
  GeoComboCostList() {};
  ~GeoComboCostList() {};
  std::vector<GeoMergeCombo> list;
  void sortByCost() { std::sort(list.begin(), list.end(), SmallerThanComboCost()); };
};
struct SingleGeoMergeEntry
{
  int mergeIdx;
  double cost;
  SingleGeoMergeEntry() : mergeIdx(0), cost(MAX_DOUBLE) {};
  SingleGeoMergeEntry(int _mergeIdx, double _cost) : mergeIdx(_mergeIdx), cost(_cost) {};
};
class FastGeoCostList
{
public:
  FastGeoCostList() { numGeoTemplatesInitialized = 0; };
  ~FastGeoCostList()
  {
    for (int partIdx = 0; partIdx < 2; partIdx++)
    {
      for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; splitDir++)
      {
        delete[] singleDistList[partIdx][splitDir];
      }
      delete[] singleDistList[partIdx];
      singleDistList[partIdx] = nullptr;
    }
  };
  SingleGeoMergeEntry** singleDistList[2];
  void init(int numTemplates, int maxNumGeoCand)
  {
    if (numGeoTemplatesInitialized == 0 || numGeoTemplatesInitialized < numTemplates)
    {
      for (int partIdx = 0; partIdx < 2; partIdx++)
      {
        singleDistList[partIdx] = new SingleGeoMergeEntry*[numTemplates];
        for (int splitDir = 0; splitDir < numTemplates; splitDir++)
        {
          singleDistList[partIdx][splitDir] = new SingleGeoMergeEntry[maxNumGeoCand];
        }
      }
      numGeoTemplatesInitialized = numTemplates;
    }
  }
  void insert(int geoIdx, int partIdx, int mergeIdx, double cost)
  {
    assert(geoIdx < numGeoTemplatesInitialized);
    singleDistList[partIdx][geoIdx][mergeIdx] = SingleGeoMergeEntry(mergeIdx, cost);
  }
  int numGeoTemplatesInitialized;
};

class EncCu
  : DecCu
{
private:
  bool m_bestModeUpdated;
  struct CtxPair
  {
    Ctx start;
    Ctx best;
  };

  std::vector<CtxPair>  m_CtxBuffer;
  CtxPair*              m_CurrCtx;
  CtxCache*             m_CtxCache;

#if ENABLE_SPLIT_PARALLELISM
  int                   m_dataId;
#endif

  //  Data : encoder control
  int                   m_cuChromaQpOffsetIdxPlus1; // if 0, then cu_chroma_qp_offset_flag will be 0, otherwise cu_chroma_qp_offset_flag will be 1.

  XUCache               m_unitCache;

  CodingStructure    ***m_pTempCS;
  CodingStructure    ***m_pBestCS;
  CodingStructure    ***m_pTempCS2;
  CodingStructure    ***m_pBestCS2;
  //  Access channel
  EncCfg*               m_pcEncCfg;
  IntraSearch*          m_pcIntraSearch;
  InterSearch*          m_pcInterSearch;
  TrQuant*              m_pcTrQuant;
  RdCost*               m_pcRdCost;
  EncSlice*             m_pcSliceEncoder;
  LoopFilter*           m_pcLoopFilter;

  CABACWriter*          m_CABACEstimator;
  RateCtrl*             m_pcRateCtrl;
  IbcHashMap            m_ibcHashMap;
  EncModeCtrl          *m_modeCtrl;

  PelStorage            m_acMergeBuffer[MMVD_MRG_MAX_RD_BUF_NUM];
  PelStorage            m_acRealMergeBuffer[MRG_MAX_NUM_CANDS];
  PelStorage            m_acMergeTmpBuffer[MRG_MAX_NUM_CANDS];
  PelStorage            m_acGeoWeightedBuffer[GEO_MAX_TRY_WEIGHTED_SAD]; // to store weighted prediction pixles
  FastGeoCostList       m_GeoCostList;
  double                m_AFFBestSATDCost;
  double                m_mergeBestSATDCost;
  MotionInfo            m_SubPuMiBuf      [( MAX_CU_SIZE * MAX_CU_SIZE ) >> ( MIN_CU_LOG2 << 1 )];

  int                   m_ctuIbcSearchRangeX;
  int                   m_ctuIbcSearchRangeY;
#if ENABLE_SPLIT_PARALLELISM
  EncLib*               m_pcEncLib;
#endif
  int                   m_bestBcwIdx[2];
  double                m_bestBcwCost[2];
  GeoMotionInfo         m_GeoModeTest[GEO_MAX_NUM_CANDS];
#if SHARP_LUMA_DELTA_QP || ENABLE_QPA_SUB_CTU
  void    updateLambda      ( Slice* slice, const int dQP,
 #if WCG_EXT && ER_CHROMA_QP_WCG_PPS
                              const bool useWCGChromaControl,
 #endif
                              const bool updateRdCostLambda );
#endif
  double                m_sbtCostSave[2];
public:
  /// copy parameters from encoder class
  void  init                ( EncLib* pcEncLib, const SPS& sps PARL_PARAM( const int jId = 0 ) );

  void setDecCuReshaperInEncCU(EncReshape* pcReshape, ChromaFormat chromaFormatIDC) { initDecCuReshaper((Reshape*) pcReshape, chromaFormatIDC); }
  /// create internal buffers
  void  create              ( EncCfg* encCfg );

  /// destroy internal buffers
  void  destroy             ();

  /// CTU analysis function
  void  compressCtu         ( CodingStructure& cs, const UnitArea& area, const unsigned ctuRsAddr, const int prevQP[], const int currQP[] );
  /// CTU encoding function
  int   updateCtuDataISlice ( const CPelBuf buf );

  EncModeCtrl* getModeCtrl  () { return m_modeCtrl; }


  void   setMergeBestSATDCost(double cost) { m_mergeBestSATDCost = cost; }
  double getMergeBestSATDCost()            { return m_mergeBestSATDCost; }
  void   setAFFBestSATDCost(double cost)   { m_AFFBestSATDCost = cost; }
  double getAFFBestSATDCost()              { return m_AFFBestSATDCost; }
  IbcHashMap& getIbcHashMap()              { return m_ibcHashMap;        }
  EncCfg*     getEncCfg()            const { return m_pcEncCfg;          }

  EncCu();
  ~EncCu();

protected:

  void xCalDebCost            ( CodingStructure &cs, Partitioner &partitioner, bool calDist = false );
  Distortion getDistortionDb  ( CodingStructure &cs, CPelBuf org, CPelBuf reco, ComponentID compID, const CompArea& compArea, bool afterDb );

  void xCompressCU            ( CodingStructure*& tempCS, CodingStructure*& bestCS, Partitioner& pm, double maxCostAllowed = MAX_DOUBLE );
#if ENABLE_SPLIT_PARALLELISM
  void xCompressCUParallel    ( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &pm );
  void copyState              ( EncCu* other, Partitioner& pm, const UnitArea& currArea, const bool isDist );
#endif

  bool
    xCheckBestMode         ( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &pm, const EncTestMode& encTestmode );

  void xCheckModeSplit        ( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &pm, const EncTestMode& encTestMode, const ModeType modeTypeParent, bool &skipInterPass );

  bool xCheckRDCostIntra(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &pm, const EncTestMode& encTestMode, bool adaptiveColorTrans);

  void xCheckDQP              ( CodingStructure& cs, Partitioner& partitioner, bool bKeepCtx = false);
  void xCheckChromaQPOffset   ( CodingStructure& cs, Partitioner& partitioner);
  void xFillPCMBuffer         ( CodingUnit &cu);

  void xCheckRDCostHashInter  ( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &pm, const EncTestMode& encTestMode );
  void xCheckRDCostAffineMerge2Nx2N
                              ( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode );
  void xCheckRDCostInter      ( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &pm, const EncTestMode& encTestMode );
  bool xCheckRDCostInterIMV(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &pm, const EncTestMode& encTestMode, double &bestIntPelCost);
  void xEncodeDontSplit       ( CodingStructure &cs, Partitioner &partitioner);

  void xCheckRDCostMerge2Nx2N ( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &pm, const EncTestMode& encTestMode );

  void xCheckRDCostMergeGeo2Nx2N(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &pm, const EncTestMode& encTestMode);

  void xEncodeInterResidual(   CodingStructure *&tempCS
                             , CodingStructure *&bestCS
                             , Partitioner &partitioner
                             , const EncTestMode& encTestMode
                             , int residualPass       = 0
                             , bool* bestHasNonResi   = NULL
                             , double* equBcwCost     = NULL
                           );
#if REUSE_CU_RESULTS
  void xReuseCachedResult     ( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &Partitioner );
#endif
  bool xIsBcwSkip(const CodingUnit& cu)
  {
    if (cu.slice->getSliceType() != B_SLICE)
    {
      return true;
    }
    return((m_pcEncCfg->getBaseQP() > 32) && ((cu.slice->getTLayer() >= 4)
       || ((cu.refIdxBi[0] >= 0 && cu.refIdxBi[1] >= 0)
       && (abs(cu.slice->getPOC() - cu.slice->getRefPOC(REF_PIC_LIST_0, cu.refIdxBi[0])) == 1
       ||  abs(cu.slice->getPOC() - cu.slice->getRefPOC(REF_PIC_LIST_1, cu.refIdxBi[1])) == 1))));
  }
  void xCheckRDCostIBCMode    ( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &pm, const EncTestMode& encTestMode );
  void xCheckRDCostIBCModeMerge2Nx2N( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode );

  void xCheckPLT              ( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode );

#if pre_ana
  public:
    RdCost*     getRDcost() { return      m_pcRdCost; }
#endif

};

//! \}



#if pre_ana
const int g_GOPSizeRA = 16;
const int g_GOPSizeLD = 8;
const int g_presizeLD = 80;
int getRPLIdxLDB(int poc);
int getRPLIdxRA(int poc);


class pre_analysis
{
public:
  
  int cfgctusize=0;
  int CTUsize = 0;
  int framew;
  int frameh;
  RdCost *m_pcRdCost;
  int TotalCTUNum;
  int hieStruct; // 0:ai,1:LD,2:RA
  DFunc costfun = DF_HAD;

  int m_size;
  int curidx;
  int swEndIdx;
  vector<PelBuf*> pre_ana_buf;
  vector<vector<vector<uint64_t>>> CUSATD;
  vector<uint64_t> FrameSATD;
  int IBCRef[4][3][2];
  vector<int> encodingorder;

#if yang2019content
  vector<int> scenechange;
  vector<uint64_t> fd;
  vector<uint64_t> pastD;
  vector<vector<int>> ctuTypeNum;
  vector<vector<int>> ctuType;
  vector<vector<double>> typeBAfactor;
  vector<vector<uint64_t>> regionalD;
  
  //0,text,1,SI,2,NI
  vector < vector<int>> cuflag;
  bool inPic(int x, int y, int xoffset, int yoffset)
  {
    if (x + xoffset >= 0 && x + xoffset < framew && y + yoffset >= 0 && y + yoffset < frameh)
      return true;
    else
      return false;
  }
  
  void transfrom(const CPelBuf &resi, CoeffBuf &dstCoeff, const int width, const int height);
//  {
//    int bitDepth = 10;
//    int maxLog2TrDynamicRange = 15;
//    const int      TRANSFORM_MATRIX_SHIFT = g_transformMatrixShift[TRANSFORM_FORWARD];
//    const uint32_t transformWidthIndex = floorLog2(width) - 1;  // nLog2WidthMinus1, since transform start from 2-point
//    const uint32_t transformHeightIndex = floorLog2(height) - 1;  // nLog2HeightMinus1, since transform start from 2-point
//
//    int trTypeHor = DCT2;
//    int trTypeVer = DCT2;
//    int  skipWidth = (trTypeHor != DCT2 && width == 32) ? 16 : width > JVET_C0024_ZERO_OUT_TH ? width - JVET_C0024_ZERO_OUT_TH : 0;
//    int  skipHeight = (trTypeVer != DCT2 && height == 32) ? 16 : height > JVET_C0024_ZERO_OUT_TH ? height - JVET_C0024_ZERO_OUT_TH : 0;
//    if (1)
//    {
//      if ((width == 4 && height > 4) || (width > 4 && height == 4))
//      {
//        skipWidth = width - 4;
//        skipHeight = height - 4;
//      }
//      else if ((width >= 8 && height >= 8))
//      {
//        skipWidth = width - 8;
//        skipHeight = height - 8;
//      }
//    }
//
//#if RExt__DECODER_DEBUG_TOOL_STATISTICS
//    if (trTypeHor != DCT2)
//    {
//      CodingStatistics::IncrementStatisticTool(CodingStatisticsClassType{ STATS__TOOL_EMT, uint32_t(width), uint32_t(height), compID });
//    }
//#endif
//
//    ALIGN_DATA(MEMORY_ALIGN_DEF_SIZE, TCoeff block[MAX_TB_SIZEY * MAX_TB_SIZEY]);
//
//    const Pel *resiBuf = resi.buf;
//    const int  resiStride = resi.stride;
//
//    for (int y = 0; y < height; y++)
//    {
//      for (int x = 0; x < width; x++)
//      {
//        block[(y * width) + x] = resiBuf[(y * resiStride) + x];
//      }
//    }
//
//    if (width > 1 && height > 1) // 2-D transform
//    {
//      const int      shift_1st = ((floorLog2(width)) + bitDepth + TRANSFORM_MATRIX_SHIFT) - maxLog2TrDynamicRange + COM16_C806_TRANS_PREC;
//      const int      shift_2nd = (floorLog2(height)) + TRANSFORM_MATRIX_SHIFT + COM16_C806_TRANS_PREC;
//      CHECK(shift_1st < 0, "Negative shift");
//      CHECK(shift_2nd < 0, "Negative shift");
//      TCoeff *tmp = (TCoeff *)alloca(width * height * sizeof(TCoeff));
//
//      fastFwdTrans1[trTypeHor][transformWidthIndex](block, tmp, shift_1st, height, 0, skipWidth);
//      fastFwdTrans1[trTypeVer][transformHeightIndex](tmp, dstCoeff.buf, shift_2nd, width, skipWidth, skipHeight);
//    }
//    else if (height == 1) //1-D horizontal transform
//    {
//      const int      shift = ((floorLog2(width)) + bitDepth + TRANSFORM_MATRIX_SHIFT) - maxLog2TrDynamicRange + COM16_C806_TRANS_PREC;
//      CHECK(shift < 0, "Negative shift");
//      CHECKD((transformWidthIndex < 0), "There is a problem with the width.");
//      fastFwdTrans1[trTypeHor][transformWidthIndex](block, dstCoeff.buf, shift, 1, 0, skipWidth);
//    }
//    else //if (iWidth == 1) //1-D vertical transform
//    {
//      int shift = ((floorLog2(height)) + bitDepth + TRANSFORM_MATRIX_SHIFT) - maxLog2TrDynamicRange + COM16_C806_TRANS_PREC;
//      CHECK(shift < 0, "Negative shift");
//      CHECKD((transformHeightIndex < 0), "There is a problem with the height.");
//      fastFwdTrans1[trTypeVer][transformHeightIndex](block, dstCoeff.buf, shift, 1, 0, skipHeight);
//    }
// 
//  }
  
  uint64_t  CalFD(int fidx, RPLEntry rpl1, RPLEntry rpl2)
  {



    
    vector<uint64_t> d = { 0,0,0,0 };
    // obtain ref list
    int numRefPics0 = rpl1.m_numRefPicsActive;
    int deltaRefPics0[MAX_NUM_REF_PICS] = { 0 };
    memcpy(deltaRefPics0, rpl1.m_deltaRefPics, numRefPics0 * sizeof(int));

    int numRefPics1 = rpl2.m_numRefPicsActive;
    int deltaRefPics1[MAX_NUM_REF_PICS] = { 0 };
    memcpy(deltaRefPics1, rpl2.m_deltaRefPics, numRefPics1 * sizeof(int));
    //printf("POC:%d\tL0:  ", fidx);
    //for (int i = 0; i < numRefPics0; i++)
    //{
    //  printf("%d  ", fidx- deltaRefPics0[i]);
    //}
    //printf("L1:  ");
    //for (int i = 0; i < numRefPics1; i++)
    //{
    //  printf("%d  ", fidx - deltaRefPics1[i]);
    //}
    //printf("\n");

    
        uint64_t dist = MAX_INT;
        uint64_t dist0 = 0;
        uint64_t dist1 = 0;


        int refidx = 0;
        // for L0
        if (numRefPics0 > 0)
        {
          
          int refpoc = fidx - deltaRefPics0[refidx];

          dist0 = m_pcRdCost->getDistPart(pre_ana_buf[fidx]->subBuf(0,0,framew,frameh), pre_ana_buf[refpoc]->subBuf(0, 0, framew, frameh), 10, COMPONENT_Y, costfun);
        }
        
        // for L1 if not LDP
        if (numRefPics1>0)
        {
          int refpoc = fidx - deltaRefPics1[refidx];

          dist1 = m_pcRdCost->getDistPart(pre_ana_buf[fidx]->subBuf(0, 0, framew, frameh), pre_ana_buf[refpoc]->subBuf(0, 0, framew, frameh), 10, COMPONENT_Y, costfun);
          
        }
        
        int ref0 = fidx - deltaRefPics0[refidx];
        int ref1 = fidx - deltaRefPics1[refidx];
        if (numRefPics0 > 0 && numRefPics1 > 0)
        {
          dist = dist0 / 2 + dist1 / 2;
          for (int ct = 0; ct < TOTAL; ct++)
          {
            
            double dd0 = regionalD[ref0][ct] / (double)(regionalD[ref0][TOTAL]+0.0001);
            double dd1 = regionalD[ref1][ct] / (double)(regionalD[ref1][TOTAL] + 0.0001);
            typeBAfactor[fidx][ct] = 0.5*dd0 + 0.5*dd1;
          }
        }
        else if (numRefPics0 == 0)
        {
          dist = dist1;
          for (int ct = 0; ct < TOTAL; ct++)
          {

            
            double dd1 = regionalD[ref1][ct] / (double)(regionalD[ref1][TOTAL] + 0.0001);
            typeBAfactor[fidx][ct] = dd1;
          }
        }
        else if (numRefPics1 == 0)
        {
          dist = dist0;
          for (int ct = 0; ct < TOTAL; ct++)
          {

            double dd0 = regionalD[ref0][ct] / (double)(regionalD[ref0][TOTAL] + 0.0001);
            
            typeBAfactor[fidx][ct] = dd0 ;
          }
        }
        

        


        return dist;




  }

  
  int check_scenechange(int fidx)
  {
    if (fd[fidx] / framew/frameh > 2500 * 16)
    {
      return 1;
    }
    if (pastD.size() != 0)
    {
      uint64_t avgD = 0;
      for (uint64_t x : pastD)
      {
        avgD += x;
      }
      
      if (fd[fidx] > 10 * avgD)
      {
        pastD.clear();
        return 1;
      }

    }
    //printf("pastD.size():%d\n", (int)pastD.size());
    return 0;
  }

  void CTU_classification(int fidx,int ctuidx)
  {
    

    
    Position CTUPos = getCTUPos(ctuidx);
    int xctu = CTUPos.x;
    int yctu = CTUPos.y;
    int IBCIdx = 0;
    
    ////////// T-CTU
    int w = min(cfgctusize, framew - xctu);
    int h = min(cfgctusize, frameh - yctu);

    //double a1 = 0;
    double a = 0;
    for (int tx = xctu; tx < xctu + w; tx++)
    {
      for (int ty = yctu; ty < yctu + h; ty++)
      {
        double ta1 = 0;
        double ta2 = 0;
        int a1idx = 0;
        int a2idx = 0;
        if (inPic(tx, ty, -1, -1))
        {
          ta1 += pow(pre_ana_buf[fidx]->at(tx, ty) - pre_ana_buf[fidx]->at(tx - 1, ty - 1), 2);
          a1idx++;
        }
        if (inPic(tx, ty, 1, -1))
        {
          ta1 += pow(pre_ana_buf[fidx]->at(tx, ty) - pre_ana_buf[fidx]->at(tx + 1, ty - 1), 2);
          a1idx++;
        }
        if (inPic(tx, ty, -1, 1))
        {
          ta1 += pow(pre_ana_buf[fidx]->at(tx, ty) - pre_ana_buf[fidx]->at(tx - 1, ty + 1), 2);
          a1idx++;
        }
        if (inPic(tx, ty, 1, 1))
        {
          ta1 += pow(pre_ana_buf[fidx]->at(tx, ty) - pre_ana_buf[fidx]->at(tx + 1, ty + 1), 2);
          a1idx++;
        }
        if (a1idx != 0)
        {
          ta1 = sqrt(ta1) / a1idx;
        }
        if (inPic(tx, ty, -1, -1) && inPic(tx, ty, 1, 1))
        {
          ta2 += pow(pre_ana_buf[fidx]->at(tx - 1, ty - 1) - pre_ana_buf[fidx]->at(tx + 1, ty + 1), 2);
          a2idx++;
        }
        if (inPic(tx, ty, 1, -1) && inPic(tx, ty, -1, 1))
        {
          ta2 += pow(pre_ana_buf[fidx]->at(tx + 1, ty - 1) - pre_ana_buf[fidx]->at(tx - 1, ty + 1), 2);
          a2idx++;
        }
        if (a2idx != 0)
        {
          ta2 = sqrt(ta2) / a2idx;
        }

        a = a + 0.5*ta1 + 0.5*ta2;
      }
    }
    a = a / w / h;
    if (a > 35)
    {
      cuflag[fidx][ctuidx] = TCTU;
      ctuTypeNum[fidx][TCTU] += 1;
      ctuType[fidx].push_back(TCTU);
      return;
    }
    uint64_t c10 = 0;
    uint64_t c64 = 0;
    for (int y = yctu; y < min(yctu + cfgctusize, frameh); y += CTUsize)
    {
      for (int x = xctu; x < min(framew, xctu + cfgctusize); x += CTUsize)
      {

        
        
        //////////
        
        int w1 = min(CTUsize, framew - x);
        int h1 = min(CTUsize, frameh - y);
        TCoeff* tbuf = (TCoeff*)xMalloc(TCoeff, w1*h1);
        memset(tbuf, 0, sizeof(TCoeff)*w1*h1);
        CoeffBuf tmp(tbuf, w1, h1);
        transfrom(pre_ana_buf[fidx]->subBuf(x, y, w1, h1), tmp, w1, h1);
        for (int ty = 0; ty < h1; ty++)
        {
          for (int tx = 0; tx < w1; tx++)
          {
            if (tx < 10 && ty < 10)
            {
              c10 += abs(tmp.at(tx, ty));
            }
            c64+= abs(tmp.at(tx, ty));
            //if(tmp.at(tx, ty)!=0)
              //printf("%d\t%d\t%d\n", tx, ty, tmp.at(tx, ty));
          }
        }
        
        xFree(tbuf);
      }
    }

    double s = c10 / double(c64);
    if (s <= 0.065)
    {
      cuflag[fidx][ctuidx] = SICTU;
      ctuTypeNum[fidx][SICTU] += 1;
      ctuType[fidx].push_back(SICTU);
    }
    else
    {
      cuflag[fidx][ctuidx] = NICTU;
      ctuTypeNum[fidx][NICTU] += 1;
      ctuType[fidx].push_back(NICTU);
    }
    
 
  }
#endif

#if wang2018frame
  vector<double> ifc;
  vector<int> iskey;
  vector<double> avgSATD;
  void setKey(int fidx)
  {
    if (fidx == 0)
    {
      ifc[fidx] = 0;
      iskey[fidx] = 1;
      return;
    }
    const int curCUSize = 16;
    const int totalCU = ((frameh + curCUSize/2) / curCUSize) * ((framew + curCUSize/2) / curCUSize);
    int sb = 0;
    for (int y = 0; y < frameh; y += curCUSize)
    {
    
      for (int x = 0; x < framew; x += curCUSize)
      {
        int w = min(curCUSize, framew - x);
        int h = min(curCUSize, frameh - y);
        Pel* tbuf = (Pel*)xMalloc(Pel, w*h);
        memset(tbuf, 0, sizeof(Pel)*w*h);
        PelBuf tmp(tbuf, w, h);
     
        auto dist = m_pcRdCost->getDistPart(pre_ana_buf[fidx]->subBuf(x, y, w, h), pre_ana_buf[fidx-1]->subBuf(x, y, w, h), 10, COMPONENT_Y, DF_SAD);
        if (dist < 2.5*w*h * 4)
        {
          sb += 1;
        }
      }
    }
    ifc[fidx] = sb / (double)totalCU;
    const double tao = 0.99;
    iskey[fidx] = ifc[fidx] < 0.99 ? 1 : 0;
  }
  void cptSATD(int fidx)
  {
    ///// only for LDB
    
    const int curCUSize = 32;
    const int totalCU = ((frameh + curCUSize / 2) / curCUSize) * ((framew + curCUSize / 2) / curCUSize);
    double satd = 0;
    for (int y = 0; y < frameh; y += curCUSize)
    {

      for (int x = 0; x < framew; x += curCUSize)
      {
        int w = min(curCUSize, framew - x);
        int h = min(curCUSize, frameh - y);
        Distortion dist = 0;
#if SATDwithCol
        if (fidx > 0)
        {        
          dist = m_pcRdCost->getDistPart(pre_ana_buf[fidx]->subBuf(x, y, w, h), pre_ana_buf[fidx - 1]->subBuf(x, y, w, h), 10, COMPONENT_Y, DF_HAD);
        }
        else 
        {
          Pel* tbuf = (Pel*)xMalloc(Pel, w*h);
          memset(tbuf, 0, sizeof(Pel)*w*h);
          PelBuf tmp(tbuf, w, h);
          tmp.fill(pre_ana_buf[fidx]->subBuf(x, y, w, h).computeAvg());
          dist = m_pcRdCost->getDistPart(pre_ana_buf[fidx]->subBuf(x, y, w, h), tmp, 10, COMPONENT_Y, DF_HAD);
        }
#endif
        satd += dist;
      }
    }
    avgSATD[fidx] = satd / frameh / framew;
  }

#endif

  pre_analysis() {
    m_size = 0; curidx = 0;
  }
  pre_analysis(int iframew, int iframeh) {
    m_size = 0; curidx = 0; framew = iframew; frameh = iframeh; 
  }
  //~pre_analysis();

  void init()
  {
    
    CTUsize = cfgctusize / 2;
    int numx = framew / cfgctusize + ((framew % cfgctusize) == 0 ? 0 : 1);
    int numy = frameh / cfgctusize + ((frameh % cfgctusize) == 0 ? 0 : 1);
    TotalCTUNum = numx * numy;
    
    


    IBCRef[0][0][0] = -cfgctusize/2;
    IBCRef[0][0][1] = 0;
    IBCRef[0][1][0] = -cfgctusize / 2;
    IBCRef[0][1][1] = cfgctusize / 2;
    IBCRef[0][2][0] = -cfgctusize;
    IBCRef[0][2][1] = cfgctusize / 2;

    IBCRef[1][0][0] = -cfgctusize / 2;
    IBCRef[1][0][1] = 0;
    IBCRef[1][1][0] = -cfgctusize;
    IBCRef[1][1][1] = cfgctusize / 2;
    IBCRef[1][2][0] = -3*cfgctusize / 2;
    IBCRef[1][2][1] = cfgctusize / 2;

    IBCRef[2][0][0] = -cfgctusize / 2;
    IBCRef[2][0][1] = 0;
    IBCRef[2][1][0] = 0;
    IBCRef[2][1][1] = -cfgctusize / 2;
    IBCRef[2][2][0] = cfgctusize / 2;
    IBCRef[2][2][1] = -cfgctusize / 2;

    IBCRef[3][0][0] = -cfgctusize / 2;
    IBCRef[3][0][1] = -cfgctusize / 2;
    IBCRef[3][1][0] = 0;
    IBCRef[3][1][1] = -cfgctusize / 2;
    IBCRef[3][2][0] = -cfgctusize / 2;
    IBCRef[3][2][1] = 0;
  }
  void updateCUSATD()
  {
    int cursize = (int)CUSATD.size();
    //FrameSATD.resize(cursize + m_size);
    CUSATD.resize(cursize + m_size);
    for (int i = 0; i < m_size; i++)
      CUSATD[i + cursize].resize(TotalCTUNum);
  }

  Position getCTUPos(int ctuidx)
  {
    int numx = framew / cfgctusize + ((framew % cfgctusize) == 0 ? 0 : 1);
    int numy = frameh / cfgctusize + ((frameh % cfgctusize) == 0 ? 0 : 1);
    int x = ctuidx % numx*cfgctusize;
    int y = ctuidx / numx * cfgctusize;
    return Position(x, y);
  }
  void createbuf(int w, int h)
  {
    for (int i = 0; i < m_size; i++)
    {
      //pre_ana_buf.push_back( new PelBuf);
      Pel *tbuf = new Pel[w*h];

      //PelBuf temp((Pel*)xMalloc(Pel, w*h), w, h);
      PelBuf temp(tbuf, w, h);
      pre_ana_buf.push_back(&temp);

    }
  }

  void addbuf(int n)
  {
    if (n > 0)
    {
      for (int i = 0; i < n; i++)
      {
        pre_ana_buf.push_back(new PelBuf);
      }
    }
  }

  template <class T>
  T Sum(const T x) {
    return x;
  }

  template <class T>
  T Sum(const std::vector<T> &v) {
    decltype(Sum(v[0])) sum = 0;
    for (const auto &x : v)
      sum += Sum(x);
    return sum;
  }


  void destroy()
  {
    while (pre_ana_buf.size() != 0)
    {
      delete pre_ana_buf.back()->buf;
      pre_ana_buf.back()->buf = NULL;
      delete pre_ana_buf.back();
      pre_ana_buf.back() = NULL;
      pre_ana_buf.pop_back();
    }

  }

  void updateFrameSATD(int startpos,int endpos)
  {
    for(int idx=startpos;idx<endpos;idx++)
    {
      for (int i = 0; i < CUSATD[idx].size(); i++)
      {
        for (int j = 0; j < CUSATD[idx][i].size(); j++)
        {
          FrameSATD[idx] += CUSATD[idx][i][j];
        }
      }
      
    }
  }

  void clearSATD()
  {
    CUSATD.clear();
    FrameSATD.clear();

  }


  vector<uint64_t> CalIntraSATD(int fidx,int ctuidx)
  {
    Pel* tbuf = (Pel*)xMalloc(Pel, CTUsize*CTUsize);

    uint64_t dist;
    Position CTUPos = getCTUPos(ctuidx);
    int xctu = CTUPos.x;
    int yctu = CTUPos.y;
    int IBCIdx = 0;
    vector<uint64_t> d = { 0,0,0,0 };
    for (int y = yctu; y < min(yctu + cfgctusize, frameh); y += CTUsize)
    {
      for (int x = xctu; x < min(framew, xctu + cfgctusize); x += CTUsize)
      {
        int w = min(CTUsize, framew - x);
        int h = min(CTUsize, frameh - y);
        Area CurCU = Area(x, y, w, h);
        if (x == xctu && y == yctu)
          IBCIdx = 0;
        else if (x == xctu + CTUsize && y == yctu)
          IBCIdx = 1;
        else if (x == xctu && y == yctu + CTUsize)
          IBCIdx = 2;
        else if (x == xctu + CTUsize && y == yctu + CTUsize)
          IBCIdx = 3;
        PelBuf tmp(tbuf, w, h);
        tmp.fill(pre_ana_buf[fidx]->subBuf(x, y, w, h).computeAvg());
        dist = m_pcRdCost->getDistPart(pre_ana_buf[fidx]->subBuf(x, y, w, h), tmp, 10, COMPONENT_Y, costfun);
        d[IBCIdx] = dist;
        //CUSATD[fidx][ctuidx].push_back(dist);
        //auto dist = RdCost::getDistPart(pre_ana_buf[idx]->subBuf(x, y, CTUsize, CTUsize), tmp, 10, COMPONENT_Y, DF_SSE);
      }
    }
    
    xFree(tbuf);
    return d;
  }

  vector<uint64_t>  CalIBCSATD(int fidx,int ctuidx)
  {
    Pel* tbuf = (Pel*)xMalloc(Pel, CTUsize*CTUsize);


    int IBCIdx = 0;
    Position CTUPos = getCTUPos(ctuidx);
    int xctu = CTUPos.x;
    int yctu = CTUPos.y;
    vector<uint64_t> d={0,0,0,0};
        for (int y = yctu; y < min(yctu + cfgctusize, frameh); y += CTUsize)
        {
          for (int x = xctu; x < min(framew, xctu + cfgctusize); x += CTUsize)
          {

            int w = min(CTUsize, framew - x);
            int h = min(CTUsize, frameh - y);
            Area CurCU = Area(x, y, w, h);
            PelBuf tmp(tbuf, w, h);
            if (x == xctu && y == yctu)
              IBCIdx = 0;
            else if (x == xctu + CTUsize && y == yctu)
              IBCIdx = 1;
            else if (x == xctu && y == yctu + CTUsize)
              IBCIdx = 2;
            else if (x == xctu + CTUsize && y == yctu + CTUsize)
              IBCIdx = 3;

            if (xctu == 384 && yctu == 0)
              int x = 0;

            uint64_t dist = MAX_INT;
            for (int refidx = 0; refidx < 3; refidx++)
            {
              if (x + IBCRef[IBCIdx][refidx][0] >= 0 && y + IBCRef[IBCIdx][refidx][1] >= 0)
              {
                int refx = x + IBCRef[IBCIdx][refidx][0];
                int refy = y + IBCRef[IBCIdx][refidx][1];
                int refw= min(min(CTUsize, framew - refx), w);
                int refh= min(min(CTUsize, frameh - refy),h);
                if (refw >= w && refh >= h)
                {
                  Area RefCU = Area(refx, refy, refw, refh);


                  //tmp.copyFrom(pre_ana_buf[fidx]->subBuf(x + IBCRef[IBCIdx][refidx][0], y + IBCRef[IBCIdx][refidx][1], w, h));
                  //dist = min(dist, m_pcRdCost->getDistPart(pre_ana_buf[fidx]->subBuf(x, y, w, h), tmp, 10, COMPONENT_Y, costfun));
                  tmp.copyFrom(pre_ana_buf[fidx]->subBuf(refx, refy, refw, refh));
                  dist = min(dist, m_pcRdCost->getDistPart(pre_ana_buf[fidx]->subBuf(x, y, refw, refh), tmp, 10, COMPONENT_Y, costfun));
                }
                
              }
            }
            d[IBCIdx] = dist;
            //if (CUSATD[fidx][ctuidx].size() < IBCIdx)
            //{
            //  CUSATD[fidx][ctuidx].push_back(0);
            //}
            //CUSATD[fidx][ctuidx].push_back(dist);
          }
        }

        xFree(tbuf);
    return d;
  }

  vector<uint64_t>  CalInterSATD(int fidx,int ctuidx, RPLEntry rpl1, RPLEntry rpl2)
  {
    

    
    Position CTUPos = getCTUPos(ctuidx);
    int xctu = CTUPos.x;
    int yctu = CTUPos.y;
    int IBCIdx = 0;
    vector<uint64_t> d = { 0,0,0,0 };
    // obtain ref list
    int numRefPics0 = rpl1.m_numRefPicsActive;
    int deltaRefPics0[MAX_NUM_REF_PICS] = {0};
    memcpy(deltaRefPics0, rpl1.m_deltaRefPics, numRefPics0 * sizeof(int));

    int numRefPics1 = rpl2.m_numRefPicsActive;
    int deltaRefPics1[MAX_NUM_REF_PICS] = { 0 };
    memcpy(deltaRefPics1, rpl2.m_deltaRefPics, numRefPics1 * sizeof(int));
    //printf("POC:%d\tL0:  ", fidx);
    //for (int i = 0; i < numRefPics0; i++)
    //{
    //  printf("%d  ", fidx- deltaRefPics0[i]);
    //}
    //printf("L1:  ");
    //for (int i = 0; i < numRefPics1; i++)
    //{
    //  printf("%d  ", fidx - deltaRefPics1[i]);
    //}
    //printf("\n");
    
    for (int y = yctu; y < min(yctu + cfgctusize, frameh); y += CTUsize)
    {
      for (int x = xctu; x < min(framew, xctu + cfgctusize); x += CTUsize)
      {
        uint64_t dist = MAX_INT;
        int w = min(CTUsize, framew - x);
        int h = min(CTUsize, frameh - y);
        Area CurCU = Area(x, y, w, h);
        if (x == xctu && y == yctu)
          IBCIdx = 0;
        else if (x == xctu + CTUsize && y == yctu)
          IBCIdx = 1;
        else if (x == xctu && y == yctu + CTUsize)
          IBCIdx = 2;
        else if (x == xctu + CTUsize && y == yctu + CTUsize)
          IBCIdx = 3;
        

        // for L0
        for (int refidx = 0; refidx < numRefPics0; refidx++)
        {
          int refpoc = fidx - deltaRefPics0[refidx];
          //if(refpoc>=0)
          //{
            dist = min(dist, m_pcRdCost->getDistPart(pre_ana_buf[fidx]->subBuf(x, y, w, h), pre_ana_buf[refpoc]->subBuf(x, y, w, h), 10, COMPONENT_Y, costfun));
          //}
          //else
          //{
            //dist = min(dist, m_pcRdCost->getDistPart(pre_ana_buf[fidx]->subBuf(x, y, w, h), pre_ana_buf[refpoc]->subBuf(x, y, w, h), 10, COMPONENT_Y, costfun));
          //}
            //if (IBCIdx == 0 && x == 0 && y == 0 && refidx==0)
            //{
            //  printf("%d->%d:  %llu\n", fidx, refpoc,dist);
            //}
        }
        // for L1 if not LDP
        for (int refidx = 0; refidx < numRefPics1; refidx++)
        {
          int refpoc = fidx - deltaRefPics1[refidx];

          dist = min(dist, m_pcRdCost->getDistPart(pre_ana_buf[fidx]->subBuf(x, y, w, h), pre_ana_buf[refpoc]->subBuf(x, y, w, h), 10, COMPONENT_Y, costfun));

        }
        // for combine if not LDP}
        if (numRefPics1 != 0 && numRefPics0 != 0)
        {
          int refpoc0 = fidx - deltaRefPics0[0];
          int refpoc1 = fidx - deltaRefPics1[0];
          if (refpoc0 != refpoc1)
          {
            //Pel* tbuf = (Pel*)xMalloc(Pel, CTUsize*CTUsize);
            Pel* tbuf = (Pel*)xMalloc(Pel, w*h);
            PelBuf tmp(tbuf, w, h);
            int BD = 10;
            ClpRng t=ClpRng();
            t.min = 0;
            t.max = (1 << BD) - 1;
            t.bd = BD;
            t.n = 0; 
           
            tmp.addAvg(pre_ana_buf[refpoc0]->subBuf(x, y, w, h), pre_ana_buf[refpoc1]->subBuf(x, y, w, h),t);
          }
        }



        d[IBCIdx] = dist;
      }
    }

    
    
    return d;
  }


};



#endif

#if UsePipe
#if iswindows
const double paraC = 4.164 * 16;
const double paraK = 0.429;
const int SlideWindowSize = 20;
const double r_beta = 5;
const double C0 = 1;
class usingpipe
{
public:
  string pip_name;
  char ReadBuf[PBUFSIZE];
  char SendBuf[PBUFSIZE];
  HANDLE h_Pipe;
  vector<double> action = { 0,0 }; // 0, QP, 1, lambda
  vector<double> state_and_reward = { 0,0 };
  double reward = 0;

  double lastLambda=0;
  double lastbpp=0;
  double lastpsnr=0;
  double ExpFrameBufLevel = 0;
  double ExpGOPBufLevel = 0;
  double sumbits=0;
  double lastmse = 0;
  double refmse = 0;
  double GOPid=0;
  usingpipe() {};

  void init(string pipdir) { pip_name.empty(); pip_name += pipdir; };
  bool isexist() { return WaitNamedPipe(TEXT(pip_name.c_str()), NMPWAIT_WAIT_FOREVER); }
  bool create_handle()
  {
    h_Pipe = CreateFile(						//管道属于一种特殊的文件
      TEXT(pip_name.c_str()),				//文件名字
      GENERIC_READ | GENERIC_WRITE,				//文件模式
      //0,											//是否共享
      FILE_SHARE_READ | FILE_SHARE_WRITE,
      NULL,
      OPEN_EXISTING,
      FILE_ATTRIBUTE_NORMAL,						//文件属性(只读，默认...)NORMAL 为默认属性
      //FILE_FLAG_OVERLAPPED,
      NULL);
    return h_Pipe != INVALID_HANDLE_VALUE;
  }
  bool write() { return WriteFile(h_Pipe, SendBuf, sizeof(SendBuf), 0, 0); }
  bool read() { return ReadFile(h_Pipe, ReadBuf, sizeof(ReadBuf), 0, 0); }
  void encode()
  {
    std::stringstream buf;
    buf.precision(4);//覆盖默认精度
    //buf.setf(std::ios::fixed);//保留小数位
    buf.setf(std::ios::scientific);
    for (int i = 0; i < state_and_reward.size(); i++)
    {
      buf << state_and_reward[i] << " ";
    }

    std::string str1;
    str1 = buf.str();
    strcpy(SendBuf, str1.c_str());
  }
  void decode()
  {
    char* d = " ";
    char *p = strtok(ReadBuf, d);
    int i = 0;
    while (p) {
      string s = p; //分割得到的字符串转换为string类型
      action[i] = atoi(s.c_str()); //存入结果数组
      p = strtok(NULL, d);
      i += 1;
    }
  }

};
#endif
#endif



#endif // __ENCMB__
