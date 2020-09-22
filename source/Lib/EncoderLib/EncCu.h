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
int getRPLIdxLDB(int poc);
int getRPLIdxRA(int poc);
class pre_analysis
{
public:
  int m_size;
  vector<PelBuf*> pre_ana_buf;
  vector<vector<vector<uint64_t>>> CUSATD;
  vector<uint64_t> FrameSATD;
  int curidx;
  int CTUsize = 64;
  int framew;
  int frameh;
  RdCost *m_pcRdCost;
  //EncCu *m_pcEncCu;
  //EncLib *m_pclib;
  int TotalCTUNum;
  //RPLEntry RefList0[MAX_GOP];
  int IBCRef[4][3][2];

  pre_analysis(int iframew, int iframeh) {
    m_size = 0; curidx = 0; framew = iframew; frameh = iframeh; 
  }
  //~pre_analysis();

  void init()
  {
    FrameSATD.resize(m_size);
    CUSATD.resize(m_size);

    int numx = framew / 128 + ((framew % 128) == 0 ? 0 : 1);
    int numy = frameh / 128 + ((frameh % 128) == 0 ? 0 : 1);
    TotalCTUNum = numx * numy;
    
    for(int i=0;i< m_size;i++)
      CUSATD[i].resize(TotalCTUNum);


    IBCRef[0][0][0] = -64;
    IBCRef[0][0][1] = 0;
    IBCRef[0][1][0] = -64;
    IBCRef[0][1][1] = 64;
    IBCRef[0][2][0] = -128;
    IBCRef[0][2][1] = 64;

    IBCRef[1][0][0] = -64;
    IBCRef[1][0][1] = 0;
    IBCRef[1][1][0] = -128;
    IBCRef[1][1][1] = 64;
    IBCRef[1][2][0] = -192;
    IBCRef[1][2][1] = 64;

    IBCRef[2][0][0] = -64;
    IBCRef[2][0][1] = 0;
    IBCRef[2][1][0] = 0;
    IBCRef[2][1][1] = 64;
    IBCRef[2][2][0] = 64;
    IBCRef[2][2][1] = -64;

    IBCRef[3][0][0] = -64;
    IBCRef[3][0][1] = -64;
    IBCRef[3][1][0] = 0;
    IBCRef[3][1][1] = -64;
    IBCRef[3][2][0] = -64;
    IBCRef[3][2][1] = 0;
  }

  Position getCTUPos(int ctuidx)
  {
    int numx = framew / 128 + ((framew % 128) == 0 ? 0 : 1);
    int numy = frameh / 128 + ((frameh % 128) == 0 ? 0 : 1);
    int x = ctuidx % numx*128;
    int y = ctuidx / numx * 128;
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

  


  vector<uint64_t> CalIntraSATD(int fidx,int ctuidx)
  {
    Pel* tbuf = (Pel*)xMalloc(Pel, CTUsize*CTUsize);

    uint64_t dist;
    Position CTUPos = getCTUPos(ctuidx);
    int xctu = CTUPos.x;
    int yctu = CTUPos.y;
    int IBCIdx = 0;
    vector<uint64_t> d = { 0,0,0,0 };
        for (int y = yctu; y < min(yctu + 128, frameh); y += CTUsize)
        {
          for (int x = xctu; x < min(framew, xctu + 128); x += CTUsize)
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
            dist = m_pcRdCost->getDistPart(pre_ana_buf[fidx]->subBuf(x, y, w, h), tmp, 10, COMPONENT_Y, DF_HAD);
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
        for (int y = yctu; y < min(yctu + 128, frameh); y += CTUsize)
        {
          for (int x = xctu; x < min(framew, xctu + 128); x += CTUsize)
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
                Area RefCU = Area(x + IBCRef[IBCIdx][refidx][0], y + IBCRef[IBCIdx][refidx][1], min(CTUsize, framew - x), min(CTUsize, frameh - y));


                tmp.copyFrom(pre_ana_buf[fidx]->subBuf(x + IBCRef[IBCIdx][refidx][0], y + IBCRef[IBCIdx][refidx][1], w, h));
                dist = min(dist, m_pcRdCost->getDistPart(pre_ana_buf[fidx]->subBuf(x, y, w, h), tmp, 10, COMPONENT_Y, DF_HAD));

                
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

    
    for (int y = yctu; y < min(yctu + 128, frameh); y += CTUsize)
    {
      for (int x = xctu; x < min(framew, xctu + 128); x += CTUsize)
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

          dist = min(dist, m_pcRdCost->getDistPart(pre_ana_buf[fidx]->subBuf(x, y, w, h), pre_ana_buf[refpoc]->subBuf(x, y, w, h), 10, COMPONENT_Y, DF_HAD));
          
        }
        // for L1 if not LDP
        for (int refidx = 0; refidx < numRefPics1; refidx++)
        {
          int refpoc = fidx - deltaRefPics1[refidx];

          dist = min(dist, m_pcRdCost->getDistPart(pre_ana_buf[fidx]->subBuf(x, y, w, h), pre_ana_buf[refpoc]->subBuf(x, y, w, h), 10, COMPONENT_Y, DF_HAD));

        }
        // for combine if not LDP}
        if (numRefPics1 != 0 && numRefPics0 != 0)
        {
          int refpoc0 = fidx - deltaRefPics0[0];
          int refpoc1 = fidx - deltaRefPics1[0];
          if (refpoc0 != refpoc1)
          {
            Pel* tbuf = (Pel*)xMalloc(Pel, CTUsize*CTUsize);
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

  //int CalRPLIdx(int POCCurr)
  //{
  //  
  //  int m_iGOPSize = m_pclib->getGOPSize();
  //  int fullListNum = m_iGOPSize;
  //  int partialListNum = m_pclib->getRPLCandidateSize(0) - m_iGOPSize;
  //  int extraNum = fullListNum;

  //  int rplPeriod = m_pclib->getIntraPeriod();
  //  if (rplPeriod < 0)  //Need to check if it is low delay or RA but with no RAP
  //  {
  //    //if (slice->getSPS()->getRPLList0()->getReferencePictureList(1)->getRefPicIdentifier(0) * slice->getSPS()->getRPLList1()->getReferencePictureList(1)->getRefPicIdentifier(0) < 0)
  //    if (m_pclib->getSPS(0)->getRPLList0()->getReferencePictureList(1)->getRefPicIdentifier(0) * m_pclib->getSPS(0)->getRPLList1()->getReferencePictureList(1)->getRefPicIdentifier(0) < 0)
  //    {
  //      rplPeriod = m_iGOPSize * 2;
  //      
  //    }
  //  }

  //  if (rplPeriod < 0)
  //  {
  //    // first 18 frames in seq
  //    if (POCCurr < (2 * m_iGOPSize + 2))
  //    {
  //      //int candidateIdx = (POCCurr + m_iGOPSize - 1 >= fullListNum + partialListNum) ? GOPid : POCCurr + m_iGOPSize - 1;
  //      //return candidateIdx;
  //      return POCCurr+7;
  //    }
  //    else
  //    {
  //      return ((POCCurr%m_iGOPSize == 0) ? m_iGOPSize - 1 : POCCurr % m_iGOPSize - 1);
  //      
  //    }
  //    extraNum = fullListNum + partialListNum;
  //  }
  //  for (; extraNum < fullListNum + partialListNum; extraNum++)
  //  {
  //    if (rplPeriod > 0)
  //    {
  //      int POCIndex = POCCurr % rplPeriod;
  //      if (POCIndex == 0)
  //      {
  //        POCIndex = rplPeriod;
  //      }
  //      if (POCIndex == m_pclib->getRPLEntry(1,extraNum).m_POC)
  //      {
  //        
  //        return extraNum;
  //      }
  //    }
  //  }

  //  
  //}
};



#endif

#endif // __ENCMB__
