//:>----------------------------------------------------------------------
//: FILE:      gl3TriggerDecider.h
//: HISTORY:
//:>----------------------------------------------------------------------
#ifndef GL3TRIGGERDECIDER
#define GL3TRIGGERDECIDER

#include <stdio.h>	// for FILE ops
#include "HLTFormats.h"


class online_tracking_collector;
class l3_algos_t;

class gl3TriggerDecider {

 public:
   gl3TriggerDecider(online_tracking_collector* _event, char* HLTparameters = "HLTparameters");
   ~gl3TriggerDecider () {};
   void setTriggers(char* HLTparameters = "HLTparameters");
   void setQA(char* outputFileName);
   void readHLTparameters(char* HLTparameters = "HLTparameters");
   void setTriggersFromRunControl(l3_algos_t* alg);
   int decide(int eventId = 0);
   void writeQA();
   void writeScalers();
   void flushQA() ;	// tonko added
   void closeQA() ;	// tonko added

   double bField, scalerCount;

   int triggerOnHighPt, triggerOnDiElectron, triggerOnHeavyFragment, triggerOnAllEvents, triggerOnRandomEvents, triggerOnBesGoodEvents;

   struct HLT_DIEP hlt_diEP;
   struct HLT_HIPT hlt_hiPt;
   struct HLT_HF hlt_hF;

 private:
   FILE *f1 ;	// tonko: moved here from static in cxx
   online_tracking_collector* event;

   int triggerBitHighPt, triggerBitDiElectron, triggerBitHeavyFragment, triggerBitAllEvents, triggerBitRandomEvents, triggerBitBesGoodEvents;
   double innerSectorGain, outerSectorGain, spaceChargeP0, spaceChargeP1;
   double paraXVertex, paraYVertex;
   double sigmaDedx1, dPOverP2, sigmaTof;

   //high pt trigger parameters
   int nHitsCut_highPt;
   double ptCut_highPt;
   
   //diElectron trigger parameters
   int nHitsCut_diElectron;
   int nDedxCut_diElectron;
   double towerEnergyCut_diElectron;
   double pCutLowForTof_diElectron;
   double pCutHighForTof_diElectron;
   double pCutLowForEmc_diElectron;
   double pCutHighForEmc_diElectron;
   double pCutLowForTofAndEmc_diElectron;
   double pCutHighForTofAndEmc_diElectron;
   double nSigmaDedxCutLowForTof_diElectron;
   double nSigmaDedxCutHighForTof_diElectron;
   double nSigmaDedxCutLowForEmc_diElectron;
   double nSigmaDedxCutHighForEmc_diElectron;
   double nSigmaDedxCutLowForTofAndEmc_diElectron;
   double nSigmaDedxCutHighForTofAndEmc_diElectron;
   double p1Cut_diElectron;
   double p2Cut_diElectron;
   double pOverECutLow_diElectron;
   double pOverECutHigh_diElectron;
   double invariantMassCutLow_diElectron;
   double invariantMassCutHigh_diElectron;
   double vertexZDiffCut_diElectron;
   double oneOverBetaCut_diElectron;
   int onlyOppPairs_diElectron;

   //heavyFragment trigger parameters
   int nHitsCut_heavyFragment;
   int nDedxCut_heavyFragment;
   int useTofMatchedGlobalTrack_heavyfragMent;
   double nSigmaDedxHe3Cut_heavyFragment;
   double nSigmaDedxTritonCut_heavyFragment;
   double nSigmaMassTritonCut_heavyFragment;
   double dcaCut_heavyFragment;

   //random trigger parameter
   int sampleScale_randomEvents;

   //debug
   int Debug_allTracks;
   int Debug_tracks;
   int Debug_towers;
   int Debug_matchs;
   int Debug_vertex;

   int nHitsCut_debug;
   int nDedxCut_debug;
   double matchZEdgeCut_debug;
   double matchPhiDiffCut_debug;
   double dedxCutLow_debug;
   double dedxCutHigh_debug;
   double pCut_debug;
   double towerEnergyCut_debug;
   
};
#endif
