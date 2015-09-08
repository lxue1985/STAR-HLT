/*:>-------------------------------------------------------------------
**: FILE:     gl3TriggerDecider.cxx
**: HISTORY:
**:<------------------------------------------------------------------*/

#include <iostream>
#include <string.h>
#include <fstream>
#include "gl3TriggerDecider.h"
#include <stdio.h>	// for FILE ops
#include "online_tracking_collector.h"
#include "gl3Bischel.h"
#include "gl3TOF.h"
#include "gl3Node.h"
#include "../L3_SUPPORT/l3_support.h"

const int maxNElectronTracks = 1000;
const int maxNTowers = 4800;

  

//***********************************************************************
// constructor
//***********************************************************************
gl3TriggerDecider::gl3TriggerDecider(online_tracking_collector* _event, char* HLTparameters)
{
  event = _event;
  setTriggers(HLTparameters);
  f1=0;
}

void gl3TriggerDecider::flushQA()
{
	if(f1) fflush(f1) ;
}

void gl3TriggerDecider::closeQA()
{
	if(f1) {
		fclose(f1) ;
		f1 = 0 ;
	}
}



void gl3TriggerDecider::setTriggers(char* HLTparameters)
{
  //set which triggers to enable here
  triggerOnHighPt = 0;
  triggerOnDiElectron = 0;
  triggerOnHeavyFragment = 0;
  triggerOnAllEvents = 0;
  triggerOnRandomEvents = 0;
  triggerOnBesGoodEvents = 0;
  
  triggerBitHighPt = 0x10000;
  triggerBitDiElectron = 0x20000;
  triggerBitHeavyFragment = 0x40000;
  triggerBitAllEvents = 0x80000;
  triggerBitRandomEvents = 0x100000;
  triggerBitBesGoodEvents = 0x200000;

   //high pt trigger parameters
   nHitsCut_highPt = 20;
   ptCut_highPt = 5.;
   
   //diElectron trigger parameters
   nHitsCut_diElectron = 15;
   nDedxCut_diElectron = 5;
   towerEnergyCut_diElectron = 0.5;
   pCutLowForTof_diElectron = 0;
   pCutHighForTof_diElectron = 0;
   pCutLowForEmc_diElectron = 0.5;
   pCutHighForEmc_diElectron = 1.e10;
   pCutLowForTofAndEmc_diElectron = 0;
   pCutHighForTofAndEmc_diElectron = 1.e10;
   nSigmaDedxCutLowForTof_diElectron = 0.;
   nSigmaDedxCutHighForTof_diElectron = 3.;
   nSigmaDedxCutLowForEmc_diElectron = -1.;
   nSigmaDedxCutHighForEmc_diElectron = 3.;
   nSigmaDedxCutLowForTofAndEmc_diElectron = -3.;
   nSigmaDedxCutHighForTofAndEmc_diElectron = 3.;
   p1Cut_diElectron = 2.;
   p2Cut_diElectron = 0.5;
   pOverECutLow_diElectron = 0.3;
   pOverECutHigh_diElectron = 1.5;
   invariantMassCutLow_diElectron = 0.;
   invariantMassCutHigh_diElectron = 1.e6;
   vertexZDiffCut_diElectron = 20.;
   oneOverBetaCut_diElectron = 0.03;
   onlyOppPairs_diElectron = 0;

   //heavyFragment trigger parameters
   nHitsCut_heavyFragment = 0;
   nDedxCut_heavyFragment = 0;
   useTofMatchedGlobalTrack_heavyfragMent = 0;
   nSigmaDedxHe3Cut_heavyFragment = 0.;
   nSigmaDedxTritonCut_heavyFragment = 0.;
   nSigmaMassTritonCut_heavyFragment = 0.;
   dcaCut_heavyFragment = 0.;

   sampleScale_randomEvents = 997;

   //debug parameters
   Debug_allTracks = 0;
   Debug_tracks = 0;
   Debug_towers = 0;
   Debug_matchs = 0;
   Debug_vertex = 0;
   
   nHitsCut_debug = 10;
   nDedxCut_debug = 0;
   dedxCutLow_debug = 0.e-6;
   dedxCutHigh_debug = 1.;
   pCut_debug = 0.;
   towerEnergyCut_debug = 0.5;
   matchPhiDiffCut_debug = 0.2;
   matchZEdgeCut_debug = 10;
   
   //others
//   innerSectorGain = 0.;
//   outerSectorGain = 0.;
   spaceChargeP0 = 0.;
   spaceChargeP1 = 0.;
   paraXVertex = 0.;
   paraYVertex = 0.;
   bField = 0.;
   scalerCount = 0.;
   sigmaDedx1 = 0.;
   dPOverP2 = 0.;
   sigmaTof = 0.;

   readHLTparameters(HLTparameters);
}

void gl3TriggerDecider::setTriggersFromRunControl(l3_algos_t* alg)
{
  if(strcmp(alg->name, "highPt") == 0)
    {
      nHitsCut_highPt = alg->rc.lx.userInt[0];
      ptCut_highPt = alg->rc.lx.userFloat[0];
    }
  else if(strcmp(alg->name, "diElectron") == 0)
    {
      nHitsCut_diElectron = alg->rc.lx.userInt[0];
      nDedxCut_diElectron = alg->rc.lx.userInt[1];
      nSigmaDedxCutLowForTof_diElectron = alg->rc.lx.userFloat[0];
      nSigmaDedxCutLowForEmc_diElectron = alg->rc.lx.userFloat[1];
      nSigmaDedxCutLowForTofAndEmc_diElectron = alg->rc.lx.userFloat[2];
      p2Cut_diElectron = alg->rc.lx.userFloat[3];
      oneOverBetaCut_diElectron = alg->rc.lx.userFloat[4];
    }
  else if(strcmp(alg->name, "heavyFragment") == 0)
    {
      nHitsCut_heavyFragment = alg->rc.lx.userInt[0];
      nDedxCut_heavyFragment = alg->rc.lx.userInt[1];
      nSigmaDedxHe3Cut_heavyFragment = alg->rc.lx.userFloat[0];
      nSigmaDedxTritonCut_heavyFragment = alg->rc.lx.userFloat[1];
      nSigmaMassTritonCut_heavyFragment = alg->rc.lx.userFloat[2];
      dcaCut_heavyFragment = alg->rc.lx.userFloat[3];
    }
  else if(strcmp(alg->name, "allEvents_with_HLT") == 0)
    {
    }
  else if(strcmp(alg->name, "randomEvents") == 0)
    {
      sampleScale_randomEvents = alg->rc.lx.userInt[0];
    }
  else if(strcmp(alg->name, "besGoodEvents") == 0)
    {
    }
  else if(strcmp(alg->name, "diV0") == 0)
    {
    }
  else if(strcmp(alg->name, "allEvents_no_HLT") == 0)
    {
    }
  else
    {
      LOG(ERR,"HLT algorithm %s not found", alg->name);      
    }
}

void gl3TriggerDecider::setQA(char* outputFileName)
{
  f1 = fopen(outputFileName, "w");

  fprintf(f1, "Settings and cuts:\n");
  fprintf(f1, "triggerOnHighPt = %i\n", triggerOnHighPt);
  fprintf(f1, "triggerOnDiElectron = %i\n", triggerOnDiElectron);
  fprintf(f1, "triggerOnHeavyFragment = %i\n", triggerOnHeavyFragment);
  fprintf(f1, "triggerOnAllEvents = %i\n", triggerOnAllEvents);
  fprintf(f1, "triggerOnRandomEvents = %i\n", triggerOnRandomEvents);
  fprintf(f1, "triggerOnBesGoodEvents = %i\n", triggerOnBesGoodEvents);
  fprintf(f1, "\n");
  fprintf(f1, "triggerBitHighPt = 0x%x\n", triggerBitHighPt);
  fprintf(f1, "triggerBitDiElectron = 0x%x\n", triggerBitDiElectron);
  fprintf(f1, "triggerBitHeavyFragment = 0x%x\n", triggerBitHeavyFragment);
  fprintf(f1, "triggerBitAllEvents = 0x%x\n", triggerBitAllEvents);
  fprintf(f1, "triggerBitRandomEvents = 0x%x\n", triggerBitRandomEvents);
  fprintf(f1, "triggerBitBesGoodEvents = 0x%x\n", triggerBitBesGoodEvents);
  fprintf(f1, "\n");
  fprintf(f1, "nHitsCut_highPt = %i\n", nHitsCut_highPt);
  fprintf(f1, "ptCut_highPt = %e\n", ptCut_highPt);
  fprintf(f1, "\n");
  fprintf(f1, "nHitsCut_heavyFragment = %i\n", nHitsCut_heavyFragment);
  fprintf(f1, "nDedxCut_heavyFragment = %i\n", nDedxCut_heavyFragment);
  fprintf(f1, "useTofMatchedGlobalTrack_heavyfragMent = %i\n", useTofMatchedGlobalTrack_heavyfragMent);
  fprintf(f1, "nSigmaDedxHe3Cut_heavyFragment = %e\n", nSigmaDedxHe3Cut_heavyFragment);
  fprintf(f1, "nSigmaDedxTritonCut_heavyFragment = %e\n", nSigmaDedxTritonCut_heavyFragment);
  fprintf(f1, "nSigmaMassTritonCut_heavyFragment = %e\n", nSigmaMassTritonCut_heavyFragment);
  fprintf(f1, "dcaCut_heavyFragment = %e\n", dcaCut_heavyFragment);
  fprintf(f1, "\n");
  fprintf(f1, "nHitsCut_diElectron = %i\n", nHitsCut_diElectron);
  fprintf(f1, "nDedxCut_diElectron = %i\n", nDedxCut_diElectron);
  fprintf(f1, "towerEnergyCut_diElectron = %e\n", towerEnergyCut_diElectron);
  fprintf(f1, "pCutLowForTof_diElectron = %e\n", pCutLowForTof_diElectron);
  fprintf(f1, "pCutHighForTof_diElectron = %e\n", pCutHighForTof_diElectron);
  fprintf(f1, "pCutLowForEmc_diElectron = %e\n", pCutLowForEmc_diElectron);
  fprintf(f1, "pCutHighForEmc_diElectron = %e\n", pCutHighForEmc_diElectron);
  fprintf(f1, "pCutLowForTofAndEmc_diElectron = %e\n", pCutLowForTofAndEmc_diElectron);
  fprintf(f1, "pCutHighForTofAndEmc_diElectron = %e\n", pCutHighForTofAndEmc_diElectron);
  fprintf(f1, "nSigmaDedxCutLowForTof_diElectron = %e\n", nSigmaDedxCutLowForTof_diElectron);
  fprintf(f1, "nSigmaDedxCutHighForTof_diElectron = %e\n", nSigmaDedxCutHighForTof_diElectron);
  fprintf(f1, "nSigmaDedxCutLowForEmc_diElectron = %e\n", nSigmaDedxCutLowForEmc_diElectron);
  fprintf(f1, "nSigmaDedxCutHighForEmc_diElectron = %e\n", nSigmaDedxCutHighForEmc_diElectron);
  fprintf(f1, "nSigmaDedxCutLowForTofAndEmc_diElectron = %e\n", nSigmaDedxCutLowForTofAndEmc_diElectron);
  fprintf(f1, "nSigmaDedxCutHighForTofAndEmc_diElectron = %e\n", nSigmaDedxCutHighForTofAndEmc_diElectron);
  fprintf(f1, "p1Cut_diElectron = %e\n", p1Cut_diElectron);
  fprintf(f1, "p2Cut_diElectron = %e\n", p2Cut_diElectron);
  fprintf(f1, "pOverECutLow_diElectron = %e\n", pOverECutLow_diElectron);
  fprintf(f1, "pOverECutHigh_diElectron = %e\n", pOverECutHigh_diElectron);
  fprintf(f1, "invariantMassCutLow_diElectron = %e\n", invariantMassCutLow_diElectron);
  fprintf(f1, "invariantMassCutHigh_diElectron = %e\n", invariantMassCutHigh_diElectron);
  fprintf(f1, "vertexZDiffCut_diElectron = %e\n", vertexZDiffCut_diElectron);
  fprintf(f1, "oneOverBetaCut_diElectron = %e\n", oneOverBetaCut_diElectron);
  fprintf(f1, "onlyOppPairs_diElectron = %i\n", onlyOppPairs_diElectron);
  fprintf(f1, "\n");
  fprintf(f1, "sampleScale_randomEvents = %e\n", sampleScale_randomEvents);
  fprintf(f1, "\n");
//  fprintf(f1, "dedx innerSectorGain = %e\n", innerSectorGain);
//  fprintf(f1, "dedx outerSectorGain = %e\n", outerSectorGain);
  fprintf(f1, "\n");
  fprintf(f1, "spaceChargeP0 = %e\n", spaceChargeP0);
  fprintf(f1, "spaceChargeP1 = %e\n", spaceChargeP1);
  fprintf(f1, "\n");
  fprintf(f1, "xVertex = %e\n", paraXVertex);
  fprintf(f1, "yVertex = %e\n", paraYVertex);
  fprintf(f1, "\n");
  fprintf(f1, "bField = %e\n", bField);
  fprintf(f1, "scalerCount = %e\n", scalerCount);
  fprintf(f1, "\n");
  fprintf(f1, "sigmaDedx1 = %e\n", sigmaDedx1);
  fprintf(f1, "dPOverP2 = %e\n", dPOverP2);
  fprintf(f1, "sigmaTof = %e\n", sigmaTof);
  fprintf(f1, "\n");

}

void gl3TriggerDecider::readHLTparameters(char* HLTparameters)
{
  string parameterName;
  ifstream ifs(HLTparameters);
  if(!ifs.fail())
    while(!ifs.eof())
      {
	ifs>>parameterName;
	if(parameterName == "triggerOnHighPt") ifs>>triggerOnHighPt;
	if(parameterName == "triggerOnDiElectron") ifs>>triggerOnDiElectron;
	if(parameterName == "triggerOnHeavyFragment") ifs>>triggerOnHeavyFragment;
	if(parameterName == "triggerOnAllEvents") ifs>>triggerOnAllEvents;
	if(parameterName == "triggerOnRandomEvents") ifs>>triggerOnRandomEvents;
	if(parameterName == "triggerOnBesGoodEvents") ifs>>triggerOnBesGoodEvents;
	if(parameterName == "triggerBitHighPt") ifs>>hex>>triggerBitHighPt;
	if(parameterName == "triggerBitDiElectron") ifs>>hex>>triggerBitDiElectron;
	if(parameterName == "triggerBitHeavyFragment") ifs>>hex>>triggerBitHeavyFragment;
	if(parameterName == "triggerBitAllEvents") ifs>>hex>>triggerBitAllEvents;
	if(parameterName == "triggerBitRandomEvents") ifs>>hex>>triggerBitRandomEvents;
	if(parameterName == "triggerBitBesGoodEvents") ifs>>hex>>triggerBitBesGoodEvents;
	if(parameterName == "nHitsCut_highPt") ifs>>dec>>nHitsCut_highPt;
	if(parameterName == "ptCut_highPt") ifs>>ptCut_highPt;
	if(parameterName == "nHitsCut_heavyFragment") ifs>>nHitsCut_heavyFragment;
	if(parameterName == "nDedxCut_heavyFragment") ifs>>nDedxCut_heavyFragment;
	if(parameterName == "useTofMatchedGlobalTrack_heavyfragMent") ifs>>useTofMatchedGlobalTrack_heavyfragMent;
	if(parameterName == "nSigmaDedxHe3Cut_heavyFragment") ifs>>nSigmaDedxHe3Cut_heavyFragment;
	if(parameterName == "nSigmaDedxTritonCut_heavyFragment") ifs>>nSigmaDedxTritonCut_heavyFragment;
	if(parameterName == "nSigmaMassTritonCut_heavyFragment") ifs>>nSigmaMassTritonCut_heavyFragment;
	if(parameterName == "dcaCut_heavyFragment") ifs>>dcaCut_heavyFragment;
	if(parameterName == "nHitsCut_diElectron") ifs>>nHitsCut_diElectron;
	if(parameterName == "nDedxCut_diElectron") ifs>>nDedxCut_diElectron;
	if(parameterName == "towerEnergyCut_diElectron") ifs>>towerEnergyCut_diElectron;
	if(parameterName == "pCutLowForTof_diElectron") ifs>>pCutLowForTof_diElectron;
	if(parameterName == "pCutHighForTof_diElectron") ifs>>pCutHighForTof_diElectron;
	if(parameterName == "pCutLowForEmc_diElectron") ifs>>pCutLowForEmc_diElectron;
	if(parameterName == "pCutHighForEmc_diElectron") ifs>>pCutHighForEmc_diElectron;
	if(parameterName == "pCutLowForTofAndEmc_diElectron") ifs>>pCutLowForTofAndEmc_diElectron;
	if(parameterName == "pCutHighForTofAndEmc_diElectron") ifs>>pCutHighForTofAndEmc_diElectron;
	if(parameterName == "nSigmaDedxCutLowForTof_diElectron") ifs>>nSigmaDedxCutLowForTof_diElectron;
	if(parameterName == "nSigmaDedxCutHighForTof_diElectron") ifs>>nSigmaDedxCutHighForTof_diElectron;
	if(parameterName == "nSigmaDedxCutLowForEmc_diElectron") ifs>>nSigmaDedxCutLowForEmc_diElectron;
	if(parameterName == "nSigmaDedxCutHighForEmc_diElectron") ifs>>nSigmaDedxCutHighForEmc_diElectron;
	if(parameterName == "nSigmaDedxCutLowForTofAndEmc_diElectron") ifs>>nSigmaDedxCutLowForTofAndEmc_diElectron;
	if(parameterName == "nSigmaDedxCutHighForTofAndEmc_diElectron") ifs>>nSigmaDedxCutHighForTofAndEmc_diElectron;
	if(parameterName == "p1Cut_diElectron") ifs>>p1Cut_diElectron;
	if(parameterName == "p2Cut_diElectron") ifs>>p2Cut_diElectron;
	if(parameterName == "pOverECutLow_diElectron") ifs>>pOverECutLow_diElectron;
	if(parameterName == "pOverECutHigh_diElectron") ifs>>pOverECutHigh_diElectron;
	if(parameterName == "invariantMassCutLow_diElectron") ifs>>invariantMassCutLow_diElectron;
	if(parameterName == "invariantMassCutHigh_diElectron") ifs>>invariantMassCutHigh_diElectron;
	if(parameterName == "vertexZDiffCut_diElectron") ifs>>vertexZDiffCut_diElectron;
	if(parameterName == "oneOverBetaCut_diElectron") ifs>>oneOverBetaCut_diElectron;
	if(parameterName == "onlyOppPairs_diElectron") ifs>>onlyOppPairs_diElectron;

	if(parameterName == "sampleScale_randomEvents") ifs>>sampleScale_randomEvents;

	if(parameterName == "Debug_allTracks") ifs>>Debug_allTracks;
	if(parameterName == "Debug_tracks") ifs>>Debug_tracks;
	if(parameterName == "Debug_towers") ifs>>Debug_towers;
	if(parameterName == "Debug_matchs") ifs>>Debug_matchs;
	if(parameterName == "Debug_vertex") ifs>>Debug_vertex;
	if(parameterName == "nHitsCut_debug") ifs>>nHitsCut_debug;
	if(parameterName == "nDedxCut_debug") ifs>>nDedxCut_debug;
	if(parameterName == "dedxCutLow_debug") ifs>>dedxCutLow_debug;
	if(parameterName == "dedxCutHigh_debug") ifs>>dedxCutHigh_debug;
	if(parameterName == "pCut_debug") ifs>>pCut_debug;
	if(parameterName == "towerEnergyCut_debug") ifs>>towerEnergyCut_debug;
	if(parameterName == "matchPhiDiffCut_debug") ifs>>matchPhiDiffCut_debug;
	if(parameterName == "matchZEdgeCut_debug") ifs>>matchZEdgeCut_debug;
//	if(parameterName == "innerSectorGain") ifs>>innerSectorGain;
//	if(parameterName == "outerSectorGain") ifs>>outerSectorGain;
	if(parameterName == "spaceChargeP0") ifs>>spaceChargeP0;
	if(parameterName == "spaceChargeP1") ifs>>spaceChargeP1;
	if(parameterName == "xVertex") ifs>>paraXVertex;
	if(parameterName == "yVertex") ifs>>paraYVertex;
	if(parameterName == "sigmaDedx1") ifs>>sigmaDedx1;
	if(parameterName == "dPOverP2") ifs>>dPOverP2;
	if(parameterName == "sigmaTof") ifs>>sigmaTof;
      }
}

int gl3TriggerDecider::decide(int eventId)
{
 if (f1)  fprintf(f1, "event:  %i\n", eventId);

  int highPtTriggered = 0;
  int diElectronTriggered = 0;
  int heavyFragmentTriggered = 0;
  int allEventsTriggered = 0;
  int randomEventsTriggered = 0;
  int besGoodEventsTriggered = 0;

  hlt_diEP.nEPairs = 0;
  hlt_hiPt.nHighPt = 0;
  hlt_hF.nHeavyFragments = 0;

  int nTracks = event->getNGlobalTracks();

  //  LOG(TERR,"GL3: %d tracks",nTracks) ;

   if(triggerOnHighPt) {
     for(int i=0; i<nTracks; i++)
       {
	 gl3Node* node = event->getNode(i);
	 gl3Track* pTrack = event->getNode(i)->primaryTrack;
	 if(!pTrack) continue;

//	LOG(TERR,"hight Pt: %d/%d : %d vs %d, %f vs %f",i,nTracks,
//	    pTrack->nHits,nHitsCut_highPt,
//	    pTrack->pt, ptCut_highPt) ;

	 if(pTrack->nHits < nHitsCut_highPt) continue;
	 if(pTrack->pt < ptCut_highPt) continue;

	 highPtTriggered = triggerBitHighPt;

	 if(!node->primarytofCell && node->primaryTrack)
	   event->tof->matchPrimaryTracks(node,event->getLmVertex());
	 double matchTowerEnergy = 0.;
	 for(int j=0; j<event->emc->getNBarrelTowers(); j++)
	   {
	     gl3EmcTower* tower = event->emc->getBarrelTower(j);
	     double phiDiff, zEdge;
	     if(!tower->matchTrack(pTrack, phiDiff, zEdge)) continue;
	     if(tower->getEnergy() < matchTowerEnergy) continue;
	     matchTowerEnergy = tower->getEnergy();
	     node->emcTower = tower;
	     node->emcMatchPhiDiff = phiDiff;
	     node->emcMatchZEdge = zEdge;
	   }
   
	 hlt_hiPt.highPtNodeSN[hlt_hiPt.nHighPt] = i;
	 hlt_hiPt.nHighPt ++;
       }
   }

   if(triggerOnHeavyFragment) {
     for(int i=0; i<nTracks; i++)
       {
	 int triggered = 0;
	 gl3Node* node = event->getNode(i);
	 gl3Track* gTrack = event->getGlobalTrack(i);
	 gl3Track* pTrack = node->primaryTrack;
	 if(gTrack->nHits < nHitsCut_heavyFragment) continue;
	 if(gTrack->nDedx < nDedxCut_heavyFragment) continue;

	 if( useTofMatchedGlobalTrack_heavyfragMent && !node->globaltofCell) continue;  ///< require tof match for global tracks
	 
	 double x0 = gTrack->r0*cos(gTrack->phi0);
	 double y0 = gTrack->r0*sin(gTrack->phi0);
	 double dca2 = sqrt(pow(x0 - event->beamlineX, 2)+pow(y0 - event->beamlineY, 2));
	 
	 double dcaToUse = 0;
	 if(event->useBeamlineToMakePrimarys) dcaToUse = dca2 ; ///< using 2D dca in pp 500GeV
	 else dcaToUse = gTrack->dca ;  ///< using 3D dca in AuAu 200 GeV

	 if(dcaToUse > dcaCut_heavyFragment || dcaToUse < 0) continue;

	 if(gTrack->nSigmaDedx(event->bischel, "He3", sigmaDedx1) > nSigmaDedxHe3Cut_heavyFragment)
	   {
	     triggered = 1;
	     if(!node->primarytofCell && node->primaryTrack)
	       event->tof->matchPrimaryTracks(node,event->getLmVertex());
	     double matchTowerEnergy = 0.;
	     for(int i=0; i<event->emc->getNBarrelTowers(); i++)
	       {
		 gl3EmcTower* tower = event->emc->getBarrelTower(i);
		 double phiDiff, zEdge;
		 if(!tower->matchTrack(gTrack, phiDiff, zEdge)) continue;
		 if(tower->getEnergy() < matchTowerEnergy) continue;
		 matchTowerEnergy = tower->getEnergy();
		 node->emcTower = tower;
		 node->emcMatchPhiDiff = phiDiff;
		 node->emcMatchZEdge = zEdge;
	       }
	   }

	 if(fabs(gTrack->nSigmaDedx(event->bischel, "Triton", sigmaDedx1)) < nSigmaDedxTritonCut_heavyFragment)
	   {
	     if(!node->primarytofCell && node->primaryTrack)
	       event->tof->matchPrimaryTracks(node,event->getLmVertex());
	     if(node->primarytofCell)
	       {
		 //		 double c = 29.979; //cm/ns 
		 double p = pTrack->pt*sqrt(1.+pow(pTrack->tanl,2));
		 double beta = node->beta;
		 double tof = node->tof;
		 double l = beta*tof;
		 double tritonMass = 2.80925;
		 double sigmaMass1 = pow(p,2)*2./l*sqrt(1.+pow(tritonMass/p,2))*sigmaTof;
		 double sigmaMass2 = pow(tritonMass,2)*2.*p*dPOverP2;
		 double sigmaMass = sqrt(pow(sigmaMass1,2)+pow(sigmaMass2,2));
		 //LOG(INFO,"sigmaMass1:%f, sigmaMass2:%f, sigmaMass:%f,p:%f,tof:%f",sigmaMass1,sigmaMass2,sigmaMass,p,tof);
		 double m2 = pow(p,2)*(pow(1./beta,2)-1);
		 if(fabs(m2-pow(tritonMass,2))<nSigmaMassTritonCut_heavyFragment*sigmaMass)
		   {
		     triggered = 1;
		     double matchTowerEnergy = 0.;
		     for(int i=0; i<event->emc->getNBarrelTowers(); i++)
		       {
			 gl3EmcTower* tower = event->emc->getBarrelTower(i);
			 double phiDiff, zEdge;
			 if(!tower->matchTrack(gTrack, phiDiff, zEdge)) continue;
			 if(tower->getEnergy() < matchTowerEnergy) continue;
			 matchTowerEnergy = tower->getEnergy();
			 node->emcTower = tower;
			 node->emcMatchPhiDiff = phiDiff;
			 node->emcMatchZEdge = zEdge;
		       }
		   }
		 
		 //		 if(mass>tritonMass-nSigmaMassTritonCut_heavyFragment*sigmaMass)
		 //		   cout<<"triton trigger: "<<gTrack->nHits<<"  "<<gTrack->nDedx<<"  "<<gTrack->pt*sqrt(1.+pow(gTrack->tanl,2))<<"  "<<gTrack->dedx<<"  "<<gTrack->nSigmaDedx(event->bischel, "Triton", sigmaDedx1)<<"  "<<mass<<"+-"<<sigmaMass<<endl;
	       }
	   }
	 if(triggered)
	   {
	     heavyFragmentTriggered = triggerBitHeavyFragment;
	     hlt_hF.heavyFragmentSN[hlt_hF.nHeavyFragments] = i;
	     hlt_hF.nHeavyFragments ++;

	   }
       }
   }


   if(triggerOnDiElectron) { 
       const double electronMass = 0.511e-3;
       int nElectronNodes = 0, nTowers = 0;
       gl3Node* electronNode[maxNElectronTracks];
       gl3EmcTower* emcTower[maxNTowers];
       int electronSelectionBits[maxNElectronTracks];
       int selectionBitTof = 0x1;
       int selectionBitEmc = 0x2;
       int selectionBitTofAndEmc = 0x4;

       for(int i=0; i<event->emc->getNBarrelTowers(); i++)
	 {
	   gl3EmcTower* tower = event->emc->getBarrelTower(i);
	   double towerEnergy = tower->getEnergy();
	   if(towerEnergy < 0.5) continue;
	   if(nTowers >= maxNTowers)
	     {
	       cout<<"gl3TriggerDecider::decide()  WARN! towerArray Full!"<<endl;
	       break;
	     }
	   emcTower[nTowers] = tower;
	   nTowers ++;
	 }

       for(int i=0; i<nTracks; i++)
	 {
	   int isElectron = 0;
	   int passTofCut = 0;
	   int passEmcCut = 0;
	   gl3Node* node = event->getNode(i);
	   gl3Track* track = node->primaryTrack;
	   if(!track) continue;
	   if(track->nHits < nHitsCut_diElectron) continue;
	   if(track->nDedx < nDedxCut_diElectron) continue;
	   double p = track->pt*sqrt(1.+pow(track->tanl, 2));
	   if(p < p2Cut_diElectron) continue;
	   double nSigmaDedxE = track->nSigmaDedx(event->bischel, "Electron", sigmaDedx1);
	   if(nSigmaDedxE > min(nSigmaDedxCutLowForTof_diElectron, nSigmaDedxCutLowForTofAndEmc_diElectron) 
	      && nSigmaDedxE < max(nSigmaDedxCutHighForTof_diElectron, nSigmaDedxCutHighForTofAndEmc_diElectron)
	      && p > min(pCutLowForTof_diElectron, pCutLowForTofAndEmc_diElectron)
	      && p < max(pCutHighForTof_diElectron, pCutHighForTofAndEmc_diElectron))
	     {
	       if(!node->primarytofCell && node->primaryTrack)
		 event->tof->matchPrimaryTracks(node,event->getLmVertex());
	       if(fabs(1./node->beta - 1.) < oneOverBetaCut_diElectron) 
		   passTofCut = 1;
	     }

	   if(nSigmaDedxE > min(nSigmaDedxCutLowForEmc_diElectron, nSigmaDedxCutLowForTofAndEmc_diElectron)
	      && nSigmaDedxE < max(nSigmaDedxCutHighForEmc_diElectron, nSigmaDedxCutHighForTofAndEmc_diElectron)
	      && p > min(pCutLowForEmc_diElectron, pCutLowForTofAndEmc_diElectron)
	      && p < max(pCutHighForEmc_diElectron, pCutHighForTofAndEmc_diElectron))
	     {
	       double matchTowerEnergy = 0.;
	       for(int j=0; j<nTowers; j++)
		 {
		   gl3EmcTower* tower = emcTower[j];
		   double phiDiff, zEdge;
		   if(!tower->matchTrack(track, phiDiff, zEdge)) continue;
		   if(tower->getEnergy() < matchTowerEnergy) continue;
		   matchTowerEnergy = tower->getEnergy();
		   node->emcTower = tower;
		   node->emcMatchPhiDiff = phiDiff;
		   node->emcMatchZEdge = zEdge;
		 }
	       if(matchTowerEnergy > towerEnergyCut_diElectron &&
		  p/matchTowerEnergy > pOverECutLow_diElectron &&
		  p/matchTowerEnergy < pOverECutHigh_diElectron) 
		 passEmcCut = 1;
	     }	   
	   if(nSigmaDedxE > nSigmaDedxCutLowForTof_diElectron 
	      && nSigmaDedxE < nSigmaDedxCutHighForTof_diElectron 
	      && p > pCutLowForTof_diElectron 
	      && p < pCutHighForTof_diElectron && passTofCut)
	     isElectron |= selectionBitTof;
	   if(nSigmaDedxE > nSigmaDedxCutLowForEmc_diElectron 
	      && nSigmaDedxE < nSigmaDedxCutHighForEmc_diElectron 
	      && p > pCutLowForEmc_diElectron 
	      && p < pCutHighForEmc_diElectron && passEmcCut)
	     isElectron |= selectionBitEmc;
	   if(nSigmaDedxE > nSigmaDedxCutLowForTofAndEmc_diElectron 
	      && nSigmaDedxE < nSigmaDedxCutHighForTofAndEmc_diElectron 
	      && p > pCutLowForTofAndEmc_diElectron 
	      && p < pCutHighForTofAndEmc_diElectron 
	      && passTofCut && passEmcCut)
	     isElectron |= selectionBitTofAndEmc;
	   
	   if(isElectron)
	     {
	       if(nElectronNodes >= maxNElectronTracks)
		 {
		   cout<<"gl3TriggerDecider::decide()  WARN! electronTrackArray Full!"<<endl;
		   break;
		 }
	       electronNode[nElectronNodes] = node;	
	       electronSelectionBits[nElectronNodes] = isElectron;
	       nElectronNodes ++;
	     }
	 }

       for(int i=0; i<nElectronNodes; i++)
	 {
	   for(int j=i+1; j<nElectronNodes; j++)
	     {
	       int i1=i;
	       int i2=j;

	       gl3Node* node1 = electronNode[i];
	       gl3Track* track1 = node1->primaryTrack;	     
	       double p1 = track1->pt*sqrt(1.+pow(track1->tanl, 2));
	       gl3Node* node2 = electronNode[j];
	       gl3Track* track2 = node2->primaryTrack;
	       double p2 = track2->pt*sqrt(1.+pow(track2->tanl, 2));
	       if(p1<p2) 
		 {
		   swap(p1, p2);
		   swap(i1, i2);
		   swap(node1, node2);
		   swap(track1,track2);
		 }
	       if(p1 < p1Cut_diElectron) continue;
	       if(track1->q*track2->q > 0 && onlyOppPairs_diElectron) continue;
	       if(fabs(track1->z0 - track2->z0) > vertexZDiffCut_diElectron) continue;
	   
	       double pt1 = track1->pt;
	       double px1 = pt1*cos(track1->psi);
	       double py1 = pt1*sin(track1->psi);
	       double pz1 = pt1*track1->tanl;
	       double E1 = sqrt(electronMass*electronMass+p1*p1);
	       double pt2 = track2->pt;
	       double px2 = pt2*cos(track2->psi);
	       double py2 = pt2*sin(track2->psi);
	       double pz2 = pt2*track2->tanl;
	       double E2 = sqrt(electronMass*electronMass+p2*p2);
	       double px = px1+px2;
	       double py = py1+py2;
	       double pz = pz1+pz2;
	       double psi = atan2(py,px);
	       double pt = sqrt(pow(px,2)+pow(py,2));
	       double tanl = atan2(pz,pt);
	       double invariantMass = sqrt(pow(E1+E2,2)-pow(pt,2)-pow(pz,2));
	       //	       cout<<"px1:"<<px1<<"  py1:"<<py1<<"  pz1:"<<pz1<<"  E1:"<<E1<<"  px2:"<<px2<<"  py2:"<<py2<<"  pz2:"<<pz2<<"  E2:"<<E2<<"  invarMass:"<<invariantMass<<endl;
	       //	       cout<<"p1:"<<p1<<"  E1:"<<node1->emcTower<<" beta1:"<<node1->beta<<"  p2:"<<p2<<"  E2:"<<node2->emcTower<<" beta2:"<<node2->beta<<"  invarMass:"<<invariantMass<<endl;
	       if(invariantMass < invariantMassCutLow_diElectron) continue;
	       if(invariantMass > invariantMassCutHigh_diElectron) continue;

	       diElectronTriggered = triggerBitDiElectron;

	       hlt_diElectronPair diElePair;
	       diElePair.dau1NodeSN = (int)(node1-event->getNode(0));
	       diElePair.dau2NodeSN = (int)(node2-event->getNode(0));
	       diElePair.invariantMass = invariantMass;
	       diElePair.pt = pt;
	       diElePair.psi = psi;
	       diElePair.tanl = tanl;
	       diElePair.dau1SelectionBit = electronSelectionBits[i1];
	       diElePair.dau2SelectionBit = electronSelectionBits[i2];
	       hlt_diEP.ePair[hlt_diEP.nEPairs] = diElePair;
	       hlt_diEP.nEPairs ++;
	       	      
	       /*
	       l3EmcTowerInfo* towerInfo1 = 0;
	       if(matchTower[i1]) towerInfo1 = matchTower[i1]->getTowerInfo();		    
	       l3EmcTowerInfo* towerInfo2 = 0;
	       if(matchTower[i2]) towerInfo2 = matchTower[i2]->getTowerInfo();	       
	       if (f1){ //write out QA messages
	       fprintf(f1, "pair:  \n");
	       fprintf(f1, "%e  ", invariantMass);
	       if(matchTower[i1])
		 fprintf(f1, "%i  %i  %i  %i  %i  %i  %i  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %i  %i  %i  %i  %e  %e  %e  %e  %e  %e  %e  %e  %e  %i  %e  %e  %e  ", track1->id, track1->flag, track1->innerMostRow, track1->outerMostRow, track1->nHits, track1->nDedx, track1->q, track1->chi2[0], track1->chi2[1], track1->dedx, track1->pt, track1->phi0, track1->psi, track1->r0, track1->tanl, track1->z0, track1->length, track1->dpt, track1->dpsi, track1->dz0, track1->eta, track1->dtanl, p1, px1, py1, pz1, towerInfo1->getDaqID(), towerInfo1->getSoftID(), towerInfo1->getCrate(), towerInfo1->getCrateSeq(), towerInfo1->getPhi(), towerInfo1->getEta(), towerInfo1->getEtaMin(), towerInfo1->getEtaMax(), towerInfo1->getZ(), towerInfo1->getZmin(), towerInfo1->getZmax(), towerInfo1->getGain(), towerInfo1->getPedestal(), matchTower[i1]->getADC(), matchTower[i1]->getEnergy(), matchPhiDiff[i1], matchZEdge[i1]);
	       else
		 fprintf(f1, "%i  %i  %i  %i  %i  %i  %i  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %i  %i  %i  %i  %e  %e  %e  %e  %e  %e  %e  %e  %e  %i  %e  %e  %e  ", track1->id, track1->flag, track1->innerMostRow, track1->outerMostRow, track1->nHits, track1->nDedx, track1->q, track1->chi2[0], track1->chi2[1], track1->dedx, track1->pt, track1->phi0, track1->psi, track1->r0, track1->tanl, track1->z0, track1->length, track1->dpt, track1->dpsi, track1->dz0, track1->eta, track1->dtanl, p1, px1, py1, pz1, 0, 0, 0, 0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0, 0., 0., 0.);
	       if(matchTower[i2])
		 fprintf(f1, "%i  %i  %i  %i  %i  %i  %i  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %i  %i  %i  %i  %e  %e  %e  %e  %e  %e  %e  %e  %e  %i  %e  %e  %e\n", track2->id, track2->flag, track2->innerMostRow, track2->outerMostRow, track2->nHits, track2->nDedx, track2->q, track2->chi2[0], track2->chi2[1], track2->dedx, track2->pt, track2->phi0, track2->psi, track2->r0, track2->tanl, track2->z0, track2->length, track2->dpt, track2->dpsi, track2->dz0, track2->eta, track2->dtanl, p2, px2, py2, pz2, towerInfo2->getDaqID(), towerInfo2->getSoftID(), towerInfo2->getCrate(), towerInfo2->getCrateSeq(), towerInfo2->getPhi(), towerInfo2->getEta(), towerInfo2->getEtaMin(), towerInfo2->getEtaMax(), towerInfo2->getZ(), towerInfo2->getZmin(), towerInfo2->getZmax(), towerInfo2->getGain(), towerInfo2->getPedestal(), matchTower[i2]->getADC(), matchTower[i2]->getEnergy(), matchPhiDiff[i2], matchZEdge[i2]);
	       else
		 fprintf(f1, "%i  %i  %i  %i  %i  %i  %i  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %i  %i  %i  %i  %e  %e  %e  %e  %e  %e  %e  %e  %e  %i  %e  %e  %e\n", track2->id, track2->flag, track2->innerMostRow, track2->outerMostRow, track2->nHits, track2->nDedx, track2->q, track2->chi2[0], track2->chi2[1], track2->dedx, track2->pt, track2->phi0, track2->psi, track2->r0, track2->tanl, track2->z0, track2->length, track2->dpt, track2->dpsi, track2->dz0, track2->eta, track2->dtanl, p2, px2, py2, pz2, 0, 0, 0, 0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0, 0., 0., 0.);
               }
	       */
	     }
	 }
   }

   if(triggerOnAllEvents) allEventsTriggered = triggerBitAllEvents;

   if(triggerOnRandomEvents && eventId%sampleScale_randomEvents == 0) randomEventsTriggered = triggerBitRandomEvents;

   if(triggerOnBesGoodEvents) {
       if(fabs(event->getVertex().Getz())<30 && sqrt(pow(event->getVertex().Getx(),2)+pow(event->getVertex().Gety(),2))<2. && event->getNPrimaryTracks()>5)
	 besGoodEventsTriggered = triggerBitBesGoodEvents;
   }

   if(f1)
     fprintf(f1, "trigger: 0x%x\n", event->trigger | highPtTriggered | diElectronTriggered | heavyFragmentTriggered | allEventsTriggered | randomEventsTriggered | besGoodEventsTriggered);

   // don't write asc files if only minbias triggred
   if(highPtTriggered || diElectronTriggered || heavyFragmentTriggered) {
     writeQA();
   }
   
  return highPtTriggered | diElectronTriggered | heavyFragmentTriggered | allEventsTriggered | randomEventsTriggered | besGoodEventsTriggered;
}
	       
void gl3TriggerDecider::writeQA()
{
  if (!f1) return;

  int nElectronTracks = 0, nTowers = 0, nMatch = 0;
  gl3Track* electronTrack[maxNElectronTracks];
  gl3EmcTower* emcTower[maxNTowers];
  gl3Track* matchTrack[maxNElectronTracks];
  gl3EmcTower* matchTower[maxNTowers];
  double matchPhiDiff[maxNTowers];
  double matchZEdge[maxNTowers];
  int nTracks = event->getNGlobalTracks();
  if(Debug_vertex)
    {
      fprintf(f1, "vertex: %e %e %e\n", event->getVertex().Getx(), event->getVertex().Gety(), event->getVertex().Getz());
      fprintf(f1, "lmVertex: %e %e %e\n", event->getLmVertex().Getx(), event->getLmVertex().Gety(), event->getLmVertex().Getz());
    }
  if(Debug_tracks || Debug_matchs)
    for(int i=0; i<nTracks; i++)
      {
	gl3Track* track = event->getGlobalTrack(i);
	if(track->nHits < nHitsCut_debug) continue;
	if(track->nDedx < nDedxCut_debug) continue;
	//	   if(!track->flag) continue; // reject secondaries
	double p = track->pt*sqrt(1.+pow(track->tanl, 2));
	if(p < pCut_debug) continue;
	if(track->dedx < dedxCutLow_debug) continue;
	if(track->dedx > dedxCutHigh_debug) continue;
	if(nElectronTracks >= maxNElectronTracks)
	  {
	    cout<<"gl3TriggerDecider::debug()  WARN! electronTrackArray Full!"<<endl;
	    break;
	  }
	double xVertex = track->getPara()->xVertex;
	double yVertex = track->getPara()->yVertex;
	track->updateToClosestApproach(xVertex, yVertex);
	electronTrack[nElectronTracks] = track;
	nElectronTracks ++;
      }
  if(Debug_towers || Debug_matchs)
    for(int i=0; i<event->emc->getNBarrelTowers(); i++)
      {
	gl3EmcTower* tower = event->emc->getBarrelTower(i);
	double towerEnergy = tower->getEnergy();
	if(towerEnergy < towerEnergyCut_debug) continue;
	if(nTowers >= maxNTowers)
	  {
	    cout<<"gl3TriggerDecider::debug()  WARN! towerArray Full!"<<endl;
	    break;
	  }
	emcTower[nTowers] = tower;
	nTowers ++;
    }
  if(Debug_matchs)
    for(int i=0; i<nElectronTracks; i++)
      for(int j=0; j<nTowers; j++)
	{
	  gl3EmcTower* tower = emcTower[j];
	  gl3Track* track = electronTrack[i];	     
	  double phiDiff, zEdge;
	  if(!tower->matchTrack(track, phiDiff, zEdge, matchPhiDiffCut_debug, matchZEdgeCut_debug)) continue;
	  matchTrack[nMatch] = electronTrack[i];
	  matchTower[nMatch] = emcTower[j];
	  matchPhiDiff[nMatch] = phiDiff;
	  matchZEdge[nMatch] = zEdge;
	  nMatch ++;
	}

  if(Debug_tracks)
    {
      fprintf(f1, "tracks:  %i\n", nElectronTracks);      
      for(int i=0; i<nElectronTracks; i++)
	{
	  gl3Track* track = electronTrack[i];	     
	  double pt = track->pt;
	  double px = pt*cos(track->psi);
	  double py = pt*sin(track->psi);
	  double pz = pt*track->tanl;
	  double p = sqrt(pow(pt,2) + pow(pz,2));
	  fprintf(f1, "%i  %i  %i  %i  %i  %i  %i  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n", track->id, track->flag, track->innerMostRow, track->outerMostRow, track->nHits, track->nDedx, track->q, track->chi2[0], track->chi2[1], track->dedx, track->pt, track->phi0, track->psi, track->r0, track->tanl, track->z0, track->length, track->dpt, track->dpsi, track->dz0, track->eta, track->dtanl, p, px, py, pz);
	}
    }
  if(Debug_towers)
    {
      fprintf(f1, "towers:    %i\n", nTowers);
      for(int i=0; i<nTowers; i++)
	{
	  gl3EmcTower* tower = emcTower[i];
	  l3EmcTowerInfo* towerInfo = tower->getTowerInfo();
	  fprintf(f1, "%i  %i  %i  %i  %e  %e  %e  %e  %e  %e  %e  %e  %e  %i  %e\n", towerInfo->getDaqID(), towerInfo->getSoftID(), towerInfo->getCrate(), towerInfo->getCrateSeq(), towerInfo->getPhi(), towerInfo->getEta(), towerInfo->getEtaMin(), towerInfo->getEtaMax(), towerInfo->getZ(), towerInfo->getZmin(), towerInfo->getZmax(), towerInfo->getGain(), towerInfo->getPedestal(), tower->getADC(), tower->getEnergy());
	}
    }
  if(Debug_matchs)
    {
      fprintf(f1, "matchs:     %i\n", nMatch);
      for(int i=0; i<nMatch; i++)
	{
	  gl3Track* track = matchTrack[i];
	  gl3EmcTower* tower = matchTower[i];
	  l3EmcTowerInfo* towerInfo = tower->getTowerInfo();
	  double pt = track->pt;
	  double px = pt*cos(track->psi);
	  double py = pt*sin(track->psi);
	  double pz = pt*track->tanl;
	  double p = sqrt(pow(pt,2) + pow(pz,2));
	  fprintf(f1, "%i  %i  %i  %i  %i  %i  %i  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %i  %i  %i  %i  %e  %e  %e  %e  %e  %e  %e  %e  %e  %i  %e  %e  %e\n", track->id, track->flag, track->innerMostRow, track->outerMostRow, track->nHits, track->nDedx, track->q, track->chi2[0], track->chi2[1], track->dedx, track->pt, track->phi0, track->psi, track->r0, track->tanl, track->z0, track->length, track->dpt, track->dpsi, track->dz0, track->eta, track->dtanl, p, px, py, pz, towerInfo->getDaqID(), towerInfo->getSoftID(), towerInfo->getCrate(), towerInfo->getCrateSeq(), towerInfo->getPhi(), towerInfo->getEta(), towerInfo->getEtaMin(), towerInfo->getEtaMax(), towerInfo->getZ(), towerInfo->getZmin(), towerInfo->getZmax(), towerInfo->getGain(), towerInfo->getPedestal(), tower->getADC(), tower->getEnergy(), matchPhiDiff[i], matchZEdge[i]);
	}
    }
}
       
void gl3TriggerDecider::writeScalers()
{
  if(f1)
    {
      fprintf(f1, "bField = %e\n", bField);
      fprintf(f1, "scalerCount = %e\n", scalerCount);
    }
}
