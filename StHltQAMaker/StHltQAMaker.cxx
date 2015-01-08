/****************************************************************************
 **
 **    Author: Liang Xue
 **    Description: HLT QA Maker to do the HLT information propaganda from 
 **    online RTS_Reader to offline StEvent
 **
 *****************************************************************************/

#include "StHltQAMaker.h"
#include "StEvent/StEvent.h"
#include "StEvent/StHltEvent.h"
#include "StEvent/StHltTrack.h"
#include "StEvent/StHltTrackNode.h"
#include "StEvent/StHltBTofHit.h"
#include "StEvent/StHltVpdHit.h"
#include "StEvent/StHltTriggerReason.h"
#include "StEvent/StHltHeavyFragment.h"
#include "StEvent/StHltHighPt.h"
#include "StEvent/StHltDiElectron.h"
#include "StEvent/StHltBEmcTowerHit.h"

#include "StThreeVectorF.hh"
#include "TVector3.h"

ClassImp(StHltQAMaker)

//________________________________________________________________________________________
StHltQAMaker::StHltQAMaker():StMaker(),mStEvent(0)
{
	mOutFile = 0;

	hvertexX = 0;
	hvertexY = 0;
	hvertexZ = 0;
	hlmvertexX = 0;
	hlmvertexY = 0;
	hlmvertexZ = 0;
	hPvpdVzvslmVz = 0;

	hnHits = 0;
	hndedx = 0;
	hgPt = 0;
	hgEta = 0;
	hgPhi = 0;
	hgChi2Xy = 0;
	hgChi2Z = 0;
	hgMult = 0;
	hDcaXy = 0;
	hDcaZ = 0;
	hgdEdxvsP = 0;

	hpPt = 0;
	hpEta = 0;
	hpPhi = 0;
	hpChi2Xy = 0;
	hpChi2Z = 0;
	hpMult = 0;
	hpdEdxvsP = 0;

	hEmcPhiDiff = 0;
	hEmczEdge = 0;
	hEmcEnergyfromNode = 0;
	hEmcEtaPhifromNode = 0;
	hEmcSoftIdfromNode = 0;
	hEmcdaqIdfromNode = 0;
	hEmcEnergy = 0;
	hEmcEtaPhi = 0;
	hEmcSoftId = 0;
	hEmcdaqId = 0;


    hLocalZ = 0;
	hLocalY = 0;
	hInBetaVsPrimPt = 0;
	hInBetaVsGlobPt = 0;
	hMatchIdvsFireId = 0;
	hTrayIDvsTrgTime = 0;
	hChannelID = 0;

	hHFdEdxvsP = 0;
	hHighptP = 0;
	hHFdEdxvsPFromHltEvent = 0;
	hHighptPFromHltEvent = 0;

	hDau1dEdxvsP = 0;
	hDau1Energy = 0;
	hDiElecInvMass = 0;
	hDau1Beta = 0;
	hDau2dEdxvsP = 0;
	hDau2Energy = 0;
	hDau2Beta = 0;
	hDiElecInvMassBG = 0;

	hDau1dEdxvsPFromHltEvent = 0;
	hDau1EnergyFromHltEvent = 0;
	hDiElecInvMassFromHltEvent = 0;
	hDau1BetaFromHltEvent = 0;
	hDau2dEdxvsPFromHltEvent = 0;
	hDau2EnergyFromHltEvent = 0;
	hDau2BetaFromHltEvent = 0;
	hDiElecInvMassBGFromHltEvent = 0;
}

//________________________________________________________________________________________
StHltQAMaker::~StHltQAMaker()
{
}

//________________________________________________________________________________________
Int_t StHltQAMaker::Init()
{
    
	char currentdir[256];
	char tmp[256];
	getcwd(tmp,256);
	strcpy(currentdir,tmp);
	printf("Current Dir : %s \n",currentdir);

	char rootfile[256];
	sprintf(rootfile,"%s/hltQA.root",currentdir);

	mOutFile = new TFile(rootfile,"RECREATE") ;

	hvertexX = new TH1F("hvertexX","vertexX",100,-10,10);
	hvertexY = new TH1F("hvertexY","vertexY",100,-10,10);
	hvertexZ = new TH1F("hvertexZ","vertexZ",400,-200,200);
	hlmvertexX = new TH1F("hlmvertexX","lmvertexX",100,-10,10);
	hlmvertexY = new TH1F("hlmvertexY","lmvertexY",100,-10,10);
	hlmvertexZ = new TH1F("hlmvertexZ","lmvertexZ",400,-200,200);
	hPvpdVzvslmVz = new TH2F("hPvpdVzvslmVz","hPvpdVzvslmVz",400,-100,100,400,-100,100);

	hnHits = new TH1I("hnHits","nHits",50,0,50);
	hndedx = new TH1I("hndedx","ndedx",50,0,50);
	hgPt  = new TH1F("hgPt","global tracks Pt distribution",150,0,15.);
	hgEta = new TH1F("hgEta","global tracks Eta distribution",120,-3,3);
	hgPhi = new TH1F("hgPhi","global tracks Phi distribution",300,0,2*PI);
	hgChi2Xy = new TH1F("hgChi2Xy","global tracks Chi2Xy distribution",100,0,10);
	hgChi2Z = new TH1F("hgChi2Z","global tracks Chi2Z distribution",100,0,10);
	hgMult = new TH1I("hgMult","global tracks Mult distribution",4000,0,4000);
	hDcaXy = new TH1F("hDcaXy","DcaXy dsitribution",300,-6,6);
	hDcaZ = new TH1F("hDcaZ","DcaZ dsitribution",300,-6,6);
	hgdEdxvsP = new TH2F("hgdEdxvsP","global tracks dEdxvsP distribution",500,-5,5,100,0,10.e-06);

	hpPt  = new TH1F("hpPt","primary tracks Pt distribution",150,0,15.);
	hpEta = new TH1F("hpEta","primary tracks Eta distribution",120,-3,3);
	hpPhi = new TH1F("hpPhi","primary tracks Phi distribution",300,0,2*PI);
	hpChi2Xy = new TH1F("hpChi2Xy","primary tracks Chi2Xy distribution",100,0,10);
	hpChi2Z = new TH1F("hpChi2Z","primary tracks Chi2Z distribution",100,0,10);
	hpMult = new TH1I("hpMult","primary tracks Mult distribution",4000,0,4000);
	hpdEdxvsP = new TH2F("hpdEdxvsP","prmary tracks dEdxvsP distribution",500,-5,5,100,0,100.e-06);

    hEmcPhiDiff = new TH1D("hEmcPhiDiff","Emc PhiDiff distribution",100,0,0.1);
	hEmczEdge   = new TH1D("hEmczEdge","Emc zEdge distribution",100,0,5.);
	hEmcEnergyfromNode = new TH1F("hEmcEnergyfromNode","Emc energy obtain from Node",200,0,20);
	hEmcEtaPhifromNode = new TH2F("hEmcEtaPhifromNode","Emc eta phi obtain from Node",120,0,2*PI,100,-1,1);
	hEmcSoftIdfromNode = new TH1I("hEmcSoftIdfromNode","Emc softId obtain from Node",5000,0,5000);
	hEmcdaqIdfromNode = new TH1I("hEmcdaqIdfromNode","Emc daqID obtain from Node",5000,0,5000);
	hEmcEnergy = new TH1F("hEmcEnergy","Emc energy obtain from collection",200,0,20);
	hEmcEtaPhi = new TH2F("hEmcEtaPhi","Emc eta phi obtain from collection",120,0,2*PI,100,-1,1);
	hEmcSoftId = new TH1I("hEmcSoftId","Emc softId obtain from collection",5000,0,5000);
	hEmcdaqId = new TH1I("hEmcdaqId","Emc daqID obtain from collection",5000,0,5000);

    hLocalZ = new TH1F("hLocalZ","tof LocalZ",200,-10,10);
	hLocalY = new TH1F("hLocalY","tof LocalY",200,-10,10);
	hInBetaVsPrimPt = new TH2F("hInBetaVsPrimPt","Inverse Beta vs Primary Pt",500,0,5,500,0,5);
    hInBetaVsGlobPt = new TH2F("hInBetaVsGlobPt","Inverse Beta vs Global Pt",500,0,5,500,0,5);
	hMatchIdvsFireId = new TH2F("hMatchIdvsFireId","tof matched ID vs fire ID",200,0,200,200,0,200);
	hTrayIDvsTrgTime = new TH2F("hTrayIDvsTrgTime","tof trayId vs trigger time",124,0,124,400,2700,3100);
	hChannelID = new TH1F("hChannelID","tof channel Id",200,0,200);

	hHFdEdxvsP = new TH2F("hHFdEdxvsP","prmary tracks dEdxvsP for Heavy Fragment distribution",500,-5,5,100,0,100.e-06);
	hHighptP  = new TH1F("hHighptP","high pt tracks Pt distribution",150,0,15.);
	hHFdEdxvsPFromHltEvent = new TH2F("hHFdEdxvsPFromHltEvent","prmary tracks dEdxvsP for Heavy Fragment distribution From Hltevent",500,-5,5,100,0,100.e-06);
	hHighptPFromHltEvent  = new TH1F("hHighptPFromHltEvent","high pt tracks Pt distribution From HltEvent",150,0,15.);
	
	hDau1dEdxvsP = new TH2F("hDau1dEdxvsP"," daughter1 dEdxvsP distribution",500,0,5,100,0,10.e-06);
	hDau2dEdxvsP = new TH2F("hDau2dEdxvsP"," daughter2 dEdxvsP distribution",500,0,5,100,0,10.e-06);
	hDau1Energy = new TH1F("hDau1Energy","daughter1 Energy",100,0,5);
	hDau2Energy = new TH1F("hDau2Energy","daughter2 Energy",100,0,5);
	hDau1Beta = new TH1F("hDau1Beta","hDau1Beta",100,0,5);
	hDau2Beta = new TH1F("hDau2Beta","hDau2Beta",100,0,5);
	hDiElecInvMass = new TH1F("hDiElecInvMass","hDiElecInvMass",120,1.,13.);
	hDiElecInvMassBG = new TH1F("hDiElecInvMassBG","hDiElecInvMassBG",120,1.,13.);

	hDau1dEdxvsPFromHltEvent = new TH2F("hDau1dEdxvsPFromHltEvent"," daughter1 dEdxvsP distribution From HltEvent",500,0,5,100,0,10.e-06);
	hDau2dEdxvsPFromHltEvent = new TH2F("hDau2dEdxvsPFromHltEvent"," daughter2 dEdxvsP distribution From HltEvent",500,0,5,100,0,10.e-06);
	hDau1EnergyFromHltEvent = new TH1F("hDau1EnergyFromHltEvent","daughter1 Energy From HltEvent",100,0,5);
	hDau2EnergyFromHltEvent = new TH1F("hDau2EnergyFromHltEvent","daughter2 Energy From HltEvent",100,0,5);
	hDau1BetaFromHltEvent = new TH1F("hDau1BetaFromHltEvent","hDau1BetaFromHltEvent",100,0,5);
	hDau2BetaFromHltEvent = new TH1F("hDau2BetaFromHltEvent","hDau2BetaFromHltEvent",100,0,5);
	hDiElecInvMassFromHltEvent = new TH1F("hDiElecInvMassFromHltEvent","hDiElecInvMassFromHltEvent",120,1.,13.);
	hDiElecInvMassBGFromHltEvent = new TH1F("hDiElecInvMassBGFromHltEvent","hDiElecInvMassBGFromHltEvent",120,1.,13.);

	return StMaker::Init();
}

//________________________________________________________________________________________
Int_t StHltQAMaker::InitRun(int runumber)
{
return StMaker::Init();
}

//________________________________________________________________________________________
void  StHltQAMaker::Clear(Option_t *)
{
	  StMaker::Clear(); // perform the basic clear (mandatory)
}

//________________________________________________________________________________________
Int_t StHltQAMaker::Make()
{

	StHltEvent *HltEvent = 0;
    
	mStEvent = dynamic_cast<StEvent*> (GetInputDS("StEvent"));

	if(mStEvent) {

		HltEvent = mStEvent->hltEvent();
		if(HltEvent){
			FillEvent(HltEvent);
			FillGlobalTrack(HltEvent);
			FillPrimaryTrack(HltEvent);
			FillNode(HltEvent);
			FillEmc(HltEvent);
			FillTof(HltEvent);
			FillVpd(HltEvent);
			FillTriggerReason(HltEvent);
			FillTriggerParticleFromHltEvent(HltEvent);
		}
		else {
			LOG_WARN << "No StHltEvent found" << endm;
		}

	}
	else {LOG_WARN << "No StEvent found" << endm; }
	
	return kStOK;

}

//________________________________________________________________________________________
void StHltQAMaker::FillEvent(StHltEvent* event)
{
  
  StThreeVectorF vertexPos   =  event->vertex();
  hvertexX->Fill(vertexPos.x());
  hvertexY->Fill(vertexPos.y());
  hvertexZ->Fill(vertexPos.z());

  StThreeVectorF lmvertexPos =  event->lowMultVertex();
  hlmvertexX->Fill(lmvertexPos.x());
  hlmvertexY->Fill(lmvertexPos.y());
  hlmvertexZ->Fill(lmvertexPos.z());

  hPvpdVzvslmVz->Fill(event->vpdVertexZ(),lmvertexPos.z());

  float t0 = event->t0();
  float innerGain = event->innerSecGain();
  float outerGain = event->outerSecGain();
  unsigned int version = event->version();
  unsigned int triggerreason = event->triggerReasonBitOred();

  printf("t0 = %f, \t innerGain = %e, \t outerGain = %e \n",t0,innerGain,outerGain);
  printf("triggerreason = 0x%X, \t version = 0x%X \n ",triggerreason,version);

}

//________________________________________________________________________________________
void StHltQAMaker::FillGlobalTrack(StHltEvent* event)
{

	StSPtrVecHltTrack& VecgTrack = event->globalTrack();

	hgMult->Fill(VecgTrack.size());

	for(unsigned int i=0; i< VecgTrack.size(); i++){

		StHltTrack *gTrack = (StHltTrack*)VecgTrack.at(i);
		if(gTrack){
			float pt  =  gTrack->pt();
			float px  =  cos(gTrack->psi())*gTrack->pt();
			float py  =  sin(gTrack->psi())*gTrack->pt();
			float pz  =  gTrack->tanl()*gTrack->pt();

			TVector3 mom(px,py,pz);
			float eta = mom.PseudoRapidity();
			float phi = mom.Phi();
			if(phi<0.) phi = phi + 2*PI;
			float p = mom.Mag();

			float dedx  =  gTrack->dedx();
			float nhits =  gTrack->nHits();
			float ndedx =  gTrack->ndedx();
			float q     =  gTrack->q();
			hnHits->Fill(nhits);
			hndedx->Fill(ndedx);

			float chi2xy = gTrack->chi2(0);
			float chi2z = gTrack->chi2(1);
			hgChi2Xy->Fill(chi2xy);
			hgChi2Z->Fill(chi2z);

			if(nhits >= 25 && fabs(eta) < 1.){
				hgPt->Fill(pt);
				hgPhi->Fill(phi);
				hgEta->Fill(eta);
			}
			if(nhits >= 20 && ndedx >= 15) hgdEdxvsP->Fill(p*q,dedx);

			float dcaX = gTrack->r0()*cos(gTrack->phi0()) - event->lowMultVertex().x();
			float dcaY = gTrack->r0()*sin(gTrack->phi0()) - event->lowMultVertex().y();
			float cross = dcaX*sin(gTrack->psi()) - dcaY*cos(gTrack->psi());
			float theSign = (cross>=0) ? 1. : -1.;
			float dcaXy = theSign*sqrt(pow(dcaX,2)+pow(dcaY,2));
			float dcaZ = gTrack->z0() - event->lowMultVertex().z();

			hDcaXy->Fill(dcaXy);
			hDcaZ->Fill(dcaZ);

			//if(gTrack->trackNode()) printf(" i = %d, globalSN = %d \n",i, gTrack->trackNode()->globalTrackSN());
		}

	}

}

//________________________________________________________________________________________
void StHltQAMaker::FillPrimaryTrack(StHltEvent* event)
{

	StSPtrVecHltTrack& VecpTrack = event->primaryTrack();

	hpMult->Fill(VecpTrack.size());

	for(unsigned int i=0; i< VecpTrack.size(); i++){

		StHltTrack *pTrack = (StHltTrack*)VecpTrack.at(i);
		if(pTrack){
			float pt  =  pTrack->pt();
			float px  =  cos(pTrack->psi())*pTrack->pt();
			float py  =  sin(pTrack->psi())*pTrack->pt();
			float pz  =  pTrack->tanl()*pTrack->pt();

			TVector3 mom(px,py,pz);
			float eta = mom.PseudoRapidity();
			float phi = mom.Phi();
			float p = mom.Mag();

			float dedx  =  pTrack->dedx();
			float nhits =  pTrack->nHits();
			float ndedx =  pTrack->ndedx();
			float q     =  pTrack->q();

			float chi2xy = pTrack->chi2(0);
			float chi2z = pTrack->chi2(1);
			hpChi2Xy->Fill(chi2xy);
			hpChi2Z->Fill(chi2z);

			if(nhits >= 25 && fabs(eta) < 1.){
				hpPt->Fill(pt);
				hpPhi->Fill(phi);
				hpEta->Fill(eta);
			}
			if(nhits >= 20 && ndedx >= 15) hpdEdxvsP->Fill(p*q,dedx);

			if(pTrack->trackNode()) printf(" i = %d, primaryTrackSN = %d \n",i, pTrack->trackNode()->primaryTrackSN());
		}

	}

}

void StHltQAMaker::FillNode(StHltEvent* event)
{

	StSPtrVecHltTrackNode& VecTrackNode = event->trackNode();
	for(unsigned int i=0; i< VecTrackNode.size(); i++){

		StHltTrackNode *trackNode = (StHltTrackNode*)VecTrackNode.at(i);
		if(trackNode){
			StHltTrack*        gTrack   =  (StHltTrack*)   trackNode->globalTrack();
			StHltTrack*        pTrack   =  (StHltTrack*)   trackNode->primaryTrack();
			StHltBTofHit*      bTofHit  =  (StHltBTofHit*) trackNode->bTofHit();
			StHltBEmcTowerHit* bEmcHit  =  (StHltBEmcTowerHit*) trackNode->bEmcTowerHit();

			double emcPhiDiff = trackNode->bEmcMatchPhiDiff();
			double emczEdge   = trackNode->bEmcMatchZEdge();
			hEmcPhiDiff->Fill(emcPhiDiff);
			if( emczEdge>0. ) hEmczEdge->Fill(emczEdge);

			float tofLocalY = trackNode->bTofCellLocalY();
			float tofLocalZ = trackNode->bTofCellLocalZ();
			float tofBeta = trackNode->beta();
			float tof     = trackNode->tof();
			float tofProjChannel = trackNode->bTofProjChannel();
			float tofPathLength = trackNode->bTofPathLength();
			if(bTofHit){
				hLocalZ->Fill(tofLocalY);
				hLocalY->Fill(tofLocalZ);
				hMatchIdvsFireId->Fill(tofProjChannel,bTofHit->channel());
				if(pTrack) hInBetaVsPrimPt->Fill(pTrack->pt(),1./tofBeta);
				if(gTrack) hInBetaVsGlobPt->Fill(gTrack->pt(),1./tofBeta);
			}
			if(bEmcHit) {
				float energy  =  bEmcHit->energy(); 
				float phi     =  bEmcHit->phi();
				float eta     =  bEmcHit->eta();
				int softId   =  bEmcHit->softId();
				int daqId    =  bEmcHit->daqId();
				hEmcEnergyfromNode->Fill(energy);
				hEmcEtaPhifromNode->Fill(phi,eta);
				hEmcSoftIdfromNode->Fill(softId);
				hEmcdaqIdfromNode->Fill(daqId);
			}
		}

	}

}


//________________________________________________________________________________________
void StHltQAMaker::FillEmc(StHltEvent* event)
{

	StSPtrVecHltBEmcTowerHit& VecbEmcHit = event->bEmcTowerHits();
	for(unsigned int i=0; i< VecbEmcHit.size(); i++){

		StHltBEmcTowerHit *bEmcHit = (StHltBEmcTowerHit*)VecbEmcHit.at(i);
		if(bEmcHit){
			float energy = bEmcHit->energy();
			float phi = bEmcHit->phi();
			float eta = bEmcHit->eta();
			int softID = bEmcHit->softId();
			int daqID = bEmcHit->daqId();
			hEmcEnergy->Fill(energy);
			hEmcEtaPhi->Fill(phi,eta);
			hEmcSoftId->Fill(softID);
			hEmcdaqId->Fill(daqID);

			if(bEmcHit->trackNode()) printf(" i = %d, BEmcSN = %d \n",i, bEmcHit->trackNode()->emcTowerSN());
		}

	}

}

//________________________________________________________________________________________
void StHltQAMaker::FillTof(StHltEvent* event)
{

	StSPtrVecHltBTofHit& VecbTofHit = event->bTofHit();
	for(unsigned int i=0; i< VecbTofHit.size(); i++){

		StHltBTofHit *bTofHit = (StHltBTofHit*)VecbTofHit.at(i);
		if(bTofHit){
			short trayId = bTofHit->trayId();
			short channelId = bTofHit->channel();
			float tdc = bTofHit->tdc();
			float triggerTime = bTofHit->triggerTime();

			hTrayIDvsTrgTime->Fill(trayId,tdc-triggerTime);
			hChannelID->Fill(channelId);

			if(bTofHit->trackNode()) printf(" i = %d, BTofSN = %d \n",i, bTofHit->trackNode()->tofHitSN());
		}

	}

}

//________________________________________________________________________________________
void StHltQAMaker::FillVpd(StHltEvent* event)
{

	StSPtrVecHltVpdHit& VecpVpdHit = event->vpdHit();
	for(unsigned int i=0; i< VecpVpdHit.size(); i++){

		StHltVpdHit *pVpdHit = (StHltVpdHit*)VecpVpdHit.at(i);
		if(pVpdHit){
			StBeamDirection direction = pVpdHit->direction();
			short channelId = pVpdHit->channel();
			float tdc = pVpdHit->tdc();
			float triggerTime = pVpdHit->triggerTime();

			//	hTrayIDvsTrgTime->Fill(trayId,tdc-triggerTime);
		}

	}

}

//________________________________________________________________________________________
void StHltQAMaker::FillTriggerParticleFromHltEvent(StHltEvent* event)
{

	StSPtrVecHltDiElectron& VecDiElectron = event->diElectron();
	for(unsigned int i=0; i<VecDiElectron.size(); i++){
		StHltDiElectron *dielectron = (StHltDiElectron*)VecDiElectron.at(i);
		if(dielectron) FillDiep(dielectron);
	}
	
	StSPtrVecHltHighPt& VecHighPt = event->highPt();
	for(unsigned int i=0; i<VecHighPt.size(); i++){
		StHltHighPt *highpt = (StHltHighPt*)VecHighPt.at(i);
		if(highpt) FillHighPt(highpt);
	}

	StSPtrVecHltHeavyFragment& VecHeavyFragment = event->heavyFragment();
	for(unsigned int i=0; i<VecHeavyFragment.size(); i++){
		StHltHeavyFragment *heavyfragment = (StHltHeavyFragment*)VecHeavyFragment.at(i);
		if(heavyfragment) FillHF(heavyfragment);
	}

}

//________________________________________________________________________________________
void StHltQAMaker::FillHighPt(StHltHighPt* highpt)
{

	StHltTrack pTrack = highpt->primaryTrack();
	if(pTrack.nHits()>25) hHighptPFromHltEvent->Fill(pTrack.pt()); 

}

//________________________________________________________________________________________
void StHltQAMaker::FillHF(StHltHeavyFragment* heavyfragment)
{

//	StHltTrack        gTrack   =  (StHltTrack&)   heavyfragment->globalTrack();
//	StHltTrack        pTrack   =  (StHltTrack&)   heavyfragment->primaryTrack();
//	StHltBTofHit      bTofHit  =  (StHltBTofHit&) heavyfragment->bTofHit();
//	StHltBEmcTowerHit bEmcHit  =  (StHltBEmcTowerHit&) heavyfragment->bEmcTowerHit();

	StHltTrack        pTrack   =  heavyfragment->primaryTrack();
	int nHits = pTrack.nHits();
	int ndedx = pTrack.ndedx();
	int q = pTrack.q();
	float dedx = pTrack.dedx();

	float px  =  cos(pTrack.psi())*pTrack.pt();
	float py  =  sin(pTrack.psi())*pTrack.pt();
	float pz  =  pTrack.tanl()*pTrack.pt();
	TVector3 mom(px,py,pz);
	float p = mom.Mag();
	if(nHits >= 20 && ndedx >= 15) hHFdEdxvsPFromHltEvent->Fill(p*q,dedx);

}

//________________________________________________________________________________________
void StHltQAMaker::FillDiep(StHltDiElectron* diep)
{
    
//	if(diep->daughter1primaryTrackSN()<0 || diep->daughter2primaryTrackSN()<0) return 0;

	int  dau1q       =  diep->daughter1primaryTrack().q();
	int  dau1nHits   =  diep->daughter1primaryTrack().nHits();
	int  dau1ndedx   =  diep->daughter1primaryTrack().ndedx();
	float dau1dedx    =  diep->daughter1primaryTrack().dedx();
	float dau1pt      =  diep->daughter1primaryTrack().pt();
	float dau1px      =  diep->daughter1primaryTrack().pt()*cos(diep->daughter1primaryTrack().psi());
	float dau1py      =  diep->daughter1primaryTrack().pt()*sin(diep->daughter1primaryTrack().psi());
	float dau1pz      =  diep->daughter1primaryTrack().pt()*diep->daughter1primaryTrack().tanl();
	float dau1phidiff =  diep->daughter1bEmcMatchPhiDiff();

	TVector3 dau1(dau1px,dau1py,dau1pz);
	float dau1eta = dau1.PseudoRapidity();
	float dau1phi = dau1.Phi();
	if(dau1phi < 0.) dau1phi += 2*PI;
	float dau1p = dau1.Mag();

	float dau1beta = diep->daughter1beta();
	float dau1Energy = -999.;
	float dau1PEratio = -999.;
	float dau1EPratio = -999.;
	if(diep->daughter1emcTowerSN()>=0){
		dau1Energy = diep->daughter1bEmcTowerHit().energy();
		dau1PEratio = dau1p/dau1Energy; 
		dau1EPratio = dau1Energy/dau1p;
	}

	hDau1BetaFromHltEvent->Fill(1/dau1beta);
	hDau1EnergyFromHltEvent->Fill(dau1Energy);
	hDau1dEdxvsPFromHltEvent->Fill(dau1p,dau1dedx);
	////////////////////////////////////////////////////////////////////

	int  dau2q       =  diep->daughter2primaryTrack().q();
	int  dau2nHits   =  diep->daughter2primaryTrack().nHits();
	int  dau2ndedx   =  diep->daughter2primaryTrack().ndedx();
	float dau2dedx    =  diep->daughter2primaryTrack().dedx();
	float dau2pt      =  diep->daughter2primaryTrack().pt();
	float dau2px      =  diep->daughter2primaryTrack().pt()*cos(diep->daughter2primaryTrack().psi());
	float dau2py      =  diep->daughter2primaryTrack().pt()*sin(diep->daughter2primaryTrack().psi());
	float dau2pz      =  diep->daughter2primaryTrack().pt()*diep->daughter2primaryTrack().tanl();
	float dau2phidiff =  diep->daughter2bEmcMatchPhiDiff();

	TVector3 dau2(dau2px,dau2py,dau2pz);
	float dau2eta = dau2.PseudoRapidity();
	float dau2phi = dau2.Phi();
	if(dau2phi < 0.) dau2phi += 2*PI;
	float dau2p = dau2.Mag();

	float dau2beta = diep->daughter2beta();
	float dau2Energy = -999.;
	float dau2PEratio = -999.;
	float dau2EPratio = -999.;
	if(diep->daughter2emcTowerSN()>=0){
		dau2Energy = diep->daughter2bEmcTowerHit().energy();
		dau2PEratio = dau2p/dau2Energy;
		dau2EPratio = dau2Energy/dau2p;
	}

	hDau2BetaFromHltEvent->Fill(1/dau2beta);
	hDau2EnergyFromHltEvent->Fill(dau2Energy);
	hDau2dEdxvsPFromHltEvent->Fill(dau2p,dau2dedx);

	float pt = diep->pt();
	float m =  diep->invariantMass();
//	if(dau1p>2.3 && dau2p>1.5 && dau1ndedx>16 && dau2ndedx>16 && dau1PEratio<1.5 && dau1PEratio>0.5 && dau2PEratio<1.5 && dau2PEratio>0.5 && dau1phidiff>0. && dau1phidiff<0.05 && dau2phidiff > 0.&& dau2phidiff<0.05 && fabs(1/dau1beta - 1)<0.04 && fabs(1/dau2beta - 1)<0.04) {
		if(dau1q*dau2q<0) hDiElecInvMassFromHltEvent->Fill(m);
		else  hDiElecInvMassBGFromHltEvent->Fill(m);
//	}

//	return 0;

}

//________________________________________________________________________________________
void StHltQAMaker::FillTriggerReason(StHltEvent* event)
{
	
	StSPtrVecHltTriggerReason& VecTriggerReason = event->triggerReason();
	for(unsigned int i=0; i<VecTriggerReason.size(); i++){

		StHltTriggerReason *triggerReason = (StHltTriggerReason*)VecTriggerReason.at(i);
		int reasonBit = triggerReason->reasonBit();

		if(reasonBit==0x10000){
		//	StHltHighPt aHighPt = (StHltHighPt)triggerReason->reason();
			StHltHighPt aHighPt = *(StHltHighPt*)triggerReason->reason();
		    FillHighPt(aHighPt);
		}
		else if(reasonBit==0x20000){
			//StHltDiElectron aDiep = (StHltDiElectron)triggerReason->reason();
			StHltDiElectron aDiep = *(StHltDiElectron*)triggerReason->reason();
			FillDiep(aDiep);
		}
		else if(reasonBit==0x40000){
			//StHltHeavyFragment aHeavyFragment = (StHltHeavyFragment)triggerReason->reason();
			StHltHeavyFragment aHeavyFragment = *(StHltHeavyFragment*)triggerReason->reason();
			FillHF(aHeavyFragment);
		}
		else continue;

	}

}


//________________________________________________________________________________________
void StHltQAMaker::FillHighPt(StHltHighPt& highpt)
{

	StHltTrack pTrack = highpt.primaryTrack();
		if(pTrack.nHits()>25) hHighptP->Fill(pTrack.pt()); 

}

//________________________________________________________________________________________
void StHltQAMaker::FillDiep(StHltDiElectron& diep)
{
	
//	if(diep.daughter1primaryTrackSN()<0 || diep.daughter2primaryTrackSN()<0) return 0;

	int  dau1q       =  diep.daughter1primaryTrack().q();
	int  dau1nHits   =  diep.daughter1primaryTrack().nHits();
	int  dau1ndedx   =  diep.daughter1primaryTrack().ndedx();
	float dau1dedx    =  diep.daughter1primaryTrack().dedx();
	float dau1pt      =  diep.daughter1primaryTrack().pt();
	float dau1px      =  diep.daughter1primaryTrack().pt()*cos(diep.daughter1primaryTrack().psi());
	float dau1py      =  diep.daughter1primaryTrack().pt()*sin(diep.daughter1primaryTrack().psi());
	float dau1pz      =  diep.daughter1primaryTrack().pt()*diep.daughter1primaryTrack().tanl();
	float dau1phidiff =  diep.daughter1bEmcMatchPhiDiff();

	TVector3 dau1(dau1px,dau1py,dau1pz);
	float dau1eta = dau1.PseudoRapidity();
	float dau1phi = dau1.Phi();
	if(dau1phi < 0.) dau1phi += 2*PI;
	float dau1p = dau1.Mag();

	float dau1beta = diep.daughter1beta();
	float dau1Energy = -999.;
	float dau1PEratio = -999.;
	float dau1EPratio = -999.;
	if(diep.daughter1emcTowerSN()>=0){
		dau1Energy = diep.daughter1bEmcTowerHit().energy();
		dau1PEratio = dau1p/dau1Energy; 
		dau1EPratio = dau1Energy/dau1p;
	}

    hDau1Beta->Fill(1/dau1beta);
	hDau1Energy->Fill(dau1Energy);
	hDau1dEdxvsP->Fill(dau1p,dau1dedx);
	////////////////////////////////////////////////////////////////////

	int  dau2q       =  diep.daughter2primaryTrack().q();
	int  dau2nHits   =  diep.daughter2primaryTrack().nHits();
	int  dau2ndedx   =  diep.daughter2primaryTrack().ndedx();
	float dau2dedx    =  diep.daughter2primaryTrack().dedx();
	float dau2pt      =  diep.daughter2primaryTrack().pt();
	float dau2px      =  diep.daughter2primaryTrack().pt()*cos(diep.daughter2primaryTrack().psi());
	float dau2py      =  diep.daughter2primaryTrack().pt()*sin(diep.daughter2primaryTrack().psi());
	float dau2pz      =  diep.daughter2primaryTrack().pt()*diep.daughter2primaryTrack().tanl();
	float dau2phidiff =  diep.daughter2bEmcMatchPhiDiff();

	TVector3 dau2(dau2px,dau2py,dau2pz);
	float dau2eta = dau2.PseudoRapidity();
	float dau2phi = dau2.Phi();
	if(dau2phi < 0.) dau2phi += 2*PI;
	float dau2p = dau2.Mag();

	float dau2beta = diep.daughter2beta();
	float dau2Energy = -999.;
	float dau2PEratio = -999.;
	float dau2EPratio = -999.;
	if(diep.daughter2emcTowerSN()>=0){
		dau2Energy = diep.daughter2bEmcTowerHit().energy();
		dau2PEratio = dau2p/dau2Energy;
		dau2EPratio = dau2Energy/dau2p;
	}

    hDau2Beta->Fill(1/dau2beta);
	hDau2Energy->Fill(dau2Energy);
	hDau2dEdxvsP->Fill(dau2p,dau2dedx);

	float pt = diep.pt();
	float m =  diep.invariantMass();
//    if(dau1p>2.3 && dau2p>1.5 && dau1ndedx>16 && dau2ndedx>16 && dau1PEratio<1.5 && dau1PEratio>0.5 && dau2PEratio<1.5 && dau2PEratio>0.5 && dau1phidiff>0. && dau1phidiff<0.05 && dau2phidiff > 0.&& dau2phidiff<0.05 && fabs(1/dau1beta - 1)<0.04 && fabs(1/dau2beta - 1)<0.04) {
		if(dau1q*dau2q<0) hDiElecInvMass->Fill(m);
        else  hDiElecInvMassBG->Fill(m);
//	}

//	return 0;

}

//________________________________________________________________________________________
void StHltQAMaker::FillHF(StHltHeavyFragment& heavyfragment)
{

//	StHltTrack        gTrack   =  (StHltTrack&)   heavyfragment.globalTrack();
//	StHltTrack        pTrack   =  (StHltTrack&)   heavyfragment.primaryTrack();
//	StHltBTofHit      bTofHit  =  (StHltBTofHit&) heavyfragment.bTofHit();
//	StHltBEmcTowerHit bEmcHit  =  (StHltBEmcTowerHit&) heavyfragment.bEmcTowerHit();

	int nHits = heavyfragment.primaryTrack().nHits();
	int ndedx = heavyfragment.primaryTrack().ndedx();
	int q = heavyfragment.primaryTrack().q();
	float dedx = heavyfragment.primaryTrack().dedx();

	float px  =  cos(heavyfragment.primaryTrack().psi())*heavyfragment.primaryTrack().pt();
	float py  =  sin(heavyfragment.primaryTrack().psi())*heavyfragment.primaryTrack().pt();
	float pz  =  heavyfragment.primaryTrack().tanl()*heavyfragment.primaryTrack().pt();
	TVector3 mom(px,py,pz);
	float p = mom.Mag();
	if(nHits >= 20 && ndedx >= 15) hHFdEdxvsP->Fill(p*q,dedx);

}

//________________________________________________________________________________________
Int_t StHltQAMaker::Finish() 
{
	mOutFile->cd();

	hvertexX->Write();
	hvertexY->Write();
	hvertexZ->Write();
	hlmvertexX->Write();
	hlmvertexY->Write();
	hlmvertexZ->Write();
	hPvpdVzvslmVz->Write();

	hnHits->Write();
	hndedx->Write();
	hgPt->Write();
	hgEta->Write();
	hgPhi->Write();
	hgChi2Xy->Write();
	hgChi2Z->Write();
	hgMult->Write();
	hDcaXy->Write();
	hDcaZ->Write();
	hgdEdxvsP->Write();

	hpPt->Write();
	hpEta->Write();
	hpPhi->Write();
	hpChi2Xy->Write();
	hpChi2Z->Write();
	hpMult->Write();
	hpdEdxvsP->Write();

	hEmcPhiDiff->Write();
	hEmczEdge->Write();
    hEmcEnergyfromNode->Write();
	hEmcEtaPhifromNode->Write();
	hEmcSoftIdfromNode->Write();
	hEmcdaqIdfromNode->Write();
	hEmcEnergy->Write();
	hEmcEtaPhi->Write();
	hEmcSoftId->Write();
	hEmcdaqId->Write();

	hLocalZ->Write();
	hLocalY->Write();
	hInBetaVsPrimPt->Write();
	hInBetaVsGlobPt->Write();
	hMatchIdvsFireId->Write();
	hTrayIDvsTrgTime->Write();
	hChannelID->Write();

	hHFdEdxvsP->Write();
	hHighptP->Write();
	hHFdEdxvsPFromHltEvent->Write();
	hHighptPFromHltEvent->Write();

    hDau1Energy->Write();
	hDau2Energy->Write();
    hDau1dEdxvsP->Write();
	hDau2dEdxvsP->Write();

	hDau1Beta->Write();
	hDau2Beta->Write();
	hDiElecInvMass->Write();
	hDiElecInvMassBG->Write();

    hDau1EnergyFromHltEvent->Write();
	hDau2EnergyFromHltEvent->Write();
    hDau1dEdxvsPFromHltEvent->Write();
	hDau2dEdxvsPFromHltEvent->Write();

	hDau1BetaFromHltEvent->Write();
	hDau2BetaFromHltEvent->Write();
	hDiElecInvMassFromHltEvent->Write();
	hDiElecInvMassBGFromHltEvent->Write();

	mOutFile->Close();
	return StMaker::Finish();
}

