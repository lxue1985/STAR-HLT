/****************************************************************************
**
**    Author: Liang Xue
**    Description: HLT QA Maker to do the HLT information propaganda from 
**    online RTS_Reader to offline StEvent
**
*****************************************************************************/

#ifndef STAR_StHltQAMaker
#define STAR_StHltQAMaker

#define PI 3.1415926

#include "StMaker.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"


class StEvent;
class StHltEvent;
class StHltHighPt;
class StHltHeavyFragment;
class StHltDiElectron;

class StHltQAMaker : public StMaker
{
	private:
		StEvent*                  mStEvent;                   //! pointer to StEvent
	    
        
		TFile *mOutFile;

        /////// QA histograms////////
		//----event level------//
		TH1F*           hvertexX;
		TH1F*           hvertexY;
		TH1F*           hvertexZ;
		TH1F*           hlmvertexX;
		TH1F*           hlmvertexY;
		TH1F*           hlmvertexZ;
		TH2F*           hPvpdVzvslmVz;

        //----global track-----//
		TH1I*           hnHits;
		TH1I*           hndedx;
		TH1F*           hgPt;
		TH1F*           hgEta;
		TH1F*           hgPhi;
		TH1F*           hgChi2Xy;
		TH1F*           hgChi2Z;
		TH1I*           hgMult;
		TH1F*           hDcaXy;
		TH1F*           hDcaZ;
		TH2F*           hgdEdxvsP;

        //----primary track-----//
		TH1F*           hpPt;
		TH1F*           hpEta;
		TH1F*           hpPhi;
		TH1F*           hpChi2Xy;
		TH1F*           hpChi2Z;
		TH1I*           hpMult;
		TH2F*           hpdEdxvsP;

		//------------Emc -----------//
		TH1D*           hEmcPhiDiff;
		TH1D*           hEmczEdge;
		TH1F*           hEmcEnergyfromNode;
		TH2F*           hEmcEtaPhifromNode;
		TH1I*           hEmcSoftIdfromNode;
		TH1I*           hEmcdaqIdfromNode;
		TH1F*           hEmcEnergy;
		TH2F*           hEmcEtaPhi;
		TH1I*           hEmcSoftId;
		TH1I*           hEmcdaqId;
		
		//------------Tof -----------//
		TH1F*           hLocalZ;
		TH1F*           hLocalY;
		TH2F*           hInBetaVsPrimPt;
		TH2F*           hInBetaVsGlobPt;
		TH2F*           hMatchIdvsFireId;
		TH2F*           hTrayIDvsTrgTime;
		TH1F*           hChannelID;

		TH2F*           hHFdEdxvsP;
		TH2F*           hHFdEdxvsPFromHltEvent;

		TH1F*           hHighptP;
		TH1F*           hHighptPFromHltEvent;

		TH2F*           hDau1dEdxvsP;
		TH1F*           hDau1Energy;
		TH2F*           hDau2dEdxvsP;
		TH1F*           hDau2Energy;
		TH1F*           hDau1Beta;
		TH1F*           hDau2Beta;

		TH1F*           hDiElecInvMass;
		TH1F*           hDiElecInvMassBG;

		TH2F*           hDau1dEdxvsPFromHltEvent;
		TH1F*           hDau1EnergyFromHltEvent;
		TH2F*           hDau2dEdxvsPFromHltEvent;
		TH1F*           hDau2EnergyFromHltEvent;
		TH1F*           hDau1BetaFromHltEvent;
		TH1F*           hDau2BetaFromHltEvent;

		TH1F*           hDiElecInvMassFromHltEvent;
		TH1F*           hDiElecInvMassBGFromHltEvent;

	public:
		StHltQAMaker();
		virtual ~StHltQAMaker();

		virtual void  Clear(Option_t *option="" );
		virtual Int_t Init(); // called once at the beginning of your job
		virtual Int_t InitRun(int runnumber);
		virtual Int_t Make(); // invoked for every event
		virtual Int_t Finish(); // called once at the end

		virtual const char *GetCVS() const {
			static const char cvs[]="Tag $Name:  $ $Id: StHltQAMaker.h,v 1.1 2010/09/23 16:05:58 Liang Xue Exp $ built "__DATE__" "__TIME__ ;
			return cvs;
		}

		void FillEvent(StHltEvent* event);
		void FillGlobalTrack(StHltEvent* event);
		void FillPrimaryTrack(StHltEvent* event);
		void FillNode(StHltEvent* event);
		void FillEmc(StHltEvent* event);
		void FillTof(StHltEvent* event);
		void FillVpd(StHltEvent* event);

		void FillTriggerParticleFromHltEvent(StHltEvent* event);
		void FillHighPt(StHltHighPt* highpt);
		void FillHF(StHltHeavyFragment* heavyfrag);
		void FillDiep(StHltDiElectron* diep);

		void FillTriggerReason(StHltEvent* event);
		void FillHighPt(StHltHighPt& highpt);
		void FillHF(StHltHeavyFragment& heavyfrag);
		void FillDiep(StHltDiElectron& diep);

		ClassDef(StHltQAMaker,1)

};

#endif
