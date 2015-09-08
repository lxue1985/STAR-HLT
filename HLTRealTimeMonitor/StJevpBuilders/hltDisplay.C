#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "string"
#include <unistd.h>
#include <sys/stat.h> // for current dir
#include <math.h>

#include "DAQ_READER/daq_dta.h"
#include "DAQ_L3/daq_l3.h"

#include <TStyle.h> 
#include "TVector3.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <TProfile.h>
#include "TF1.h"

#include <DAQ_HLT/daq_hlt.h>
#include <DAQ_SC/daq_sc.h>
#include "DAQ_TRG/daq_trg.h"
#include "HLTFormats.h"
#include "DAQ_READER/daq_det.h"
#include "daqFormats.h"

#include "hltDisplay.h"

using namespace std;


void StHltDisplay::startrun(daqReader *rdr)
{
	sprintf(Currentrun,"%s/run10_hlt_%d_current_hist.root",Destindir,rdr->run);
	sprintf(CurrentNtuple,"%s/run10_hlt_%d_current_ntuple.root",Destindir,rdr->run);
	runnumber = rdr->run;

	for(int i=0;i<55;i++){
		getPlotByIndex(i)->getHisto(0)->histo->Reset();
	}
	hDiElectronInvMassTpxEmcBG->Reset();
	hDiElectronInvMassFullRangeBG->Reset();
	hDiElectronInvMassCutBG->Reset();
	hHFM_dEdx->Reset();
	//	hMatchannel3D->Reset();

	printf("Start run : %d  \n", rdr->run) ;
};

void StHltDisplay::stoprun(daqReader *rdr)
{
	gStyle->SetOptStat(000000);
	gStyle->SetOptFit(111);

	hDiElectronInvMassTpxEmc->SetLineColor(4);
	hDiElectronInvMassFullRange->SetLineColor(4);
	hDiElectronInvMassCut->SetLineColor(4);

	//----------lndedx pion----------//
	float low = -13.12;
	float high = -12.8;
	TF1 *fit = new TF1("fit","gaus",low,high);
	fit->SetParName(0,"Apt");
	fit->SetParName(1,"Pos");
	fit->SetParName(2,"Sig");
	//	  fit->SetParameter(0,10000);
	fit->SetParameter(1,-12.92);
	fit->SetParameter(2,0.08);
	hLn_dEdx->Fit(fit,"EMR") ;

	if(iBin%100==0){
		hBeamX->Reset();
		hBeamY->Reset();
		hMeanDcaXy->Reset();
		hInnerGain->Reset();
		hOuterGain->Reset();
		iBin = 0;
	}
	iBin ++;

	TF1 *func = new TF1("func","gaus",-6.,6.); 
	func->SetParName(0,"Apt");
	func->SetParName(1,"Mean");
	func->SetParName(2,"Sigma");
	func->SetParameter(1,0.);
	func->SetParameter(2,0.4);

	int maxBin = hDcaXy->GetMaximumBin();
	double maxVal = -6. + 0.1*maxBin; 
	hDcaXy->Fit(func,"EMR","",maxVal-1.8,maxVal+1.8) ;

	double meanpar = func->GetParameter(1);
	double errpar = func->GetParError(1);

	//for beamX beamY and gain parameters mean DcaXy//
	char num[60];
	sprintf(num,"%i",runnumber);
	/*	
		string strday(num);
		string strrun(num);
		int lenday = strday.length();
		strday.replace(lenday-3,3,"");
		strday.replace(0,2,"");
		int lenrun = strrun.length();
		strrun.replace(0,5,"");
		int tmp_daynum = atoi(strday.c_str());
		int tmp_runnum = atoi(strrun.c_str()); 
	 */
	printf("run number :%d iBin :%d BeamX :%f BeamY :%f innerGainSector :%f outerGainSector :%f \n",runnumber,iBin,BeamX,BeamY,innerGainPara,outerGainPara);

//for(int i=1;i<100;i++){

	if((iBin-1)%3==0){
		hBeamX->GetXaxis()->SetBinLabel(iBin,num);
		hBeamY->GetXaxis()->SetBinLabel(iBin,num);
		hMeanDcaXy->GetXaxis()->SetBinLabel(iBin,num);
		hInnerGain->GetXaxis()->SetBinLabel(iBin,num);
		hOuterGain->GetXaxis()->SetBinLabel(iBin,num);
	}

/*  if(i%3==0){
		hBeamX->GetXaxis()->SetBinLabel(i,num);
		hBeamY->GetXaxis()->SetBinLabel(i,num);
		hMeanDcaXy->GetXaxis()->SetBinLabel(i,num);
		hInnerGain->GetXaxis()->SetBinLabel(i,num);
		hOuterGain->GetXaxis()->SetBinLabel(i,num);
	}
}
*/
	hBeamX->SetBinContent(iBin,BeamX);
	hBeamY->SetBinContent(iBin,BeamY);
	if(hDcaXy->GetEntries()>100){
		hMeanDcaXy->SetBinContent(iBin,meanpar);
		//hMeanDcaXy->SetBinError(iBin,errpar);
	}
	if(innerGainPara>0)hInnerGain->SetBinContent(iBin,innerGainPara);
	if(outerGainPara>0)hOuterGain->SetBinContent(iBin,outerGainPara);

    char OutParas[256];
	sprintf(OutParas,"%s/%d.dat",Destindir,runnumber);
	ofstream outstream;
	outstream.open(OutParas);
	outstream<<"beamX"<<"    "<<BeamX<<endl;
	outstream<<"beamY"<<"    "<<BeamY<<endl;
	outstream<<"innerGain"<<"    "<<innerGainPara<<endl;
	outstream<<"outerGain"<<"    "<< outerGainPara<<endl;
	outstream<<"dcaXy"<<"    "<<meanpar <<endl;
	outstream.close();

	WriteHistogram(Currentrun,CurrentNtuple);

	printf("Stop run = %d \n ", rdr->run) ;
};

///---------------fill QA plots-------------///
void StHltDisplay::event(daqReader *rdr)
{
	int triggerBitHighPt = 0x10000;
	int triggerBitDiElectron = 0x20000;
	int triggerBitHeavyFragment = 0x40000;
	int triggerBitAllEvents = 0x80000;
	int triggerBitRandomEvents = 0x100000;
	int triggerBitBesgoodEvents = 0x200000;

	daq_dta *dd ;
	dd = rdr->det("hlt")->get("gl3") ;
	if(!dd) return ;

	HLT_EVE  *hlt_eve; HLT_TOF  *hlt_tof; HLT_PVPD *hlt_pvpd; HLT_EMC  *hlt_emc; HLT_GT   *hlt_gt;
	HLT_PT   *hlt_pt;  HLT_NODE *hlt_node; HLT_HIPT *hlt_hipt;HLT_DIEP *hlt_diep;HLT_HF *hlt_hf;

	while(dd && dd->iterate()){
		hlt_gl3_t *hlt = (hlt_gl3_t *) dd->Void ;

		if(strcmp(hlt->name,"HLT_EVE")==0) hlt_eve = (HLT_EVE *)hlt->data;
		else if(strcmp(hlt->name,"HLT_TOF")==0) hlt_tof = (HLT_TOF *)hlt->data;
		else if(strcmp(hlt->name,"HLT_PVPD")==0) hlt_pvpd = (HLT_PVPD *)hlt->data;
		else if(strcmp(hlt->name,"HLT_EMC")==0) hlt_emc = (HLT_EMC *)hlt->data;
		else if(strcmp(hlt->name,"HLT_GT")==0) hlt_gt = (HLT_GT *)hlt->data;
		else if(strcmp(hlt->name,"HLT_PT")==0) hlt_pt = (HLT_PT *)hlt->data;
		else if(strcmp(hlt->name,"HLT_NODE")==0) hlt_node = (HLT_NODE *)hlt->data;
		else if(strcmp(hlt->name,"HLT_HIPT")==0) hlt_hipt = (HLT_HIPT *)hlt->data;
		else if(strcmp(hlt->name,"HLT_DIEP")==0) hlt_diep = (HLT_DIEP *)hlt->data;
		else if(strcmp(hlt->name,"HLT_HF")==0) hlt_hf = (HLT_HF *)hlt->data;
	}

	printf("sequence %d: gTracks %d: tofHits: %d: pvpdHits %d: eePairs %d: emcTowers %d: highPt %d : heavyFrag %d \n",rdr->seq, hlt_gt->nGlobalTracks, hlt_tof->nTofHits, hlt_pvpd->nPvpdHits, hlt_diep->nEPairs, hlt_emc->nEmcTowers, hlt_hipt->nHighPt, hlt_hf->nHeavyFragments) ;

	int decision = hlt_eve->hltDecision;
	int upc = 0x1; int mtd = 0x4000; int npehtnozdc = 0x40; int vpdmb = 0x20; int central = 0x4; 
	int npeht = 0x40000; int npe11 = 0x800; int npe15 = 0x10000; int npe18 = 0x20000; int atom = 0x8;

	// if( !(decision & triggerBitBesgoodEvents) ) return;  //only monitor HLT good event in BES

	if(rdr->daqbits & upc)
	{
		if(decision & triggerBitHighPt) hUpc->Fill(0.);
		if(decision & triggerBitDiElectron) hUpc->Fill(1.);
		if(decision & triggerBitHeavyFragment) hUpc->Fill(2.);
	}
	if(rdr->daqbits & mtd)
	{
		if(decision & triggerBitHighPt) hMtd->Fill(0.);
		if(decision & triggerBitDiElectron) hMtd->Fill(1.);
		if(decision & triggerBitHeavyFragment) hMtd->Fill(2.);
	}
	if(rdr->daqbits & npehtnozdc)
	{
		if(decision & triggerBitHighPt) hNpeHt_25_NoZdc->Fill(0.);
		if(decision & triggerBitDiElectron) hNpeHt_25_NoZdc->Fill(1.);
		if(decision & triggerBitHeavyFragment) hNpeHt_25_NoZdc->Fill(2.);
	}
	if(rdr->daqbits & vpdmb)
	{
		if(decision & triggerBitHighPt) hVpdMb->Fill(0.);
		if(decision & triggerBitDiElectron) hVpdMb->Fill(1.);
		if(decision & triggerBitHeavyFragment) hVpdMb->Fill(2.);
	}
	if(rdr->daqbits & central)
	{
		if(decision & triggerBitHighPt) hCentral->Fill(0.);
		if(decision & triggerBitDiElectron) hCentral->Fill(1.);
		if(decision & triggerBitHeavyFragment) hCentral->Fill(2.);
	}
	if(rdr->daqbits & npeht)
	{
		if(decision & triggerBitHighPt) hNpeHt_25->Fill(0.);
		if(decision & triggerBitDiElectron) hNpeHt_25->Fill(1.);
		if(decision & triggerBitHeavyFragment) hNpeHt_25->Fill(2.);
	}
	if(rdr->daqbits & npe11)
	{
		if(decision & triggerBitHighPt) hNpe->Fill(0.);
		if(decision & triggerBitDiElectron) hNpe->Fill(1.);
		if(decision & triggerBitHeavyFragment) hNpe->Fill(2.);
	}
	if(rdr->daqbits & npe15)
	{
		if(decision & triggerBitHighPt) hNpe->Fill(3.);
		if(decision & triggerBitDiElectron) hNpe->Fill(4.);
		if(decision & triggerBitHeavyFragment) hNpe->Fill(5.);
	}
	if(rdr->daqbits & npe18)
	{
		if(decision & triggerBitHighPt) hNpe->Fill(6.);
		if(decision & triggerBitDiElectron) hNpe->Fill(7.);
		if(decision & triggerBitHeavyFragment) hNpe->Fill(8.);
	}
	if(rdr->daqbits & atom)
	{
		if(decision & triggerBitHighPt) hAtomcule->Fill(0.);
		if(decision & triggerBitDiElectron) hAtomcule->Fill(1.);
		if(decision & triggerBitHeavyFragment) hAtomcule->Fill(2.);
	}


	if(decision) hEvtsAccpt->Fill(0.);  // event statistic
	//          if(decision & triggerBitBesgoodEvents) EvtsAccpt->Fill(1.);
	//          if(decision & triggerBitRandomEvents) EvtsAccpt->Fill(2.);

	if(decision & triggerBitHighPt) hEvtsAccpt->Fill(1.);
	if(decision & triggerBitDiElectron) hEvtsAccpt->Fill(2.);
	if(decision & triggerBitHeavyFragment) hEvtsAccpt->Fill(3.);
	if(decision & triggerBitRandomEvents) hEvtsAccpt->Fill(4.); 

	//----------fill events-----------//
	float vertX = hlt_eve->vertexX ;
	float vertY = hlt_eve->vertexY ;
	float vertZ = hlt_eve->vertexZ ;
	float lmvertX = hlt_eve->lmVertexX ;
	float lmvertY = hlt_eve->lmVertexY ;
	float lmvertZ = hlt_eve->lmVertexZ ;
	if( hlt_eve->version >= 0x20100216){
		innerGainPara = hlt_eve->innerSectorGain ;
		outerGainPara = hlt_eve->outerSectorGain ;
	}

	float VzVpd =  hlt_eve->vpdVertexZ ;

	BeamX = lmvertX ;
	BeamY = lmvertY ;

	hVertexX->Fill(vertX);
	hVertexY->Fill(vertY);
	hVertexZ->Fill(vertZ);
	hLm_VertexX->Fill(lmvertX);
	hLm_VertexY->Fill(lmvertY);
	hLm_VertexZ->Fill(lmvertZ); 
	hVzvpd_Vz->Fill(VzVpd,  lmvertZ);
	hVzDiff->Fill(VzVpd - lmvertZ);


	//------------------fill tof hits--------------//
	for(u_int i=0 ; i < hlt_tof->nTofHits ; i++){
		short trayId     = hlt_tof->tofHit[i].trayId;
		short channel    = hlt_tof->tofHit[i].channel;
		float tdc         = hlt_tof->tofHit[i].tdc;
		//	  float tof         = hlt_tof->tofHit[i].tof; 
		float triggertime = hlt_tof->tofHit[i].triggertime;
		hTrayID_TrgTime->Fill(trayId, tdc - triggertime);
		hchannelID->Fill(channel); 
	}

	//-----------------fill pvpd hit------------//
	for(u_int i=0 ; i < hlt_pvpd->nPvpdHits ; i++){
		short trayId     = hlt_pvpd->pvpdHit[i].trayId;
		//		  short channel    = hlt_pvpd->pvpdHit[i].channel;
		float tdc         = hlt_pvpd->pvpdHit[i].tdc;
		//		  float tof         = hlt_pvpd->pvpdHit[i].tof;
		float triggertime = hlt_pvpd->pvpdHit[i].triggertime;
		hTrayID_TrgTime->Fill(trayId, tdc - triggertime);
	}

	//--------------fill Emc--------------//
	for(u_int i=0 ; i < hlt_emc->nEmcTowers ; i++){
		//		  int adc = hlt_emc->emcTower[i].adc;
		float energy     = hlt_emc->emcTower[i].energy;
		float phi   = hlt_emc->emcTower[i].phi;
		float  eta   = hlt_emc->emcTower[i].eta;
		//		  float  z     = hlt_emc->emcTower[i].z; 
		int softId  = hlt_emc->emcTower[i].softId;
		int daqId   = hlt_emc->emcTower[i].daqId;
		hTowerEnergy->Fill(energy);//run
		hTowerDaqId->Fill(daqId);  //run
		hTowerSoftId->Fill(softId);  //run
		hTowerEtaPhi->Fill(phi,eta);  //run
	}


	//---------------global track-------------//

	for(u_int i=0 ; i < hlt_gt->nGlobalTracks ; i++)
	{
		int nHits = hlt_gt->globalTrack[i].nHits;
		int ndedx = hlt_gt->globalTrack[i].ndedx;
		hnhits->Fill(nHits);
		hnDedx->Fill(ndedx);

		if(hlt_gt->globalTrack[i].flag < 0.) continue;
		//		  if(hlt_gt->globalTrack[i].q != +1 && hlt_gt->globalTrack[i].q != -1) continue;
		float pt = hlt_gt->globalTrack[i].pt;
		float px = cos(hlt_gt->globalTrack[i].psi)*hlt_gt->globalTrack[i].pt ;
		float py = sin(hlt_gt->globalTrack[i].psi)*hlt_gt->globalTrack[i].pt ;
		float pz = hlt_gt->globalTrack[i].tanl*hlt_gt->globalTrack[i].pt ;
		int  q  = hlt_gt->globalTrack[i].q  ;

		TVector3 mom(px,py,pz);
		float eta=mom.PseudoRapidity();
		float phi = mom.Phi();
		if(phi<0.0) phi += twopi;
		float p = mom.Mag();
		float dedx = hlt_gt->globalTrack[i].dedx;

		hGlob_Eta->Fill(eta);
		if(nHits >= 25 && fabs(eta) < 1.){
			hGlob_Pt->Fill(pt);
			hGlob_Phi->Fill(phi);
		}
		if(nHits >= 20 && ndedx >= 15) hGlob_dEdx->Fill(p*q,dedx);
		if(nHits >= 20 && ndedx >= 20) {
			hdEdx->Fill(p*q,dedx) ; // fhor HF reference
		}
	}//loop of tracks


	//--------------primary tracks---------------//
	int count = 0 ;
	for(u_int i=0 ; i < hlt_node->nNodes ; i++){    // primary track loop
		int globalTrackSN  =  hlt_node->node[i].globalTrackSN;
		int primaryTrackSN =  hlt_node->node[i].primaryTrackSN;

		hlt_track GTrack = hlt_gt->globalTrack[globalTrackSN] ;
		double dcaX = GTrack.r0*cos(GTrack.phi0) - hlt_eve->lmVertexX;
		double dcaY = GTrack.r0*sin(GTrack.phi0) - hlt_eve->lmVertexY;
		double cross = dcaX*sin(GTrack.psi) - dcaY*cos(GTrack.psi);
		double theSign = (cross>=0) ? 1. : -1.;
		double dcaXy = theSign*sqrt(pow(dcaX,2)+pow(dcaY,2));
		double dcaZ = GTrack.z0 - hlt_eve->lmVertexZ;
		hDcaXy->Fill(dcaXy); 
		hDcaZ->Fill(dcaZ); 
		if(primaryTrackSN < 0 ) continue ;
		count++ ;

		hlt_track PTrack = hlt_pt->primaryTrack[primaryTrackSN] ;
		if(PTrack.flag < 0.) continue;
		//		  if(PTrack.q != +1 && PTrack.q != -1) continue;

		int nHits = PTrack.nHits;
		int ndedx = PTrack.ndedx;
		int  q  = PTrack.q  ;
		float pt = PTrack.pt;
		float px = cos(PTrack.psi)*PTrack.pt ;
		float py = sin(PTrack.psi)*PTrack.pt ;
		float pz = PTrack.tanl*PTrack.pt ;

		TVector3 mom(px,py,pz);
		float eta=mom.PseudoRapidity();
		float phi = mom.Phi();
		if(phi<0.0) phi += twopi;
		float p = mom.Mag();
		float dedx = PTrack.dedx;

		hPrim_Eta->Fill(eta);
		if(nHits >= 25 && fabs(eta) < 1. ){
			hPrim_Pt->Fill(pt);
			hPrim_Phi->Fill(phi);
		}
		if(nHits >= 20 && ndedx >= 15) hPrim_dEdx->Fill(p*q,dedx);
		if(nHits >= 20 && ndedx >= 15 && p >= 0.5 && p <= 0.6 ) hLn_dEdx->Fill(log(dedx));
	}

	primaryTracks = count;
	hglobalMult->Fill(hlt_gt->nGlobalTracks);
	hprimaryMult->Fill(count);


	//------------------fille nodes-----------------//
	for(u_int i=0 ; i < hlt_node->nNodes ; i++){
		int primaryTrackSN =  hlt_node->node[i].primaryTrackSN;
		int tofHitSN       =  hlt_node->node[i].tofHitSN;
		int emcTowerSN     =  hlt_node->node[i].emcTowerSN;

		hlt_track NTrack =  hlt_pt->primaryTrack[primaryTrackSN] ;

		float px = NTrack.pt*cos(NTrack.psi);
		float py = NTrack.pt*sin(NTrack.psi);
		float pz = NTrack.tanl*NTrack.pt;
		float p = sqrt(px*px + py*py +pz*pz) ;

		if(tofHitSN >= 0) {
			float localY = hlt_node->node[i].localY;
			float localZ = hlt_node->node[i].localZ;
			float beta   = hlt_node->node[i].beta;
			//		  float tof    = hlt_node->node[i].tof;
			int  projChannel = hlt_node->node[i].projChannel;
			hLocalZ->Fill(localZ);
			hLocalY->Fill(localY);
			if(primaryTrackSN >= 0) {
				hInverseBeta->Fill(p, 1/beta);
				for(u_int j=0 ; j < hlt_tof->nTofHits ; j++){
					int Proj_trayId = hlt_tof->tofHit[tofHitSN].trayId ;
					int fire_trayId  = hlt_tof->tofHit[j].trayId ;
					if(Proj_trayId == fire_trayId ){
						hMatchId_fiberId->Fill(projChannel, hlt_tof->tofHit[j].channel);
						//	hMatchannel3D->Fill(projChannel, hlt_tof->tofHit[j].channel,Proj_trayId);
					}
				}
			}
		}

		if(emcTowerSN >= 0 && NTrack.nHits>20 && NTrack.ndedx>15) {
			double emcMatchPhiDiff = hlt_node->node[i].emcMatchPhiDiff;
			double emcMatchZEdge   = hlt_node->node[i].emcMatchZEdge;
			hMatchPhi_Diff->Fill(emcMatchPhiDiff);
			if(emcMatchZEdge > 0.) hzEdge->Fill(emcMatchZEdge);
		}
	}

	//-------------------Heavy Fragment-----------------//
	for(u_int i=0 ; i < hlt_hf->nHeavyFragments ; i++){
		int heavyFrag_NodeSN = hlt_hf->heavyFragmentSN[i] ;
		int heavyFragmentglobSN  = hlt_node->node[heavyFrag_NodeSN].globalTrackSN ; 
		hlt_track HFtrack = hlt_gt->globalTrack[heavyFragmentglobSN] ;
		int nHits =  HFtrack.nHits ;
		int ndedx =  HFtrack.ndedx ;
		int q     =  HFtrack.q ;
		float hfpx    = HFtrack.pt*cos(HFtrack.psi) ;
		float hfpy    = HFtrack.pt*sin(HFtrack.psi) ;
		float hfpz    = HFtrack.pt*HFtrack.tanl ;
		float hfp     = sqrt(hfpx*hfpx + hfpy*hfpy + hfpz*hfpz) ;
		float hfdedx  =  HFtrack.dedx ;

		if(nHits >= 20 && ndedx >= 15) {
			hHFM_dEdx->Fill(hfp*q , hfdedx) ;
		}
	}

	//---------------  Di electrons ------------------// 
	for(u_int i=0 ; i < hlt_diep->nEPairs ; i++) {
		int Daughter1NodeSN = hlt_diep->ePair[i].dau1NodeSN ;
		int Daughter2NodeSN = hlt_diep->ePair[i].dau2NodeSN ;
		int Daughter1TrackSN = hlt_node->node[Daughter1NodeSN].primaryTrackSN ;
		int Daughter2TrackSN = hlt_node->node[Daughter2NodeSN].primaryTrackSN ;
		int Daughter1EmcSN = hlt_node->node[Daughter1NodeSN].emcTowerSN ;
		int Daughter2EmcSN = hlt_node->node[Daughter2NodeSN].emcTowerSN ;
		int Daughter1TofSN = hlt_node->node[Daughter1NodeSN].tofHitSN ;
		int Daughter2TofSN = hlt_node->node[Daughter2NodeSN].tofHitSN ;

		/////////////// ---Daughter 1-- ////////////
		if(Daughter1TrackSN < 0) continue;
		hlt_track Daughter1Track =  hlt_pt->primaryTrack[Daughter1TrackSN] ;

		float Daughter1beta = -999. ; float Daughter1phidiff = -999.;
		float Daughter1_PE_ratio = -999. ; float Daughter1_EP_ratio = -999.;

		float Daughter1q     = Daughter1Track.q ;
		float Daughter1pt    = Daughter1Track.pt ;
		float Daughter1px    = Daughter1Track.pt*cos(Daughter1Track.psi) ; 
		float Daughter1py    = Daughter1Track.pt*sin(Daughter1Track.psi) ;
		float Daughter1pz    = Daughter1Track.pt*Daughter1Track.tanl ;
		float Daughter1nHits = Daughter1Track.nHits ;
		float Daughter1dedx  = Daughter1Track.dedx ;
		int Daughter1ndedx = Daughter1Track.ndedx ;


		TVector3 Daughter1( Daughter1px, Daughter1py, Daughter1pz ) ;
		float Daughter1p = Daughter1.Mag() ;

		double dedx1E = getDedxMeanElectron(Daughter1p) ;
		float nSigma1 = log(Daughter1dedx/dedx1E)/A*sqrt(Daughter1ndedx) ;

		float Daughter1eta = Daughter1.PseudoRapidity() ;
		//  if( fabs(Daughter1eta) > 1.) continue ;
		float Daughter1phi = Daughter1.Phi() ;
		if( Daughter1phi < 0.) Daughter1phi += twopi ;

		hdEdx_P1->Fill(Daughter1p , Daughter1dedx ) ;
		if(Daughter1EmcSN >= 0){
			float Daughter1TowerEnergy = hlt_emc->emcTower[Daughter1EmcSN].energy;
			Daughter1_PE_ratio = Daughter1p/Daughter1TowerEnergy ;
			Daughter1_EP_ratio = Daughter1TowerEnergy/Daughter1p ;
			hDaughter1P_TowerEnergy->Fill(Daughter1_EP_ratio) ;
			Daughter1phidiff = hlt_node->node[Daughter1NodeSN].emcMatchPhiDiff ;
		}
		if(Daughter1TofSN >= 0.){
			Daughter1beta = hlt_node->node[Daughter1NodeSN].beta ;
		}

		//////////////// ---daughter 2--- /////////////////
		if(Daughter2TrackSN < 0.) continue;
		hlt_track Daughter2Track =  hlt_pt->primaryTrack[Daughter2TrackSN] ;
		float Daughter2phidiff = -999.;  float Daughter2beta = -999. ;
		float Daughter2_PE_ratio = -999.; float Daughter2_EP_ratio = -999.; 

		float Daughter2q     =  Daughter2Track.q  ;
		float Daughter2pt    =  Daughter2Track.pt ;
		float Daughter2px    =  Daughter2Track.pt*cos(Daughter2Track.psi) ;
		float Daughter2py    =  Daughter2Track.pt*sin(Daughter2Track.psi) ;
		float Daughter2pz    =  Daughter2Track.pt*Daughter2Track.tanl ;
		float Daughter2nHits =  Daughter2Track.nHits ;
		float Daughter2dedx  =  Daughter2Track.dedx ; 
		int Daughter2ndedx = Daughter2Track.ndedx ;

		TVector3 Daughter2( Daughter2px, Daughter2py, Daughter2pz ) ;
		float Daughter2p = Daughter2.Mag() ;

		double dedx2E = getDedxMeanElectron(Daughter2p) ;
		float nSigma2 = log(Daughter2dedx/dedx2E)/A*sqrt(Daughter2ndedx) ;

		float Daughter2eta = Daughter2.PseudoRapidity();
		//	  if( fabs( Daughter2eta ) > 1. ) continue ;
		float Daughter2phi = Daughter2.Phi() ;
		if( Daughter2phi < 0.0 ) Daughter2phi += twopi ;

		hdEdx_P2->Fill(Daughter2p , Daughter2dedx ) ;
		if(Daughter2EmcSN >= 0) {
			float Daughter2TowerEnergy = hlt_emc->emcTower[Daughter2EmcSN].energy; 
			Daughter2_PE_ratio = Daughter2p/Daughter2TowerEnergy ; 
			Daughter2_EP_ratio = Daughter2TowerEnergy/Daughter2p  ;
			hDaughter2P_TowerEnergy->Fill(Daughter2_EP_ratio) ;
			Daughter2phidiff = hlt_node->node[Daughter2NodeSN].emcMatchPhiDiff ;
		}
		if(Daughter2TofSN >= 0.){
			Daughter2beta = hlt_node->node[Daughter2NodeSN].beta ;
		}

		//////////////////  ----jpsi ----  ///////////////////
		float pt = hlt_diep->ePair[i].pt ;
		float px = cos(hlt_diep->ePair[i].psi)*hlt_diep->ePair[i].pt ;
		float py = sin(hlt_diep->ePair[i].psi)*hlt_diep->ePair[i].pt ;
		float pz = hlt_diep->ePair[i].tanl*hlt_diep->ePair[i].pt ;
		float m = hlt_diep->ePair[i].invariantMass ;

		if(Daughter1q*Daughter2q < 0.) {
			hDiElectronInvMassFullRange->Fill(m);
		}
		else{
			hDiElectronInvMassFullRangeBG->Fill(m);
		}

		if( nSigma1>-0.9 && nSigma2>-0.9 && Daughter1p > 2.3 && Daughter2p > 1.5 && Daughter1ndedx>16 && Daughter2ndedx>16 && Daughter1_PE_ratio<1.5 && Daughter1_PE_ratio>0.5 && Daughter2_PE_ratio<1.5 && Daughter2_PE_ratio>0.5 && Daughter1phidiff>0. && Daughter1phidiff<0.05 && Daughter2phidiff > 0.&& Daughter2phidiff<0.05 && fabs(1/Daughter1beta - 1)<0.04 && fabs(1/Daughter2beta - 1)<0.04){
			if( Daughter1q*Daughter2q < 0.){
				hDiElectronInvMassCut->Fill(m);
			}
			else{
				hDiElectronInvMassCutBG->Fill(m);
			}
		}

		if( nSigma1>-0.9 && nSigma2>-0.9 && Daughter1p > 2.3 && Daughter2p > 1.5 && Daughter1ndedx>16 && Daughter2ndedx>16 && Daughter1_PE_ratio<1.5 && Daughter1_PE_ratio>0.5 && Daughter2_PE_ratio<1.5 && Daughter2_PE_ratio>0.5 && Daughter1phidiff>0. && Daughter1phidiff<0.05 && Daughter2phidiff > 0.&& Daughter2phidiff<0.05){
			if( Daughter1q*Daughter2q < 0.){
				hDiElectronInvMassTpxEmc->Fill(m);
			}
			else{
				hDiElectronInvMassTpxEmcBG->Fill(m);
			}

			if(Daughter1TofSN >= 0.) hDaughter1TpxEmcInverseBeta->Fill(1/Daughter1beta); 
			if(Daughter2TofSN >= 0.) hDaughter2TpxEmcInverseBeta->Fill(1/Daughter2beta);
		}

		TLorentzVector jpsi(0,0,0,0) ;

		jpsi.SetXYZM(px,py,pz,m) ;
		float rapidity = jpsi.Rapidity() ;
		hDiLeptonRapidity->Fill(rapidity);
	}
} ;

/////////////////----fill ntuple file----//////////////////////////

void StHltDisplay::ntuple(daqReader *rdr)
{
	daq_dta *dd ;
	dd = rdr->det("hlt")->get("gl3") ;
	if(!dd) return ;

	HLT_EVE *hlt_eve; HLT_TOF *hlt_tof; HLT_PVPD *hlt_pvpd; HLT_EMC *hlt_emc; HLT_GT *hlt_gt;
	HLT_PT *hlt_pt; HLT_NODE *hlt_node; HLT_HIPT *hlt_hipt; HLT_DIEP *hlt_diep; HLT_HF *hlt_hf;

	while(dd && dd->iterate()){
		hlt_gl3_t *hlt = (hlt_gl3_t *) dd->Void ;

		if(strcmp(hlt->name,"HLT_EVE")==0) hlt_eve = (HLT_EVE *)hlt->data;
		else if(strcmp(hlt->name,"HLT_TOF")==0) hlt_tof = (HLT_TOF *)hlt->data;
		else if(strcmp(hlt->name,"HLT_PVPD")==0) hlt_pvpd = (HLT_PVPD *)hlt->data;
		else if(strcmp(hlt->name,"HLT_EMC")==0) hlt_emc = (HLT_EMC *)hlt->data;
		else if(strcmp(hlt->name,"HLT_GT")==0) hlt_gt = (HLT_GT *)hlt->data;
		else if(strcmp(hlt->name,"HLT_PT")==0) hlt_pt = (HLT_PT *)hlt->data;
		else if(strcmp(hlt->name,"HLT_NODE")==0) hlt_node = (HLT_NODE *)hlt->data;
		else if(strcmp(hlt->name,"HLT_HIPT")==0) hlt_hipt = (HLT_HIPT *)hlt->data;
		else if(strcmp(hlt->name,"HLT_DIEP")==0) hlt_diep = (HLT_DIEP *)hlt->data;
		else if(strcmp(hlt->name,"HLT_HF")==0) hlt_hf = (HLT_HF *)hlt->data;
	}


	int runID = rdr->run ;
	int eventID = rdr->seq ;
	int daqBits = rdr->daqbits ; 

	//----------fill events-----------//
	float vertX = hlt_eve->vertexX ;
	float vertY = hlt_eve->vertexY ;
	float vertZ = hlt_eve->vertexZ ;
	float lmvertX = hlt_eve->lmVertexX ;
	float lmvertY = hlt_eve->lmVertexY ;

	int globMult = hlt_gt->nGlobalTracks ;
	int primMult = primaryTracks ;

	//-------------------Heavy Fragment-----------------//
	for(u_int i=0 ; i < hlt_hf->nHeavyFragments ; i++){
		int heavyFragNodeSN = hlt_hf->heavyFragmentSN[i] ;
		int heavyFragmentglobSN  = hlt_node->node[heavyFragNodeSN].globalTrackSN ; 
		int heavyFragmentprimSN  = hlt_node->node[heavyFragNodeSN].primaryTrackSN ;
		int heavyFragEmcSN = hlt_node->node[heavyFragNodeSN].emcTowerSN ;
		int heavyFragTofSN = hlt_node->node[heavyFragNodeSN].tofHitSN ; 

		float beta = -999.; float tof = -999.; float localY = -999.; float localZ = -999.;
		float Energy = -999. ;float PhiDiff = -999.; float Zedge = 999.;

		hlt_track gHFtrack = hlt_gt->globalTrack[heavyFragmentglobSN] ;
		hlt_track pHFtrack = hlt_pt->primaryTrack[heavyFragmentprimSN] ;

		if(heavyFragTofSN >= 0){
			beta  = hlt_node->node[heavyFragNodeSN].beta;
			tof   = hlt_node->node[heavyFragNodeSN].tof ;
			localZ = hlt_node->node[heavyFragNodeSN].localZ;
			localY = hlt_node->node[heavyFragNodeSN].localY;
		}

		if(heavyFragEmcSN >= 0){
			Energy = hlt_emc->emcTower[heavyFragEmcSN].energy ; 
			PhiDiff = hlt_node->node[heavyFragNodeSN].emcMatchPhiDiff;
			Zedge  = hlt_node->node[heavyFragNodeSN].emcMatchZEdge ;
		}
		float primpx =  pHFtrack.pt*cos(pHFtrack.psi) ;
		float primpy =  pHFtrack.pt*sin(pHFtrack.psi) ;
		float primpz =  pHFtrack.pt*pHFtrack.tanl ;
		float primpt    =  pHFtrack.pt ;
		float primp = sqrt(primpx*primpx + primpy*primpy + primpz*primpz) ;

		int nHits =  gHFtrack.nHits ;
		int ndedx =  gHFtrack.ndedx ;
		int q     =  gHFtrack.q ;
		float px   = gHFtrack.pt*cos(gHFtrack.psi) ;
		float py    = gHFtrack.pt*sin(gHFtrack.psi) ;
		float pz    = gHFtrack.pt*gHFtrack.tanl ;
		float pt    =  gHFtrack.pt ;
		float dedx  = gHFtrack.dedx ;
		float tanl  =  gHFtrack.tanl ;
		float psi  =   gHFtrack.psi ;
		float chi2XY = gHFtrack.chi2[0];
		float chi2Z =  gHFtrack.chi2[1];
		//	  cout<<"beta =" << beta<<endl;

		float dcaX = gHFtrack.r0*cos(gHFtrack.phi0) - hlt_eve->lmVertexX;
		float dcaY = gHFtrack.r0*sin(gHFtrack.phi0) - hlt_eve->lmVertexY;
		float cross = dcaX*sin(gHFtrack.psi) - dcaY*cos(gHFtrack.psi);
		float theSign = (cross>=0) ? 1. : -1.;
		float globDcaXY =theSign*sqrt(pow(dcaX,2)+pow(dcaY,2)); 
		float globDcaZ =  gHFtrack.z0 - hlt_eve->lmVertexZ;

		TVector3 mom(px,py,pz) ;
		float p = mom.Mag();
		float eta = mom.PseudoRapidity() ;

		float Msquare = 0 ;
		if(heavyFragTofSN >= 0){
			Msquare = pow(p,2)*(1/beta/beta -1) ; 
		}

		//	  cout << "Msquare = " << Msquare <<endl;

		float dedxTri = getDedxMeanTriton(p) ;
		float dedxRatioTri =  log(dedx/dedxTri) ;
		float nSigmaTriton =  log(dedx/dedxTri)/A*sqrt(ndedx) ;

		float dedxHe3 = getDedxMeanHe3(p) ;
		float dedxRatioHe3 =  log(dedx/dedxHe3) ;
		float nSigmaHe3 =  log(dedx/dedxHe3)/A*sqrt(ndedx) ;

		float dedxHe4 = getDedxMeanHe4(p) ; 
		float dedxRatioHe4 =  log(dedx/dedxHe4) ;
		float nSigmaHe4 = log(dedx/dedxHe4)/A*sqrt(ndedx) ;

		float hfrag[39] = { runID, eventID, daqBits, lmvertX, lmvertY, globMult, primMult, vertX, vertY, vertZ, chi2XY, chi2Z, globDcaXY, globDcaZ, p, pt, primp, primpt, q, tanl, psi, eta, ndedx, nHits, dedx, Energy, PhiDiff, Zedge, beta, tof, localZ, localY , Msquare, dedxRatioTri, dedxRatioHe3, dedxRatioHe4, nSigmaTriton, nSigmaHe3, nSigmaHe4} ;

		heavyfrag->Fill(hfrag);
	}

	//---------------  Di electrons ------------------// 
	for(u_int i=0 ; i < hlt_diep->nEPairs ; i++) {
		int dau1NodeSN = hlt_diep->ePair[i].dau1NodeSN ;
		int dau2NodeSN = hlt_diep->ePair[i].dau2NodeSN ;
		int dau1TrackSN = hlt_node->node[dau1NodeSN].primaryTrackSN ;
		int dau2TrackSN = hlt_node->node[dau2NodeSN].primaryTrackSN ;
		int dau1EmcSN = hlt_node->node[dau1NodeSN].emcTowerSN ;
		int dau2EmcSN = hlt_node->node[dau2NodeSN].emcTowerSN ;
		int dau1TofSN = hlt_node->node[dau1NodeSN].tofHitSN ;
		int dau2TofSN = hlt_node->node[dau2NodeSN].tofHitSN ;

		/////////////////////daughter 1  /////////////////// 
		if(dau1TrackSN < 0.) continue;
		hlt_track dau1Track =  hlt_pt->primaryTrack[dau1TrackSN] ;

		float dau1localY = -999.; float dau1localZ = -999.; float dau1beta = -999.; float dau1tof= -999.;
		float dau1TowerEnergy  = -999.; float dau1Zedge = -999.; float dau1phidiff  = -999.;

		float dau1q     = dau1Track.q ;
		float dau1pt    = dau1Track.pt ;
		float dau1px    = dau1Track.pt*cos(dau1Track.psi) ;
		float dau1py    = dau1Track.pt*sin(dau1Track.psi) ;
		float dau1pz    = dau1Track.pt*dau1Track.tanl ;
		float dau1nHits = dau1Track.nHits ;
		float dau1dedx  = dau1Track.dedx ;
		int dau1ndedx = dau1Track.ndedx ;
		float dau1chi2XY = dau1Track.chi2[0] ;
		float dau1chi2Z =  dau1Track.chi2[1] ;
		float dau1tanl = dau1Track.tanl ;

		float dau1dcaX = dau1Track.r0*cos(dau1Track.phi0) - hlt_eve->lmVertexX;
		float dau1dcaY = dau1Track.r0*sin(dau1Track.phi0) - hlt_eve->lmVertexY;
		float dau1cross = dau1dcaX*sin(dau1Track.psi) - dau1dcaY*cos(dau1Track.psi);
		float dau1theSign = (dau1cross>=0) ? 1. : -1.;
		float dau1DcaXY = dau1theSign*sqrt(pow(dau1dcaX,2)+pow(dau1dcaY,2)); 
		float dau1DcaZ =  dau1Track.z0 - hlt_eve->lmVertexZ;

		if(dau1EmcSN >= 0){
			dau1TowerEnergy = hlt_emc->emcTower[dau1EmcSN].energy ;
			dau1Zedge = hlt_node->node[dau1NodeSN].emcMatchZEdge ;
			dau1phidiff =  hlt_node->node[dau1NodeSN].emcMatchPhiDiff ;
		}
		if(dau1TofSN >= 0.){
			dau1localY =  hlt_node->node[dau1NodeSN].localY ;
			dau1localZ =  hlt_node->node[dau1NodeSN].localZ;
			dau1beta   =  hlt_node->node[dau1NodeSN].beta ;
			dau1tof    =  hlt_node->node[dau1NodeSN].tof ;
		}


		TVector3 dau1( dau1px, dau1py, dau1pz ) ;

		double dau1p = dau1.Mag() ;
		double dedx1E = getDedxMeanElectron(dau1p) ;
		float nSigma1 = log(dau1dedx/dedx1E)/A*sqrt(dau1ndedx) ;
		float dau1eta = dau1.PseudoRapidity() ;
		float dau1phi = dau1.Phi() ; 


		///////////////////// daughter 2 ////////////////////
		if(dau2TrackSN < 0.) continue;
		hlt_track dau2Track =  hlt_pt->primaryTrack[dau2TrackSN] ;

		float dau2localY = -999.; float dau2localZ = -999.; float dau2beta = -999.;float dau2tof = -999.;
		float dau2TowerEnergy = -999.; float dau2Zedge = -999.; float dau2phidiff = -999.;

		float dau2q     =  dau2Track.q  ;
		float dau2pt    =  dau2Track.pt ;
		float dau2px    =  dau2Track.pt*cos(dau2Track.psi) ;
		float dau2py    =  dau2Track.pt*sin(dau2Track.psi) ;
		float dau2pz    =  dau2Track.pt*dau2Track.tanl ;
		float dau2nHits =  dau2Track.nHits ;
		float dau2dedx  =  dau2Track.dedx ;
		int dau2ndedx = dau2Track.ndedx ;
		float dau2chi2XY = dau2Track.chi2[0] ;
		float dau2chi2Z =  dau2Track.chi2[1] ;
		float dau2tanl =  dau2Track.tanl ;

		float dau2dcaX = dau2Track.r0*cos(dau2Track.phi0) - hlt_eve->lmVertexX;
		float dau2dcaY = dau2Track.r0*sin(dau2Track.phi0) - hlt_eve->lmVertexY;
		float dau2cross = dau2dcaX*sin(dau2Track.psi) - dau2dcaY*cos(dau2Track.psi);
		float dau2theSign = (dau2cross>=0) ? 1. : -1.;
		float dau2DcaXY = dau2theSign*sqrt(pow(dau2dcaX,2)+pow(dau2dcaY,2));
		float dau2DcaZ =  dau2Track.z0 - hlt_eve->lmVertexZ ;

		if(dau2EmcSN >= 0 ){
			dau2TowerEnergy = hlt_emc->emcTower[dau2EmcSN].energy ;
			dau2Zedge = hlt_node->node[dau2NodeSN].emcMatchZEdge ;
			dau2phidiff =  hlt_node->node[dau2NodeSN].emcMatchPhiDiff ;
		}
		if( dau2TofSN >= 0){
			dau2localY =  hlt_node->node[dau2NodeSN].localY ;
			dau2localZ =  hlt_node->node[dau2NodeSN].localZ;
			dau2beta   =  hlt_node->node[dau2NodeSN].beta ;
			dau2tof    =  hlt_node->node[dau2NodeSN].tof ;
		}

		TVector3 dau2( dau2px, dau2py, dau2pz ) ;

		double dau2p = dau2.Mag() ;
		double dedx2E = getDedxMeanElectron(dau2p) ;
		float nSigma2 = log(dau2dedx/dedx2E)/A*sqrt(dau2ndedx) ;
		float dau2eta = dau2.PseudoRapidity();
		float dau2phi = dau2.Phi() ; 


		/////////////////-------- jpsi -----/////////////////
		float pt = hlt_diep->ePair[i].pt ;
		float px = cos(hlt_diep->ePair[i].psi)*hlt_diep->ePair[i].pt ;
		float py = sin(hlt_diep->ePair[i].psi)*hlt_diep->ePair[i].pt ;
		float pz = hlt_diep->ePair[i].tanl*hlt_diep->ePair[i].pt ;
		float m = hlt_diep->ePair[i].invariantMass ;

		TVector3 tt(px,py,pz) ;
		float p = tt.Mag() ;
		float eta = tt.PseudoRapidity() ;

		float array[54] = {runID, eventID, daqBits, lmvertX, lmvertY, globMult, primMult, vertX, vertY, vertZ, p, pt, eta, m, dau1chi2XY, dau1chi2Z, dau1DcaXY, dau1DcaZ, dau1q, dau1p, dau1phi, dau1tanl, dau1eta, dau1nHits, dau1ndedx, dau1dedx, nSigma1, dau1TowerEnergy, dau1phidiff, dau1Zedge, dau1beta, dau1tof, dau1localY, dau1localZ, dau2chi2XY, dau2chi2Z, dau2DcaXY, dau2DcaZ, dau2q, dau2p, dau2phi, dau2tanl, dau2eta, dau2nHits, dau2ndedx, dau2dedx, nSigma2, dau2TowerEnergy, dau2phidiff, dau2Zedge, dau2beta, dau2tof, dau2localY, dau2localZ} ;

		jpsi->Fill(array) ;
	}

};


///////// MAIN FUNCTION
int main(int argc, char *argv[]) 
{
	StHltDisplay hltdis;

	//begin  added by Maxim Naglis
	//records the current directory
	getcwd(hltdis.Currentdir,256);
	char tmp[256];
	getcwd(tmp,256);
	strcpy(hltdis.Currentdir,tmp);

	struct stat64 st;
	sprintf(hltdis.Destindir,"%s/current",hltdis.Currentdir);         
	if(stat64(hltdis.Destindir,&st) == 0)
	{
		printf("%s exist.\n",hltdis.Destindir);
	}
	else
	{       
		printf("%s does not exist. Create.\n",hltdis.Destindir);
		if (mkdir(hltdis.Destindir, S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH) != 0)
			perror("mkdir() error");
	}
	//end  added by Maxim Naglis

	hltdis.Main(argc, argv);
}

//-----------write histograms------------//
void StHltDisplay::WriteHistogram(char *outFile, char *outNtuple)
{	
	char histfile[256];
	sprintf(histfile,"%s",outFile);
	TFile file1(histfile,"RECREATE");
	for(int i=0;i<55;i++){
		getPlotByIndex(i)->getHisto(0)->histo->Write();
	}
	hHFM_dEdx->Write();
	hDiElectronInvMassTpxEmcBG->Write();
	hDiElectronInvMassFullRangeBG->Write();
	hDiElectronInvMassCutBG->Write();
	//	hMatchannel3D->Write();
	file1.Close();

	char ntuplefile[256];
	sprintf(ntuplefile,"%s",outNtuple);
	TFile file2(ntuplefile,"RECREATE");
	jpsi->Write();
	heavyfrag->Write();
	file2.Close();
}

//-----------nSigma calculation----------///
double StHltDisplay::getDedxMeanTriton(double p)
{
	int n = 11900 ;
	double Min = 0.1 ;
	double Max = 12. ;
	if(p<0.1) return 0;
	if(p>12.) return dedxMean_Tri[n];

	double Bin = (Max - Min)/n ;
	int lndex = (int)((p-Min)/Bin) ;
	double dp = p - Min - lndex*Bin ;

	return (1. - dp/Bin)*dedxMean_Tri[lndex] + dp/Bin*dedxMean_Tri[lndex+1];
}

//-----------input dEdx theoretical value----------///
void StHltDisplay::inputTritondEdx()
{
	string tem;
	ifstream ifs(dEdxMeanFiles[6]);
	getline(ifs,tem);

	for(int j=0;j<11901;j++)
	{
		ifs>>tem>>dedxMean_Tri[j];
		dedxMean_Tri[j] *= 1.e-06;

	}
}

double StHltDisplay::getDedxMeanHe3(double p)
{
	int n = 11900 ;
	double Min = 0.1 ;
	double Max = 12. ;
	if(p<0.1) return 0;
	if(p>12.) return dedxMean_He3[n];

	double Bin = (Max - Min)/n ;
	int lndex = (int)((p-Min)/Bin) ;
	double dp = p - Min - lndex*Bin ;

	return (1. - dp/Bin)*dedxMean_He3[lndex] + dp/Bin*dedxMean_He3[lndex+1];
}

void StHltDisplay::inputHe3dEdx()
{
	string tem;
	ifstream ifs(dEdxMeanFiles[5]);
	getline(ifs,tem);

	for(int j=0;j<11901;j++)
	{
		ifs>>tem>>dedxMean_He3[j];
		dedxMean_He3[j] *= 1.e-06;

	}
}

double StHltDisplay::getDedxMeanHe4(double p)
{
	int n = 11900 ;
	double Min = 0.1 ;
	double Max = 12. ;
	if(p<0.1) return 0;
	if(p>12.) return dedxMean_He4[n];

	double Bin = (Max - Min)/n ;
	int lndex = (int)((p-Min)/Bin) ;
	double dp = p - Min - lndex*Bin ;

	return (1. - dp/Bin)*dedxMean_He4[lndex] + dp/Bin*dedxMean_He4[lndex+1];
}

void StHltDisplay::inputHe4dEdx()
{
	string tem;
	ifstream ifs(dEdxMeanFiles[7]);
	getline(ifs,tem);

	for(int j=0;j<11901;j++)
	{
		ifs>>tem>>dedxMean_He4[j];
		dedxMean_He4[j] *= 1.e-06;
	}
}

double StHltDisplay::getDedxMeanElectron(double p)
{
	int n = 11900 ;
	double Min = 0.1 ;
	double Max = 12. ;
	if(p<0.1) return 0;
	if(p>12.) return dedxMean_Elec[n];
	double Bin = (Max - Min)/n ;
	int lndex = (int)((p-Min)/Bin) ;
	double dp = p - Min - lndex*Bin ;

	return (1. - dp/Bin)*dedxMean_Elec[lndex] + dp/Bin*dedxMean_Elec[lndex+1];
}

void StHltDisplay::inputElectrondEdx()
{
	ifstream ifs(dEdxMeanFiles[3]);
	string tem;
	getline(ifs,tem);

	for(int j=0;j<11901;j++)
	{
		ifs>>tem>>dedxMean_Elec[j];
		dedxMean_Elec[j] *= 1.e-06;
	}
}

//-----------initialization----------------//
void StHltDisplay::initialize(int argc, char *argv[])
{	
	gStyle->SetPalette(1);
	gStyle->SetOptLogz(1);

    cout<< "Initialization ..." <<endl; 

	index = 0;
	runnumber = 0;
	iBin = 0;
    cout<< "iBin in initial function ..." << iBin <<endl; 
	twopi = 6.2831854;
	pi = 3.1415927;
	A = 0.3736 ;
	BeamX = 0;
	BeamY = 0;
	innerGainPara = -999;
	outerGainPara = -999;

	sprintf(dEdxMeanFiles[0],"%s/StRoot/StJevpPool/StJevpBuilders/dedx_mean_Pion",Currentdir);
	sprintf(dEdxMeanFiles[1],"%s/StRoot/StJevpPool/StJevpBuilders/dedx_mean_Kaon",Currentdir);
	sprintf(dEdxMeanFiles[2],"%s/StRoot/StJevpPool/StJevpBuilders/dedx_mean_Proton",Currentdir);
	sprintf(dEdxMeanFiles[3],"%s/StRoot/StJevpPool/StJevpBuilders/dedx_mean_Electron",Currentdir);
	sprintf(dEdxMeanFiles[4],"%s/StRoot/StJevpPool/StJevpBuilders/dedx_mean_Deuteron",Currentdir);
	sprintf(dEdxMeanFiles[5],"%s/StRoot/StJevpPool/StJevpBuilders/dedx_mean_He3",Currentdir);
	sprintf(dEdxMeanFiles[6],"%s/StRoot/StJevpPool/StJevpBuilders/dedx_mean_Triton",Currentdir);
	sprintf(dEdxMeanFiles[7],"%s/StRoot/StJevpPool/StJevpBuilders/dedx_mean_He4",Currentdir);

	for(int i=0;i<60;i++) { HLTPlots[i] = new JevpPlot();}


	jpsi = new TNtuple("jpsi", "jpsi", "runID:eventID:daqBits:beamPosX:beamPosY:globMult:primMult:vtxX:vtxY:vtxZ:P:Pt:Eta:M:dghtr1Chi2XY:dghtr1Chi2Z:dghtr1GlobalDcaXY:dghtr1GlobalDcaZ:dghtr1q:dghtr1p:dghtr1phi:dghtr1tanl:dghtr1eta:dghtr1nHits:dghtr1ndedx:dghtr1dedx:dghtr1nSigma:dghtr1Energy:dghtr1PhiDiff:dghtr1Zedge:dghtr1beta:dghtr1tof:dghtr1LocalY:dghtr1LocalZ:dghtr2Chi2XY:dghtr2Chi2Z:dghtr2GlobalDcaXY:dghtr2GlobalDcaZ:dghtr2q:dghtr2p:dghtr2phi:dghtr2tanl:dghtr2eta:dghtr2nHits:dghtr2ndedx:dghtr2dedx:dghtr2nSigma:dghtr2Energy:dghtr2PhiDiff:dghtr2Zedge:dghtr2beta:dghtr2tof:dghtr2LocalY:dghtr2LocalZ");

	heavyfrag = new TNtuple("heavyfrag","heavyfrag","runID:eventID:daqBits:beamPosX:beamPosY:globMult:primMult:vtxX:vtxY:vtxZ:Chi2XY:Chi2Z:GlobalDcaXY:GlobalDcaZ:globP:globPt:primP:primPt:q:tanl:psi:eta:nDedx:nHits:dedx:Energy:phidiff:zedge:beta:tof:localZ:localY:Msquare:dedxRatioTri:dedxRatioHe3:dedxRatioHe4:nSigmaTriton:nSigmaHe3:nSigmaHe4") ;

	//	hMatchannel3D = new TH3D("matchannel3D","matchannel3D",200,0,200,200,0,200,120,0,120) ;

	inputElectrondEdx() ;
	inputTritondEdx() ;
	inputHe3dEdx() ;
	inputHe4dEdx() ;

	//-------- hlt plots every run---------//
	HLTPlots[index]->logy = 1 ; //0
	hEvtsAccpt = new TH1I("EvtsAccpt","EvtsAccpt",10,0.,10);
	PlotHisto *ph = new PlotHisto();
	ph->histo = hEvtsAccpt;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	//---------Tracks--------//
	index++ ;    //1
	hnhits = new TH1I("nHits","nHits",50,0,50);
	ph = new PlotHisto();
	ph->histo = hnhits;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ; //2
	hnDedx = new TH1I("nDedx","nDedx",50,0,50);
	ph = new PlotHisto();
	ph->histo = hnDedx;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ; //3
	hDcaXy = new TH1D("DcaXy","DcaXy",120,-6.,6.);
	ph = new PlotHisto();
	ph->histo = hDcaXy;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ; //4
	hDcaZ = new TH1D("DcaZ","DcaZ",120,-6.,6.);
	ph = new PlotHisto();
	ph->histo = hDcaZ;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ; //5
	hLn_dEdx = new TH1D("Ln_dEdx","Ln_dEdx",500,-13.3,-12.3);  //
	ph = new PlotHisto();
	ph->histo = hLn_dEdx;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);


	//----------Glob Tracks---------//
	index++ ;//6
	HLTPlots[index]->logy=1;
	hGlob_Pt = new TH1D("Glob_Pt","Glob_Pt",150,0.,15.);
	ph = new PlotHisto();
	ph->histo = hGlob_Pt;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//7
	hGlob_Phi = new TH1D("Glob_Phi","Glob_Phi",360,0.,twopi);
	ph = new PlotHisto();
	ph->histo = hGlob_Phi;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;  //8
	hGlob_Eta = new TH1D("Glob_Eta","Glob_Eta",120,-3,3);
	ph = new PlotHisto();
	ph->histo = hGlob_Eta;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//9
	HLTPlots[index]->optstat = 0;
	HLTPlots[index]->setDrawOpts("colz");
	hGlob_dEdx = new TH2F("Glob_dEdx","Glob_dEdx",200,-5,5,100,0,1.e-5);
	ph = new PlotHisto();
	ph->histo = hGlob_dEdx;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	//------Prim Tracks-----//

	index++ ;//10
	HLTPlots[index]->logy=1;
	hPrim_Pt = new TH1D("Prim_Pt","Prim_Pt",150,0.,15.);
	ph = new PlotHisto();
	ph->histo = hPrim_Pt;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//11
	hPrim_Phi = new TH1D("Prim_Phi","Prim_Phi",360,0.,twopi);
	ph = new PlotHisto();
	ph->histo = hPrim_Phi;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//12
	hPrim_Eta = new TH1D("Prim_Eta","Prim_Eta",120,-3,3);
	ph = new PlotHisto();
	ph->histo = hPrim_Eta;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//13
	HLTPlots[index]->optstat = 0;
	HLTPlots[index]->setDrawOpts("colz");
	hPrim_dEdx = new TH2F("Prim_dEdx","Prim_dEdx",200,-5,5,100,0,1.e-5);
	ph = new PlotHisto();
	ph->histo = hPrim_dEdx;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	//------Event------//
	index++ ;//14
	hVertexX = new TH1D("VertexX","VertexX",200,-2.,2.);
	ph = new PlotHisto();
	ph->histo = hVertexX;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//15
	hVertexY = new TH1D("VertexY","VertexY",200,-2.,2.);
	ph = new PlotHisto();
	ph->histo = hVertexY;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//16
	hVertexZ = new TH1D("VertexZ","VertexZ",1000,-200.,200.);
	ph = new PlotHisto();
	ph->histo = hVertexZ;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//17
	hLm_VertexX = new TH1D("Lm_VertexX","Lm_VertexX",200,-2.,2.);
	ph = new PlotHisto();
	ph->histo = hLm_VertexX;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//18
	hLm_VertexY = new TH1D("Lm_VertexY","Lm_VertexY",200,-2.,2.);
	ph = new PlotHisto();
	ph->histo = hLm_VertexY;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//19
	hLm_VertexZ = new TH1D("Lm_VertexZ","Lm_VertexZ",1000,-200.,200.);
	ph = new PlotHisto();
	ph->histo = hLm_VertexZ;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//20
	HLTPlots[index]->logy=1;
	hglobalMult = new TH1I("globalMult", "globalMult",4200,0,4200);
	ph = new PlotHisto();
	ph->histo = hglobalMult;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//21
	HLTPlots[index]->logy=1;
	hprimaryMult = new TH1I("primaryMult", "primaryMult",2200,0,2200);
	ph = new PlotHisto();
	ph->histo = hprimaryMult;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	//--------Emc---------//
	index++ ;//22
	hMatchPhi_Diff = new TH1D("Emc_matchPhiDiff","Emc_matchPhiDiff",50,0.,0.1);
	ph = new PlotHisto();
	ph->histo = hMatchPhi_Diff;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//23
	hTowerEnergy = new TH1D("Emc_towerEnergy","Emc_towerEnergy",200,0.,20.);
	ph = new PlotHisto();
	ph->histo = hTowerEnergy;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//24
	hTowerDaqId = new TH1I("Emc_towerDaqId","Emc_towerDaqId",5000,0.,5000.);
	ph = new PlotHisto();
	ph->histo = hTowerDaqId;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//25
	hTowerSoftId = new TH1I("Emc_towerSoftId","Emc_towerSoftId",5000,0.,5000.);
	ph = new PlotHisto();
	ph->histo = hTowerSoftId;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//26
	hzEdge = new TH1D("Emc_zEdge","Emc_zEdge",100,0.,5.);
	ph = new PlotHisto();
	ph->histo = hzEdge;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//27
	HLTPlots[index]->optstat = 0;
	HLTPlots[index]->setDrawOpts("colz");
	hTowerEtaPhi = new TH2F("Emc_towerEtaPhi","Emc_towerEtaPhi",120,-pi,pi,40,-1,1);
	ph = new PlotHisto();
	ph->histo = hTowerEtaPhi;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	//-------jspi invariant mass----//
	index++ ;//28
	hDiElectronInvMassTpxEmc = new TH1D("DiElectronInvMassTpxEmc","DiElectronInvMassTpxEmc",120,1.,13.);
	ph = new PlotHisto();
	ph->histo = hDiElectronInvMassTpxEmc;
	HLTPlots[index]->addHisto(ph);

	hDiElectronInvMassTpxEmcBG = new TH1D("DiElectronInvMassTpxEmcBG","DiElectronInvMassTpxEmcBG",120,1.,13.);
	ph = new PlotHisto();
	ph->histo = hDiElectronInvMassTpxEmcBG;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//29
	hDiElectronInvMassFullRange = new TH1D("DiElectronInvMassFullRange ","DiElectronInvMassFullRange",130,0.,13.);
	ph = new PlotHisto();
	ph->histo = hDiElectronInvMassFullRange ;
	HLTPlots[index]->addHisto(ph);

	hDiElectronInvMassFullRangeBG = new TH1D ("DiElectronInvMassFullRangeBG","DiElectronInvMassFullRangeBG",130,0.,13.);
	ph = new PlotHisto();
	ph->histo = hDiElectronInvMassFullRangeBG;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//30
	hDiElectronInvMassCut = new TH1D("DiElectronInvMassCut","DiElectronInvMassCut",120,1.,13.);
	ph = new PlotHisto();
	ph->histo = hDiElectronInvMassCut;
	HLTPlots[index]->addHisto(ph);

	hDiElectronInvMassCutBG = new TH1D("DiElectronInvMassCutBG","DiElectronInvMassCutBG",120,1.,13.);
	ph = new PlotHisto();
	ph->histo = hDiElectronInvMassCutBG;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	//--------daug e1------//
	index++ ;//31
	HLTPlots[index]->optstat = 0;
	HLTPlots[index]->setDrawOpts("colz");
	hdEdx_P1 = new TH2F("dEdx_P1","dEdx_P1",200,0.,10.,55,0.,5.5e-06);
	ph = new PlotHisto();
	ph->histo = hdEdx_P1;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//32
	hDaughter1P_TowerEnergy = new TH1D("Daughter1P_TowerEnergy","Daughter1P_TowerEnergy",100,0,5);
	ph = new PlotHisto();
	ph->histo = hDaughter1P_TowerEnergy;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ; //33
	hDaughter1TpxEmcInverseBeta = new TH1D("Daughter1TpxEmcInverseBeta","Daughter1TpxEmcInverseBeta",100,0,5);
	ph = new PlotHisto();
	ph->histo =  hDaughter1TpxEmcInverseBeta;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	//------daug e2-----//
	index++ ;//34
	HLTPlots[index]->optstat = 0;
	HLTPlots[index]->setDrawOpts("colz");
	hdEdx_P2 = new TH2F("dEdx_P2","dEdx_P2",200,0.,10.,55,0.,5.5e-06);
	ph = new PlotHisto();
	ph->histo = hdEdx_P2;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//35
	hDaughter2P_TowerEnergy = new TH1D("Daughter2P_TowerEnergy","Daughter2P_TowerEnergy",100,0,5);
	ph = new PlotHisto();
	ph->histo = hDaughter2P_TowerEnergy;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ; //36
	hDaughter2TpxEmcInverseBeta = new TH1D("Daughter2TpxEmcInverseBeta","Daughter2TpxEmcInverseBeta",100,0,5);
	ph = new PlotHisto();
	ph->histo =  hDaughter2TpxEmcInverseBeta;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//37
	hDiLeptonRapidity = new TH1D("DiLeptonRapidity","DiLeptonRapidity",150,-7.5,7.5);
	ph = new PlotHisto();
	ph->histo = hDiLeptonRapidity;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	//------Tof histograms-------//
	index++ ;//38
	hLocalZ = new TH1D("Tof_LocalZ","Tof_LocalZ",100,-5.,5.0);
	ph = new PlotHisto();
	ph->histo = hLocalZ;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//39
	hLocalY = new TH1D("Tof_LocalY","Tof_LocalY",300,-15.,15.);
	ph = new PlotHisto();
	ph->histo = hLocalY;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//40
	HLTPlots[index]->optstat = 0;
	HLTPlots[index]->setDrawOpts("colz");
	hInverseBeta = new TH2F("Tof_InverseBeta","Tof_InverseBeta",500,0,5,500,0.0,5.);
	ph = new PlotHisto();
	ph->histo = hInverseBeta;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//41
	HLTPlots[index]->optstat = 0;
	HLTPlots[index]->setDrawOpts("colz");
	hMatchId_fiberId = new TH2F("Tof_matchId_fireId","Tof_matchId_fireId",200,0,200,200,0,200); 
	ph = new PlotHisto();
	ph->histo = hMatchId_fiberId;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//42
	HLTPlots[index]->optstat = 0;
	HLTPlots[index]->setDrawOpts("colz");
	hTrayID_TrgTime = new TH2F("Tof_TrayID_TrgTime","Tof_TrayID_TrgTime",124,0.,124,400,2700,3100);
	ph = new PlotHisto();
	ph->histo = hTrayID_TrgTime;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//43
	hchannelID = new TH1D("Tof_channelID","Tof_channelID",200,0,200);
	ph = new PlotHisto();
	ph->histo = hchannelID;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//44
	HLTPlots[index]->optstat = 0; 
	HLTPlots[index]->setDrawOpts("colz");
	hVzvpd_Vz = new TH2F("Vzvpd_Vz","Vzvpd_Vz",400,-100,100,400,-100,100);
	ph = new PlotHisto();
	ph->histo = hVzvpd_Vz;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//45
	hVzDiff = new TH1D("VzDiff","VzDiff",200,-20,20);
	ph = new PlotHisto();
	ph->histo = hVzDiff;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++ ;//46
	HLTPlots[index]->optstat = 0; 
	HLTPlots[index]->setDrawOpts("colz");
	hdEdx = new TH2F("dEdx","dEdx",500,-5,5,300,0,3.e-5);
	ph = new PlotHisto(); 
	ph->histo = hdEdx;
	HLTPlots[index]->addHisto(ph);

	hHFM_dEdx = new TH2F("HFM_dEdx","HFM_dEdx",500,-5,5,300,0,3.e-5);
	ph = new PlotHisto();
	ph->histo = hHFM_dEdx;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);


	//--------------------------------------------------------------------
	index++;//47
	hUpc = new TH1I("UPC","UPC",10,0.,10);
	ph = new PlotHisto(); 
	ph->histo = hUpc;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++;//48
	hMtd = new TH1I("MTD","MTD",10,0.,10);
	ph = new PlotHisto();
	ph->histo = hMtd;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++;  //49
	hNpeHt_25_NoZdc = new TH1I("NPEHt_25_NOZDC","NPEHT_25_NOZDC",10,0.,10);
	ph = new PlotHisto();
	ph->histo = hNpeHt_25_NoZdc;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++;//50
	hVpdMb = new TH1I("VPDMB","VPDMB",10,0.,10);
	ph = new PlotHisto();
	ph->histo = hVpdMb;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++;//51
	hCentral = new TH1I("Central","Central",10,0.,10);
	ph = new PlotHisto();
	ph->histo = hCentral;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++;//52
	hNpeHt_25 = new TH1I("NPEHt_25","NPEHt_25",10,0.,10);
	ph = new PlotHisto();
	ph->histo = hNpeHt_25;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++;//53
	hNpe = new TH1I("NPE11_NPE15_NPE18","NPE11_NPE15_NPE18",10,0.,10);
	ph = new PlotHisto();
	ph->histo = hNpe;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++;//54
	hAtomcule = new TH1I("Atomcule","Atomcule",10,0.,10);
	ph = new PlotHisto();
	ph->histo = hAtomcule;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);


	//----------------BeamX BeamY Position-----------------//
	index++;//55
	hBeamX = new TH1D("BeamX","BeamX",105,0.,105);
	ph = new PlotHisto();
	ph->histo = hBeamX;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++;//56
	hBeamY = new TH1D("BeamY","BeamY",105,0.,105);
	ph = new PlotHisto();
	ph->histo = hBeamY;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++; //57
	hInnerGain = new TH1D("innerGain","innerGain",105,0.,105);
	ph = new PlotHisto();
	ph->histo = hInnerGain;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++; //58
	hOuterGain = new TH1D("outerGain","outerGain",105,0.,105);
	ph = new PlotHisto();
	ph->histo = hOuterGain;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	index++;//59
	hMeanDcaXy = new TH1D("meanDcaXy","meanDcaXy",105,0.,105);
	ph = new PlotHisto();
	ph->histo = hMeanDcaXy;
	HLTPlots[index]->addHisto(ph);
	addPlot(HLTPlots[index]);

	/////////////////////////////////////////
	hGlob_Phi->SetMinimum(0.) ;
	hPrim_Phi->SetMinimum(0.) ;
	hBeamX->SetLineColor(kBlue);
	hBeamY->SetLineColor(kRed);
	hInnerGain->SetLineColor(kBlue); 
	hOuterGain->SetLineColor(kRed);
	hHFM_dEdx->SetMarkerStyle(30);
	hHFM_dEdx->SetMarkerSize(0.9);
	hHFM_dEdx->SetMarkerColor(2);


	//---------------set X Y titles---------------//

	hnhits->GetXaxis()->SetTitle("nHits") ;
	hnDedx->GetXaxis()->SetTitle("ndedx") ;
	hDcaXy->GetXaxis()->SetTitle("DcaXY") ;
	hDcaZ->GetXaxis()->SetTitle("DcaZ") ;
	hLn_dEdx->GetXaxis()->SetTitle("log(dEdx) GeV/cm)") ;
	hGlob_Pt->GetXaxis()->SetTitle("Pt") ;
	hGlob_Phi->GetXaxis()->SetTitle("#phi") ;
	hGlob_Eta->GetXaxis()->SetTitle("#eta") ;
	hGlob_dEdx->GetXaxis()->SetTitle("Global Momentum") ;
	hGlob_dEdx->GetYaxis()->SetTitle("dEdx in GeV/cm") ;
	hPrim_Pt->GetXaxis()->SetTitle("Pt") ;
	hPrim_Phi->GetXaxis()->SetTitle("#phi") ;
	hPrim_Eta->GetXaxis()->SetTitle("#eta") ;
	hPrim_dEdx->GetXaxis()->SetTitle("Primary Mommentum") ;
	hPrim_dEdx->GetYaxis()->SetTitle("dEdx in GeV/cm") ;

	hVertexX->GetXaxis()->SetTitle("VertexX") ;
	hVertexY->GetXaxis()->SetTitle("VertexY") ;
	hVertexZ->GetXaxis()->SetTitle("VertexZ") ;
	hLm_VertexX->GetXaxis()->SetTitle("LmVertexX") ;
	hLm_VertexY->GetXaxis()->SetTitle("LmVertexY") ;
	hLm_VertexZ->GetXaxis()->SetTitle("LmVertexZ") ;
	hglobalMult->GetXaxis()->SetTitle("Multiplicity") ;
	hprimaryMult->GetXaxis()->SetTitle("Multiplicity") ;

	hMatchPhi_Diff->GetXaxis()->SetTitle("matchPhiDiff") ;
	hTowerEnergy->GetXaxis()->SetTitle("TowerEnergy") ;
	hTowerDaqId->GetXaxis()->SetTitle("TowerDaqId") ;
	hTowerSoftId->GetXaxis()->SetTitle("TowerSoftId") ;
	hzEdge->GetXaxis()->SetTitle("zEdge") ;
	hTowerEtaPhi->GetXaxis()->SetTitle("#phi") ;
	hTowerEtaPhi->GetYaxis()->SetTitle("#eta") ;

	hDiElectronInvMassTpxEmc->GetXaxis()->SetTitle("M_{inv}(ee) GeV/c^{2}") ;
	hDiElectronInvMassTpxEmcBG->GetXaxis()->SetTitle("M_{inv}(ee) GeV/c^{2}") ;
	hDiElectronInvMassFullRange->GetXaxis()->SetTitle("M_{inv}(ee) GeV/c^{2}") ;
	hDiElectronInvMassFullRangeBG->GetXaxis()->SetTitle("M_{inv}(ee) GeV/c^{2}") ;
	hDiElectronInvMassCut->GetXaxis()->SetTitle("M_{inv}(ee) GeV/c^{2}") ;
	hDiElectronInvMassCutBG->GetXaxis()->SetTitle("M_{inv}(ee) GeV/c^{2}") ;

	hdEdx_P1->GetXaxis()->SetTitle("Daughter1 Momentum") ;
	hdEdx_P1->GetYaxis()->SetTitle("dEdx in GeV/cm") ;
	hDaughter1P_TowerEnergy->GetXaxis()->SetTitle("TowerEnergy/P") ;
	hDaughter1TpxEmcInverseBeta->GetXaxis()->SetTitle("1/#beta") ;
	hdEdx_P2->GetXaxis()->SetTitle("Daughter2 Momentum") ;
	hdEdx_P2->GetYaxis()->SetTitle("dEdx in GeV/cm") ;
	hDaughter2P_TowerEnergy->GetXaxis()->SetTitle("TowerEnergy/P") ;
	hDaughter2TpxEmcInverseBeta->GetXaxis()->SetTitle("1/#beta") ;
	hDiLeptonRapidity->GetXaxis()->SetTitle("Rapidity") ;

	hLocalZ->GetXaxis()->SetTitle("LocalZ") ;
	hLocalY->GetXaxis()->SetTitle("LocalY") ;
	hInverseBeta->GetXaxis()->SetTitle("Momentum") ;
	hInverseBeta->GetYaxis()->SetTitle("1/#beta") ;
	hMatchId_fiberId->GetXaxis()->SetTitle("matchId") ;
	hMatchId_fiberId->GetYaxis()->SetTitle("fiberId") ;
	hTrayID_TrgTime->GetXaxis()->SetTitle("TrayId") ;
	hTrayID_TrgTime->GetYaxis()->SetTitle("TriggerTime") ;
	hchannelID->GetXaxis()->SetTitle("ChannelId") ;
	hVzvpd_Vz->GetXaxis()->SetTitle("pvpd VertexZ") ;
	hVzvpd_Vz->GetYaxis()->SetTitle("LmVertexZ") ;
	hVzDiff->GetXaxis()->SetTitle("Vzvpd - LmVertexZ") ;
	hdEdx->GetXaxis()->SetTitle("Primary Momentum") ;
	hdEdx->GetYaxis()->SetTitle("dEdx in GeV/cm") ;
	hHFM_dEdx->GetXaxis()->SetTitle("Primary Momentum") ;
	hHFM_dEdx->GetYaxis()->SetTitle("dEdx in GeV/cm") ;
	hBeamX->GetXaxis()->SetTitle("run number");
	hBeamY->GetXaxis()->SetTitle("run number");
	hInnerGain->GetXaxis()->SetTitle("run number");
	hOuterGain->GetXaxis()->SetTitle("run number");
	hMeanDcaXy->GetXaxis()->SetTitle("run number");


	/////////////////////////////////////////////////////////
	//-------------Bin lebal-------------//
	hEvtsAccpt->GetXaxis()->SetBinLabel(1,"All");
	//        hEvtsAccpt->GetXaxis()->SetBinLabel(2,"BES good Events");
	//        hEvtsAccpt->GetXaxis()->SetBinLabel(3,"HLT zerobias");
	hEvtsAccpt->GetXaxis()->SetBinLabel(2,"High Pt");
	hEvtsAccpt->GetXaxis()->SetBinLabel(3,"J/psi");
	hEvtsAccpt->GetXaxis()->SetBinLabel(4,"Heavy Fragment");
	hEvtsAccpt->GetXaxis()->SetBinLabel(5,"HLT zerobias");

	hUpc->GetXaxis()->SetBinLabel(1,"High Pt");
	hUpc->GetXaxis()->SetBinLabel(2,"J/psi");
	hUpc->GetXaxis()->SetBinLabel(3,"Heavy Fragment");
	hMtd->GetXaxis()->SetBinLabel(1,"High Pt");
	hMtd->GetXaxis()->SetBinLabel(2,"J/psi");
	hMtd->GetXaxis()->SetBinLabel(3,"Heavy Fragment");
	hNpeHt_25_NoZdc->GetXaxis()->SetBinLabel(1,"High Pt");
	hNpeHt_25_NoZdc->GetXaxis()->SetBinLabel(2,"J/psi");
	hNpeHt_25_NoZdc->GetXaxis()->SetBinLabel(3,"Heavy Fragment");
	hVpdMb->GetXaxis()->SetBinLabel(1,"High Pt");
	hVpdMb->GetXaxis()->SetBinLabel(2,"J/psi");
	hVpdMb->GetXaxis()->SetBinLabel(3,"Heavy Fragment");
	hCentral->GetXaxis()->SetBinLabel(1,"High Pt");
	hCentral->GetXaxis()->SetBinLabel(2,"J/psi");
	hCentral->GetXaxis()->SetBinLabel(3,"Heavy Fragment");
	hNpeHt_25->GetXaxis()->SetBinLabel(1,"High Pt");
	hNpeHt_25->GetXaxis()->SetBinLabel(2,"J/psi");
	hNpeHt_25->GetXaxis()->SetBinLabel(3,"Heavy Fragment");
	hNpe->GetXaxis()->SetBinLabel(1,"High Pt NPE11");
	hNpe->GetXaxis()->SetBinLabel(2,"J/psi NPE11");
	hNpe->GetXaxis()->SetBinLabel(3,"Heavy Fragment NPE11");
	hNpe->GetXaxis()->SetBinLabel(4,"High Pt NPE15");
	hNpe->GetXaxis()->SetBinLabel(5,"J/psi NPE15");
	hNpe->GetXaxis()->SetBinLabel(6,"Heavy Fragment NPE15");
	hNpe->GetXaxis()->SetBinLabel(7,"High Pt NPE18");
	hNpe->GetXaxis()->SetBinLabel(8,"J/psi NPE18");
	hNpe->GetXaxis()->SetBinLabel(9,"Heavy Fragment NPE18");
	hAtomcule->GetXaxis()->SetBinLabel(1,"High Pt");
	hAtomcule->GetXaxis()->SetBinLabel(2,"J/psi");
	hAtomcule->GetXaxis()->SetBinLabel(3,"Heavy Fragment");

}
