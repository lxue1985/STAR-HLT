#include "StJevpPool/StJevpPlot/JevpPlotSet.h"
#include "DAQ_READER/daq_dta.h"
#include "DAQ_L3/daq_l3.h"

#include "TNtuple.h"
#include <TH1I.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3D.h>


class StHltDisplay : public JevpPlotSet {

	public:
		
		void WriteHistogram(char *outFile, char *ourNtuple); //kinds of histograms of current run
		void inputElectrondEdx() ;
		double getDedxMeanElectron(double p);
		void inputTritondEdx();
		double getDedxMeanTriton(double p);
		void inputHe3dEdx();
		double getDedxMeanHe3(double p);
		void inputHe4dEdx() ;
		double getDedxMeanHe4(double p);

		void initialize(int argc, char *argv[]) ;
		
		void event(daqReader *rdr) ;   ///< event by event fill QA plots
		void ntuple(daqReader *rdr) ;  ///< event by event fill ntuple files
		int selectEvent(daqReader *rdr) { return 1; };
		int selectRun(daqReader *rdr) { return 1; };

		void startrun(daqReader *rdr);
		void stoprun(daqReader *rdr);
		
		char Currentrun[256];
		char CurrentNtuple[256];
		char Currentdir[256];
		char Destindir[256];
		char dEdxMeanFiles[8][256];

		private:
		
		JevpPlot *HLTPlots[60];

		TNtuple* jpsi ;
		TNtuple* heavyfrag ;

		int index;
		int primaryTracks ;
		int runnumber;
		int iBin;
		double innerGainPara;
		double outerGainPara;
		double BeamX;
		double BeamY;
		double twopi;
		double pi;
		double A;
		double dedxMean_Elec[11901];
		double dedxMean_Tri[11901];
		double dedxMean_He3[11901];
		double dedxMean_He4[11901];


		TH1I *hEvtsAccpt;
		
		///< track level QA
		TH1I *hnhits;
		TH1I *hnDedx; 
		TH1D *hDcaXy;
		TH1D *hDcaZ ;
		TH1D *hLn_dEdx;
		TH1D *hGlob_Pt;
		TH1D *hGlob_Phi;
		TH1D *hGlob_Eta;
		TH2F *hGlob_dEdx;
		TH1D *hPrim_Pt;
		TH1D *hPrim_Phi;
		TH1D *hPrim_Eta;
		TH2F *hPrim_dEdx;
		
		///< event level QA
		TH1D *hVertexX; 
		TH1D *hVertexY;
		TH1D *hVertexZ;
		TH1D *hLm_VertexX;
		TH1D *hLm_VertexY;
		TH1D *hLm_VertexZ;
		TH1I *hglobalMult;
		TH1I *hprimaryMult;
		
		///< emc QA
		TH1D *hMatchPhi_Diff;
		TH1D *hTowerEnergy ;
		TH1I *hTowerDaqId; 
		TH1I *hTowerSoftId;
		TH1D *hzEdge;
		TH2F *hTowerEtaPhi;

		///< di-electron QA
		TH1D *hDiElectronInvMassTpxEmc;
		TH1D *hDiElectronInvMassTpxEmcBG;
		TH1D *hDiElectronInvMassFullRange;
		TH1D *hDiElectronInvMassFullRangeBG;
		TH1D *hDiElectronInvMassCut;
		TH1D *hDiElectronInvMassCutBG;
		TH2F *hdEdx_P1; 
		TH1D *hDaughter1P_TowerEnergy;
		TH1D *hDaughter1TpxEmcInverseBeta;
		TH2F *hdEdx_P2;
		TH1D *hDaughter2P_TowerEnergy;
		TH1D *hDaughter2TpxEmcInverseBeta; 
		TH1D *hDiLeptonRapidity;

		///< tof qualities 
		TH1D *hLocalZ;
		TH1D *hLocalY;
		TH2F *hInverseBeta;
		TH2F *hMatchId_fiberId;
		TH2F *hTrayID_TrgTime;
		TH1D *hchannelID;
		TH2F *hVzvpd_Vz ;
		TH1D *hVzDiff;
//		TH3D *hMatchannel3D ;
		
		///< heavy fragments
		TH2F *hdEdx;
		TH2F *hHFM_dEdx;
		
		///< different event tpye
		TH1I *hUpc ;
		TH1I *hMtd ;
		TH1I *hNpeHt_25_NoZdc;
		TH1I *hVpdMb; 
		TH1I *hCentral;
		TH1I *hNpeHt_25;
		TH1I *hNpe;
		TH1I *hAtomcule ;
		
		///< run by run display
		TH1D *hBeamX;
		TH1D *hBeamY; 
		TH1D *hInnerGain;
		TH1D *hOuterGain;
		TH1D *hMeanDcaXy;

};
