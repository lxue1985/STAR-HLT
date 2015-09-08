#ifndef _HLT_GL3_H_
#define _HLT_GL3_H_


#include <sys/types.h>
#include <string.h>
#include <getopt.h>
#include <stdlib.h>

#include <rtsLog.h>
#include <DAQ_READER/daqReader.h>
#include <DAQ_READER/daq_dta.h>
#include <DAQ_TPX/daq_tpx.h>
#include <DAQ_TRG/daq_trg.h>
#include <DAQ_TOF/daq_tof.h>
#include <DAQ_SC/daq_sc.h>


#include <GL3/hlt_entity.h>
#include <L3_SUPPORT/l3_support.h>


#include "online_tracking_collector.h"
#include "HLTFormats.h"

//#define TIMING_HLT
#ifdef TIMING_HLT
#include <time.h>
#endif

class hlt_gl3 : public hlt_entity {
public:

#ifdef TIMING_HLT
  double timeNow, timeOld;
  double timeResetEvent, timeReadEmc, timeReadTof, timeReadTpc, timeFinalize, timeDecide, timeFill;
  int nTracks, nEvents;
#endif

	hlt_gl3() {
		name = "gl3" ;
		rts_id = L3_ID ;
	}

	~hlt_gl3() {;} ;


	online_tracking_collector *gl3 ;

	char parametersDirectory[64] ;
	char beamlineFile[128] ;
	char GainParameters[128] ;

	u_int algorithms_requested ;

	int run_config() {

		if(gl3) delete(gl3) ;


		LOG(INFO,"  directory %s, field %f",hlt_config.config_directory, hlt_config.l3_sc.mag_field) ;



		sprintf(parametersDirectory, "%s/", hlt_config.config_directory);
		sprintf(beamlineFile, "%sbeamline",parametersDirectory);
		sprintf(GainParameters, "%sGainParameters",parametersDirectory);

		gl3 = new online_tracking_collector(hlt_config.l3_sc.mag_field, 
						    szGL3_default_mHits, 
						    szGL3_default_mTracks, 
						    parametersDirectory,  
						    beamlineFile) ;


		gl3->readGainPara(GainParameters);

		// check triggers...
		for(int i=0;i<32;i++) {
			l3_algos_t *alg = &(hlt_config.hlt_algo[i]) ;

			if(alg->rc.lx.id) ;	
			else continue ;	// nixed

			// only the following and _no_ other triggers will be necessary
			// in this run...

			LOG(INFO,"  hlt %4d: \"%s\" [RC %s]: HLT mask 0x%08X [RC 0x%016llX]: params: %d ...; %f ...; %d",
			    alg->id, alg->name, alg->rc.name,
			    alg->hlt_mask, alg->rc.daq_TrgId,
			    alg->rc.lx.userInt[0],
			    alg->rc.lx.userFloat[0],
			    alg->rc.lx.specialProcessing
			    ) ;
			
			gl3->triggerDecider->setTriggersFromRunControl(alg);
		}
			    

		return 0 ;
	}

	int run_start() {
		if(hlt_config.standalone) gl3->triggerDecider->setQA("qa") ;



#ifdef TIMING_HLT
		timeResetEvent = 0.;
		timeReadEmc = 0.;
		timeReadTof = 0.;
		timeReadTpc = 0.;
		timeFinalize = 0.;
		timeDecide = 0.;
		timeFill = 0.;
		nTracks = 0;
		nEvents = 0;
#endif

		return 0 ;
	}

	int run_stop() {
		gl3->triggerDecider->flushQA() ;
		gl3->triggerDecider->closeQA() ;
		gl3->writeBeamline(beamlineFile,run_number);

		gl3->CalidEdx(GainParameters,run_number);

#ifdef TIMING_HLT		
		printf("nEvents:     %i \n", nEvents);
		printf("nTracks:     %i \n", nTracks);
		printf("\t\t time \t\t per evt (ms) \t per track (micro s) \n");
		printf("timeResetEvent:\t %f\t %f\t %f\n", timeResetEvent, timeResetEvent/nEvents*1000., timeResetEvent/nTracks*1000000.);
		printf("timeReadEmc:\t %f\t %f\t %f\n", timeReadEmc, timeReadEmc/nEvents*1000., timeReadEmc/nTracks*1000000.);
		printf("timeReadTof:\t %f\t %f\t %f\n", timeReadTof, timeReadTof/nEvents*1000., timeReadTof/nTracks*1000000.);
		printf("timeReadTpc:\t %f\t %f\t %f\n", timeReadTpc, timeReadTpc/nEvents*1000., timeReadTpc/nTracks*1000000.);
		printf("timeFinalize:\t %f\t %f\t %f\n", timeFinalize, timeFinalize/nEvents*1000., timeFinalize/nTracks*1000000.);
		printf("timeDecide:\t %f\t %f\t %f\n", timeDecide, timeDecide/nEvents*1000., timeDecide/nTracks*1000000.);
		printf("timeFill:\t %f\t %f\t %f\n", timeFill, timeFill/nEvents*1000., timeFill/nTracks*1000000.);
		double totalTime = timeResetEvent+timeReadEmc+timeReadTof+timeReadTpc+timeFinalize+timeDecide+timeFill;
		printf("totalTime:\t %f\t %f\t %f\n", totalTime, totalTime/nEvents*1000., totalTime/nTracks*1000000.);
#ifdef TIMING_GL3
		printf("---------------------\n");
		printf("timeMakeNodes:\t %f\t %f\t %f\n", gl3->timeMakeNodes, gl3->timeMakeNodes/nEvents*1000., gl3->timeMakeNodes/nTracks*1000000.);
		printf("timeMatchEMC:\t %f\t %f\t %f\n", gl3->timeMatchEMC, gl3->timeMatchEMC/nEvents*1000., gl3->timeMatchEMC/nTracks*1000000.);
		printf("timeMatchTOFg:\t %f\t %f\t %f\n", gl3->timeMatchTOFg, gl3->timeMatchTOFg/nEvents*1000., gl3->timeMatchTOFg/nTracks*1000000.);
		printf("timeMakeVertex:\t %f\t %f\t %f\n", gl3->timeMakeVertex, gl3->timeMakeVertex/nEvents*1000., gl3->timeMakeVertex/nTracks*1000000.);
		printf("timeMakePrimaries:\t %f\t %f\t %f\n", gl3->timeMakePrimaries, gl3->timeMakePrimaries/nEvents*1000., gl3->timeMakePrimaries/nTracks*1000000.);
		printf("timeMatchTOFp:\t %f\t %f\t %f\n", gl3->timeMatchTOFp, gl3->timeMatchTOFp/nEvents*1000., gl3->timeMatchTOFp/nTracks*1000000.);
		printf("timeMeanDedx:\t %f\t %f\t %f\n", gl3->timeMeanDedx, gl3->timeMeanDedx/nEvents*1000., gl3->timeMeanDedx/nTracks*1000000.);
#endif
#endif

		return 0 ;
	}

	int event_start(void *hlt_mask_p=0) {
		int pval = PROFILER(0) ;

#ifdef TIMING_HLT
		timeOld = (double)(clock())/CLOCKS_PER_SEC;
#endif
		gl3->resetEvent() ;

#ifdef TIMING_HLT
		timeNow = (double)(clock())/CLOCKS_PER_SEC;
		timeResetEvent += timeNow - timeOld;
#endif

		// this is the mask of requested algorithms!
		if(hlt_mask_p == 0) {
			algorithms_requested = 0xFFFFFFFF ;	// ALL, knowing nothing better
		}
		else {
			algorithms_requested = *((u_int *)hlt_mask_p) ;
		}

		PROFILER(pval) ;	// 27 us

		return 0 ;
	}

	int event_data(int rts_id, int sector, void *buff, int bytes) {
		int ret ;

		LOG(DBG,"%d: %d %d %d",evt_counter,rts_id,sector,bytes) ;
		if(!buff) return -1;

		int pval = PROFILER(0) ;

		switch(rts_id) {
		case TPX_ID :
#ifdef TIMING_HLT
		  timeOld = (double)(clock())/CLOCKS_PER_SEC;
#endif
		  ret = gl3->Tracking_readSectorTracks((char *)buff) ;

#ifdef TIMING_HLT
		  timeNow = (double)(clock())/CLOCKS_PER_SEC;
		  timeReadTpc += timeNow - timeOld;
#endif		
		  PROFILER(pval) ;	// 33 us per sector; 800 us total
		  break ;
		case TRG_ID :
#ifdef TIMING_HLT
		  timeOld = (double)(clock())/CLOCKS_PER_SEC;
#endif
		  gl3->trigger = ((gl3_trg_send_t *)buff)->daq_TrgId ;
		  gl3->eventNumber = ((gl3_trg_send_t *)buff)->eventNumber ;

		//Tonko: this is gone now! -- use raw data only
		//  gl3->emc->readFromGl3Trg((gl3_trg_send_t *)buff);

#ifdef TIMING_HLT
         	  timeNow = (double)(clock())/CLOCKS_PER_SEC;
		  timeReadEmc += timeNow - timeOld;
#endif
		  PROFILER(pval) ;	// 4.7 us
		  break ;
		case BTOW_ID :
#ifdef TIMING_HLT
		  timeOld = (double)(clock())/CLOCKS_PER_SEC;
#endif

		  gl3->emc->readRawBtowers((btow_t *)buff);

#ifdef TIMING_HLT
		  timeNow = (double)(clock())/CLOCKS_PER_SEC;
		  timeReadEmc += timeNow - timeOld;
#endif
		  PROFILER(pval) ;	// 490 us
		  break ;
		case TOF_ID :

#ifdef TIMING_HLT
		  timeOld = (double)(clock())/CLOCKS_PER_SEC;
#endif
		  gl3->tof->readFromTofMachine((char*)buff);

#ifdef TIMING_HLT
		  timeNow = (double)(clock())/CLOCKS_PER_SEC;
		  timeReadTof += timeNow - timeOld;
#endif

		  PROFILER(pval) ;	// 212 us
		  break ;
		default :
			LOG(ERR,"Unknown detector %d",rts_id) ;
			break ;
		}

		return 0 ;
	}

	int event_do(hlt_blob *out, int *blobs) {
		int decision ;
		int bytes ;

		int pval = PROFILER(0) ;

#ifdef TIMING_HLT
		timeOld = (double)(clock())/CLOCKS_PER_SEC;
#endif

		gl3->finalizeReconstruction() ;

#ifdef TIMING_HLT
		timeNow = (double)(clock())/CLOCKS_PER_SEC;
		timeFinalize += timeNow - timeOld;
		timeOld = timeNow;
#endif

		pval = PROFILER(pval) ;	// 108 us

		// USE algoritms_requested here!
		decision = gl3->triggerDecider->decide(gl3->eventNumber) ;

		pval = PROFILER(pval) ;	// 25 us(?)

#ifdef TIMING_HLT
		timeNow = (double)(clock())/CLOCKS_PER_SEC;
		timeDecide += timeNow - timeOld;
		timeOld = timeNow;
#endif
  
		//--------fill global track--------//
		bytes = gl3->FillGlobTracks(sizeof(store_GT), store_GT) ;

		LOG(DBG,"bytes %d, globle tracks %d",bytes,gl3->getNGlobalTracks()) ;
// cout<<"bytes_GT ="<<bytes<<endl;

		out[0].buff = store_GT ;
		out[0].name = "HLT_GT" ;
		out[0].bytes = bytes ;

                //---------fill primary tracks-------//
		bytes = gl3->FillPrimTracks(sizeof(store_PT), store_PT) ;
		LOG(DBG,"bytes %d, primary tracks %d",bytes,gl3->getNPrimaryTracks()) ;
		out[1].buff = store_PT ;
		out[1].name = "HLT_PT" ;
		out[1].bytes = bytes ;

		//--------fill event---------//
		bytes = gl3->FillEvent(sizeof(store_EVE), store_EVE, decision) ;
// cout<<"decision = "<<decision<<endl;
// cout<<"bytes_Eve = "<<bytes<<endl;
		LOG(DBG,"bytes %d",bytes) ;

		out[2].buff = store_EVE ;
		out[2].name = "HLT_EVE" ;
		out[2].bytes = bytes ;

		//----------fill Tof hits------------//
		bytes = gl3->FillTofHits(sizeof(store_TOF), store_TOF) ;
		LOG(DBG,"bytes %d",bytes) ;
//cout<<"bytes_Tof ="<<bytes<<endl;
		out[3].buff = store_TOF ;
		out[3].name = "HLT_TOF" ;
		out[3].bytes = bytes ;

		//----------fill pvpd Hits------------//
		bytes = gl3->FillPvpdHits(sizeof(store_PVPD), store_PVPD) ;
		LOG(DBG,"bytes %d",bytes) ;
//cout<<"bytes_Pvpd ="<<bytes<<endl;
		out[4].buff = store_PVPD ;
		out[4].name = "HLT_PVPD" ;
		out[4].bytes = bytes ;

		//-----------fill emc towers -----------//
		bytes = gl3->FillEmc(sizeof(store_EMC), store_EMC) ;
		LOG(DBG,"bytes %d",bytes) ;
//cout<<"bytes_Emc ="<<bytes<<endl;
		out[5].buff = store_EMC ;
		out[5].name = "HLT_EMC" ;
		out[5].bytes = bytes ;

                //----------fill Nodes------------//
                bytes = gl3->FillNodes(sizeof(store_NODE),store_NODE);
                LOG(DBG,"bytes %d",bytes) ;
                out[6].buff = store_NODE ;
                out[6].name = "HLT_NODE" ;
                out[6].bytes = bytes ;

                //-----------fill di electrons---------//
		bytes = sizeof(HLT_DIEP)+(gl3->triggerDecider->hlt_diEP.nEPairs-1000)*sizeof(hlt_diElectronPair);
		out[7].buff = (char *)&gl3->triggerDecider->hlt_diEP;
		out[7].name = "HLT_DIEP";
		out[7].bytes = bytes;

                //-----------fill high pt---------//
                bytes = sizeof(HLT_HIPT)+(gl3->triggerDecider->hlt_hiPt.nHighPt-1000)*sizeof(int);
                out[8].buff = (char *)&gl3->triggerDecider->hlt_hiPt;
                out[8].name = "HLT_HIPT";
                out[8].bytes = bytes;

                //-----------fill high pt---------//
                bytes = sizeof(HLT_HF)+(gl3->triggerDecider->hlt_hF.nHeavyFragments-1000)*sizeof(int);
                out[9].buff = (char *)&gl3->triggerDecider->hlt_hF;
                out[9].name = "HLT_HF";
                out[9].bytes = bytes;

		*blobs = 10 ;
#ifdef TIMING_HLT
		timeNow = (double)(clock())/CLOCKS_PER_SEC;
		timeFill += timeNow - timeOld;
		nTracks += gl3->getNGlobalTracks();
		if(gl3->getNGlobalTracks() > 0) nEvents ++;
#endif

		pval = PROFILER(pval) ;	// 60 us

		return decision ;
	}

private:
	char store_GT[sizeof(HLT_GT)] ;	// NOTE: this determines the max amount of data in the current
					// scheme -- look at "fillTracks" above
        char store_PT[sizeof(HLT_PT)] ;
        char store_EVE[sizeof(HLT_EVE)] ;
        char store_TOF[sizeof(HLT_TOF)] ;
        char store_PVPD[sizeof(HLT_PVPD)] ;
        char store_EMC[sizeof(HLT_EMC)] ;
        char store_NODE[sizeof(HLT_NODE)] ;

} ;



#endif
