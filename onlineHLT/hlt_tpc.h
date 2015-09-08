#ifndef _HLT_tpc_H_
#define _HLT_tpc_H_

#include <GL3/hlt_entity.h>


#include "online_tracking_sector.h"
#include <DAQ_TPC/daq_tpc.h>

//#define TIMING_HLT
#ifdef TIMING_HLT
#include <time.h>
#endif

class hlt_tpc : public hlt_entity {
public:
#ifdef TIMING_HLT
  double timeNow, timeOld;
  double timeStart, timeData, timeDo;
  int nTracks;
#endif

	hlt_tpc() { 
		name = "tpc" ;
		tpc_m = 0 ;
		time_scaler_check = 0 ;
	} ;

	virtual ~hlt_tpc() {
		if(tpc_m) delete tpc_m ;
	} ;


	char beamlineFile[256];
	char GainParameters[256] ;

	int run_config() {
		if(tpc_m) delete tpc_m ;
		char mapName[256] ;

		if(standalone) {
			sprintf(mapName,"../HLT/tpcHitMap_sector%d.bin",sector) ;
			sprintf(beamlineFile, "../HLT/beamline");
			sprintf(GainParameters, "../HLT/GainParameters");
			tpc_m = new online_tracking_sector(SL3_BFIELD,sector-1,mapName,"../HLT/HLTparameters","../HLT/sc.asc","../HLT/beamline") ;
		}
		else {
			sprintf(mapName,"/RTS/conf/hlt/tpcHitMap_sector%d.bin",sector) ;
			sprintf(beamlineFile, "/RTS/conf/hlt_config/beamline");
			sprintf(GainParameters, "/RTS/conf/hlt_config/GainParameters");
			tpc_m = new online_tracking_sector(SL3_BFIELD,sector-1,mapName,"/RTS/conf/hlt/HLTparameters","/RTS/conf/sc/sc.asc", "/RTS/conf/hlt_config/beamline") ;
		}

		return 0 ;
	}

	int run_start() { 

#ifdef TIMING_HLT
		timeStart = 0.;
		timeData = 0.;
		timeDo = 0.;
		nTracks = 0;
#endif

		if(!standalone) {
			tpc_m->readScalers("/RTS/conf/sc/sc.asc") ;
		}
		tpc_m->readBeamline(beamlineFile);
		tpc_m->fDedx->ReadHLTparameters(GainParameters);

		time_scaler_check = time(0) + L3_SCALER_CHECK ;
		return 0 ; 
	} ;

	int run_stop() {
#ifdef TIMING_HLT		
		printf("evt_counter:     %i \n", evt_counter);
		printf("nTracks:     %i \n", nTracks);
		printf("\t\t time \t\t per evt (ms) \t per track (micro s) \n");
		printf("timeStart:\t %f\t %f\t %f\n", timeStart, timeStart/evt_counter*1000., timeStart/nTracks*1000000.);
		printf("timeData:\t %f\t %f\t %f\n", timeData, timeData/evt_counter*1000., timeData/nTracks*1000000.);
		printf("timeDo:\t\t %f\t %f\t %f\n", timeDo, timeDo/evt_counter*1000., timeDo/nTracks*1000000.);
		double totalTime = timeStart+timeData+timeDo;
		printf("totalTime:\t %f\t %f\t %f\n\n", totalTime, totalTime/evt_counter*1000., totalTime/nTracks*1000000.);
#endif

		return 0 ;
	}

	int event_start(void *data) {
#ifdef TIMING_HLT
	  timeOld = (double)(clock())/CLOCKS_PER_SEC;
#endif

		time_t now = time(0) ;
		if(standalone) {
			//int rich_scaler = (int) data ;	// by agreement
		  tpc_m->tpcHitMap->setSpaceChargeFromScalerCounts((int*)data);
		}
		else if(!standalone && (time_scaler_check <= now)) {
			tpc_m->readScalers("/RTS/conf/sc/sc.asc") ;
			time_scaler_check = now + L3_SCALER_CHECK ;
		}
		
		tpc_m->nHits = 0 ;

#ifdef TIMING_HLT
		timeNow = (double)(clock())/CLOCKS_PER_SEC;
		timeStart += timeNow - timeOld;
#endif

		return 0 ;
	}

	int event_data(int rts_id, int rdo, void *start, int bytes) {
#ifdef TIMING_HLT
	  timeOld = (double)(clock())/CLOCKS_PER_SEC;
#endif

		tpc_t *pTPC = (tpc_t*)start;
		tpc_m->Tracking_load_cluster(pTPC);
			

#ifdef TIMING_HLT
		timeNow = (double)(clock())/CLOCKS_PER_SEC;
		timeData += timeNow - timeOld;
#endif

		return 0 ;
	}

	int event_do(hlt_blob *out, int *blobs) {

#ifdef TIMING_HLT
		  timeOld = (double)(clock())/CLOCKS_PER_SEC;
#endif

		int bytes_used ;

		tpc_m->Tracking_track_sector((char *)output_store, bytes_used) ;
		//		printf("sector size: %i \n", bytes_used);
		if((u_int)bytes_used > sizeof(output_store)) {
			LOG(ERR,"Overwrote memory!") ;
		}

		LOG(NOTE,"tpc: s %2d: nhits %d",sector,tpc_m->nHits) ;

		out->buff = output_store ;
		out->bytes = bytes_used ;
		out->name = "sl3" ;	// always
		*blobs = 1 ;	// always for SL3 tpc

#ifdef TIMING_HLT
		  timeNow = (double)(clock())/CLOCKS_PER_SEC;
		  timeDo += timeNow - timeOld;
#endif		

		return 1 ;
	}

private:
	char output_store[2*1024*1024] ;	// internal buffer which keeps results
	time_t time_scaler_check ;
	online_tracking_sector *tpc_m ;
} ;


#endif
