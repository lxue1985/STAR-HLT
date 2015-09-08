#ifndef _HLT_TPX_H_
#define _HLT_TPX_H_

#include <GL3/hlt_entity.h>


#include "online_tracking_sector.h"

//#define TIMING_HLT
#ifdef TIMING_HLT
#include <time.h>
#endif


class hlt_tpx : public hlt_entity {
public:
#ifdef TIMING_HLT
  double timeNow, timeOld;
  double timeStart, timeData, timeDo;
  int nTracks;
#endif

	hlt_tpx() { 
		name = "tpx" ;
		rts_id = TPX_ID ;
		tpx_t = 0 ;
	} ;

	virtual ~hlt_tpx() {
		if(tpx_t) delete tpx_t ;
	} ;


	char beamlineFile[256];
	char GainParameters[256] ;
	char hltParameters[256] ;
	char mapName[256] ;

	int run_config() {
		if(tpx_t) delete tpx_t ;


		LOG(INFO,"Using config_directory %s, mag field %f",
		    hlt_config.config_directory,
		    hlt_config.l3_sc.mag_field) ;

		sprintf(mapName,"%s/tpcHitMap_sector%d.bin",hlt_config.config_directory,sector) ;
		sprintf(beamlineFile, "%s/beamline",hlt_config.config_directory);
		sprintf(GainParameters, "%s/GainParameters",hlt_config.config_directory);
		sprintf(hltParameters,"%s/HLTparameters",hlt_config.config_directory) ;
		
		tpx_t = new online_tracking_sector(hlt_config.l3_sc.mag_field,
						   sector-1,	// "sector" starts from 1
						   mapName,
						   hltParameters,
						   beamlineFile) ;
                tpx_t->tpcHitMap->setDriftVelocity(hlt_config.l3_sc.drift_velocity_east);
		// this is done in the constructor... tpx_t->readBeamline(beamlineFile);
		tpx_t->fDedx->ReadGainparameters(GainParameters);
		tpx_t->put_bField(hlt_config.l3_sc.mag_field) ;

		return 0 ;
	}

	int run_start() { 

#ifdef TIMING_HLT
		timeStart = 0.;
		timeData = 0.;
		timeDo = 0.;

		nTracks = 0;
#endif


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

		if(data) {
			tpx_t->tpcHitMap->setSpaceChargeFromScalerCounts((int*)data);
		}
		else {
			// no RICH scaler -- do nothing!
			LOG(WARN,"No RICH scaler data!") ;
		}

		tpx_t->nHits = 0 ;

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


		int nhits = tpx_t->readSectorFromESB(sector, (char *)start, bytes/4) ;


#ifdef TIMING_HLT
		timeNow = (double)(clock())/CLOCKS_PER_SEC;
		timeData += timeNow - timeOld;
#endif



		return nhits ;
	}

	int event_do(hlt_blob *out, int *blobs) {
		int ret = 0 ;	// assume 0 tracks



#ifdef TIMING_HLT
		timeOld = (double)(clock())/CLOCKS_PER_SEC;
#endif

		int bytes_used = sizeof(output_store);

		// returns ntracks...
		ret = tpx_t->Tracking_track_sector((char *)output_store, bytes_used) ;

		if((u_int)bytes_used > sizeof(output_store)) {
			LOG(ERR,"Overwrote memory: used %d, avail %d!",bytes_used,sizeof(output_store)) ;
			ret = -1 ;
		}

		LOG(NOTE,"tpx: s %2d: bytes %d/%d, nhits %d",sector,bytes_used,sizeof(output_store),tpx_t->nHits) ;

		out->buff = output_store ;
		out->bytes = bytes_used ;
		out->name = "tpx_sl3" ;	// always

		*blobs = 1 ;	// always for SL3 tpx

#ifdef TIMING_HLT
		timeNow = (double)(clock())/CLOCKS_PER_SEC;
		timeDo += timeNow - timeOld;
#endif		


		return ret ;
	}

private:
	char output_store[1024*1024] ;	// internal buffer which keeps tracks; 1 MB is a lot
	online_tracking_sector *tpx_t ;
} ;


#endif
