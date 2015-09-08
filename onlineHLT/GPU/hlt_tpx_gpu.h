#ifndef _HLT_TPX_GPU_H_
#define _HLT_TPX_GPU_H_

#include <GL3/hlt_entity.h>


#include "online_GPU_tracking.h"
#include "../online_tracking_sector.h"


//#define TIMING_HLT
#ifdef TIMING_HLT
#include <time.h>
#endif

class hlt_tpx_gpu : public hlt_entity {
public:
#define GPUEVENTBUFFER 6
	int cores;
#ifdef TIMING_HLT
  double timeNow, timeOld;
  double timeStart, timeData, timeDo;
  int nTracks;
#endif

	hlt_tpx_gpu() { 
		name = "tpx_gpu" ;
		tpx_gpu_t = 0 ;
		time_scaler_check = 0 ;
		cores=GPUEVENTBUFFER;
	} ;

	virtual ~hlt_tpx_gpu() {
		if(tpx_gpu_t) delete tpx_gpu_t ;
	} ;




	int run_config() {
		if(tpx_gpu_t) delete tpx_gpu_t ;

		if(standalone) {
			tpx_gpu_t = new online_GPU_tracking() ;
		}
		else {
			tpx_gpu_t = new online_GPU_tracking() ;
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

		tpx_gpu_t->setup_memory(cores);
		time_scaler_check = time(0) + L3_SCALER_CHECK ;
		return 0 ; 
	} ;

	int run_stop() {
	tpx_gpu_t->free_memory();
		return 0 ;
	}

	int event_start(void *data) {


		return 0 ;
	}

	int event_data(int rts_id, int rdo, void *start, int bytes) {

	int i;
	for(i=0;i<tpx_gpu_t->dup;i++){
	tpx_gpu_t->hit_counter[i]=0;
	}

		return 0 ;
	}
	
	
	int load_hit_from_CPU_tracker(int position, online_tracking_sector* tracker) {

		tpx_gpu_t->load_hit_to_buffer(position,tracker);
		return 0 ;
	}	

	int copy_GPU_back( ) {
		tpx_gpu_t->copy_GPU_back();
		return 0 ;
	}	
	
	int load_track_to_CPU_tracker(int position) {
		tpx_gpu_t->build_track(position);
		return 0 ;
	}		

	int GPU_tracking() {
		tpx_gpu_t->set_balance(30);
		tpx_gpu_t->fire_GPU(30);
		return 0 ;
	}	
	
	int event_do(hlt_blob *out, int *blobs) {


		return 1 ;
	}

public:
	//char output_store[2*1024*1024] ;	// internal buffer which keeps results
	time_t time_scaler_check ;
	online_GPU_tracking *tpx_gpu_t ;
} ;


#endif
