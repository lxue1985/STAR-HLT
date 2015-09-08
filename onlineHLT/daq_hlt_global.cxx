#include <sys/types.h>
#include <errno.h>
#include <getopt.h>
#include <fcntl.h>
#include <poll.h>
#include <assert.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <signal.h>

#include <RC_Config.h>
#include <rtsMonitor.h>
#include <rtsLog.h>
#include <tasks.h>
#include <cmds.h>
#include <iccp.h>
#include <SUNRT/ipcQLib.hh>
#include <SFS/sfs_index.h>

#include <DAQ_READER/daqReader.h>
#include <DAQ_READER/daq_dta.h>
#include <DAQ_TPX/daq_tpx.h>
#include <DAQ_TRG/daq_trg.h>
#include <DAQ_SC/daq_sc.h>

#include <L3_SUPPORT/l3_support.h>

#include "daq_hlt.h"
#include "hlt_io_class.h"
#include "../HLT_TOF/tof_hlt_sector.h"

#include "online_tracking_collector.h"
#include "online_tracking_sector.h"
#include "lxAlloc.hh"
#include <iomanip>

extern const char *BuildDate ;

#define TOKENS_MAX_IX	26

#define DET_MAX_MOVERS	12

//#define RUN8_OR_BEFORE

#define TIMING_HLT
#ifdef TIMING_HLT
#include <time.h>
double timeNow, timeOld;
double timeResetEvent, timeReadEmc, timeReadTof, timeReadTpc, timeFinalize, timeDecide, timeFill;
int nTracks;
#endif

char* daqFileName;

static struct gl3_mon {
	double kbEvb ;	// output, accepted
	double kbL3 ;	// input
	int busy ;
	rtsMonStruct rts ;
	rtsMonL1 trg ;
} Mon ;

static struct tokens {
	time_t time_in ;
	u_int status ;

	u_int eventNumber ;	// from EVB/Trigger data
	u_int daq_TrgId ;	// from EVB/Trigger data

	char *data_mem[TOKENS_MAX_IX] ;	// TPX: 0..23, EVB 24, TOF 25
	int data_bytes[TOKENS_MAX_IX] ;


} tokens[4096] ;

static EvbChooser evb_chooser ;

static struct evtSender_t {
	sfs_index sfs ;

	u_char is_used ;

	u_int dstNode ;
	u_int dstTask ;

} evtSender[DET_MAX_MOVERS] ;

static int global_run;
static int global_eventNumber;
static tof_hlt_sector *tof_tracker;

static online_tracking_collector *gl3 ;
static int gl3_evt_cou ;
static int gl3_accept_cou ;
static int gl3_id ;
static int gl3_run_number ;
static int gl3_mode ;
static int gl3_force_rts ;
static char gl3_wdir[1024] ;
static time_t time_to_scalers ;

static ipcQClass *gl3_que ;
static lxAlloc *lx_mem ;
static hlt_io *hlt_recv ;
static sfs_index sfs_debug ;

static void *sender_thread(void *arg) ;
static void *checker_thread(void *arg)  ;
static void *monitor_thread(void *arg) ;

static int do_save(int t, char *buff, int bytes, int accepted) ;
static int clr_token(int t, int bad) ;
static int full_read(int desc, char *buff, int bytes) ;

static char scalerFile[128] ;

static int run_stop() 
{
#ifdef TIMING_HLT
  cout<<"gl3_evt_cou:     "<<gl3_evt_cou<<endl;
  cout<<"nTracks:         "<<nTracks<<endl;
  cout<<"                 "<<setw(16)<<"time"<<setw(16)<<"per event (ms)"<<setw(16)<<"per track (micro s)"<<endl;
  cout<<"timeResetEvent:  "<<setw(16)<<timeResetEvent<<setw(16)<<timeResetEvent/gl3_evt_cou*1000.<<setw(16)<<timeResetEvent/nTracks*1000000.<<endl;
  cout<<"timeReadEmc:     "<<setw(16)<<timeReadEmc<<setw(16)<<timeReadEmc/gl3_evt_cou*1000.<<setw(16)<<timeReadEmc/nTracks*1000000.<<endl;
  cout<<"timeReadTof:     "<<setw(16)<<timeReadTof<<setw(16)<<timeReadTof/gl3_evt_cou*1000.<<setw(16)<<timeReadTof/nTracks*1000000.<<endl;
  cout<<"timeReadTpc:     "<<setw(16)<<timeReadTpc<<setw(16)<<timeReadTpc/gl3_evt_cou*1000.<<setw(16)<<timeReadTpc/nTracks*1000000.<<endl;
  cout<<"timeFinalize:    "<<setw(16)<<timeFinalize<<setw(16)<<timeFinalize/gl3_evt_cou*1000.<<setw(16)<<timeFinalize/nTracks*1000000.<<endl;
  cout<<"timeDecide:      "<<setw(16)<<timeDecide<<setw(16)<<timeDecide/gl3_evt_cou*1000.<<setw(16)<<timeDecide/nTracks*1000000.<<endl;
  cout<<"timeFill:        "<<setw(16)<<timeFill<<setw(16)<<timeFill/gl3_evt_cou*1000.<<setw(16)<<timeFill/nTracks*1000000.<<endl;
  double totalTime = timeResetEvent+timeReadEmc+timeReadTof+timeReadTpc+timeFinalize+timeDecide+timeFill;
  cout<<"  total time:    "<<setw(16)<<totalTime<<setw(16)<<totalTime/gl3_evt_cou*1000.<<setw(16)<<totalTime/nTracks*1000000.<<endl;
#endif

	LOG(INFO,"GL3 %d: stopping run %08u",gl3_id,gl3_run_number) ;

	sfs_debug.umount() ;

	for(int i=0;i<DET_MAX_MOVERS;i++) {
		if(evtSender[i].is_used == 0) continue ;

		LOG(TERR,"Unmounting sfs %d",i) ;
		evtSender[i].sfs.umount() ;
		evtSender[i].is_used = 0 ;

	}

	// flush and close the QA file!
	gl3->triggerDecider->flushQA() ;
	gl3->triggerDecider->closeQA() ;
	
	return 0 ;
}


static int run_start(int run)
{
#ifdef TIMING_HLT
  timeResetEvent = 0.;
  timeReadEmc = 0.;
  timeReadTof = 0.;
  timeReadTpc = 0.;
  timeFinalize = 0.;
  timeDecide = 0.;
  timeFill = 0.;
#endif
	char s_name[128] ;


	gl3_run_number = run ;


	for(int i=0;i<4096;i++) {
		clr_token(i,0) ;
	}


	Mon.rts.tknIn = -1 ;
	Mon.rts.tknOut = -1 ;
	Mon.rts.tknBad = -1 ;
	Mon.rts.couEvtsRun = 0 ;
	Mon.rts.couEvtsIn = 0 ;
	Mon.rts.couEvtsBad = 0 ;
	Mon.rts.busy = 0 ;
	Mon.rts.evtsSec = 0 ;
	Mon.rts.kbSecEvb = 0 ;
	Mon.rts.kbSecAux = 0 ;

	Mon.kbEvb = 0.0 ;
	Mon.kbL3 = 0.0 ;

	memset(&Mon.trg,0,sizeof(Mon.trg)) ;


	
	if(gl3_mode) {	// local from file
	  string fileName(daqFileName);
	  fileName.replace(0, fileName.find("/st_"), "");
	  fileName.replace(fileName.find(".daq"), 4, ".sfs");
	  sprintf(s_name,"%s/%s",gl3_wdir,fileName.data()) ;	// in local dir...
	}
	else {	// from net
		sprintf(s_name,"/RTScache/data/gl3_%08d_%02d.sfs",gl3_run_number,gl3_id) ;
	}

	errno =0 ;
	if(sfs_debug.mount(s_name, O_WRONLY|O_CREAT|O_TRUNC,0666)< 0) {
		LOG(ERR,"GL3 %d: sfs_debug.mount(%s): [%s]",gl3_id, s_name, strerror(errno)) ;
	}


	gl3_evt_cou = 0 ;
	gl3_accept_cou = 0 ;

#ifdef TIMING_HLT
	nTracks = 0;
#endif

	evb_chooser.nservers = 0 ;

	if(gl3_mode) {	// local rom file
	  string fileName(daqFileName);
	  fileName.replace(0, fileName.find("/st_"), "");
	  fileName.replace(fileName.find(".daq"), 4, ".asc");
	  sprintf(s_name,"%s/%s",gl3_wdir,fileName.data()) ;	// in local dir...
	}
	else {	// from net
		sprintf(s_name,"/RTScache/data/gl3_%08d_%02d.asc",gl3_run_number,gl3_id) ;


		// open run configuration file
		
		STAR_CFG *conf = (STAR_CFG *)malloc(sizeof(STAR_CFG)) ;
				
		int fd = open("/RTS/conf/handler/cfg_04000",O_RDONLY,0666) ;
		if(fd < 0) {
			LOG(ERR,"Can't open RC configuration file [%s]",strerror(errno)) ;
		}
		else {
			extern void swapSTAR_CFG(STAR_CFG *blah) ;

			int ret = read(fd,(char *)conf,sizeof(STAR_CFG)) ;
			if(ret != sizeof(STAR_CFG)) {
				LOG(ERR,"error reading RC configuration [%s]",strerror(errno)) ;
			}
			else {
				swapSTAR_CFG(conf);
				evb_chooser.configure(conf,0) ;
			}
		}

		free(conf) ;


	}

#if 0
	LOG(TERR,"Would have seen %d EVBs but zapping to 0!",evb_chooser.nservers) ;
	evb_chooser.nservers = 0 ;
#endif
	for(int i=0;i<DET_MAX_MOVERS;i++) evtSender[i].is_used = 0 ;


	for(int i=0;i<evb_chooser.nservers;i++) {
		int ret = evtSender[i].sfs.mount(evb_chooser.servers[i].ip, evb_chooser.servers[i].port) ;
		if(ret < 0) {
			LOG(ERR,"Can't SFS-mount EVB server %d",i) ;
			continue ;
		}
		else {
			LOG(TERR,"Mounted EVB %d: node 0x%04X:%d",i,evb_chooser.servers[i].node,
			    evb_chooser.servers[i].task) ;
		}
		evtSender[i].is_used = 1 ;
		evtSender[i].dstNode = evb_chooser.servers[i].node ;
		evtSender[i].dstTask = evb_chooser.servers[i].task ;
	}


	time_to_scalers = time(0) + L3_SCALER_CHECK ;
	int rich7 = gl3->readScalers(scalerFile);
	//let's open this QA file
	gl3->triggerDecider->setQA(s_name) ;

	int all = gl3->triggerDecider->triggerOnAllEvents ;
	LOG(INFO,"GL3 %d: starting run %08u, EVBs %d, current dir %s, scaler %u, all events %d",
	    gl3_id, gl3_run_number, evb_chooser.nservers, gl3_wdir,rich7, all) ;

	return 0 ;
}


static int do_config()
{
	if(gl3) {
		delete gl3 ;
		gl3 = 0 ;
	}

	if(tof_tracker) {
	  delete tof_tracker ;
	  tof_tracker = 0 ;
	}
	tof_tracker = new tof_hlt_sector();

	// need to reload emc stuff from /RTS/conf...
	char parametersDirectory[64] ;
	
	if(gl3_mode && !gl3_force_rts) 
	  sprintf(parametersDirectory, "%s", "../HLT/");
	else 
	  sprintf(parametersDirectory, "%s", "/RTS/conf/tpx/");

	gl3 = new online_tracking_collector(szGL3_default_mHits, szGL3_default_mTracks, parametersDirectory) ;
	
	return 0 ;

}

static int clr_token(int t, int bad)
{

	for(int i=0;i<TOKENS_MAX_IX;i++) {
		if(tokens[t].data_mem[i]) {
			lx_mem->free(tokens[t].data_mem[i]) ;
			tokens[t].data_mem[i] = 0 ;
		}
	}

	Mon.rts.couEvtsIn-- ;
	Mon.rts.tknOut = t ;

	if(bad) {
		Mon.rts.tknBad = t ;
		Mon.rts.couEvtsBad++ ;
	}

	memset((void *)(tokens+t),0,sizeof(tokens[0])) ;

	return 0 ;
}

static int do_track(int t, char **out, int *ret_bytes)
{
	int ret ;
	static char *buffer ;
	static int b_alloced ;
	int decision = 0 ;
	u_int evt_seq ;

	// allocate a global L3_GTD buffer on first run...
	if(!buffer) {
		b_alloced = 4*1024*1024 ;	//4 MB?
		buffer = (char *) valloc(b_alloced) ;
	}

	*out = buffer ;	// will return this buffer...

	gl3_evt_cou++ ;
	evt_seq = tokens[t].eventNumber ;

	LOG(NOTE,"GL3 %d: resetting event %d",gl3_id,gl3_evt_cou) ;

#ifdef TIMING_HLT
	timeNow = (double)(clock())/CLOCKS_PER_SEC;
#endif

	gl3->resetEvent() ;	// get ready for new event...

#ifdef TIMING_HLT
	timeOld = timeNow;
	timeNow = (double)(clock())/CLOCKS_PER_SEC;
	timeResetEvent += timeNow - timeOld;
#endif

	// grab the Trigger data from EVBs bank and decode
	if(tokens[t].data_mem[L3_SUPPORT_EVB_IX]) {
	  //		gl3->setToken(tokens[t].daq_TrgId);
	  gl3->trigger = tokens[t].daq_TrgId;
	  gl3->emc->readFromGl3Trg((gl3_trg_send_t *)tokens[t].data_mem[L3_SUPPORT_EVB_IX]) ;
	}

#ifdef TIMING_HLT
	timeOld = timeNow;
	timeNow = (double)(clock())/CLOCKS_PER_SEC;
	timeReadEmc += timeNow - timeOld;
#endif

	//--------------add by Bingchu -----------
	if(tokens[t].data_mem[L3_SUPPORT_TOF_IX]) {
		//TofSend *ts = (TofSend *)tokens[t].data_mem[L3_SUPPORT_TOF_IX];
		gl3->tof->readFromTofMachine(tokens[t].data_mem[L3_SUPPORT_TOF_IX]);
		//LOG(TERR," trace 2003, receiving tof hits: %d, 1st hit tdc: %f, last hit tdc: %f", ts->nHits,ts->cell[0].tdc,ts->cell[ts->nHits-1].tdc);
	}

#ifdef TIMING_HLT
	timeOld = timeNow;
	timeNow = (double)(clock())/CLOCKS_PER_SEC;
	timeReadTof += timeNow - timeOld;
#endif

	// grab sector tracks, if present, and readin
	for(int i=0;i<24;i++) {
		if(tokens[t].data_mem[i] == 0) continue ;

		LOG(DBG,"GL3 %d: adding sector %d",gl3_id,i+1) ;

		// load tracks from sector's L3_SECTP
		ret = gl3->Tracking_readSectorTracks(tokens[t].data_mem[i]) ;

		LOG(NOTE,"GL3 %d: Token %d: sector %d: tracks read %d",gl3_id,t,i+1,ret) ;
	}

#ifdef TIMING_HLT
	timeOld = timeNow;
	timeNow = (double)(clock())/CLOCKS_PER_SEC;
	timeReadTpc += timeNow - timeOld;
#endif

	gl3->finalizeReconstruction() ;	// I hope this is what is necessary...

#ifdef TIMING_HLT
	timeOld = timeNow;
	timeNow = (double)(clock())/CLOCKS_PER_SEC;
	timeFinalize += timeNow - timeOld;
#endif


	// do the decision
	decision = gl3->triggerDecider->decide(evt_seq) ;	// and now what?
#ifdef TIMING_HLT
	timeOld = timeNow;
	timeNow = (double)(clock())/CLOCKS_PER_SEC;
	timeDecide += timeNow - timeOld;
#endif

	if(decision) {
		gl3_accept_cou++ ;
	}

	LOG(TERR,"GL3 %d: evt %d, seq %d: T %4d: tracks %d: Trg 0x%08X: decision 0x%X",gl3_id,gl3_evt_cou,evt_seq,t,gl3->getNGlobalTracks(),tokens[t].daq_TrgId,decision) ;

#ifdef TIMING_HLT
	timeOld = timeNow;
	timeNow = (double)(clock())/CLOCKS_PER_SEC;
#endif
	// fill the output banks...
	ret = gl3->fillTracks(b_alloced,buffer,t) ;	// I hope this fills L3_GTD correctly...

#ifdef TIMING_HLT
	timeOld = timeNow;
	timeNow = (double)(clock())/CLOCKS_PER_SEC;
	timeFill += timeNow - timeOld;
#endif

	((struct L3_GTD *)buffer)->bh.w9 = evt_seq ;	// this is where I put this unique number -- HACK!
	((struct L3_GTD *)buffer)->bh.crc = decision ;	// this is where I put this decision -- HACK!

	LOG(DBG,"GL3 %d: bank length is %d",gl3_id,ret) ;
	*ret_bytes = ret ;

	// flush QA every now and then
	if((gl3_accept_cou % 10)==0) gl3->triggerDecider->flushQA() ;

#ifdef TIMING_HLT
	nTracks += gl3->getNGlobalTracks();
#endif
	return decision ;	
}


int process_data(int desc, u_int sector, int t, char *from, int bytes)
{
	int run ;
	int s ;
	char *mem ;
	int ret ;
	int err ;
	int eventNumber ;
	u_int daq_TrgId = 0 ;

	err = 0 ;
	eventNumber = -1 ;	// mark as unknown...



	if(sector<TOKENS_MAX_IX) s = sector ;
	else {
		return -100 ;	// will close socket
	}


	if((t<0) || (t>=4096)) {
		return -101 ;	// will close socket
	}

	// check for overflow!
	mem = tokens[t].data_mem[s] ;

	if(mem) {	// overflow!
		tokens[t].status |= 0x10 ;		// mark as bad
		err = 1 ;				// set error
		mem = (char *) malloc(bytes) ;		// need temp alloc to read in socket
	}
	else {
		mem = lx_mem->alloc(bytes) ;		// fresh allocation

	}

	if(desc<0) {	// this is in-memory, so copy
		memcpy(mem, from, bytes) ;
	}
	else {
		ret = full_read(desc,mem,bytes) ;
		if(ret != bytes) {
			if(err) free(mem) ;
			else lx_mem->free(mem) ;

			if(ret <= 0) return ret ;	// return the same code if negative
			return -102 ;
		}
	}


	if(err) {
		free(mem) ;	// relase the temporary allocation
		return 1 ;	// don't return <= 0 because I don't want to close the socket for an overflow
	}

	switch(s) {
	case L3_SUPPORT_EVB_IX :
		// extract necessary things
		run = ((gl3_trg_send_t *)mem)->run_number ;
		eventNumber= ((gl3_trg_send_t *)mem)->eventNumber ;
		daq_TrgId = ((gl3_trg_send_t *)mem)->daq_TrgId ;
		break ;
	case L3_SUPPORT_TOF_IX :
		run = global_run;
		eventNumber = global_eventNumber;
		//LOG(TERR,"\t run:%d event:%d ",run, eventNumber);
		break;
	default :	// TPX
		run = ((L3_SECTP *)mem)->bh.w9 ;
		break ;
	}

	if(run != gl3_run_number) {
		// this will zap all structures so I need to temporary copy this data!

		char *tmp = (char *)malloc(bytes) ;	// temp alloc
		memcpy(tmp, mem, bytes) ;		// temp copy

		lx_mem->free(mem) ;			// free what we had

		run_stop() ;
		run_start(run) ;	// this frees all allocation!


		mem = lx_mem->alloc(bytes) ;		// redo allocation
		memcpy(mem, tmp, bytes) ;		// recopy back
		free(tmp) ;				// free the temp alloc
	}


	tokens[t].data_mem[s] = mem ;
	tokens[t].data_bytes[s] = bytes ;

	if(tokens[t].status == 0) {	// first time...
		tokens[t].time_in = time(0) ;
		tokens[t].status = 1 ;	// first

		Mon.rts.tknIn = t ;
		Mon.rts.couEvtsIn++ ;
		Mon.rts.couEvtsRun++ ;
		Mon.rts.couEvtsAll++ ;
	}


	if(eventNumber >= 0) {	// protects from TPX...
		tokens[t].eventNumber= eventNumber ;
		tokens[t].daq_TrgId = daq_TrgId ;
	}

	// figure out if all contribs are done?
	u_int mask = 0 ;
	for(int i=0;i<TOKENS_MAX_IX;i++) {
		if(tokens[t].data_mem[i]) mask |= (1<<i) ;
	}

	if(mask == 0x001FFFFFF) {	// all sectors and Evb..
		ic_msg ic ;

		tokens[t].status |= 2 ;	// all elements are in!

		ic.head.token = t ;
		ic.head.daq_cmd = GL3_ANNOUNCE_TOKEN ;

		LOG(NOTE,"Sending token %d",t) ;
		gl3_que->send(&ic,sizeof(ic),1) ;
		

	}


	return 1 ;
}

/*
	This runs only when emulating from a file
*/
void *sender_thread(void *arg)
{
	char *fname = (char *) arg ;
	daqFileName = fname;
	online_tracking_sector *sl3[24] ;
	char *tmp_mem = (char *)malloc(SL3_TRK_BYTES) ;

	for(int i=0;i<24;i++) {
		char mapName[64] ;
		char HLTparameters[64], scalerFile[64];

		if(gl3_mode && !gl3_force_rts) {
			sprintf(mapName,"tpcHitMap_sector%d.bin",i+1) ;
			sprintf(HLTparameters,"HLTparameters");
			sprintf(scalerFile, "sc.asc");
		}
		else {
			sprintf(mapName,"/RTS/conf/tpx/tpcHitMap_sector%d.bin",i+1) ;
			sprintf(HLTparameters,"/RTS/conf/tpx/HLTparameters");
			sprintf(scalerFile, "/RTS/conf/sc/sc.asc");
		}

		sl3[i] = new online_tracking_sector(SL3_BFIELD,i,mapName,HLTparameters, scalerFile) ;
	}


	daqReader *dr = new daqReader(fname) ;

	while(dr->get(0,0)) {
		int got_some = 0 ;
		daq_dta *dd ;
		daq_dta *dd2 = dr->det("sc")->get() ;
		dd2->iterate();
		sc_t *sc_p = (sc_t *) dd2->Void ;
		//		cout<<"scaler count: "<<sc_p->rich_scalers[7]<<endl;
		tof_tracker->new_event();

		for(int s=0;s<24;s++) {	// note: we count from 0 here...
			int got_sector ;

			LOG(DBG,"Reading sector %d",s+1) ;

			sl3[s]->tpcHitMap->setSpaceChargeFromScalerCount(sc_p->rich_scalers[7]);

#ifndef RUN8_OR_BEFORE
			sl3[s]->nHits = 0 ;
			got_sector = 0 ;
			for(int r=1;r<=6;r++) {
				dd = dr->det("tpx")->get("cld_raw",s+1,r) ;
				while(dd && dd->iterate()) {

					got_sector = 1 ;

					if(dd->ncontent == 0) continue ;

					LOG(DBG,"Loading sector %d, %d bytes",s+1,dd->ncontent) ;
					sl3[s]->readSectorFromESB(s+1, (char *)dd->Void, dd->ncontent/4) ;

				}
			}
#else
			tpc_t *  pTPC;
			dd = dr->det("tpx")->get("legacy",s+1);
			if(dd) 
			  {
			    dd->iterate();
			    pTPC = (tpc_t *)dd->Void;        
			  }

			else 
			  {
			    dd = dr->det("tpc")->get("legacy",s+1);
			    if(dd) 
			      {				
				dd->iterate();
				pTPC = (tpc_t *)dd->Void;
			      }			    
			    else 
			      {
				pTPC = NULL;
			      }
			  }
			got_sector = 0;
			if(pTPC)
			  {
			    got_sector = 1;
			    if(!pTPC->has_clusters) 
			      got_sector = 0;
			  }
			  
			sl3[s]->nHits = 0 ;
			if(got_sector == 0) continue ;
			sl3[s]->Tracking_load_cluster(pTPC);
#endif

			if(got_sector == 0) continue ;	// this sector was not in the data at all..


			got_some = 1 ;


			int bytes_used ;
			sl3[s]->Tracking_track_sector(tmp_mem, bytes_used) ;

			((L3_SECP *)tmp_mem)->bh.token = dr->token ;
			((L3_SECP *)tmp_mem)->bh.w9 = dr->run ;
			
			process_data(-1, s, dr->token, tmp_mem, bytes_used) ;

		}

		LOG(TERR,"from file: seq %d, token %d, run %d: had TPX %d",dr->seq,dr->token,dr->run,got_some) ;

		if(got_some) {	// if any sector, do the trigger as well
			// now do trigger!
			dd = dr->det("trg")->get("raw") ;
			if(dd && dd->iterate()) {


				gl3_trg_parse((char *)dd->Void, dd->ncontent, (gl3_trg_send_t *)tmp_mem) ;

				gl3_trg_send_t *g = (gl3_trg_send_t *) tmp_mem ;

				LOG(NOTE,"token %d/%d, seq %d",dr->token,g->token,g->eventNumber) ;

				((gl3_trg_send_t *)tmp_mem)->token = dr->token;
				((gl3_trg_send_t *)tmp_mem)->run_number = dr->run;
				
				process_data(-1, L3_SUPPORT_EVB_IX, dr->token, tmp_mem, sizeof(gl3_trg_send_t)) ;


			}

		}
		int sector_from_daq = 1;
		int bytes_used = 0;
		int st = 0;	
		char *buff_s = (char *)malloc(SL3_TRK_BYTES);
		//int ret;
		u_int len;
		//ret = createShmSegment(SL3_SEGMENT_ID(0), SL3_TRK_BYTES) ; 
		//tofsend_t *ts = (tofsend_t *)tmp_mem;
		//for(int j=0;j<4;j++) ts->tmp_mem[j] = getShmPtr(SL3_SEGMENT_ID(0),0);
		//int size =0;	
		//TofSend *ts = (TofSend *) tmp_mem;	
		for(int rdo=0;rdo<4;rdo++) {
			//LOG(TERR,"\t trace 10001" );
		  	dd = dr->det("tof")->get("raw",sector_from_daq,rdo+1) ;
		  while(dd && dd->iterate()) {
			global_run = dr->run;
			memcpy(buff_s,dd->Byte, dd->ncontent);
			len = dd->ncontent/4;	
			st += dd->ncontent;	
			tof_tracker->read_rdo_event(rdo+1,buff_s,len);
			
			/*memcpy(buff_s, dd->Byte, dd->ncontent);
			uint *ui2 = (u_int *)buff_s;
			for(int j=0;j<40;j++) LOG(TERR," trace 10003 buff_s\t %2d:0x%08X",j,ui2[j]);
			buff_s += dd->ncontent;*/
			//memcpy(ts->tmp_mem[rdo], dd->Byte, dd->ncontent);
			//tmp_mem -= st*4;
			//dword[2*rdo] = ((u_int)(buff_s - shared_mem))/4 ; // offset
            //dword[2*rdo+1] = dd->ncontent/4 ; // words, by usual DAQ standard
			//ts->tmp_mem[rdo] += dd->ncontent ;
			//uint *uii = (u_int *)ts->tmp_mem[rdo];
			//`for(int j=0;j<20;j++) LOG(TERR," trace 10002 tofsend\t %2d:0x%08X",j,uii[j]);
			//LOG(TERR,"\t trace 10001" );
			//LOG(TERR,"\t tof hits: %d, bytes_used: %d, tdc 0: %f", ts->nHits,bytes_used, ts->cell[0].tdc);
			//buff_s += dd->ncontent;
			//LOG(TERR,"\t trace 10002" );
		  }			
		}
	    bytes_used = tof_tracker->do_event(tmp_mem,SL3_TRK_BYTES);
	    //if(tokens[dr->token].data_mem[1]) LOG(TERR,"10002.1 TPC mem found");
	    process_data(-1, L3_SUPPORT_TOF_IX, dr->token, tmp_mem, bytes_used);
        //LOG(TERR,"\t trace 10003" );
		free (buff_s);

		
	}

	LOG(WARN,"End Of File reached...") ;
	sleep(3) ;	// give tracking time...
	run_stop();
	exit(0) ;

	for(;;) sleep(1000) ;	// never return...
	
	return 0 ;
}

int full_read(int desc, char *buff, int bytes)
{
	int ret_bytes = bytes ;
	int read_bytes = 0 ;
	int ret ;

	do {
		errno = 0 ;
		ret = read(desc, buff, bytes) ;
		if(ret <= 0) {
			if(errno) LOG(ERR,"GL3 %d: read failed, expect %d, got %d [%s]",gl3_id, ret_bytes, read_bytes, strerror(errno)) ;
			return ret ;
		}
		buff += ret ;
		bytes -= ret ;
		read_bytes += ret ;

	} while(bytes) ;

	return ret_bytes ;
}

int read_cback(int desc)
{
	iccp2k i2c ;
	int ret ;


	ret = full_read(desc,(char *)&i2c,sizeof(i2c)) ;
	if(ret <= 0) return ret ;

	int t = i2c.token ;
	int bytes = i2c.words * 4 ;
	int s = i2c.srcTask ;

	if(s==100) s = L3_SUPPORT_EVB_IX ;	// rewrite old style...
	LOG(NOTE,"GL3 %d: got T %4d: component %2d: bytes %d (%d with i2c)",gl3_id,t,s,bytes,bytes+ret) ;


	Mon.kbL3 += bytes + sizeof(i2c) ;

	ret = process_data(desc, s,  t, 0, bytes) ;

	return ret ;
}

	
	
int main(int argc, char *argv[])
{
	int c ;
	char fname[128] ;

	extern char *optarg;

	// defaults
	fname[0] = 0 ;
	gl3_id = 0 ;
	gl3_force_rts = 0 ;

	rtsLogLevel(WARN) ;
	
	// parse line args...

	while((c = getopt(argc, argv, "i:d:f:R")) != EOF) {
		switch(c) {
		case 'i' :
			gl3_id = atoi(optarg) ;
			break ;
		case 'd' :
			rtsLogLevel(optarg) ;
			break ;
		case 'f' :
			strcpy(fname,optarg) ;
			break ;
		case 'R' :
			gl3_force_rts = 1 ;
			break ;
		}
	}
	

	if(fname[0]) {	// log to stderr if running standalone from a file
		rtsLogOutput(RTS_LOG_STDERR) ;
		gl3_mode = 1 ;	// local
	}
	else {
		gl3_mode = 0 ;	// network...
		rtsLogAddDest(RTS_DAQMAN, RTS_LOG_PORT_L3) ;		
	}


	getcwd(gl3_wdir,sizeof(gl3_wdir)) ;
	signal(SIGPIPE,SIG_IGN) ;	// just in case...

	LOG(INFO,"GL3 %d starting from %s, current dir %s, build %s",gl3_id,gl3_mode?"file":"network",gl3_wdir,BuildDate) ;



	
	// main queue...
	gl3_que = new ipcQClass(GL3_TASK + gl3_id,1,0) ;


	char *arena = (char *)valloc(24*1024*4096) ;

	lx_mem = new lxAlloc(arena, 24*1024, 4096,1) ;

	do_config() ;	// create GL3 tracker and setup...


	pthread_t thr_id ;


	if(gl3_mode) {	// local
		// data reader thread, pass the filename
		pthread_create(&thr_id, 0, sender_thread, fname) ;

	}
	else {	// online!

		// checker thread
		pthread_create(&thr_id, 0, checker_thread, 0) ;

		// checker thread
		pthread_create(&thr_id, 0, monitor_thread, 0) ;

		// network reader thread
		hlt_recv = new hlt_io ;
		hlt_recv->use_read_cback(read_cback) ;
		hlt_recv->start_recv(GL3_PORT_START+gl3_id,24*6) ;	// large number of connections...
		
	}
	



	
	for(;;) {	// forever...
	int ret ;
	ic_msg ic ;
	int t ;


	ret = gl3_que->receive(&ic, sizeof(ic), 1) ;


	t = ic.head.token ;

	LOG(DBG,"Got token %d...",t) ;
	time_t now = time(0) ;

	switch(ic.head.daq_cmd) {
	case GL3_ANNOUNCE_TOKEN :
		char *data ;
		int decision ;
		int bytes ;
		int accepted ;

		accepted = 0 ;
		bytes = 0 ;
		data = 0 ;

		Mon.busy = 1 ;

		if(gl3_evt_cou>10000) 
		  {
		    run_stop();
		    return 0;	
		  }
		
		if(tokens[t].eventNumber >= 0x7FFFFFFF) {
			LOG(WARN,"GL3 %d: T %4d: eventNumber 0x%08X - discarding",gl3_id,t,tokens[t].eventNumber) ;
		}
		else {
			decision = do_track(t, &data, &bytes) ;

			if(decision < 0) {
				LOG(ERR,"Problem with token %d -- discarding",t) ;
			}
			else if(decision) {
				accepted = 1 ;
				//do_save(t, data, bytes) ;
				Mon.trg.trgs[32].accepted++ ;
			}
			else {
				Mon.trg.trgs[32].aborted++ ;
			}
		}

#if 0	// used for testing
		if(evb_chooser.nservers) {
			if((t % 11)==0) {
				LOG(TERR,"T %4d: overiding decision!",t) ;
				accepted = 1 ;
			}
		}
#endif
		do_save(t, data, bytes, accepted) ;
		clr_token(t,0) ;


		if(gl3_mode == 0) {	// only when online...
			if(time_to_scalers <= now) {
				int rich7 = gl3->readScalers(scalerFile) ;
				time_to_scalers = now + L3_SCALER_CHECK ;
				LOG(TERR,"GL3 %d: scaler %u",gl3_id,rich7) ;
				gl3->triggerDecider->writeScalers();
			}
		}

		Mon.busy = 0 ;

		break ;
	case DAQ_TIMEOUT :

		u_int sta = 0 ;
		for(int i=0;i<25;i++) {
			if(tokens[t].data_mem[i]) sta |= (1<<i) ;
		}

		LOG(ERR,"GL3 %d: expiring T %4d: status 0x%02X: 0x%07X",gl3_id,t,tokens[t].status,sta) ;

		do_save(t, data, bytes, 0) ;
		clr_token(t, 1) ;


		break ;
	}


	}

	return 0 ;
	
}

int do_save(int t, char *buff, int bytes, int accepted)
{
	iccp2k i2c ;
	fs_iovec fsiovec[4] ;
	int n_iovec ;
	char sfs_buff[256] ;
	sfs_index sfs ;
	char name[128] ;
	int len ;
	int tot_bytes ;
	int ret ;
	int evb_ix = -1 ;


	ret = 0 ;

	n_iovec = 0 ;
	tot_bytes = 0 ;

	if(evb_chooser.nservers) {
		evb_ix = evb_chooser.chooseIdx(t) ;
	}

	if(accepted && (gl3_accept_cou<100)) {	// local dump first 100
		gbPayload gb ;	// God...
		char gb_buff[256] ;

		timeval tv ;
		gettimeofday(&tv,0) ;

		memset(&gb,0,sizeof(gb)) ;

	
		gb.gbPayloadVersion = GB_PAYLOAD_VERSION ;
		gb.token = t ;
		gb.sec = tv.tv_sec ;
		gb.usec = tv.tv_usec ;
		gb.rtsDetMask |= (1<<L3_SYSTEM) ;
		gb.eventNumber = tokens[t].eventNumber ;


		sprintf(name,"/#%d/EventSummary",tokens[t].eventNumber) ;
		len = sfs.putfileheader(gb_buff,name,sizeof(gb),SFS_ATTR_CD) ;
		fsiovec[n_iovec].filename = 0 ;
		fsiovec[n_iovec].buff = gb_buff ;
		fsiovec[n_iovec].len = len  ;

		tot_bytes += fsiovec[n_iovec].len ;
		n_iovec++ ;

	
		fsiovec[n_iovec].filename = 0 ;
		fsiovec[n_iovec].buff = &gb ;
		fsiovec[n_iovec].len = sizeof(gb) ;

		tot_bytes += fsiovec[n_iovec].len ;
		n_iovec++ ;
	

		sprintf(name,"l3/l3_gtd") ;
		len = sfs.putfileheader(sfs_buff,name,bytes) ;
		fsiovec[n_iovec].filename = 0 ;
		fsiovec[n_iovec].buff = sfs_buff ;
		fsiovec[n_iovec].len = len ;

		tot_bytes += fsiovec[n_iovec].len ;
		n_iovec++ ;



		fsiovec[n_iovec].filename = 0 ;
		fsiovec[n_iovec].buff = buff ;
		fsiovec[n_iovec].len = bytes ;

		tot_bytes += fsiovec[n_iovec].len ;
		n_iovec++ ;



		errno = 0 ;

		ret = sfs_debug.writev(fsiovec,n_iovec) ;
		if(ret != tot_bytes) {
			LOG(ERR,"GL3 %d: error in writev, expect %d, got %d [%s]",gl3_id, tot_bytes, ret, strerror(errno)) ;
			return -1 ;
		}
		else {
			if(evb_ix < 0) {
				LOG(TERR,"GL3 %d: seq %d, T %4d: saved, %d bytes (no EVBs)",gl3_id, tokens[t].eventNumber,t,ret) ;
				return 0 ;	// no evbs in this run... OK...
			}
		}

	}	// accepted

	if(gl3_mode) return 0 ;	// not online so don't connect to evbs!
	if(evb_ix < 0) return 0 ;


	// and now to the EVBs
	n_iovec = 0 ;
	tot_bytes = 0 ;


	i2c.token = t ;
	i2c.srcNode = Mon.rts.node ;
	i2c.srcTask = Mon.rts.task ;
	i2c.dstNode = evtSender[evb_ix].dstNode ;
	i2c.dstTask = evtSender[evb_ix].dstTask ;
	i2c.words = 0 ;

	if(accepted) {
		i2c.cmd = CMD2_DATA ;
	}
	else {
		i2c.cmd = CMD2_L3_RELEASE ;
	}

	fsiovec[n_iovec].filename = 0 ;
	fsiovec[n_iovec].buff = &i2c ;
	fsiovec[n_iovec].len = sizeof(i2c) ;

	tot_bytes += fsiovec[n_iovec].len ;
	n_iovec++ ;


	if(accepted) {
		sprintf(name,"l3_gtd") ;
		len = sfs.putfileheader(sfs_buff,name,bytes) ;
		fsiovec[n_iovec].filename = 0 ;
		fsiovec[n_iovec].buff = sfs_buff ;
		fsiovec[n_iovec].len = len ;

		tot_bytes += fsiovec[n_iovec].len ;
		i2c.words += fsiovec[n_iovec].len ;
		n_iovec++ ;

		fsiovec[n_iovec].filename = 0 ;
		fsiovec[n_iovec].buff = buff ;
		fsiovec[n_iovec].len = bytes ;

		tot_bytes += fsiovec[n_iovec].len ;
		i2c.words += fsiovec[n_iovec].len ;
		n_iovec++ ;

		i2c.words /= 4 ;	// to get to words...
	} 

	errno = 0 ;

	ret = evtSender[evb_ix].sfs.writev(fsiovec,n_iovec) ;
	if(ret != tot_bytes) {
		LOG(ERR,"GL3 %d: EVB %d: error in writev, expect %d, got %d [%s]",gl3_id, evb_ix,tot_bytes, ret, strerror(errno)) ;
		return -1 ;
	}
	else {
		LOG(TERR,"GL3 %d: EVB %d: seq %d, T %4d: %d bytes, accepted %d",gl3_id, evb_ix,tokens[t].eventNumber,t,ret,accepted) ;
	}

	Mon.kbEvb += tot_bytes ;
	
	return ret ;

}

void *checker_thread(void *arg)
{
	time_t now, delta ;
	ic_msg ic ;

	ic.head.daq_cmd = DAQ_TIMEOUT ;

	for(;;) {

	sleep(2) ;

	now = time(0) ;

	for(int t=0;t<4096;t++) {
		int expire ;
		u_int st = tokens[t].status ;

		if(tokens[t].time_in == 0) continue ;

		expire = 0 ;

		delta = now - tokens[t].time_in ;

		if((st==1) && (delta > 10)) {
			expire = 1 ;
		}

		if((st & 0xFFF0) && (delta>2)) {	// some error, but still wait 2 seconds...
			expire = 1 ;
		}

		if(expire) {

			tokens[t].status |= 0x100 ;	// expired

			LOG(NOTE,"GL3 %d: expiring T %4d: delta %d sec, status 0x%02X",gl3_id,t,delta,tokens[t].status) ;

			ic.head.token = t ;
			gl3_que->send(&ic,sizeof(ic),1) ;

		}



	}
		
	}	// forever

	return 0 ;
}


void *monitor_thread(void *arg)
{
	struct sockaddr_in serverAddr ;
	struct timeval ts_old, ts_new ;
	

	errno = 0 ;
	int odesc = socket(AF_INET, SOCK_DGRAM, 0) ;
	
	int socksize = sizeof(struct sockaddr_in) ;
	memset(&serverAddr, 0, socksize) ;

	serverAddr.sin_family = AF_INET ;
	serverAddr.sin_port = htons(RTS_MON_PORT) ;

	memset(&Mon,0,sizeof(Mon)) ;

	Mon.rts.node = GL3_NODES(gl3_id) ;
	Mon.rts.task = GL3_TASK ;
	Mon.rts.couEvtsAll = 0 ;
	Mon.rts.version = RTS_MON_VERSION_X ;
	Mon.rts.size = sizeof(rtsMonStruct) + sizeof(rtsMonL1) ;


        if(((serverAddr.sin_addr.s_addr = inet_addr(RTS_MON_HOST)) == (in_addr_t)-1)) {
		LOG(DBG,"Unknown server %s, trying daqman.daq.bnl.local!",(int)RTS_MON_HOST,0,0,0,0) ;
                errno = 0 ;
                struct hostent *hent = gethostbyname("daqman.daq.bnl.local") ;
                if(hent == NULL) {
                        LOG(WARN,"Can't get daqman.daq.bnl.local IP -- trying localhost",(int)strerror(errno),0,0,0,0) ;

			errno = 0 ;
			hent = gethostbyname("localhost.localdomain") ;

		}
		memcpy(&serverAddr.sin_addr.s_addr,*(hent->h_addr_list),sizeof(serverAddr.sin_addr.s_addr)) ;
	}


	gettimeofday(&ts_old, 0) ;

	double old_kbEvb = 0.0 ;
	double old_kbL3 = 0.0 ;
	u_int old_evts = 0 ;

	for(;;) {



	int busy ;

	busy = 0 ;

#define LOOPS_P_SEC	20
	for(int i=0;i<LOOPS_P_SEC;i++) {
		usleep(1000000/LOOPS_P_SEC) ;

		if(Mon.busy) busy += (100/LOOPS_P_SEC) ;
	}

	Mon.rts.busy = busy ;

	gettimeofday(&ts_new, 0) ;

	Mon.rts.tim = ts_new.tv_sec ;

	double usec = ((double)(ts_new.tv_sec * 1000000.0) + (double)ts_new.tv_usec) ;
	usec -= ((double)(ts_old.tv_sec * 1000000.0) + (double)ts_old.tv_usec) ;


	ts_old = ts_new ;



	if(usec == 0.0) {
		LOG(WARN,"usec 0?") ;
		usec = 1.0 ;
	}


	if(Mon.rts.couEvtsRun < old_evts) old_evts = 0 ;
	Mon.rts.evtsSec = (u_int)(((double)(Mon.rts.couEvtsRun - old_evts)/usec)*1000000.0) ;
	old_evts = Mon.rts.couEvtsRun ;


	if((Mon.kbEvb - old_kbEvb) < 0.0) old_kbEvb = 0.0 ;
	Mon.rts.kbSecEvb = (u_int)(((1000000.0 / 1024.0) * (Mon.kbEvb - old_kbEvb)) / usec) ;
	old_kbEvb = Mon.kbEvb ;


	if((Mon.kbL3 - old_kbL3) < 0.0) old_kbL3 = 0.0 ;
	Mon.rts.kbSecAux = (u_int)((1000000.0 * (Mon.kbL3 - old_kbL3)) / usec / 1024.0) ;
	old_kbL3 = Mon.kbL3 ;


//	Mon.rts.couEvtsIn = atomic_read(&atomic_couEvtsIn) ;

	// HACK -- swap so that DAQ monitoring sees the input rate
	u_int tmp = Mon.rts.kbSecEvb ;
	Mon.rts.kbSecEvb = Mon.rts.kbSecAux ;
	Mon.rts.kbSecAux = tmp ;

	char cbuff[4096] ;
	char *buff = cbuff ;

	memcpy(buff,&Mon.rts,17*4) ;
	buff += 17*4 ;
	memcpy(buff,&Mon.trg,sizeof(rtsMonL1)) ;
	buff += sizeof(rtsMonL1) ;

	errno = 0 ;
#if 0
	int ret = sendto(odesc,cbuff, buff-cbuff,
			 0, (struct sockaddr *) &serverAddr, socksize) ;
		
	if(ret <= 0) {
		LOG(ERR,"sendto: %s",strerror(errno)) ;
	}
#endif


	}

	return 0 ;
}
