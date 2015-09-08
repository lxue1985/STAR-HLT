#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <string.h>
#include<sys/time.h>
#include<unistd.h>

#include <daqFormats.h>
#include <rts.h>
#include "rtsLog.h"
#include <DAQ_READER/daqReader.h>
#include <DAQ_READER/daq_det.h>
#include <DAQ_READER/daq_dta.h>
#include <DAQ_L3/daq_l3.h>
#include <DAQ_TRG/daq_trg.h>
#include "eventTrackerLib.hh"
#include "online_tracking_collector.h"
#include "online_tracking_sector.h"
#include "online_tracking_TpcHitMap.h"
#include "l3EmcCalibration.h"

#include <rtsLog.h>
#include "gl3Histo.h"


//#include "lib.cuh"
#include "clib_kernel.h"
#include "online_GPU_tracking.h"

#define  maxTracks 1000
void Hcalc_hash(int number,float *fpara,float *fhit,int * volumeC,int * rowC);
void Hclear(int * dataC,int number);
void HSclass_hash(int number,float *fpara,float *fhit, int * volume,int * row);
void Htracking(int number,float *para,float *fhit,int * volumeC,int * rowC,int * trackC,float *track);
void setParameters (clib_FtfPara & para) ;



unsigned int utime(unsigned int* dp,int i);
void udiff(unsigned int* dp,int i);

clib_FtfPara * getParaV;


#include "gl3Track.h"

#ifdef OLD_DAQ_READER
#include <evpReader.hh>

#else /* OLD_DAQ_READER */
#include <DAQ_READER/daqReader.h>

#endif /* OLD_DAQ_READER */
#include "FtfSl3.h"

#ifndef OLD_DAQ_READER
#include <DAQ_READER/daq_dta.h>
#include <rtsLog.h>	// for my LOG() call
// this needs to be always included
// only the detectors we will use need to be included
// for their structure definitions...
#include <DAQ_BSMD/daq_bsmd.h>
#include <DAQ_BTOW/daq_btow.h>
#include <DAQ_EMC/daq_emc.h>
#include <DAQ_ESMD/daq_esmd.h>
#include <DAQ_ETOW/daq_etow.h>
#include <DAQ_FPD/daq_fpd.h>
#include <DAQ_FTP/daq_ftp.h>
#include <DAQ_L3/daq_l3.h>
#include <DAQ_PMD/daq_pmd.h>
#include <DAQ_PP2PP/daq_pp2pp.h>
#include <DAQ_RIC/daq_ric.h>
#include <DAQ_SC/daq_sc.h>
#include <DAQ_SSD/daq_ssd.h>
#include <DAQ_SVT/daq_svt.h>
#include <DAQ_TOF/daq_tof.h>
#include <DAQ_TPC/daq_tpc.h>
#include <DAQ_TPX/daq_tpx.h>
#include <DAQ_TRG/daq_trg.h>


#endif /* OLD_DAQ_READER */
int what=GL3_READ_TPC_TRACKS;
inline float fswap(float swapped)

{

    unsigned int* uintptr = (unsigned int*) &swapped;
    unsigned int uintvar = l2h32(*uintptr);
    float* floatvar = (float*)&uintvar;
    return *floatvar;
    
}


#ifdef UNIX_LITTLE_ENDIAN
#define l2hfloat(x) (x)
#define b2hfloat(x) (fswap(x))

#else
#define l2hfloat(x) (fswap(x))
#define b2hfloat(x) (x)

#endif
// control....
char g_fn[255];    
// filename

uint g_seq = 0;        
// only write this event if != 0

float g_bfield = 1000;
uint g_ptracks = 0;
uint g_nftracks = 0;
uint g_nctracks = 0;
uint g_pause = 0;
uint g_vertex = 0;
int copyl3_t(l3_t &l3, L3_P *l3p);
void printL3Info(l3_t& l3)

{

    //printf("%d tracks, %d clusters: Vertex = (%f, %f %f)\n",l3.tracks_num, l3.cluster_num, l3.xVertex,l3.yVertex,l3.zVertex);
    if(g_ptracks) 
    {

        for(u_int i=0;i<l3.tracks_num;i++) 
        {

            global_track *track = &l3.track[i];
            //printf("%5d: pt=%5.3f z0=%7.2f q=%2d nHits=%2d ndedx=%2d ",
            //  i, track->pt, track->z0, track->q,
            //  track->nHits, track->ndedx);
            //printf("flag=0x%04x iRow=%2d oRow=%2d\n",
            //   track->flag, track->innerMostRow,
            //  track->outerMostRow);
            
        }

        
    }

    
}


unsigned int utime(unsigned int* dp,int i){
struct timeval tv;
struct timezone tz;
gettimeofday (&tv , &tz);
dp[i]=(tv.tv_sec%1000)*1000000+tv.tv_usec;
return dp[i];
}
void udiff(unsigned int* dp,int i){

int l;
for(l=1;l<=i;l++){
printf("time diff from %d to %d is  %d micro-second\n",l-1,l,dp[l]-dp[l-1] );
}
}

void parseArgs(int argc, char *argv[])

{

    g_fn[0] = 0;
    for(int i=1;i<argc;i++) 
    {

        if(argv[i][0] != '-') 
        {   
                        // should be filename...

            if(g_fn[0] != 0) goto badargs;
            strcpy(g_fn, argv[i]);
            
        }

        else if(strcmp(argv[i], "-event") == 0) 
        {

            i++;
            g_seq = atoi(argv[i]);
            ////printf("g_seq = %d\n",g_seq);
            
        }

        else if(strcmp(argv[i], "-B") == 0) 
        {

            i++;
            g_bfield = atof(argv[i]);
            //printf("g_bfield = %f\n",g_bfield);
            
        }

        else if (strcmp(argv[i], "-p") == 0) 
        {

            g_ptracks = 1;
            
        }

        else if (strcmp(argv[i], "-nft") == 0) 
        {

            g_nftracks = 1;
            
        }

        else if (strcmp(argv[i], "-nct") == 0) 
        {

            g_nctracks = 1;
            
        }

        else if (strcmp(argv[i], "-pause") == 0) 
        {

            g_pause = 1;
            
        }

        else if (strcmp(argv[i], "-vtx") == 0) 
        {

            rtsLogLevel(CRIT);
            g_vertex = 1;
            
        }

        else 
        {

            goto badargs;
            
        }

        
    }

    if(g_fn[0] == 0) goto badargs;
    return;
    badargs:
    //printf("eventTracker <-event #> <-B> <-p> <-nft> <-nct> filename\n");
    //printf("\t-event    -> event number\n");
    //printf("\t-B        -> mag field\n");
    //printf("\t-p        -> print track params\n");
    //printf("\t-nft      -> don't read tracks from file\n");
    //printf("\t-nct      -> don't calculate tracks\n");
    //printf("\t-pause    -> pause before and after\n");
    //printf("\t-vtx      -> produce output for vtx scan\n");
    exit(0);
    
}

float getfield(daqReader *rdr);
float getfield(daqReader *rdr)
{

    float bField=g_bfield;
    if(bField == 1000) 
    {  
                // try to read from file

        bField =-0.5;       
                // use default

        daq_dta *dd;
        dd = rdr->det("sc")->get("legacy");
        if(dd) 
        {

            dd->iterate();
            sc_t *sc = (sc_t *)dd->Void;
            if(sc->valid) bField = sc->mag_field;
            
        }

        
    }       

    if(fabs(bField) < .1) bField = .1;
    return bField;
    
}

tpc_t *getTPC(daqReader *rdr,int i);
tpc_t *getTPC(daqReader *rdr,int i)
{

    tpc_t *  sxm_pTPC;
    daq_dta *dd;
    dd = rdr->det("tpx")->get("legacy",i+1);
    if(dd) 
    {

        LOG(NOTE, "There is tpx data...");
        dd->iterate();
        sxm_pTPC = (tpc_t *)dd->Void;
        
    }

    else 
    {

        LOG(NOTE, "No tpx data for sector %d check for TPC",i);
        dd = rdr->det("tpc")->get("legacy",i+1);
        if(dd) 
        {

            dd->iterate();
            sxm_pTPC = (tpc_t *)dd->Void;
            int cl_found = 0;
            for(int pr=0;pr<45;pr++) 
            {

                cl_found += sxm_pTPC->cl_counts[pr];
                
            }

            LOG(NOTE, "Found tpc data for sector %d... %d clusters found",i,cl_found);
            
        }

        else 
        {

            sxm_pTPC = NULL;
            
        }

        
    }

    if(!sxm_pTPC) 
    {

        LOG(WARN, "No data for TPC sector %d",i+1,0,0,0,0);
        return NULL;
        
    }

    //LOG(NOTE, "Tpc reader done");
    if(!sxm_pTPC->has_clusters) 
    {

        LOG(WARN, "TPC sector %d has no clusters",i);
        return NULL;
        
    }

    // read in clusters...
    if(what & GL3_READ_TPC_CLUSTERS) 
    {

        LOG(DBG, "Reading clusters");
        // gl3->readClustersFromEvpReader(i+1);
        int nnn=0;
        for(int i=0;i<45;i++) 
        {

            nnn += sxm_pTPC->cl_counts[i];
            
        }

        LOG(DBG, "clusters done %d",nnn);
        
    }

    return sxm_pTPC;
    
}


int load_TPC(daqReader *rdr, int &sector,online_GPU_tracking *co_tracker){

	int i=0;
    daq_dta *dd;
    tpc_t *Tracking_pTPC;

	for(i=0;i<co_tracker->dup;i++){
		if(sector==0){

			
        char *mem = rdr->get(0,EVP_TYPE_PHYS);
        if(!mem) 
        {

            if(rdr->status == EVP_STAT_EOR) 
            {

                if(!g_vertex) 
                {

                    //printf("End of run...\n");
                    
                }

                break;
                
            }

            else 
            {

                //printf("Error reading an event\n");
                break;
                
            }

            
        }

        // We have an event of some kind...
        //
        // do filtering...
        //
        if(g_seq != 0) 
        {                
                        // event number

            if(rdr->seq != g_seq) continue;
            
        }

        if(rdr->token == 0) continue;
        if(!g_vertex) 
        { 

            //printf("**** Event %d (seq = %d): %d bytes, token %d, triggers = 0x%x\n",
            //  rdr->event_number, rdr->seq, rdr->bytes, rdr->token, rdr->daqbits);
            
        }

        DATAP *datap = (DATAP *)mem;
        if(!datap) 
        {

            //printf("Error reading datap:\n");
            break;
            
        }

        // First Use old L3 reader to read L3 from datafile if its there...
        if(!g_nftracks) 
        {

            dd = rdr->det("l3")->get("legacy");
            //LOG("JEFF", "blih");
            if(!dd) 
            {

                ////printf("No L3 banks in data file %d\n",ret);
                
            }

            else 
            {

                dd->iterate();
                l3_t *pL3 = (l3_t *)dd->Void;
                //printf("This comes from the datafile L3 banks...------ len=%d\n",ret);
                printL3Info(*pL3);
                //printf("End Datafile L3 banks-------------------------\n");
                
            }

            
        }		
		
		};


            Tracking_pTPC=getTPC(rdr,sector);
            int size;

			int ipsix=i;

 
			//tracker->nHits = 0;
			co_tracker->ptracker[ipsix]->nHits=0;
			//tracker->setTrackingAngles(i+1);   
			co_tracker->ptracker[ipsix]->setTrackingAngles(sector+1);   
			//tracker->sector_ID=i;
			co_tracker->ptracker[ipsix]->sector_ID=sector;
			int nhits= co_tracker->ptracker[ipsix]->Tracking_load_cluster(Tracking_pTPC);
			//cout<<"nhitsnhitsnhits  "<<nhits<<endl;
			co_tracker->ptracker[ipsix]->setPointers();
			co_tracker->load_hit_to_buffer(ipsix);
			co_tracker->set_bfield(getfield(rdr),ipsix); 
			sector++;
			if(sector==24)sector=0;
		

	
	}

return i;
}


// Dumps tracks....
int main(int argc, char *argv[])

{

//cudaInit(argc, argv);

unsigned int ts[100];





//clib_FtfPara * m_Device_para;
clib_FtfPara * m_Host_para;



int max_track=2000;

//small_FtfHit *smhit=new small_FtfHit[max_hit];




int Di;
//for(Di=0;Di<dup;Di++){
//map[Di]=0;
//manager[Di]=0;
//}






//small_FtfTrack *smtrack=new small_FtfTrack[max_track];
//small_hash* hit_hash=new small_hash[41*11*46];

//small_track_manager track_manager;


	clib_FtfPara bm_Host_para;
	m_Host_para=&bm_Host_para;
    //allocateArray((void**)&m_Device_para, sizeof(clib_FtfPara));



getParaV=&bm_Host_para;



   online_tracking_TpcHitMap* tpcHitMap = new online_tracking_TpcHitMap();
    tpcHitMap->loadMap();

	l3EmcCalibration *bemcCalib = new l3EmcCalibration(4800);
    bemcCalib->loadTextMap("/home/sunxm/backup/usr/emcmap.txt");
    l3EmcCalibration * eemcCalib = new l3EmcCalibration(720);


    int ret = 0;
    l3_t *l3_datafile;
    l3_t l3_legacy;
    rtsLogOutput(RTS_LOG_STDERR);
    parseArgs(argc, argv);
    rtsLogLevel(WARN);
    // //printf("sizeof tpc %d\n",sizeof(tpc));
    daqReader *rdr = new daqReader(g_fn);
    if(!rdr) 
    {

        //printf("Error getting daqReader\n");
        return 0;
        
    }

    // Buffer for event storage...
    L3_P *l3p = (L3_P *)malloc(szL3_max);
    EventTracker *evtTracker;
    if(g_bfield == 1000)
    evtTracker = new EventTracker();
    else
    evtTracker = new EventTracker(g_bfield);
    char tmp[200];
    if(g_pause) 
    {

        //printf("Enter something: ");
        scanf("%s", tmp);
        
    }

    int eve=0;


        online_tracking_collector *eventbuilder=new online_tracking_collector(NULL,bemcCalib,eemcCalib);

        //online_tracking_sector *tracker = new online_tracking_sector(-0.5,"/RTS/conf/L3/map.bin",0);
		online_GPU_tracking *co_tracker=new online_GPU_tracking();
		online_GPU_tracking *ALT_co_tracker=new online_GPU_tracking();

int trackers=6;
		co_tracker->setup_memory(trackers);
		ALT_co_tracker->setup_memory(trackers);

		co_tracker->set_tpcmap(tpcHitMap);
		ALT_co_tracker->set_tpcmap(tpcHitMap);

		//tracker->setTpcHitMap(tpcHitMap);	
	        L3_SECTP *sectp = (L3_SECTP *)malloc(szSECP_max);
	


        float bField=-0.5;
        int what=GL3_READ_TPC_TRACKS;
      
		
		
        eventbuilder->resetEvent();





        // need temporary track memory...

        int i;
       
float max=-10000;
float min=10000;
int sector=0;

////////////////////////////////////////////prepare data for GPU
				load_TPC(rdr, sector,co_tracker);
				co_tracker->set_balance(30);
				co_tracker->fire_GPU(30);

				
//////////////////////////////////////////// GPU is working ,CPU is idle, so CPU prepare data for next GPU run
				load_TPC(rdr, sector,ALT_co_tracker);
				ALT_co_tracker->set_balance(30);


//////////////////////////////////////////// copy the first GPU result, start the prepared next GPU run
				 ret=co_tracker->copy_GPU_back();
				ALT_co_tracker->fire_GPU(30);

//////////////////////////////////////////// GPU is working ,CPU is idle, so CPU finalize GPU data and prepare data for co_tracker
				
				for(Di=0;Di<co_tracker->dup;Di++){
				co_tracker->build_track(Di);
				co_tracker->ptracker[Di]->fillTracks(szSECP_max, (char *)sectp, 0);
								//tracker->Tracking_track_sector((char *)sectp,size);

				int n = eventbuilder->Tracking_readSectorTracks((char *)sectp);
				}

				load_TPC(rdr, sector,co_tracker);
				co_tracker->set_balance(30);

////////////////////////////////////////////// CPU finalize GPU data for the second GPU run, 
				ret=ALT_co_tracker->copy_GPU_back();
				co_tracker->fire_GPU(30);

//////////////////////////////////////////// GPU is working ,CPU is idle, so CPU finalize GPU data and prepare data for ALT_co_tracker
				
				for(Di=0;Di<ALT_co_tracker->dup;Di++){
				ALT_co_tracker->build_track(Di);
				ALT_co_tracker->ptracker[Di]->fillTracks(szSECP_max, (char *)sectp, 0);
								//tracker->Tracking_track_sector((char *)sectp,size);

				int n = eventbuilder->Tracking_readSectorTracks((char *)sectp);
				}

				load_TPC(rdr, sector,ALT_co_tracker);
				ALT_co_tracker->set_balance(30); 

/////////////////////////////////////////////

        eventbuilder->finalizeReconstruction();
        printf("getNTracks %d \n",eventbuilder->getNTracks());

		
		//int ntra=eventbuilder->getNTracks();


        
    

        free(sectp);
        //delete tracker;
        delete eventbuilder;
		
		
    if(g_pause) 
    {

        //printf("Enter something: ");
        scanf("%s", tmp);
        
    }

    free(l3p);



	


}

// Copy from type L3_P --> l3_t
//
int copyl3_t(l3_t &l3, L3_P *l3p)

{

    int len = l3p->bh.length;
    l3.max_channels = 1280000 ;
    l3.channels = 0 ;
    l3.mode = 1 ;
    // Tonko, zero this out and make non-sensical in case of further
    // problems down the road
    l3.tracks_num = 0 ;
    l3.cluster_num = 0 ;
    l3.xVertex = -1000.0 ;
    l3.yVertex = -1000.0 ;
    l3.zVertex = -1000.0 ;
    LOG(DBG,"L3_P bytes %d",len,0,0,0) ;
    if(checkBank(l3p->bh.bank_type,"L3_P") != 0) 
    {

        return -1 ;
        
    }

    if(l3p->tracks.len && l3p->tracks.off)
    {

        struct L3_GTD* l3gtd = 
        (struct L3_GTD*)((char*)l3p + l2h32(l3p->tracks.off)*4) ;
        // Tonko, sanity check
        if(checkBank(l3gtd->bh.bank_type,"L3_GTD") < 0) 
        {

            return -1 ;
            
        }

        LOG("JEFF", "l3gtd->nTracks=%d nHits=%d", l2h32(l3gtd->nTracks), l2h32(l3gtd->nHits));
        l3.tracks_num = l2h32(l3gtd->nTracks);
        l3.cluster_num = l2h32(l3gtd->nHits);
        l3.xVertex = l2hfloat(l3gtd->xVert);
        l3.yVertex = l2hfloat(l3gtd->yVert);
        l3.zVertex = l2hfloat(l3gtd->zVert);
        // Tonko, sanity check
        if(l3.tracks_num >= L3_MAX_NR_TRACKS) 
        {

            LOG(ERR,"L3 track number %d > %d!",l3.tracks_num,L3_MAX_NR_TRACKS ,0,0,0) ;
            return -1 ;
            
        }

        for (unsigned int i=0; i<l3.tracks_num; i++) 
        {

            global_track *tr = &(l3gtd->track[i]);
            l3.track[i].id = l2h32(tr->id);
            
            #ifndef UNIX_LITTLE_ENIDAN
            l3.track[i].flag         = tr->flag;
            l3.track[i].innerMostRow = tr->innerMostRow;
            l3.track[i].outerMostRow = tr->outerMostRow;
            l3.track[i].nHits = tr->nHits;
            l3.track[i].ndedx = tr->ndedx;
            l3.track[i].q     = tr->q;
            
            #else
            l3.track[i].flag = ( ((unsigned short)tr->innerMostRow) |
            ((unsigned short)tr->outerMostRow)<<8);
            l3.track[i].innerMostRow = (char)( tr->flag & 0x00ff );
            l3.track[i].outerMostRow = (char)((tr->flag & 0xff00)>>8);
            l3.track[i].nHits    = (unsigned char)tr->q;
            l3.track[i].reserved = (char)tr->ndedx;
            l3.track[i].ndedx    = (unsigned char)tr->reserved;
            l3.track[i].q        = (char)tr->nHits;
            
            #endif
            l3.track[i].chi2[0] = l2hfloat(tr->chi2[0]);
            l3.track[i].chi2[1] = l2hfloat(tr->chi2[1]);
            l3.track[i].dedx    = l2hfloat(tr->dedx);
            l3.track[i].pt      = l2hfloat(tr->pt);
            l3.track[i].phi0    = l2hfloat(tr->phi0);
            l3.track[i].psi     = l2hfloat(tr->psi);
            l3.track[i].r0      = l2hfloat(tr->r0);
            l3.track[i].tanl    = l2hfloat(tr->tanl);
            l3.track[i].z0      = l2hfloat(tr->z0);
            l3.track[i].length  = l2hfloat(tr->length);
            l3.track[i].dpt     = l2hfloat(tr->dpt);
            l3.track[i].dpsi    = l2hfloat(tr->dpsi);
            l3.track[i].dz0     = l2hfloat(tr->dz0);
            l3.track[i].dtanl   = l2hfloat(tr->dtanl);
            
        }

        
    }

    
    #ifdef SHOW_DEBUG_INFO
    
    #ifdef UNIX_LITTLE_ENDIAN
    //printf("Running on LITTLE endian machine\n");
    
    #else
    //printf("Running on BIG endian machine\n");
    
    #endif
    //printf("\nVertex: (%6.2f/%6.2f/%6.2f)\n",
    l3.xVertex, l3.yVertex, l3.zVertex);
    //printf("Tracks: %5d   Clusters %7d\n",
    l3.tracks_num, l3.cluster_num);
    for (unsigned int i=0; i<l3.tracks_num; i++) 
    {

        //printf("%5d: pt=%5.3f z0=%7.2f q=%2d nHits=%2d ndedx=%2d ",
        i, l3.track[i].pt, l3.track[i].z0, l3.track[i].q, 
        l3.track[i].nHits, l3.track[i].ndedx); 
        //printf("flag=0x%04x iRow=%2d oRow=%2d\n",
        l3.track[i].flag, l3.track[i].innerMostRow, 
        l3.track[i].outerMostRow); 
        
    }

    
    #endif
    return len ;
    
}



void Hcalc_hash(int number,float *fpara,float *fhit,int * volumeC,int * rowC){



// int index = blockIdx.x*blockDim.x + threadIdx.x;




clib_FtfContainer * vc=(clib_FtfContainer *)volumeC;
clib_FtfContainer * rc=(clib_FtfContainer *)rowC;
clib_FtfHit *fh=(clib_FtfHit *)fhit;
clib_FtfPara *para=(clib_FtfPara *)fpara;


int index;
for(index=0;index<number;index++){

float thisHit_x=fh[index].x;
float thisHit_y=fh[index].y;
float thisHit_z=fh[index].z;
float thisHit_row=fh[index].row;


 double       r2            = thisHit_x * thisHit_x + thisHit_y * thisHit_y ;
 float       r             = (float)sqrt ( r2 ) ;
 float       phi           = (float)atan2(thisHit_y,thisHit_x) + para->phiShift ;
        if ( phi < 0 ) phi = phi + twoPi ;

 float       eta           = (float)seta(r,thisHit_z) ;
 float       phiSlice=(para->phiMax-para->phiMin)/para->nPhi;
 int phiIndex = (int)( (phi-para->phiMin)/phiSlice + 1.);
 float etaSlice =(para->etaMax-para->etaMin)/para->nEta;
 int etaIndex = (int)((eta - para->etaMin)/etaSlice + 1.);
 
 

 
 fh[index].r=r;
 fh[index].eta=eta;
 fh[index].phi=phi;
 fh[index].track=0;


 fh[index].phiIndex=phiIndex;
 fh[index].etaIndex=etaIndex;
 fh[index].nextVolumeHit=0;
 fh[index].nextRowHit=0;
 fh[index].nextTrackHit=0;


 
 float        x            = thisHit_x - para->xVertex ;
 float       y            = thisHit_y - para->yVertex ;
        r2           = x * x + y * y ;
 double       invR2        = 1.F / r2 ;
        fh[index].xp    =     x * invR2 ;
        fh[index].yp    =   - y * invR2 ;
		fh[index].wxy   =   r2 * r2 /  ( square(para->xyErrorScale)
			* ( square(fh[index].dx) + square(fh[index].dy) ) ) ;


		fh[index].s=0;
		fh[index].wz=(double)(1./ square ( para->szErrorScale * fh[index].dz ));




}


}
void Hclear(int * dataC,int number){




clib_FtfContainer * fc=(clib_FtfContainer *)dataC;

int index;
for(index=0;index<number;index++){
fc[index].first=NULL;
fc[index].last=NULL;
}



}
void HSclass_hash(int number,float *fpara,float *fhit, int * volume,int * row){

clib_FtfContainer * volumeC=(clib_FtfContainer *)volume;
clib_FtfContainer * rowC=(clib_FtfContainer *)row;
clib_FtfHit *fh=(clib_FtfHit *)fhit;
clib_FtfPara *para=(clib_FtfPara *)fpara;

     int nPhiPlusOne    = para->nPhi + 1 ;
    int nEtaPlusOne    = para->nEta + 1 ;
    int nPhiEtaPlusOne = nPhiPlusOne * nEtaPlusOne ;
 

 
 


int il;

for(il=0;il<number;il++){


// if(ihit[iitem*il+0]>df_nPhi||ihit[iitem*il+0]<1) continue;
// if(ihit[iitem*il+1]>df_nEta||ihit[iitem*il+1]<1) continue;

clib_FtfHit *thisHit=&fh[il];

 if(thisHit->phiIndex>para->nPhi||thisHit->phiIndex<1) continue;
 if(thisHit->etaIndex>para->nEta||thisHit->etaIndex<1) continue;


int localRow = int(thisHit->row) - para->rowInnerMost+1 ;

int	    volumeIndex = localRow  * nPhiEtaPlusOne + 
        thisHit->phiIndex * nEtaPlusOne + thisHit->etaIndex ;

		
	
        if (volumeC[volumeIndex].first == 0 ) 
        volumeC[volumeIndex].first = (void *)thisHit ;
        else{
        ((clib_FtfHit *)(volumeC[volumeIndex].last))->nextVolumeHit = thisHit ;
		}
        volumeC[volumeIndex].last = (void *)thisHit ;
        /*-------------------------------------------------------------------------
        Set row pointers
        -------------------------------------------------------------------------*/
        if ( rowC[localRow].first == NULL )
        
        {

            rowC [localRow].first = (void *)thisHit ;
            
        }

        else{
        ((clib_FtfHit *)(rowC[localRow].last))->nextRowHit = thisHit ;
		}
        rowC[localRow].last = (void *)thisHit ;

}



}
void Htracking(int number,float *fpara,float *fhit,int * ivolumeC,int * irowC,int * itrackC,float *ftrack){




clib_FtfContainer * volumeC=(clib_FtfContainer *)ivolumeC;
clib_FtfContainer * rowC=(clib_FtfContainer *)irowC;
clib_FtfContainer * trackC=(clib_FtfContainer *)itrackC;
clib_FtfTrack * track=(clib_FtfTrack *)ftrack;
clib_FtfHit *hit=(clib_FtfHit *)fhit;
clib_FtfPara *para=(clib_FtfPara *)fpara;

//clib_FtfPara bpara;
//clib_FtfPara *para=&bpara;

//
//para->setDefaults();
//setParameters(*para);


int nTracks=0;



  para->init = 0 ;
//----------------------------------------------------------------------------
//     Allocate volumes 
//---------------------------------------------------------------------------*/
   para->nRowsPlusOne = ( para->rowOuterMost - para->rowInnerMost ) / para->modRow + 2 ;
   if ( para->nRowsPlusOne < 1 ) {
   	
     return  ;
   }
   para->nPhiPlusOne    = para->nPhi + 1 ;
   para->nEtaPlusOne    = para->nEta + 1 ;
   para->nPhiEtaPlusOne = para->nPhiPlusOne * para->nEtaPlusOne ;
   if ( para->mergePrimaries ) {
      para->nPhiTrackPlusOne = para->nPhiTrack + 1 ;
      para->nEtaTrackPlusOne = para->nEtaTrack + 1 ;
   }





    //
    //     Set conformal coordinates if we are working with primaries
    //
    int nHitsSegment   = (short)para->nHitsForSegment;  
    if ( para->primaries ) 
    
    {

        
        para->minHitsForFit = 1 ;
        para->nHitsForSegment = max(2,nHitsSegment);
        
    }

    else 
    
    {

        para->minHitsForFit = 2 ;
        para->nHitsForSegment = max(3,nHitsSegment);
        
    }



		
		
    for ( int ir = para->nRowsPlusOne - 1 ; ir>=para->minHitsPerTrack ; ir--) 
    
	{

	
			
		
        if ( rowC[ir].first &&  (((clib_FtfHit *)rowC[ir].first)->row) < para->rowEnd ) 
        break ;

        for ( clib_FtfHit *firstHit = (clib_FtfHit *)rowC[ir].first ;
        firstHit != 0 ;
        firstHit = (clib_FtfHit *)(firstHit->nextRowHit) ) 
        
		{
		
			
			
		   
            //
            //     Check hit was not used before
            //
            if ( firstHit->track != 0  ) continue ;
            //
            //     One more track
            //
            nTracks++ ;
			            //
            //
            if ( nTracks > maxTracks )
            
            {

                
                nTracks = maxTracks  ;
                return  ;
                
            }

            //
            //     Initialize variables before going into track hit loop
            //
            clib_FtfTrack *thisTrack = &track[nTracks-1];

            thisTrack->para     = para ;

            thisTrack->id       = nTracks ;
            thisTrack->firstHit = thisTrack->lastHit = firstHit ;
            thisTrack->innerMostRow = thisTrack->outerMostRow = firstHit->row ;
 
            thisTrack->xRefHit  = firstHit->x ;
            thisTrack->yRefHit  = firstHit->y ;
            thisTrack->xLastHit = firstHit->x ;
            thisTrack->yLastHit = firstHit->y ;
            
            //
            //              Set fit parameters to zero
            //
            thisTrack->reset ( ) ;
      
			//
            //      Go into hit looking loop
            //
            if ( thisTrack->buildTrack ( firstHit, volumeC ) ) 
            
            {
				

                //
                //    Merge Tracks if requested
                //
                if ( para->primaries &&
                para->mergePrimaries == 1 &&
                para->fillTracks &&
                thisTrack->mergePrimary( trackC )  ) 
                
                {
					

                    nTracks-- ;
                    thisTrack->deleteCandidate ( ) ;
                    
                }

                
            }

            else
            
            {

                //
                //      If track was not built delete candidate
                //
                thisTrack->deleteCandidate ( ) ;
                nTracks-- ;
                
            }

            //
            //       End loop over hits inside row
            //
            
       


             
        }

        
    }
	
	
	
	

    cout<<"HnTracks   "<< nTracks<<endl ;  



// if(rowC[30].first!=NULL){
// clib_FtfHit *firstHit = (clib_FtfHit *)rowC[30].first;
// fpara[0]=firstHit->x;
// }else{
// fpara[0]=0;
// }



}

void setParameters (clib_FtfPara & para) 

{

    // FtfPara* para = &(para) ;
    para.hitChi2Cut        =  50.  ;
    para.goodHitChi2       =  20.  ;
    para.trackChi2Cut      =  10.  ;
    para.maxChi2Primary    = 50.   ; 
    // track with vertex if refit and minPtHelix is reasonable.
    para.segmentRowSearchRange = 2   ;
    para.trackRowSearchRange = 3    ;
    para.dphi              = 0.08F  ;
    para.deta              = 0.08F  ;
    para.dphiMerge        = 0.01F  ;
    para.detaMerge        = 0.02F  ;
    para.etaMinTrack      = -2.2F  ;
    para.etaMaxTrack      =  2.2F  ;
    para.dEdx             =  1     ;
    // get Errors = 1 makes FPE's
    para.getErrors        =  0     ;
    para.goBackwards      =  1     ;
    para.goodDistance     =  5.F   ;
    para.mergePrimaries   =  0     ;
	para.primaries   =  1     ;
    para.maxDistanceSegment = 50.F ;
    para.minHitsPerTrack  = 5      ;
    para.nHitsForSegment  = 2      ;
    para.nEta             = 40     ;
    para.nEtaTrack        = 40     ;
    para.nPhi             = 10     ;
    para.nPhiTrack        = 40     ;
    para.nPrimaryPasses   = 1      ;
    para.nSecondaryPasses = 0      ;
    para.xyErrorScale     = 1.0F   ;
    para.szErrorScale     = 1.0F   ;
    para.phiClosed        = 0      ;
    //   para.ptMinHelixFit    = 100.F  ;
    para.ptMinHelixFit    = 0.1  ; 
    // enables refit
    para.rVertex          = 0.F    ;
    para.xVertex          = 0.F    ;
    para.yVertex          = 0.F    ;
    para.zVertex          = 0.F    ;
    para.dxVertex         = 0.05F ;
    para.dyVertex         = 0.05F ;
    para.phiVertex        = 0.F    ;
    para.zMax             = 205. ;
    
    
}