#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <string.h>
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
#include "gl3TriggerDecider.h"

#include <rtsLog.h>
#include "gl3Histo.h"

#include <L3_SUPPORT/l3_support.h>





#include "mTrackEvent.h"


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
#include <DAQ_PP2PP/daq_pp2pp.h>
#include <DAQ_RIC/daq_ric.h>
#include <DAQ_SC/daq_sc.h>
#include <DAQ_SSD/daq_ssd.h>
#include <DAQ_SVT/daq_svt.h>
#include <DAQ_TOF/daq_tof.h>
#include <DAQ_TPC/daq_tpc.h>
#include <DAQ_TPX/daq_tpx.h>
#include <DAQ_TRG/daq_trg.h>


//====================QIU Hao add to do QA
#include <stdio.h>
//====================QIU Hao add to do QA



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
        for(int j=0;j<45;j++) 
        {

            nnn += sxm_pTPC->cl_counts[j];
            
        }

        LOG(DBG, "clusters done %d",nnn);
        
    }

    return sxm_pTPC;
    
}

// Dumps tracks....
int main(int argc, char *argv[])

{  
  daq_dta *dd;
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
      
      printf("Error getting daqReader\n");
      return 0;
      
    }
  
  // Buffer for event storage...
  L3_P *l3p = (L3_P *)malloc(szL3_max);
  char tmp[200];
  if(g_pause) 
    {

      printf("Enter something: ");
      scanf("%s", tmp);
      
    }
  
  int eve=0;
  
  online_tracking_sector *tracker[24];
  for(int i=0; i<24; i++)
    {
      char mapName[256];
      sprintf(mapName, "tpcHitMap_sector%d.bin", i+1); 
      
      tracker[i] = new online_tracking_sector(-0.5,i,mapName);
    }
  
  online_tracking_collector *eventbuilder=new online_tracking_collector(NULL,NULL,NULL);
    
  L3_SECTP *sectp = (L3_SECTP *)malloc(szSECP_max);
    
  
  eventbuilder->triggerDecider->setQA("../../../../HLTQA.data");
  
  
  
  for(;;){
    
    //    printf("event %d \n", eve);
    eve++;
    if(eve>500) break;
    char *mem = rdr->get(0,EVP_TYPE_PHYS);
    if(!mem) 
      {
	
	if(rdr->status == EVP_STAT_EOR) 
	  {
	    
	    if(!g_vertex) 
	      {
		
		printf("End of run...\n");
                    
	      }
	    
	    break;
	    
	  }
	
	else 
	  {
	    
                printf("Error reading an event\n");
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
    
    const int trigger[4] = {0x10, 0x20, 0x80, 0x100};
    int triggered = 0;
    for(int i=0; i<4; i++)
      {
	if(rdr->daqbits & trigger[i])
	  triggered = 1;
      }
    //    if(!triggered) continue;
    
    if(!g_vertex) 
      { 
	
	printf("**** Event %d (seq = %d): %d bytes, token %d, triggers = 0x%x\n",
	  rdr->event_number, rdr->seq, rdr->bytes, rdr->token, rdr->daqbits);
	
      }
    
    DATAP *datap = (DATAP *)mem;
    if(!datap) 
      {
	    
	printf("Error reading datap:\n");
	break;
	
      }
	
    // First Use old L3 reader to read L3 from datafile if its there...
    if(!g_nftracks) 
      {
	
	dd = rdr->det("l3")->get("legacy");
	//	LOG("JEFF", "blih");
	if(!dd) 
	  {
	    
	    //	    printf("No L3 banks in data file %d\n",ret);
	    
	  }
	
	else 
	  {
	    
	    dd->iterate();
	    l3_t *pL3 = (l3_t *)dd->Void;
	    //	    printf("This comes from the datafile L3 banks...------ len=%d\n",ret);
	    //	    printL3Info(*pL3);
	    //	    printf("End Datafile L3 banks-------------------------\n");
	    
	  }
	    
	
      }
    
    //////////////////////////////////////////////////////////////////////////////////////
    tpc_t *Tracking_pTPC;
    float bField=-0.5;
    int what=GL3_READ_TPC_TRACKS;
    daq_dta *dd;
    
    
    eventbuilder->resetEvent();
    bField=getfield(rdr);
    eventbuilder->setBField(bField);    
    
    dd = rdr->det("trg")->get("raw") ;
    if(dd && dd->iterate()) {
      gl3_trg_send_t gl3_data ;	


      gl3_trg_parse((char *)dd->Void, dd->ncontent, &gl3_data) ;
      
      eventbuilder->emc->readFromGl3Trg(&gl3_data);    
    }
    //    eventbuilder->emc->readFromDaqReader(rdr);              //Run 8
    //    eventbuilder->emc->getTrackingSeeds();
	
    //    printf("field=%f",bField);
    // need temporary track memory...
    
    int i;
    int ret;
    float max=-10000;
    float min=10000;
    
    for(i=0;i<24;i++) 
      {	
	tracker[i]->para.bField = fabs(bField); 
	// +1 , -1 or 0 if zero field...
	tracker[i]->para.bFieldPolarity = (bField>0) ? 1 : -1;		

	/*
	//for Run 9
	tracker[i]->nHits = 0 ;
	tracker[i]->setTrackingAngles(i+1);   	
	tracker[i]->sector_ID=i;

	for(int r=1;r<=6;r++) {
	  dd = rdr->det("tpx")->get("cld_raw",i+1,r) ;
	  while(dd && dd->iterate()) {
	    
	    
	    if(dd->ncontent == 0) continue ;
	    
	    LOG(DBG,"Loading sector %d, %d bytes",i+1,dd->ncontent) ;
	    tracker[i]->readSectorFromESB(i+1, (char *)dd->Void, dd->ncontent/4) ;
	    
	  }
	}
	int size;
	*/

	//for Run 8
	Tracking_pTPC=getTPC(rdr,i);
	if(!Tracking_pTPC) continue;
	int size;
	
	
	tracker[i]->nHits=0;	
	tracker[i]->setTrackingAngles(i+1);   	
	tracker[i]->sector_ID=i;
	tracker[i]->Tracking_load_cluster(Tracking_pTPC);
	

	//cout<<"load_time  "<<tracker[i]->diftime<<endl;	
	

	//full tracking 
	tracker[i]->Tracking_track_sector((char *)sectp,size);
	//tracking with BTOW as seeds
	//	tracker[i]->Tracking_track_sector((float*)eventbuilder->emc->trackingSeeds[i],eventbuilder->emc->nTrackingSeeds[i],(char *)sectp,size);
	

      

	LOG(DBG, "SECP size = %d",sectp->bh.length*4 + sectp->banks[0].len*4);
	
	int n = eventbuilder->Tracking_readSectorTracks((char *)sectp);
	//	cout<<n<<" tracks in sector "<<i+1<<endl;
	
	      
	//cout<<"track_time  "<<((float)tracker[i]->diftime)/((float)n)<<endl;
	
      }

    
    
    eventbuilder->finalizeReconstruction();
    printf("%d  tracks \n",eventbuilder->getNGlobalTracks());
    
    eventbuilder->triggerDecider->decide(rdr->seq);
    
  }
    
  free(sectp);
  //  delete []tracker;
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


