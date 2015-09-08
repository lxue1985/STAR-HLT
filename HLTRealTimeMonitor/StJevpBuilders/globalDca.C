#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/types.h>
#include <stdlib.h>

// this needs to be always included
#include <DAQ_READER/daqReader.h>
#include <DAQ_READER/daq_dta.h>

// only the detectors we will use need to be included
// for their structure definitions...
#include <DAQ_TPX/daq_tpx.h>
#include <DAQ_HLT/daq_hlt.h>
#include "HLTFormats.h"

#include "math.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"

int main(int argc, char *argv[])
{
  TFile* f1 = new TFile("dca.root", "RECREATE");
  TNtuple* ntp = new TNtuple("dca", "dca", "dcaXy:dcaZ:nHits:p");

	char *mountpoint = NULL;
	char *inputFile = argv[argc-1];

	class daqReader *evp ;			// tha main guy
	printf("inputFile: %s\n", inputFile);
	evp = new daqReader(inputFile) ;	// create it with the filename argument..
	if(mountpoint) {
	  evp->setEvpDisk(mountpoint);
	}

	int good=0;
	int bad=0;
	
	for(;;) {
	        char *ret = evp->get(0,EVP_TYPE_ANY);

		printf("ret: %s \n", ret);
       
		if(ret) {
		  if(evp->status) {
		    printf("evp status is non-null [0x08X, %d dec] \n",evp->status) ;
			continue ;
		  }
		  good++;
		}
		else {    // something other than an event, check what.
		  switch(evp->status) {
		  case EVP_STAT_OK:   // just a burp!
		    continue;
		  case EVP_STAT_EOR:
		    printf("End of Run/File \n");
		    if(evp->IsEvp()) {   // but event pool, keep trying...
		      printf("Wait a second... \n");
		      sleep(1);
		      continue;
		    }
		    break;        // file, we're done...
		  case EVP_STAT_EVT:
		    bad++;
		    printf("Problem getting event - skipping [good %d, bad %d] \n",good,bad);
		    sleep(1);
		    continue;
		  case EVP_STAT_CRIT:
		    printf("evp->status CRITICAL (?) \n") ;
		    return -1;
		  }
		}

		if(evp->status == EVP_STAT_EOR) {
		  printf("End of File reached... \n") ;
			break ;	// of the for() loop...
		}

	       
		daq_dta *dd ;	// generic data pointer; reused all the time


		printf("sequence %d: token %4d, trgcmd %2d, daqcmd %2d, time %u, detectors 0x%08X (status 0x%X) \n",evp->seq, evp->token, evp->trgcmd, evp->daqcmd,
		       evp->evt_time, evp->detectors, evp->status) ;

		printf("trginfo: seq = #%d  token = %d detectors = 0x%x triggers = 0x%x \n",
			   evp->seq,
			   evp->token,
			   evp->detectors,
			   evp->daqbits);

		/*************************** HLT10 **************************/
		int found = 0 ;
		dd = evp->det("hlt")->get("gl3") ;

		HLT_EVE  *hlt_eve ;
		HLT_TOF  *hlt_tof ;
		HLT_PVPD *hlt_pvpd ;
		HLT_EMC  *hlt_emc ;
		HLT_GT   *hlt_gt ;
		HLT_PT   *hlt_pt ;
		HLT_NODE *hlt_node ;
		HLT_HIPT *hlt_hipt ;
		HLT_DIEP *hlt_diep ;
		HLT_HF *hlt_hf ;
		
		
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

		for(int i=0; i<hlt_gt->nGlobalTracks; i++)
		  {
		    hlt_track track = hlt_gt->globalTrack[i];
		    double p = track.pt*sqrt(1.+pow(track.tanl,2));
		    double dcaX = track.r0*cos(track.phi0) - hlt_eve->lmVertexX;
		    double dcaY = track.r0*sin(track.phi0) - hlt_eve->lmVertexY;
		    double cross = dcaX*sin(track.psi) - dcaY*cos(track.psi);
		    double theSign = (cross>=0) ? 1. : -1.;
		    double dcaXy = theSign*sqrt(pow(dcaX,2)+pow(dcaY,2));
		    double dcaZ = track.z0 - hlt_eve->lmVertexZ;
		    ntp->Fill(dcaXy, dcaZ, track.nHits, p);
		
		  }
	}
	TCanvas* c1 = new TCanvas();
	ntp->Draw("dcaXy>>hh(100, -10., 10.)");
	hh->Fit("gaus", "", "", -3., 3.);
	c1->SaveAs("dcaXy.gif");
	ntp->Draw("dcaZ>>hh(100, -50., 50.)");
	hh->Fit("gaus", "", "", -5., 5.);
	c1->SaveAs("dcaZ.gif");
	ntp->Write();
	f1->Close();
	//	delete evp ;	// cleanup i.e. if running through a set of files.

	return 0 ;
}

