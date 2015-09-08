//:>------------------------------------------------------------------
//: FILE:       online_tracking_subsector.h
//: HISTORY:
//:             28oct1996 version 1.00
//:             23aug1999 ppy printVols and printRows deleted
//:             26jan2000 ppy destructor added
//:             27jan2000 ppy VOLUME, ROW and AREA classes replaced by
//:                           FtfContainer
//:             11feb2000 ppy timeout added, variables initialCpuTime and
//:                           initialRealTime added
//:<------------------------------------------------------------------

#ifndef online_tracking_subsector_
#define online_tracking_subsector_
#include <string.h>
#include "FtfGeneral.h"
#include "FtfPara.h"
#include "FtfHit.h"
#include "FtfTrack.h"
#include "FtfMcTrack.h"
#include   <stdio.h>   
#include   <stdlib.h>                                                
#include   <time.h>   

#define MILLION 1000000L
#define THOUSAND 1000  

class online_tracking_subsector 

{

    public:

 long diftime;
 struct itimerspec nvalue, ovalue;
 timer_t timeid;

int timer_run();
int timer_stop();
void timer_print(char *name);

    online_tracking_subsector( ) ;
    ~online_tracking_subsector( ) ;
    friend class FtfTrack ;
    int    dEdx_s                    ( int a ) ;
    int     getTracks               ( ) ;
    int     HybridgetTracks               ( ) ;
    int     getTracks               ( float* seed,int seednum,int seedlayer=4) ;
	int     getTracks (float* seed,int seednum,int seedlayer,float seed_eta_phi_cut ) ;
    void    mergePrimaryTracks      ( ) ;
    double  process ( ) ;
    double  process (float* seed,int seednum ) ;
    double  process (float* seed,int seednum,int seedlayer,float seed_eta_phi_cut ) ;
    int     reset                   ( ) ;
    int     setConformalCoordinates ( ) ;
    int     setPointers             ( ) ;
    double  CpuTime                 ( ) ;
    double  RealTime                ( ) ;

    //
    //
    int           nHits      ;  
    int           nHitsOutOfRange ;
    int           maxHits    ;  
    FtfHit        *hit       ;  
    int           nTracks    ; 
    FtfTrack      *track     ;  
    FtfPara       para       ;
    int           maxTracks  ;
    int           nMcTracks  ;
    FtfMcTrack    *mcTrack    ;
    FtfContainer  *volumeC ;
    FtfContainer  *rowC    ;
    FtfContainer  *trackC  ;
    double        initialCpuTime ;
    double        initialRealTime ;
    double        cpuTime ;
    double        realTime ;
    private: 
    FtfTrack      *currentTrack ;
    
} 

;

#endif

