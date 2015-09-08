/*:>-------------------------------------------------------------------
**: FILE:     online_tracking_sector.h 
**: HISTORY:  may 29, 1999  Frame for sl3 tracking
**:           aug 23, 1999  add setup with number of hits and tracks
**:           sep 28, 1999  cs: new input functions (pointer to
**:                         bank instead of FILE) using
**:                         daqFormats.h
**:           oct 25, 1999  ppy: use sl3CoordianteTransform header
**:           11/24/99      cle: include <L3Formats.h> intead daqFormats.h
**:                         commented out fillUSTracks
**:                         added public sectorNr variable
**:           12/03/99      ppy: sectorGeometry class added
**:                              variable sectorGeo added
**:                              clean extra includes    
**:           12/06/99      ppy: method added to check whether track can be merged
**:           01/26/00      ppy: delete rawToGlobal declaration
**:           feb 10, 2000  ppy: add xyError and zError
**:           mar 21  2000  ppy: add fill Hits
**:           apr  8  2000  ppy: add variables to cut hit based on charge and time
**:           apr 18  2000  ppy: modify readMezzannine to include rb and mz as arguments
**:           apr 19  2000  ppy: dEdx from Christof
**:           apr 19  2000  cs ppy: include FtfDedx added
**:<------------------------------------------------------------------*/

#ifndef online_tracking_sector_
#define online_tracking_sector_
#include "l3CoordinateTransformer.h"
#include "online_tracking_subsector.h"
#include "online_tracking_TpcHitMap.h"
#include "FtfDedx.h"
#include "daqFormats.h"
#include "sizes.h"
#include "FtfSl3.h"
#include "mTrackEvent.h"

#include <rtsLog.h>

#ifdef OLD_DAQ_READER
#include <evpReader.hh>

#else /* OLD_DAQ_READER */
#include <DAQ_READER/daqReader.h>

#endif /* OLD_DAQ_READER */
#include <rtsSystems.h>
#include <daqFormats.h>

#ifndef OLD_DAQ_READER
#include <DAQ_TPC/daq_tpc.h>

#endif /* OLD_DAQ_READER */
class online_tracking_sector: public online_tracking_subsector 
{
public:
    int                        sector_ID;
    tpc_t*                     Tracking_pTPC;
private:   
    int                        sectorNr; 
    short                      debugLevel  ;
    double                     xyError ;
    double                     zError ;
    sectorGeometry*            sectorGeo ;
public:
    int                        minTimeBin ;
    int                        maxTimeBin ;
    int                        minClusterCharge ;
    int                        maxClusterCharge ;
    int                        embedded;
    FtfDedx*                   fDedx;
    online_tracking_TpcHitMap* tpcHitMap;
    // flag to show that the current event being read in is embedded
    //
    //     Constructor
    //
    online_tracking_sector ( float bField,int in_sector_ID, char *map_fname=0, char* HLTparameters="HLTparameters", char* beamlineFile = "beamline")    
    { 

	LOG(DBG,"field %f, sector %d, map %s",bField,in_sector_ID,map_fname) ;

        debugLevel            = 0 ;
        xyError               = 0.3 ;
        zError                = 1.0 ;
        minTimeBin            = 0 ;
        maxTimeBin            = 380 ;
        minClusterCharge      = 80 ;
        maxClusterCharge      = 100000 ;
        coordinateTransformer = new l3CoordinateTransformer();
        coordinateTransformer->Set_parameters_by_hand(0.581, 200.668, 201.138 );
        sectorGeo             = new sectorGeometry[24];
        fDedx                 = new FtfDedx((l3CoordinateTransformer*)coordinateTransformer, 0, 0.7, 0, HLTparameters);
        //////////////////////////////////////////
        setup();
        para.bField           = fabs(bField); 
        // +1 , -1 or 0 if zero field...
        para.bFieldPolarity   = (bField>0) ? 1 : -1;
        ////// check parms  ///////////////
        setXyError(.12) ; 
        setZError(.24) ;
        para.ptMinHelixFit    = 0.;
        para.maxChi2Primary   = 0.;  
        para.trackChi2Cut     = 10 ; 
        para.hitChi2Cut       = 50 ;
        para.goodHitChi2      = 20 ; 
        reset();
        setTrackingAngles(in_sector_ID+1); 
        sector_ID             = in_sector_ID;
        //////////////////////////////////////////

	LOG(DBG,"Doing hit map...") ;

	tpcHitMap = new online_tracking_TpcHitMap(HLTparameters, in_sector_ID);
	LOG(DBG,"hitMap allocated... Loading... [%p]",tpcHitMap) ;

	if(!map_fname) {
            tpcHitMap->loadMap();
	}
	else {
            tpcHitMap->loadMap(map_fname);
	}

	readHLTparameters(HLTparameters);
	readBeamline(beamlineFile);        
    }

    ;
    ~online_tracking_sector ( ) 
    
    {

        if ( track != 0 ) delete []track ;
        if ( hit   != 0 ) delete []hit   ;
        delete []sectorGeo ;
        delete fDedx;
        delete tpcHitMap;
	delete coordinateTransformer;
    }

    ;

    void put_bField(double bField) ;

    void readHLTparameters(char* HLTparameters = "HLTparameters");

    void readBeamline(char* beamlineFile = "beamline");

    int readScalers(char* scalerFile = "sc.asc");

    //
    //  Coordinate Transformer classes
    //
    inline l3CoordinateTransformer* getCoordinateTransformer()
    
    { 

        return (l3CoordinateTransformer*)coordinateTransformer ; 
        
    } 

    ;

    inline online_tracking_TpcHitMap* getTpcHitMap()
    { return (online_tracking_TpcHitMap*) tpcHitMap; };

    //
    //  Sector phase space
    //
    sectorGeometry* getSectorGeometry ( int n ) 
    
    {

        if ( n < 1 && n > 24 ) 
        
            { 

                printf ( "online_tracking_sector::getSectorGeometry: wrong sector %d, returning sector 1 \n", n ) ;
                n = 0 ;
            
            }

        return &(sectorGeo[n]);
        
    }

    //
    int   fillTracks      ( int maxBytes, char* buff, unsigned int token ) ;

    int   Tracking_track_sector ( char *out_buf,int &out_buf_used ) ;
    int   Tracking_track_sector (float* seed,int seednum , char *out_buf,int &out_buf_used ) ;
    int   Tracking_track_sector (float* detadphi,int region_num, int region_type, char *out_buf,int &out_buf_used ) ;
    int   Tracking_load_cluster ( tpc_t *in_pTPC );

    // Tonko:
    int readSectorFromESB(int sector, char *mem, int len) ;

    int   Tracking_load_cluster ( mTrackEvent *in_mTrackEvent,int secID);
    int   Tracking_load_cluster (int MCtrack_number,float* vertex, mTrackEvent *in_mTrackEvent,int secID,int *ID_list,int* hit_number_list);

    int   fillHits        ( unsigned int maxBytes, char* buff, unsigned int token ) ;
    float getXyError      ( ) 
    
    { 

        return xyError ; 
        
    } 

    ;
    float getZError       ( ) 
    
    { 

        return zError  ; 
        
    } 

    ;
    int   getNrTracks     ( ) ;
    int   canItBeMerged   ( FtfTrack* thisTrack ) ;
    int   dEdx ( ) ;
    int   processSector   ( ) ;
    int   processSector   (float* seed,int seednum ) ;
    int   processSector   (float* seed,int seednum,int seedlayer,float seed_eta_phi_cut ) ;
    //int   processData     (TPCSECLP* seclp, char* trackBuffer, int maxTrackBytes, int& nTrackBytes,
    //char* hitBuffer, int maxHitBytes, int& nHitBytes);
    //int   readMezzanine   ( int sector, int readOutBoard,
    //int MezzanninneNr, struct TPCMZCLD_local *mzcld );
    //int   readSector      ( struct bankHeader *bank ) ; 
    //int   readSector (DATAP *datap, int sector);
    //int   readSector      ( struct TPCSECLP *seclp1, struct TPCSECLP *seclp2 ); 
    int setTrackingAngles(int hypersector);
    int readSectorFromEvpReader(int sector);
    int   setParameters   ( ) ;
    void  setCoordinateTransformer ( l3CoordinateTransformer* in )
    
    { 

        coordinateTransformer = in ; 
        
    } 

    ;

    void setTpcHitMap ( online_tracking_TpcHitMap* in)
    {tpcHitMap = (online_tracking_TpcHitMap*) in; };


    void  setDebugLevel   ( int _debugLevel ) 
    
    { 

        debugLevel = _debugLevel ; 
        
    } 

    ;
    void  setSector       ( int _sector     ) 
    
    { 

        sectorNr   = _sector     ; 
        
    } 

    ;
    void  setXyError      ( float _xyError  ) 
    
    { 

        xyError    = _xyError    ; 
        
    } 

    ;
    void  setZError       ( float _zError   ) 
    
    { 

        zError     = _zError     ; 
        
    } 

    ;
    int   setup           ( int maxHitsIn=szSL3_default_mHits, int maxTracksIn=szSL3_default_mTracks ) ;
    void  Print() const
    
    {

        
    }

    ;
private:
    l3CoordinateTransformer*     coordinateTransformer ;
    
}

    ;

#endif

