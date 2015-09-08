//:>------------------------------------------------------------------
//: FILE:       online_tracking_subsector.cxx
//: HISTORY:
//:             28oct1996 version 1.00
//:             03jun1999 ppy para.fillTracks included. Merging only when tracks filled
//:             03jun1999 ppy add a function for real time, clock gives cpu time
//:             06may1999 ppy getTracks returns 1 for error
//:             10aug1999 ppy nHitsForSegment set to at least 3 for
//:                           secondary search.
//:             23aug1999 ppy ClassImp added with ROOT flag
//:             21dec1999 ppy printf replaced by fprintf(stderr,...
//:             26jan2000 ppy malloc replaced with new, destructor function added
//:             27jan2000 ppy refHit replaced by xRefHit and yRefHit
//:             28jan2000 ppy track id from 1 to N
//:              1feb2000 ppy track id starting at 1
//:             11feb2000 ppy timeout added, variables initialCpuTime and
//:                           initialRealTime added
//:             10apr2000 ppy deleteCandidate added when track is merged
//:             21aug2000 ppy replace prints with l3Log
//:<------------------------------------------------------------------
//:>------------------------------------------------------------------
//: CLASS:       online_tracking_subsector, steers track finding
//: DESCRIPTION: Functions associated with this class
//: AUTHOR:      ppy - Pablo Yepes, yepes@physics.rice.edu
//:>------------------------------------------------------------------
#include "online_tracking_subsector.h"
#include   <stdio.h>   
#include   <stdlib.h>
#include <signal.h>
#include   <time.h>
#include   <sys/time.h> 

#include "../GL3/profiler.hh" 

PROFILER_DECLARE ;

//*********************************************************************
//      Initializes the package
//*********************************************************************
online_tracking_subsector::online_tracking_subsector ( ) 

{

    //
    hit        = 0 ;  
    track      = 0 ;
    mcTrack    = 0 ;
    trackC     = 0 ;
    volumeC    = 0 ;
    rowC       = 0 ;
    nHitsOutOfRange = 0 ;

    if (timer_create(CLOCK_REALTIME, NULL, &timeid) == -1) {
        perror("Failed to create a timer based on CLOCK_REALTIME");
        return ;
    }
    
    
}

//*********************************************************************
//      Initializes the package
//*********************************************************************
online_tracking_subsector::~online_tracking_subsector ( ) 

{

    //
    if ( mcTrack ) delete[] mcTrack ;
    if ( volumeC ) delete[] volumeC  ;
    if ( rowC    ) delete[] rowC ;
    if ( trackC    ) delete[] trackC ;
    
}

//*********************************************************************
//      Steers the tracking
//*********************************************************************

double online_tracking_subsector::process ( float* seed,int seednum,int seedlayer,float seed_eta_phi_cut ) 

{  

    //-----------------------------------------------------------------
    //        Make sure there is something to work with
    //------------------------------------------------------------------ 
    if ( nHits <= 0 ) 
    
    {

        if ( para.infoLevel > 2 ) LOG(ERR, "fft: Hit structure is empty \n " ) ;
        return 1 ;
        
    }

    //
    initialCpuTime  = CpuTime ( );
    initialRealTime = RealTime ( );
    //
    //        General initialization
    //
    if ( para.init == 0 ) 
    
    {

        if ( reset ( ) ) return 1 ;
        
    }

    //
    //      Event reset and set pointers
    //
    if ( para.eventReset  && setPointers ( ) ) return 1 ;   
    //
    //      Build primary tracks now
    //
    short i ;
    para.primaries = 1 ;
    for ( i = 0 ; i < para.nPrimaryPasses ; i++ )
    if ( getTracks (seed,seednum,seedlayer,seed_eta_phi_cut ) ) break ;
    para.primaries = 0 ;
    for ( i = 0 ; i < para.nSecondaryPasses ; i++ )
    if ( getTracks (seed,seednum,seedlayer,seed_eta_phi_cut ) ) break ;
    //       if ( para.dEdx ) dEdx ( ) ;
    cpuTime  = CpuTime  ( ) - initialCpuTime  ;
    realTime = RealTime ( ) - initialRealTime ;
    
    #ifdef DEBUG
    if ( para.infoLevel > 0 )
    LOG(ERR, "online_tracking_subsector::process: cpu %7.3f real %f7.2 \n", cpuTime, realTime ) ;
    
    #endif
    //printf("online_tracking_subsector::process: cpu %7.3f real %7.2f \n", cpuTime, realTime );
    return cpuTime ;
    
} 




double online_tracking_subsector::process ( float* seed,int seednum ) 

{  

    //-----------------------------------------------------------------
    //        Make sure there is something to work with
    //------------------------------------------------------------------ 
    if ( nHits <= 0 ) 
    
    {

        if ( para.infoLevel > 2 ) LOG(ERR, "fft: Hit structure is empty \n " ) ;
        return 1 ;
        
    }

    //
    initialCpuTime  = CpuTime ( );
    initialRealTime = RealTime ( );
    //
    //        General initialization
    //
    if ( para.init == 0 ) 
    
    {

        if ( reset ( ) ) return 1 ;
        
    }

    //
    //      Event reset and set pointers
    //
    if ( para.eventReset  && setPointers ( ) ) return 1 ;   
    //
    //      Build primary tracks now
    //
    short i ;
    para.primaries = 1 ;
    for ( i = 0 ; i < para.nPrimaryPasses ; i++ )
    if ( getTracks (seed,seednum ) ) break ;
    para.primaries = 0 ;
    for ( i = 0 ; i < para.nSecondaryPasses ; i++ )
    if ( getTracks (seed,seednum ) ) break ;
    //       if ( para.dEdx ) dEdx ( ) ;
    cpuTime  = CpuTime  ( ) - initialCpuTime  ;
    realTime = RealTime ( ) - initialRealTime ;
    
    #ifdef DEBUG
    if ( para.infoLevel > 0 )
    LOG(ERR, "online_tracking_subsector::process: cpu %7.3f real %f7.2 \n", cpuTime, realTime ) ;
    
    #endif
    //printf("online_tracking_subsector::process: cpu %7.3f real %7.2f \n", cpuTime, realTime );
    return cpuTime ;
    
} 



    


#ifdef HYBRID_TRACKER

#include "hyb/HybridTrack.h"
#include "hyb/HybridHit.h"
#include "hyb/HybridPara.h"

int online_tracking_subsector::HybridgetTracks ( )

{

    //
    //     Set conformal coordinates if we are working with primaries
    //
    int nHitsSegment   = (short)para.nHitsForSegment;
    if ( para.primaries )

    {

        setConformalCoordinates ( ) ;
        para.minHitsForFit = 1 ;
        para.nHitsForSegment = max(2,nHitsSegment);

    }

    else

    {

        para.minHitsForFit = 2 ;
        para.nHitsForSegment = max(3,nHitsSegment);

    }
        float maxchi=para.maxChi2Primary;
        para.maxChi2Primary=3.0;

    //
    //               Loop over rows
    //
    for ( int ir = para.nRowsPlusOne - 1 ; ir>=para.minHitsPerTrack ; ir--)

        {

        //
        //           Loop over hits in this particular row
        //
        if ( rowC[ir].first &&  (((HybridHit *)rowC[ir].first)->row) < para.rowEnd )
        break ;
        //if(rowC[ir].first!=0)
        //printf("%d  %f %f %f\n",((HybridHit *)rowC[ir].first)->x,((HybridHit *)rowC[ir].first)->y,((HybridHit *)rowC[ir].first)->z,ir);
        //printf("  %f %f %f in %d\n",((HybridHit *)rowC[ir].first)->x,((HybridHit *)rowC[ir].first)->y,((HybridHit *)rowC[ir].first)->z,ir);
        //    if ( (((HybridHit *)rowC[ir].first)->row) < para.rowEnd ) break ;
        for ( HybridHit *firstHit = (HybridHit *)rowC[ir].first ;
        firstHit != 0 ;
        firstHit = (HybridHit *)(firstHit->nextRowHit) )

                {

            //
            //     Check hit was not used before
            //
            if ( firstHit->track != 0  ) continue ;
            //
            //     One more track
            nTracks++ ;
            //
            //
            if ( nTracks > maxTracks )

            {

                LOG(ERR, "\n online_tracking_subsector::getTracks: Max nr tracks reached !") ;
                nTracks = maxTracks  ;
                                para.maxChi2Primary=maxchi;
                return 1 ;

            }

            //
            //     Initialize variables before going into track hit loop
            //
            HybridTrack *thisTrack = (HybridTrack *)(&track[nTracks-1]);
            thisTrack->para     = &para ;
            thisTrack->id       = nTracks ;
            thisTrack->firstHit = thisTrack->lastHit = firstHit ;
            thisTrack->innerMostRow = thisTrack->outerMostRow = firstHit->row ;
            thisTrack->xRefHit  = firstHit->x ;
            thisTrack->yRefHit  = firstHit->y ;
            thisTrack->xLastHit = firstHit->x ;
            thisTrack->yLastHit = firstHit->y ;

            #ifdef TRDEBUG
            thisTrack->debugNew ( ) ;

            #endif
            //
            //              Set fit parameters to zero
            //
            thisTrack->reset ( ) ;
            //
            //      Go into hit looking loop
            //
            if ( thisTrack->buildTrack ( firstHit, (HybridContainer*)volumeC ) )

            {

                //
                //    Merge Tracks if requested
                //
                if ( para.primaries &&
                para.mergePrimaries == 1 &&
                para.fillTracks &&
                thisTrack->mergePrimary( (HybridContainer*)trackC )  )

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

        //       End loop over rows
        //
        //    Check time
        //
        if ( CpuTime() - initialCpuTime > para.maxTime )

        {

            LOG(ERR, "online_tracking_subsector::getTracks: tracker time out after %f\n",
            CpuTime() - initialCpuTime ) ;
            break ;

        }


    }

    //
    para.nHitsForSegment = nHitsSegment ;
        para.maxChi2Primary=maxchi;
    //
    return 0 ;

}



double online_tracking_subsector::process (  )

{

    //-----------------------------------------------------------------
    //        Make sure there is something to work with
    //------------------------------------------------------------------ 
    if ( nHits <= 0 )

    {

        if ( para.infoLevel > 2 ) LOG(ERR, "fft: Hit structure is empty \n " ) ;
        return 1 ;

    }

    //
    initialCpuTime  = CpuTime ( );
    initialRealTime = RealTime ( );
    //
    //        General initialization
    //
    if ( para.init == 0 )

    {

        if ( reset ( ) ) return 1 ;

    }

    //
    //      Event reset and set pointers
    //
    if ( para.eventReset  && setPointers ( ) ) return 1 ;
    //
    //      Build primary tracks now
    //
    short i ;
    para.primaries = 1 ;
    for ( i = 0 ; i < para.nPrimaryPasses ; i++ )
    {if ( HybridgetTracks ( ) ) break ;
    if ( getTracks ( ) ) break ;}
    //
    //      Look for secondaries
    //
    para.primaries = 0 ;
    for ( i = 0 ; i < para.nSecondaryPasses ; i++ )
    {if ( HybridgetTracks ( ) ) break ;
    if ( getTracks ( ) ) break ;}
    //      if ( para.dEdx ) dEdx ( ) ;
    cpuTime  = CpuTime  ( ) - initialCpuTime  ;
    realTime = RealTime ( ) - initialRealTime ;
    #ifdef DEBUG
    if ( para.infoLevel > 0 )
    LOG(ERR, "online_tracking_subsector::process: cpu %7.3f real %f7.2 \n", cpuTime, realTime ) ;

    #endif
    //printf("online_tracking_subsector::process: cpu %7.3f  \n", cpuTime );
    return cpuTime ;

}


#else



double online_tracking_subsector::process (  )

{

    //-----------------------------------------------------------------
    //        Make sure there is something to work with
    //------------------------------------------------------------------ 
    if ( nHits <= 0 )

    {

        if ( para.infoLevel > 2 ) LOG(ERR, "fft: Hit structure is empty \n " ) ;
        return 1 ;

    }

    //
    initialCpuTime  = CpuTime ( );
    initialRealTime = RealTime ( );
    //
    //        General initialization
    //
    int pval = PROFILER(0) ;
    if ( para.init == 0 )

    {

        if ( reset ( ) ) return 1 ;

    }
    pval = PROFILER(pval) ;	// 1 us
    //
    //      Event reset and set pointers
    //
    if ( para.eventReset  && setPointers ( ) ) return 1 ;

    pval = PROFILER(pval) ;	// 620 us

    //
    //      Build primary tracks now
    //
    int i ;
    para.primaries = 1 ;
    for ( i = 0 ; i < para.nPrimaryPasses ; i++ )
    if ( getTracks ( ) ) break ;
    //
    //      Look for secondaries
    //
    pval = PROFILER(pval) ;	// 3300 us
    para.primaries = 0 ;
    for ( i = 0 ; i < para.nSecondaryPasses ; i++ )
    if ( getTracks ( ) ) break ;
    //      if ( para.dEdx ) dEdx ( ) ;
    cpuTime  = CpuTime  ( ) - initialCpuTime  ;
    realTime = RealTime ( ) - initialRealTime ;
    #ifdef DEBUG
    if ( para.infoLevel > 0 )
    LOG(ERR, "online_tracking_subsector::process: cpu %7.3f real %f7.2 \n", cpuTime, realTime ) ;

    #endif
    //printf("online_tracking_subsector::process: cpu %7.3f  \n", cpuTime );

    pval = PROFILER(pval) ;	// 1.7
    return cpuTime ;

}


#endif


//********************************************************************
//     Calculates deposited Energy
//********************************************************************
int online_tracking_subsector::dEdx_s (int a ) 

{

    for ( int i = 0 ; i<nTracks ; i++ )
    
    {

        track[i].dEdx( ) ;
        
    }
return 0;
    
}

//**********************************************************************
//	Recontruct primary tracks
//**********************************************************************
int online_tracking_subsector::getTracks (float* seed,int seednum,int seedlayer,float seed_eta_phi_cut ) 

{

    //
    //     Set conformal coordinates if we are working with primaries
    //
    int nHitsSegment   = (short)para.nHitsForSegment;  
    if ( para.primaries ) 
    
    {

        setConformalCoordinates ( ) ;
        para.minHitsForFit = 1 ;
        para.nHitsForSegment = max(2,nHitsSegment);
        
    }

    else 
    
    {

        para.minHitsForFit = 2 ;
        para.nHitsForSegment = max(3,nHitsSegment);
        
    }

    //
    //               Loop over rows
    //


    int si;
    for(si=0;si<seednum;si++)
    {		
		
		int ir=int(seed[5*si]);
        //
        //           Loop over hits in this particular row
        //
        if ( rowC[ir].first &&  (((FtfHit *)rowC[ir].first)->row) < para.rowEnd ) 
        break ;
        //if(rowC[ir].first!=0)
        //printf("%d  %f %f %f\n",((FtfHit *)rowC[ir].first)->x,((FtfHit *)rowC[ir].first)->y,((FtfHit *)rowC[ir].first)->z,ir);
        //printf("  %f %f %f in %d\n",((FtfHit *)rowC[ir].first)->x,((FtfHit *)rowC[ir].first)->y,((FtfHit *)rowC[ir].first)->z,ir);
        //    if ( (((FtfHit *)rowC[ir].first)->row) < para.rowEnd ) break ;
        for ( FtfHit *firstHit = (FtfHit *)rowC[ir].first ;
        firstHit != 0 ;
        firstHit = (FtfHit *)(firstHit->nextRowHit) ) 
        
		{

            //
            //     Check hit was not used before
            //
            if ( firstHit->track != 0  ) continue ;
			if ( firstHit->eta<seed[5*si+1]||firstHit->eta>seed[5*si+2]  ) continue ;
			if ( firstHit->phi<seed[5*si+3]||firstHit->phi>seed[5*si+4]  ) continue ;
		
            //
            //     One more track
            //
            nTracks++ ;
            //
            //
            if ( nTracks > maxTracks )
            
            {

                LOG(ERR, "\n online_tracking_subsector::getTracks: Max nr tracks reached !") ;
                nTracks = maxTracks  ;
                return 1 ;
                
            }

            //
            //     Initialize variables before going into track hit loop
            //
            FtfTrack *thisTrack = &track[nTracks-1];
            thisTrack->para     = &para ;
            thisTrack->id       = nTracks ;
            thisTrack->firstHit = thisTrack->lastHit = firstHit ;
            thisTrack->innerMostRow = thisTrack->outerMostRow = firstHit->row ;
            thisTrack->xRefHit  = firstHit->x ;
            thisTrack->yRefHit  = firstHit->y ;
            thisTrack->xLastHit = firstHit->x ;
            thisTrack->yLastHit = firstHit->y ;
            
            #ifdef TRDEBUG
            thisTrack->debugNew ( ) ;
            
            #endif
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
                if ( para.primaries &&
                para.mergePrimaries == 1 &&
                para.fillTracks &&
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

        //       End loop over rows
        //
        //    Check time
        //
        if ( CpuTime() - initialCpuTime > para.maxTime ) 
        
        {

            LOG(ERR, "online_tracking_subsector::getTracks: tracker time out after %f\n", 
            CpuTime() - initialCpuTime ) ;
            break ;
            
        }

        
    }

    //
    para.nHitsForSegment = nHitsSegment ;  
    //
    return 0 ;

}




int online_tracking_subsector::getTracks (float* seed,int seednum,int seedlayer ) 

{

    //
    //     Set conformal coordinates if we are working with primaries
    //
    int nHitsSegment   = (short)para.nHitsForSegment;  
    if ( para.primaries ) 
    
    {

        setConformalCoordinates ( ) ;
        para.minHitsForFit = 1 ;
        para.nHitsForSegment = max(2,nHitsSegment);
        
    }

    else 
    
    {

        para.minHitsForFit = 2 ;
        para.nHitsForSegment = max(3,nHitsSegment);
        
    }

    //
    //               Loop over rows
    //
    int si;
    for(si=0;si<seednum;si++)
    {

        float dis=1000;
        float cdis=1000;
        FtfHit *firstHit;


        for ( int ir = seed[5*si] ; ir>=seed[5*si+1] ; ir--) 
        
        {

            //
            //           Loop over hits in this particular row
            //
            if ( rowC[ir].first &&  (((FtfHit *)rowC[ir].first)->row) < para.rowEnd ) 
            break ;
            //if(rowC[ir].first!=0)
            // printf("  %f %f %f in %d\n",((FtfHit *)rowC[ir].first)->x,((FtfHit *)rowC[ir].first)->y,((FtfHit *)rowC[ir].first)->z,ir);
            //    if ( (((FtfHit *)rowC[ir].first)->row) < para.rowEnd ) break ;
            //    if ( (((FtfHit *)rowC[ir].first)->row) < para.rowEnd ) break ;
            for ( FtfHit *sfirstHit = (FtfHit *)rowC[ir].first ;
            sfirstHit != 0 ;
            sfirstHit = (FtfHit *)(sfirstHit->nextRowHit) ) 
            
            {

				 if ( sfirstHit->track != 0  ) continue ;

                // cdis=sqrt(pow(sfirstHit->x-seed[5*si+2],2)+pow(sfirstHit->y-seed[5*si+3],2)+pow(sfirstHit->z-seed[5*si+4],2));
                // if(cdis<dis)
                // {

                    // dis=cdis;
                    // firstHit=sfirstHit;
                    
                // }
				 /*
				 float sx=seed[5*si+2];
				 float sy=seed[5*si+3];
				 float sz=seed[5*si+4];
				 float sr=sqrt(sx*sx+sy*sy);
				 float angle=atan(sy/sx);
				 if(angle<0)angle=angle+3.14159;
				 if(angle<0&&sx>0)angle=angle+2*3.14159;
				 */

				 float rcut=333.;
				 float beamlineZRange=100.;
				 float towerZRange=8.;
				 float towerLocalXRange=10.;

				 float sr=seed[5*si+2];
				 float angle=seed[5*si+3];
				 float sz=seed[5*si+4];
				 angle=angle-3.14159/2;
				 float px=sfirstHit->x*cos(angle)+sfirstHit->y*sin(angle);
				 float py=-sfirstHit->x*sin(angle)+sfirstHit->y*cos(angle);
				 float cr=sqrt(rcut*rcut-sr*sr/4);// x pos of center
				 float pcr;// p to center
				 if(px>0){pcr=sqrt((px+cr)*(px+cr)+(py-sr/2)*(py-sr/2));}
			     else{pcr=sqrt((px-cr)*(px-cr)+(py-sr/2)*(py-sr/2));}
				
				 if(pcr>rcut+towerLocalXRange) continue;
				 
				 float pr=sqrt(px*px+py*py);
				 if(fabs(sfirstHit->z-sz*pr/sr)>beamlineZRange*(sr-pr)/sr+towerZRange) continue;
				firstHit=sfirstHit;

				
				
				
				
				
				
        
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

					LOG(ERR, "\n online_tracking_subsector::getTracks: Max nr tracks reached !") ;
					nTracks = maxTracks  ;
					return 1 ;
					
				}

				//
				//     Initialize variables before going into track hit loop
				//
				FtfTrack *thisTrack = &track[nTracks-1];
				thisTrack->para     = &para ;
				thisTrack->id       = nTracks ;
				thisTrack->firstHit = thisTrack->lastHit = firstHit ;
				thisTrack->innerMostRow = thisTrack->outerMostRow = firstHit->row ;
				thisTrack->xRefHit  = firstHit->x ;
				thisTrack->yRefHit  = firstHit->y ;
				thisTrack->xLastHit = firstHit->x ;
				thisTrack->yLastHit = firstHit->y ;
				
				#ifdef TRDEBUG
				thisTrack->debugNew ( ) ;
				
				#endif
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
					if ( para.primaries &&
					para.mergePrimaries == 1 &&
					para.fillTracks &&
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
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
                
            }

            
        }



        //
        //       End loop over hits inside row
        //
        //       End loop over rows
        //
        //    Check time
        //
        // if ( CpuTime() - initialCpuTime > para.maxTime )
        // {
        // LOG(ERR, "online_tracking_subsector::getTracks: tracker time out after %f\n",
        // CpuTime() - initialCpuTime ) ;
        // break ;
        // }
        
    }      

    //
    para.nHitsForSegment = nHitsSegment ;  
    //
    return 0 ;
    
}

int online_tracking_subsector::getTracks ( ) 

{


    int pval = PROFILER(0) ;

    //
    //     Set conformal coordinates if we are working with primaries
    //
    int nHitsSegment   = (short)para.nHitsForSegment;  
    if ( para.primaries ) 
    
    {

        setConformalCoordinates ( ) ;
        para.minHitsForFit = 1 ;
        para.nHitsForSegment = max(2,nHitsSegment);


    }

    else 
    
    {

        para.minHitsForFit = 2 ;
        para.nHitsForSegment = max(3,nHitsSegment);
        
    }

    pval = PROFILER(pval) ;	// 70 us

    int build_yes = 0 ;
    int delete_1_yes = 0 ;
    int delete_2_yes = 0 ;

    //
    //               Loop over rows
    //
    for ( int ir = para.nRowsPlusOne - 1 ; ir>=para.minHitsPerTrack ; ir--) 
    
	{

        //
        //           Loop over hits in this particular row
        //
        if ( rowC[ir].first &&  (((FtfHit *)rowC[ir].first)->row) < para.rowEnd ) 
        break ;
        //if(rowC[ir].first!=0)
        //printf("%d  %f %f %f\n",((FtfHit *)rowC[ir].first)->x,((FtfHit *)rowC[ir].first)->y,((FtfHit *)rowC[ir].first)->z,ir);
        //printf("  %f %f %f in %d\n",((FtfHit *)rowC[ir].first)->x,((FtfHit *)rowC[ir].first)->y,((FtfHit *)rowC[ir].first)->z,ir);
        //    if ( (((FtfHit *)rowC[ir].first)->row) < para.rowEnd ) break ;
        for ( FtfHit *firstHit = (FtfHit *)rowC[ir].first ;
        firstHit != 0 ;
        firstHit = (FtfHit *)(firstHit->nextRowHit) ) 
        
		{

		int ival= PROFILER(0) ;

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

                LOG(ERR, "\n online_tracking_subsector::getTracks: Max nr tracks reached !") ;
                nTracks = maxTracks  ;
                return 1 ;
                
            }

            //
            //     Initialize variables before going into track hit loop
            //
            FtfTrack *thisTrack = &track[nTracks-1];
            thisTrack->para     = &para ;
            thisTrack->id       = nTracks ;
            thisTrack->firstHit = thisTrack->lastHit = firstHit ;
            thisTrack->innerMostRow = thisTrack->outerMostRow = firstHit->row ;
            thisTrack->xRefHit  = firstHit->x ;
            thisTrack->yRefHit  = firstHit->y ;
            thisTrack->xLastHit = firstHit->x ;
            thisTrack->yLastHit = firstHit->y ;
            
            #ifdef TRDEBUG
            thisTrack->debugNew ( ) ;
            
            #endif
            //
            //              Set fit parameters to zero
            //
            thisTrack->reset ( ) ;
            //
            //      Go into hit looking loop
            //
	    ival = PROFILER(ival) ;	// 1us


	    int iret = thisTrack->buildTrack ( firstHit, volumeC ) ;

	    ival = PROFILER(ival) ;

//            if ( thisTrack->buildTrack ( firstHit, volumeC ) ) 
	    if(iret) 
            {
		build_yes++ ;
                //
                //    Merge Tracks if requested
                //
                if ( para.primaries &&
                para.mergePrimaries == 1 &&
                para.fillTracks &&
                thisTrack->mergePrimary( trackC )  ) 
                
                {

		    delete_1_yes++ ;
                    nTracks-- ;
                    thisTrack->deleteCandidate ( ) ;
                    
                }

                
            }

            else
            
            {	
		delete_2_yes++ ;

                //
                //      If track was not built delete candidate
                //
		int jval = PROFILER(0) ;
                thisTrack->deleteCandidate ( ) ;
                nTracks-- ;
		PROFILER(jval) ;
                
            }

	    ival = PROFILER(ival) ;	// 4 us
            //
            //       End loop over hits inside row
            //
            
        }

        //       End loop over rows
        //
        //    Check time
        //
        if ( CpuTime() - initialCpuTime > para.maxTime ) 
        
        {

            LOG(ERR, "online_tracking_subsector::getTracks: tracker time out after %f\n", 
            CpuTime() - initialCpuTime ) ;
            break ;
            
        }

        
    }

    pval = PROFILER(pval) ;
//    LOG(TERR,"nTracks %d: %d %d %d",nTracks, build_yes, delete_1_yes, delete_2_yes) ;
    //
    para.nHitsForSegment = nHitsSegment ;  
    //
    return 0 ;
    
}

//********************************************************************
//
void online_tracking_subsector::mergePrimaryTracks ( ) 

{

    //
    //   Reset area keeping track pointers
    //
    memset ( trackC, 0, para.nPhiTrackPlusOne*para.nEtaTrackPlusOne*sizeof(FtfContainer) ) ;  
    //
    //    Loop over tracks
    //
    for ( int i = 0 ; i < nTracks ; i++ ) 
    
    {

        currentTrack = &(track[i]);
        if ( currentTrack->flag < 0 ) continue ;
        //
        //  reset link to following track
        //
        currentTrack->nxatrk = 0 ;
        //
        //    Try to merge this track
        //    if track is not merged is added
        //    to the track volume (area)
        //
        if ( currentTrack->mergePrimary ( trackC ) ) 
        
        {

            currentTrack->flag = -1 ;
            
        }

        
    }

    
}

//********************************************************************
//      Resets program
//*********************************************************************
int online_tracking_subsector::reset (void) 

{

    double phiDiff ;
    //
    //   Initialization flag in principle assume failure
    //
    para.init = 0 ;
    //----------------------------------------------------------------------------
    //     Allocate volumes
    //---------------------------------------------------------------------------*/
    para.nRowsPlusOne = ( para.rowOuterMost - para.rowInnerMost ) / para.modRow + 2 ;
    if ( para.nRowsPlusOne < 1 ) 
    
    {

        LOG(ERR, " Rows: Outer Most Inner Most %d % d \n", 
        para.rowOuterMost,  para.rowInnerMost ) ;
        return 1 ;
        
    }

    para.nPhiPlusOne    = para.nPhi + 1 ;
    para.nEtaPlusOne    = para.nEta + 1 ;
    para.nPhiEtaPlusOne = para.nPhiPlusOne * para.nEtaPlusOne ;
    if ( para.mergePrimaries ) 
    
    {

        para.nPhiTrackPlusOne = para.nPhiTrack + 1 ;
        para.nEtaTrackPlusOne = para.nEtaTrack + 1 ;
        
    }

    //
    //-->    Allocate volume memory
    //
    if (volumeC != NULL) delete []volumeC; 
    
    #ifdef TRDEBUG
    LOG(ERR, "Allocating %d bytes of memory for volume\n",
    para.nRowsPlusOne*
    para.nPhiPlusOne*
    para.nEtaPlusOne*sizeof(VOLUME));
    
    #endif
    int nVolumes = para.nRowsPlusOne*para.nPhiPlusOne *
    para.nEtaPlusOne ;
    volumeC = new FtfContainer[nVolumes];
    if(volumeC == NULL) 
    
    {

        LOG(ERR, "Problem with memory allocation... exiting\n" ) ;
        return 1 ;
        
    }

    //
    //      Allocate row memory
    //
    if ( rowC != NULL ) delete[] rowC ;
    
    #ifdef TRDEBUG
    LOG(ERR, "Allocating %d bytes of memory for rowC\n",
    para.nRowsPlusOne*sizeof(ROW));
    
    #endif
    rowC = new FtfContainer[para.nRowsPlusOne];
    if ( rowC == NULL) 
    
    {

        LOG(ERR, "Problem with memory allocation... exiting\n" ) ;
        exit(0);
        
    }

    //
    //       Allocate track area memory
    //
    if ( para.mergePrimaries ) 
    
    {

        if (trackC    != NULL) delete []trackC     ;
        
        #ifdef TRDEBUG
        LOG(ERR, "Allocating %d bytes of memory for track_area\n",
        para.nPhiTrackPlusOne*
        para.nEtaTrackPlusOne*sizeof(AREA));
        
        #endif
        int nTrackVolumes = para.nPhiTrackPlusOne*
        para.nEtaTrackPlusOne ;
        trackC    = new FtfContainer[nTrackVolumes];
        if(trackC == NULL) 
        
        {

            LOG(ERR,  "Problem with memory allocation... exiting\n" ) ;
            return 1 ;
            
        }

        else
        
        {

            //
            //   Check there is some memory allocated
            //
            if ( trackC == 0 )
            
            {

                LOG(ERR,"online_tracking_subsector::reset: Merging option not available \n " ) ; 
                LOG(ERR, " when option was not used the first time         \n " ) ; 
                return 1 ;
                
            }

            
        }

        
    }

    /*--------------------------------------------------------------------------
    Look whether the phi range is closed (< 5 degrees )
    -------------------------------------------------------------------------- */
    phiDiff = para.phiMax - para.phiMin ;
    if ( phiDiff > 2. * pi + 0.1 ) 
    
    {

        LOG(ERR, "online_tracking_subsector::reset: Wrong phi range %f, %f ", 
        para.phiMin*toDeg, para.phiMax*toDeg ) ;
        return 1 ;
        
    }

    if ( fabs(phiDiff-twoPi ) < pi / 36. ) para.phiClosed = 1 ;
    else 
    para.phiClosed = 0 ;
    /*--------------------------------------------------------------------------
    Calculate volume dimensions
    -------------------------------------------------------------------------- */
    para.phiSlice   = (para.phiMax - para.phiMin)/para.nPhi ;
    para.etaSlice   = (para.etaMax - para.etaMin)/para.nEta ;
    /*--------------------------------------------------------------------------
    If needed calculate track area dimensions
    -------------------------------------------------------------------------- */
    para.phiSliceTrack   = (para.phiMaxTrack - para.phiMinTrack)/para.nPhiTrack ;
    para.etaSliceTrack   = (para.etaMaxTrack - para.etaMinTrack)/para.nEtaTrack ;
    //
    //    Set vertex parameters
    //
    if ( para.xVertex != 0 || para.yVertex != 0 ) 
    
    { 

        para.rVertex   = (double)sqrt (para.xVertex*para.xVertex +
        para.yVertex*para.yVertex) ;
        para.phiVertex = (double)atan2(para.yVertex,para.xVertex);
        
    }

    else 
    
    {

        para.rVertex   = 0.F ;
        para.phiVertex = 0.F ;
        
    }

    if ( para.dxVertex != 0 || para.dyVertex != 0 )
    para.xyWeightVertex = 1.F / ((double)sqrt(para.dxVertex*para.dxVertex+
    para.dyVertex*para.dyVertex) ) ;
    else para.xyWeightVertex = 1.0F ;
    //
    //   Set # hits & tracks to zero
    //
    //nHits   = 0 ;
    // nTracks = 0 ;
    //
    //    Set initialization flag to true
    //
    para.init = 1 ;
    return 0 ;
    
}

//*********************************************************************
//	Set hit pointers
//*********************************************************************
int online_tracking_subsector::setConformalCoordinates ( )

{

    /*-------------------------------------------------------------------------
    Loop over hits 
    -------------------------------------------------------------------------*/
    FtfHit* thisHit ;
    double x, y, r2, invR2 ;
    for ( int ihit = 0 ; ihit<nHits ; ihit++ )
    
    {

        /*-------------------------------------------------------------------------
        Transform coordinates
        -------------------------------------------------------------------------*/
        thisHit = &(hit[ihit]) ;
#if 0	// Tonko: bad code

        uint v1 = (uint)volumeC;
        uint v2 = (uint)&volumeC[20746];
        uint h = (uint)thisHit;
        if((h>v1) && (h<v2)) 
        
        {

            printf("hit: 0x%x v1=0x%x v2=0x%x ihit=%d nHits=%d hit=0x%x\n",
            h,v1,v2,ihit,nHits,(uint)hit);
            
        }
#endif

        x            = thisHit->x - para.xVertex ;
        y            = thisHit->y - para.yVertex ;
        r2           = x * x + y * y ;
        invR2        = 1.F / r2 ;
        thisHit->xp    =     x * invR2 ;
        thisHit->yp    =   - y * invR2 ;
        thisHit->wxy   =   r2 * r2 /  ( square(para.xyErrorScale)
        * ( square(thisHit->dx) + square(thisHit->dy) ) ) ;
        
    } 

    return 0 ;
    
} 

//********************************************************************
//	Set hit pointers
//********************************************************************
int online_tracking_subsector::setPointers ( )

{

    int ihit, localRow ;
    register int volumeIndex;
    double r, r2, phi, eta ;
    FtfHit *thisHit ;
    //
    nHitsOutOfRange = 0 ;
    //
    //   Set volumes to zero
    //
    memset ( rowC,   0, para.nRowsPlusOne*sizeof(FtfContainer) ) ;
    int n = para.nRowsPlusOne*para.nEtaPlusOne*para.nPhiPlusOne ;
    memset ( volumeC, 0, n*sizeof(FtfContainer) ) ;
    if ( para.mergePrimaries )
    
    { 

        memset ( trackC, 0, para.nPhiTrackPlusOne*para.nEtaTrackPlusOne*sizeof(FtfContainer) ) ;  
        
    }

    /*-------------------------------------------------------------------------
    Loop over hits 
    -------------------------------------------------------------------------*/
    for ( ihit = 0 ; ihit<nHits ; ihit++ )
    
    {

        /*-------------------------------------------------------------------------
        Check whether row is to be considered
        -------------------------------------------------------------------------*/
        thisHit = &(hit[ihit]) ;
        localRow = thisHit->row - para.rowInnerMost ;
        if ( fmod((double)localRow,(double)para.modRow) != 0 ) continue ;
        if ( thisHit->row < para.rowInnerMost ) continue ;
        if ( thisHit->row > para.rowOuterMost ) continue ;
        /*
        *->    Get "renormalized" index
        */
        localRow = localRow / para.modRow + 1 ;
        /*-------------------------------------------------------------------------
        Transform coordinates
        -------------------------------------------------------------------------*/
        r2            = thisHit->x * thisHit->x + thisHit->y * thisHit->y ;
        r             = (double)sqrt ( r2 ) ;
        phi           = (double)atan2(thisHit->y,thisHit->x) + para.phiShift ;
        if ( phi < 0 ) phi = phi + twoPi ;
        // l3Log("r: %f, z: %f\n",r, thisHit->z);
        eta           = (double)seta(r,thisHit->z) ;
        if ( para.szFitFlag ) 
        
        {

            thisHit->s  = 0.F ;
            thisHit->wz = (double)(1./ square ( para.szErrorScale * thisHit->dz ));
            
        }

        thisHit->r   = r   ;
        //if(thisHit->row==20)printf("Xiangming Sun r=  %f\n",r);
        thisHit->phi = phi ;
        thisHit->eta = eta ;
        /*-------------------------------------------------------------------------
        Set pointers
        -------------------------------------------------------------------------*/
        thisHit->nextVolumeHit  = 
        thisHit->nextRowHit     = 0 ;
        /*-------------------------------------------------------------------------
        Get phi index for hit
        -------------------------------------------------------------------------*/
        thisHit->phiIndex = (int)( (thisHit->phi-para.phiMin)/para.phiSlice + 1.);
        if ( thisHit->phiIndex < 1 || thisHit->phiIndex > para.nPhi ) 
        
        {

            if ( para.infoLevel > 2 ) 
            
            {

                LOG(ERR, " === > Hit %d has Phi = %f \n", (int)thisHit->id,    
                thisHit->phi*toDeg ) ;
                LOG(ERR, " Phi index %d out of range \n", thisHit->phiIndex ) ;
                
            }

            nHitsOutOfRange++ ;
            continue ;
            
        } 

        /*-------------------------------------------------------------------------
        Get eta index for hit
        -------------------------------------------------------------------------*/
        thisHit->etaIndex = (int)((thisHit->eta - para.etaMin)/para.etaSlice + 1.);
        if ( thisHit->etaIndex < 1 || thisHit->etaIndex > para.nEta ) 
        
        {

            if ( para.infoLevel > 2 ) 
            
            {

                LOG(ERR, " \n === > Hit %d has Eta = %f  z %f ", (int)thisHit->id, 
                thisHit->eta, thisHit->z ) ;
                LOG(ERR, " \n Eta min/max %f %f ", para.etaMin, para.etaMax ) ;
                LOG(ERR, " \n Eta slice   %f    ", para.etaSlice ) ;
                LOG(ERR, " \n Eta index %d out of range ", thisHit->etaIndex ) ;	  
                
            }

            nHitsOutOfRange++ ;
            continue ;
            
        }

        //
        //    Reset track assigment
        //
        thisHit->nextTrackHit  = 0 ;
        thisHit->track         = 0 ;
        /* ------------------------------------------------------------------------- 
        Increase nr. of hits in small volume  WARNING! C-arrays go from 0
        Set Id of first hit in this vol. and link to next hit in previous
        hit in the same volume
        -------------------------------------------------------------------------*/
        volumeIndex = localRow  * para.nPhiEtaPlusOne + 
        thisHit->phiIndex * para.nEtaPlusOne + thisHit->etaIndex ;
        if (volumeC[volumeIndex].first == 0 ) 
        volumeC[volumeIndex].first = (void *)thisHit ;
        else
        ((FtfHit *)(volumeC[volumeIndex].last))->nextVolumeHit = thisHit ;
        volumeC[volumeIndex].last = (void *)thisHit ;
        /*-------------------------------------------------------------------------
        Set row pointers
        -------------------------------------------------------------------------*/
        if ( rowC[localRow].first == NULL )
        
        {

            rowC [localRow].first = (void *)thisHit ;
            
        }

        else
        ((FtfHit *)(rowC[localRow].last))->nextRowHit = thisHit ;
        rowC[localRow].last = (void *)thisHit ;
        
    }

    return 0 ;
    
} 

//***********************************************************************
//     For timing
//***********************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
double online_tracking_subsector::CpuTime( void ) 

{

    return (double)(clock())/CLOCKS_PER_SEC;
    
}
int online_tracking_subsector::timer_run()
{

    ovalue.it_interval.tv_sec = 0;
    ovalue.it_interval.tv_nsec = 0;
    ovalue.it_value.tv_sec = MILLION; /* a large number */
    ovalue.it_value.tv_nsec = 0;
    
    if (timer_settime(timeid, 0, &ovalue, NULL) == -1) {
        perror("Failed to set interval timer"); 
        return -1;
    }
    
    return 0;
}

int online_tracking_subsector::timer_stop()
{
    if (timer_gettime(timeid, &nvalue) == -1) {
        perror("Failed to get interval timer value");
        return -1;
    }
    diftime = MILLION*(ovalue.it_value.tv_sec - nvalue.it_value.tv_sec) +
            (ovalue.it_value.tv_nsec - nvalue.it_value.tv_nsec)/THOUSAND;    
    return 0;
}

void online_tracking_subsector::timer_print(char *name)
{

    printf("The %s took %ld microseconds or %f seconds.\n", name,
            diftime, diftime/(double)MILLION);
}
//

#ifdef I386
double online_tracking_subsector::RealTime (void) 

{

    const long nClicks = 400000000 ;
    unsigned long eax, edx;
    asm volatile("rdtsc":"=a" (eax), "=d" (edx));
    double realTime = (double)(eax)/ nClicks;
    return realTime;
    
}


#else
double online_tracking_subsector::RealTime (void) 

{

    return 1. ;
    
}


#endif

