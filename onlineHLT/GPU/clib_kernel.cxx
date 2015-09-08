
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>



#include "clib_kernel.h"

#define fitem 6
#define iitem 4
#define titem 20





blocks blockIdx;
blocks blockDim;
blocks threadIdx;
int __mul24(int a,int b){
return a*b;
}





/*__global__*/ void tracking(int number,float *para,float *fhit,int *ihit,int * volumeC,int * rowC,float *track){



int nTracks=0;

	for ( int ir = df_nRow ; ir>=df_minHitsPerTrack ; ir--) 
    
	{


		 

		

        for ( int first_index=rowC[2*ir]; 
        first_index != -1 ;
        first_index=ihit[iitem*first_index+1] ) 
        
		{

		

		
		
			//clib_FtfTrack  body_thisTrack ;
			//clib_FtfTrack *thisTrack=&body_thisTrack;
			//thisTrack->fhit=fhit;
			//thisTrack->ihit=ihit;
			//thisTrack->track=track;
			//thisTrack->id= nTracks;
			//
			//clib_FtfHit body_firstHit;
			//clib_FtfHit* firstHit=&body_firstHit;
			//firstHit->setup(first_index,fhit);
			//
			//thisTrack->firstHit = thisTrack->lastHit = (void *)first_index ;
			//thisTrack->innerMostRow = thisTrack->outerMostRow = ihit[iitem*first_index+2] ;
   //         thisTrack->xRefHit  = firstHit->x ;
   //         thisTrack->yRefHit  = firstHit->y ;
   //         thisTrack->xLastHit = firstHit->x ;
   //         thisTrack->yLastHit = firstHit->y ;	
			//thisTrack->reset ( ) ;
			
            // if ( thisTrack->buildTrack ( firstHit, volumeC ) ) 
            
            // {


                
            // }

            // else
            
            // {

                // thisTrack->deleteCandidate ( ) ;
                // nTracks-- ;
                
            // }
			
			
			



		}
	}







}


/*__global__*/ void class_hash(int *ihit,int number,int * volumeC,int * rowC){

int il;

for(il=0;il<number;il++){


// if(ihit[iitem*il+0]>df_nPhi||ihit[iitem*il+0]<1) continue;
// if(ihit[iitem*il+1]>df_nEta||ihit[iitem*il+1]<1) continue;


int vi=ihit[iitem*il+2];
 if(vi<0) continue;

int localrow=ihit[iitem*il+2]/(41*11);



int pindex=(vi/41)%11;
int eindex=vi%41;






if(volumeC[2*vi]==-1)
volumeC[2*vi]=il;
else
ihit[iitem*(volumeC[2*vi+1])+0]=il;
volumeC[2*vi+1]=il;

if(rowC[2*localrow]==-1)
rowC[2*localrow]=il;
else
ihit[iitem*(rowC[2*localrow+1])+1]=il;
rowC[2*localrow+1]=il;
}

}




/*__global__*/ void
clear(int * dataC,int number){
int index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
if(index>=number) return;
dataC[2*index+0]=-1;
dataC[2*index+1]=-1;
}

/*__global__*/ void
calc_hash(float *fhit,int *ihit,float *para,int number,int * volumeC,int * rowC)
{
// int index = blockIdx.x*blockDim.x + threadIdx.x;
int index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

if(index>=number) return;





float thisHit_x=fhit[fitem*index+0];
float thisHit_y=fhit[fitem*index+1];
float thisHit_z=fhit[fitem*index+2];
float thisHit_row=fhit[fitem*index+3];


 float       r2            = thisHit_x * thisHit_x + thisHit_y * thisHit_y ;
 float       r             = (float)sqrt ( r2 ) ;
 float       phi           = (float)atan2(thisHit_y,thisHit_x) + para[3] ;
        if ( phi < 0 ) phi = phi + twoPi ;

 float       eta           = (float)seta(r,thisHit_z) ;
 float       phiSlice=(para[2]-para[1])/df_nPhi;
 int phiIndex = (int)( (phi-para[1])/phiSlice + 1.);
 float etaSlice =(para[5]-para[4])/df_nEta;
 int etaIndex = (int)((eta - para[4])/etaSlice + 1.);
 
 
     int nPhiPlusOne    = df_nPhi + 1 ;
    int nEtaPlusOne    = df_nEta + 1 ;
    int nPhiEtaPlusOne = nPhiPlusOne * nEtaPlusOne ;
 
 int localRow = int(thisHit_row) - df_rowInnerMost+1 ;
 
 
 int volumeIndex = localRow  * nPhiEtaPlusOne + phiIndex * nEtaPlusOne + etaIndex ;
 
 
 
 
 if(phiIndex>df_nPhi||phiIndex<1) volumeIndex=-1;
 if(etaIndex>df_nEta||etaIndex<1) volumeIndex=-1;
 
 
 ihit[iitem*index+0]=-1;
 ihit[iitem*index+1]=-1;
 ihit[iitem*index+2]=volumeIndex;
 ihit[iitem*index+3]=-1;
 
 float        x            = thisHit_x - df_xVertex ;
 float       y            = thisHit_y - df_yVertex ;
        r2           = x * x + y * y ;
 float       invR2        = 1.F / r2 ;
        fhit[fitem*index+3]    =     x * invR2 ;
        fhit[fitem*index+4]    =   - y * invR2 ;
        fhit[fitem*index+5]   =   r2 * r2 /  ( square(df_xyErrorScale)
        * ( square(x * invR2) + square(- y * invR2) ) ) ;
 
 
 
 
}



/*__global__*/ void
test(float* data)
{
 int index = blockIdx.x*blockDim.x + threadIdx.x;


 /*
 clib_FtfPara para;
 para.setDefaults();
 if(index==0){data[index]=para.modRow;}
 if(index==1){data[index]=para.infoLevel;}
 if(index==2){data[index]=para.hitChi2Cut;}
 if(index==3){data[index]=para.goodHitChi2;}
 if(index==4){data[index]=para.trackChi2Cut;}
 if(index==5){data[index]=para.maxChi2Primary;}
 if(index==6){data[index]=para.segmentRowSearchRange;}
 if(index==7){data[index]=para.trackRowSearchRange;}
 if(index==8){data[index]=para.dEdx;}
 if(index==9){data[index]=para.dEdxNTruncate;}
 */

 
data[index]=dtest(data);
}


/*__device__*/ float
dtest(float *data)
{
float re=0;
int i;
 int length=1000;
 float buf[1000];
 for(i=0;i<length;i++){
buf[i]=2;
 }
 for(i=0;i<length;i++){
re+=buf[i];
 }




return re;
}

void clib_FtfPara::prt(){

  cout<<"       infoLevel    "<< infoLevel<<endl;       // Level of information prcout<<"ed about progress
  cout<<"       segmentRowSearchRange    "<< segmentRowSearchRange<<endl;       // Row search range for segments 
  cout<<"       trackRowSearchRange    "<<trackRowSearchRange <<endl;         // Row search range for tracks 
  cout<<"       dEdx      "<<dEdx <<endl;          // dEdx switch
  cout<<"       dEdxNTruncate     "<<dEdxNTruncate <<endl;  // # pocout<<"s to truncate in dEdx
  cout<<"       minHitsForDedx    "<<minHitsForDedx <<endl;  // cs: min number of hits for dEdx calculation
  cout<<"       eventReset       "<<eventReset <<endl;   // Flag to reset event in fft 
  cout<<"       getErrors        "<< getErrors<<endl;   // Flag to switch error calculation
  cout<<"       fillTracks       "<< fillTracks<<endl;   // Flag to switch FtfTrack class filling
  cout<<"       ghostFlag        "<< ghostFlag<<endl;   // =1 when there are ghost hits
  cout<<"       goBackwards      "<< goBackwards<<endl;   // Flag to go backwards at the end of track reco
  cout<<"       init    "<< init<<endl;            // Control initialization 
  cout<<"       mergePrimaries     "<< mergePrimaries<<endl; // Switch to control primary merging 
  cout<<"       minHitsPerTrack    "<<minHitsPerTrack <<endl; // Minimum # hits per track 
  cout<<"       modRow    "<< modRow<<endl;          // Modulo pad row number to use 
  cout<<"       nHitsForSegment    "<<nHitsForSegment <<endl; // # hits in initial segments 
  cout<<"       minHitsForFit    "<<minHitsForFit <<endl;
  cout<<"       nEta    "<< nEta<<endl;            // # volumes in eta 
  cout<<"       nEtaTrack    "<< nEtaTrack<<endl;       // # Track areas in eta 
  cout<<"       nPhi    "<< nPhi<<endl;            // # volumes in nphi 
  cout<<"       nPhiTrack    "<<nPhiTrack <<endl;       // # Track areas in nphi 
  cout<<"       nPrimaryPasses    "<<nPrimaryPasses <<endl;  // # iterations looking for primaries
  cout<<"       nSecondaryPasses    "<< nSecondaryPasses<<endl;// # iterations looking for secondaries
  cout<<"       vertexConstrainedFit    "<<vertexConstrainedFit <<endl; // 
  cout<<"       parameterLocation    "<< parameterLocation<<endl; // 1=inner most pocout<<", 0=closest approach
  cout<<"     maxChi2Primary     "<<maxChi2Primary <<endl; // maximum chi2 to be considered primary 
  cout<<"       rowInnerMost    "<<rowInnerMost <<endl;    // Row where end track search 
  cout<<"       rowOuterMost    "<<rowOuterMost <<endl;    // Outer most row to consider tin tracking
  cout<<"       rowStart    "<<rowStart <<endl;        // Row where start track search
  cout<<"       rowEnd      "<<rowEnd <<endl;        // Row where end   track search
  cout<<"       szFitFlag    "<<szFitFlag <<endl;       // Switch for sz fit 
  cout<<"     bField          "<< bField<<endl;    // Magnetic field  (magnitude = >0)
  cout<<"       bFieldPolarity    "<<bFieldPolarity <<endl;  // polarity of field (1 or -1)
  cout<<"     hitChi2Cut    "<<hitChi2Cut <<endl;      // Maximum hit chi2 
  cout<<"     goodHitChi2    "<<goodHitChi2 <<endl;     // Chi2 to stop looking for next hit 
  cout<<"     trackChi2Cut    "<<trackChi2Cut <<endl;    // Maximum track chi2 
  cout<<"     deta    "<< deta<<endl;            // Eta search range 
  cout<<"     dphi    "<<dphi <<endl;            // Phi search range 
  cout<<"     detaMerge     "<<detaMerge <<endl;      // Eta difference for track merge 
  cout<<"     dphiMerge     "<<dphiMerge <<endl;      // Phi difference for track merge
  cout<<"     distanceMerge     "<<distanceMerge <<endl;  // Maximum distance for reference pocout<<" to merge secondaries
  cout<<"     etaMin    "<<etaMin <<endl;          // Min eta to consider 
  cout<<"     etaMinTrack    "<<etaMinTrack<<endl;    // Track min eta to consider 
  cout<<"     etaMax    "<< etaMax<<endl;          // Max eta to consider 
  cout<<"     etaMaxTrack     "<< etaMaxTrack<<endl;    // Track max eta to consider 
  cout<<"     goodDistance     "<< goodDistance<<endl;   // In segment building
  // distance consider good enough 
  cout<<"     phiMin    "<<phiMin <<endl;          // Min phi to consider 
  cout<<"     phiMinTrack    "<<phiMinTrack<<endl;    // Track min phi to consider 
  cout<<"     phiMax    "<<phiMax <<endl;          // Max phi to consider 
  cout<<"     phiMaxTrack     "<<phiMaxTrack <<endl;    // Track max phi to consider 
  cout<<"     phiShift          "<<phiShift <<endl;  // Shift in phi when calculating phi
  cout<<"     ptMinHelixFit     "<<ptMinHelixFit <<endl;  // Minimum pt to apply helix fit
  cout<<"     maxDistanceSegment    "<<maxDistanceSegment <<endl; // Maximum distance for segments 
  cout<<"     segmentMaxAngle    "<< segmentMaxAngle<<endl; // Maximum angle between to consecutive track pieces 
  // when forming segments. A piece is the connection 
  // two hits
  cout<<"     szErrorScale    "<< szErrorScale<<endl;    // sz error scale 
  cout<<"     xyErrorScale    "<< xyErrorScale<<endl;    // xy error scale 
  cout<<"     xVertex          "<< xVertex<<endl;   // x position primary vertex 
  cout<<"     yVertex          "<<yVertex <<endl;   // y position primary vertex 
  cout<<"     dxVertex         "<< dxVertex<<endl;
  cout<<"     dyVertex         "<<dyVertex <<endl;
  cout<<"     zVertex          "<<zVertex <<endl;
  cout<<"     xyWeightVertex    "<<xyWeightVertex <<endl;  // Weight vertex in x-y
  cout<<"     phiVertex          "<<phiVertex <<endl;
  cout<<"     rVertex            "<< rVertex<<endl;
  cout<<"     maxTime            "<<maxTime <<endl; // maxTime tracker can run
  cout<<"       phiClosed     "<< phiClosed<<endl;
  cout<<"       primaries      "<< primaries<<endl;
  cout<<"       nRowsPlusOne       "<< nRowsPlusOne<<endl; // Number volumes + 1
  cout<<"       nPhiPlusOne       "<< nPhiPlusOne<<endl; // Number volumes + 1
  cout<<"       nEtaPlusOne     "<< nEtaPlusOne<<endl; // Number volumes + 1 
  cout<<"       nPhiEtaPlusOne     "<<nPhiEtaPlusOne <<endl; // Number volumes + 1 
  cout<<"       nPhiTrackPlusOne     "<< nPhiTrackPlusOne<<endl;                
  cout<<"       nEtaTrackPlusOne     "<< nEtaTrackPlusOne<<endl;                
  cout<<"     phiSlice     "<<phiSlice <<endl;
  cout<<"     etaSlice     "<<etaSlice <<endl;
  cout<<"     phiSliceTrack     "<<phiSliceTrack <<endl;
  cout<<"     etaSliceTrack     "<<etaSliceTrack <<endl;
  cout<<"     zMax     "<< zMax<<endl;
  cout<<"		current_ntrack    "<<current_ntrack <<endl;

}

/*__device__*/ void clib_FtfPara::setDefaults (void)
{
/*  Define cuts - this should be obsolete */

   modRow          = 1    ;
   infoLevel       = 0 ;
   hitChi2Cut      = 500.F  ;
   goodHitChi2     = 100.F ;
   trackChi2Cut    = 250.F ;
   maxChi2Primary  = 0. ;
   segmentRowSearchRange = 1 ;
   trackRowSearchRange   = 3 ;
   dEdx              = 0     ;
   dEdxNTruncate     = 20    ;
   minHitsForDedx    = 15    ;
   dphi              = 0.10F * modRow ;
   deta              = 0.10F * modRow ;
   dphiMerge         = 0.02F  ;
   detaMerge         = 0.02F  ;
   distanceMerge     = 2. ;
   etaMin            = -2.5F  ;
   etaMinTrack       = -2.2F  ;
   etaMax            =  2.5F  ;
   etaMaxTrack       =  2.2F  ;
   eventReset        =  1     ;
   getErrors         =  0     ;
   fillTracks        =  1     ;
   ghostFlag         =  0     ;
   goBackwards       =  0     ;
   goodDistance      =  1.F * modRow ;
   init              =  0 ;
   mergePrimaries    =  1    ;
   parameterLocation =  1    ;
   phiMin            =  (float)(-0.000001/toDeg)  ;
   phiMinTrack       =  (float)(-0.000001/toDeg)  ;
   phiMax            = (float)(360.2/toDeg)  ;
   phiMaxTrack       = (float)(360.2/toDeg)  ;
   maxDistanceSegment = 100.F * modRow ;
   minHitsPerTrack   = 5      ;
   nHitsForSegment   = 2      ;
   nEta              = 60     ;
   nEtaTrack         = 60     ;
   nPhi              = 20     ;
   nPhiTrack         = 60     ;
   nPrimaryPasses    = 1      ;
   nSecondaryPasses  = 0      ;
   vertexConstrainedFit = 0 ;
   rowInnerMost      = 1      ;
   rowOuterMost      = 45     ;
   rowStart          = 45     ;
   rowEnd            =  1     ;
   segmentMaxAngle   = 10.F/toDeg ;
   szFitFlag         = 1      ;
   xyErrorScale      = 1.0F   ;
   szErrorScale      = 1.0F   ;
   bField            = 0.0F   ;
   phiShift          = 0.0    ;
   zMax              = 210.   ;
   
   ptMinHelixFit     = 0.F  ;
   rVertex           = 0.F    ;
   xVertex           = 0.F    ;
   yVertex           = 0.F    ;
   zVertex           = 0.F    ;
   dxVertex          = 0.005F ;
   dyVertex          = 0.005F ;
   phiVertex         = 0.F    ;
   maxTime           = 1.e18 ; // by default tracker can run as long as the age of the Universe

   return  ;
}

/*__device__*/ void DftfMatrixDiagonal ( float *h, float &h11, float &h22, float &h33 ){
	float f1, f2, f3 ;

	f1 = h[5]*h[6]-h[8]*h[1] ;
	f2 = h[4]*h[8]-h[5]*h[5] ;
	f3 = h[8]*h[0]-h[2]*h[2] ;
	h11 = float(h[8] / ( f3 - f1 * f1 / f2 )) ;

	f1 = h[2]*h[1]-h[0]*h[5] ;
	f2 = h[8]*h[0]-h[2]*h[2] ;
	f3 = h[0]*h[4]-h[1]*h[1] ;
	h22 = float(h[0] / ( f3 - f1 * f1 / f2 )) ;

	f1 = h[1]*h[5]-h[4]*h[2] ;
	f2 = h[0]*h[4]-h[1]*h[1] ;
	f3 = h[4]*h[8]-h[7]*h[7] ;
	h33 = float(h[4] / ( f3 - f1 * f1 / f2 )) ;
}



/*__device__*/ void clib_FtfHit::setStatus ( clib_FtfTrack* this_track ) {

   track = (clib_FtfBaseTrack *)this_track ;

}


/*__device__*/   void   clib_FtfHit::setup ( int i, float*fhit ) {
this->x=fhit[fitem*i];
this->y=fhit[fitem*i+1];
this->z=fhit[fitem*i+2];
this->xp=fhit[fitem*i+3];
this->yp=fhit[fitem*i+4];
this->wxy=fhit[fitem*i+5];
this->id=i;

}


clib_FtfPara * getPara(){

	return 0;

}


int small_segmentHitSelection ( small_FtfHit *baseHit, small_FtfHit *candidateHit,small_FtfTrack *tra ){


   float dx, dy, dr, d3, dangle ;
   float dphi, deta ;
   float   angle ;
  
//
//   select hit with the
//   the smallest value of d3 (defined below)
//
   dphi  = (float)fabs((baseHit->phi) - (candidateHit->phi)) ; 
   if ( dphi > pi ) dphi = (float)fabs( twoPi - dphi ) ;
   if ( dphi > getPara()->dphi && dphi < twoPi -getPara()->dphi ) return 0 ;
//
//    Make sure we want to look at the difference in eta
//
   if ( baseHit->dz < 1000. && candidateHit->dz < 1000. ){
        deta  = (float)fabs((baseHit->eta) - (candidateHit->eta)) ; 
        if ( deta > getPara()->deta ) return 0 ;
   }
   else deta = 0.F ;
  
   dr    = (float)fabs((float)(baseHit->row - candidateHit->row));
   d3    = (float)(toDeg * dr * ( dphi  + deta ) ) ;
//
//     If initial segment is longer than 2 store angle info in 
//     a1Xy and a1_sz
//

  
   if ( getPara()->nHitsForSegment > 2 && tra->nHits-1 < getPara()->nHitsForSegment ) {
	  ;
      dx = candidateHit->x - baseHit->x ;
      dy = candidateHit->y - baseHit->y ;
      angle = (float)atan2 ( dy, dx ) ;
	 
      if ( angle < 0  ) angle = angle + twoPi ;
      tra->lastXyAngle = angle ;
   }
#ifdef TRDEBUG
   if ( getPara()->trackDebug && getPara()->debugLevel >= 3 ) {
      LOG(ERR, "FtfTrack::segmentHitSelection:\n");
      LOG(ERR, "dr,dphi,deta,distance, Min distance  %7.2f %7.2f %7.2f %7.2f %7.2f\n",
                dr,dphi,deta,d3,chi2[0] ) ;
      if ( d3 < chi2[0] )
	LOG(ERR, "Best point, keep it !!!\n" );  
      else{
         LOG(ERR, "Worse than previous, reject !!\n" );
//       candidateHit->Show ( getPara()->color_transparent );
      }
      debugAsk() ;
   }
#endif
   if ( d3 < tra->chi2[0] ) {
//
//   For second hit onwards check the difference in angle 
//   between the last two track segments
//
    
      if ( tra->nHits > 1 ) {
	 dx     = candidateHit->x - baseHit->x ;
         dy     = candidateHit->y - baseHit->y ;
         angle  = (float)atan2 ( dy, dx ) ;
         if ( angle < 0  ) angle = angle + twoPi ;
	    dangle = (float)fabs ( tra->lastXyAngle - angle );
		  
	    tra->lastXyAngle = angle ;
         if ( dangle > getPara()->segmentMaxAngle ) return 0 ;
      }
//
//    Check whether this is the "closest" hit
//
      tra->chi2[0]          = d3 ;
      if ( d3 < getPara()->goodDistance ) return 2 ;
	  return 1 ;
   }
//
//    If hit does not fulfill criterai return 0
//
   return 0 ;
}


int small_followHitSelection ( small_FtfHit *baseHit, small_FtfHit *candidateHit,  int way, small_FtfTrack *tra ){
//
   float lszChi2 = 0 ;
   float lchi2 ;
   float slocal=0, deta, dphi ;
   float dx, dy, dxy, dsz, temp ;
//
//           Check delta eta 
//
//   if ( baseHit->dz < 1000. && candidateHit->dz < 1000 ){
      deta = fabs((baseHit->eta)-(candidateHit->eta)) ;
      if ( deta > getPara()->deta ) return 0 ; 
//   }
//   else deta = 0.F ;
//
//           Check delta phi
//
  dphi = fabs((baseHit->phi)-(candidateHit->phi)) ;
  if ( dphi > getPara()->dphi && dphi < twoPi-getPara()->dphi ) return 0 ;
//
//      If looking for secondaries calculate conformal coordinates
//

//
//      Calculate distance in x and y
//
   temp = (tra->a2Xy * candidateHit->xp - candidateHit->yp + tra->a1Xy) ;
   dxy  = temp * temp / ( tra->a2Xy * tra->a2Xy + 1.F ) ;
//
//    Calculate chi2
//
   lchi2    = (dxy * candidateHit->wxy) ;

   if ( lchi2 > tra->chi2[0] ) return 0 ;
//
//      Now in the sz plane
//
   if ( getPara()->szFitFlag ){
//
//        Get "s" and calculate distance hit-line
//
      dx     = baseHit->x - candidateHit->x ;
      dy     = baseHit->y - candidateHit->y ;
      slocal = baseHit->s - way * sqrt ( dx * dx + dy * dy ) ;

      temp = (tra->a2Sz * slocal - candidateHit->z + tra->a1Sz) ;
      dsz  = temp * temp / ( tra->a2Sz * tra->a2Sz + 1 ) ;
//
//              Calculate chi2
//
      lszChi2 = dsz * candidateHit->wz ;
      lchi2 += lszChi2 ;
   } 
   else {
      lszChi2 = 0.F ;
      //slocal = 0;
   }
//
//         Check whether the chi2 square is better than previous one
//
   if ( lchi2 < tra->chi2[0] ) {
      tra->chi2[0]       = (float)lchi2    ;
      tra->chi2[1]       = (float)lszChi2 ;
      
      if ( getPara()->szFitFlag  ) candidateHit->s = (float)slocal ;
//
//       if a good chi2 is found let's stop here
//
      if ( lchi2 < getPara()->goodHitChi2 ) return 2 ;

      return 1 ;
   }
//
//     Return the selected hit
//
   return 0 ;
}






void small_add ( small_FtfHit *thisHit, int way, small_FtfTrack *tra  )
{

//	cout<<"small_add  "<<thisHit->id<<endl;
//
//      Increment # hits in this track
//
  if(tra->nHits>=50) return;
  tra->hitmap[tra->nHits]=thisHit->id;

	tra->nHits++ ; 
	tra->lasthit.phi=	thisHit->phi;
	tra->lasthit.dz=	thisHit->dz;
	tra->lasthit.eta=	thisHit->eta;
	tra->lasthit.row=	thisHit->row;
	tra->lasthit.x=	thisHit->x;
	tra->lasthit.y=	thisHit->y;
	tra->lasthit.s=	thisHit->s;

	

//
//    Check whether a fit update is needed
//
  if ( tra->nHits < getPara()->minHitsForFit ) return ;
//
//    Include hit in xy fit parameter calculation
//
  
  //if(thisHit->id==3584){
  //cout<<tra->s11Xy<<endl;
  //cout<<tra->s12Xy<<endl;
  //cout<<tra->s22Xy<<endl;
  //cout<<tra->g1Xy<<endl;
  //cout<<tra->g2Xy<<endl;
  //cout<<tra->a1Xy<<endl;
  //cout<<tra->a2Xy<<endl;
  //cout<<tra->a1Sz<<endl;
  //cout<<tra->a2Sz<<endl;
  //cout<<thisHit->wxy<<endl;
  //cout<<thisHit->xp<<endl;
  //cout<<thisHit->yp<<endl;
  //cout<<thisHit->s<<endl;
  //cout<<thisHit->z<<endl;
  //cout<<tra->nHits<<endl;

  //}


  tra->s11Xy = tra->s11Xy + thisHit->wxy ;
  tra->s12Xy = tra->s12Xy + thisHit->wxy * thisHit->xp ;
  tra->s22Xy = tra->s22Xy + thisHit->wxy * square(thisHit->xp) ;
  tra->g1Xy  = tra->g1Xy  + thisHit->wxy * thisHit->yp ;
  tra->g2Xy  = tra->g2Xy  + thisHit->wxy * thisHit->xp * thisHit->yp ;
  
 
  if ( tra->nHits > getPara()->minHitsForFit  )
  {
     tra->ddXy  = tra->s11Xy * tra->s22Xy - square ( tra->s12Xy ) ;
     if ( tra->ddXy != 0 ) {
        tra->a1Xy  = ( tra->g1Xy * tra->s22Xy - tra->g2Xy * tra->s12Xy ) / tra->ddXy ;
        tra->a2Xy  = ( tra->g2Xy * tra->s11Xy - tra->g1Xy * tra->s12Xy ) / tra->ddXy ;
     }
     else {
		 //LOG(ERR, "FtfTrack:add: ddSz = 0 \n" ) ;
     }
  }
//
//     Now in the sz plane
//
  if ( getPara()->szFitFlag ) {
     tra->s11Sz = tra->s11Sz + thisHit->wz ;
     tra->s12Sz = tra->s12Sz + thisHit->wz * thisHit->s ;
     tra->s22Sz = tra->s22Sz + thisHit->wz * thisHit->s * thisHit->s ;
     tra->g1Sz  = tra->g1Sz  + thisHit->wz * thisHit->z ;
     tra->g2Sz  = tra->g2Sz  + thisHit->wz * thisHit->s * thisHit->z ;
  
     if ( tra->nHits > getPara()->minHitsForFit ) {
		
        tra->ddSz  = tra->s11Sz * tra->s22Sz -  tra->s12Sz * tra->s12Sz ;
	if ( tra->ddSz != 0 ) {
           tra->a1Sz  = ( tra->g1Sz * tra->s22Sz - tra->g2Sz * tra->s12Sz ) / tra->ddSz ;
           tra->a2Sz  = ( tra->g2Sz * tra->s11Sz - tra->g1Sz * tra->s12Sz ) / tra->ddSz ;
         }
         else
         {
            if ( getPara()->infoLevel > 0 ) {
               //LOG(ERR, "FtfTrack:add: ddSz = 0 \n" ) ;
            }
         }
      }
   }
}




void small_reset ( small_FtfTrack *tra)
{
/*----------------------------------------------------------------------
                Set fit parameters to zero
----------------------------------------------------------------------*/
 
  //tra->flag     = getPara()->primaries ;


  tra->nHits    = 0 ;
  tra->s11Xy   = 
  tra->s12Xy   = 
  tra->s22Xy   = 
  tra->g1Xy    = 
  tra->g2Xy    = 
  tra->chi2[0]  = 0.F ;
 
  //tra->nxatrk   = 0 ;
  if ( getPara()->szFitFlag ) 
  {
     tra->s11Sz =
     tra->s12Sz =
     tra->s22Sz =
     tra->g1Sz  =
     tra->g2Sz  =
     tra->chi2[1]  = 
     tra->length         = 0.F ;
  }


}


void small_track_manager::ini(){
	for(int i=0;i<2000;i++){
		this->track_id[i]=i-1;
	}
	track_id[0]=0;
	this->maxused=0;
	this->number=0;

}

int small_track_manager::pop(){

	number++;
	if(maxused<number) maxused=number;
	return track_id[number-1];

}
int small_track_manager::push(int id){
number--;
track_id[number]=id;
return number;

}



//void small_track_assign::ini(){
//
//	this->number=0;
//
//}
//
//
//int small_track_assign::push(int id){
//
//track_id[number]=id;
//number++;
//return number;
//
//}




int element_setbit(int & data,int position){
if(position>31) return -1;

unsigned short mask=0;
mask=1<<position;
data=data|mask;
return 1;

}


int element_getbit(int data,int position){

int sd=data>>position;
return sd%2;

}



int small_segmentHitGroup ( small_FtfHit *candidateHit,int num,small_FtfHit *selected_hit,int &selected_pos, small_FtfTrack *tra ){
		int result=0;
	for(int i=0;i<num;i++){
				result= small_segmentHitSelection(&tra->lasthit,&candidateHit[i],tra);
		             if ( result > 0 ) {
							selected_pos=i;
						   selected_hit->id=candidateHit[i].id;
						   selected_hit->row=candidateHit[i].row;

						   selected_hit->x=candidateHit[i].x;
						   selected_hit->y=candidateHit[i].y;
						   selected_hit->z=candidateHit[i].z;
						   selected_hit->xp=candidateHit[i].xp;
						   selected_hit->yp=candidateHit[i].yp;
						   selected_hit->eta=candidateHit[i].eta;
						   selected_hit->phi=candidateHit[i].phi;
						   selected_hit->wxy=candidateHit[i].wxy;
						   selected_hit->wz=candidateHit[i].wz;

						   selected_hit->dx=candidateHit[i].dx;
						   selected_hit->dy=candidateHit[i].dy;
						   selected_hit->dz=candidateHit[i].dz;
						   selected_hit->s=candidateHit[i].s;


                if ( result ==2  ) return 2 ; 
             }
	
	}

return result;

}


int small_followHitGroup (  small_FtfHit *candidateHit,  int way, int num,small_FtfHit *selected_hit,int &selected_pos, small_FtfTrack *tra ){


		int result=0;
	for(int i=0;i<num;i++){
				result= small_followHitSelection(&tra->lasthit,&candidateHit[i],way,tra);
		             if ( result > 0 ) {
							selected_pos=i;
						   selected_hit->id=candidateHit[i].id;
						   selected_hit->row=candidateHit[i].row;

						   selected_hit->x=candidateHit[i].x;
						   selected_hit->y=candidateHit[i].y;
						   selected_hit->z=candidateHit[i].z;
						   selected_hit->xp=candidateHit[i].xp;
						   selected_hit->yp=candidateHit[i].yp;
						   selected_hit->eta=candidateHit[i].eta;
						   selected_hit->phi=candidateHit[i].phi;
						   selected_hit->wxy=candidateHit[i].wxy;
						   selected_hit->wz=candidateHit[i].wz;

						   selected_hit->dx=candidateHit[i].dx;
						   selected_hit->dy=candidateHit[i].dy;
						   selected_hit->dz=candidateHit[i].dz;
						   selected_hit->s=candidateHit[i].s;


                if ( result ==2  ) return 2 ; 
             }
	
	}

return result;


}