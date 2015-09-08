
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "clib_kernel.h"


/*__device__*/ void clib_FtfTrack::add ( clib_FtfHit *thisHit, int way )
{
//
//      Increment # hits in this track
//
  nHits++ ; 
//
//         Update pointers
//
  if ( way < 0 || nHits == 1 ) {
     if ( nHits > 1 ) ((clib_FtfBaseHit *)lastHit)->nextTrackHit = thisHit ;
     lastHit = thisHit ;
     innerMostRow = ((clib_FtfBaseHit *)lastHit)->row ;
     xLastHit = ((clib_FtfBaseHit *)lastHit)->x ;
     yLastHit = ((clib_FtfBaseHit *)lastHit)->y ;
  }
  else {
     ((clib_FtfBaseHit *)thisHit)->nextTrackHit = firstHit ; 
     firstHit = thisHit ;
     outerMostRow = ((clib_FtfBaseHit *)firstHit)->row ;
  }
//
//        Declare hit as used and fill chi2
//
  thisHit->setStatus ( this ) ;
//
//    Check whether a fit update is needed
//

  
  if ( nHits < getPara()->minHitsForFit ) return ;
//
//    Include hit in xy fit parameter calculation
//
  
  s11Xy = s11Xy + thisHit->wxy ;
  s12Xy = s12Xy + thisHit->wxy * thisHit->xp ;
  s22Xy = s22Xy + thisHit->wxy * square(thisHit->xp) ;
  g1Xy  = g1Xy  + thisHit->wxy * thisHit->yp ;
  g2Xy  = g2Xy  + thisHit->wxy * thisHit->xp * thisHit->yp ;
  
 

  if ( nHits > getPara()->minHitsForFit  )
  {
     ddXy  = s11Xy * s22Xy - square ( s12Xy ) ;
     if ( ddXy != 0 ) {
        a1Xy  = ( g1Xy * s22Xy - g2Xy * s12Xy ) / ddXy ;
        a2Xy  = ( g2Xy * s11Xy - g1Xy * s12Xy ) / ddXy ;
		
     }
     else {
        if ( getPara()->infoLevel > 0 ) {
           //LOG(ERR, "clib_FtfTrack:add: ddXy = 0 \n" ) ;
        }
     }
  }
//
//     Now in the sz plane
//
  if ( getPara()->szFitFlag ) {
     s11Sz = s11Sz + thisHit->wz ;
     s12Sz = s12Sz + thisHit->wz * thisHit->s ;
     s22Sz = s22Sz + thisHit->wz * thisHit->s * thisHit->s ;
     g1Sz  = g1Sz  + thisHit->wz * thisHit->z ;
     g2Sz  = g2Sz  + thisHit->wz * thisHit->s * thisHit->z ;
  
     if ( nHits > getPara()->minHitsForFit ) {
		
        ddSz  = s11Sz * s22Sz -  s12Sz * s12Sz ;
	if ( ddSz != 0 ) {
           a1Sz  = ( g1Sz * s22Sz - g2Sz * s12Sz ) / ddSz ;
           a2Sz  = ( g2Sz * s11Sz - g1Sz * s12Sz ) / ddSz ;
         }
         else
         {
            if ( getPara()->infoLevel > 0 ) {
               //LOG(ERR, "clib_FtfTrack:add: ddSz = 0 \n" ) ;
            }
         }
      }
   }
}
//****************************************************************************
//   Fill track information tables
//****************************************************************************
/*__device__*/ void clib_FtfTrack::add ( clib_FtfTrack *piece ) 
{
//
//   Get circle parameters
//
  s11Xy += piece->s11Xy  ;
  s12Xy += piece->s12Xy  ;
  s22Xy += piece->s22Xy  ;
  g1Xy  += piece->g1Xy   ;
  g2Xy  += piece->g2Xy   ;

  ddXy  =   s11Xy * s22Xy - square ( s12Xy ) ;
  a1Xy  = ( g1Xy * s22Xy - g2Xy * s12Xy ) / ddXy ;
  a2Xy  = ( g2Xy * s11Xy - g1Xy * s12Xy ) / ddXy ;
//
//     Now in the sz plane
//
  if ( getPara()->szFitFlag ) {
     float det1 = s11Sz * s22Sz - s12Sz * s12Sz ;
     dtanl = (float) ( s11Sz / det1 );
     dz0   = (float) ( s22Sz / det1 );

     float det2 = piece->s11Sz * piece->s22Sz - piece->s12Sz * piece->s12Sz ;
     piece->dtanl = (float) ( piece->s11Sz / det2 );
     piece->dz0   = (float) ( piece->s22Sz / det2 );
      
     float weight1 = 1./(dtanl*dtanl);
     float weight2 = 1./(piece->dtanl*piece->dtanl);
     float weight  = (weight1+weight2);
     tanl = ( weight1 * tanl + weight2 * piece->tanl ) / weight ; 

     weight1 = 1./(dz0*dz0);
     weight2 = 1./(piece->dz0*piece->dz0);
     weight  = (weight1+weight2);
     z0   = ( weight1 * z0 + weight2 * piece->z0 ) / weight ; 
   }

//
//  Add space points to first track
//
    int counter ;
    if ( piece->outerMostRow < outerMostRow ){
      if ( lastHit != NULL ) {
         counter = 0 ;
	 for ( currentHit   = piece->firstHit ; 
	       currentHit != 0 && counter < piece->nHits ;
	       currentHit  = ((clib_FtfBaseHit *)currentHit)->nextTrackHit  ) {
	    ((clib_FtfBaseHit *)currentHit)->track = this   ;
	    counter++ ;
	 }
	 ((clib_FtfBaseHit *)lastHit)->nextTrackHit = piece->firstHit ;
	 lastHit         = piece->lastHit ;
      }
      piece->firstHit = 0 ;
      innerMostRow = piece->innerMostRow ;
      xLastHit     = piece->xLastHit ;
      yLastHit     = piece->yLastHit ;
    }
    else {
       if ( piece->lastHit != NULL ) {
	  counter = 0 ;
	  for ( currentHit   = (clib_FtfHit *)piece->firstHit ; 
		currentHit != 0 && counter < piece->nHits ;
		currentHit  = ((clib_FtfBaseHit *)currentHit)->nextTrackHit  ) {
	     ((clib_FtfBaseHit *)currentHit)->track = this   ;
	     counter++;
	  }
	  ((clib_FtfBaseHit *)piece->lastHit)->nextTrackHit = firstHit ;
	  firstHit               = piece->firstHit ;
       }
       outerMostRow = piece->outerMostRow ;
       piece->firstHit = 0 ;
    }
//
//
   nHits  += piece->nHits ;
   chi2[0] += piece->chi2[0] ;
   chi2[1] += piece->chi2[1] ;
//
//   Update track parameters
//
//
   getPara()->szFitFlag = 0 ;
   if ( getPara()->fillTracks ) fill ( ) ;
   getPara()->szFitFlag = 1 ;
//
//
//   Declare track 2 not to be used
//
   piece->flag    = -1 ;
}

//****************************************************************************
//   Control how the track gets built
//****************************************************************************
/*__device__*/ int clib_FtfTrack::buildTrack ( clib_FtfHit *frstHit, clib_FtfContainer *volume ) {
//
//   Add first hit to track
//
   add ( frstHit, GO_DOWN ) ;
//
//    Try to build a segment first
//
   if ( !segment ( volume, GO_DOWN ) ) return 0 ;
//
//    If segment build go for a real track with a fit
//
   int rowToStop = getPara()->rowInnerMost ;
   if ( !follow ( volume, GO_DOWN, rowToStop ) ) return 0 ;
//
//    Now to extent track the other direction if requested
//
   if ( getPara()->goBackwards ) {
      follow ( volume, GO_UP, getPara()->rowOuterMost ) ;
   }
//
//  Fill tracks
//
    if ( getPara()->fillTracks ) 
      fill ( ) ;


   return 1 ;
}
//***************************************************************************
//   Calculates dEdx
//***************************************************************************
/*__device__*/ void clib_FtfTrack::dEdx (  ){
   int i, j ;
   clib_FtfHit *nextHit ;
   int nTruncate = max(1,
	           getPara()->dEdxNTruncate*nHits/100) ;
   nTruncate = min(nHits/2,nTruncate) ;
//
//   Define array to keep largest de's
//
   float *de = new float[nTruncate] ;
//
//    Reset
//
   dedx = 0.F ;
   memset ( de, 0, nTruncate*sizeof(float) ) ;
//
//
//
   for  ( nextHit = (clib_FtfHit *)firstHit ; 
          nextHit != 0 ;
          nextHit = (clib_FtfHit *)nextHit->nextTrackHit) { 
    
      dedx += nextHit->q ;
	 
      if ( nextHit->q < de[0] ) continue ;

      for ( i = nTruncate-1 ; i>=0 ; i-- ){
         if ( nextHit->q > de[i] ){
            for ( j=0 ; j<i ; j++ ) de[j] = de[j+1] ;
            de[i] = nextHit->q ;
            break ;
	 }
      }
   }
//
//    Subtract largest de
//
   for ( i=0 ; i<nTruncate ; i++ ) dedx -= de[i] ;
   dedx = dedx / length ;
/*   End track in required volume condition */
      
}
//*********************************************************************** 
//   Delete track candidate 
//***********************************************************************
/*__device__*/ void clib_FtfTrack::deleteCandidate(void)
{
  clib_FtfHit *curentHit = (clib_FtfHit *)firstHit ;
  clib_FtfHit *nextHit ;

  while ( curentHit != 0 )
  {
    nextHit            = (clib_FtfHit *)curentHit->nextTrackHit;
    curentHit->nextTrackHit     =  0 ;
    curentHit->xyChi2   =
    curentHit->szChi2   =  
    curentHit->s        =  0.F ;

    curentHit->setStatus ( 0 ) ;
    curentHit = nextHit;
  }
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//    Fills track variables with or without fit
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*__device__*/ void clib_FtfTrack::fill (  ) {
//
//   Get circle parameters
//
 
	
   float xc, yc ;
   float rc   = sqrt ( a2Xy * a2Xy + 1 ) / ( 2 * fabs(a1Xy) ) ;
   pt          = (float)fabs(bFactor * getPara()->bField * rc );
   float xParameters = 0.; 
   float yParameters = 0.;
//
   if ( pt > getPara()->ptMinHelixFit ) {
      float combinedChi2 = 0.5*(chi2[0]+chi2[1])/nHits ;
      if ( getPara()->primaries && combinedChi2 < getPara()->maxChi2Primary ) 
         getPara()->vertexConstrainedFit = 1 ;
      else 
         getPara()->vertexConstrainedFit = 0 ;

      fitHelix ( ) ;

//    if ( id == 303) Print(31);
      if ( getPara()->vertexConstrainedFit && getPara()->parameterLocation ){ 
	 updateToRadius ( sqrt(xLastHit*xLastHit+yLastHit*yLastHit) ) ;
      }
      else if ( !getPara()->vertexConstrainedFit && !getPara()->parameterLocation ) {
         updateToClosestApproach ( getPara()->xVertex, getPara()->yVertex, 2000. ) ;
      }

//    if ( id == 1 ) Print(30);
   }
   else{

      if ( getPara()->primaries ){ 
	    //cout << "." << flush;
         fillPrimary ( xc, yc, rc, getPara()->xVertex, getPara()->yVertex ) ;
         if ( getPara()->parameterLocation ) {// give track parameters at inner most point
	    updateToRadius ( sqrt(xLastHit*xLastHit+yLastHit*yLastHit) ) ;
         }
      }
      else { // Secondaries now
         xc = - a2Xy / ( 2. * a1Xy ) + xRefHit ;
         yc = - 1.   / ( 2. * a1Xy ) + yRefHit ;
         if ( getPara()->parameterLocation ) { // give track parameters at inner most point
	   xParameters = xLastHit ; 
	   yParameters = yLastHit ; 
	 }
	 else { // give parameters at point of closest approach
	    getClosest ( getPara()->xVertex, getPara()->yVertex,
	                 rc, xc, yc, xParameters, yParameters ) ;
	 }
         fillSecondary ( xc, yc, xParameters, yParameters ) ;
      }
//
//    Get Errors
//
      if ( getPara()->getErrors ) {
         getErrorsCircleFit (  (float)xc, (float)yc, (float)rc ) ;
         float det = s11Sz * s22Sz - s12Sz * s12Sz ;
         dtanl = (float) ( s11Sz / det );
         dz0   = (float) ( s22Sz / det );
      }
   }
}
//****************************************************************************     
//     Fill track information variables
//****************************************************************************
/*__device__*/ void clib_FtfTrack::fillPrimary (  float &xc, float &yc, float &rc,
                              float xPar, float yPar ) {
//
//   Get circle parameters
//
   xc = getPara()->xVertex - a2Xy / ( 2. * a1Xy ) ;
   yc = getPara()->yVertex - 1.   / ( 2. * a1Xy ) ;
//
//   Get track parameters
//
   float angle_vertex  = atan2 ( yPar-yc, xPar-xc ) ;
   if ( angle_vertex < 0. ) angle_vertex = angle_vertex + twoPi ;

   float dx_last    = xLastHit - xc ;
   float dy_last    = yLastHit - yc ;
   float angle_last = atan2 ( dy_last, dx_last ) ;
   if ( angle_last < 0. ) angle_last = angle_last + twoPi ;
//
//       Get the rotation
//
   float d_angle = angle_last - angle_vertex ;

// if ( d_angle >  pi ) d_angle -= twoPi  ;
   if ( d_angle < -pi ) d_angle += twoPi  ;

   q = getPara()->bFieldPolarity * ( ( d_angle < 0 ) ? 1 : -1 ) ;
   r0   = sqrt(xPar*xPar+yPar*yPar) ;
   phi0 = atan2(yPar,xPar) ;
   if ( phi0 < 0 ) phi0 += 2. * M_PI ;
   psi  = (float)(angle_vertex - getPara()->bFieldPolarity*q * 0.5F * pi) ;
   if ( psi < 0     )  psi = (float)(psi + twoPi );
   if ( psi > twoPi )  psi = (float)(psi - twoPi );
//
//      Get z parameters if needed       
//
   if ( getPara()->szFitFlag == 1 ){
      tanl = -(float)a2Sz ;
      z0   =  (float)(a1Sz + a2Sz * ( length - rc * d_angle * getPara()->bFieldPolarity*q ) );
   }
   else if ( getPara()->szFitFlag == 2 ) {
      tanl = ((clib_FtfBaseHit *)firstHit)->z /
          (float)sqrt ( ((clib_FtfBaseHit *)firstHit)->x*((clib_FtfBaseHit *)firstHit)->x + 
	                 ((clib_FtfBaseHit *)firstHit)->y*((clib_FtfBaseHit *)firstHit)->y ) ;
      z0      = 0.F ;
   }
//
//    Store some more track info
//
   eta     = seta(1.,tanl )   ;
//
//   Set primary track
//
   flag = 1 ;

}
//****************************************************************************
//   
//   Fill track information tables
//
//****************************************************************************
/*__device__*/ void clib_FtfTrack::fillSecondary ( float &xc, float &yc,
                               float xPar, float yPar )
{
/*--------------------------------------------------------------------------
     Get angles for initial and final points
------------------------------------------------------------------------------*/
   float dx1    = ((clib_FtfBaseHit *)firstHit)->x - xc ;
   float dy1    = ((clib_FtfBaseHit *)firstHit)->y - yc ;
   float angle1 = atan2 ( dy1, dx1 ) ;
   if ( angle1 < 0. ) angle1 = angle1 + twoPi ;

   float dx2    = xLastHit - xc ;
   float dy2    = yLastHit - yc ;
   float angle2 = atan2 ( dy2, dx2 ) ;
   if ( angle2 < 0. ) angle2 = angle2 + twoPi ;
/*--------------------------------------------------------------------------
     Get the rotation
------------------------------------------------------------------------------*/
   float dangle = angle2 - angle1 ;
 //  if ( dangle >  pi ) dangle =   dangle - twoPi  ;
   if ( dangle < -pi ) dangle =   dangle + twoPi  ;

   q    =  getPara()->bFieldPolarity * ( ( dangle > 0 ) ? 1 : -1 ) ;
   r0   = ((clib_FtfHit *)lastHit)->r   ;
   phi0 = ((clib_FtfHit *)lastHit)->phi ;
   psi  = (float)(angle2 - getPara()->bFieldPolarity * q * piHalf );
   if ( psi < 0     ) psi = (float)(psi + twoPi );
//
//      Get z parameters if needed       
//
   if ( getPara()->szFitFlag ){
      tanl = -(float)a2Sz ;
      z0   =  (float)(a1Sz + a2Sz * length  );
   }
   else{
      tanl = ((clib_FtfBaseHit *)firstHit)->z /
           (float)sqrt ( ((clib_FtfBaseHit *)firstHit)->x*((clib_FtfBaseHit *)firstHit)->x + 
	                  ((clib_FtfBaseHit *)firstHit)->y*((clib_FtfBaseHit *)firstHit)->y ) ;
      z0      = 0.F ;
   }
//
//-->    Store some more track info
//
   eta     = seta(1., tanl )   ;
//
//    Set primary track flag
//
   flag = 0 ;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//        Adds hits to a track chosing the closest to fit
// Arguments:
//              volume:	      volume pointer
//              way   :       which way to procede in r (negative or positive)
//              row_to_stop:  row index where to stop
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*__device__*/ int clib_FtfTrack::follow ( clib_FtfContainer *volume, int way, int ir_stop ) {

   clib_FtfHit *nextHit ;
   
   if ( way < 0 )
      nextHit = (clib_FtfHit *)lastHit ;
   else
      nextHit = (clib_FtfHit *)firstHit ; 

//
//     Define variables to keep total chi2
//
   float xyChi2 = chi2[0] ;
   float szChi2 = chi2[1] ;
   
//
//    Loop as long a a hit is found and the segment
//    is shorter than n_hit_segm
//
   while ( way * nextHit->row < way * ir_stop ) {
//
//      Select next hit
//
      chi2[0] = getPara()->hitChi2Cut ;

      nextHit = seekNextHit ( volume, nextHit, way*getPara()->trackRowSearchRange, USE_FOLLOW ) ;


//
//    Stop if nothing found
//
      if ( nextHit == 0 ) break ;
//
//   Keep total chi2
//
      float lxyChi2 = chi2[0]-chi2[1] ;
      xyChi2 += lxyChi2 ;
      nextHit->xyChi2 = lxyChi2 ;
//
//   if sz fit update track length
//
      if ( getPara()->szFitFlag  ) {
         length = nextHit->s ;
         szChi2 += chi2[1]  ;
         nextHit->szChi2 = chi2[1] ;

      }
//
//     Add hit to track
//
      add ( nextHit, way ) ;

   } // End while
//
//    Check # hits
//
   if ( nHits < getPara()->minHitsPerTrack ) return 0 ; 
//
//   Store track chi2
//
   chi2[0] = xyChi2 ;
   chi2[1] = szChi2 ;
//
//        Check total chi2
//
   float normalized_chi2 = (chi2[0]+chi2[1])/nHits ;
   if ( normalized_chi2 > getPara()->trackChi2Cut ) return 0 ;
//
   return 1 ;
}
/*******************************************************************************
        Reconstructs tracks
*********************************************************************************/
/*__device__*/ int clib_FtfTrack::followHitSelection ( clib_FtfHit *baseHit, clib_FtfHit *candidateHit, int way ){

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
   if ( getPara()->primaries == 0 ){
      float xx = candidateHit->x - xRefHit ;
      float yy = candidateHit->y - yRefHit ;
      float rr = xx * xx + yy * yy ;
      candidateHit->xp =   xx / rr ;
      candidateHit->yp = - yy / rr ;

      candidateHit->wxy  = rr * rr /
                        ( square(getPara()->xyErrorScale)  *
                        ( square(candidateHit->dx) + square(candidateHit->dy) ) ) ;
   }
//
//      Calculate distance in x and y
//
   temp = (a2Xy * candidateHit->xp - candidateHit->yp + a1Xy) ;
   dxy  = temp * temp / ( a2Xy * a2Xy + 1.F ) ;

//
//    Calculate chi2
//
   lchi2    = (dxy * candidateHit->wxy) ;

   if ( lchi2 > chi2[0] ) return 0 ;
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

      temp = (a2Sz * slocal - candidateHit->z + a1Sz) ;
      dsz  = temp * temp / ( a2Sz * a2Sz + 1 ) ;
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
   if ( lchi2 < chi2[0] ) {
      chi2[0]       = (float)lchi2    ;
      chi2[1]       = (float)lszChi2 ;
      
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
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//    Merges tracks
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*__device__*/ int clib_FtfTrack::mergePrimary ( clib_FtfContainer *trackArea ){
   short  track_merged ;
   register int  areaIndex ;
   int    i_phi, i_eta ;
   clib_FtfTrack *i_track = 0 ;
   int    ip, ie ;
   float  delta_psi ;
//
//   Check Track is primary
//
   if ( flag != 1 ) return 0 ;
//-
//   Get track area       
//
   i_phi = (int)(( psi - getPara()->phiMinTrack ) / getPara()->phiSliceTrack + 1 );
   if ( i_phi < 0 ) {
       //LOG(ERR, " Track phi index too low  %d \n", i_phi ) ;
       i_phi = 1 ;
   }
   if ( i_phi >= getPara()->nPhiTrackPlusOne ) {
       //LOG(ERR, " Track phi index too high %d \n", i_phi ) ;
       i_phi = getPara()->nPhiTrack ;
   }
//
//     Now eta
//
   i_eta = (int)(( eta - getPara()->etaMinTrack ) / getPara()->etaSliceTrack + 1 );
   if ( i_eta <= 0 ) {
       //LOG(ERR, " Track eta index too low  %d \n", i_eta ) ;
       i_eta = 1 ;
   }
   if ( i_eta >= getPara()->nEtaTrackPlusOne ) {
       //LOG(ERR, " Track eta index too high %d \n", i_eta ) ;
       i_eta = getPara()->nEtaTrack ;
   }
//
//     Loop around selected area
//
   track_merged = 0 ;
   for ( ip = max(i_phi-1,1) ; ip < min(i_phi+2,getPara()->nPhiTrackPlusOne) ; ip++ ) {
      for ( ie = max(i_eta-1,1) ; ie < min(i_eta+2,getPara()->nEtaTrackPlusOne) ; ie++ ) {
         areaIndex = ip * getPara()->nEtaTrackPlusOne + ie ;
//
//    Loop over tracks
//
         for ( i_track = (clib_FtfTrack *)trackArea[areaIndex].first ; 
               i_track != 0 ;
               i_track = i_track->getNextTrack()  ) {
//
//    Reject track if it is not good
//
         if ( i_track->flag < 0 ) continue ; 
//
// Compare both tracks
//
//   No overlapping tracks
			short delta1 = i_track->outerMostRow - outerMostRow ;
			short delta2 = i_track->innerMostRow - innerMostRow ;
			if ( delta1 * delta2 <= 0 ) continue ;
//
//    Tracks close enough
//
            if ( fabs(eta-i_track->eta) > getPara()->detaMerge ) continue ;
            delta_psi = (float)fabs(psi - i_track->psi) ;
            if ( delta_psi > getPara()->dphiMerge && delta_psi < twoPi - getPara()->dphiMerge ) continue ;

            i_track->add ( this ) ;

            track_merged = 1 ;
            break ;
         }
      }
   }
//
//->  If track not matched add it
//
   if ( track_merged == 0 ) {
      areaIndex = i_phi * getPara()->nEtaTrackPlusOne + i_eta ;
      if ( trackArea[areaIndex].first == 0 )
         trackArea[areaIndex].first = 
         trackArea[areaIndex].last = (void *)this  ;
      else {
         ((clib_FtfTrack *)trackArea[areaIndex].last)->nxatrk = this ; 
	 trackArea[areaIndex].last = (void *)this ;
      }
   }
   return track_merged ;
}
/************************************************************************* 
	Recontruct primary tracks 
*************************************************************************/
/*__device__*/ void clib_FtfTrack::reset (void)
{
/*----------------------------------------------------------------------
                Set fit parameters to zero
----------------------------------------------------------------------*/

  flag     = getPara()->primaries ;
  nHits    = 0 ;
  s11Xy   = 
  s12Xy   = 
  s22Xy   = 
  g1Xy    = 
  g2Xy    = 
  chi2[0]  = 0.F ;
  nxatrk   = 0 ;
  if ( getPara()->szFitFlag ) 
  {
     s11Sz =
     s12Sz =
     s22Sz =
     g1Sz  =
     g2Sz  =
     chi2[1]  = 
     length         = 0.F ;
  }
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     Function to look for next hit
// Input:	volume:         Volume pointer
//          baseHit:       Last point in track
//          n_r_steps:      How many rows search and which way (up or down)
//		    which_function: Function to be used to decide whether the hit is good
// Returns:	Selected hit
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*__device__*/ clib_FtfHit *clib_FtfTrack::seekNextHit ( clib_FtfContainer  *volume, 
			        clib_FtfHit *baseHit,
				int     n_r_steps,
				int which_function ) {
#define N_LOOP 9 
   int loop_eta[N_LOOP] = { 0, 0, 0,-1,-1,-1, 1, 1, 1 } ;
   int loop_phi[N_LOOP] = { 0,-1, 1, 0,-1, 1, 0,-1, 1 };


   int ir, irp, ipp, itp, k;
   register int areaIndex ; 
   int result ;
   
//-------------------------------------------------------------------------------
//     Calculate limits on the volume loop
//-----------------------------------------------------------------------------*/
   int initialRow, way ;
   if ( n_r_steps < 0 ) {
      initialRow = max(1, (baseHit->row - getPara()->rowInnerMost)/getPara()->modRow);
      n_r_steps  = min(initialRow,-n_r_steps ) ;
      way        = -1 ;
   }
   else {
      initialRow = max(1, (baseHit->row - getPara()->rowInnerMost + 2)/getPara()->modRow);
      n_r_steps  = min((getPara()->rowOuterMost-initialRow+1),n_r_steps) ;
      way = 1 ;
   }
   
   clib_FtfHit *selected_hit  = 0 ;
//
//      Loop over modules
//
   for ( ir = 0 ; ir < n_r_steps ; ir++ ){
      irp = initialRow + way * ir ;
      for ( k=0; k< N_LOOP; k++){ 
         ipp = baseHit->phiIndex + loop_phi[k];
//
//--   Gymnastics if phi is closed
//
         if ( ipp < 1 ) {
            if ( getPara()->phiClosed )
               ipp = getPara()->nPhi + ipp ;
            else
               continue ;
         }
         else if ( ipp > getPara()->nPhi ) {
            if ( getPara()->phiClosed )
               ipp = ipp - getPara()->nPhi ;
            else
               continue ;
         }
//
//     Now get eta index
//
         itp = baseHit->etaIndex + loop_eta[k];
         if ( itp <     1      ) continue ;
         if ( itp > getPara()->nEta ) continue ;
//

//
//       Now loop over hits in each volume 
//
         areaIndex = irp   * getPara()->nPhiEtaPlusOne + ipp * getPara()->nEtaPlusOne + itp ;
         for ( clib_FtfHit *candidateHit = (clib_FtfHit *)volume[areaIndex].first ; 
             candidateHit != 0 ;
             candidateHit = (clib_FtfHit *)candidateHit->nextVolumeHit ){

//----------------------------------------------------------------------------
//         Check whether the hit was used before
//--------------------------------------------------------------------------*/
             if ( candidateHit->track != 0 ) continue ;
//--------------------------------------------------------------------------
//         If first points, just choose the closest hit
//-------------------------------------------------------------------------- */
             if ( which_function == USE_SEGMENT ) 
	        result = segmentHitSelection ( baseHit, candidateHit ) ;
	     else 
                result = followHitSelection  ( baseHit, candidateHit, way ) ;

		

//
//     Check result
//
             if ( result > 0 ) {
	        selected_hit = candidateHit ;
                if ( result ==2  ) goto found ; 
             }
//
//       End hit loop  
//
         }
//
//     End row loop      
//
      }
//
//   End volume loop inside cone      
//
   }
found: ;

   return selected_hit ;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   Forms segments
//   Arguments:
//             volume     :    volume pointer
//             way        :    whether to go to negative or positive ir
//             row_to_stop:    row index where to stop
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*__device__*/ int clib_FtfTrack::segment( clib_FtfContainer *volume, int way ){
//
//   Define some variables
//
	float dx, dy, rr ;
	clib_FtfHit* nextHit ;
//
//   Check which way to go
//
	if ( way < 0 ) 
	   nextHit = (clib_FtfHit *)lastHit ;
	else
	   nextHit = (clib_FtfHit *)firstHit ;

//
//    Loop as long a a hit is found and the segment
//    is shorter than n_hit_segm
//
   while ( nextHit != 0 && nHits < getPara()->nHitsForSegment ) {
      chi2[0] = getPara()->maxDistanceSegment ; ;
      nextHit = seekNextHit ( volume, nextHit, way*getPara()->segmentRowSearchRange, 
                              USE_SEGMENT ) ;

//
//     If sz fit update s
//
      if ( nextHit != 0 ){
//
//   Calculate track length if sz plane considered
//
         if ( getPara()->szFitFlag  ){
            dx = ((clib_FtfBaseHit *)nextHit)->x - ((clib_FtfBaseHit *)lastHit)->x ;
            dy = ((clib_FtfBaseHit *)nextHit)->y - ((clib_FtfBaseHit *)lastHit)->y ;
            length    += (float)sqrt ( dx * dx + dy * dy ) ;
            nextHit->s      = length ;
         }
//
//   Calculate conformal coordinates
//
         if ( getPara()->primaries == 0 ){
            rr = square ( xRefHit - nextHit->x ) +
                 square ( yRefHit - nextHit->y ) ;


            nextHit->xp    =   ( nextHit->x - xRefHit ) / rr ;
            nextHit->yp    = - ( nextHit->y - yRefHit ) / rr ;
            nextHit->wxy   = rr * rr / ( square(getPara()->xyErrorScale)  *
                                         square(nextHit->dx) + square(nextHit->dy) ) ;
         }
//
//     Add hit to track
//
	 add ( nextHit, way ) ;
      }
   } // End while ( lastHit ...
//
//    If number of hits is as expected return 1 
//
   if ( nHits == getPara()->nHitsForSegment )
      return 1 ;
   else
      return 0 ;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     Routine to look for segments.
//	 Arguments:
//	 baseHit:       Hit from which track is being extrapolated
//   candidateHit:  Hit being examined as a candidate to which extend track
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*__device__*/ int clib_FtfTrack::segmentHitSelection ( clib_FtfHit *baseHit, clib_FtfHit *candidateHit ){


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
   if ( getPara()->nHitsForSegment > 2 && nHits-1 < getPara()->nHitsForSegment ) {
      dx = candidateHit->x - baseHit->x ;
      dy = candidateHit->y - baseHit->y ;
      angle = (float)atan2 ( dy, dx ) ;
      if ( angle < 0  ) angle = angle + twoPi ;
      lastXyAngle = angle ;
   }

   if ( d3 < chi2[0] ) {
//
//   For second hit onwards check the difference in angle 
//   between the last two track segments
//
    
      if ( nHits > 1 ) {
	 dx     = candidateHit->x - baseHit->x ;
         dy     = candidateHit->y - baseHit->y ;
         angle  = (float)atan2 ( dy, dx ) ;
         if ( angle < 0  ) angle = angle + twoPi ;
	    dangle = (float)fabs ( lastXyAngle - angle );
		  
	    lastXyAngle = angle ;
         if ( dangle > getPara()->segmentMaxAngle ) return 0 ;
      }
//
//    Check whether this is the "closest" hit
//
      chi2[0]          = d3 ;
      if ( d3 < getPara()->goodDistance ) return 2 ;
	  return 1 ;
   }
//
//    If hit does not fulfill criterai return 0
//
   return 0 ;
}

