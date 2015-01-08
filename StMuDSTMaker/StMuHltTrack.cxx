/*!
 * \class StMuHltTrack 
 * \author Liang Xue , Nov 2010
 */

#include "StMuHltTrack.h"


#include "StHltTrackNode.h"
#include "StHltTrack.h"

 
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
StMuHltTrack::StMuHltTrack(const StHltTrack* hltTrack) : mId(-999), mFlag(-999),mInnerMostRow(-999),
  mOuterMostRow(-999),mReserved(-999),mNHits(-999), mNHitsDedx(-999),mQ(0),mChiSqXY(-999.),
  mChiSqZ(-999.),mdEdx(-999.),mPt(-999.),mPhi0(-999.),mPsi(-999.),mR0(-999.),mTanl(-999.),mZ0(-999.),
  mLength(-999.),mDpt(-999.),mDpsi(-999.),mDz0(-999.),mDtanl(-999.),mP(0.,0.,0.),
  mIndex2Global(-999),mIndex2Primary(-999),mIndex2BEmcTowerHit(-999),mIndex2BTofHit(-999)
{
	if(hltTrack){

		StHltTrackNode* trackNode = (StHltTrackNode*)hltTrack->trackNode();
		if(trackNode){
			mIndex2Global = trackNode->globalTrackSN();
			mIndex2Primary = trackNode->primaryTrackSN();
			mIndex2BEmcTowerHit = trackNode->emcTowerSN();
			mIndex2BTofHit = trackNode->tofHitSN();
		}

		mP = momentumAtPrimaryVertex(hltTrack); 

		mId = hltTrack->id();
		mFlag = hltTrack->flag();
		mInnerMostRow = hltTrack->innerMostRow();
		mOuterMostRow = hltTrack->outerMostRow();
		mReserved = hltTrack->reserved();
		mNHits = hltTrack->nHits();
		mNHitsDedx = hltTrack->ndedx();
		mQ = hltTrack->q();
		mChiSqXY = hltTrack->chi2(0);
		mChiSqZ = hltTrack->chi2(1);
		mdEdx = hltTrack->dedx();
		mPt = hltTrack->pt();
		mPhi0 = hltTrack->phi0();
		mPsi = hltTrack->psi();
		mR0 = hltTrack->r0();
		mTanl = hltTrack->tanl();
		mZ0 = hltTrack->z0();
		mLength = hltTrack->length();
		mDpt = hltTrack->dpt();
		mDpsi = hltTrack->dpsi();
		mDz0 = hltTrack->dz0();
		mDtanl = hltTrack->dtanl();

	}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void StMuHltTrack::Print(Option_t *option) const {
  //
  // Print out some of the essential HLT Track info. 
  // Note: more is stored in the HLT Track; may expand printout
  // or specify options some time
  //
  cout << "HLT Track, id " << mId << ", flag " << mFlag << " (>0 is OK)" << endl;
  cout << "momentum " << mP << endl;
  cout << "phi " << mPhi0 << ", pt " << mPt << endl;
  cout << "Total hits: " << mNHits << ", nHitsDedx " << mNHitsDedx <<endl;

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
StThreeVectorD StMuHltTrack::momentumAtPrimaryVertex(const StHltTrack* track) const
{
	float px = cos(track->psi())*track->pt();
	float py = sin(track->psi())*track->pt();
	float pz = track->tanl()*track->pt();
	StThreeVectorD momentum(px,py,pz);
	return  momentum;
}



ClassImp(StMuHltTrack)
