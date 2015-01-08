/*!
 * \class StMuHltHighPt 
 * \author Liang Xue , Nov 2010
 */

#include "StMuHltHighPt.h"

#include "StHltHighPt.h"
#include "StHltTrack.h"

 
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
StMuHltHighPt::StMuHltHighPt(const StHltHighPt* highpt) : mBEmcMatchPhiDiff(-999.),mBEmcMatchZEdge(-999.),
	mBTofProjChannel(-999.),mBTofCellLocalY(-999.),mBTofCellLocalZ(-999.),mBTofPathLength(-999.),
	mBeta(-999.),mTof(-999.),mIndex2Global(-999),mIndex2Primary(-999),mIndex2BTofHit(-999),
	mIndex2BEmcTowerHit(-999)
{

	if(highpt){
		
		StHltTrack gTrack = highpt->globalTrack();
        mGlobalTrack = StMuHltTrack(&gTrack);

		StHltTrack pTrack = highpt->primaryTrack();
        mPrimaryTrack = StMuHltTrack(&pTrack);

		StHltBTofHit bTofHit = highpt->bTofHit();
		mBTofHit = StMuHltBTofHit(&bTofHit);

        StHltBEmcTowerHit bEmcTowerHit = highpt->bEmcTowerHit();
		mBEmcTowerHit = StMuHltBEmcTowerHit(&bEmcTowerHit);

		mIndex2Global     =   highpt->globalTrackSN();
		mIndex2Primary    =   highpt->primaryTrackSN();
		mIndex2BTofHit       =   highpt->tofHitSN();
		mIndex2BEmcTowerHit       =   highpt->emcTowerSN();

		mBEmcMatchPhiDiff =   highpt->bEmcMatchPhiDiff();
		mBEmcMatchZEdge   =   highpt->bEmcMatchZEdge();
		mBTofProjChannel  =   highpt->bTofProjChannel();
		mBTofCellLocalY   =   highpt->bTofCellLocalY();
		mBTofCellLocalZ   =   highpt->bTofCellLocalZ();
		mBTofPathLength   =   highpt->bTofPathLength();
		mBeta             =   highpt->beta();
		mTof              =   highpt->tof();

	}

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void StMuHltHighPt::Print(Option_t *option) const {
	//
	// Print out some of the essential HighPt HLT Track info. 
	//
	cout << " Index2Global Track : " << mIndex2Global << " Index2Primary Track : " << mIndex2Primary << endl;
	cout << " Index2BTof Hit : " << mIndex2BTofHit << " Index2BEmcTower Hit : " << mIndex2BEmcTowerHit << endl;
	cout << " BEMC Matched PhiDiff : " << mBEmcMatchPhiDiff << " BEMC Matched ZEdge : " << mBEmcMatchZEdge <<endl;
	cout << " BTOF Projected Channel ID : " << mBTofProjChannel << " BTOF Cell LocalY : " << mBTofCellLocalY << endl;
	cout << " BTOF Cell LocalZ : " << mBTofCellLocalZ << " BTOF PathLength : " <<  mBTofPathLength <<endl;
	cout << " BTOF beta : " << mBeta << " BTOF tof : " <<  mTof <<endl;

}



ClassImp(StMuHltHighPt)
