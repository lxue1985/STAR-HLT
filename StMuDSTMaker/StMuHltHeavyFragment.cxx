/*!
 * \class StMuHltHeavyFragment 
 * \author Liang Xue , Nov 2010
 */

#include "StMuHltHeavyFragment.h"


#include "StHltHeavyFragment.h"
#include "StHltTrack.h"

 
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
StMuHltHeavyFragment::StMuHltHeavyFragment(const StHltHeavyFragment* heavyfragment) :
	mBEmcMatchPhiDiff(-999.),mBEmcMatchZEdge(-999.),mBTofProjChannel(-999.),mBTofCellLocalY(-999.),
	mBTofCellLocalZ(-999.),mBTofPathLength(-999.),mBeta(-999.),mTof(-999.),mIndex2Global(-999),
	mIndex2Primary(-999),mIndex2BTofHit(-999),mIndex2BEmcTowerHit(-999)
{

	if(heavyfragment){
		
		StHltTrack gTrack = heavyfragment->globalTrack();
        mGlobalTrack = StMuHltTrack(&gTrack);

		StHltTrack pTrack = heavyfragment->primaryTrack();
        mPrimaryTrack = StMuHltTrack(&pTrack);

		StHltBTofHit bTofHit = heavyfragment->bTofHit();
		mBTofHit = StMuHltBTofHit(&bTofHit);

        StHltBEmcTowerHit bEmcTowerHit = heavyfragment->bEmcTowerHit();
		mBEmcTowerHit = StMuHltBEmcTowerHit(&bEmcTowerHit);


		mIndex2Global     =   heavyfragment->globalTrackSN();
		mIndex2Primary    =   heavyfragment->primaryTrackSN();
		mIndex2BTofHit       =   heavyfragment->tofHitSN();
		mIndex2BEmcTowerHit       =   heavyfragment->emcTowerSN();
		mBEmcMatchPhiDiff =   heavyfragment->bEmcMatchPhiDiff();
		mBEmcMatchZEdge   =   heavyfragment->bEmcMatchZEdge();
		mBTofProjChannel  =   heavyfragment->bTofProjChannel();
		mBTofCellLocalY   =   heavyfragment->bTofCellLocalY();
		mBTofCellLocalZ   =   heavyfragment->bTofCellLocalZ();
		mBTofPathLength   =   heavyfragment->bTofPathLength();
		mBeta             =   heavyfragment->beta();
		mTof              =   heavyfragment->tof();

	}

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void StMuHltHeavyFragment::Print(Option_t *option) const {
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


ClassImp(StMuHltHeavyFragment)
