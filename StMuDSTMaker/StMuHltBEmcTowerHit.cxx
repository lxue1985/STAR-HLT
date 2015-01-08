/*!
 * \class StMuHltBEmcTowerHit 
 * \author Liang Xue , Nov 2010
 */

/** @class StMuHltBEmcTowerHit
 * Class holding the hlt online BEmc properties
 * All unitis are in standard STAR units: GeV,cm
 **/


#include "StMuHltBEmcTowerHit.h"
#include "StHltBEmcTowerHit.h"
#include "StHltTrackNode.h"

 
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
StMuHltBEmcTowerHit::StMuHltBEmcTowerHit(const StHltBEmcTowerHit* hltBEmcTowerHit) :
	mAdc(-999),mSoftId(-999),mDaqId(-999),mEnergy(-999.),mPhi(-999.),mEta(-999.),mZ(-999.),
	mIndex2Global(-999),mIndex2Primary(-999),mIndex2BEmcTowerHit(-999),mIndex2BTofHit(-999)
{

	if(hltBEmcTowerHit){

		StHltTrackNode* trackNode = (StHltTrackNode*) hltBEmcTowerHit->trackNode();
		if(trackNode){
			mIndex2Global     =   trackNode->globalTrackSN();
			mIndex2Primary    =   trackNode->primaryTrackSN();
			mIndex2BEmcTowerHit = trackNode->emcTowerSN();
			mIndex2BTofHit = trackNode->tofHitSN();
		}

		mAdc      =   hltBEmcTowerHit->adc();
		mSoftId   =   hltBEmcTowerHit->softId();
		mDaqId    =   hltBEmcTowerHit->daqId();
		mEnergy   =   hltBEmcTowerHit->energy();
		mPhi      =   hltBEmcTowerHit->phi();
		mEta      =   hltBEmcTowerHit->eta();
		mZ        =   hltBEmcTowerHit->z();

	}

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void StMuHltBEmcTowerHit::Print(Option_t *option) const {
  //
  // Print out some of the essential Emc info. 
  //
  cout << " Adc : " << mAdc << " soft ID : " << mSoftId << endl;
  cout << " daq ID : " << mDaqId << " Energy : " << mEnergy << endl;
  cout << " phi : " << mPhi << " eta : " << mEta <<endl;
  cout << " z : " << mZ <<endl;

}


ClassImp(StMuHltBEmcTowerHit)
