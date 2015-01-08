/*!
 * \class StMuHltBTofHit 
 * \author Liang Xue , Nov 2010
 */

/** @class StMuHltBTofHit
 * Class holding the hlt btof hit properties
 * All unitis are in standard STAR units: GeV,cm
 **/
#include "StMuHltBTofHit.h"

#include "StHltBTofHit.h"
#include "StHltTrackNode.h"

 
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
StMuHltBTofHit::StMuHltBTofHit(const StHltBTofHit* hltBTofHit) :
   mTrayId(-999),mChannel(-999),mTdc(-999.),mTot(-999.),mTof(-999.),mTriggerTime(-999.),
   mIndex2Global(-999),mIndex2Primary(-999),mIndex2BEmcTowerHit(-999),mIndex2BTofHit(-999)
{

	if(hltBTofHit){	

		StHltTrackNode* trackNode = (StHltTrackNode*) hltBTofHit->trackNode();
		if(trackNode){
			mIndex2Global     =   trackNode->globalTrackSN();
			mIndex2Primary    =   trackNode->primaryTrackSN();
			mIndex2BEmcTowerHit = trackNode->emcTowerSN();
			mIndex2BTofHit = trackNode->tofHitSN();
		}

		mTrayId       =   hltBTofHit->trayId();
		mChannel      =   hltBTofHit->channel();
		mTdc          =   hltBTofHit->tdc();
		mTot          =   hltBTofHit->tot();
		mTof          =   hltBTofHit->tof();
		mTriggerTime  =   hltBTofHit->triggerTime();

	}

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void StMuHltBTofHit::Print(Option_t *option) const {
  //
  // Print out some of the essential HighPt HLT Track info. 
  //
  cout << " Tray ID : " << mTrayId << " Channel ID : " << mChannel << endl;
  cout << " tdc : " << mTdc << " tot : " << mTot << endl;
  cout << " tof : " << mTof << " trigger Time : " << mTriggerTime <<endl;

}


ClassImp(StMuHltBTofHit)


