/*!
 * \class StMuHltVpdHit
 * \author Liang Xue , Nov 2010
 */

#include "StMuHltVpdHit.h"


#include "StHltVpdHit.h"

 
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
StMuHltVpdHit::StMuHltVpdHit(const StHltVpdHit* hltVpdHit) : mDirection(east),
mChannel(-999),mTdc(-999.),mTot(-999.),mTof(-999.),mTriggerTime(-999.)
{

	if(hltVpdHit){

		mDirection    =   hltVpdHit->direction();
		mChannel      =   hltVpdHit->channel();
		mTdc          =   hltVpdHit->tdc();
		mTot          =   hltVpdHit->tot();
		mTof          =   hltVpdHit->tof();
		mTriggerTime  =   hltVpdHit->triggerTime();

	}

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void StMuHltVpdHit::Print(Option_t *option) const {
  //
  // Print out some of the essential Hlt Vpd info. 
  //
  cout << " direction : " << mDirection << " Channel ID : " << mChannel << endl;
  cout << " tdc : " << mTdc << " tot : " << mTot << endl;
  cout << " tof : " << mTof << " trigger Time : " << mTriggerTime <<endl;

}


ClassImp(StMuHltVpdHit)
