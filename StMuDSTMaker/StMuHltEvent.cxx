/*!
 * \class StMuHltEvent 
 * \author Liang Xue , Nov 2010
 */

#include "StEvent/StEvent.h"
#include "StEvent/StHltEvent.h"

#include "StEvent/StHltTrack.h" 
#include "StEvent/StHltBTofHit.h"
#include "StEvent/StHltVpdHit.h"
#include "StEvent/StHltTriggerReason.h"
#include "StEvent/StHltHighPt.h"
#include "StEvent/StHltHeavyFragment.h"
#include "StEvent/StHltDiElectron.h"
#include "StEvent/StHltBEmcTowerHit.h"

#include "StMuHltTrack.h"
#include "StMuHltBTofHit.h"
#include "StMuHltBEmcTowerHit.h"
#include "StMuHltVpdHit.h"
#include "StMuHltHighPt.h"
#include "StMuHltHeavyFragment.h"
#include "StMuHltDiElectron.h"
#include "StMuHltEvent.h"

ClassImp(StMuHltEvent)

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
StMuHltEvent::StMuHltEvent() : mVersion(-999),mTriggerReasonBitOred(-999),mVpdVertexZ(-999.),mT0(-999.),
mInnerSecGain(-999.),mOuterSecGain(-999.),mVertex(-999.,-999.,-999.),mLowMultVertex(-999.,-999.,-999.)
{ 
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
StMuHltEvent::StMuHltEvent(const StHltEvent &hlteventOb) : mVersion(-999),mTriggerReasonBitOred(-999),
mVpdVertexZ(-999.),mT0(-999.),mInnerSecGain(-999.),mOuterSecGain(-999.),mVertex(-999.,-999.,-999.),
mLowMultVertex(-999.,-999.,-999.)
{ 
	fillEvent(hlteventOb);  /// fill hlt event info
	cout<<"In StMuHltEvent Contructor"<<endl;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
StMuHltEvent::~StMuHltEvent()
{
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void StMuHltEvent::clear()
{
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void StMuHltEvent::fillEvent(const StHltEvent &hlteventOb)
{

	const StHltEvent *hltevent = &hlteventOb;
	
	/// classes that we just copy from StHltEvent
	if(hltevent){

		cout<<"Before filling event stuffs"<<endl;
		mVersion = hltevent->version();
		mTriggerReasonBitOred = hltevent->triggerReasonBitOred();
		mVpdVertexZ =  hltevent->vpdVertexZ();
		mT0 =  hltevent->t0();
		mInnerSecGain =  hltevent->innerSecGain();
		mOuterSecGain =  hltevent->outerSecGain();
		mVertex = hltevent->vertex();
		mLowMultVertex = hltevent->lowMultVertex();
	
		cout<<"Before filling tracks etc"<<endl;
		
		fillPrimaryTracks(hltevent);  // fill global tracks
		fillGlobalTracks(hltevent);   // fill primary tracks
		fillBTofHit(hltevent);        // fill btof info
		fillBEmcTowerHit(hltevent);   // fill bemc info
		fillVpdHit(hltevent);         // fill vpd info
		
		cout<<"Before filling HighPt HeavyFragment Di-electron pairs"<<endl;
		fillTriggerReason(hltevent);  // fill highpt heavyfragment di-electron
		
	}
	else cout<< "no StHltEvent in this event !" << endl;

} 

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void StMuHltEvent::fillPrimaryTracks(const StHltEvent *hltevent){
	
	const StSPtrVecHltTrack& VecpTrack = hltevent->primaryTrack();
	
	for(u_int i=0; i< VecpTrack.size(); i++){

		StHltTrack *pTrack = (StHltTrack*)VecpTrack.at(i);
		mHltPrimaryTracks.push_back(StMuHltTrack(pTrack));
		
	}	
	
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void StMuHltEvent::fillGlobalTracks(const StHltEvent* hltevent){

	const StSPtrVecHltTrack& VecgTrack = hltevent->globalTrack();

	for(u_int i=0; i< VecgTrack.size(); i++){

		StHltTrack *gTrack = (StHltTrack*)VecgTrack.at(i);
		mHltGlobalTracks.push_back(StMuHltTrack(gTrack));

	}

}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void StMuHltEvent::fillBTofHit(const StHltEvent* hltevent){

	const StSPtrVecHltBTofHit& VecbTofHit = hltevent->bTofHit();

	for(u_int i=0; i< VecbTofHit.size(); i++){

		StHltBTofHit *bTofHit = (StHltBTofHit*)VecbTofHit.at(i);
		mHltBTofHit.push_back(StMuHltBTofHit(bTofHit));

	}

}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void StMuHltEvent::fillBEmcTowerHit(const StHltEvent* hltevent){

	const StSPtrVecHltBEmcTowerHit& VecbEmcHit = hltevent->bEmcTowerHits();

	for(u_int i=0; i< VecbEmcHit.size(); i++){

		StHltBEmcTowerHit *bEmcHit = (StHltBEmcTowerHit*)VecbEmcHit.at(i);
		mHltBEmcTowerHit.push_back(StMuHltBEmcTowerHit(bEmcHit));

	}

}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void StMuHltEvent::fillVpdHit(const StHltEvent* hltevent){

	const StSPtrVecHltVpdHit& VecpVpdHit = hltevent->vpdHit();

	for(u_int i=0; i< VecpVpdHit.size(); i++){

		StHltVpdHit *pVpdHit = (StHltVpdHit*)VecpVpdHit.at(i);
		mHltVpdHit.push_back(StMuHltVpdHit(pVpdHit));

	}

}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void StMuHltEvent::fillTriggerReason(const StHltEvent* hltevent){
	
	const StSPtrVecHltTriggerReason& VecTriggerReason = hltevent->triggerReason();

	for(u_int i=0; i< VecTriggerReason.size(); i++){

		StHltTriggerReason *triggerReason = (StHltTriggerReason*)VecTriggerReason.at(i);
		int reasonBit = triggerReason->reasonBit();
		if(reasonBit==highPtBit){
		   StHltHighPt* highpt = (StHltHighPt*)triggerReason->reason();
		   mHltHighPt.push_back(StMuHltHighPt(highpt));
		}
		else if(reasonBit==heavyFragmentBit){
			StHltHeavyFragment* heavyfragment = (StHltHeavyFragment*)triggerReason->reason();
			mHltHeavyFragment.push_back(StMuHltHeavyFragment(heavyfragment));
		}
		else if(reasonBit==diElectronBit){
            StHltDiElectron* Diep = (StHltDiElectron*)triggerReason->reason();
			mHltDielectron.push_back(StMuHltDiElectron(Diep));
		}
		else continue;

	}

}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
StMuHltTrackPrimaryVector StMuHltEvent::primaryTracks() {return mHltPrimaryTracks;}

StMuHltTrackGlobalVector StMuHltEvent::globalTracks() {return mHltGlobalTracks;}

StMuHltBTofHitVector StMuHltEvent::bTofHit() {return mHltBTofHit;}

StMuHltBEmcTowerHitVector StMuHltEvent::bEmcTowerHit() {return mHltBEmcTowerHit;}

StMuHltVpdHitVector StMuHltEvent::vpdHit() {return mHltVpdHit;}

StMuHltHighPtVector StMuHltEvent::highPt() {return mHltHighPt;}

StMuHltHeavyFragmentVector StMuHltEvent::heavyFragment() {return mHltHeavyFragment;}

StMuHltDiElectronVector StMuHltEvent::diElectron() {return mHltDielectron;}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
StMuHltTrack* StMuHltEvent::primaryTracks(int i) {return &mHltPrimaryTracks[i];}

StMuHltTrack* StMuHltEvent::globalTracks(int i) {return &mHltGlobalTracks[i];}

StMuHltBTofHit* StMuHltEvent::bTofHit(int i) {return &mHltBTofHit[i];}

StMuHltBEmcTowerHit* StMuHltEvent::bEmcTowerHit(int i) {return &mHltBEmcTowerHit[i];}

StMuHltVpdHit* StMuHltEvent::vpdHit(int i) {return &mHltVpdHit[i];}

StMuHltHighPt* StMuHltEvent::highPt(int i) {return &mHltHighPt[i];}

StMuHltHeavyFragment* StMuHltEvent::heavyFragment(int i) {return &mHltHeavyFragment[i];}

StMuHltDiElectron* StMuHltEvent::diElectron(int i) {return &mHltDielectron[i];}





