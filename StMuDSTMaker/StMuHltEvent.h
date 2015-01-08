/*!
 * \class StMuHltEvent 
 * \author Liang Xue , Nov 2010
 */


/**
  @class StMuHltEvent
  The StMuHltEvent class holds the hlt online event-wise information of the STAR's common muDst.  
  Most of its data members are classes from the StHltEvent package.
 **/

#ifndef StMuHltEvent_h
#define StMuHltEvent_h

#define highPtBit 0x10000
#define diElectronBit 0x20000
#define heavyFragmentBit 0x40000
#define allEvent 0x80000
#define randomEvents 0x100000
#define BESGoodEvents 0x200000

#include <vector>
#include "TObject.h"
#include "StThreeVectorF.hh"


class StHltEvent;

class StMuHltTrack;
class StMuHltBTofHit;
class StMuHltBEmcTowerHit;
class StMuHltVpdHit;
class StMuHltHighPt;
class StMuHltHeavyFragment;
class StMuHltDiElectron;


typedef vector<StMuHltTrack> StMuHltTrackGlobalVector;
typedef vector<StMuHltTrack> StMuHltTrackPrimaryVector;
typedef vector<StMuHltBTofHit> StMuHltBTofHitVector;
typedef vector<StMuHltBEmcTowerHit> StMuHltBEmcTowerHitVector;
typedef vector<StMuHltVpdHit> StMuHltVpdHitVector;
typedef vector<StMuHltHighPt> StMuHltHighPtVector;
typedef vector<StMuHltHeavyFragment> StMuHltHeavyFragmentVector;
typedef vector<StMuHltDiElectron> StMuHltDiElectronVector;


class StMuHltEvent : public TObject {
	public:
		StMuHltEvent();
		StMuHltEvent(const StHltEvent &hlteventOb); 
		~StMuHltEvent();

		int  version();
		int  triggerReasonBitOred();
		float vpdVertexZ();
		float t0();
		float innerSecGain();
		float outerSecGain();
		StThreeVectorF& vertex();
		StThreeVectorF& lowMultVertex();
        
		StMuHltTrackPrimaryVector primaryTracks();   /// return vector of primary tracks
		StMuHltTrackGlobalVector globalTracks();     /// return vector of global tracks
		StMuHltBTofHitVector bTofHit();              /// return vector of btof hits
		StMuHltBEmcTowerHitVector bEmcTowerHit();    /// return vector of bemc hits
		StMuHltVpdHitVector vpdHit();                /// return vector of vpd hits
		StMuHltHighPtVector highPt();                /// return vector of highpt
		StMuHltHeavyFragmentVector heavyFragment();  /// return vector of heavy fragment
		StMuHltDiElectronVector diElectron();        /// return vector of di-electron

	    StMuHltTrack* primaryTracks(int i);          /// return pointer of particular primary track
		StMuHltTrack* globalTracks(int i);           /// return pointer of particular global track
		StMuHltBTofHit* bTofHit(int i);              /// return pointer of particular btof hit
		StMuHltBEmcTowerHit* bEmcTowerHit(int i);    /// return pointer of particular bemc hit
		StMuHltVpdHit* vpdHit(int i);                /// return pointer of particular vpd hit
		StMuHltHighPt* highPt(int i);                /// return pointer of particular highpt 
        StMuHltHeavyFragment* heavyFragment(int i);  /// return pointer of particular heavy fragment
		StMuHltDiElectron* diElectron(int i);        /// return pointer of particular di-electron

	protected:

		Int_t mVersion;
		Int_t mTriggerReasonBitOred;    /// reason with "OR" operator
		Float_t mVpdVertexZ;
		Float_t mT0;
		Float_t mInnerSecGain; //dedx gain
		Float_t mOuterSecGain;
		StThreeVectorF mVertex;
		StThreeVectorF mLowMultVertex;

		void clear();
		void fillEvent(const StHltEvent &hlteventOb);
		void fillPrimaryTracks(const StHltEvent *hltevent);
		void fillGlobalTracks(const StHltEvent *hltevent);
		void fillBTofHit(const StHltEvent *hltevent);
		void fillBEmcTowerHit(const StHltEvent *hltevent);
		void fillVpdHit(const StHltEvent *hltevent);
		void fillTriggerReason(const StHltEvent *hltevent);
        
     public:

		StMuHltTrackGlobalVector mHltGlobalTracks;     /// vector of global tracks
		StMuHltTrackPrimaryVector mHltPrimaryTracks;   /// vector of primary tracks
		StMuHltBTofHitVector mHltBTofHit;              /// vector of btof hits
		StMuHltBEmcTowerHitVector mHltBEmcTowerHit;    /// vector of bemc hits
		StMuHltVpdHitVector mHltVpdHit;                /// vector of vpd hits
		StMuHltHighPtVector mHltHighPt;                /// vector of highpt
		StMuHltHeavyFragmentVector mHltHeavyFragment;  /// vector of heavy fragment
		StMuHltDiElectronVector mHltDielectron;        /// vectr of di-electron
		

		ClassDef(StMuHltEvent,1)
};

inline int StMuHltEvent::version() { return mVersion;}
inline int StMuHltEvent::triggerReasonBitOred() { return mTriggerReasonBitOred;}
inline float StMuHltEvent::vpdVertexZ() { return mVpdVertexZ;}
inline float StMuHltEvent::t0() {return mT0;}
inline float StMuHltEvent::innerSecGain() { return mInnerSecGain; }
inline float StMuHltEvent::outerSecGain() { return mOuterSecGain; }
inline StThreeVectorF& StMuHltEvent::vertex() {return mVertex;}
inline StThreeVectorF& StMuHltEvent::lowMultVertex() {return mLowMultVertex;}


#endif
