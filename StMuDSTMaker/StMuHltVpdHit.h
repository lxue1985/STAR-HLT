/*!
 * \class StMuHltVpdHit
 * \author Liang Xue , Nov 2010
 */



/** @class StMuHltVpdHit
 * Class holding the hlt online vpd properties
 * All unitis are in standard STAR units: GeV,cm
 **/

#ifndef StMuHltVpdHit_h
#define StMuHltVpdHit_h

#include "TObject.h"
#include "StEnumerations.h"

class StMuHltTrack;
class StHltVpdHit;

class StMuHltVpdHit : public TObject {
	public:
		StMuHltVpdHit() {/*no-op*/};          // Default constructor
		StMuHltVpdHit(const StHltVpdHit*);   // Constructor from StHltVpdHit
		~StMuHltVpdHit() {}; 

		int    channel() const;                 // Returns Vpd channel ID  ( channelID = nModule*6 + nCell)
		int    module() const;                  // Returns Vpd module ID  ( moduleID = nModule/6)
		int    cell() const;                    // Returns Vpd cell ID  ( cellID = nModule%6)
		StBeamDirection    direction() const;   // Returns direction of VPD
		float   tdc() const;                     // Return  Vpd tdc
		float   tot() const;                     // Returns Vpd tot
		float   tof() const;                     // Returns Vpd  tof
		float   triggerTime() const;             // Returns Vpd trigger time

		virtual void Print(Option_t* option = "") const;  // Print vpd info

	protected:
	    
		StBeamDirection mDirection;
	    Int_t   mChannel;
		Float_t mTdc;
		Float_t mTot;
		Float_t mTof;
	    Float_t mTriggerTime;

		ClassDef(StMuHltVpdHit,1)
};

inline int    StMuHltVpdHit::channel() const {return mChannel;}
inline int    StMuHltVpdHit::module() const {return mChannel/6;}
inline int    StMuHltVpdHit::cell() const {return mChannel%6;}
inline StBeamDirection   StMuHltVpdHit::direction() const {return mDirection;}
inline float   StMuHltVpdHit::tdc() const {return mTdc;}
inline float   StMuHltVpdHit::tot() const {return mTot;}
inline float   StMuHltVpdHit::tof() const {return mTof;}
inline float   StMuHltVpdHit::triggerTime() const {return mTriggerTime;}

#endif
