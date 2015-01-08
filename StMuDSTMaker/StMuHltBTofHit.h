/*!
 * \class StMuHltBTofHit 
 * \author Liang Xue , Nov 2010
 */

/** @class StMuHltBTofHit
 * Class holding the hlt btof hit properties
 * All unitis are in standard STAR units: GeV,cm
 **/

#ifndef StMuHltBTofHit_h
#define StMuHltBTofHit_h

#include "TObject.h"


class StHltBTofHit;

class StMuHltBTofHit : public TObject {
	public:
		StMuHltBTofHit() {/*no-op*/};          // Default constructor
		StMuHltBTofHit(const StHltBTofHit*);   // Constructor from StHltBTofHit
		~StMuHltBTofHit() {};

		int    channel() const;       // Returns bTof channel ID  ( channelID = nModule*6 + nCell)
		int    module() const;        // Returns bTof module ID  ( moduleID = nModule/6)
		int    cell() const;          // Returns bTof cell ID  ( cellID = nModule%6)
		int    tray() const;          // Returns bTof tray ID
		float   tdc() const;           // Returns bTof tdc
		float   tot() const;           // Returns bTof tot
		float   tof() const;           // Returns bTof  tof
		float   triggerTime() const;   // Returns bTof trigger time

		int index2Global() const;                        // Returns index of associated global track.
		int index2Primary() const;                       // Returns index of associated primary track.
		int index2BEmcTowerHit() const;                  // Returns index of associated bemc tower hit.
		int index2BTofHit() const;                       // Returns index of associated btof hit.

		virtual void Print(Option_t* option = "") const;  // Print track info

	protected:
	    
	    Int_t   mTrayId;
		Int_t   mChannel;
		Float_t mTdc;
		Float_t mTot;
		Float_t mTof;
	    Float_t mTriggerTime;
		
		Int_t mIndex2Global;
		Int_t mIndex2Primary;
		Int_t mIndex2BEmcTowerHit;
		Int_t mIndex2BTofHit;
		void setIndex2Global(int i) {mIndex2Global=i;}             // Set index of associated global track.
		void setIndex2Primary(int i) {mIndex2Primary=i;}           // Set index of associated primary track.
		void setIndex2BEmcTowerHit(int i) {mIndex2BEmcTowerHit=i;} // Set index of associated bemc tower hit.
		void setIndex2BTofHit(int i) {mIndex2BTofHit=i;}           // Set index of associated btof hit.

		ClassDef(StMuHltBTofHit,1)
};

inline int    StMuHltBTofHit::channel() const {return mChannel;}
inline int    StMuHltBTofHit::module() const {return mChannel/6;}
inline int    StMuHltBTofHit::cell() const {return mChannel%6;}
inline int    StMuHltBTofHit::tray() const {return mTrayId;}
inline float   StMuHltBTofHit::tdc() const {return mTdc;}
inline float   StMuHltBTofHit::tot() const {return mTot;}
inline float   StMuHltBTofHit::tof() const {return mTof;}
inline float   StMuHltBTofHit::triggerTime() const {return mTriggerTime;}

inline int StMuHltBTofHit::index2Global() const {return mIndex2Global;}
inline int StMuHltBTofHit::index2Primary() const {return mIndex2Primary;}
inline int StMuHltBTofHit::index2BEmcTowerHit() const {return mIndex2BEmcTowerHit;}
inline int StMuHltBTofHit::index2BTofHit() const {return mIndex2BTofHit;}

#endif
