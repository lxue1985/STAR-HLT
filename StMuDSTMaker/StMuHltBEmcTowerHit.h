/*!
 * \class StMuHltBEmcTowerHit 
 * \author Liang Xue , Nov 2010
 */

/** @class StMuHltBEmcTowerHit
 * Class holding the hlt online BEmc properties
 * All unitis are in standard STAR units: GeV,cm
 **/

#ifndef StMuHltBEmcTowerHit_h
#define StMuHltBEmcTowerHit_h

#include "TObject.h"



class StHltBEmcTowerHit;

class StMuHltBEmcTowerHit : public TObject {
	public:
		StMuHltBEmcTowerHit() {/*no-op*/};               // Default constructor
		StMuHltBEmcTowerHit(const StHltBEmcTowerHit*);   // Constructor from StHltBEmcTowerHit
		~StMuHltBEmcTowerHit() {};

		int    adc() const;          // Returns bEmc adc
		int    softId() const;       // Returns bEmc soft ID
		int    daqId() const;        // Returns bEmc daq ID
		float   energy() const;       // Returns bEmc energy
		float   phi() const;          // Returns bEmc phi
		float   eta() const;          // Returns bEmc eta
		float   z() const;            // Returns bEmc z
        
		int index2Global() const;                        // Returns index of associated global track.
		int index2Primary() const;                       // Returns index of associated primary track.
		int index2BEmcTowerHit() const;                  // Returns index of associated bemc tower hit.
		int index2BTofHit() const;                       // Returns index of associated btof hit.

		virtual void Print(Option_t* option = "") const; // Print track info

	protected:
	    
	    Int_t      mAdc;
		Int_t      mSoftId;
		Int_t      mDaqId;
		Float_t    mEnergy;   //with online calibration
		Float_t    mPhi;
		Float_t    mEta;
		Float_t    mZ;
		Int_t mIndex2Global;
		Int_t mIndex2Primary;
		Int_t mIndex2BEmcTowerHit;
		Int_t mIndex2BTofHit;
		void setIndex2Global(int i) {mIndex2Global=i;}             // Set index of associated global track.
		void setIndex2Primary(int i) {mIndex2Primary=i;}           // Set index of associated primary track.
		void setIndex2BEmcTowerHit(int i) {mIndex2BEmcTowerHit=i;} // Set index of associated bemc tower hit.
		void setIndex2BTofHit(int i) {mIndex2BTofHit=i;}           // Set index of associated btof hit.

		ClassDef(StMuHltBEmcTowerHit,1)
};

inline int    StMuHltBEmcTowerHit::adc() const {return mAdc;}
inline int    StMuHltBEmcTowerHit::softId() const {return mSoftId;}
inline int    StMuHltBEmcTowerHit::daqId() const {return mDaqId;}
inline float   StMuHltBEmcTowerHit::energy() const {return mEnergy;}
inline float   StMuHltBEmcTowerHit::phi() const {return mPhi;}
inline float   StMuHltBEmcTowerHit::eta() const {return mEta;}
inline float   StMuHltBEmcTowerHit::z() const {return mZ;}

inline int StMuHltBEmcTowerHit::index2Global() const {return mIndex2Global;}
inline int StMuHltBEmcTowerHit::index2Primary() const {return mIndex2Primary;}
inline int StMuHltBEmcTowerHit::index2BEmcTowerHit() const {return mIndex2BEmcTowerHit;}
inline int StMuHltBEmcTowerHit::index2BTofHit() const {return mIndex2BTofHit;}

#endif
