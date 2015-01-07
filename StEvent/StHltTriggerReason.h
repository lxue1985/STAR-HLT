/*!
 * \class StHltTriggerReason 
 * \author Liang Xue & Aihong Tang, Nov 2010
 */



//need to add enumeration for reason names.
#ifndef StHltTriggerReason_hh
#define StHltTriggerReason_hh

#include <Stiostream.h>
#include "StObject.h"
#include "StArray.h"

class StHltTriggerReasonCapable;

class StHltTriggerReason : public StObject {
	public:
		StHltTriggerReason();
		~StHltTriggerReason();

		int reasonBit() const;
		StHltTriggerReasonCapable* reason();
		const StHltTriggerReasonCapable* reason() const;

		void setReasonBit(int);
		void setReason(StHltTriggerReasonCapable*);

	private:
		int mReasonBit;  
		/* 
		   high pt 0x10000
		   diElectron 0x20000
		   heavyFragment 0x40000
		   allEvents 0x80000
		   randomEvents 0x100000
		   BESGoodEvents 0x200000
		 */
#ifdef __CINT__
		StObjLink mReason;
#else
		StLink<StHltTriggerReasonCapable> mReason;
#endif //__CINT__

		// can be an pointer to StHltHighPt StHltHeavyFragment or StHltDielectron

		ClassDef(StHltTriggerReason,1)

};

inline int StHltTriggerReason::reasonBit() const {return mReasonBit;}
inline StHltTriggerReasonCapable* StHltTriggerReason::reason() { return mReason; }
inline const StHltTriggerReasonCapable* StHltTriggerReason::reason() const { return mReason; }

ostream& operator<<(ostream&, const StHltTriggerReason&); //print operator

#endif
