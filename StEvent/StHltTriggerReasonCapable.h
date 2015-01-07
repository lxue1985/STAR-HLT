/*!
 * \class StHltTriggerReasonCapable 
 * \author Liang Xue & Aihong Tang, Nov 2010
 */


/// a virtual class to make StHltHeavyFragment StHltHighPt StHltDiElectron capable
#ifndef StHltTriggerReasonCapable_hh
#define StHltTriggerReasonCapable_hh

#include <Stiostream.h>
#include "StObject.h"
#include "StArray.h"

class StHltTriggerReasonCapable : public StObject {
public:
	StHltTriggerReasonCapable();
	~StHltTriggerReasonCapable();

	ClassDef(StHltTriggerReasonCapable,1)

};


ostream& operator<<(ostream&, const StHltTriggerReasonCapable&); //print operator

#endif
