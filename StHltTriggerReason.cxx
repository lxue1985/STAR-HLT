/*!
 * \class StHltTriggerReason 
 * \author Liang Xue & Aihong Tang, Nov 2010
 */

#include "StHltTriggerReason.h"
#include "StHltTriggerReasonCapable.h"

ClassImp(StHltTriggerReason)

StHltTriggerReason::StHltTriggerReason()
{
	mReasonBit = 0;
	mReason = 0;
}

StHltTriggerReason::~StHltTriggerReason(){ /* noop */}

void StHltTriggerReason::setReasonBit(int val) { mReasonBit = val; }

void StHltTriggerReason::setReason(StHltTriggerReasonCapable* val) { mReason = val; }

ostream&
operator<<(ostream &os, const StHltTriggerReason& r)
{
	os << " reasonBit "<<r.reasonBit()<<endl;
	return os;
}
