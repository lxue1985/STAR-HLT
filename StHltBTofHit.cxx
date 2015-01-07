/*!
 * \class StHltBTofHit 
 * \author Liang Xue & Aihong Tang, Nov 2010
 */

#include "StHltBTofHit.h"
#include "StHltTrackNode.h"

ClassImp(StHltBTofHit)

StHltBTofHit::StHltBTofHit()
{
        mTrayId = 0;
        mChannel = 0; // = nModule*6 + nCell
        mTdc = 0;
        mTot = 0;
        mTof = 0;
        mTriggerTime = 0;
	mTrackNode = 0;
}

StHltBTofHit::~StHltBTofHit(){/* noop */}

void
StHltBTofHit::setTrayId(short val)
{
	mTrayId = val;
}

void
StHltBTofHit::setChannel(short val)
{
	mChannel = val;
}

void
StHltBTofHit::setTdc(float val)
{
	mTdc = val;
}

void
StHltBTofHit::setTot(float val)
{
	mTot = val;
}

void
StHltBTofHit::setTof(float val)
{
	mTof = val;
}

void
StHltBTofHit::setTriggerTime(float val)
{
	mTriggerTime = val;
}

void
StHltBTofHit::setTrackNode(StHltTrackNode* val)
{
	mTrackNode = val;
}


ostream&
operator<<(ostream &os, const StHltBTofHit& hit)
{
	os << " trayId "<<hit.trayId()<<" channel "<<hit.channel()
	   << " module "<<hit.module()<<" cell "<<hit.cell()<<endl
	   << " tdc "<<hit.tdc()<<" tot "<<hit.tot()<<" tof "<<hit.tof()<<" triggerTime "<<hit.triggerTime()<<endl;
	return os;
}

