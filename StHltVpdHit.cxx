/*!
 * \class StHltVpdHit 
 * \author Liang Xue & Aihong Tang, Nov 2010
 */

#include "StHltVpdHit.h"
#include "StHltTrackNode.h"

ClassImp(StHltVpdHit)

StHltVpdHit::StHltVpdHit()
{
        mDirection = east; // should there be a "unknown" entry in StBeamDirection ?
        mChannel = 0; // total 19 channels for each side (east/west)
        mTdc = 0;
        mTot = 0;
        mTof = 0;
        mTriggerTime = 0;
}

StHltVpdHit::~StHltVpdHit(){/* noop */}

void
StHltVpdHit::setDirection(StBeamDirection val)
{
	mDirection = val;
}

void
StHltVpdHit::setChannel(short val)
{
	mChannel = val;
}

void
StHltVpdHit::setTdc(float val)
{
	mTdc = val;
}

void
StHltVpdHit::setTot(float val)
{
	mTot = val;
}

void
StHltVpdHit::setTof(float val)
{
	mTof = val;
}

void
StHltVpdHit::setTriggerTime(float val)
{
	mTriggerTime = val;
}

ostream&
operator<<(ostream &os, const StHltVpdHit& hit)
{
	os << " direction "<<hit.direction()<<" channel "<<hit.channel()
	   << " module "<<hit.module()<<" cell "<<hit.cell()<<endl
	   << " tdc "<<hit.tdc()<<" tot "<<hit.tot()<<" tof "<<hit.tof()<<" triggerTime "<<hit.triggerTime()<<endl;
	return os;
}

