/*!
 * \class StHltBEmcTowerHit 
 * \author  Liang Xue & Aihong Tang,  Nov 2010
 */

#include "StHltBEmcTowerHit.h"
#include "StHltTrackNode.h"

ClassImp(StHltBEmcTowerHit)

StHltBEmcTowerHit::StHltBEmcTowerHit()
{
	mAdc = 0;
	mEnergy = 0;
	mPhi = 0; 
	mEta = 0; 
	mZ = 0;
 	mSoftId = 0;
	mDaqId = 0;
	mTrackNode = 0; 
}

StHltBEmcTowerHit::~StHltBEmcTowerHit(){/* noop */}

void 
StHltBEmcTowerHit::setAdc(int val)
{
	mAdc = val;
}

void 
StHltBEmcTowerHit::setEnergy(float val)
{
	mEnergy = val;
}

void
StHltBEmcTowerHit::setPhi(float val)
{
	mPhi = val;
}

void
StHltBEmcTowerHit::setEta(float val)
{
	mEta = val;
}

void
StHltBEmcTowerHit::setZ(float val)
{
	mZ = val;
}

void
StHltBEmcTowerHit::setSoftId(int val)
{
	mSoftId = val;
}

void
StHltBEmcTowerHit::setDaqId(int val)
{
	mDaqId = val;
}

void
StHltBEmcTowerHit::setTrackNode(StHltTrackNode* val)
{
	    mTrackNode = val;
}

ostream&
operator<<(ostream &os, const StHltBEmcTowerHit& hit)
{
	os << " adc "<<hit.adc()<<" energy "<<hit.energy()<<" phi "<<hit.phi()<<" eta "<<hit.eta()<<" z "<<hit.z()<<" softId "<<hit.softId()<<" daqId "<<hit.daqId()<<endl;
	return os;
}
