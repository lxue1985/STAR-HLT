/*!
 * \class StMuHltDiElectron 
 * \author Liang Xue , Nov 2010
 */

#include "StMuHltDiElectron.h"

#include "StHltDiElectron.h"
#include "StHltTrack.h"

 
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
StMuHltDiElectron::StMuHltDiElectron(const StHltDiElectron* hltDiElectron) :
	mDaughter1BEmcMatchPhiDiff(-999.),mDaughter1BEmcMatchZEdge(-999.),mDaughter1BTofProjChannel(-999.),
	mDaughter1BTofCellLocalY(-999.),mDaughter1BTofCellLocalZ(-999.),mDaughter1BTofPathLength(-999.),
	mDaughter1Beta(-999.),mDaughter1Tof(-999.),mDaughter1Index2Global(-999),mDaughter1Index2Primary(-999),
	mDaughter1Index2BTof(-999),mDaughter1Index2BEmc(-999),
	mDaughter2BEmcMatchPhiDiff(-999.),mDaughter2BEmcMatchZEdge(-999.),mDaughter2BTofProjChannel(-999.),
	mDaughter2BTofCellLocalY(-999.),mDaughter2BTofCellLocalZ(-999.),mDaughter2BTofPathLength(-999.),
	mDaughter2Beta(-999.),mDaughter2Tof(-999.),mDaughter2Index2Global(-999),mDaughter2Index2Primary(-999),
	mDaughter2Index2BTof(-999),mDaughter2Index2BEmc(-999),
	mPt(-999.),mPsi(-999.),mTanl(-999.),mInvariantMass(-999.)
{

	// daughter 1 information

	StHltTrack daughter1gTrack = hltDiElectron->daughter1globalTrack();
	mDaughter1GlobalTrack = StMuHltTrack(&daughter1gTrack);

	StHltTrack daughter1pTrack = hltDiElectron->daughter1primaryTrack();
	mDaughter1PrimaryTrack = StMuHltTrack(&daughter1pTrack);

	StHltBTofHit daughter1bTofHit = hltDiElectron->daughter1bTofHit();
	mDaughter1BTofHit = StMuHltBTofHit(&daughter1bTofHit);

	StHltBEmcTowerHit daughter1bEmcTowerHit = hltDiElectron->daughter1bEmcTowerHit();
	mDaughter1BEmcTowerHit = StMuHltBEmcTowerHit(&daughter1bEmcTowerHit);

	mDaughter1Index2Global     =   hltDiElectron->daughter1globalTrackSN();
	mDaughter1Index2Primary    =   hltDiElectron->daughter1primaryTrackSN();
	mDaughter1Index2BTof       =   hltDiElectron->daughter1tofHitSN();
	mDaughter1Index2BEmc       =   hltDiElectron->daughter1emcTowerSN();

	mDaughter1BEmcMatchPhiDiff =   hltDiElectron->daughter1bEmcMatchPhiDiff();
	mDaughter1BEmcMatchZEdge   =   hltDiElectron->daughter1bEmcMatchZEdge();
	mDaughter1BTofProjChannel  =   hltDiElectron->daughter1bTofProjChannel();
	mDaughter1BTofCellLocalY   =   hltDiElectron->daughter1bTofCellLocalY();
	mDaughter1BTofCellLocalZ   =   hltDiElectron->daughter1bTofCellLocalZ();
	mDaughter1BTofPathLength   =   hltDiElectron->daughter1bTofPathLength();
	mDaughter1Beta             =   hltDiElectron->daughter1beta();
	mDaughter1Tof              =   hltDiElectron->daughter1tof();


	// daughter 2 information

	StHltTrack daughter2gTrack = hltDiElectron->daughter2globalTrack();
	mDaughter2GlobalTrack = StMuHltTrack(&daughter2gTrack);

	StHltTrack daughter2pTrack = hltDiElectron->daughter2primaryTrack();
	mDaughter2PrimaryTrack = StMuHltTrack(&daughter2pTrack);

	StHltBTofHit daughter2bTofHit = hltDiElectron->daughter2bTofHit();
	mDaughter2BTofHit = StMuHltBTofHit(&daughter2bTofHit);

	StHltBEmcTowerHit daughter2bEmcTowerHit = hltDiElectron->daughter2bEmcTowerHit();
	mDaughter2BEmcTowerHit = StMuHltBEmcTowerHit(&daughter2bEmcTowerHit);

	mDaughter2Index2Global     =   hltDiElectron->daughter2globalTrackSN();
	mDaughter2Index2Primary    =   hltDiElectron->daughter2primaryTrackSN();
	mDaughter2Index2BTof       =   hltDiElectron->daughter2tofHitSN();
	mDaughter2Index2BEmc       =   hltDiElectron->daughter2emcTowerSN();

	mDaughter2BEmcMatchPhiDiff =   hltDiElectron->daughter2bEmcMatchPhiDiff();
	mDaughter2BEmcMatchZEdge   =   hltDiElectron->daughter2bEmcMatchZEdge();
	mDaughter2BTofProjChannel  =   hltDiElectron->daughter2bTofProjChannel();
	mDaughter2BTofPathLength   =   hltDiElectron->daughter2bTofPathLength();
	mDaughter2Beta             =   hltDiElectron->daughter2beta();
	mDaughter2Tof              =   hltDiElectron->daughter2tof();

    /// di-electron pair information
	mPt                        =   hltDiElectron->pt();
	mPsi                       =   hltDiElectron->psi();
	mTanl                      =   hltDiElectron->tanl();
	mInvariantMass             =   hltDiElectron->invariantMass();


}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void StMuHltDiElectron::Print(Option_t *option) const {
  //
  // Print out some of the essential HLT Di-Electron info. 
  //
  cout << " Daughter1 Index2Global Track : " << mDaughter1Index2Global << " Daughter1 Index2Primary Track : " << mDaughter1Index2Primary << endl;
  cout << " Daughter1 Index2BTof Hit : " << mDaughter1Index2BTof << " Daughter1 Index2BEmcTower Hit : " << mDaughter1Index2BEmc << endl;
  cout << " Daughter1 BEMC Matched PhiDiff : " << mDaughter1BEmcMatchPhiDiff << " Daughter1 BEMC Matched ZEdge : " << mDaughter1BEmcMatchZEdge <<endl;
  cout << " Daughter1 BTOF Projected Channel ID : " << mDaughter1BTofProjChannel << " Daughter1 BTOF Cell LocalY : " << mDaughter1BTofCellLocalY << endl;
  cout << " Daughter1 BTOF Cell LocalZ : " << mDaughter1BTofCellLocalZ << " Daughter1 BTOF PathLength : " <<  mDaughter1BTofPathLength <<endl;
  cout << " Daughter1 BTOF beta : " << mDaughter1Beta << " Daughter1 BTOF tof : " <<  mDaughter1Tof <<endl;

  cout << " Daughter2 Index2Global Track : " << mDaughter2Index2Global << " Daughter2 Index2Primary Track : " << mDaughter2Index2Primary << endl;
  cout << " Daughter2 Index2BTof Hit : " << mDaughter2Index2BTof << " Daughter2 Index2BEmcTower Hit : " << mDaughter2Index2BEmc << endl;
  cout << " Daughter2 BEMC Matched PhiDiff : " << mDaughter2BEmcMatchPhiDiff << " Daughter2 BEMC Matched ZEdge : " << mDaughter2BEmcMatchZEdge <<endl;
  cout << " Daughter2 BTOF Projected Channel ID : " << mDaughter2BTofProjChannel << " Daughter2 BTOF Cell LocalY : " << mDaughter2BTofCellLocalY << endl;
  cout << " Daughter2 BTOF Cell LocalZ : " << mDaughter2BTofCellLocalZ << " Daughter2 BTOF PathLength : " <<  mDaughter2BTofPathLength <<endl;
  cout << " Daughter2 BTOF beta : " << mDaughter2Beta << " Daughter2 BTOF tof : " <<  mDaughter2Tof <<endl;

}


ClassImp(StMuHltDiElectron)
