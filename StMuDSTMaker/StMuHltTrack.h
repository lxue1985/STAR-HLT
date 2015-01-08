/*!
 * \class StMuHltTrack 
 * \author Liang Xue , Nov 2010
 */



/** @class StMuHltTrack
 * Class holding the hlt online track properties
 * All unitis are in standard STAR units: GeV,cm
 **/

#ifndef StMuHltTrack_h
#define StMuHltTrack_h


#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StThreeVectorD.hh"

#include "TObject.h"


class StHltTrack;

class StMuHltTrack : public TObject {


	public:
	
	    StMuHltTrack() {/*no-op*/};       // Default constructor
		StMuHltTrack(const StHltTrack*);  // Constructor from StHltTrack
    	~StMuHltTrack(){};
		int id() const;                   // Returns the track id
		int flag() const;         // Returns flag
		int innerMostRow() const;         // Return most inner row.
		int outerMostRow() const;         // Return most outer row.
		int reserved() const;
		int nHits() const;       // Return total number of hits on track.
		int nHitsDedx() const;   // Return number of hits used for dEdx. 
		int q() const;                    // Return charge.
		double chi2xy() const;            // Returns chi2xy of fit.      
		double chi2z() const;             // Returns chi2z of fit.      
		double dEdx() const;              // Returns measured dE/dx value.
		double pt() const;                // Returns pT at point (r,phi,z).
		double phi0() const;               // Returns phi.
		double psi() const;               // Returns psi.
		double r0() const;                // Returns r0.
		double tanl() const;              // Returns tanl.
		double z0() const;                // Returns z0.
		double length() const;            // Returns length of track (cm).
		double dpt() const;               // Returns delta pt.
		double dpsi() const;              // Returns delta psi.
		double dz0() const;               // Returns delta z0.
		double dtanl() const;             // Returns delta tanl.

		const StThreeVectorF &p() const;                 // Returns 3-momentum.
		int index2Global() const;                        // Returns index of associated global track.
		int index2Primary() const;                       // Returns index of associated primary track.
		int index2BEmcTowerHit() const;                  // Returns index of associated bemc tower hit.
		int index2BTofHit() const;                       // Returns index of associated btof hit.

		virtual void Print(Option_t* option = "") const; // Print track info



	protected:

        Int_t mId;            //primary key
		Int_t mFlag;
		Int_t mInnerMostRow;
		Int_t mOuterMostRow;
		Int_t mReserved;
		Int_t mNHits;
		Int_t mNHitsDedx;       
		Int_t mQ;            //charge
		Float_t mChiSqXY;     //chi squared of the momentum fit
		Float_t mChiSqZ;
		Float_t mdEdx;
		Float_t mPt;          //pt at (r, phi, z)
		Float_t mPhi0;        //azimuthal angle of point where parameters are given
		Float_t mPsi;         //azimuthal angle of the momentum at (r,...
		Float_t mR0;          //r (in cyl. coord.) for point where parameters are given
		Float_t mTanl;        //tg of the dip angle at (r,phi,z)
		Float_t mZ0;          //z coordinate of point where parameters are given
		Float_t mLength;
		Float_t mDpt;
		Float_t mDpsi;
		Float_t mDz0;
		Float_t mDtanl;

		StThreeVectorF mP;
		Int_t mIndex2Global;
		Int_t mIndex2Primary;
		Int_t mIndex2BEmcTowerHit;
		Int_t mIndex2BTofHit;

		StThreeVectorD momentumAtPrimaryVertex(const StHltTrack*) const;
		void setIndex2Global(int i) {mIndex2Global=i;}             // Set index of associated global track.
		void setIndex2Primary(int i) {mIndex2Primary=i;}           // Set index of associated primary track.
		void setIndex2BEmcTowerHit(int i) {mIndex2BEmcTowerHit=i;} // Set index of associated bemc tower hit.
		void setIndex2BTofHit(int i) {mIndex2BTofHit=i;}           // Set index of associated btof hit.

		ClassDef(StMuHltTrack,1)
};

inline int StMuHltTrack::id() const {return mId;}
inline int StMuHltTrack::flag() const {return mFlag;}
inline int StMuHltTrack::innerMostRow() const {return mInnerMostRow;}
inline int StMuHltTrack::outerMostRow() const {return mOuterMostRow;}
inline int StMuHltTrack::reserved() const {return mReserved;}
inline int StMuHltTrack::nHits() const {return mNHits;}
inline int StMuHltTrack::nHitsDedx() const {return mNHitsDedx;}
inline int StMuHltTrack::q() const {return mQ;}
inline double StMuHltTrack::chi2xy() const {return mChiSqXY;}
inline double StMuHltTrack::chi2z() const {return mChiSqZ;}
inline double StMuHltTrack::dEdx() const {return mdEdx;}
inline double StMuHltTrack::pt() const {return mPt;}
inline double StMuHltTrack::phi0() const {return mPhi0;}
inline double StMuHltTrack::psi() const {return mPsi;}
inline double StMuHltTrack::r0() const {return mR0;}
inline double StMuHltTrack::tanl() const {return mTanl;}
inline double StMuHltTrack::z0() const {return mZ0;}
inline double StMuHltTrack::length() const {return mLength;}
inline double StMuHltTrack::dpt() const {return mDpt;}
inline double StMuHltTrack::dpsi() const {return mDpsi;}
inline double StMuHltTrack::dz0() const {return mDz0;}
inline double StMuHltTrack::dtanl() const {return mDtanl;}

inline const StThreeVectorF &StMuHltTrack::p() const {return mP;}
inline int StMuHltTrack::index2Global() const {return mIndex2Global;}
inline int StMuHltTrack::index2Primary() const {return mIndex2Primary;}
inline int StMuHltTrack::index2BEmcTowerHit() const {return mIndex2BEmcTowerHit;}
inline int StMuHltTrack::index2BTofHit() const {return mIndex2BTofHit;}


#endif

