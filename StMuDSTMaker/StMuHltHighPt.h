/*!
 * \class StMuHltHighPt 
 * \author Liang Xue , Nov 2010
 */


/** @class StMuHltHighPt
 * Class holding the hlt online high pt  properties
 * All unitis are in standard STAR units: GeV,cm
 **/

#ifndef StMuHltHighPt_h
#define StMuHltHighPt_h

#include "TObject.h"

#include "StMuHltTrack.h"
#include "StMuHltBTofHit.h"
#include "StMuHltBEmcTowerHit.h"


class StHltHighPt;

class StMuHltHighPt : public TObject {


	public:

		StMuHltHighPt() {/*no-op*/};       // Default constructor
		StMuHltHighPt(const StHltHighPt*); // Constructor from StHltHighPt
		~StMuHltHighPt() {};

		double bEmcMatchPhiDiff() const;     // Returns bEmc Phi Difference 
		double bEmcMatchZEdge() const;       // Returns bEmc Matched z Edge
		float   bTofProjChannel() const;      // Returns bTof Project Channel ID
		float   bTofCellLocalY() const;       // Returns bTof Cell LocalY
		float   bTofCellLocalZ() const;       // Returns bTof Cell LocalZ
		float   bTofPathLength() const;       // Returns bTof Path Length
		float   beta() const;                 // Returns bTof beta
		float   tof() const;                  // Returns bTof Time of Flight

		int    index2Global() const;         // Returns index of associated global track.
		int    index2Primary() const;        // Returns index of associated primary track.
		int    index2BTof() const;           // Returns index of associated btof hit.
		int    index2BEmc() const;           // Returns index of associated bemc hit.

		const StMuHltTrack& primaryTrack() const;        // Returns object of associated global track.
		const StMuHltTrack& globalTrack() const;         // Returns object of associated primary track.
		const StMuHltBTofHit& bTofHit() const;           // Returns object of associated btof hit.
		const StMuHltBEmcTowerHit& bEmcTowerHit() const; // Returns object of associated bemc hit.
		
		virtual void Print(Option_t* option = "") const;  // Print hlt high pt info
		
		

	protected:
		
		// HighPt information from Node
		Double_t mBEmcMatchPhiDiff;
		Double_t mBEmcMatchZEdge;
		Float_t mBTofProjChannel;
		Float_t mBTofCellLocalY;
		Float_t mBTofCellLocalZ;
		Float_t mBTofPathLength;
		Float_t mBeta;
	    Float_t mTof;
		Int_t mIndex2Global;
		Int_t mIndex2Primary;
		Int_t mIndex2BTofHit;
		Int_t mIndex2BEmcTowerHit;

		StMuHltTrack mGlobalTrack;          ///< associat object to global track
		StMuHltTrack mPrimaryTrack;         ///< associat object to primary track
		StMuHltBTofHit mBTofHit;            ///< associat object to bTofHit
		StMuHltBEmcTowerHit mBEmcTowerHit;  ///< associat object to bEmcTowerHit

		void setIndex2Global(int i) {mIndex2Global=i;} ///< Set index of associated global track.
		void setIndex2Primary(int i) { mIndex2Primary=i;} ///< Set index of associated primary track.
		void setIndex2BTof(int i) { mIndex2BTofHit=i;} ///< Set index of associated btof hit.
		void setIndex2BEmc(int i) { mIndex2BEmcTowerHit=i;} ///< Set index of associated bemc tower hit.

		ClassDef(StMuHltHighPt,1)
};

inline double StMuHltHighPt::bEmcMatchPhiDiff() const {return mBEmcMatchPhiDiff;}
inline double StMuHltHighPt::bEmcMatchZEdge() const {return mBEmcMatchZEdge;}
inline float   StMuHltHighPt::bTofProjChannel() const {return mBTofProjChannel;}
inline float   StMuHltHighPt::bTofCellLocalY() const {return mBTofCellLocalY;}
inline float   StMuHltHighPt::bTofCellLocalZ() const {return mBTofCellLocalZ;}
inline float   StMuHltHighPt::bTofPathLength() const {return mBTofPathLength;}
inline float   StMuHltHighPt::beta() const {return mBeta;}
inline float   StMuHltHighPt::tof() const {return mTof;}

inline int StMuHltHighPt::index2Global() const {return mIndex2Global;}
inline int StMuHltHighPt::index2Primary() const {return mIndex2Primary;}
inline int StMuHltHighPt::index2BTof() const {return mIndex2BTofHit;}
inline int StMuHltHighPt::index2BEmc() const {return mIndex2BEmcTowerHit;}

inline const StMuHltTrack& StMuHltHighPt::globalTrack() const { return mGlobalTrack; }
inline const StMuHltTrack& StMuHltHighPt::primaryTrack() const { return mPrimaryTrack; }
inline const StMuHltBTofHit& StMuHltHighPt::bTofHit() const { return mBTofHit; }
inline const StMuHltBEmcTowerHit& StMuHltHighPt::bEmcTowerHit() const { return mBEmcTowerHit; }


#endif

