/*!
 * \class StMuHltDiElectron 
 * \author Liang Xue , Nov 2010
 */




/** @class StMuHltDiElectron
 * Class holding the hlt online di-electron properties
 * All unitis are in standard STAR units: GeV,cm
 **/

#ifndef StMuHltDiElectron_h
#define StMuHltDiElectron_h

#include "TObject.h"

#include "StMuHltTrack.h"
#include "StMuHltBTofHit.h"
#include "StMuHltBEmcTowerHit.h"

class StHltDiElectron;

class StMuHltDiElectron : public TObject {

	public:

		StMuHltDiElectron() {/*no-op*/};       // Default constructor
		StMuHltDiElectron(const StHltDiElectron*); // Constructor from StHltDiElectron
		~StMuHltDiElectron() {};


		// Daughter1 methods
		
		double daughter1bEmcMatchPhiDiff() const;     // Returns daughter1 bEmc Phi Difference 
		double daughter1bEmcMatchZEdge() const;       // Returns daughter1 bEmc Matched z Edge
		float   daughter1bTofProjChannel() const;      // Returns daughter1 bTof Project Channel ID
		float   daughter1bTofCellLocalY() const;       // Returns daughter1 bTof Cell LocalY
		float   daughter1bTofCellLocalZ() const;       // Returns daughter1 bTof Cell LocalZ
		float   daughter1bTofPathLength() const;       // Returns daughter1 bTof Path Length
		float   daughter1beta() const;                 // Returns daughter1 bTof beta
		float   daughter1tof() const;                  // Returns daughter1 bTof Time of Flight

		int    daughter1index2Global() const;         // Returns daughter1 index of associated global track.
		int    daughter1index2Primary() const;        // Returns daughter1 index of associated primary track.
		int    daughter1index2BTof() const;           // Returns daughter1 index of associated btof hit.
		int    daughter1index2BEmc() const;           // Returns daughter1 index of associated bemc tower hit.
		
		const StMuHltTrack& daughter1primaryTrack() const;        // Returns object of daughter1 associated global track.
		const StMuHltTrack& daughter1globalTrack() const;         // Returns object of daughter1 associated primary track.
		const StMuHltBTofHit& daughter1bTofHit() const;           // Returns object of daughter1 associated btof hit.
		const StMuHltBEmcTowerHit& daughter1bEmcTowerHit() const; // Returns object of daughter1 associated bemc hit.

		
		
		// Daughter2 methods
		
		double daughter2bEmcMatchPhiDiff() const;     // Returns daughter2 bEmc Phi Difference 
		double daughter2bEmcMatchZEdge() const;       // Returns daughter2 bEmc Matched z Edge
		float   daughter2bTofProjChannel() const;      // Returns daughter2 bTof Project Channel ID
		float   daughter2bTofCellLocalY() const;       // Returns daughter2 bTof Cell LocalY
		float   daughter2bTofCellLocalZ() const;       // Returns daughter2 bTof Cell LocalZ
		float   daughter2bTofPathLength() const;       // Returns daughter2 bTof Path Length
		float   daughter2beta() const;                 // Returns daughter2 bTof beta
		float   daughter2tof() const;                  // Returns daughter2 bTof Time of Flight

		int    daughter2index2Global() const;         // Returns idaughter2 ndex of associated global track.
		int    daughter2index2Primary() const;        // Returns idaughter2 ndex of associated primary track.
		int    daughter2index2BTof() const;           // Returns idaughter2 ndex of associated btof hit.
		int    daughter2index2BEmc() const;           // Returns daughter2 index of associated bemc tower hit.
		
		const StMuHltTrack& daughter2primaryTrack() const;        // Returns object of daughter2 associated global track.
		const StMuHltTrack& daughter2globalTrack() const;         // Returns object of daughter2 associated primary track.
		const StMuHltBTofHit& daughter2bTofHit() const;           // Returns object of daughter2 associated btof hit.
		const StMuHltBEmcTowerHit& daughter2bEmcTowerHit() const; // Returns object of daughter2 associated bemc hit.

		virtual void Print(Option_t* option = "") const;  // Print track info


	protected:
	    
		// Daughter1 hlt information from node

		Double_t mDaughter1BEmcMatchPhiDiff;
		Double_t mDaughter1BEmcMatchZEdge;
		Float_t  mDaughter1BTofProjChannel;
		Float_t  mDaughter1BTofCellLocalY;
		Float_t  mDaughter1BTofCellLocalZ;
		Float_t  mDaughter1BTofPathLength;
		Float_t  mDaughter1Beta;
	    Float_t  mDaughter1Tof;

		Int_t mDaughter1Index2Global;
		Int_t mDaughter1Index2Primary;
		Int_t mDaughter1Index2BTof;
		Int_t mDaughter1Index2BEmc;

		StMuHltTrack mDaughter1GlobalTrack;          ///< associat object to daughter1 global track
		StMuHltTrack mDaughter1PrimaryTrack;         ///< associat object to daughter1 primary track
		StMuHltBTofHit mDaughter1BTofHit;            ///< associat object to daughter1 bTofHit
		StMuHltBEmcTowerHit mDaughter1BEmcTowerHit;  ///< associat object to daughter1 bEmcTowerHit


		// Daughter2 hlt information from node

		Double_t mDaughter2BEmcMatchPhiDiff;
		Double_t mDaughter2BEmcMatchZEdge;
		Float_t  mDaughter2BTofProjChannel;
		Float_t  mDaughter2BTofCellLocalY;
		Float_t  mDaughter2BTofCellLocalZ;
		Float_t  mDaughter2BTofPathLength;
		Float_t  mDaughter2Beta;
	    Float_t  mDaughter2Tof;

		Int_t mDaughter2Index2Global;
		Int_t mDaughter2Index2Primary;
		Int_t mDaughter2Index2BTof;
		Int_t mDaughter2Index2BEmc;

		StMuHltTrack mDaughter2GlobalTrack;          ///< associat object to daughter2 global track
		StMuHltTrack mDaughter2PrimaryTrack;         ///< associat object to daughter2 primary track
		StMuHltBTofHit mDaughter2BTofHit;            ///< associat object to daughter2 bTofHit
		StMuHltBEmcTowerHit mDaughter2BEmcTowerHit;  ///< associat object to daughter2 bEmcTowerHit

		// Di-Electron information
		Float_t mPt;
		Float_t mPsi;
		Float_t mTanl;
		Float_t mInvariantMass;

		void setDaughter1Index2Global(int i) {mDaughter1Index2Global=i;}// Set index of associated global track.
		void setDaughter1Index2Primary(int i) { mDaughter1Index2Primary=i;}// Set index of associated primary track.
		void setDaughter1Index2BTof(int i) { mDaughter1Index2BTof=i;} // Set index of associated btof hit.
		void setDaughter1Index2BEmc(int i) { mDaughter1Index2BEmc=i;} // Set index of associated bemc tower hit.

		void setDaughter2Index2Global(int i) {mDaughter2Index2Global=i;}// Set index of associated global track.
		void setDaughter2Index2Primary(int i) { mDaughter2Index2Primary=i;}// Set index of associated primary track.
		void setDaughter2Index2BTof(int i) { mDaughter2Index2BTof=i;} // Set index of associated btof hit.
		void setDaughter2Index2BEmc(int i) { mDaughter2Index2BEmc=i;} // Set index of associated bemc tower hit.

		ClassDef(StMuHltDiElectron,1)
};


inline double StMuHltDiElectron::daughter1bEmcMatchPhiDiff() const {return mDaughter1BEmcMatchPhiDiff;}
inline double StMuHltDiElectron::daughter1bEmcMatchZEdge() const {return mDaughter1BEmcMatchZEdge;}
inline float   StMuHltDiElectron::daughter1bTofProjChannel() const {return mDaughter1BTofProjChannel;}
inline float   StMuHltDiElectron::daughter1bTofCellLocalY() const {return mDaughter1BTofCellLocalY;}
inline float   StMuHltDiElectron::daughter1bTofCellLocalZ() const {return mDaughter1BTofCellLocalZ;}
inline float   StMuHltDiElectron::daughter1bTofPathLength() const {return mDaughter1BTofPathLength;}
inline float   StMuHltDiElectron::daughter1beta() const {return mDaughter1Beta;}
inline float   StMuHltDiElectron::daughter1tof() const {return mDaughter1Tof;}

inline int StMuHltDiElectron::daughter1index2Global() const {return mDaughter1Index2Global;}
inline int StMuHltDiElectron::daughter1index2Primary() const {return mDaughter1Index2Primary;}
inline int StMuHltDiElectron::daughter1index2BTof() const {return mDaughter1Index2BTof;}
inline int StMuHltDiElectron::daughter1index2BEmc() const {return mDaughter1Index2BEmc;}

inline const StMuHltTrack& StMuHltDiElectron::daughter1globalTrack() const { return mDaughter1GlobalTrack; }
inline const StMuHltTrack& StMuHltDiElectron::daughter1primaryTrack() const { return mDaughter1PrimaryTrack; }
inline const StMuHltBTofHit& StMuHltDiElectron::daughter1bTofHit() const { return mDaughter1BTofHit; }
inline const StMuHltBEmcTowerHit& StMuHltDiElectron::daughter1bEmcTowerHit() const { return mDaughter1BEmcTowerHit; }


inline double StMuHltDiElectron::daughter2bEmcMatchPhiDiff() const {return mDaughter2BEmcMatchPhiDiff;}
inline double StMuHltDiElectron::daughter2bEmcMatchZEdge() const {return mDaughter2BEmcMatchZEdge;}
inline float   StMuHltDiElectron::daughter2bTofProjChannel() const {return mDaughter2BTofProjChannel;}
inline float   StMuHltDiElectron::daughter2bTofCellLocalY() const {return mDaughter2BTofCellLocalY;}
inline float   StMuHltDiElectron::daughter2bTofCellLocalZ() const {return mDaughter2BTofCellLocalZ;}
inline float   StMuHltDiElectron::daughter2bTofPathLength() const {return mDaughter2BTofPathLength;}
inline float   StMuHltDiElectron::daughter2beta() const {return mDaughter2Beta;}
inline float   StMuHltDiElectron::daughter2tof() const {return mDaughter2Tof;}

inline int StMuHltDiElectron::daughter2index2Global() const {return mDaughter2Index2Global;}
inline int StMuHltDiElectron::daughter2index2Primary() const {return mDaughter2Index2Primary;}
inline int StMuHltDiElectron::daughter2index2BTof() const {return mDaughter2Index2BTof;}
inline int StMuHltDiElectron::daughter2index2BEmc() const {return mDaughter2Index2BEmc;}

inline const StMuHltTrack& StMuHltDiElectron::daughter2globalTrack() const { return mDaughter2GlobalTrack; }
inline const StMuHltTrack& StMuHltDiElectron::daughter2primaryTrack() const { return mDaughter2PrimaryTrack; }
inline const StMuHltBTofHit& StMuHltDiElectron::daughter2bTofHit() const { return mDaughter2BTofHit; }
inline const StMuHltBEmcTowerHit& StMuHltDiElectron::daughter2bEmcTowerHit() const { return mDaughter2BEmcTowerHit; }

#endif
