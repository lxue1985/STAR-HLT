/*!
 * \class StHltMaker 
 * \author Liang Xue , Nov 2010
 */


#ifndef STAR_StHltMaker
#define STAR_StHltMaker


#include "StRTSBaseMaker.h"

namespace star {
namespace rts {
namespace hlt {
class HLT_EVE;
class HLT_TOF;
class HLT_PVPD;
class HLT_EMC;
class HLT_GT;
class HLT_PT;
class HLT_NODE;
class HLT_HIPT;
class HLT_DIEP;
class HLT_HF;
}}}

using namespace std;
using namespace star;
using namespace star::rts;
using namespace star::rts::hlt;

class StEvent;
class StHltEvent;

class StHltMaker : public StRTSBaseMaker {

 private:

    StEvent*                         mStEvent;                   //! pointer to StEvent
    StHltEvent*                      mStHltEvent;                //! pointer to StHltEvent

	unsigned int mnHighPt;                 // number of HighPt tracks
	unsigned int mnHeavyFragment;          // number of HeavyFragment tracks
	unsigned int mnDielectron;             // number of di-electron pairs
	
	int mHighPtNodeSN[1000];               // series number of HighPt tracks
	int mHeavyFragmentNodeSN[1000];        // series number of HeavyFragment tracks
	int mDaughter1NodeSN[1000];            // series number of daughter1 electron
	int mDaughter2NodeSN[1000];            // series number of daughter2 electron

 protected:
 
    virtual StRtsTable*              GetNextGl3();

	StHltEvent*                      GetHltEvent();                  //! create StHltEvent point

    virtual void ProcessBank(const HLT_EVE  *bank);
    virtual void ProcessBank(const HLT_TOF  *bank);
    virtual void ProcessBank(const HLT_PVPD *bank);
    virtual void ProcessBank(const HLT_EMC  *bank);
    virtual void ProcessBank(const HLT_GT   *bank);
    virtual void ProcessBank(const HLT_PT   *bank);
    virtual void ProcessBank(const HLT_NODE *bank);
    virtual void ProcessBank(const HLT_HIPT *bank);
    virtual void ProcessBank(const HLT_DIEP *bank);
    virtual void ProcessBank(const HLT_HF   *bank);
    
 public:

  StHltMaker(const char *name="HLT") ;

  virtual       ~StHltMaker();
  virtual void  Clear(Option_t *option="");
  virtual Int_t Init();
  virtual Int_t Make();
  virtual Int_t Finish();

  virtual Int_t InitRun  (int runumber); /// Overload empty StMaker::InitRun 

  void  FillNodePointer(StHltEvent*);    /// create pointers in track node
  void  FillHighPt(StHltEvent*);         /// fill high pt information 
  void  FillHeavyFragment(StHltEvent*);  /// fill heavy fragment information`
  void  FillDielectron(StHltEvent*);     /// fill di-electron information
  void  FillTriggerReason(StHltEvent*);  /// fill triggered prticles (highpt heavyfragment di-electron) to trigger reason

  /// Displayed on session exit, leave it as-is please ...
  virtual const char *GetCVS() const {
    static const char cvs[]="Tag $Name:  $ $Id: StHltMaker.h,v 1.1 2009/12/23 21:05:58 fine Exp $ built "__DATE__" "__TIME__ ; 
    return cvs;
  }

  // obtain the whole list of leading edge hits
  // to obtain the published result use StMaker::GetDataSet("pp2ppRawHits");

  ClassDef(StHltMaker,0)   //StAF chain virtual base class for Makers

};


#endif
