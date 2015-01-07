/*!
 * \class StHltVpdHit 
 * \author Liang Xue & Aihong Tang, Nov 2010
 */


#ifndef StHltVpdHit_hh
#define StHltVpdHit_hh

#include <Stiostream.h>
#include "StObject.h"
#include "StArray.h"
#include "StEnumerations.h"

class StHltVpdHit : public StObject {
public:
    StHltVpdHit();
    ~StHltVpdHit();

	StBeamDirection direction() const; //east or west
	short channel() const; 
	float tdc() const;
	float tot() const;
	float tof() const;
	float triggerTime() const;
	short module() const;
	short cell() const;


	void setDirection(StBeamDirection);
	void setChannel(short);
	void setTdc(float);
	void setTot(float);
	void setTof(float);
	void setTriggerTime(float);


private:
	StBeamDirection mDirection;
	short mChannel; // total 19 channels. 
	float mTdc;
	float mTot;
	float mTof;
	float mTriggerTime;

        ClassDef(StHltVpdHit,1)
};

inline StBeamDirection StHltVpdHit::direction() const {return mDirection;}
inline short StHltVpdHit::channel() const {return mChannel;}
inline float StHltVpdHit::tdc() const {return mTdc;}
inline float StHltVpdHit::tot() const {return mTot;}
inline float StHltVpdHit::tof() const {return mTof;}
inline float StHltVpdHit::triggerTime() const {return mTriggerTime;}
inline short StHltVpdHit::module() const {return mChannel/6;}
inline short StHltVpdHit::cell() const {return mChannel%6;}


ostream& operator<<(ostream&, const StHltVpdHit&); //Printting operator

#endif

