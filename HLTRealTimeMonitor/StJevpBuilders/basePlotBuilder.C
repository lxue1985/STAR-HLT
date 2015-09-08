#include <stdio.h>
#include <stdlib.h>

#include "StJevpPool/StJevpPlot/JevpPlotSet.h"
#include "DAQ_READER/daqReader.h"
#include "StJevpPool/StJevpServer/RunStatus.h"

#include <TH1I.h>
#include <TH2F.h>

#include <math.h>

// This is the one PlotSet that is guarenteed to always exist
// It's main purpose is to provide the run information 
// To the server...
//
// It has no plots (currently)
//

class BasePlotSet : public JevpPlotSet {
public:
  RunStatus status;
  int run;

  void initialize(int argc, char *argv[]) {
    run = 0;
    status.setNumPlotSets(1);

    for(int i=1;i<argc;i++) {
      if(memcmp(argv[i],"-numplotsets", 11) == 0) {
	i++;
	status.setNumPlotSets(atoi(argv[i+1]));
      }
    }
  };
  
  void startrun(daqReader *rdr) {
    run = rdr->run;
    status.setRunNumber(rdr->run);
    // status.setEventNumber();   // What is this, what is it for?
    // status.setEventCounter();  // Doesn't make much sense anymore...
    status.setEndOfRun(0);        // What are the values?
    // status.setLiveSource(0);   // not done from here...
    // status.setToken();         // doesn't make sense anymore
    // status.setEmpty();
    // status.setType();
    // status. setStatus()        // the reader status, but why?
    status.setTime(rdr->evt_time);
    // status.setTriggerBits();   // by plot, not here...
    // status.setDetectorBits();  // by plot, not here...
    status.setTriggerBitsRun(rdr->evpgroupsinrun);   // which triggers are in the run?
    status.setDetectorBitsRun(rdr->detsinrun);  // which detectors are in the run?
    // setEndOfRunActionPerformed(); // not here...

    send((TObject *)&status);
  };

  void stoprun(daqReader *rdr) {
    printf("Stopping run #%d\n",run);
    status.setEndOfRun(1);
    send((TObject *)&status);
  };
};


int main(int argc, char *argv[])
{
  BasePlotSet me;
  
  me.Main(argc, argv);
}
