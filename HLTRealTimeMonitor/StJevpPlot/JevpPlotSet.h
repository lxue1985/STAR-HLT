#ifndef _JEVPPLOTSET_H_
#define _JEVPPLOTSET_H_

// This is the user class for creating plot's
// 

#include <TROOT.h>
#include <TSocket.h>
#include "JevpPlot.h"
#include "DAQ_READER/daqReader.h"

class JevpPlotSet : public TObject {
  
 public:
 
  JevpPlotSet();
  // Plot management
  //
  // Generally, one adds plots only at initialization
  // and accesses the plots using getPlot()
  //
  // Importantly, the memory still owned by the calling class
  int addPlot(JevpPlot *plot);
  JevpPlot *getPlot(char *name);

  void removePlot(char *name);
  int getNumberOfPlots();
  JevpPlot *getPlotByIndex(int i);
  void dump();

  // initialize() is called once at the begininning of the program
  // the arguments are the command line arguments of the program
  virtual void initialize(int argc, char *argv[]);

  // startrun is called once at the begining of the run
  virtual void startrun(daqReader *rdr);
  
  // stoprun is called once at the end of the run
  virtual void stoprun(daqReader *rdr);
  void _stoprun(daqReader *rdr);

  // event is called once for each event
  virtual void event(daqReader *rdr);
  //ntuple is called onec for each event
  virtual void ntuple(daqReader *rdr);

  // this is the selection criteria for events
  virtual int selectEvent(daqReader *rdr);

  // this is the selection criteria for runs
  virtual int selectRun(daqReader *rdr);

  void Main(int argc, char *argv[]);

  int send(TObject *msg); // used only for special cases...such as base server run info
 private:
  TSocket *socket;

  TList plots;    // The plots built
  
  unsigned int current_run;

  daqReader *reader;
  char *diska;      // event pool path
  char *daqfile;    // data file / null for live
  char *server;     // server
  int serverport;   // server port
  char *pdf;        // direct pdf file output
  char *loglevel;
  int update_time;

  int parseArgs(int argc, char *argv[]);   
  int connect(char *host, int port);
  int updatePlots();
  void writePdfFile();

  ClassDef(JevpPlotSet, 1);
};

#endif
