#include <TROOT.h>
#include <TServerSocket.h>
#include <TSocket.h>
#include <TMessage.h>
#include <TMonitor.h>
#include <TClass.h>

#include "EvpMessage.h"
#include "StJevpPool/StJevpPlot/JevpPlot.h"
#include "DisplayDefs.h"
#include "RunStatus.h"
#include "EvpConstants.h"

class JevpServer {
 public:
  TServerSocket *ssocket;
  TMonitor *mon;
  char *refplotdir;
  DisplayDefinition *display;

  int run;  // will be run info...

  TList plots;
  RunStatus *status;

  JevpServer() {
    ssocket = NULL;
    mon = NULL;
    status = NULL;
    refplotdir = EVP_REFPLOT_DIR;
  };
  
  int init(int port, char *displayDefFile);
  void getMessage();
  void shiftRefPlotsUp(char *name, int idx);
  void shiftRefPlotsDown(char *name, int idx);
  int getMaxRef(char *name);
  void deleteReferencePlot(char *name, int idx);
  void handleEvpMessage(TSocket *s, EvpMessage *msg);
  void handleEvpPlot(TSocket *s, JevpPlot *plot);
  void saveReferencePlot(JevpPlot *plot);
  JevpPlot *getPlot(char *name);
  void handleGetPlot(TSocket *s, char *argstring);
  void handleSwapRefs(char *args);
  void clearForNewRun();
  void dump();
  int updateDisplayDefs(char *fn);

  char *getParamFromString(char *dest, char *source, char *param=NULL);
  void writePdf(char *fn);
};
