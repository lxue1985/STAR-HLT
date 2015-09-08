#include <ctype.h>
#include <TROOT.h>
#include <TServerSocket.h>
#include <TSocket.h>
#include <TMessage.h>
#include <TMonitor.h>
#include <TClass.h>
#include <TCanvas.h>
#include <StPDFUtilities/PdfIndex.hh>
#include <dirent.h>
#include <TFile.h>

#include "EvpConstants.h"
#include "JevpServer.h"
#include "EvpMessage.h"  
#include "RunStatus.h"

int JevpServer::init(int port, char *displayDefFile) {
  ssocket = new TServerSocket(port,kTRUE,100);
  mon = new TMonitor();

  mon->Add(ssocket);
  
  display = NULL;
  updateDisplayDefs(displayDefFile);

  return 0;
}  

// Parse a string of the form
// "defaultparam param1=x param2=y"
char *JevpServer::getParamFromString(char *dest, char *source, char *param)
{
  char *tmp = dest;
  char *str = source;

  // Find the "param=" and position directly after it...
  if(param != NULL) {
    str = strstr(source, param);
    if(!str) return NULL;

    str += strlen(param);
    if(*str != '=') {
      return NULL;
    }
    str++;
  }
  
  // copy till whitespace or end...
  while((*str != '\0') && !isspace(*str)) {
    *tmp = *str;
    tmp++;
    str++;
  }
  *tmp = '\0';
  return dest;  
}

void JevpServer::getMessage() {
  TMessage *mess;
  TSocket *s;

  printf("Got a message\n");

  s = mon->Select();

  if(s == ssocket) {
    TSocket *nsock = ssocket->Accept();
    TInetAddress adr = nsock->GetInetAddress();
    mon->Add(nsock);
  }
  else {
    // read...
  
    int ret = s->Recv(mess);
    if(ret == 0) {    // Got a disconnection...
      mon->Remove(s);
      delete s;
      delete mess;
      return;
    }

    // Handle control messages...
    if(strcmp(mess->GetClass()->GetName(),"EvpMessage")==0) {

      EvpMessage *msg = (EvpMessage *)mess->ReadObject(mess->GetClass());
	
      handleEvpMessage(s, msg);
	
      delete msg;
    }
    else if (strcmp(mess->GetClass()->GetName(), "RunStatus") == 0) {
      RunStatus *newstat = (RunStatus *)mess->ReadObject(mess->GetClass());

      if(status) {
	delete status;
      }
      
      status = newstat;

      status->dump(); 
    }
    else if (strcmp(mess->GetClass()->GetName(), "JevpPlot")==0) {
      JevpPlot *plot = (JevpPlot *)mess->ReadObject(mess->GetClass());
      handleEvpPlot(s, plot);
    }
    else {
      printf("Got invalid message type: %s\n",mess->GetClass()->GetName());
    }
	      
    delete(mess);
  }    
}

JevpPlot *JevpServer::getPlot(char *name) {
  JevpPlot *curr = (JevpPlot *)plots.First();

  while(curr) {
    if(strcmp(curr->GetPlotName(), name) == 0) {
      return curr;
    }
    curr = (JevpPlot *)plots.After(curr);
  }

  return NULL;
}

void JevpServer::handleSwapRefs(char *name)
{
  char name1[256];
  char name2[256];
  char tmp[256];
  char base[256];
  int idx1, idx2;
  sscanf(name, "%s %d %d", base, &idx1, &idx2);
  
  printf("Swapping %s (%d <--> %d)\n",base,idx1,idx2);
  sprintf(name1, "%s/REF.%s.%d.root",refplotdir, base, idx1);
  sprintf(name2, "%s/REF.%s.%d.root",refplotdir, base, idx2);

  sprintf(tmp, "%s/REF.%s.root.tmp",refplotdir, base);
  rename(name1, tmp);
  rename(name2, name1);
  rename(tmp, name2);
}

int JevpServer::getMaxRef(char *name)
{
  int maxid = 0;
  struct dirent *dirent;

  DIR *dir = opendir(refplotdir);
  if(dir == NULL) {
    printf("Error opening dir %s\n", refplotdir);
    return -1;
  }

  char basename[256];
  sprintf(basename, "REF.%s.", name);
  
  while((dirent = readdir(dir)) != NULL) {
    if(memcmp(basename, dirent->d_name, strlen(basename)) == 0) {
      char *tmp = dirent->d_name;
      tmp += strlen(basename);

      int id = atoi(tmp);     
      if(id > maxid) maxid = id;
    }
  }
  
  closedir(dir);
  return maxid;
}

void JevpServer::shiftRefPlotsUp(char *name, int first)
{
  int max = getMaxRef(name);
  for(int i=max;i>=first;i--) {
    char dst[256];
    char src[256];
    sprintf(dst, "%s/REF.%s.%d.root",refplotdir,name,i+1);
    sprintf(src, "%s/REF.%s.%d.root",refplotdir,name,i);   
    printf("Renaming file: %s --> %s\n", src, dst);
    rename(src, dst);
  }
}

void JevpServer::shiftRefPlotsDown(char *name, int first)
{
  int max = getMaxRef(name);
  for(int i=first;i<=max;i++) {
    char dst[256];
    char src[256];
    sprintf(dst, "%s/REF.%s.%d.root",refplotdir,name,i-1);
    sprintf(src, "%s/REF.%s.%d.root",refplotdir,name,i);

    printf("Renaming file: %s --> %s\n", src, dst);
    rename(src, dst);
  }
}

void JevpServer::deleteReferencePlot(char *name, int refid) {
  char filename[256];

  sprintf(filename, "%s/REF.%s.%d.root",refplotdir,name,refid);
  printf("Deleting file: %s\n", filename);
  unlink(filename);
  shiftRefPlotsDown(name,refid+1);
}

void JevpServer::saveReferencePlot(JevpPlot *plot) {

  char plotname[256];

  if(plot->refid > 0) {
    shiftRefPlotsUp(plot->GetPlotName(), plot->refid);
  }
  
  sprintf(plotname, "%s/REF.%s.%d.root",refplotdir,plot->GetPlotName(), plot->refid);

  // Now actually save plot to the file plotname...
  TFile f(plotname, "new");
  plot->Write();
  f.Close();
}

void JevpServer::handleEvpPlot(TSocket *s, JevpPlot *plot) {

  if(plot->refid != 0) {     // If this is a reference plot 
    saveReferencePlot(plot); // it goes to a disk file...
    return;
  }

  // Otherwise, the plot stays in memory...
  // Remove previous version...
  JevpPlot *curr = (JevpPlot *)plots.First();

  while(curr) {

    if(strcmp(curr->GetPlotName(), plot->GetPlotName()) == 0) {
      plots.Remove(curr);
      delete curr;
      break;
    }

    curr = (JevpPlot *)plots.After(curr);
  }
  
  plots.Add(plot);
}

void JevpServer::handleGetPlot(TSocket *s, char *argstring) 
{
  JevpPlot *plot=NULL;
  char refidstr[20];
  char runidstr[20];
  char plotname[80];

  printf("argstring is %s\n",argstring);
  if(!getParamFromString(plotname, argstring)) {
    printf("No plot indicated in getplot?\n");
    return;
  }
 
  printf("Plotname is %s\n",plotname);

  if(getParamFromString(refidstr, argstring, "refid")) {
    char fn[256];
    sprintf(fn, "%s/REF.%s.%d.root", EVP_REFPLOT_DIR, plotname, atoi(refidstr));
      
    printf("Reading [%s] from file %s\n",plotname, fn);

    TFile *f1 = new TFile(fn);
    if(!f1) {
      printf("Error opening file: %s",fn);
      plot = NULL;
    }
    else {
      //f1->GetObject(plotname, plot);
      f1->GetObject("JevpPlot",plot);
      f1->Close();

      if(plot) {
	printf("Got plot.....xxx\n");
	plot->refid = atoi(refidstr);
      }
      else {
	printf("Didn't get plot %s\n",plotname);
      }
    }
  }
  else if (getParamFromString(runidstr, argstring, "run")) {    
    char fn[256];
    sprintf(fn, "%s/%d.root",EVP_SAVEPLOT_DIR, atoi(runidstr));

    TFile *f1 = new TFile(fn);
    if(!f1) {
      printf("Error opening file: %s",fn);
      plot = NULL;
    } 
    else {
      f1->GetObject(plotname, plot);
      f1->Close();
    }
  }
  else {
    printf("getplot..\n");
    plot = getPlot(plotname);
  }

    
  if(!plot) {
    char tmp[100];
    sprintf(tmp, "No plot %s",plotname);
    EvpMessage m;
    m.setCmd("noplot");
    m.setArgs(tmp);
    TMessage mess(kMESS_OBJECT);
    mess.WriteObject(&m);
    s->Send(mess);
  } else {
    TMessage mess(kMESS_OBJECT);
    mess.WriteObject(plot);
    s->Send(mess);
  }
}

void JevpServer::handleEvpMessage(TSocket *s, EvpMessage *msg)
{
  if(strcmp(msg->cmd, "newrun") == 0) {
    clearForNewRun();
  }
  else if(strcmp(msg->cmd, "dump") == 0) {
    dump();
  }
  else if(strcmp(msg->cmd, "display_desc") == 0) {  // Display Descriptor
    EvpMessage m;
    m.setCmd("xml");
    m.setArgs(display->textBuff);

    TMessage mess(kMESS_OBJECT);
    mess.WriteObject(&m);
    s->Send(mess);
  }
  else if(strcmp(msg->cmd, "stoprun") == 0) {
    writePdf("pdf.pdf");
  }
  else if(strcmp(msg->cmd, "getplot") == 0) {
    handleGetPlot(s,msg->args);
  }
  else if(strcmp(msg->cmd, "swaprefs") == 0) {
    handleSwapRefs(msg->args);
  }
  else if(strcmp(msg->cmd, "deleteplot") == 0) {
    char str[256];
    int idx;
    sscanf(msg->args, "%s %d", str, &idx);
    deleteReferencePlot(str,idx);
  }
  else {
    printf("Unknown command: %s\n",msg->cmd);
  }
}


void JevpServer::clearForNewRun()
{
  // Delete all from histogram list
  // First free the actual histo, then remove the link...
  JevpPlot *h = (JevpPlot *)plots.First();

  while(h) {
    delete h;
    plots.Remove(plots.FirstLink());
    h = (JevpPlot *)plots.First();
  }
}

void JevpServer::dump()
{
  JevpPlot *curr = (JevpPlot *)plots.First();

  int i=0;
  while(curr) {

    printf("Histogram[%d]: %s\n",i, curr->GetPlotName());
    
    i++;
    curr = (JevpPlot *)plots.After(curr);
  } 
}

int JevpServer::updateDisplayDefs(char *fn)
{
  if(display) delete display;
  display = new DisplayDefinition();
  display->Read(fn);
  return 0;
}

void JevpServer::writePdf(char *fn)
{
  int firstdraw=1;

  // PDF file index support...
  PdfIndex index;   
  index_entry *index_tab = NULL;
  int page = 1;

  // TAB's and descriptors...
  int tab;
  int subtab;
  int hist;

  TabDescriptor *td;
  TabDescriptor *std;
  HistogramDescriptor *hd;

  // The ROOT Canvas and pdf names for the drawings...
  TCanvas canvas("c1");
  char firstname[256];
  char lastname[256];

  strcpy(firstname, fn);
  strcat(firstname, "(");
  strcpy(lastname, fn);
  strcat(lastname, ")");

  // Build the structure of the file...
  for(tab=0;;tab++) {
    td = display->GetTab(tab);
    if(td == NULL) break;
    index_tab = index.add(NULL, td->name, page, 0);
    
    for(subtab=0;;subtab++) {
      std = display->GetTab(tab,subtab);
      if(std == NULL) {   // Do we move on?
	break;
	continue;
      }

      index.add(index_tab, std->name, page, 0);
      page++;

      // Build the canvas...
      CanvasDescriptor *cd = std->canvas;

      canvas.Divide(cd->wide,cd->deep);

      for(int x=0;x<cd->wide;x++) {
	for(int y=0;y<cd->deep;y++) {
	  int i=x*cd->deep + y;

	  canvas.cd(i+1);

	  // Find the name for the current subcanvas...
 	  hd = display->GetHist(tab,subtab,i);
	  if(hd == NULL) continue;

	  // Have the histogram name, now find the histogram
	  JevpPlot *plot = getPlot(hd->name);

	  if(plot) {
	    plot->draw();
	  }
	  else {
	    //printf("Can't find plot %s\n",hd->name);
	  }
	}
      }

      if(firstdraw) 
	canvas.Print(firstname,"pdf,Portrait");
      else
	canvas.Print(fn,"pdf,Portrait");
    }
  }
  TCanvas summary("c2");
  summary.Print(lastname, "pdf,Portrait");

  index.CreateIndexedFile(fn,fn);
}
