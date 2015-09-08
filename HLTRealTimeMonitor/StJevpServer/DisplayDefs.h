#ifndef _DISPLAY_DEFS_H_
#define _DISPLAY_DEFS_H_

#include <TObject.h>
#include <libxml/xmlreader.h>

// The parameters here can override the actual histograms settings
// These are set from the display file...
class HistogramDescriptor {
 public:
  char *name;
  int isLogX;
  int isLogY;
  HistogramDescriptor *next;
  
  HistogramDescriptor(char *name);
  ~HistogramDescriptor();
};

// The canvas has a list of histograms
class CanvasDescriptor {
 public:
  int wide;
  int deep;
  HistogramDescriptor *first_histo;
  
  int nhistos() {
    int n=0;

    printf("cd->wide=%d deepd=%d\n",wide,deep);
    HistogramDescriptor *hd = first_histo;
    while(hd) {
      n++;
      printf("(%d) hd->name=%s\n",n,hd->name);
      hd = hd->next;
    }
    return n;
  }
      
  CanvasDescriptor() { wide = 1; deep = 1; first_histo=NULL; };
  ~CanvasDescriptor() { if(first_histo) delete first_histo; };
};

// The tabs are a tree of depth 2...
class TabDescriptor {
 public:
  char *name;
  TabDescriptor *child;
  TabDescriptor *next;
  CanvasDescriptor *canvas;
  
  TabDescriptor(char *name);
  ~TabDescriptor();
};


// The DisplayDefinition class has a split personality
//
//  * tabs contains the parsed tree
//  * textBuff contains the raw xml
//
//  The raw xml is sent to the presenter from the server
//  truely, this is a presenter oriented structure,
//  however the server needs it because the server creates the
//  pdf file for the database.
//  
//  given this fact, I want to make the presenter truly independent
//  of the computer, so I need at least one (the server's)
//  display definition available ***without*** and presenter viewable
//  config files, so I let the presenter download the display defs
//  this is better as text than as a tree...
//
//  Of course, the presenter can have an optional override with
//  its own display def...
//
class DisplayDefinition {
 public:
  TabDescriptor *tabs;

  char *textBuff;
  int textBuffLen;
  
  DisplayDefinition() { tabs = NULL; textBuff = NULL;}
  void chomp(char *to, char *from, int max);
  int Read(char *fn);
  int ReadBuff(char *buff, int len);
  HistogramDescriptor *ReadHistoDescriptor(xmlTextReaderPtr reader);
  CanvasDescriptor *ReadCanvasDescriptor(xmlTextReaderPtr reader);
  TabDescriptor *ReadTabDescriptor(xmlTextReaderPtr reader);
  
  TabDescriptor *GetTab(int i, int j=-1);
  HistogramDescriptor *GetHist(int i, int j, int k);
 
  void dumptab(TabDescriptor *tab, int indent);
  void dump() { dumptab(tabs,0); };
  ~DisplayDefinition() { if(tabs) delete tabs; if(textBuff) free(textBuff); };

  static u_int getTabBase();
  static u_int getTabDepthMult(u_int idx);
  static u_int getTabNextIdx(u_int idx);
  static u_int getTabChildIdx(u_int idx);
  static u_int getTabIdxAtDepth(u_int idx, u_int depth);
  void *getTab(u_int combo_index, int *isCanvas);
};
  

  
#endif
