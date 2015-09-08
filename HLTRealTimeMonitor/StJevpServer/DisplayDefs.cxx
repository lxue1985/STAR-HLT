#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "DisplayDefs.h"
#include <libxml/xmlreader.h>

HistogramDescriptor::HistogramDescriptor(char *name) {
  next = NULL;
  isLogX=-1;
  isLogY=-1;
  this->name = (char *)malloc(strlen(name)+1);
  strcpy(this->name, name);
}

HistogramDescriptor::~HistogramDescriptor() {
  if(next) delete next;
  free(name);
}

TabDescriptor::TabDescriptor(char *name) {  
  next = NULL;
  child = NULL;
  next = NULL;
  this->name = (char *)malloc(strlen(name)+1);
  strcpy(this->name, name);
}

TabDescriptor::~TabDescriptor() {
  if(child) delete child;
  if(next) delete next;
  free(name);
}

void DisplayDefinition::chomp(char *to, char *str, int max)
{
  while(isspace(*str)) {
    if(*str == '\0') {
      to[0] = '\0';
      return;
    }
    str++;
  }

  strncpy(to,str,max);
  to[max-1] = '\0';

  int n=strlen(to);
  for(int i=n-1;i>=0;i--) {
    if(isspace(to[i])) to[i] = '\0';
    else return;
  }
}

// have already read type=ELEMENT, name=histogram
HistogramDescriptor *DisplayDefinition::ReadHistoDescriptor(xmlTextReaderPtr reader)
{
  char name[100];
  int xlog=-1;
  int ylog=-1;
  name[0] = '\0';

  int startdepth = xmlTextReaderDepth(reader);

  int ret = xmlTextReaderRead(reader);
  while(ret == 1) {
    char *element = (char *)xmlTextReaderConstName(reader);
    int element_len = strlen(element);
    int type = xmlTextReaderNodeType(reader);
    int depth = xmlTextReaderDepth(reader);

    if((type == XML_READER_TYPE_TEXT) && (depth == startdepth+1)) {  // name
      chomp(name, (char *)xmlTextReaderConstValue(reader),100);
    }

    if(type == XML_READER_TYPE_ELEMENT) {
  
      if((element_len == 9) && (memcmp("xlogscale", element, 9)==0)) {
	xlog = 1;
      }
      
      if((element_len == 9) && (memcmp("ylogscale", element, 9)==0)) {
	ylog = 1;
      }

    }
    
    if(type == XML_READER_TYPE_END_ELEMENT) {
      if((element_len == 9) && (memcmp(element,"histogram",9) == 0)) {
	break;
      }
    }

    ret = xmlTextReaderRead(reader);
  }

  // Create histo...
  if(ret != 1) {
    return NULL;
  }

  HistogramDescriptor *hd = new HistogramDescriptor(name);
  hd->isLogX = xlog;
  hd->isLogY = ylog;

  return hd;
}

// have already read type=ELEMENT, name=Canvas
CanvasDescriptor *DisplayDefinition::ReadCanvasDescriptor(xmlTextReaderPtr reader)
{
  char tmp[100];
  int wide=1;
  int deep=1;
  HistogramDescriptor *first = NULL;
  HistogramDescriptor *curr = NULL;

  int ret = xmlTextReaderRead(reader);
  while(ret == 1) {
    char *element = (char *)xmlTextReaderConstName(reader);
    int element_len = strlen(element);
    int type = xmlTextReaderNodeType(reader);

    if(type == XML_READER_TYPE_ELEMENT) {
  
      if((element_len == 4) && (memcmp("wide", element, 4)==0)) {
	ret = xmlTextReaderRead(reader);
	type = xmlTextReaderNodeType(reader);
	//printf("checking wide %d-%d\n",type,XML_READER_TYPE_ELEMENT);
	if(type == XML_READER_TYPE_TEXT) {
	  chomp(tmp, (char *)xmlTextReaderConstValue(reader),100);
	  wide = atoi(tmp);
	  //printf("wide = %d\n",wide);
	}
      }
      
      if((element_len == 4) && (memcmp("deep", element, 4)==0)) {
	ret = xmlTextReaderRead(reader);
	type = xmlTextReaderNodeType(reader);
	if(type == XML_READER_TYPE_TEXT) {
	  chomp(tmp, (char *)xmlTextReaderConstValue(reader),100);
	  deep = atoi(tmp);
	}
      }

      if((element_len == 9) && (memcmp("histogram", element, 9)==0)) {
	HistogramDescriptor *hd = ReadHistoDescriptor(reader);
	if(first == NULL) {
	  first = hd;
	  curr = hd;
	}
	else {
	  curr->next = hd;
	  curr = hd;
	}
      }
    }
    
    if(type == XML_READER_TYPE_END_ELEMENT) {
      if((element_len == 6) && (memcmp("canvas",element,6) == 0)) {
	break;
      }
    }

    ret = xmlTextReaderRead(reader);
  }

  // Create Canvas...
  if(ret != 1) {
    return NULL;
  }

  CanvasDescriptor *cd = new CanvasDescriptor();
  cd->wide = wide;
  cd->deep = deep;
  cd->first_histo = first;

  return cd;
}


// have already read type=ELEMENT, name=tab
TabDescriptor *DisplayDefinition::ReadTabDescriptor(xmlTextReaderPtr reader)
{
  char name[100];
  TabDescriptor *child = NULL;
  TabDescriptor *curr = NULL;
  CanvasDescriptor *cd = NULL;

  name[0] = '\0';

  // Get "tab" depth
  int startdepth = xmlTextReaderDepth(reader);

  int ret = xmlTextReaderRead(reader);
  while(ret == 1) {
    char *element = (char *)xmlTextReaderConstName(reader);
    int element_len = strlen(element);
    int type = xmlTextReaderNodeType(reader);

    int depth = xmlTextReaderDepth(reader);

    if((type == XML_READER_TYPE_TEXT) && (depth == startdepth+1)) {  // name
      chomp(name, (char *)xmlTextReaderConstValue(reader),100);
    }

    if(type == XML_READER_TYPE_ELEMENT) {
  
      if((element_len == 3) && (memcmp("tab", element, 3)==0)) {
	TabDescriptor *td = ReadTabDescriptor(reader);
	if(!child) {
	  child = td;
	  curr = td;
	}
	else {
	  curr->next = td;
	  curr = td;
	}
      }
      
      if((element_len == 6) && (memcmp("canvas", element, 6)==0)) {
	cd = ReadCanvasDescriptor(reader);
      }
    }
    
    if(type == XML_READER_TYPE_END_ELEMENT) {
      if((element_len == 3) && (memcmp("tab",element, 3) == 0)) {
	break;
      }
    }

    ret = xmlTextReaderRead(reader);
  }

  // Create histo...
  if(ret != 1) {
    return NULL;
  }

  TabDescriptor *td = new TabDescriptor(name);
  td->child = child;
  td->canvas = cd;
  
  return td;
}

int DisplayDefinition::ReadBuff(char *buff, int len)
{
  char templat[100];
  strcpy(templat, "/tmp/evpreader.XXXXXX");

  int fd = mkstemp(templat);
  int ret;

  while(len > 0) {
    ret = write(fd, buff, len);
    if(ret <= 0) { return -1; };
    buff += ret;
    len -= ret;
  }

  close(fd);
  
  ret = Read(templat);
  
  unlink(templat);
  
  return ret;
}

int DisplayDefinition::Read(char *fn)
{
  LIBXML_TEST_VERSION

  xmlTextReaderPtr reader;
  int err;
  TabDescriptor *curr=NULL;
  int ret;

  reader = xmlReaderForFile(fn, NULL, XML_PARSE_NOBLANKS);
  if(!reader) {
    printf("No reader\n");
    err = -1;
  }
  else {
    ret = xmlTextReaderRead(reader);
    while(ret == 1) {
      char *element = (char *)xmlTextReaderConstName(reader);
      int element_len = strlen(element);
      int type = xmlTextReaderNodeType(reader);

      if(type == XML_READER_TYPE_ELEMENT) {
	
	if((element_len == 3) && (memcmp("tab", element, 3)==0)) {
	  if(!tabs) {
	    tabs = ReadTabDescriptor(reader);
	    curr = tabs;
	  }
	  else {
	    curr->next = ReadTabDescriptor(reader);
	    curr = curr->next;
	  }
	}
      }
      ret = xmlTextReaderRead(reader);
    }
  }
  xmlCleanupParser();


  // Now cache the old xml file text...
  int fd = open(fn, O_RDONLY);
  if(fd < 0) {
    char tmp[100];
    printf("pwd = %s\n",getwd(tmp));
    printf("can't open %s\n",fn);
    return -1;
  }

  struct stat sbuff;
  ret = fstat(fd, &sbuff);
  int len = sbuff.st_size;

  if(textBuff) free(textBuff);
  textBuff = (char *)malloc(len);
  textBuffLen = len;

  char *buff = textBuff;
  while(len) {
    ret = read(fd, buff, len);
    len -= ret;
    buff += ret;
  }
  
  close(fd);
  //

  return 0;
}

  
void DisplayDefinition::dumptab(TabDescriptor *tab, int indent)
{
  for(int i=0;i<indent;i++) printf(" ");
  
  printf("%s\n",tab->name);
  if(tab->canvas) {
    HistogramDescriptor *hd = tab->canvas->first_histo;

    while(hd) {
      for(int i=0;i<indent+4;i++) printf(" ");
      printf("%s\n",hd->name);
      hd = hd->next;
    }
  }

  if(tab->child) {
    dumptab(tab->child, indent+8);
  }
  
  if(tab->next) {
    dumptab(tab->next, indent);
  }
}


#define TAB_BASE 100

u_int DisplayDefinition::getTabBase()
{
  return TAB_BASE;
}

// Gets the multiplier to access the final
u_int DisplayDefinition::getTabDepthMult(u_int idx)
{
  u_int depth = 1;
  idx /= TAB_BASE;
  while(idx) {
    idx = idx / TAB_BASE;
    depth = depth * TAB_BASE;
  }
  return depth;
}

u_int DisplayDefinition::getTabNextIdx(u_int idx)
{
  u_int m = getTabDepthMult(idx);
  idx += m;
  return idx;
}

u_int DisplayDefinition::getTabChildIdx(u_int idx)
{
  u_int m = getTabDepthMult(idx);
  m *= TAB_BASE;
  idx += m;
  return idx;
}

u_int DisplayDefinition::getTabIdxAtDepth(u_int idx, u_int depth)
{
  while(depth) {
    idx /= TAB_BASE;
    depth--;
  }

  return idx % TAB_BASE;
}

// Returns either a tab or a canvas descriptor....
//
// TabDescriptor *tab;     // if *isCanvas=0
// CanvasDescriptor *can;  // if *isCanvas=1
//
void *DisplayDefinition::getTab(u_int combo_index, int *isCanvas)
{

  printf("-------getTabCall idx = %d------------\n",combo_index);
  TabDescriptor *tab = tabs;
  
  *isCanvas = 0;



  for(int depth=0;;depth++) {
    int idx = getTabIdxAtDepth(combo_index, depth);
    int next_idx = getTabIdxAtDepth(combo_index, depth+1);

    printf("idx(%d @ %d)=%d   next=%d\n",combo_index,depth,idx,next_idx);
    
    int x;
    for(x=1;x<idx;x++) {
      printf("---->horizontal search[%d / %d] looking at: %s\n",x,idx,tab ? tab->name : "null");
      if(tab == NULL) return NULL;
      tab = tab->next;
    }
    
    if(tab == NULL) {
      return NULL;
    }

    printf("---->horizonatl search[%d / %d] found %s\n",x,idx,tab->name);
    
    if(next_idx == 0) {
      printf("Got a tab...%s\n",tab->name);
      return tab;
    }

    if(tab->canvas) {   // we are done... one way or the other....

      printf("Got a canvas\n");
      // check next index...
      int idxp1 = getTabIdxAtDepth(combo_index, depth+1);
      int idxp2 = getTabIdxAtDepth(combo_index, depth+2);
      
      if((idxp1 == 1) && (idxp2 == 0)) {
	*isCanvas = 1;
	return tab->canvas;
      }
      else {
	return NULL;
      }
    }
    
    tab = tab->child;  
    if(!tab) {
      printf("a tab with no canvas\n");
      return NULL;    // So this tree had no canvas?
    }
  }
}


TabDescriptor *DisplayDefinition::GetTab(int i, int j)
{
  TabDescriptor *tab = tabs;
  
  for(int x=0;x<i;x++) {
    if(tab == NULL) return NULL;
    tab = tab->next;
  }

  if(j == -1) return tab;

  tab = tab->child;
  for(int x=0;x<j;x++) {
    if(tab == NULL) return NULL;
    tab = tab->next;
  }
  
  return tab;
}

HistogramDescriptor *DisplayDefinition::GetHist(int i, int j, int k)
{
  TabDescriptor *tab = GetTab(i,j);
  if(!tab) return NULL;
  
  if(!tab->canvas) return NULL;
  
  HistogramDescriptor *hd = tab->canvas->first_histo;
  
  for(int x=0;x<k;x++) {
    if(hd == NULL) return NULL;
    hd = hd->next;
  }

  return hd;
}
    
    
