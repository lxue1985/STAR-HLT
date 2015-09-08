#include <stdio.h>
#include <RTS/src/SUNRT/EthernetServer.hh>
#include <StPDFUtilities/PdfIndex.hh>
#include "DisplayDefs.h"
#include "JevpServer.h"

void create_hist();

class myserver :  public EthernetServer 
{
public:
  myserver() : EthernetServer(3900, 4, 20) {};

  void processPacket(int address, char *packet) {};

  void processPacket(int address, char *packet, int socket) {
    printf("hello\n");
  }
};


int main(int argc, char *argv[])
{
  //printf("Hi Jeff\n");

  //create_hist();

  //printf("Done\n");

  //myserver serv;
  //serv.run();

//   DisplayDefinition display;
  
//   display.Read("/home/jml/cvs2/offline/users/jml/StJevpPool/StJevpData/HistoDefs.txt");
//   display.dumptab(display.tabs, 0);

//   HistogramDescriptor *hd = display.GetHist(1,2,1);
//   printf("\n\n----> [1,2,1] --> %s\n",hd->name);

  JevpServer serv;
  serv.init(EVP_PORT,"/home/jml/cvs2/offline/users/jml/StJevpPool/StJevpData/HistoDefs.txt");
  for(;;) {
    serv.getMessage();
  }
}
