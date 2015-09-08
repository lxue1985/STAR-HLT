#include <stdio.h>
#include <RTS/src/SUNRT/EthernetServer.hh>
#include <StPDFUtilities/PdfIndex.hh>
#include <TError.h>
#include "DisplayDefs.h"
#include "JevpServer.h"
#include "EvpConstants.h"

int main(int argc, char *argv[])
{
  gErrorIgnoreLevel = kBreak;   // suppress root messages...

  JevpServer serv;
  serv.init(EVP_PORT,"/home/jml/cvs2/offline/users/jml/StJevpPool/StJevpData/HistoDefs.txt");

  //char *pdfFile="pdf.pdf";
  // serv.writePdf(pdfFile);
  for(;;) {
    serv.getMessage();
  }

}
