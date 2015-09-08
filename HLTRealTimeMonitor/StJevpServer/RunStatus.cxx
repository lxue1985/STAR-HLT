#include "RunStatus.h"

ClassImp(RunStatus) ;

void RunStatus::dump()
{
  printf("Run Status----------------\n");
  printf("run=%d\n",mRun);
  printf("End=%d\n",mEnd);
  printf("Trigger Bits=0x%x\n",mTriggerBitsRun);
  printf("Detector Bits=0x%x\n",mDetectorBitsRun);
  printf("--------------------------\n");
}



/***************************************************************************
 *
 * $Id: RunStatus.cxx,v 1.1 2009/10/20 20:38:09 jml Exp $
 *
 * Author: Frank Laue, laue@bnl.gov
 ***************************************************************************
 *
 * Description:
 *
 ***************************************************************************
 *
 * $Log: RunStatus.cxx,v $
 * Revision 1.1  2009/10/20 20:38:09  jml
 * getting up to date...
 *
 * Revision 1.1  2009/01/23 16:11:06  jeromel
 * Import from online/RTS/src/
 *
 * Revision 1.1  2007/02/27 15:23:39  laue
 * Initial version
 *
 * Revision 1.1  2006/10/04 20:31:16  laue
 * Initial Version
 *
 *
 ***************************************************************************/

