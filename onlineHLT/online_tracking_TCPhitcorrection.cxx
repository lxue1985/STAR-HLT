/***********************************************************************
 *
 * $Id: online_tracking_TCPhitcorrection.cxx,v 1.2 2009/11/15 21:19:54 qiuh Exp $
 *
 * Author: QIU Hao   12/3/2008
 *
 ***********************************************************************
 *
 * Description: a map to correct tpc hit distortion
 *
 ***********************************************************************
 */

#include <iostream>
#include <fstream>
#include "stdio.h"
#include "FtfGeneral.h"
#include "online_tracking_TCPhitcorrection.h"
#include   <stdio.h>   
#include   <stdlib.h>                                                




  
//   StTpcDb* tpcDbIn = gStTpcDb;
online_tracking_TCPhitcorrection::online_tracking_TCPhitcorrection()
{

}
  


void online_tracking_TCPhitcorrection::loadMap(const char* fileName)
{
  /*
  ifstream ifstr(fileName);
  for(int iSector=0; iSector<24; iSector++)
    for(int iPadRow=0; iPadRow<45; iPadRow++)
      for(int iLocalXGrid=0; iLocalXGrid<nLocalXGrids+1; iLocalXGrid++)
	for(int iLocalZGrid=0; iLocalZGrid<nLocalZGrids+1; iLocalZGrid++)
	  for(int iComponent=0; iComponent<4; iComponent++)
	    ifstr>>tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid][iLocalZGrid][iComponent];
  */
  FILE* f1 = fopen(fileName, "r");
  fread(tpcHitCorrectionMap, sizeof(tpcHitCorrectionMap), 1, f1);
  fclose(f1);

}



void online_tracking_TCPhitcorrection::printMap()
{

}

void online_tracking_TCPhitcorrection::getCorrection(int sector, int padRow, const double* position, double spaceCharge, double* correctionXy)
{

      int iSector = sector-1;
      int iPadRow = padRow-1;
      double localXGrid = position[0]/localXGridSize[iPadRow]+(double)nLocalXGrids/2.;
      double localZGrid = position[2]/localZGridSize;
      if(localXGrid<0 || localXGrid>=(double)nLocalXGrids+1. || localZGrid<0)
	{
	  cout<<"online_tracking_TCPhitcorrection:WARN  - online_tracking_TCPhitcorrection::getCorrection(): out of grids. localXGrid:"<<localXGrid<<" localZGrid:"<<localZGrid<<endl;
	  correctionXy[0] = 0.;
	  correctionXy[1] = 0.;
	  return;
	}	  
      int iLocalXGrid = (int)localXGrid;
      int iLocalZGrid = (int)localZGrid;
      if(iLocalZGrid <= nLocalZGrids)
	{
	  double dLocalX = localXGrid-iLocalXGrid;
	  double dLocalZ = localZGrid-iLocalZGrid;
	  double correctionX00 = tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid][iLocalZGrid][0] + tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid][iLocalZGrid][2]*spaceCharge;
	  double correctionY00 = tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid][iLocalZGrid][1] + tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid][iLocalZGrid][3]*spaceCharge;
	  double correctionX01 = tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid][iLocalZGrid+1][0] + tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid][iLocalZGrid+1][2]*spaceCharge;
	  double correctionY01 = tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid][iLocalZGrid+1][1] + tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid][iLocalZGrid+1][3]*spaceCharge;
	  double correctionX10 = tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid+1][iLocalZGrid][0] + tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid+1][iLocalZGrid][2]*spaceCharge;
	  double correctionY10 = tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid+1][iLocalZGrid][1] + tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid+1][iLocalZGrid][3]*spaceCharge;
	  double correctionX11 = tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid+1][iLocalZGrid+1][0] + tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid+1][iLocalZGrid+1][2]*spaceCharge;
	  double correctionY11 = tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid+1][iLocalZGrid+1][1] + tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid+1][iLocalZGrid+1][3]*spaceCharge;
	  correctionXy[0] = correctionX00*(1-dLocalX)*(1-dLocalZ) + correctionX01*(1-dLocalX)*dLocalZ + correctionX10*dLocalX*(1-dLocalZ) + correctionX11*dLocalX*dLocalZ;
	  correctionXy[1] = correctionY00*(1-dLocalX)*(1-dLocalZ) + correctionY01*(1-dLocalX)*dLocalZ + correctionY10*dLocalX*(1-dLocalZ) + correctionY11*dLocalX*dLocalZ;
	}
      else     // if(iLocalZGrid > nLocalZGrids)
	{
	  iLocalZGrid = nLocalZGrids+1;
	  double dLocalX = localXGrid-iLocalXGrid;
	  double correctionX0 = tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid][iLocalZGrid][0] + tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid][iLocalZGrid][2]*spaceCharge;
	  double correctionY0 = tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid][iLocalZGrid][1] + tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid][iLocalZGrid][3]*spaceCharge;
	  double correctionX1 = tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid+1][iLocalZGrid][0] + tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid+1][iLocalZGrid][2]*spaceCharge;
	  double correctionY1 = tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid+1][iLocalZGrid][1] + tpcHitCorrectionMap[iSector][iPadRow][iLocalXGrid+1][iLocalZGrid][3]*spaceCharge;
	  correctionXy[0] = correctionX0*(1-dLocalX) + correctionX1*dLocalX;
	  correctionXy[1] = correctionY0*(1-dLocalX) + correctionY1*dLocalX;
	}


  

}
