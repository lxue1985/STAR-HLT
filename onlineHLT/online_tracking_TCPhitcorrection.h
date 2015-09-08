/***********************************************************************
 *
 * $Id: online_tracking_TCPhitcorrection.h,v 1.1 2009/03/06 18:32:09 tonko Exp $
 *
 * Author: QIU Hao   12/3/2008
 *
 ***********************************************************************
 *
 * Description: a map to correct tpc hit distortion
 *
 ***********************************************************************
 */
#ifndef online_tracking_TCPhitcorrection_H
#define online_tracking_TCPhitcorrection_H

const int nLocalXGrids = 30;
const int nLocalZGrids = 50;
const double minLocalZ = 0.;
const double maxLocalZ = 208.7;
const double localZGridSize = (maxLocalZ-minLocalZ)/nLocalZGrids;



class online_tracking_TCPhitcorrection 
{
  
 public:
  online_tracking_TCPhitcorrection();
  virtual ~online_tracking_TCPhitcorrection () {};


  void loadMap(const char* fileName = "tpcHitCorrectionMap.bin");
  void printMap();
  void getCorrection(int sector, int padRow, const double* position, double spaceCharge, double* correctionXy);

  double tpcHitCorrectionMap[24][45][nLocalXGrids+1][nLocalZGrids+1][4];
  //[sector-1][padrow-1][localSecotrXGrid][localSectorZGrid][componentIndex]
  //componentIndex = 0, 1 : independent on luminosity
  //componentIndex = 2, 3 : propotinal to luminosity
  //componentIndex = 0, 2 : in global x direction
  //componentIndex = 1, 3 : in global y direction

  double padRowLocalY[45], localXGridSize[45];


    
};
    
#endif
