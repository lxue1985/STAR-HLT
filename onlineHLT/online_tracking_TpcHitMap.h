/***********************************************************************
 *
 * $Id: online_tracking_TpcHitMap.h,v 1.6 2011/01/13 20:38:47 qiuh Exp $
 *
 * Author: QIU Hao   12/3/2008
 * modified by Hongwei Ke 1/10/2011
 *
 ***********************************************************************
 *
 * Description: a map to correct tpc hit distortion
 *
 ***********************************************************************
 */
#ifndef online_tracking_TpcHitMap_H
#define online_tracking_TpcHitMap_H

const int nPadGrids = 30;
const int nTimeBucketGrids = 40;
const int nTimeBuckets = 400;
const double timeBucketGridSize = (double)nTimeBuckets/(double)nTimeBucketGrids;


class online_tracking_TpcHitMap 
{
  
public:
    online_tracking_TpcHitMap(char* HLTparameters="HLTparameters", int in_sectorID=0);
    virtual ~online_tracking_TpcHitMap () {};
    void readHLTparameters(char* HLTparameters="HLTparameters");
    void loadMap(const char* fileName = "tpcHitMap.bin");
    void mapping(double* xyz, int padRow, double pad, double timeBucket);
    void setSpaceChargeFromScalerCount(double scalerCount);  //use BBC or ZDC count
    void setSpaceChargeFromScalerCounts(int *scalerCounts);  //see below
    void setDriftVelocity(double inDriftVelocity);
    double tpcHitMap[45][nPadGrids+1][nTimeBucketGrids+1][3][2];
    //[padRow-1][padGrid][timeBucketGrid][x/y/z][independent or propotional to luminosity]

private:

    double padGridSize[45];
    double spaceChargeP0, spaceChargeP1;
    double spaceCharge;

    int scalerToUse;  // 0=RICH, 1=BBCx, 2=ZDCx, 3=ZDCsum, 4=BBCsum ;

    int    sectorID;            // 0 - 23, assigned by online_tracking_sector when construct
    double mapDriftVelocity;
    double driftVelocity;
};
    
#endif
/*
  scalerCounts[i]
  i    count
  0    bbce
  1    bbcw
  2    bbc and
  3    yellow beam background
  4    blue beam background
  5    zdce
  6    zdcw
  7    zdc and
  ...  the rest keep changing year by year
  http://online.star.bnl.gov/cgi-bin/db_scripts/cgi/database_scaler.pl
*/
