#ifndef CLIB_KERNEL_H
#define CLIB_KERNEL_H


struct blocks{
int x;
int y;
};
int __mul24(int a,int b);


#include "FtfGeneral.h"
#include "FtfTrack.h"





/*__global__*/ void test(float* data);
/*__device__*/ float dtest(float *data);







#define min(a,b)    ( ( (a) < (b) ) ? (a) : (b) )
#define max(a,b)    ( ( (a) > (b) ) ? (a) : (b) )
#define seta(r,z)   (float)(3.0F * (z) / (fabs(z)+2.0F*(r)))
#define reta(eta,r) ((2.F*(r)*eta / ( 3 - fabs(eta)) )) 
#define sgn(a)      (float)( ( (a) > 0   ) ? (1) :(-1) )
#define square(a)   (float)( (a) * (a) )




    #define df_xyError                 0.12 
    #define df_zError                  0.24 
    #define df_hitChi2Cut          50.  
    #define df_goodHitChi2         20.  
    #define df_trackChi2Cut        10.  
    #define df_maxChi2Primary     50.    
    #define df_segmentRowSearchRange  2   
    #define df_trackRowSearchRange  3    
    #define df_dphi               0.08F  
    #define df_deta               0.08F  
    #define df_dphiMerge         0.01F  
    #define df_detaMerge         0.02F  
    #define df_etaMinTrack       -2.2F  
    #define df_etaMaxTrack        2.2F  
    #define df_etaMin       -2.5F  
    #define df_etaMax        2.5F  
    #define df_dEdx               1     
    #define df_getErrors          0     
    #define df_goBackwards        1     
    #define df_goodDistance       5.F   
    #define df_mergePrimaries     0     
    #define df_maxDistanceSegment  50.F 
    #define df_minHitsPerTrack   5      
    #define df_nHitsForSegment   2      
    #define df_nEta              40     
    #define df_nEtaTrack         40     
    #define df_nPhi              10     
    #define df_nPhiTrack         40     
    #define df_nPrimaryPasses    1      
    #define df_nSecondaryPasses  0      
    #define df_xyErrorScale      1.0F   
    #define df_szErrorScale      1.0F   
    #define df_phiClosed         0      
    #define df_ptMinHelixFit     0.1   
    #define df_rVertex           0.F    
    #define df_xVertex           0.F    
    #define df_yVertex           0.F    
    #define df_zVertex           0.F    
    #define df_dxVertex          0.05F 
    #define df_dyVertex          0.05F 
    #define df_phiVertex         0.F    
    #define df_zMax              205. 
    #define df_rowInnerMost      1 
    #define df_nRow              45 








/*__constant__*/   // float toDeg =57.29577951F ;
/*__constant__*/   // float pi    =3.141592654F ;
/*__constant__*/   // float twoPi =2.F*3.141592654F ; 
/*__constant__*/   // float piHalf = 0.5*3.141592654F ;
/*__constant__*/   // float bFactor = 0.0029979 ;

/*__constant__*/   // int  USE_SEGMENT= 1 ;
/*__constant__*/   // int  USE_FOLLOW = 2 ;
/*__constant__*/   // int  GO_DOWN    =-1 ;
/*__constant__*/   // int  GO_UP      = 1 ;





/*__global__*/ void tracking(int number,float *para,float *fhit,int *ihit,int * volumeC,int * rowC,float *track);

/*__global__*/ void class_hash(int *ihit,int number,int * volumeC,int * rowC);



/*__global__*/ void
clear(int * dataC,int number);

/*__global__*/ void
calc_hash(float *fhit,int *ihit,float *para,int number,int * volumeC,int * rowC);


/*__device__*/ class clib_FtfTrack;

/*__device__*/ class clib_FtfPara {          
 public:
	 void prt();
/*__device__*/  void      setDefaults ( ) ;

  int       infoLevel;       // Level of information printed about progress
  int       segmentRowSearchRange;       // Row search range for segments 
  int       trackRowSearchRange;         // Row search range for tracks 
  int       dEdx  ;          // dEdx switch
  int       dEdxNTruncate ;  // # points to truncate in dEdx
  int       minHitsForDedx;  // cs: min number of hits for dEdx calculation
  int       eventReset   ;   // Flag to reset event in fft 
  int       getErrors    ;   // Flag to switch error calculation
  int       fillTracks   ;   // Flag to switch FtfTrack class filling
  int       ghostFlag    ;   // =1 when there are ghost hits
  int       goBackwards  ;   // Flag to go backwards at the end of track reco
  int       init;            // Control initialization 
  int       mergePrimaries ; // Switch to control primary merging 
  int       minHitsPerTrack; // Minimum # hits per track 
  int       modRow;          // Modulo pad row number to use 
  int       nHitsForSegment; // # hits in initial segments 
  int       minHitsForFit;
  int       nEta;            // # volumes in eta 
  int       nEtaTrack;       // # Track areas in eta 
  int       nPhi;            // # volumes in nphi 
  int       nPhiTrack;       // # Track areas in nphi 
  int       nPrimaryPasses;  // # iterations looking for primaries
  int       nSecondaryPasses;// # iterations looking for secondaries
  int       vertexConstrainedFit; // 
  int       parameterLocation; // 1=inner most point, 0=closest approach
  float     maxChi2Primary ; // maximum chi2 to be considered primary 
  int       rowInnerMost;    // Row where end track search 
  int       rowOuterMost;    // Outer most row to consider tin tracking
  int       rowStart;        // Row where start track search
  int       rowEnd  ;        // Row where end   track search
  int       szFitFlag;       // Switch for sz fit 
  float     bField      ;    // Magnetic field  (magnitude = >0)
  int       bFieldPolarity;  // polarity of field (1 or -1)
  float     hitChi2Cut;      // Maximum hit chi2 
  float     goodHitChi2;     // Chi2 to stop looking for next hit 
  float     trackChi2Cut;    // Maximum track chi2 
  float     deta;            // Eta search range 
  float     dphi;            // Phi search range 
  float     detaMerge ;      // Eta difference for track merge 
  float     dphiMerge ;      // Phi difference for track merge
  float     distanceMerge ;  // Maximum distance for reference point to merge secondaries
  float     etaMin;          // Min eta to consider 
  float     etaMinTrack ;    // Track min eta to consider 
  float     etaMax;          // Max eta to consider 
  float     etaMaxTrack ;    // Track max eta to consider 
  float     goodDistance ;   // In segment building
  // distance consider good enough 
  float     phiMin;          // Min phi to consider 
  float     phiMinTrack ;    // Track min phi to consider 
  float     phiMax;          // Max phi to consider 
  float     phiMaxTrack ;    // Track max phi to consider 
  float     phiShift      ;  // Shift in phi when calculating phi
  float     ptMinHelixFit ;  // Minimum pt to apply helix fit
  float     maxDistanceSegment; // Maximum distance for segments 
  float     segmentMaxAngle; // Maximum angle between to consecutive track pieces 
  // when forming segments. A piece is the connection 
  // two hits
  float     szErrorScale;    // sz error scale 
  float     xyErrorScale;    // xy error scale 
  float     xVertex      ;   // x position primary vertex 
  float     yVertex      ;   // y position primary vertex 
  float     dxVertex     ;
  float     dyVertex     ;
  float     zVertex      ;
  float     xyWeightVertex;  // Weight vertex in x-y
  float     phiVertex      ;
  float     rVertex        ;
  float     maxTime        ; // maxTime tracker can run
  int       phiClosed ;
  int       primaries  ;
  int       nRowsPlusOne, nPhiPlusOne   ; // Number volumes + 1
  int       nEtaPlusOne, nPhiEtaPlusOne ; // Number volumes + 1 
  int       nPhiTrackPlusOne, nEtaTrackPlusOne ;                
  float     phiSlice, etaSlice ;
  float     phiSliceTrack, etaSliceTrack ;
  float     zMax ;
  int		current_ntrack;
  float     TC[1] ;

} ;



/*__device__*/ class clib_Ftf3DHit
{
 public:
    /*__device__*/ void set ( float _x, float _y, float _z )
	{ x = _x ; y = _y ; z = _z ; } ; 

    float x ;
    float y ;
    float z ;
};

/*__device__*/ class clib_FtfBaseHit
{         
   public:


      int          id ;
      short        row   ;         // Row #     
      void         *track;         // track to which the pnt was assigned
      void         *nextTrackHit  ;// Next track hit             
      void         *nextMcTrackHit;// Next MC track hit
      float        xyChi2 ;        // Chi2 in xy                 
      float        szChi2 ;        // Chi2 in sz                 
      float        x    ; 
      float        y    ;
      float        z    ;
      float        dx   ;          // error on the x coordinate  
      float        dy   ;          // error on the y coordinate  
      float        dz   ;          // error on the z coordinate  
      float        q    ;          // total charge assigned to this point 
      float        wxy  ;          // x-y weight x-y
      float        wz   ;          // z weight on z
      float        s    ;          // Track trajectory
} ;




/*__device__*/ class clib_FtfBaseTrack 
{ 
      
 public:
    void        *firstHit;// First hit belonging to track
    void        *lastHit ;// Last  hit belonging to track
    void        *currentHit ;
    /*__device__*/ int         fitHelix   (  ) ;
    /*__device__*/ int         refitHelix ( int mode, int modEqual, int rowMin, int rowMax ) ;
    /*__device__*/ int         fitCircle  (  ) ;
    /*__device__*/ int         fitLine    (  ) ;
    /*__device__*/ clib_FtfBaseHit* getCurrentHit ( ) { return (clib_FtfBaseHit *)currentHit ; } ;
    /*__device__*/ clib_FtfPara*    getPara()         { return (clib_FtfPara *)para ; } ;
    /*__device__*/ int         getErrorsCircleFit ( float a, float b, float r ) ;


    /*__device__*/ float   arcLength       ( float x1, float y1, float x2, float y2 ) ;
    /*__device__*/ clib_Ftf3DHit closestApproach ( float xBeam, float yBeam ) ;
    /*__device__*/ clib_Ftf3DHit extraRadius     ( float r ) ;
    /*__device__*/ int      extraRCyl       ( float &r, float &phi, float &z,
                               float &rc, float &xc, float &yc ) ;
    /*__device__*/ int      intersectorZLine    ( float a, float b, 
				   clib_Ftf3DHit& cross1, clib_Ftf3DHit& cross2 ) ;
    /*__device__*/ int      intersectorZLine    ( float a, float b, clib_Ftf3DHit& cross ) ;
    /*__device__*/ int      intersectorYCteLine ( float a, clib_Ftf3DHit& cross ) ;
    /*__device__*/ clib_Ftf3DHit getClosest      ( float xBeam, float yBeam,
	                       float &rc, float &xc, float &yc ) ;
    /*__device__*/ int      getClosest      ( float xBeam, float yBeam,
	                       float rc, float xc, float yc,
	                       float &xClosest, float &yClosest ) ;

    /*__device__*/ void     updateToRadius  ( float r ) ;
    /*__device__*/ void     updateToClosestApproach ( float xBeam, float yBeam, float rMax=10000. ) ;
    /*__device__*/ int      phiRotate       ( float deltaPhi ) ;
    // JB
    /*__device__*/ clib_Ftf3DHit extrapolate2PathLength(float pathlength);
    /*__device__*/ float getRadius();
    /*__device__*/ float getXCenter();
    /*__device__*/ float getYCenter();
    /*__device__*/ float pathLength(float Rx, float Ry, float Rz, float Nx, float Ny, float Nz );
    
    /*__device__*/ inline     void startLoop( ){ currentHit = firstHit ; } ;
	/*__device__*/    inline     void nextHit (){ currentHit = ((clib_FtfBaseHit *)currentHit)->nextTrackHit ; } ;
    /*__device__*/ inline     int  done     ( ) { return currentHit != 0 ; } ;
    /*__device__*/ void       Print       ( int level ) ;
	  

    int       id     ;  // primary key 
    short     flag   ;  // Primaries flag=1, Secondaries flag=0      
    char      innerMostRow ;
    char      outerMostRow ;
    short     nHits  ;  // Number of points assigned to that track
    short     nDedx  ;  // Number of points used for dEdx
    short     q  ;      // charge 
    float    chi2[2];  // chi squared of the momentum fit 
    float    dedx;     // dE/dx information 
    float    pt  ;     // pt (transverse momentum) at (r,phi,z) 
    float    phi0;     // azimuthal angle of point where parameters are given 
    float    psi ;     // azimuthal angle of the momentum at (r,.. 
    float    r0  ;     // r (in cyl. coord.) for point where parameters given 
    float    tanl;     // tg of the dip angle at (r,phi,z) 
    float    z0  ;     // z coordinate of point where parameters are given 
    float    length ;
    float    dpt ;
    float    dpsi;
    float    dz0 ;
    float    eta ;
    float    dtanl ;

    void     *para  ;    // Parameters pointer     
} ;

/*__device__*/ void DftfMatrixDiagonal ( float *h, float &h11, float &h22, float &h33 );

/*__device__*/ class clib_FtfContainer{
   public:
     void *first; 
     void *last;  
} ;



/*__device__*/ class clib_FtfHit: public clib_FtfBaseHit {         
 public:
/*__device__*/   void         setStatus ( clib_FtfTrack* this_track ) ;
/*__device__*/   void         setup ( int i, float*fhit ) ;
  
  long         id ;
  short        phiIndex ;        // Phi index    
  short        etaIndex ;        // Eta index    
  short        flags    ;        // various flags      
  short        sector   ;        // various flags      

  void*        nextVolumeHit ;  // Next volume hit            
  void*        nextRowHit    ;  // Next row hit               
  float        r    ;            // radius                     
  float        phi  ;            // azimuthal angle            
  float        dphi ;            // Error in phi               
  float        eta  ;            // hit pseudorapidity         
  float        xp   ;            // x conformal coordinate 
  float        yp   ;            // y conformal coordinate 
  short        buffer1 ;          //
  short        buffer2 ;         
  int          Mc_key;
  int          Mc_track_key;
  unsigned short hardwareId ;
} ;




/*__device__*/ class clib_FtfTrack : public clib_FtfBaseTrack {
	  
public:

/*__device__*/    void      add                   ( clib_FtfHit   *thisHit, int way ) ;
/*__device__*/    void      add                   ( clib_FtfTrack *thisTrack ) ;
/*__device__*/    int       buildTrack            ( clib_FtfHit *firstHit, clib_FtfContainer *volume ) ;
/*__device__*/    void      dEdx                  ( ) ;
/*__device__*/    void      deleteCandidate       ( ) ;
/*__device__*/    void      fill                  ( ) ;
/*__device__*/    void      fillPrimary           ( float &xc, float &yc, float &rc,
                                     float xPar, float yPar ) ;
/*__device__*/    void      fillSecondary         ( float &xc, float &yc, float xPar, float yPar ) ;
/*__device__*/    int       follow                ( clib_FtfContainer *volume, int way, int rowToStop ) ;
/*__device__*/    int       followHitSelection    ( clib_FtfHit *baseHit, clib_FtfHit *candidateHit, int way ) ;
/*__device__*/    clib_FtfTrack* getNextTrack ( )      { return (clib_FtfTrack *)nxatrk ; } ; 
/*__device__*/    int       mergePrimary          ( clib_FtfContainer   *trackArea ) ;
/*__device__*/    void      reset                 ( ) ;
/*__device__*/    clib_FtfHit    *seekNextHit          ( clib_FtfContainer  *volume, clib_FtfHit *baseHit,
	                             int nradiusSteps, int whichFunction ) ;
/*__device__*/    int     segment               ( clib_FtfContainer *volume, int way ) ;
/*__device__*/    int     segmentHitSelection ( clib_FtfHit *baseHit, clib_FtfHit *candidateHit ) ;


float *fhit;
int *ihit;
float *track;




   void    *nxatrk  ;      
        

		
   float   lastXyAngle ;    // Angle in the xy plane of line connecting to last hits        
   typedef float vfit ;
		
   vfit    xRefHit ;
   vfit    yRefHit ;
   vfit    xLastHit ;
   vfit    yLastHit ;

   vfit    s11Xy  ;       // Fit Parameters
   vfit    s12Xy  ;
   vfit    s22Xy  ;
   vfit    g1Xy   ;
   vfit    g2Xy   ;       
   vfit    s11Sz  ;
   vfit    s12Sz  ;
   vfit    s22Sz  ;
   vfit    g1Sz   ;
   vfit    g2Sz   ; 

   vfit    ddXy, a1Xy, a2Xy ;    /*fit par in xy */
   vfit    ddSz, a1Sz, a2Sz ;    /*fit par in sz */

//private:

	   
} ;

 class small_FtfHit {         
 public:

  long         id ;
         
  float        phi  ;            // azimuthal angle            
  float        eta  ;            // hit pseudorapidity         
  float        xp   ;            // x conformal coordinate 
  float        yp   ;            // y conformal coordinate 

      short        row   ;         // Row #     
      float        xyChi2 ;        // Chi2 in xy                 
      float        szChi2 ;        // Chi2 in sz                 
      float        x    ; 
      float        y    ;
      float        z    ;
      float        dx   ;          // error on the x coordinate  
      float        dy   ;          // error on the y coordinate  
      float        dz   ;          // error on the z coordinate  
      float        wxy  ;          // x-y weight x-y
      float        wz   ;          // z weight on z
      float        s    ;          // Track trajectory

	  


} ;


class small_FtfTrack{

public:
   float   lastXyAngle ;    // Angle in the xy plane of line connecting to last hits        
		


   //float    xRefHit ;
   //float    yRefHit ;
   //float    xLastHit ;
   //float    yLastHit ;

   float    s11Xy  ;       // Fit Parameters
   float    s12Xy  ;
   float    s22Xy  ;
   float    g1Xy   ;
   float    g2Xy   ;       
   float    s11Sz  ;
   float    s12Sz  ;
   float    s22Sz  ;
   float    g1Sz   ;
   float    g2Sz   ; 

   float    ddXy, a1Xy, a2Xy ;    /*fit par in xy */
   float    ddSz, a1Sz, a2Sz ;    /*fit par in sz */

    int       id     ;  // primary key 
    //short     flag   ;  // Primaries flag=1, Secondaries flag=0      
    //char      innerMostRow ;
    //char      outerMostRow ;
    short     nHits  ;  // Number of points assigned to that track
    //short     nDedx  ;  // Number of points used for dEdx
    //short     q  ;      // charge 
    float    chi2[2];  // chi squared of the momentum fit 
    //float    dedx;     // dE/dx information 
    //float    pt  ;     // pt (transverse momentum) at (r,phi,z) 
    //float    phi0;     // azimuthal angle of point where parameters are given 
    //float    psi ;     // azimuthal angle of the momentum at (r,.. 
    //float    r0  ;     // r (in cyl. coord.) for point where parameters given 
    //float    tanl;     // tg of the dip angle at (r,phi,z) 
    //float    z0  ;     // z coordinate of point where parameters are given 
    float    length ;
    //float    dpt ;
    //float    dpsi;
    //float    dz0 ;
    //float    eta ;
    //float    dtanl ;
	int hitmap[50];
	small_FtfHit lasthit;

};




 class small_hash{
 public:
	 int start;
	 int end;
	 int map;
 
 };

 //class small_track_assign{
 //public:

	// void ini();

	// int push(int id);

	// int number;
	// int track_id[20];
 //
 //};


 class small_track_manager{
 public:
	 void ini();
	 int pop();
	 int push(int id);
	 int number;
	 int maxused;
	 int track_id[2000];
 
 };


// clib_FtfPara * getParaV;
clib_FtfPara * getPara();

int small_segmentHitSelection ( small_FtfHit *baseHit, small_FtfHit *candidateHit,small_FtfTrack *tra );

int small_segmentHitGroup ( small_FtfHit *candidateHit,int num,small_FtfHit *selected_hit,int &selected_pos, small_FtfTrack *tra );


int small_followHitSelection ( small_FtfHit *baseHit, small_FtfHit *candidateHit,  int way, small_FtfTrack *tra );

int small_followHitGroup (  small_FtfHit *candidateHit,  int way, int num,small_FtfHit *selected_hit, int &selected_pos, small_FtfTrack *tra );


void small_add ( small_FtfHit *thisHit, int way, small_FtfTrack *tra  );
void small_reset ( small_FtfTrack *tra);


int element_setbit(int & data,int position);

int element_getbit(int data,int position);








#endif

