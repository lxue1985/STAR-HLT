#ifndef LIB_KERNEL_H
#define LIB_KERNEL_H


#define FETCH(t, i) tex1Dfetch(t##Tex, i)
#define HASHSIZE 20746
#define max_hit 8000


#define    getPara_infoLevel    0
#define    getPara_segmentRowSearchRange    2
#define    getPara_trackRowSearchRange    3
#define    getPara_dEdx      1
#define    getPara_dEdxNTruncate     20
#define    getPara_minHitsForDedx    15
#define    getPara_eventReset        1
#define    getPara_getErrors         0
#define    getPara_fillTracks        1
#define    getPara_ghostFlag        0
#define    getPara_goBackwards      1
#define    getPara_init    0
#define    getPara_mergePrimaries     0
#define    getPara_minHitsPerTrack    5
#define    getPara_modRow    1
#define    getPara_nHitsForSegment    2
#define    getPara_minHitsForFit    1
#define    getPara_nEta    40
#define    getPara_nEtaTrack    40
#define    getPara_nPhi    10
#define    getPara_nPhiTrack    40
#define    getPara_nPrimaryPasses    1
#define    getPara_nSecondaryPasses    0
#define    getPara_vertexConstrainedFit    0
#define    getPara_parameterLocation    1
#define    getPara_ maxChi2Primary     0
#define    getPara_rowInnerMost    1
#define    getPara_rowOuterMost    45
#define    getPara_rowStart    45
#define    getPara_rowEnd      1
#define    getPara_szFitFlag    1
#define    getPara_bField        0.501
#define    getPara_bFieldPolarity    -1
#define    getPara_hitChi2Cut    50
#define    getPara_goodHitChi2    20
#define    getPara_trackChi2Cut    10
#define    getPara_deta    0.08
#define    getPara_dphi    0.08
#define    getPara_detaMerge     0.02
#define    getPara_dphiMerge     0.01
#define    getPara_distanceMerge     2
#define    getPara_etaMin    -0.4
#define    getPara_etaMinTrack    -2.2
#define    getPara_etaMax    2.3
#define    getPara_etaMaxTrack     2.2
#define    getPara_goodDistance     5
#define    getPara_phiMin    0.785398
#define    getPara_phiMinTrack    -1.74533e-08
#define    getPara_phiMax    1.309
#define    getPara_phiMaxTrack     6.28668
#define    getPara_phiShift       0
#define    getPara_ptMinHelixFit     0
#define    getPara_maxDistanceSegment    50
#define    getPara_segmentMaxAngle    0.174533
#define    getPara_szErrorScale    1
#define    getPara_xyErrorScale    1
#define    getPara_xVertex       0
#define    getPara_yVertex       0
#define    getPara_dxVertex      0.05
#define    getPara_dyVertex      0.05
#define    getPara_zVertex       0
#define    getPara_xyWeightVertex    -1.77759e-06
#define    getPara_phiVertex       0
#define    getPara_rVertex         0
#define    getPara_maxTime         1e+18
#define    getPara_phiClosed     0
#define    getPara_primaries      1
#define    getPara_nRowsPlusOne    -1073975984
#define    getPara_nPhiPlusOne    -1208504332
#define    getPara_nEtaPlusOne     48
#define    getPara_nPhiEtaPlusOne     -1245493376
#define    getPara_nPhiTrackPlusOne     -1244130928
#define    getPara_nEtaTrackPlusOne     -1208578100
#define    getPara_phiSlice     0
#define    getPara_etaSlice     -1.64829e-06
#define    getPara_phiSliceTrack     -1.60932e-06
#define    getPara_etaSliceTrack     0
#define    getPara_zMax     205
#define    getPara_current_ntrack    -1073975976




#define toDeg 57.29577951F 
#define pi    3.141592654F 
#define twoPi 2.F*pi  
#define GO_DOWN -1 
#define N_LOOP 9 
//
//-->   Functions
//
#define min(a,b)    ( ( (a) < (b) ) ? (a) : (b) )
#define max(a,b)    ( ( (a) > (b) ) ? (a) : (b) )
#define seta(r,z)   (float)(3.0F * (z) / (fabs(z)+2.0F*(r)))
#define reta(eta,r) ((2.F*(r)*eta / ( 3 - fabs(eta)) )) 
#define sgn(a)      (float)( ( (a) > 0   ) ? (1) :(-1) )
#define square(a)   (float)( (a) * (a) )




texture<float4, 1, cudaReadModeElementType> positionTex;
texture<int2, 1, cudaReadModeElementType> hashTex;

__global__ void test(float4* position, int2*hash,int*map,int* manager, float *out,int* number);
__device__ float dtest(float *data);




__device__ class small_FtfHit {         
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


__device__ class small_FtfTrack{

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




__device__ class small_hash{
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


__device__ struct small_track_manager{
 public:


	 int number;
	 int maxused;
	 int track_id[2000];
 
 };




__device__ int small_segmentHitSelection ( small_FtfHit *baseHit, small_FtfHit *candidateHit,small_FtfTrack *tra );

__device__ int small_segmentHitGroup ( small_FtfHit *candidateHit,int num,small_FtfHit *selected_hit,int &selected_pos, small_FtfTrack *tra );


__device__ int small_followHitSelection ( small_FtfHit *baseHit, small_FtfHit *candidateHit,  int way, small_FtfTrack *tra );

__device__ int small_followHitGroup (  small_FtfHit *candidateHit,  int way, int num,small_FtfHit *selected_hit, int &selected_pos, small_FtfTrack *tra );


__device__ void small_add ( small_FtfHit *thisHit, int way, small_FtfTrack *tra  );
__device__ void small_reset ( small_FtfTrack *tra);


__device__ int element_setbit(int & data,int position);

__device__ int element_getbit(int data,int position);







#endif

