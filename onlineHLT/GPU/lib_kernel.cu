#include <cstdlib>
#include <cstdio>
#include <string.h>


#include <cuda_gl_interop.h>
#include "lib_kernel.cuh"


__global__ void test(float4* position, int2*hash,int*map,int* manager, float *out, int *number){


#define DIM 24
#define MAXTRACK 500

 int index;

 index = blockIdx.x*blockDim.x + threadIdx.x;
 

__shared__ int sindex[DIM];
__shared__ int sblockIdx[DIM];

sindex[threadIdx.x]=index;
sblockIdx[threadIdx.x]=manager[blockIdx.x+30];

if(sindex[threadIdx.x]>=number[sblockIdx[threadIdx.x]]) return;

small_FtfTrack* out_track=  (small_FtfTrack*)(out+sblockIdx[threadIdx.x]*sizeof(small_FtfTrack)*MAXTRACK);

// __shared__   int  track_index[DIM];


__shared__  small_FtfTrack localtrack[DIM];
__shared__  small_FtfTrack* track[DIM];


__shared__   int  hashIndex[DIM];


__shared__   int  candidate[DIM];
__shared__   int  cnum[DIM];


	__shared__   int  cei[DIM];
	__shared__   int  cpi[DIM];
	__shared__   int  cri[DIM];

   __shared__   int  loop_eta[9];
   __shared__   int  loop_phi[9];
__shared__   int  rhoi[DIM];
__shared__   int  phii[DIM];
__shared__   int  etai[DIM];
__shared__   int  miss[DIM];
__shared__   int  ccri[DIM];

__shared__   int  result[DIM];	

__shared__   int  curmap[DIM];
__shared__   int  curpos[DIM];	
__shared__   small_FtfHit selected_hit[DIM];
__shared__   int  map_index[DIM];
__shared__   int  pos_index[DIM];
__shared__   int  ccei[DIM];
__shared__   int  ccpi[DIM];
__shared__   int  cchashIndex[DIM];
__shared__   int  hit_num[DIM];
__shared__   int  start[DIM];
__shared__   int  hl[DIM];

__shared__   float  dx[DIM];
__shared__   float  dy[DIM];

 __shared__   float  xyChi2[DIM];
 __shared__   float  szChi2[DIM];
__shared__   float  lxyChi2[DIM];
// __shared__   float  lszChi2[DIM];
// __shared__   float  slocal[DIM];


		

__shared__   float  r2[DIM];

__shared__   float  invR2[DIM];


   
   loop_eta[0]=0;
   loop_eta[1]=0;
   loop_eta[2]=0;
   loop_eta[3]=-1;
   loop_eta[4]=-1;
   loop_eta[5]=-1;
   loop_eta[6]=1;
   loop_eta[7]=1;
   loop_eta[8]=1;
   

	loop_phi[0]=0;   
	loop_phi[1]=-1;   
	loop_phi[2]=1;   
	loop_phi[3]=0;   
	loop_phi[4]=-1;   
	loop_phi[5]=1;   
	loop_phi[6]=0;   
	loop_phi[7]=-1;   
	loop_phi[8]=1;   
   
   
   phii[threadIdx.x]=sindex[threadIdx.x]/41;
   etai[threadIdx.x]=sindex[threadIdx.x]%41;
   
   
   do{
   
   
   cnum[threadIdx.x]=atomicAdd(&map[sblockIdx[threadIdx.x]],1);
   
   
   if(cnum[threadIdx.x]>number[sblockIdx[threadIdx.x]]) break;
   if(position[(sblockIdx[threadIdx.x])*max_hit+cnum[threadIdx.x]].w<0) continue;


   candidate[threadIdx.x]=cnum[threadIdx.x];

   
   cnum[threadIdx.x]=position[(sblockIdx[threadIdx.x])*max_hit+cnum[threadIdx.x]].w;
   
   rhoi[threadIdx.x]=cnum[threadIdx.x]/451;
   
   if(rhoi[threadIdx.x]<7) break;
   
   cnum[threadIdx.x]=cnum[threadIdx.x]%451;

   phii[threadIdx.x]=cnum[threadIdx.x]/41;
   etai[threadIdx.x]=cnum[threadIdx.x]%41;   
 
				// for( rhoi[threadIdx.x]=45;rhoi[threadIdx.x]>0;rhoi[threadIdx.x]--){
				// for( phii[threadIdx.x]=10;phii[threadIdx.x]>0;phii[threadIdx.x]--){	
				// for( etai[threadIdx.x]=40;etai[threadIdx.x]>0;etai[threadIdx.x]--){



 	     hashIndex[threadIdx.x] =  (rhoi[threadIdx.x])  * (41*11) + 
        (phii[threadIdx.x]) * 41 + (etai[threadIdx.x]) ;





 //gmap[threadIdx.x]=map[hashIndex[threadIdx.x]];
 //ghit_num[threadIdx.x]=hash[(sblockIdx[threadIdx.x])*HASHSIZE+hashIndex[threadIdx.x]].y-hash[(sblockIdx[threadIdx.x])*HASHSIZE+hashIndex[threadIdx.x]].x;



// for( candidate[threadIdx.x]=gstart[threadIdx.x]; candidate[threadIdx.x]<gstart[threadIdx.x]+ghit_num[threadIdx.x];candidate[threadIdx.x]++){


//cout<<"candidate[threadIdx.x]   "<<candidate[threadIdx.x]<<" hashIndex[threadIdx.x] "<<hashIndex[threadIdx.x]<<" hit_hash[hashIndex[threadIdx.x]].end "<<hit_hash[hashIndex[threadIdx.x]].end<<endl;







	



						 
		// track_index[threadIdx.x]=manager[cnum[threadIdx.x]+1];

	track[threadIdx.x]=&localtrack[threadIdx.x];


	

	small_reset(track[threadIdx.x]);





					selected_hit[threadIdx.x].x=position[(sblockIdx[threadIdx.x])*max_hit+candidate[threadIdx.x]].x;
					selected_hit[threadIdx.x].y=position[(sblockIdx[threadIdx.x])*max_hit+candidate[threadIdx.x]].y;
					selected_hit[threadIdx.x].z=position[(sblockIdx[threadIdx.x])*max_hit+candidate[threadIdx.x]].z;


					
					

					  r2[threadIdx.x]           =square(selected_hit[threadIdx.x].x)+square(selected_hit[threadIdx.x].y) ;


					 
				   // r[threadIdx.x]             = (float)sqrt ( r2[threadIdx.x] ) ;
				  selected_hit[threadIdx.x].phi           = (float)atan2(selected_hit[threadIdx.x].y,selected_hit[threadIdx.x].x) + getPara_phiShift ;
				  if (selected_hit[threadIdx.x].phi < 0 ) selected_hit[threadIdx.x].phi = selected_hit[threadIdx.x].phi + twoPi ;
				  // l3Log("r: %f, z: %f\n",r, selected_hit[threadIdx.x].z);
				  selected_hit[threadIdx.x].eta           = (float)seta((float)sqrt ( r2[threadIdx.x] ),selected_hit[threadIdx.x].z) ;

				  if ( getPara_szFitFlag ) {
					selected_hit[threadIdx.x].s  = 0.F ;
					selected_hit[threadIdx.x].wz = (float)(1./ square ( getPara_szErrorScale *0.24 ));
				  }


				  
				  
				  
				  

					 invR2[threadIdx.x]        = 1.F / r2[threadIdx.x] ;
					selected_hit[threadIdx.x].xp    =     selected_hit[threadIdx.x].x * invR2[threadIdx.x] ;
					selected_hit[threadIdx.x].yp    =   - selected_hit[threadIdx.x].y * invR2[threadIdx.x] ;
					selected_hit[threadIdx.x].wxy   =   r2[threadIdx.x] * r2[threadIdx.x] /  ( square(getPara_xyErrorScale)
					* ( square(0.12) + square(0.12) ) ) ;
					selected_hit[threadIdx.x].id= candidate[threadIdx.x];
					selected_hit[threadIdx.x].row=rhoi[threadIdx.x];








	small_add(&selected_hit[threadIdx.x],GO_DOWN,track[threadIdx.x]);


//position[(sblockIdx[threadIdx.x])*max_hit+candidate[threadIdx.x]].w=-1;


candidate[threadIdx.x]=0;





//		cout<<" smhit[candidate[threadIdx.x]] "<<smhit[candidate[threadIdx.x]].id<<endl;


	 cei[threadIdx.x]=etai[threadIdx.x];
	 cpi[threadIdx.x]=phii[threadIdx.x];
	 cri[threadIdx.x]=rhoi[threadIdx.x];
//cout<<" bug1 "<<hashIndex[threadIdx.x]<<"  "<<etai[threadIdx.x]<<"  "<<phii[threadIdx.x]<<"  "<<rhoi[threadIdx.x]<<"  "<<track[threadIdx.x]->nHits<<"  "<<candidate[threadIdx.x]<<endl;


 miss[threadIdx.x]=0;
	for( ccri[threadIdx.x]=cri[threadIdx.x]-1;ccri[threadIdx.x]>0;ccri[threadIdx.x]--){

	








	if(track[threadIdx.x]->nHits >= getPara_nHitsForSegment)break;
		track[threadIdx.x]->chi2[0] = getPara_maxDistanceSegment ;
		
		

		 result[threadIdx.x]=0;



		for( map_index[threadIdx.x]=0;map_index[threadIdx.x]<9;map_index[threadIdx.x]++){
			
				
				


				 ccei[threadIdx.x]=cei[threadIdx.x]+loop_eta[map_index[threadIdx.x]];
				 ccpi[threadIdx.x]=cpi[threadIdx.x]+loop_phi[map_index[threadIdx.x]];

				if(ccei[threadIdx.x]>40||ccei[threadIdx.x]<1) continue;
				if(ccpi[threadIdx.x]>10||ccpi[threadIdx.x]<1) continue;

 				 cchashIndex[threadIdx.x] =  (ccri[threadIdx.x])  * (41*11) + 
				(ccpi[threadIdx.x]) * 41 + (ccei[threadIdx.x]) ;
//		 cout<<" ccri[threadIdx.x] "<<ccri[threadIdx.x]<<" ccpi[threadIdx.x] "<<ccpi[threadIdx.x]<<"ccei[threadIdx.x] "<<ccei[threadIdx.x]<<endl;

				//maplist[threadIdx.x][map_index[threadIdx.x]]=map[cchashIndex[threadIdx.x]];

				
				
				 hit_num[threadIdx.x]=hash[(sblockIdx[threadIdx.x])*HASHSIZE+cchashIndex[threadIdx.x]].y-hash[(sblockIdx[threadIdx.x])*HASHSIZE+cchashIndex[threadIdx.x]].x;
				 start[threadIdx.x]=hash[(sblockIdx[threadIdx.x])*HASHSIZE+cchashIndex[threadIdx.x]].x;
				
				

				for( hl[threadIdx.x]=0;hl[threadIdx.x]<hit_num[threadIdx.x];hl[threadIdx.x]++){





					selected_hit[threadIdx.x].x=position[(sblockIdx[threadIdx.x])*max_hit+start[threadIdx.x]+hl[threadIdx.x]].x;
					selected_hit[threadIdx.x].y=position[(sblockIdx[threadIdx.x])*max_hit+start[threadIdx.x]+hl[threadIdx.x]].y;
					selected_hit[threadIdx.x].z=position[(sblockIdx[threadIdx.x])*max_hit+start[threadIdx.x]+hl[threadIdx.x]].z;



					  r2[threadIdx.x]           =square(selected_hit[threadIdx.x].x)+square(selected_hit[threadIdx.x].y) ;


					 
				   // r[threadIdx.x]             = (float)sqrt ( r2[threadIdx.x] ) ;
				  selected_hit[threadIdx.x].phi           = (float)atan2(selected_hit[threadIdx.x].y,selected_hit[threadIdx.x].x) + getPara_phiShift ;
				  if (selected_hit[threadIdx.x].phi < 0 ) selected_hit[threadIdx.x].phi = selected_hit[threadIdx.x].phi + twoPi ;
				  // l3Log("r: %f, z: %f\n",r, selected_hit[threadIdx.x].z);
				  selected_hit[threadIdx.x].eta           = (float)seta((float)sqrt ( r2[threadIdx.x] ),selected_hit[threadIdx.x].z) ;

				  if ( getPara_szFitFlag ) {
					selected_hit[threadIdx.x].s  = 0.F ;
					selected_hit[threadIdx.x].wz = (float)(1./ square ( getPara_szErrorScale *0.24 ));
				  }

				  
				  
				  
				  

					 invR2[threadIdx.x]        = 1.F / r2[threadIdx.x] ;
					selected_hit[threadIdx.x].xp    =     selected_hit[threadIdx.x].x * invR2[threadIdx.x] ;
					selected_hit[threadIdx.x].yp    =   - selected_hit[threadIdx.x].y * invR2[threadIdx.x] ;
					selected_hit[threadIdx.x].wxy   =   r2[threadIdx.x] * r2[threadIdx.x] /  ( square(getPara_xyErrorScale)
					* ( square(0.12) + square(0.12) ) ) ;
					selected_hit[threadIdx.x].id=start[threadIdx.x]+hl[threadIdx.x];
					selected_hit[threadIdx.x].row=ccri[threadIdx.x];


				
					 
				
							result[threadIdx.x]= small_segmentHitSelection(&track[threadIdx.x]->lasthit,&selected_hit[threadIdx.x],track[threadIdx.x]);

				
							
							
							
								 if ( result[threadIdx.x] > 0 ) {
										pos_index[threadIdx.x]=start[threadIdx.x]+hl[threadIdx.x];

							if ( result[threadIdx.x] ==2  ) break ; 
						 }				
				}



	
	


             if ( result[threadIdx.x] > 0 ) {
				curmap[threadIdx.x]=map_index[threadIdx.x];
				curpos[threadIdx.x]=pos_index[threadIdx.x];
                if ( result[threadIdx.x] ==2  ) break; ; 
             }


		}
//cout<<" result[threadIdx.x] "<<result[threadIdx.x]<<endl;
			

		if(result[threadIdx.x]>0){
			miss[threadIdx.x]=0;
			//element_setbit(maplist[threadIdx.x][curmap[threadIdx.x]],curpos[threadIdx.x]);
			if ( getPara_szFitFlag  ){
				 dx[threadIdx.x] = selected_hit[threadIdx.x].x - track[threadIdx.x]->lasthit.x ;
				 dy[threadIdx.x] = selected_hit[threadIdx.x].y - track[threadIdx.x]->lasthit.y ;
            track[threadIdx.x]->length    += (float)sqrt ( dx[threadIdx.x] * dx[threadIdx.x] + dy[threadIdx.x] * dy[threadIdx.x] ) ;
            selected_hit[threadIdx.x].s      = track[threadIdx.x]->length ;
         }
			small_add(&selected_hit[threadIdx.x],GO_DOWN,track[threadIdx.x]);

			// if(position[(sblockIdx[threadIdx.x])*max_hit+curpos[threadIdx.x]].w<0){candidate[threadIdx.x]++;}else{
			// candidate[threadIdx.x]=0;
			// }
			// if(candidate[threadIdx.x]>2) break;

			
			cei[threadIdx.x]=cei[threadIdx.x]+loop_eta[curmap[threadIdx.x]];
			cpi[threadIdx.x]=cpi[threadIdx.x]+loop_phi[curmap[threadIdx.x]];
			cri[threadIdx.x]=ccri[threadIdx.x];
 				 cchashIndex[threadIdx.x] =  (ccri[threadIdx.x])  * (41*11) + 
				(cpi[threadIdx.x]) * 41 + (cei[threadIdx.x]) ;
				//map[cchashIndex[threadIdx.x]]=maplist[threadIdx.x][curmap[threadIdx.x]];
				position[(sblockIdx[threadIdx.x])*max_hit+curpos[threadIdx.x]].w=-1;

		}else{
			miss[threadIdx.x]++;
			if(miss[threadIdx.x]>=2)	break;		
		}
		
		
	
	
	}




	if(track[threadIdx.x]->nHits < getPara_nHitsForSegment||candidate[threadIdx.x]>2){ 
	// cnum[threadIdx.x]=atomicSub(&manager[0],1);
	// manager[cnum[threadIdx.x]]=track_index[threadIdx.x];
	continue;}

//cout<<" track[threadIdx.x]->nHits "<<track[threadIdx.x]->nHits<<endl;


    xyChi2[threadIdx.x] = track[threadIdx.x]->chi2[0] ;
    szChi2[threadIdx.x] = track[threadIdx.x]->chi2[1] ;


	for( ccri[threadIdx.x]=cri[threadIdx.x]-1;ccri[threadIdx.x]>getPara_rowInnerMost;ccri[threadIdx.x]--){

		
	
		track[threadIdx.x]->chi2[0] = getPara_hitChi2Cut ;
		
		 result[threadIdx.x]=0;

//cout<<" bug2 "<<track_index[threadIdx.x]<<endl;

		for( map_index[threadIdx.x]=0;map_index[threadIdx.x]<9;map_index[threadIdx.x]++){
			
				
				


				 ccei[threadIdx.x]=cei[threadIdx.x]+loop_eta[map_index[threadIdx.x]];
				 ccpi[threadIdx.x]=cpi[threadIdx.x]+loop_phi[map_index[threadIdx.x]];

				if(ccei[threadIdx.x]>40||ccei[threadIdx.x]<1) continue;
				if(ccpi[threadIdx.x]>10||ccpi[threadIdx.x]<1) continue;
//		 cout<<"follow  ccri[threadIdx.x] "<<ccri[threadIdx.x]<<" ccpi[threadIdx.x] "<<ccpi[threadIdx.x]<<" ccei[threadIdx.x] "<<ccei[threadIdx.x]<<endl;



 				 cchashIndex[threadIdx.x] =  (ccri[threadIdx.x])  * (41*11) + 
				(ccpi[threadIdx.x]) * 41 + (ccei[threadIdx.x]) ;

				//maplist[threadIdx.x][map_index[threadIdx.x]]=map[cchashIndex[threadIdx.x]];

				 hit_num[threadIdx.x]=hash[(sblockIdx[threadIdx.x])*HASHSIZE+cchashIndex[threadIdx.x]].y-hash[(sblockIdx[threadIdx.x])*HASHSIZE+cchashIndex[threadIdx.x]].x;
				 start[threadIdx.x]=hash[(sblockIdx[threadIdx.x])*HASHSIZE+cchashIndex[threadIdx.x]].x;



				

				for( hl[threadIdx.x]=0;hl[threadIdx.x]<hit_num[threadIdx.x];hl[threadIdx.x]++){








					selected_hit[threadIdx.x].x=position[(sblockIdx[threadIdx.x])*max_hit+start[threadIdx.x]+hl[threadIdx.x]].x;
					selected_hit[threadIdx.x].y=position[(sblockIdx[threadIdx.x])*max_hit+start[threadIdx.x]+hl[threadIdx.x]].y;
					selected_hit[threadIdx.x].z=position[(sblockIdx[threadIdx.x])*max_hit+start[threadIdx.x]+hl[threadIdx.x]].z;



					 // x[threadIdx.x]            = selected_hit[threadIdx.x].x - getPara_xVertex ;
					 // y[threadIdx.x]            = selected_hit[threadIdx.x].y - getPara_yVertex ;
					 // r2[threadIdx.x]           = x[threadIdx.x] * x[threadIdx.x] + y[threadIdx.x] * y[threadIdx.x] ;
					  r2[threadIdx.x]           =square(selected_hit[threadIdx.x].x)+square(selected_hit[threadIdx.x].y) ;


					 
				   // r[threadIdx.x]             = (float)sqrt ( r2[threadIdx.x] ) ;
				  selected_hit[threadIdx.x].phi           = (float)atan2(selected_hit[threadIdx.x].y,selected_hit[threadIdx.x].x) + getPara_phiShift ;
				  if (selected_hit[threadIdx.x].phi < 0 ) selected_hit[threadIdx.x].phi = selected_hit[threadIdx.x].phi + twoPi ;
				  // l3Log("r: %f, z: %f\n",r, selected_hit[threadIdx.x].z);
				  selected_hit[threadIdx.x].eta           = (float)seta((float)sqrt ( r2[threadIdx.x] ),selected_hit[threadIdx.x].z) ;

				  if ( getPara_szFitFlag ) {
					selected_hit[threadIdx.x].s  = 0.F ;
					selected_hit[threadIdx.x].wz = (float)(1./ square ( getPara_szErrorScale *0.24 ));
				  }


				  
				  
				  
				  

					 invR2[threadIdx.x]        = 1.F / r2[threadIdx.x] ;
					selected_hit[threadIdx.x].xp    =     selected_hit[threadIdx.x].x * invR2[threadIdx.x] ;
					selected_hit[threadIdx.x].yp    =   - selected_hit[threadIdx.x].y * invR2[threadIdx.x] ;
					selected_hit[threadIdx.x].wxy   =   r2[threadIdx.x] * r2[threadIdx.x] /  ( square(getPara_xyErrorScale)
					* ( square(0.12) + square(0.12) ) ) ;
					selected_hit[threadIdx.x].id=start[threadIdx.x]+hl[threadIdx.x];
					selected_hit[threadIdx.x].row=ccri[threadIdx.x];





result[threadIdx.x]= small_followHitSelection(&track[threadIdx.x]->lasthit,&selected_hit[threadIdx.x],GO_DOWN,track[threadIdx.x]);
		         





				 if ( result[threadIdx.x] > 0 ) {
							pos_index[threadIdx.x]=start[threadIdx.x]+hl[threadIdx.x];;


                if ( result[threadIdx.x] ==2  ) break ; 
             }
	
					
				
				}







				
				
				
				
				
				
				
             if ( result[threadIdx.x] > 0 ) {
				curmap[threadIdx.x]=map_index[threadIdx.x];
				curpos[threadIdx.x]=pos_index[threadIdx.x];
                if ( result[threadIdx.x] ==2  ) break; ; 
             }


		}
//cout<<" bug3 "<<track_index[threadIdx.x]<<endl;
		if(result[threadIdx.x]>0){
			miss[threadIdx.x]=0;
			//element_setbit(maplist[threadIdx.x][curmap[threadIdx.x]],curpos[threadIdx.x]);
			



       lxyChi2[threadIdx.x] = track[threadIdx.x]->chi2[0]-track[threadIdx.x]->chi2[1] ;
      xyChi2[threadIdx.x] += lxyChi2[threadIdx.x] ;
      selected_hit[threadIdx.x].xyChi2 = lxyChi2[threadIdx.x] ;
//
//   if sz fit update track[threadIdx.x] length
//
      if ( getPara_szFitFlag  ) {
         track[threadIdx.x]->length = selected_hit[threadIdx.x].s ;
         szChi2[threadIdx.x] += track[threadIdx.x]->chi2[1]  ;
         selected_hit[threadIdx.x].szChi2 = track[threadIdx.x]->chi2[1] ;

      }
//
//     Add hit to track[threadIdx.x]
//
			small_add(&selected_hit[threadIdx.x],GO_DOWN,track[threadIdx.x]);

			
			// if(position[(sblockIdx[threadIdx.x])*max_hit+curpos[threadIdx.x]].w<0){candidate[threadIdx.x]++;}else{
			// candidate[threadIdx.x]=0;
			// }
			// if(candidate[threadIdx.x]>2) break;
			
			
//cout<<" bug4 "<<track[threadIdx.x]->nHits<<endl;
	
			
			// if(track[threadIdx.x]->nHits == getPara_minHitsPerTrack){
			
		// for( map_index[threadIdx.x]=0;map_index[threadIdx.x]<getPara_minHitsPerTrack;map_index[threadIdx.x]++){
// localtrack[threadIdx.x].hitmap[map_index[threadIdx.x]];	
// position[(sblockIdx[threadIdx.x])*max_hit+localtrack[threadIdx.x].hitmap[map_index[threadIdx.x]]].w=-1;		
	// }
	// }else if(track[threadIdx.x]->nHits > getPara_minHitsPerTrack){
	
// position[(sblockIdx[threadIdx.x])*max_hit+curpos[threadIdx.x]].w=-1;	
	// }			
			
			// position[(sblockIdx[threadIdx.x])*max_hit+curpos[threadIdx.x]].w=-1;	
			


		
			
			
			
			cei[threadIdx.x]=cei[threadIdx.x]+loop_eta[curmap[threadIdx.x]];
			cpi[threadIdx.x]=cpi[threadIdx.x]+loop_phi[curmap[threadIdx.x]];
			cri[threadIdx.x]=ccri[threadIdx.x];
 				 cchashIndex[threadIdx.x] =  (ccri[threadIdx.x])  * (41*11) + 
				(cpi[threadIdx.x]) * 41 + (cei[threadIdx.x]) ;
				// map[cchashIndex[threadIdx.x]]=maplist[threadIdx.x][curmap[threadIdx.x]];
				position[(sblockIdx[threadIdx.x])*max_hit+curpos[threadIdx.x]].w=-1;

		}else{
			miss[threadIdx.x]++;
			if(miss[threadIdx.x]>=3)	break;		
		}
		
		
	
	
	}


   if ( track[threadIdx.x]->nHits < getPara_minHitsPerTrack||candidate[threadIdx.x]>2 ) { 
   // cnum[threadIdx.x]=atomicSub(&manager[0],1);
   // manager[cnum[threadIdx.x]]=track_index[threadIdx.x];
   continue ;} 
//
//   Store track[threadIdx.x] chi2
//

	
   track[threadIdx.x]->chi2[0] = xyChi2[threadIdx.x] ;
   track[threadIdx.x]->chi2[1] = szChi2[threadIdx.x] ;
//
//        Check total chi2
//
   float normalized_chi2 = (track[threadIdx.x]->chi2[0]+track[threadIdx.x]->chi2[1])/track[threadIdx.x]->nHits ;
   if ( normalized_chi2 > getPara_trackChi2Cut ) {
   // cnum[threadIdx.x]=atomicSub(&manager[0],1);
   // manager[cnum[threadIdx.x]]=track_index[threadIdx.x];
   continue ;}	

   
cnum[threadIdx.x]=atomicAdd(&manager[sblockIdx[threadIdx.x]],1);   
   
   
// out[track_index[threadIdx.x]]=localtrack[threadIdx.x].nHits;	


 out_track[cnum[threadIdx.x]].nHits=localtrack[threadIdx.x].nHits;
 out_track[cnum[threadIdx.x]].lastXyAngle=localtrack[threadIdx.x].lastXyAngle;	
 out_track[cnum[threadIdx.x]].s11Xy=localtrack[threadIdx.x].s11Xy;	
 out_track[cnum[threadIdx.x]].s12Xy=localtrack[threadIdx.x].s12Xy;	
 out_track[cnum[threadIdx.x]].s22Xy=localtrack[threadIdx.x].s22Xy;	
 out_track[cnum[threadIdx.x]].g1Xy=localtrack[threadIdx.x].g1Xy;	
 out_track[cnum[threadIdx.x]].g2Xy=localtrack[threadIdx.x].g2Xy;	
 out_track[cnum[threadIdx.x]].s11Sz=localtrack[threadIdx.x].s11Sz;	
 out_track[cnum[threadIdx.x]].s12Sz=localtrack[threadIdx.x].s12Sz;	
 out_track[cnum[threadIdx.x]].s22Sz=localtrack[threadIdx.x].s22Sz;	
 out_track[cnum[threadIdx.x]].g1Sz=localtrack[threadIdx.x].g1Sz;	
 out_track[cnum[threadIdx.x]].g2Sz=localtrack[threadIdx.x].g2Sz;	
 out_track[cnum[threadIdx.x]].ddXy=localtrack[threadIdx.x].ddXy;	
 out_track[cnum[threadIdx.x]].a1Xy=localtrack[threadIdx.x].a1Xy;	
 out_track[cnum[threadIdx.x]].a2Xy=localtrack[threadIdx.x].a2Xy;	
 out_track[cnum[threadIdx.x]].a2Sz=localtrack[threadIdx.x].a2Sz;	
 out_track[cnum[threadIdx.x]].a1Sz=localtrack[threadIdx.x].a1Sz;	
 out_track[cnum[threadIdx.x]].chi2[0]=localtrack[threadIdx.x].chi2[0];	
 out_track[cnum[threadIdx.x]].chi2[1]=localtrack[threadIdx.x].chi2[1];	
 out_track[cnum[threadIdx.x]].length=localtrack[threadIdx.x].length;	
 

		


		for( map_index[threadIdx.x]=0;map_index[threadIdx.x]<localtrack[threadIdx.x].nHits;map_index[threadIdx.x]++){
	out_track[cnum[threadIdx.x]].hitmap[map_index[threadIdx.x]]=localtrack[threadIdx.x].hitmap[map_index[threadIdx.x]];	
	}
	

 }while(1);




 


				// }
				// }}
				
 //__syncthreads();				
 //out[sindex[threadIdx.x]]=manager[sblockIdx[threadIdx.x]];				
			
			









	




}






__device__ float dtest(float *data){


return 0;
}





__device__ int small_segmentHitSelection ( small_FtfHit *baseHit, small_FtfHit *candidateHit,small_FtfTrack *tra ){


   float dx, dy, dr, d3, dangle ;
   float dphi, deta ;
   float   angle ;
  
//
//   select hit with the
//   the smallest value of d3 (defined below)
//
   dphi  = (float)fabs((baseHit->phi) - (candidateHit->phi)) ; 
   if ( dphi > pi ) dphi = (float)fabs( twoPi - dphi ) ;
   if ( dphi > getPara_dphi && dphi < twoPi -getPara_dphi ) return 0 ;
//
//    Make sure we want to look at the difference in eta
//
   if ( baseHit->dz < 1000. && candidateHit->dz < 1000. ){
        deta  = (float)fabs((baseHit->eta) - (candidateHit->eta)) ; 
        if ( deta > getPara_deta ) return 0 ;
   }
   else deta = 0.F ;
  
   dr    = (float)fabs((float)(baseHit->row - candidateHit->row));
   d3    = (float)(toDeg * dr * ( dphi  + deta ) ) ;
//
//     If initial segment is longer than 2 store angle info in 
//     a1Xy and a1_sz
//

  
   if ( getPara_nHitsForSegment > 2 && tra->nHits-1 < getPara_nHitsForSegment ) {
	  ;
      dx = candidateHit->x - baseHit->x ;
      dy = candidateHit->y - baseHit->y ;
      angle = (float)atan2 ( dy, dx ) ;
	 
      if ( angle < 0  ) angle = angle + twoPi ;
      tra->lastXyAngle = angle ;
   }

   if ( d3 < tra->chi2[0] ) {
//
//   For second hit onwards check the difference in angle 
//   between the last two track segments
//
    
      if ( tra->nHits > 1 ) {
	 dx     = candidateHit->x - baseHit->x ;
         dy     = candidateHit->y - baseHit->y ;
         angle  = (float)atan2 ( dy, dx ) ;
         if ( angle < 0  ) angle = angle + twoPi ;
	    dangle = (float)fabs ( tra->lastXyAngle - angle );
		  
	    tra->lastXyAngle = angle ;
         if ( dangle > getPara_segmentMaxAngle ) return 0 ;
      }
//
//    Check whether this is the "closest" hit
//
      tra->chi2[0]          = d3 ;
      if ( d3 < getPara_goodDistance ) return 2 ;
	  return 1 ;
   }
//
//    If hit does not fulfill criterai return 0
//
   return 0 ;
}


__device__ int small_followHitSelection ( small_FtfHit *baseHit, small_FtfHit *candidateHit,  int way, small_FtfTrack *tra ){
//
   float lszChi2 = 0 ;
   float lchi2 ;
   float slocal=0, deta, dphi ;
   float dx, dy, dxy, dsz, temp ;
//
//           Check delta eta 
//
//   if ( baseHit->dz < 1000. && candidateHit->dz < 1000 ){
      deta = fabs((baseHit->eta)-(candidateHit->eta)) ;
      if ( deta > getPara_deta ) return 0 ; 
//   }
//   else deta = 0.F ;
//
//           Check delta phi
//
  dphi = fabs((baseHit->phi)-(candidateHit->phi)) ;
  if ( dphi > getPara_dphi && dphi < twoPi-getPara_dphi ) return 0 ;
//
//      If looking for secondaries calculate conformal coordinates
//

//
//      Calculate distance in x and y
//
   temp = (tra->a2Xy * candidateHit->xp - candidateHit->yp + tra->a1Xy) ;
   dxy  = temp * temp / ( tra->a2Xy * tra->a2Xy + 1.F ) ;
//
//    Calculate chi2
//
   lchi2    = (dxy * candidateHit->wxy) ;

   if ( lchi2 > tra->chi2[0] ) return 0 ;
//
//      Now in the sz plane
//
   if ( getPara_szFitFlag ){
//
//        Get "s" and calculate distance hit-line
//
      dx     = baseHit->x - candidateHit->x ;
      dy     = baseHit->y - candidateHit->y ;
      slocal = baseHit->s - way * sqrt ( dx * dx + dy * dy ) ;

      temp = (tra->a2Sz * slocal - candidateHit->z + tra->a1Sz) ;
      dsz  = temp * temp / ( tra->a2Sz * tra->a2Sz + 1 ) ;
//
//              Calculate chi2
//
      lszChi2 = dsz * candidateHit->wz ;
      lchi2 += lszChi2 ;
   } 
   else {
      lszChi2 = 0.F ;
      //slocal = 0;
   }
//
//         Check whether the chi2 square is better than previous one
//
   if ( lchi2 < tra->chi2[0] ) {
      tra->chi2[0]       = (float)lchi2    ;
      tra->chi2[1]       = (float)lszChi2 ;
      
      if ( getPara_szFitFlag  ) candidateHit->s = (float)slocal ;
//
//       if a good chi2 is found let's stop here
//
      if ( lchi2 < getPara_goodHitChi2 ) return 2 ;

      return 1 ;
   }
//
//     Return the selected hit
//
   return 0 ;
}






__device__ void small_add ( small_FtfHit *thisHit, int way, small_FtfTrack *tra  )
{

//	cout<<"small_add  "<<thisHit->id<<endl;
//
//      Increment # hits in this track
//
  if(tra->nHits>=50) return;
  tra->hitmap[tra->nHits]=thisHit->id;


	tra->nHits++ ; 
	tra->lasthit.phi=	thisHit->phi;
	tra->lasthit.dz=	thisHit->dz;
	tra->lasthit.eta=	thisHit->eta;
	tra->lasthit.row=	thisHit->row;
	tra->lasthit.x=	thisHit->x;
	tra->lasthit.y=	thisHit->y;
	tra->lasthit.s=	thisHit->s;

	

//
//    Check whether a fit update is needed
//
  if ( tra->nHits < getPara_minHitsForFit ) return ;
//
//    Include hit in xy fit parameter calculation
//
  
  //if(thisHit->id==3584){
  //cout<<tra->s11Xy<<endl;
  //cout<<tra->s12Xy<<endl;
  //cout<<tra->s22Xy<<endl;
  //cout<<tra->g1Xy<<endl;
  //cout<<tra->g2Xy<<endl;
  //cout<<tra->a1Xy<<endl;
  //cout<<tra->a2Xy<<endl;
  //cout<<tra->a1Sz<<endl;
  //cout<<tra->a2Sz<<endl;
  //cout<<thisHit->wxy<<endl;
  //cout<<thisHit->xp<<endl;
  //cout<<thisHit->yp<<endl;
  //cout<<thisHit->s<<endl;
  //cout<<thisHit->z<<endl;
  //cout<<tra->nHits<<endl;

  //}


  tra->s11Xy = tra->s11Xy + thisHit->wxy ;
  tra->s12Xy = tra->s12Xy + thisHit->wxy * thisHit->xp ;
  tra->s22Xy = tra->s22Xy + thisHit->wxy * square(thisHit->xp) ;
  tra->g1Xy  = tra->g1Xy  + thisHit->wxy * thisHit->yp ;
  tra->g2Xy  = tra->g2Xy  + thisHit->wxy * thisHit->xp * thisHit->yp ;
  
 
  if ( tra->nHits > getPara_minHitsForFit  )
  {
     tra->ddXy  = tra->s11Xy * tra->s22Xy - square ( tra->s12Xy ) ;
     if ( tra->ddXy != 0 ) {
        tra->a1Xy  = ( tra->g1Xy * tra->s22Xy - tra->g2Xy * tra->s12Xy ) / tra->ddXy ;
        tra->a2Xy  = ( tra->g2Xy * tra->s11Xy - tra->g1Xy * tra->s12Xy ) / tra->ddXy ;
     }
     else {
		 //LOG(ERR, "FtfTrack:add: ddSz = 0 \n" ) ;
     }
  }
//
//     Now in the sz plane
//
  if ( getPara_szFitFlag ) {
     tra->s11Sz = tra->s11Sz + thisHit->wz ;
     tra->s12Sz = tra->s12Sz + thisHit->wz * thisHit->s ;
     tra->s22Sz = tra->s22Sz + thisHit->wz * thisHit->s * thisHit->s ;
     tra->g1Sz  = tra->g1Sz  + thisHit->wz * thisHit->z ;
     tra->g2Sz  = tra->g2Sz  + thisHit->wz * thisHit->s * thisHit->z ;
  
     if ( tra->nHits > getPara_minHitsForFit ) {
		
        tra->ddSz  = tra->s11Sz * tra->s22Sz -  tra->s12Sz * tra->s12Sz ;
	if ( tra->ddSz != 0 ) {
           tra->a1Sz  = ( tra->g1Sz * tra->s22Sz - tra->g2Sz * tra->s12Sz ) / tra->ddSz ;
           tra->a2Sz  = ( tra->g2Sz * tra->s11Sz - tra->g1Sz * tra->s12Sz ) / tra->ddSz ;
         }
         else
         {
            if ( getPara_infoLevel > 0 ) {
               //LOG(ERR, "FtfTrack:add: ddSz = 0 \n" ) ;
            }
         }
      }
   }
}




__device__ void small_reset ( small_FtfTrack *tra)
{
/*----------------------------------------------------------------------
                Set fit parameters to zero
----------------------------------------------------------------------*/
 
  //tra->flag     = getPara_primaries ;


  tra->nHits    = 0 ;
  tra->s11Xy   = 
  tra->s12Xy   = 
  tra->s22Xy   = 
  tra->g1Xy    = 
  tra->g2Xy    = 
  tra->chi2[0]  = 0.F ;
 
  //tra->nxatrk   = 0 ;
  if ( getPara_szFitFlag ) 
  {
     tra->s11Sz =
     tra->s12Sz =
     tra->s22Sz =
     tra->g1Sz  =
     tra->g2Sz  =
     tra->chi2[1]  = 
     tra->length         = 0.F ;
  }


}


// __device__ void small_track_manager::ini(){
	// for(int i=0;i<2000;i++){
		// this->track_id[i]=i;
	// }
	// this->maxused=0;
	// this->number=0;

// }




//void small_track_assign::ini(){
//
//	this->number=0;
//
//}
//
//
//int small_track_assign::push(int id){
//
//track_id[number]=id;
//number++;
//return number;
//
//}




__device__ int element_setbit(int & data,int position){
if(position>31) return -1;

unsigned short mask=0;
mask=1<<position;
data=data|mask;
return 1;

}


__device__ int element_getbit(int data,int position){

int sd=data>>position;
return sd%2;

}



__device__ int small_segmentHitGroup ( small_FtfHit *candidateHit,int num,small_FtfHit *selected_hit,int &selected_pos, small_FtfTrack *tra ){
		int result=0;
	for(int i=0;i<num;i++){
				result= small_segmentHitSelection(&tra->lasthit,&candidateHit[i],tra);
		             if ( result > 0 ) {
							selected_pos=i;
						   selected_hit->id=candidateHit[i].id;
						   selected_hit->row=candidateHit[i].row;

						   selected_hit->x=candidateHit[i].x;
						   selected_hit->y=candidateHit[i].y;
						   selected_hit->z=candidateHit[i].z;
						   selected_hit->xp=candidateHit[i].xp;
						   selected_hit->yp=candidateHit[i].yp;
						   selected_hit->eta=candidateHit[i].eta;
						   selected_hit->phi=candidateHit[i].phi;
						   selected_hit->wxy=candidateHit[i].wxy;
						   selected_hit->wz=candidateHit[i].wz;

						   selected_hit->dx=candidateHit[i].dx;
						   selected_hit->dy=candidateHit[i].dy;
						   selected_hit->dz=candidateHit[i].dz;
						   selected_hit->s=candidateHit[i].s;


                if ( result ==2  ) return 2 ; 
             }
	
	}

return result;

}


__device__ int small_followHitGroup (  small_FtfHit *candidateHit,  int way, int num,small_FtfHit *selected_hit,int &selected_pos, small_FtfTrack *tra ){


		int result=0;
	for(int i=0;i<num;i++){
				result= small_followHitSelection(&tra->lasthit,&candidateHit[i],way,tra);
		             if ( result > 0 ) {
							selected_pos=i;
						   selected_hit->id=candidateHit[i].id;
						   selected_hit->row=candidateHit[i].row;

						   selected_hit->x=candidateHit[i].x;
						   selected_hit->y=candidateHit[i].y;
						   selected_hit->z=candidateHit[i].z;
						   selected_hit->xp=candidateHit[i].xp;
						   selected_hit->yp=candidateHit[i].yp;
						   selected_hit->eta=candidateHit[i].eta;
						   selected_hit->phi=candidateHit[i].phi;
						   selected_hit->wxy=candidateHit[i].wxy;
						   selected_hit->wz=candidateHit[i].wz;

						   selected_hit->dx=candidateHit[i].dx;
						   selected_hit->dy=candidateHit[i].dy;
						   selected_hit->dz=candidateHit[i].dz;
						   selected_hit->s=candidateHit[i].s;


                if ( result ==2  ) return 2 ; 
             }
	
	}

return result;


}







