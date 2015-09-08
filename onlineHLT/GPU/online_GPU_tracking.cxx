#include "online_GPU_tracking.h"
#include "lib.cuh";
#include "online_tracking_sector.h"
#include "clib_kernel.h"
#include "online_tracking_TpcHitMap.h"

#define MAXTRACK 500
#define MAXHIT 8000


online_GPU_tracking::online_GPU_tracking(void)
{
}

online_GPU_tracking::~online_GPU_tracking(void)
{
}

int online_GPU_tracking::setup_memory(int dup_in ){

		  dup=dup_in;

		  max_hit=MAXHIT;


		hit_address=new FtfHit*[max_hit*dup];
		position=new float[max_hit*4*dup] ;
		hash=new int[41*11*46*2*dup];
		map=new int[1*dup];
		manager=new int[30+30];
		balance=manager+30;


		out=new float[sizeof(small_FtfTrack)*MAXTRACK*dup] ;
		hit_counter=new int[dup];
		track_number=new int[dup];

		ptracker=new online_tracking_sector*[dup];





		allocateArray((void**)&Dposition, sizeof(float)*max_hit*4*dup);
		allocateArray((void**)&Dhash, sizeof(int)*41*11*46*2*dup);
		allocateArray((void**)&Dmap, sizeof(int)*1*dup);
		allocateArray((void**)&Dout, sizeof(float)*sizeof(small_FtfTrack)*MAXTRACK*dup);
		allocateArray((void**)&Dmanager, sizeof(int)*1*dup);
		allocateArray((void**)&Dnumber, sizeof(int)*1*dup);
	return 0;
}
int	online_GPU_tracking::free_memory(){
		delete []hit_address;
		delete []position;
		delete []hash;
		delete []map;
		delete []manager;
		delete []out;
		delete []hit_counter;
		delete []track_number;
		freeArray(Dposition);
		freeArray(Dhash);
		freeArray(Dmap);
		freeArray(Dout);
		freeArray(Dmanager);
		freeArray(Dnumber);
		// int i=0;
		// for(i=0;i<dup;i++){
		// delete ptracker[i];
		// }
		delete []ptracker;

}



int   online_GPU_tracking::load_hit_to_buffer( int pos,online_tracking_sector *tracker  ){
// cout<<"load_hit_to_buffer  "<<pos<<"  nhits  "<<tracker->nHits<<endl;
this->ptracker[pos]=tracker;
 int idup=pos;
 hit_counter[idup]=0;
 

				for(int rhoi=45;rhoi>=0;rhoi--){

				for(int etai=0;etai<=40;etai++){
				for(int phii=0;phii<=10;phii++){

				




 
 	    int areaIndex = (rhoi)  * (41*11) + 
        (phii) * 41 + (etai) ;
 	    int hashIndex =  (rhoi)  * (41*11) + 
        (phii) * 41 + (etai) ;

		//hit_hash[hashIndex].start=hit_counter;

		
				//for(int idup=0;idup<dup;idup++){
			   
				hash[41*11*46*2*idup+hashIndex*2]=hit_counter[idup];
			   //}


         for ( FtfHit *candidateHit = (FtfHit *)tracker->volumeC[areaIndex].first ; 
             candidateHit != 0 ;
             candidateHit = (FtfHit *)candidateHit->nextVolumeHit ){

			   //smhit[hit_counter].x=candidateHit->x;
			   //smhit[hit_counter].y=candidateHit->y;
			   //smhit[hit_counter].z=candidateHit->z;


			   //for(int idup=0;idup<dup;idup++){
			   
			   hit_address[max_hit*idup+hit_counter[idup]]=candidateHit;

			   position[max_hit*4*idup+4*hit_counter[idup]+0]=candidateHit->x;
			   position[max_hit*4*idup+4*hit_counter[idup]+1]=candidateHit->y;
			   position[max_hit*4*idup+4*hit_counter[idup]+2]=candidateHit->z;
			   position[max_hit*4*idup+4*hit_counter[idup]+3]=(float)hashIndex;

			   
			   //}



	
			 hit_counter[idup]++;

		 }



		 //hit_hash[hashIndex].end=hit_counter;
		 //hit_hash[hashIndex].map=0;
			   //for(int idup=0;idup<dup;idup++){
			   
			    hash[41*11*46*2*idup+hashIndex*2+1]=hit_counter[idup];
				 //map[41*11*46*idup+hashIndex]=0;
			   //}



				}}}








return hit_counter[idup];
}

int online_GPU_tracking::fire_GPU(int cores){
		int i=0;
		for(i=0;i<dup;i++){
		//ptracker[i]=new online_tracking_sector(-0.5,"/RTS/conf/L3/map.bin",0);
		manager[i]=0;
		map[i]=0;

		}

copyArrayToDevice(Dmanager, manager, 0, sizeof(int)*(30+30));

copyArrayToDevice(Dposition, position, 0, sizeof(float)*max_hit*4*dup);

copyArrayToDevice(Dhash, hash, 0, sizeof(int)*41*11*46*2*dup);

copyArrayToDevice(Dmap, map, 0, sizeof(int)*1*dup);
copyArrayToDevice(Dnumber, hit_counter, 0, sizeof(int)*1*dup);

Gtest(Dposition,Dhash,Dmap,Dmanager,Dout,Dnumber,cores);


return 0;

}

int online_GPU_tracking::copy_GPU_back(){



copyArrayFromDevice(out, Dout,  sizeof(float)*sizeof(small_FtfTrack)*MAXTRACK*dup);
copyArrayFromDevice(manager, Dmanager,  sizeof(int)*1*dup);

int idup=0;
int icore=0;
for(idup=0;idup<dup;idup++){
track_number[idup]=manager[idup];

// cout<<idup<<" track number "<<track_number[idup]<<endl;
}


return dup;
}


int online_GPU_tracking::build_track( int pos ){

	online_tracking_sector *tracker=this->ptracker[pos];
 int idup=pos;


int NPI=4;
int Di=0;
small_FtfTrack *pout=(small_FtfTrack *)(out+idup*sizeof(small_FtfTrack)*MAXTRACK);

for(NPI=0;NPI<this->track_number[idup];NPI++){

small_FtfTrack *ptrack=&pout[NPI];

//cout<<idup<<" build_track "<<track_number[idup]<<endl;

for( Di=0;Di<ptrack->nHits;Di++){

	if(Di<ptrack->nHits-1)
	(hit_address[max_hit*idup+ptrack->hitmap[Di]])->nextTrackHit=hit_address[max_hit*idup+ptrack->hitmap[Di+1]];
	else
	(hit_address[max_hit*idup+ptrack->hitmap[Di]])->nextTrackHit=NULL;
	
	}

 tracker->track[NPI].para     = &tracker->para ;
 tracker->track[NPI].reset();
 tracker->track[NPI].nHits=ptrack->nHits;
 tracker->track[NPI].lastXyAngle=ptrack->lastXyAngle;	
 tracker->track[NPI].s11Xy=ptrack->s11Xy;	
 tracker->track[NPI].s12Xy=ptrack->s12Xy;	
 tracker->track[NPI].s22Xy=ptrack->s22Xy;	
 tracker->track[NPI].g1Xy=ptrack->g1Xy;	
 tracker->track[NPI].g2Xy=ptrack->g2Xy;	
 tracker->track[NPI].s11Sz=ptrack->s11Sz;	
 tracker->track[NPI].s12Sz=ptrack->s12Sz;	
 tracker->track[NPI].s22Sz=ptrack->s22Sz;	
 tracker->track[NPI].g1Sz=ptrack->g1Sz;	
 tracker->track[NPI].g2Sz=ptrack->g2Sz;	
 tracker->track[NPI].ddXy=ptrack->ddXy;	
 tracker->track[NPI].a1Xy=ptrack->a1Xy;	
 tracker->track[NPI].a2Xy=ptrack->a2Xy;	
 tracker->track[NPI].a2Sz=ptrack->a2Sz;	
 tracker->track[NPI].a1Sz=ptrack->a1Sz;	
 tracker->track[NPI].chi2[0]=ptrack->chi2[0];	
 tracker->track[NPI].chi2[1]=ptrack->chi2[1];	
 tracker->track[NPI].length=ptrack->length;


  tracker->track[NPI].firstHit=hit_address[max_hit*idup+ptrack->hitmap[0]];
  tracker->track[NPI].lastHit=hit_address[max_hit*idup+ptrack->hitmap[0]];
  tracker->track[NPI].currentHit=hit_address[max_hit*idup+ptrack->hitmap[0]];


            





 tracker->track[NPI].fill();
 tracker->fDedx->TruncatedMean(&tracker->track[NPI]);

//cout<<"tracker->track[NPI].pt  "<<tracker->track[NPI].pt<<endl;
//cout<<"tracker->track[NPI].dedx  "<<tracker->track[NPI].dedx<<endl;
}
tracker->nTracks=this->track_number[idup];
return 0;
}





int online_GPU_tracking::set_balance(int cores){
if(cores>30) cores=30;
if(cores<this->dup) cores=this->dup;
int total_hits=0;
int core_used=this->dup;
int core_rest=cores-this->dup;




int i=0;
for(i=0;i<this->dup;i++){
	total_hits=total_hits+this->ptracker[i]->nHits;
	balance[i]=i;
}

for(i=0;i<this->dup;i++){

	int core_for_this=int(float(core_rest*this->ptracker[i]->nHits)/float(total_hits)+0.5);

	if(core_for_this==0)core_for_this=1;
	//cout<<"balance "<<i<<" nhits  "<<this->ptracker[i]->nHits<<" = "<<core_for_this+1<<endl;
	int j=0;

	for(j=core_used;j<core_used+core_for_this;j++){
		balance[j]=i;
	}
	core_rest=core_rest-core_for_this;
	core_used=core_used+core_for_this;
	total_hits=total_hits-this->ptracker[i]->nHits;
}



return total_hits+core_rest;
}