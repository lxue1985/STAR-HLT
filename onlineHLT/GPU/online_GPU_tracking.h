#ifndef online_GPU_tracking_
#define online_GPU_tracking_
class  FtfHit;
class online_tracking_sector;



class online_GPU_tracking
{
public:
	online_GPU_tracking(void);
	~online_GPU_tracking(void);
int	setup_memory(int dup_in);
int	free_memory();
int   load_hit_to_buffer( int pos,online_tracking_sector *tracker  );

int copy_GPU_back();
int fire_GPU(int cores);
int build_track( int pos );
int set_balance(int cores);

		int dup;

		int max_hit;

		FtfHit **hit_address;
		float  *position;
		int *hash;
		int *map;
		int *manager;
		int *balance;
		float  *out;

		float * Dposition;
		int * Dhash;
		int * Dmap;
		float *Dout;
		int *Dmanager;
		int *Dnumber;


		int *hit_counter;
		int *track_number;
		online_tracking_sector **ptracker;



};



#endif