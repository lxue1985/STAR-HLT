// ana.h: interface for the ana class.
//
//////////////////////////////////////////////////////////////////////

#ifndef ana_H
#define ana_H

class ana{
public:
	void end();
	
	float fbuf[100];
	float info[100];
	int buf[100];
	int *offset;
	char ofile[200];  //  open file
	char wfile[200];  //  write file
//	bool write;
//	int lengthf;
	int maxinum;       // ndata selected
	int maxfnum;int infolength;
	int maxlength,nowlength;// ndata selected
	int ndataw;       // write ndata selected
//	bool write;        //write or not
	int nevent;        // nevent
//	int num;   
//	int numw; 
	int *idata;
	int **idata_block;
	int *i1;
	int *i2;
	int *i3;
	int *i4;
	float *fdata; 
	float **fdata_block;
	int *idataw;
	float *fdataw; 


//	void ana();
	int	write(int inum, int fnum,int infonum, int length,int maxlength, int *idata, float *fdata);
	int	start(int ddd=0);   //run
//	void    analysis();
    int	get(int);
 ~ana();
};


#endif // !defined(AFX_ANA_H__C3B702CF_A9BD_4D07_8ED3_F44C9DC9B104__INCLUDED_)
