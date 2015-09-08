// ana.cpp: implementation of the ana class.
//
//////////////////////////////////////////////////////////////////////

#include "ana.h"
#include <stdio.h>
#include <iostream.h>
#include <fstream.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////



ana::~ana()
{

}
int ana::start(int ddd){
int	bre;
	int length,inum,fnum,infol;
maxinum=0;
maxfnum=0;
maxlength=0;
FILE *p;
p=fopen(ofile,"rb");
if(p==0)return 0;
//	for(int i=0;i<lengthf;i++){
int i=0;do{
		if(i!=0){	
bre=fread(&inum,4,1,p);
nevent=i;
if(bre==0) break;

	fread(&fnum,4,1,p);


	fread(&length,4,1,p);
	fread(&infol,4,1,p);
	if(inum>maxinum)maxinum=inum;
	if(fnum>maxfnum)maxfnum=fnum;
	if(length>maxlength)maxlength=length;
	fseek(p,4*((inum+fnum)*length+infol),1);
//		offset[i]=offset[i-1]+4*(inum+fnum)*length+16+4*infol;
		}

//		else{offset[0]=0;}
i++;	}while(ddd==0||i<ddd+1);

offset=new int[nevent];

fclose(p);
p=fopen(ofile,"rb");

i=0;do{
		if(i!=0){	
bre=fread(&inum,4,1,p);
nevent=i;
if(bre==0) break;

	fread(&fnum,4,1,p);


	fread(&length,4,1,p);
	fread(&infol,4,1,p);
	if(inum>maxinum)maxinum=inum;
	if(fnum>maxfnum)maxfnum=fnum;
	if(length>maxlength)maxlength=length;
	fseek(p,4*((inum+fnum)*length+infol),1);
		offset[i]=offset[i-1]+4*(inum+fnum)*length+16+4*infol;}

		else{offset[0]=0;}
i++;	}while(ddd==0||i<ddd+1);
fclose(p);

idata=new int[maxinum*maxlength];
idata_block=new int*[maxinum];
fdata=new float[maxfnum*maxlength];
fdata_block=new float*[maxfnum];

nevent--;
return 1;
}


int ana::get(int n){
FILE *p;
int i;
p=fopen(ofile,"rb");

int bbb=fseek(p,offset[n],0);
	fread(&maxinum,4,1,p);
	fread(&maxfnum,4,1,p);
	fread(&nowlength,4,1,p);
	fread(&infolength,4,1,p);
if(bbb!=0) return  0;
if(infolength!=0) fread(info,4*infolength,1,p);
if(maxinum!=0)	{fread(idata,4*maxinum*nowlength,1,p);

for(i=0;i<maxinum;i++){
idata_block[i]=idata+i*nowlength;
}


}
if(maxfnum!=0)	{fread(fdata,4*maxfnum*nowlength,1,p);

for(i=0;i<maxfnum;i++){
fdata_block[i]=fdata+i*nowlength;
}
}
	//num=maxlength;
//	cout<<n<<"n=="<<nowlength<<endl;
fclose(p);
return  nowlength;
}

int ana::write(int inum, int fnum,int infonum, int length,int maxlength, int *idata, float *fdata)
{
	int i;
FILE* p=fopen(wfile,"ab");
//else FILE* p=fopen(wfile,"w");
	fwrite( &inum,4,1,p);
	fwrite( &fnum,4,1,p);
	fwrite( &length,4,1,p);
	fwrite( &infonum,4,1,p);
	
if(infonum!=0)	fwrite( info,4*infonum,1,p);
	
	for( i=0;i<inum;i++)
	{	
	fwrite( idata+i*maxlength,4*length,1,p);
	}
	for( i=0;i<fnum;i++)
	{

		fwrite( fdata+i*maxlength,4*length,1,p);
	}
		fclose(p);

return 0;
}

void ana::end()
{
delete []offset;
delete []idata;
delete []idata_block;
delete []fdata;
delete []fdata_block;
}
