#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/wait.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <netdb.h>
#include "data.h"
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
data::data()

{

    max_length=11000;
    max_number=255;
    midlength=0;
    number=0;
    length=3;
    source=new char[max_length];
    result=new char[max_length];
    
}

data::data(int aa)

{

    max_length=aa;
    max_number=255;
    result=new char[max_length];
    number=0;
    length=3;
    source=new char[max_length];
    
}

data::~data()

{

    delete []source;
    delete []result;
    
}

void data::start()

{

    number=0;
    length=3;
    
}

void data::end()

{

    result[0]=number;
    ((short*)(result+1))[0]=length;
    
}

int data::put(char * msource,int mlength,int mtype)

{

    if(length+midlength+mlength+4>max_length) return -1;
    if(number==255) return -2;
    number++;
    ((short*)(result+length+0))[0]=mtype;
    ((short*)(result+length+2))[0]=mlength+midlength;
    //((short *)result[length+2])[0]=mlength;
    for(int i=0;i<mlength;i++)
    {

        result[length+midlength+4+i]=255-msource[i];
        
    }

    length=length+midlength+mlength+4;
    midlength=0;
    return length;
    
}

int data::get()

{

    if(number==result[0]) return -1;
    for(int i=0;i<((short*)(result+length+2))[0];
    i++)
    {

        source[i]=255-result[length+4+i];
        
    }

    if(length+((short*)(result+length+2))[0]+4>max_length) return -2;
    type=((short*)(result+length+0))[0];
    mlength=((short*)(result+length+2))[0];
    length=length+mlength+4;
    number++;
    return result[0]-number;
    
}

int data::put(char * msource1,int mlength1,char * msource2,int mlength2,int mtype)

{

    if(length+mlength1+mlength2+4>max_length) return -1;
    if(number==255) return -2;
    number++;
    ((short*)(result+length+0))[0]=mtype;
    ((short*)(result+length+2))[0]=mlength1+mlength2;
    //((short *)result[length+2])[0]=mlength;
    int i;
    for(i=0;i<mlength1;i++)
    {

        result[length+4+i]=255-msource1[i];
        
    }

    for( i=0;i<mlength2;i++)
    {

        result[length+4+i+mlength1]=255-msource2[i];
        
    }

    length=length+mlength1+mlength2+4;
    return length;
    
}

int data::put(char * msource1,int mlength1,char * msource2,int mlength2,char * msource3,int mlength3,int mtype)

{

    if(length+mlength1+mlength2+mlength3+4>max_length) return -1;
    if(number==255) return -2;
    number++;
    ((short*)(result+length+0))[0]=mtype;
    ((short*)(result+length+2))[0]=mlength1+mlength2+mlength3;
    //((short *)result[length+2])[0]=mlength;
    int i;
    for(i=0;i<mlength1;i++)
    {

        result[length+4+i]=255-msource1[i];
        
    }

    for( i=0;i<mlength2;i++)
    {

        result[length+4+i+mlength1]=255-msource2[i];
        
    }

    for( i=0;i<mlength3;i++)
    {

        result[length+4+i+mlength1+mlength2]=255-msource3[i];
        
    }

    length=length+mlength1+mlength2+mlength3+4;
    return length;
    
}

int data::put(char * msource1,int mlength1,char * msource2,int mlength2,char * msource3,int mlength3,char * msource4,int mlength4,int mtype)

{

    if(length+mlength1+mlength2+mlength3+mlength4+4>max_length) return -1;
    if(number==255) return -2;
    number++;
    ((short*)(result+length+0))[0]=mtype;
    ((short*)(result+length+2))[0]=mlength1+mlength2+mlength3+mlength4;
    //((short *)result[length+2])[0]=mlength;
    int i;
    for(i=0;i<mlength1;i++)
    {

        result[length+4+i]=255-msource1[i];
        
    }

    for( i=0;i<mlength2;i++)
    {

        result[length+4+i+mlength1]=255-msource2[i];
        
    }

    for( i=0;i<mlength3;i++)
    {

        result[length+4+i+mlength1+mlength2]=255-msource3[i];
        
    }

    for( i=0;i<mlength4;i++)
    {

        result[length+4+i+mlength1+mlength2+mlength3]=255-msource4[i];
        
    }

    length=length+mlength1+mlength2+mlength3+mlength4+4;
    return length;
    
}

int data::put(char *msource, int mlength)

{

    if(length+mlength+4>max_length) return -1;
    if(number==255) return -2;
    int i;
    for(i=0;i<mlength;i++)
    {

        result[length+midlength+4+i]=255-msource[i];
        
    }

    midlength+=mlength;
    return length;
    
}


