// tcp_client.cpp: implementation of the tcp_client class.
//
//////////////////////////////////////////////////////////////////////
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
#include "tcp_client.h"
using namespace std;
#define MAX_MSG 100
#define LINE_ARRAY_SIZE (MAX_MSG+1)
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
tcp_client::tcp_client()

{

    MAX_LENGTH=5000;
    SEG=10000;
    client_ID=0;
    data_valid=false;
    
}

tcp_client::~tcp_client()

{

    
}

int tcp_client::buildserver(char * IP,unsigned short port)
{

    hostInfo = gethostbyname(IP);
    if (hostInfo == NULL) 
    {

        cout << "problem interpreting host: " <<"\n";
        exit(1);
        
    }

    socketDescriptor = socket(AF_INET, SOCK_STREAM, 0);
    if (socketDescriptor < 0) 
    {

        cerr << "cannot create socket\n";
        exit(1);
        
    }		   

    int optval = 1 ;
    setsockopt(socketDescriptor,SOL_SOCKET, SO_KEEPALIVE, (char *)&optval, sizeof(optval)) ;
    setsockopt(socketDescriptor,SOL_SOCKET, SO_REUSEADDR, (char *)&optval, sizeof(optval)) ;	  
    serverAddress.sin_family = hostInfo->h_addrtype;
    memcpy((char *) &serverAddress.sin_addr.s_addr,
    hostInfo->h_addr_list[0], hostInfo->h_length);
    serverAddress.sin_family = AF_INET;
    serverAddress.sin_port = htons(port);
    serverAddress.sin_addr.s_addr = inet_addr("131.243.31.213");
    if (connect(socketDescriptor,
    (struct sockaddr *) &serverAddress,
    sizeof(serverAddress)) < 0) 
    {

        cerr << "cannot connect\n";
        exit(1);
        
    }		   

    data recv_data;
    while (1) 
    {

        int ds=recv(socketDescriptor, recv_data.result, recv_data.max_length, 0);
        if(ds<=0)
        {

            cout<<"recv error "<<ds<<endl;
            break;
            
        }

        on_receive(recv_data);	
        // cout << "  --  " << recv_data.result << "\n";
        // Convert line to upper case.
        
    }			   

    return 0;
    
}

void tcp_client::buildserver()

{

    pthread_t bthreads;
    pthread_attr_t bpthread_custom_attr;
    pthread_attr_init(&bpthread_custom_attr);
    pthread_create(&bthreads, &bpthread_custom_attr, tcp_client_build_function, (void *)(this));					   
    
}

int tcp_client::on_receive(data &dat)
{

    // int i;
    // for (i = 0; dat.result[i] != '\0'; i++)
    // dat.result[i] = toupper(dat.result[i]);
    // if (send(p[num].connectSocket, dat.result, strlen(dat.result) + 1, 0) < 0)
    // cerr << "Error: cannot send modified data";
    dat.start();
    while(dat.get()>=0)
    {

        //	cout<<"type "<<dat.type<<endl;
        if(dat.type==1)
        {

            cout<<"receive  "<<dat.source<<endl;
            
        }

        if(dat.type==10)
        {

            int *plength=(int *)dat.source;
            int i;
            for(i=0;i<dat.mlength-4;i++)
            {

                data_buf[i+plength[0]*SEG]=dat.source[4+i];
                
            }

            this->data_valid=false;
            
        }

        if(dat.type==11)
        {

            int *plength=(int *)dat.source;
            int i;
            for(i=0;i<dat.mlength-4;i++)
            {

                data_buf[i+plength[0]*SEG]=dat.source[4+i];
                
            }

            this->data_valid=true;
            this->data_length=plength[0]*SEG+dat.mlength-4;
            cout<<"total "<<this->data_length<<endl;
            cout<<"content "<<data_buf<<endl;
            
        }

        if(dat.type==12)
        {

            data sd;
            sd.start();
            sd.put((char*)&this->client_ID,4,13);
            sd.end();
            this->on_send(sd);
            
        }

        if(dat.type==14)
        {

            send_flag=true;
            
        }

        //cout<<"next type "<<dat.type<<endl;
        
    }	     	

    return 0;
    
}

int tcp_client::on_send(data &dat)
{

    if (send(socketDescriptor, dat.result, ((short*)(dat.result+1))[0], 0) < 0)
    {

        cerr << "Error: cannot send modified data";
        return -1;
        
    }

    //cout<<"send"<<((short*)(dat.result+1))[0]<<endl;
    return 0;
    
}

void* tcp_client_server_function(void *arg)

{

    return (NULL);
    
}

void* tcp_client_build_function(void *arg)

{

    tcp_client *p=(tcp_client *)arg;
    p->buildserver("131.243.31.213",1600);
    return (NULL);
    
}

int tcp_client::send_data(char* buf,int length)
{

    int integer=length/SEG;
    int res=length%SEG;
    int i;
    data dat;
    for(i=0;i<integer;i++)
    {

        dat.start();
        dat.put((char*)&i,4,(buf+SEG*i),SEG,10);
        dat.end();
        this->on_send(dat);
        send_flag=false;
        int ii=0;
        do
        {

            ii++;
            if(ii>100000000||send_flag)break;
            
        }
        while(1);
        
    }

    dat.start();
    dat.put((char*)&i,4,(buf+SEG*i),res,11);
    dat.end();
    this->on_send(dat);
    int ii=0;
    do
    {

        ii++;
        if(ii>100000000||send_flag)break;
        
    }
    while(1);
    return 0;
    
}


