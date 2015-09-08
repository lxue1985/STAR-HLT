// tcp_server.cpp: implementation of the tcp_server class.
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
#include <netinet/tcp.h>
#include "data.h"
#include "tcp_server.h"
using namespace std;
#define MAX_MSG 100
#define LINE_ARRAY_SIZE (MAX_MSG+1)
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
tcp_server::tcp_server()

{

    SEG=10000;
    MAX_LENGTH=5000;
    connections=24;
    p=(parm *)malloc(sizeof(parm)*connections);
    threads=(pthread_t *)malloc(connections*sizeof(*threads));
    client_ID=0;
    data_buf=new char*[connections];
    ID2index=new int[connections];
    int i;
    for(i=0;i<connections;i++)
    {

        ID2index[i]=-1;
        
    }

    
}

tcp_server::~tcp_server()

{

    free(p);
    free(threads);
    delete []data_buf;
    delete []ID2index;
    
}

int tcp_server::buildserver(unsigned short port)
{

    pthread_attr_init(&pthread_custom_attr);	
    Primarytcp_server = socket(AF_INET, SOCK_STREAM, 0);
    int optval = 1 ;
    setsockopt(Primarytcp_server,SOL_SOCKET, SO_KEEPALIVE, (char *)&optval, sizeof(optval)) ;
    setsockopt(Primarytcp_server,SOL_SOCKET, SO_REUSEADDR, (char *)&optval, sizeof(optval)) ;	  
    if (Primarytcp_server < 0) 
    {

        cerr << "cannot create listen socket";
        exit(1);
        
    }

    local.sin_family = AF_INET;
    local.sin_addr.s_addr = htonl(INADDR_ANY);
    local.sin_port = htons(port);
    ServerPort=htons(port);
    if (bind(Primarytcp_server,
    (struct sockaddr *) &local,
    sizeof(local)) < 0) 
    {

        cerr << "cannot bind socket";
        exit(1);
        
    }  

    listen(Primarytcp_server, connections);
    while (1) 
    {

        cout << "Waiting for TCP connection on port " << Primarytcp_server << endl;
        // Accept a connection with a client that is requesting one.  The
        // accept() call is a blocking call; i.e., this thread of
        // execution stops until a connection comes in.
        // connectSocket is a new socket that the system provides,
        // separate from listenSocket.  We *could* accept more
        // connections on listenSocket, before connectSocket is closed,
        // but this program doesn't do that.
        int clientAddressLength = sizeof(sockaddr_in);
        int free_client;
        free_client=this->find_free_socket();
        cout<<"free_client  "<<free_client<<endl;
        if(free_client>=connections) break;
        if(free_client<0) break;
        p[free_client].index=free_client;
        p[free_client].tcp=this;
        p[free_client].connectSocket = accept(Primarytcp_server,(struct sockaddr *) &p[free_client].clientAddress,(socklen_t*)&clientAddressLength);
        if(p[free_client].connectSocket!=-1)
        {					   

            pthread_create(&threads[free_client], &pthread_custom_attr, tcp_server_server_function, (void *)(p+free_client));			   
            
        }
        else
        {

            cout <<"accept error"<< " ...\n";
            
        }

        
    }	

    return 0;
    
}

void tcp_server::buildserver()

{

    parm bp;
    pthread_t bthreads;
    pthread_attr_t bpthread_custom_attr;
    pthread_attr_init(&bpthread_custom_attr);
    pthread_create(&bthreads, &bpthread_custom_attr, tcp_server_build_function, (void *)(this));					   
    
}

int tcp_server::find_free_socket()
{

    int keepAlive = 1;
        //设定KeepAlive

    int keepIdle = 5;
        //开始首次KeepAlive探测前的TCP空闭时间   

    int keepInterval = 5;
        //两次KeepAlive探测间的时间间隔   

    int keepCount = 3;
        //判定断开前的KeepAlive探测次数   

    int i=0;
    for(i=0;i<connections;i++)
    {

        if(setsockopt(this->p[i].connectSocket,SOL_SOCKET,SO_KEEPALIVE,(void*)&keepAlive,sizeof(keepAlive)) == -1)   
        
        {   

            return i;
            
        }   

        if(setsockopt(this->p[i].connectSocket,SOL_TCP,TCP_KEEPIDLE,(void *)&keepIdle,sizeof(keepIdle)) == -1)   
        
        {   

            return i;
            
        }   

        if(setsockopt(this->p[i].connectSocket,SOL_TCP,TCP_KEEPINTVL,(void *)&keepInterval,sizeof(keepInterval)) == -1)   
        
        {   

            return i;
            
        }   

        if(setsockopt(this->p[i].connectSocket,SOL_TCP,TCP_KEEPCNT,(void *)&keepCount,sizeof(keepCount)) == -1)   
        
        {   

            return i;
            
        }  

        
    }

    return -1;
    
}

int tcp_server::on_receive(int num,data &dat)
{

    dat.start();
    while(dat.get()>=0)
    {

        if(dat.type==1)
        {

            on_send(num,dat );
            cout<<dat.source<<endl;
            
        }

        if(dat.type==10)
        {

            int *plength=(int *)dat.source;
            int i;
            for(i=0;i<dat.mlength-4;i++)
            {

                //data_buf[p[num].client_ID][i+plength[0]*SEG]=dat.source[4+i];
                
            }

            this->data_valid=false;
            data flg;
            flg.start();
            flg.put(0,0,14);
            flg.end();
            on_send(num,flg);
            //cout<<"parts "<<endl;
            
        }

        if(dat.type==11)
        {

            data flg;
            flg.start();
            flg.put(0,0,14);
            flg.end();
            on_send(num,flg);	
            int *plength=(int *)dat.source;
            int i;
            for(i=0;i<dat.mlength-4;i++)
            {

                //data_buf[p[num].client_ID][i+plength[0]*SEG]=dat.source[4+i];
                
            }

            this->data_valid=true;
            this->data_length=plength[0]*SEG+dat.mlength-4;
            cout<<"total "<<this->data_length<<endl;
            //cout<<"content "<<data_buf[p[num].client_ID]<<endl;
            send_data(p[num].client_ID,"server got it",strlen("server got it")+1);
            
        }

        if(dat.type==13)
        {

            int *plength=(int *)dat.source;
            this->p[num].client_ID=plength[0];
            ID2index[this->p[num].client_ID]=num;
            cout<<"in "<<num<<" client_ID "<<plength[0]<<endl;
            
        }

        if(dat.type==14)
        {

            p[num].send_flag=true;
            
        }

        
    }	

    // int i;
    // for (i = 0; dat.result[i] != '\0'; i++)
    // dat.result[i] = toupper(dat.result[i]);
    // if (send(p[num].connectSocket, dat.result, strlen(dat.result) + 1, 0) < 0)
    // cerr << "Error: cannot send modified data";
    return 0;
    
}

int tcp_server::on_send(int num,data &dat)
{

    if (send(p[num].connectSocket, dat.result, ((short*)(dat.result+1))[0], 0) < 0)
    {

        cerr << "Error: cannot send modified data";
        return -1;
        
    }

    return 0;
    
}

void* tcp_server_server_function(void *arg)

{

    parm *p=(parm *)arg;
    if (p->connectSocket < 0) 
    {

        cerr << "cannot accept connection ";
        exit(1);
        
    }

    int i;
    // Show the IP address of the client.
    // inet_ntoa() converts an IP address from binary form to the
    // standard "numbers and dots" notation.
    cout << "  connected to " << inet_ntoa(p->clientAddress.sin_addr);
    // Show the client's port number.
    // ntohs() converts a short int from network byte order (which is
    // Most Significant Byte first) to host byte order (which on x86,
    // for example, is Least Significant Byte first).
    cout << ":" << ntohs(p->clientAddress.sin_port) << "\n";
    // Read lines from socket, using recv(), storing them in the line
    // array.  If no messages are currently available, recv() blocks
    // until one arrives.
    // First set line to all zeroes, so we'll know where the end of
    // the string is.
    data recv_data;
    recv_data.start();
    recv_data.put(0,0,12);
    recv_data.end();
    p->tcp->on_send(p->index,recv_data);
    while (1) 
    {

        int ds=recv(p->connectSocket, recv_data.result, recv_data.max_length, 0);
        if(ds<=0)
        {

            p->connectSocket=0;
            cout<<"recv error "<<ds<<endl;
            p->tcp->ID2index[p->client_ID]=-1;
            break;
            
        }

        p->tcp->on_receive(p->index,recv_data);	
        //cout << "  --  " << recv_data.result << "\n";
        // Convert line to upper case.
        
    }	

    return (NULL);
    
}

void* tcp_server_build_function(void *arg)

{

    tcp_server *p=(tcp_server *)arg;
    p->buildserver(1600);
    return (NULL);
    
}

int tcp_server::send_data(int num, char* buf,int length)
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
        this->on_send(num,dat);
        p[num].send_flag=false;
        int ii=0;
        do
        {

            ii++;
            if(ii>100000000||p[num].send_flag)break;
            
        }
        while(1);
        
    }

    dat.start();
    dat.put((char*)&i,4,(buf+SEG*i),res,11);
    dat.end();
    this->on_send(num,dat);
    int ii=0;
    do
    {

        ii++;
        if(ii>100000000||p[num].send_flag)break;
        
    }
    while(1);
    return 0;
    
}


