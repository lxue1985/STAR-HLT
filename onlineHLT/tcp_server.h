// tcp_server.h: interface for the tcp_server class.
//
//////////////////////////////////////////////////////////////////////
#include <arpa/inet.h>
#include <netdb.h>
#include <netinet/in.h>
#include <unistd.h>
#include <iostream>
#include <pthread.h>

#if !defined(AFX_tcp_server_H__ABC7954F_9495_426A_A3F2_05DD3DE34998__INCLUDED_)
#define AFX_tcp_server_H__ABC7954F_9495_426A_A3F2_05DD3DE34998__INCLUDED_
#include "data.h"
void *tcp_server_server_function(void *arg);
void* tcp_server_build_function(void *arg);
class tcp_server;
typedef struct 
{

    sockaddr_in clientAddress;
    int connectSocket;
    int index;
    int client_ID;
    tcp_server * tcp;
    bool send_flag;
    
} 
parm;
class tcp_server  

{

    public:
    int * ID2index;
    int SEG;
    parm *p;
    pthread_t *threads;
    pthread_attr_t pthread_custom_attr;
    int connections;
    int client_ID;
    char ** data_buf;
    bool data_valid;
    int data_length;
    int Primarytcp_server;
    sockaddr_in local;	
    int ServerIP;
    unsigned short ServerPort; 
    void buildserver();
    tcp_server();
    virtual ~tcp_server();
    int buildserver(unsigned short port);
    int on_receive(int,data &dat);
    int find_free_socket();
    int on_send(int,data &dat);
    int send_data(int num, char* buf,int length);
    public:
    int MAX_LENGTH;
    
}
;

#endif // !defined(AFX_tcp_server_H__ABC7954F_9495_426A_A3F2_05DD3DE34998__INCLUDED_)

