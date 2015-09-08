// tcp_client.h: interface for the tcp_client class.
//
//////////////////////////////////////////////////////////////////////
#include <arpa/inet.h>
#include <netdb.h>
#include <netinet/in.h>
#include <unistd.h>
#include <iostream>
#include <pthread.h>

#if !defined(AFX_tcp_client_H__ABC7954F_9495_426A_A3F2_05DD3DE34998__INCLUDED_)
#define AFX_tcp_client_H__ABC7954F_9495_426A_A3F2_05DD3DE34998__INCLUDED_
#include "data.h"
void *tcp_client_server_function(void *arg);
void* tcp_client_build_function(void *arg);
class tcp_client  

{

    public:
    bool send_flag;
    int SEG;
    int socketDescriptor;
    unsigned short int serverPort;
    struct sockaddr_in serverAddress;
    struct hostent *hostInfo;
    int send_data(char*,int);
    char * data_buf;
    bool data_valid;
    int data_length;	
    pthread_t threads;
    pthread_attr_t pthread_custom_attr;
    int client_ID;
    int Primarytcp_client;
    unsigned short ServerPort; 
    void buildserver();
    tcp_client();
    virtual ~tcp_client();
    int buildserver(char *,unsigned short port);
    int on_receive(data &dat);
    int on_send(data &dat);
    public:
    int MAX_LENGTH;
    
}
;

#endif // !defined(AFX_tcp_client_H__ABC7954F_9495_426A_A3F2_05DD3DE34998__INCLUDED_)

