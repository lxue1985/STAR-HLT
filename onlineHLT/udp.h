// udp.h: interface for the udp class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_UDP_H__ABC7954F_9495_426A_A3F2_05DD3DE34998__INCLUDED_)
#define AFX_UDP_H__ABC7954F_9495_426A_A3F2_05DD3DE34998__INCLUDED_
#include <arpa/inet.h>
#include <netdb.h>
#include <netinet/in.h>
#include <unistd.h>
#include <iostream>
#include "data.h"
class udp  

{

    public:
    udp();
    virtual ~udp();
    int buildserver(unsigned short  port);
    int buildclient();
    int send(data & dat);
    int receive();
    public:
    int reply();
    void set(char * ServerIP,int ServerPort);
    int PrimaryUDP;
        //change

    int MAX_LENGTH;
    int ServerIP;
    unsigned short ServerPort; 
    sockaddr_in local;
    sockaddr_in send_remote;
    sockaddr_in recv_remote;
    int send_length;
    char * send_pointer;
    int recv_length;
    int remote_length;
    char * recv_pointer;
    
}
;

#endif // !defined(AFX_UDP_H__ABC7954F_9495_426A_A3F2_05DD3DE34998__INCLUDED_)

