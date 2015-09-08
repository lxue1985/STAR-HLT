// udp.cpp: implementation of the udp class.
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
#include "udp.h"
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
udp::udp()

{

    MAX_LENGTH=5000;
    
}

udp::~udp()

{

    
}

int udp::buildserver(unsigned short  port)
{

    //statemachine.server_port=port;
    //WSADATA wsaData;
    //if (WSAStartup(MAKEWORD(2, 2), &wsaData) != 0)
    //{
    ////	//Beep(10000,10);
    //}
    //else{
    //
    //}
    PrimaryUDP = socket(AF_INET, SOCK_DGRAM, 0);
    if (PrimaryUDP < 0)
    
    {

        /*
        //Beep(10000,10);
        //Sleep(1000);
        //Beep(10000,10);
        */
        
    }

    local.sin_family=AF_INET;
    local.sin_port= htons(port); 
    local.sin_addr.s_addr = htonl(INADDR_ANY);
    int nResult=bind(PrimaryUDP,(sockaddr*)&local,sizeof(sockaddr));
    if(nResult==-1)
    {

        buildserver(port+1);
        /*
        //Beep(10000,10);
        //Sleep(1000);
        //Beep(10000,10);
        //Sleep(1000);
        //Beep(10000,10);
        */	
        
    }

    return 0;
    
}

int udp::buildclient()
{

    //WSADATA wsaData;
    //if (WSAStartup(MAKEWORD(2, 2), &wsaData) != 0)
    //{
    //	//Beep(10000,10);
    //}
    //else{
    //
    //}
    PrimaryUDP = socket(AF_INET, SOCK_DGRAM, 0);
    if (PrimaryUDP < 0)
    
    {

        //Beep(10000,10);
        //Sleep(1000);
        //Beep(10000,10);
        
    }

    local.sin_addr.s_addr = htonl(INADDR_ANY);
    local.sin_family = AF_INET;
    local.sin_port = 0;
    if (bind(PrimaryUDP, (struct sockaddr*)&local, sizeof(local)) < 0)
    
    {

        //Beep(10000,10);
        //Sleep(1000);
        //Beep(10000,10);
        //Sleep(1000);
        //Beep(10000,10);
        
    }

    //	getsockname(PrimaryUDP, (struct sockaddr*)&local, &length);
    return 0;
    
}

int udp::send(data & dat)
{

    if(((short*)(dat.result+1))[0]==3) return 0;
    //	aaa<<"send IP is  "<<send_remote.sin_addr.S_un.S_addr<<endl;
    //aaa<<" send port is  "<<send_remote.sin_port<<endl;
    //	for(int i=0;i<10;i++){
    if	(sendto(PrimaryUDP, (const char*)dat.result, ((short*)(dat.result+1))[0], 0, (const sockaddr*)&send_remote, sizeof(send_remote))==((short*)(dat.result+1))[0])
    {

        return ((short*)(dat.result+1))[0];
        //	};
        //	//Beep(10000,20);
        
    }

    return -1;
    
}

int udp::receive()
{

    remote_length=sizeof(recv_remote);
    recv_length=recvfrom(PrimaryUDP, (char *)recv_pointer, MAX_LENGTH, 0, (sockaddr*)&recv_remote, (socklen_t*)&remote_length);
    return recv_length;
    
}

void udp::set(char* ServerIP,int ServerPort)

{

    inet_aton(ServerIP,&send_remote.sin_addr);
    send_remote.sin_family = AF_INET;
    send_remote.sin_port = htons(ServerPort);
    
}

int udp::reply()

{

    send_remote.sin_addr.s_addr = recv_remote.sin_addr.s_addr;
    send_remote.sin_family = AF_INET;
    send_remote.sin_port = recv_remote.sin_port;
    return 0;
    
}


