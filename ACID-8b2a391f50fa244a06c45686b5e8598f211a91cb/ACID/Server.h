#pragma once

#include <iostream>
#include <Ws2tcpip.h>
#include <winsock2.h>
#include <intrin.h>
#include <stdio.h>
#include <string>
#include <vector>

#pragma comment(lib, "ws2_32.lib")     //winsock library file

#define BUFFER_SIZE   1024


class Server
{
private:
	SOCKET       in;
	SOCKET       BindSocket(u_short port); // Create socket
	std::string  m_ipAddr;
	int          m_port;

	void Cleanup(); // cleanup Socket

public:
	Server();
	Server(std::string ipAddr, int port );

	~Server();
	
	int Init();		// Initialize Winsock;
	int Receive(std::vector< std::vector<double> > &dataHash, int iteration = 1000, int id=NULL);		// Receive message from client
};
