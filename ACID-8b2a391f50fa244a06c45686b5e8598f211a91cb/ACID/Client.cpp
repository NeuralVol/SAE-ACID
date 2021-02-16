#include "Client.h"
#include <string.h>


// return of -1 indicates error. 

Client::Client()
{
}


Client::Client(const std::string& addr, int xport, int port)
	:m_ipAddr(addr), x_port(xport), m_port(port), out(), serverHint(), clientHint()
{
	Init();				//Initialize Winsock 
	CreateHint(x_port);	//Create hint structure for Server
}


Client::~Client()
{
	Cleanup();
}


//Start winsock
int Client::Init()
{
	WSAData data;
	WORD version = MAKEWORD(2, 2);
	int ws_ok = WSAStartup(version, &data);
	
	if (ws_ok == SOCKET_ERROR)
	{
		std::cout << "Client can't start Winsock" << WSAGetLastError()<<std::endl;
		return ws_ok;
	}
	std::cout << "Winsock successfully started" << std::endl;
	return ws_ok ;
}


// Create hint structure for Server
int Client::CreateHint(u_short port)
{	
	int sts;
	serverHint.sin_family = AF_INET;
	serverHint.sin_port = htons(port);  // little to big endian
	sts = inet_pton(AF_INET, m_ipAddr.c_str(), &serverHint.sin_addr);

	if (sts == SOCKET_ERROR || sts == 0)
	{
		std::cout << "Error converting IP address.Error: " << WSAGetLastError() << std::endl;
	}
	std::cout << "IP address converted" << std::endl;
	return sts;
}


//Bind Socket to IP address and Port
int Client::Bind(u_short port)
{
	//Client hint structure	
	clientHint.sin_family = AF_INET;
	//Client port you want to send from
	clientHint.sin_port = htons(port);  
	clientHint.sin_addr.S_un.S_addr = ADDR_ANY;

	out = socket(AF_INET, SOCK_DGRAM, 0); //Create Socket
	if (out == INVALID_SOCKET)
	{
		std::cout << "Socket not created.Error: " << WSAGetLastError() << std::endl;
	}

	if (bind(out, (sockaddr*)&clientHint, sizeof(clientHint)) == SOCKET_ERROR) 
	{
		std::cout << "Can't bind Socket " << WSAGetLastError() << std::endl;
		return -1;
	}
	return out;
}


//Send message, automatic bind and unbind socket
int Client::Send(const char *msg, size_t size)
{
	//Bind Socket
	Bind(m_port);
	int sts;
	int send_ok = sendto(out, msg, size, 0, (sockaddr*)&serverHint, sizeof(serverHint));
	if (send_ok == SOCKET_ERROR)
	{
		std::cout << "Message not sent to X-plane. Error: " << WSAGetLastError() << std::endl;
		return send_ok;
	}
	//std::cout << "Message Sent to X-Plane" << std::endl;

	//Close Socket
	sts = closesocket(out);
	if (sts == SOCKET_ERROR)
	{
		std::cout << "Socket was not closed. Error: " << WSAGetLastError() << std::endl;
	}
	return 0;

}


void Client::Cleanup()
{
	int sts;
	sts = WSACleanup();

	if (sts == SOCKET_ERROR)
	{
		std::cout << "Could not stop Winsock.Error: " << WSAGetLastError() << std::endl;
	}
}


void Client::CmndGen(char* command, const char* arg) {
	char type[] = "CMND/";
	for (int i = 0; i <= 4; i++)
		command[i] = type[i];
	for (int i = 5; i <= strlen(arg)+5; i++)
		command[i] = arg[i - 5];
}


void Client::DrefGen(char* command, const char* arg, int val) {
	char type[] = "DREF\00";
	float_converter converter_object{ val };
	for (int i = 0; i <= 4; i++)
		command[i] = type[i];
	for (int i = 5; i <= 8; i++)
		command[i] = converter_object.hex[i - 5];
	for (int i = 9; i <= strlen(arg)+9; i++)
		command[i] = arg[i - 9];
}


void Client::RrefGen(char* command, const char* arg, int freq, int val) {
	char type[] = "RREF\00";
	int_converter converter_object{ freq };
	for (int i = 0; i <= 4; i++)
		command[i] = type[i];
	for (int i = 5; i <= 8; i++)
		command[i] = converter_object.hex[i - 5];
	converter_object.val = val;
	for (int i = 9; i <= 12; i++)
		command[i] = converter_object.hex[i - 9];
	for (int i = 13; i <= strlen(arg)+13; i++)
		command[i] = arg[i - 13];
}

