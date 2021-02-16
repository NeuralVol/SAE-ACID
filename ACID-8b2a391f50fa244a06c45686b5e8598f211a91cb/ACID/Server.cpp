#include "Server.h"



Server::Server()
{
}


Server::Server(std::string ipAddr, int port)
	: m_ipAddr(ipAddr), m_port(port), in()
{
	if (Init() < 1)
		quick_exit(0);
}


Server::~Server()
{
	Cleanup();	
}


// Initialize winsock;
int Server::Init()
{
	WSAData data;
	WORD ver = MAKEWORD(2, 2);
	int wsInit = WSAStartup(ver, &data);

	if (wsInit == SOCKET_ERROR)
	{
		std::cout << "Client can't start Winsock" << WSAGetLastError() << std::endl;
		return wsInit;
	}
	std::cout << "Winsock successfully started" << std::endl;
	return 1;
}


// Create server hint structure and bind socket to ip address
SOCKET Server::BindSocket(u_short port)
{
	sockaddr_in serverHint;
	serverHint.sin_addr.S_un.S_addr = ADDR_ANY;
	serverHint.sin_family = AF_INET;
	serverHint.sin_port = htons(port);  // little to big endian

	in = socket(AF_INET, SOCK_DGRAM, 0);
	if (in == INVALID_SOCKET)
	{
		std::cout << "Socket not created.Error: " << WSAGetLastError() << std::endl;
	}

	if (bind(in, (sockaddr*)&serverHint, sizeof(serverHint)) == SOCKET_ERROR)
	{
		std::cout << "Can't bind Socket" << WSAGetLastError() << std::endl;
		return in;
	}
	return in;
}


// Receive message automatic bind and unbind
int Server::Receive(std::vector< std::vector<double> > &dataHash, int iteration, int id)
{
	sockaddr_in client;	
	int client_len = sizeof(client);
	memset(&client, 0, client_len);

	char buf[BUFFER_SIZE];
	memset(buf, 0, BUFFER_SIZE);
	//Bind Socket
	in = BindSocket(m_port);

	float delta = 0.1;	//range of acceptance
	bool cont = true;	//flag stop if press <esc>
	long count = 0;		//frequency reduction ratio

	while (cont && count < iteration)
	{
		//wait for message
		int bytes_in = recvfrom(in, buf, 509, 0, (sockaddr*)&client, &client_len);
		if (bytes_in == SOCKET_ERROR)
		{
			std::cout << "Error receiving from client: " << WSAGetLastError << std::endl;
			continue;
		}

		//reduce radio by increase modulus
		//if (count % 2 == 0) 
		//problem doesn't copy null
		//*data = (char *)malloc(strlen(buf) + 1);
		//strcpy(*data, buf);
		//for (int i = 0; i < 50; i++)
		//	printf("%x ", buf[i]);
		//std::cout << std::endl;

		if((count+1)%100 == 0)
			std::cout << "Message received from Port: " << m_port  << " . Iteration: " << count+1 << std::endl;

		int mId_next = 0;
		float mFlt_next = 0;
		int ptr = 13;
		int receive = 0;

		memcpy(&mId_next, buf + 5, 4);
		memcpy(&mFlt_next, buf + 9, 4);

		//multiple end point
		if (id != NULL) {
			while (mId_next != 0 && mFlt_next != 0) {
				//std::cout << "Id: " << mId_next << " Float: " << mFlt_next << std::endl;
				if (mId_next == id) {
					receive = 1;
					break;
				}						
				memcpy(&mId_next, buf + ptr, 4);
				ptr += 4;
				memcpy(&mFlt_next, buf + ptr, 4);
				ptr += 4;
			}
			std::cout << "Id: " << id << "receive: " << receive << std::endl;
			closesocket(in);
			return receive;
		}

		//Loop in case combined message
		while (mId_next !=0 || mFlt_next !=0) {
			dataHash[mId_next].push_back(mFlt_next);
			//std::cout << "Id: " << mId_next << " Float: " << mFlt_next << std::endl;
			memcpy(&mId_next, buf + ptr, 4);
			ptr += 4;
			memcpy(&mFlt_next, buf + ptr, 4);
			ptr += 4;
		}

		//Display client info
		char client_ip[256];
		memset(client_ip, 0, 256);
		inet_ntop(AF_INET, &client.sin_addr, client_ip, 256);

		count++;
		if (GetAsyncKeyState(VK_ESCAPE))
			cont = false;
	}

	//Close Socket
	closesocket(in);
	wprintf(L"Exiting.\n\n");
	return 1;
}


// Cleanup socket
void Server::Cleanup()
{
	WSACleanup();
}

