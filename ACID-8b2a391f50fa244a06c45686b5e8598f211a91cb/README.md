# ACID

ACID UDP communication.

## Installation

1. Configure X-plane
```
General/Data/Output network data to Log.txt
Network/UDP port/Port we receive on 49010
Data Output/General/Send network data output
Data Output/Dataref/Network computer/Port 55555
```

2. Setup Ip address
```
ACID.h
#define XPLANE_PORT  49000  		//Port to send message to X-plane
#define SOCKET_PORT  55555  		//Port to bind socket to.
#define IP_ADRRESS  "192.168.1.3"	//Change this to one of ips in X-plane network tab
```

3. Retarget solution
```
Retarget solution vs2019
Right click on the project in the Solution Explorer
Properties/Configuration properties/General
Change "Target Platform" to Visual Studio 2019 (v141)
Change "Target Platform Version" to 10.0.17763.0
```

4. Add firewall exception
```
Allow local network/private
```

## Usage
Reference mapping: Client.h
```bash

pip X-plane.exe
pip ACID.exe
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## Authority
NeuralVol property.

## License
[MIT](https://choosealicense.com/licenses/mit/)
