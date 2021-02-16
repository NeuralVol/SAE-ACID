#include <iostream>
#include <string>
#include "ACID.h"


int main()
{
	using namespace std;
	cout << "\t\t\t\t\tStart Main" << endl << endl;

	ACID obj;

	/********************************************/
	//obj.SendDref(override_control_surfaces_r, 1);
	//system("pause");
	//obj.SendDref(rudder1_deg(10), 20);

	/*******************************************/
	/*obj.SendRref(latitude_r, 1);
	obj.SendRref(longitude_r, 2);
	obj.SendRref(elevation_r, 3);*/


	obj.RequestRref();
	obj.Receive(300);
	obj.StopAll();
	obj.Process(1);

	system("pause");
	return 0;
}
