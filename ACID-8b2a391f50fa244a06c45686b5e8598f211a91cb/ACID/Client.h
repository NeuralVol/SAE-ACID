#pragma once

#include <iostream>
#include <string>
#include <WS2tcpip.h>

#pragma comment(lib, "ws2_32.lib")

#define on_r "DREF\00\00\00\x80\x3f"
#define off_r "DREF\00\00\00\00\00"
#define quit_r "CMND/sim/operation/quit"

/********************************************DREF mapping**************************************************/
#define night_vision_on_r "sim/cockpit/electrical/night_vision_on\00"
#define override_throttles_r "sim/operation/override/override_throttles\00"
#define override_mixture_r "sim/operation/override/override_mixture\00"
#define override_control_surfaces_r "sim/operation/override/override_control_surfaces\00"
#define override_toe_brakes_r "sim/operation/override/override_toe_brakes\00"
#define ENGN_thro_use(...) "sim/flightmodel/engine/ENGN_thro_use["#__VA_ARGS__"]\00"
#define ENGN_mixt(...) "sim/flightmodel/engine/ENGN_mixt["#__VA_ARGS__"]\00"
#define aileron1_deg(...) "sim/flightmodel2/wing/aileron1_deg["#__VA_ARGS__"]\00"
#define aileron2_deg(...) "sim/flightmodel2/wing/aileron2_deg["#__VA_ARGS__"]"
#define elevator1_deg(...) "sim/flightmodel2/wing/elevator1_deg["#__VA_ARGS__"]\00"
#define elevator2_deg(...) "sim/flightmodel2/wing/elevator2_deg["#__VA_ARGS__"]\00"
#define rudder1_deg(...) "sim/flightmodel2/wing/rudder1_deg["#__VA_ARGS__"]\00"
#define rudder2_deg(...) "sim/flightmodel2/wing/rudder2_deg["#__VA_ARGS__"]\00"
#define flap1_deg(...) "sim/flightmodel2/wing/flap1_deg["#__VA_ARGS__"]\00" 
#define flap2_deg(...) "sim/flightmodel2/wing/flap2_deg["#__VA_ARGS__"]\00" 
#define speedbrake1_deg(...) "sim/flightmodel2/wing/speedbrake1_deg["#__VA_ARGS__"]\00"  
#define speedbrake2_deg(...) "sim/flightmodel2/wing/speedbrake2_deg["#__VA_ARGS__"]\00" 
#define spoiler1_deg(...) "sim/flightmodel2/wing/spoiler1_deg["#__VA_ARGS__"]\00" 
#define spoiler2_deg(...) "sim/flightmodel2/wing/spoiler2_deg["#__VA_ARGS__"]\00" 
#define deploy_ration(...) "sim/flightmodel2/gear/deploy_ration["#__VA_ARGS__"]\00"   
#define left_brake_ratio_r "sim/cockpit2/controls/left_brake_ratio\00"
#define right_brake_ratio_r "sim/cockpit2/controls/right_brake_ratio\00"
#define gear_handle_down_r "sim/cockpit2/controls/gear_handle_down\00"
#define parking_brake_ratio_r "sim/cockpit2/controls/parking_brake_ratio\00"


/********************************************RREF mapping**************************************************/
#define latitude_r "sim/flightmodel/position/latitude\00"
#define longitude_r "sim/flightmodel/position/longitude\00"
#define elevation_r "sim/flightmodel/position/elevation\00"
#define theta_r "sim/flightmodel/position/theta\00"
#define phi_r "sim/flightmodel/position/phi\00" //roll angle (deg)
#define psi_r "sim/flightmodel/position/psi\00"
#define alpha_r "sim/flightmodel/position/alpha\00" //angle of attack (deg)
#define beta_r "sim/flightmodel/position/beta\00" //sideslip angle (deg)
#define vpath_r "sim/flightmodel/position/vpath\00"
#define hpath_r "sim/flightmodel/position/hpath\00"
#define P_r "sim/flightmodel/position/P\00" //roll rate (deg/s)
#define Q_r "sim/flightmodel/position/Q\00" //pitch rate (deg/s)
#define R_r "sim/flightmodel/position/R\00" //yaw rate (deg/s)
#define P_dot_r "sim/flightmodel/position/P_dot\00" //deg/s^2
#define Q_dot_r "sim/flightmodel/position/Q_dot\00" //deg/s^2
#define R_dot_r "sim/flightmodel/position/R_dot\00" //deg/s^2
#define groundspeed_r "sim/flightmodel/position/groundspeed\00"
#define true_airspeed_r "sim/flightmodel/position/true_airspeed\00" //airspeed (m/s)
#define local_ax_r "sim/flightmodel/position/local_ax\00" //m/s^2
#define local_ay_r "sim/flightmodel/position/local_ay\00" //m/s^2
#define local_az_r "sim/flightmodel/position/local_az\00" //m/s^2
#define radio_altimeter_height_ft_pilot_r "sim/cockpit2/guages/indicators/radio_altimeter_height_ft_pilot\00"
#define rho_r "sim/weather/rho\00" //kg/m^3
#define acf_cgY_r "sim/aircraft/weight/acf_cgY\00"
#define acf_cgZ_r "sim/aircraft/weight/acf_cgZ\00"
#define acf_cgY_original_r "sim/aircraft/weight/acf_cgY_original\00"
#define acf_cgZ_original_r "sim/aircraft/weight/acf_cgZ_original\00"
#define acf_Jxx_unitmass_r "sim/aircraft/weight/acf_Jxx_unitmass\00" //m^2
#define acf_Jyy_unitmass_r "sim/aircraft/weight/acf_Jyy_unitmass\00" //m^2
#define acf_Jzz_unitmass_r "sim/aircraft/weight/acf_Jzz_unitmass\00" //m^2
#define m_total_r "sim/flightmodel/weight/m_total\00" //kg
#define POINT_thrust(...) "sim/flightmodel/engine/POINT_thrust["#__VA_ARGS__"]\00" //thrust (N)
#define ail1_def(...) "sim/flightmodel/controls/ail1_def["#__VA_ARGS__"]\00" //aileron deflection 1
#define ail2_def(...) "sim/flightmodel/controls/ail2_def["#__VA_ARGS__"]\00" //aileron deflection 2
#define elv1_def(...) "sim/flightmodel/controls/elv1_def["#__VA_ARGS__"]\00" //elevator deflection 1
#define elv2_def(...) "sim/flightmodel/controls/elv2_def["#__VA_ARGS__"]\00" //elevator deflection 2
#define rudd_def(...) "sim/flightmodel/controls/rudd_def["#__VA_ARGS__"]\00" //rudder deflection 1
#define rudd2_def(...) "sim/flightmodel/controls/rudd2_def["#__VA_ARGS__"]\00" //rudder deflection 2
#define fla1_def(...) "sim/flightmodel/controls/fla1_def["#__VA_ARGS__"]\00" //flaps deflection 1
#define fla2_def(...) "sim/flightmodel/controls/fla2_def["#__VA_ARGS__"]\00"// flaps deflection 2
#define splr_def(...) "sim/flightmodel/controls/splr_def["#__VA_ARGS__"]\00"
#define splr2_def(...) "sim/flightmodel/controls/splr2_def["#__VA_ARGS__"]\00"
#define SF_r "sim/flightmodel/forces/fside_aero" //N
#define Lift_r "sim/flightmodel/forces/fnrml_aero" //N
#define D_r "sim/flightmodel/forces/faxil_aero" //N
#define L_r "sim/flightmodel/forces/L_aero" //Nm
#define M_r "sim/flightmodel/forces/M_aero" //Nm
#define N_r "sim/flightmodel/forces/N_aero" //Nm


typedef union {
	float val;
	struct {
		// 32 bits IEEE, order from the LSB to the MSB. 
		unsigned int mantissa : 23;
		unsigned int exponent : 8;
		unsigned int sign : 1;
	} raw;
	char hex[sizeof(float)];
} float_converter;

typedef union {
	// 32 bits hex int, order from the LSB to the MSB.
	int val;
	char hex[sizeof(int)];
} int_converter;


class Client
{
private:
	std::string m_ipAddr;	//PC IP address
	int x_port;				//Port number to send message to X-plane
	int m_port;				//Socket port writing to
	SOCKET out;				//Socket status
	sockaddr_in serverHint;
	sockaddr_in clientHint;

	int CreateHint(u_short port);   //Create hint structure for Server
	int Bind(u_short port);			//Bind Socket to IP address and Port
	void Cleanup();					//Terminates of the Winsock

public:
	Client();
	Client(const std::string& addr, int xport, int port);

	~Client();

	int Init();  //Initialize winsock
	int Send(const char* msg, size_t size);  //Send message(DREF,RREF,CMND) to X-plane

	void CmndGen(char* command, const char* arg);						//Generate command from cmd agr
	void DrefGen(char* command, const char* arg, int val);				//Generate dref from cmd agr, value
	void RrefGen(char* command, const char* arg, int freq, int val);	//Generate dref from cmd agr, value and frequency
};
