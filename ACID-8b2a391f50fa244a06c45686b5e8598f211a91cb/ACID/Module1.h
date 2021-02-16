#pragma once

#include "../Vol_Client/Client.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

#define pi 3.1415926535 //I need this for the Compute functions

struct DataIn //this is the object type of data inputs from XPlane
{
	vector<double> value; //value of the data at a point in time
	int type; //the ID of the data input that the user defines in main.cpp
	void degtorad(); //function that converts degrees to radians
	void AssignData(vector<vector<double>>& dataHash); //function that takes the data from UDP and assigns it to the proper variable 
};


struct DataOut //this is the object type of data outputs
{
	vector<float> value;
	int datasource;
	float ComputeC_D(float S, float v, float rho, float alpha, float beta, float m, float a_x, float T, float  a_y, float a_z, float& Dcalc, int datasource, float D); //function that calculates C_D
	float ComputeC_SF(float S, float v, float rho, float alpha, float beta, float m, float a_x, float T, float  a_y, float a_z, float& SFcalc, int datasource, float SF); //function that calculates C_SF
	float ComputeC_L(float S, float v, float rho, float alpha, float beta, float m, float a_x, float T, float  a_y, float a_z, float& Liftcalc, int datasource, float Lift); //function that calculates C_L
	float ComputeC_l(float m, float S, float v, float rho, float b, float p_dot, float p, float q, float r_dot, float r, float I_xx, float I_yy, float I_zz, float& Lcalc, int datasource, float L); //function that calculates C_l
	float ComputeC_m(float m, float S, float v, float rho, float c, float q_dot, float p, float r, float I_xx, float I_yy, float I_zz, float& Mcalc, int datasource, float M); //function that calculates C_m
	float ComputeC_n(float m, float S, float v, float rho, float b, float p_dot, float p, float q, float r_dot, float r, float I_xx, float I_yy, float I_zz, float& Ncalc, int datasource, float N); //function that calculates C_n
};


class Dataset //this is the overall object type of module 1; it contains all the functions and objects that module 1 uses
{
	DataIn I_xx; //moment of inertia about the local x axis (kg m^2)
	DataIn I_yy; //moment of inertia about the local y axis (kg m^2)
	DataIn I_zz; //moment of inertia about the local z axis (kg m^2)
	DataIn W; //mass (kg)
	DataIn a_x; //acceleration along the local x axis (m/s^2)
	DataIn a_y; //acceleration along the local y axis (m/s^2)
	DataIn a_z; //acceleration along the local z axis (m/s^2)
	DataIn alpha; //angle of attack (rad)
	DataIn beta; //sideslip angle (rad)
	DataIn rho; //air density (kg/m^3)
	DataIn v; // airspeed (m/s)
	DataIn S; //wing area (m^2)
	DataIn b; //wingspan (m)
	DataIn c; //mean chord (m)
	DataIn p; //roll rate (rad/s)
	DataIn q; //pitch rate (rad/s)
	DataIn r; //yaw rate (rad/s)
	DataIn p_dot; //roll acceleration (rad/s^2)
	DataIn q_dot; //pitch acceleration (rad/s^2)
	DataIn r_dot; //yaw acceleration (rad/s^2)
	DataIn T; //thrust (N)
	DataIn Lift; //lift force (N)
	DataIn D; //drag force (N)
	DataIn SF; //sideforce (N)
	DataIn L; //rolling moment (Nm)
	DataIn M; //pitching moment (Nm)
	DataIn N; //yawing moment (Nm)
	
public:
	DataOut C_D; //coefficient of drag
	DataOut C_SF; //coefficient of sideforce
	DataOut C_L; //coefficient of lift
	DataOut C_l; // coefficient of rolling moment
	DataOut C_m; //coefficient of pitching moment
	DataOut C_n; //coefficient of yawing moment

	int counter; //counts the number of data inputs assigned from UDP to module 1. If this value is not exactly what the program expects, it will spit out an error when calculating coefficients.
	int inputcount; //nb of inputs we're capturing. This is initialized in the constructor but should be modified manually if we decide we need more inputs.

	void AssignVar(const char* arg, int id); //assigns the id of each input to the correct variable and increments "counter" each time this is done succesfully

	Dataset();// initializes module 1
	void Store(vector<vector<double>> &dataHash); //takes the data captured from UDP and assigns it to the variables

	void Calculate(bool debug = false); //calculates all the coefficients
};
