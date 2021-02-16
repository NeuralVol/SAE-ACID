#include "Module1.h"



Dataset::Dataset() 
{ //all types are initialized to 0 and inputcount is initialized based on the number of required inputs
    I_xx.type = 0;
    I_yy.type = 0;
    I_zz.type = 0;
    W.type = 0;
    a_x.type = 0;
    a_y.type = 0;
    a_z.type = 0;
    alpha.type = 0;
    beta.type = 0;
    rho.type = 0;
    v.type = 0;
    S.type = 0;
    b.type = 0;
    c.type = 0;
    p.type = 0;
    q.type = 0;
    r.type = 0;
    p_dot.type = 0;
    q_dot.type = 0;
    r_dot.type = 0;
    T.type = 0;
    Lift.type = 0;
    D.type=0;
    SF.type = 0;
    L.type = 0;
    M.type = 0;
    N.type=0; 
    counter = 0;
    inputcount = 27; //MANUALLY CHANGE THIS IF YOU ADD EXTRA INPUTS
    float x;
    int y;
    init:
    cout << "INPUT SOURCE OF DATA: 0 FOR AIRCRAFT STATE VARIABLES / 1 FOR XPLANE FORCES" <<'\n'; //populates coefficient outputs with either data directly from X-Plane or data calculated from aircraft state variables from X-Plane
    cin >> y;
    C_D.datasource = y;
    C_L.datasource = y;
    C_SF.datasource = y;
    C_l.datasource = y;
    C_m.datasource = y;
    C_n.datasource = y;
    if (!(y == 0 || y == 1))
    {
        cout << "INPUT NOT RECOGNIZED" << '\n';
            goto init; //if the user doesn't input either a 0 or 1, it will go back to the input and ask for another input
    }
    cout << "INPUT WING AREA IN m^2 (16.17 for Cessna Skyhawk)" << '\n'; //initializes wing area based on user input
    cin >> x;
    S.value.push_back(x);
    cout << "INPUT WINGSPAN IN m (11 for Cessna Skyhawk)" << '\n'; //initializes wingspan based on user input
    cin >> x;
    b.value.push_back(x);
    cout << "INPUT MEAN CHORD IN m (1.47 for Cessna Skyhawk)" << '\n'; //initializes mean chord based on user input
    cin >> x;
    c.value.push_back(x);
}

void DataIn::degtorad()
{

    for (int i = 0; i < value.size(); i++)
    {
        value[i] = value[i] * pi / 180; 
    }

}

void DataIn::AssignData(vector<vector<double>>& dataHash)
{
    if (type != 0) //type is 0 when the input required from UDP hasn't been assigned properly
    {
        for (int i = 0; i < dataHash[type].size(); i++) 
        {
            value.push_back(dataHash[type][i]);
        }
    }
}

float DataOut::ComputeC_D(float S, float v, float rho, float alpha, float beta, float m, float a_x, float T, float  a_y, float a_z, float& Dcalc, int datasource, float D)
{
    float dyn;
    float result;
    dyn = 0.5 * rho * v * v;
    if (datasource == 0)
    {
        result = 1 / (dyn * S) * ((cos(alpha) * cos(-beta)) * (m * a_x - T) - sin(-beta) * m * a_y + cos(-beta) * sin(alpha) * m * (a_z+9.81));
        Dcalc = dyn * S * result;
        return result;
    }
    else 
    {
        result = 1 / (dyn * S) * D;
        Dcalc = D;
        return result;
    }
}

float DataOut::ComputeC_SF(float S, float v, float rho, float alpha, float beta, float m, float a_x, float T, float  a_y, float a_z, float& SFcalc, int datasource, float SF)
{
    float dyn;
    float result;
    dyn = 0.5 * rho * v * v;
    if (datasource == 0)
    {
        result = 1 / (dyn * S) * ((cos(alpha) * sin(-beta)) * (m * a_x - T) + cos(-beta) * m * a_y + sin(-beta) * sin(alpha) * m * (a_z+9.81));
        SFcalc = dyn * S * result;
        return result;
    }
    else
    {
        result = 1 / (dyn * S) * SF;
        SFcalc = SF;
        return result;
    }
}
float DataOut::ComputeC_L(float S, float v, float rho, float alpha, float beta, float m, float a_x, float T, float  a_y, float a_z, float& Liftcalc, int datasource, float Lift)
{
    float dyn;
    float result;
    dyn = 0.5 * rho * v * v;
    if (datasource == 0)
    {
    result= 1 / (dyn * S) * (-sin(alpha) * (m * a_x - T)  + cos(alpha) * m * (a_z+9.81));
    Liftcalc = dyn * S * result;
    return result;
    }
    else
    {
    result = 1 / (dyn * S) * Lift;
    Liftcalc = Lift;
    return result;
    }
}

float DataOut::ComputeC_l(float m, float S, float v, float rho, float b, float p_dot, float p, float q, float r_dot, float r, float I_xx, float I_yy, float I_zz, float& Lcalc, int datasource, float L)
{
    float dyn;
    float I_xz;
    float result;
    I_xz = 0; //we can change this later if we want a different value for I_xz
    dyn = 0.5 * rho * v * v;
    if (datasource == 0)
    { 
   result= m / (dyn * S*b) * (I_xx*p_dot-I_xz*(p*q+r_dot)+(I_zz-I_yy)*q*r); //the m infront of the eqn is because the moments of inertia are per unit mass of the a/c
   Lcalc = dyn * S *b* result;
   return result;
    }
    else
    {
        result = 1 / (dyn * S*b) * L;
        Lcalc = L;
        return result;
    }
}

float DataOut::ComputeC_m(float m, float S, float v, float rho, float c, float q_dot, float p, float r, float I_xx, float I_yy, float I_zz, float& Mcalc, int datasource, float M)
{
    float dyn;
    float I_xz;
    float result;
    I_xz = 0; //we can change this later if we want a different value for I_xz
    dyn = 0.5 * rho * v * v;
    if (datasource == 0)
    { 
		result= m / (dyn * S * c) * (I_yy * q_dot + I_xz * (p * p - r*r) + (I_xx - I_zz) * p * r); //the m infront of the eqn is because the moments of inertia are per unit mass of the a/c
		Mcalc = dyn * S *c* result;
		return result;
    }
    else
    {
        result = 1 / (dyn * S*c) * M;
        Mcalc = M;
        return result;
    }
}

float DataOut::ComputeC_n(float m, float S, float v, float rho, float b, float p_dot, float p, float q, float r_dot, float r, float I_xx, float I_yy, float I_zz, float& Ncalc, int datasource, float N)
{
    float dyn;
    float I_xz;
    float result;
    I_xz = 0; //we can change this later if we want a different value for I_xz
    dyn = 0.5 * rho * v * v;
    if (datasource == 0)
    { 
    result= m / (dyn * S * b) * (I_zz * r_dot - I_xz * (p_dot - q*r) + (I_yy - I_xx) * p * q); //the m infront of the eqn is because the moments of inertia are per unit mass of the a/c
    Ncalc = dyn * S *b* result;
    return result;
    }
    else
    {
        result = 1 / (dyn * S*b) * N;
        Ncalc = N;
        return result;
    }
}

void Dataset::Store(vector<vector<double>> &dataHash)
{
    I_xx.AssignData(dataHash);
    I_yy.AssignData(dataHash);
    I_zz.AssignData(dataHash);
    W.AssignData(dataHash);
    a_x.AssignData(dataHash);
    a_y.AssignData(dataHash);
    a_z.AssignData(dataHash);
    Lift.AssignData(dataHash);
    D.AssignData(dataHash);
    SF.AssignData(dataHash);
    L.AssignData(dataHash);
    M.AssignData(dataHash);
    N.AssignData(dataHash);
    alpha.AssignData(dataHash);
    alpha.degtorad();
    beta.AssignData(dataHash);
    beta.degtorad();
    rho.AssignData(dataHash);
    v.AssignData(dataHash);
    p.AssignData(dataHash);
    q.AssignData(dataHash);
    r.AssignData(dataHash);
    p_dot.AssignData(dataHash);
    q_dot.AssignData(dataHash);
    r_dot.AssignData(dataHash);
    T.AssignData(dataHash);
    S.value.resize(p.value.size()); //making S the same size as the other vectors
    for (int i = 1; i < S.value.size(); i++)
    {
        S.value[i] = S.value[0];
    }
    b.value.resize(p.value.size()); //making b the same size as the other vectors
    for (int i = 1; i < b.value.size(); i++)
    {
        b.value[i] = b.value[0];
    }
    c.value.resize(p.value.size()); //making c the same size as the other vectors
    for (int i = 1; i < c.value.size(); i++)
    {
        c.value[i] = c.value[0];
    }
}

void Dataset::AssignVar(const char* arg, int id)
{
    using namespace std;
    if (arg == acf_Jxx_unitmass_r)
    {
        I_yy.type = id;
        counter++;
    }
    if (arg == acf_Jyy_unitmass_r)
    {
        I_zz.type = id;
        counter++;
    }
    if (arg == acf_Jzz_unitmass_r)
    {
        I_xx.type = id;
        counter++;
    }
    if (arg == m_total_r)
    {
        W.type = id;
        counter++;
    }
    if (arg == local_ax_r)
    {
        a_y.type = id;
        counter++;
    }
    if (arg == local_ay_r)
    {
        a_z.type = id;
        counter++;
    }
    if (arg == local_az_r)
    {
        a_x.type = id;
        counter++;
    }
    if (arg == alpha_r)
    {
        alpha.type = id;
        counter++;
    }
    if (arg == beta_r)
    {
        beta.type = id;
        counter++;
    }
    if (arg == rho_r)
    {
        rho.type = id;
        counter++;
    }
    if (arg == true_airspeed_r)
    {
        v.type = id;
        counter++;
    }
    if (arg == P_r)
    {
        p.type = id;
        counter++;
    }
    if (arg == Q_r)
    {
        q.type = id;
        counter++;
    }
    if (arg == R_r)
    {
        r.type = id;
        counter++;
    }
    if (arg == P_dot_r)
    {
        p_dot.type = id;
        counter++;
    }
    if (arg == Q_dot_r)
    {
        q_dot.type = id;
        counter++;
    }
    if (arg == R_dot_r)
    {
        r_dot.type = id;
        counter++;
    }
    if (arg == POINT_thrust(...))
    {
        T.type = id;
        counter++;
    }
    if (arg == SF_r)
    {
        SF.type = id;
        counter++;
    }
    if (arg == Lift_r)
    {
        Lift.type = id;
        counter++;
    }
    if (arg == D_r)
    {
        D.type = id;
        counter++;
    }
    if (arg == L_r)
    {
        L.type = id;
        counter++;
    }
    if (arg == M_r)
    {
        M.type = id;
        counter++;
    }
    if (arg == N_r)
    {
        N.type = id;
        counter++;
    }

}

void Dataset::Calculate(bool debug)
{
	//Checks to see if all the input variables have been assigned data from UDP. 
	//The minus 3 is because S, b and c are inputted manually
    if (counter == inputcount-3) { 
        vector<float> Dcalc;
        Dcalc.resize(p.value.size());
        vector<float> Liftcalc;
        Liftcalc.resize(p.value.size());
        vector<float> SFcalc;
        SFcalc.resize(p.value.size());
        vector<float> Lcalc;
        Lcalc.resize(p.value.size());
        vector<float> Mcalc;
        Mcalc.resize(p.value.size());
        vector<float> Ncalc;
        Ncalc.resize(p.value.size());

        cout << "COEFFICIENTS:" << '\n';
        for (int i = 0; i < p.value.size(); i++)
        {
            C_D.value.push_back(C_D.ComputeC_D(S.value[i], v.value[i], rho.value[i], alpha.value[i], beta.value[i], W.value[i], a_x.value[i], T.value[i], a_y.value[i], a_z.value[i], Dcalc[i], C_D.datasource, D.value[i]));
            C_SF.value.push_back(C_SF.ComputeC_SF(S.value[i], v.value[i], rho.value[i], alpha.value[i], beta.value[i], W.value[i], a_x.value[i], T.value[i], a_y.value[i], a_z.value[i], SFcalc[i], C_SF.datasource, SF.value[i]));
            C_L.value.push_back(C_L.ComputeC_L(S.value[i], v.value[i], rho.value[i], alpha.value[i], beta.value[i], W.value[i], a_x.value[i], T.value[i], a_y.value[i], a_z.value[i], Liftcalc[i], C_L.datasource, Lift.value[i]));
            C_l.value.push_back(C_l.ComputeC_l(W.value[i], S.value[i], v.value[i], rho.value[i], b.value[i], p_dot.value[i], p.value[i], q.value[i], r_dot.value[i], r.value[i], I_xx.value[i], I_yy.value[i], I_zz.value[i],Lcalc[i], C_l.datasource, L.value[i]));
            C_m.value.push_back(C_m.ComputeC_m(W.value[i], S.value[i], v.value[i], rho.value[i], c.value[i], q_dot.value[i], p.value[i], r.value[i], I_xx.value[i], I_yy.value[i], I_zz.value[i], Mcalc[i], C_m.datasource, M.value[i]));
            C_n.value.push_back(C_n.ComputeC_n(W.value[i], S.value[i], v.value[i], rho.value[i], b.value[i], p_dot.value[i], p.value[i], q.value[i], r_dot.value[i], r.value[i], I_xx.value[i], I_yy.value[i], I_zz.value[i], Ncalc[i], C_n.datasource, N.value[i]));
			if (debug) {
				cout << i << ": " << "C_D: " << C_D.value[i] << " C_SF: " << C_SF.value[i] << " C_L: " << C_L.value[i] << " C_l: " << C_l.value[i] << " C_m: " << C_m.value[i] << " C_n: " << C_n.value[i] << '\n';
				cout << "    Force values from XPlane: D: " << D.value[i] << " SF: " << SF.value[i] << " Lift: " << Lift.value[i] << " L: " << L.value[i] << " M: " << M.value[i] << " N: " << N.value[i] << '\n';
				cout << "    Calculated force values:  D: " << Dcalc[i] << " SF: " << SFcalc[i] << " Lift: " << Liftcalc[i] << " L: " << Lcalc[i] << " M: " << Mcalc[i] << " N: " << Ncalc[i] << '\n';
			}
		}
    }
    else
    {
        cout << "ERROR GETTING ALL THE DATA FROM XPLANE" << '\n';
    }
}