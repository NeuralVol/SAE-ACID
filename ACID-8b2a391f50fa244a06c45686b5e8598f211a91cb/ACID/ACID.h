#pragma once

#include "../Vol_Client/Client.h"
#include "../Vol_Server/Server.h"
#include "Module1.h"
#include <time.h>
#include <fstream>
#include <algorithm>
#include <vector>

#define XPLANE_PORT  49000  //Port to send message to X-plane
#define SOCKET_PORT  55555  //Port to bind socket to.
#define IP_ADRRESS  "192.168.1.4"

using matrix = std::vector<std::vector<double>>;
static const double TOL = 10e-7;
static int ID = 0;


class ACID
{
private:
	Client nclient;
	Server nserver;
	Dataset Input;
	std::ofstream myfile;
	std::vector< std::vector<double> > dataHash;

	char* data;
	int index;
	static const int CMD_LIMIT = 100;
	char cmdArr[CMD_LIMIT][413];
	static const int ID_RANGE = 200;
	const char* cmdName[ID_RANGE];
	bool startDebug;
	int thisId;

	template<class T, size_t row, size_t col>
	void printArray(T const (&arr)[row][col]);
	template<class T>
	void printVector(std::vector<T> const &vector);
	template<class T>
	void printMatrix(std::vector<std::vector<T>> const &mat);
	template<class T1, class T2>
	void printMatrix(std::vector<std::vector<T1>> const &mat1, std::vector<std::vector<T2>> const &mat2, int stop = -999999);
	template <class T, size_t row, size_t col>
	void transposeArray(T(&base)[row][col], T(&target)[col][row]);
	template<class T>
	void transposeVector(std::vector<std::vector<T>> &base);

	double norm(std::vector<double> const &vec);
	void scalar_div(std::vector<double> &in, double r, std::vector<double> &out);
	void scalar_sub(std::vector<double> &in, double r, std::vector<double> &out);
	double dot_product(std::vector<double> const &x, std::vector<double> const &y);
	void cross_product(double vector_a[], double vector_b[], double temp[]);
	double r_square(matrix const& z, matrix const& z_);
	
	matrix transpose(matrix const &mat);
	matrix inverse(matrix const &mat);
	void complement(matrix &mat);
	double variance(std::vector<double> const &vec);
	double least_square_cost(matrix const &base, matrix const &z);	
	double rscore(matrix const &z, matrix const &z_);

	void Polynomial(int variable, int degree, std::vector< std::vector<int> > &poly, std::vector<int> data = {});
	void TestOrthogonal(matrix const &base);
	void ModifiedGS(matrix &base, matrix &target);
	void EliminateError(matrix &base, matrix &base_temp, matrix &X, matrix &X_temp);
	void EliminateError(matrix &base, matrix &X, matrix &mapping, int count);
	void SortError(matrix const &base, matrix const &z, matrix &mapping);

public:
	ACID();
	~ACID();

	void Process(int iteration, bool debug = false);
	void EliminateDependence(matrix &origin, std::vector<std::vector<int>> &poly);
	void InitializeData(matrix &origin, matrix &z, std::vector<std::vector<int>> const &poly);
	void UpdateData(matrix &origin, matrix &z, std::vector<std::vector<int>> const &poly);
	void InitializeRandomData(matrix &origin, matrix &z, std::vector<std::vector<int>> const &poly);
	void UpdateRandomData(matrix &origin, matrix &z, std::vector<std::vector<int>> const &poly, bool noise = true);
	void PolyGen(int variable, int degree, std::vector< std::vector<int> > &poly);
	void QRGen(matrix &base, matrix &target, matrix &z);
	void MSFE(matrix &base, matrix& X, matrix const &z);
	
	void SendCmnd(const char *arg);
	void SendDref(const char *arg, int val);
	void SendRref(const char *arg, int id, int freq=20);
	void StopRref(const char *arg, int id);

	void RequestRref();
	void StopAll();
	void Receive(int iteration);
	void Record();
	void CalculateCoeff();//module1
};
