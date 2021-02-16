#include "ACID.h"


/* ---------------- operator overloading for matrix calculation ---------------- */
matrix operator*(matrix const &m1, matrix const &m2)
{
	if (m1.at(0).size() != m2.size()) {
		std::cout << "ERROR MATRIX MULTIPLICATION" << std::endl;
		system("pause");
		std::exit(0);
	}
		
	matrix result(m1.size(), std::vector<double>(m2.at(0).size()));

	for (std::size_t row = 0; row < result.size(); ++row)
		for (std::size_t col = 0; col < result.at(0).size(); ++col)
			for (std::size_t inner = 0; inner < m2.size(); ++inner)
				result.at(row).at(col) += m1.at(row).at(inner) * m2.at(inner).at(col);

	return result;
}


matrix operator*(double r, matrix const &m1)
{
	matrix result(m1.size(), std::vector<double>(m1.at(0).size()));

	for (std::size_t row = 0; row < result.size(); ++row)
		for (std::size_t col = 0; col < result.at(0).size(); ++col)
			result.at(row).at(col) = m1.at(row).at(col)*r;

	return result;
}


matrix operator-(matrix const &m1, matrix const &m2)
{
	matrix result(m1.size(), std::vector<double>(m1.at(0).size()));

	for (std::size_t row = 0; row < result.size(); ++row)
		for (std::size_t col = 0; col < result.at(0).size(); ++col)
			result.at(row).at(col) = m1.at(row).at(col) - m2.at(row).at(col);

	return result;
}


ACID::ACID()
	: nclient(IP_ADRRESS, XPLANE_PORT, SOCKET_PORT), nserver(IP_ADRRESS, SOCKET_PORT)
{
	thisId = ++ID;
	startDebug = false;
	srand(1);//or 1, 2, time(NULL)
	index = 0;	
	dataHash.resize(ID_RANGE);
	std::cout << std::endl;
}


ACID::~ACID()
{
}


/* ---------------- printing method for ADT ---------------- */
template<class T, size_t row, size_t col>
void ACID::printArray(T const (&arr)[row][col]) {
	for (std::size_t i(0); i < row; ++i) {
		for (std::size_t j(0); j < col; ++j)
			if (abs(arr[i][j])<TOL)
				std::cout << 0 << " ";
			else
				std::cout << arr[i][j] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


template<class T>
void ACID::printVector(std::vector<T> const & vector) {
	for (auto const& e : vector)
		if (abs(e)<TOL)
			std::cout << 0 << " ";
		else
			std::cout << e << " ";
	std::cout << std::endl << std::endl;
}


template<class T>
void ACID::printMatrix(std::vector<std::vector<T>> const &mat) {
	for (std::vector<T> row : mat) {
		for (T val : row)
			if (abs(val)<TOL)
				std::cout << 0 << " ";
			else
				std::cout << val << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


template<class T1, class T2>
void ACID::printMatrix(std::vector<std::vector<T1>> const &mat1, std::vector<std::vector<T2>> const &mat2, int stop) {
	for (int i = 0; i < mat1.size(); ++i) {
		if (mat1[i][0] == stop)
			continue;
		std::cout << "Row: " << i << "\t";
		for (int j = 0; j < mat1[0].size(); ++j)
			std::cout << mat1[i][j] << " ";
		std::cout << "\t";
		for (int j = 0; j < ((mat2[0].size()>20) ? 20 : mat2[0].size()); ++j)
			std::cout << mat2[i][j] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


template <class T, size_t row, size_t col>
void ACID::transposeArray(T (&base)[row][col], T (&target)[col][row]) {
	for (int i = 0; i < row; ++i)
		for (int j = 0; j < col; ++j)
			target[j][i] = base[i][j];
}


// Inplace transpose vector
template<class T>
void ACID::transposeVector(std::vector<std::vector<T>>& base) {
	//problem here resize just fill new element with default value the old one remain the same
	//std::vector< std::vector<double> > target(base.size(), vector<double>(base.size(), 0));		

	std::vector<std::vector<T>> temp = base;
	base.clear();
	base.resize(temp[0].size(), std::vector<T>(temp.size(), 0));

	for (int i = 0; i < temp.size(); ++i)
		for (int j = 0; j < temp[0].size(); ++j)
			base[j][i] = temp[i][j];
}


// Function to calculate the 2-norm of the vectors
double ACID::norm(std::vector<double> const &vec) {
	double sum = 0;

	for (int i = 0; i < vec.size(); ++i)
		sum += vec[i] * vec[i];

	return sqrt(sum);
}


// Function to divide 1d vector by a scalar
void ACID::scalar_div(std::vector<double> &in, double r, std::vector<double> &out) {
	for (int i = 0; i < in.size(); i++)
		out[i] = in[i] / r;
}


// Function to subtract 1d vector by a scalar
void ACID::scalar_sub(std::vector<double> &in, double r, std::vector<double> &out) {
	for (int i = 0; i < in.size(); i++)
		out[i] -= r * in[i];
}


// Function to calculate dot product of two vectors
double ACID::dot_product(std::vector<double> const &x, std::vector<double> const &y) {
	double sum = 0;

	for (int i = 0; i < x.size(); i++)
		sum += x[i] * y[i];

	return sum;
}


// Function to calculate cross product of two vectors
void ACID::cross_product(double vector_a[], double vector_b[], double temp[]) {
	temp[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
	temp[1] = vector_a[0] * vector_b[2] - vector_a[2] * vector_b[0];
	temp[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];
}


// Create transpose matrix
matrix ACID::transpose(matrix const &mat) {
	matrix transpose;
	transpose.resize(mat[0].size(), std::vector<double>(mat.size(), 0));

	for (int i = 0; i < mat.size(); ++i)
		for (int j = 0; j < mat[0].size(); ++j)
			transpose[j][i] = mat[i][j];

	return transpose;
}


// Inverse matrix method
// Gaussian Elimination combined indentity matrix
// Determinant is very slow and complicated for large matrix
matrix ACID::inverse(matrix const &mat)
{
	//problem here doesn't change the size of current vector
	//A.resize(temp-1, std::vector<double>(temp+1, 0));
	//vector<double> line(2 * n, 0);
	//vector< vector<double> > A(n, line);

	std::vector< std::vector<double> > A = mat;
	std::vector< std::vector<double> > inv = mat;
	int size = mat.size();
	for (int i = 0; i < size; i++)
		A[i].resize(size * 2, 0);
	for (int i = 0; i < size; i++)
		A[i][size + i] = 1;
	int n = A.size();

	for (int i = 0; i < n; i++) {
		// Search for maximum in this column
		double maxEl = abs(A[i][i]);
		int maxRow = i;
		for (int k = i + 1; k < n; k++) {
			if (abs(A[k][i]) > maxEl) {
				maxEl = A[k][i];
				maxRow = k;
			}
		}

		// Swap maximum row with current row (column by column)
		for (int k = i; k < 2 * n; k++) {
			double tmp = A[maxRow][k];
			A[maxRow][k] = A[i][k];
			A[i][k] = tmp;
		}

		// Make all rows below this one 0 in current column
		for (int k = i + 1; k < n; k++) {
			double c = -A[k][i] / A[i][i];
			for (int j = i; j < 2 * n; j++) {
				if (i == j) {
					A[k][j] = 0;
				}
				else {
					A[k][j] += c * A[i][j];
				}
			}
		}
	}

	// Solve equation Ax=b for an upper triangular matrix A
	for (int i = n - 1; i >= 0; i--) {
		for (int k = n; k < 2 * n; k++) {
			A[i][k] /= A[i][i];
		}
		A[i][i] = 1;//not necessary

		for (int rowModify = i - 1; rowModify >= 0; rowModify--) {
			for (int columModify = n; columModify < 2 * n; columModify++) {
				A[rowModify][columModify] -= A[i][columModify] * A[rowModify][i];
			}
			A[rowModify][i] = 0;//not necessary
		}
	}

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			inv[i][j] = A[i][j + size];

	return inv;
}


void ACID::complement(matrix &mat)
{
	for (int i = 0; i < mat.size(); ++i)
		for (int j = 0; j < mat[0].size(); ++j)
			if(mat[i][j] > 0)
				mat[i][j] = -mat[i][j];
}


// Get variance of zi
double ACID::variance(std::vector<double> const &vec)
{
	double var = 0;
	double N = vec.size();
	double mean = 0;

	for (int i = 0; i < N; i++)
		mean += vec[i];
	mean = mean / N;

	for (int i = 0; i < N; i++)
		var += (vec[i] - mean)*(vec[i] - mean);
	var = var / (N - 1);

	return var;
}


// Get least_square_cost support for a_ calculation
double ACID::least_square_cost(matrix const &base, matrix const &z)
{
	double cost = 0;
	double element = 0;
	double numerator = 0;
	double denominator = 0;
	matrix mat = transpose(base);
	matrix trans;
	for (int i = 0; i < mat.size(); i++) {
		trans.clear();
		trans.resize(1, mat[i]);
		numerator = (trans*z)[0][0];
		denominator = (trans*transpose(trans))[0][0];
		element = numerator * numerator / denominator;
		cost += element;
	}		
	return cost;
}


// Coefficient of determination R^2
double ACID::rscore(matrix const &z, matrix const &z_)
{
	double rscore = 0;
	double mean = 0;
	double numerator = 0;
	double denominator = 0;
	double N = z.size();
	double temp;
	matrix dif = z;

	//mean of z'
	for (int i = 0; i < N; i++)
		mean += z_[i][0];
	mean = mean / N;

	for (int i = 0; i < N; i++) {
		temp = z[i][0] - z_[i][0];
		numerator += temp * temp;
	}

	for (int i = 0; i < N; i++) {
		temp = z[i][0] - mean;
		denominator += temp * temp;
	}

	//numerator = (transpose(dif)*dif)[0][0];
	//denominator = (transpose(z)*z)[0][0] - N*mean*mean;
	std::cout << "numerator: " << numerator << "denominator: " << denominator << std::endl << std::endl;

	rscore = 1 - (numerator / denominator);
	return rscore;
}


// Coefficient of determination R^2
double ACID::r_square(matrix const &z, matrix const &z_)
{
	double rscore = 0;
	double mean = 0;
	double numerator = 0;
	double denominator = 0;
	double N = z.size();
	matrix dif = z;

	for (int i = 0; i < N; i++)
		dif[i][0] -= z_[i][0];

	for (int i = 0; i < N; i++)
		mean += z_[i][0];
	mean = mean / N;

	numerator = (transpose(dif)*dif)[0][0];
	denominator = (transpose(z)*z)[0][0] - N*mean*mean;

	rscore = 1 - (numerator / denominator);
	return rscore;
}


/* ---------------- generate all possible Monomials ---------------- */
/*  Solve k-combination with repetitions, or k-multicombination, or 
multisubset of size k from a set S is given by a sequence of k not 
necessarily distinct elements of S, where order is not taken into 
account. Resolve by twelve-fold way technique. 
Result	: from n choose n + k -1 combinations
curDeg	: current degreee of monomial

Input variables:
poly	: vector contain result monomials
data	: container for testing contain one path of recursion
variable	: number of columns in Base.
degree	: number of rows in Base.									*/
void ACID::Polynomial(int variable, int degree, std::vector< std::vector<int> > &poly, std::vector<int> data)
{
	int curDeg;
	if (variable == 1)	{
		data.push_back(degree);
		poly.push_back(data);
	}
	else
		for (curDeg = 0; curDeg <= degree; curDeg++) {
			std::vector<int> vec;
			vec.insert(vec.end(), data.begin(), data.end());
			vec.push_back(degree - curDeg);
			Polynomial(variable - 1, curDeg, poly, vec);
		}
}


// Test orthonormal set
void ACID::TestOrthogonal(matrix const &base) {
	printf("Numerical verification that {q_1, ..., q_%i} is an "
		"orthonormal set:\n", base.size());
	for (int i = 0; i < base.size(); i++) {
		for (int j = i; j < base.size(); j++) {
			double x = dot_product(base[i], base[j]);
			printf("q_%i * q_%i = %lg\n", i + 1, j + 1, x);
		}
	}
	printf("\n");
}


/* ------------------ modified gramSchmidt --------------------	*/
/*  Given a matrix A of dimension m by n, this algorithm
computes a QR decomposition of A, where Q is a unitary
m by n matrix and R is a n by n upper triangular matrix
and A = QR -> QR factorization.
This implement base on the nature of C++ array, each 
vector is stored as sub-array for the ease of access
Transpose target automatic after finish transformation

Input variables:
Base	: pointer to array of arrays, the ith array of
which should correspond to the ith column of the
matrix Base. During the algorithm, the columns of Q
will replace the columns of Base.
Target	: pointer to array of arrays in which the ith
column of the upper triangular matrix R will be
stored in the ith subarray of Target.							*/
void ACID::ModifiedGS(matrix &base, matrix &target) 
{
	if (base.size() > base.at(0).size()) {
		std::cout << "ERROR No sample < No variable" << std::endl;
		system("pause");
		std::exit(0);
	}

	target.clear();
	target.resize(base.size(), std::vector<double>(base.size(), 0));
	int i, j;
	double anorm = 10e-7;

	for (i = 0; i < base.size(); i++) {
		target[i][i] = norm(base[i]);                  // r_ii = ||a_i||

		if (target[i][i] > TOL) {
			scalar_div(base[i], target[i][i], base[i]);   // a_i = a_i/r_ii
		}
		else if (i == 0) { // set base[0] = [1 0 0 ... 0]^T
			base[i][0] = 1;
			for (j = 1; j < base[0].size(); j++) {
				base[i][j] = 0;
			}
		}
		else { // need to choose a_i orthogonal to < a_1, ... a_{i-1} >
			for (j = 0; j < base[0].size(); j++) {
				base[i][j] = -base[0][i] * base[0][j];
			}
			base[i][i] += 1;

			for (j = 1; j < i; j++) {
				scalar_sub(base[j], base[j][i], base[i]);
			}

			anorm = norm(base[i]);
			scalar_div(base[i], anorm, base[i]);
		}

		for (j = i + 1; j < base.size(); j++) {
			target[j][i] = dot_product(base[i], base[j]); // r_ij = a_i*a_j
			scalar_sub(base[i], target[j][i], base[j]);   // a_j -= r_ij a_i
		}
	}

	bool done = true;
	for (int i = 0; i < target.size(); i++)
		if (abs(target[i][i]) < TOL)
			done = false;

	if (done) {
		for (int a = 0; a < base.size(); a++)
			for (int b = 0; b < a; b++) {
				double tmp = target[a][b];
				target[a][b] = target[b][a];
				target[b][a] = tmp;
			}
	}		

	//Remarks: Problem if A is nearly rank-deficient, P can't be orthogonal
	/*for (; i < base[0].size(); i++) {
		for (j = 0; j < base[0].size(); j++) {
			base[i][j] = -base[0][i] * base[0][j];
		}
		base[i][i] += 1;

		for (j = 1; j < i; j++) {
			scalar_sub(base[j], base[j][i], base[i]);
		}

		anorm = norm(base[i]);
		scalar_div(base[i], anorm, base[i]);
	}*/
}


// Eliminate pi to minimize the sum error in case there is no order
void ACID::EliminateError(matrix &base, matrix &base_temp, matrix &X, matrix &X_temp)
{
	transposeVector(base);
	transposeVector(X);

	//remove insignificant p
	bool remove;
	for (int i = 0; i < base.size(); i++) {
		remove = true;
		for (int j = 0; j < base_temp.size(); j++)
			if (base[i] == base_temp[j]) {
				remove = false;
				break;
			}
		if (remove)
			base[i] = std::vector<double>(base[0].size(), 0);
	}

	//remove insignificant x
	for (int i = 0; i < X.size(); i++) {
		remove = true;
		for (int j = 0; j < X_temp.size(); j++)
			if (X[i] == X_temp[j]) {
				remove = false;
				break;
			}
		if (remove)
			X[i] = std::vector<double>(X[0].size(), 0);
	}

	transposeVector(base);
	transposeVector(X);

	if (startDebug) {
		printMatrix(base);
		printMatrix(X);
	}
}


// Eliminate pi to minimize the sum error bottom up - order set
void ACID::EliminateError(matrix &base, matrix &X, matrix &mapping, int count)
{
	transposeVector(base);
	transposeVector(X);

	//base.resize(base.size() - count);
	//base.resize(base.size() + count, std::vector<double>(base[0].size(), 0));
	//X.resize(X.size() - count);
	//X.resize(X.size() + count, std::vector<double>(X[0].size(), 0));

	double last_index = mapping[1].size() - 1;
	for (int i = 0; i < count; i++) {
		base[mapping[1][last_index - i]] = std::vector<double>(base[0].size(), 0);
		X[mapping[1][last_index - i]] = std::vector<double>(X[0].size(), 0);
	}

	//printMatrix(X);

	transposeVector(base);
	transposeVector(X);

	if (startDebug) {
		printMatrix(base);
		printMatrix(X);
	}
}


// Sort Q base on the signification of error in J
//std::sort(rank.begin(), rank.end());
//problem even sort after QR, re QR don't produce the same result
//transform R is difficult as it need to be square matrix, need to re-orthogonalize
void ACID::SortError(matrix const &base, matrix const &z, matrix &mapping)
{
	double ele;
	double cost;
	matrix temp = base;
	transposeVector(temp);
	matrix pj_trans(1, temp[0]);	
	mapping.clear();
	mapping.resize(3, base[0]);
	transposeVector(mapping);
	
	for (int i = 0; i < temp.size(); i++) {
		pj_trans.clear();
		pj_trans.push_back(temp[i]);
		ele = (pj_trans*z)[0][0];
		cost = ele * ele;
		temp[i].push_back(cost);
		mapping[i][0] = cost;
		mapping[i][1] = i;
	}

	//sort with customized iterator
	std::sort(mapping.begin(), mapping.end(),
		[](const std::vector<double> &lhs, const std::vector<double> &rhs)
			{return lhs[0] > rhs[0];});

	for (int i = 0; i < temp.size(); i++)
		mapping[i][2] = i;

	cout << "\t\t\t\tSortError" << endl << endl;
	transposeVector(mapping);
	printMatrix(mapping);
	//remove row contain ranking
	transposeVector(temp);
	temp.resize(temp.size() - 1);
}


/*------------------ Regression Process ------------------*/
/* Important: n ≤ m ≤ nc (4)
This proccess only find m, need another filter to get n
Feature:
Esc			: exit the loop
Enter		: update data for another iteration
Debug		: turn on debug mode if true
Rscore		: >0.9 give reliable interpolation
Correlation	: the comparison between z and z_.			*/
//calculate theta before MSFE create loss of precision
//base = X;
//QRGen(base, target, X);
//theta = inverse(target)*transpose(base)*z; 
//MSFE(base, z, X);
void ACID::Process(int iteration, bool debug)
{
	using namespace std;
	startDebug = debug;
	matrix origin;
	matrix X;
	matrix base;
	matrix target;
	matrix theta;
	matrix z;
	matrix z_;
	vector<vector<int>> poly;					//10 & 4 -> 1001 combinations
	int variable = 10;							//parameter in 12
	int degree = 3;								//maximum degree 
	PolyGen(variable, degree, poly);			//generate monomials
	InitializeData(origin, z, poly);			//construct initial database
	EliminateDependence(origin, poly);			//remove dependent vector

	int count = 0;
	double r_score = 0;
	int limit = iteration;
	while (count < limit)
	{
		UpdateData(origin, z, poly);			//append incoming sample
		X = origin;								//regressor matrix
		base = origin;							//X=base*target
		QRGen(base, target, z);					//QR factorization
		MSFE(base, z, X);						//minimize predicted squared error

		theta = inverse(target)*transpose(base)*z;
		z_ = X*theta;
		r_score = r_square(z, z_);
		if(r_score>0.85 && r_score < 1)
			count = limit - 1;

		if (count == (limit-1)) {
			cout << "***					Coefficient					***" << endl << endl;
			printMatrix(theta, transpose(X), 0);
			cout << "***					Correlation					***" << endl << endl;
			printMatrix(z, z_);
		}

		cout << "R score: " << r_score << endl << endl;
		if (GetAsyncKeyState(VK_ESCAPE) && (count != limit - 1))
			count = limit-1;
		else 
			count++;
	}
}


// Eliminate some p that is dependent vector
void ACID::EliminateDependence(matrix &origin, std::vector<std::vector<int>> &poly)
{
	cout << "origin: " << origin.size() << " " << origin[0].size() << endl;
	cout << "poly: " << poly.size() << " " << poly[0].size() << endl;
	transposeVector(origin);

	int dim = 0;
	matrix target;
	matrix temp = origin;

	ModifiedGS(temp, target);

	for (int i = 0; i < target.size(); i++)
		if (abs(target[i][i]) < TOL)
		{
			origin.erase(origin.begin() + dim);
			poly.erase(poly.begin() + dim);
		}
		else
			dim++;

	if (startDebug)
		TestOrthogonal(temp);

	transposeVector(origin);
	cout << "origin: " << origin.size() << " " << origin[0].size() << endl;
	cout << "poly: " << poly.size() << " " << poly[0].size() << endl << endl;
}


/* airspeed 6, pitch rate 4, angle of attack 1, sideslip 2, roll rate 3, yaw rate 5, 
thrust 25, roll angle 26, elevator angle 27, rudder angle 28, aileron angle 27, flap angle 30,
Matrix dataX containt 0-5 and 24-29 rows of dataHash (remove the first unused row)*/
void ACID::InitializeData(matrix &origin, matrix &z, std::vector<std::vector<int>> const &poly)
{
	using namespace std;
	//matrix contain parameters for polynomial
	matrix dataX = dataHash;
	dataX.erase(dataX.begin());
	dataX.resize(30);
	for(int i = 0; i < 18; i++)
		dataX.erase(dataX.begin()+6);
	transposeVector(dataX);
	complement(dataX);

	int temp = 0;
	//cout << "Display data set? (0-no; 1-dataX; 2-monomials)\n";
	//cin >> temp;
	//if (temp == 1) {
	//	cout << "***\t\t\tdataX***" << endl << endl;
	//	printMatrix(dataX);
	//}

	vector<double> sample;
	double factor;
	//generate random monomial sample
	for (int i = 0; i < dataX.size(); i++) {
		sample.clear();
		for (int j = 0; j < poly.size(); j++) {
			factor = 1;
			for (int k = 0; k < poly[0].size(); k++)
				factor *= pow(dataX[i][k], poly[j][k]);
			sample.push_back(factor);
		}
		origin.push_back(sample);
	}

	//sample vector z
	z.resize(origin.size(), std::vector<double>(1));
	//generate random z
	for (int i = 0; i < z.size(); i++)
		z[i][0] = Input.C_SF.value[i];

	cout << "***					Sample Data					***" << endl << endl;
	if (temp == 2)
		printMatrix(origin);	
	cout << "Initial Sample: " << origin.size() << endl << endl;
}


void ACID::UpdateData(matrix &origin, matrix &z, std::vector<std::vector<int>> const &poly)
{
	std::vector<double> vecX;
	for (int i = 0; i < 30; i++)
		vecX.push_back(dataHash[i + 1][dataHash[1].size() - 1]);
	for (int i = 0; i < 18; i++)
		vecX.erase(vecX.begin() + 6);

	std::vector<double> sample;
	double factor;
	for (int j = 0; j < poly.size(); j++) {
		factor = 1;
		for (int k = 0; k < poly[0].size(); k++)
			factor *= pow(vecX[k], poly[j][k]);
		sample.push_back(factor);
	}
	origin.push_back(sample);

	// Update last sample
	//int r = origin.size() - 1;
	int r = dataHash[1].size() - 1;
	z.push_back(std::vector<double>(1, 0));
	z[z.size() - 1][0] = Input.C_SF.value[r] * (1 + (rand() % 10) / 100);

	cout << "Parameter: AOA-Sideslip-Roll-Pitch-Yaw-Airspeed-Thrust-Rollangle-Elevatorangle-Rudderangle-Aileronangle-Flagangle" << endl;
	printVector(vecX);
	cout << "Origin: " << origin.size() << "x" << origin[0].size() << endl << endl;
}


void ACID::InitializeRandomData(matrix &origin, matrix &z, std::vector<std::vector<int>> const &poly)
{
	using namespace std;
	matrix dataX = {
		/*{ 1, 	3, 	8, 	15, 24 },
		{ 3, 	4, 	5, 	12, 21 },
		{ 8, 	5, 	9, 	7, 	16 },
		{ 15, 	12,	7, 	16, 9 },
		{ 24, 	21, 16, 9, 	25 },
		{ 35,	24, 19, 13, 8 },
		{ 32,	45, 11, 17, 5 },
		{ 27,	40, 55, 15, 9 },
		{ 20,	33, 48, 26,	7 },
		{ 11,	24, 39, 28, 16 },
		{ 36,	13, 28, 45, 27 },
		{ 13,	49, 15, 32, 41 },
		{ 28,	15, 23, 17, 55 },
		{ 31,	19, 33, 59, 47 },
		{ 22,	32, 17, 41, 61 },
		{ 31,	51, 36, 19, 57 },
		{ 29,	19, 43, 51, 63 },
		{ 35, 	32, 27, 20, 11 },
		{ 25, 	45, 40, 33, 24 },
		{ 13, 	27, 21, 48, 39 },
		{ 41, 	35, 7,  53, 47 },
		{ 45, 	42, 36, 27, 29 },
		{ 36, 	41, 32, 43, 49 },
		{ 31, 	29, 13, 37, 23 },
		{ 53, 	17, 39, 21, 41 } };*/
		{ 1, 	3,  31, 29, 13, 37, 23, 53, 17, 39, 21, 41 } };

	vector<double> sample;
	double factor;
	//generate random monomial sample
	for (int i = 0; i < dataX.size(); i++) {
		sample.clear();
		for (int j = 0; j < poly.size(); j++) {
			factor = 1;
			for (int k = 0; k < poly[0].size(); k++)
				factor *= pow(dataX[i][k], poly[j][k]);
			sample.push_back(factor);
		}
		origin.push_back(sample);
	}

	//sample vector z
	z.resize(origin.size(), std::vector<double>(1)); 
	//generate random z
	for (int i = 0; i < z.size(); i++) 
		z[i][0] = origin[i][0]*2 + origin[i][1]*3 + origin[i][2]*5 
		+ origin[i][3]*7 + origin[i][5]*(rand()%3);
	
	//at least 200 initial sample
	int limit = (poly.size() > 200) ? (ceil(poly.size()/100.0) * 100.0) : 200;
	int iteration = limit - origin.size();
	for (int i = 0; i < iteration; i++)
		UpdateRandomData(origin, z, poly, false);

	cout << "***					Sample Data					***" << endl << endl;
	printMatrix(origin);
	cout << "Initial Sample: " <<origin.size() << endl << endl;
}


void ACID::UpdateRandomData(matrix &origin, matrix &z, std::vector<std::vector<int>> const &poly, bool noise)
{
	// Generate random numbers by lambda func and fill it in vector
	std::vector<double> vecX(origin[0].size());
	std::generate(vecX.begin(), vecX.end(), []() {
		return (rand() % 25);
	});

	std::vector<double> sample;
	double factor;
	for (int j = 0; j < poly.size(); j++) {
		factor = 1;
		for (int k = 0; k < poly[0].size(); k++)
			factor *= pow(vecX[k], poly[j][k]);
		sample.push_back(factor);
	}
	origin.push_back(sample);

	// Update last sample
	int r = origin.size() - 1;
	double ratio = (noise) ? (rand() % 2 + 0.5) : 1;
	z.push_back(std::vector<double>(1, origin[r][0] * 2 + origin[r][1] * 3 
		+ origin[r][2] * 5 + origin[r][3] * 7 + origin[r][5]*ratio));
}


// Generate all possible monomial (decrease degree)
void ACID::PolyGen(int variable, int degree, std::vector< std::vector<int> > &poly)
{
	using namespace std;
	
	for (int i = 0; i < degree; i++)
		Polynomial(variable, degree - i, poly);
	poly.push_back(std::vector<int>(poly[0].size(), 0));

	cout << "***					Monomials					***" << endl << endl;
	for (const auto &row : poly) {
		cout << " ";
		for (const auto &element : row) {
			cout << element << " ";
		}
		cout << endl;
	}
	cout << " Total: " << poly.size() << endl << endl;
}


// QR factorization
//col =< row and rank should be equal col to avoid rank-deficient
void ACID::QRGen(matrix &base, matrix &target, matrix &z)
{
	matrix target_temp;
	/*matrix X = {
		{ -1,	-1,	1, -1 },
		{ 1,	3,	3, 1 },
		{ -1,	-1, 5, 3 },
		{ 1,	3,	7, 5 },
		{ -1,	-1,	9, 15 } };
		{ 1, 	3, 	8, 	15, 24, 35, 48, 63 },
		{ 3, 	4, 	5, 	12, 21, 32, 45, 60 },
		{ 8, 	5, 	9, 	7, 	16, 27, 40, 55 },
		{ 15, 	12,	7, 	16, 9, 	20, 33, 48 },
		{ 24, 	21, 16, 9, 	25, 11, 24, 39 },
		{ 35, 	32, 27, 20, 11, 36, 13, 28 },
		{ 48, 	45, 40, 33, 24, 13, 49, 15 },
		{ 63, 	60, 55, 48, 39, 28, 15, 64 },
		{ 80, 	77, 72, 65, 56, 45, 32, 17 },
		{ 99, 	96, 91, 84, 75, 64, 51, 36 } };
		{ 1, 	3, 	8, 	15, 24, 35, 48, 63, 80, 8 },
		{ 3, 	4, 	5, 	12, 21, 32, 45, 60, 77, 5 },
		{ 8, 	5, 	9, 	7, 	16, 27, 40, 55, 72, 9 },
		{ 15, 	12,	7, 	16, 9, 	20, 33, 48, 65,	7 },
		{ 24, 	21, 16, 9, 	25, 11, 24, 39, 56, 16 },
		{ 35, 	32, 27, 20, 11, 36, 13, 28, 45, 27 },
		{ 48, 	45, 40, 33, 24, 13, 49, 15, 32, 41 },
		{ 63, 	60, 55, 48, 39, 28, 15, 64, 17, 55 },
		{ 71, 	67, 69, 53, 47, 31, 19, 33, 59, 47 },
		{ 80, 	77, 72, 65, 56, 45, 32, 17, 81, 72 },
		{ 99, 	96, 91, 84, 75, 64, 51, 36, 19, 91 },
		{ 95, 	29, 13, 97, 61, 59, 19, 43, 71, 83 } };*/	

	transposeVector(base);
	//QR factorization
	ModifiedGS(base, target); 
	transposeVector(base);

	if (startDebug) {
		printMatrix(base);
		printMatrix(target);
	}
}


// Mean squared fit error
//lost precision due to backward subtistution only eliminate in the back
//o_max estimate by divide 25
//matrix base_temp = base;
//matrix X_temp = X;
//SortError(base_temp, z, X_temp);
void ACID::MSFE(matrix &base, matrix &X, matrix const &z)
{
	using namespace std;
	matrix a_; //coefficient estimation vector â
	matrix esp; //error vector espsilon
	double J = 0; //theory vector J based on exact a Eq.12
	double J_ = 0; //estimation vector J Eq.20

	matrix trans; //transpose of base
	double N; //sample times
	double n; //number of qj
	double var = 0;
	double PSE = 0;

	cout << "***					MSFE					***" << endl << endl;
	//sort qj to achieve global minimun
	matrix mapping;
	SortError(base, z, mapping);
	matrix base_temp = base;
	matrix X_temp = X;

	trans = transpose(base_temp);
	N = base_temp.size();
	n = base_temp[0].size();

	a_ = trans*z; 
	//a_ = inverse(trans*base)*trans*z;  //aj=pj-t.z
	//J = 0.5*(transpose(z)*z - transpose(a_)*a_)[0][0];
	esp = z - transpose(trans) * a_;
	J = 0.5 * (transpose(esp) * esp)[0][0];

	var = variance(transpose(z)[0]);
	PSE = (2*J + var*n) / N;
	cout << "PSE: " << PSE << " No: " << n << " N: " << N << endl << endl;

	int count = 0;
	double prePSE = PSE;	
	int total = base_temp[0].size();
	// Set limit removal
	int limit = total*0.95;
	transposeVector(base_temp);
	transposeVector(X_temp);
	
	for(int i = 1; i < limit; i++)
	{	
		double last_index = mapping[1].size();
		base_temp[mapping[1][last_index - i]] = std::vector<double>(base_temp[0].size(), 0);
		X_temp[mapping[1][last_index - i]] = std::vector<double>(X_temp[0].size(), 0);
		//matrix test(base_temp.size(), std::vector<double>(1, 0));
		//printMatrix(base_temp, test);

		trans = base_temp;
		N = base_temp[0].size();
		n = base_temp.size() - i;

		a_ = trans*z;
		esp = z - transpose(trans)*a_;
		J = 0.5*(transpose(esp)*esp)[0][0];

		PSE = (2*J + var*n) / N;
		cout << "J: " << J << " var*n: " << var*n << " No: " << n << endl << endl;
		//if (PSE < prePSE) { // incorrect result
		if (abs(PSE) < abs(prePSE)) {
			count++;
			prePSE = PSE;
		}
		else
			break;
	}

	/***************************************************************/
	//J_ = 0.5*((transpose(z)*z)[0][0] - least_square_cost(base, z));
	//cout << "J_: " << J_ << endl << endl;
	//transposeVector(base_temp);
	//transposeVector(X_temp);

	//Wrong way eliminate exceed 1 vector
	//EliminateError(base, base_temp, X, X_temp);
	//count = (count > limit) ? limit : count;

	EliminateError(base, X, mapping, count);
	cout << "PSE min: " << prePSE << endl << "Pj remove: " << count << endl << endl;
}


void ACID::SendCmnd(const char* arg) {
	char command[509];
	nclient.CmndGen(command, arg);
	nclient.Send(command, sizeof(command));
}


void ACID::SendDref(const char* arg, int val) {
	char command[509];
	nclient.DrefGen(command, arg, val);
	nclient.Send(command, sizeof(command));
}


void ACID::SendRref(const char* arg, int id, int freq) {
	char command[413];
	nclient.RrefGen(command, arg, freq, id);
	nclient.Send(command, sizeof(command));
	//while (nserver.Receive(dataHash, id) == 0)
	//	nclient.Send(command, sizeof(command));
	nclient.RrefGen(cmdArr[index++], arg, 0, id);
	cmdName[id] = arg;
	Input.AssignVar(arg, id); //module1
}


void ACID::StopRref(const char* arg, int id) {
	char command[413];
	nclient.RrefGen(command, arg, 0, id);
	nclient.Send(command, sizeof(command));
}


void ACID::RequestRref()
{
	SendRref(alpha_r, 1); //angle of attack
	SendRref(beta_r, 2); //sideslip angle
	SendRref(P_r, 3); //roll rate
	SendRref(Q_r, 4); //pitch rate
	SendRref(R_r, 5); //yaw rate
	SendRref(true_airspeed_r, 6); //airspeed
	SendRref(acf_Jxx_unitmass_r, 7); // 2nd moment of inertia about the x axis
	SendRref(acf_Jyy_unitmass_r, 8); // 2nd moment of inertia about the y axis
	SendRref(acf_Jzz_unitmass_r, 9); // 2nd moment of inertia about the z axis
	SendRref(P_dot_r, 10); //roll acceleration
	SendRref(Q_dot_r, 11); //pitch acceleration
	SendRref(R_dot_r, 12); //yaw acceleration
	SendRref(local_ax_r, 13); //acceleration in the x direction
	SendRref(local_ay_r, 14); //acceleration in the y direction
	SendRref(local_az_r, 15); //acceleration in the z direction
	SendRref(POINT_thrust(...), 16); //thrust
	SendRref(m_total_r, 17); //total aircraft mass
	SendRref(Lift_r, 18); //lift force
	SendRref(D_r, 19); // drag force
	SendRref(SF_r, 20); //sideforce
	SendRref(L_r, 21); //moment about the x axis / rolling moment
	SendRref(M_r, 22); //moment about the y axis / pitching moment
	SendRref(N_r, 23); //moment about the z axis / yawing moment
	SendRref(rho_r, 24); //air density
	SendRref(POINT_thrust(0), 25); //thrust
	SendRref(phi_r, 26); //roll angle	
	SendRref(elv1_def(8), 27); //elevator angle
	SendRref(rudd_def(10), 28); //rudder angle
	SendRref(ail1_def(0), 29); //aileron angle
	SendRref(fla1_def(0), 30); //flap angle
}


//for (int j = 0; j < 50; j++)
//	printf("%x ", cmdArr[i][j]);
void ACID::StopAll() {
	std::cout << "\t\t\t\t\tStop All" << std::endl << std::endl;
	for (int i = 0; i < index; i++) {
		nclient.Send(cmdArr[i], sizeof(cmdArr[i]));
	}
}


void ACID::Receive(int iteration) {
	myfile.open("cmdRecord.bin", std::ios::binary);
	for (int i = 0; i < index; i++) {
		myfile << cmdArr[i] << std::endl;
	}
	myfile.close();

	std::cout << "\t\t\t\t\tReceive" << std::endl << std::endl;
	nserver.Receive(dataHash, iteration);
	CalculateCoeff();//module1
	Record();
}


void ACID::Record() {
	using namespace std;
	int row_val = 0;
	myfile.open("dataRecord.txt");
	for (const auto &row : dataHash) {
		if (row.size() != 0) {
			myfile << "Id: " << row_val << "\t| " << cmdName[row_val] << endl;
			for (const auto &element : row) {
				myfile << element << ", ";
			}
			myfile << endl << endl;
		}
		row_val++;
	}
	myfile.close();
}


void ACID::CalculateCoeff() {
	Input.Store(dataHash);
	Input.Calculate();
}
