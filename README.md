//Written by Khang Pham at LSU

///////////////////////////////////////////////////////////////
//This program can be used to solve a 2D Schrodinger Equation//
//with uncoupled potential using the Thomas algorithm for    //
//periodic boudary and alternate direction implicit method   //
///////////////////////////////////////////////////////////////

#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>
#include <string.h>
#include <assert.h>
#include <array>
#include <iostream>
#include <complex>
#include <set>
#include <algorithm>
#include <future>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <stdio.h>
#include <string>
#include <math.h>

#define PI M_PI
#define h_bar 1. 
#define m1 1.
#define m2 1.
#define L 10. //Length of axis
#define T 10  //Total run time
#define grid_point 100. //number of points on an axis

typedef std::complex<double> compx;
#define I compx(0.,1.)
#define one compx(1.,0.)
#define two compx(2.,0.)

class wave_function {
public:
	wave_function();
	wave_function(bool);
	double dt = .01;
	double dx = L/grid_point;
	std::vector<std::vector<compx>> value;
	void solve_triag();
	double potential(int, int);
	double real_space(int);
	compx sec_space_deriv(char, int, int);
	void normalize();
};

//default constructor that initializes all value of psi to 0
wave_function::wave_function() {
	value.resize(grid_point);
	for (int l = 0; l < grid_point; l++) {
		value.resize(grid_point);
		for (int m = 0; m < grid_point; m++)
			value.at(m).push_back(0.0);
	}
}

//constructor that initializes psi to a desired initial condition
wave_function::wave_function(bool init) {
	if (init) {
		compx psi00, psi10, psi12;
		value.resize(grid_point);
		for (int l = 0; l < grid_point; l++) {
			value[l].resize(grid_point);
			for (int m = 0; m < grid_point; m++) {
				psi00 = pow(m1 / (PI*h_bar), .25)*exp(-m1*pow(real_space(l), 2.) / (2.*h_bar))*pow(m1*3.0 / (PI*h_bar), .25)*exp(-m1*pow(real_space(m), 2.) / (2.*h_bar)); 
				psi10 = pow(m1 / (PI*h_bar), .25)*pow(2., .5)*pow(m1 / (h_bar), .5)*real_space(l)*exp(-m1*pow(real_space(l), 2.) / (2.*h_bar))*pow(m1 / (PI*h_bar), .25)*exp(-m1*pow(real_space(m), 2.) / (2.*h_bar));
				psi12 = pow(m1 / (PI*h_bar), .25)*pow(2., .5)*pow(m1 / (h_bar), .5)*real_space(l)*exp(-m1*pow(real_space(l), 2.) / (2.*h_bar))*pow(m1 / (PI*h_bar*4.), .25)*(2.*(m1 / h_bar)*pow(real_space(m),2.) - 1.)*exp(-m1*pow(real_space(m), 2.) / (2.*h_bar));
				value[l][m] = psi00 +psi10 + psi12;
			}
		}
		normalize();
	}
	else if (!init) {
		wave_function();
	}
}



//Solving for values using Thomas algorithm and ADI method
void wave_function::solve_triag() {

	compx a = -h_bar / (two*m1), b = -h_bar / (two*m2);
	double r = dt / (dx*dx);
	compx A = (I*r*a / two), B = (I*r*b / two), C = I*dt / two;
	compx mid = one - two*A;
	compx x_N, q_N, tvalue;
	std::vector<compx> alpha(grid_point - 2);
	std::vector<compx> beta1(grid_point - 1);
	std::vector<compx> beta2(grid_point - 1);
	wave_function tmp; 
	wave_function x_1, x_2; 

	//Solving for the x1 direction first
	for (int x2 = 0; x2 < grid_point; x2++) {

		alpha[0] = A / mid;
		beta1[0] = ((one - C*potential(0, x2))*value[0][x2] - B*sec_space_deriv('y', 0, x2)) / mid;
		beta2[0] = -A / mid;

		//Forward run
		for (int l = 1; l < grid_point - 2; l++) {
			alpha[l] = A / (mid - A*alpha[l - 1]);

			beta1[l] = ((one - C*potential(l, x2))*value[l][x2] - B*sec_space_deriv('y', l, x2)
				- A*beta1[l - 1]) / (mid - A*alpha[l - 1]);
			beta2[l] = -A*beta2[l - 1] / (mid - A*alpha[l - 1]);
		}

		beta1[grid_point - 2] = ((one - C*potential(grid_point - 2, x2))*value[grid_point - 2][x2] - B*sec_space_deriv('y', grid_point - 2, x2)
			- A*beta1[grid_point - 3]) / (mid - A*alpha[grid_point - 3]);
		beta2[grid_point - 2] = (-A - A*beta2[grid_point - 3]) / (mid - A*alpha[grid_point - 3]);

		//Backward run
		x_1.value[grid_point - 2][x2] = beta1[grid_point - 2];
		x_2.value[grid_point - 2][x2] = beta2[grid_point - 2];
		for (int l = grid_point - 3; l >= 0; l--) {
			x_1.value[l][x2] = beta1[l] - alpha[l] * x_1.value[l + 1][x2];
			x_2.value[l][x2] = beta2[l] - alpha[l] * x_2.value[l + 1][x2];
		}

		compx q_N = (one - C*potential(grid_point - 1, x2))*value[grid_point - 1][x2] - B*sec_space_deriv('y', grid_point - 1, x2);
		compx x_N = (q_N - A*x_1.value[0][x2] - A*x_1.value[grid_point - 2][x2])
			/ (mid + A*x_2.value[0][x2] + A*x_2.value[grid_point - 2][x2]);

		//Solve for tmp values (aka n*)
		tmp.value[grid_point - 1][x2] = x_N;
		for (int l = 0; l < grid_point - 1; l++) {
			tmp.value[l][x2] = x_1.value[l][x2] + x_2.value[l][x2] * x_N;
			tvalue = tmp.value[l][x2];
		}
	}

	tmp.normalize();

	//Solving for the x2 direction
	for (int x1 = 0; x1 < grid_point; x1++) {

		alpha[0] = B / (one - two*B + C*potential(x1, 0));
		beta1[0] = (tmp.value[x1][0] - A*tmp.sec_space_deriv('x', x1, 0)) / (one - two*B + C*potential(x1, 0));
		beta2[0] = -B / (one - two*B + C*potential(x1, 0));

		//Forward run
		for (int m = 1; m < grid_point - 2; m++) {
			alpha[m] = B / (one - two*B + C*potential(x1, m) - B*alpha[m - 1]);

			beta1[m] = (tmp.value[x1][m] - A*tmp.sec_space_deriv('x', x1, m) - B*beta1[m - 1])
				/ (one - two*B + C*potential(x1, m) - B*alpha[m - 1]);
			beta2[m] = -B*beta2[m - 1] / (one - two*B + C*potential(x1, m) - B*alpha[m - 1]);
		}

		beta1[grid_point - 2] = (tmp.value[x1][grid_point - 2] - A*tmp.sec_space_deriv('x', x1, grid_point - 2) - B*beta1[grid_point - 3])
			/ (one - two*B + C*potential(x1, grid_point - 2) - B*alpha[grid_point - 3]);
		beta2[grid_point - 2] = (-B - B*beta2[grid_point - 3]) / (one - two*B + C*potential(x1, grid_point - 2) - B*alpha[grid_point - 3]);

		//Backward run
		x_1.value[x1][grid_point - 2] = beta1[grid_point - 2];
		x_2.value[x1][grid_point - 2] = beta2[grid_point - 2];
		for (double m = grid_point - 3; m >= 0; m--) { //-2 because of the same reason above
			x_1.value[x1][m] = beta1[m] - alpha[m] * x_1.value[x1][m + 1];
			x_2.value[x1][m] = beta2[m] - alpha[m] * x_2.value[x1][m + 1];
		}

		q_N = tmp.value[x1][grid_point - 1] - A*tmp.sec_space_deriv('x', x1, grid_point - 1);
		x_N = (q_N - B*x_1.value[x1][0] - B*x_1.value[x1][grid_point - 2])
			/ ((one - two*B + C*potential(x1, grid_point - 1)) + B*x_2.value[x1][0] + B*x_2.value[x1][grid_point - 2]);

		value[x1][grid_point - 1] = x_N;
		for (int m = 0; m < grid_point - 1; m++) {
			value[x1][m] = x_1.value[x1][m] + x_2.value[x1][m] * x_N;
		}

	}

	normalize();
}

double wave_function::potential(int x1, int x2) {
	//2d harmonic oscillator
	return .5*m1*real_space(x1)*real_space(x1) + .5*m2*real_space(x2)*real_space(x2);
}

double wave_function::real_space(int index) {
	return (-L / 2.0) + (index*dx);
}

//x stands for x1 or particle 1 and y stands for x2 or particle 2
compx wave_function::sec_space_deriv(char x1_or_x2, int l, int m) {
	if (x1_or_x2 == 'x') {
		if (l == 0)
			return value[l + 1][m] - two * value[l][m] + value[grid_point - 1][m];
		else if (l == grid_point - 1)
			return value[0][m] - two * value[l][m] + value[l - 1][m];
		else
			return value[l + 1][m] - two * value[l][m] + value[l - 1][m];
	}
	else if (x1_or_x2 == 'y') {
		if (m == 0)
			return value[l][m + 1] - two * value[l][m] + value[l][grid_point - 1];
		else if (m == grid_point - 1)
			return value[l][0] - two * value[l][m] + value[l][m - 1];
		else
			return value[l][m + 1] - two * value[l][m] + value[l][m - 1];
	}
}

//Copied from stack overflow to deal with string for now 
template <typename NEW>
std::string to_string_with_precision(const NEW a_value, const int n = 6)
{
	std::ostringstream out;
	out << std::fixed << std::setprecision(n) << a_value;
	return out.str();
}


//normalization
void wave_function::normalize() {
	compx sum = 0;
	for (int i1 = 0; i1 < grid_point; i1++) {
		for (int i2 = 0; i2 < grid_point; i2++) {
			sum += pow(abs(value[i1][i2]), 2.0);
		}
	}
	sum *= dx * dx;
	compx amplitude = one / pow(sum, .5);
	for (int i1 = 0; i1 < grid_point; i1++) {
		for (int i2 = 0; i2 < grid_point; i2++) {
			value[i1][i2] *= amplitude;
		}
	}
}

int main() {
	wave_function v(true);
	std::ofstream file1; //contains values
	std::ofstream file2; //potential data
	file2.open("potential.dat");
	for (int i = 0; i < grid_point; i++)
		for (int j = 0; j < grid_point; j++)
			file2 << v.real_space(i) << "\t" << v.real_space(j) << "\t" << v.potential(i, j) << std::endl;

	int index = 0; //this is to space out data files
	for (double k = v.dt; k <= T; k += v.dt) {
		std::cout << index << std::endl;
		if (index % 10 == 0) {
			file1.open("data_" + to_string_with_precision(index, 0) + ".dat");
			file1 << "Time   " << k - v.dt << std::endl
				<< "x" << "\t" << "y" << "\t" << "imag" << "\t" << "real" << "\t" << "abs" << std::endl;
			for (int i = 0; i < grid_point; i++) {
				for (int j = 0; j < grid_point; j++) {
					file1 << v.real_space(i) << "\t"
						<< v.real_space(j) << "\t"
						<< imag(v.value[i][j]) << "\t"
						<< real(v.value[i][j]) << "\t"
						<< abs(v.value[i][j]) << std::endl;
				}
			}
			file1.close();
		}
		v.solve_triag();
		index++;
	}
	getchar();
}
