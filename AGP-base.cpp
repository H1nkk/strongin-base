#include <iostream>
#include <cmath>
#include <algorithm>
#include <map>
#include <numbers>
#include <random>
#include <chrono>
#include <vector>

using namespace std;

double r = 2.0; // method parameter
const double E = 1e-5; // epsilon
const double a = 2, b = 7; // left and right bounds
const int ITERMAX = 5000;

map<double (*)(double), double> extremums;
map<double (*)(double), double> leftBound;
map<double (*)(double), double> rightBound;
vector<double (*)(double)> funcs;

struct info {
	double extremumArg; // значение точки экстремума
	double extremumVal; // значение функции в точке экстремума
	int iterCount; // число совершенных итераций
	info(double extremumArg, double  extremumVal, int iterCount) : extremumArg(extremumArg), extremumVal(extremumVal), iterCount(iterCount) {}
};

double f1(double x) {
	if (x < (a + b) / 2) return 7 - x;
	double intersec = 7 - (a + b) / 2;
	return intersec - (a + b) / 2 + x;
}

double f2(double x) {
	return x * sin(x);
}

double f3(double x) {
	return x * cos(x);
}

double f4(double x) {
	return 354318.0 * x * cos(x);
}

double f5(double x) {
	return 3 * x * cos(4 * x);
}

double f6(double x) {
	return x * cos(x * x);
}

double f7(double x) {
	return -x - a;
}

double becnhFunc1(double x) {
	return sin(x) + sin(10.0 * x / 3.0);
}

double becnhFunc2(double x) {
	double res = 0;
	for (int k = 1; k <= 5; k++) {
		res += k * sin((k + 1) * x + k);
	}
	return -res;
}

double becnhFunc3(double x) {
	return (3.0 * x - 1.4) * sin(18.0 * x);
}

double becnhFunc4(double x) {
	double res = -(x + sin(x));
	res *= exp(-(x * x));
	return res;
}

double becnhFunc5(double x) {
	double res = sin(x) + sin(10.0 * x / 3.0) + log(x) - 0.84 * x + 3.0;
	return res;
}

double becnhFunc6(double x) {
	double res = sin(x) + sin(10.0 * x / 3.0) + log(x) - 0.84 * x + 3.0;
	return res;
}

double becnhFunc7(double x) {
	double res = -sin(2 * numbers::pi_v<double> *x) * exp(-x);
	return res;
}

double becnhFunc8(double x) {
	double res = (x * x - 5.0 * x + 6.0);
	res /= (x * x + 1.0);
	return res;
}

double becnhFunc9(double x) {
	double res = -x + sin(3.0 * x) - 1;
	return res;
}

double becnhFunc10(double x) {
	double res = 2 * (x - 3.0) * (x - 3.0) + exp(x * x * 0.5);
	return res;
}

void initMaps() {
	extremums[f1] = (a + b) / 2;

	// для a = 2, b = 7:
	extremums[f2] = 4.91318;
	extremums[f3] = 3.42562;
	extremums[f4] = 3.42562;
	extremums[f5] = 7;
	extremums[f6] = 6.86546;

	extremums[becnhFunc1] = 5.145735;
	leftBound[becnhFunc1] = 2.7;
	rightBound[becnhFunc1] = 7.5;
	funcs.push_back(becnhFunc1);

	extremums[becnhFunc2] = 5.791785; // тут несколько экстремумов
	leftBound[becnhFunc2] = 0.0;
	rightBound[becnhFunc2] = 10.0;
	funcs.push_back(becnhFunc2);

	extremums[becnhFunc3] = 0.96609;
	leftBound[becnhFunc3] = 0;
	rightBound[becnhFunc3] = 1.2;
	funcs.push_back(becnhFunc3);

	extremums[becnhFunc4] = 0.679560;
	leftBound[becnhFunc4] = -10;
	rightBound[becnhFunc4] = 10;
	funcs.push_back(becnhFunc4);

	extremums[becnhFunc5] = 5.19978;
	leftBound[becnhFunc5] = 2.7;
	rightBound[becnhFunc5] = 7.5;
	funcs.push_back(becnhFunc5);

	extremums[becnhFunc6] = 5.19978;
	leftBound[becnhFunc6] = 2.7;
	rightBound[becnhFunc6] = 7.5;
	funcs.push_back(becnhFunc6);

	extremums[becnhFunc7] = 0.224885;
	leftBound[becnhFunc7] = 0;
	rightBound[becnhFunc7] = 4;
	funcs.push_back(becnhFunc7);

	extremums[becnhFunc8] = 2.41420;
	leftBound[becnhFunc8] = -5;
	rightBound[becnhFunc8] = 5;
	funcs.push_back(becnhFunc8);

	extremums[becnhFunc9] = 5.877287;
	leftBound[becnhFunc9] = 0;
	rightBound[becnhFunc9] = 6.5;
	funcs.push_back(becnhFunc9);

	extremums[becnhFunc10] = 1.590700;
	leftBound[becnhFunc10] = -3;
	rightBound[becnhFunc10] = 3;
	funcs.push_back(becnhFunc10);
}

info AGP(double a, double b, double (*func)(double x)) {
	vector<double> dots = { a,b };
	vector<double> value = { func(a), func(b) }; // value[i] = значение функции в точке i 
	vector<double> R(1);
	double M;
	double m;

	sort(dots.begin(), dots.end());

	int iteration;
	for (iteration = 0; iteration < ITERMAX; iteration++) {

		int dotsCount = dots.size();

		for (int i = 0; i < dots.size(); i++) {
			value[i] = func(dots[i]);
		}

		M = fabs((value[dotsCount - 1] - value[dotsCount - 2]) / (dots[dotsCount - 1] - dots[dotsCount - 2]));

		for (int i = 1; i < dots.size(); i++) {
			M = max(M, fabs((value[i] - value[i - 1]) / (dots[i] - dots[i - 1])));
		}



		if (M > 0) {
			m = r * M;
		}
		else {
			m = 1;
		}

		double Rmax = -INFINITY;
		int maxInd = -1;
		for (int i = 1; i < dots.size(); i++) {
			R[i - 1] = m * (dots[i] - dots[i - 1])
				+ (value[i] - value[i - 1]) * (value[i] - value[i - 1]) / (m * (dots[i] - dots[i - 1]))
				- 2 * (value[i] - value[i - 1]);

			if (R[i - 1] > Rmax) {
				Rmax = R[i - 1];
				maxInd = i - 1;
			}
		}
		R.push_back(0);
		double newDot = 0.5 * (dots[maxInd + 1] + dots[maxInd]) - (value[maxInd + 1] - value[maxInd]) * 0.5 / m;
		dots.push_back(newDot);
		value.push_back(0);

		sort(dots.begin(), dots.end());
		double minLength = dots[1] - dots[0];
		for (int i = 2; i < dots.size(); i++) {
			minLength = min(minLength, dots[i] - dots[i - 1]);
		}
		if (minLength <= E) break;

	}

	double extrArg = dots[0], funcMin = func(dots[0]);
	for (auto dot : dots) {
		if (func(dot) < funcMin) {
			funcMin = func(dot);
			extrArg = dot;
		}
	}
	info res = { extrArg, funcMin, iteration };
	return res;
}

double findBestR(double (*testingFuncion)(double)) {
	cout << fixed;

	double extrPoint = extremums[testingFuncion];

	double best = r, diff = fabs(AGP(a, b, testingFuncion).extremumArg - extrPoint);
	for (r = 1.1; r < 100.0; r += 0.1) {
		double curdiff = fabs(AGP(a, b, testingFuncion).extremumArg - extrPoint);
		if (curdiff <= diff) {
			best = r;
			diff = curdiff;
			cout << best << ' ' << diff << '\n';
			if (curdiff == 0) break;
		}
	}

	cout << best << ' ' << diff << '\n';
	return best;
}

void timeMeasure() {
	double (*testingFunction)(double) = f2;

	r = 2.0;
	info res{ 0,0,0 };
	double avgTimeSpent = 0;
	int iters = 100;
	double foo = 4;
	for (int i = 0; i < iters; i++) { // усредняем время, т.к. оно колеблется - думаю, что из-за частого выделения памяти
		double (*testingFunction)(double) = funcs[rand() % 10];

		auto start = chrono::high_resolution_clock::now();
		res = AGP(a, b, testingFunction);
		auto stop = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
		double timeSpent = duration.count() / 1000000.0;
		avgTimeSpent += timeSpent;
		foo += res.extremumArg;
	}
	avgTimeSpent /= iters;

	cout << "Extremum argument: " << res.extremumArg << '\n';
	cout << "Iteration count: " << res.iterCount << '\n';
	cout << "Extremum value: " << res.extremumVal << '\n';
	cout << "Calculated in " << avgTimeSpent << " seconds\n";
	cout << "\nfoo for release:" << foo << "\n##########################\n";
}

void benchTests() {
	for (int i = funcs.size() - 1; i >= 0; i--) {
		info res = AGP(leftBound[funcs[i]], rightBound[funcs[i]], funcs[i]);
		if (res.extremumArg != extremums[funcs[i]]) {
			cout << "Func " << i + 1 << ". AGP result: " << res.extremumArg << ", actual result: " << extremums[funcs[i]] << '\n';
			cout << "Difference in results: " << fabs(res.extremumArg - extremums[funcs[i]]);
			cout << "\nIterations count : " << res.iterCount << "\n\n";
			cout << flush;

		}
	}
}

int main() {
	initMaps();

	cout << fixed;

	benchTests();

	timeMeasure();

	return 0;
}