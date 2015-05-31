#include <iostream>
#include <iomanip>
#include <vector>
#include "src/keplerian_toolbox.h"
#include "src/exposin/lambert_exposin.h"

using namespace std;
using namespace kep_toolbox;

int main() {
	array3D r1 = {149.0e9, 0.0, 0.0};
	array3D r2 = {200.0e9, -100.0e9, 0.0};
	double tof = 86400 * 565.0;
	double mu = 1.32e20;
	int lw = 1;
	int N = 0;
	double k2 = 0.175;
	lambert_exposin problem(r1, r2, tof, mu, lw, N, k2);

	array3D v1, v2;
	problem.get_v1(v1);
	problem.get_v2(v2);
	cout << v1 << endl;
	cout << v2 << endl;
	return 0;
}
