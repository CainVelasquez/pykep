#include <iostream>
#include <iomanip>
#include <vector>
#include "src/keplerian_toolbox.h"

using namespace std;
using namespace kep_toolbox;

int main() {
    array3D r1 = {149.0e9, 0.0, 0.0};
    array3D r2 = {200.0e9, -100.0e9, 0.0};
    double tof = 86400 * 565.0;
    double mu = 1.32e20;
    bool lw = true;
    lambert_exposin problem(r1, r2, tof, mu, lw);
    cout << "Revs: 0 to " << problem.get_revs() << endl;
    cout << "Done" << endl;
    return 0;
}
