//
// Created by jackcook on 7/15/20.
//

#include <iostream>
#include <vector>

using namespace std;

#ifndef CPPGFUNCTION_INTERPOLATION_H
#define CPPGFUNCTION_INTERPOLATION_H

namespace jcc::interpolation {

    double linterp(double xp, double x0, double y0, double x1, double y1);
    void interp1d(vector<double>& xp, vector<double>& yp, vector<double>& x, vector<double>& y, int ypSize=1);

} // jcc::interpolation

#endif //CPPGFUNCTION_INTERPOLATION_H
