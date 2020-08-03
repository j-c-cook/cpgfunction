// -*- lsst-c++ -*-

//
// Created by jackcook on 7/11/20.
//

#ifndef CPPGFUNCTION_GFUNCTION_H
#define CPPGFUNCTION_GFUNCTION_H

#include <iostream>
#include <vector>
#include "boreholes.h"
#include "heat_transfer.h"

using namespace std;

/** The functions will be listed here in the order which they appear in the associated cpp file **/

namespace gt::gfunction {
        /**
         * Uniform borehole wall temperature (UBWHT) g-function calculation method
         *
         * This function superimposes the finite line source (FLS) solution to estimate the g-function of a
         * geothermal bore field. Each borehole is modeled as a series of finite line source segments, as proposed
         * in [CITE: CimminoBernier2014].
         *
         * @param gfunction
         * @param boreholes
         * @param time
         * @param alpha
         * @param nSegments
         * @param use_similarities
         * @param disp
         */
        void uniform_temperature(vector<double>& gfunction, vector<gt::boreholes::Borehole> boreholes,
                vector<double>& time, double alpha, int nSegments=12, bool use_similarities=true ,bool disp=false);
        void _borehole_segments(vector<gt::boreholes::Borehole>& boreSegments,
                                vector<gt::boreholes::Borehole>& boreholes, int nSegments);
        void load_history_reconstruction(vector<vector<double>>& q_reconstructed, vector<double>& time, vector<double>& _time,
                vector<vector<double> >& Q, vector<double>& dt, const int p);
        void _temporal_superposition(vector<double>& Tb_0, vector<vector<vector<double> > >& h_ij,
                                     vector<vector<vector<double> > >& dh_ij,
                std::vector<std::vector<double>>& q_reconstructed, int p);
        void _solve_eqn(vector<double>& x, vector<vector<double>>& A, vector<double>& b);

} // namespace gt::gfunction

#endif //CPPGFUNCTION_GFUNCTION_H
