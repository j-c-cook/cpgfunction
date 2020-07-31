//
// Created by jackcook on 7/11/20.
//

#ifndef CPPGFUNCTION_HEAT_TRANSFER_H
#define CPPGFUNCTION_HEAT_TRANSFER_H

#include <iostream>
#include <vector>
#include "gfunction.h"
#include "boreholes.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/asio.hpp>

using namespace std;


namespace gt::heat_transfer {
    double finite_line_source(double time_, double alpha, gt::boreholes::Borehole& b1, gt::boreholes::Borehole& b2,
            bool reaSource=true, bool imgSource=true);
    void thermal_response_factors(std::vector< std::vector< std::vector<double> > >& h_ij,
            std::vector<gt::boreholes::Borehole>& boreSegments, std::vector<double>& time,
            double alpha, bool use_similaries, bool disp=false);

    struct SimilaritiesType {
        ~SimilaritiesType() {} // destructor

        int nSim = 0;
        vector< vector <tuple<int, int> > > Sim;
        vector<tuple<double, double>> HSim;
        vector<tuple<double, double>> DSim;
        vector<double> disSim;

        SimilaritiesType() {} // constructor
    };

    struct Similarity {
        ~Similarity() {} // destructor

        // delcare variables
        // positive similarities
        vector<int> nSimPos;  // number of positive similarities
        vector< vector <tuple<int, int> > > simPos;  // number of sim positions
        vector<double> disSimPos;  // the distances between each position
        vector<tuple<int, int>> HSimPos;  // the heights of each borehole
        vector<tuple<int, int>> DSimPos;  // the burial depth of each borehole
        // negative similarities
        vector<int> nSimNeg;
        vector< vector <tuple<int, int> > > simNeg;
        vector<double> disSimNeg;
        vector<tuple<int, int>> HSimNeg;
        vector<tuple<int, int>> DSimNeg;

        Similarity(){} // constructor

        void similarities(SimilaritiesType &SimReal, SimilaritiesType &SimImage,
                vector<gt::boreholes::Borehole> &boreSegments,
                bool splitRealAndImage=true, double disTol=0.1, double tol=1.0e-6);

        void _similarities_group_by_distance(vector<gt::boreholes::Borehole> &boreSegments,
                vector< vector < tuple <int, int> > > &Pairs,
                vector<int> &nPairs, vector<double> &disPairs, int &nDis,
                double disTol=0.1);
        void _similarities_one_distance(SimilaritiesType &SimT, vector< tuple <int, int> >& pairs,
                vector<gt::boreholes::Borehole> &boreSegments, const string& kind, double tol=1.0e-6);
    };

} // namespace gt::heat_transfer

#endif //CPPGFUNCTION_HEAT_TRANSFER_H
