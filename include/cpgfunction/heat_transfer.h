//
// Created by jackcook on 7/11/20.
//

#ifndef CPPGFUNCTION_HEAT_TRANSFER_H
#define CPPGFUNCTION_HEAT_TRANSFER_H

#include <iostream>
#include <vector>
#include <cpgfunction/boreholes.h>

using namespace std;


namespace gt { namespace heat_transfer {

    struct SegmentResponse {
        ~SegmentResponse() {} // destructor

        int nSources;
        int nSum;
//        vector < vector < double > > h_ij;
        vector<double> H_ij;
        vector<gt::boreholes::Borehole> boreSegments;

        SegmentResponse(int nSources, int nSum, int nt) : nSources(nSources), boreSegments(nSources),
        H_ij(nSum * nt, 0), nSum(nSum)
        {} // constructor

        // storage_mode = 1 is the reduced segment response vector
        int storage_mode = 1;

//        void ReSizeContainers(int n, int nt);
        void get_h_value(double &h, int i, int j, int k);
        void get_index_value(int &index, int i, int j);
    };  // struct SegmentResponse();

    double finite_line_source(double time_, double alpha, gt::boreholes::Borehole& b1, gt::boreholes::Borehole& b2,
            bool reaSource=true, bool imgSource=true);
    void thermal_response_factors(SegmentResponse &SegRes, std::vector<gt::boreholes::Borehole>& boreSegments, std::vector<double>& time,
            double alpha, bool use_similaries, bool disp=false);

} } // namespace gt::heat_transfer

#endif //CPPGFUNCTION_HEAT_TRANSFER_H
