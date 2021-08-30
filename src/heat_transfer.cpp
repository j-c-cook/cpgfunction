//
// Created by jackcook on 7/11/20.
//

#include <cpgfunction/heat_transfer.h>
#include <stdexcept>
#include <thread>
#include <boost/asio.hpp>
#include <cpgfunction/boreholes.h>
#include <qdt.h>

using namespace boost::math::quadrature;

namespace gt { namespace heat_transfer {
    double finite_line_source(const double time_, const double alpha, gt::boreholes::Borehole &b1,
                              gt::boreholes::Borehole &b2, bool reaSource, bool imgSource) {

        auto _Ils = [&b1, &b2, reaSource, imgSource](const double s) {
            auto _erfint = [](const double x) {
                return x * std::erf(x) - (1 / sqrt(M_PI)) * (1 - exp(-pow(x, 2)));
            };
            double r = b1.distance(b2);
            double func = 0.;
            // function to integrate
            if (reaSource) {
                // Real part of the FLS solution
                func += _erfint(double(b2.D - b1.D + b2.H) * s);
                func += -_erfint(double(b2.D - b1.D) * s);
                func += _erfint(double(b2.D - b1.D - b1.H) * s);
                func += -_erfint(double(b2.D - b1.D + b2.H - b1.H) * s);
            } // fi reaSource
            if (imgSource) {
                // Image part of the FLS solution
                func += _erfint(double(b2.D + b1.D + b2.H) * s);
                func += -_erfint(double(b2.D + b1.D) * s);
                func += _erfint(double(b2.D + b1.D + b1.H) * s);
                func += -_erfint(double(b2.D + b1.D + b2.H + b1.H) * s);
            } // fi imgSource
            double a = 0.5 / (b2.H * pow(s, 2)) * func * exp(-pow(r, 2) * pow(s, 2));
            return a;
        }; // auto _Ils

        // lower bound of integration
        double a = double(1.) / sqrt(double(4.) * alpha * time_);
        // Evaluate the integral using Gauss-Kronrod
        double Q;
        auto method = qdt::adaptive(qdt::gauss_kronrod());
        Q = method.integrate(_Ils, a, qdt::INF);
        return Q;

    } // void finite_line_source

    void
    thermal_response_factors(SegmentResponse &SegRes, std::vector<gt::boreholes::Borehole> &boreSegments,
                             std::vector<double> &time,
                             const double alpha, bool use_similaries, bool disp) {
        // total number of line sources
        int nSources = boreSegments.size();
        // number of time values
        int nt = time.size();

        // Open up processes here
        // Create a vector of threads
        //may return 0 when not able to detect
        const auto processor_count = thread::hardware_concurrency();
        // Launch the pool with n threads.
        boost::asio::thread_pool pool(processor_count);
        if (disp) {
            cout << "\tDetected " << processor_count << " as the number of available threads" << endl;
        }

        gt:boreholes::SimilaritiesType SimReal; // positive
        gt::boreholes::SimilaritiesType SimImage; // negative

        int COUNT=0;
//        auto check_key = [&h_map](Key &sim_key) {
//            // Key is not present
//            if (h_map.find(sim_key) == h_map.end())
//                return false;
//            return true;
//        }; // check_key

        auto sum_to_n = [](const int n) {
            return n * (n + 1) / 2;
        };

        if (use_similaries) {
            auto start = std::chrono::steady_clock::now();
            // Calculations with similarities
            if (disp) {
                cout << "Identifying similarities..." << endl;
            }
            bool splitRealAndImage = true;
            double disTol = 0.1;
            double tol = 1.0e-6;
            gt::boreholes::Similarity sim;
            sim.similarities(SimReal, SimImage, boreSegments, splitRealAndImage, disTol, tol);

            //---
            // Adaptive hashing scheme if statement
            // Determine the Segment Response storing mode here
            int Ntot = sum_to_n(nSources);

            // lambda function for calculating h at each time step
            auto _calculate_h = [&boreSegments, &splitRealAndImage, &time, &alpha, &nt, &SegRes, &Ntot](boreholes::SimilaritiesType &SimReal,
                    int s, bool reaSource, bool imgSource) {
                // begin function
                int n1;
                int n2;
                gt::boreholes::Borehole b1;
                gt::boreholes::Borehole b2;
                // begin thread
                n1 = get<0>(SimReal.Sim[s][0]);
                n2 = get<1>(SimReal.Sim[s][0]);
                b1 = boreSegments[n1];
                b2 = boreSegments[n2];
                vector<double> hPos(nt);
                if (splitRealAndImage) {
                    for (int k=0; k<nt; k++) {
                        hPos[k] = finite_line_source(time[k], alpha, b1, b2, reaSource, imgSource);
                    }  // next k
                    int i;
                    int j;
                    if (SegRes.storage_mode==0) {
                        // not looping through every (i, j), real and image stored separately
                        // TODO: make the option to store entire matrix
                        throw std::invalid_argument("This storage mode does not currently exist.");
                    } else if (SegRes.storage_mode==1) {
                        // will loop through every (i, j), will combine real+image
                        int index;
                        for (int k=0; k<SimReal.Sim[s].size(); k++) {
                            i = get<0>(SimReal.Sim[s][k]);
                            j = get<1>(SimReal.Sim[s][k]);
                            for (int t=0; t<time.size(); t++){
                                // must consider real and image source separate when combining
                                if (i <= j) {
                                    // we want to store n2, n1
                                    SegRes.get_index_value(index, i, j);
                                    SegRes.h_ij[index][t] += b2.H / b1.H * hPos[t]; // non-critical race condition
                                } else {
                                    SegRes.get_index_value(index, j, i);
                                    SegRes.h_ij[index][t] += hPos[t]; // non-critical race condition
                                }  // else ()
                            }  // next t
                        }  // next k
                    }  // else if(SegRes.storage_mode==1)
                } else {
                    throw std::invalid_argument( "Not yet written yet.");
                }
            };
            auto end = std::chrono::steady_clock::now();
            if (disp) {
                double milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                double seconds = milli / 1000;
                std::cout << "Elapsed time in seconds : "
                          << seconds
                          << " sec" << std::endl;
                std::cout << "Calculating segment to segment response factors ..." << std::endl;
            } // end if

            // inputs
            bool reaSource;
            bool imgSource;
            for (int s=0; s<SimReal.nSim; s++) {
                reaSource = true;
                imgSource = false;
                boost::asio::post(pool, [&_calculate_h, &SimReal, s, reaSource, imgSource]
                { _calculate_h(SimReal, s, reaSource, imgSource); });
//                _calculate_h(SimReal, s, reaSource, imgSource, hash_mode);
            } // next s
            if (splitRealAndImage) {
                reaSource = false;
                imgSource = true;
                for (int s=0; s<SimImage.nSim; s++) {
                    boost::asio::post(pool, [&_calculate_h, &SimImage, s, reaSource, imgSource]
                    { _calculate_h(SimImage, s, reaSource, imgSource); });
//                    _calculate_h(SimImage, s, reaSource, imgSource, hash_mode);
                }
            }
            pool.join();
            auto end2 = std::chrono::steady_clock::now();
            if (disp) {
                double milli = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - end).count();
                double seconds = milli / 1000;
                std::cout << "Elapsed time in seconds : "
                          << seconds
                          << " sec" << std::endl;
            }
        } else {
            if (disp) {
                std::cout << "Calculating segment to segment response factors ..." << std::endl;
            } // end if
            auto start = std::chrono::steady_clock::now();
            bool sameSegment;
            bool otherSegment;

//            auto _fill_line = [&h_ij, &time, &boreSegments](const int i, const int j, const double alpha,
//                    bool sameSegment, bool otherSegment) {
//                auto _dot_product = [&h_ij, &time](const int i, const int j, const double constant) {
//                    for (int k=0; k < time.size(); k++) {
//                        h_ij[j][i][k+1] = constant * h_ij[i][j][k+1];
//                    } // end for
//                };
//                double h;
//                double constant;
//                gt::boreholes::Borehole b1;
//                gt::boreholes::Borehole b2;
//                b2 = boreSegments[i];
//                for (int k = 0; k < time.size(); k++) {
//                    double t = time[k];
//                    if (sameSegment && not otherSegment){
//                        b1 = boreSegments[i];
//                        h = finite_line_source(t, alpha, b2, b2);
//                    } else if (otherSegment && not sameSegment) {
//                        b1 = boreSegments[j];
//                        h = finite_line_source(t, alpha, b1, b2);
//                    } else {
//                        throw std::invalid_argument( "sameSegment and otherSegment cannot both be true" );
//                    } // end if
//                    h_ij[i][j][k+1] = h;
//                    if (otherSegment && not sameSegment) {
//                        constant = double(b2.H / b1.H);
//                        _dot_product(i, j, constant);
//                    } // end if
//                }; // end for
//            }; // auto _fill_line
//
//            for (int i = 0; i < nSources; i++) {
//                // Segment to same-segment thermal response factor
//                // FLS solution for combined real and image sources
//                sameSegment = true;
//                otherSegment = false;
//                boost::asio::post(pool, [i, alpha, sameSegment, otherSegment, &_fill_line]
//                { _fill_line(i, i, alpha, true, false); });
////                _fill_line(i, i, alpha, sameSegment, otherSegment); // could call with no threading during debugging
//
//                // Segment to other segment thermal response factor
//                for (int j = i + 1; j<nSources; j++) {
//                    sameSegment = false;
//                    otherSegment = true;
//                    boost::asio::post(pool, [i, j, alpha, sameSegment, otherSegment, &_fill_line]
//                    { _fill_line(i, j, alpha, sameSegment, otherSegment); });
////                    _fill_line(i, j, alpha, sameSegment, otherSegment);  // could call with no threading during debugging
//                } // end for
//            } // fi (end if)

            /** Wait for all the threads in vector to join **/
            pool.join();
            auto end = std::chrono::steady_clock::now();
            if (disp) {
                double milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                double seconds = milli / 1000;
                std::cout << "Elapsed time in seconds : "
                          << seconds
                          << " sec" << std::endl;
            }
            // Iterate over the thread vector
        } // fi similarity

//        for (int i = 0; i < nSources; i++) {
//            for (int j=0; j<nSources; j++) {
//                for (int k=0; k<nt; k++) {
//                    std::cout << h_ij[i][j][k+1] << "\t";
//                }
//                std::cout << std::endl;
//            }
//            int a = 1;
//        }
//        gt::boreholes::Borehole b1;
//        gt::boreholes::Borehole b2;
//        int Ntot = sum_to_n(nSources);
//        for (int i = 0; i < nSources; i++) {
//            for (int j=0; j<nSources; j++) {
//                for (int k=0; k<nt; k++) {
//                    double v;
//                    double h;
//                    double hij;
//                    double ratio;
//                    if (i <= j) {
//                        int index = SegRes.get_index_value(i, j);
//                        h = SegRes.h_ij[index][k];
//                    } else {
//                        int index = SegRes.get_index_value(j, i);
//                        b1 = boreSegments[i];
//                        b2 = boreSegments[j];
//                        h = b2.H/b1.H * SegRes.h_ij[index][k];
//                    }
////                    hij = h_ij[i][j][k+1];
//                    SegRes.get_h_value(v, i, j, k);
//                    h_ij[i][j][k+1] = h;
//                    if (h - v > 1.0e-6) {
//                        ratio = h / hij;
//                        int a = 1;
//                    }
//                }
//            }
//        }
        int a =1;
    } // void thermal_response_factors

//    void SegmentResponse::ReSizeContainers(const int n, const int nt) {
//        switch(storage_mode) {
//            case 0 :
//                cout << "Case 0 not written yet" << endl;
//                break;
//            case 1 :
//                h_ij.resize(n, vector<vector<double>>(n, vector<double>(nt, 0)));
//                break;
//            default:
//                throw invalid_argument("The case selected is not currently implemented.");
//        }
//
//    }

    void SegmentResponse::get_h_value(double &h, const int i, const int j, const int k) {
        int index;
        switch (storage_mode) {
            case 0 :
                cout << "Case 0 not written yet" << endl;
                break;
            case 1 :
                if (i <= j) {
                    get_index_value(index, i, j);
                    h = h_ij[index][k];
                } else {
                    get_index_value(index, j, i);
                    h = boreSegments[j].H/boreSegments[i].H * h_ij[index][k];
                }
                break;
            default:
                throw invalid_argument("The case selected is not currently implemented.");
        }  // switch();
    }  // SegmentResponse::get_h_value();

    void SegmentResponse::get_index_value(int &index, const int i, const int j) {
        index = i * (2*nSources - i - 1) / 2 + j;
    }  // SegmentResponse::get_index_value();


} } // namespace gt::heat_transfer