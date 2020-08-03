// -*- lsst-c++ -*-

//
// Created by jackcook on 7/11/20.
//

#include "gfunction.h"
#include <chrono>
#include "spline.h"
#include "gauss_jacobi.h"
#include "../jcc/interpolation.h"
#include <thread>
#include <boost/asio.hpp>


extern "C" void dgesv_( int *n, int *nrhs, double  *a, int *lda, int *ipiv, double *b, int *lbd, int *info  );

using namespace std;  // lots of vectors, only namespace to be used

namespace gt::gfunction {
    // The uniform borehole wall temperature (UBWHT) g-function calculation. Originally presented in
    // Cimmino and Bernier (2014) and a later paper on speed improvements by Cimmino (2018)
    void uniform_temperature(vector<double>& gfunction, vector<gt::boreholes::Borehole> boreholes,
                             vector<double>& time, const double alpha, const int nSegments,
                             const bool use_similarities, const bool disp) {
        // TODO: place this resizing into a "general" namespace
        int len_actual = time.size(); // check to see if there is enough space in the vector
        int len_g = gfunction.size();
        // if need be, resize the vector to be the same size as the number of boreholes needed
        if (len_actual != len_g) {
            gfunction.resize(len_actual);
        } else ; // else do nothing

        if (disp) {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "Calculating g-function for uniform borehole wall temperature" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        auto startall = std::chrono::steady_clock::now();
        // Open up processes here
        // Create a vector of threads
        //may return 0 when not able to detect
        const auto processor_count = thread::hardware_concurrency();
//        // Launch the pool with n threads.
//        boost::asio::thread_pool pool(processor_count);
//        cout << "\tDetected " << processor_count << " as the number of available threads" << endl;

        // Number of boreholes
        int nbh = boreholes.size();
        // Total number of line sources
        int nSources = nSegments * nbh;
        // Number of time values
        int nt = time.size();

        // Split boreholes into segments
        vector<gt::boreholes::Borehole> boreSegments(nSources);
        _borehole_segments(boreSegments, boreholes, nSegments);

        // Initialize segment-to-segment response factors (https://slaystudy.com/initialize-3d-vector-in-c/)
        // NOTE: (nt + 1), the first row will be full of zeros for later interpolation
        vector< vector< vector<double> > > h_ij(nSources ,
                vector< vector<double> > (nSources, vector<double> (nt+1, 0.0)) );
        // Calculate segment to segment thermal response factors
        auto start = std::chrono::steady_clock::now();
        gt::heat_transfer::thermal_response_factors(h_ij, boreSegments, time, alpha, use_similarities, disp);
        auto end = std::chrono::steady_clock::now();

        if (disp) {
            std::cout << "Building and solving system of equations ..." << std::endl;
        }
        // -------------------------------------------------------------------------
        // Build a system of equation [A]*[X] = [B] for the evaluation of the
        // g-function. [A] is a coefficient matrix, [X] = [Qb,Tb] is a state
        // space vector of the borehole heat extraction rates and borehole wall
        // temperature (equal for all segments), [B] is a coefficient vector.
        // -------------------------------------------------------------------------

        // -------- timings for debug
        double milli = 0;
        double segment_length_time = 0;
        double time_vector_time = 0;
        double segment_h_values_time = 0;
        double fill_A_time = 0;
        double load_history_reconstruction_time = 0;
        double temporal_superposition_time = 0;
        double fill_gsl_matrices_time = 0;
        double LU_decomposition_time = 0;

        auto start2 = std::chrono::steady_clock::now();

        // ------ Segment lengths -------
        start = std::chrono::steady_clock::now();
        std::vector<float> Hb(nSources);
        auto _segmentlengths = [&boreSegments, &Hb](const int nSources) {
            for (int b=0; b<nSources; b++) {
                Hb[b] = boreSegments[b].H;
            } // next b
        }; // auto _segmentlengths
//        boost::asio::post(pool, [nSources, &boreSegments, &Hb, &_segmentlengths]{ _segmentlengths(nSources); });
        _segmentlengths(nSources);
        end = std::chrono::steady_clock::now();
        milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        segment_length_time += milli;

        // ------ time vectors ---------
        start = std::chrono::steady_clock::now();
        // create new time vector that starts at 0
        std::vector<double> _time_untouched(time.size()+1);
        std::vector<double> _time(time.size()+1);
        std::vector<double> dt(time.size());

        auto _fill_time = [&_time, &time, &dt, &_time_untouched]() {
            for (int i=0; i<_time.size(); i++) {
                if (i==0) {
                    _time_untouched[0] = 0;
                    _time[0] = 0;
                    dt[i] = time[i];
                } else {
                    _time_untouched[i] = time[i-1];
                    _time[i] = time[i-1];
                    dt[i] = time[i] - time[i-1];
                } // fi
            } // next i
        }; // auto _fill_time
//        boost::asio::post(pool, [&_fill_time, &time, &_time]{ _fill_time() ;});
        _fill_time();

        end = std::chrono::steady_clock::now();
        milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        time_vector_time += milli;

//         TODO: come back to thinking about how to thread this or another (depracated)

//        auto _fill_dt = [&dt, &time]() {
//            for (int i=0; i<time.size(); i++) {
//                if (i==0) {
//                    dt[i] = time[i];
//                } else {
//                    dt[i] = time[i] - time[i-1];
//                } // fi
//            } // end i
//        };
//        _fill_dt();

//        pool.join(); // starting up a new idea after this, pool will close here


//        cout << "\tDetected " << processor_count << " as the number of available threads" << endl;

        // ---------- segment h values -------------
        /** Starting up pool2 here **/
        // Launch the pool with n threads.
        auto tic = std::chrono::steady_clock::now();
        boost::asio::thread_pool pool2(processor_count);
        auto toc = std::chrono::steady_clock::now();
        if (disp) {
            double milli = std::chrono::duration_cast<std::chrono::milliseconds>(tic - toc).count();
            double seconds = milli;
            std::cout << "Time to open a pool : "
                      << seconds
                      << " sec" << std::endl;
        }

        start = std::chrono::steady_clock::now();
        // h_dt : an interpolation through the k direction of the cube
//        vector< vector< vector<double> > > h_dt(nSources ,
//                vector< vector<double> > (nSources, vector<double> (nt)) );
        vector< vector< vector<double> > > h_dt(nSources ,vector< vector<double> > (nSources, vector<double> (nt)) );
        // dh_ij is the difference of h_ij, ie. dh_[i][j][k] = h_ij[i][j][k]-h_ij[i][j][k-1], ie. \delta h_ij
        std::vector< std::vector< std::vector<double> > > dh_ij(nSources ,
                std::vector< std::vector<double> > (nSources, std::vector<double> (nt)) );

        // Thermal response factors evaluated at t=dt (h_dt)
        auto _interpolate = [&h_ij, &h_dt, &_time, &dt, &dh_ij](const int i) {
            std::vector<double> y(_time.size());
            for (int j=0; j <h_ij[i].size(); j++) {
                for (int k=0; k<h_ij[i][j].size(); k++) {
                    if (k==1) {
                        dh_ij[i][j][k-1] = h_ij[i][j][k];  // start at time "1"
                    } else if (k>1) {
                        dh_ij[i][j][k-1] = h_ij[i][j][k] - h_ij[i][j][k-1];

                    } // fi
                    y[k] = h_ij[i][j][k];
                } // end k

                std::vector<double> yp(dt.size());
                jcc::interpolation::interp1d(dt, yp, _time, y, yp.size() );

                for (int k=0; k<h_dt[i][j].size(); k++) {
                    h_dt[i][j][k] = yp[k];
                } // end k
            } // next j

        }; // _interpolate
        // h_dt for loop
        for (int i=0; i<h_ij.size(); i++) {
//            for (int j=0; j<h_ij[i].size(); j++) {
//                boost::asio::post(pool2, [&_interpolate, i]{ _interpolate(i) ;});
                 _interpolate(i);
//            } // end j
        } // end i

        // if _interpolate threaded, join() pools here.
        pool2.join();  // need interpolated values moving forward
        /** Starting up pool3 here **/
        // Launch the pool with n threads.
//        boost::asio::thread_pool pool3(processor_count);
//        cout << "\tDetected " << processor_count << " as the number of available threads" << endl;

//        auto _dh_ij = [&h_ij, &dh_ij](const int i, const int j) {
//            for (int k=1; k<h_ij[i][j].size(); k++) {
//                if (k==1) {
//                    dh_ij[i][j][k-1] = h_ij[i][j][k];  // start at time "1"
//                } else {
//                    dh_ij[i][j][k-1] = h_ij[i][j][k] - h_ij[i][j][k-1];
//                } // fi
//            } // next k
//        }; // +dh_ij
//        // Thermal response factor increments
//        for (int i=0; i<h_ij.size(); i++) {
//            for (int j=0; j<h_ij[0].size(); j++) {
////                boost::asio::post(pool3, [&_dh_ij, i, j]{ _dh_ij(i, j) ;});
//                _dh_ij(i, j);
//            } // next j
//        } // next i

        // after threaded, join() pool here
//        pool3.join();

//        h_ij = vector<vector<vector<double> >> ();  // destruct h_ij; no longer needed

        end = std::chrono::steady_clock::now();
        milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        segment_h_values_time += milli;

        // after interpolation scheme, get rid of h_ij first
        // Initialize segment heat extraction rates
        vector<vector<double> > Q(nSources, vector<double> (nt));

        // Define A and b for utitilizing Ax=b
        /**
         * A = [    [   ],
         *          [hb[0:len(hb), 0]
         *     ]
         * b = [    [  ],
         *          [sum(hb)]
         *     ]
         * **/

        int SIZE = nSources + 1;
        // _gesv initializiation
        int nrhs = 1; // number of columns in the b Matrix
        int lda = SIZE;
        int ldb = SIZE;
        std::vector<int> _ipiv(SIZE);
        int info;
        vector<double> A_ (SIZE * SIZE);
        vector<double> b_ (SIZE);

        std::vector<std::vector <double>> A(nSources +1, std::vector<double> (nSources +1));
        std::vector<double> b(nSources + 1);

        // Fill A
        int n = SIZE - 1;
//        for (int j=0; j<SIZE; j++) {
//            if (j==n) {
//                gsl_matrix_set (_A, n, n, 0);
//                A[n][n] = 0;
//            } else {
//                gsl_matrix_set (_A, n, j, Hb[j]);
//                A[n][j] = Hb[j];
//            } // fi
//        } // next j

        // Energy conservation: sum([Qb*Hb)] = sum([Hb])
//        int n = A.size() - 1;
//        int n = SIZE -1;
//        _fill_A(n, SIZE);
        n = b.size() - 1;
        double Hb_sum=0;
        for (auto & _hb : Hb) {
            Hb_sum += _hb;
        }
        b[n] = Hb_sum;

        // Build and solve the system of equations at all times

        // the loop p=n depends on what occured at p=n-1, so this will be be in series
        // however threading will be interspersed throughout to make use of as many threads as possible
        // TODO: since this is in series, move the variable declarations here
//        vector<vector<double>> h_ij_dt (nSources, vector<double> (nSources));
        std::vector<double> Tb_0 (nSources);
        // Restructured load history
        // create interpolation object for accumulated heat extraction
        std::vector<std::vector<double>> q_reconstructed (nSources, std::vector<double> (nt));

        for (int p=0; p<nt; p++) {
            // current thermal response factor matrix
            auto _fill_h_ij_dt = [&h_dt, &A] (const int i, const int p) {
                int m = h_dt[0].size();
                for (int j=0; j<A[i].size(); j++) {
                    if (j==A[i].size()-1) {
                        A[i][j] = -1;
                    } else {
                        A[j][i] = h_dt[i][j][p];
                    }
//                    h_ij_dt[j][i] = h_dt[i][j][p];

                } // next j
                ;
            }; // _fill_h_ij_dt
            // _fill_A
            //        auto _fill_A = [&A, &Hb, &_A](const int i, const int SIZE) { // TODO: keep in mind this function can make use of threading

            // ------------- fill A ------------
            start = std::chrono::steady_clock::now();
            auto _fillA = [&Hb, &h_dt, &A, &A_, &h_ij, &dt, &_time_untouched](int i, int p, int SIZE) {
                int n = SIZE - 1;
                vector<double> yp(1);
                vector<double> xp(1);
                for (int j=0; j<SIZE; j++) {
                    if (i == n) { // then we are referring to Hb
                        if (j==n) {
                            A_[i+j*SIZE] = 0;
//                            A[i][n] = 0;
                        } else {
                            A_[i+j*SIZE] = Hb[j];
//                            A[i][j] = Hb[j];
                        } // fi
                    } else {
                        if (j==A[i].size()-1) {
                            A_[i+j*SIZE] = -1;
//                            A[i][j] = -1;
                        } else {
                            if (i==7 && j==7 &&p==4) {
                                int a = 1;
                            }
                            xp[0] = dt[p];
                            vector<double> x = _time_untouched;
                            vector<double> y = h_ij[i][j];
//                            cout << "i: " << i << " j: " << j << " p: " << p << endl;

                            if (yp.size() > 1) {
                                int a = 1;
                            }
                            jcc::interpolation::interp1d(xp, yp, _time_untouched, y );
                            double a = h_dt[i][j][p];
//                            A_[j+i*SIZE] = yp[0];
                            A_[j+i*SIZE] = h_dt[i][j][p];
//                            A[j][i] = h_dt[i][j][p];
                        } // fi
                    } // fi
                } // next k
            };
            boost::asio::thread_pool pool3(processor_count);
            // A needs filled each loop because the _gsl partial pivot decomposition modifies the matrix
            for (int i=0; i<SIZE; i++) {
//                boost::asio::post(pool3, [&_fillA, i, p, SIZE]{ _fillA(i, p, SIZE) ;});
                _fillA(i, p, SIZE);
            }
            pool3.join();

            end = std::chrono::steady_clock::now();  // _fill_A
            milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            fill_A_time += milli;

            // ----- load history reconstruction -------
            start = std::chrono::steady_clock::now();
            load_history_reconstruction(q_reconstructed,time, _time, Q, dt, p);
            end = std::chrono::steady_clock::now();
            milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            load_history_reconstruction_time += milli;

            // ----- temporal superposition
            start = std::chrono::steady_clock::now();
            _temporal_superposition(Tb_0, dh_ij, h_ij, q_reconstructed, p);
            // fill b with -Tb
            for (int i=0; i<Tb_0.size(); i++) {
                b[i] = -Tb_0[i];
            }
            end = std::chrono::steady_clock::now();
            milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            temporal_superposition_time += milli;

            int m = A.size();
            int n = A[0].size();
            vector<double> x(b.size());
//            _solve_eqn(x, A, b);
            /** was _solve_eqn **/

            // ---- fill gsl matrix A and b -----
            start = std::chrono::steady_clock::now();

            for (int i=0; i<b.size(); i++) {
                b_[i] = b[i];
            } // next i

            end = std::chrono::steady_clock::now();
            milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            fill_gsl_matrices_time += milli;

            // ----- LU decomposition -----
            start = std::chrono::steady_clock::now();
            dgesv_( &n, &nrhs, &*A_.begin(), &lda, &*_ipiv.begin(), &*b_.begin(), &ldb, &info );

            for (int i=0; i<SIZE; i++) {
                x[i] = b_[i];
            } // next i
            end = std::chrono::steady_clock::now();
            milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            LU_decomposition_time += milli;

            // ---- Save Q's for next p ---
            for (int j=0; j<Q.size(); j++) {
                Q[j][p] = x[j];
            } // next j
            // the borehole wall temperatures are equal for all segments
            double Tb = x[x.size()-1];
            gfunction[p] = Tb;
        } // next p
        segment_length_time /= 1000;
        time_vector_time /= 1000;
        segment_h_values_time /= 1000;
        segment_length_time /= 1000;
        fill_A_time /= 1000;
        load_history_reconstruction_time /= 1000;
        temporal_superposition_time /= 1000;
        fill_gsl_matrices_time /= 1000;
        LU_decomposition_time /= 1000;

        cout << "------ timings report -------" << endl;
        cout << " t\t " << " t/p\t" << "name" << endl;
        cout << segment_length_time << "\t" << segment_length_time << "\t" << "segment length time" << endl;
        cout << time_vector_time << "\t" << time_vector_time << "\t" << "time vector time" << endl;
        cout << segment_h_values_time << "\t" << segment_h_values_time << "\t" << "segment h values time" << endl;
        cout << fill_A_time << "\t" << fill_A_time / double(nt) << "\t" << "time to fill vector A" << endl;
        cout << load_history_reconstruction_time << "\t" << load_history_reconstruction_time / double(nt)
        << "\t" << "load hist reconstruction" << endl;
        cout << temporal_superposition_time << "\t" << temporal_superposition_time / double(nt)
             << "\t" << "temporal superposition time:" << endl;
        cout << fill_gsl_matrices_time << "\t" << fill_gsl_matrices_time / double(nt)
             << "\t" << "gsl fill matrices time" << endl;
        cout << LU_decomposition_time << "\t" << LU_decomposition_time/double(nt)
             << "\t" << "LU decomp time" << endl;

        auto end2 = std::chrono::steady_clock::now();
        if (disp) {
            double milli1 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
            double seconds1 = milli1 / 1000;
            double milli2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - startall).count();
            double seconds2 = milli2 / 1000;
            std::cout << "Elapsed time in seconds : "
                      << seconds1
                      << " sec" << std::endl;
            std::cout << "Total time for g-function evaluation : "
                    << seconds2
                    << " sec" << std::endl;
        }
        int a = 1;
    } // void uniform_temperature

    void _borehole_segments(std::vector<gt::boreholes::Borehole>& boreSegments,
            std::vector<gt::boreholes::Borehole>& boreholes, const int nSegments) {
        double H;
        double D;
        int count = 0;
        // Split boreholes into segments
        for(auto& b : boreholes) {
            // TODO: maybe thread this later on
            for (int i=0; i<nSegments; i++) {
                H = b.H / double(nSegments);
                D = b.D + double(i) * b.H / double(nSegments);
                boreSegments[count] = gt::boreholes::Borehole(H, D, b.r_b, b.x, b.y);
                count++;
            }  // end for
        } // end for
    } // void _borehole_segments

    void load_history_reconstruction(std::vector<std::vector<double>>& q_reconstructed,
            vector<double>& time, vector<double>& _time, vector<vector<double> >& Q,
            vector<double>& dt, const int p) {
        // for s in range p+1
        int nSources = Q.size();

        // Inverted time steps
        std::vector<double> dt_reconstructed (p+1);
        for (int i=p; i>=0; i--) {
            dt_reconstructed[p-i] = dt[i];  // reverse the dt
        }
        // t_restructured is [0, cumsum(dt_reversed)]
        std::vector<double> t_reconstructed(p+2);  // will start at 0
        for (int i=1; i<=p+1; i++) {
            t_reconstructed[i] = dt_reconstructed[i-1];
        }
        for (int i=1; i<=p+1; i++) {
            t_reconstructed[i] = t_reconstructed[i] + t_reconstructed[i-1];
        }
         // local time vector
        std::vector<double> t(p+3);
        for (int i=0; i<t.size(); i++) {
            if (i==t.size()-1) {
                t[i] = _time[i-1] + _time[1];
            } else {
                t[i] = _time[i];
            }
        }
        int _tsize = t.size();
        // Q*dt
        std::vector<std::vector <double> > Q_dt (nSources, std::vector<double> (t.size()));
        auto _Q_dot_dt = [&Q_dt, &Q, &dt, &p, &_tsize, &t](const int i) {
            for (int j = 1; j<_tsize; j++) {
                if (j>=p+1) {
                    Q_dt[i][j] = Q_dt[i][j-1];
                } else {
                    Q_dt[i][j] = Q[i][j-1] * dt[j-1] + Q_dt[i][j-1];
                }  // fi
            } // next j
        };
        for (int i=0; i<nSources; i++) {
            _Q_dot_dt(i);  // could be threaded here, if timings ever prove necessary
        } // next i

        auto _interpolate = [&Q_dt, &q_reconstructed, &t, &t_reconstructed, &dt_reconstructed, &p](const int i) {
            int n = t.size();
            std::vector<double> y(n);
            for (int j=0; j<n; j++) {
                y[j] = Q_dt[i][j];
            }

            std::vector<double> yp(n);
            jcc::interpolation::interp1d(t_reconstructed, yp, t, y);

            for (int j=0; j<p; j++) {
                double c = yp[j];
                double d = yp[j+1];
                double e = dt_reconstructed[j];
//                q_reconstructed[i][j] = (s(t_reconstructed[j+1]) - s(t_reconstructed[j])) / dt_reconstructed[j];
                q_reconstructed[i][j] = (d - c) / e;
            }
        }; // _interpolate
        for (int i=0; i<nSources; i++) {
            _interpolate(i); // could be threaded here, but this function doesn't take long at all
        }
    } // load_history_reconstruction

    void _temporal_superposition(vector<double>& Tb_0, vector<vector<vector<double> > >& h_ij,
                                 vector<vector<vector<double> > >& dh_ij,
                                 std::vector<std::vector<double>>& q_reconstructed, int p)
            {
        const auto processor_count = thread::hardware_concurrency();
        // Launch the pool with n threads.
        boost::asio::thread_pool pool(processor_count);

        std::fill(Tb_0.begin(), Tb_0.end(), 0);
        // Number of heat sources
        int nSources = q_reconstructed.size();
        // Number of time steps
        int nt = p + 1;

        auto _borehole_wall_temp = [&h_ij, &q_reconstructed, &Tb_0, dh_ij](const int i, const int nSources, const int nt){
            for (int j =0; j<nSources; j++) {
                for (int k=0; k<nt; k++) {
//                    if (k==0) {
//                        Tb_0[i] += h_ij[i][j][k+1] * q_reconstructed[j][nt-k-1] ;
//                    } else if (k>0) {
//                        Tb_0[i] += (h_ij[i][j][k+1]- h_ij[i][j][k]) * q_reconstructed[j][nt-k-1] ;
//                    }
                    Tb_0[i] += dh_ij[i][j][k] * q_reconstructed[j][nt-k-1] ;

                }
            }
        };
        for (int i=0; i<nSources; i++) {
//            boost::asio::post(pool, [&_borehole_wall_temp, i, nSources, nt]
//            { _borehole_wall_temp(i, nSources, nt); });
            _borehole_wall_temp(i, nSources, nt);
        }

        pool.join();
    }  // _temporal_superposition

} // namespace gt::gfunction
