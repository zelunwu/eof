#include "eof/eof.hpp"

#include <cmath>
#include <iostream>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

int g_failed = 0;
int g_passed = 0;

void check(bool cond, const char* msg) {
    if (cond) {
        ++g_passed;
    } else {
        ++g_failed;
        std::cerr << "FAIL: " << msg << "\n";
    }
}

}  // namespace

int main() {
    using eof::Matrix;

    {
        const int n_loc = 10;
        const int n_time = 40;
        Matrix X(n_loc, n_time);
        for (int i = 0; i < n_loc; ++i) {
            for (int t = 0; t < n_time; ++t) {
                X(i, t) = std::sin(2.0 * M_PI * i / n_loc) *
                              std::cos(2.0 * M_PI * t / n_time) +
                          0.25 * std::cos(2.0 * M_PI * i / n_loc) *
                              std::sin(4.0 * M_PI * t / n_time);
            }
        }
        eof::PackedEof solver(X, 2);
        Matrix Xa = X;
        Xa.colwise() -= X.rowwise().mean();
        const double rel = (solver.reconstruct() - Xa).norm() / (Xa.norm() + 1e-15);
        check(rel < 1e-10, "rank-2 reconstruction");
        Matrix gram = solver.modes().transpose() * solver.modes();
        check((gram - Matrix::Identity(2, 2)).norm() < 1e-8, "orthonormal EOFs");
        check(solver.expvar()(0) > solver.expvar()(1), "ordered expvar");
    }

    {
        eof::PackedEof solver(Matrix::Random(6, 15));
        check(std::abs(solver.expvar().sum() - 100.0) < 1e-6, "expvar sums to 100");
    }

    {
        const int n_lon = 4, n_lat = 3, n_time = 50;
        std::vector<double> data(static_cast<size_t>(n_lon * n_lat * n_time), 0.0);
        const int nsp = n_lon * n_lat;
        for (int t = 0; t < n_time; ++t) {
            for (int j = 0; j < n_lat; ++j) {
                for (int i = 0; i < n_lon; ++i) {
                    const int s = i + j * n_lon;
                    double v = 0.0;
                    if (j == 1) {
                        v = std::sin(2.0 * M_PI * t / 20.0);
                    } else if (j == 2) {
                        v = std::sin(2.0 * M_PI * t / 13.0 + 0.4);
                    }
                    data[static_cast<size_t>(s + t * nsp)] = v;
                }
            }
        }
        for (int t = 0; t < n_time; ++t) {
            data[static_cast<size_t>(t * nsp)] = std::nan("");
        }
        double lat[3] = {-75.0, 0.0, 75.0};
        eof::GridOptions opt;
        opt.n_eof = 2;
        eof::GridEof weighted(data.data(), n_lon, n_lat, n_time, lat, 3, opt);
        eof::GridEof unweighted(data.data(), n_lon, n_lat, n_time, nullptr, 0, opt);
        auto maps = weighted.maps();
        check(std::isnan(maps[0]), "NaN mask preserved");
        check(std::isfinite(maps[1]), "ocean cell finite");
        const double ratio_w = weighted.expvar()(0) / weighted.expvar()(1);
        const double ratio_u = unweighted.expvar()(0) / unweighted.expvar()(1);
        check(ratio_w > ratio_u * 1.05, "area weight prefers equator");
    }

    {
        const int n_lon = 3, n_lat = 2, n_time = 40;
        const int nsp = n_lon * n_lat;
        std::vector<double> d1(static_cast<size_t>(nsp * n_time));
        std::vector<double> d2(static_cast<size_t>(nsp * n_time));
        for (int t = 0; t < n_time; ++t) {
            const double a = std::sin(2.0 * M_PI * t / 10.0);
            const double b = std::cos(2.0 * M_PI * t / 10.0);
            for (int p = 0; p < nsp; ++p) {
                d1[static_cast<size_t>(p + t * nsp)] = a;
                d2[static_cast<size_t>(p + t * nsp)] = 100.0 * b;
            }
        }
        eof::SpatialWeights ones(n_lon, n_lat);
        eof::MultivariateEof raw(d1.data(), n_lon, n_lat, d2.data(), n_lon, n_lat, n_time,
                                 ones, ones, false, 2);
        eof::MultivariateEof norm(d1.data(), n_lon, n_lat, d2.data(), n_lon, n_lat, n_time,
                                  ones, ones, true, 2);
        check(raw.expvar()(0) / raw.expvar()(1) > 10.0, "MEOF raw variance imbalance");
        check(std::abs(norm.expvar()(0) - norm.expvar()(1)) < 5.0,
              "MEOF normalize balances variance");
    }

    std::cout << "passed=" << g_passed << " failed=" << g_failed << "\n";
    return g_failed == 0 ? 0 : 1;
}
