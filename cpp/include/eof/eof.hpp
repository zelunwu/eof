#pragma once

/// Object-oriented EOF analysis matching the MATLAB / Python toolbox.
///
/// Layout is column-major (Eigen / MATLAB):
///   packed X is n_locations x n_timesteps
///   a grid is (n_lon, n_lat, n_time) with space index = lon + lat * n_lon.

#include <Eigen/Dense>
#include <string>
#include <vector>

namespace eof {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

/// Multiplicative weights on an (n_lon x n_lat) grid.
class SpatialWeights {
public:
    SpatialWeights(int n_lon, int n_lat);
    SpatialWeights(int n_lon, const Vector& lat_deg, const std::string& method = "sqrtcos");
    SpatialWeights(int n_lon, int n_lat, const double* weights, int n_weights);

    int n_lon() const { return n_lon_; }
    int n_lat() const { return n_lat_; }
    const Matrix& array() const { return W_; }
    double at(int lon, int lat) const { return W_(lon, lat); }

private:
    int n_lon_;
    int n_lat_;
    Matrix W_;
};

/// EOF of a packed (n_locations x n_timesteps) matrix.
class PackedEof {
public:
    explicit PackedEof(Matrix X, int n_eof = -1);

    int n_locations() const { return static_cast<int>(L_.rows()); }
    int n_timesteps() const { return static_cast<int>(Y_.cols()); }
    int n_eof() const { return static_cast<int>(L_.cols()); }

    const Matrix& modes() const { return L_; }
    const Matrix& pcs() const { return Y_; }
    const Vector& eigenvalues() const { return eig_; }
    const Vector& expvar() const { return expvar_; }
    Matrix reconstruct(int n_eof = -1) const;

private:
    Matrix L_;
    Matrix Y_;
    Vector eig_;
    Vector expvar_;
};

struct GridOptions {
    int n_eof = -1;
    bool unweight = true;
    bool standardize = false;
    std::string weighting = "sqrtcos";
};

/// EOF of a (lon, lat, time) field stored column-major.
class GridEof {
public:
    GridEof(const double* data, int n_lon, int n_lat, int n_time,
            const SpatialWeights& weights, int n_eof = -1);

    GridEof(const double* data, int n_lon, int n_lat, int n_time,
            const double* lat, int n_lat_len,
            const GridOptions& opt = {});

    int n_lon() const { return n_lon_; }
    int n_lat() const { return n_lat_; }
    int n_time() const { return n_time_; }
    int n_eof() const { return solver_.n_eof(); }

    const PackedEof& packed() const { return solver_; }
    const SpatialWeights& weights() const { return weights_; }
    const Vector& eigenvalues() const { return solver_.eigenvalues(); }
    const Vector& expvar() const { return solver_.expvar(); }

    /// n_lon * n_lat * n_eof, NaN on masked cells. lon fastest, then lat, then mode.
    std::vector<double> maps(bool unweight = true, bool standardize = false) const;
    Matrix pcs(bool standardize = false) const;

private:
    void build(const double* data, int n_eof);

    int n_lon_;
    int n_lat_;
    int n_time_;
    SpatialWeights weights_;
    std::vector<int> keep_;
    PackedEof solver_;
};

/// Combined EOF of two grids that share n_time.
class MultivariateEof {
public:
    MultivariateEof(const double* data1, int n_lon1, int n_lat1,
                    const double* data2, int n_lon2, int n_lat2, int n_time,
                    const SpatialWeights& w1, const SpatialWeights& w2,
                    bool normalize = true, int n_eof = -1);

    int n_eof() const { return solver_.n_eof(); }
    const Vector& eigenvalues() const { return solver_.eigenvalues(); }
    const Vector& expvar() const { return solver_.expvar(); }
    const PackedEof& packed() const { return solver_; }

    std::vector<double> maps1(bool unweight = true) const;
    std::vector<double> maps2(bool unweight = true) const;
    const Matrix& pcs() const { return solver_.pcs(); }

private:
    static PackedEof solve(const double* data1, int n_lon1, int n_lat1,
                           const double* data2, int n_lon2, int n_lat2, int n_time,
                           const SpatialWeights& w1, const SpatialWeights& w2,
                           bool normalize, int n_eof,
                           std::vector<int>& keep1, std::vector<int>& keep2,
                           double& s1, double& s2);

    int n_lon1_, n_lat1_, n_lon2_, n_lat2_, n_time_;
    SpatialWeights w1_;
    SpatialWeights w2_;
    bool normalize_;
    double s1_ = 1.0;
    double s2_ = 1.0;
    std::vector<int> keep1_;
    std::vector<int> keep2_;
    PackedEof solver_;
};

}  // namespace eof
