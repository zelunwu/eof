#include "eof/eof.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace eof {
namespace {

void pin_signs(Matrix& L, Matrix& Y) {
    for (int k = 0; k < L.cols(); ++k) {
        Eigen::Index idx = 0;
        L.col(k).cwiseAbs().maxCoeff(&idx);
        if (L(idx, k) < 0.0) {
            L.col(k) *= -1.0;
            Y.row(k) *= -1.0;
        }
    }
}

std::pair<Vector, Matrix> eigensolve(const Matrix& S_in, int n_eof) {
    Matrix S = 0.5 * (S_in + S_in.transpose());
    Eigen::SelfAdjointEigenSolver<Matrix> es(S);
    if (es.info() != Eigen::Success) {
        throw std::runtime_error("eof: eigendecomposition failed");
    }
    const Vector evals = es.eigenvalues();
    const Matrix evecs = es.eigenvectors();
    const int n = static_cast<int>(evals.size());
    n_eof = std::min(n_eof, n);
    Vector out_val(n_eof);
    Matrix out_vec(n, n_eof);
    for (int i = 0; i < n_eof; ++i) {
        const int src = n - 1 - i;
        out_val(i) = std::max(evals(src), 0.0);
        out_vec.col(i) = evecs.col(src);
    }
    return {out_val, out_vec};
}

std::vector<int> finite_locations(const double* data, int n_space, int n_time) {
    std::vector<int> keep;
    keep.reserve(n_space);
    for (int s = 0; s < n_space; ++s) {
        bool ok = true;
        for (int t = 0; t < n_time; ++t) {
            if (!std::isfinite(data[s + t * n_space])) {
                ok = false;
                break;
            }
        }
        if (ok) {
            keep.push_back(s);
        }
    }
    return keep;
}

Matrix pack_weighted(const double* data, int n_lon, int n_lat, int n_time,
                     const SpatialWeights& W, const std::vector<int>& keep) {
    const int n_space = n_lon * n_lat;
    Matrix X(static_cast<int>(keep.size()), n_time);
    for (int t = 0; t < n_time; ++t) {
        for (int r = 0; r < static_cast<int>(keep.size()); ++r) {
            const int s = keep[static_cast<size_t>(r)];
            X(r, t) = data[s + t * n_space] * W.at(s % n_lon, s / n_lon);
        }
    }
    return X;
}

double field_std(const Matrix& d) {
    Matrix a = d;
    a.colwise() -= a.rowwise().mean();
    const double n = static_cast<double>(a.size());
    if (n <= 1.0) {
        return 0.0;
    }
    return std::sqrt(a.squaredNorm() / (n - 1.0));
}

void unweight_rows(Matrix& L, const SpatialWeights& W, const std::vector<int>& keep) {
    const int n_lon = W.n_lon();
    for (int r = 0; r < static_cast<int>(keep.size()); ++r) {
        const int s = keep[static_cast<size_t>(r)];
        const double w = W.at(s % n_lon, s / n_lon);
        const double den = (w == 0.0) ? std::numeric_limits<double>::infinity() : w;
        for (int k = 0; k < L.cols(); ++k) {
            double v = L(r, k) / den;
            if (!std::isfinite(v)) {
                v = 0.0;
            }
            L(r, k) = v;
        }
    }
}

std::vector<double> scatter_maps(const Matrix& L, const std::vector<int>& keep,
                                 int n_lon, int n_lat) {
    const int n_space = n_lon * n_lat;
    const int n_eof = static_cast<int>(L.cols());
    std::vector<double> maps(static_cast<size_t>(n_space * n_eof),
                             std::numeric_limits<double>::quiet_NaN());
    for (int k = 0; k < n_eof; ++k) {
        for (int r = 0; r < static_cast<int>(keep.size()); ++r) {
            maps[static_cast<size_t>(keep[static_cast<size_t>(r)] + k * n_space)] = L(r, k);
        }
    }
    return maps;
}

}  // namespace

SpatialWeights::SpatialWeights(int n_lon, int n_lat)
    : n_lon_(n_lon), n_lat_(n_lat), W_(Matrix::Ones(n_lon, n_lat)) {
    if (n_lon < 1 || n_lat < 1) {
        throw std::invalid_argument("SpatialWeights: invalid grid size");
    }
}

SpatialWeights::SpatialWeights(int n_lon, const Vector& lat_deg, const std::string& method)
    : n_lon_(n_lon), n_lat_(static_cast<int>(lat_deg.size())), W_(n_lon, lat_deg.size()) {
    if (n_lon < 1 || n_lat_ < 1) {
        throw std::invalid_argument("SpatialWeights: invalid grid size");
    }
    for (int j = 0; j < n_lat_; ++j) {
        if (lat_deg(j) > 90.0 || lat_deg(j) < -90.0) {
            throw std::invalid_argument("latitude must lie in [-90, 90]");
        }
        const double c = std::cos(lat_deg(j) * M_PI / 180.0);
        double w = 1.0;
        if (method == "sqrtcos") {
            w = std::sqrt(std::max(c, 0.0));
        } else if (method == "cos") {
            w = std::max(c, 0.0);
        } else if (method == "none") {
            w = 1.0;
        } else {
            throw std::invalid_argument("unknown weighting method");
        }
        W_.col(j).setConstant(w);
    }
}

SpatialWeights::SpatialWeights(int n_lon, int n_lat, const double* weights, int n_weights)
    : n_lon_(n_lon), n_lat_(n_lat), W_(n_lon, n_lat) {
    if (weights == nullptr) {
        throw std::invalid_argument("SpatialWeights: null weights");
    }
    const int n_space = n_lon * n_lat;
    if (n_weights == n_lat) {
        for (int j = 0; j < n_lat; ++j) {
            W_.col(j).setConstant(weights[j]);
        }
    } else if (n_weights == n_space) {
        for (int j = 0; j < n_lat; ++j) {
            for (int i = 0; i < n_lon; ++i) {
                W_(i, j) = weights[i + j * n_lon];
            }
        }
    } else {
        throw std::invalid_argument("SpatialWeights: bad weights length");
    }
}

PackedEof::PackedEof(Matrix X, int n_eof) {
    if (X.rows() < 1 || X.cols() < 1) {
        throw std::invalid_argument("PackedEof: X must be non-empty");
    }
    const int n_locations = static_cast<int>(X.rows());
    const int n_timesteps = static_cast<int>(X.cols());
    X.colwise() -= X.rowwise().mean();
    if (n_eof < 0) {
        n_eof = std::min(n_locations, n_timesteps);
    }
    n_eof = std::min(n_eof, std::min(n_locations, n_timesteps));
    if (n_eof < 1) {
        throw std::invalid_argument("PackedEof: n_eof must be >= 1");
    }
    const double total_var = X.squaredNorm() / static_cast<double>(n_timesteps);
    if (n_locations <= n_timesteps) {
        auto pair = eigensolve(X * X.transpose(), n_eof);
        eig_ = pair.first / static_cast<double>(n_timesteps);
        L_ = std::move(pair.second);
        Y_ = (X.transpose() * L_).transpose();
    } else {
        auto pair = eigensolve(X.transpose() * X, n_eof);
        eig_ = pair.first / static_cast<double>(n_timesteps);
        Vector scale(n_eof);
        const double tiny = std::numeric_limits<double>::min();
        for (int k = 0; k < n_eof; ++k) {
            scale(k) = 1.0 / std::sqrt(std::max(eig_(k) * n_timesteps, tiny));
        }
        L_ = X * (pair.second * scale.asDiagonal());
        Y_ = L_.transpose() * X;
    }
    pin_signs(L_, Y_);
    expvar_ = Vector::Zero(n_eof);
    if (total_var > 0.0) {
        expvar_ = eig_ / total_var * 100.0;
    }
}

Matrix PackedEof::reconstruct(int n_eof) const {
    if (n_eof < 0 || n_eof > L_.cols()) {
        n_eof = static_cast<int>(L_.cols());
    }
    return L_.leftCols(n_eof) * Y_.topRows(n_eof);
}

GridEof::GridEof(const double* data, int n_lon, int n_lat, int n_time,
                 const SpatialWeights& weights, int n_eof)
    : n_lon_(n_lon),
      n_lat_(n_lat),
      n_time_(n_time),
      weights_(weights),
      solver_(Matrix::Zero(1, 1), 1) {
    if (weights.n_lon() != n_lon || weights.n_lat() != n_lat) {
        throw std::invalid_argument("GridEof: weight grid mismatch");
    }
    build(data, n_eof);
}

GridEof::GridEof(const double* data, int n_lon, int n_lat, int n_time,
                 const double* lat, int n_lat_len, const GridOptions& opt)
    : n_lon_(n_lon),
      n_lat_(n_lat),
      n_time_(n_time),
      weights_(n_lon, n_lat),
      solver_(Matrix::Zero(1, 1), 1) {
    if (lat != nullptr && n_lat_len > 0) {
        if (n_lat_len != n_lat) {
            throw std::invalid_argument("GridEof: len(lat) must equal n_lat");
        }
        Vector latv = Eigen::Map<const Vector>(lat, n_lat);
        weights_ = SpatialWeights(n_lon, latv, opt.weighting);
    }
    build(data, opt.n_eof);
}

void GridEof::build(const double* data, int n_eof) {
    if (data == nullptr || n_lon_ < 1 || n_lat_ < 1 || n_time_ < 1) {
        throw std::invalid_argument("GridEof: invalid input");
    }
    keep_ = finite_locations(data, n_lon_ * n_lat_, n_time_);
    if (keep_.empty()) {
        throw std::invalid_argument("GridEof: no finite spatial locations");
    }
    Matrix X = pack_weighted(data, n_lon_, n_lat_, n_time_, weights_, keep_);
    solver_ = PackedEof(std::move(X), n_eof);
}

std::vector<double> GridEof::maps(bool unweight, bool standardize) const {
    Matrix L = solver_.modes();
    if (unweight) {
        unweight_rows(L, weights_, keep_);
    }
    if (standardize) {
        for (int k = 0; k < L.cols(); ++k) {
            const double stdn = std::sqrt(std::abs(solver_.eigenvalues()(k)));
            if (stdn > 0.0) {
                L.col(k) *= stdn;
            }
        }
    }
    return scatter_maps(L, keep_, n_lon_, n_lat_);
}

Matrix GridEof::pcs(bool standardize) const {
    Matrix Y = solver_.pcs();
    if (standardize) {
        for (int k = 0; k < Y.rows(); ++k) {
            const double stdn = std::sqrt(std::abs(solver_.eigenvalues()(k)));
            if (stdn > 0.0) {
                Y.row(k) /= stdn;
            }
        }
    }
    return Y;
}

PackedEof MultivariateEof::solve(const double* data1, int n_lon1, int n_lat1,
                                 const double* data2, int n_lon2, int n_lat2, int n_time,
                                 const SpatialWeights& w1, const SpatialWeights& w2,
                                 bool normalize, int n_eof,
                                 std::vector<int>& keep1, std::vector<int>& keep2,
                                 double& s1, double& s2) {
    keep1 = finite_locations(data1, n_lon1 * n_lat1, n_time);
    keep2 = finite_locations(data2, n_lon2 * n_lat2, n_time);
    Matrix X1 = pack_weighted(data1, n_lon1, n_lat1, n_time, w1, keep1);
    Matrix X2 = pack_weighted(data2, n_lon2, n_lat2, n_time, w2, keep2);
    s1 = s2 = 1.0;
    if (normalize) {
        s1 = field_std(X1);
        s2 = field_std(X2);
        if (s1 > 0.0) {
            X1 /= s1;
        }
        if (s2 > 0.0) {
            X2 /= s2;
        }
    }
    Matrix X(X1.rows() + X2.rows(), n_time);
    X.topRows(X1.rows()) = X1;
    X.bottomRows(X2.rows()) = X2;
    return PackedEof(std::move(X), n_eof);
}

MultivariateEof::MultivariateEof(const double* data1, int n_lon1, int n_lat1,
                                 const double* data2, int n_lon2, int n_lat2, int n_time,
                                 const SpatialWeights& w1, const SpatialWeights& w2,
                                 bool normalize, int n_eof)
    : n_lon1_(n_lon1),
      n_lat1_(n_lat1),
      n_lon2_(n_lon2),
      n_lat2_(n_lat2),
      n_time_(n_time),
      w1_(w1),
      w2_(w2),
      normalize_(normalize),
      solver_(solve(data1, n_lon1, n_lat1, data2, n_lon2, n_lat2, n_time, w1, w2,
                    normalize, n_eof, keep1_, keep2_, s1_, s2_)) {}

std::vector<double> MultivariateEof::maps1(bool unweight) const {
    Matrix L = solver_.modes().topRows(static_cast<int>(keep1_.size()));
    if (normalize_) {
        L *= s1_;
    }
    if (unweight) {
        unweight_rows(L, w1_, keep1_);
    }
    return scatter_maps(L, keep1_, n_lon1_, n_lat1_);
}

std::vector<double> MultivariateEof::maps2(bool unweight) const {
    Matrix L = solver_.modes().bottomRows(static_cast<int>(keep2_.size()));
    if (normalize_) {
        L *= s2_;
    }
    if (unweight) {
        unweight_rows(L, w2_, keep2_);
    }
    return scatter_maps(L, keep2_, n_lon2_, n_lat2_);
}

}  // namespace eof
