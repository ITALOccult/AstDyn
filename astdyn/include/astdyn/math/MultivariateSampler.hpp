#ifndef ASTDYN_MATH_MULTIVARIATE_SAMPLER_HPP
#define ASTDYN_MATH_MULTIVARIATE_SAMPLER_HPP

#include <Eigen/Dense>
#include <random>
#include <vector>

namespace astdyn::math {

/**
 * @brief Sampler for multivariate normal distributions using Cholesky decomposition.
 */
class MultivariateNormalSampler {
public:
    using MatrixXd = Eigen::MatrixXd;
    using VectorXd = Eigen::VectorXd;

    /**
     * @brief Initialize with mean and covariance matrix.
     * @param mean Vector representing the mean values (size N)
     * @param cov Covariance matrix (size N x N)
     */
    MultivariateNormalSampler(const VectorXd& mean, const MatrixXd& cov)
        : mean_(mean), gen_(std::random_device{}())
    {
        // LL^T = P (Cholesky decomposition)
        Eigen::LLT<MatrixXd> llt(cov);
        if (llt.info() != Eigen::Success) {
            throw std::runtime_error("MultivariateNormalSampler: covariance matrix is not positive-definite");
        }
        L_ = llt.matrixL();
        dist_ = std::normal_distribution<double>(0.0, 1.0);
    }

    /**
     * @brief Samples a single point from the distribution.
     * x = mu + L*z, where z ~ N(0, I)
     */
    VectorXd sample() {
        VectorXd z(mean_.size());
        for (int i = 0; i < z.size(); ++i) {
            z(i) = dist_(gen_);
        }
        return mean_ + L_ * z;
    }

    /**
     * @brief Generates N samples.
     */
    std::vector<VectorXd> samples(int n) {
        std::vector<VectorXd> result;
        result.reserve(n);
        for (int i = 0; i < n; ++i) {
            result.push_back(sample());
        }
        return result;
    }

private:
    VectorXd mean_;
    MatrixXd L_;
    std::mt19937 gen_;
    std::normal_distribution<double> dist_;
};

} // namespace astdyn::math

#endif // ASTDYN_MATH_MULTIVARIATE_SAMPLER_HPP
