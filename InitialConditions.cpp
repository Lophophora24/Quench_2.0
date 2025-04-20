#include "InitialConditions.h"

double phi_0[SIZE_X];
double pi_0[SIZE_X];

double phi_averaged_with_eta = 0;

/*
  We need a functor that can pretend it's const,
  but to be a good random number generator
  it needs mutable state.
*/
namespace Eigen {
    namespace internal {
        template<typename Scalar>
        struct scalar_normal_dist_op
        {
            static boost::mt19937 rng;    // The uniform pseudo-random algorithm
            mutable boost::normal_distribution<Scalar> norm;  // The gaussian combinator

            EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

                template<typename Index>
            inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
        };

        template<typename Scalar> boost::mt19937 scalar_normal_dist_op<Scalar>::rng;

        template<typename Scalar>
        struct functor_traits<scalar_normal_dist_op<Scalar> >
        {
            enum { Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false };
        };
    } // end namespace internal
} // end namespace Eigen

int size = SIZE_X - 1; // Dimensionality (rows)
int nn = M;     // How many samples (columns) to draw
// Define mean and covariance of the distribution
Eigen::VectorXd mean_phi(size);
Eigen::MatrixXd covar_phi(size, size);
// Define mean and covariance of the distribution
Eigen::VectorXd mean_pi(size);
Eigen::MatrixXd covar_pi(size, size);

Eigen::MatrixXd A(size, size);
Eigen::MatrixXd B(size, size);

void init_()
{
    for (int i = 0; i < SIZE_X-1; ++i) {
        mean_phi(i) = 0;
        mean_pi(i) = 0;
        for (int j = 0; j < SIZE_X-1; ++j) {
            covar_phi(i, j) = 0;
            covar_pi(i, j) = 0;
        }
    }

    for (int i = 0; i < SIZE_X-1; ++i) {
        for (int j = 0; j < SIZE_X-1; ++j) {
            A(i, j) = 0;
            B(i, j) = 0;

            if (i == j) {
                A(i, j) = h / TMP;
                B(i, j) = (m * m * h) / TMP + 2. / (h * TMP);
            }
            if (abs(i - j) == 1) {
                B(i, j) = -1. / (h * TMP);
            }
        }
    }

    B(0, SIZE_X - 2) = -1. / (h * TMP);
    B(SIZE_X - 2, 0) = -1. / (h * TMP);

    covar_pi = A.inverse();
    covar_phi = B.inverse();
}

Eigen::MatrixXd normTransform_phi(size, size);

Eigen::LLT<Eigen::MatrixXd> cholSolver_phi(covar_phi);

Eigen::MatrixXd normTransform_pi(size, size);

Eigen::LLT<Eigen::MatrixXd> cholSolver(covar_pi);

Eigen::MatrixXd samples_phi(size, nn);

Eigen::MatrixXd samples_pi(size, nn);




void init_for_generating()
{
    
 

    // We can only use the cholesky decomposition if 
    // the covariance matrix is symmetric, pos-definite.
    // But a covariance matrix might be pos-semi-definite.
    // In that case, we'll go to an EigenSolver
    if (cholSolver_phi.info() == Eigen::Success) {
        // Use cholesky solver
        normTransform_phi = cholSolver_phi.matrixL();
    }
    else {
        // Use eigen solver
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar_phi);
        normTransform_phi = eigenSolver.eigenvectors()
            * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }

    

    // We can only use the cholesky decomposition if 
    // the covariance matrix is symmetric, pos-definite.
    // But a covariance matrix might be pos-semi-definite.
    // In that case, we'll go to an EigenSolver
    if (cholSolver.info() == Eigen::Success) {
        // Use cholesky solver
        normTransform_pi = cholSolver.matrixL();
    }
    else {
        // Use eigen solver
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar_pi);
        normTransform_pi = eigenSolver.eigenvectors()
            * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }

    Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
    Eigen::internal::scalar_normal_dist_op<double>::rng.seed(std::chrono::system_clock::now().time_since_epoch().count()); // Seed the rng

    samples_phi = (normTransform_phi
        * Eigen::MatrixXd::NullaryExpr(size, nn, randN)).colwise()
        + mean_phi;

    samples_pi = (normTransform_pi
        * Eigen::MatrixXd::NullaryExpr(size, nn, randN)).colwise()
        + mean_pi;
}



/*
  Draw nn samples from a size-dimensional normal distribution
  with a specified mean and covariance
*/
void generate_initial_conditions(int iter_num)
{
    phi_averaged_with_eta = 0;

    for (int i = 0; i < SIZE_X-1; ++i) {
		phi_averaged_with_eta += eta(x[i]) * samples_phi(i, iter_num) * h;
	}

    for (int i = 0; i < SIZE_X - 1; ++i) {
        phi_0[i] = samples_phi(i, iter_num);
        pi_0[i] = samples_pi(i, iter_num) - 2. * al * eta(x[i]) * pow(phi_averaged_with_eta, 1); //phi^2
    }

    phi_0[SIZE_X - 1] = phi_0[0];
    pi_0[SIZE_X - 1] = pi_0[0];
}