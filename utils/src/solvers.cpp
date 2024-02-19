#include "solvers.h"
#include "geometrycentral/numerical/linear_solvers.h"

/*
 * Compute the inverse of a sparse diagonal matrix.
 *
 * Input: A sparse diagonal matrix <M>.
 * Returns: The inverse of M, which is also a sparse diagonal matrix.
 */
SparseMatrix<double> sparseInverseDiagonal(SparseMatrix<double>& M) {

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    SparseMatrix<double> inv(M.rows(), M.cols());
    for (int i = 0; i < M.rows(); i++) {
        tripletList.push_back(T(i, i, 1.0 / M.coeffRef(i, i)));
    }
    inv.setFromTriplets(tripletList.begin(), tripletList.end());
    return inv;
}

/*
 * Computes the residual of Ax - λx, where x has unit norm and λ = x.Ax.
 *
 * Input: <A>, the complex sparse matrix whose eigendecomposition is being computed; and <x>, the current guess for the
 * smallest eigenvector
 * Returns: The residual
 */
double residual(const SparseMatrix<std::complex<double>>& A, const Vector<std::complex<double>>& x) {

    // TODO
    Vector<std::complex<double>> Ax = A * x;
    std::complex<double> lambda = x.dot(Ax);
    Vector<std::complex<double>> residual = Ax - lambda * x;
    return residual.norm();
    // return 0; // placeholder
}

/*
 * Solves Ax = λx, where λ is the smallest nonzero eigenvalue of A, and x is the corresponding eigenvector.
 *
 * Input: <A>, the complex positive definite sparse matrix whose eigendecomposition is being computed.
 * Returns: The smallest eigenvector of A.
 */
Vector<std::complex<double>> solveInversePowerMethod(const SparseMatrix<std::complex<double>>& A) {

    // TODO
    SparseMatrix<std::complex<double>> localA = A;
    size_t N = A.rows();
    Vector<std::complex<double>> y = Vector<std::complex<double>>::Zero(N);
    y(0) = 1.0;
    PositiveDefiniteSolver<std::complex<double>> solver(localA);
    while (residual(A, y) > 1.0e-10) {
        y = solver.solve(y);
        // center around the origin
        std::complex<double> avg = 0.0;
        for (size_t i = 0; i < N; ++i) {
            avg += y(i);
        }
        avg /= (double)N;
        for (size_t i = 0; i < N; ++i) {
            y(i) -= avg;
        }
        // normalize y
        y.normalize();
    }
    return y;
    // return Vector<std::complex<double>>::Zero(1);
}