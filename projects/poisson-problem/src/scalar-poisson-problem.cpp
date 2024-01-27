// Implement member functions for ScalarPoissonProblem class.
#include "scalar-poisson-problem.h"
#include "geometrycentral/numerical/linear_solvers.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
ScalarPoissonProblem::ScalarPoissonProblem(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;

    // TODO: Build member variables A (Laplace matrix), M (mass matrix), total area
    this->A = inputGeo->laplaceMatrix();
    this->M = inputGeo->massMatrix();
    this->totalArea = inputGeo->totalArea();
    // this->A = identityMatrix<double>(1); // placeholder
    // this->M = identityMatrix<double>(1); // placeholder
    // this->totalArea = 0;                 // placeholder
}

/*
 * Computes the solution of the poisson problem Ax = -M(rho - rhoBar), where A is the POSITIVE DEFINITE Laplace matrix
 * and M is the mass matrix.
 *
 * Input: <rho>, the density of vertices in the mesh.
 * Returns: The solution vector.
 */
Vector<double> ScalarPoissonProblem::solve(const Vector<double>& rho) const {

    // TODO
    // Note: Geometry Central has linear solvers: https://geometry-central.net/numerical/linear_solvers/
    size_t n = rho.rows();
    double rhoBar = 0.0;
    for (size_t i = 0; i < n; ++i)
       rhoBar += rho(i);
    rhoBar /= (double)n;
    // compute -M(rho - rhoBar)
    Vector<double> rhs = Vector<double>::Zero(n);
    for (size_t i = 0; i < n; ++i) {
       rhs(i) = rhoBar - rho(i);
    }
    rhs = M * rhs;
    SparseMatrix<double> lhs = A;
    Vector<double> result = solvePositiveDefinite(lhs, rhs);
    double avg = 0.0;
    for (size_t i = 0; i < n; ++i)
        avg += result(i);
    avg /= (double)n;
    for (size_t i = 0; i < n; ++i)
        result(i) -= avg;
    return result;
    // return Vector<double>::Zero(rho.rows()); // placeholder
}