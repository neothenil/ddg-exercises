// Implement member functions for MeanCurvatureFlow class.
#include "mean-curvature-flow.h"
#include "geometrycentral/numerical/linear_solvers.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
MeanCurvatureFlow::MeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
SparseMatrix<double> MeanCurvatureFlow::buildFlowOperator(const SparseMatrix<double>& M, double h) const {

    // TODO
    return M + h * geometry->laplaceMatrix();
    // return identityMatrix<double>(1); // placeholder
}

/*
 * Performs mean curvature flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void MeanCurvatureFlow::integrate(double h) {

    // TODO
    // Note: Geometry Central has linear solvers: https://geometry-central.net/numerical/linear_solvers/
    // Note: Update positions via geometry->inputVertexPositions
    SparseMatrix<double> M = geometry->massMatrix();
    SparseMatrix<double> A = buildFlowOperator(M, h);
    size_t N = mesh->nVertices();
    Vector<double> bx = Vector<double>::Zero(N);
    Vector<double> by = Vector<double>::Zero(N);
    Vector<double> bz = Vector<double>::Zero(N);
    for (Vertex v : mesh->vertices()) {
        size_t idx = v.getIndex();
        Vector3 b = M.coeff(idx, idx) * geometry->inputVertexPositions[v];
        bx(idx) = b.x;
        by(idx) = b.y;
        bz(idx) = b.z;
    }
    PositiveDefiniteSolver<double> solver(A);
    Vector<double> fx = Vector<double>::Zero(N);
    Vector<double> fy = Vector<double>::Zero(N);
    Vector<double> fz = Vector<double>::Zero(N);
    solver.solve(fx, bx);
    solver.solve(fy, by);
    solver.solve(fz, bz);
    for (Vertex v : mesh->vertices()) {
        size_t idx = v.getIndex();
        geometry->inputVertexPositions[v] = Vector3{fx(idx), fy(idx), fz(idx)};
        // geometry->inputVertexPositions[v] = geometry->inputVertexPositions[v]; // placeholder
    }
}