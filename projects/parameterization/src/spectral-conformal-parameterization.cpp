// Implement member functions for SpectralConformalParameterization class.
#include "spectral-conformal-parameterization.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
SpectralConformalParameterization::SpectralConformalParameterization(ManifoldSurfaceMesh* inputMesh,
                                                                     VertexPositionGeometry* inputGeo) {

    this->mesh = inputMesh;
    this->geometry = inputGeo;
}

/*
 * Builds the complex conformal energy matrix EC = ED - A.
 *
 * Input:
 * Returns: A complex sparse matrix representing the conformal energy
 */
SparseMatrix<std::complex<double>> SpectralConformalParameterization::buildConformalEnergy() const {

    // TODO
    // Construct the matrix in A.
    const std::complex<double> constant(0.0, -0.25);
    size_t nVertices = mesh->nVertices();
    SparseMatrix<std::complex<double>> matrixA(nVertices, nVertices);
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;
    for (BoundaryLoop bl : mesh->boundaryLoops()) {
        for (Halfedge he : bl.adjacentHalfedges()) {
            // Half edges in boundary loop have opposite direction.
            size_t col = he.tailVertex().getIndex();
            size_t row = he.tipVertex().getIndex();
            triplets.push_back(Eigen::Triplet<std::complex<double>>(row, col, constant));
            triplets.push_back(Eigen::Triplet<std::complex<double>>(col, row, -constant));
        }
    }
    matrixA.setFromTriplets(triplets.begin(), triplets.end());
    return 0.5 * geometry->complexLaplaceMatrix() - matrixA;
    // return identityMatrix<std::complex<double>>(1); // placeholder
}


/*
 * Flattens the input surface mesh with 1 or more boundaries conformally.
 *
 * Input:
 * Returns: A MeshData container mapping each vertex to a vector of planar coordinates.
 */
VertexData<Vector2> SpectralConformalParameterization::flatten() const {

    // TODO
    auto vertexData = VertexData<Vector2>(*mesh);
    SparseMatrix<std::complex<double>> conformalEnergy = buildConformalEnergy();
    Vector<std::complex<double>> eigenVector = solveInversePowerMethod(conformalEnergy);
    for (Vertex v : mesh->vertices()) {
        size_t idx = v.getIndex();
        vertexData[v] = Vector2{eigenVector(idx).real(), eigenVector(idx).imag()};
    }
    return vertexData;
    // return VertexData<Vector2>(*mesh); // placeholder
}