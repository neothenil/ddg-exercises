// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {

    // TODO
    size_t nVertices = mesh.nVertices();
    SparseMatrix<double> h0(nVertices, nVertices);
    std::vector<Eigen::Triplet<double>> triplets;

    for (Vertex v : mesh.vertices()) {
        size_t idx = v.getIndex();
        double dualArea = barycentricDualArea(v);
        triplets.push_back(Eigen::Triplet<double>(idx, idx, dualArea));
    }

    h0.setFromTriplets(triplets.begin(), triplets.end());
    return h0;
    // return identityMatrix<double>(1); // placeholder
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    // TODO
    size_t nEdges = mesh.nEdges();
    SparseMatrix<double> h1(nEdges, nEdges);
    std::vector<Eigen::Triplet<double>> triplets;

    for (Edge e : mesh.edges()) {
        size_t idx = e.getIndex();
        Halfedge he1, he2;
        he1 = e.halfedge();
        he2 = he1.twin();
        double ratio = (cotan(he1) + cotan(he2)) / 2.0;
        triplets.push_back(Eigen::Triplet<double>(idx, idx, ratio));
    }

    h1.setFromTriplets(triplets.begin(), triplets.end());
    return h1;
    // return identityMatrix<double>(1); // placeholder
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

    // TODO
    size_t nFaces = mesh.nFaces();
    SparseMatrix<double> h2(nFaces, nFaces);
    std::vector<Eigen::Triplet<double>> triplets;

    for (Face f : mesh.faces()) {
        size_t idx = f.getIndex();
        double area = faceArea(f);
        triplets.push_back(Eigen::Triplet<double>(idx, idx, 1.0 / area));
    }

    h2.setFromTriplets(triplets.begin(), triplets.end());
    return h2;
    // return identityMatrix<double>(1); // placeholder
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    // TODO
    size_t nEdges = mesh.nEdges();
    size_t nVertices = mesh.nVertices();
    SparseMatrix<double> d0(nEdges, nVertices);
    std::vector<Eigen::Triplet<double>> triplets;

    for (Edge e : mesh.edges()) {
        size_t idxEdge = e.getIndex();
        Vertex first = e.firstVertex();
        Vertex second = e.secondVertex();
        size_t idxFirstVert = first.getIndex();
        size_t idxSecondVert = second.getIndex();
        triplets.push_back(Eigen::Triplet<double>(idxEdge, idxFirstVert, -1.0));
        triplets.push_back(Eigen::Triplet<double>(idxEdge, idxSecondVert, 1.0));
    }

    d0.setFromTriplets(triplets.begin(), triplets.end());
    return d0;
    // return identityMatrix<double>(1); // placeholder
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    // TODO
    size_t nEdges = mesh.nEdges();
    size_t nFaces = mesh.nFaces();
    SparseMatrix<double> d1(nFaces, nEdges);
    std::vector<Eigen::Triplet<double>> triplets;

    for (Face f : mesh.faces()) {
        size_t idxFace = f.getIndex();
        for (Halfedge he : f.adjacentHalfedges()) {
            size_t idxEdge = he.edge().getIndex();
            // If the half edge has same orientation with edge, set v to 1, else -1.
            double v = he.orientation() ? 1.0 : -1.0;
            triplets.push_back(Eigen::Triplet<double>(idxFace, idxEdge, v));
        }
    }

    d1.setFromTriplets(triplets.begin(), triplets.end());
    return d1;
    // return identityMatrix<double>(1); // placeholder
}

} // namespace surface
} // namespace geometrycentral