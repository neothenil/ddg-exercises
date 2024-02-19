// PLEASE READ:
//
// This file implements additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because
// we are "inside" the class, we no longer have to call
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
// Functions in this file can be called from other projects simply by using geometry->cotan(he),
// geometry->barycentricDualArea(v), etc. where "geometry" is a pointer to a VertexPositionGeometry. This avoids having
// to declare a GeometryRoutines object in every project, and also mimics the way that geometry routines are normally
// called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.


#include "geometrycentral/surface/vertex_position_geometry.h"
#include <complex>

namespace geometrycentral {
namespace surface {

/*
 * Compute the Euler characteristic of the mesh.
 */
int VertexPositionGeometry::eulerCharacteristic() const {
    return (int)mesh.nVertices() - (int)mesh.nEdges() + (int)mesh.nFaces();
}

/*
 * Compute the mean length of all the edges in the mesh.
 *
 * Input:
 * Returns: The mean edge length.
 */
double VertexPositionGeometry::meanEdgeLength() const {

    double total = 0.0;
    for (Edge e : mesh.edges()) {
        total += edgeLength(e);
    }
    return total / mesh.nEdges();
}

/*
 * Compute the total surface area of the mesh.
 *
 * Input:
 * Returns: The surface area of the mesh.
 */
double VertexPositionGeometry::totalArea() const {

    double total = 0.0;
    for (Face f : mesh.faces()) {
        total += faceArea(f);
    }
    return total;
}

/*
 * Computes the cotangent of the angle opposite to a halfedge. (Do NOT use built-in function for this)
 *
 * Input: The halfedge whose cotan weight is to be computed.
 * Returns: The cotan of the angle opposite the given halfedge.
 */
double VertexPositionGeometry::cotan(Halfedge he) const {

    // TODO
    if (he.isInterior()) {
        Vector3 a = inputVertexPositions[he.next().tipVertex()];
        Vector3 b = inputVertexPositions[he.tailVertex()];
        Vector3 c = inputVertexPositions[he.tipVertex()];
        Vector3 ab = b - a;
        Vector3 ac = c - a;
        return dot(ab, ac) / norm(cross(ab, ac));
    }
    return 0.0;
    // return 0; // placeholder
}

/*
 * Computes the barycentric dual area of a vertex.
 *
 * Input: The vertex whose barycentric dual area is to be computed.
 * Returns: The barycentric dual area of the given vertex.
 */
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {

    // TODO
    double sum = 0.0;
    for (Face f : v.adjacentFaces()) {
        sum += faceArea(f);
    }
    return sum / 3.0;
    // return 0; // placeholder
}

/*
 * Computes the angle (in radians) at a given corner. (Do NOT use built-in function for this)
 *
 *
 * Input: The corner at which the angle needs to be computed.
 * Returns: The angle clamped between 0 and π.
 */
double VertexPositionGeometry::angle(Corner c) const {

    // TODO
    Halfedge outgoing = c.halfedge();
    Halfedge income = outgoing.next().next();
    Vector3 i = inputVertexPositions[c.vertex()];
    Vector3 j = inputVertexPositions[outgoing.tipVertex()];
    Vector3 k = inputVertexPositions[income.tailVertex()];
    Vector3 ij = j - i;
    Vector3 ik = k - i;
    return acos(dot(ij, ik) / norm(ij) / norm(ik));
    // return 0; // placeholder
}

/*
 * Computes the signed angle (in radians) between two adjacent faces. (Do NOT use built-in function for this)
 *
 * Input: The halfedge (shared by the two adjacent faces) on which the dihedral angle is computed.
 * Returns: The dihedral angle.
 */
double VertexPositionGeometry::dihedralAngle(Halfedge he) const {

    // TODO
    if (he.edge().isBoundary())
        return 0.0;
    Vector3 i = inputVertexPositions[he.tailVertex()];
    Vector3 j = inputVertexPositions[he.tipVertex()];
    Vector3 k = inputVertexPositions[he.next().tipVertex()];
    Vector3 l = inputVertexPositions[he.twin().next().tipVertex()];
    Vector3 ij = j - i, ik = k - i, il = l - i;
    Vector3 Nijk = cross(ij, ik);
    Vector3 Njil = cross(il, ij);
    Nijk = Nijk / norm(Nijk);
    Njil = Njil / norm(Njil);
    return atan2(dot(ij, cross(Nijk, Njil))/norm(ij), dot(Nijk, Njil));
    // return 0; // placeholder
}

/*
 * Computes the normal at a vertex using the "equally weighted" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "equally weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalEquallyWeighted(Vertex v) const {

    // TODO
    Vector3 normal = Vector3::zero();
    for (Face f : v.adjacentFaces()) {
        Halfedge he1 = f.halfedge();
        Halfedge he2 = he1.next();
        Vector3 ab = inputVertexPositions[he1.tipVertex()] -
                     inputVertexPositions[he1.tailVertex()];
        Vector3 bc = inputVertexPositions[he2.tipVertex()] -
                     inputVertexPositions[he2.tailVertex()];
        Vector3 faceNormal = cross(ab, bc);
        normal += faceNormal / norm(faceNormal);
    }
    return normal / norm(normal);
    // return {0, 0, 0}; // placeholder
}

/*
 * Computes the normal at a vertex using the "tip angle weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "tip angle weights" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAngleWeighted(Vertex v) const {

    // TODO
    Vector3 normal = Vector3::zero();
    for (Halfedge he : v.outgoingHalfedges()) {
        Corner c = he.corner();
        double phi = angle(c);
        Halfedge he1 = he.next();
        Vector3 ij = inputVertexPositions[he.tipVertex()] -
                     inputVertexPositions[he.tailVertex()];
        Vector3 jk = inputVertexPositions[he1.tipVertex()] -
                     inputVertexPositions[he1.tailVertex()];
        Vector3 faceNormal = cross(ij, jk);
        normal += phi * faceNormal / norm(faceNormal);
    }
    return normal / norm(normal);
    // return {0, 0, 0}; // placeholder
}

/*
 * Computes the normal at a vertex using the "inscribed sphere" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "inscribed sphere" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalSphereInscribed(Vertex v) const {

    // TODO
    Vector3 normal = Vector3::zero();
    for (Halfedge outgoing : v.outgoingHalfedges()) {
        Halfedge incoming = outgoing.next().next();
        Vector3 ij = inputVertexPositions[outgoing.tipVertex()] -
                     inputVertexPositions[outgoing.tailVertex()];
        Vector3 ik = inputVertexPositions[incoming.tailVertex()] -
                     inputVertexPositions[incoming.tipVertex()];
        Vector3 faceNormal = cross(ij, ik);
        double factor = norm(ij) * norm(ij) * norm(ik) * norm(ik);
        normal += faceNormal / factor;
    }
    return normal / norm(normal);
    // return {0, 0, 0}; // placeholder
}

/*
 * Computes the normal at a vertex using the "face area weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "face area weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAreaWeighted(Vertex v) const {

    // TODO
    Vector3 normal = Vector3::zero();
    for (Halfedge outgoing : v.outgoingHalfedges()) {
        Halfedge incoming = outgoing.next().next();
        Vector3 ij = inputVertexPositions[outgoing.tipVertex()] -
                     inputVertexPositions[outgoing.tailVertex()];
        Vector3 ik = inputVertexPositions[incoming.tailVertex()] -
                     inputVertexPositions[incoming.tipVertex()];
        Vector3 faceNormal = cross(ij, ik);
        normal += faceNormal * 0.5;
    }
    return normal / norm(normal);
    // return {0, 0, 0}; // placeholder
}

/*
 * Computes the normal at a vertex using the "Gauss curvature" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "Gauss curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalGaussianCurvature(Vertex v) const {

    // TODO
    Vector3 normal = Vector3::zero();
    for (Halfedge outgoing : v.outgoingHalfedges()) {
        Vector3 ij = inputVertexPositions[outgoing.tipVertex()] -
                     inputVertexPositions[outgoing.tailVertex()];
        double theta = dihedralAngle(outgoing);
        normal += theta * ij / norm(ij);
    }
    return normal / norm(normal);
    // return {0, 0, 0}; // placeholder
}

/*
 * Computes the normal at a vertex using the "mean curvature" method (equivalent to the "area gradient" method).
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "mean curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalMeanCurvature(Vertex v) const {

    // TODO
    Vector3 normal = Vector3::zero();
    for (Halfedge outgoing : v.outgoingHalfedges()) {
        Vector3 ij = inputVertexPositions[outgoing.tipVertex()] -
                     inputVertexPositions[outgoing.tailVertex()];
        double cotan1 = cotan(outgoing);
        double cotan2 = cotan(outgoing.twin());
        normal += (cotan1 + cotan2) * ij;
    }
    return normal / norm(normal);
    // return {0, 0, 0}; // placeholder
}

/*
 * Computes the angle defect at a vertex.
 *
 * Input: The vertex whose angle defect is to be computed.
 * Returns: The angle defect of the given vertex.
 */
double VertexPositionGeometry::angleDefect(Vertex v) const {

    // TODO
    double total = 0.0;
    for (Halfedge he : v.outgoingHalfedges()) {
        Corner c = he.corner();
        total += angle(c);
    }
    return 2 * PI - total;
    // return 0; // placeholder
}

/*
 * Computes the total angle defect of the mesh.
 *
 * Input:
 * Returns: The total angle defect
 */
double VertexPositionGeometry::totalAngleDefect() const {

    // TODO
    double total = 0.0;
    for (Vertex v : mesh.vertices()) {
        total += angleDefect(v);
    }
    return total; // placeholder
}

/*
 * Computes the (integrated) scalar mean curvature at a vertex.
 *
 * Input: The vertex whose mean curvature is to be computed.
 * Returns: The mean curvature at the given vertex.
 */
double VertexPositionGeometry::scalarMeanCurvature(Vertex v) const {

    // TODO
    double ans = 0.0;
    for (Halfedge outgoing : v.outgoingHalfedges()) {
        Vector3 ij = inputVertexPositions[outgoing.tipVertex()] -
                     inputVertexPositions[outgoing.tailVertex()];
        double theta = dihedralAngle(outgoing);
        ans += theta * norm(ij);
    }
    return 0.5 * ans;
    // return 0; // placeholder
}

/*
 * Computes the circumcentric dual area of a vertex.
 *
 * Input: The vertex whose circumcentric dual area is to be computed.
 * Returns: The circumcentric dual area of the given vertex.
 */
double VertexPositionGeometry::circumcentricDualArea(Vertex v) const {

    // TODO
    double ans = 0.0;
    for (Halfedge outgoing : v.outgoingHalfedges()) {
        Vector3 ij = inputVertexPositions[outgoing.tipVertex()] -
                     inputVertexPositions[outgoing.tailVertex()];
        double cotan1 = cotan(outgoing);
        double cotan2 = cotan(outgoing.twin());
        ans += norm(ij) * norm(ij) * (cotan1 + cotan2);
    }
    return 1.0 / 8.0 * ans;
    // return 0; // placeholder
}

/*
 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
 *
 * Input: The vertex on which the principal curvatures need to be computed.
 * Returns: A std::pair containing the minimum and maximum principal curvature values at a vertex.
 */
std::pair<double, double> VertexPositionGeometry::principalCurvatures(Vertex v) const {

    // TODO
    std::pair<double, double>  ans(0, 0);
    double area = circumcentricDualArea(v);
    double gaussian = angleDefect(v) / area;
    double mean = scalarMeanCurvature(v) / area;
    ans.first = mean - sqrt(mean * mean - gaussian);
    ans.second = mean + sqrt(mean * mean - gaussian);
    return ans;
    // return std::make_pair(0, 0); // placeholder
}


/*
 * Builds the sparse POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace matrix,
 * multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse positive definite Laplace matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrix() const {

    // TODO
    size_t nVertices = mesh.nVertices();
    SparseMatrix<double> matrix(nVertices, nVertices);
    std::vector<Eigen::Triplet<double>> triplets;

    for (Vertex v : mesh.vertices()) {
        size_t row = v.getIndex();
        double diag = 0.0;
        for (Halfedge outgoing : v.outgoingHalfedges()) {
            size_t col = outgoing.tipVertex().getIndex();
            double cotan1 = cotan(outgoing);
            double cotan2 = cotan(outgoing.twin());
            double elementVal = 0.5 * (cotan1 + cotan2);
            diag += elementVal;
            triplets.push_back(Eigen::Triplet<double>(row, col, -elementVal));
        }
        diag += 1e-8;  // shift the diagonal element by a small constant
        triplets.push_back(Eigen::Triplet<double>(row, row, diag));
    }

    matrix.setFromTriplets(triplets.begin(), triplets.end());
    return matrix;
    // return identityMatrix<double>(1); // placeholder
}

/*
 * Builds the sparse diagonal mass matrix containing the barycentric dual area of each vertex.
 *
 * Input:
 * Returns: Sparse mass matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {

    // TODO
    size_t nVertices = mesh.nVertices();
    SparseMatrix<double> matrix(nVertices, nVertices);
    std::vector<Eigen::Triplet<double>> triplets;

    for (Vertex v : mesh.vertices()) {
        size_t idx = v.getIndex();
        double area = barycentricDualArea(v);
        triplets.push_back(Eigen::Triplet<double>(idx, idx, area));
    }

    matrix.setFromTriplets(triplets.begin(), triplets.end());
    return matrix;
    // return identityMatrix<double>(1); // placeholder
}

/*
 * Builds the sparse complex POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace
 * matrix, multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse complex positive definite Laplace matrix for the mesh.
 */
SparseMatrix<std::complex<double>> VertexPositionGeometry::complexLaplaceMatrix() const {

    // TODO
    size_t nVertices = mesh.nVertices();
    SparseMatrix<std::complex<double>> matrix(nVertices, nVertices);
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;

    for (Vertex v : mesh.vertices()) {
        size_t row = v.getIndex();
        std::complex<double> diag = 0.0;
        for (Halfedge outgoing : v.outgoingHalfedges()) {
            size_t col = outgoing.tipVertex().getIndex();
            double cotan1 = cotan(outgoing);
            double cotan2 = cotan(outgoing.twin());
            std::complex<double> elementVal = 0.5 * (cotan1 + cotan2);
            diag += elementVal;
            triplets.push_back(Eigen::Triplet<std::complex<double>>(row, col, -elementVal));
        }
        diag += 1e-8; // shift the diagonal element by a small constant
        triplets.push_back(Eigen::Triplet<std::complex<double>>(row, row, diag));
    }

    matrix.setFromTriplets(triplets.begin(), triplets.end());
    return matrix;
    // return identityMatrix<std::complex<double>>(1); // placeholder
}

/*
 * Compute the center of mass of a mesh.
 */
Vector3 VertexPositionGeometry::centerOfMass() const {

    // Compute center of mass.
    Vector3 center = {0.0, 0.0, 0.0};
    for (Vertex v : mesh.vertices()) {
        center += inputVertexPositions[v];
    }
    center /= mesh.nVertices();

    return center;
}

/*
 * Centers a mesh about the origin.
 * Also rescales the mesh to unit radius if <rescale> == true.
 */
void VertexPositionGeometry::normalize(const Vector3& origin, bool rescale) {

    // Compute center of mass.
    Vector3 center = centerOfMass();

    // Translate to origin [of original mesh].
    double radius = 0;
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] -= center;
        radius = std::max(radius, inputVertexPositions[v].norm());
    }

    // Rescale.
    if (rescale) {
        for (Vertex v : mesh.vertices()) {
            inputVertexPositions[v] /= radius;
        }
    }

    // Translate to origin [of original mesh].
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] += origin;
    }
}

} // namespace surface
} // namespace geometrycentral