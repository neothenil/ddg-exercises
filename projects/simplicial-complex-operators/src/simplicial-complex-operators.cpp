// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

static size_t countOfAInB(const std::set<size_t>& A, const std::set<size_t>& B);

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    
    size_t nEdges = mesh->nEdges();
    size_t nVertices = mesh->nVertices();
    SparseMatrix<size_t> adjMatrix(nEdges, nVertices);
    std::vector<Eigen::Triplet<size_t>> triplets;
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();

    for (Edge e : mesh->edges()) {
        size_t idxEdge = geometry->edgeIndices[e];
        Vertex first = e.firstVertex();
        Vertex second = e.secondVertex();
        size_t idxFirstVert = geometry->vertexIndices[first];
        size_t idxSecondVert = geometry->vertexIndices[second];
        triplets.push_back(Eigen::Triplet<size_t>(idxEdge, idxFirstVert, 1));
        triplets.push_back(Eigen::Triplet<size_t>(idxEdge, idxSecondVert, 1));
    }

    adjMatrix.setFromTriplets(triplets.begin(), triplets.end());
    return adjMatrix;

    // return identityMatrix<size_t>(1); // placeholder
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    // TODO

    size_t nEdges = mesh->nEdges();
    size_t nFaces = mesh->nFaces();
    SparseMatrix<size_t> adjMatrix(nFaces, nEdges);
    std::vector<Eigen::Triplet<size_t>> triplets;
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    for (Face f : mesh->faces()) {
        size_t idxFace = geometry->faceIndices[f];
        for (Edge e : f.adjacentEdges()) {
            size_t idxEdge = geometry->edgeIndices[e];
            triplets.push_back(Eigen::Triplet<size_t>(idxFace, idxEdge, 1));
        }
    }

    adjMatrix.setFromTriplets(triplets.begin(), triplets.end());
    return adjMatrix;

    //return identityMatrix<size_t>(1); // placeholder
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    // TODO
    size_t num = mesh->nVertices();
    Vector<size_t> v = Vector<size_t>::Zero(num);

    for (const auto& idx : subset.vertices) {
        v(idx) = 1;
    }

    return v;

    //return Vector<size_t>::Zero(1);
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    // TODO
    size_t num = mesh->nEdges();
    Vector<size_t> v = Vector<size_t>::Zero(num);

    for (const auto& idx : subset.edges) {
        v(idx) = 1;
    }

    return v;

    //return Vector<size_t>::Zero(1);
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO
    size_t num = mesh->nFaces();
    Vector<size_t> v = Vector<size_t>::Zero(num);

    for (const auto& idx : subset.faces) {
        v(idx) = 1;
    }

    return v;

    //return Vector<size_t>::Zero(1);
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // TODO
    MeshSubset starSubset = subset.deepCopy();
    for (const auto& idx : starSubset.vertices) {
        std::set<size_t> adjEdges;
        for (SparseMatrix<size_t>::InnerIterator it(A0, idx); it; ++it) {
            if (it.value()) {
                adjEdges.insert(it.row());
            }
        }
        starSubset.addEdges(adjEdges);
    }
    for (const auto& idx : starSubset.edges) {
        std::set<size_t> adjFaces;
        for (SparseMatrix<size_t>::InnerIterator it(A1, idx); it; ++it) {
            if (it.value()) {
                adjFaces.insert(it.row());
            }
        }
        starSubset.addFaces(adjFaces);
    }
    return starSubset;

    //return subset; // placeholder
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
    MeshSubset closureSubset = subset.deepCopy();
    Eigen::SparseMatrix<size_t, Eigen::RowMajor> rowMajorA0 = A0;
    Eigen::SparseMatrix<size_t, Eigen::RowMajor> rowMajorA1 = A1;

    for (const auto& idxFace : closureSubset.faces) {
        std::set<size_t> edges;
        for (Eigen::SparseMatrix<size_t, Eigen::RowMajor>::InnerIterator it(rowMajorA1, idxFace); it; ++it) {
            if (it.value()) {
                edges.insert(it.col());
            }
        }
        closureSubset.addEdges(edges);
    }
    for (const auto& idxEdge : closureSubset.edges) {
        std::set<size_t> vertices;
        for (Eigen::SparseMatrix<size_t, Eigen::RowMajor>::InnerIterator it(rowMajorA0, idxEdge); it; ++it) {
            if (it.value()) {
                vertices.insert(it.col());
            }
        }
        closureSubset.addVertices(vertices);
    }

    return closureSubset;

    //return subset; // placeholder
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    MeshSubset closureOfStar = closure(star(subset));
    MeshSubset starOfClosure = star(closure(subset));
    closureOfStar.deleteSubset(starOfClosure);
    return closureOfStar;

    //return subset; // placeholder
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    Eigen::SparseMatrix<size_t, Eigen::RowMajor> rowMajorA0 = A0;
    Eigen::SparseMatrix<size_t, Eigen::RowMajor> rowMajorA1 = A1;

    // All edges of each face are in this subset.
    for (const auto& idxFace : subset.faces) {
        std::set<size_t> edges;
        for (Eigen::SparseMatrix<size_t, Eigen::RowMajor>::InnerIterator it(rowMajorA1, idxFace); it; ++it) {
            if (it.value()) {
                edges.insert(it.col());
            }
        }
        for (const auto& idxEdge : edges) {
            if (subset.edges.find(idxEdge) == subset.edges.cend()) return false;
        }
    }
    // All vertices of each edge are in this subset.
    for (const auto& idxEdge : subset.edges) {
        std::set<size_t> vertices;
        for (Eigen::SparseMatrix<size_t, Eigen::RowMajor>::InnerIterator it(rowMajorA0, idxEdge); it; ++it) {
            if (it.value()) {
                vertices.insert(it.col());
            }
        }
        for (const auto& idxVertex : vertices) {
            if (subset.vertices.find(idxVertex) == subset.vertices.cend()) return false;
        }
    }
    return true; // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    size_t nFaces, nEdges, nVertices;
    nFaces = subset.faces.size();
    nEdges = subset.edges.size();
    nVertices = subset.vertices.size();

    if (nFaces == 0 && nEdges == 0) return 0;
    if (nFaces == 0) {
        MeshSubset comparedSubset({}, subset.edges, {});
        comparedSubset = closure(comparedSubset);
        if (subset.equals(comparedSubset)) {
            return 1;
        } else {
            return -1;
        }
    }
    MeshSubset comparedSubset({}, {}, subset.faces);
    comparedSubset = closure(comparedSubset);
    if (subset.equals(comparedSubset)) {
        return 2;
    } else {
        return -1;
    }

    //return -1; // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    MeshSubset result;
    for (const auto& idx : subset.edges) {
        std::set<size_t> adjFaces;
        for (SparseMatrix<size_t>::InnerIterator it(A1, idx); it; ++it) {
            if (it.value()) {
                adjFaces.insert(it.row());
            }
        }
        if (countOfAInB(adjFaces, subset.faces) == 1) result.addEdge(idx);
    }
    for (const auto& idx : subset.vertices) {
        std::set<size_t> adjEdges;
        for (SparseMatrix<size_t>::InnerIterator it(A0, idx); it; ++it) {
            if (it.value()) {
                adjEdges.insert(it.row());
            }
        }
        if (countOfAInB(adjEdges, subset.edges) == 1) result.addVertex(idx);
    }
    return closure(result);

    //return subset; // placeholder
}

static size_t countOfAInB(const std::set<size_t>& A, const std::set<size_t>& B) {
    size_t result = 0;
    for (const auto& e : A) {
        result += B.count(e);
    }
    return result;
}