#ifndef TOROIDAL2DVERTEXMESHWITHMUTABLESIZE_HPP_
#define TOROIDAL2DVERTEXMESHWITHMUTABLESIZE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellPopulation.hpp"
#include "MutableVertexMesh.hpp"
#include "FarhadifarForce.hpp"

/**
 * A subclass of MutableVertexMesh<2,2> for a rectangular mesh with
 * periodic left and right boundaries and top and bottom boundaries,
 * representing a toroidal geometry.
 *
 * The class works by overriding calls such as ReMesh() and
 * GetVectorFromAtoB() so that simulation classes can treat this
 * class in exactly the same way as a MutableMesh<2,2>.
 */
class Toroidal2dVertexMeshWithMutableSize : public MutableVertexMesh<2,2>
{
    friend class TestToroidal2dVertexMeshWithMutableSize;

private:

    /** The width of the mesh, taking account of left-right periodicity. */
    double mWidth;

    /** The height of the mesh, taking account of top-bottom periodicity. */
    double mHeight;

    /** The box coords of the mesh. */
    double mBoxLowerY;
    double mBoxUpperY;
    double mBoxLowerX;
    double mBoxUpperX;

    /** Reference stress tensor. */
    c_matrix<double, 2,2> mReferenceStress;

    /**
     * Auxiliary mesh pointer, created/updated when GetMeshForVtk() is called
     * and stored so that it may be deleted by the destructor.
     */
    VertexMesh<2,2>* mpMeshForVtk;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archives the member variables of the object which
     * have to be preserved during its lifetime.
     *
     * The remaining member variables are re-initialised before being used
     * by each ReMesh() call so they do not need to be archived.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<MutableVertexMesh<2,2> >(*this);
        archive & mWidth;
        archive & mHeight;
        archive & mBoxLowerY;
        archive & mBoxUpperY;
        archive & mBoxLowerX;
        archive & mBoxUpperX;
        mpMeshForVtk = nullptr;
    }

public:

    /**
     * Default constructor.
     *
     * @param width the width of the mesh
     * @param height the height of the mesh
     * @param nodes vector of pointers to nodes
     * @param vertexElements vector of pointers to VertexElements
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangement (defaults to 0.01)
     * @param t2Threshold the maximum threshold distance for Type 2 swaps (defaults to 0.001)
     */
    Toroidal2dVertexMeshWithMutableSize(double width,
                         double height,
                         std::vector<Node<2>*> nodes,
                         std::vector<VertexElement<2,2>*> vertexElements,
                         double cellRearrangementThreshold=0.01,
                         double t2Threshold=0.001);

    /**
     * Constructor.
     */
    Toroidal2dVertexMeshWithMutableSize();

    /**
     * Destructor.
     */
    ~Toroidal2dVertexMeshWithMutableSize();

    /**
     * Overridden GetVectorFromAtoB() method.
     *
     * This method evaluates the (surface) distance between
     * two points in a 2D toroidal geometry.
     *
     * @param rLocation1 the x and y co-ordinates of point 1
     * @param rLocation2 the x and y co-ordinates of point 2
     *
     * @return the vector from location1 to location2
     */
    c_vector<double, 2> GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2);

    /**
     * Overridden SetNode() method.
     *
     * If the location should be set outside a toroidal boundary
     * move it back onto the cylinder.
     *
     * @param nodeIndex is the index of the node to be moved
     * @param point is the new target location of the node
     */
    void SetNode(unsigned nodeIndex, ChastePoint<2> point);

    /**
     * Overridden GetWidth() method.
     *
     * Calculate the 'width' of any dimension of the mesh, taking periodicity
     * into account.
     *
     * @param rDimension a dimension (0 or 1)
     * @return The maximum distance between any nodes in this dimension.
     */
    double GetWidth(const unsigned& rDimension) const;

    /**
     * New SetWidth() method.
     *
     * Calculate the 'width' of any dimension of the mesh, taking periodicity
     * into account.
     *
     * @param rDimension a dimension (0 or 1)
     * @return The maximum distance between any nodes in this dimension.
     */
    void SetWidth(const unsigned& rDimension, double width);

    /**
     * New GetNearestNodeIndexIgnorePeriodicity() method.
     *
     * Calculates the index of the node closest to the given vector, but doesn't account for periodicity.
     * @param rTestPoint reference to the point
     * @return node index
     */
     unsigned GetNearestNodeIndexIgnorePeriodicity(const ChastePoint<2>& rTestPoint);

    /**
     * New GetBoxCoords method.
     *
     * Get position the corners of the periodic box.
     *
     * @param rDimension. References which coord is being changed according to
     * 0, 1, 2, 3 = mBoxLowerX, mBoxUpperX, mBoxLowerY, mBoxUpperY
     * @return value of desired component
     */
    double GetBoxCoords(const unsigned& rDimension);

    /**
     * New SetBoxCoords method.
     *
     * Change position the corners of the periodic box.
     *
     * @param rDimension. References which coord is being changed according to
     * 0, 1, 2, 3 = mBoxLowerX, mBoxUpperX, mBoxLowerY, mBoxUpperY
     * @param newVal new position that corners will be moved to
     */
    void SetBoxCoords(const unsigned& rDimension, double newVal);

    /**
     * New RefitPeriodicBox method.
     *
     * Change re-find the corners of the periodic box, given the width and
     * height.
     *
     * @param rDimension. References which coord is being changed according to
     * 0, 1, 2, 3 = mBoxLowerX, mBoxUpperX, mBoxLowerY, mBoxUpperY
     * @param newVal new position that corners will be moved to
     */
    void RefitPeriodicBox();

    /**
     * New GetBoundaryNodes method.
     *
     * Return boundary nodes in the tissue.
     *
     */
     std::set<unsigned> GetBoundaryNodes();

     /**
      * New SetReferenceStress() method.
      *
      * Sets the stress tensor of the inital tissue. This stress is then used as
      * a base state, about which relaxation is done.
      *
      * @param rCellPopulation The reference cell population to get stress from
      * @param useZero. Boolean check to use zero stress as base.
      */
     void SetReferenceStress(AbstractCellPopulation<2,2>& rCellPopulation, FarhadifarForce<2>* p_force, bool useZero = false);

     /**
      * New RelaxPeriodicBox() method.
      *
      * Calculates the x&y componets of tissue stress and relaxes the box in
      * response via contractions/expansions.
      *
      * @param rCellPopulation The reference cell population to get stress from
      */
     void RelaxPeriodicBox(AbstractCellPopulation<2,2>& rCellPopulation, FarhadifarForce<2>* p_force, double stiffness = 1000);

    /**
     * Overridden AddNode() method.
     *
     * @param pNewNode the node to be added to the mesh
     *
     * @return the global index of the new node
     */
    unsigned AddNode(Node<2>* pNewNode);

    /**
     * Overridden GetMeshForVtk() method.
     *
     * Return a pointer to an extended mesh that is a 'non-periodic'
     * version of our mesh. This can then be used when writing to
     * VTK.
     *
     * @return a non-periodic vertex mesh
     */
     VertexMesh<2,2>* GetMeshForVtk();

     /**
      * Construct the mesh using a MeshReader.
      *
      * @param rMeshReader the mesh reader
      * @param width the mesh width
      * @param height the mesh height
      */
     void ConstructFromMeshReader(AbstractMeshReader<2,2>& rMeshReader, double width, double height);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Toroidal2dVertexMeshWithMutableSize)

#endif /*TOROIDAL2DVERTEXMESHWITHMUTABLESIZE_HPP_*/
