#ifndef TOROIDALHONEYCOMBVERTEXMESHGENERATORMUTABLE_HPP_
#define TOROIDALHONEYCOMBVERTEXMESHGENERATORMUTABLE_HPP_

#include <cmath>
#include <vector>

#include "HoneycombVertexMeshGenerator.hpp"
#include "Toroidal2dVertexMeshWithMutableSize.hpp"

/**
 * Honeycomb mesh generator that creates a 2D "toroidal" mesh (one in which
 * periodicity is imposed on the left and right and top and bottom boundaries)
 * for use with vertex-based simulations.
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 */
class ToroidalHoneycombVertexMeshGeneratorMutable : HoneycombVertexMeshGenerator
{
public:

    /**
     * Constructor.
     *
     * @param numElementsAcross  The number of columns of elements in the mesh.  This MUST be an even number.
     * @param numElementsUp  The number of rows of elements in the mesh.   This MUST be an even number.
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangement (defaults to 0.01)
     * @param t2Threshold the maximum threshold distance for Type 2 swaps (defaults to 0.001)
     */
    ToroidalHoneycombVertexMeshGeneratorMutable(unsigned numElementsAcross,
                                         unsigned numElementsUp,
                                         double cellRearrangementThreshold=0.01,
                                         double t2Threshold=0.001,
                                        double elementArea=0);
    /**
     * @return a 2D honeycomb mesh
     */
    MutableVertexMesh<2,2>* GetMesh();

    /**
     * @return a 2D honeycomb mesh with periodic left/right and top/bottom boundaries
     */
    Toroidal2dVertexMeshWithMutableSize* GetMutableToroidalMesh();
};

#endif /*TOROIDALHONEYCOMBVERTEXMESHGENERATORMUTABLE_HPP_*/
