#include "SidekickBoundaryConditionToroidal.hpp"
#include "Toroidal2dVertexMeshWithMutableSize.hpp"

SidekickBoundaryConditionToroidal::SidekickBoundaryConditionToroidal(AbstractCellPopulation<2>* pCellPopulation)
    : AbstractCellPopulationBoundaryCondition<2>(pCellPopulation)
{
}

void SidekickBoundaryConditionToroidal::ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
{

    // Pointer to mesh --
    AbstractMesh<2, 2>& r_mesh = this->mpCellPopulation->rGetMesh();
    Toroidal2dVertexMeshWithMutableSize* p_static_cast_mesh_toroidal = static_cast<Toroidal2dVertexMeshWithMutableSize*>(&r_mesh);

    // Coords of box
    double currentXLower = p_static_cast_mesh_toroidal->GetBoxCoords(0);
    double currentXUpper = p_static_cast_mesh_toroidal->GetBoxCoords(1);
    double currentYLower = p_static_cast_mesh_toroidal->GetBoxCoords(2);
    double currentYUpper = p_static_cast_mesh_toroidal->GetBoxCoords(3);
    // Centroid of box
    double boxCentroidX = (currentXLower + currentXUpper)/2;
    double boxCentroidY = (currentYLower + currentYUpper)/2;
    // Angles to corners of box. Used to see which side boundary vertices
    // are on
    double anglex0y0 = std::atan2(currentYLower-boxCentroidY,
                                currentXLower-boxCentroidX);
    double anglex0y1 = std::atan2(currentYUpper-boxCentroidY,
                                currentXLower-boxCentroidX);
    // double anglex1y1 = std::atan2(currentYUpper-boxCentroidY,
                                // currentXUpper-boxCentroidX);
    // double anglex1y0 = std::atan2(currentYLower-boxCentroidY,
                                // currentXUpper-boxCentroidX);

    // Set to store boundary nodes
    std::set<unsigned> boundaryNodes = p_static_cast_mesh_toroidal->GetBoundaryNodes();
    // If it was a LHS boundary node, check where it is and fix it
    for(auto n_index : boundaryNodes)
    {
        // Get the node
        Node<2>* p_node = this->mpCellPopulation->GetNode(n_index);
        double nodeX = p_node->rGetLocation()[0];
        double nodeY = p_node->rGetLocation()[1];
        // angle to node
        double angleNode = std::atan2(nodeY-boxCentroidY,
                                    nodeX-boxCentroidX);
        // If on LHS
         if( angleNode <= anglex0y0 || angleNode >= anglex0y1 )
         {
             // Set its location to the old location
             c_vector<double, 2> old_node_location;
             old_node_location = rOldLocations.find(p_node)->second;

             p_node->rGetModifiableLocation()[0] = old_node_location[0];
         }
    }
    // p_static_cast_mesh_toroidal->SetBoxCoords(0, 0);
    // p_static_cast_mesh_toroidal->RefitPeriodicBox();


}

bool SidekickBoundaryConditionToroidal::VerifyBoundaryCondition()
{
    ///\todo Check the boundary condition is actually satisfied
    bool condition_satisfied = true;
    return condition_satisfied;
}

void SidekickBoundaryConditionToroidal::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    AbstractCellPopulationBoundaryCondition<2>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
