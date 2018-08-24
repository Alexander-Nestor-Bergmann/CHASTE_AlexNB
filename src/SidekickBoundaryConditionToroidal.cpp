#include "SidekickBoundaryConditionToroidal.hpp"
#include "Toroidal2dVertexMeshWithMutableSize.hpp"

SidekickBoundaryConditionToroidal::SidekickBoundaryConditionToroidal(AbstractCellPopulation<2>* pCellPopulation)
    : AbstractCellPopulationBoundaryCondition<2>(pCellPopulation)
{
}

void SidekickBoundaryConditionToroidal::ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
{

    // Pointer to mesh
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


    // std::cout << "/* message */" << '\n';
    // std::cout << "/* message */" << '\n';
    // Set to store boundary
    std::set<unsigned> boundaryNodes;
    // Loop over every node
    for (unsigned node_index=0; node_index<this->mpCellPopulation->GetNumNodes(); node_index++)
    {
        // Get the node
        Node<2>* p_node = this->mpCellPopulation->GetNode(node_index);

        // Get neighbours
        std::set<unsigned> neighbourNodes = this->mpCellPopulation->GetNeighbouringNodeIndices(node_index);

        // // Iterate over connected nodes
        for(auto neighbourIndex : neighbourNodes)
        {
             // Get dist between nodes
             double nodeY = p_node->rGetLocation()[1];
             double nodeX = p_node->rGetLocation()[0];
             Node<2>* neighbour_p_node = this->mpCellPopulation->GetNode(neighbourIndex);
             double neighbourY = neighbour_p_node->rGetLocation()[1];
             double neighbourX = neighbour_p_node->rGetLocation()[0];

             // angle to node
             double angleNode = std::atan2(nodeY-boxCentroidY,
                                         nodeX-boxCentroidX);

             double width = p_static_cast_mesh_toroidal->GetWidth(0);
             double height = p_static_cast_mesh_toroidal->GetWidth(1);
             // If on other side of mesh, mark as boundary node.
             if( (fabs(nodeX - neighbourX) > 0.5*width || fabs(nodeY - neighbourY) > 0.5*height) && (angleNode <= anglex0y0 || angleNode >= anglex0y1) )
             {
                 boundaryNodes.insert(node_index);
                 for(auto nabIndex : neighbourNodes)
                 {
                     Node<2>* temp_neighbour_p_node = this->mpCellPopulation->GetNode(nabIndex);
                     // Get position
                     double temp_neighbourY = temp_neighbour_p_node->rGetLocation()[1];
                     double temp_neighbourX = temp_neighbour_p_node->rGetLocation()[0];
                     // angle to node
                     double angleNeighbour = std::atan2(temp_neighbourY-boxCentroidY,
                                                 temp_neighbourX-boxCentroidX);
                     if( (angleNeighbour <= anglex0y0 || angleNeighbour >= anglex0y1) )
                     {
                         boundaryNodes.insert(nabIndex);
                     }

                 }
                 break;
             }
        }

    }
    // If it was a boundary node, check where it is and pull it
    for(auto n_index : boundaryNodes)
    {
        // Get the node
        Node<2>* p_node = this->mpCellPopulation->GetNode(n_index);
        // Set its location to the old location
        c_vector<double, 2> old_node_location;
        old_node_location = rOldLocations.find(p_node)->second;
        p_node->rGetModifiableLocation()[0] = old_node_location[0];

    }

        // // Iterate over connected nodes
        // for (Node<2>::ContainingElementIterator it = p_node->ContainingElementsBegin();
        //      it != p_node->ContainingElementsEnd();
        //      ++it)
        //  {
        //      // Get dist between nodes
        //      double separation = this->mpCellPopulation->rGetMesh().GetDistanceBetweenNodes(node_index, *it);
        //      double width = p_static_cast_mesh_toroidal->GetWidth(0);
        //      double height = p_static_cast_mesh_toroidal->GetWidth(1);
        //      // If on other side of mesh, mark as boundary node.
        //      if( separation > 0.5*width || separation > 0.5*height)
        //      {
        //          IsBoundaryNode = true;
        //      }
        //  }
        //
        //  // If it was a boundary node, check where it is and pull it
        //  if(IsBoundaryNode)
        //  {
        //      double nodeX = p_node->rGetLocation()[0];
        //      double nodeY = p_node->rGetLocation()[1];
        //      // angle to node
        //      double angleNode = std::atan2(nodeY-boxCentroidY,
        //                                  nodeX-boxCentroidX);
        //
        //      // If we are at the end of the box, fix it.
        //      if ( angleNode <= anglex0y0 || angleNode >= anglex0y1 )
        //      {
        //          // std::cout << p_node->rGetLocation()[0] << '\n';
        //
        //          c_vector<double, 2> old_node_location;
        //          old_node_location = rOldLocations.find(p_node)->second;
        //          // std::cout << old_node_location[0] << '\n';
        //          p_node->rGetModifiableLocation()[0] = old_node_location[0];
        //
        //          // std::cout << p_node->rGetModifiableLocation()[1] << ", ";
        //          // std::cout << "/* message */" << '\n';
        //      }
        //  }
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
