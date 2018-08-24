#include "ExtrinsicPullModifierToroidal.hpp"
#include "Toroidal2dVertexMeshWithMutableSize.hpp"

ExtrinsicPullModifierToroidal::ExtrinsicPullModifierToroidal()
    : AbstractCellBasedSimulationModifier<2>(),
      mApplyExtrinsicPullToAllNodes(true),
      mSpeed(1.0)
{
}

ExtrinsicPullModifierToroidal::~ExtrinsicPullModifierToroidal()
{
}

void ExtrinsicPullModifierToroidal::UpdateAtEndOfTimeStep(AbstractCellPopulation<2,2>& rCellPopulation)
{
    double epsilon = 0.8;

    double dt = SimulationTime::Instance()->GetTimeStep();
    unsigned num_nodes = rCellPopulation.GetNumNodes();
    ChasteCuboid<2> bounds = rCellPopulation.rGetMesh().CalculateBoundingBox();
    double x_min = bounds.rGetLowerCorner()[0];
    double x_max = bounds.rGetUpperCorner()[0];
    // double y_min = bounds.rGetLowerCorner()[1];
    // double y_max = bounds.rGetUpperCorner()[1];

    if (mApplyExtrinsicPullToAllNodes)
    {
        // Pull on all nodes, with a constant strain rate
        double width = x_max - x_min;
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            Node<2>* p_node = rCellPopulation.GetNode(node_index);
            double speed = mSpeed * (p_node->rGetLocation()[0] - x_min) / width;

            // Respect the SidekickBoundaryCondition...
            if (p_node->rGetLocation()[0] > x_min + epsilon)
            {
//                if (p_node->rGetLocation()[0] < x_max - epsilon)
//                {
                    p_node->rGetModifiableLocation()[0] += speed*dt;
//                }
            }
        }
    }
    else
    {
        // double pullSpeed;
        // pullSpeed = mSpeed/40.;
        double YshrinkSpeed = mSpeed/2;

        // Pointer to mesh
        AbstractMesh<2, 2>& r_mesh = rCellPopulation.rGetMesh();
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
        double anglex1y1 = std::atan2(currentYUpper-boxCentroidY,
                                    currentXUpper-boxCentroidX);
        double anglex1y0 = std::atan2(currentYLower-boxCentroidY,
                                    currentXUpper-boxCentroidX);

        // Find the boundary nodes
        unsigned num_nodes = rCellPopulation.GetNumNodes();
        // Set to store boundary
        std::set<unsigned> boundaryNodes;
        // Iterate over nodes
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            // Get the node
            Node<2>* p_node = rCellPopulation.GetNode(node_index);

            // Get neighbours
            std::set<unsigned> neighbourNodes = rCellPopulation.GetNeighbouringNodeIndices(node_index);

            // // Iterate over connected nodes
            for(auto neighbourIndex : neighbourNodes)
            {
                 // Get dist between nodes
                 double nodeY = p_node->rGetLocation()[1];
                 double nodeX = p_node->rGetLocation()[0];
                 Node<2>* neighbour_p_node = rCellPopulation.GetNode(neighbourIndex);
                 double neighbourY = neighbour_p_node->rGetLocation()[1];
                 double neighbourX = neighbour_p_node->rGetLocation()[0];

                 double width = p_static_cast_mesh_toroidal->GetWidth(0);
                 double height = p_static_cast_mesh_toroidal->GetWidth(1);
                 // If on other side of mesh, mark as boundary node.
                 if( fabs(nodeX - neighbourX) > 0.5*width || fabs(nodeY - neighbourY) > 0.5*height)
                 {
                     boundaryNodes.insert(node_index);
                     for(auto nabIndex : neighbourNodes)
                     {
                         boundaryNodes.insert(nabIndex);
                     }
                     break;
                 }
            }
        }

        // If it was a boundary node, check where it is and pull it
        for(auto n_index : boundaryNodes)
        {
            // Get the node
            Node<2>* p_node = rCellPopulation.GetNode(n_index);
            double nodeX = p_node->rGetLocation()[0];
            double nodeY = p_node->rGetLocation()[1];

            // angle to node
            double angleNode = std::atan2(nodeY-boxCentroidY,
                                        nodeX-boxCentroidX);
            // Stretch/compress at the correct boundaries
            if (anglex1y0 <= angleNode && anglex1y1 >= angleNode)
            {
                p_node->rGetModifiableLocation()[0] += mSpeed*dt;
            }
            if (anglex1y1 <= angleNode && anglex0y1 >= angleNode)
            {
                p_node->rGetModifiableLocation()[1] -= YshrinkSpeed*dt;
            }
            if (anglex0y0 <= angleNode && anglex1y0 >= angleNode)
            {
                p_node->rGetModifiableLocation()[1] += YshrinkSpeed*dt;
            }
        }

        // Shrink box
        p_static_cast_mesh_toroidal->SetBoxCoords(1, currentXUpper+mSpeed*dt);
        p_static_cast_mesh_toroidal->SetBoxCoords(2,                                                                     currentYLower+YshrinkSpeed*dt);
        p_static_cast_mesh_toroidal->SetBoxCoords(3,                                                                     currentYUpper-YshrinkSpeed*dt);
    }
}

void ExtrinsicPullModifierToroidal::SetupSolve(AbstractCellPopulation<2,2>& rCellPopulation, std::string outputDirectory)
{
}

void ExtrinsicPullModifierToroidal::ApplyExtrinsicPullToAllNodes(bool applyExtrinsicPullToAllNodes)
{
    mApplyExtrinsicPullToAllNodes = applyExtrinsicPullToAllNodes;
}

void ExtrinsicPullModifierToroidal::SetSpeed(double speed)
{
    mSpeed = speed;
}

void ExtrinsicPullModifierToroidal::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ApplyExtrinsicPullToAllNodes>" << mApplyExtrinsicPullToAllNodes << "</ApplyExtrinsicPullToAllNodes>\n";
    *rParamsFile << "\t\t\t<Speed>" << mSpeed << "</Speed>\n";

    // Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<2>::OutputSimulationModifierParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ExtrinsicPullModifierToroidal)
