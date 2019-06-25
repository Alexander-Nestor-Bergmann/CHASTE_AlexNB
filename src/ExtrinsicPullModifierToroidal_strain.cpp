#include "ExtrinsicPullModifierToroidal.hpp"
#include "Toroidal2dVertexMeshWithMutableSize.hpp"

ExtrinsicPullModifierToroidal::ExtrinsicPullModifierToroidal()
    : AbstractCellBasedSimulationModifier<2>(),
      mApplyExtrinsicPullToAllNodes(true),
      mSpeed(1.0),
      mApplyExtrinsicPull(true)
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

    if (mApplyExtrinsicPullToAllNodes && mApplyExtrinsicPull)
    {
        // Pull on all nodes, with a constant strain rate
        double width = x_max - x_min;
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            Node<2>* p_node = rCellPopulation.GetNode(node_index);
            double speed = mSpeed * (p_node->rGetLocation()[0] - x_min) / width;

            p_node->rGetModifiableLocation()[0] *= (1+mSpeed*dt);
            p_node->rGetModifiableLocation()[1] /= (1+mSpeed*dt);

//             // Respect the SidekickBoundaryCondition...
//             if (p_node->rGetLocation()[0] > x_min + epsilon)
//             {
// //                if (p_node->rGetLocation()[0] < x_max - epsilon)
// //                {
//                     p_node->rGetModifiableLocation()[0] += speed*dt;
// //                }
//             }
        }
        // Pointer to mesh
        AbstractMesh<2, 2>& r_mesh = rCellPopulation.rGetMesh();
        Toroidal2dVertexMeshWithMutableSize* p_static_cast_mesh_toroidal = static_cast<Toroidal2dVertexMeshWithMutableSize*>(&r_mesh);

        double currentXLower = p_static_cast_mesh_toroidal->GetBoxCoords(0);
        double currentXUpper = p_static_cast_mesh_toroidal->GetBoxCoords(1);
        double currentYLower = p_static_cast_mesh_toroidal->GetBoxCoords(2);
        double currentYUpper = p_static_cast_mesh_toroidal->GetBoxCoords(3);

        p_static_cast_mesh_toroidal->SetBoxCoords(0, currentXLower*(1+mSpeed*dt));
        p_static_cast_mesh_toroidal->SetBoxCoords(1, currentXUpper*(1+mSpeed*dt));
        p_static_cast_mesh_toroidal->SetBoxCoords(2,                                                                     currentYLower/(1+mSpeed*dt));
        p_static_cast_mesh_toroidal->SetBoxCoords(3,                                                                     currentYUpper/(1+mSpeed*dt));

        currentXLower = p_static_cast_mesh_toroidal->GetBoxCoords(0);
        currentXUpper = p_static_cast_mesh_toroidal->GetBoxCoords(1);
        currentYLower = p_static_cast_mesh_toroidal->GetBoxCoords(2);
        currentYUpper = p_static_cast_mesh_toroidal->GetBoxCoords(3);
        p_static_cast_mesh_toroidal->SetWidth(0, currentXUpper - currentXLower);
        p_static_cast_mesh_toroidal->SetWidth(1, currentYUpper - currentYLower);
    }
    else if (mApplyExtrinsicPull)
    {
        // double YshrinkSpeed = mSpeed/2;

        // Pointer to mesh
        AbstractMesh<2, 2>& r_mesh = rCellPopulation.rGetMesh();
        Toroidal2dVertexMeshWithMutableSize* p_static_cast_mesh_toroidal = static_cast<Toroidal2dVertexMeshWithMutableSize*>(&r_mesh);

        // STRETCH LHS BECAUSE THAT'S HOW STRESS RELAXATION RELAXES.
        // Coords of box
        // double currentXLower = p_static_cast_mesh_toroidal->GetBoxCoords(0);
        double currentXUpper = p_static_cast_mesh_toroidal->GetBoxCoords(1);
        // double currentYLower = p_static_cast_mesh_toroidal->GetBoxCoords(2);
        // double currentYUpper = p_static_cast_mesh_toroidal->GetBoxCoords(3);

        // // Centroid of box
        // double boxCentroidX = (currentXLower + currentXUpper)/2;
        // double boxCentroidY = (currentYLower + currentYUpper)/2;
        // // Angles to corners of box. Used to see which side boundary vertices
        // // are on
        // double anglex0y0 = std::atan2(currentYLower-boxCentroidY,
        //                             currentXLower-boxCentroidX);
        // double anglex0y1 = std::atan2(currentYUpper-boxCentroidY,
        //                             currentXLower-boxCentroidX);
        //
        // // Find the boundary nodes
        // std::set<unsigned> boundaryNodes = p_static_cast_mesh_toroidal->GetBoundaryNodes();
        //
        // // If it was a boundary node, check where it is and pull it appropriately
        // for(auto n_index : boundaryNodes)
        // {
        //     // Get the node
        //     Node<2>* p_node = p_static_cast_mesh_toroidal->GetNode(n_index);
        //     double nodeX = p_node->rGetLocation()[0];
        //     double nodeY = p_node->rGetLocation()[1];
        //
        //     // angle to node
        //     double angleNode = std::atan2(nodeY-boxCentroidY,
        //                                 nodeX-boxCentroidX);
        //     // Stretch/compress at the correct boundaries. Note, don't need to stretch upper and RHS because they are deermined by height and width of box, so just stretch the box after.
        //     // LHS
        //     if ( angleNode <= anglex0y0 || angleNode >= anglex0y1 )
        //     {
        //         p_node->rGetModifiableLocation()[0] -= 0.5*mSpeed*dt;
        //     }
        // }

        // Reset the size of the box
        // p_static_cast_mesh_toroidal->SetBoxCoords(0, currentXLower - 0.5*mSpeed*dt);
        p_static_cast_mesh_toroidal->SetBoxCoords(1, currentXUpper + mSpeed*dt);

        // // // Coords of box
        // double currentXLower = p_static_cast_mesh_toroidal->GetBoxCoords(0);
        // double currentXUpper = p_static_cast_mesh_toroidal->GetBoxCoords(1);
        // // double currentYLower = p_static_cast_mesh_toroidal->GetBoxCoords(2);
        // // double currentYUpper = p_static_cast_mesh_toroidal->GetBoxCoords(3);
        // //
// //
//         // // Don't need to move upper and RHS. just change box size.
//         // // Centroid of box
//         // double boxCentroidX = (currentXLower + currentXUpper)/2;
//         // double boxCentroidY = (currentYLower + currentYUpper)/2;
//         // // Angles to corners of box. Used to see which side boundary vertices
//         // // are on
//         // double anglex0y0 = std::atan2(currentYLower-boxCentroidY,
//         //                             currentXLower-boxCentroidX);
//         // double anglex0y1 = std::atan2(currentYUpper-boxCentroidY,
//         //                             currentXLower-boxCentroidX);
//         // // double anglex1y1 = std::atan2(currentYUpper-boxCentroidY,
//         //                             // currentXUpper-boxCentroidX);
//         // double anglex1y0 = std::atan2(currentYLower-boxCentroidY,
//         //                             currentXUpper-boxCentroidX);
//         //
//         // // \TODO Don't need to pull every node, just stretch the RHS of the box, because the RHS vertices are defined by the width.
//         // // Find the boundary nodes
//         // unsigned num_nodes = rCellPopulation.GetNumNodes();
//         // // Set to store boundary
//         // std::set<unsigned> boundaryNodes = p_static_cast_mesh_toroidal->GetBoundaryNodes();
//         //
//         // // If it was a boundary node, check where it is and pull it appropriately
//         // for(auto n_index : boundaryNodes)
//         // {
//         //     // Get the node
//         //     Node<2>* p_node = rCellPopulation.GetNode(n_index);
//         //     double nodeX = p_node->rGetLocation()[0];
//         //     double nodeY = p_node->rGetLocation()[1];
//         //
//         //     // angle to node
//         //     double angleNode = std::atan2(nodeY-boxCentroidY,
//         //                                 nodeX-boxCentroidX);
//         //     // Stretch/compress at the correct boundaries
//         //     // if (anglex1y0 <= angleNode && anglex1y1 >= angleNode)
//         //     // {
//         //     //     p_node->rGetModifiableLocation()[0] += mSpeed*dt;
//         //     // }
//         //     // if (anglex1y1 <= angleNode && anglex0y1 >= angleNode)
//         //     // {
//         //     //     p_node->rGetModifiableLocation()[1] -= YshrinkSpeed*dt;
//         //     // }
//         //     if (anglex0y0 <= angleNode && anglex1y0 >= angleNode)
//         //     {
//         //         p_node->rGetModifiableLocation()[1] += YshrinkSpeed*dt;
//         //     }
//         // }
//
//         // Shrink box
//         // p_static_cast_mesh_toroidal->SetBoxCoords(0, 0);
//         // p_static_cast_mesh_toroidal->SetBoxCoords(1, currentXUpper+mSpeed*dt);
//         p_static_cast_mesh_toroidal->SetBoxCoords(0, currentXLower-mSpeed*dt);
//         // p_static_cast_mesh_toroidal->SetBoxCoords(2,                                                                     currentYLower+YshrinkSpeed*dt);
//         // p_static_cast_mesh_toroidal->SetBoxCoords(3,                                                                     currentYUpper-YshrinkSpeed*dt);
//         // p_static_cast_mesh_toroidal->RefitPeriodicBox();
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

void ExtrinsicPullModifierToroidal::ApplyExtrinsicPull(bool applyExtrinsicPull)
{
    mApplyExtrinsicPull = applyExtrinsicPull;
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
