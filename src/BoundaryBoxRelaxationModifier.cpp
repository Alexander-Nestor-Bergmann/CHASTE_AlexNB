#include "BoundaryBoxRelaxationModifier.hpp"
#include "Toroidal2dVertexMeshWithMutableSize.hpp"
#include "StressTensor.hpp"

BoundaryBoxRelaxationModifier::BoundaryBoxRelaxationModifier()
    : AbstractCellBasedSimulationModifier<2>(),
      mStiffness(1000),
      mRelaxPeriodicBox(true),
      mpForce(nullptr)
{
}

BoundaryBoxRelaxationModifier::~BoundaryBoxRelaxationModifier()
{
}

void BoundaryBoxRelaxationModifier::UpdateAtEndOfTimeStep(AbstractCellPopulation<2,2>& rCellPopulation)
{
    if (mRelaxPeriodicBox)
    {
        // Pointer to mesh
        AbstractMesh<2, 2>& r_mesh = rCellPopulation.rGetMesh();
        Toroidal2dVertexMeshWithMutableSize* p_mesh = static_cast<Toroidal2dVertexMeshWithMutableSize*>(&r_mesh);

        // Coords of box
        double currentXLower = p_mesh->GetBoxCoords(0);
        double currentXUpper = p_mesh->GetBoxCoords(1);
        double currentYLower = p_mesh->GetBoxCoords(2);
        double currentYUpper = p_mesh->GetBoxCoords(3);

        // StressTensor
        c_matrix<double, 2,2> stressTensor2d = GetTissueStressTensor(rCellPopulation, mpForce.get());

        // Force in x (stress * box height)
        double xForce = - stressTensor2d(0,0) * p_mesh->GetWidth(1);
        // Force in y (stress * box width)
        double yForce = - stressTensor2d(1,1) * p_mesh->GetWidth(0);

        // Displacement of box is force/stiffness
        double deltaX = xForce / mStiffness;
        double deltaY = yForce / mStiffness;

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
        //                             currentXUpper-boxCentroidX);
        double anglex1y0 = std::atan2(currentYLower-boxCentroidY,
                                    currentXUpper-boxCentroidX);

        // Find the boundary nodes
        std::set<unsigned> boundaryNodes = p_mesh->GetBoundaryNodes();

        // If it was a boundary node, check where it is and pull it appropriately
        for(auto n_index : boundaryNodes)
        {
            // Get the node
            Node<2>* p_node = p_mesh->GetNode(n_index);
            double nodeX = p_node->rGetLocation()[0];
            double nodeY = p_node->rGetLocation()[1];

            // angle to node
            double angleNode = std::atan2(nodeY-boxCentroidY,
                                        nodeX-boxCentroidX);
            // Stretch/compress at the correct boundaries. Note, don't need to stretch upper and RHS because they are deermined by height and width of box, so just stretch the box after.
            // Bottom
            if (anglex0y0 <= angleNode && anglex1y0 >= angleNode)
            {
                p_node->rGetModifiableLocation()[1] += deltaY;
            }
            // LHS
            if ( angleNode <= anglex0y0 || angleNode >= anglex0y1 )
            {
                p_node->rGetModifiableLocation()[0] += deltaX;
            }
        }
        // Reset the size of the box
        p_mesh->SetBoxCoords(0, currentXLower + deltaX);
        p_mesh->SetBoxCoords(1, currentXUpper - deltaX);
        p_mesh->SetBoxCoords(2, currentYLower + deltaY);
        p_mesh->SetBoxCoords(3, currentYUpper - deltaY);
    }
}

void BoundaryBoxRelaxationModifier::SetupSolve(AbstractCellPopulation<2,2>& rCellPopulation, std::string outputDirectory)
{
}

void BoundaryBoxRelaxationModifier::RelaxPeriodicBox(bool relaxPeriodicBox)
{
    mRelaxPeriodicBox = relaxPeriodicBox;
}

void BoundaryBoxRelaxationModifier::SetStiffness(double stiffness)
{
    mStiffness = stiffness;
}

void BoundaryBoxRelaxationModifier::SetForcePointer(boost::shared_ptr<FarhadifarForce<2> > p_force)
{
    mpForce = p_force;
}

void BoundaryBoxRelaxationModifier::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<RelaxPeriodicBox>" << mRelaxPeriodicBox << "</RelaxPeriodicBox>\n";
    *rParamsFile << "\t\t\t<Stiffness>" << mStiffness << "</Stiffness>\n";

    // Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<2>::OutputSimulationModifierParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(BoundaryBoxRelaxationModifier)
