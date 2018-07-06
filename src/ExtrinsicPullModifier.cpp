#include "ExtrinsicPullModifier.hpp"

template<unsigned DIM>
ExtrinsicPullModifier<DIM>::ExtrinsicPullModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mApplyExtrinsicPullToAllNodes(true),
      mSpeed(1.0)
{
}

template<unsigned DIM>
ExtrinsicPullModifier<DIM>::~ExtrinsicPullModifier()
{
}

template<unsigned DIM>
void ExtrinsicPullModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    double epsilon = 0.8;

    double dt = SimulationTime::Instance()->GetTimeStep();
    unsigned num_nodes = rCellPopulation.GetNumNodes();
    ChasteCuboid<DIM> bounds = rCellPopulation.rGetMesh().CalculateBoundingBox();
    double x_min = bounds.rGetLowerCorner()[0];
    double x_max = bounds.rGetUpperCorner()[0];

    if (mApplyExtrinsicPullToAllNodes)
    {
        // Pull on all nodes, with a constant strain rate
        double width = x_max - x_min;
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            Node<DIM>* p_node = rCellPopulation.GetNode(node_index);
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
        // Pull on the right-most nodes only, with a constant speed
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            Node<DIM>* p_node = rCellPopulation.GetNode(node_index);
            if (fabs(p_node->rGetLocation()[0] - x_max) < 0.1)
            {
                p_node->rGetModifiableLocation()[0] += mSpeed*dt;
            }
        }
    }
}

template<unsigned DIM>
void ExtrinsicPullModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
}

template<unsigned DIM>
void ExtrinsicPullModifier<DIM>::ApplyExtrinsicPullToAllNodes(bool applyExtrinsicPullToAllNodes)
{
    mApplyExtrinsicPullToAllNodes = applyExtrinsicPullToAllNodes;
}

template<unsigned DIM>
void ExtrinsicPullModifier<DIM>::SetSpeed(double speed)
{
    mSpeed = speed;
}

template<unsigned DIM>
void ExtrinsicPullModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ApplyExtrinsicPullToAllNodes>" << mApplyExtrinsicPullToAllNodes << "</ApplyExtrinsicPullToAllNodes>\n";
    *rParamsFile << "\t\t\t<Speed>" << mSpeed << "</Speed>\n";

    // Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ExtrinsicPullModifier<1>;
template class ExtrinsicPullModifier<2>;
template class ExtrinsicPullModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ExtrinsicPullModifier)
