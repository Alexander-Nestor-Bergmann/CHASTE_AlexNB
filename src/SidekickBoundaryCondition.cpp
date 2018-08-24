#include "SidekickBoundaryCondition.hpp"

template<unsigned DIM>
SidekickBoundaryCondition<DIM>::SidekickBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation)
    : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation)
{
}

template<unsigned DIM>
void SidekickBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations)
{
    double epsilon = 0.8;

    ChasteCuboid<DIM> bounds = this->mpCellPopulation->rGetMesh().CalculateBoundingBox();
    double x_min = bounds.rGetLowerCorner()[0];
    double x_max = bounds.rGetUpperCorner()[0];
    // double y_max = bounds.rGetUpperCorner()[1];

    // Loop over every node
    for (unsigned node_index=0; node_index<this->mpCellPopulation->GetNumNodes(); node_index++)
    {
        Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

        if (p_node->IsBoundaryNode())
        {
            c_vector<double, DIM> old_node_location;
            old_node_location = rOldLocations.find(p_node)->second;

            // If the node lies on the left, then revert its x coordinate
            if (p_node->rGetLocation()[0] < x_min + epsilon)
            {
                p_node->rGetModifiableLocation()[0] = old_node_location[0];
            }

           // // If the node lies on the right, then revert its x coordinate
           // if (p_node->rGetLocation()[0] > x_max - epsilon)
           // {
           //     p_node->rGetModifiableLocation()[0] = old_node_location[0];
           // }

           // // Try a y flagpole condition on the pulled end.
           // if (p_node->rGetLocation()[1] > y_max - epsilon)
           // {
           //     p_node->rGetModifiableLocation()[1] = old_node_location[1];
           // }
        }
    }}

template<unsigned DIM>
bool SidekickBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
    ///\todo Check the boundary condition is actually satisfied
    bool condition_satisfied = true;
    return condition_satisfied;
}

template<unsigned DIM>
void SidekickBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class SidekickBoundaryCondition<1>;
template class SidekickBoundaryCondition<2>;
template class SidekickBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SidekickBoundaryCondition)
