#include "ConstantTargetAreaModifier.hpp"

template<unsigned DIM>
ConstantTargetAreaModifier<DIM>::ConstantTargetAreaModifier()
    : AbstractTargetAreaModifier<DIM>()
{
}

template<unsigned DIM>
ConstantTargetAreaModifier<DIM>::~ConstantTargetAreaModifier()
{
}

template<unsigned DIM>
void ConstantTargetAreaModifier<DIM>::UpdateTargetAreaOfCell(CellPtr pCell)
{
    // Get target area A of a healthy cell in S, G2 or M phase
    double cell_target_area = this->mReferenceTargetArea;

    // Set cell data
    pCell->GetCellData()->SetItem("target area", cell_target_area);
}

template<unsigned DIM>
void ConstantTargetAreaModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractTargetAreaModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ConstantTargetAreaModifier<1>;
template class ConstantTargetAreaModifier<2>;
template class ConstantTargetAreaModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ConstantTargetAreaModifier)
