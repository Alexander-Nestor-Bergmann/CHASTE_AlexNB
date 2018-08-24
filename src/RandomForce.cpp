#include "RandomForce.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned DIM>
RandomForce<DIM>::RandomForce(double diffusionConstant)
    : AbstractForce<DIM>(),
      mDiffusionConstant(diffusionConstant)
{
}

template<unsigned DIM>
RandomForce<DIM>::~RandomForce()
{
}

template<unsigned DIM>
void RandomForce<DIM>::SetDiffusionConstant(double diffusionConstant)
{
    mDiffusionConstant = diffusionConstant;
}

template<unsigned DIM>
void RandomForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("RandomForce is to be used with a VertexBasedCellPopulation only");
    }

    // Helper variable that is a static cast of the cell population
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    double dt = SimulationTime::Instance()->GetTimeStep();

    // Iterate over vertices in the cell population
    for (unsigned node_index=0; node_index<p_cell_population->GetNumNodes(); node_index++)
    {
        c_vector<double, DIM> force_contribution;
        for (unsigned i=0; i<DIM; i++)
        {
            double nu = p_cell_population->GetDampingConstant(node_index);
            double xi = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
            force_contribution[i] = nu*sqrt(2.0*mDiffusionConstant/dt)*xi;
        }
        rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force_contribution);
    }
}

template<unsigned DIM>
void RandomForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Output member variable then call method on direct parent class
    *rParamsFile << "\t\t\t<DiffusionConstant>" << mDiffusionConstant << "</DiffusionConstant> \n";

    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class RandomForce<1>;
template class RandomForce<2>;
template class RandomForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomForce)
