#ifndef SidekickBoundaryCondition_HPP_
#define SidekickBoundaryCondition_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

template<unsigned DIM>
class SidekickBoundaryCondition : public AbstractCellPopulationBoundaryCondition<DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<DIM> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     */
    SidekickBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation);

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SidekickBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a SidekickBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const SidekickBoundaryCondition<DIM>* t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
}

/**
 * De-serialize constructor parameters and initialize a SidekickBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, SidekickBoundaryCondition<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)SidekickBoundaryCondition<DIM>(p_cell_population);
}
}
} // namespace ...

#endif /*SidekickBoundaryCondition_HPP_*/
