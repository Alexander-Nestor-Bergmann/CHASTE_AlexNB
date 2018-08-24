#ifndef SidekickBoundaryConditionToroidal_HPP_
#define SidekickBoundaryConditionToroidal_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>


class SidekickBoundaryConditionToroidal : public AbstractCellPopulationBoundaryCondition<2>
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
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<2> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     */
    SidekickBoundaryConditionToroidal(AbstractCellPopulation<2>* pCellPopulation);

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations);

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

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a SidekickBoundaryConditionToroidal.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const SidekickBoundaryConditionToroidal* t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<2>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
}

/**
 * De-serialize constructor parameters and initialize a SidekickBoundaryConditionToroidal.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, SidekickBoundaryConditionToroidal* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)SidekickBoundaryConditionToroidal(p_cell_population);
}
}
} // namespace ...

#endif /*SidekickBoundaryConditionToroidal_HPP_*/
