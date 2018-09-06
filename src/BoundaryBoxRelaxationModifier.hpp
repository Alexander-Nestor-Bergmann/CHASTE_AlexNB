#ifndef BOUNDARYBOXRELAXATIONMODIFIER_HPP_
#define BOUNDARYBOXRELAXATIONMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"
#include "FarhadifarForce.hpp"

class BoundaryBoxRelaxationModifier : public AbstractCellBasedSimulationModifier<2,2>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<2,2> >(*this);
        archive & mStiffness;
        archive & mRelaxPeriodicBox;
    }

    double mStiffness;
    bool mRelaxPeriodicBox;
    boost::shared_ptr<FarhadifarForce<2> > mpForce;

public:

    /**
     * Default constructor.
     */
    BoundaryBoxRelaxationModifier();

    /**
     * Destructor.
     */
    virtual ~BoundaryBoxRelaxationModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specifies what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<2,2>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<2,2>& rCellPopulation, std::string outputDirectory);

        /**
     * @param relaxPeriodicBox: whether to allow the box to change shape to dissipate stress.
     */
    void RelaxPeriodicBox(bool relaxPeriodicBox);

    /**
     * @param stiffness the stiffness of the surrounding tissue.
     */
    void SetStiffness(double stiffness);

    /**
     * @param p_force. Pointer to simulation farhadifar force.
    */
    void SetForcePointer(boost::shared_ptr<FarhadifarForce<2> > p_force);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(BoundaryBoxRelaxationModifier)

#endif /*BOUNDARYBOXRELAXATIONMODIFIER_HPP_*/
