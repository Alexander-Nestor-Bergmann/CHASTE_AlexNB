#ifndef EXTRINSICPULLMODIFIERTOROIDAL_HPP_
#define EXTRINSICPULLMODIFIERTOROIDAL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"

class ExtrinsicPullModifierToroidal : public AbstractCellBasedSimulationModifier<2,2>
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
        archive & mApplyExtrinsicPullToAllNodes;
        archive & mPinAnteriorMostCells;
        archive & mSpeed;
        archive & mApplyExtrinsicPull;
    }

    bool mApplyExtrinsicPullToAllNodes;
    bool mPinAnteriorMostCells;
    double mSpeed;
    bool mApplyExtrinsicPull;

public:

    /**
     * Default constructor.
     */
    ExtrinsicPullModifierToroidal();

    /**
     * Destructor.
     */
    virtual ~ExtrinsicPullModifierToroidal();

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
     * @param applyExtrinsicPullToAllNodes whether to apply the extrinsic pull to all nodes in the tissue
     */
    void ApplyExtrinsicPullToAllNodes(bool applyExtrinsicPullToAllNodes);

    /**
     * @param speed the speed of the extrinsic pull
     */
    void SetSpeed(double speed);

    /**
     * @param applyExtrinsicPull whether to apply the extrinsic pull at all.
     */
    void ApplyExtrinsicPull(bool applyExtrinsicPullToAllNodes);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(ExtrinsicPullModifierToroidal)

#endif /*EXTRINSICPULLMODIFIERTOROIDAL_HPP_*/
