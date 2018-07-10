#ifndef CONSTANTTARGETAREAMODIFIER_HPP_
#define CONSTANTTARGETAREAMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include "AbstractTargetAreaModifier.hpp"

/**
 * A modifier class in which the target area property of each cell is held constant over time.
 */
template<unsigned DIM>
class ConstantTargetAreaModifier : public AbstractTargetAreaModifier<DIM>
{
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
        archive & boost::serialization::base_object<AbstractTargetAreaModifier<DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    ConstantTargetAreaModifier();

    /**
     * Destructor.
     */
    virtual ~ConstantTargetAreaModifier();

    /**
     * Overridden UpdateTargetAreaOfCell() method.
     *
     * @param pCell pointer to the cell
     */
    virtual void UpdateTargetAreaOfCell(const CellPtr pCell);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ConstantTargetAreaModifier)

#endif /*CONSTANTTARGETAREAMODIFIER_HPP_*/
