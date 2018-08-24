#ifndef RANDOMFORCE_HPP_
#define RANDOMFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * A force class to model the random motion of vertices.
 * For use with a VertexBasedCellPopulation only.
 */
template<unsigned DIM>
class RandomForce : public AbstractForce<DIM>
{
private:

    /** Diffusion constant */
    double mDiffusionConstant;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mDiffusionConstant;
    }

public :

    /**
     * Constructor.
     *
     * @params diffusionConstant the value to assign to the diffusion constant (defaults to 0.01)
     */
    RandomForce(double diffusionConstant=0.01);

    /**
     * Destructor.
     */
    ~RandomForce();

    /*
     * Method to set the diffusion constant.
     *
     * @params diffusionConstant the value to assign to the diffusion constant.
     */
    void SetDiffusionConstant(double diffusionConstant);

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the tissue
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomForce)

#endif /*RANDOMFORCE_HPP_*/
