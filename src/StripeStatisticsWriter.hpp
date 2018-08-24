#ifndef STRIPESTATISTICSWRITER_HPP_
#define STRIPESTATISTICSWRITER_HPP_

#include "AbstractCellPopulationWriter.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class StripeStatisticsWriter : public AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    StripeStatisticsWriter();

    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {}

    virtual void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
    {}

    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
    {}

    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
    {}

    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(StripeStatisticsWriter)

#endif /* STRIPESTATISTICSWRITER_HPP_ */