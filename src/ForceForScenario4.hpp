#ifndef FORCEFORSCENARIO4_HPP_
#define FORCEFORSCENARIO4_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "FarhadifarForce.hpp"
#include <iostream>

template<unsigned DIM>
class ForceForScenario4: public FarhadifarForce<DIM>
{
private:

    double mCellBoundaryAdhesionParameter;
    double mHomotypicLineTensionParameter;
    double mHeterotypicLineTensionParameter;
    double mSupercontractileLineTensionParameter;

    unsigned mNumStripes;

    bool mUseCombinedInterfacesForLineTension;
    bool mUseDistinctStripeMismatchesForCombinedInterfaces;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<FarhadifarForce<DIM> >(*this);
        archive & mHomotypicLineTensionParameter;
        archive & mHeterotypicLineTensionParameter;
        archive & mSupercontractileLineTensionParameter;
        archive & mNumStripes;
        archive & mUseCombinedInterfacesForLineTension;
        archive & mUseDistinctStripeMismatchesForCombinedInterfaces;
    }

    double GetCombinedInterfaceLength(Node<DIM>* pNode,
                                      unsigned elemIndex,
                                      unsigned cell1StripeIdentity,
                                      unsigned cell2StripeIdentity,
                                      VertexBasedCellPopulation<DIM>& rVertexCellPopulation);

    double GetCombinedInterfaceScaleFactor(Node<DIM>* pNodeA,
                                           Node<DIM>* pNodeB,
                                           unsigned element1Index,
                                           unsigned element2Index,
                                           unsigned cell1StripeIdentity,
                                           unsigned cell2StripeIdentity,
                                           VertexBasedCellPopulation<DIM>& rVertexCellPopulation);

public:

    ForceForScenario4();

    ~ForceForScenario4()
    {}

    double GetLineTensionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation);

    void SetHomotypicLineTensionParameter(double homotypicLineTensionParameter);
    void SetHeterotypicLineTensionParameter(double heterotypicLineTensionParameter);
    void SetSupercontractileLineTensionParameter(double supercontractileLineTensionParameter);
    void SetNumStripes(unsigned numStripes);
    void SetUseCombinedInterfacesForLineTension(bool useCombinedInterfaceLineTension);
    void SetUseDistinctStripeMismatchesForCombinedInterfaces(bool useDistinctStripeMismatchesForCombinedInterfaces);
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ForceForScenario4)

#endif /*FORCEFORSCENARIO4_HPP_*/
