#ifndef FORCEFORSCENARIO4_HPP_
#define FORCEFORSCENARIO4_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "FarhadifarForce.hpp"
#include <iostream>


class ForceForScenario4: public FarhadifarForce<2>
{
private:

    double mCellBoundaryAdhesionParameter;
    double mHomotypicLineTensionParameter;
    double mHeterotypicLineTensionParameter;
    double mSupercontractileLineTensionParameter;

    unsigned mNumStripes;

    bool mUseCombinedInterfacesForLineTension;
    bool mUseDistinctStripeMismatchesForCombinedInterfaces;

    bool mAddPosteriorPull;
    double mPosteriorPull;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<FarhadifarForce<2> >(*this);
        archive & mHomotypicLineTensionParameter;
        archive & mHeterotypicLineTensionParameter;
        archive & mSupercontractileLineTensionParameter;
        archive & mNumStripes;
        archive & mUseCombinedInterfacesForLineTension;
        archive & mUseDistinctStripeMismatchesForCombinedInterfaces;

        archive & mAddPosteriorPull;
        archive & mPosteriorPull;
    }

    double GetCombinedInterfaceLength(Node<2>* pNode,
                                      unsigned elemIndex,
                                      unsigned cell1StripeIdentity,
                                      unsigned cell2StripeIdentity,
                                      VertexBasedCellPopulation<2>& rVertexCellPopulation);

    double GetCombinedInterfaceScaleFactor(Node<2>* pNodeA,
                                           Node<2>* pNodeB,
                                           unsigned element1Index,
                                           unsigned element2Index,
                                           unsigned cell1StripeIdentity,
                                           unsigned cell2StripeIdentity,
                                           VertexBasedCellPopulation<2>& rVertexCellPopulation);

public:

    ForceForScenario4();

    ~ForceForScenario4()
    {}

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the force on each node in the vertex-based cell population based on the energy function
     * Farhadifar's model.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<2>& rCellPopulation);
    void SetAddPosteriorPull(bool addPosteriorPull);
    void SetPosteriorPull(double posteriorPull);



    double GetLineTensionParameter(Node<2>* pNodeA, Node<2>* pNodeB, VertexBasedCellPopulation<2>& rVertexCellPopulation);

    void SetHomotypicLineTensionParameter(double homotypicLineTensionParameter);
    void SetHeterotypicLineTensionParameter(double heterotypicLineTensionParameter);
    void SetSupercontractileLineTensionParameter(double supercontractileLineTensionParameter);
    void SetNumStripes(unsigned numStripes);
    void SetUseCombinedInterfacesForLineTension(bool useCombinedInterfaceLineTension);
    void SetUseDistinctStripeMismatchesForCombinedInterfaces(bool useDistinctStripeMismatchesForCombinedInterfaces);
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(ForceForScenario4)

#endif /*FORCEFORSCENARIO4_HPP_*/
