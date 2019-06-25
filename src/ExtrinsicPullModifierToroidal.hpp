#ifndef EXTRINSICPULLMODIFIERTOROIDAL_HPP_
#define EXTRINSICPULLMODIFIERTOROIDAL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"

class ExtrinsicPullModifierToroidal : public AbstractCellBasedSimulationModifier<2,2>
{
private:

    bool mApplyExtrinsicPullToAllNodes;
    bool mPinAnteriorMostCells;
    double mSpeed;
    bool mApplyExtrinsicPull;

    double mCellBoundaryAdhesionParameter;
    double mHomotypicLineTensionParameter;
    double mHeterotypicLineTensionParameter;
    double mSupercontractileLineTensionParameter;

    unsigned mNumStripes;

    bool mUseCombinedInterfacesForLineTension;
    bool mUseDistinctStripeMismatchesForCombinedInterfaces;
    double mAreaElasticityParameter;
    double mPerimeterContractilityParameter;
    double mBoundaryLineTensionParameter;
    double mPullingForce;
    double mGradualPullIncreaseTime;
    double mSimulationEndtime;

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

        archive & mHomotypicLineTensionParameter;
        archive & mHeterotypicLineTensionParameter;
        archive & mSupercontractileLineTensionParameter;
        archive & mNumStripes;
        archive & mUseCombinedInterfacesForLineTension;
        archive & mUseDistinctStripeMismatchesForCombinedInterfaces;
        archive & mAreaElasticityParameter;
        archive & mPerimeterContractilityParameter;
        archive & mBoundaryLineTensionParameter;
        archive & mPullingForce;
        archive & mGradualPullIncreaseTime;
        archive & mSimulationEndtime;
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

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void SetGradualIncreasePullTime(double gradualPullIncreaseTime);
    void SetSimulationEndTime(double simulationEndtime);

    void SetAreaElasticityParameter(double areaElasticityParameter);
    double GetAreaElasticityParameter();
    void SetPerimeterContractilityParameter(double perimeterContractilityParameter);
    double GetPerimeterContractilityParameter();
    void SetBoundaryLineTensionParameter(double boundaryLineTensionParameter);
    double GetBoundaryLineTensionParameter();
    void SetPullingForce(double pullingForce);

    c_matrix<double, 2, 2> GetTissueStressTensor(AbstractCellPopulation<2,2>& rCellPopulation);
    c_matrix<double, 2,2> GetSingleCellStressTensor(AbstractCellPopulation<2,2>& rCellPopulation, unsigned elementIndex);

    c_vector<double, 2> GetForceContribution(AbstractCellPopulation<2>& rCellPopulation);

    double GetLineTensionParameter(Node<2>* pNodeA, Node<2>* pNodeB, VertexBasedCellPopulation<2>& rVertexCellPopulation);

    void SetHomotypicLineTensionParameter(double homotypicLineTensionParameter);

    void SetHeterotypicLineTensionParameter(double heterotypicLineTensionParameter);

    void SetSupercontractileLineTensionParameter(double supercontractileLineTensionParameter);

    void SetNumStripes(unsigned numStripes);

    void SetUseCombinedInterfacesForLineTension(bool useCombinedInterfaceLineTension);

    void SetUseDistinctStripeMismatchesForCombinedInterfaces(bool useDistinctStripeMismatchesForCombinedInterfaces);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void OutputForceParameters(out_stream& rParamsFile);

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
