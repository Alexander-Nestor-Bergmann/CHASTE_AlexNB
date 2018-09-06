#ifndef STRESSTENSOR_HPP_
#define STRESSTENSOR_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "MutableVertexMesh.hpp"
#include "AbstractCellPopulation.hpp"
#include "FarhadifarForce.hpp"

//
// class StressTensor
// {
// private:
//
//     /** Needed for serialization. */
//     friend class boost::serialization::access;
//     /**
//      * Boost Serialization method for archiving/checkpointing.
//      * Archives the object and its member variables.
//      *
//      * @param archive  The boost archive.
//      * @param version  The current version of this class.
//      */
//     // template<class Archive>
//     // void serialize(Archive & archive, const unsigned int version)
//     // {
//     //     archive & boost::serialization::base_object<MutableVertexMesh<2,2> >(*this);
//     // }
//
// public:
//
//     /**
//      * Default constructor.
//      */
//     StressTensor();
//
//     /**
//      * Destructor.
//      */
//     virtual ~StressTensor();

    /**
     * GetSingleCellStressTensor() method.
     *
     * Gets the stress tensor for a specified cell.
     *
     c_vector<double, 2> applied_force_0;
     applied_force_0 = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
     * @param rCellPopulation reference to the cell population
     * @param elementIndex reference to the cell population
     */
    c_matrix<double, 2,2> GetSingleCellStressTensor(AbstractCellPopulation<2,2>& rCellPopulation, FarhadifarForce<2>* p_force, unsigned elementIndex);

    /**
     * GetTissueStressTensor() method.
     *
     * Gets the stress tensor for the tissue.
     *
     * @param rCellPopulation reference to the cell population
     */
    c_matrix<double, 2,2> GetTissueStressTensor(AbstractCellPopulation<2,2>& rCellPopulation, FarhadifarForce<2>* p_force);

// };

// #include "SerializationExportWrapper.hpp"
// CHASTE_CLASS_EXPORT(StressTensor)

#endif /*STRESSTENSOR_HPP_*/
