#include "StripeStatisticsWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
StripeStatisticsWriter<ELEMENT_DIM, SPACE_DIM>::StripeStatisticsWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("stripestatistics.dat")
{
}

// We neglect boundary edges here
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StripeStatisticsWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    // Make sure the cell population is updated
    pCellPopulation->Update();

    // Initialise helper variables
    double total_num_edges = 0.0;
    double total_edges_length = 0.0;
    double mismatch_one_num_edges = 0.0;
    double mismatch_one_boundary_length = 0.0;
    double mismatch_two_num_edges = 0.0;
    double mismatch_two_boundary_length = 0.0;

    // Iterate over cells
    for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
         cell_iter != pCellPopulation->End();
         ++cell_iter)
    {
        // Find this cell's stripe identity
        unsigned cell_stripe_identity = cell_iter->GetCellData()->GetItem("stripe");

        // Get the set of neighbouring element indices
        unsigned elem_index = pCellPopulation->GetLocationIndexUsingCell(*cell_iter);
        std::set<unsigned> neighbour_elem_indices = pCellPopulation->rGetMesh().GetNeighbouringElementIndices(elem_index);

        // Iterate over these neighbours
        for (std::set<unsigned>::iterator neighbour_iter = neighbour_elem_indices.begin();
             neighbour_iter != neighbour_elem_indices.end();
             ++neighbour_iter)
        {
            // Get the length of the edge shared with this neighbour
            unsigned neighbour_index = *neighbour_iter;
            double edge_length = pCellPopulation->rGetMesh().GetEdgeLength(elem_index, neighbour_index);

            total_edges_length += edge_length;
            total_num_edges += 1.0;

            // Find this neighbour's stripe identity
            CellPtr p_neighbour = pCellPopulation->GetCellUsingLocationIndex(*neighbour_iter);
            unsigned neighbour_stripe_identity = p_neighbour->GetCellData()->GetItem("stripe");

            unsigned num_stripes = 4; ///\todo remove hardcoding
            unsigned mismatch = abs(int(cell_stripe_identity - neighbour_stripe_identity));
            if (mismatch > num_stripes/2)
            {
                mismatch = num_stripes - mismatch;
            }

            if (mismatch == 1)
            {
                mismatch_one_num_edges += 1.0;
                mismatch_one_boundary_length += edge_length;
            }
            else if (mismatch == 2)
            {
                mismatch_two_num_edges += 1.0;
                mismatch_two_boundary_length += edge_length;
            }
            else
            {
                assert(mismatch == 0);
            }
        }
    }

    // We have counted each cell-cell edge twice
    total_num_edges *= 0.5;
    total_edges_length *= 0.5;
    mismatch_one_num_edges *= 0.5;
    mismatch_one_boundary_length *= 0.5;
    mismatch_two_num_edges *= 0.5;
    mismatch_two_boundary_length *= 0.5;

    *this->mpOutStream << total_num_edges << "\t" << total_edges_length
                       << "\t" << mismatch_one_num_edges << "\t" << mismatch_one_boundary_length
                       << "\t" << mismatch_two_num_edges << "\t" << mismatch_two_boundary_length;
}

// Explicit instantiation
template class StripeStatisticsWriter<1,1>;
template class StripeStatisticsWriter<1,2>;
template class StripeStatisticsWriter<2,2>;
template class StripeStatisticsWriter<1,3>;
template class StripeStatisticsWriter<2,3>;
template class StripeStatisticsWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(StripeStatisticsWriter)
