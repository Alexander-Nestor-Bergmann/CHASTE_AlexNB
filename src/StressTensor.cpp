#include "StressTensor.hpp"
#include "MutableVertexMesh.hpp"

// StressTensor::StressTensor()
// {
// }
//
//
// StressTensor::~StressTensor()
// {
// }

c_matrix<double, 2,2> GetSingleCellStressTensor(AbstractCellPopulation<2,2>& rCellPopulation, unsigned elementIndex)
{
    AbstractMesh<2, 2>& r_mesh = rCellPopulation.rGetMesh();
    // Pointer to mesh
    MutableVertexMesh<2,2>* p_mesh = static_cast<MutableVertexMesh<2,2>*>(&r_mesh);

    // Initialise stress tensor
    c_matrix<double, 2,2> stressTensor;
    stressTensor(0,0) = 0;
    stressTensor(0,1) = 0;
    stressTensor(1,0) = 0;
    stressTensor(1,1) = 0;

    // Get a pointer to the cell
    CellPtr cell = rCellPopulation.GetCellUsingLocationIndex(elementIndex);
    // Get cell element
    VertexElement<2, 2>* p_element = p_mesh->GetElement(elementIndex);

    // Cell centroid:
    c_vector<double, 2> centroid = rCellPopulation.GetLocationOfCellCentre(cell);

    // Loop over nodes owned by this element
    c_vector<double, 2> force;
    c_vector<double, 2> nodeLocation;
    for (unsigned local_index=0; local_index<p_element->GetNumNodes(); local_index++)
    {
        // Get the node pointer
        Node<2>* p_node = rCellPopulation.GetNode(local_index);
        // Get the location vector from centroid.
        nodeLocation = p_node->rGetLocation();
        // Get the force
        force = p_node->rGetAppliedForce();
        // Add to stress tensor
        stressTensor(0,0) = nodeLocation[0]*force[0];
        stressTensor(0,1) = nodeLocation[0]*force[1];
        stressTensor(1,0) = nodeLocation[1]*force[0];
        stressTensor(1,1) = nodeLocation[1]*force[1];
    }

    return stressTensor;
}

c_matrix<double, 2, 2> GetTissueStressTensor(AbstractCellPopulation<2,2>& rCellPopulation)
{
    AbstractMesh<2, 2>& r_mesh = rCellPopulation.rGetMesh();

    // Pointer to mesh
    MutableVertexMesh<2,2>* p_mesh = static_cast<MutableVertexMesh<2,2>*>(&r_mesh);

    // Initialise stress tensor
    c_matrix<double, 2,2> tissue_stress_tensor = zero_matrix<double>(2,2);

    // Get num cells
    unsigned num_cells = p_mesh->GetNumElements();

    // Iterate over cells and add individual tensors
    c_matrix<double, 2,2> cell_stress_tensor;
    double total_area = 0.0;
    for (unsigned elem_idx = 0 ; elem_idx < num_cells ; elem_idx++)
    {
        // Get area
        double cell_area = p_mesh->GetSurfaceAreaOfElement(elem_idx);
        total_area += cell_area;

        // Add stress tensors
        cell_stress_tensor = GetSingleCellStressTensor(rCellPopulation, elem_idx);
        tissue_stress_tensor += cell_stress_tensor*cell_area;
    }

    tissue_stress_tensor /= total_area;

    return tissue_stress_tensor;
}

// Serialization for Boost >= 1.36
// #include "SerializationExportWrapperForCpp.hpp"
// CHASTE_CLASS_EXPORT(StressTensor)
