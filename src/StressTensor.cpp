#include "StressTensor.hpp"
#include "MutableVertexMesh.hpp"
#include "Toroidal2dVertexMeshWithMutableSize.hpp"
#include "FarhadifarForce.hpp"
#include "SmartPointers.hpp"

// StressTensor::StressTensor()
// {
// }
//
//
// StressTensor::~StressTensor()
// {
// }

c_matrix<double, 2,2> GetSingleCellStressTensor(AbstractCellPopulation<2,2>& rCellPopulation, FarhadifarForce<2>* p_force, unsigned elementIndex)
{
    AbstractMesh<2, 2>& r_mesh = rCellPopulation.rGetMesh();
    // Pointer to mesh
    // MutableVertexMesh<2,2>* p_mesh = static_cast<MutableVertexMesh<2,2>*>(&r_mesh);
    Toroidal2dVertexMeshWithMutableSize* p_mesh = static_cast<Toroidal2dVertexMeshWithMutableSize*>(&r_mesh);

    // Cast to vertex-based cell population
    VertexBasedCellPopulation<2>* vertex_population = static_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation);

    // Initialise stress tensor
    c_matrix<double, 2,2> stressTensor = zero_matrix<double>(2,2);

    // Get a pointer to the cell
    CellPtr cell = rCellPopulation.GetCellUsingLocationIndex(elementIndex);
    // Get cell element
    VertexElement<2, 2>* p_element = p_mesh->GetElement(elementIndex);

    // Get area and perim:
    double element_area = p_mesh->GetVolumeOfElement(elementIndex);
    double element_perimeter =
                        p_mesh->GetSurfaceAreaOfElement(elementIndex);
    // Target area
    double target_area = rCellPopulation.GetCellUsingLocationIndex(elementIndex)->GetCellData()->GetItem("target area");

    // Cell centroid:
    c_vector<double, 2> centroid = rCellPopulation.GetLocationOfCellCentre(cell);

    // Num vertices belonging cell
    unsigned numVertices = p_element->GetNumNodes();
    // Loop over nodes owned by this element
    for (unsigned local_index=0; local_index<numVertices; local_index++)
    {
        // Get the global node index
        unsigned global_index = p_element->GetNodeGlobalIndex(local_index);
        // Get the node pointer
        Node<2>* p_this_node = p_mesh->GetNode(global_index);

        // Get the location of the vertex
        c_vector<double, 2> nodeLocation = p_element->GetNodeLocation(local_index);

        // Get the previous and next nodes in this element
        unsigned previous_node_local_index = (numVertices+local_index-1)%numVertices;
        // Next node
        unsigned next_node_local_index = (local_index+1)%numVertices;
        // Pointers to the nodes
        Node<2>* p_previous_node = p_element->GetNode(previous_node_local_index);
        // Next node
        Node<2>* p_next_node = p_element->GetNode(next_node_local_index);
        // Compute the line tension parameter (NOTE this is actually lambda/2) for each of these edges
        double previous_edge_line_tension_parameter = p_force->GetLineTensionParameter(p_previous_node, p_this_node, *vertex_population);
        double next_edge_line_tension_parameter = p_force->GetLineTensionParameter(p_this_node, p_next_node, *vertex_population);

        // Location of neighbours.
        c_vector<double, 2> previous_node_location = p_element->GetNodeLocation(previous_node_local_index);
        c_vector<double, 2> next_node_location = p_element->GetNodeLocation(next_node_local_index);
        c_vector<double, 2> difference_vector = p_mesh->GetVectorFromAtoB(previous_node_location, next_node_location);

        // Cross product
        c_vector<double, 2> element_area_gradient;
        element_area_gradient[0] = 0.5*difference_vector[1];
        element_area_gradient[1] = -0.5*difference_vector[0];
        // Elasticity term
        double bulk_elasticity_parameter = p_force->GetAreaElasticityParameter();
        c_vector<double, 2> area_elasticity_contribution = - bulk_elasticity_parameter*(element_area -
                target_area) * element_area_gradient;

        // length of edges and unit tangents:
        // previous node to this one
        double previous_edge_length = p_mesh->GetDistanceBetweenNodes(global_index, p_element->GetNodeGlobalIndex(previous_node_local_index));
        c_vector<double, 2> previous_edge_gradient = p_mesh->GetVectorFromAtoB(previous_node_location, nodeLocation)/previous_edge_length;
        // this node to next
        double next_edge_length = p_mesh->GetDistanceBetweenNodes(global_index, p_element->GetNodeGlobalIndex(next_node_local_index));
        c_vector<double, 2> next_edge_gradient = p_mesh->GetVectorFromAtoB(nodeLocation, next_node_location)/next_edge_length;

        // Add the force contribution from cell-cell line tension (this has already been halved so no need for the factor of 1/2)
        c_vector<double, 2> line_tension_contribution =                next_edge_line_tension_parameter*next_edge_gradient - previous_edge_line_tension_parameter*previous_edge_gradient;

        // Add the force contribution from this cell's perimeter contractility
        c_vector<double, 2> element_perimeter_gradient = next_edge_gradient - previous_edge_gradient;
        double gamma = p_force->GetPerimeterContractilityParameter();
        c_vector<double, 2> perimeter_contractility_contribution = gamma* element_perimeter*element_perimeter_gradient;

        // Add the force contributions
        c_vector<double, 2> force_on_node;
        force_on_node = area_elasticity_contribution + perimeter_contractility_contribution + line_tension_contribution;

        // Vector to centroid, accounting for periodicity.
        c_vector<double, 2> nodeToCentroid = p_mesh->GetVectorFromAtoB(centroid, nodeLocation);

        // Add to stress tensor
        stressTensor(0,0) += nodeToCentroid[0]*force_on_node[0];
        stressTensor(0,1) += nodeToCentroid[0]*force_on_node[1];
        stressTensor(1,0) += nodeToCentroid[1]*force_on_node[0];
        stressTensor(1,1) += nodeToCentroid[1]*force_on_node[1];
    }

    stressTensor /= element_area;

    // // Test of the stress tensor to compare with P_eff
    // std::cerr << "/* error message */" << '\n';
    // cout << stressTensor(0,0) << ", " << stressTensor(0,1) << ",\n" << stressTensor(1,0) << ", " << stressTensor(1,1) << "\n";
    //
    // cout << "L< G (" << lambda1 << ", " << gamma1 << ")" << '\n';
    //
    // double L0 = -(2*lambda1)/(2*gamma1);
    // double T = gamma1*(element_perimeter - L0);
    // double pEff = element_area - 1 + 0.5*T*(element_perimeter/element_area);
    //
    // std::cout << "PEFF " << pEff << '\n';
    // std::cerr << "/* error message */" << '\n';

    return stressTensor;
}

c_matrix<double, 2, 2> GetTissueStressTensor(AbstractCellPopulation<2,2>& rCellPopulation, FarhadifarForce<2>* p_force)
{
    AbstractMesh<2, 2>& r_mesh = rCellPopulation.rGetMesh();

    // Pointer to mesh
    // MutableVertexMesh<2,2>* p_mesh = static_cast<MutableVertexMesh<2,2>*>(&r_mesh);
    Toroidal2dVertexMeshWithMutableSize* p_mesh = static_cast<Toroidal2dVertexMeshWithMutableSize*>(&r_mesh);

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
        double cell_area = p_mesh->GetVolumeOfElement(elem_idx);
        total_area += cell_area;

        // Add stress tensors
        cell_stress_tensor = GetSingleCellStressTensor(rCellPopulation, p_force, elem_idx);
        tissue_stress_tensor += cell_stress_tensor*cell_area;
    }

    tissue_stress_tensor /= total_area;

    return tissue_stress_tensor;
}

// Serialization for Boost >= 1.36
// #include "SerializationExportWrapperForCpp.hpp"
// CHASTE_CLASS_EXPORT(StressTensor)


// // ALTERNATIVE METHOD TO GET FORCE ON VERTEX
// // Add the force contribution from the area elasticity (note the minus sign)
// double bulk_elasticity_parameter = 1;
// c_vector<double, 2> element_area_gradient =
//         p_mesh->GetAreaGradientOfElementAtNode(p_element, local_index);
// c_vector<double, 2> area_elasticity_contribution =  bulk_elasticity_parameter*(element_area -
//         target_area) * element_area_gradient;
//
// // Get the previous and next nodes in this element
// unsigned previous_node_local_index = (numVertices+local_index-1)%numVertices;
// // Node<2>* p_previous_node = p_element->GetNode(previous_node_local_index);
// // Next node
// unsigned next_node_local_index = (local_index+1)%numVertices;
// // Node<2>* p_next_node = p_element->GetNode(next_node_local_index);
//
// // Compute the line tension parameter for each of these edges
// double previous_edge_line_tension_parameter = lambda;
// double next_edge_line_tension_parameter = lambda;
//
// // Compute the gradient of each these edges, computed at the present node
// c_vector<double, 2> previous_edge_gradient =
//         -p_mesh->GetNextEdgeGradientOfElementAtNode(p_element, previous_node_local_index);
// c_vector<double, 2> next_edge_gradient = p_mesh->GetNextEdgeGradientOfElementAtNode(p_element, local_index);
//
// // Add the force contribution from cell-cell and cell-boundary line tension (note the minus sign)
// c_vector<double, 2> line_tension_contribution = 0.5 *  (previous_edge_line_tension_parameter*previous_edge_gradient +
//         next_edge_line_tension_parameter*next_edge_gradient);
//
// // Add the force contribution from this cell's perimeter contractility (note the minus sign)
// c_vector<double, 2> element_perimeter_gradient = previous_edge_gradient + next_edge_gradient;
// c_vector<double, 2> perimeter_contractility_contribution =  gamma* element_perimeter*element_perimeter_gradient;
