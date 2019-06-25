#include "ExtrinsicPullModifierToroidal.hpp"
#include "Toroidal2dVertexMeshWithMutableSize.hpp"
#include <math.h>

ExtrinsicPullModifierToroidal::ExtrinsicPullModifierToroidal()
    : AbstractCellBasedSimulationModifier<2>(),
      mApplyExtrinsicPullToAllNodes(true),
      mSpeed(1.0),
      mApplyExtrinsicPull(true),

      // For calculating stress
      mHomotypicLineTensionParameter(1.0),
      mHeterotypicLineTensionParameter(1.0),
      mSupercontractileLineTensionParameter(1.0),
      mNumStripes(4),
      mUseCombinedInterfacesForLineTension(false),
      mUseDistinctStripeMismatchesForCombinedInterfaces(false),

      mAreaElasticityParameter(1.0),
      mPerimeterContractilityParameter(0.04),
      mBoundaryLineTensionParameter(1),
      mPullingForce(0.0),
      mGradualPullIncreaseTime(1.0),
      mSimulationEndtime(1.0)

{
}

ExtrinsicPullModifierToroidal::~ExtrinsicPullModifierToroidal()
{
}

void ExtrinsicPullModifierToroidal::UpdateAtEndOfTimeStep(AbstractCellPopulation<2,2>& rCellPopulation)
{
    double dt = SimulationTime::Instance()->GetTimeStep();

    if (mApplyExtrinsicPull)
    {

        // Pointer to mesh
        AbstractMesh<2, 2>& r_mesh = rCellPopulation.rGetMesh();
        Toroidal2dVertexMeshWithMutableSize* p_static_cast_mesh_toroidal = static_cast<Toroidal2dVertexMeshWithMutableSize*>(&r_mesh);

        // // STRETCH LHS BECAUSE THAT'S HOW STRESS RELAXATION RELAXES.
        // // Coords of box
        // double currentXUpper = p_static_cast_mesh_toroidal->GetBoxCoords(1);
        //
        // // Reset the size of the box
        // p_static_cast_mesh_toroidal->SetBoxCoords(1, currentXUpper + mSpeed*dt);

        // Coords of box
        // double currentXLower = p_static_cast_mesh_toroidal->GetBoxCoords(0);
        double currentXUpper = p_static_cast_mesh_toroidal->GetBoxCoords(1);
        // double currentYLower = p_static_cast_mesh_toroidal->GetBoxCoords(2);
        // double currentYUpper = p_static_cast_mesh_toroidal->GetBoxCoords(3);
        // // Centroid of box
        // double boxCentroidX = (currentXLower + currentXUpper)/2;
        // double boxCentroidY = (currentYLower + currentYUpper)/2;
        // // Angles to corners of box. Used to see which side boundary vertices
        // // are on
        // double anglex1y1 = std::atan2(currentYUpper-boxCentroidY,
        //                             currentXUpper-boxCentroidX);
        // double anglex1y0 = std::atan2(currentYLower-boxCentroidY,
        //                             currentXUpper-boxCentroidX);
        //
        // // Get boxwidth
        // double box_height = currentYUpper - currentYLower;

        // Strain that we will pul to find stiffness
        double small_strain = 1 + 0.01;

        // Get stress based on force
        double stress_from_posterior = dt*(mPullingForce); //* box_height;

        // Get original tissue stress
        c_matrix<double, 2, 2> original_tissue_stress = GetTissueStressTensor(rCellPopulation);

        p_static_cast_mesh_toroidal->SetBoxCoords(1, currentXUpper*small_strain);

        // // Get the posterior nodes
        // std::vector<float> posterior_nodes;
        // std::set<unsigned> boundaryNodes = p_static_cast_mesh_toroidal->GetBoundaryNodes();
        // // If it was a RHS boundary node, check where it is and fix it
        // for(auto n_index : boundaryNodes)
        // {
        //     // Get the node
        //     Node<2>* p_node = rCellPopulation.GetNode(n_index);
        //     double nodeX = p_node->rGetLocation()[0];
        //     double nodeY = p_node->rGetLocation()[1];
        //
        //     // Angle to node
        //     double angleNode = std::atan2(nodeY-boxCentroidY, nodeX-boxCentroidX);
        //
        //     // If on RHS
        //     if (angleNode <= anglex1y1 && angleNode >= anglex1y0)
        //     {
        //         // Set its location to the old location
        //         posterior_nodes.push_back(n_index);
        //
        //         p_node->rGetModifiableLocation()[0] = nodeX*small_strain;
        //         p_node->rGetModifiableLocation()[1] = nodeY*small_strain;
        //     }
        // }

        // Get tissue stress after deformation
        // Get original tissue stress
        c_matrix<double, 2, 2> deformed_tissue_stress = GetTissueStressTensor(rCellPopulation);

        // And stiffness
        double stiffness = fabs((deformed_tissue_stress(0,0) - original_tissue_stress(0,0))/small_strain);

        // Calculate strain
        double posterior_strain = stress_from_posterior / stiffness;

        // Calculate posterior pull delay
        double sim_time = SimulationTime::Instance()->GetTime();
        double gradual_scaling = 1;
        if (sim_time < mGradualPullIncreaseTime)
        {
            gradual_scaling = sim_time / mGradualPullIncreaseTime;
        }
        else
        {
            gradual_scaling = 1 - (sim_time - mGradualPullIncreaseTime) / (1.5*mSimulationEndtime - mGradualPullIncreaseTime);
        }


        // Pull the RHS
        p_static_cast_mesh_toroidal->SetBoxCoords(1, currentXUpper*(1+posterior_strain*gradual_scaling));



        // // Move all the nodes appropriately
        // for(auto n_index : posterior_nodes)
        // {
        //     // Get the node
        //     Node<2>* p_node = rCellPopulation.GetNode(n_index);
        //     double nodeX = p_node->rGetLocation()[0];
        //     double nodeY = p_node->rGetLocation()[1];
        //
        //     p_node->rGetModifiableLocation()[0] = (nodeX/small_strain)*(1+posterior_strain);
        //     p_node->rGetModifiableLocation()[1] = (nodeY/small_strain)*(1+posterior_strain);
        // }
    }
}

void ExtrinsicPullModifierToroidal::SetGradualIncreasePullTime(double gradualPullIncreaseTime)
{
    mGradualPullIncreaseTime = gradualPullIncreaseTime;
}

void ExtrinsicPullModifierToroidal::SetSimulationEndTime(double simulationEndtime)
{
    mSimulationEndtime = simulationEndtime;
}

void ExtrinsicPullModifierToroidal::SetPullingForce(double pullingForce)
{
    mPullingForce = pullingForce;
}

void ExtrinsicPullModifierToroidal::SetupSolve(AbstractCellPopulation<2,2>& rCellPopulation, std::string outputDirectory)
{
}

void ExtrinsicPullModifierToroidal::ApplyExtrinsicPullToAllNodes(bool applyExtrinsicPullToAllNodes)
{
    mApplyExtrinsicPullToAllNodes = applyExtrinsicPullToAllNodes;
}

void ExtrinsicPullModifierToroidal::SetSpeed(double speed)
{
    mSpeed = speed;
}

void ExtrinsicPullModifierToroidal::ApplyExtrinsicPull(bool applyExtrinsicPull)
{
    mApplyExtrinsicPull = applyExtrinsicPull;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ExtrinsicPullModifierToroidal::SetAreaElasticityParameter(double areaElasticityParameter)
{
    mAreaElasticityParameter = areaElasticityParameter;
}

double ExtrinsicPullModifierToroidal::GetAreaElasticityParameter()
{
    return mAreaElasticityParameter;
}

void ExtrinsicPullModifierToroidal::SetPerimeterContractilityParameter(double perimeterContractilityParameter)
{
    mPerimeterContractilityParameter = perimeterContractilityParameter;
}

double ExtrinsicPullModifierToroidal::GetPerimeterContractilityParameter()
{
    return mPerimeterContractilityParameter;
}

void ExtrinsicPullModifierToroidal::SetBoundaryLineTensionParameter(double boundaryLineTensionParameter)
{
    mBoundaryLineTensionParameter = boundaryLineTensionParameter;
}

double ExtrinsicPullModifierToroidal::GetBoundaryLineTensionParameter()
{
    return mBoundaryLineTensionParameter;
}




c_matrix<double, 2,2> ExtrinsicPullModifierToroidal::GetSingleCellStressTensor(AbstractCellPopulation<2,2>& rCellPopulation, unsigned elementIndex)
{
    AbstractMesh<2, 2>& r_mesh = rCellPopulation.rGetMesh();
    // Pointer to mesh
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
        double previous_edge_line_tension_parameter = GetLineTensionParameter(p_previous_node, p_this_node, *vertex_population);
        double next_edge_line_tension_parameter = GetLineTensionParameter(p_this_node, p_next_node, *vertex_population);

        // Location of neighbours.
        c_vector<double, 2> previous_node_location = p_element->GetNodeLocation(previous_node_local_index);
        c_vector<double, 2> next_node_location = p_element->GetNodeLocation(next_node_local_index);
        c_vector<double, 2> difference_vector = p_mesh->GetVectorFromAtoB(previous_node_location, next_node_location);

        // Cross product
        c_vector<double, 2> element_area_gradient;
        element_area_gradient[0] = 0.5*difference_vector[1];
        element_area_gradient[1] = -0.5*difference_vector[0];
        // Elasticity term
        double bulk_elasticity_parameter = GetAreaElasticityParameter();
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
        double gamma = GetPerimeterContractilityParameter();
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

    return stressTensor;
}

c_matrix<double, 2, 2> ExtrinsicPullModifierToroidal::GetTissueStressTensor(AbstractCellPopulation<2,2>& rCellPopulation)
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
        cell_stress_tensor = GetSingleCellStressTensor(rCellPopulation, elem_idx);
        tissue_stress_tensor += cell_stress_tensor*cell_area;
    }

    tissue_stress_tensor /= total_area;

    return tissue_stress_tensor;
}


// c_vector<double, 2> ExtrinsicPullModifierToroidal::GetForceContribution(AbstractCellPopulation<2>& rCellPopulation)
// {
//     // Initialise stress tensor
//     c_matrix<double, 2,2> TissueStressTensor = zero_matrix<double>(2,2);
//     // Initialise area
//     double tissue_area = 0;
//
//     // Throw an exception message if not using a VertexBasedCellPopulation
//     ///\todo: check whether this line influences profiling tests - if so, we should remove it.
//     if (dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation) == nullptr)
//     {
//         EXCEPTION("ExtrinsicPullModifierToroidal is to be used with a VertexBasedCellPopulation only");
//     }
//
//     // Define some helper variables
//     VertexBasedCellPopulation<2>* p_cell_population = static_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation);
//     unsigned num_nodes = p_cell_population->GetNumNodes();
//     unsigned num_elements = p_cell_population->GetNumElements();
//
//     // Begin by computing the area and perimeter of each element in the mesh, to avoid having to do this multiple times
//     std::vector<double> element_areas(num_elements);
//     std::vector<double> element_perimeters(num_elements);
//     std::vector<double> target_areas(num_elements);
//     for (typename VertexMesh<2,2>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
//          elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
//          ++elem_iter)
//     {
//         unsigned elem_index = elem_iter->GetIndex();
//         element_areas[elem_index] = p_cell_population->rGetMesh().GetVolumeOfElement(elem_index);
//         element_perimeters[elem_index] = p_cell_population->rGetMesh().GetSurfaceAreaOfElement(elem_index);
//         try
//         {
//             // If we haven't specified a growth modifier, there won't be any target areas in the CellData array and CellData
//             // will throw an exception that it doesn't have "target area" entries.  We add this piece of code to give a more
//             // understandable message. There is a slight chance that the exception is thrown although the error is not about the
//             // target areas.
//             target_areas[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("target area");
//         }
//         catch (Exception&)
//         {
//             EXCEPTION("You need to add an AbstractTargetAreaModifier to the simulation in order to use a FarhadifarForce");
//         }
//
//
//         // Add cell area
//         tissue_area += p_mesh->GetVolumeOfElement(elementIndex);
//         // Cell centroid:
//         c_vector<double, 2> centroid = rCellPopulation.GetLocationOfCellCentre(cell);
//
//
//         // Iterate over vertices in the cell population
//         for (unsigned node_index=0; node_index<num_nodes; node_index++)
//         {
//             Node<2>* p_this_node = p_cell_population->GetNode(node_index);
//
//             /*
//              * The force on this Node is given by the gradient of the total free
//              * energy of the CellPopulation, evaluated at the position of the vertex. This
//              * free energy is the sum of the free energies of all CellPtrs in
//              * the cell population. The free energy of each CellPtr is comprised of three
//              * terms - an area deformation energy, a perimeter deformation energy
//              * and line tension energy.
//              *
//              * Note that since the movement of this Node only affects the free energy
//              * of the CellPtrs containing it, we can just consider the contributions
//              * to the free energy gradient from each of these CellPtrs.
//              */
//             c_vector<double, 2> area_elasticity_contribution = zero_vector<double>(2);
//             c_vector<double, 2> perimeter_contractility_contribution = zero_vector<double>(2);
//             c_vector<double, 2> line_tension_contribution = zero_vector<double>(2);
//
//             // Find the indices of the elements owned by this node
//             std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();
//
//             // Iterate over these elements
//             for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
//                  iter != containing_elem_indices.end();
//                  ++iter)
//             {
//                 // Get this element, its index and its number of nodes
//                 VertexElement<2, 2>* p_element = p_cell_population->GetElement(*iter);
//                 unsigned elem_index = p_element->GetIndex();
//                 unsigned num_nodes_elem = p_element->GetNumNodes();
//
//                 // Find the local index of this node in this element
//                 unsigned local_index = p_element->GetNodeLocalIndex(node_index);
//
//                 // Add the force contribution from this cell's area elasticity (note the minus sign)
//                 c_vector<double, 2> element_area_gradient =
//                         p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);
//                 area_elasticity_contribution -= GetAreaElasticityParameter()*(element_areas[elem_index] -
//                         target_areas[elem_index])*element_area_gradient;
//
//                 // Get the previous and next nodes in this element
//                 unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
//                 Node<2>* p_previous_node = p_element->GetNode(previous_node_local_index);
//
//                 unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
//                 Node<2>* p_next_node = p_element->GetNode(next_node_local_index);
//
//                 // Compute the line tension parameter for each of these edges - be aware that this is half of the actual
//                 // value for internal edges since we are looping over each of the internal edges twice
//                 double previous_edge_line_tension_parameter = GetLineTensionParameter(p_previous_node, p_this_node, *p_cell_population);
//                 double next_edge_line_tension_parameter = GetLineTensionParameter(p_this_node, p_next_node, *p_cell_population);
//
//                 // Compute the gradient of each these edges, computed at the present node
//                 c_vector<double, 2> previous_edge_gradient =
//                         -p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, previous_node_local_index);
//                 c_vector<double, 2> next_edge_gradient = p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);
//
//                 // Add the force contribution from cell-cell and cell-boundary line tension (note the minus sign)
//                 line_tension_contribution -= previous_edge_line_tension_parameter*previous_edge_gradient +
//                         next_edge_line_tension_parameter*next_edge_gradient;
//
//                 // Add the force contribution from this cell's perimeter contractility (note the minus sign)
//                 c_vector<double, 2> element_perimeter_gradient = previous_edge_gradient + next_edge_gradient;
//                 perimeter_contractility_contribution -= GetPerimeterContractilityParameter()* element_perimeters[elem_index]*
//                                                                                                          element_perimeter_gradient;
//             }
//
//             c_vector<double, 2> force_on_node = area_elasticity_contribution + perimeter_contractility_contribution + line_tension_contribution;
//             // p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force_on_node);
//             // return force_on_node;
//
//             // Add to stress tensor
//             stressTensor(0,0) += nodeToCentroid[0]*force_on_node[0];
//             stressTensor(0,1) += nodeToCentroid[0]*force_on_node[1];
//             stressTensor(1,0) += nodeToCentroid[1]*force_on_node[0];
//             stressTensor(1,1) += nodeToCentroid[1]*force_on_node[1];
//         }
//     }
//
// }


double ExtrinsicPullModifierToroidal::GetCombinedInterfaceLength(Node<2>* pNode,
                                                      unsigned elemIndex,
                                                      unsigned cell1StripeIdentity,
                                                      unsigned cell2StripeIdentity,
                                                      VertexBasedCellPopulation<2>& rVertexCellPopulation)
{
    VertexElement<2, 2>* p_elem = rVertexCellPopulation.GetElement(elemIndex);
    unsigned num_nodes = p_elem->GetNumNodes();

    Node<2>* p_prev_node = pNode;
    std::set<unsigned> prev_node_elems = p_prev_node->rGetContainingElementIndices();
    unsigned prev_local_idx = p_elem->GetNodeLocalIndex(pNode->GetIndex());

    double combined_interface_length = 0.0;
    bool part_of_combined_interface = true;

    // We include this condition in the while loop to avoid an infinite loop in
    // the case of a cell surrounded by cells of different stripe identities
    unsigned num_prev_nodes_considered = 0;

    while (part_of_combined_interface && (num_prev_nodes_considered < num_nodes))
    {
        unsigned this_node_index = p_prev_node->GetIndex();

        std::set<unsigned> this_node_elems = prev_node_elems;
        prev_local_idx = (prev_local_idx - 1 + num_nodes)%num_nodes;

        p_prev_node = p_elem->GetNode(prev_local_idx);
        prev_node_elems = p_prev_node->rGetContainingElementIndices();

        part_of_combined_interface = false;

        std::set<unsigned> shared_elems;
        std::set_intersection(prev_node_elems.begin(), prev_node_elems.end(), this_node_elems.begin(), this_node_elems.end(), std::inserter(shared_elems, shared_elems.begin()));

        if (shared_elems.size() == 2)
        {
            CellPtr p_cell_1 = rVertexCellPopulation.GetCellUsingLocationIndex(*(shared_elems.begin()));
            CellPtr p_cell_2 = rVertexCellPopulation.GetCellUsingLocationIndex(*(++(shared_elems.begin())));
            unsigned cell_1_stripe_identity = p_cell_1->GetCellData()->GetItem("stripe");
            unsigned cell_2_stripe_identity = p_cell_2->GetCellData()->GetItem("stripe");

            if (mUseDistinctStripeMismatchesForCombinedInterfaces)
            {
                // Only continue if the neighbouring cells have the same stripes (e.g. 1 and 2) as the interface for which we are calling this method
                if (   ((cell_1_stripe_identity == cell1StripeIdentity) && (cell_2_stripe_identity == cell2StripeIdentity))
                    || ((cell_1_stripe_identity == cell2StripeIdentity) && (cell_2_stripe_identity == cell1StripeIdentity)) )
                {
                    part_of_combined_interface = true;
                    combined_interface_length += rVertexCellPopulation.rGetMesh().GetDistanceBetweenNodes(p_prev_node->GetIndex(), this_node_index);
                }
            }
            else if (cell_1_stripe_identity != cell_2_stripe_identity)
            {
                part_of_combined_interface = true;
                combined_interface_length += rVertexCellPopulation.rGetMesh().GetDistanceBetweenNodes(p_prev_node->GetIndex(), this_node_index);
            }
        }

        num_prev_nodes_considered++;
    }

    Node<2>* p_next_node = pNode;
    std::set<unsigned> next_node_elems = p_next_node->rGetContainingElementIndices();
    unsigned next_local_idx = p_elem->GetNodeLocalIndex(pNode->GetIndex());

    part_of_combined_interface = true;

    // We include this condition in the while loop to avoid an infinite loop in
    // the case of a cell surrounded by cells of different stripe identities
    unsigned num_next_nodes_considered = 0;

    while (part_of_combined_interface && (num_next_nodes_considered < num_nodes))
    {
        unsigned this_node_index = p_next_node->GetIndex();
        std::set<unsigned> this_node_elems = next_node_elems;

        next_local_idx = (next_local_idx + 1)%num_nodes;

        p_next_node = p_elem->GetNode(next_local_idx);
        next_node_elems = p_next_node->rGetContainingElementIndices();

        part_of_combined_interface = false;

        std::set<unsigned> shared_elems;
        std::set_intersection(this_node_elems.begin(), this_node_elems.end(), next_node_elems.begin(), next_node_elems.end(), std::inserter(shared_elems, shared_elems.begin()));

        if (shared_elems.size() == 2)
        {
            CellPtr p_cell_1 = rVertexCellPopulation.GetCellUsingLocationIndex(*(shared_elems.begin()));
            CellPtr p_cell_2 = rVertexCellPopulation.GetCellUsingLocationIndex(*(++(shared_elems.begin())));
            unsigned cell_1_stripe_identity = p_cell_1->GetCellData()->GetItem("stripe");
            unsigned cell_2_stripe_identity = p_cell_2->GetCellData()->GetItem("stripe");

            if (mUseDistinctStripeMismatchesForCombinedInterfaces)
            {
                // Only continue if the neighbouring cells have the same stripes (e.g. 1 and 2) as the interface for which we are calling this method
                if (   ((cell_1_stripe_identity == cell1StripeIdentity) && (cell_2_stripe_identity == cell2StripeIdentity))
                    || ((cell_1_stripe_identity == cell2StripeIdentity) && (cell_2_stripe_identity == cell1StripeIdentity)) )
                {
                    part_of_combined_interface = true;
                    combined_interface_length += rVertexCellPopulation.rGetMesh().GetDistanceBetweenNodes(this_node_index, p_next_node->GetIndex());
                }
            }
            else if (cell_1_stripe_identity != cell_2_stripe_identity)
            {
                part_of_combined_interface = true;
                combined_interface_length += rVertexCellPopulation.rGetMesh().GetDistanceBetweenNodes(this_node_index, p_next_node->GetIndex());
            }
        }
        num_next_nodes_considered++;
    }

    return combined_interface_length;
}


double ExtrinsicPullModifierToroidal::GetCombinedInterfaceScaleFactor(Node<2>* pNodeA,
                                                            Node<2>* pNodeB,
                                                            unsigned element1Index,
                                                            unsigned element2Index,
                                                            unsigned cell1StripeIdentity,
                                                            unsigned cell2StripeIdentity,
                                                            VertexBasedCellPopulation<2>& rVertexCellPopulation)
{
    // Compute the total length of the 'boundary' interface that contains this edge, from each cell's perspective

    // Cell 1
    double interface_length_in_cell_1 = GetCombinedInterfaceLength(pNodeA, element1Index, cell1StripeIdentity, cell2StripeIdentity, rVertexCellPopulation);

    // Cell 2
    double interface_length_in_cell_2 = GetCombinedInterfaceLength(pNodeA, element2Index, cell1StripeIdentity, cell2StripeIdentity, rVertexCellPopulation);

    assert(interface_length_in_cell_1 > 0);
    assert(interface_length_in_cell_2 > 0);

    // Now compute the ratio of the longer to the short 'boundary' interfaces
    double shorter_interface = std::min(interface_length_in_cell_1, interface_length_in_cell_2);
    double scale_factor = 1.0/shorter_interface;

    return scale_factor;
}


double ExtrinsicPullModifierToroidal::GetLineTensionParameter(Node<2>* pNodeA,
                                                    Node<2>* pNodeB,
                                                    VertexBasedCellPopulation<2>& rVertexCellPopulation)
{
    double line_tension = 0.0;

    // If either node is a boundary node, then use mBoundaryLineTensionParameter...
    if (pNodeA->IsBoundaryNode() || pNodeB->IsBoundaryNode())
    {
        line_tension = GetBoundaryLineTensionParameter();
    }
    else
    {
        // ...otherwise find the elements shared by these nodes
        std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
        std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();
        std::set<unsigned> shared_elements;
        std::set_intersection(elements_containing_nodeA.begin(), elements_containing_nodeA.end(),
                              elements_containing_nodeB.begin(), elements_containing_nodeB.end(),
                              std::inserter(shared_elements, shared_elements.begin()));
        assert(!shared_elements.empty());

        // If the nodes share a single edge, then it must be on the boundary, so use mBoundaryLineTensionParameter...
        if (shared_elements.size() == 1)
        {
            line_tension = GetBoundaryLineTensionParameter();
        }
        else
        {
            // ...otherwise proceed
            assert(shared_elements.size() == 2);
            unsigned elem_1_index = *(shared_elements.begin());
            unsigned elem_2_index = *(++(shared_elements.begin()));

            // Find the stripe identities on either side of the edge
            CellPtr p_cell_1 = rVertexCellPopulation.GetCellUsingLocationIndex(elem_1_index);
            CellPtr p_cell_2 = rVertexCellPopulation.GetCellUsingLocationIndex(elem_2_index);

            unsigned cell_1_stripe_identity = p_cell_1->GetCellData()->GetItem("stripe");
            unsigned cell_2_stripe_identity = p_cell_2->GetCellData()->GetItem("stripe");

            // If the edge is within a stripe (an 'internal' interface), then use mHomotypicLineTensionParameter
            if (cell_1_stripe_identity == cell_2_stripe_identity)
            {
                // The factor of 0.5 arises because each internal interface is visited twice in our iteration
                line_tension = 0.5*mHomotypicLineTensionParameter;
            }
            else
            {
                // Label numbers wrap around, so check to find smallest difference in stripe identities
                unsigned mismatch = abs(int(cell_1_stripe_identity - cell_2_stripe_identity));
                if (mismatch > mNumStripes/2)
                {
                    mismatch = mNumStripes - mismatch;
                }

                // If the edge is between stripes whose identities differ by 1, then use mHeterotypicLineTensionParameter
                if (mismatch == 1)
                {
                    // The factor of 0.5 arises because each internal interface is visited twice in our iteration
                    line_tension = 0.5*mHeterotypicLineTensionParameter;
                }
                else
                {
                    // The factor of 0.5 arises because each internal interface is visited twice in our iteration
                    line_tension = 0.5*mSupercontractileLineTensionParameter;
                }

                if (mUseCombinedInterfacesForLineTension)
                {
                    double scale_factor = GetCombinedInterfaceScaleFactor(pNodeA, pNodeB, elem_1_index, elem_2_index, cell_1_stripe_identity, cell_2_stripe_identity, rVertexCellPopulation);
                    line_tension *= scale_factor;
                }
            }
        }
    }

    return line_tension;
}


void ExtrinsicPullModifierToroidal::SetHomotypicLineTensionParameter(double labelledCellLabelledCellAdhesionEnergyParameter)
{
    mHomotypicLineTensionParameter = labelledCellLabelledCellAdhesionEnergyParameter;
}


void ExtrinsicPullModifierToroidal::SetHeterotypicLineTensionParameter(double labelledCellCellAdhesionEnergyParameter)
{
    mHeterotypicLineTensionParameter = labelledCellCellAdhesionEnergyParameter;
}


void ExtrinsicPullModifierToroidal::SetSupercontractileLineTensionParameter(double supercontractileLineTensionParameter)
{
    mSupercontractileLineTensionParameter = supercontractileLineTensionParameter;
}


void ExtrinsicPullModifierToroidal::SetNumStripes(unsigned numStripes)
{
    mNumStripes = numStripes;
}


void ExtrinsicPullModifierToroidal::SetUseCombinedInterfacesForLineTension(bool useCombinedInterfaceLineTension)
{
    mUseCombinedInterfacesForLineTension = useCombinedInterfaceLineTension;
}


void ExtrinsicPullModifierToroidal::SetUseDistinctStripeMismatchesForCombinedInterfaces(bool useDistinctStripeMismatchesForCombinedInterfaces)
{
    mUseDistinctStripeMismatchesForCombinedInterfaces = useDistinctStripeMismatchesForCombinedInterfaces;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void ExtrinsicPullModifierToroidal::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ApplyExtrinsicPullToAllNodes>" << mApplyExtrinsicPullToAllNodes << "</ApplyExtrinsicPullToAllNodes>\n";
    *rParamsFile << "\t\t\t<Speed>" << mSpeed << "</Speed>\n";

    // Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<2>::OutputSimulationModifierParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ExtrinsicPullModifierToroidal)
