#include "Toroidal2dVertexMeshWithMutableSize.hpp"
#include "ForceForScenario4.hpp"
#include "StressTensor.hpp"


ForceForScenario4::ForceForScenario4()
    : FarhadifarForce<2>(),
      mHomotypicLineTensionParameter(1.0),
      mHeterotypicLineTensionParameter(1.0),
      mSupercontractileLineTensionParameter(1.0),
      mNumStripes(4),
      mUseCombinedInterfacesForLineTension(false),
      mUseDistinctStripeMismatchesForCombinedInterfaces(false),
      mAddPosteriorPull(false),
      mPosteriorPull(0.0)
{
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void ForceForScenario4::SetAddPosteriorPull(bool addPosteriorPull)
{
    mAddPosteriorPull = addPosteriorPull;
}


void ForceForScenario4::SetPosteriorPull(double posteriorPull)
{
    mPosteriorPull = posteriorPull;
}



void ForceForScenario4::AddForceContribution(AbstractCellPopulation<2>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    ///\todo: check whether this line influences profiling tests - if so, we should remove it.
    if (dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("ForceForScenario4 is to be used with a VertexBasedCellPopulation only");
    }

    // Define some helper variables
    VertexBasedCellPopulation<2>* p_cell_population = static_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation);
    unsigned num_nodes = p_cell_population->GetNumNodes();
    unsigned num_elements = p_cell_population->GetNumElements();

    // Begin by computing the area and perimeter of each element in the mesh, to avoid having to do this multiple times
    std::vector<double> element_areas(num_elements);
    std::vector<double> element_perimeters(num_elements);
    std::vector<double> target_areas(num_elements);
    for (typename VertexMesh<2,2>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        element_areas[elem_index] = p_cell_population->rGetMesh().GetVolumeOfElement(elem_index);
        element_perimeters[elem_index] = p_cell_population->rGetMesh().GetSurfaceAreaOfElement(elem_index);
        try
        {
            // If we haven't specified a growth modifier, there won't be any target areas in the CellData array and CellData
            // will throw an exception that it doesn't have "target area" entries.  We add this piece of code to give a more
            // understandable message. There is a slight chance that the exception is thrown although the error is not about the
            // target areas.
            target_areas[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("target area");
        }
        catch (Exception&)
        {
            EXCEPTION("You need to add an AbstractTargetAreaModifier to the simulation in order to use a FarhadifarForce");
        }
    }

    // Iterate over vertices in the cell population
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<2>* p_this_node = p_cell_population->GetNode(node_index);

        /*
         * The force on this Node is given by the gradient of the total free
         * energy of the CellPopulation, evaluated at the position of the vertex. This
         * free energy is the sum of the free energies of all CellPtrs in
         * the cell population. The free energy of each CellPtr is comprised of three
         * terms - an area deformation energy, a perimeter deformation energy
         * and line tension energy.
         *
         * Note that since the movement of this Node only affects the free energy
         * of the CellPtrs containing it, we can just consider the contributions
         * to the free energy gradient from each of these CellPtrs.
         */
        c_vector<double, 2> area_elasticity_contribution = zero_vector<double>(2);
        c_vector<double, 2> perimeter_contractility_contribution = zero_vector<double>(2);
        c_vector<double, 2> line_tension_contribution = zero_vector<double>(2);

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<2, 2>* p_element = p_cell_population->GetElement(*iter);
            unsigned elem_index = p_element->GetIndex();
            unsigned num_nodes_elem = p_element->GetNumNodes();

            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);

            // Add the force contribution from this cell's area elasticity (note the minus sign)
            c_vector<double, 2> element_area_gradient =
                    p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);
            area_elasticity_contribution -= this->GetAreaElasticityParameter()*(element_areas[elem_index] -
                    target_areas[elem_index])*element_area_gradient;

            // Get the previous and next nodes in this element
            unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
            Node<2>* p_previous_node = p_element->GetNode(previous_node_local_index);

            unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
            Node<2>* p_next_node = p_element->GetNode(next_node_local_index);

            // Compute the line tension parameter for each of these edges - be aware that this is half of the actual
            // value for internal edges since we are looping over each of the internal edges twice
            double previous_edge_line_tension_parameter = GetLineTensionParameter(p_previous_node, p_this_node, *p_cell_population);
            double next_edge_line_tension_parameter = GetLineTensionParameter(p_this_node, p_next_node, *p_cell_population);

            // Compute the gradient of each these edges, computed at the present node
            c_vector<double, 2> previous_edge_gradient =
                    -p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, previous_node_local_index);
            c_vector<double, 2> next_edge_gradient = p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);

            // Add the force contribution from cell-cell and cell-boundary line tension (note the minus sign)
            line_tension_contribution -= previous_edge_line_tension_parameter*previous_edge_gradient +
                    next_edge_line_tension_parameter*next_edge_gradient;

            // Add the force contribution from this cell's perimeter contractility (note the minus sign)
            c_vector<double, 2> element_perimeter_gradient = previous_edge_gradient + next_edge_gradient;
            perimeter_contractility_contribution -= this->GetPerimeterContractilityParameter()* element_perimeters[elem_index]*
                                                                                                     element_perimeter_gradient;
        }

        c_vector<double, 2> force_on_node = area_elasticity_contribution + perimeter_contractility_contribution + line_tension_contribution;
        p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force_on_node);
    }

    // If posterior pull is active, find boundary and pull
    if (mAddPosteriorPull)
    {

        // Get the cell population and pointer to mesh
        AbstractMesh<2, 2>& r_mesh = rCellPopulation.rGetMesh();
        Toroidal2dVertexMeshWithMutableSize* p_static_cast_mesh_toroidal = static_cast<Toroidal2dVertexMeshWithMutableSize*>(&r_mesh);

        // Coords of box
        double currentXLower = p_static_cast_mesh_toroidal->GetBoxCoords(0);
        double currentXUpper = p_static_cast_mesh_toroidal->GetBoxCoords(1);
        double currentYLower = p_static_cast_mesh_toroidal->GetBoxCoords(2);
        double currentYUpper = p_static_cast_mesh_toroidal->GetBoxCoords(3);

        // Centroid of box
        double boxCentroidX = (currentXLower + currentXUpper)/2;
        double boxCentroidY = (currentYLower + currentYUpper)/2;
        // Angles to corners of box. Used to see which side boundary vertices
        // are on
        double anglex0y0 = std::atan2(currentYLower-boxCentroidY,
                                    currentXLower-boxCentroidX);
        double anglex0y1 = std::atan2(currentYUpper-boxCentroidY,
                                    currentXLower-boxCentroidX);

        // Find the boundary nodes
        std::set<unsigned> boundaryNodes = p_static_cast_mesh_toroidal->GetBoundaryNodes();

        // If it was a boundary node, check where it is and pull it appropriately
        for(auto n_index : boundaryNodes)
        {
            // Get the node
            Node<2>* p_node = p_static_cast_mesh_toroidal->GetNode(n_index);
            double nodeX = p_node->rGetLocation()[0];
            double nodeY = p_node->rGetLocation()[1];

            // angle to node
            double angleNode = std::atan2(nodeY-boxCentroidY,
                                        nodeX-boxCentroidX);
            // Stretch/compress at the correct boundaries. Note, don't need to stretch upper and RHS because they are deermined by height and width of box, so just stretch the box after.
            // LHS
            if ( angleNode <= anglex0y0 || angleNode >= anglex0y1 )
            {
                c_vector<double, 2> force_on_node;
                force_on_node(0) = mPosteriorPull;
                force_on_node(1) = 0;
                p_cell_population->GetNode(n_index)->AddAppliedForceContribution(force_on_node);
            }
        }
    }
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double ForceForScenario4::GetCombinedInterfaceLength(Node<2>* pNode,
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


double ForceForScenario4::GetCombinedInterfaceScaleFactor(Node<2>* pNodeA,
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


double ForceForScenario4::GetLineTensionParameter(Node<2>* pNodeA,
                                                    Node<2>* pNodeB,
                                                    VertexBasedCellPopulation<2>& rVertexCellPopulation)
{
    double line_tension = 0.0;

    // If either node is a boundary node, then use mBoundaryLineTensionParameter...
    if (pNodeA->IsBoundaryNode() || pNodeB->IsBoundaryNode())
    {
        line_tension = this->mBoundaryLineTensionParameter;
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
            line_tension = this->mBoundaryLineTensionParameter;
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


void ForceForScenario4::SetHomotypicLineTensionParameter(double labelledCellLabelledCellAdhesionEnergyParameter)
{
    mHomotypicLineTensionParameter = labelledCellLabelledCellAdhesionEnergyParameter;
}


void ForceForScenario4::SetHeterotypicLineTensionParameter(double labelledCellCellAdhesionEnergyParameter)
{
    mHeterotypicLineTensionParameter = labelledCellCellAdhesionEnergyParameter;
}


void ForceForScenario4::SetSupercontractileLineTensionParameter(double supercontractileLineTensionParameter)
{
    mSupercontractileLineTensionParameter = supercontractileLineTensionParameter;
}


void ForceForScenario4::SetNumStripes(unsigned numStripes)
{
    mNumStripes = numStripes;
}


void ForceForScenario4::SetUseCombinedInterfacesForLineTension(bool useCombinedInterfaceLineTension)
{
    mUseCombinedInterfacesForLineTension = useCombinedInterfaceLineTension;
}


void ForceForScenario4::SetUseDistinctStripeMismatchesForCombinedInterfaces(bool useDistinctStripeMismatchesForCombinedInterfaces)
{
    mUseDistinctStripeMismatchesForCombinedInterfaces = useDistinctStripeMismatchesForCombinedInterfaces;
}


void ForceForScenario4::OutputForceParameters(out_stream& rParamsFile)
{
    // Output member variables then call method on direct parent class
    *rParamsFile << "\t\t\t<HomotypicLineTensionParameter>" << mHomotypicLineTensionParameter << "</HomotypicLineTensionParameter> \n";
    *rParamsFile << "\t\t\t<HeterotypicLineTensionParameter>" << mHeterotypicLineTensionParameter << "</HeterotypicLineTensionParameter> \n";
    *rParamsFile << "\t\t\t<NumStripes>" << mNumStripes << "</NumStripes> \n";
    *rParamsFile << "\t\t\t<UseCombinedInterfacesForLineTension>" << mUseCombinedInterfacesForLineTension << "</UseCombinedInterfacesForLineTension> \n";

    FarhadifarForce<2>::OutputForceParameters(rParamsFile);
}

// // Explicit instantiation
// template class ForceForScenario4<1>;
// template class ForceForScenario4<2>;
// template class ForceForScenario4<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ForceForScenario4)
