#include "Toroidal2dVertexMeshWithMutableSize.hpp"
#include "AbstractCellPopulation.hpp"
#include "StressTensor.hpp"

Toroidal2dVertexMeshWithMutableSize::Toroidal2dVertexMeshWithMutableSize
                                           (double width,
                                           double height,
                                           std::vector<Node<2>*> nodes,
                                           std::vector<VertexElement<2, 2>*> vertexElements,
                                           double cellRearrangementThreshold,
                                           double t2Threshold)
    : MutableVertexMesh<2,2>(nodes, vertexElements, cellRearrangementThreshold, t2Threshold),
      mWidth(width),
      mHeight(height),
      mBoxLowerY(0.0),
      mBoxUpperY(height),
      mBoxLowerX(0.0),
      mBoxUpperX(width),
      mReferenceStress(zero_matrix<double>(2,2)),
      mpMeshForVtk(nullptr)
{
    // Call ReMesh() to remove any deleted nodes and relabel
    ReMesh();
}

Toroidal2dVertexMeshWithMutableSize::Toroidal2dVertexMeshWithMutableSize()
{
}

Toroidal2dVertexMeshWithMutableSize::~Toroidal2dVertexMeshWithMutableSize()
{
    delete mpMeshForVtk;
}

c_vector<double, 2> Toroidal2dVertexMeshWithMutableSize::GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2)
{
    assert(mWidth > 0.0);
    assert(mHeight > 0.0);
    // assert(mBoxUpperX - mBoxLowerX == mWidth)
    // assert(mBoxUpperY - mBoxLowerY == mHeight)

    c_vector<double, 2> vector = rLocation2 - rLocation1;
    vector[0] = fmod(vector[0], mWidth);
    vector[1] = fmod(vector[1], mHeight);

    // If the points are more than halfway across the domain, measure the other way
    if (vector[0] > 0.5*mWidth)
    {
        vector[0] -= mWidth;
    }
    else if (vector[0] < -0.5*mWidth)
    {
        vector[0] += mWidth;
    }

    // If the points are more than halfway up the domain, measure the other way
    if (vector[1] > 0.5*mHeight)
    {
        vector[1] -= mHeight;
    }
    else if (vector[1] < -0.5*mHeight)
    {
        vector[1] += mHeight;
    }
    return vector;
}

void Toroidal2dVertexMeshWithMutableSize::SetNode(unsigned nodeIndex, ChastePoint<2> point)
{
    double x_coord = point.rGetLocation()[0];
    double y_coord = point.rGetLocation()[1];

    // Perform a periodic movement if necessary
    if (x_coord > mBoxUpperX)
    {
        // Move point left
        // point.SetCoordinate(0, x_coord - mWidth);
    }
    else if (x_coord < mBoxLowerX)
    {
        // Move point right
        // point.SetCoordinate(0, x_coord + mWidth);

        // Fix the left-hand side at within box.
        // point.SetCoordinate(0, mBoxLowerX);
    }
    if (y_coord > mBoxUpperY)
    {
        // Move point down
        point.SetCoordinate(1, y_coord - mHeight);
    }
    else if (y_coord < mBoxLowerY)
    {
        // // Move point up
        point.SetCoordinate(1, y_coord + mHeight);

        // Fix the bottom at box.
        // point.SetCoordinate(1, mBoxLowerY);
    }

    // Update the node's location
    MutableVertexMesh<2,2>::SetNode(nodeIndex, point);
}

double Toroidal2dVertexMeshWithMutableSize::GetWidth(const unsigned& rDimension) const
{
    assert(rDimension==0 || rDimension==1);

    double width = mWidth;
    if (rDimension == 1)
    {
        width = mHeight;
    }

    return width;
}

void Toroidal2dVertexMeshWithMutableSize::SetWidth(const unsigned& rDimension, double width)
{
    // This function assumes that the box shrinks from top or RHS only.
    assert(rDimension==0 || rDimension==1);

    if (rDimension == 0)
    {
        mWidth = width;
        // mBoxUpperX = width;
    }
    else
    {
        mHeight = width;
        // mBoxUpperY = width;
    }

}

unsigned Toroidal2dVertexMeshWithMutableSize::GetNearestNodeIndexIgnorePeriodicity(const ChastePoint<2>& rTestPoint)
{
    unsigned num_nodes = GetNumNodes();
    // Hold the best distance from node to point found so far
    // and the (local) node at which this was recorded
    unsigned best_node_index = 0u;
    double best_node_point_distance = DBL_MAX;

    const c_vector<double, 2>& test_location = rTestPoint.rGetLocation();
    // Now loop through the nodes, calculating the distance and updating best_node_point_distance
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        // Calculate the distance from the chosen point to the current node
        double node_point_distance = norm_2(mNodes[node_index]->rGetLocation() - test_location );
        // Update the "best" distance and node index if necessary
        if (node_point_distance < best_node_point_distance)
        {
            best_node_index = node_index;
            best_node_point_distance = node_point_distance;
        }
    }

    return best_node_index;
}

double Toroidal2dVertexMeshWithMutableSize::GetBoxCoords(const unsigned& rDimension)
{
    assert(rDimension==0 || rDimension==1 || rDimension==2 || rDimension==3);

    double returnVal=0;

    if (rDimension==0)
    {
        returnVal = mBoxLowerX;
    }
    if (rDimension==1)
    {
        returnVal = mBoxUpperX;
    }
    if (rDimension==2)
    {
        returnVal = mBoxLowerY;
    }
    if (rDimension==3)
    {
        returnVal = mBoxUpperY;
    }

    return returnVal;

}

void Toroidal2dVertexMeshWithMutableSize::RefitPeriodicBox()
{
    // // \TODO can't auto find the boxwidth because this toroidal framework uses that as an input to get the position of the boundary.
    // // Bottom left corner:
    // c_vector<double, 2> boxCoord;
    // boxCoord(0) = mBoxUpperX;
    // boxCoord(1) = mBoxUpperY;
    // // Find the node closest to bottom left corner
    // unsigned top_node_index = GetNearestNodeIndexIgnorePeriodicity(boxCoord);
    // c_vector<double, 2> top_node_location;
    // top_node_location = GetNode(top_node_index)->rGetLocation();
    // // Upper Corner
    // boxCoord(0) = top_node_location[0];
    // boxCoord(1) = top_node_location[1] - mHeight;
    // // Find the node closest to bottom left corner
    // unsigned bottom_node_index = GetNearestNodeIndexIgnorePeriodicity(boxCoord);
    //
    // // Top right corner
    // boxCoord(0) = top_node_location[0] - mWidth;
    // boxCoord(1) = top_node_location[1];
    // // Find the node closest to bottom left corner
    // unsigned left_node_index = GetNearestNodeIndexIgnorePeriodicity(boxCoord);
    //
    // std::cout << left_node_index << ", " <<  top_node_index << ", " << bottom_node_index << '\n';
    // std::cout << top_node_location[0] << ", " << top_node_location[1] << '\n';
    // std::cout << GetNode(bottom_node_index)->rGetLocation()[0] << ", " << GetNode(bottom_node_index)->rGetLocation()[1] << '\n';
    // std::cout << GetNode(left_node_index)->rGetLocation()[0] << ", " << GetNode(left_node_index)->rGetLocation()[1] << '\n';
    // std::cout << "/* message */" << '\n';
    // abort();
    //
    // // Set the new box size.
    // SetBoxCoords(0, GetNode(bottom_node_index)->rGetLocation()[0]);
    // SetBoxCoords(1, top_node_location[0]);
    // SetBoxCoords(2, GetNode(left_node_index)->rGetLocation()[1]);
    // SetBoxCoords(3, top_node_location[1]);

    // Initialise minimum values as
    double minX = 100000;
    double minY = 100000;

    unsigned num_nodes = GetNumNodes();
    for (unsigned index=0; index<num_nodes; index++)
    {
        c_vector<double, 2> location;
        location = GetNode(index)->rGetLocation();
        // If the location is lower than the current min, update.
        if( location[0] < minX )
        {
            minX = location[0];
        }
        if( location[1] < minY )
        {
            minY = location[1];
        }
    }
    double oldHeight = GetWidth(1);
    double oldWidth = GetWidth(0);
    SetBoxCoords(0, minX);
    SetBoxCoords(1, minX + oldWidth);
    SetBoxCoords(2, minY);
    SetBoxCoords(3, minY + oldHeight);

}

void Toroidal2dVertexMeshWithMutableSize::SetBoxCoords(const unsigned& rDimension, double newVal)
{
    assert(rDimension==0 || rDimension==1 || rDimension==2 || rDimension==3);

    if (rDimension==0)
    {
        mBoxLowerX = newVal;
        mWidth = mBoxUpperX - mBoxLowerX;
    }
    if (rDimension==1)
    {
        mBoxUpperX = newVal;
        mWidth = mBoxUpperX - mBoxLowerX;
    }
    if (rDimension==2)
    {
        mBoxLowerY = newVal;
        mHeight = mBoxUpperY - mBoxLowerY;
    }
    if (rDimension==3)
    {
        mBoxUpperY = newVal;
        mHeight = mBoxUpperY - mBoxLowerY;
    }

}


std::set<unsigned> Toroidal2dVertexMeshWithMutableSize::GetBoundaryNodes()
{
    // Set to store boundary
    std::set<unsigned> boundaryNodes;

    unsigned num_nodes = this->GetNumNodes();
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        // Get the node
        Node<2>* p_node = this->GetNode(node_index);

        // Get neighbours
        std::set<unsigned> neighbourNodes = GetNeighbouringNodeIndices(node_index);

        // // Iterate over connected nodes
        for(auto neighbourIndex : neighbourNodes)
        {
             // Get dist between nodes
             double nodeY = p_node->rGetLocation()[1];
             double nodeX = p_node->rGetLocation()[0];
             Node<2>* neighbour_p_node = this->GetNode(neighbourIndex);
             double neighbourY = neighbour_p_node->rGetLocation()[1];
             double neighbourX = neighbour_p_node->rGetLocation()[0];

             double width = this->GetWidth(0);
             double height = this->GetWidth(1);
             // If on other side of mesh, mark as boundary node.
             if( fabs(nodeX - neighbourX) > 0.5*width || fabs(nodeY - neighbourY) > 0.5*height)
             {
                 boundaryNodes.insert(node_index);
                 for(auto nabIndex : neighbourNodes)
                 {
                     boundaryNodes.insert(nabIndex);
                 }
                 break;
             }
        }
    }

    return boundaryNodes;
}


void Toroidal2dVertexMeshWithMutableSize::SetReferenceStress(AbstractCellPopulation<2,2>& rCellPopulation, FarhadifarForce<2>* p_force, bool useZero)
{
    if (useZero)
    {
        mReferenceStress(0,0) = 0;
        mReferenceStress(0,1) = 0;
        mReferenceStress(1,0) = 0;
        mReferenceStress(1,1) = 0;
    }
    else
    {
        // Get StressTensor
        c_matrix<double, 2,2> stressTensor2d = GetTissueStressTensor(rCellPopulation, p_force);
        // Set Reference stress
        mReferenceStress(0,0) = stressTensor2d(0,0);
        mReferenceStress(0,1) = stressTensor2d(0,1);
        mReferenceStress(1,0) = stressTensor2d(1,0);
        mReferenceStress(1,1) = stressTensor2d(1,1);
    }
}

void Toroidal2dVertexMeshWithMutableSize::RelaxPeriodicBox(AbstractCellPopulation<2,2>& rCellPopulation, FarhadifarForce<2>* p_force, double stiffness)
{
    // StressTensor
    c_matrix<double, 2,2> stressTensor2d = GetTissueStressTensor(rCellPopulation, p_force);
    // Force in x
    double xForce = (mReferenceStress(0,0)-stressTensor2d(0,0)) * mHeight;
    // Force in y
    double yForce = (mReferenceStress(1,1) - stressTensor2d(1,1)) * mWidth;

    // Displacement of box is force/stiffness
    double deltaX = xForce / stiffness;
    double deltaY = yForce / stiffness;

    // Centroid of box
    double boxCentroidX = (mBoxLowerX + mBoxUpperX)/2;
    double boxCentroidY = (mBoxLowerY + mBoxUpperY)/2;
    // Angles to corners of box. Used to see which side boundary vertices
    // are on
    double anglex0y0 = std::atan2(mBoxLowerY-boxCentroidY,
                                mBoxLowerX-boxCentroidX);
    double anglex0y1 = std::atan2(mBoxUpperY-boxCentroidY,
                                mBoxLowerX-boxCentroidX);
    // double anglex1y1 = std::atan2(mBoxUpperY-boxCentroidY,
    //                             mBoxUpperX-boxCentroidX);
    double anglex1y0 = std::atan2(mBoxLowerY-boxCentroidY,
                                mBoxUpperX-boxCentroidX);

    // Find the boundary nodes
    std::set<unsigned> boundaryNodes = this->GetBoundaryNodes();

    // If it was a boundary node, check where it is and pull it appropriately
    for(auto n_index : boundaryNodes)
    {
        // Get the node
        Node<2>* p_node = this->GetNode(n_index);
        double nodeX = p_node->rGetLocation()[0];
        double nodeY = p_node->rGetLocation()[1];

        // angle to node
        double angleNode = std::atan2(nodeY-boxCentroidY,
                                    nodeX-boxCentroidX);
        // Stretch/compress at the correct boundaries. Note, don't need to stretch upper and RHS because they are deermined by height and width of box, so just stretch the box after.
        // Bottom
        if (anglex0y0 <= angleNode && anglex1y0 >= angleNode)
        {
            p_node->rGetModifiableLocation()[1] += deltaY;
        }
        // LHS
        if ( angleNode <= anglex0y0 || angleNode >= anglex0y1 )
        {
            p_node->rGetModifiableLocation()[0] += deltaX;
        }
    }

    // Reset the size of the box
    SetBoxCoords(0, mBoxLowerX + deltaX);
    SetBoxCoords(1, mBoxUpperX - deltaX);
    SetBoxCoords(2, mBoxLowerY + deltaY);
    SetBoxCoords(3, mBoxUpperY - deltaY);

}

unsigned Toroidal2dVertexMeshWithMutableSize::AddNode(Node<2>* pNewNode)
{
    unsigned node_index = MutableVertexMesh<2,2>::AddNode(pNewNode);

    // If necessary move it to be back onto the torus
    ChastePoint<2> new_node_point = pNewNode->GetPoint();
    SetNode(node_index, new_node_point);

    return node_index;
}

VertexMesh<2, 2>* Toroidal2dVertexMeshWithMutableSize::GetMeshForVtk()
{
    unsigned num_nodes = GetNumNodes();

    std::vector<Node<2>*> temp_nodes(4*num_nodes);
    std::vector<VertexElement<2, 2>*> elements;

    // Create four copies of each node
    for (unsigned index=0; index<num_nodes; index++)
    {
        c_vector<double, 2> location;
        location = GetNode(index)->rGetLocation();

        // Node copy at original location
        Node<2>* p_node = new Node<2>(index, false, location[0], location[1]);
        temp_nodes[index] = p_node;

        // Node copy shifted right
        p_node = new Node<2>(num_nodes + index, false, location[0] + mWidth, location[1]);
        temp_nodes[num_nodes + index] = p_node;

        // Node copy shifted up
        p_node = new Node<2>(2*num_nodes + index, false, location[0], location[1] + mHeight);
        temp_nodes[2*num_nodes + index] = p_node;

        // Node copy shifted right and up
        p_node = new Node<2>(3*num_nodes + index, false, location[0] + mWidth, location[1] + mHeight);
        temp_nodes[3*num_nodes + index] = p_node;
    }

    // Iterate over elements
    for (VertexMesh<2,2>::VertexElementIterator elem_iter = GetElementIteratorBegin();
         elem_iter != GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        unsigned num_nodes_in_elem = elem_iter->GetNumNodes();

        std::vector<Node<2>*> elem_nodes;

        // Compute whether the element straddles either periodic boundary
        bool element_straddles_left_right_boundary = false;
        bool element_straddles_top_bottom_boundary = false;

        const c_vector<double, 2>& r_this_node_location = elem_iter->GetNode(0)->rGetLocation();
        for (unsigned local_index=0; local_index<num_nodes_in_elem; local_index++)
        {
            const c_vector<double, 2>& r_next_node_location = elem_iter->GetNode((local_index+1)%num_nodes_in_elem)->rGetLocation();
            c_vector<double, 2> vector;
            vector = r_next_node_location - r_this_node_location;

            if (fabs(vector[0]) > 0.5*mWidth)
            {
                element_straddles_left_right_boundary = true;
            }
            if (fabs(vector[1]) > 0.5*mHeight)
            {
                element_straddles_top_bottom_boundary = true;
            }
        }

        // Use the above information when duplicating the element
        for (unsigned local_index=0; local_index<num_nodes_in_elem; local_index++)
        {
            unsigned this_node_index = elem_iter->GetNodeGlobalIndex(local_index);

            // If the element straddles the left/right periodic boundary...
            if (element_straddles_left_right_boundary)
            {
                double boxMeanX = (mBoxUpperX + mBoxLowerX)/2;
                // ...and this node is located to the left of the centre of the mesh...
                bool node_is_right_of_centre = (elem_iter->GetNode(local_index)->rGetLocation()[0] > boxMeanX);
                if (!node_is_right_of_centre)
                {
                    // ...then choose the equivalent node to the right
                    this_node_index += num_nodes;
                }
            }

            // If the element straddles the top/bottom periodic boundary...
            if (element_straddles_top_bottom_boundary)
            {
                double boxMeanY = (mBoxUpperY + mBoxLowerY)/2;
                // ...and this node is located below the centre of the mesh...
                bool node_is_above_centre = (elem_iter->GetNode(local_index)->rGetLocation()[1] > boxMeanY);
                if (!node_is_above_centre)
                {
                    // ...then choose the equivalent node above
                    this_node_index += 2*num_nodes;
                }
            }

            elem_nodes.push_back(temp_nodes[this_node_index]);
        }

        VertexElement<2,2>* p_element = new VertexElement<2,2>(elem_index, elem_nodes);
        elements.push_back(p_element);
    }

    // Now delete any nodes from the mesh for VTK that are not contained in any elements
    std::vector<Node<2>*> nodes;
    unsigned count = 0;
    for (unsigned index=0; index<temp_nodes.size(); index++)
    {
        unsigned num_elems_containing_this_node = temp_nodes[index]->rGetContainingElementIndices().size();

        if (num_elems_containing_this_node == 0)
        {
            // Avoid memory leak
            delete temp_nodes[index];
        }
        else
        {
            temp_nodes[index]->SetIndex(count);
            nodes.push_back(temp_nodes[index]);
            count++;
        }
    }

    mpMeshForVtk = new VertexMesh<2,2>(nodes, elements);
    return mpMeshForVtk;
}

void Toroidal2dVertexMeshWithMutableSize::ConstructFromMeshReader(AbstractMeshReader<2,2>& rMeshReader, double width, double height)
{
    assert(rMeshReader.HasNodePermutation() == false);

    // Store numbers of nodes and elements
    unsigned num_nodes = rMeshReader.GetNumNodes();
    unsigned num_elements = rMeshReader.GetNumElements();

    // Reserve memory for nodes
    this->mNodes.reserve(num_nodes);

    rMeshReader.Reset();

    // Add nodes
    std::vector<double> node_data;
    for (unsigned node_idx = 0 ; node_idx < num_nodes ; node_idx++)
    {
        node_data = rMeshReader.GetNextNode();
        node_data.pop_back();
        this->mNodes.push_back(new Node<2>(node_idx, node_data, false));
    }

    rMeshReader.Reset();

    // Reserve memory for elements
    mElements.reserve(rMeshReader.GetNumElements());

    // Add elements
    for (unsigned elem_idx = 0 ; elem_idx < num_elements ; elem_idx++)
    {
        // Get the data for this element
        ElementData element_data = rMeshReader.GetNextElementData();

        // Get the nodes owned by this element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes_in_element = element_data.NodeIndices.size();
        for (unsigned j=0; j<num_nodes_in_element; j++)
        {
            assert(element_data.NodeIndices[j] < this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        // Use nodes and index to construct this element
        VertexElement<2,2>* p_element = new VertexElement<2,2>(elem_idx, nodes);
        mElements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = (unsigned) element_data.AttributeValue;
            p_element->SetAttribute(attribute_value);
        }
    }

    /*
     * Set width and height from function arguments, and validate by checking area is correct
     */
    this->mWidth = width;
    this->mHeight = height;

    double total_surface_area = 0.0;
    for (unsigned elem_idx = 0 ; elem_idx < num_elements ; elem_idx++)
    {
        total_surface_area += this->GetVolumeOfElement(elem_idx);
    }

    if (fabs(mWidth * mHeight - total_surface_area) > 1e-6)
    {
        EXCEPTION("Mesh width and height do not match sheet surface area.");
    }

    // Set default parameter values
    this->mCellRearrangementRatio = 1.5;
    this->mCellRearrangementThreshold = 0.01;
    this->mT2Threshold = 0.001;
    this->mMeshChangesDuringSimulation = true;

    this->mBoxLowerY = 0;
    this->mBoxUpperY = height;
    this->mBoxLowerX = 0;
    this->mBoxUpperX = width;
    this->mpMeshForVtk = nullptr;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Toroidal2dVertexMeshWithMutableSize)
