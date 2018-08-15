#include "ExtrinsicPullModifier.hpp"

template<unsigned DIM>
ExtrinsicPullModifier<DIM>::ExtrinsicPullModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mApplyExtrinsicPullToAllNodes(true),
      mSpeed(1.0)
{
}

template<unsigned DIM>
ExtrinsicPullModifier<DIM>::~ExtrinsicPullModifier()
{
}

template<unsigned DIM>
void ExtrinsicPullModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    double epsilon = 0.8;

    double dt = SimulationTime::Instance()->GetTimeStep();
    unsigned num_nodes = rCellPopulation.GetNumNodes();
    ChasteCuboid<DIM> bounds = rCellPopulation.rGetMesh().CalculateBoundingBox();
    double x_min = bounds.rGetLowerCorner()[0];
    double x_max = bounds.rGetUpperCorner()[0];

    if (mApplyExtrinsicPullToAllNodes)
    {
        // Pull on all nodes, with a constant strain rate
        double width = x_max - x_min;
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            Node<DIM>* p_node = rCellPopulation.GetNode(node_index);
            double speed = mSpeed * (p_node->rGetLocation()[0] - x_min) / width;

            // Respect the SidekickBoundaryCondition...
            if (p_node->rGetLocation()[0] > x_min + epsilon)
            {
//                if (p_node->rGetLocation()[0] < x_max - epsilon)
//                {
                    p_node->rGetModifiableLocation()[0] += speed*dt;
//                }
            }
        }
    }
    else
    {
        // Get the centroid of the boundary nodes
        // Define variables for centroid position
        double centroidX = 0;
        double centroidY = 0;
        double counter = 0;
        // Sum the X,Y coords for all points
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            Node<DIM>* tempNode = rCellPopulation.GetNode(node_index);

            std::set<unsigned> neighbourNodes = rCellPopulation.GetNeighbouringNodeIndices(node_index);

            bool isBoundaryNode = false;
            int numNeighbours = neighbourNodes.size();
            if( numNeighbours < 3 )
            {
                isBoundaryNode = true;
            }
            else
            {
                // for (int idx=0; idx<numNeighbours; idx++)
                for(auto tempNeighbour : neighbourNodes)
                {
                    // Node<DIM>* tempNeighbour = rCellPopulation.GetNode(neighbourNodes[idx]);
                    std::set<unsigned> tempNeighbourNodes = rCellPopulation.GetNeighbouringNodeIndices(tempNeighbour);

                    if (tempNeighbourNodes.size() < 3)
                    {
                        isBoundaryNode = true;
                    }
                }
            }

            // if (tempNode->IsBoundaryNode()) // GetNeighbouringElementIndices
            // // GetNeighbouringNodeIndices
            if (isBoundaryNode)
            {
                centroidX += tempNode->rGetLocation()[0];
                centroidY += tempNode->rGetLocation()[1];
                counter += 1;
            }
        }
        // Divide by total to get centroid.
        centroidX /= counter;
        centroidY /= counter;

        // Make a vector of node indices:
        std::vector< int > boundaryNodes;
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            // Get the node at this index
            Node<DIM>* p_node = rCellPopulation.GetNode(node_index);

            std::set<unsigned> neighbourNodes = rCellPopulation.GetNeighbouringNodeIndices(node_index);

            bool isBoundaryNode = false;
            int numNeighbours = neighbourNodes.size();
            if( numNeighbours < 3 )
            {
                isBoundaryNode = true;
            }
            else
            {
                // for (int idx=0; idx<numNeighbours; idx++)
                for(auto tempNeighbour : neighbourNodes)
                {
                    // Node<DIM>* tempNeighbour = rCellPopulation.GetNode(neighbourNodes[idx]);
                    std::set<unsigned> tempNeighbourNodes = rCellPopulation.GetNeighbouringNodeIndices(tempNeighbour);

                    if (tempNeighbourNodes.size() < 3)
                    {
                        isBoundaryNode = true;
                    }
                }
            }

            // if (tempNode->IsBoundaryNode()) // GetNeighbouringElementIndices
            // // GetNeighbouringNodeIndices
            if (isBoundaryNode)
            {
                // Get the X and Y coords
                double currentX = p_node->rGetLocation()[0];
                double currentY = p_node->rGetLocation()[1];

                // If this is the first boundary node, just add it to the list.
                if (boundaryNodes.size() == 0)
                {
                    boundaryNodes.push_back(node_index);
                }
                // else, compare it to the current nodes; put it infront of
                // the first node that makes a larger angle, relative to
                // centroid.
                else
                {
                    int i = 0; // counter
                    bool done = false; // Breaks when we have placed the vertex
                    while( not done )
                    {
                        // If we have reached the end of the list, break
                        if (i == boundaryNodes.size())
                        {
                            boundaryNodes.push_back(node_index);
                            done = true;
                        }

                        // Get the x and y coords from the next node in
                        // boundaryNodes to compare against.
                        Node<DIM>* tempNode = rCellPopulation.GetNode(boundaryNodes[i]);
                        // Get the coords of the node at the current index in
                        // the ordered list of boundary nodes.
                        double previousX = tempNode->rGetLocation()[0];
                        double previousY = tempNode->rGetLocation()[1];

                        // Calculate the angle that the two nodes make relative
                        // to the centroid
                        double prevAngle;
                        double dY = previousY - centroidY;
                        double dX = previousX - centroidX;
                        prevAngle = std::atan2( dY, dX );
                        double currentAngle;
                        dX = currentX - centroidX;
                        dY = currentY - centroidY;
                        currentAngle = std::atan2(dY, dX);
                        // If the angle of new node is smaller, it comes first
                        // so append it at this location.
                        if ( currentAngle < prevAngle && not done )
                        {
                            boundaryNodes.insert(boundaryNodes.begin() + i, node_index);
                            done = true;
                        }

                        // update count
                        i += 1;
                    }
                }
            }
        }

        // Get the top and bottom rightmost corners of bounding box.
        double lowerX = bounds.rGetUpperCorner()[0];
        double lowerY = bounds.rGetLowerCorner()[1];
        double upperX = bounds.rGetUpperCorner()[0];
        double upperY = bounds.rGetUpperCorner()[1];

        // Make variables to hold indices of nodes closest to bounding box.
        int upperNode = 0;
        int lowerNode = 0;
        // Euclid distane variables
        double dist2;
        double upperDist = 10000000.; // Initialise as huge
        double lowerDist = 10000000.; // Initialise as huge
        double tempDist;
        double xDiff;
        double yDiff;
        for (unsigned i=0; i < boundaryNodes.size(); i++)
        {
            // Get the node pointer
            Node<DIM>* tempNode = rCellPopulation.GetNode(boundaryNodes[i]);

            // calc dist to upper corner
            xDiff = upperX - tempNode->rGetLocation()[0];
            yDiff = upperY - tempNode->rGetLocation()[1];
            dist2 = xDiff*xDiff + yDiff*yDiff;
            tempDist = std::sqrt(dist2);
            // If smaller, mark this as the closest node.
            if( tempDist < upperDist )
            {
                upperDist = tempDist;
                upperNode = i;
            }

            // calc dist to lower corner
            xDiff = lowerX - tempNode->rGetLocation()[0];
            yDiff = lowerY - tempNode->rGetLocation()[1];
            dist2 = xDiff*xDiff + yDiff*yDiff;
            tempDist = std::sqrt(dist2);
            // If smaller, mark this as the closest node.
            if( tempDist < lowerDist )
            {
                lowerDist = tempDist;
                lowerNode = i;
            }
        }

        std::cout << "\n" << '\n';
        std::cout << "Boundary\n" << '\n';
        for (unsigned i=0; i < boundaryNodes.size(); i++)
        {
            Node<DIM>* tempNode = rCellPopulation.GetNode(boundaryNodes[i]);
            cout << tempNode->rGetLocation()[0] << ", ";
        }
        std::cout << "\n" << '\n';
        for (unsigned i=0; i < boundaryNodes.size(); i++)
        {
            Node<DIM>* tempNode = rCellPopulation.GetNode(boundaryNodes[i]);
            cout << tempNode->rGetLocation()[1] << ", ";
        }
        std::cout << "\n" << '\n';
        std::cout << "PULL\n" << '\n';

        // Iterate over the rightmost nodes and pull them.
        for (unsigned i=0; i < boundaryNodes.size(); i++)
        {
            // If it lies between the top and bottom:
            if( i >= lowerNode && i <= upperNode)
            {
                Node<DIM>* p_node = rCellPopulation.GetNode(boundaryNodes[i]);
                p_node->rGetModifiableLocation()[0] += mSpeed*dt;

                cout << p_node->rGetLocation()[0] << ", " << p_node->rGetLocation()[1] << ", ";

            }
            else if( i > upperNode)
            {
                break;
            }
        }
        // std::cout << "\n" << '\n';


        // // Pull on the right-most nodes only, with a constant speed
        // for (unsigned node_index=0; node_index<num_nodes; node_index++)
        // {
        //     Node<DIM>* p_node = rCellPopulation.GetNode(node_index);
        //     if (fabs(p_node->rGetLocation()[0] - x_max) < 0.1)
        //     {
        //         p_node->rGetModifiableLocation()[0] += mSpeed*dt;
        //     }
        // }




    }
}

template<unsigned DIM>
void ExtrinsicPullModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
}

template<unsigned DIM>
void ExtrinsicPullModifier<DIM>::ApplyExtrinsicPullToAllNodes(bool applyExtrinsicPullToAllNodes)
{
    mApplyExtrinsicPullToAllNodes = applyExtrinsicPullToAllNodes;
}

template<unsigned DIM>
void ExtrinsicPullModifier<DIM>::SetSpeed(double speed)
{
    mSpeed = speed;
}

template<unsigned DIM>
void ExtrinsicPullModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ApplyExtrinsicPullToAllNodes>" << mApplyExtrinsicPullToAllNodes << "</ApplyExtrinsicPullToAllNodes>\n";
    *rParamsFile << "\t\t\t<Speed>" << mSpeed << "</Speed>\n";

    // Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ExtrinsicPullModifier<1>;
template class ExtrinsicPullModifier<2>;
template class ExtrinsicPullModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ExtrinsicPullModifier)
