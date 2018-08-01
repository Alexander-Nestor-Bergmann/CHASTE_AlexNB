#ifndef TESTSDKSIMULATIONS_HPP_
#define TESTSDKSIMULATIONS_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "ConstantTargetAreaModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "FakePetscSetup.hpp"
#include "ExtrinsicPullModifier.hpp"
#include "SidekickBoundaryCondition.hpp"
#include "ForceForScenario4.hpp"

static const double M_DT = 0.01;
static const double M_RELAXATION_TIME = 2;
static const double M_EXTENSION_TIME = 3;
static const double M_VIS_TIME_STEP = 1;
static const unsigned M_NUM_CELLS_WIDE = 14;
static const unsigned M_NUM_CELLS_HIGH = 20;

bool myfunction (int i,int j) { return (i<j); }

class TestSdkSimulations : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void xTestAllInOnego() throw (Exception)
    {

        // \TODO Define boundary vertices at beginning, that need to have the extrinsic pull applied, then only update after T1s on those vertices. Looks like epsilon needs to be increased slightly because some of the boundary vells are not being pulled
        // Looks like some rosettes are being formed -- what stops them from being resovled?
        // IS it possible to use cout to print current time?

        // Create box. Find vertices, v1,v2 closest to top and bottom corners. Order all boundary vertices. Mark all vertices between v1 and v2 to have external stretch condition applied.


        // Specify simulation rules
        bool check_internal_intersections = false;

        // Specify mechanical parameter values
        double k = 1.0;
        double lambda_bar = 0.05;
        double gamma_bar = 0.04;
        double heterotypic_line_tension_multiplier = 2.0;
        double supercontractile_line_tension_multiplier = 4.0;

        // Initialise various singletons
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        CellPropertyRegistry::Instance()->Clear();
        CellId::ResetMaxCellId();

        // Generate a vertex mesh
        HoneycombVertexMeshGenerator honeycomb_generator(M_NUM_CELLS_WIDE, M_NUM_CELLS_HIGH);
        MutableVertexMesh<2,2>* p_mesh = honeycomb_generator.GetMesh();
        p_mesh->SetCheckForInternalIntersections(check_internal_intersections);

        // Create some non-proliferating cells
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Bestow cell stripe identities
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->GetCellData()->SetItem("stripe", 1);
        }

        // Create a cell population that associates the cells with the vertex mesh
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputResultsForChasteVisualizer(true);
        cell_population.SetOutputCellRearrangementLocations(false);

        // Create a simulation using the cell population
        OffLatticeSimulation<2> simulation(cell_population);
        simulation.SetOutputDirectory("TestAllInOnego");
        simulation.SetEndTime(M_RELAXATION_TIME);

        simulation.SetDt(M_DT);
        unsigned output_time_step_multiple = (unsigned) (0.1*M_RELAXATION_TIME/M_DT);
        simulation.SetSamplingTimestepMultiple(output_time_step_multiple);

        // Create the appropriate force law(s) for the specified geometry
//        MAKE_PTR(SidekickForce<2>, p_force);
        MAKE_PTR(ForceForScenario4<2>, p_force);
        p_force->SetNumStripes(4);
        p_force->SetAreaElasticityParameter(k);
        p_force->SetPerimeterContractilityParameter(gamma_bar*k);
        p_force->SetHomotypicLineTensionParameter(lambda_bar*pow(k,1.5));
        p_force->SetHeterotypicLineTensionParameter(heterotypic_line_tension_multiplier*lambda_bar*pow(k,1.5));
        p_force->SetSupercontractileLineTensionParameter(supercontractile_line_tension_multiplier*lambda_bar*pow(k,1.5));
        p_force->SetBoundaryLineTensionParameter(lambda_bar*lambda_bar*pow(k,1.5));

        simulation.AddForce(p_force);

        // Pass in a target area modifier (needed, but not used)
        MAKE_PTR(ConstantTargetAreaModifier<2>, p_growth_modifier);
        simulation.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulation.Solve();

        simulation.SetEndTime(M_RELAXATION_TIME + M_EXTENSION_TIME);

        // Bestow cell stripe identities
        for (unsigned i=0; i<simulation.rGetCellPopulation().GetNumRealCells(); i++)
        {
            unsigned row = i/M_NUM_CELLS_WIDE;
            unsigned col = i%M_NUM_CELLS_WIDE;

            CellPtr p_cell = simulation.rGetCellPopulation().GetCellUsingLocationIndex(i);
            if (row%4 == 0)
            {
                if ((col%7 == 0) || (col%7 == 4))      { p_cell->GetCellData()->SetItem("stripe", 1); }
                else if ((col%7 == 1) || (col%7 == 5)) { p_cell->GetCellData()->SetItem("stripe", 2); }
                else if ((col%7 == 2) || (col%7 == 6)) { p_cell->GetCellData()->SetItem("stripe", 3); }
                else                                   { p_cell->GetCellData()->SetItem("stripe", 4); }
            }
            else if (row%4 == 1)
            {
                if ((col%7 == 0) || (col%7 == 3))      { p_cell->GetCellData()->SetItem("stripe", 1); }
                else if ((col%7 == 1) || (col%7 == 4)) { p_cell->GetCellData()->SetItem("stripe", 2); }
                else if (col%7 == 5)                   { p_cell->GetCellData()->SetItem("stripe", 3); }
                else                                   { p_cell->GetCellData()->SetItem("stripe", 4); }
            }
            else if (row%4 == 2)
            {
                if ((col%7 == 0) || (col%7 == 4))      { p_cell->GetCellData()->SetItem("stripe", 1); }
                else if ((col%7 == 1) || (col%7 == 5)) { p_cell->GetCellData()->SetItem("stripe", 2); }
                else if (col%7 == 2)                   { p_cell->GetCellData()->SetItem("stripe", 3); }
                else                                   { p_cell->GetCellData()->SetItem("stripe", 4); }
            }
            else
            {
                if ((col%7 == 0) || (col%7 == 3))      { p_cell->GetCellData()->SetItem("stripe", 1); }
                else if ((col%7 == 1) || (col%7 == 4)) { p_cell->GetCellData()->SetItem("stripe", 2); }
                else if ((col%7 == 2) || (col%7 == 5)) { p_cell->GetCellData()->SetItem("stripe", 3); }
                else                                   { p_cell->GetCellData()->SetItem("stripe", 4); }
            }
        }

        p_force->SetUseCombinedInterfacesForLineTension(true);
        p_force->SetUseDistinctStripeMismatchesForCombinedInterfaces(true);

        // Impose a sliding condition at each boundary
        MAKE_PTR_ARGS(SidekickBoundaryCondition<2>, p_bc, (&(simulation.rGetCellPopulation())));
        simulation.AddCellPopulationBoundaryCondition(p_bc);


        // Get the centroid of the boundary nodes
        /// Get total  umber of nodes
        unsigned num_nodes = simulation.rGetCellPopulation().GetNumNodes();
        // Define variable
        double centroidX = 0;
        double centroidY = 0;
        double counter = 0;
        // Sum the X,Y coords for all points
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            Node<2>* tempNode = simulation.rGetCellPopulation().GetNode(node_index);

            if (tempNode->IsBoundaryNode())
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
            Node<2>* p_node = simulation.rGetCellPopulation().GetNode(node_index);

            // If the node is a boundary node, order and store it.
            if (p_node->IsBoundaryNode())
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
                        Node<2>* tempNode = simulation.rGetCellPopulation().GetNode(boundaryNodes[i]);
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

        ChasteCuboid<2> bounds = simulation.rGetCellPopulation().rGetMesh().CalculateBoundingBox();
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
            Node<2>* tempNode = simulation.rGetCellPopulation().GetNode(boundaryNodes[i]);

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

        // Loop over the boundary nodes and return only the nodes between the
        // upper and lower
        for (int i=0; i < boundaryNodes.size(); i++)
        {
            if( i >= lowerNode && i <= upperNode)
            {
                Node<2>* tempNode = simulation.rGetCellPopulation().GetNode(boundaryNodes[i]);
                cout << tempNode->rGetLocation()[0] << ", ";
            }
            else if( i > upperNode)
            {
                break;
            }
        }
        std::cout << "\n" << '\n';
        for (int i=0; i < boundaryNodes.size(); i++)
        {
            if( i >= lowerNode && i <= upperNode)
            {
                Node<2>* tempNode = simulation.rGetCellPopulation().GetNode(boundaryNodes[i]);
                cout << tempNode->rGetLocation()[1] << ", ";
            }
            else if( i > upperNode)
            {
                break;
            }
        }

        Node<2>* upNode = simulation.rGetCellPopulation().GetNode(boundaryNodes[upperNode]);
        Node<2>* lowNode = simulation.rGetCellPopulation().GetNode(boundaryNodes[lowerNode]);

        cout << upNode->rGetLocation()[0] << ", " << upNode->rGetLocation()[1] << "\n";
        cout << lowNode->rGetLocation()[0] << ", " << lowNode->rGetLocation()[1] << "\n";

        // for (unsigned i=0; i < boundaryNodes.size(); i++)
        // {
        //     Node<2>* tempNode = simulation.rGetCellPopulation().GetNode(boundaryNodes[i]);
        //     cout << tempNode->rGetLocation()[0] << ", ";
        // }
        // std::cout << "\n" << '\n';
        // for (unsigned i=0; i < boundaryNodes.size(); i++)
        // {
        //     Node<2>* tempNode = simulation.rGetCellPopulation().GetNode(boundaryNodes[i]);
        //     cout << tempNode->rGetLocation()[1] << ", ";
        // }


        ///\todo work out whether we need to impose the sliding BC here too (!)
        MAKE_PTR(ExtrinsicPullModifier<2>, p_modifier);
        p_modifier->ApplyExtrinsicPullToAllNodes(false);
        p_modifier->SetSpeed(0.1);
        simulation.AddSimulationModifier(p_modifier);

        simulation.Solve();
    }
};

#endif /* TESTSDKSIMULATIONS_HPP_*/
