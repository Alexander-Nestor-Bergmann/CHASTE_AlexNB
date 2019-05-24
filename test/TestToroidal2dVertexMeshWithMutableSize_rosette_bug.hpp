#ifndef TESTSDKSIMULATIONSWITHTOROIDALMESH_HPP_
#define TESTSDKSIMULATIONSWITHTOROIDALMESH_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellsGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "ConstantTargetAreaModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "FakePetscSetup.hpp"
#include "ExtrinsicPullModifierToroidal.hpp"
#include "SidekickBoundaryConditionToroidal.hpp"
#include "ForceForScenario4.hpp"
#include "Toroidal2dVertexMeshWithMutableSize.hpp"
#include "ToroidalHoneycombVertexMeshGeneratorMutable.hpp"
#include "BoundaryBoxRelaxationModifier.hpp"
#include "StressTensor.hpp"
#include <Eigen/Dense>
#include "Debug.hpp"

static const double M_DT = 0.5; // 0.1
static const double M_RELAXATION_TIME = 3;
static const double M_EXTENSION_TIME = 1500;
static const double M_VIS_TIME_STEP = 1;
static const double M_PULL = 0.01; // 0.005 good 0.04 max
static const bool M_APPLY_EXTRINSIC_PULL = true;
static const bool M_RELAX_PERIODIC_BOX = true;
static const bool M_APPLY_FLAGPOLE_CONDITION = true;
static const double M_TISSUE_STIFFNESS = 500; // 500 for a 14x20 tissue, 2000 for 28x40
static const unsigned M_NUM_CELLS_WIDE = 14; // 14
static const unsigned M_NUM_CELLS_HIGH = 20; // 20
static const double M_ROSETTE_PROBABILITY = 0.5;

// NOTE Mesh will only flip cells from top to bottom, not left to right. To
// change this turn on the "SetNode" functions for the x-axis.

class TestSdkSimulationsWithToroidalMeshRosetteBug : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestAllInOnego()
    {
        // Specify simulation rules
        bool check_internal_intersections = false;
        bool use_combined_interfaces_for_line_tension = true;
        bool use_distinct_stripe_mismatches_for_combined_interfaces = false;
        std::string output_name("TestRosetteError");

        // Specify mechanical parameter values
        // -0.259,0.172
        double k = 1.0;
        double lambda_bar = 0.05; // Area = 0.69496417 for L,G = 0.05, 0.04
        double gamma_bar = 0.04;
        // double lambda_bar = -0.569;
        // double gamma_bar = 0.145;
        double heterotypic_line_tension_multiplier = 2.0;
        double supercontractile_line_tension_multiplier = 8.0; // 8.0 for scenario 4, 2.0 for normal

        // Initialise various singletons
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        CellPropertyRegistry::Instance()->Clear();
        CellId::ResetMaxCellId();

        // Create mesh
        ToroidalHoneycombVertexMeshGeneratorMutable generator(M_NUM_CELLS_WIDE, M_NUM_CELLS_HIGH, 0.01, 0.001, 0.66584922);
        Toroidal2dVertexMeshWithMutableSize* p_mesh = generator.GetMutableToroidalMesh();
        p_mesh->SetCheckForInternalIntersections(check_internal_intersections);

        // Set the T1 threshold to be very small so there are no exchanges.
        p_mesh->SetCellRearrangementThreshold(.01);
        // p_mesh->SetT2Threshold(0);

        // Enforce rosettes rather than doing T1s
        p_mesh->SetProtorosetteFormationProbability(M_ROSETTE_PROBABILITY);

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
        // simulation.SetOutputDirectory("TestAllInOnego");
        simulation.SetOutputDirectory(output_name);
        simulation.SetEndTime(M_RELAXATION_TIME);

        simulation.SetDt(M_DT);
        unsigned output_time_step_multiple = (unsigned) (M_VIS_TIME_STEP);
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
        // p_force->SetBoundaryLineTensionParameter(lambda_bar*lambda_bar*pow(k,1.5));
        p_force->SetBoundaryLineTensionParameter(lambda_bar*pow(k,1.5));
        p_force->SetUseCombinedInterfacesForLineTension(false);

        simulation.AddForce(p_force); // Can also add this after initial solve.

        // Pass in a target area modifier (needed, but not used)
        MAKE_PTR(ConstantTargetAreaModifier<2>, p_growth_modifier);
        simulation.AddSimulationModifier(p_growth_modifier);

        MARK;
        // Run simulation
        simulation.Solve();
        MARK;

        // Before bestowing stripes, set the reference stress about which to relax as the starting stress.
        p_mesh->SetReferenceStress(cell_population, p_force.get(), false);
        // // StressTensor
        // c_matrix<double, 2,2> stressTensorPre = GetTissueStressTensor(cell_population, p_force.get());
        // // Output the tensor.
        // std::cout << stressTensorPre(0,0) << ", " << stressTensorPre(0,1) << '\n';
        // std::cout << stressTensorPre(1,0) << ", " << stressTensorPre(1,1) << '\n';

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
            // // Make the last stripe dark blue
            // if( col == M_NUM_CELLS_WIDE-1 )
            // {
            //     p_cell->GetCellData()->SetItem("stripe", 1);
            // }
        }

        p_force->SetUseCombinedInterfacesForLineTension(use_combined_interfaces_for_line_tension);
        p_force->SetUseDistinctStripeMismatchesForCombinedInterfaces(use_distinct_stripe_mismatches_for_combined_interfaces);

        // // Set the reference stress about which to relax as the starting stress. Here is with stripes.
        // p_mesh->SetReferenceStress(cell_population);

        // Extrinsic pull
        MAKE_PTR(ExtrinsicPullModifierToroidal, p_modifier);
        p_modifier->ApplyExtrinsicPullToAllNodes(false);
        p_modifier->ApplyExtrinsicPull(M_APPLY_EXTRINSIC_PULL);
        // p_modifier->RelaxPeriodicBox(M_RELAX_PERIODIC_BOX);
        p_modifier->SetSpeed(M_PULL);
        simulation.AddSimulationModifier(p_modifier);

        // Relaxing the boundary
        MAKE_PTR(BoundaryBoxRelaxationModifier, p_boundary_modifier);
        p_boundary_modifier->RelaxPeriodicBox(M_RELAX_PERIODIC_BOX);
        p_boundary_modifier->SetStiffness(M_TISSUE_STIFFNESS);
        p_boundary_modifier->SetForcePointer(p_force);
        simulation.AddSimulationModifier(p_boundary_modifier);

        // Impose a sliding condition at each boundary
        if (M_APPLY_FLAGPOLE_CONDITION)
        {
            MAKE_PTR_ARGS(SidekickBoundaryConditionToroidal, p_bc, (&(simulation.rGetCellPopulation())));
            simulation.AddCellPopulationBoundaryCondition(p_bc);
        }

        // Run the simulation
        simulation.SetEndTime(M_RELAXATION_TIME + M_EXTENSION_TIME);

        MARK;
        simulation.Solve();
        MARK;
    }
};

#endif /* TESTSDKSIMULATIONSWITHTOROIDALMESH_HPP_*/
