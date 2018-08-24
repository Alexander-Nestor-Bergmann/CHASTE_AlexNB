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

static const double M_DT = 0.1;
static const double M_RELAXATION_TIME = 9;
static const double M_EXTENSION_TIME = 300;
static const double M_VIS_TIME_STEP = 10;
static const unsigned M_NUM_CELLS_WIDE = 12;
static const unsigned M_NUM_CELLS_HIGH = 20;

class TestSdkSimulationsWithToroidalMesh : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestAllInOnego() throw (Exception)
    {
        // Specify simulation rules
        bool check_internal_intersections = false;
        bool use_combined_interfaces_for_line_tension = false;
        bool use_distinct_stripe_mismatches_for_combined_interfaces = false;
        std::string output_name("TestToroidalConstPull_Scenario2");

        // Specify mechanical parameter values
        // -0.259,0.172
        double k = 1.0;
        double lambda_bar = 0.05;
        double gamma_bar = 0.04;
        // double lambda_bar = -0.569;
        // double gamma_bar = 0.145;
        double heterotypic_line_tension_multiplier = 2.0;
        double supercontractile_line_tension_multiplier = 2.0;
        // double heterotypic_line_tension_multiplier = 0.1;
        // double supercontractile_line_tension_multiplier = 0.1;

        // Initialise various singletons
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        CellPropertyRegistry::Instance()->Clear();
        CellId::ResetMaxCellId();

        // Create mesh
        ToroidalHoneycombVertexMeshGeneratorMutable generator(M_NUM_CELLS_WIDE, M_NUM_CELLS_HIGH);
        Toroidal2dVertexMeshWithMutableSize* p_mesh = generator.GetMutableToroidalMesh();
        p_mesh->SetCheckForInternalIntersections(check_internal_intersections);

        // Set the T1 threshold to be very small so there are no exchanges.
        p_mesh->SetCellRearrangementThreshold(.01);

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
        p_force->SetUseCombinedInterfacesForLineTension(use_combined_interfaces_for_line_tension);


        simulation.AddForce(p_force); // Can also add this after initial solve.

        // Pass in a target area modifier (needed, but not used)
        MAKE_PTR(ConstantTargetAreaModifier<2>, p_growth_modifier);
        simulation.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulation.Solve();

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
            // Make the last stripe dark blue
            if( col == M_NUM_CELLS_WIDE-1 )
            {
                p_cell->GetCellData()->SetItem("stripe", 1);
            }
        }

        p_force->SetUseCombinedInterfacesForLineTension(use_combined_interfaces_for_line_tension);
        p_force->SetUseDistinctStripeMismatchesForCombinedInterfaces(use_distinct_stripe_mismatches_for_combined_interfaces);

        ///\todo work out whether we need to impose the sliding BC here too (!)
        MAKE_PTR(ExtrinsicPullModifierToroidal, p_modifier);
        p_modifier->ApplyExtrinsicPullToAllNodes(false);
        p_modifier->SetSpeed(0.04);
        simulation.AddSimulationModifier(p_modifier);

        // Impose a sliding condition at each boundary
        MAKE_PTR_ARGS(SidekickBoundaryConditionToroidal, p_bc, (&(simulation.rGetCellPopulation())));
        simulation.AddCellPopulationBoundaryCondition(p_bc);

        simulation.SetEndTime(M_RELAXATION_TIME + M_EXTENSION_TIME);
        simulation.Solve();
    }
};

#endif /* TESTSDKSIMULATIONSWITHTOROIDALMESH_HPP_*/
