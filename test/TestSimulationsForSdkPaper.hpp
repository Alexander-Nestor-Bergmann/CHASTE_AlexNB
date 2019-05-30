#ifndef TESTSIMULATIONSFORSDKPAPER_HPP_
#define TESTSIMULATIONSFORSDKPAPER_HPP_

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

static const double M_DT = 0.5; // 0.1
static const double M_RELAXATION_TIME = 3;
static const double M_EXTENSION_TIME = 1500;
static const double M_VIS_TIME_STEP = 10;
static const double M_PULL = 0.01; // 0.005 good 0.04 max
static const double M_TISSUE_STIFFNESS = 500; // 500 for a 14x20 tissue, 2000 for 28x40
static const unsigned M_NUM_CELLS_WIDE = 14; // 14
static const unsigned M_NUM_CELLS_HIGH = 20; // 20

// NOTE Mesh will only flip cells from top to bottom, not left to right. To
// change this turn on the "SetNode" functions for the x-axis.

// Before running this test suite, edit MutableVertexMesh::HandleHighOrderJunctions() to
// call PerformRosetteRankIncrease() then RemoveDeletedNodes() regardless of the rosette rank
// of each node. We have checked and this works fine for us.

class TestSimulationsForSdkPaper : public AbstractCellBasedWithTimingsTestSuite
{
public:

    // Tetley et al "scenario 4" with extrinsic pull
    void xTestSimulation1()
    {
        // Specify simulation rules
        bool use_combined_interfaces_for_line_tension = true;
        bool use_distinct_stripe_mismatches_for_combined_interfaces = false;
        std::string output_name("TestSimulation1");

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
        p_mesh->SetCheckForInternalIntersections(false);

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
        simulation.SetOutputDirectory(output_name);
        simulation.SetEndTime(M_RELAXATION_TIME);

        simulation.SetDt(M_DT);
        unsigned output_time_step_multiple = (unsigned) (M_VIS_TIME_STEP);
        simulation.SetSamplingTimestepMultiple(output_time_step_multiple);

        // Create the appropriate force law(s) for the specified geometry
        MAKE_PTR(ForceForScenario4<2>, p_force);
        p_force->SetNumStripes(4);
        p_force->SetAreaElasticityParameter(k);
        p_force->SetPerimeterContractilityParameter(gamma_bar*k);
        p_force->SetHomotypicLineTensionParameter(lambda_bar*pow(k,1.5));
        p_force->SetHeterotypicLineTensionParameter(heterotypic_line_tension_multiplier*lambda_bar*pow(k,1.5));
        p_force->SetSupercontractileLineTensionParameter(supercontractile_line_tension_multiplier*lambda_bar*pow(k,1.5));
        p_force->SetBoundaryLineTensionParameter(lambda_bar*pow(k,1.5));
        p_force->SetUseCombinedInterfacesForLineTension(false);

        simulation.AddForce(p_force); // Can also add this after initial solve.

        // Pass in a target area modifier (needed, but not used)
        MAKE_PTR(ConstantTargetAreaModifier<2>, p_growth_modifier);
        simulation.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulation.Solve();

        // Before bestowing stripes, set the reference stress about which to relax as the starting stress.
        p_mesh->SetReferenceStress(cell_population, p_force.get(), false);

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

        p_force->SetUseCombinedInterfacesForLineTension(use_combined_interfaces_for_line_tension);
        p_force->SetUseDistinctStripeMismatchesForCombinedInterfaces(use_distinct_stripe_mismatches_for_combined_interfaces);

        // Extrinsic pull
        MAKE_PTR(ExtrinsicPullModifierToroidal, p_modifier);
        p_modifier->ApplyExtrinsicPullToAllNodes(false);
        p_modifier->ApplyExtrinsicPull(true);
        p_modifier->SetSpeed(M_PULL);
        simulation.AddSimulationModifier(p_modifier);

        // Relaxing the boundary
        MAKE_PTR(BoundaryBoxRelaxationModifier, p_boundary_modifier);
        p_boundary_modifier->SetStiffness(M_TISSUE_STIFFNESS);
        p_boundary_modifier->SetForcePointer(p_force);
        simulation.AddSimulationModifier(p_boundary_modifier);

        // Impose a sliding condition at each boundary
        MAKE_PTR_ARGS(SidekickBoundaryConditionToroidal, p_bc, (&(simulation.rGetCellPopulation())));
        simulation.AddCellPopulationBoundaryCondition(p_bc);

        // Run the simulation
        simulation.SetEndTime(M_RELAXATION_TIME + M_EXTENSION_TIME);
        simulation.Solve();
    }

    // Same simulation as 1, but prevent T1s
    void xTestSimulation2()
    {
        // Specify simulation rules
        bool use_combined_interfaces_for_line_tension = true;
        bool use_distinct_stripe_mismatches_for_combined_interfaces = false;
        std::string output_name("TestSimulation2");

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
        p_mesh->SetCheckForInternalIntersections(false);

        // Enforce rosettes rather than doing T1s
        p_mesh->SetProtorosetteFormationProbability(1.0);

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
        simulation.SetOutputDirectory(output_name);
        simulation.SetEndTime(M_RELAXATION_TIME);

        simulation.SetDt(M_DT);
        unsigned output_time_step_multiple = (unsigned) (M_VIS_TIME_STEP);
        simulation.SetSamplingTimestepMultiple(output_time_step_multiple);

        // Create the appropriate force law(s) for the specified geometry
        MAKE_PTR(ForceForScenario4<2>, p_force);
        p_force->SetNumStripes(4);
        p_force->SetAreaElasticityParameter(k);
        p_force->SetPerimeterContractilityParameter(gamma_bar*k);
        p_force->SetHomotypicLineTensionParameter(lambda_bar*pow(k,1.5));
        p_force->SetHeterotypicLineTensionParameter(heterotypic_line_tension_multiplier*lambda_bar*pow(k,1.5));
        p_force->SetSupercontractileLineTensionParameter(supercontractile_line_tension_multiplier*lambda_bar*pow(k,1.5));
        p_force->SetBoundaryLineTensionParameter(lambda_bar*pow(k,1.5));
        p_force->SetUseCombinedInterfacesForLineTension(false);

        simulation.AddForce(p_force); // Can also add this after initial solve.

        // Pass in a target area modifier (needed, but not used)
        MAKE_PTR(ConstantTargetAreaModifier<2>, p_growth_modifier);
        simulation.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulation.Solve();

        // Before bestowing stripes, set the reference stress about which to relax as the starting stress.
        p_mesh->SetReferenceStress(cell_population, p_force.get(), false);

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

        p_force->SetUseCombinedInterfacesForLineTension(use_combined_interfaces_for_line_tension);
        p_force->SetUseDistinctStripeMismatchesForCombinedInterfaces(use_distinct_stripe_mismatches_for_combined_interfaces);

        // Extrinsic pull
        MAKE_PTR(ExtrinsicPullModifierToroidal, p_modifier);
        p_modifier->ApplyExtrinsicPullToAllNodes(false);
        p_modifier->ApplyExtrinsicPull(true);
        p_modifier->SetSpeed(M_PULL);
        simulation.AddSimulationModifier(p_modifier);

        // Relaxing the boundary
        MAKE_PTR(BoundaryBoxRelaxationModifier, p_boundary_modifier);
        p_boundary_modifier->SetStiffness(M_TISSUE_STIFFNESS);
        p_boundary_modifier->SetForcePointer(p_force);
        simulation.AddSimulationModifier(p_boundary_modifier);

        // Impose a sliding condition at each boundary
        MAKE_PTR_ARGS(SidekickBoundaryConditionToroidal, p_bc, (&(simulation.rGetCellPopulation())));
        simulation.AddCellPopulationBoundaryCondition(p_bc);

        // Run the simulation
        simulation.SetEndTime(M_RELAXATION_TIME + M_EXTENSION_TIME);
        simulation.Solve();
    }

    // Same as simulation 2, but with 50% chance of rosettes instead of T1
    void xTestSimulation3()
    {
        // Specify simulation rules
        bool use_combined_interfaces_for_line_tension = true;
        bool use_distinct_stripe_mismatches_for_combined_interfaces = false;
        std::string output_name("TestSimulation3");

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
        p_mesh->SetCheckForInternalIntersections(false);

        // Enforce rosettes rather than T1s with some probability
        p_mesh->SetProtorosetteFormationProbability(0.5);

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
        simulation.SetOutputDirectory(output_name);
        simulation.SetEndTime(M_RELAXATION_TIME);

        simulation.SetDt(M_DT);
        unsigned output_time_step_multiple = (unsigned) (M_VIS_TIME_STEP);
        simulation.SetSamplingTimestepMultiple(output_time_step_multiple);

        // Create the appropriate force law(s) for the specified geometry
        MAKE_PTR(ForceForScenario4<2>, p_force);
        p_force->SetNumStripes(4);
        p_force->SetAreaElasticityParameter(k);
        p_force->SetPerimeterContractilityParameter(gamma_bar*k);
        p_force->SetHomotypicLineTensionParameter(lambda_bar*pow(k,1.5));
        p_force->SetHeterotypicLineTensionParameter(heterotypic_line_tension_multiplier*lambda_bar*pow(k,1.5));
        p_force->SetSupercontractileLineTensionParameter(supercontractile_line_tension_multiplier*lambda_bar*pow(k,1.5));
        p_force->SetBoundaryLineTensionParameter(lambda_bar*pow(k,1.5));
        p_force->SetUseCombinedInterfacesForLineTension(false);

        simulation.AddForce(p_force); // Can also add this after initial solve.

        // Pass in a target area modifier (needed, but not used)
        MAKE_PTR(ConstantTargetAreaModifier<2>, p_growth_modifier);
        simulation.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulation.Solve();

        // Before bestowing stripes, set the reference stress about which to relax as the starting stress.
        p_mesh->SetReferenceStress(cell_population, p_force.get(), false);

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

        p_force->SetUseCombinedInterfacesForLineTension(use_combined_interfaces_for_line_tension);
        p_force->SetUseDistinctStripeMismatchesForCombinedInterfaces(use_distinct_stripe_mismatches_for_combined_interfaces);

        // Extrinsic pull
        MAKE_PTR(ExtrinsicPullModifierToroidal, p_modifier);
        p_modifier->ApplyExtrinsicPullToAllNodes(false);
        p_modifier->ApplyExtrinsicPull(true);
        p_modifier->SetSpeed(M_PULL);
        simulation.AddSimulationModifier(p_modifier);

        // Relaxing the boundary
        MAKE_PTR(BoundaryBoxRelaxationModifier, p_boundary_modifier);
        p_boundary_modifier->SetStiffness(M_TISSUE_STIFFNESS);
        p_boundary_modifier->SetForcePointer(p_force);
        simulation.AddSimulationModifier(p_boundary_modifier);

        // Impose a sliding condition at each boundary
        MAKE_PTR_ARGS(SidekickBoundaryConditionToroidal, p_bc, (&(simulation.rGetCellPopulation())));
        simulation.AddCellPopulationBoundaryCondition(p_bc);

        // Run the simulation
        simulation.SetEndTime(M_RELAXATION_TIME + M_EXTENSION_TIME);
        simulation.Solve();
    }

    ///\todo Higher order rosettes are more stable... note that MutableVertexMesh
    // already has two distinct member variables for setting the probability per
    // unit time of resolving a 'protorosette' (4-way junction) and a 'rosette'
    // (>4 way junction). This may be sufficient for our purposes; if not, we can
    // edit MutableVertexMesh::CheckForRosettes() as we see fit.
    void TestSimulation4()
    {
        // Specify simulation rules
        bool use_combined_interfaces_for_line_tension = true;
        bool use_distinct_stripe_mismatches_for_combined_interfaces = false;
        std::string output_name("TestSimulation4");

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
        p_mesh->SetCheckForInternalIntersections(false);

        // Enforce rosettes with some probability
        p_mesh->SetProtorosetteFormationProbability(0.5);
        p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(0.0);
        p_mesh->SetRosetteResolutionProbabilityPerTimestep(0.0);

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
        simulation.SetOutputDirectory(output_name);
        simulation.SetEndTime(M_RELAXATION_TIME);

        simulation.SetDt(M_DT);
        unsigned output_time_step_multiple = (unsigned) (M_VIS_TIME_STEP);
        simulation.SetSamplingTimestepMultiple(output_time_step_multiple);

        // Create the appropriate force law(s) for the specified geometry
        MAKE_PTR(ForceForScenario4<2>, p_force);
        p_force->SetNumStripes(4);
        p_force->SetAreaElasticityParameter(k);
        p_force->SetPerimeterContractilityParameter(gamma_bar*k);
        p_force->SetHomotypicLineTensionParameter(lambda_bar*pow(k,1.5));
        p_force->SetHeterotypicLineTensionParameter(heterotypic_line_tension_multiplier*lambda_bar*pow(k,1.5));
        p_force->SetSupercontractileLineTensionParameter(supercontractile_line_tension_multiplier*lambda_bar*pow(k,1.5));
        p_force->SetBoundaryLineTensionParameter(lambda_bar*pow(k,1.5));
        p_force->SetUseCombinedInterfacesForLineTension(false);

        simulation.AddForce(p_force); // Can also add this after initial solve.

        // Pass in a target area modifier (needed, but not used)
        MAKE_PTR(ConstantTargetAreaModifier<2>, p_growth_modifier);
        simulation.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulation.Solve();

        // Before bestowing stripes, set the reference stress about which to relax as the starting stress.
        p_mesh->SetReferenceStress(cell_population, p_force.get(), false);

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

        p_force->SetUseCombinedInterfacesForLineTension(use_combined_interfaces_for_line_tension);
        p_force->SetUseDistinctStripeMismatchesForCombinedInterfaces(use_distinct_stripe_mismatches_for_combined_interfaces);

        // Extrinsic pull
        MAKE_PTR(ExtrinsicPullModifierToroidal, p_modifier);
        p_modifier->ApplyExtrinsicPullToAllNodes(false);
        p_modifier->ApplyExtrinsicPull(true);
        p_modifier->SetSpeed(M_PULL);
        simulation.AddSimulationModifier(p_modifier);

        // Relaxing the boundary
        MAKE_PTR(BoundaryBoxRelaxationModifier, p_boundary_modifier);
        p_boundary_modifier->SetStiffness(M_TISSUE_STIFFNESS);
        p_boundary_modifier->SetForcePointer(p_force);
        simulation.AddSimulationModifier(p_boundary_modifier);

        // Impose a sliding condition at each boundary
        MAKE_PTR_ARGS(SidekickBoundaryConditionToroidal, p_bc, (&(simulation.rGetCellPopulation())));
        simulation.AddCellPopulationBoundaryCondition(p_bc);

        // Run the simulation
        simulation.SetEndTime(M_RELAXATION_TIME + M_EXTENSION_TIME);
        simulation.Solve();
    }

    ///\todo Higher order rosettes are more stable... note that MutableVertexMesh
    // already has two distinct member variables for setting the probability per
    // unit time of resolving a 'protorosette' (4-way junction) and a 'rosette'
    // (>4 way junction). This may be sufficient for our purposes; if not, we can
    // edit MutableVertexMesh::CheckForRosettes() as we see fit.
    void xTestSimulation5()
    {
        // Specify simulation rules
        bool use_combined_interfaces_for_line_tension = true;
        bool use_distinct_stripe_mismatches_for_combined_interfaces = false;
        std::string output_name("TestSimulation5");

        // Specify mechanical parameter values
        // -0.259,0.172
        double k = 1.0;
        // double lambda_bar = 0.05; // Area = 0.69496417 for L,G = 0.05, 0.04
        // double gamma_bar = 0.04;
        double lambda_bar = 0.05; // -0.569;
        double gamma_bar = 0.08; // 0.145;
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
        p_mesh->SetCheckForInternalIntersections(false);

        // Enforce rosettes with some probability
        p_mesh->SetProtorosetteFormationProbability(0.5);
        p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(0.0);
        p_mesh->SetRosetteResolutionProbabilityPerTimestep(0.0);

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
        simulation.SetOutputDirectory(output_name);
        simulation.SetEndTime(M_RELAXATION_TIME);

        simulation.SetDt(M_DT);
        unsigned output_time_step_multiple = (unsigned) (M_VIS_TIME_STEP);
        simulation.SetSamplingTimestepMultiple(output_time_step_multiple);

        // Create the appropriate force law(s) for the specified geometry
        MAKE_PTR(ForceForScenario4<2>, p_force);
        p_force->SetNumStripes(4);
        p_force->SetAreaElasticityParameter(k);
        p_force->SetPerimeterContractilityParameter(gamma_bar*k);
        p_force->SetHomotypicLineTensionParameter(lambda_bar*pow(k,1.5));
        p_force->SetHeterotypicLineTensionParameter(heterotypic_line_tension_multiplier*lambda_bar*pow(k,1.5));
        p_force->SetSupercontractileLineTensionParameter(supercontractile_line_tension_multiplier*lambda_bar*pow(k,1.5));
        p_force->SetBoundaryLineTensionParameter(lambda_bar*pow(k,1.5));
        p_force->SetUseCombinedInterfacesForLineTension(false);

        simulation.AddForce(p_force); // Can also add this after initial solve.

        // Pass in a target area modifier (needed, but not used)
        MAKE_PTR(ConstantTargetAreaModifier<2>, p_growth_modifier);
        simulation.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulation.Solve();

        // Before bestowing stripes, set the reference stress about which to relax as the starting stress.
        p_mesh->SetReferenceStress(cell_population, p_force.get(), false);

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

        p_force->SetUseCombinedInterfacesForLineTension(use_combined_interfaces_for_line_tension);
        p_force->SetUseDistinctStripeMismatchesForCombinedInterfaces(use_distinct_stripe_mismatches_for_combined_interfaces);

        // Extrinsic pull
        MAKE_PTR(ExtrinsicPullModifierToroidal, p_modifier);
        p_modifier->ApplyExtrinsicPullToAllNodes(false);
        p_modifier->ApplyExtrinsicPull(true);
        p_modifier->SetSpeed(M_PULL);
        simulation.AddSimulationModifier(p_modifier);

        // Relaxing the boundary
        MAKE_PTR(BoundaryBoxRelaxationModifier, p_boundary_modifier);
        p_boundary_modifier->SetStiffness(M_TISSUE_STIFFNESS);
        p_boundary_modifier->SetForcePointer(p_force);
        simulation.AddSimulationModifier(p_boundary_modifier);

        // Impose a sliding condition at each boundary
        MAKE_PTR_ARGS(SidekickBoundaryConditionToroidal, p_bc, (&(simulation.rGetCellPopulation())));
        simulation.AddCellPopulationBoundaryCondition(p_bc);

        // Run the simulation
        simulation.SetEndTime(M_RELAXATION_TIME + M_EXTENSION_TIME);
        simulation.Solve();
    }
};

#endif /* TESTSIMULATIONSFORSDKPAPER_HPP_*/
