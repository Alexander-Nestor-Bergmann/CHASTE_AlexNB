#ifndef TESTTOROIDAL2DVERTEXMESHWITHMUTABLESIZE_HPP_
#define TESTTOROIDAL2DVERTEXMESHWITHMUTABLESIZE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

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
#include "ExtrinsicPullModifierToroidal.hpp"
#include "FarhadifarForce.hpp"
#include "Toroidal2dVertexMeshWithMutableSize.hpp"
#include "ToroidalHoneycombVertexMeshGeneratorMutable.hpp"

static const double M_DT = 0.1;
static const double M_RELAXATION_TIME = 50;
static const double M_VIS_TIME_STEP = 10;
static const unsigned M_NUM_CELLS_WIDE = 4;
static const unsigned M_NUM_CELLS_HIGH = 4;
static const double M_PULL = 0.04;

class TestToroidal2dVertexMeshWithMutableSize : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestPullOnVertices()
    {
        // Initialise various singletons
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        CellPropertyRegistry::Instance()->Clear();
        CellId::ResetMaxCellId();

        // Create mesh
        ToroidalHoneycombVertexMeshGeneratorMutable generator(M_NUM_CELLS_WIDE, M_NUM_CELLS_HIGH);
        Toroidal2dVertexMeshWithMutableSize* p_mesh = generator.GetMutableToroidalMesh();
        p_mesh->SetCheckForInternalIntersections(false);

        // Create some non-proliferating cells
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create a cell population that associates the cells with the vertex mesh
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputResultsForChasteVisualizer(true);
        cell_population.SetOutputCellRearrangementLocations(false);

        // Create a simulation using the cell population
        OffLatticeSimulation<2> simulation(cell_population);
        simulation.SetOutputDirectory("TestToroidalPull");
        simulation.SetEndTime(M_RELAXATION_TIME);

        simulation.SetDt(M_DT);
        unsigned output_time_step_multiple = (unsigned) (M_VIS_TIME_STEP);
        simulation.SetSamplingTimestepMultiple(output_time_step_multiple);

        // Make the Farhadifar force
        MAKE_PTR(FarhadifarForce<2>, p_force);
        // before passing the force to the simulation
        simulation.AddForce(p_force);

        // Pass in a target area modifier (needed, but not used)
        MAKE_PTR(ConstantTargetAreaModifier<2>, p_growth_modifier);
        simulation.AddSimulationModifier(p_growth_modifier);

        // Extrinsic pull on RHS nodes
        MAKE_PTR(ExtrinsicPullModifierToroidal, p_modifier);
        p_modifier->ApplyExtrinsicPullToAllNodes(false);
        p_modifier->SetSpeed(M_PULL);
        simulation.AddSimulationModifier(p_modifier);

        simulation.Solve();
    }
};

#endif /* TESTTOROIDAL2DVERTEXMESHWITHMUTABLESIZE_HPP_*/
