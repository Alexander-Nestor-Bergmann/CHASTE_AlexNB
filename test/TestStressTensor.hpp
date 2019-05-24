#ifndef TESTSTRESSTENSOR_HPP_
#define TESTSTRESSTENSOR_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Toroidal2dVertexMeshWithMutableSize.hpp"
#include "Toroidal2dVertexMesh.hpp"
#include "ToroidalHoneycombVertexMeshGeneratorMutable.hpp"
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
#include "SidekickBoundaryCondition.hpp"
#include "ForceForScenario4.hpp"
#include "StressTensor.hpp"
#include <Eigen/Dense>

static const double M_DT = 0.1;
static const double M_RELAXATION_TIME = 20;
static const double M_VIS_TIME_STEP = 1;
static const unsigned M_NUM_CELLS_WIDE = 10;
static const unsigned M_NUM_CELLS_HIGH = 10;
static const double M_PULL = 0.02;

class TestStressTensor : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestTissueStressTensor()
    {
        // Initialise various singletons
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        CellPropertyRegistry::Instance()->Clear();
        CellId::ResetMaxCellId();

        // Create mesh
        ToroidalHoneycombVertexMeshGeneratorMutable generator(M_NUM_CELLS_WIDE, M_NUM_CELLS_HIGH, 0.01, 0.001, 0.65381798);
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
        simulation.SetOutputDirectory("TestStressTensor");
        simulation.SetEndTime(M_RELAXATION_TIME);

        simulation.SetDt(M_DT);
        unsigned output_time_step_multiple = (unsigned) (M_VIS_TIME_STEP);
        std::cout << 0.1*M_RELAXATION_TIME/M_DT << '\n';
        simulation.SetSamplingTimestepMultiple(output_time_step_multiple);

        // Make the Farhadifar force
        MAKE_PTR(FarhadifarForce<2>, p_force);
        // before passing the force to the simulation
        simulation.AddForce(p_force);

        // Pass in a target area modifier (needed, but not used)
        MAKE_PTR(ConstantTargetAreaModifier<2>, p_growth_modifier);
        simulation.AddSimulationModifier(p_growth_modifier);

        // ///\todo work out whether we need to impose the sliding BC here too (!)
        // MAKE_PTR(ExtrinsicPullModifierToroidal, p_modifier);
        // p_modifier->ApplyExtrinsicPullToAllNodes(false);
        // p_modifier->SetSpeed(M_PULL);
        // simulation.AddSimulationModifier(p_modifier);

        // Solve
        simulation.Solve();

        // StressTensor
        // StressTensor<2> stressTensorObject;
        // c_matrix<double, 2,2> stressTensor2d = stressTensorObject.GetTissueStressTensor(cell_population);
        c_matrix<double, 2,2> stressTensor2d = GetTissueStressTensor(cell_population, p_force.get());

        // Convert to Eigen::Matrix type
        Eigen::Matrix2d eigenStressTensor(2,2);

        // Assign values
        // eigenStressTensor << 1, 2, 3, 4;
        eigenStressTensor(0,0) = stressTensor2d(0,0);
        eigenStressTensor(0,1) = stressTensor2d(0,1);
        eigenStressTensor(1,0) = stressTensor2d(1,0);
        eigenStressTensor(1,1) = stressTensor2d(1,1);

        // Output the tensor
        std::cout << eigenStressTensor << '\n';

        // Create the solver
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(eigenStressTensor);

        // Check for success
        if (eigensolver.info() != Eigen::Success)
        {
            abort();
        }

        // Output results
        cout << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << endl;
        cout << "Here's a matrix whose columns are eigenvectors of A \n"
             << "corresponding to these eigenvalues:\n"
             << eigensolver.eigenvectors() << endl;

        cout << stressTensor2d(0,0) << ", " << stressTensor2d(0,1) << ", " << stressTensor2d(1,0) << ", " << stressTensor2d(1,1) << ", ";
        std::cerr << "/* error message */" << '\n';

        // simulation.Solve();
    }
};

#endif /* TestStressTensor_HPP_*/
