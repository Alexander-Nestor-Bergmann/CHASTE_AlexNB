# CMake generated Testfile for 
# Source directory: /home/chaste/src/projects/AlexNB/test
# Build directory: /home/chaste/projects/AlexNB/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(TestHello "/home/chaste/projects/AlexNB/test/TestHello")
set_tests_properties(TestHello PROPERTIES  LABELS "Continuous_project_AlexNB" PROCESSORS "1" WORKING_DIRECTORY "/home/chaste/src/")
add_test(TestSdkSimulations "/home/chaste/projects/AlexNB/test/TestSdkSimulations")
set_tests_properties(TestSdkSimulations PROPERTIES  LABELS "Continuous_project_AlexNB" PROCESSORS "1" WORKING_DIRECTORY "/home/chaste/src/")
add_test(TestToroidal2dVertexMeshWithMutableSize "/home/chaste/projects/AlexNB/test/TestToroidal2dVertexMeshWithMutableSize")
set_tests_properties(TestToroidal2dVertexMeshWithMutableSize PROPERTIES  LABELS "Continuous_project_AlexNB" PROCESSORS "1" WORKING_DIRECTORY "/home/chaste/src/")
add_test(TestSdkSimulationsWithToroidalMesh "/home/chaste/projects/AlexNB/test/TestSdkSimulationsWithToroidalMesh")
set_tests_properties(TestSdkSimulationsWithToroidalMesh PROPERTIES  LABELS "Continuous_project_AlexNB" PROCESSORS "1" WORKING_DIRECTORY "/home/chaste/src/")
add_test(TestStressTensor "/home/chaste/projects/AlexNB/test/TestStressTensor")
set_tests_properties(TestStressTensor PROPERTIES  LABELS "Continuous_project_AlexNB" PROCESSORS "1" WORKING_DIRECTORY "/home/chaste/src/")
add_test(TestToroidal2dVertexMeshWithMutableSize_rosette_bug "/home/chaste/projects/AlexNB/test/TestToroidal2dVertexMeshWithMutableSize_rosette_bug")
set_tests_properties(TestToroidal2dVertexMeshWithMutableSize_rosette_bug PROPERTIES  LABELS "Continuous_project_AlexNB" PROCESSORS "1" WORKING_DIRECTORY "/home/chaste/src/")
