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
add_test(TestSdkSimulations_NB "/home/chaste/projects/AlexNB/test/TestSdkSimulations_NB")
set_tests_properties(TestSdkSimulations_NB PROPERTIES  LABELS "Continuous_project_AlexNB" PROCESSORS "1" WORKING_DIRECTORY "/home/chaste/src/")
