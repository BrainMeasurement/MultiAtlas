#-----------------------------------------------------------------------------
set(TEMP ${CMAKE_CURRENT_BINARY_DIR}/../Testing/Temporary)
set(INPUT ${CMAKE_CURRENT_SOURCE_DIR}/../Data/Input)
set(BASELINE ${CMAKE_CURRENT_SOURCE_DIR}/../Data/Baseline)

message (STATUS ${CMAKE_CURRENT_BINARY_DIR})
message (STATUS ${CMAKE_CURRENT_SOURCE_DIR})
#-----------------------------------------------------------------------------
# Training module
set(Training_CLP IGR3D_MABMIS_Training)

#-----------------------------------------------------------------------------
add_executable(${Training_CLP}Test ${Training_CLP}Test.cxx)
target_link_libraries(${Training_CLP}Test ${Training_CLP}Lib)
set_target_properties(${Training_CLP}Test PROPERTIES LABELS ${Training_CLP})


configure_file(${INPUT}/TrainingData/TrainingData.xml ${TEMP}/TrainingData.xml)

#-----------------------------------------------------------------------------
set(testname ${Training_CLP}Test)
ExternalData_add_test(${Training_CLP}Data NAME ${testname} COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${Training_CLP}Test>
  ModuleEntryPoint
   -i 5,3,2 -s 1.5
   --trainingXML ${TEMP}/TrainingData.xml --atlasTreeXML ${TEMP}/TrainedAtlas.xml
  )
set_property(TEST ${testname} PROPERTY LABELS ${Training_CLP})

set(testname compareTraingingXMLs)
ExternalData_add_test(${Training_CLP}Data NAME ${testname} COMMAND diff ${INPUT}/AtlasTree/TrainedAtlas.xml ${TEMP}/TrainedAtlas.xml)
set_property(TEST compareTraingingXMLs APPEND PROPERTY DEPENDS ${Training_CLP}Test)

#-----------------------------------------------------------------------------
# Testing module
set(Testing_CLP IGR3D_MABMIS_Testing)

#-----------------------------------------------------------------------------
add_executable(${Testing_CLP}Test ${Testing_CLP}Test.cxx)
target_link_libraries(${Testing_CLP}Test ${Testing_CLP}Lib)
set_target_properties(${Testing_CLP}Test PROPERTIES LABELS ${Testing_CLP})

#-----------------------------------------------------------------------------
configure_file(${INPUT}/AtlasTree/TrainedAtlas.xml ${TEMP}/TrainedAtlas.xml)
set(testname ${Testing_CLP}Test)
ExternalData_add_test(${Testing_CLP}Data NAME ${testname} COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${Testing_CLP}Test>
  --compareRadiusTolerance 3 --compare DATA{${BASELINE}/UA_cbq_000_seg.nii.gz}
            ${TEMP}/UA_cbq_000_seg.nii.gz
  ModuleEntryPoint
   -i 5,3,2 -s 1.5 -o ${TEMP}
   --atlasTreeXML ${TEMP}/TrainedAtlas.xml --imageListXML ${INPUT}/TestData/TestData.xml
  )
set_property(TEST ${testname} PROPERTY LABELS ${Testing_CLP})
