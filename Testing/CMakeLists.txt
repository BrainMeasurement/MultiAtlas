#-----------------------------------------------------------------------------
set(TEMP ${CMAKE_CURRENT_BINARY_DIR}/../Testing/Temporary)
set(INPUT ${CMAKE_CURRENT_SOURCE_DIR}/../Data/Input)
set(BASELINE ${CMAKE_CURRENT_SOURCE_DIR}/../Data/Baseline)

message (STATUS "MABMIS tests CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message (STATUS "MABMIS tests CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")


macro(copy_files_recursively srcDir destDir)
    message(STATUS "Configuring directory ${destDir}")
    make_directory(${destDir})

    file(GLOB templateFiles RELATIVE ${srcDir} ${srcDir}/*)
    foreach(templateFile ${templateFiles})
        set(srcTemplatePath ${srcDir}/${templateFile})
        if(NOT IS_DIRECTORY ${srcTemplatePath})
            message(STATUS "Configuring file ${templateFile}")
            configure_file(
                    ${srcTemplatePath}
                    ${destDir}/${templateFile}
                    COPYONLY)
        endif(NOT IS_DIRECTORY ${srcTemplatePath})
    endforeach(templateFile)
endmacro(copy_files_recursively)


#-----------------------------------------------------------------------------
# Training module
set(Training_CLP IGR3D_MABMIS_Training)

#-----------------------------------------------------------------------------
add_executable(${Training_CLP}Test Test.cxx)
target_link_libraries(${Training_CLP}Test ${Training_CLP}Lib ${SlicerExecutionModel_EXTRA_EXECUTABLE_TARGET_LIBRARIES})


configure_file(${INPUT}/TrainingData/TrainingData.xml ${TEMP}/TrainingData.xml COPYONLY)
copy_files_recursively(${INPUT}/TrainingData/ ${TEMP}/)

#-----------------------------------------------------------------------------
set(testname ${Training_CLP}Test)
ExternalData_add_test(${Training_CLP}Data NAME ${testname} COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${Training_CLP}Test>
  ModuleEntryPoint
   -i 5,3,2 -s 1.5
   --trainingXML ${TEMP}/TrainingData.xml --atlasTreeXML ${TEMP}/TrainedAtlas.xml
  )
set_property(TEST ${testname} PROPERTY LABELS MABMIS)

set(testname ${Training_CLP}_compareTraingingXMLs)
if (WIN32)
	ExternalData_add_test(${Training_CLP}Data NAME ${testname} COMMAND diff ${INPUT}/AtlasTree/TrainedAtlas-win.xml ${TEMP}/TrainedAtlas.xml)
else()
	ExternalData_add_test(${Training_CLP}Data NAME ${testname} COMMAND diff ${INPUT}/AtlasTree/TrainedAtlas.xml ${TEMP}/TrainedAtlas.xml)
endif()
set_property(TEST ${testname} PROPERTY LABELS MABMIS)
set_property(TEST ${testname} APPEND PROPERTY DEPENDS ${Training_CLP}Test)

#-----------------------------------------------------------------------------
# Testing module
set(Testing_CLP IGR3D_MABMIS_Testing)

#-----------------------------------------------------------------------------
add_executable(${Testing_CLP}Test Test.cxx)
target_link_libraries(${Testing_CLP}Test ${Testing_CLP}Lib ${SlicerExecutionModel_EXTRA_EXECUTABLE_TARGET_LIBRARIES})

#-----------------------------------------------------------------------------
configure_file(${INPUT}/AtlasTree/TrainedAtlas.xml ${TEMP}/TrainedAtlas.xml)
copy_files_recursively(${INPUT}/AtlasTree/TrainedAtlas ${TEMP}/TrainedAtlas/)

copy_files_recursively(${INPUT}/TestData ${TEMP}/TestData/)

set(testname ${Testing_CLP}Test)
ExternalData_add_test(${Testing_CLP}Data NAME ${testname} COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${Testing_CLP}Test>
  --compareRadiusTolerance 3 --compare DATA{${BASELINE}/UA_cbq_000_seg.nii.gz}
            ${TEMP}/TestData/UA_cbq_000_seg.nii.gz
  ModuleEntryPoint
   -i 5,3,2 -s 1.5 -o ${TEMP}
   --atlasTreeXML ${TEMP}/TrainedAtlas.xml --imageListXML ${TEMP}/TestData/TestData.xml
  )
set_property(TEST ${testname} PROPERTY LABELS MABMIS)


#-----------------------------------------------------------------------------
# Direct invoke module
set(Testing_CLP IGR3D_MABMIS_Direct_Invoke)

add_executable(${Testing_CLP}Test Test.cxx)
target_link_libraries(${Testing_CLP}Test ${Testing_CLP}Lib ${SlicerExecutionModel_EXTRA_EXECUTABLE_TARGET_LIBRARIES})

# Test creating the atlas
set(testname ${Testing_CLP}TrainingTest)
ExternalData_add_test(${Testing_CLP}Data NAME ${testname} COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${Testing_CLP}Test>
  ModuleEntryPoint
   -i 5,3,2 -s 1.5 --mode "(re)train atlas" --imageDir ${TEMP}
   --atlasTreeXML ${TEMP}/Atlas.xml --imageListXML ${TEMP}/List.xml
   --image0 ${INPUT}/TrainingData/AA_cbq_000.nii.gz --segmentation0 ${INPUT}/TrainingData/AA_seg_000.nii.gz
   --image1 ${INPUT}/TrainingData/AB_cbq_000.nii.gz --segmentation1 ${INPUT}/TrainingData/AB_seg_000.nii.gz
   --image2 ${INPUT}/TrainingData/AC_cbq_000.nii.gz --segmentation2 ${INPUT}/TrainingData/AC_seg_000.nii.gz
   --image3 ${INPUT}/TrainingData/AD_cbq_000.nii.gz --segmentation3 ${INPUT}/TrainingData/AD_seg_000.nii.gz
  )
set_property(TEST ${testname} PROPERTY LABELS MABMIS)

# Test using the atlas
make_directory(${TEMP}/SegmentationResults)
set(testname ${Testing_CLP}SegmentationTest)
ExternalData_add_test(${Testing_CLP}Data NAME ${testname} COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${Testing_CLP}Test>
  --compareRadiusTolerance 3 --compare DATA{${BASELINE}/UA_cbq_000_seg.nii.gz} ${TEMP}/SegmentationResults/UA-label.nrrd
  ModuleEntryPoint
   -i 5,3,2 -s 1.5 --mode "Direct invoke" --imageDir ${TEMP}
   --atlasTreeXML ${TEMP}/Atlas.xml --imageListXML ${TEMP}/List.xml
   --image0 ${INPUT}/TestData/UA_cbq_000.nii.gz --segmentation0 ${TEMP}/SegmentationResults/UA-label.nrrd
   --image1 ${INPUT}/TestData/UB_cbq_000.nii.gz --segmentation1 ${TEMP}/SegmentationResults/UB-label.nrrd
   --image2 ${INPUT}/TestData/VA_cbq_000.nii.gz --segmentation2 ${TEMP}/SegmentationResults/VA-label.nrrd
   --image3 ${INPUT}/TestData/VB_cbq_000.nii.gz --segmentation3 ${TEMP}/SegmentationResults/VB-label.nrrd
  )
set_property(TEST ${testname} PROPERTY LABELS MABMIS)
set_property(TEST ${testname} APPEND PROPERTY DEPENDS ${Testing_CLP}TrainingTest)

set(testname ${Testing_CLP}SegmentationTest2)
ExternalData_add_test(${Testing_CLP}Data NAME ${testname} COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${Testing_CLP}Test>
  --compareRadiusTolerance 3 --compare DATA{${BASELINE}/UA_cbq_000_seg.nii.gz} ${TEMP}/SegmentationResults/UA-label.nrrd
  ModuleEntryPoint
   -i 5,3,2 -s 1.5 --mode "Direct invoke" --imageDir ${TEMP}
   --atlasTreeXML ${TEMP}/Atlas.xml #no imageListXML parameter
   --image0 ${INPUT}/TestData/UA_cbq_000.nii.gz --segmentation0 ${TEMP}/SegmentationResults/UA-label.nrrd
  )
set_property(TEST ${testname} PROPERTY LABELS MABMIS)
set_property(TEST ${testname} APPEND PROPERTY DEPENDS ${Testing_CLP}TrainingTest)
