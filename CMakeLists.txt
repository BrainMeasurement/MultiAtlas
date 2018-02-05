cmake_minimum_required(VERSION 3.5)

#-----------------------------------------------------------------------------
project(MABMIS)
set(EXTENSION_HOMEPAGE "http://wiki.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/MABMIS")
set(EXTENSION_CATEGORY "Segmentation")
set(EXTENSION_CONTRIBUTORS "Xiaofeng Liu (GE Global Research), Minjeong Kim (UNC), Dinggang Shen (UNC), Jim Miller (GE Global Research)")
set(EXTENSION_ICONURL "http://wiki.slicer.org/slicerWiki/images/e/e2/MABMIS_Icon.png")
set(EXTENSION_DESCRIPTION "Multi-Atlas Based Group Segmentation")
set(EXTENSION_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/2/2c/MABMIS_trainning_GUI.png http://wiki.slicer.org/slicerWiki/images/6/64/MABMIS_testing_GUI.png http://wiki.slicer.org/slicerWiki/images/c/c3/MABMIS_algorithm.png")

find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

set(MABMIS_ITK_COMPONENTS
  ITKIOImageBase
  ITKTransform
  ITKIOTransformBase
  ITKIOXML
  ITKImageGrid
  ITKImageFunction
  )
find_package(ITK 4.6 COMPONENTS ${MABMIS_ITK_COMPONENTS} REQUIRED)
set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1) # See Libs/ITKFactoryRegistration/CMakeLists.txt
include(${ITK_USE_FILE})

find_package(Slicer)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
set(MODULE_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  )
message(STATUS ${MODULE_TARGET_LIBRARIES})

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME IGR3D_MABMIS_Training
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  ADDITIONAL_SRCS
	itkMABMISAtlasXMLFile.cxx
	Training.cxx
	commonMABMIS.cxx
  #EXECUTABLE_ONLY
  )

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME IGR3D_MABMIS_Testing
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  ADDITIONAL_SRCS
	itkMABMISAtlasXMLFile.cxx
	Testing.cxx
	commonMABMIS.cxx
  #EXECUTABLE_ONLY
  )

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME IGR3D_MABMIS_Direct_Invoke
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  ADDITIONAL_SRCS
	itkMABMISAtlasXMLFile.cxx
	Testing.cxx
	Training.cxx
	commonMABMIS.cxx
  #EXECUTABLE_ONLY
  )


#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
