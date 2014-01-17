cmake_minimum_required(VERSION 2.8.9)

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  set(EXTENSION_NAME MABMIS)
  set(EXTENSION_HOMEPAGE "http://wiki.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/MABMIS")
  set(EXTENSION_CATEGORY "Segmentation")
  set(EXTENSION_CONTRIBUTORS "Xiaofeng Liu (GE GLobal Research), Minjeong Kim (UNC), Dinggang Shen (UNC), Jim Miller (GE GLobal Research)")
  set(EXTENSION_ICONURL "http://wiki.slicer.org/slicerWiki/images/e/e2/MABMIS_Icon.png")
  set(EXTENSION_DESCRIPTION "Multi-Atlas Based Group Segmentation")
  set(EXTENSION_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/2/2c/MABMIS_trainning_GUI.png http://wiki.slicer.org/slicerWiki/images/6/64/MABMIS_testing_GUI.png http://wiki.slicer.org/slicerWiki/images/c/c3/MABMIS_algorithm.png")
endif()

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

#-----------------------------------------------------------------------------
set(MODULE_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
set(MODULE_NAME IGR3D_MABMIS_Training)

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  ADDITIONAL_SRCS
	itkMABMISAtlasXMLFile.cxx
  #EXECUTABLE_ONLY
  )

  message(STATUS ${ITK_LIBRARIES})
#-----------------------------------------------------------------------------
set(MODULE_NAME IGR3D_MABMIS_Testing)

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  ADDITIONAL_SRCS
	itkMABMISAtlasXMLFile.cxx
  #EXECUTABLE_ONLY
  )


#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()




################################################
