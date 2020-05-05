/*=========================================================================
  Copyright (c) IDEA LAB, UNC-Chapel Hill, 2013.

     MABMIS (Multi-Atlas-Based Multi-Image Segmentation)

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "IGR3D_MABMIS_TestingCLP.h"
#include "Testing.h"

int main( int argc, char *argv[] )
{
  PARSE_ARGS;

  // step 1: read the test image list
  itk::MABMISImageDataXMLFileReader::Pointer imageListXMLReader = itk::MABMISImageDataXMLFileReader::New();
  ImageListXML = ReplacePathSepForOS(ImageListXML); 
  AtlaseTreeXML = ReplacePathSepForOS(AtlaseTreeXML); 
  OutputFolder = ReplacePathSepForOS(OutputFolder); 

  imageListXMLReader->SetFilename(ImageListXML);

  try
    {
    imageListXMLReader->GenerateOutputInformation();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while reading file" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
  itk::MABMISImageData * inputImageData = imageListXMLReader->GetOutputObject();

  // if the data path is empty, use the path of the xml file instead
  if( inputImageData->m_DataDirectory.size() <= 1 )
    {
    const size_t dir_sep = ImageListXML.find_last_of(FILESEP);
    if( dir_sep != std::string::npos )
      {
      inputImageData->m_DataDirectory = ImageListXML.substr(0, dir_sep);
      }
    else
      {
      inputImageData->m_DataDirectory.resize(0);
      }
    }
  if( !inputImageData->m_DataDirectory.empty() )
    {
    if( !(inputImageData->m_DataDirectory[inputImageData->m_DataDirectory.size() - 1] == FILESEP) )
      {
      inputImageData->m_DataDirectory = inputImageData->m_DataDirectory + FILESEP;
      }
    }

  // output directory
  if( OutputFolder.empty() )
    {
    inputImageData->m_OutputDirectory = inputImageData->m_DataDirectory;
    }
  else
    {
    if( !(OutputFolder[OutputFolder.size() - 1] == FILESEP) )
      {
      OutputFolder = OutputFolder + FILESEP;
      }
    inputImageData->m_OutputDirectory = OutputFolder;
    }
  itksys::SystemTools::MakeDirectory(inputImageData->m_OutputDirectory);

  // load the tree-structured atlases that is generated in the training step
  itk::MABMISAtlasXMLFileReader::Pointer treeAtlasXMLReader = itk::MABMISAtlasXMLFileReader::New();
  treeAtlasXMLReader->SetFilename(AtlaseTreeXML);
  try
    {
    treeAtlasXMLReader->GenerateOutputInformation();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while reading file" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
  itk::MABMISAtlas * atlasTree = treeAtlasXMLReader->GetOutputObject();

  // set the atlas path as the same path as the xml file.
  const size_t dir_sep = AtlaseTreeXML.find_last_of(FILESEP);
  if( dir_sep != std::string::npos )
    {
    atlasTree->m_AtlasDirectory = ReplacePathSepForOS(AtlaseTreeXML.substr(0, dir_sep + 1) + atlasTree->m_AtlasDirectory);
    }

  if( atlasTree->m_AtlasDirectory.size() == 0 )
    {
    atlasTree->m_AtlasDirectory = ReplacePathSepForOS(".");
    }

  if( !(atlasTree->m_AtlasDirectory[atlasTree->m_AtlasDirectory.size() - 1] == FILESEP) )
    {
    atlasTree->m_AtlasDirectory = ReplacePathSepForOS(atlasTree->m_AtlasDirectory + FILESEP );
    }

  int retVal = Testing( inputImageData, atlasTree, iterations, SmoothingKernelSize);
  delete inputImageData;
  delete atlasTree;
  return retVal;
}

