/*=========================================================================
  Copyright (c) IDEA LAB, UNC-Chapel Hill, 2013.

     MABMIS (Multi-Atlas-Based Multi-Image Segmentation)

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "Training.h"

// for math
#include <vcl_iostream.h>
#include <vnl/vnl_random.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matlab_print.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_svd_economy.h>
#include <string>
#include <vector>

// basic itk
#include "itkVector.h"

// registration
#include "itkImageRegistrationMethod.h"
#include "itkSymmetricForcesDemonsRegistrationFilter.h"

// interpolators
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

// filter
#include "itkResampleImageFilter.h"

#include "itkCastImageFilter.h"
#include "itkWarpImageFilter.h"

#include "itkHistogramMatchingImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkWarpVectorImageFilter.h"
#include "itkInverseDisplacementFieldImageFilter.h"

// for affine transformation
#include "itkTransform.h"
#include "itkAffineTransform.h"
#include "itkImageRegistrationMethod.h"

// for Diffeomorphic Demons
#include <itkCommand.h>
#include <itkDiffeomorphicDemonsRegistrationFilter.h>
#include <itkMultiResolutionPDEDeformableRegistration.h>

#include <metaCommand.h>


DeformationFieldType::SpacingType   df_spacing;
DeformationFieldType::DirectionType df_direction;
DeformationFieldType::PointType     df_origin;

typedef itk::Vector<ShortPixelType, ImageDimension>      ShortVectorPixelType;
typedef itk::Image<ShortVectorPixelType, ImageDimension> ShortDeformationFieldType;
typedef itk::ImageFileWriter<ShortDeformationFieldType>  ShortDeformationFieldWriterType;

void SearchRootAmongAtlases(std::vector<std::string> imgfilenames, int & root);

int RegistrationBetweenRootandAtlases(int root, std::vector<std::string> imageFileNames, std::vector<int> iterations, double sigma);

int BuildStatisticalDeformationModel(int root,  std::vector<std::string> imageFileNames, int simulatedAtlasSize);

std::vector<std::string> GenerateSimulatedData(int root, std::vector<std::string> imgfilenames, int simulatedAtlasSize);

void strtrim(std::string& str)
{
  std::string::size_type pos = str.find_last_not_of(' ');

  if( pos != std::string::npos )
    {
    str.erase(pos + 1);
    pos = str.find_first_not_of(' ');
    if( pos != std::string::npos )
      {
      str.erase(0, pos);
      }
    }
  else
    {
    str.erase(str.begin(), str.end() );
    }
}

int Training( itk::MABMISImageData* trainingData, std::string outputFile,
              std::vector<int> iterations, double sigma)
{
  // sanity check
	std::cout << "m_DataDirectory: " << trainingData->m_DataDirectory << std::endl;
  if( trainingData->m_SegmentationFileNames.size() != trainingData->m_ImageFileNames.size() )
    {
    std::cerr << "The numbers of image files and segmentation files are NOT equal!!!" << std::endl;
    return -1;
    }
  if( trainingData->m_ImageFileNames.size() == 0 )
    {
    std::cerr << "No Atlas is provided. QUIT. " << std::endl;
    return -1;
    }
  if( trainingData->m_NumberImageData != trainingData->m_ImageFileNames.size() )
    {
    trainingData->m_NumberImageData = trainingData->m_ImageFileNames.size();
    }

  // step 0: files and folders
  // get output path from the output file name
  size_t      sep = outputFile.find_last_of(FILESEP);
  std::string outputFolder = "";
  if( sep != std::string::npos )
    {
    outputFolder = outputFile.substr(0, sep) + FILESEP;
    }
  else
    {
    sep = -1;
    }
  // get output file name without extension, and use it to create a directory to save trained atlas
  size_t      fnSep = outputFile.find_last_of('.');
  std::string outputxmlname = "";
  if( fnSep != std::string::npos )
    {
    outputxmlname = outputFile.substr(sep + 1, fnSep - sep - 1);
    }
  else
    {
    outputxmlname = outputFile.substr(sep + 1, std::string::npos);
    }

  std::string atlasFolder;
  if( outputFolder.empty() )
    {
    atlasFolder = outputxmlname;
    }
  else
    {
    atlasFolder = outputFolder + outputxmlname;
    }
  atlasFolder = ReplacePathSepForOS(atlasFolder); 
  if( !atlasFolder.empty() )
    {
    atlasFolder = atlasFolder + FILESEP;
    }
  itksys::SystemTools::MakeDirectory(atlasFolder.c_str() );

  // copy training data into the atlas folder
  std::vector<std::string> imageFiles(trainingData->m_NumberImageData);
  std::vector<std::string> segmentFiles(trainingData->m_NumberImageData);
  for( int i = 0; i < imageFiles.size(); ++i )
    {
    imageFiles[i] = ReplacePathSepForOS(atlasFolder + trainingData->m_ImageFileNames[i]);
    segmentFiles[i] = ReplacePathSepForOS(atlasFolder + trainingData->m_SegmentationFileNames[i]);

    // copy
    std::string fromImagefile = trainingData->m_ImageFileNames[i];
    if( !trainingData->m_DataDirectory.empty() )
      {
      fromImagefile  = ReplacePathSepForOS(trainingData->m_DataDirectory + trainingData->m_ImageFileNames[i]);
      }
	
    bool ret1=itksys::SystemTools::CopyFileAlways(fromImagefile.c_str(), imageFiles[i].c_str() );
	if (!ret1) 
      {
	  std::cerr << "ERROR: Cannot copy atlas image file to trained atlas folder!" << std::endl;
	  return -1;
      }

    // if file extension is '.nii.gz', also copy the 'img file'
    if( itksys::SystemTools::Strucmp(fromImagefile.substr(fromImagefile.size() - 3, 3).c_str(), "hdr") == 0 )
      {
      std::string fromfileImg = fromImagefile.substr(0, fromImagefile.size() - 3) + "img";
      std::string tofileImg = imageFiles[i].substr(0, imageFiles[i].size() - 3) + "img";
      itksys::SystemTools::CopyFileAlways(fromfileImg.c_str(), tofileImg.c_str() );
      }

    std::string fromSegfile = trainingData->m_DataDirectory + trainingData->m_SegmentationFileNames[i];

    bool ret2 = itksys::SystemTools::CopyFileAlways(fromSegfile.c_str(), segmentFiles[i].c_str() );
	if (!ret1) 
      {
	  std::cerr << "ERROR: Cannot copy atlas image file to trained atlas folder!" << std::endl;
	  return -1;
      }
    // if file extension is '.nii.gz', also copy the 'img file'
    if( itksys::SystemTools::Strucmp(fromSegfile.substr(fromSegfile.size() - 3, 3).c_str(), "hdr") == 0 )
      {
      std::string fromfileImg = fromSegfile.substr(0, fromSegfile.size() - 3) + "img";
      std::string tofileImg = segmentFiles[i].substr(0, segmentFiles[i].size() - 3) + "img";
      itksys::SystemTools::CopyFileAlways(fromfileImg.c_str(), tofileImg.c_str() );
      }
    }

  // step 1: find root of the atlas
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "1. Find root atlas ... " << std::endl;;
  ////////////////////////////////
  // build tree on atlases to find the root
  int root = -1;    // root node of tree
  SearchRootAmongAtlases(imageFiles, root);
  std::cout << "Done." << std::endl;;

  // step 2.1: register all atlas to the root
  //////////////////////////////
  // registration between root and other atlases
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "2. Run coarse registration ... " << std::endl;
  int val = RegistrationBetweenRootandAtlases(root, imageFiles, iterations, sigma);
  if( val != 0 )
    {
    return -1;
    }
  std::cout << "Done!" << std::endl;

  int simulatedAtlasSize = 2 * imageFiles.size();

  std::cout << "---------------------------------------" << std::endl;
  std::cout << "3. Build statistical deformation model..." << std::endl;
  if( BuildStatisticalDeformationModel(root, imageFiles, simulatedAtlasSize) != 0 )
    {
    return -1;
    }
  std::cout << "Done. " << std::endl;

  // step 2.5: generate templaet images
  /////////////////////////
  // generate simulated template with deformation field
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "4. Generate simulated templates ... " << std::endl;
  std::vector<std::string> simulatedImageFiles = GenerateSimulatedData(root, imageFiles, simulatedAtlasSize);
  std::cout << "Done. " << std::endl;

  /////////////////////////
  // delete subsampled deformation field
  // if (!isDebug)	if (false)
    {
    for( int i = 0; i < imageFiles.size(); ++i )
      {
      if( i == root )
        {
        continue;
        }
      char i_str[10], root_str[10];
      sprintf(i_str, "%03d", i);
      sprintf(root_str, "%03d", root);
      std::string subdeformationFieldFileName = std::string(i_str) + "_to_" + std::string(root_str)
        + "deform_000_sub.nii.gz";
      subdeformationFieldFileName = outputFolder + subdeformationFieldFileName;
      basicoperator->RemoveFile(subdeformationFieldFileName);
      }
    }

  //////////////////////////////
  // build the combinative tree
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "5. Build combinative tree ... " << std::endl;
  int             totalAtlasSize = imageFiles.size() + simulatedAtlasSize;
  int             tree_size = totalAtlasSize;
  vnl_vector<int> tree(tree_size);    // each tree

  // fill distance matrix
  std::vector<std::string> atlasfilenames;
  for( int i = 0; i < tree_size; ++i )
    {
    std::string atlasfilename;

    if( i < imageFiles.size() )
      {
      atlasfilename.append(imageFiles[i]);
      }
    else
      {
      atlasfilename.append(simulatedImageFiles[i - imageFiles.size()]);
      }

    // atlasfilename = outputFolder + atlasfilename;
    atlasfilenames.push_back(atlasfilename);
    }

  tree = treeoperator->BuildCombinativeTree(root, atlasfilenames, tree_size, tree); // cout << "pass:
                                                                                    // BuildCombinativeTree " << endl;

  std::cout << "Done. " << std::endl;

  // write trained atlas to xml file
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "6. Write to output file..." << std::endl;
  itk::MABMISAtlasXMLFileWriter::Pointer atlasWriter = itk::MABMISAtlasXMLFileWriter::New();
  itk::MABMISAtlas                       atlas;
  atlas.m_NumberAllAtlases = totalAtlasSize;
  atlas.m_NumberSimulatedAtlases = simulatedAtlasSize;
  atlas.m_AtlasDirectory = std::string(".") + FILESEP + outputxmlname;
  atlas.m_Tree.resize(tree_size);
  atlas.m_IsSimulatedImage.resize(tree_size);
  atlas.m_RealImageIDs.resize(totalAtlasSize - simulatedAtlasSize);
  atlas.m_SimulatedImageIDs.resize(simulatedAtlasSize);
  for( int i = 0; i < tree_size; ++i )
    {
    atlas.m_Tree[i] = tree[i];
    }
  atlas.m_AtlasFilenames.resize(atlasfilenames.size() );
  int countRealImages = 0, countSimulatedImages = 0;
  for( int i = 0; i < atlasfilenames.size(); ++i )
    {
    const size_t      sep = atlasfilenames[i].find_last_of(FILESEP);
    std::string fname = atlasfilenames[i];
    if( sep != std::string::npos )
      {
      fname = atlasfilenames[i].substr(sep + 1, std::string::npos);
/*			fname = atlasfilenames[i].substr(0, sep);
      sep = fname.find_last_of(FILESEP);
      if (sep != std::string::npos)
        fname = atlasfilenames[i].substr(sep+1, std::string::npos);
      else
        fname = atlasfilenames[i];
        */
      }
    else
      {
      fname = atlasfilenames[i];
      }

    atlas.m_AtlasFilenames[i] = fname;

    if( i < totalAtlasSize - simulatedAtlasSize )
      {
      atlas.m_IsSimulatedImage[i] = false;
      atlas.m_RealImageIDs[countRealImages++] = i;
      }
    else
      {
      atlas.m_IsSimulatedImage[i] = true;
      atlas.m_SimulatedImageIDs[countSimulatedImages++] = i;
      }
    }
  // atlas segmentation files
  atlas.m_AtlasSegmentationFilenames.resize(segmentFiles.size() );
  for( int i = 0; i < segmentFiles.size(); ++i )
    {
    const size_t      sep = segmentFiles[i].find_last_of(FILESEP);
    std::string fname = segmentFiles[i];
    if( sep != std::string::npos )
      {
      fname = segmentFiles[i].substr(sep + 1, std::string::npos);
/*			fname = segmentFiles[i].substr(0, sep);
      sep = fname.find_last_of(FILESEP);
      if (sep != std::string::npos)
        fname = segmentFiles[i].substr(sep+1, std::string::npos);
      else
        fname = segmentFiles[i];
*/
      }
    else
      {
      fname = segmentFiles[i];
      }
    atlas.m_AtlasSegmentationFilenames[i] = fname;
    }

  treeoperator->FindRoot(tree, tree_size, atlas.m_TreeRoot);
  treeoperator->GetTreeHeight(tree, tree_size, atlas.m_TreeHeight);

  atlasWriter->SetObject(&atlas);
  atlasWriter->SetFilename(outputFile);
  atlasWriter->WriteFile();

  std::cout << "Done." << std::endl;

  return EXIT_SUCCESS;
}

void SearchRootAmongAtlases(std::vector<std::string> imgfilenames, int & root)
{
  int tree_size = imgfilenames.size();

  vnl_matrix<double> distanceMatrix(tree_size, tree_size);
  for( int i = 0; i < tree_size; ++i )
    {
    for( int j = 0; j < tree_size; ++j )
      {
      distanceMatrix[i][j] = 0.0;
      }
    }

  vnl_vector<int> tree(tree_size);  // each tree

  // fill distance matrix
  imgoperator->PairwiseDistanceAmongImages(imgfilenames, tree_size, distanceMatrix);
  if( isDebug )
    {
    std::string distFileName = "dist_all_atlases.txt";
    basicoperator->SaveMatrix2File(distanceMatrix, tree_size, tree_size, distFileName);
    }

  // build MST
  tree = treeoperator->generateMSTFromMatrix(distanceMatrix, tree_size, tree);

  treeoperator->FindRoot(tree, tree_size, root);

  // save distance matrix and tree
  if( isDebug )
    {
    std::string treeFileName = "tree_all_atlases.txt";
    treeoperator->SaveTreeWithInfo(tree, tree_size, treeFileName);
    }

  std::cout << "The root is " << root << ": " << imgfilenames[root] << ". " << std::endl;
}

int RegistrationBetweenRootandAtlases(int root, std::vector<std::string> imageFileNames,
                                      std::vector<int> iterations, double sigma)
{
  // get the path from the filename
  const size_t      sep = imageFileNames[0].find_last_of(FILESEP);
  std::string outputFolder = "";

  if( sep != std::string::npos )
    {
    outputFolder = imageFileNames[0].substr(0, sep);
    }
  if( !outputFolder.empty() )
    {
    outputFolder = outputFolder + FILESEP;
    }

  const int atlas_size = imageFileNames.size();
  std::cout << "Register between root and atlases ... " << std::endl;
  for( int i = 0; i < atlas_size; ++i )
    {
    std::cout << i << ", ";
    if( i == root )
      {
      continue;
      }
    // if not the root image do registration
    std::string fixedImageFileName = imageFileNames[root];
    std::string movingImageFileName = imageFileNames[i];

    char i_str[10], root_str[10];
    sprintf(i_str, "%03d", i);
    sprintf(root_str, "%03d", root);
    std::string deformedImageFileName = std::string(i_str) + "_to_" + std::string(root_str) + "_cbq_reg.nii.gz";
    std::string      deformationFieldFileName = std::string(i_str) + "_to_" + std::string(root_str) + "_deform_000.nii.gz";

    deformedImageFileName  = outputFolder + deformedImageFileName;
    deformationFieldFileName = outputFolder + deformationFieldFileName;

    if( !itksys::SystemTools::FileExists(deformationFieldFileName.c_str(), true) )
      {
      // rough registration

      int val = regoperator->DiffeoDemonsRegistrationWithParameters(
          fixedImageFileName, movingImageFileName,
          deformedImageFileName, deformationFieldFileName,
          sigma + 0.5, doHistMatchTrue, iterations);
      if( val != 0 )
        {
        std::cout << "Cannot perform registration between :" << std::endl;
        std::cout << fixedImageFileName << " and " << std::endl;
        std::cout << movingImageFileName << "." << std::endl;
        std::cout << "Please verify the input images are correct." << std::endl;
        }
      }

    // ---------------------------------------------------------------------
    // ------------- Generate the inverse deformation field

    DeformationFieldType::Pointer deformationField = 0;
    dfoperator->ReadDeformationField(deformationFieldFileName, deformationField);

    DeformationFieldType::Pointer inversedDeformationField = DeformationFieldType::New();
    inversedDeformationField->SetRegions(deformationField->GetLargestPossibleRegion() );
    inversedDeformationField->SetSpacing(deformationField->GetSpacing() );
    inversedDeformationField->SetDirection(deformationField->GetDirection() );
    inversedDeformationField->SetOrigin(deformationField->GetOrigin() );
    inversedDeformationField->Allocate();

    dfoperator->InverseDeformationField3D(deformationField, inversedDeformationField);

    std::string invDeformedImageFileName = std::string(root_str) + "_to_" + std::string(i_str) + "_cbq_reg.nii.gz";
    std::string      invDeformationFileName = std::string(root_str) + "_to_" + std::string(i_str) + "_deform_000.nii.gz";
    invDeformedImageFileName = outputFolder + invDeformedImageFileName;
    invDeformationFileName = outputFolder + invDeformationFileName;

    std::string invDeformedSegmentFileName = std::string(root_str) + "_to_" + std::string(i_str) + "seg_000.nii.gz";
    invDeformedSegmentFileName = outputFolder + invDeformedSegmentFileName;

    dfoperator->WriteDeformationField(invDeformationFileName, inversedDeformationField);

    // update
    regoperator->DiffeoDemonsRegistrationWithInitialWithParameters(
      movingImageFileName, fixedImageFileName,
      invDeformationFileName,
      invDeformedImageFileName, invDeformationFileName,
      sigma, doHistMatch, iterations);

    if( isEvaluate )
      {
      dfoperator->ApplyDeformationFieldAndWriteWithTypeWithFileNames(
        fixedImageFileName,
        invDeformationFileName,
        invDeformedSegmentFileName, false);
      }
    }
  std::cout << std::endl;
  std::cout << "Done. " << std::endl;

  return 0;
}

int BuildStatisticalDeformationModel(int root,  std::vector<std::string> imageFileNames, int simulatedAtlasSize)
{
  ///////////////////////////////
  // do PCA simulation
  // std::cout << "Build deformation field model ... ";
  std::vector<std::string> allDeformationFieldFileNames;

  // get the path from the filename
  const size_t      sep = imageFileNames[0].find_last_of(FILESEP);
  std::string outputFolder = "";
  if( sep != std::string::npos )
    {
    outputFolder = imageFileNames[0].substr(0, sep);
    }
  if( !outputFolder.empty() )
    {
    outputFolder = outputFolder + FILESEP;
    }

  int atlas_size = imageFileNames.size();

  datasimulator->SetRoot(root);
  datasimulator->SetAtlasSize(atlas_size);
  datasimulator->SetSimulateSize(simulatedAtlasSize);
  for( int i = 0; i < atlas_size; ++i )
    {
    if( i == root )
      {
      continue;
      }

    char i_str[10], root_str[10];
    sprintf(i_str, "%03d", i);
    sprintf(root_str, "%03d", root);

    std::string deformationFieldFileName = std::string(i_str) + "_to_" + std::string(root_str) + "_deform_000.nii.gz";
    deformationFieldFileName = outputFolder + deformationFieldFileName;
    allDeformationFieldFileNames.push_back(deformationFieldFileName);
    }

  return datasimulator->DoPCATraining(allDeformationFieldFileNames, atlas_size - 1, imageFileNames, root);
}

// return: simulated template file names
std::vector<std::string> GenerateSimulatedData(int root, std::vector<std::string> imageFiles, int simulatedAtlasSize)
{
  // output folder
  const size_t      sep = imageFiles[0].find_last_of(FILESEP);
  std::string outputFolder = "";

  if( sep != std::string::npos )
    {
    outputFolder = imageFiles[0].substr(0, sep);
    }
  if( !outputFolder.empty() )
    {
    outputFolder = outputFolder + FILESEP;
    }

  std::vector<std::string> simulateDeformationFieldFileNames(0);
  std::vector<std::string> simulateTemplateFileNames(0);
  std::cout << "Generating simulated images: " << std::endl;
  for( int i = 0; i < simulatedAtlasSize; ++i )
    {
    std::cout << i << ", ";

    std::string index_string;
    basicoperator->myitoa( i, index_string, 3 );

    std::string simulateDeformationFieldFileName = "simulated_deform_" + index_string + ".nii.gz";
    std::string simulateTemplateFileName = "simulated_cbq_" + index_string + ".nii.gz";
    simulateDeformationFieldFileName = outputFolder + simulateDeformationFieldFileName;
    simulateTemplateFileName = outputFolder + simulateTemplateFileName;

    simulateDeformationFieldFileNames.push_back(simulateDeformationFieldFileName);
    simulateTemplateFileNames.push_back(simulateTemplateFileName);

    // load simulated deformation field
    DeformationFieldType::Pointer df = 0;
    dfoperator->ReadDeformationField(simulateDeformationFieldFileNames[i], df);

    DeformationFieldType::Pointer invdf = DeformationFieldType::New();
    invdf->SetRegions(df->GetLargestPossibleRegion() );
    invdf->SetSpacing(df->GetSpacing() );
    invdf->SetDirection(df->GetDirection() );
    invdf->SetOrigin(df->GetOrigin() );
    invdf->Allocate();

    dfoperator->InverseDeformationField3D(df, invdf);
    InternalImageType::Pointer rootImg = 0;
    imgoperator->ReadImage(imageFiles[root], rootImg);
    const itk::ImageIOBase::IOComponentType ioType = imgoperator->GetIOPixelType(imageFiles[root]);
    InternalImageType::Pointer simImg = 0;
    dfoperator->ApplyDeformationField(rootImg, invdf, simImg, true);
    imgoperator->WriteImage(simulateTemplateFileNames[i], simImg, ioType);
    }
  std::cout << std::endl;
  
  return simulateTemplateFileNames;
}
