/*=========================================================================
  Copyright (c) IDEA LAB, UNC-Chapel Hill, 2013.

     MABMIS (Multi-Atlas-Based Multi-Image Segmentation)

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
// #pragma warning ( disable : 4244 )
#endif

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
#include "itkImage.h"
#include "itkVector.h"

// registration
#include "itkImageRegistrationMethod.h"
#include "itkSymmetricForcesDemonsRegistrationFilter.h"

// interpolators
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

// reader / writer
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// filter
#include "itkResampleImageFilter.h"

#include "itkCastImageFilter.h"
#include "itkWarpImageFilter.h"

#include "itkHistogramMatchingImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkDivideByConstantImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkWarpVectorImageFilter.h"
#include "itkInverseDeformationFieldImageFilter.h"

// for affine transformation
#include "itkTransform.h"
#include "itkAffineTransform.h"
#include "itkImageRegistrationMethod.h"

// for Diffeomorphic Demons
#include <itkCommand.h>
#include <itkDiffeomorphicDemonsRegistrationFilter.h>
#include <itkMultiResolutionPDEDeformableRegistration.h>

// including itksys::SystemTools::MakeDirectory(char*)
#include <itksys/SystemTools.hxx>
#include <metaCommand.h>

#define DBL_MAX         1.7976931348623158e+308 /* max value */

// To include all related header files
#include "IGR3D_MABMIS_TrainingCLP.h"
#include "itkMABMISImageOperationFilter.h"
#include "itkMABMISDeformationFieldFilter.h"
#include "itkMABMISSimulateData.h"
#include "itkMABMISImageRegistrationFilter.h"
#include "itkMABMISTreeOperation.h"
#include "itkMABMISBasicOperationFilter.h"

#include "itkMABMISAtlasXMLFile.h"

#ifdef _WIN32
#define FILESEP '\\'
#else
#define  FILESEP '/'
#endif

using namespace std;

typedef double CoordinateRepType;
const   unsigned int SpaceDimension = ImageDimension;

// basic data type
typedef unsigned char                                  CharPixelType;  // for image IO usage
typedef float                                          FloatPixelType; // for
typedef int                                            IntPixelType;
typedef short                                          ShortPixelType;
typedef float                                          InternalPixelType; // for internal processing usage
typedef itk::Vector<InternalPixelType, ImageDimension> VectorPixelType;

// basic image type
typedef itk::Image<CharPixelType, ImageDimension>     CharImageType;
typedef itk::Image<IntPixelType, ImageDimension>      IntImageType;
typedef itk::Image<ShortPixelType, ImageDimension>    ShortImageType;
typedef itk::Image<FloatPixelType, ImageDimension>    FloatImageType;
typedef itk::Image<InternalPixelType, ImageDimension> InternalImageType;
typedef itk::Image<VectorPixelType, ImageDimension>   DeformationFieldType;

// basic iterator type
typedef itk::ImageRegionIterator<DeformationFieldType> DeformationFieldIteratorType;
typedef itk::ImageRegionIterator<InternalImageType>    InternalImageIteratorType;
typedef itk::ImageRegionIterator<CharImageType>        CharImageIteratorType;

// basic image reader/writer related type
typedef itk::ImageFileReader<CharImageType>     CharImageReaderType;
typedef itk::ImageFileReader<InternalImageType> InternalImageReaderType;
typedef itk::ImageFileWriter<InternalImageType> InternalImageWriterType;

typedef itk::WarpImageFilter<InternalImageType, InternalImageType, DeformationFieldType> InternalWarpFilterType;
typedef itk::ImageFileWriter<CharImageType>                                              CharImageWriterType;
typedef itk::ImageFileWriter<IntImageType>                                               IntImageWriterType;
typedef itk::ImageFileWriter<FloatImageType>                                             FloatImageWriterType;
typedef itk::ImageFileWriter<ShortImageType>                                             ShortImageWriterType;

typedef itk::ImageFileReader<DeformationFieldType> DeformationFieldReaderType;
typedef itk::ImageFileWriter<DeformationFieldType> DeformationFieldWriterType;

//////////////////////////////////////////////////////////////////////////////
// image filter type
typedef itk::ResampleImageFilter<InternalImageType, InternalImageType>          ResampleFilterType;
typedef itk::HistogramMatchingImageFilter<InternalImageType, InternalImageType> InternalHistMatchFilterType;

////////////////////////////////////////////////////////////////////////////
// operation on deformation fields
typedef itk::WarpVectorImageFilter<DeformationFieldType, DeformationFieldType,
                                   DeformationFieldType>                                      WarpVectorFilterType;
typedef itk::InverseDeformationFieldImageFilter<DeformationFieldType,
                                                DeformationFieldType>
                                                                                              InverseDeformationFieldImageFilterType;
typedef itk::AddImageFilter<DeformationFieldType, DeformationFieldType, DeformationFieldType> AddImageFilterType;

// global bool variables to adjust the  procedure

bool isEvaluate = false; // if false, we do not know the ground-truth of labels
bool isDebug = false;    // false;//true; // if true, print out more information

// -----------------------------------------------------------------------------
// global variables

// int localPatchSize = 1; //(2r+1)*(2r+1)*(2r+1) is the volume of local patch

// demons registration parameters
// int iterInResolutions[4][3]={{5,3,2},{10,5,5},{15,10,5},{20,15,10}};
// int itereach = 2; //
// int itereach0 = 0;int itereach1 = 1;int itereach2 = 2;int itereach3 = 3;
// double sigmaDef = 1.5;
// double sigmaDef10 = 1.0;double sigmaDef15 = 1.5;double sigmaDef20 = 2.0;
// double sigmaDef25 = 2.5;double sigmaDef30 = 3.0;double sigmaDef35 = 3.5;
bool doHistMatch = true;
bool doHistMatchTrue = true; bool doHistMatchFalse = false;

// -----------------------------------------------------------------------------

typedef itk::Statistics::MABMISSimulateData<InternalImageType, InternalImageType> DataSimulatorType;
DataSimulatorType::Pointer datasimulator;

typedef itk::Statistics::MABMISImageOperationFilter<CharImageType, CharImageType> ImageOperationFilterType;
ImageOperationFilterType::Pointer imgoperator;
typedef itk::Statistics::MABMISDeformationFieldFilter<InternalImageType,
                                                      InternalImageType> DeformationFieldOperationFilterType;
DeformationFieldOperationFilterType::Pointer dfoperator;
typedef itk::Statistics::MABMISImageRegistrationFilter<CharImageType, CharImageType> ImageRegistrationFilterType;
ImageRegistrationFilterType::Pointer regoperator;
typedef itk::Statistics::MABMISTreeOperation<InternalImageType, InternalImageType> TreeOperationType;
TreeOperationType::Pointer treeoperator;

typedef itk::Statistics::MABMISBasicOperationFilter<CharImageType, CharImageType> BasicOperationFilterType;
BasicOperationFilterType::Pointer basicoperator;

//
DeformationFieldType::SpacingType   df_spacing;
DeformationFieldType::DirectionType df_direction;
DeformationFieldType::PointType     df_origin;

typedef itk::Vector<ShortPixelType, ImageDimension>      ShortVectorPixelType;
typedef itk::Image<ShortVectorPixelType, ImageDimension> ShortDeformationFieldType;
typedef itk::ImageFileWriter<ShortDeformationFieldType>  ShortDeformationFieldWriterType;

int ReadImgInfo(std::vector<std::string> imgfilenames);

void SearchRootAmongAtlases(std::vector<std::string> imgfilenames, int & root);

int RegistrationBetweenRootandAtlases(int root, std::vector<std::string> imageFileNames, std::vector<int> iterations,
                                      double sigma);

int BuildStatisticalDeformationModel(int root,  std::vector<std::string> imageFileNames, int simulatedAtlasSize);

std::vector<std::string> GenerateSimulatedData(int root, std::vector<std::string> imgfilenames, int simulatedAtlasSize);

void strtrim(string& str)
{
  string::size_type pos = str.find_last_not_of(' ');

  if( pos != string::npos )
    {
    str.erase(pos + 1);
    pos = str.find_first_not_of(' ');
    if( pos != string::npos )
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
  if( !atlasFolder.empty() )
    {
    atlasFolder = atlasFolder + FILESEP;
    }
  itksys::SystemTools::MakeDirectory(atlasFolder.c_str() );

  // copy training data into the atlas folder
  std::vector<std::string> imageFiles(trainingData->m_NumberImageData);
  std::vector<std::string> segmentFiles(trainingData->m_NumberImageData);
  for( int i = 0; i < imageFiles.size(); i++ )
    {
    imageFiles[i] = atlasFolder + trainingData->m_ImageFileNames[i];
    segmentFiles[i] = atlasFolder + trainingData->m_SegmentationFileNames[i];

    // copy
    std::string fromImagefile = trainingData->m_ImageFileNames[i];
    if( !trainingData->m_DataDirectory.empty() )
      {
      fromImagefile  = trainingData->m_DataDirectory + trainingData->m_ImageFileNames[i];
      }
    itksys::SystemTools::CopyFileAlways(fromImagefile.c_str(), imageFiles[i].c_str() );
    // if file extension is '.hdr', also copy the 'img file'
    if( itksys::SystemTools::Strucmp(fromImagefile.substr(fromImagefile.size() - 3, 3).c_str(), "hdr") == 0 )
      {
      std::string fromfileImg = fromImagefile.substr(0, fromImagefile.size() - 3) + "img";
      std::string tofileImg = imageFiles[i].substr(0, imageFiles[i].size() - 3) + "img";
      itksys::SystemTools::CopyFileAlways(fromfileImg.c_str(), tofileImg.c_str() );
      }

    std::string fromSegfile = trainingData->m_DataDirectory + trainingData->m_SegmentationFileNames[i];
    itksys::SystemTools::CopyFileAlways(fromSegfile.c_str(), segmentFiles[i].c_str() );
    // if file extension is '.hdr', also copy the 'img file'
    if( itksys::SystemTools::Strucmp(fromSegfile.substr(fromSegfile.size() - 3, 3).c_str(), "hdr") == 0 )
      {
      std::string fromfileImg = fromSegfile.substr(0, fromSegfile.size() - 3) + "img";
      std::string tofileImg = segmentFiles[i].substr(0, segmentFiles[i].size() - 3) + "img";
      itksys::SystemTools::CopyFileAlways(fromfileImg.c_str(), tofileImg.c_str() );
      }
    }

  // step 1: find root of the atlas
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "Find root atlas ... " << std::endl;;
  ////////////////////////////////
  // build tree on atlases to find the root
  int root = -1;    // root node of tree
  SearchRootAmongAtlases(imageFiles, root);
  std::cout << "Done." << std::endl;;

  // step 2.1: register all atlas to the root
  //////////////////////////////
  // registration between root and other atlases
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "Run coarse registration ... " << std::endl;
  int val = RegistrationBetweenRootandAtlases(root, imageFiles, iterations, sigma);
  if( val != 0 )
    {
    return -1;
    }
  std::cout << "Done!" << std::endl;

  int simulatedAtlasSize = 2 * imageFiles.size();

  std::cout << "---------------------------------------" << std::endl;
  std::cout << "Build statistical deformation model..." << std::endl;
  if( BuildStatisticalDeformationModel(root, imageFiles, simulatedAtlasSize) != 0 )
    {
    return -1;
    }
  std::cout << "Done. " << std::endl;

  // step 2.5: generate templaet images
  /////////////////////////
  // generate simulated template with deformation field
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "Generate simulated templates ... " << std::endl;
  std::vector<std::string> simulatedImageFiles = GenerateSimulatedData(root, imageFiles, simulatedAtlasSize);
  std::cout << "Done. " << std::endl;

  /////////////////////////
  // delete subsampled deformation field
  // if (!isDebug)	if (false)
    {
    for( int i = 0; i < imageFiles.size(); i++ )
      {
      if( i == root )
        {
        continue;
        }
      char i_str[10], root_str[10];
      sprintf(i_str, "%03d", i);
      sprintf(root_str, "%03d", root);
      std::string subdeformationFieldFileName = std::string(i_str) + "_to_" + std::string(root_str)
        + "deform_000_sub.mha";
      subdeformationFieldFileName = outputFolder + subdeformationFieldFileName;
      basicoperator->RemoveFile(subdeformationFieldFileName);
      }
    }

  //////////////////////////////
  // build the combinative tree
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "Build combinative tree ... " << std::endl;
  int             totalAtlasSize = imageFiles.size() + simulatedAtlasSize;
  int             tree_size = totalAtlasSize;
  vnl_vector<int> tree(tree_size);    // each tree

  // fill distance matrix
  vector<string> atlasfilenames;
  for( int i = 0; i < tree_size; i++ )
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
  std::cout << "Write to output file..." << std::endl;
  itk::MABMISAtlasXMLFileWriter::Pointer atlasWriter = itk::MABMISAtlasXMLFileWriter::New();
  itk::MABMISAtlas                       atlas;
  atlas.m_NumberAllAtlases = totalAtlasSize;
  atlas.m_NumberSimulatedAtlases = simulatedAtlasSize;
  atlas.m_AtlasDirectory = std::string(".") + FILESEP + outputxmlname;
  atlas.m_Tree.resize(tree_size);
  atlas.m_IsSimulatedImage.resize(tree_size);
  atlas.m_RealImageIDs.resize(totalAtlasSize - simulatedAtlasSize);
  atlas.m_SimulatedImageIDs.resize(simulatedAtlasSize);
  for( int i = 0; i < tree_size; i++ )
    {
    atlas.m_Tree[i] = tree[i];
    }
  atlas.m_AtlasFilenames.resize(atlasfilenames.size() );
  int countRealImages = 0, countSimulatedImages = 0;
  for( int i = 0; i < atlasfilenames.size(); i++ )
    {
    size_t      sep = atlasfilenames[i].find_last_of(FILESEP);
    std::string fname = atlasfilenames[i];
    if( sep != string::npos )
      {
      fname = atlasfilenames[i].substr(sep + 1, string::npos);
/*			fname = atlasfilenames[i].substr(0, sep);
      sep = fname.find_last_of(FILESEP);
      if (sep != string::npos)
        fname = atlasfilenames[i].substr(sep+1, string::npos);
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
  for( int i = 0; i < segmentFiles.size(); i++ )
    {
    size_t      sep = segmentFiles[i].find_last_of(FILESEP);
    std::string fname = segmentFiles[i];
    if( sep != string::npos )
      {
      fname = segmentFiles[i].substr(sep + 1, string::npos);
/*			fname = segmentFiles[i].substr(0, sep);
      sep = fname.find_last_of(FILESEP);
      if (sep != string::npos)
        fname = segmentFiles[i].substr(sep+1, string::npos);
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

template <class T>
int DoIt( itk::MABMISImageData* trainingData, std::string outputFile,
          std::vector<int> iterations, double sigma)
{
  return Training(trainingData, outputFile, iterations, sigma);
}

int main( int argc, char *argv[] )
{
  PARSE_ARGS;

  // step 1: read the training atlases
  itk::MABMISImageDataXMLFileReader::Pointer trainingXMLReader = itk::MABMISImageDataXMLFileReader::New();
  trainingXMLReader->SetFilename(TrainingDataXML);
  try
    {
    trainingXMLReader->GenerateOutputInformation();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while reading file" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
  itk::MABMISImageData * trainingData = trainingXMLReader->GetOutputObject();

  // if the data path is empty, use the path of the xml file instead
  if( trainingData->m_DataDirectory.size() <= 1 )
    {
    size_t sep = TrainingDataXML.find_last_of(FILESEP);
    if( sep != std::string::npos )
      {
      trainingData->m_DataDirectory = TrainingDataXML.substr(0, sep);
      }
    else
      {
      trainingData->m_DataDirectory.resize(0);
      }
    }
  if( !trainingData->m_DataDirectory.empty() )
    {
    if( !(trainingData->m_DataDirectory[trainingData->m_DataDirectory.size() - 1] == FILESEP) )
      {
      trainingData->m_DataDirectory = trainingData->m_DataDirectory + FILESEP;
      }
    }

  // Will look into getting rid of these global variables later ---Xiaofeng
  datasimulator = DataSimulatorType::New();
  imgoperator = ImageOperationFilterType::New();
  dfoperator = DeformationFieldOperationFilterType::New();
  regoperator = ImageRegistrationFilterType::New();
  treeoperator = TreeOperationType::New();
  basicoperator = BasicOperationFilterType::New();

  // get number of iterations

  int retVal = DoIt<unsigned short>(trainingData, TrainingOutputFile, iterations, SmoothingKernelSize);
  delete trainingData;

  return retVal;

/*
  std::fstream pFile;
  std::string buffer;

  pFile.open (SubjectList.c_str());
  if (pFile == NULL) {
    perror ("Error opening file");
  }
  else  {
    while ( !pFile.eof() )
    {
      getline (pFile, buffer,'\n');
      strtrim( buffer );
      if (buffer.size() > 0)
      {
        sub_ids.push_back(buffer);
      }
    }
  }
  pFile.close();



  filenumber =  sub_ids.size();
  atlas_size = AtlasSize;

  if ((atlas_size<=0) || (atlas_size>=filenumber)){
    std::cerr << "atlas_size should be within [1," << filenumber-1 << "]!" << std::endl;
    return 102; //EXIT_FAILURE;
  }

  sample_size = filenumber - atlas_size;
  simulate_size = 2*atlas_size;

  int numAllCombinations = (int)pow((float)numSampleEachDirection, numEigenVector);
  if (numAllCombinations < simulate_size) {
        simulate_size = numAllCombinations;
  }

  allatlasnumber = atlas_size + simulate_size;
  allfilenumber = atlas_size + sample_size + simulate_size;

  std::string imagefilename;
  std::string segmentfilename;

  for (int i = 0; i< filenumber; i++) {
    imagefilename = sub_ids[i] + "_cbq_" + "000" + ".hdr";
    imgfilenames.push_back(imagefilename);
    segmentfilename = sub_ids[i] + "_seg_" + "000" + ".hdr";
    segfilenames.push_back(segmentfilename);
  }


  datasimulator = DataSimulatorType::New();
  imgoperator = ImageOperationFilterType::New();
  dfoperator = DeformationFieldOperationFilterType::New();
  regoperator = ImageRegistrationFilterType::New();
  treeoperator = TreeOperationType::New();
  basicoperator = BasicOperationFilterType::New();

  // read the size of images
  ReadImgInfo(imgfilenames);

  datasimulator->SetImx(imx);
  datasimulator->SetImy(imy);
  datasimulator->SetImz(imz);

  return DoIt<unsigned short>( argc, argv );
*/
  return EXIT_SUCCESS;
}

void SearchRootAmongAtlases(std::vector<std::string> imgfilenames, int & root)
{
  int tree_size = imgfilenames.size();

  vnl_matrix<double> distanceMatrix(tree_size, tree_size);
  for( int i = 0; i < tree_size; i++ )
    {
    for( int j = 0; j < tree_size; j++ )
      {
      distanceMatrix[i][j] = 0.0;
      }
    }

  vnl_vector<int> tree(tree_size);  // each tree

  // fill distance matrix
  imgoperator->PairwiseDistanceAmongImages(imgfilenames, tree_size, distanceMatrix);
  if( isDebug )
    {
    string distFileName = "dist_all_atlases.txt";
    basicoperator->SaveMatrix2File(distanceMatrix, tree_size, tree_size, distFileName);
    }

  // build MST
  tree = treeoperator->generateMSTFromMatrix(distanceMatrix, tree_size, tree);

  treeoperator->FindRoot(tree, tree_size, root);

  // save distance matrix and tree
  if( isDebug )
    {
    string treeFileName = "tree_all_atlases.txt";
    treeoperator->SaveTreeWithInfo(tree, tree_size, treeFileName);
    }

  std::cout << "The root is " << root << ": " << imgfilenames[root] << ". " << std::endl;
}

int RegistrationBetweenRootandAtlases(int root, std::vector<std::string> imageFileNames,
                                      std::vector<int> iterations, double sigma)
{
  // get the path from the filename
  size_t      sep = imageFileNames[0].find_last_of(FILESEP);
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
  std::cout << "Register between root and atlases ... ";
  for( int i = 0; i < atlas_size; i++ )
    {
    std::cout << i << ", ";
    if( i == root )
      {
      continue;
      }
    // if not the root image do registration
    string fixedImageFileName = imageFileNames[root];
    string movingImageFileName = imageFileNames[i];

    char i_str[10], root_str[10];
    sprintf(i_str, "%03d", i);
    sprintf(root_str, "%03d", root);
    std::string deformedImageFileName = std::string(i_str) + "_to_" + std::string(root_str) + "_cbq_reg.hdr";
    string      deformationFieldFileName = std::string(i_str) + "_to_" + std::string(root_str) + "_deform_000.mha";

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

    std::string invDeformedImageFileName = std::string(root_str) + "_to_" + std::string(i_str) + "_cbq_reg.hdr";
    string      invDeformationFileName = std::string(root_str) + "_to_" + std::string(i_str) + "_deform_000.mha";
    invDeformedImageFileName = outputFolder + invDeformedImageFileName;
    invDeformationFileName = outputFolder + invDeformationFileName;

    std::string invDeformedSegmentFileName = std::string(root_str) + "_to_" + std::string(i_str) + "seg_000.hdr";
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
      dfoperator->ApplyDeformationFieldAndWriteWithFileNames(
        fixedImageFileName,
        invDeformationFileName,
        invDeformedSegmentFileName, false);
      }
    }
  std::cout << "Done. " << std::endl;

  return 0;
}

int BuildStatisticalDeformationModel(int root,  std::vector<std::string> imageFileNames, int simulatedAtlasSize)
{
  ///////////////////////////////
  // do PCA simulation
  // std::cout << "Build deformation field model ... ";
  vector<string> allDeformationFieldFileNames;

  // get the path from the filename
  size_t      sep = imageFileNames[0].find_last_of(FILESEP);
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
  for( int i = 0; i < atlas_size; i++ )
    {
    string allDeformationFieldFileName;
    if( i == root )
      {
      continue;
      }

    char i_str[10], root_str[10];
    sprintf(i_str, "%03d", i);
    sprintf(root_str, "%03d", root);

    std::string deformationFieldFileName = std::string(i_str) + "_to_" + std::string(root_str) + "_deform_000.mha";
    deformationFieldFileName = outputFolder + deformationFieldFileName;
    allDeformationFieldFileNames.push_back(deformationFieldFileName);
    }

  return datasimulator->DoPCATraining(allDeformationFieldFileNames, atlas_size - 1, imageFileNames, root);
}

// return: simulated template file names
std::vector<std::string> GenerateSimulatedData(int root, std::vector<std::string> imageFiles, int simulatedAtlasSize)
{
  // output folder
  size_t      sep = imageFiles[0].find_last_of(FILESEP);
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
  std::cout << "The simulated images: " << std::endl;
  for( int i = 0; i < simulatedAtlasSize; i++ )
    {
    std::cout << i << ", ";

    string index_string;
    basicoperator->myitoa( i, index_string, 3 );

    string simulateDeformationFieldFileName = "simulated_deform_" + index_string + ".mha";
    string simulateTemplateFileName = "simulated_cbq_" + index_string + ".hdr";
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
    InternalImageType::Pointer simImg = 0;
    dfoperator->ApplyDeformationField(rootImg, invdf, simImg, true);
    imgoperator->WriteImage(simulateTemplateFileNames[i], simImg);
    }

  return simulateTemplateFileNames;
}

int ReadImgInfo(std::vector<std::string> imgfilenames)
{
  InternalImageType::Pointer inputImage1 = 0;

  return imgoperator->ReadImage(imgfilenames[0], inputImage1);
  // InternalImageType::SizeType input_size = inputImage1->GetLargestPossibleRegion().GetSize();
  // imx = input_size[0]; imy = input_size[1]; imz = input_size[2];
}
