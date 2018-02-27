#include "IGR3D_MABMIS_Direct_InvokeCLP.h"

#include "Testing.h"
#include "Training.h"

#include "itkAffineTransform.h"
#include "itkTransformFileReader.h"
#include "itkResampleImageFilter.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLabelImageGaussianInterpolateImageFunction.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkPluginUtilities.h"

const unsigned int Dimension = 3;
typedef itk::Image<short, Dimension> ShortImageType;
typedef itk::AffineTransform<double, Dimension> AffineType;
typedef itk::Transform<double, Dimension, Dimension> TransformType;
AffineType::Pointer identity = AffineType::New();

// interpolationQuality:
// -1 = linear (default and fall-back)
//  0 = nearest neighbor
//  1 = LabelImageGaussianInterpolate (smoother alternative to nearest neighbor)
//  3 = cubic BSpline
//  5 = WindowedSinc with radius of 5 voxels
template<typename PixelType>
void resampleAndWrite(const std::string & inFile, const std::string & outFile,
    ShortImageType::RegionType region,
    ShortImageType::PointType origin,
    ShortImageType::SpacingType spacing,
    ShortImageType::DirectionType direction,
    TransformType::Pointer transform,
    int interpolationQuality)
{
  typedef itk::Image<PixelType, Dimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inFile);
  reader->Update();

  typedef itk::MinimumMaximumImageCalculator<ImageType> MinMaxType;
  typename MinMaxType::Pointer minMax = MinMaxType::New();
  minMax->SetImage(reader->GetOutput());
  minMax->ComputeMinimum();

  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  typename ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
  switch (interpolationQuality)
    {
    case 0:
      resampleFilter->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<ImageType>::New());
      break;
    case 1:
      resampleFilter->SetInterpolator(itk::LabelImageGaussianInterpolateImageFunction<ImageType>::New());
      break;
    case 3:
      resampleFilter->SetInterpolator(itk::BSplineInterpolateImageFunction<ImageType>::New()); //order 3 (=cubic) by default
      break;
    case 5:
      resampleFilter->SetInterpolator(itk::WindowedSincInterpolateImageFunction<ImageType, 5>::New());
      break;
    case -1:
    default:
      // ResampleImageFilter uses linear by default
      break;
    }
  resampleFilter->SetInput(reader->GetOutput());
  resampleFilter->SetSize(region.GetSize());
  resampleFilter->SetOutputOrigin(origin);
  resampleFilter->SetOutputSpacing(spacing);
  resampleFilter->SetOutputDirection(direction);
  resampleFilter->SetTransform(transform);
  resampleFilter->SetDefaultPixelValue(minMax->GetMinimum());
  resampleFilter->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outFile);
  writer->SetInput(resampleFilter->GetOutput());
  writer->SetUseCompression(true);
  writer->Update();
}

//does the template-based dispatch
void resampleAndWrite(const std::string & inFile, const std::string & outFile,
    ShortImageType::RegionType region,
    ShortImageType::PointType origin,
    ShortImageType::SpacingType spacing,
    ShortImageType::DirectionType direction,
    TransformType::Pointer transform,
    int interpolationQuality,
    itk::ImageIOBase::IOComponentType componentType)
{
switch (componentType)
  {
  case itk::ImageIOBase::UCHAR:
    resampleAndWrite<unsigned char>(inFile, outFile, region, origin, spacing, direction, transform, interpolationQuality);
    break;
  case itk::ImageIOBase::CHAR:
    resampleAndWrite<char>(inFile, outFile, region, origin, spacing, direction, transform, interpolationQuality);
    break;
  case itk::ImageIOBase::USHORT:
    resampleAndWrite<unsigned short>(inFile, outFile, region, origin, spacing, direction, transform, interpolationQuality);
    break;
  case itk::ImageIOBase::SHORT:
    resampleAndWrite<short>(inFile, outFile, region, origin, spacing, direction, transform, interpolationQuality);
    break;
  case itk::ImageIOBase::UINT:
    resampleAndWrite<unsigned int>(inFile, outFile, region, origin, spacing, direction, transform, interpolationQuality);
    break;
  case itk::ImageIOBase::INT:
    resampleAndWrite<int>(inFile, outFile, region, origin, spacing, direction, transform, interpolationQuality);
    break;
  case itk::ImageIOBase::ULONG:
    resampleAndWrite<unsigned long>(inFile, outFile, region, origin, spacing, direction, transform, interpolationQuality);
    break;
  case itk::ImageIOBase::LONG:
    resampleAndWrite<long>(inFile, outFile, region, origin, spacing, direction, transform, interpolationQuality);
    break;
  //case itk::ImageIOBase::ULONGLONG:
  //  resampleAndWrite<unsigned long long>(inFile, outFile, region, origin, spacing, direction, transform, interpolationQuality);
  //  break;
  //case itk::ImageIOBase::LONGLONG:
  //  resampleAndWrite<long long>(inFile, outFile, region, origin, spacing, direction, transform, interpolationQuality);
  //  break;
  case itk::ImageIOBase::FLOAT:
    resampleAndWrite<float>(inFile, outFile, region, origin, spacing, direction, transform, interpolationQuality);
    break;
  case itk::ImageIOBase::DOUBLE:
    resampleAndWrite<double>(inFile, outFile, region, origin, spacing, direction, transform, interpolationQuality);
    break;
  case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
  default:
    itkGenericExceptionMacro("Error: unknown component type in image " << inFile);
    break;
  }
}

// resample inFile into outFile so it matches the image grid of refFile
// uses the pixel type of the inFile
void resampleOrCopy(const std::string & inFile,
    const std::string & outFile,
    const std::string & refFile,
    TransformType::Pointer transform,
    int interpolationQuality)
{
  typedef itk::ImageFileReader<ShortImageType> ShortReaderType;
  ShortReaderType::Pointer shortReader = ShortReaderType::New();
  shortReader->SetFileName(refFile);
  shortReader->UpdateOutputInformation();
  ShortImageType::Pointer image = shortReader->GetOutput();

  //reference image metadata
  ShortImageType::RegionType region = image->GetLargestPossibleRegion();
  ShortImageType::PointType origin = image->GetOrigin();
  ShortImageType::SpacingType spacing = image->GetSpacing();
  ShortImageType::DirectionType direction = image->GetDirection();

  //inImage componenet type
  itk::ImageIOBase::IOPixelType pixelType;
  itk::ImageIOBase::IOComponentType componentType;
  itk::GetImageType(inFile, pixelType, componentType);

  shortReader->SetFileName(inFile);
  shortReader->UpdateOutputInformation();
  image = shortReader->GetOutput();

  bool transformIsIdentity = false;
  AffineType * aTransform = dynamic_cast<AffineType *>(transform.GetPointer());
  if (aTransform && aTransform->Metric(identity) == 0.0)
    {
    transformIsIdentity = true;
    }

  if (region != image->GetLargestPossibleRegion()
      || origin != image->GetOrigin()
      || spacing != image->GetSpacing()
      || direction != image->GetDirection()
      || transformIsIdentity != true
      || GetExtension(inFile) != GetExtension(outFile))
    {
    std::cout << "Resampling " << inFile << std::endl;
    std::cout << "      into " << outFile << std::endl;
    resampleAndWrite(inFile, outFile,
        region, origin, spacing, direction, transform,
        interpolationQuality, componentType);
    }
  else
    {
    std::cout << "Copying " << inFile << std::endl;
    std::cout << "     to " << outFile << std::endl;
    itksys::SystemTools::CopyFileAlways(inFile, outFile);
    }
}

//transform the images to match the atlas root and write them
void resampleOrCopyN(const std::string & rootFilename, const std::string & outDir,
    const std::vector<std::string> & inFiles,
    const std::vector<std::string> & outFiles,
    std::vector<TransformType::Pointer > & transforms,
    int interpolationQuality)
{
  for (unsigned i = 0; i < inFiles.size(); i++)
    {
    resampleOrCopy(inFiles[i], outDir + outFiles[i], rootFilename, transforms[i], interpolationQuality);
    }
}

#define HANDLE_IMAGE(n)                                 \
{                                                       \
  if (image##n##Arg.isSet())                            \
    {                                                   \
    imageFileNames.push_back(image##n);                 \
    if (transform##n##Arg.isSet())                      \
      {                                                 \
      transformFileNames.push_back(transform##n);       \
      }                                                 \
    else                                                \
      {                                                 \
      transformFileNames.push_back("");                 \
      }                                                 \
    if (segmentation##n##Arg.isSet())                   \
      {                                                 \
      segmentationFileNames.push_back(segmentation##n); \
      }                                                 \
    else                                                \
      {                                                 \
      segmentationFileNames.push_back("");              \
      }                                                 \
    }                                                   \
}
//segmentationFileNames

int main( int argc, char *argv[] )
{
  PARSE_ARGS;

  int segmentationInterpolationQuality = 1;
  if (segmentationInterpolationArg.isSet() && segmentationInterpolation == "Nearest")
    {
    segmentationInterpolationQuality = 0;
    }

  int imageInterpolationQuality = 3;
  if (imageInterpolationArg.isSet() && imageInterpolation == "Linear")
    {
    imageInterpolationQuality = -1;
    }
  if (imageInterpolationArg.isSet() && imageInterpolation == "WindowedSinc")
    {
    imageInterpolationQuality = 5;
    }

  imageListXML = ReplacePathSepForOS(imageListXML);
  atlasTreeXML = ReplacePathSepForOS(atlasTreeXML);
  imageDir = ReplacePathSepForOS(imageDir);

  try
    {
    itk::MABMISImageData miData;
    const std::string extension = ".nrrd";
    std::string listDir = itksys::SystemTools::GetParentDirectory(imageListXML);
    listDir = itksys::SystemTools::GetRealPath(listDir);
    imageDir = itksys::SystemTools::GetRealPath(imageDir);
    if (!imageDirArg.isSet() || imageDir.empty() || imageDir == ".")
      {
      imageDir = listDir;
      }
    imageDir += '/';

    std::vector<TransformType::Pointer > transforms;
    std::vector<std::string> imageFileNames;
    std::vector<std::string> segmentationFileNames;
    itk::MABMISAtlas* atlas = nullptr;
    std::string rootFilename;
    if (atlasTreeXMLArg.isSet() && itksys::SystemTools::FileExists(atlasTreeXML))
      {
      //read atlas tree
      itk::MABMISAtlasXMLFileReader::Pointer atlasReader =
          itk::MABMISAtlasXMLFileReader::New();
      atlasReader->SetFilename(atlasTreeXML);
      atlasReader->GenerateOutputInformation();
      atlas = atlasReader->GetOutputObject();

      //get atlas root image
      for (unsigned i = 0; i < atlas->m_NumberAllAtlases; i++)
        {
        if (atlas->m_Tree[i] == atlas->m_TreeRoot)
          {
          rootFilename = atlas->m_AtlasFilenames[i];
          break;
          }
        }

      const std::string& atlDir = atlas->m_AtlasDirectory;
      if (atlDir.size() > 2 && atlDir[0] == '.' && (atlDir[1] == '\\' || atlDir[1] == '/'))
        {
        listDir = itksys::SystemTools::GetParentDirectory(atlasTreeXML)
              + atlas->m_AtlasDirectory.substr(1) + '/';
        }
      else
        {
        listDir = atlas->m_AtlasDirectory + '/';
        }
      atlas->m_AtlasDirectory = listDir;
      rootFilename = listDir + rootFilename;

      if (mode == "(re)train atlas")
        {
        std::cout << "Re-training the atlas: " << atlasTreeXML << std::endl;
        imageDir = listDir; //new images should be put into the old atlas directory too
        for (unsigned i = 0; i < atlas->m_NumberAllAtlases - atlas->m_NumberSimulatedAtlases; i++)
          {
          miData.m_ImageFileNames.push_back(atlas->m_AtlasFilenames[i]);
          imageFileNames.push_back(imageDir + atlas->m_AtlasFilenames[i]);
          miData.m_SegmentationFileNames.push_back(atlas->m_AtlasSegmentationFilenames[i]);
          segmentationFileNames.push_back(imageDir + atlas->m_AtlasSegmentationFilenames[i]);
          transforms.push_back(identity.GetPointer());
          }
        }
      }
    else
      {
      if (mode == "Direct invoke")
        {
        itkGenericExceptionMacro("Atlas XML file is missing " << atlasTreeXML);
        }
      else if (mode == "(re)train atlas")
        {
        std::cout << "Initial creation of the atlas: " << atlasTreeXML << std::endl;
        }
      }
    unsigned numberOfPreExistingImages = miData.m_ImageFileNames.size();
    std::cout << "Using " << imageDir << " as a temporary folder." << std::endl;
    itksys::SystemTools::MakeDirectory(imageDir);

    //check which parameters are present
    std::vector<std::string> transformFileNames;
    HANDLE_IMAGE(0);
    HANDLE_IMAGE(1);
    HANDLE_IMAGE(2);
    HANDLE_IMAGE(3);
    HANDLE_IMAGE(4);
    HANDLE_IMAGE(5);
    HANDLE_IMAGE(6);
    HANDLE_IMAGE(7);
    HANDLE_IMAGE(8);
    HANDLE_IMAGE(9);

    //fill miData and read array of transforms
    itk::TransformFileReader::Pointer transformReader = itk::TransformFileReader::New();
    for (unsigned i = numberOfPreExistingImages; i < imageFileNames.size(); i++)
      {
      miData.m_ImageFileNames.push_back("image" + std::to_string(i) + extension);
      miData.m_SegmentationFileNames.push_back("image" + std::to_string(i) + std::string("-label") + extension);
      if (transformFileNames[i - numberOfPreExistingImages].empty())
        {
        transforms.push_back(identity.GetPointer());
        }
      else
        {
        transformReader->SetFileName(transformFileNames[i - numberOfPreExistingImages]);
        transformReader->Update();
        if (transformReader->GetTransformList()->size() > 1)
          {
          itkGenericExceptionMacro("Only simple transforms are supported. "
              << transformFileNames[i] << " contains more than one transform!");
          }
        itk::TransformBaseTemplate<double>::Pointer genericTransform = transformReader->GetTransformList()->front();
        TransformType::Pointer aTransform = static_cast<TransformType *>(genericTransform.GetPointer());
        if (aTransform.IsNull())
          {
          itkGenericExceptionMacro("Could not interpret " << transformFileNames[i] << " as a supported transform!");
          }
        transforms.push_back(aTransform);
        }
      }
    miData.m_NumberImageData = miData.m_ImageFileNames.size();
    if (miData.m_NumberImageData == 0 || (atlas && atlas->m_AtlasFilenames.size() == miData.m_NumberImageData))
      {
      itkGenericExceptionMacro("Some images have to be provided");
      }
    if (rootFilename.empty()) //there was no atlas
      {
      rootFilename = imageFileNames[0];
      }

    if (mode == "(re)train atlas")
      {
      //check if all images have segmentations
      for (unsigned i = 0; i < miData.m_SegmentationFileNames.size(); i++)
        {
        if (miData.m_SegmentationFileNames[i].empty())
          {
          itkGenericExceptionMacro(<< miData.m_ImageFileNames[i] << " lacks a segmentation!");
          }
        }
      }

    resampleOrCopyN(rootFilename, imageDir, imageFileNames, miData.m_ImageFileNames, transforms, imageInterpolationQuality);
    if (mode == "Create imageXML" || imageListXMLArg.isSet())
      {
      itk::MABMISImageDataXMLFileWriter::Pointer miWriter =
          itk::MABMISImageDataXMLFileWriter::New();
      miWriter->SetObject(&miData);
      miWriter->SetFilename(imageListXML);
      miWriter->WriteFile();
      }
    //set path after the xml file is written, so xml uses relative paths
    miData.m_DataDirectory = imageDir;
    miData.m_OutputDirectory = imageDir;

    if (mode == "Direct invoke")
      {
      Testing(&miData, atlas, iterations, sigma);

      //inverse the transforms and resample the segmentations back into original image grids
      for (unsigned i = numberOfPreExistingImages; i < segmentationFileNames.size(); i++)
        {
        std::string segDir = itksys::SystemTools::GetParentDirectory(segmentationFileNames[i]);
        itksys::SystemTools::MakeDirectory(segDir);

        TransformType::Pointer invT = transforms[i]->GetInverseTransform();
        if (!invT) //not possible to invert
          {
          invT = identity;
          std::cerr << "Error: it was not possible to invert tranform " << transformFileNames[i] << std::endl;
          std::cerr << "Resampling segmentation back using identity transform" << std::endl;
          }

        resampleOrCopy(imageDir + miData.m_SegmentationFileNames[i], segmentationFileNames[i],
            imageFileNames[i], invT, segmentationInterpolationQuality);
        }
      }
    if (mode == "(re)train atlas")
      {
      //segmentations are now inputs which accompany the images, and need to be transformed the same way
      resampleOrCopyN(rootFilename, imageDir, segmentationFileNames, miData.m_SegmentationFileNames, transforms, segmentationInterpolationQuality);
      Training(&miData, atlasTreeXML, iterations, sigma);
      }
    }
  catch (itk::ExceptionObject &exc)
    {
    std::cout << exc;
    return EXIT_FAILURE;
    }
  catch (std::runtime_error &exc)
    {
    std::cout << "Runtime error has occurred: " << exc.what() << std::endl;
    return EXIT_FAILURE;
    }
  catch (...)
    {
    std::cout << "An unknown error has occurred!" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS; //no error has occurred
}
