#include "IGR3D_MABMIS_Direct_InvokeCLP.h"

#include "Testing.h"
#include "Training.h"

#include "itkAffineTransform.h"
#include "itkTransformFileReader.h"
#include "itkResampleImageFilter.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkPluginUtilities.h"

const unsigned int Dimension = 3;
typedef itk::Image<short, Dimension> ShortImageType;
typedef itk::AffineTransform<double, Dimension> AffineType;
typedef itk::Transform<double, Dimension, Dimension> TransformType;

template<typename PixelType>
void resampleAndWrite(const std::string & inFile, const std::string & outFile,
    ShortImageType::RegionType region,
    ShortImageType::PointType origin,
    ShortImageType::SpacingType spacing,
    ShortImageType::DirectionType direction,
    TransformType::Pointer transform)
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
  //resampleFilter->SetInterpolator(itk::WindowedSincInterpolateImageFunction<ImageType, 5>::New());
  resampleFilter->SetInterpolator(itk::BSplineInterpolateImageFunction<ImageType>::New()); //cubic by default
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
    itk::ImageIOBase::IOComponentType componentType)
{
switch (componentType)
  {
  case itk::ImageIOBase::UCHAR:
    resampleAndWrite<unsigned char>(inFile, outFile, region, origin, spacing, direction, transform);
    break;
  case itk::ImageIOBase::CHAR:
    resampleAndWrite<char>(inFile, outFile, region, origin, spacing, direction, transform);
    break;
  case itk::ImageIOBase::USHORT:
    resampleAndWrite<unsigned short>(inFile, outFile, region, origin, spacing, direction, transform);
    break;
  case itk::ImageIOBase::SHORT:
    resampleAndWrite<short>(inFile, outFile, region, origin, spacing, direction, transform);
    break;
  case itk::ImageIOBase::UINT:
    resampleAndWrite<unsigned int>(inFile, outFile, region, origin, spacing, direction, transform);
    break;
  case itk::ImageIOBase::INT:
    resampleAndWrite<int>(inFile, outFile, region, origin, spacing, direction, transform);
    break;
  case itk::ImageIOBase::ULONG:
    resampleAndWrite<unsigned long>(inFile, outFile, region, origin, spacing, direction, transform);
    break;
  case itk::ImageIOBase::LONG:
    resampleAndWrite<long>(inFile, outFile, region, origin, spacing, direction, transform);
    break;
  //case itk::ImageIOBase::ULONGLONG:
  //  resampleAndWrite<unsigned long long>(inFile, outFile, region, origin, spacing, direction, transform);
  //  break;
  //case itk::ImageIOBase::LONGLONG:
  //  resampleAndWrite<long long>(inFile, outFile, region, origin, spacing, direction, transform);
  //  break;
  case itk::ImageIOBase::FLOAT:
    resampleAndWrite<float>(inFile, outFile, region, origin, spacing, direction, transform);
    break;
  case itk::ImageIOBase::DOUBLE:
    resampleAndWrite<double>(inFile, outFile, region, origin, spacing, direction, transform);
    break;
  case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
  default:
    itkGenericExceptionMacro("Error: unknown component type in image " << inFile);
    break;
  }
}

#define HANDLE_IMAGE(n)                           \
{                                                 \
  if (image##n##Arg.isSet())                      \
    {                                             \
    imageFileNames.push_back(image##n);           \
    if (transform##n##Arg.isSet())                \
      {                                           \
      transformFileNames.push_back(transform##n); \
      }                                           \
    else                                          \
      {                                           \
      transformFileNames.push_back("");           \
      }                                           \
    }                                             \
}

int main( int argc, char *argv[] )
{
  PARSE_ARGS;

  imageListXML = ReplacePathSepForOS(imageListXML);
  atlasTreeXML = ReplacePathSepForOS(atlasTreeXML);
  imageDir = ReplacePathSepForOS(imageDir);

  try
    {
    itk::MABMISImageData miData;
    const std::string extension = ".nrrd";
    std::cout << "Examining inputs" << std::endl;
    if (!imageListXMLArg.isSet())
      {
      imageListXML = ".";
      }
    std::string listDir = itksys::SystemTools::GetParentDirectory(imageListXML);
    listDir = itksys::SystemTools::GetRealPath(listDir);
    imageDir = itksys::SystemTools::GetRealPath(imageDir);
    if (imageDirArg.isSet() && !imageDir.empty() && imageDir != "." && imageDir != listDir)
      {
      miData.m_DataDirectory = imageDir;
      miData.m_OutputDirectory = imageDir;
      }
    else
      {
      imageDir = listDir;
      }

    //check which parameters are present
    std::vector<std::string> imageFileNames;
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
    std::vector<TransformType::Pointer > transforms;
    AffineType::Pointer identity = AffineType::New();
    itk::TransformFileReader::Pointer transformReader = itk::TransformFileReader::New();
    for (unsigned i = 0; i < imageFileNames.size(); i++)
      {
      miData.m_ImageFileNames.push_back("image" + std::to_string(i) + extension);
      miData.m_SegmentationFileNames.push_back("image" + std::to_string(i) + "-label" + extension);
      if (transformFileNames[i].empty())
        {
        transforms.push_back(identity.GetPointer());
        }
      else
        {
        transformReader->SetFileName(transformFileNames[i]);
        transformReader->Update();
        if (transformReader->GetTransformList()->size() > 1)
          {
          itkGenericExceptionMacro("Only simple affine transforms are supported. "
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
    miData.m_NumberImageData = imageFileNames.size();
    if (miData.m_NumberImageData == 0)
      {
      itkGenericExceptionMacro("Some images have to be provided");
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

    itk::MABMISAtlas* atlas = nullptr;
    std::string rootFilename;
    if (atlasTreeXMLArg.isSet())
      {
      //read atlas tree
      itk::MABMISAtlasXMLFileReader::Pointer atlasReader =
          itk::MABMISAtlasXMLFileReader::New();
      atlasReader->SetFilename(atlasTreeXML);
      atlasReader->GenerateOutputInformation();
      atlas = atlasReader->GetOutputObject();
      
      //get atlas root image
      unsigned i = 0;
      for (unsigned i = 0; i < atlas->m_NumberAllAtlases; i++)
        {
        if (atlas->m_Tree[i] == atlas->m_TreeRoot)
          {
          rootFilename = atlas->m_AtlasFilenames[i];
          break;
          }
        }
      
      if (atlas->m_AtlasDirectory[0] == '.') //relative path
        {
        rootFilename = itksys::SystemTools::GetParentDirectory(atlasTreeXML)
            + atlas->m_AtlasDirectory.substr(1) + '/' + rootFilename;
        }
      else
        {
        rootFilename = atlas->m_AtlasDirectory + '/' + rootFilename;
        }
      }
    else
      {
      rootFilename = imageFileNames[0];
      }

    typedef itk::ImageFileReader<ShortImageType> ShortReaderType;
    ShortReaderType::Pointer shortReader = ShortReaderType::New();
    shortReader->SetFileName(rootFilename);
    shortReader->UpdateOutputInformation();

    //atlas root metadata
    ShortImageType::RegionType region = shortReader->GetOutput()->GetLargestPossibleRegion();
    ShortImageType::PointType origin = shortReader->GetOutput()->GetOrigin();
    ShortImageType::SpacingType spacing = shortReader->GetOutput()->GetSpacing();
    ShortImageType::DirectionType direction = shortReader->GetOutput()->GetDirection();

    //transform the images to match the atlas root and write them
    for (unsigned i = 0; i < imageFileNames.size(); i++)
      {
      itk::ImageIOBase::IOPixelType pixelType;
      itk::ImageIOBase::IOComponentType componentType;
      itk::GetImageType(imageFileNames[i], pixelType, componentType);

      shortReader->SetFileName(imageFileNames[i]);
      shortReader->UpdateOutputInformation();

      if (region != shortReader->GetOutput()->GetLargestPossibleRegion()
          || origin != shortReader->GetOutput()->GetOrigin()
          || spacing != shortReader->GetOutput()->GetSpacing()
          || direction != shortReader->GetOutput()->GetDirection()
          || GetExtension(imageFileNames[i]) != GetExtension(miData.m_ImageFileNames[i]))
        {
        std::cout << "Resampling " << imageFileNames[i] << std::endl;
        resampleAndWrite(imageFileNames[i], imageDir + '/' + miData.m_ImageFileNames[i],
            region, origin, spacing, direction, transforms[i], componentType);      
        }
      else
        {
        std::cout << "Copying " << imageFileNames[i] << std::endl;
        itksys::SystemTools::CopyFileAlways(imageFileNames[i], imageDir + '/' + miData.m_ImageFileNames[i]);
        }
      }

    //write the XML or invoke testing
    if (mode == "Create imageXML")
      {
      itk::MABMISImageDataXMLFileWriter::Pointer miWriter =
          itk::MABMISImageDataXMLFileWriter::New();
      miWriter->SetObject(&miData);
      miWriter->SetFilename(imageListXML);
      miWriter->WriteFile();
      }
    else if (mode == "Direct invoke")
      {
      Testing(&miData, atlas, iterations, sigma);
      }
    else if (mode == "(re)train atlas")
      {
      if (atlas) //merge images from the atlas into the list
        {
        for (unsigned i = 0; i < atlas->m_NumberAllAtlases; i++)
          {
          miData.m_ImageFileNames.push_back(atlas->m_AtlasFilenames[i]);
          miData.m_SegmentationFileNames.push_back(atlas->m_AtlasSegmentationFilenames[i]);
          }
        miData.m_NumberImageData = miData.m_ImageFileNames.size();
        miData.m_DataDirectory = atlas->m_AtlasDirectory; // needed?
        }

      Training(&miData, atlasTreeXML, iterations, sigma);
      }
    else
      {
      itkGenericExceptionMacro("Unknown processing mode");
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
