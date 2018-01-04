#include "itkMABMISImageOperationFilter.h"
#ifndef __itkMABMISImageOperationFilter_hxx
#define __itkMABMISImageOperationFilter_hxx

#include "itkMABMISImageOperationFilter.h"

namespace itk
{
namespace Statistics
{
template <class TInputImage, class TOutputImage>
MABMISImageOperationFilter<TInputImage, TOutputImage>
::MABMISImageOperationFilter()
{
}

template <class TInputImage, class TOutputImage>
MABMISImageOperationFilter<TInputImage, TOutputImage>
::~MABMISImageOperationFilter()
{
}

template <class TInputImage, class TOutputImage>
void
MABMISImageOperationFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}

template <class TInputImage, class TOutputImage>
int
MABMISImageOperationFilter<TInputImage, TOutputImage>
::ReadImage(const std::string filename, InternalImageType::Pointer& image)
{
  // std::cout << "Reading "<< filename << std:: endl;
  InternalImageReaderType::Pointer internalImageReader = InternalImageReaderType::New();
  internalImageReader->SetFileName(filename);
  try
    {
    internalImageReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << std::endl;
    return -1;
    }
  image = internalImageReader->GetOutput();
  return 0;
}

template<class TInputImage, class TOutputImage>
itk::ImageIOBase::IOComponentType
MABMISImageOperationFilter<TInputImage, TOutputImage>
::GetIOPixelType(const std::string & filename)
{
  itk::ImageIOBase::Pointer imageIO;

  try
    {
    imageIO = itk::ImageIOFactory::CreateImageIO(filename.c_str(), itk::ImageIOFactory::ReadMode);
    if( imageIO )
      {
      imageIO->SetFileName(filename);
      imageIO->ReadImageInformation();
      return imageIO->GetComponentType();
      }
    else
      {
      std::cout << "Could not create the imageIO for file " << filename << "." << std::endl;
      exit( EXIT_FAILURE );
      }
    }
  catch( itk::ExceptionObject& err )
    {
    std::cout << "Could not read the image information of " << filename << "." << std::endl;
    std::cout << err << std::endl;
    exit( EXIT_FAILURE );
    }

  return itk::ImageIOBase::UCHAR; //this should never be reached
}

template <class TInputImage, class TOutputImage>
template <class PixelType>
int
MABMISImageOperationFilter<TInputImage, TOutputImage>
::WriteImage(const std::string filename, InternalImageType::Pointer image)
{
  typedef itk::Image<PixelType, ImageDimension> GivenImageType;
  typedef itk::CastImageFilter<InternalImageType, GivenImageType> CasterType;
  typename CasterType::Pointer caster = CasterType::New();
  typedef itk::ImageFileWriter<GivenImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  caster->SetInput(image);
  writer->SetInput(caster->GetOutput());
  writer->SetFileName(filename);
  writer->SetUseCompression(false);
  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject & err)
    {
    std::cerr << err << std::endl;
    return -1;
    }
  return 0;
}

template<class TInputImage, class TOutputImage>
int
MABMISImageOperationFilter<TInputImage, TOutputImage>
::WriteImage(const std::string filename, InternalImageType::Pointer image, itk::ImageIOBase::IOComponentType ioType)
{
  if( ioType == itk::ImageIOBase::UCHAR )
    {
    return this->WriteImage<unsigned char>(filename, image);
    }
  else if( ioType == itk::ImageIOBase::CHAR )
    {
    return this->WriteImage<char>(filename, image);
    }
  else if( ioType == itk::ImageIOBase::USHORT )
    {
    return this->WriteImage<unsigned short>(filename, image);
    }
  else if( ioType == itk::ImageIOBase::SHORT )
    {
    return this->WriteImage<short>(filename, image);
    }
  else if( ioType == itk::ImageIOBase::UINT )
    {
    return this->WriteImage<unsigned int>(filename, image);
    }
  else if( ioType == itk::ImageIOBase::INT )
    {
    return this->WriteImage<int>(filename, image);
    }
  else if( ioType == itk::ImageIOBase::ULONG )
    {
    return this->WriteImage<unsigned long>(filename, image);
    }
  else if( ioType == itk::ImageIOBase::LONG )
    {
    return this->WriteImage<long>(filename, image);
    }
  //else if( ioType == itk::ImageIOBase::ULONGLONG )
  //  {
  //  return this->WriteImage<unsigned long long>(filename, image);
  //  }
  //else if( ioType == itk::ImageIOBase::LONGLONG )
  //  {
  //  return this->WriteImage<long long>(filename, image);
  //  }
  else if(ioType == itk::ImageIOBase::FLOAT)
    {
    return this->WriteImage<float>(filename, image);
    }
  else if(ioType == itk::ImageIOBase::DOUBLE)
    {
    return this->WriteImage<InternalPixelType>(filename, image);
    }
  else
    {
    std::cerr << itk::ImageIOBase::GetComponentTypeAsString(ioType) << " not supported as I/O pixel type" << std::endl;
    return -1;
    }
  return -2;
}

template <class TInputImage, class TOutputImage>
void
MABMISImageOperationFilter<TInputImage, TOutputImage>
::WriteImage(const std::string filename, InternalImageType::Pointer image, char* outputType)
{
  itk::ImageIOBase::IOComponentType ioType = itk::ImageIOBase::GetComponentTypeFromString(outputType);
  return this->WriteImage(filename, image, ioType);
}

// calculate distance between two images by mean squared difference
template <class TInputImage, class TOutputImage>
double
MABMISImageOperationFilter<TInputImage, TOutputImage>
::calculateDistanceMSD(const std::string& imageName1, const std::string&  imageName2)
{
  double dist;

  InternalImageReaderType::Pointer imageReader1  = InternalImageReaderType::New();
  InternalImageReaderType::Pointer imageReader2 = InternalImageReaderType::New();

  imageReader1->SetFileName( imageName1 );
  imageReader2->SetFileName( imageName2 );
  try
    {
    imageReader1->Update();
    imageReader2->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

  InternalImageType::Pointer image1 = imageReader1->GetOutput();
  InternalImageType::Pointer image2 = imageReader2->GetOutput();

  InternalImageType::SizeType size;
  size = image1->GetLargestPossibleRegion().GetSize();
  InternalImageIteratorType it1(image1, image1->GetLargestPossibleRegion() );
  InternalImageIteratorType it2(image2, image2->GetLargestPossibleRegion() );
  dist = 0.0;
  for( it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2 )
    {
    const float value1 = it1.Get();
    const float value2 = it2.Get();
    dist += (value1 - value2) * (value1 - value2);
    }
  dist = dist / (size[0] * size[1] * size[2]);

  return dist;
}

// calculate pair-wise distance in an image population
template <class TInputImage, class TOutputImage>
void
MABMISImageOperationFilter<TInputImage, TOutputImage>
::PairwiseDistanceAmongImages(const std::vector<std::string>& imageFileNames, int totalNumber, vnl_matrix<double>& distanceMatrix)

{
  for( int i = 0; i < totalNumber; ++i )
    {
    // std::cout << imageFileNames[i] << std::endl;
    for( int j = i; j < totalNumber; ++j )
      {
      if( j == i )
        {
        distanceMatrix[i][j] = 0.0;
        }
      else
        {
        const double dist = calculateDistanceMSD(imageFileNames[i], imageFileNames[j]);
        // cout << "dist:: " << dist << endl;
        distanceMatrix[i][j] = sqrt(dist);
        distanceMatrix[j][i] = sqrt(dist);
        }
      }
    }
  // std::cout << "done!" << endl;
  return;
}
} // namespace Statistics
} // namespace itk

#endif
