#ifndef __itkMABMISImageOperationFilter_txx
#define __itkMABMISImageOperationFilter_txx

#include "itkMABMISImageOperationFilter.h"

using namespace std;

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
::ReadImage(string filename, InternalImageType::Pointer& image)
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

template <class TInputImage, class TOutputImage>
void
MABMISImageOperationFilter<TInputImage, TOutputImage>
::WriteImageUCHAR(string filename, InternalImageType::Pointer image)
{
  Internal2CharCastFilterType::Pointer caster = Internal2CharCastFilterType::New();

  caster->SetInput(image);
  CharImageWriterType::Pointer writer = CharImageWriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(caster->GetOutput() );
  writer->SetUseCompression( false );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << std::endl;
    return;
    }
  return;
}

template <class TInputImage, class TOutputImage>
void
MABMISImageOperationFilter<TInputImage, TOutputImage>
::WriteImageINT(string filename, InternalImageType::Pointer image)
{
  Internal2IntCastFilterType::Pointer caster = Internal2IntCastFilterType::New();

  caster->SetInput(image);
  IntImageWriterType::Pointer writer = IntImageWriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(caster->GetOutput() );
  writer->SetUseCompression( false );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << std::endl;
    return;
    }
  return;
}

template <class TInputImage, class TOutputImage>
void
MABMISImageOperationFilter<TInputImage, TOutputImage>
::WriteImageSHORT(string filename, InternalImageType::Pointer image)
{
  Internal2ShortCastFilterType::Pointer caster = Internal2ShortCastFilterType::New();

  caster->SetInput(image);
  ShortImageWriterType::Pointer writer = ShortImageWriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(caster->GetOutput() );
  writer->SetUseCompression( false );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << std::endl;
    return;
    }
  return;
}

template <class TInputImage, class TOutputImage>
void
MABMISImageOperationFilter<TInputImage, TOutputImage>
::WriteImageFLOAT(string filename, InternalImageType::Pointer image)
{
  Internal2FloatCastFilterType::Pointer caster = Internal2FloatCastFilterType::New();

  caster->SetInput(image);
  FloatImageWriterType::Pointer writer = FloatImageWriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(caster->GetOutput() );
  writer->SetUseCompression( false );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << std::endl;
    return;
    }
  return;
}

template <class TInputImage, class TOutputImage>
void
MABMISImageOperationFilter<TInputImage, TOutputImage>
::WriteImage(string filename, InternalImageType::Pointer image)
{
  Internal2CharCastFilterType::Pointer caster = Internal2CharCastFilterType::New();

  caster->SetInput(image);
  CharImageWriterType::Pointer writer = CharImageWriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(caster->GetOutput() );
  writer->SetUseCompression( false );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << std::endl;
    return;
    }
  return;
}

template <class TInputImage, class TOutputImage>
void
MABMISImageOperationFilter<TInputImage, TOutputImage>
::WriteImage(string filename, InternalImageType::Pointer image, char* outputType)
{
  //
  if( strcmp(outputType, "uchar") == 0 )
    {
    WriteImageUCHAR(filename, image);
    }
  else if( strcmp(outputType, "short") == 0 )
    {
    WriteImageSHORT(filename, image);
    }
  else if( strcmp(outputType, "int") == 0 )
    {
    WriteImageINT(filename, image);
    }
  else if( strcmp(outputType, "float") == 0 )
    {
    WriteImageFLOAT(filename, image);
    }
  else // default
    {
    WriteImage(filename, image);
    }
}

// calculate distance between two images by mean squared difference
template <class TInputImage, class TOutputImage>
double
MABMISImageOperationFilter<TInputImage, TOutputImage>
::calculateDistanceMSD(string& imageName1, string&  imageName2)
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
  float value1, value2;
  for( it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2 )
    {
    value1 = it1.Get();
    value2 = it2.Get();
    dist += (value1 - value2) * (value1 - value2);
    }
  dist = dist / (size[0] * size[1] * size[2]);

  return dist;
}

// calculate pair-wise distance in an image population
template <class TInputImage, class TOutputImage>
void
MABMISImageOperationFilter<TInputImage, TOutputImage>
::PairwiseDistanceAmongImages(vector<string> imageFileNames, int totalNumber, vnl_matrix<double>& distanceMatrix)

{
  for( int i = 0; i < totalNumber; i++ )
    {
    // std::cout << imageFileNames[i] << std::endl;
    for( int j = i; j < totalNumber; j++ )
      {
      if( j == i )
        {
        distanceMatrix[i][j] = 0.0;
        }
      else
        {
        double dist = 0.0;
        // dist = calculateDistance(imageFileNames[i], imageFileNames[j]);
        dist = calculateDistanceMSD(imageFileNames[i], imageFileNames[j]);
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
