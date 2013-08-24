#ifndef __itkMABMISImageOperationFilter_h
#define __itkMABMISImageOperationFilter_h

#include <itkImage.h>
#include <itkImageToImageFilter.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkCastImageFilter.h"

#include "itkImage.h"

#define ImageDimension 3

using namespace std;

namespace itk
{
namespace Statistics
{
template <class TInputImage, class TOutputImage>
class MABMISImageOperationFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:

  typedef MABMISImageOperationFilter                    Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  typedef TInputImage                 ImageType;
  typedef typename ImageType::Pointer ImagePointerType;

  typedef itk::ImageFileReader<ImageType>     ImageReaderType;
  typedef itk::ImageRegionIterator<ImageType> ImageIteratorType;

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
  typedef itk::ImageFileReader<CharImageType>                     CharImageReaderType;
  typedef itk::ImageFileReader<InternalImageType>                 InternalImageReaderType;
  typedef itk::ImageFileWriter<InternalImageType>                 InternalImageWriterType;
  typedef itk::CastImageFilter<InternalImageType, CharImageType>  Internal2CharCastFilterType;
  typedef itk::CastImageFilter<InternalImageType, FloatImageType> Internal2FloatCastFilterType;
  typedef itk::CastImageFilter<InternalImageType, IntImageType>   Internal2IntCastFilterType;
  typedef itk::CastImageFilter<InternalImageType, ShortImageType> Internal2ShortCastFilterType;

  typedef itk::WarpImageFilter<InternalImageType, InternalImageType, DeformationFieldType> InternalWarpFilterType;
  typedef itk::ImageFileWriter<CharImageType>                                              CharImageWriterType;
  typedef itk::ImageFileWriter<IntImageType>                                               IntImageWriterType;
  typedef itk::ImageFileWriter<FloatImageType>                                             FloatImageWriterType;
  typedef itk::ImageFileWriter<ShortImageType>                                             ShortImageWriterType;

  typedef itk::ImageFileReader<DeformationFieldType> DeformationFieldReaderType;
  typedef itk::ImageFileWriter<DeformationFieldType> DeformationFieldWriterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MABMISImageOperationFilter, ImageToImageFilter);

  void PairwiseDistanceAmongImages(vector<std::string> imageFileNames, int totalNumber,
                                   vnl_matrix<double>& distanceMatrix);

  double calculateDistanceMSD(string& filename1, string& filename2);

  int ReadImage(string filename, InternalImageType::Pointer& image);

  void WriteImageUCHAR(string filename, InternalImageType::Pointer image);

  void WriteImageINT(string filename, InternalImageType::Pointer image);

  void WriteImageSHORT(string filename, InternalImageType::Pointer image);

  void WriteImageFLOAT(string filename, InternalImageType::Pointer image);

  void WriteImage(string filename, InternalImageType::Pointer image);

  void WriteImage(string filename, InternalImageType::Pointer image, char* outputType);

private:
  MABMISImageOperationFilter(const Self &); // purposely not implemented
  void operator=(const Self &);             // purposely not implemented

  int m_Root;
protected:
  MABMISImageOperationFilter();
  ~MABMISImageOperationFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;
};
} // namespace itk
} // namespace Statistics

#include "itkMABMISImageOperationFilter.hxx"

#endif
