#ifndef __itkMABMISImageOperationFilter_h
#define __itkMABMISImageOperationFilter_h

#include "commonMABMIS.h"

#include "itkCastImageFilter.h"
#include "itkWarpImageFilter.h"

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
  typedef itk::Image<InternalPixelType, ImageDimension> InternalImageType;
  typedef itk::Image<VectorPixelType, ImageDimension>   DeformationFieldType;

  // basic iterator type
  typedef itk::ImageRegionIterator<DeformationFieldType> DeformationFieldIteratorType;
  typedef itk::ImageRegionIterator<InternalImageType>    InternalImageIteratorType;

  // basic image reader/writer related type
  typedef itk::ImageFileReader<InternalImageType>                 InternalImageReaderType;
  typedef itk::ImageFileWriter<InternalImageType>                 InternalImageWriterType;

  typedef itk::WarpImageFilter<InternalImageType, InternalImageType, DeformationFieldType> InternalWarpFilterType;

  typedef itk::ImageFileReader<DeformationFieldType> DeformationFieldReaderType;
  typedef itk::ImageFileWriter<DeformationFieldType> DeformationFieldWriterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MABMISImageOperationFilter, ImageToImageFilter);

  void PairwiseDistanceAmongImages(const std::vector<std::string>& imageFileNames, int totalNumber,
                                   vnl_matrix<double>& distanceMatrix);

  double calculateDistanceMSD(const std::string& filename1, const std::string& filename2);

  int ReadImage(const std::string filename, InternalImageType::Pointer& image);

  itk::ImageIOBase::IOComponentType GetIOPixelType(const std::string& filename);

  /** Write the image to disk, after casting it the the provided pixel type. */
  template <class PixelType>
  int WriteImage(const std::string filename, InternalImageType::Pointer image);

  int WriteImage(const std::string filename, InternalImageType::Pointer image, itk::ImageIOBase::IOComponentType ioType);

  void WriteImage(const std::string filename, InternalImageType::Pointer image, char* outputType);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MABMISImageOperationFilter);

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
