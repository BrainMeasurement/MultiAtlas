#ifndef __itkMABMISDeformationFieldFilter_h
#define __itkMABMISDeformationFieldFilter_h

#include "commonMABMIS.h"

#include "itkWarpImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkMultiplyImageFilter.h"

// interpolator
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkMABMISImageOperationFilter.h"
#include "itkWarpVectorImageFilter.h"

namespace itk
{
namespace Statistics
{
template <class TInputImage, class TOutputImage>
class MABMISDeformationFieldFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:

  typedef MABMISDeformationFieldFilter                  Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  typedef TInputImage                 ImageType;
  typedef typename ImageType::Pointer ImagePointerType;

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

  typedef itk::ImageFileReader<DeformationFieldType> DeformationFieldReaderType;
  typedef itk::ImageFileWriter<DeformationFieldType> DeformationFieldWriterType;

  typedef itk::WarpVectorImageFilter<DeformationFieldType, DeformationFieldType,
                                     DeformationFieldType>         WarpVectorFilterType;
  typedef itk::MultiplyImageFilter<DeformationFieldType, float,
                                             DeformationFieldType> MultiplyDeformationFieldFilterType;
  typedef itk::DivideImageFilter<DeformationFieldType, float,
                                           DeformationFieldType>   DivideDeformationFieldFilterType;
  typedef itk::AddImageFilter<DeformationFieldType, DeformationFieldType,
                              DeformationFieldType>                AddImageFilterType;

// basic iterator type
  typedef itk::ImageRegionIterator<DeformationFieldType> DeformationFieldIteratorType;
  typedef itk::ImageRegionIterator<InternalImageType>    InternalImageIteratorType;
  typedef itk::ImageRegionIterator<CharImageType>        CharImageIteratorType;

// interpolator type
  typedef itk::LinearInterpolateImageFunction<InternalImageType, double>          InternalLinearInterpolatorType;
  typedef itk::NearestNeighborInterpolateImageFunction<InternalImageType, double> InternalNNInterpolatorType;

  typedef itk::WarpImageFilter<InternalImageType, InternalImageType, DeformationFieldType> InternalWarpFilterType;

  typedef itk::Statistics::MABMISImageOperationFilter<ImageType, ImageType> ImageOperationType;
  typename ImageOperationType::Pointer imgoperator;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MABMISDeformationFieldFilter, ImageToImageFilter);

  DeformationFieldType::SpacingType   df_spacing;
  DeformationFieldType::DirectionType df_direction;
  DeformationFieldType::PointType     df_origin;

  int ReadDeformationField(std::string filename, DeformationFieldType::Pointer & deformationfield);

  void WriteDeformationField(std::string  filename, DeformationFieldType::Pointer deformationfield);

  void ComposeDeformationFieldsAndSave(std::string inputDeformationFieldFileName, std::string deformationFieldFileName,
                                       std::string composedDeformationFieldFileName);

  void ComposeDeformationFields(DeformationFieldType::Pointer input, DeformationFieldType::Pointer deformationField,
                                DeformationFieldType::Pointer & composedDeformationField);

  void InverseDeformationField3D(DeformationFieldType::Pointer deformationField,
                                 DeformationFieldType::Pointer & deformationFieldInverse);

  void ApplyDeformationField(InternalImageType::Pointer movingImage, DeformationFieldType::Pointer deformationField,
                             InternalImageType::Pointer & deformedImage, bool isLinearInterpolator);

  void ApplyDeformationFieldAndWriteWithTypeWithFileNames(std::string  movingImageFileName,
                                                          std::string deformationFieldFileName,
                                                          std::string deformedImageFileName, bool isLinear);

  void DownResampleDeformationField(std::string deformationFieldFileName, std::string resampledDeformationFieldFileName,
                                    int sampleRate);

  void UpResampleDeformationField(std::string deformationFieldFileName, std::string  resampledDeformationFieldFileName,
                                  int sampleRate);

  itkSetMacro(Imx, int);
  itkSetMacro(Imy, int);
  itkSetMacro(Imz, int);
private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MABMISDeformationFieldFilter);

  int m_Imx;
  int m_Imy;
  int m_Imz;
protected:
  MABMISDeformationFieldFilter();
  ~MABMISDeformationFieldFilter();
  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;
};
} // namespace itk
} // namespace Statistics

#include "itkMABMISDeformationFieldFilter.hxx"

#endif
