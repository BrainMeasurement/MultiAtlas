#ifndef __itkMABMISSimulateData_h
#define __itkMABMISSimulateData_h

#include "commonMABMIS.h"

#include "itkMABMISDeformationFieldFilter.h"
#include "itkMABMISImageOperationFilter.h"
#include "itkMABMISBasicOperationFilter.h"

namespace itk
{
namespace Statistics
{
template <class TInputImage, class TOutputImage>
class MABMISSimulateData : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:

  typedef MABMISSimulateData                            Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  typedef TInputImage                 ImageType;
  typedef typename ImageType::Pointer ImagePointerType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MABMISSimulateData, ImageToImageFilter);

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

  DeformationFieldType::SpacingType   df_spacing;
  DeformationFieldType::DirectionType df_direction;
  DeformationFieldType::PointType     df_origin;

  typedef itk::Statistics::MABMISDeformationFieldFilter<ImageType, ImageType> DeformationFieldOperationType;
  typedef itk::Statistics::MABMISImageOperationFilter<ImageType, ImageType>   ImageOperationType;
  typedef itk::Statistics::MABMISBasicOperationFilter<ImageType, ImageType>   BasicOperationType;

  typename DeformationFieldOperationType::Pointer dfoperator;
  typename ImageOperationType::Pointer imgoperator;
  typename BasicOperationType::Pointer basicoperator;

  int DoPCATraining(const std::vector<std::string>& deformationFieldFileNames, int numFiles,
                    const std::vector<std::string>& allImgFileName, int root);

  void LoadIntoArray(const std::string& resampledDeformationFieldFileName, float* df_vector);

  void SaveFromArray(const std::string& deformationFieldFileName, float* df_vector, int sx, int sy, int sz);

  itkSetMacro(Root, int);
  itkSetMacro(Imx, int);
  itkSetMacro(Imy, int);
  itkSetMacro(Imz, int);
  itkSetMacro(AtlasSize, int);
  itkSetMacro(SimulateSize, int);
private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MABMISSimulateData);

  int m_Root;

  int m_Imx;
  int m_Imy;
  int m_Imz;

  int m_AtlasSize;
  int m_SimulateSize;
protected:
  MABMISSimulateData();
  ~MABMISSimulateData();
  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;
};
} // namespace itk
} // namespace Statistics

#include "itkMABMISSimulateData.hxx"

#endif
