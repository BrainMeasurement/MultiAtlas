#ifndef __itkMABMISSimulateData_h
#define __itkMABMISSimulateData_h

#include <itkImage.h>
#include <itkImageToImageFilter.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include "itkImage.h"

#include "itkMABMISDeformationFieldFilter.h"
#include "itkMABMISImageOperationFilter.h"
#include "itkMABMISBasicOperationFilter.h"

#define ImageDimension 3

namespace itk
{
namespace Statistics
{
template <class TInputImage, class TOutputImage>
class ITK_EXPORT MABMISSimulateData : public ImageToImageFilter<TInputImage, TOutputImage>
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

  int  DoPCATraining(std::vector<std::string> deformationFieldFileNames, int numFiles,
                     std::vector<std::string> allImgFileName, int root);

  void LoadIntoArray(string resampledDeformationFieldFileName, float* df_vector);

  void SaveFromArray(string  deformationFieldFileName, float* df_vector, int sx, int sy, int sz);

  itkSetMacro(Root, int);
  itkSetMacro(Imx, int);
  itkSetMacro(Imy, int);
  itkSetMacro(Imz, int);
  itkSetMacro(AtlasSize, int);
  itkSetMacro(SimulateSize, int);
private:
  MABMISSimulateData(const Self &); // purposely not implemented
  void operator=(const Self &);     // purposely not implemented

  int m_Root;

  int m_Imx;
  int m_Imy;
  int m_Imz;

  int m_AtlasSize;
  int m_SimulateSize;
protected:
  MABMISSimulateData();
  ~MABMISSimulateData();
  void PrintSelf(std::ostream& os, Indent indent) const;
};
} // namespace itk
} // namespace Statistics

#include "itkMABMISSimulateData.hxx"

#endif
