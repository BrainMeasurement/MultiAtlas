#ifndef __itkMABMISBasicOperationFilter_h
#define __itkMABMISBasicOperationFilter_h

#include <itkImage.h>
#include <itkImageToImageFilter.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#define ImageDimension 3

namespace itk
{
namespace Statistics
{
template <class TInputImage, class TOutputImage>
class ITK_EXPORT MABMISBasicOperationFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:

  typedef MABMISBasicOperationFilter                    Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  typedef TInputImage                 ImageType;
  typedef typename ImageType::Pointer ImagePointerType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MABMISBasicOperationFilter, ImageToImageFilter);

  void RemoveFile(string filename);

  void bubbleSort(double* arr, int* index, int n);

  void myitoa(int num, string& str, int digit);

  void SaveMatrix2File(vnl_matrix<double> matrix, int iSize, int jSize, string martixFileName);

private:
  MABMISBasicOperationFilter(const Self &); // purposely not implemented
  void operator=(const Self &);             // purposely not implemented

protected:
  MABMISBasicOperationFilter();
  ~MABMISBasicOperationFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;
};
} // namespace itk
} // namespace Statistics

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMABMISBasicOperationFilter.txx"
#endif

#endif
