#ifndef __itkMABMISBasicOperationFilter_h
#define __itkMABMISBasicOperationFilter_h

#include "commonMABMIS.h"

namespace itk
{
namespace Statistics
{
template <class TInputImage, class TOutputImage>
class MABMISBasicOperationFilter : public ImageToImageFilter<TInputImage, TOutputImage>
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

  void RemoveFile(std::string filename);

  void bubbleSort(double* arr, int* index, int n);

  void myitoa(int num, std::string& str, int digit);

  void SaveMatrix2File(vnl_matrix<double> matrix, int iSize, int jSize, std::string martixFileName);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MABMISBasicOperationFilter);

protected:
  MABMISBasicOperationFilter();
  ~MABMISBasicOperationFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;
};
} // namespace itk
} // namespace Statistics

#include "itkMABMISBasicOperationFilter.hxx"
#endif
