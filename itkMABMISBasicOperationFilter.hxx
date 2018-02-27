#ifndef __itkMABMISBasicOperationFilter_hxx
#define __itkMABMISBasicOperationFilter_hxx

#include "itkMABMISBasicOperationFilter.h"
#include <iomanip>

namespace itk
{
namespace Statistics
{
template <class TInputImage, class TOutputImage>
MABMISBasicOperationFilter<TInputImage, TOutputImage>
::MABMISBasicOperationFilter()
{
}

template <class TInputImage, class TOutputImage>
MABMISBasicOperationFilter<TInputImage, TOutputImage>
::~MABMISBasicOperationFilter()
{
}

template <class TInputImage, class TOutputImage>
void
MABMISBasicOperationFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}

template <class TInputImage, class TOutputImage>
void
MABMISBasicOperationFilter<TInputImage, TOutputImage>
::RemoveFile(std::string filename)
{
  int retVal;
  if( itksys::SystemTools::FileExists(filename.c_str(), true) )
    {
    if (!itksys::SystemTools::RemoveFile(filename))
      {
      std::cout << " Error deleting " << filename << std::endl;
      }
    }
  return;
}

// bubble sort on double array
template <class TInputImage, class TOutputImage>
void
MABMISBasicOperationFilter<TInputImage, TOutputImage>
::bubbleSort(double* arr, int* index, int n)
{
  bool   swapped = true;
  int    j = 0;
  double tmp;
  int    index_tmp;

  while( swapped )
    {
    swapped = false;
    ++j;
    for( int i = 0; i < n - j; ++i )
      {
      if( arr[i] > arr[i + 1] )
        {
        tmp = arr[i];
        arr[i] = arr[i + 1];
        arr[i + 1] = tmp;

        index_tmp = index[i];
        index[i] = index[i + 1];
        index[i + 1] = index_tmp;

        swapped = true;
        }
      }
    }
}

template <class TInputImage, class TOutputImage>
std::string
MABMISBasicOperationFilter<TInputImage, TOutputImage>
::myitoa(int number, int digits)
{
  std::stringstream ss;
  ss << std::setw(digits) << std::setfill('0') << number;
  return ss.str();
}

template <class TInputImage, class TOutputImage>
void
MABMISBasicOperationFilter<TInputImage, TOutputImage>
::SaveMatrix2File(vnl_matrix<double> matrix, int iSize, int jSize, std::string martixFileName)
{
  std::ofstream outfile;

  outfile.open(martixFileName.c_str() );
  for( int i = 0; i < iSize; ++i )
    {
    for( int j = 0; j < jSize; ++j )
      {
      outfile <<  ' ';
      outfile << std::right << std::fixed << std::setw(8) << std::setprecision(2) << matrix[i][j];
      }
    outfile << std::endl;
    }
  outfile.close();
  return;
}
} // namespace Statistics
} // namespace itk

#endif
