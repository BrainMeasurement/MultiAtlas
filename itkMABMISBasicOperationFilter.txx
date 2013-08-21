#ifndef __itkMABMISBasicOperationFilter_txx
#define __itkMABMISBasicOperationFilter_txx

#include "itkMABMISBasicOperationFilter.h"

using namespace std;

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
  if (itksys::SystemTools::FileExists(filename.c_str(), true))
	{
#if defined(_WIN32) && !defined(__CYGWIN__)
		std::ostringstream oss (std::ostringstream::out);
		oss << "del /F " << filename << std::endl;
		//std::cerr << oss.str().c_str();
		system(oss.str().c_str());

#else
		std::ostringstream oss (std::ostringstream::out);
		oss << "rm -f " << filename << std::endl;
		//std::cerr << oss.str().c_str();
		system(oss.str().c_str());
#endif
	}
	return;
}


// bubble sort on double array
template <class TInputImage, class TOutputImage>
void
MABMISBasicOperationFilter<TInputImage, TOutputImage>
::bubbleSort(double* arr, int* index, int n) 
{
	bool swapped = true;
	int j = 0;
	double tmp;
	int index_tmp;
	while (swapped) 
	{
		swapped = false;
		j++;
		for (int i = 0; i < n - j; i++) 
		{
			if (arr[i] > arr[i + 1]) 
			{
				tmp = arr[i];
				arr[i] = arr[i + 1];
				arr[i + 1] = tmp;

				index_tmp = index[i];
				index[i] = index[i+1];
				index[i+1] = index_tmp;

				swapped = true;
			}
		}
	}
}

template <class TInputImage, class TOutputImage>
void
MABMISBasicOperationFilter<TInputImage, TOutputImage>
::myitoa(int num, string& str, int digit)
{
  str = "";

  if (num<10) str.append("00");
  else if (num<100) str.append("0");

  stringstream st;
  st << num;
  str.append(st.str());

	return;   
}

template <class TInputImage, class TOutputImage>
void
MABMISBasicOperationFilter<TInputImage, TOutputImage>
::SaveMatrix2File(vnl_matrix<double> matrix, int iSize, int jSize, std::string martixFileName)
{

	std::ofstream outfile;
	outfile.open (martixFileName.c_str()); 

	for(int i = 0; i < iSize; i++)
	{
		for(int j = 0; j < jSize; j++)
		{

				outfile <<  ' ' ;
				outfile << std::right << std::fixed << std::setw(8) << std::setprecision(2) << matrix[i][j];
		}
		outfile << std::endl;
	}
	outfile.close();
	return;
}



} // namespace Statistics
}  // namespace itk

#endif
