#ifndef __commonMABMIS_h
#define __commonMABMIS_h

// for math
#include <vnl/vnl_random.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matlab_print.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_svd_economy.h>
#include <string>
#include <vector>

#ifdef WIN32
constexpr char FILESEP = '\\';
#else
constexpr char FILESEP = '/';
#endif

// including itksys::SystemTools::MakeDirectory(char*)
#include <itksys/SystemTools.hxx>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageToImageFilter.h"
#include "itkImageRegionIterator.h"

constexpr unsigned int ImageDimension = 3;

// basic itk
#include "itkVector.h"

// registration
#include "itkImageRegistrationMethod.h"
#include "itkSymmetricForcesDemonsRegistrationFilter.h"

// interpolator
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkNearestNeighborExtrapolateImageFunction.h"

// filter
#include "itkResampleImageFilter.h"

#include "itkCastImageFilter.h"
#include "itkWarpImageFilter.h"

#include "itkHistogramMatchingImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkWarpVectorImageFilter.h"
#include "itkInverseDisplacementFieldImageFilter.h"

// for affine transformation
#include "itkTransform.h"
#include "itkAffineTransform.h"
#include "itkImageRegistrationMethod.h"

// for Diffeomorphic Demons
#include <itkCommand.h>
#include <itkDiffeomorphicDemonsRegistrationFilter.h>
#include <itkMultiResolutionPDEDeformableRegistration.h>

// To include all related header files
#include "itkMABMISBasicOperationFilter.h"
#include "itkMABMISImageOperationFilter.h"
#include "itkMABMISDeformationFieldFilter.h"
#include "itkMABMISSimulateData.h"
#include "itkMABMISImageRegistrationFilter.h"
#include "itkMABMISTreeOperation.h"
#include "itkMABMISAtlasXMLFile.h"

typedef double CoordinateRepType;

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
typedef itk::ImageFileReader<CharImageType>     CharImageReaderType;
typedef itk::ImageFileReader<InternalImageType> InternalImageReaderType;
typedef itk::ImageFileWriter<InternalImageType> InternalImageWriterType;

typedef itk::WarpImageFilter<InternalImageType, InternalImageType, DeformationFieldType> InternalWarpFilterType;
typedef itk::ImageFileWriter<CharImageType>                                              CharImageWriterType;
typedef itk::ImageFileWriter<IntImageType>                                               IntImageWriterType;
typedef itk::ImageFileWriter<FloatImageType>                                             FloatImageWriterType;
typedef itk::ImageFileWriter<ShortImageType>                                             ShortImageWriterType;

typedef itk::ImageFileReader<DeformationFieldType> DeformationFieldReaderType;
typedef itk::ImageFileWriter<DeformationFieldType> DeformationFieldWriterType;

//////////////////////////////////////////////////////////////////////////////
// image filter type
typedef itk::HistogramMatchingImageFilter<InternalImageType, InternalImageType> InternalHistMatchFilterType;

////////////////////////////////////////////////////////////////////////////
// operation on deformation fields
typedef itk::WarpVectorImageFilter<DeformationFieldType, DeformationFieldType, DeformationFieldType> WarpVectorFilterType;
typedef itk::AddImageFilter<DeformationFieldType, DeformationFieldType, DeformationFieldType>        AddImageFilterType;

typedef itk::Vector<ShortPixelType, ImageDimension>      ShortVectorPixelType;
typedef itk::Image<ShortVectorPixelType, ImageDimension> ShortDeformationFieldType;
typedef itk::ImageFileWriter<ShortDeformationFieldType>  ShortDeformationFieldWriterType;

std::string ReplacePathSepForOS(const std::string & input);

std::string GetRootName(const std::string & filename);

std::string GetExtension(const std::string & filename);


// global bool variables to adjust the  procedure
constexpr bool isEvaluate = false; // if false, we do not know the ground-truth of labels
constexpr bool isDebug = false;    // false;//true; // if true, print out more information
constexpr int localPatchSize = 1;  // (2r+1)*(2r+1)*(2r+1) is the volume of local patch
constexpr bool doHistMatch = true;
constexpr bool doHistMatchTrue = true;
constexpr bool doHistMatchFalse = false;

//demons registration parameters
//extern int iterInResolutions[4][3];
//extern int itereach;
//extern int itereach0;
//extern int itereach1;
//extern int itereach2;
//extern int itereach3;
//extern double sigmaDef;
//extern double sigmaDef10;
//extern double sigmaDef15;
//extern double sigmaDef20;
//extern double sigmaDef25;
//extern double sigmaDef30;
//extern double sigmaDef35;

typedef itk::Statistics::MABMISSimulateData<InternalImageType, InternalImageType> DataSimulatorType;
extern DataSimulatorType::Pointer datasimulator;
typedef itk::Statistics::MABMISImageOperationFilter<CharImageType, CharImageType> ImageOperationFilterType;
extern ImageOperationFilterType::Pointer imgoperator;
typedef itk::Statistics::MABMISDeformationFieldFilter<InternalImageType,
    InternalImageType> DeformationFieldOperationFilterType;
extern DeformationFieldOperationFilterType::Pointer dfoperator;
typedef itk::Statistics::MABMISImageRegistrationFilter<CharImageType, CharImageType> ImageRegistrationFilterType;
extern ImageRegistrationFilterType::Pointer regoperator;
typedef itk::Statistics::MABMISTreeOperation<InternalImageType, InternalImageType> TreeOperationType;
extern TreeOperationType::Pointer treeoperator;
typedef itk::Statistics::MABMISBasicOperationFilter<CharImageType, CharImageType> BasicOperationFilterType;
extern BasicOperationFilterType::Pointer basicoperator;

constexpr char ImageSuffix[] = "_cbq_000.nii.gz";
constexpr char DeformationSuffix[] = "_deform_000.nii.gz";
constexpr char SegmentationSuffix[] = "_seg_000.nii.gz";
constexpr char RegisteredSuffix[] = "_cbq_reg.nii.gz";

#endif //__commonMABMIS_h
