#ifndef __itkMABMISImageRegistrationFilter_hxx
#define __itkMABMISImageRegistrationFilter_hxx

#include "itkMABMISImageRegistrationFilter.h"

namespace itk
{
namespace Statistics
{
template <class TInputImage, class TOutputImage>
MABMISImageRegistrationFilter<TInputImage, TOutputImage>
::MABMISImageRegistrationFilter()
{
  dfoperator = DeformationFieldOperationType::New();
  imgoperator = ImageOperationType::New();
}

template <class TInputImage, class TOutputImage>
MABMISImageRegistrationFilter<TInputImage, TOutputImage>
::~MABMISImageRegistrationFilter()
{
}

template <class TInputImage, class TOutputImage>
void
MABMISImageRegistrationFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}

template <class TInputImage, class TOutputImage>
int
MABMISImageRegistrationFilter<TInputImage, TOutputImage>
::DiffeoDemonsRegistrationWithParameters(std::string fixedImageFileName, std::string  movingImageFileName,
                                         std::string deformedImageFileName, std::string deformationFieldFileName,
                                         double sigmaDef, bool doHistMatch, std::vector<int> iterInResolutions)
{
  // for debugging
  //std::cerr << "DiffeoDemonsRegistrationWithParameters" << std::endl;
  //std::cerr << fixedImageFileName << std::endl;
  //std::cerr << movingImageFileName << std::endl;
  //std::cerr << deformedImageFileName << std::endl;
  //std::cerr << deformationFieldFileName << std::endl;

  int res = iterInResolutions.size();
  // read fixed and moving images
  InternalImageReaderType::Pointer fixedImageReader   = InternalImageReaderType::New();
  InternalImageReaderType::Pointer movingImageReader  = InternalImageReaderType::New();

  fixedImageReader->SetFileName( fixedImageFileName );
  movingImageReader->SetFileName( movingImageFileName );

  // do histogram matching and some parameters need to set manually
  InternalHistMatchFilterType::Pointer matcher = InternalHistMatchFilterType::New();
  if( doHistMatch )
    {
    matcher->SetInput( movingImageReader->GetOutput() );
    matcher->SetReferenceImage( fixedImageReader->GetOutput() );
    matcher->SetNumberOfHistogramLevels( 1024 ); // 1024
    matcher->SetNumberOfMatchPoints( 7 );
    matcher->ThresholdAtMeanIntensityOn();
    }

  // Set up the demons filter
  typedef itk::PDEDeformableRegistrationFilter
    <InternalImageType, InternalImageType, DeformationFieldType>   BaseRegistrationFilterType;
  BaseRegistrationFilterType::Pointer filter;

  // s <- s o exp(u) (Diffeomorphic demons)
  typedef itk::DiffeomorphicDemonsRegistrationFilter
    <InternalImageType, InternalImageType, DeformationFieldType>
    ActualRegistrationFilterType;
  typedef ActualRegistrationFilterType::GradientType GradientType;

  ActualRegistrationFilterType::Pointer actualfilter = ActualRegistrationFilterType::New();

  float maxStepLength = 2.0;
  actualfilter->SetMaximumUpdateStepLength( maxStepLength );
  unsigned int gradientType = 0;
  actualfilter->SetUseGradientType(static_cast<GradientType>(gradientType) );
  filter = actualfilter;

  // set up smoothing kernel for deformation field
  // float sigmaDef = 1.5;
  if( sigmaDef > 0.1 )
    {
    filter->SmoothDisplacementFieldOn();
    filter->SetStandardDeviations( sigmaDef );
    }
  else
    {
    filter->SmoothDisplacementFieldOff();
    }

  // set up smoothing kernel for update field
  float sigmaUp = 0.0;
  if( sigmaUp > 0.1 )
    {
    filter->SmoothUpdateFieldOn();
    filter->SetUpdateFieldStandardDeviations( sigmaUp );
    }
  else
    {
    filter->SmoothUpdateFieldOff();
    }

  // filter->SetIntensityDifferenceThreshold( 0.001 );

  // Set up the multi-resolution filter
  typedef itk::MultiResolutionPDEDeformableRegistration<
      InternalImageType, InternalImageType, DeformationFieldType, InternalPixelType>   MultiResRegistrationFilterType;
  MultiResRegistrationFilterType::Pointer multires = MultiResRegistrationFilterType::New();

  typedef itk::VectorLinearInterpolateNearestNeighborExtrapolateImageFunction<DeformationFieldType,
                                                                              double> FieldInterpolatorType;

  FieldInterpolatorType::Pointer VectorInterpolator = FieldInterpolatorType::New();

#if ( ITK_VERSION_MAJOR > 3 ) || ( ITK_VERSION_MAJOR == 3 && ITK_VERSION_MINOR > 8 )
  multires->GetFieldExpander()->SetInterpolator(VectorInterpolator);
#endif
  std::vector<unsigned int> curNumIterations;
  // unsigned int curNumOfIterations[] = {15,10,5};
  for( int i = 0; i < res; ++i )
    {
    curNumIterations.push_back(iterInResolutions[i]);
    }

  multires->SetRegistrationFilter( filter );
  multires->SetNumberOfLevels( curNumIterations.size() );
  multires->SetNumberOfIterations( &curNumIterations[0] );
  multires->SetFixedImage( fixedImageReader->GetOutput() );
  if( doHistMatch )
    {
    multires->SetMovingImage( matcher->GetOutput() );
    }
  else
    {
    multires->SetMovingImage( movingImageReader->GetOutput() );
    }

  // Compute the deformation field
  try
    {
    multires->UpdateLargestPossibleRegion();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught at demons registration between input images!" << std::endl;
    std::cerr << excep << std::endl;
    return -1;
    }

  // write deformation field into a file
  DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
  fieldWriter->SetFileName( deformationFieldFileName );
  fieldWriter->SetInput( multires->GetOutput() );
  fieldWriter->Update();

  dfoperator->ApplyDeformationFieldAndWriteWithTypeWithFileNames(movingImageFileName,
                                                                 deformationFieldFileName, deformedImageFileName, true);

  //// write deformed image
  // InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
  // InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

  // warper->SetInput( movingImageReader->GetOutput() );
  // warper->SetInterpolator( interpolator );
  // warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
//	warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
// warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
// warper->SetDeformationField( multires->GetOutput() );
// CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
//	Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
//	writer->SetFileName( deformedImageFileName );
//	caster->SetInput( warper->GetOutput() );
// writer->SetInput( caster->GetOutput()   );
//	writer->Update();

  return 0;
}

template <class TInputImage, class TOutputImage>
int
MABMISImageRegistrationFilter<TInputImage, TOutputImage>
::DiffeoDemonsRegistrationWithInitialWithParameters(std::string  fixedImageFileName, std::string movingImageFileName,
                                                    std::string initDeformationFieldFileName, std::string deformedImageFileName,
                                                    std::string deformationFieldFileName, double sigmaDef, bool doHistMatch,
                                                    std::vector<int> iterInResolutions)
{
  int res = iterInResolutions.size();
  // read initial deformation field file
  DeformationFieldType::Pointer initDeformationField = 0;

  dfoperator->ReadDeformationField(initDeformationFieldFileName, initDeformationField);

  // read fixed and moving images
  InternalImageReaderType::Pointer fixedImageReader   = InternalImageReaderType::New();
  InternalImageReaderType::Pointer movingImageReader  = InternalImageReaderType::New();

  fixedImageReader->SetFileName( fixedImageFileName );
  movingImageReader->SetFileName( movingImageFileName );

  // do histogram matching and some parameters need to set manually
  InternalHistMatchFilterType::Pointer matcher = InternalHistMatchFilterType::New();
  if( doHistMatch )
    {
    matcher->SetInput( movingImageReader->GetOutput() );
    matcher->SetReferenceImage( fixedImageReader->GetOutput() );
    matcher->SetNumberOfHistogramLevels( 1024 ); // 1024
    matcher->SetNumberOfMatchPoints( 7 );
    matcher->ThresholdAtMeanIntensityOn();
    }

  // Set up the demons filter
  typedef itk::PDEDeformableRegistrationFilter
    <InternalImageType, InternalImageType, DeformationFieldType>   BaseRegistrationFilterType;
  BaseRegistrationFilterType::Pointer filter;

  // s <- s o exp(u) (Diffeomorphic demons)
  typedef itk::DiffeomorphicDemonsRegistrationFilter
    <InternalImageType, InternalImageType, DeformationFieldType>
    ActualRegistrationFilterType;
  typedef ActualRegistrationFilterType::GradientType GradientType;

  ActualRegistrationFilterType::Pointer actualfilter = ActualRegistrationFilterType::New();

  float maxStepLength = 2.0;
  actualfilter->SetMaximumUpdateStepLength( maxStepLength );
  unsigned int gradientType = 0;
  actualfilter->SetUseGradientType(static_cast<GradientType>(gradientType) );
  filter = actualfilter;

  // set up smoothing kernel for deformation field
  // float sigmaDef = 1.5;
  if( sigmaDef > 0.1 )
    {
    filter->SmoothDisplacementFieldOn();
    filter->SetStandardDeviations( sigmaDef );
    }
  else
    {
    filter->SmoothDisplacementFieldOff();
    }

  // set up smoothing kernel for update field
  float sigmaUp = 0.0;
  if( sigmaUp > 0.1 )
    {
    filter->SmoothUpdateFieldOn();
    filter->SetUpdateFieldStandardDeviations( sigmaUp );
    }
  else
    {
    filter->SmoothUpdateFieldOff();
    }

  // filter->SetIntensityDifferenceThreshold( 0.001 );

  // Set up the multi-resolution filter
  typedef itk::MultiResolutionPDEDeformableRegistration<
      InternalImageType, InternalImageType, DeformationFieldType, InternalPixelType>   MultiResRegistrationFilterType;
  MultiResRegistrationFilterType::Pointer multires = MultiResRegistrationFilterType::New();

  typedef itk::VectorLinearInterpolateNearestNeighborExtrapolateImageFunction<DeformationFieldType,
                                                                              double> FieldInterpolatorType;

  FieldInterpolatorType::Pointer VectorInterpolator = FieldInterpolatorType::New();

#if ( ITK_VERSION_MAJOR > 3 ) || ( ITK_VERSION_MAJOR == 3 && ITK_VERSION_MINOR > 8 )
  multires->GetFieldExpander()->SetInterpolator(VectorInterpolator);
#endif
  std::vector<unsigned int> curNumIterations;
  // unsigned int curNumOfIterations[] = {15,10,5};
  for( int i = 0; i < res; ++i )
    {
    curNumIterations.push_back(iterInResolutions[i]);
    }

  multires->SetRegistrationFilter( filter );
  multires->SetNumberOfLevels( curNumIterations.size() );
  multires->SetNumberOfIterations( &curNumIterations[0] );
  multires->SetFixedImage( fixedImageReader->GetOutput() );
  if( doHistMatch )
    {
    multires->SetMovingImage( matcher->GetOutput() );
    }
  else
    {
    multires->SetMovingImage( movingImageReader->GetOutput() );
    }

  // set initial
  multires->SetArbitraryInitialDisplacementField( initDeformationField );

  // Compute the deformation field
  try
    {
    multires->UpdateLargestPossibleRegion();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught at demons registration between input images!" << std::endl;
    std::cerr << excep << std::endl;
    return -1;
    }

  // write deformation field into a file
  DeformationFieldWriterType::Pointer fieldWriter = DeformationFieldWriterType::New();
  fieldWriter->SetFileName( deformationFieldFileName );
  fieldWriter->SetInput( multires->GetOutput() );
  fieldWriter->Update();

  dfoperator->ApplyDeformationFieldAndWriteWithTypeWithFileNames(movingImageFileName,
                                                                 deformationFieldFileName, deformedImageFileName, true);

  //// write deformed image
  // InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();
  // InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();

  // warper->SetInput( movingImageReader->GetOutput() );
  // warper->SetInterpolator( interpolator );
  // warper->SetOutputSpacing( fixedImageReader->GetOutput()->GetSpacing() );
  // warper->SetOutputOrigin( fixedImageReader->GetOutput()->GetOrigin() );
  // warper->SetOutputDirection( fixedImageReader->GetOutput()->GetDirection() );
  // warper->SetDeformationField( multires->GetOutput() );
  // CharImageWriterType::Pointer      writer =  CharImageWriterType::New();
  // Internal2CharCastFilterType::Pointer  caster =  Internal2CharCastFilterType::New();
  // writer->SetFileName( deformedImageFileName );
  // caster->SetInput( warper->GetOutput() );
  // writer->SetInput( caster->GetOutput()   );
  // writer->Update();

  // if (isCompressed)
  //	CompressDeformationField2Short(deformationFieldFileName);
  return 0;
}
} // namespace Statistics
} // namespace itk

#endif
