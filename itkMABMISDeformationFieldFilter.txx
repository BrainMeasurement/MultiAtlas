#ifndef __itkMABMISDeformationFieldFilter_txx
#define __itkMABMISDeformationFieldFilter_txx

#include "itkMABMISDeformationFieldFilter.h"

using namespace std;

namespace itk
{
namespace Statistics
{

template <class TInputImage, class TOutputImage>
MABMISDeformationFieldFilter<TInputImage, TOutputImage>
::MABMISDeformationFieldFilter()
{
	imgoperator = ImageOperationType::New();	
}
template <class TInputImage, class TOutputImage>
MABMISDeformationFieldFilter<TInputImage, TOutputImage>
::~MABMISDeformationFieldFilter()
{
}

template <class TInputImage, class TOutputImage>
void
MABMISDeformationFieldFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
   Superclass::PrintSelf( os, indent );
}


template <class TInputImage, class TOutputImage>
int
MABMISDeformationFieldFilter<TInputImage, TOutputImage>
::ReadDeformationField(std::string filename, DeformationFieldType::Pointer &deformationfield)
{
	DeformationFieldReaderType::Pointer deformationFieldReader = DeformationFieldReaderType::New();
	deformationFieldReader->SetFileName( filename.c_str() );
	try
	{
		deformationFieldReader->Update();
	}
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << err << std::endl; 
		return -1;
	} 
	deformationfield = deformationFieldReader->GetOutput();
	return 0;
}

template <class TInputImage, class TOutputImage>
void
MABMISDeformationFieldFilter<TInputImage, TOutputImage>
::ComposeDeformationFieldsAndSave(std::string inputDeformationFieldFileName, std::string deformationFieldFileName, std::string composedDeformationFieldFileName)
{
	DeformationFieldType::Pointer inputDeformationField = DeformationFieldType::New();
	ReadDeformationField(inputDeformationFieldFileName, inputDeformationField);

	DeformationFieldType::Pointer deformationField = DeformationFieldType::New();
	ReadDeformationField(deformationFieldFileName, deformationField);

	DeformationFieldType::Pointer composedDeformationField = DeformationFieldType::New();

	ComposeDeformationFields(inputDeformationField, deformationField, composedDeformationField);

	WriteDeformationField(composedDeformationFieldFileName, composedDeformationField);

	return;
}


// compose deformation fields
template <class TInputImage, class TOutputImage>
void
MABMISDeformationFieldFilter<TInputImage, TOutputImage>
::ComposeDeformationFields(DeformationFieldType::Pointer input, DeformationFieldType::Pointer deformationField, DeformationFieldType::Pointer &composedDeformationField)
{
	WarpVectorFilterType::Pointer  vectorWarper = WarpVectorFilterType::New();
	vectorWarper->SetInput( input );
	vectorWarper->SetDeformationField( deformationField );
	vectorWarper->SetOutputOrigin(deformationField->GetOrigin());
	vectorWarper->SetOutputSpacing(deformationField->GetSpacing());
	vectorWarper->SetOutputDirection(deformationField->GetDirection());

	AddImageFilterType::Pointer addImageSum = AddImageFilterType::New();
	addImageSum->SetInput1(vectorWarper->GetOutput());
	addImageSum->SetInput2(deformationField);
	try
	{
		addImageSum->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr << err << std::endl;
		return;
	}
	composedDeformationField = addImageSum->GetOutput();
	composedDeformationField->DisconnectPipeline();
	return;

}

template <class TInputImage, class TOutputImage>
void
MABMISDeformationFieldFilter<TInputImage, TOutputImage>
::ApplyDeformationField(InternalImageType::Pointer movingImage, DeformationFieldType::Pointer deformationField,InternalImageType::Pointer &deformedImage, bool isLinearInterpolator)
{
	if (isLinearInterpolator)
	{
		InternalLinearInterpolatorType::Pointer interpolator = InternalLinearInterpolatorType::New();
		InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();

		warper->SetInput( movingImage );
		warper->SetInterpolator( interpolator );
		warper->SetOutputSpacing( movingImage->GetSpacing() );
		warper->SetOutputOrigin( movingImage->GetOrigin() );
		warper->SetOutputDirection( movingImage->GetDirection() );
		warper->SetDeformationField( deformationField );
		try
		{
			warper->Update();
		}
		catch( itk::ExceptionObject & err ) 
		{ 
			std::cerr << "ExceptionObject caught !" << std::endl; 
			std::cerr << err << std::endl; 
		} 
		deformedImage = warper->GetOutput();
		deformedImage->DisconnectPipeline();
	}
	else
	{
		InternalNNInterpolatorType::Pointer interpolator = InternalNNInterpolatorType::New();
		InternalWarpFilterType::Pointer warper = InternalWarpFilterType::New();

		warper->SetInput( movingImage );
		warper->SetInterpolator( interpolator );
		warper->SetOutputSpacing( movingImage->GetSpacing() );
		warper->SetOutputOrigin( movingImage->GetOrigin() );
		warper->SetOutputDirection( movingImage->GetDirection() );
		warper->SetDeformationField( deformationField );
		try
		{
			warper->Update();
		}
		catch( itk::ExceptionObject & err ) 
		{ 
			std::cerr << "ExceptionObject caught !" << std::endl; 
			std::cerr << err << std::endl; 
		} 
		deformedImage = warper->GetOutput();
		deformedImage->DisconnectPipeline();
	}

}


// apply deformation field on image and write deformed image
template <class TInputImage, class TOutputImage>
void
MABMISDeformationFieldFilter<TInputImage, TOutputImage>
::ApplyDeformationFieldAndWriteWithFileNames(string movingImageName, string deformationFieldFileName,  string deformedImageName, bool isLinearInterpolator)
{

	typename ImageOperationType::Pointer imageoperator = ImageOperationType::New();

	DeformationFieldType::Pointer deformationField = DeformationFieldType::New();
	ReadDeformationField(deformationFieldFileName, deformationField);

	InternalImageType::Pointer movingImage = InternalImageType::New();
	imageoperator->ReadImage(movingImageName,movingImage);

	InternalImageType::Pointer deformedImage = InternalImageType::New();

	ApplyDeformationField(movingImage, deformationField, deformedImage, isLinearInterpolator);

	imageoperator->WriteImage(deformedImageName, deformedImage);

}

template <class TInputImage, class TOutputImage>
void
MABMISDeformationFieldFilter<TInputImage, TOutputImage>
::ApplyDeformationFieldAndWriteWithTypeWithFileNames(std::string  movingImageFileName, std::string deformationFieldFileName, std::string deformedImageFileName, bool isLinear)
{
  
	itk::ImageIOBase::Pointer imageIO;
	try
	{
		imageIO = itk::ImageIOFactory::CreateImageIO(movingImageFileName.c_str(), itk::ImageIOFactory::ReadMode);
		if ( imageIO )
		{
			imageIO->SetFileName(movingImageFileName);
			imageIO->ReadImageInformation();
		}
		else
		{
			std::cout << "Could not read the image information of "<< movingImageFileName <<"." << std::endl;
			exit( EXIT_FAILURE );
		}
	}
	catch( itk::ExceptionObject& err )
	{
		std::cout << "Could not read the image information of "<< movingImageFileName <<"." << std::endl;
		std::cout << err << std::endl;
		exit( EXIT_FAILURE );
	}
	//{UNKNOWNCOMPONENTTYPE,UCHAR,CHAR,USHORT,SHORT,UINT,INT,ULONG,LONG,FLOAT,DOUBLE} IOComponentType;
	//                    0     1    2      3     4    5   6     7    8     9     10
	int input_type = imageIO->GetComponentType(); // 9:float, 10:double

	//
	DeformationFieldType::Pointer deformationField = DeformationFieldType::New();
	ReadDeformationField(deformationFieldFileName, deformationField);

	InternalImageType::Pointer movingImage = InternalImageType::New();
	imgoperator->ReadImage(movingImageFileName,movingImage);

	InternalImageType::Pointer deformedImage = InternalImageType::New();

	ApplyDeformationField(movingImage, deformationField, deformedImage, isLinear);

	if (input_type == 1) // UCHAR
	{		
		imgoperator->WriteImageUCHAR(deformedImageFileName, deformedImage);
	}
	else if (input_type == 4) // SHORT
	{
		imgoperator->WriteImageSHORT(deformedImageFileName, deformedImage);
	}
	else if (input_type == 6) // INT
	{
		imgoperator->WriteImageINT(deformedImageFileName, deformedImage);
	}
	else if (input_type == 9) // FLOAT
	{
		imgoperator->WriteImageFLOAT(deformedImageFileName, deformedImage);
	}
	else if (input_type == 10) // DOUBLE
	{
		imgoperator->WriteImageFLOAT(deformedImageFileName, deformedImage);
	}
	else
	{
		imgoperator->WriteImageFLOAT(deformedImageFileName, deformedImage);
	}
  
}
// end
// 


template <class TInputImage, class TOutputImage>
void
MABMISDeformationFieldFilter<TInputImage, TOutputImage>
::WriteDeformationField(std::string  filename, DeformationFieldType::Pointer deformationfield)
{
	DeformationFieldWriterType::Pointer deformationFieldWriter = DeformationFieldWriterType::New();
	deformationFieldWriter->SetFileName( filename );
	deformationFieldWriter->SetInput (deformationfield);
	deformationFieldWriter->SetUseCompression( false );
	try
	{
		deformationFieldWriter->Update();
	}
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cerr << err << std::endl; 
		return;
	} 
	return;
}

template <class TInputImage, class TOutputImage>
void
MABMISDeformationFieldFilter<TInputImage, TOutputImage>
::DownResampleDeformationField(std::string deformationFieldFileName, std::string resampledDeformationFieldFileName, int sampleRate)
{
  

	int rx, ry, rz;
	rx = sampleRate;
	ry = sampleRate;
	rz = sampleRate;

	DeformationFieldType::Pointer originImage = 0;
	
	ReadDeformationField(deformationFieldFileName, originImage);

	int im_x, im_y, im_z;
	int im_xn, im_yn, im_zn;
	DeformationFieldType::SizeType im_size = originImage->GetLargestPossibleRegion().GetSize();
	im_x = im_size[0]; im_xn = (im_x-1)/rx+1;
	im_y = im_size[1]; im_yn = (im_y-1)/ry+1;
	im_z = im_size[2]; im_zn = (im_z-1)/rz+1;



	//cerr << "rx: " << rx << "ry: " << ry << "rz: " << rz << endl;
	//cerr << "im_x: " << im_x << "im_y: " << im_y << "im_z: " << im_z << endl;
	

	// create an array to store original image
	float*** originImageRawX = new float**[im_z];
	float*** originImageRawY = new float**[im_z];
	float*** originImageRawZ = new float**[im_z];
	for (int k = 0; k < im_z; k++) 
	{
		originImageRawX[k] = new float*[im_y];
		originImageRawY[k] = new float*[im_y];
		originImageRawZ[k] = new float*[im_y];
		for (int j = 0; j < im_y; j++)
		{
			originImageRawX[k][j] = new float[im_x];
			originImageRawY[k][j] = new float[im_x];
			originImageRawZ[k][j] = new float[im_x];
			for (int i = 0; i < im_x; i++)
			{
				originImageRawX[k][j][i] = 0.0;
				originImageRawY[k][j][i] = 0.0;
				originImageRawZ[k][j][i] = 0.0;
			}
		}
	}

	// create an array to store sampled image
	float*** sampledImageRawX = new float**[im_zn];
	float*** sampledImageRawY = new float**[im_zn];
	float*** sampledImageRawZ = new float**[im_zn];
	for (int k = 0; k < im_zn; k++) 
	{
		sampledImageRawX[k] = new float*[im_yn];
		sampledImageRawY[k] = new float*[im_yn];
		sampledImageRawZ[k] = new float*[im_yn];
		for (int j = 0; j < im_yn; j++)
		{
			sampledImageRawX[k][j] = new float[im_xn];
			sampledImageRawY[k][j] = new float[im_xn];
			sampledImageRawZ[k][j] = new float[im_xn];
			for (int i = 0; i < im_xn; i++)
			{
				sampledImageRawX[k][j][i] = 0.0;
				sampledImageRawY[k][j][i] = 0.0;
				sampledImageRawZ[k][j][i] = 0.0;
			}
		}
	}

	// load original image
	DeformationFieldIteratorType itOrigin(originImage, originImage->GetLargestPossibleRegion() );
	VectorPixelType vectorPixel;
	DeformationFieldType::IndexType idx;
	int idx_x, idx_y, idx_z;

	for (itOrigin.GoToBegin(); !itOrigin.IsAtEnd(); ++itOrigin)
	{
		vectorPixel = itOrigin.Get();
		idx = itOrigin.GetIndex();
		idx_x = idx.GetElement(0);
		idx_y = idx.GetElement(1);
		idx_z = idx.GetElement(2);

		//pixel = itOrigin.Get();
		//idx = itOrigin.GetIndex();
		originImageRawX[idx[2]][idx[1]][idx[0]] = vectorPixel.GetElement(0);
		originImageRawY[idx[2]][idx[1]][idx[0]] = vectorPixel.GetElement(1);
		originImageRawZ[idx[2]][idx[1]][idx[0]] = vectorPixel.GetElement(2);
	}


	// resample image
	for (int k = 0; k < im_zn; k++) 
	{
		for (int j = 0; j < im_yn; j++)
		{
			for (int i = 0; i < im_xn; i++)
			{
				sampledImageRawX[k][j][i] = originImageRawX[k*rz][j*ry][i*rx]/rx;
				sampledImageRawY[k][j][i] = originImageRawY[k*rz][j*ry][i*rx]/ry;
				sampledImageRawZ[k][j][i] = originImageRawZ[k*rz][j*ry][i*rx]/rz;
			}
		}
	}

	// create the resampled image
	DeformationFieldType::Pointer sampledImage = DeformationFieldType::New();
	DeformationFieldType::IndexType start;
	start[0] = 0;start[1] = 0;start[2] = 0;
	DeformationFieldType::SizeType size;
	size[0] = im_xn;size[1] = im_yn;size[2] = im_zn;
	DeformationFieldType::RegionType region;
	region.SetSize (size);
	region.SetIndex (start);

	sampledImage -> SetRegions (region);
	// sampledImage -> SetRegions (originImage->GetLargestPossibleRegion());
	sampledImage -> SetSpacing (originImage->GetSpacing());
	sampledImage -> SetDirection (originImage->GetDirection());
	sampledImage -> SetOrigin (originImage->GetOrigin());

	sampledImage -> Allocate();

	DeformationFieldIteratorType itSampled(sampledImage, sampledImage->GetLargestPossibleRegion() );
	for (itSampled.GoToBegin(); !itSampled.IsAtEnd(); ++itSampled)
	{
		idx = itSampled.GetIndex();
		vectorPixel.SetElement(0, sampledImageRawX[idx[2]][idx[1]][idx[0]]);
		vectorPixel.SetElement(1, sampledImageRawY[idx[2]][idx[1]][idx[0]]);
		vectorPixel.SetElement(2, sampledImageRawZ[idx[2]][idx[1]][idx[0]]);
		itSampled.Set(vectorPixel);
	}

	WriteDeformationField(resampledDeformationFieldFileName, sampledImage);


	// delete newed variables
	for (int k = 0; k < im_z; k++) 
	{
		for (int j = 0; j < im_y; j++)
		{
			delete[] originImageRawX[k][j];
			delete[] originImageRawY[k][j];
			delete[] originImageRawZ[k][j];
		}
		delete[] originImageRawX[k];
		delete[] originImageRawY[k];
		delete[] originImageRawZ[k];
	}
	delete[] originImageRawX;
	delete[] originImageRawY;
	delete[] originImageRawZ;
	for (int k = 0; k < im_zn; k++) 
	{
		for (int j = 0; j < im_yn; j++)
		{
			delete[] sampledImageRawX[k][j];
			delete[] sampledImageRawY[k][j];
			delete[] sampledImageRawZ[k][j];
		}
		delete[] sampledImageRawX[k];
		delete[] sampledImageRawY[k];
		delete[] sampledImageRawZ[k];
	}
	delete[] sampledImageRawX;
	delete[] sampledImageRawY;
	delete[] sampledImageRawZ;
  
	return;
}

template <class TInputImage, class TOutputImage>
void
MABMISDeformationFieldFilter<TInputImage, TOutputImage>
::UpResampleDeformationField(std::string deformationFieldFileName, std::string  resampledDeformationFieldFileName, int sampleRate)
{
  
	// upsample
	int rx, ry, rz;
	rx = sampleRate;
	ry = sampleRate;
	rz = sampleRate;

	DeformationFieldType::Pointer originImage = 0;

	ReadDeformationField(deformationFieldFileName, originImage);


	int im_x, im_y, im_z;
	int im_xn, im_yn, im_zn;
	DeformationFieldType::SizeType im_size = originImage->GetLargestPossibleRegion().GetSize();
	im_x = im_size[0];
	im_y = im_size[1];
	im_z = im_size[2];
	//im_xn = im_x*rx; im_yn = im_y*ry; im_zn = im_z*rz;
	im_xn = m_Imx; im_yn = m_Imy; im_zn = m_Imz;

	// create an array to store original image
	float*** originImageRawX = new float**[im_z];
	float*** originImageRawY = new float**[im_z];
	float*** originImageRawZ = new float**[im_z];
	for (int k = 0; k < im_z; k++) 
	{
		originImageRawX[k] = new float*[im_y];
		originImageRawY[k] = new float*[im_y];
		originImageRawZ[k] = new float*[im_y];
		for (int j = 0; j < im_y; j++)
		{
			originImageRawX[k][j] = new float[im_x];
			originImageRawY[k][j] = new float[im_x];
			originImageRawZ[k][j] = new float[im_x];
			for (int i = 0; i < im_x; i++)
			{
				originImageRawX[k][j][i] = 0.0;
				originImageRawY[k][j][i] = 0.0;
				originImageRawZ[k][j][i] = 0.0;
			}
		}
	}

	// create an array to store sampled image
	float*** sampledImageRawX = new float**[im_zn];
	float*** sampledImageRawY = new float**[im_zn];
	float*** sampledImageRawZ = new float**[im_zn];
	for (int k = 0; k < im_zn; k++) 
	{
		sampledImageRawX[k] = new float*[im_yn];
		sampledImageRawY[k] = new float*[im_yn];
		sampledImageRawZ[k] = new float*[im_yn];
		for (int j = 0; j < im_yn; j++)
		{
			sampledImageRawX[k][j] = new float[im_xn];
			sampledImageRawY[k][j] = new float[im_xn];
			sampledImageRawZ[k][j] = new float[im_xn];
			for (int i = 0; i < im_xn; i++)
			{
				sampledImageRawX[k][j][i] = 0.0;
				sampledImageRawY[k][j][i] = 0.0;
				sampledImageRawZ[k][j][i] = 0.0;
			}
		}
	}

	// load original image
	DeformationFieldIteratorType itOrigin(originImage, originImage->GetLargestPossibleRegion() );
	VectorPixelType vectorPixel;
	DeformationFieldType::IndexType idx;
	int idx_x, idx_y, idx_z;

	for (itOrigin.GoToBegin(); !itOrigin.IsAtEnd(); ++itOrigin)
	{
		vectorPixel = itOrigin.Get();
		idx = itOrigin.GetIndex();
		idx_x = idx.GetElement(0);
		idx_y = idx.GetElement(1);
		idx_z = idx.GetElement(2);

		//pixel = itOrigin.Get();
		//idx = itOrigin.GetIndex();
		originImageRawX[idx[2]][idx[1]][idx[0]] = vectorPixel.GetElement(0);
		originImageRawY[idx[2]][idx[1]][idx[0]] = vectorPixel.GetElement(1);
		originImageRawZ[idx[2]][idx[1]][idx[0]] = vectorPixel.GetElement(2);
	}


	// resample image
	for (int k = 0; k < im_zn; k++) 
	{
		for (int j = 0; j < im_yn; j++)
		{
			for (int i = 0; i < im_xn; i++)
			{
				sampledImageRawX[k][j][i] = originImageRawX[k/rz][j/ry][i/rx]*rx;
				sampledImageRawY[k][j][i] = originImageRawY[k/rz][j/ry][i/rx]*ry;
				sampledImageRawZ[k][j][i] = originImageRawZ[k/rz][j/ry][i/rx]*rz;
			}
		}
	}

	// create the resampled image
	DeformationFieldType::Pointer sampledImage = DeformationFieldType::New();
	DeformationFieldType::IndexType start;
	start[0] = 0;start[1] = 0;start[2] = 0;
	DeformationFieldType::SizeType size;
	size[0] = im_xn;size[1] = im_yn;size[2] = im_zn;
	DeformationFieldType::RegionType region;
	region.SetSize (size);
	region.SetIndex (start);

	sampledImage -> SetRegions (region);
	// sampledImage -> SetRegions (originImage->GetLargestPossibleRegion());
	sampledImage -> SetSpacing (originImage->GetSpacing());
	sampledImage -> SetDirection (originImage->GetDirection());
	sampledImage -> SetOrigin (originImage->GetOrigin());

	sampledImage -> Allocate();

	DeformationFieldIteratorType itSampled(sampledImage, sampledImage->GetLargestPossibleRegion() );
	for (itSampled.GoToBegin(); !itSampled.IsAtEnd(); ++itSampled)
	{
		idx = itSampled.GetIndex();
		vectorPixel.SetElement(0, sampledImageRawX[idx[2]][idx[1]][idx[0]]);
		vectorPixel.SetElement(1, sampledImageRawY[idx[2]][idx[1]][idx[0]]);
		vectorPixel.SetElement(2, sampledImageRawZ[idx[2]][idx[1]][idx[0]]);
		itSampled.Set(vectorPixel);
	}

	WriteDeformationField(resampledDeformationFieldFileName, sampledImage);


	// delete newed variables
	for (int k = 0; k < im_z; k++) 
	{
		for (int j = 0; j < im_y; j++)
		{
			delete[] originImageRawX[k][j];
			delete[] originImageRawY[k][j];
			delete[] originImageRawZ[k][j];
		}
		delete[] originImageRawX[k];
		delete[] originImageRawY[k];
		delete[] originImageRawZ[k];
	}
	delete[] originImageRawX;
	delete[] originImageRawY;
	delete[] originImageRawZ;
	for (int k = 0; k < im_zn; k++) 
	{
		for (int j = 0; j < im_yn; j++)
		{
			delete[] sampledImageRawX[k][j];
			delete[] sampledImageRawY[k][j];
			delete[] sampledImageRawZ[k][j];
		}
		delete[] sampledImageRawX[k];
		delete[] sampledImageRawY[k];
		delete[] sampledImageRawZ[k];
	}
	delete[] sampledImageRawX;
	delete[] sampledImageRawY;
	delete[] sampledImageRawZ;

  
	return;
}






template <class TInputImage, class TOutputImage>
void
MABMISDeformationFieldFilter<TInputImage, TOutputImage>
::InverseDeformationField3D(DeformationFieldType::Pointer deformationField, 
							   DeformationFieldType::Pointer &deformationFieldInverse)
{
	int SHIFT = 2;
	float OUTSIDE = 0;//100000.0;//0;
	int samplenum = 1;
	float interval = 1.0/(2*samplenum+1);
	unsigned int i, j, k;
	int x, y, z;
	float ii, jj, kk;
	float mdl_subvoxelx, mdl_subvoxely, mdl_subvoxelz;
	float disp_subvoxelx, disp_subvoxely, disp_subvoxelz;
	mdl_subvoxelx = 0.0; mdl_subvoxely = 0.0; mdl_subvoxelz = 0.0; disp_subvoxelx = 0.0; disp_subvoxely = 0.0; disp_subvoxelz = 0.0;
	DeformationFieldIteratorType dfNewIt ( deformationFieldInverse, deformationFieldInverse->GetRequestedRegion() );
	DeformationFieldIteratorType dfIt ( deformationField, deformationField->GetRequestedRegion() );

	// initialize three 3D float matrix to store deformation field
	//unsigned int image_size = deformationField->GetRequestedRegion().GetSize()[0];
	unsigned int x_size = deformationField->GetRequestedRegion().GetSize()[0];
	unsigned int y_size = deformationField->GetRequestedRegion().GetSize()[1];
	unsigned int z_size = deformationField->GetRequestedRegion().GetSize()[2];
	
	float*** dfx = new float**[x_size];
	float*** dfy = new float**[x_size];
	float*** dfz = new float**[x_size];
	for(i = 0; i < x_size; i++) 
	{
		dfx[i] = new float*[y_size];
		dfy[i] = new float*[y_size];
		dfz[i] = new float*[y_size];
	}
	for(i = 0; i < x_size; i++) 
	{
		for (j = 0; j < y_size; j++)
		{
			dfx[i][j] = new float[z_size];
			dfy[i][j] = new float[z_size];
			dfz[i][j] = new float[z_size];
		}
	}

	// load deformationFieldBA into 3 3D matrix dfx and dfy and dfz
	VectorPixelType vectorPixel;
	DeformationFieldType::IndexType idx;
	//for (i = 0, j = 0, k = 0, dfIt.GoToBegin(); !dfIt.IsAtEnd(); ++dfIt)
	for (dfIt.GoToBegin(); !dfIt.IsAtEnd(); ++dfIt)
	{
		vectorPixel = dfIt.Get();
		idx = dfIt.GetIndex();
		dfx[idx[0]][idx[1]][idx[2]] = vectorPixel.GetElement(0);
		dfy[idx[0]][idx[1]][idx[2]] = vectorPixel.GetElement(1);
		dfz[idx[0]][idx[1]][idx[2]] = vectorPixel.GetElement(2);
		//dfx[i][j][k] = vectorPixel.GetElement(0);
		//dfy[i][j][k] = vectorPixel.GetElement(1);
		//dfz[i][j][k] = vectorPixel.GetElement(2);
		//k++;
		//if (k == z_size)
		//{
		//	j++;
		//	k = 0;
		//	if (j == image_size)
		//	{
		//		i++;
		//		j = 0;
		//	}
		//}
		//i++;
		//if (i == image_size)
		//{
		//	j++;
		//	i = 0;
		//	if (j == image_size)
		//	{
		//		k++;
		//		j = 0;
		//	}
		//}
		//if ((i == image_size-1) && (j == image_size-1) && (k == z_size-1))
		//{
		//	std::cerr << "read end!" << std::endl;
		//}
	}

	// allocate some internal data matrices, weights matrix and	enlarged inverse df matrix with borders
	float*** totalweights = new float**[x_size+2*SHIFT];
	float***	    rdfbx = new float**[x_size+2*SHIFT];
	float***	    rdfby = new float**[x_size+2*SHIFT];
	float***	    rdfbz = new float**[x_size+2*SHIFT];
	for(i = 0; i < x_size+2*SHIFT; i++)
	{
		totalweights[i] = new float*[y_size+2*SHIFT];
			   rdfbx[i] = new float*[y_size+2*SHIFT];
			   rdfby[i] = new float*[y_size+2*SHIFT];
			   rdfbz[i] = new float*[y_size+2*SHIFT];
	}
	for(i = 0; i < x_size+2*SHIFT; i++) 
	{
		for (j = 0; j < y_size+2*SHIFT; j++)
		{
			totalweights[i][j] = new float[z_size+2*SHIFT];
				   rdfbx[i][j] = new float[z_size+2*SHIFT];
				   rdfby[i][j] = new float[z_size+2*SHIFT];
				   rdfbz[i][j] = new float[z_size+2*SHIFT];
		}
	}

	// initialize these matrices
	for (i = 0; i < x_size+2*SHIFT; i++)
	{	
		for (j = 0; j < y_size+2*SHIFT; j++)
		{	
			for (k = 0; k < z_size+2*SHIFT; k++)
			{
				totalweights[i][j][k] = 0.0;	
					   rdfbx[i][j][k] = 0.0;
					   rdfby[i][j][k] = 0.0;
					   rdfbz[i][j][k] = 0.0;
			}
		}
	}

	// estimating
	for (i = 0; i < x_size; i++)
	{
		// std::cerr << i << ", ";
		for (j = 0; j < y_size; j++)
		{
			// std::cerr << "(" << i << ", " << j << "), " ;
			for (k = 0; k < z_size; k++)
			{
				for (x = -samplenum; x <= samplenum; x++)
				{
					for (y = -samplenum; y <= samplenum; y++)
					{
						for (z = -samplenum; z <= samplenum; z++)
						{
							mdl_subvoxelx = x*interval + i;
							mdl_subvoxely = y*interval + j;
							mdl_subvoxelz = z*interval + k;

							// call interpolateDisplacement
							{
								int ni, nj, nk, nip1, njp1, nkp1;
								float b, c, d, b1, c1, d1;
								ni = (int)mdl_subvoxelx; nip1 = ni + 1;
								nj = (int)mdl_subvoxely; njp1 = nj + 1;
								nk = (int)mdl_subvoxelz; nkp1 = nk + 1;

								if((ni>=0)&&(ni<(int)x_size-1)&&(nj>=0)&&(nj<(int)y_size-1)&&(nk>=0)&&(nk<(int)z_size-1))
								{
									b = mdl_subvoxelx - ni; b1 = 1.0-b;
									c = mdl_subvoxely - nj; c1 = 1.0-c;
									d = mdl_subvoxelz - nk; d1 = 1.0-d;

									disp_subvoxelx = (d1*( dfx[ni][nj][nk]*(b1*c1) + dfx[nip1][nj][nk]*(b*c1)+ 
										dfx[ni][njp1][nk]*(b1*c)+ dfx[nip1][njp1][nk]*(b*c) ) +
										d*(dfx[ni][nj][nkp1]*(b1*c1) + dfx[nip1][nj][nkp1]*(b*c1)+ 
										dfx[ni][njp1][nkp1]*(b1*c)+ dfx[nip1][njp1][nkp1]*(b*c)))
										/(d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c))) ; 
									disp_subvoxely = (d1*( dfy[ni][nj][nk]*(b1*c1) + dfy[nip1][nj][nk]*(b*c1)+ 
										dfy[ni][njp1][nk]*(b1*c)+ dfy[nip1][njp1][nk]*(b*c) ) +
										d*(dfy[ni][nj][nkp1]*(b1*c1) + dfy[nip1][nj][nkp1]*(b*c1)+ 
										dfy[ni][njp1][nkp1]*(b1*c)+ dfy[nip1][njp1][nkp1]*(b*c)))
										/(d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c))) ; 
									disp_subvoxelz = (d1*( dfz[ni][nj][nk]*(b1*c1) + dfz[nip1][nj][nk]*(b*c1)+ 
										dfz[ni][njp1][nk]*(b1*c)+ dfz[nip1][njp1][nk]*(b*c) ) +
										d*(dfz[ni][nj][nkp1]*(b1*c1) + dfz[nip1][nj][nkp1]*(b*c1)+ 
										dfz[ni][njp1][nkp1]*(b1*c)+ dfz[nip1][njp1][nkp1]*(b*c)))
										/(d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c))) ; 
									//disp_subvoxely = ( dfy[ni][nj]*(b1*c1) + 
									//	dfy[nip1][nj]*(b*c1)+ 
									//	dfy[ni][njp1]*(b1*c)+ 
									//	dfy[nip1][njp1]*(b*c) )/( (b1*c1)+(b*c1)+(b1*c)+(b*c) ) ; 

								}
								else if (((ni==(int)x_size-1)&&(nj>=0)&&(nj<(int)y_size-1)  &&(nk>=0)&&(nk<(int)z_size-1))
									   ||((ni>=0)&&(ni<(int)x_size-1) &&(nj==(int)y_size-1) &&(nk>=0)&&(nk<(int)z_size-1))
									   ||((ni>=0)&&(ni<(int)x_size-1) &&(nj>=0)&&(nj<(int)y_size-1)  &&(nk=(int)z_size-1)))
								{
									disp_subvoxelx = dfx[ni][nj][nk];
									disp_subvoxely = dfy[ni][nj][nk];
									disp_subvoxelz = dfz[ni][nj][nk];
								}

							}

							ii = mdl_subvoxelx + disp_subvoxelx;
							jj = mdl_subvoxely + disp_subvoxely;
							kk = mdl_subvoxelz + disp_subvoxelz;

							// call iterativeEstimate
							{
								int ni,nj,nk,nip1,njp1,nkp1;
								float b,c,d,b1,c1,d1,combined,weight;
								ii += SHIFT; jj += SHIFT; kk += SHIFT;
								ni = (int)ii; nip1 = ni+1;
								nj = (int)jj; njp1 = nj+1;
								nk = (int)kk; nkp1 = nk+1;
								if(ni>=0 && ni<(int)x_size+2*SHIFT-1  &&  nj>=0 && nj<(int)y_size+2*SHIFT-1 &&  nk>=0 && nk<(int)z_size+2*SHIFT-1)
								{
									b = ii-ni ;        b1 = 1.0-b ;
									c = jj-nj ;        c1 = 1.0-c ;
									d = kk-nk ;        d1 = 1.0-d ;
									combined = d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c));

									weight = d1*(b1*c1)/combined ;
									totalweights[ni][nj][nk] += weight ;
										   rdfbx[ni][nj][nk] += disp_subvoxelx*weight ;
										   rdfby[ni][nj][nk] += disp_subvoxely*weight ;
										   rdfbz[ni][nj][nk] += disp_subvoxelz*weight ;

									weight = d1*(b*c1)/combined ;
									totalweights[nip1][nj][nk] += weight ;
										   rdfbx[nip1][nj][nk] += disp_subvoxelx*weight ;
										   rdfby[nip1][nj][nk] += disp_subvoxely*weight ;
										   rdfbz[nip1][nj][nk] += disp_subvoxelz*weight ;

									weight = d1*(b1*c)/combined ;
									totalweights[ni][njp1][nk] += weight ;
										   rdfbx[ni][njp1][nk] += disp_subvoxelx*weight ;
										   rdfby[ni][njp1][nk] += disp_subvoxely*weight ;
										   rdfbz[ni][njp1][nk] += disp_subvoxelz*weight ;

									weight = d1*(b*c)/combined ;
									totalweights[nip1][njp1][nk]  += weight ;
										   rdfbx[nip1][njp1][nk] += disp_subvoxelx*weight ;
										   rdfby[nip1][njp1][nk] += disp_subvoxely*weight ;
										   rdfbz[nip1][njp1][nk] += disp_subvoxelz*weight ;

									weight = d*(b1*c1)/combined ;
									totalweights[ni][nj][nkp1] += weight ;
										   rdfbx[ni][nj][nkp1] += disp_subvoxelx*weight ;
										   rdfby[ni][nj][nkp1] += disp_subvoxely*weight ;
										   rdfbz[ni][nj][nkp1] += disp_subvoxelz*weight ;

									weight = d*(b*c1)/combined ;
									totalweights[nip1][nj][nkp1] += weight ;
										   rdfbx[nip1][nj][nkp1] += disp_subvoxelx*weight ;
										   rdfby[nip1][nj][nkp1] += disp_subvoxely*weight ;
										   rdfbz[nip1][nj][nkp1] += disp_subvoxelz*weight ;

									weight = d*(b1*c)/combined ;
									totalweights[ni][njp1][nkp1] += weight ;
										   rdfbx[ni][njp1][nkp1] += disp_subvoxelx*weight ;
										   rdfby[ni][njp1][nkp1] += disp_subvoxely*weight ;
										   rdfbz[ni][njp1][nkp1] += disp_subvoxelz*weight ;

									weight = d*(b*c)/combined ;
									totalweights[nip1][njp1][nkp1] += weight ;
										   rdfbx[nip1][njp1][nkp1] += disp_subvoxelx*weight ;
										   rdfby[nip1][njp1][nkp1] += disp_subvoxely*weight ;
										   rdfbz[nip1][njp1][nkp1] += disp_subvoxelz*weight ;
								}

							}
						}// end for z
					}// end for y
				}// end for z
			}// end for k
		}// end for j
	}// end for i

	// allocate inverse deformation field
	float*** rdfx = new float**[x_size];
	float*** rdfy = new float**[x_size];
	float*** rdfz = new float**[x_size];
	for(i = 0; i<x_size; i++)
	{
		rdfx[i] = new float*[y_size];
		rdfy[i] = new float*[y_size];
		rdfz[i] = new float*[y_size];
	}
	for(i = 0; i<x_size; i++) 
	{
		for (j = 0; j < y_size; j++)
		{
			rdfx[i][j] = new float[z_size];
			rdfy[i][j] = new float[z_size];
			rdfz[i][j] = new float[z_size];
		}
	}
	int count_outside = 0;
	// normalize the enlarged rdfb to rdf
	for (i = 0; i< x_size; i++)
	{
		for (j = 0; j < y_size; j++)
		{
			for (k = 0; k < z_size; k++)
			{
				if (totalweights[i+SHIFT][j+SHIFT][k+SHIFT]>0)
				{
					rdfx[i][j][k] = rdfbx[i+SHIFT][j+SHIFT][k+SHIFT]/(-totalweights[i+SHIFT][j+SHIFT][k+SHIFT]);
					rdfy[i][j][k] = rdfby[i+SHIFT][j+SHIFT][k+SHIFT]/(-totalweights[i+SHIFT][j+SHIFT][k+SHIFT]);
					rdfz[i][j][k] = rdfbz[i+SHIFT][j+SHIFT][k+SHIFT]/(-totalweights[i+SHIFT][j+SHIFT][k+SHIFT]);
				}
				else
				{
					rdfx[i][j][k] = OUTSIDE;
					rdfy[i][j][k] = OUTSIDE;
					rdfz[i][j][k] = OUTSIDE;
					count_outside++;
				}
			}
		}
	}
	//for (i = 0, j = 0, k = 0, dfNewIt.GoToBegin(); !dfNewIt.IsAtEnd(); ++dfNewIt)
	for (dfNewIt.GoToBegin(); !dfNewIt.IsAtEnd(); ++dfNewIt)
	{
		idx = dfNewIt.GetIndex();
		vectorPixel.SetElement(0,rdfx[idx[0]][idx[1]][idx[2]]);
		vectorPixel.SetElement(1,rdfy[idx[0]][idx[1]][idx[2]]);
		vectorPixel.SetElement(2,rdfz[idx[0]][idx[1]][idx[2]]);
		dfNewIt.Set(vectorPixel);
		//k++;
		//if (k == z_size)
		//{
		//	j++;
		//	k = 0;
		//	if (j == image_size)
		//	{
		//		i++;
		//		j = 0;
		//	}
		//}
		//i++;
		//if (i == image_size)
		//{
		//	j++;
		//	i = 0;
		//	if (j == image_size)
		//	{
		//		k++;
		//		j = 0;
		//	}
		//}
		//j++;
		//if (j==image_size)
		//{
		//	i++;
		//	j=0;
		//}
	}

	// free memory
	for(i = 0; i<x_size; i++)
	{
		for(j = 0; j<y_size; j++)
		{
			delete[] dfx[i][j];			delete[] dfy[i][j];			delete[] dfz[i][j];
			delete[] rdfx[i][j];		delete[] rdfy[i][j];		delete[] rdfz[i][j];
		}
		delete[] dfx[i];		delete[] dfy[i];		delete[] dfz[i];
		delete[] rdfx[i];		delete[] rdfy[i];		delete[] rdfz[i];
	}
	for(i = 0; i<x_size+2*SHIFT; i++)
	{
		for(j = 0; j<y_size+2*SHIFT; j++)
		{
			delete[] rdfbx[i][j];			delete[] rdfby[i][j];			delete[] rdfbz[i][j];
			delete[] totalweights[i][j];
		}
		delete[] rdfbx[i];		delete[] rdfby[i];		delete[] rdfbz[i];
		delete[] totalweights[i]; 
	}
	delete[] dfx;	delete[] dfy;	delete[] dfz;
	delete[] rdfx;	delete[] rdfy;	delete[] rdfz;
	delete[] rdfbx;	delete[] rdfby;	delete[] rdfbz;
	delete[] totalweights;

}




} // namespace Statistics
}  // namespace itk

#endif
