#ifndef __itkMABMISSimulateData_hxx
#define __itkMABMISSimulateData_hxx

#include "itkMABMISSimulateData.h"

//For definition of FILESEP
#include "itkMABMISAtlasXMLFile.h"

int numEigenVector = 4;         // t
int numSampleEachDirection = 4; // n  n^t total number of intermediate templates

namespace itk
{
namespace Statistics
{
template <class TInputImage, class TOutputImage>
MABMISSimulateData<TInputImage, TOutputImage>
::MABMISSimulateData()
{
  dfoperator = DeformationFieldOperationType::New();
  imgoperator = ImageOperationType::New();
  basicoperator = BasicOperationType::New();
}

template <class TInputImage, class TOutputImage>
MABMISSimulateData<TInputImage, TOutputImage>
::~MABMISSimulateData()
{
}

template <class TInputImage, class TOutputImage>
void
MABMISSimulateData<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}

template <class TInputImage, class TOutputImage>
int
MABMISSimulateData<TInputImage, TOutputImage>
::DoPCATraining(std::vector<std::string> deformationFieldFileNames, int numFiles,
                std::vector<std::string> allImgFileName, int root)
{
  bool doPCATraining = true;
  int  sampleRate = 2;
  bool doResampleDeformationField = true;
  int  size_x = 0; int size_y = 0; int size_z = 0;
  int  size_xn = 0; int size_yn = 0; int size_zn = 0;
  int  size_dfn = 0; // size of sub-sampled deformation field

  float c4[] = {-0.8416, -0.2533, 0.2533, 0.8416};

  float* c = NULL;

  if( numSampleEachDirection >= 3 && numSampleEachDirection <= 5 )
    {
    c = new float[numSampleEachDirection];
    for( int i = 0; i < numSampleEachDirection; ++i )
      {
      c[i] = c4[i];
      }
    }
  else
    {
    std::cerr << "The number of samples on each eigen vector direction should be between 3 and 5!" << std::endl;
    exit(1);
    }
  // some global parameters
  // SVD   (Vt == V^T,transpose; Si = inv(S), inverse)
  // D = U S Vt   =>  DV = US => DVSi = U

  ///////////////////////////////////////
  // initialize

  DeformationFieldType::Pointer curDeformationField = 0;

  dfoperator->ReadDeformationField(deformationFieldFileNames[0], curDeformationField);

  DeformationFieldType::SizeType im_size = curDeformationField->GetLargestPossibleRegion().GetSize();
  size_x = im_size[0]; size_y = im_size[1]; size_z = im_size[2];
  size_xn =
    (size_x - 1) / sampleRate + 1; size_yn = (size_y - 1) / sampleRate + 1; size_zn = (size_z - 1) / sampleRate + 1;
  size_dfn = (size_xn) * (size_yn) * (size_zn) * 3;

  df_spacing = curDeformationField->GetSpacing();
  df_direction = curDeformationField->GetDirection();
  df_origin = curDeformationField->GetOrigin();

  // initialize
  ////////////////////////////////////////
  // do PCA training
  ////////////////////////////////////////
  // std::cerr << "Resample deformation fields ..." << std::endl;
  vnl_matrix<float> df_eigenvector(size_dfn, numEigenVector);
  vnl_vector<float> df_eigenvalues(numEigenVector);
  ////////////////////////////////////////
  // down sample
  for( int i = 0; i < numFiles; ++i )
    {
    std::string deformationFieldFileName;
    deformationFieldFileName = deformationFieldFileNames[i];

    std::string resampledDeformationFieldFileName;
    std::string tempstring;
    tempstring = deformationFieldFileNames[i];

    tempstring.erase(tempstring.end() - 4, tempstring.end() );
    resampledDeformationFieldFileName.append(tempstring);
    resampledDeformationFieldFileName.append("_sub.nii.gz");

    if( doResampleDeformationField )
      {
      dfoperator->DownResampleDeformationField(deformationFieldFileName, resampledDeformationFieldFileName, sampleRate);
      }
    }

  vnl_matrix<float> df_mat(size_dfn, numFiles, 0.0);
  vnl_vector<float> df_mean(size_dfn, 0.0);

  if( doPCATraining )
    {
    // build all vectors from deformation fields
    for( int i = 0; i < numFiles; ++i )
      {
      std::string resampledDeformationFieldFileName;
      std::string tempstring;
      tempstring = deformationFieldFileNames[i];

      tempstring.erase(tempstring.end() - 4, tempstring.end() );
      resampledDeformationFieldFileName.append(tempstring);
      resampledDeformationFieldFileName.append("_sub.nii.gz");

      float* df_arr = new float[size_dfn];
      LoadIntoArray(resampledDeformationFieldFileName, df_arr);
      vnl_vector<float> df_cur(size_dfn);
      df_cur.copy_in(df_arr);

      df_mean += df_cur;
      df_mat.set_column(i, df_cur);
      delete[] df_arr;
      }

    // calculate mean vector
    df_mean /= (float)(numFiles);
    // subtract mean vector from each vector
    for( int i = 0; i < numFiles; ++i )
      {
      df_mat.set_column(i, df_mat.get_column(i) - df_mean);
      }

    // do SVD
    vnl_svd_economy<float> svd_e(df_mat);
    // build eigen-vector matrix of original deformation field matrix
    for( int i = 0; i < numEigenVector; ++i )
      {
      // U = D*V*Sinv
      df_eigenvalues.data_block()[i] = svd_e.lambdas().data_block()[i];
      df_eigenvector.set_column(i, df_mat * svd_e.V().get_column(i) );

      // df_eigenvector.scale_column(i, 1/df_eigenvalues.data_block()[i]);
      }
    // calculate coefficients and generate all intermediate templates
    // char intermediateTemplateListName[MAX_FILE_NAME_LENGTH];
    // strcpy(intermediateTemplateListName, "simulated_image_list.txt");
    // std::ofstream intermediateFileNamesFile;
    // intermediateFileNamesFile.open (intermediateTemplateListName);
    } // end if doPCA

  // for debugging
  // cerr << "Pass: PCA" << std::endl;

  // make a temporary folder to store the intermediate files
  std::string tempFolder = "temp_PCATraining";
  char        numStr[10];
  int         num = rand() % 10000;
  sprintf(numStr, "%04d", num);
  tempFolder = tempFolder + numStr;

  itksys::SystemTools::MakeDirectory(tempFolder.c_str() );

  // output folder, after selection of simulated atlases
  std::string outputFolder = "";
  const size_t      sep = allImgFileName[0].find_last_of(FILESEP);
  if( sep != std::string::npos )
    {
    outputFolder = allImgFileName[0].substr(0, sep);
    }

  if( 1 ) // if (0) // if (1)
    {
    std::vector<int> coeff(numEigenVector,0);
    const int numAllCombinations = (int)pow( (float)numSampleEachDirection, numEigenVector);

    // template image: the root of the tree
    InternalImageType::Pointer templateImage = 0;
    if( imgoperator->ReadImage(allImgFileName[root], templateImage) != 0 )
      {
      std::cerr << " Cannot read image: " << allImgFileName[root] << std::endl;
      std::cerr << "Please verify the file exists. " << std::endl;
      return -1;
      }

    std::cout << "Generate intermediate deformations from PCA results... " << std::endl;
    for( int i = 0; i < numAllCombinations; ++i )
      {
      std::cerr << i << ", ";
      for( int j = 0; j < numEigenVector; ++j )
        {
        coeff[j] = 0;
        }
      int num = i;
      int index = 0;
      while( num > 0 )
        {
        coeff[index] = (num % numSampleEachDirection);
        num = (num - coeff[index]) / numSampleEachDirection;
        index++;
        }

      // generate subsampled intermediate deformation field
      // create a new file list
      vnl_vector<float> df_intermediate_sub(size_dfn);
      df_intermediate_sub = df_mean;
      for( int j = 0; j < numEigenVector; ++j )
        {
        // df_intermediate_sub += c[coeff[j]]*df_eigenvector.get_column(j)*df_eigenvalues.data_block()[j];
        df_intermediate_sub += c[coeff[j]] * df_eigenvector.get_column(j);
        }

      std::string index_string;    basicoperator->myitoa( i, index_string, 3 );
      std::string intermediateSubDeformationFieldFileName = "inter_deform_sub_000.nii.gz";
      intermediateSubDeformationFieldFileName.erase(
        intermediateSubDeformationFieldFileName.end() - 7, intermediateSubDeformationFieldFileName.end() );
      intermediateSubDeformationFieldFileName.append(index_string);
      intermediateSubDeformationFieldFileName.append(".nii.gz");
      intermediateSubDeformationFieldFileName = tempFolder + FILESEP + intermediateSubDeformationFieldFileName;

      std::string intermediateDeformationFieldFileName = "inter_deform_000.nii.gz";
      intermediateDeformationFieldFileName.erase(
        intermediateDeformationFieldFileName.end() - 7, intermediateDeformationFieldFileName.end() );
      intermediateDeformationFieldFileName.append(index_string);
      intermediateDeformationFieldFileName.append(".nii.gz");
      intermediateDeformationFieldFileName = tempFolder + FILESEP + intermediateDeformationFieldFileName;

      SaveFromArray(intermediateSubDeformationFieldFileName, df_intermediate_sub.data_block(), size_xn, size_yn,
                    size_zn);

      // upSample the above deformation field and save
      if( doResampleDeformationField )
        {
        dfoperator->SetImx(size_x);
        dfoperator->SetImy(size_y);
        dfoperator->SetImz(size_z);

        dfoperator->UpResampleDeformationField(intermediateSubDeformationFieldFileName,
                                               intermediateDeformationFieldFileName, sampleRate);
        }

      std::string intermediateTemplateFileName = "inter_template_000.nii.gz";
      intermediateTemplateFileName.erase(intermediateTemplateFileName.end() - 7, intermediateTemplateFileName.end() );
      intermediateTemplateFileName.append(index_string);
      intermediateTemplateFileName.append(".nii.gz");
      intermediateTemplateFileName = tempFolder + FILESEP + intermediateTemplateFileName;

      // reverse deformation field
      std::string intermediateReversedDeformationFieldFileName = "inter_deform_reverse_000.nii.gz";
      intermediateReversedDeformationFieldFileName.erase(
        intermediateReversedDeformationFieldFileName.end() - 7,
        intermediateReversedDeformationFieldFileName.end() );
      intermediateReversedDeformationFieldFileName.append(index_string);
      intermediateReversedDeformationFieldFileName.append(".nii.gz");
      intermediateReversedDeformationFieldFileName = tempFolder + FILESEP
        + intermediateReversedDeformationFieldFileName;

      DeformationFieldType::Pointer deformationField = 0;
      if( dfoperator->ReadDeformationField(intermediateDeformationFieldFileName, deformationField) != 0 )
        {
        std::cerr << " Cannot read deformation field: " << intermediateDeformationFieldFileName << std::endl;
        std::cerr << "Please verify the file exists. " << std::endl;
        return -1;
        }

      DeformationFieldType::Pointer reversedDeformationField = DeformationFieldType::New();
      reversedDeformationField->SetRegions(deformationField->GetLargestPossibleRegion() );
      reversedDeformationField->SetSpacing(deformationField->GetSpacing() );
      reversedDeformationField->SetDirection(deformationField->GetDirection() );
      reversedDeformationField->SetOrigin(deformationField->GetOrigin() );
      reversedDeformationField->Allocate();
      dfoperator->InverseDeformationField3D(deformationField, reversedDeformationField);
      dfoperator->WriteDeformationField(intermediateReversedDeformationFieldFileName, reversedDeformationField);

      // apply intermediate deformation field to get intermediate template
      InternalImageType::Pointer intermediateTemplateImage = 0;
      dfoperator->ApplyDeformationField(templateImage, reversedDeformationField, intermediateTemplateImage, true);
      // write image
      imgoperator->WriteImage(intermediateTemplateFileName, intermediateTemplateImage);

      // intermediateFileNamesFile << intermediateTemplateFileName << std::endl;
      } // end for numOfAllCombinations
        // intermediateFileNamesFile.close();

    // do selection
      {
      std::cout << "Select the best simulated images... " << std::endl;
      // calculate pairwise distance between each atlas image and each simulated image
      double* * dist = new double *[m_AtlasSize];
      for( int i = 0; i < m_AtlasSize; ++i )
        {
        dist[i] = new double[numAllCombinations];
        for( int j = 0; j < numAllCombinations; ++j )
          {
          dist[i][j] = 0.0;
          }
        }
      double distmax = 0.0;
      for( int i = 0; i < m_AtlasSize; ++i )
        {
        for( int j = 0; j < numAllCombinations; ++j )
          {
          std::string index_string;     basicoperator->myitoa( j, index_string, 3 );
          std::string curInterTempFileName = "inter_template_000.nii.gz";
          curInterTempFileName.erase(curInterTempFileName.end() - 7, curInterTempFileName.end() );
          curInterTempFileName.append(index_string);
          curInterTempFileName.append(".nii.gz");
          curInterTempFileName = tempFolder + FILESEP + curInterTempFileName;

          dist[i][j] = imgoperator->calculateDistanceMSD(allImgFileName[i], curInterTempFileName);
          if( dist[i][j] > distmax )
            {
            distmax = dist[i][j];
            }
          }
        }

      // calculate the minimum distance between each simulated image and all atlas images
      double* dist_min_each_inter = new double[numAllCombinations];
      int*    index = new int[numAllCombinations];
      for( int j = 0; j < numAllCombinations; ++j )
        {
        index[j] = j;
        dist_min_each_inter[j] = distmax;
        for( int i = 0; i < m_AtlasSize; ++i )
          {
          if( dist[i][j] < dist_min_each_inter[j] )
            {
            dist_min_each_inter[j] = dist[i][j];
            }
          }
        }
      // sort
      basicoperator->bubbleSort(dist_min_each_inter, index, numAllCombinations);
      // copy files to new names after selection
      for( int i = 0; i < m_SimulateSize; ++i )
        {
        int index_inter = index[i + (numAllCombinations - m_SimulateSize)];

        std::string index_string2;      basicoperator->myitoa( index_inter, index_string2, 3 );
        std::string curInterTempFileName = "inter_deform_000.nii.gz";
        curInterTempFileName.erase(curInterTempFileName.end() - 7, curInterTempFileName.end() );
        curInterTempFileName.append(index_string2);
        curInterTempFileName.append(".nii.gz");
        curInterTempFileName = tempFolder + FILESEP + curInterTempFileName;

        std::string index_string3;    basicoperator->myitoa( i, index_string3, 3 );
        std::string outDFName = "simulated_deform_" + index_string3 + ".nii.gz";
        outDFName = outputFolder + FILESEP + outDFName;

        itksys::SystemTools::CopyFileAlways(curInterTempFileName.c_str(), outDFName.c_str() );
        }

      delete[] index;
      for( int i = 0; i < m_AtlasSize; ++i )
        {
        delete[] dist[i];
        }
      delete[] dist;
      delete[] dist_min_each_inter;
      } // do selection
        // remove irrelevant files
    for( int i = 0; i < numAllCombinations; ++i )
      {
      std::string index_string;     basicoperator->myitoa( i, index_string, 3 );

      std::string curInterTempFileName = "inter_template_" + index_string + ".nii.gz";
      std::string curInterTempDeformFileName = "inter_deform_" + index_string + ".nii.gz";
      std::string curInterTempRevFileName = "inter_deform_reverse_" + index_string + ".nii.gz";
      std::string curInterTempSubFileName = "inter_deform_sub_" + index_string + ".nii.gz";

      curInterTempFileName = tempFolder + FILESEP + curInterTempFileName;
      curInterTempDeformFileName = tempFolder + FILESEP + curInterTempDeformFileName;
      curInterTempRevFileName = tempFolder + FILESEP + curInterTempRevFileName;
      curInterTempSubFileName = tempFolder + FILESEP + curInterTempSubFileName;
      // for debugging
      basicoperator->RemoveFile(curInterTempFileName);
      basicoperator->RemoveFile(curInterTempDeformFileName);
      basicoperator->RemoveFile(curInterTempRevFileName);
      basicoperator->RemoveFile(curInterTempSubFileName);
      }
    itksys::SystemTools::RemoveADirectory(tempFolder.c_str() );
    std::cerr << "Done!" << std::endl;
    // clear newed variables
    } // end if (1)

  // do PCA training
  ////////////////////////////////////////
  delete[] c;

  return 0;
}

template <class TInputImage, class TOutputImage>
void
MABMISSimulateData<TInputImage, TOutputImage>
::SaveFromArray(std::string  deformationFieldFileName, float* df_vector, int sx, int sy, int sz)
{
  std::cerr << "deformationFieldFileName: " << deformationFieldFileName << std::endl;

  // save deformation field array into file
  // create the resampled image
  DeformationFieldType::Pointer   deformationField = DeformationFieldType::New();
  DeformationFieldType::IndexType start;
  start[0] = 0; start[1] = 0; start[2] = 0;
  DeformationFieldType::SizeType size;
  size[0] = sx; size[1] = sy; size[2] = sz;
  DeformationFieldType::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);

  deformationField->SetRegions(region);
  // sampledImage -> SetRegions (originImage->GetLargestPossibleRegion());
  // deformationField -> SetSpacing (originImage->GetSpacing());
  // deformationField -> SetDirection (originImage->GetDirection());
  // deformationField -> SetOrigin (originImage->GetOrigin());
  deformationField->SetSpacing(df_spacing);
  deformationField->SetDirection(df_direction);
  deformationField->SetOrigin(df_origin);

  deformationField->Allocate();

  DeformationFieldIteratorType itDF(deformationField, deformationField->GetLargestPossibleRegion() );
  VectorPixelType              vectorPixel;
  int                          index = 0;
  for( itDF.GoToBegin(); !itDF.IsAtEnd(); ++itDF )
    {
    vectorPixel.SetElement(0, df_vector[index]);
    vectorPixel.SetElement(1, df_vector[index + 1]);
    vectorPixel.SetElement(2, df_vector[index + 2]);
    itDF.Set(vectorPixel);
    index += 3;
    }

  dfoperator->WriteDeformationField(deformationFieldFileName, deformationField);

  return;
}

template <class TInputImage, class TOutputImage>
void
MABMISSimulateData<TInputImage, TOutputImage>
::LoadIntoArray(std::string resampledDeformationFieldFileName, float* df_vector)
{
  // read deformation field and load it into array
  DeformationFieldType::Pointer dfImage = 0;

  dfoperator->ReadDeformationField(resampledDeformationFieldFileName, dfImage);

  DeformationFieldType::SizeType im_size = dfImage->GetLargestPossibleRegion().GetSize();

  // load original image
  DeformationFieldIteratorType itOrigin(dfImage, dfImage->GetLargestPossibleRegion() );
  int index = 0;
  for( itOrigin.GoToBegin(); !itOrigin.IsAtEnd(); ++itOrigin )
    {
    const VectorPixelType & vectorPixel = itOrigin.Get();

    df_vector[index] = vectorPixel.GetElement(0);
    df_vector[index + 1] = vectorPixel.GetElement(1);
    df_vector[index + 2] = vectorPixel.GetElement(2);
    index += 3;
    }

  return;
}
} // namespace Statistics
} // namespace itk

#endif
