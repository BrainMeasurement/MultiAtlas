/*=========================================================================
  Copyright (c) IDEA LAB, UNC-Chapel Hill, 2013.

     MABMIS (Multi-Atlas-Based Multi-Image Segmentation)

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "IGR3D_MABMIS_TrainingCLP.h"
#include "Training.h"

int main( int argc, char *argv[] )
{
  PARSE_ARGS;

  // step 1: read the training atlases
  itk::MABMISImageDataXMLFileReader::Pointer trainingXMLReader = itk::MABMISImageDataXMLFileReader::New();
  TrainingDataXML = ReplacePathSepForOS(TrainingDataXML); 
  TrainingOutputFile = ReplacePathSepForOS(TrainingOutputFile); 
  trainingXMLReader->SetFilename(TrainingDataXML);
  try
    {
    trainingXMLReader->GenerateOutputInformation();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while reading file" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
  itk::MABMISImageData * trainingData = trainingXMLReader->GetOutputObject();

  // if the data path is empty, use the path of the xml file instead
  if( trainingData->m_DataDirectory.size() <= 1 )
    {
    const size_t sep = TrainingDataXML.find_last_of(FILESEP);
	//std::cout <<"sep=" << sep << std::endl;
    if( sep != std::string::npos )
      {
      trainingData->m_DataDirectory = TrainingDataXML.substr(0, sep);
      }
    else
      {
      trainingData->m_DataDirectory.resize(0);
      }
    }
  if( !trainingData->m_DataDirectory.empty() )
    {
    if( !(trainingData->m_DataDirectory[trainingData->m_DataDirectory.size() - 1] == FILESEP) )
      {
      trainingData->m_DataDirectory = trainingData->m_DataDirectory + FILESEP;
      }
    }

  int retVal = Training(trainingData, TrainingOutputFile, iterations, SmoothingKernelSize);
  delete trainingData;
  return retVal;
}
