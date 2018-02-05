#include "commonMABMIS.h"

#include <algorithm>

std::string ReplacePathSepForOS( const std::string & input )
{
  std::string output = input;
#ifdef _WIN32
  std::replace(output.begin(), output.end(), '/', FILESEP);
#else
  std::replace(output.begin(), output.end(), '\\', FILESEP);
#endif
  return output;
}

std::string
GetExtension(const std::string & filename)
{
  std::string fileExt( itksys::SystemTools::GetFilenameLastExtension(filename) );
  //If the last extension is .gz, then need to pull off 2 extensions.
  //.gz is the only valid compression extension.
  if ( fileExt == std::string(".gz") )
    {
    fileExt = itksys::SystemTools::GetFilenameLastExtension( 
              itksys::SystemTools::GetFilenameWithoutLastExtension(filename) );
    fileExt += ".gz";
    }
  return ( fileExt );
}

std::string
GetRootName(const std::string & filename)
{
  const std::string fileExt = GetExtension(filename);

  // Create a base filename
  // i.e Image.hdr --> Image
  if ( fileExt.length() > 0                    //Ensure that an extension was found
       && filename.length() > fileExt.length() //Ensure that the filename does
                                               // not contain only the extension
       )
    {
    const std::string::size_type it = filename.find_last_of(fileExt);
    const std::string            baseName( filename, 0, it - ( fileExt.length() - 1 ) );
    return ( baseName );
    }
  //Default to return same as input when the extension is nothing (Analyze)
  return ( filename );
}

// demons registration parameters
// int iterInResolutions[4][3]={{5,3,2},{10,5,5},{15,10,5},{20,15,10}};
// int itereach = 2;
// int itereach0 = 0;
// int itereach1 = 1;
// int itereach2 = 2;
// int itereach3 = 3;
// double sigmaDef = 1.5;
// double sigmaDef10 = 1.0;
// double sigmaDef15 = 1.5;
// double sigmaDef20 = 2.0;
// double sigmaDef25 = 2.5;
// double sigmaDef30 = 3.0;
// double sigmaDef35 = 3.5;

DataSimulatorType::Pointer datasimulator = DataSimulatorType::New();
ImageOperationFilterType::Pointer imgoperator = ImageOperationFilterType::New();
DeformationFieldOperationFilterType::Pointer dfoperator = DeformationFieldOperationFilterType::New();
ImageRegistrationFilterType::Pointer regoperator = ImageRegistrationFilterType::New();
TreeOperationType::Pointer treeoperator = TreeOperationType::New();
BasicOperationFilterType::Pointer basicoperator = BasicOperationFilterType::New();
