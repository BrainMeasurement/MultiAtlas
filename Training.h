#ifndef __Training_h
#define __Training_h

#include <vector>
#include "itkMABMISAtlasXMLFile.h"

std::string ReplacePathSepForOS(const std::string & input);

int Training(itk::MABMISImageData* trainingData, std::string outputFile,
             std::vector<int> iterations, double sigma);

#endif //__Training_h
