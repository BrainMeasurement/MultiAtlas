#ifndef __Training_h
#define __Training_h

#include "commonMABMIS.h"

int Training(itk::MABMISImageData* trainingData, std::string outputFile,
             std::vector<int> iterations, double sigma);

#endif //__Training_h
