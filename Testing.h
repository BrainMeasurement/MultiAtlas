#ifndef __Testing_h
#define __Testing_h

#include <vector>
#include "itkMABMISAtlasXMLFile.h"

std::string ReplacePathSepForOS(const std::string & input);

int Testing(itk::MABMISImageData* imageData, itk::MABMISAtlas* atlasTree,
            std::vector<int> iterations, double sigma);

#endif //__Testing_h
