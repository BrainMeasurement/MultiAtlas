#ifndef __Testing_h
#define __Testing_h

#include "commonMABMIS.h"

int Testing(itk::MABMISImageData* imageData, itk::MABMISAtlas* atlasTree,
            std::vector<int> iterations, double sigma);

#endif //__Testing_h
