#ifndef __itkMABMISTreeOperation_h
#define __itkMABMISTreeOperation_h

#include <itkImage.h>
#include <itkImageToImageFilter.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include <itkImage.h>

#include "itkMABMISDeformationFieldFilter.h"
#include "itkMABMISImageOperationFilter.h"
#include "itkMABMISBasicOperationFilter.h"

#define ImageDimension 3

namespace itk
{
namespace Statistics
{
template <class TInputImage, class TOutputImage>
class MABMISTreeOperation : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:

  typedef MABMISTreeOperation                           Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  typedef TInputImage                 ImageType;
  typedef typename ImageType::Pointer ImagePointerType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MABMISTreeOperation, ImageToImageFilter);

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

  DeformationFieldType::SpacingType   df_spacing;
  DeformationFieldType::DirectionType df_direction;
  DeformationFieldType::PointType     df_origin;

  typedef itk::Statistics::MABMISDeformationFieldFilter<ImageType, ImageType> DeformationFieldOperationType;
  typedef itk::Statistics::MABMISImageOperationFilter<ImageType, ImageType>   ImageOperationType;
  typedef itk::Statistics::MABMISBasicOperationFilter<ImageType, ImageType>   BasicOperationType;

  typename DeformationFieldOperationType::Pointer dfoperator;
  typename ImageOperationType::Pointer imgoperator;
  typename BasicOperationType::Pointer basicoperator;

  vnl_vector<int> BuildCombinativeTree(int root, std::vector<std::string> atlasfilenames, int tree_size,
                                       vnl_vector<int> tree);

  vnl_vector<int>  BuildIncrementalTree(int sample_size, int tree_size, int allfilenumber,  int* sample_cur,
                                        int* atlas_cur, vnl_matrix<double> cross_dist, vnl_vector<int> itree);

  vnl_vector<int> generateMSTFromMatrix(vnl_matrix<double> curDistance, int nnode, vnl_vector<int> treeMST );

  vnl_vector<int> generateMSTFromMatrixWithFixedRoot(vnl_matrix<double> curDistance, int nnode, vnl_vector<int> treeMST,
                                                     int root );

  void SaveTreeWithInfo(vnl_vector<int> tree, int treeSize, std::string filename);

  void FindCenterBasedOnShortestTotalDistance(vnl_matrix<double> distanceMatrix, int matrixSize, int & center);

  vnl_vector<int> FromEdgeSegment2Tree(vnl_vector<int> tree, bool* treeEdgeUsed, int* v1, int* v2, int curRoot,
                                       int matrixSize);

  void FindRoot(vnl_vector<int> tree, int treeSize, int & root);

  void GetTreeHeight(vnl_vector<int> tree, int treeSize, int & treeHeight);

  void CalculateNodeSegmentsToRoot(vnl_vector<int> tree, int treeSize, vnl_vector<int>& height);

  bool ValidateTree(vnl_vector<int> tree, int tree_size);

  itkSetMacro(Root, int);
  itkSetMacro(Imx, int);
  itkSetMacro(Imy, int);
  itkSetMacro(Imz, int);
  itkSetMacro(AtlasSize, int);
  itkSetMacro(SimulateSize, int);
  itkSetMacro(AllAtlasNumber, int);
private:
  MABMISTreeOperation(const Self &); // purposely not implemented
  void operator=(const Self &);      // purposely not implemented

  int m_Root;

  int m_Imx;
  int m_Imy;
  int m_Imz;

  int m_AtlasSize;
  int m_SimulateSize;
  int m_AllAtlasNumber;

  bool isDebug;
protected:
  MABMISTreeOperation();
  ~MABMISTreeOperation();
  void PrintSelf(std::ostream& os, Indent indent) const;
};
} // namespace itk
} // namespace Statistics

#include "itkMABMISTreeOperation.hxx"
#endif
