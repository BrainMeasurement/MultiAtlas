#ifndef __itkMABMISTreeOperation_txx
#define __itkMABMISTreeOperation_txx

#include "itkMABMISTreeOperation.h"

using namespace std;

namespace itk
{
namespace Statistics
{
template <class TInputImage, class TOutputImage>
MABMISTreeOperation<TInputImage, TOutputImage>
::MABMISTreeOperation()
{
  isDebug = true;

  dfoperator = DeformationFieldOperationType::New();
  imgoperator = ImageOperationType::New();
  basicoperator = BasicOperationType::New();
}

template <class TInputImage, class TOutputImage>
MABMISTreeOperation<TInputImage, TOutputImage>
::~MABMISTreeOperation()
{
}

template <class TInputImage, class TOutputImage>
void
MABMISTreeOperation<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}

template <class TInputImage, class TOutputImage>
vnl_vector<int>
MABMISTreeOperation<TInputImage, TOutputImage>
::BuildCombinativeTree(int root, std::vector<std::string> atlasfilenames, int tree_size, vnl_vector<int> tree)
{
  vnl_matrix<double> distanceMatrix(tree_size, tree_size);
  for( int i = 0; i < tree_size; i++ )
    {
    for( int j = 0; j < tree_size; j++ )
      {
      distanceMatrix[i][j] = 0.0;
      }
    }

  imgoperator->PairwiseDistanceAmongImages(atlasfilenames, tree_size, distanceMatrix);

  // build MST
  tree = generateMSTFromMatrixWithFixedRoot(distanceMatrix, tree_size, tree, root);

  // save distance matrix and tree
  if( isDebug )
    {
    string distFileName = "dist_all_atlases_and_simulated.txt";
    basicoperator->SaveMatrix2File(distanceMatrix, tree_size, tree_size, distFileName);
    string treeFileName = "tree_all_atlases_and_simulated.txt";
    SaveTreeWithInfo(tree, tree_size, treeFileName);
    }

  return tree;
}

template <class TInputImage, class TOutputImage>
vnl_vector<int>
MABMISTreeOperation<TInputImage, TOutputImage>
::BuildIncrementalTree(int sample_size, int tree_size, int allfilenumber,  int* sample_cur, int* atlas_cur,
                       vnl_matrix<double> cross_dist, vnl_vector<int> itree)
{
  int itree_size_cur = tree_size; // allfilenumber
  int sample_left_size = allfilenumber - itree_size_cur;

  for( int ii = 0; ii < sample_size; ii++ )
    {
    // find the min for each row
    double* test2train_min_dist = new double[sample_left_size];
    int*    test2train_min_index = new int[sample_left_size];
      {
      int*    index_train = new int[itree_size_cur];
      double* dist_temp1 = new double[itree_size_cur];
      for( int j = 0; j < sample_left_size; j++ )
        {
        for( int i = 0; i < itree_size_cur; i++ )
          {
          index_train[i] = i;
          dist_temp1[i] = cross_dist[sample_cur[j] - m_AllAtlasNumber][atlas_cur[i]];
          }
        basicoperator->bubbleSort(dist_temp1, index_train, itree_size_cur);
        test2train_min_index[j] = atlas_cur[index_train[0]];
        test2train_min_dist[j] = dist_temp1[0];
        }
      delete[] dist_temp1;
      delete[] index_train;
      }
    // find the sample index to be attached
    int*    index_test = new int[sample_left_size];
    double* dist_temp2 = new double[sample_left_size];
    for( int j = 0; j < sample_left_size; j++ )
      {
      index_test[j] = j;
      dist_temp2[j] = test2train_min_dist[j];
      }
    basicoperator->bubbleSort(dist_temp2, index_test, sample_left_size);

    // attach it to the tree
    int sample_to_attach = sample_cur[index_test[0]];
    int best_match_atlas = test2train_min_index[index_test[0]]; // parent node
    itree[sample_to_attach] = best_match_atlas;

    // include this sample to the list of current atlas
    atlas_cur[itree_size_cur] = sample_to_attach;
    // switch this sample to  the end of the list of sample_cur
    for( int i = index_test[0]; i < sample_size; i++ )
      {
      sample_cur[i] = sample_cur[i + 1];
      }
    sample_cur[sample_size - 1] = sample_to_attach;

    delete[] dist_temp2;
    delete[] index_test;

    // the tree increases and sample left decreases
    itree_size_cur++;
    sample_left_size--;
    // delete
    // for (int i = 0; i < sample_left_size; i++) delete[] dist_matrix_cur[i];
    // delete[] dist_matrix_cur;
    delete[] test2train_min_dist;
    delete[] test2train_min_index;
    }

  // for debugging
  // cerr << "itree: " << endl;
  // for (int i = 0; i < tree_size; i++)
  //  cerr <<"* " << itree[i] << endl;

  return itree;
}

template <class TInputImage, class TOutputImage>
vnl_vector<int>
MABMISTreeOperation<TInputImage, TOutputImage>
::generateMSTFromMatrix(vnl_matrix<double> curDistance, int nnode, vnl_vector<int> treeMST )
{
  // for debugging

  string distFileName;

  distFileName = "distance_matrix_generateMSTFromMatrix.txt";
  basicoperator->SaveMatrix2File(curDistance, nnode, nnode, distFileName);

  // implement Prim algorithm
  int*    tree_v1 = new int[nnode - 1];
  int*    tree_v2 = new int[nnode - 1];
  double* tree_el = new double[nnode - 1];

  double* close_lowcost = new double[nnode];
  int*    close_vec = new int[nnode];
  for( int i = 1; i < nnode; i++ )
    {
    close_vec[i] = 0;
    close_lowcost[i] = curDistance[i][0];
    }
  close_lowcost[0] = -1;
  for( int i = 0; i < nnode - 1; i++ )
    {
    // find node in graph
    int k;
    for( k = 1; k < nnode; k++ )
      {
      if( close_lowcost[k] != -1 )
        {
        break;
        }
      }
    for( int j = k + 1; j < nnode; j++ )
      {
      if( close_lowcost[j] != -1 )
        {
        if( close_lowcost[j] < close_lowcost[k] )
          {
          k = j;
          }
        }
      }

    // add to tree
    tree_el[i] = close_lowcost[k];
    tree_v1[i] = k;
    tree_v2[i] = close_vec[k];
    close_lowcost[k] = -1;
    // adjust the candidate node set
    for( int j = 1; j < nnode; j++ )
      {
      if( close_lowcost[j] > curDistance[j][k] )
        {
        close_lowcost[j] = curDistance[j][k];
        close_vec[j] = k;
        }
      }
    }

  // implement Prim algorithm
  /////////////////////////////////////
  // output tree with node and edge

  // for (int i = 0; i < nnode-1; i++)
  //	std::cout << "(" << tree_v1[i] << ", "	<< tree_v2[i] << ") = " << tree_el[i] << std::endl;

  // output tree with node and edge
  /////////////////////////////////////////////
  // rebuild the distance matrix based on connectivity on MST
  vnl_matrix<double> curDistanceTemp(nnode, nnode);
  for( int i = 0; i < nnode; i++ )
    {
    for( int j = 0; j < nnode; j++ )
      {
      if( i != j )
        {
        curDistanceTemp.put(i, j, DBL_MAX);
        }
      else
        {
        curDistanceTemp.put(i, j, 0);
        }
      }
    }
  for( int i = 0; i < nnode - 1; i++ )
    {
    curDistanceTemp[tree_v1[i]][tree_v2[i]] = tree_el[i];
    curDistanceTemp[tree_v2[i]][tree_v1[i]] = tree_el[i];
    }

  // rebuild the distance matrix based on connectivity on MST
  /////////////////////////////////////////////
  // set root to the default on in this particular scenario
  // otherwise, find the root on MST based on the shortest total distance to all other nodes

  int root = nnode - 1;

  // find root based on the tree path
  FindCenterBasedOnShortestTotalDistance(curDistanceTemp, nnode, root);

  // for debugging
  // cerr << "root: " << root <<  endl;

  // find root based on the original graph
  // FindCenterBasedOnShortestTotalDistance(curDistance, nnode, root);

  // set root to the default on in this particular scenario
  // otherwise, find the root on MST based on the shortest total distance to all other nodes
  /////////////////////////////////////
  // generate tree based on MST and root
  // int* treeMST = new int[nnode];
  bool* treeEdgeUsed = new bool[nnode - 1];
  for( int i = 0; i < nnode - 1; i++ )
    {
    treeEdgeUsed[i] = false;
    }
  treeMST[root] = -1;

  treeMST = FromEdgeSegment2Tree(treeMST, treeEdgeUsed, tree_v1, tree_v2, root, nnode);

  // for debugging
  // cerr << "treeMST: " << endl;
  // for (int i = 0; i < nnode; i++)
  //  cerr << "!! " << treeMST[i] << endl;

  // check if the generated tree is valid
  bool isAllEdgeUsed = true;
  for( int i = 0; i < nnode - 1; i++ )
    {
    if( !treeEdgeUsed[i] )
      {
      isAllEdgeUsed = false;
      }
    }
  if( isAllEdgeUsed )
    {
    // //SaveTreeWithInfo(treeMST, nnode, currentTreeFileName);
    // std::cout << "Tree is generated!@" << std::endl;
    }
  else
    {
    std::cout << "Tree cannot be generated!@" << std::endl;
    }
  // generate tree based on MST and root
  /////////////////////////////////////
  // clean up
  // delete[] treeMST;
  delete[] close_vec;
  delete[] close_lowcost;
  delete[] tree_el;
  delete[] tree_v1;
  delete[] tree_v2;

  return treeMST;
}

template <class TInputImage, class TOutputImage>
vnl_vector<int>
MABMISTreeOperation<TInputImage, TOutputImage>
::generateMSTFromMatrixWithFixedRoot(vnl_matrix<double> curDistance, int nnode, vnl_vector<int> treeMST, int root )
{
  string distFileName;

  distFileName = "distance_matrix_generateMSTFromMatrixWithFixedRoot.txt";
  basicoperator->SaveMatrix2File(curDistance, nnode, nnode, distFileName);

  // implement Prim algorithm
  int*    tree_v1 = new int[nnode - 1];
  int*    tree_v2 = new int[nnode - 1];
  double* tree_el = new double[nnode - 1];

  double* close_lowcost = new double[nnode];
  int*    close_vec = new int[nnode];
  for( int i = 1; i < nnode; i++ )
    {
    close_vec[i] = 0;
    close_lowcost[i] = curDistance[i][0];
    }
  close_lowcost[0] = -1;
  for( int i = 0; i < nnode - 1; i++ )
    {
    // find node in graph
    int k;
    for( k = 1; k < nnode; k++ )
      {
      if( close_lowcost[k] != -1 )
        {
        break;
        }
      }
    for( int j = k + 1; j < nnode; j++ )
      {
      if( close_lowcost[j] != -1 )
        {
        if( close_lowcost[j] < close_lowcost[k] )
          {
          k = j;
          }
        }
      }

    // add to tree
    tree_el[i] = close_lowcost[k];
    tree_v1[i] = k;
    tree_v2[i] = close_vec[k];
    close_lowcost[k] = -1;
    // adjust the candidate node set
    for( int j = 1; j < nnode; j++ )
      {
      if( close_lowcost[j] > curDistance[j][k] )
        {
        close_lowcost[j] = curDistance[j][k];
        close_vec[j] = k;
        }
      }
    }

  // implement Prim algorithm
  /////////////////////////////////////
  // output tree with node and edge

  // for (int i = 0; i < nnode-1; i++)
  //	std::cout << "(" << tree_v1[i] << ", "	<< tree_v2[i] << ") = " << tree_el[i] << std::endl;

  // output tree with node and edge
  /////////////////////////////////////////////
  // rebuild the distance matrix based on connectivity on MST

  double* * curDistanceTemp = new double *[nnode];
  for( int i = 0; i < nnode; i++ )
    {
    curDistanceTemp[i] = new double[nnode];
    for( int j = 0; j < nnode; j++ )
      {
      if( i != j )
        {
        curDistanceTemp[i][j] = DBL_MAX;
        }
      else
        {
        curDistanceTemp[i][j] = 0;
        }
      }
    }
  for( int i = 0; i < nnode - 1; i++ )
    {
    curDistanceTemp[tree_v1[i]][tree_v2[i]] = tree_el[i];
    curDistanceTemp[tree_v2[i]][tree_v1[i]] = tree_el[i];
    }

  // rebuild the distance matrix based on connectivity on MST
  /////////////////////////////////////////////
  // set root to the default on in this particular scenario
  // otherwise, find the root on MST based on the shortest total distance to all other nodes

  // int root = nnode-1;

  // FindCenterBasedOnShortestTotalDistance(curDistance, nnode, root);

  // set root to the default on in this particular scenario
  // otherwise, find the root on MST based on the shortest total distance to all other nodes
  /////////////////////////////////////
  // generate tree based on MST and root
  // int* treeMST = new int[nnode];
  bool* treeEdgeUsed = new bool[nnode - 1];
  for( int i = 0; i < nnode - 1; i++ )
    {
    treeEdgeUsed[i] = false;
    }
  treeMST[root] = -1;

  treeMST = FromEdgeSegment2Tree(treeMST, treeEdgeUsed, tree_v1, tree_v2, root, nnode);

  // check if the generated tree is valid
  bool isAllEdgeUsed = true;
  for( int i = 0; i < nnode - 1; i++ )
    {
    if( !treeEdgeUsed[i] )
      {
      isAllEdgeUsed = false;
      }
    }
  if( isAllEdgeUsed )
    {
    // //SaveTreeWithInfo(treeMST, nnode, currentTreeFileName);
    // std::cout << "Tree is generated!@" << std::endl;
    }
  else
    {
    std::cout << "Tree cannot be generated!@" << std::endl;
    }

  // for debugging

  SaveTreeWithInfo(treeMST, nnode, "treeMST.txt");

  // generate tree based on MST and root
  /////////////////////////////////////
  // clean up
  // delete[] treeMST;
  delete[] close_vec;
  delete[] close_lowcost;
  delete[] tree_el;
  delete[] tree_v1;
  delete[] tree_v2;
  for( int i = 0; i < nnode; i++ )
    {
    delete[] curDistanceTemp[i];
    }
  delete[] curDistanceTemp;
  return treeMST;
}

template <class TInputImage, class TOutputImage>
void
MABMISTreeOperation<TInputImage, TOutputImage>
::SaveTreeWithInfo(vnl_vector<int> tree, int treeSize, std::string filename)
{
  std::ofstream outfile;

  outfile.open(filename.c_str() );
  outfile << "Tree Size : " << treeSize << std::endl;
  int root;
  FindRoot(tree, treeSize, root);
  outfile << "Tree Root : " << root << std::endl;
  int height;
  GetTreeHeight(tree, treeSize, height); // cerr << "pass: GetTreeHeight" << endl;
  outfile << "Tree Height : " << height << std::endl;
  for( int i = 0; i < treeSize; i++ )
    {
    outfile <<  ' ';
    outfile << std::right << std::setw(3) << i;
    outfile << ' ';
    outfile << std::right << std::fixed << std::setw(8) << std::setprecision(2) << tree[i];
    outfile << std::endl;
    }
  outfile.close();

  return;
}

template <class TInputImage, class TOutputImage>
void
MABMISTreeOperation<TInputImage, TOutputImage>
::FindCenterBasedOnShortestTotalDistance(vnl_matrix<double> distanceMatrix, int matrixSize, int & center)
{
//   //for debugging
//	cerr << "FindCenterBasedOnShortestTotalDistance:distanceMatrix: " << endl;
//	for (int i = 0; i < matrixSize; i++)
//  for (int j = 0; j < matrixSize; j++)
//	    cerr << distanceMatrix[i][j] <<" " << endl;

  // input distance matrix should have very large value on the position defined by two un-connecting nodes

  ////////////////////////
  double* * distance = new double *[matrixSize];

  for( int i = 0; i < matrixSize; i++ )
    {
    distance[i] = new double[matrixSize];
    }

  double* * distance_temp = new double *[matrixSize];
  for( int i = 0; i < matrixSize; i++ )
    {
    distance_temp[i] = new double[matrixSize];
    }
  // make a copy of distanceMatrix into distance and find maximum
  for( int i = 0; i < matrixSize; i++ )
    {
    for( int j = 0; j < matrixSize; j++ )
      {
      distance[i][j] = distanceMatrix[i][j];
      }
    }
  // make distance matrix be symmetric
  for( int i = 0; i < matrixSize; i++ )
    {
    for( int j = i + 1; j < matrixSize; j++ )
      {
      if( distance[i][j] > distance[j][i] )
        {
        distance[i][j] = distance[j][i];
        }
      else
        {
        distance[j][i] = distance[i][j];
        }
      }
    }
  // make a copy of the symmetric distance matrix
  for( int i = 0; i < matrixSize; i++ )
    {
    for( int j = 0; j < matrixSize; j++ )
      {
      distance_temp[i][j] = distance[i][j];
      }
    }

  // find the minimum distance between two points on kNN map
  int ind = 0;
  while( ind < matrixSize )
    {
    for( int i = 0; i < matrixSize; i++ )
      {
      for( int j = 0; j < matrixSize; j++ )
        {
        if( distance[i][ind] + distance[ind][j] < distance[i][j] )
          {
          distance_temp[i][j] = distance[i][ind] + distance[ind][j];
          }
        }
      }
    for( int i = 0; i < matrixSize; i++ )
      {
      for( int j = 0; j < matrixSize; j++ )
        {
        distance[i][j] = distance_temp[i][j];
        }
      }
    ind = ind + 1;
    }

  // calculate the sum of each row

  int* index = new int[matrixSize];

  double* distanceSum = new double[matrixSize];
  for( int i = 0; i < matrixSize; i++ )
    {
    distanceSum[i] = 0.0;
    }
  for( int i = 0; i < matrixSize; i++ )
    {
    for( int j = 0; j < matrixSize; j++ )
      {
      distanceSum[i] += distance[i][j];
      }
    index[i] = i;
    }
  basicoperator->bubbleSort(distanceSum, index, matrixSize);

  // for debugging
  // cerr << "bubbleSort:" << endl;
  // for (int i = 0; i < matrixSize; i++) cerr << index[i] << " " << endl;

  // center = index[matrixSize-1];
  center = index[0];
  // clean up
  delete[] index;
  delete[] distanceSum;
  for( int i = 0; i < matrixSize; i++ )
    {
    delete[] distance[i];
    }
  delete[] distance;
  for( int i = 0; i < matrixSize; i++ )
    {
    delete[] distance_temp[i];
    }
  delete[] distance_temp;

  return;
}

template <class TInputImage, class TOutputImage>
vnl_vector<int>
MABMISTreeOperation<TInputImage, TOutputImage>
::FromEdgeSegment2Tree(vnl_vector<int> tree, bool* treeEdgeUsed, int* v1, int* v2, int curRoot, int matrixSize)
{
  int child = -1;

  for( int i = 0; i < matrixSize - 1; i++ )
    {
    if( (treeEdgeUsed[i] == false) && ( (v1[i] == curRoot) || (v2[i] == curRoot) ) )
      {
      child = v1[i] + v2[i] - curRoot;

      // cerr << "v1: " << v1[i] << "v2: " << v2[i] << "curRoot: " << curRoot << "child: " << child << endl;
      tree[child] = curRoot;
      treeEdgeUsed[i] = true;
      tree = FromEdgeSegment2Tree(tree, treeEdgeUsed, v1, v2, child, matrixSize);
      }
    }
  return tree;
}

// Find the root node in a tree
template <class TInputImage, class TOutputImage>
void
MABMISTreeOperation<TInputImage, TOutputImage>
::FindRoot(vnl_vector<int> tree, int treeSize, int & root)
{
  //  cerr << "tree: " << endl;
  // for debugging
  // for (int i = 0; i < treeSize; i++)
  //  cerr << "**" << tree[i] << " " << endl;
  for( int i = 0; i < treeSize; i++ )
    {
    if( tree[i] == -1 )
      {
      root = i;
      return;
      }
    }
  root = -1;
  return;
}

// calculate the height of each node in a tree
template <class TInputImage, class TOutputImage>
void
MABMISTreeOperation<TInputImage, TOutputImage>
::CalculateNodeSegmentsToRoot(vnl_vector<int> tree, int treeSize, vnl_vector<int>& height)
{
  for( int i = 0; i < treeSize; i++ )
    {
    height[i] = 0;
    }
  for( int i = 0; i < treeSize; i++ )
    {
    bool isAtRoot = false;
    int  curNode = i;
    while( !isAtRoot )
      {
      if( tree[curNode] >= 0 )
        {
        height[i]++;
        curNode = tree[curNode];
        }
      else
        {
        isAtRoot = true;
        }
      }
    }
  return;
}

// get the height of the tree
template <class TInputImage, class TOutputImage>
void
MABMISTreeOperation<TInputImage, TOutputImage>
::GetTreeHeight(vnl_vector<int> tree, int treeSize, int & treeHeight)
{
  // cerr << "tree: " << endl;
  // for debugging
  // for (int i = 0; i < treeSize; i++)
  //  cerr << tree[i] << " " << endl;

  vnl_vector<int> height(treeSize);
  CalculateNodeSegmentsToRoot(tree, treeSize, height); // cerr << "pass: CalculateNodeSegmentsToRoot" << endl;
  treeHeight = 0;
  for( int i = 0; i < treeSize; i++ )
    {
    if( height[i] > treeHeight )
      {
      treeHeight = height[i];
      }
    }
  return;
}

template <class TInputImage, class TOutputImage>
bool
MABMISTreeOperation<TInputImage, TOutputImage>
::ValidateTree(vnl_vector<int> tree, int tree_size)
{
  // for debugging
  // cerr << "ValidateTree: " << endl;
  // for (int i = 0; i < tree_size; i++)
  //  cerr << "t " << tree[i] << endl;

  // to check if there is one and only one root
  int root;
  int root_num = 0;

  for( int i = 0; i < tree_size; i++ )
    {
    if( tree[i] == -1 )
      {
      root = i;
      root_num++;
      }
    if( (tree[i] < -1) || (tree[i] >= tree_size) )
      {
      return false;
      }
    }
  if( root_num != 1 )
    {
    return false;
    }
  // for debugging
  // cerr << "root_num: " << root_num ;
  // to check if all node can trace back to the root
  for( int i = 0; i < tree_size; i++ )
    {
    //
    bool isRooted = false;
    int  max_step = tree_size;
    int  j = 0;
    int  cur_node = i;
    while( j < max_step )
      {
      if( cur_node == root )
        {
        isRooted = true;
        break;
        }
      else
        {
        cur_node = tree[cur_node];
        j++;
        }
      }

    if( !isRooted )
      {
      return false;
      }

    // for debugging
    // cerr << "isRooted: " << isRooted ;
    }

  return true;
}
} // namespace Statistics
} // namespace itk

#endif
