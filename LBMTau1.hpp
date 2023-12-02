#pragma once

#include <vector>
#include <iostream>
#include <cstring>
#include <bitset>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cassert>
#ifdef USE_LARGE_PAGES
  #include <sys/mman.h>
#endif

template <class DividentType>
unsigned ceilDivision (DividentType dividend, unsigned divisor)
{
  return (dividend + divisor - 1) / divisor ;
}

struct dim3
{
  int x = 1 ;
  int y = 1 ;
  int z = 1 ;
} ;

using IndexType = size_t ;

inline IndexType rowMajorIndex 
(
  unsigned x, unsigned y, unsigned z,
  dim3 latticeDimension
)
{
  return x + y * latticeDimension.x + 
    static_cast <IndexType>(z) * latticeDimension.x * latticeDimension.y ;
}

struct
NodeType
{
  using DType = uint8_t ;

  enum class Type : DType
  {
    FLUID = 0, 
    SOLID = 1
  } ;

  static constexpr DType typeMask_  = 0xF ;
  static constexpr int   typeWidth_ = 4 ; // In bits, must match typeMask_.

  DType packedType_ = static_cast<DType> (Type::FLUID) ;

  bool isFluid() const 
  {
    return Type::FLUID == static_cast<Type> (packedType_ & typeMask_) ; 
  }
  bool isSolid() const 
  {
    return Type::SOLID == static_cast<Type> (packedType_ & typeMask_) ; 
  }
  bool hasOnlyFluidNeighbours() const
  {
    return 0 == (packedType_ & (~ typeMask_)) ;
  }
  bool isFluidOnlyFragment() const
  {
    return isFluid() && hasOnlyFluidNeighbours() ;
  }
  DType packedValue() const { return packedType_ ; }
  DType type() const { return packedType_ & typeMask_ ; }

  //TODO: WARNING: setFluid() overrides information about neighbourhood
  //      but setSolid() does NOT override neighbour information.
  //      In current version, first ALL nodes must be set to fluid and, 
  //      after that, some nodes can be set to solid.
  void setFluid() 
  { 
    packedType_ = static_cast<DType> (Type::FLUID) ; 
  }
  void setSolid() 
  { 
    packedType_ |= static_cast<DType> (Type::SOLID) ; 
  }
  void setNonFluidNeighbour()
  {
    packedType_ |= (1 << typeWidth_) ;
  }
} ;

template <unsigned D, unsigned Q, class DataType> struct LatticeArrangement ;

template <class DataType>
struct LatticeArrangement <3,19,DataType>
{
  static constexpr int D =  3 ;
  static constexpr int Q = 19 ;

  static constexpr int ex (int idx) 
  { 
  constexpr int ex_ [Q] = 
// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
  {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0};
  
    return ex_ [idx] ; 
  }
  static constexpr int ey (int idx)
  {
  constexpr int ey_ [Q] = 
// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
  {0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1};
  
    return ey_ [idx] ;
  }
  static constexpr int ez (int idx)
  {
  constexpr int ez_ [Q] = 
// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
  {0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1};
  
    return ez_ [idx] ;
  }
  static constexpr int inv (int idx)
  {
  constexpr int inv_ [Q] = 
// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
  {0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17};

    return inv_ [idx] ;
  }

  static constexpr DataType a = static_cast <DataType> (1) /  3 ;
  static constexpr DataType b = static_cast <DataType> (1) / 18 ;
  static constexpr DataType c = static_cast <DataType> (1) / 36 ;

  static constexpr DataType w (int idx)
  {
  constexpr DataType w_ [Q] = 
// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
  {a, b, b, b, b, b, b, c, c, c, c, c, c, c, c, c, c, c, c};  

    return w_ [idx] ;
  }

// Original order from 
//  The Lattice Boltzmann Method Principles and Practice
//  by T. Krüger et.al.
/** D3Q19 **/
  //                    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
  //const int ex [q] = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0};
  //const int ey [q] = {0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1};
  //const int ez [q] = {0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1};
  //const int inv[q] = {0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17};
  //const float w[q] = {a, b, b, b, b, b, b, c, c, c, c, c, c, c, c, c, c, c, c};  
} ;


template <class DataType>
struct LatticeArrangement <3,27,DataType>
{
  static constexpr int D =  3 ;
  static constexpr int Q = 27 ;

  static constexpr int ex (int idx) 
  { 
  constexpr int ex_ [Q] = 
// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
  {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1};
  
    return ex_ [idx] ; 
  }
  static constexpr int ey (int idx)
  {
  constexpr int ey_ [Q] = 
// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
  {0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, 1,-1};
  
    return ey_ [idx] ;
  }
  static constexpr int ez (int idx)
  {
  constexpr int ez_ [Q] = 
// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
  {0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1};
  
    return ez_ [idx] ;
  }
  static constexpr int inv (int idx)
  {
  constexpr int inv_ [Q] = 
// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
  {0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17,20,19,22,21,24,23,26,25};

    return inv_ [idx] ;
  }

  static constexpr DataType a = static_cast <DataType> (8) /  27 ;
  static constexpr DataType b = static_cast <DataType> (2) /  27 ;
  static constexpr DataType c = static_cast <DataType> (1) /  54 ;
  static constexpr DataType d = static_cast <DataType> (1) / 216 ;

  static constexpr DataType w (int idx)
  {
  constexpr DataType w_ [Q] = 
// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
  {a, b, b, b, b, b, b, c, c, c, c, c, c, c, c, c, c, c, c, d, d, d, d, d, d, d, d};  

    return w_ [idx] ;
  }

// Original order from 
//  The Lattice Boltzmann Method Principles and Practice
//  by T. Krüger et.al.
/** D3Q27 **/
  //const float a=8.0/27.0,b=2.0/27.0,c=1.0/54.0,d=1.0/216.0;
  //                    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
  //const int ex [q] = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1};
  //const int ey [q] = {0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, 1,-1};
  //const int ez [q] = {0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1};
  //const int inv[q] = {0, 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17,20,19,22,21,24,23,26,25};
  //const float w[q] = {a, b, b, b, b, b, b, c, c, c, c, c, c, c, c, c, c, c, c, d, d, d, d, d, d, d, d};  
} ;



// Node data for nodes from a single chunk is packed into consequtive four vectors.
// The first vector contains vx, the next vectors vy, vz, and densities for all
// nodes from the chunk.
template <class L, class T> 
struct LBMTau1Vector3D
{
  using DataType = T ;
  using LatticeArrangement = L ;

  static constexpr unsigned VSize = 4 ; //FIXME: CAN NOT be easily changed now.
  static constexpr unsigned TILE_SIZE = 4 ;

  using V  __attribute__ ((aligned (64))) __attribute__ ((vector_size (VSize * sizeof (T))))  = T ;

  static constexpr unsigned N_NODES_PER_BATCH = VSize ;

   LBMTau1Vector3D (unsigned width, unsigned height, unsigned depth) ;

  ~LBMTau1Vector3D () 
  {
    free (dataPtr_) ;
    dataPtr_ = NULL ;
    nBytesAllocated_ = 0 ;
  }

  void iter() 
  {
    iter (direction_) ;
    direction_ = (direction_ + 1) % 2 ;
  }

  void iter (int direction) 
  {
#ifdef KERNEL_REFERENCE
    iterReference (direction) ;
#elif defined KERNEL_VECTORIZED_LOOP
    iterVectorizedLoop (direction) ;
#elif defined KERNEL_OPT_MANUAL
         if constexpr (3 == L::D  &&  19 == L::Q)
    {
      iterD3Q19 (direction) ;
    }
    else if constexpr (3 == L::D  &&  27 == L::Q)
    {
      iterD3Q27 (direction) ;
    }
#elif defined KERNEL_COPY_ONLY
    iterCopyOnly (direction) ;
#endif
  }

  void iterD3Q19 (int direction) ;
  void iterD3Q27 (int direction) ;
  void iterReference (int direction) ;
  void iterVectorizedLoop (int direction) ;
  void iterCopyOnly (int direction) ;

  void setF (DataType fx, DataType fy, DataType fz)
  {
    fx_ = fx ;
    fy_ = fy ; 
    fz_ = fz ;
  }

  size_t nNodes() const { return nNodes_ ; }
  size_t nBytesAllocated() const { return nBytesAllocated_ ; }

  dim3 latticeDimension() const
  {
    dim3 lDim ;
    lDim.x = w() ;
    lDim.y = h() ;
    lDim.z = d() ;
    return lDim ;
  }

  IndexType index
  (
    unsigned x, unsigned y, unsigned z, 
    dim3 latticeDimension,
    unsigned setIndex = 0 // For multiple copies of data, usually 0 or 1.
  )
  {
    return rowMajorIndex (x, y, z, latticeDimension) +
           setIndex * secondCopyOffset() ;
  }
  IndexType index
  (
    unsigned x, unsigned y, unsigned z, 
    unsigned setIndex = 0
  )
  {
    return index (x,y,z, latticeDimension(), setIndex) ;
  }
  // Since node data is unapcked into separate vectors, returns index of 
  // the first vector from chunk.
  IndexType indexOfNodeBatch 
  (
    unsigned x, unsigned y, unsigned z, 
    dim3 latticeDimension,
    unsigned setIndex = 0 // For multiple copies of data, usually 0 or 1.
  )
  {
    unsigned xBatch = x - (x % N_NODES_PER_BATCH) ;
    return rowMajorIndex (xBatch, y, z, latticeDimension) +
           setIndex * secondCopyOffset() ;
  }
  IndexType indexOfNodeBatch
  (
    unsigned x, unsigned y, unsigned z, 
    unsigned setIndex = 0
  )
  {
    return indexOfNodeBatch (x,y,z, latticeDimension(), setIndex) ;
  }

  static constexpr size_t alignmentSecondCopy = 1 * 1024 ;
  IndexType secondCopyOffset() const
  {
    // Must be aligned too.
    return (nNodes() + alignmentSecondCopy) & (0ll - alignmentSecondCopy) ;
  }

  unsigned w() const { return w_ ; }
  unsigned h() const { return h_ ; }
  unsigned d() const { return d_ ; }

  V * dataPtr() { return dataPtr_ ; }

  V * dataPtr_ ;
  size_t nBytesAllocated_ ;
  size_t componentOffset_ ;
  int direction_ ;

  // All fields of a node are packed into a single vector.
  // These methods REQUIRE that N_NODES_PER_BATCH = VSize.
  T get_vx (unsigned x, unsigned y, unsigned z, unsigned setIndex = 0) 
  { 
    unsigned xInBatch = x % N_NODES_PER_BATCH ;
    return dataPtr() [indexOfNodeBatch (x, y, z, latticeDimension(), setIndex) + 0] [xInBatch] ; 
  }
  T get_vy (unsigned x, unsigned y, unsigned z, unsigned setIndex = 0) 
  { 
    unsigned xInBatch = x % N_NODES_PER_BATCH ;
    return dataPtr() [indexOfNodeBatch (x, y, z, latticeDimension(), setIndex) + 1] [xInBatch] ; 
  }
  T get_vz (unsigned x, unsigned y, unsigned z, unsigned setIndex = 0) 
  { 
    unsigned xInBatch = x % N_NODES_PER_BATCH ;
    return dataPtr() [indexOfNodeBatch (x, y, z, latticeDimension(), setIndex) + 2] [xInBatch] ; 
  }
  T get_r  (unsigned x, unsigned y, unsigned z, unsigned setIndex = 0) 
  { 
    unsigned xInBatch = x % N_NODES_PER_BATCH ;
    return dataPtr() [indexOfNodeBatch (x, y, z, latticeDimension(), setIndex) + 3] [xInBatch] ; 
  }

  LBMTau1Vector3D & set_vx (T val, unsigned x, unsigned y, unsigned z, unsigned setIndex = 0) 
  { 
    unsigned xInBatch = x % N_NODES_PER_BATCH ;
    dataPtr() [indexOfNodeBatch (x, y, z, latticeDimension(), setIndex) + 0] [xInBatch] = val ; 
    return * this ;
  }
  LBMTau1Vector3D & set_vy (T val, unsigned x, unsigned y, unsigned z, unsigned setIndex = 0) 
  { 
    unsigned xInBatch = x % N_NODES_PER_BATCH ;
    dataPtr() [indexOfNodeBatch (x, y, z, latticeDimension(), setIndex) + 1] [xInBatch] = val ; 
    return * this ;
  }
  LBMTau1Vector3D & set_vz (T val, unsigned x, unsigned y, unsigned z, unsigned setIndex = 0) 
  { 
    unsigned xInBatch = x % N_NODES_PER_BATCH ;
    dataPtr() [indexOfNodeBatch (x, y, z, latticeDimension(), setIndex) + 2] [xInBatch] = val ; 
    return * this ;
  }
  LBMTau1Vector3D & set_r  (T val, unsigned x, unsigned y, unsigned z, unsigned setIndex = 0) 
  { 
    unsigned xInBatch = x % N_NODES_PER_BATCH ;
    dataPtr() [indexOfNodeBatch (x, y, z, latticeDimension(), setIndex) + 3] [xInBatch] = val ; 
    return * this ;
  }

  void setSolidNeighbours (unsigned x, unsigned y, unsigned z)
  {
    for (unsigned k=0; k < L::Q; k++)
    { 
      unsigned xp = (x + L::ex(k) + w()) % w() ;
      unsigned yp = (y + L::ey(k) + h()) % h() ;
      unsigned zp = (z + L::ez(k) + d()) % d() ;

      auto idxNghbr = index (xp,yp,zp, latticeDimension()) ;
      nodeTypes_ [idxNghbr].setNonFluidNeighbour() ;

      setSolidNodeMask (xp,yp,zp) ;
    }   
  }

  void setSolid (unsigned x, unsigned y, unsigned z)
  {
    const IndexType idx = index (x, y, z, latticeDimension()) ;
    nodeTypes_ [idx].setSolid() ;
    setSolidNodeMask (x,y,z) ;

    setSolidNeighbours (x,y,z) ;

    set_vx (NAN, x,y,z, 0) ;
    set_vy (NAN, x,y,z, 0) ;
    set_vz (NAN, x,y,z, 0) ;
    set_r  (NAN, x,y,z, 0) ;

    set_vx (NAN, x,y,z, 1) ;
    set_vy (NAN, x,y,z, 1) ;
    set_vz (NAN, x,y,z, 1) ;
    set_r  (NAN, x,y,z, 1) ;
  }

  void resetAllNodesToFluid()
  {
    #pragma omp parallel for
    for (IndexType i=0 ; i < nodeTypes_.size() ; i++)
    {
      nodeTypes_ [i].setFluid() ;
    }
    
    const IndexType W = w() ;
    const IndexType H = h() ;
    const IndexType D = d() ;
    constexpr unsigned TS = TILE_SIZE ;
    // Keep the same thread mapping as for LBM iteration to force 
    // page allocation on proper NUMA nodes.
    #pragma omp parallel for 
    for (IndexType zp=0 ; zp < D ; zp += TS) 
      for (IndexType yp=0 ; yp < H ; yp += TS)
        for (IndexType zo=0 ; zo < TS ; zo++) 
        for (IndexType yo=0 ; yo < TS ; yo++) 
        {
          IndexType z = zp + zo ;
          IndexType y = yp + yo ;

          for (IndexType x=0 ; x < W ; x++)
          {
            set_vx (  0, x,y,z, 0) ;
            set_vy (  0, x,y,z, 0) ;
            set_vz (  0, x,y,z, 0) ;
            set_r  (  1, x,y,z, 0) ;

            set_vx (NAN, x,y,z, 1) ;
            set_vy (NAN, x,y,z, 1) ;
            set_vz (NAN, x,y,z, 1) ;
            set_r  (NAN, x,y,z, 1) ;
          }
        }

    resetAllNodeMasksToFluid() ;
  }

  std::vector <NodeType> nodeTypes_ ;
  unsigned w_ ;
  unsigned h_ ;
  unsigned d_ ;
  size_t nNodes_ ;  
  
  NodeType & t (unsigned x, unsigned y, unsigned z) 
    { return nodeTypes_ [index (x, y, z, latticeDimension())] ; }
  NodeType & t  (IndexType index) 
    { return nodeTypes_ [index] ; }

  // NodeMask stores information if neighbouring N_NODES_PER_BATCH are fluid AND
  // have only fluid neighbours. If yes then neighbouring N_NODES_PER_BATCH nodes
  // can be processed with efficient vector code. 
  // Two bits of NodeMask corresponds to neighbouring N_NODES_PER_BATCH nodes.
  // Bits equal to 00 denote fluid-only batch that can be processed with vector 
  // code, bits equal to 01 force operations on separate nodes according
  // to their types and neighbourhood, and bits equal to 10 mark batch containing
  // solid nodes only that can be skipped.
  static constexpr int FLUID_BATCH = 0 ; // DO NOT CHANGE, hidden dependencies!!!
  static constexpr int MIXED_BATCH = 1 ;
  static constexpr int SOLID_BATCH = 2 ;
  static constexpr unsigned N_BITS_PER_MASK = 64 ;
  static constexpr unsigned N_BITS_PER_BATCH =  2 ;
  static constexpr unsigned N_BATCHES_PER_MASK = N_BITS_PER_MASK / N_BITS_PER_BATCH ;
  static constexpr unsigned N_NODES_PER_MASK = N_NODES_PER_BATCH * N_BATCHES_PER_MASK ;
  using NodeMask = std::bitset <N_BITS_PER_MASK> ;
  std::vector <NodeMask> nodeMask_ ;
  unsigned nMasksPerRow_ ;

  void allocateNodeMask()
  {
    std::cout << "N_BATCHES_PER_MASK = " << N_BATCHES_PER_MASK << "\n" ;
    std::cout << "N_NODES_PER_MASK   = " << N_NODES_PER_MASK   << "\n" ;

    nMasksPerRow_ = ceilDivision (w(), N_NODES_PER_MASK) ;

    const IndexType nRows = static_cast<IndexType> (h()) * d() ;
    const IndexType nMasks = nRows * nMasksPerRow_ ;

    std::cout << "nMasksPerRow_ = " << nMasksPerRow_
              << " nRows = " << nRows
              << " nMasks = " << nMasks
              << "\n" ;

    nodeMask_.resize (nMasks) ;
  }
  
  //TODO: Rename to setMixedNodeMask?
  void setSolidNodeMask (unsigned x, unsigned y, unsigned z)
  {
    markBatchAsMixed (x,y,z) ;

    unsigned xb0 = x / N_NODES_PER_BATCH ;

    bool areAllNodesSolid = true ;
    for (unsigned xo=0 ; xo < N_NODES_PER_BATCH ; xo++)
    {
      unsigned x = xb0 * N_NODES_PER_BATCH + xo ;

      if (not t (x,y,z).isSolid()) areAllNodesSolid = false ;
    }

    if (areAllNodesSolid) markBatchAsSolid (x,y,z) ;
  }

  // Returns index of the first word with bitmask for the given row.
  IndexType indexNodeMaskRow (IndexType nodeY, IndexType nodeZ)
  {
    return (nodeZ * h() + nodeY) * nMasksPerRow_ ;
  }
  IndexType indexNodeMask (IndexType nodeX, IndexType nodeY, IndexType nodeZ)
  {
    return indexNodeMaskRow (nodeY, nodeZ) + (nodeX / N_NODES_PER_MASK) ;
  }
  unsigned indexNodeMaskBit (IndexType nodeX)
  {
    return ((nodeX % N_NODES_PER_MASK) / N_NODES_PER_BATCH) * N_BITS_PER_BATCH ;
  }

  unsigned getBatchType (IndexType nodeX, IndexType nodeY, IndexType nodeZ)
  {
    unsigned type = 0 ;

    const IndexType iMask = indexNodeMask (nodeX, nodeY, nodeZ) ;
    const unsigned  iBit  = indexNodeMaskBit (nodeX) ;

    type = nodeMask_ [iMask][iBit] + nodeMask_ [iMask][iBit +1] * 2 ;
    return type ;
  }
  bool isBatchMixed (unsigned batchType)
  {
    return (1 == batchType) ;
  }
  bool isBatchSolid (unsigned batchType)
  {
    return (2 == batchType) ;
  }
  bool isBatchFluid (unsigned batchType)
  {
    return (0 == batchType) ;
  }

  bool isBatchMixed (IndexType nodeX, IndexType nodeY, IndexType nodeZ)
  {
    const IndexType iMask = indexNodeMask (nodeX, nodeY, nodeZ) ;
    const unsigned  iBit  = indexNodeMaskBit (nodeX) ;
    return nodeMask_ [iMask].test (iBit) ;
  }
  bool isBatchSolid (IndexType nodeX, IndexType nodeY, IndexType nodeZ)
  {
    const IndexType iMask = indexNodeMask (nodeX, nodeY, nodeZ) ;
    const unsigned  iBit  = indexNodeMaskBit (nodeX) ;
    return nodeMask_ [iMask].test (iBit + 1) ;
  }
  bool isBatchFluid (IndexType nodeX, IndexType nodeY, IndexType nodeZ)
  {
    return not isBatchMixed (nodeX,nodeY,nodeZ) 
        && not isBatchSolid (nodeX,nodeY,nodeZ) ;
  }
  void markBatchAsMixed (IndexType nodeX, IndexType nodeY, IndexType nodeZ)
  {
    const IndexType iMask = indexNodeMask (nodeX, nodeY, nodeZ) ;
    const unsigned  iBit  = indexNodeMaskBit (nodeX) ;
  
    nodeMask_ [iMask].  set (iBit) ;
    nodeMask_ [iMask].reset (iBit +1) ;
  }
  void markBatchAsSolid (IndexType nodeX, IndexType nodeY, IndexType nodeZ)
  {
    const IndexType iMask = indexNodeMask (nodeX, nodeY, nodeZ) ;
    const unsigned  iBit  = indexNodeMaskBit (nodeX) ;
  
    nodeMask_ [iMask].reset (iBit) ;
    nodeMask_ [iMask].  set (iBit +1) ;
  }

  void resetAllNodeMasksToFluid()
  {
    NodeMask solidMask ;
    solidMask.set() ; //FIXME: Sets BOTH mixed and solid types.
    std::fill (nodeMask_.begin(), nodeMask_.end(), solidMask) ;

    // Update nodeMask_ for FLUID nodes ONLY.
    for (IndexType z=0 ; z < d() ; z++)
      for (IndexType y=0 ; y < h() ; y++)
      {
        for (IndexType x=0 ; x < w() ; x++)
        {
          const IndexType iMask = indexNodeMask (x,y,z) ;
          const unsigned  iBit  = indexNodeMaskBit (x) ;
          
          nodeMask_ [iMask].reset (iBit) ;
          nodeMask_ [iMask].reset (iBit + 1) ;
        }
      }
  }

  void checkNodeMasks()
  {
    using std::cout ;
    cout << "Checking node masks...\n" ;

    const auto D = d() ;
    const auto H = h() ;
    const auto W = w() ;

    for (IndexType z=0 ; z < D ; z++)
      for (IndexType y=0 ; y < H ; y++)
        for (IndexType xb=0 ; xb < W ; xb += N_NODES_PER_BATCH)
        {
          bool areNonFluidNeighboursInBatch = false ;
          bool areFluidNodesInBatch         = false ;

          for (IndexType xo = 0 ; xo < N_NODES_PER_BATCH ; xo++)
          {
            IndexType x = xb + xo ;

            if (t (x,y,z).isFluid()) areFluidNodesInBatch = true ;

            bool hasNonFluidNeighbour = false ;

            for(int k=0; k < L::Q; k++)
            { 
              int ip = ( x + L::ex(k) + W ) % W;
              int jp = ( y + L::ey(k) + H ) % H;
              int mp = ( z + L::ez(k) + D ) % D;
              
              if (not t (ip,jp,mp).isFluid()) 
              {
                hasNonFluidNeighbour = true ;
                areNonFluidNeighboursInBatch = true ;
              }
            }

            if (t (x,y,z).isFluidOnlyFragment() && hasNonFluidNeighbour)
            {
              cout << "ERROR at " << x << "," << y << "," << z
                   << " : node marked as FLUID_ONLY but has non-fluid neighbours\n" ;
            }
            if (not t (x,y,z).isFluidOnlyFragment() && not hasNonFluidNeighbour)
            {
              cout << "ERROR at " << x << "," << y << "," << z
                   << " : node marked as NON_FLUID but has only fluid neighbours\n" ;
            }
          }

          if (isBatchFluid (xb,y,z) && areNonFluidNeighboursInBatch)
          {
            cout << "ERROR at " << xb << "," << y << "," << z
                 << " : batch marked as FLUID_ONLY but nodes have non-fluid neighbours\n" ;
          }
          if ((not isBatchFluid (xb,y,z)) && not areNonFluidNeighboursInBatch)
          {
            cout << "ERROR at " << xb << "," << y << "," << z
                 << " : batch marked as NON_FLUID but nodes have only fluid neighbours\n" ;
          }
          if (isBatchSolid (xb,y,z) && areFluidNodesInBatch)
          {
            cout << "ERROR at " << xb << "," << y << "," << z
                 << " : batch marked as SOLID_ONLY but has fluid nodes\n" ;
          }
          if ((not isBatchSolid (xb,y,z)) && (not areFluidNodesInBatch))
          {
            cout << "ERROR at " << xb << "," << y << "," << z
                 << " : batch marked as NOT SOLID_ONLY but has only solid nodes\n" ;
          }
          if (isBatchMixed (xb,y,z) && (not areFluidNodesInBatch))
          {
            cout << "ERROR at " << xb << "," << y << "," << z
                 << " : batch marked as MIXED but there are NO fluid nodes in batch\n" ;
          }
          if (isBatchMixed (xb,y,z) && not areNonFluidNeighboursInBatch)
          {
            cout << "ERROR at " << xb << "," << y << "," << z
                 << " : batch marked as MIXED but there are NO non-fluid neighbours in batch\n" ;
          }
        }
    cout << "Node masks checked.\n" ;
  }


  struct Statistics
  {
    size_t nNodes               = 0 ;
    size_t nFluidNodes          = 0 ;
    size_t nBoundaryNodes       = 0 ;
    size_t nSolidNodes          = 0 ;
    size_t nComputationalNodes  = 0 ;
    size_t nChunks              = 0 ;
    size_t nFluidChunks         = 0 ;
    size_t nSolidChunks         = 0 ;
    size_t nMixedChunks         = 0 ;
    size_t nComputationalChunks = 0 ;
    size_t nFluidChunksAfterMixed = 0 ;
    size_t nMixedChunksBetweenFluid = 0 ;
  } ;

  Statistics computeStatistics()
  {
    Statistics s ;

    const IndexType W = w() ;
    const IndexType H = h() ;
    const IndexType D = d() ;

    s.nNodes = nNodes() ;

    for (IndexType z=0 ; z < D ; z ++) 
      for (IndexType y=0 ; y < H ; y ++)
        for (IndexType x=0 ; x < W ; x++)
        {
          NodeType nt = t (x,y,z) ;
               if (nt.isFluidOnlyFragment()) s.nFluidNodes    ++ ;
          else if (nt.isFluid())             s.nBoundaryNodes ++ ;
          else if (nt.isSolid())             s.nSolidNodes    ++ ;
        }

    s.nComputationalNodes = s.nFluidNodes + s.nBoundaryNodes ;

    assert (s.nComputationalNodes + s.nSolidNodes == s.nNodes) ;


    s.nChunks = nNodes() / N_NODES_PER_BATCH ;

    for (IndexType z=0 ; z < D ; z++)
      for (IndexType y=0 ; y < H ; y++)
        for (IndexType xb=0 ; xb < W ; xb += N_NODES_PER_BATCH)
        {
          const auto ct = getBatchType (xb,y,z) ;

               if (isBatchFluid (ct))  s.nFluidChunks ++ ;
          else if (isBatchSolid (ct))  s.nSolidChunks ++ ;
          else if (isBatchMixed (ct)) 
          {
            s.nMixedChunks ++ ;

            IndexType xNext = (xb + N_NODES_PER_BATCH) % W ;
            if (isBatchFluid (xNext,y,z))
            {
              s.nFluidChunksAfterMixed ++ ;

              IndexType xPrev = (xb - N_NODES_PER_BATCH + W) % W ;
              if (isBatchFluid (xPrev,y,z)) s.nMixedChunksBetweenFluid ++ ;
            }
          }
        }

    s.nComputationalChunks = s.nFluidChunks + s.nMixedChunks ;

    assert (s.nComputationalChunks + s.nSolidChunks == s.nChunks) ;

    return s ;
  }

  private:
  
    void computeBoundaryBatch 
    (
      IndexType xp, IndexType y, IndexType z,
      IndexType W, IndexType H, IndexType D,
      size_t idx_x00_y00_z00,
      V * ptrSrc, V * ptrDst,
      int direction,
      T fx, T fy, T fz
    ) ;
    void computeFluidBatch 
    (
      IndexType xp, IndexType y, IndexType z,
      IndexType W, IndexType H, IndexType D,
      size_t idx_x00_y00_z00,
      V * ptrSrc, V * ptrDst,
      int direction,
      T fx, T fy, T fz
    ) ;

    DataType fx_, fy_, fz_ ;
} ;



template <class L, class T>
LBMTau1Vector3D <L,T>::
LBMTau1Vector3D (unsigned width, unsigned height, unsigned depth)
{
  using std::cout ;

  direction_ = 0 ;
  fx_ = fy_ = fz_ = 0 ;

  w_ = width ;
  h_ = height ;
  d_ = depth ;
  nNodes_ = (size_t)(w_) * h_ * d_ ;

  cout << "nNodes = " << nNodes() << "\n" ;

  cout << "sizeof (DataType): " << sizeof (DataType) << "\n" ;

  constexpr unsigned nValuesPerNode = 4 ; // Velocity 3D and density
  constexpr unsigned nCopiesPerNode = 2 ; // Input (read) and output (write).
  constexpr unsigned nBytesPerNode = nValuesPerNode * nCopiesPerNode * sizeof (T) ;

  cout << "nBytesPerNode = " << nBytesPerNode << "\n" ;
  const size_t nRequiredBytes = nBytesPerNode * nNodes() ;
  cout << "nRequiredBytes = " << nRequiredBytes << "\n" ;

#ifdef USE_LARGE_PAGES
  long alignSize = 4 * 1024*1024 ; // Universal large page size.
#else
  long pageSize = sysconf(_SC_PAGESIZE) ;
  if (pageSize < 0) 
  {
    cout << "ERROR in sysconf(_SC_PAGESIZE)\n" ;
    perror ("") ;
    exit (-1) ;
  }
  cout << "PAGE_SIZE = " << pageSize << "\n" ;

  long alignSize = pageSize ; // Universal large page size.
#endif

  size_t nBytesToAlloc = nRequiredBytes + alignmentSecondCopy * nBytesPerNode + alignSize ;
  int err = posix_memalign ((void**)(& dataPtr_), alignSize, nBytesToAlloc) ;
  if (0 != err)
  {
    cout << "ERROR in posix_memalign(...), can not allocate " << nRequiredBytes << " bytes.\n" ;
    cout << "Error code : " << err << " " << std::strerror (err) << "\n" ;
    exit (-2) ;
  }

#ifdef USE_LARGE_PAGES
  err = madvise (dataPtr_, nBytesToAlloc, MADV_SEQUENTIAL | MADV_HUGEPAGE) ;
  if (0 != err)
  {
    cout << "ERROR in madvise(...)\n" ;
    cout << "Error code : " << err << " " << std::strerror (err) << "\n" ;
    exit (-2) ;
  }
#endif

  nodeTypes_.resize (nNodes()) ;
  cout <<"nodeTypes ptr = " << (void*)(nodeTypes_.data()) << "\n" ;
  cout << "dataPtr = " << dataPtr() << "\n" ;
  cout << "secondCopyOffset = " << (void*)(secondCopyOffset()) << "\n" ;

  allocateNodeMask() ;

  resetAllNodesToFluid() ;
}

//
// Memory copy only kernel.
//
template <class L, class T>
void
LBMTau1Vector3D <L,T>::
iterCopyOnly (int direction)
{
  const IndexType W = w() ;
  const IndexType H = h() ;
  const IndexType D = d() ;

  const IndexType offset = secondCopyOffset() ;
  const IndexType offsetSrc =      direction  * offset ;
  const IndexType offsetDst = (1 - direction) * offset ;

  V * __restrict__ vPtr = dataPtr() ;

  constexpr unsigned TS = TILE_SIZE ;
  #pragma omp parallel for 
  for (IndexType zp=0 ; zp < D ; zp += TS) 
    for (IndexType yp=0 ; yp < H ; yp += TS)
      for (IndexType zo=0 ; zo < TS ; zo++) 
      for (IndexType yo=0 ; yo < TS ; yo++) 
      {
        IndexType z = zp + zo ;
        IndexType y = yp + yo ;

        IndexType idx_x00_y00_z00 = y * W + z * W * H ;  
        V * ptrSrc = vPtr + offsetSrc + idx_x00_y00_z00 ;
        V * ptrDst = vPtr + offsetDst + idx_x00_y00_z00 ;
        
        constexpr unsigned NN = N_NODES_PER_BATCH ;
        for (IndexType x = 0 ; x < W ; x += NN)   
        {
          V vx = ptrSrc [x + 0] ;
          V vy = ptrSrc [x + 1] ;
          V vz = ptrSrc [x + 2] ;
          V vr = ptrSrc [x + 3] ;

          __builtin_nontemporal_store (vx, ptrDst + x + 0) ;
          __builtin_nontemporal_store (vy, ptrDst + x + 1) ;
          __builtin_nontemporal_store (vz, ptrDst + x + 2) ;
          __builtin_nontemporal_store (vr, ptrDst + x + 3) ;
        }
      }
}

//
// Reference kernel, single thread.
//
template <class L, class T>
void
LBMTau1Vector3D <L,T>::
iterReference (int direction)
{
  constexpr T one = 1 ;
  constexpr T oph = 1.5 ;
  constexpr T three = 3. ;
  constexpr T fpf = 4.5 ;

  const T fx = fx_, fy = fy_, fz = fz_ ;

  const IndexType W = w() ;
  const IndexType H = h() ;
  const IndexType D = d() ;

  for (IndexType z=0 ; z < D ; z ++) 
    for (IndexType y=0 ; y < H ; y ++)
      for (IndexType x=0 ; x < W ; x++)
      {
        if (t (x,y,z).isFluid())
        {
          T rn = 0 ;
          T un = 0 ;
          T vn = 0 ;
          T xn = 0 ;  

          for(int k=0; k < L::Q ; k++)
          { 
            int ip = ( x + L::ex(k) + W ) % W;
            int jp = ( y + L::ey(k) + H ) % H;
            int mp = ( z + L::ez(k) + D ) % D;
            int ik = L::inv (k) ;

            T f ;
            if (t (ip,jp,mp).isFluid())
            {
              T r_,vx,vy,vz ;
              r_ = get_r  (ip,jp,mp, direction) ;
              vx = get_vx (ip,jp,mp, direction) + fx ;   
              vy = get_vy (ip,jp,mp, direction) + fy ;
              vz = get_vz (ip,jp,mp, direction) + fz ; 
              T vel = L::ex(ik)*vx+L::ey(ik)*vy+L::ez(ik)*vz ;
              f = L::w(ik)*r_*(one- (vx*vx+vy*vy+vz*vz)*oph+vel*three+vel*vel*fpf);
            } else
            {
              f = L::w(ik) * get_r (x,y,z, direction) ;
            }

            rn += f;
            un += L::ex(ik)*f;
            vn += L::ey(ik)*f;
            xn += L::ez(ik)*f;
          }

          const T invR = one / rn ;
          un *= invR ;
          vn *= invR ;
          xn *= invR ;

          set_r  (rn, x,y,z, (direction+1)%2) ;
          set_vx (un, x,y,z, (direction+1)%2) ;
          set_vy (vn, x,y,z, (direction+1)%2) ;
          set_vz (xn, x,y,z, (direction+1)%2) ;
        }
#ifdef FILL_SOLID_NODES_WITH_NAN
        else
        {
          set_r  (NAN, x,y,z, (direction+1)%2) ;
          set_vx (NAN, x,y,z, (direction+1)%2) ;
          set_vy (NAN, x,y,z, (direction+1)%2) ;
          set_vz (NAN, x,y,z, (direction+1)%2) ;
        }
#endif
      }
}

//
//  Automatically vectorized kernels.
//
template <class L, class T>
inline
void
LBMTau1Vector3D <L,T>::
computeBoundaryBatch 
(
  IndexType xp, IndexType y, IndexType z,
  IndexType W, IndexType H, IndexType D,
  size_t idx_x00_y00_z00,
  V * ptrSrc, V * ptrDst,
  int direction,
  T fx, T fy, T fz  
)
{
  constexpr unsigned NN = N_NODES_PER_BATCH ;
  constexpr T one = 1 ;
  constexpr T oph = 1.5 ;
  constexpr T three = 3. ;
  constexpr T fpf = 4.5 ;

  V vvx = {(T)(0)} ;
  V vvy = {(T)(0)} ;
  V vvz = {(T)(0)} ;
  V vr  = {(T)(0)} ;

  using NT = NodeType::DType ;
  using VT __attribute__ ((aligned (sizeof(NT)))) 
    __attribute__ ((vector_size (VSize * sizeof (NT)))) = NT ;
  using I = typename std::conditional <sizeof (int32_t) == sizeof (T), 
        int32_t, int64_t>::type ;
  using VI __attribute__ ((vector_size (VSize * sizeof (I)))) = I ;

  NT * pNT = (NT *)(nodeTypes_.data()) ;

  // BEWARE: Assumes that direction 0 is (0,0,0).
  VT nnt_next = *(VT *)(pNT + idx_x00_y00_z00 + xp) ;
  V vx_next   = ptrSrc [idx_x00_y00_z00 + xp + 0] ;
  V vy_next   = ptrSrc [idx_x00_y00_z00 + xp + 1] ;
  V vz_next   = ptrSrc [idx_x00_y00_z00 + xp + 2] ;
  V  r_next   = ptrSrc [idx_x00_y00_z00 + xp + 3] ;
  V r_current = r_next ;
#ifdef FILL_SOLID_NODES_WITH_NAN
  VT nnt_current = nnt_next ;
#endif

  #pragma unroll
  for(int k=0; k < (L::Q - 1) ; k++)
  {
    const int ik = L::inv (k) ;

    const int k_next = k+1 ;  
    const int jp = ( y  + L::ey(k_next) + H ) % H;
    const int mp = ( z  + L::ez(k_next) + D ) % D;

    IndexType idx_0 = jp * W + mp * W * H + xp ;

    VT nnt = nnt_next ;
    V  vx  =  vx_next ;
    V  vy  =  vy_next ;
    V  vz  =  vz_next ;
    V  r_  =   r_next ;

    nnt_next = *(VT *)(pNT + idx_0) ;
     vx_next = ptrSrc [idx_0 + 0] ;
     vy_next = ptrSrc [idx_0 + 1] ;
     vz_next = ptrSrc [idx_0 + 2] ;
      r_next = ptrSrc [idx_0 + 3] ;

    T svx, svy, svz, sr ;
    NT snnt ;

    if (-1 == L::ex(k_next))
    {
      int ip = (xp - 1 + W) % W ;
      IndexType nghbrTypeIdx = index (ip,jp,mp, 0) ;

      svx = get_vx (ip,jp,mp, direction) ;
      svy = get_vy (ip,jp,mp, direction) ;
      svz = get_vz (ip,jp,mp, direction) ;
      sr  = get_r  (ip,jp,mp, direction) ;
      snnt= t (nghbrTypeIdx).packedValue() ;
    }
    if (1 == L::ex(k_next))
    {
      int ip = (xp + NN) % W ;
      IndexType nghbrTypeIdx = index (ip,jp,mp, 0) ;

      svx = get_vx (ip,jp,mp, direction) ;
      svy = get_vy (ip,jp,mp, direction) ;
      svz = get_vz (ip,jp,mp, direction) ;
      sr  = get_r  (ip,jp,mp, direction) ;
      snnt= t (nghbrTypeIdx).packedValue() ;
    }

    vx += fx ;
    vy += fy ;
    vz += fz ;

    V vel = (T)(L::ex(ik)) * vx + (T)(L::ey(ik)) * vy + (T)(L::ez(ik)) * vz ;
    V f_nghbr = (T)(L::w(ik))*r_*(one- (vx*vx+vy*vy+vz*vz)*oph+vel*three+vel*vel*fpf);
    V f_current = (T)(L::w(ik)) * r_current ;

    nnt &= 1 ;
    VI fluidMask =  __builtin_convertvector(nnt, VI) ;
    fluidMask -= 1 ;
    V f = (V)(fluidMask & (VI)(f_nghbr)) + (V)(~ fluidMask & (VI)(f_current)) ;

    vr  += f ;
    vvx += (T)(L::ex(ik))*f;
    vvy += (T)(L::ey(ik))*f;
    vvz += (T)(L::ez(ik))*f;  
    
    if (-1 == L::ex(k_next))
    {
      vx_next  = __builtin_shufflevector ( vx_next,  vx_next, -1,0,1,2) ;
      vy_next  = __builtin_shufflevector ( vy_next,  vy_next, -1,0,1,2) ;
      vz_next  = __builtin_shufflevector ( vz_next,  vz_next, -1,0,1,2) ;
      r_next   = __builtin_shufflevector (  r_next,   r_next, -1,0,1,2) ;
      nnt_next = __builtin_shufflevector (nnt_next, nnt_next, -1,0,1,2) ;

      vx_next  [0] = svx ;
      vy_next  [0] = svy ;
      vz_next  [0] = svz ;
      r_next   [0] = sr  ;
      nnt_next [0] = snnt ;
    }
    if (1 == L::ex(k_next))
    {
      vx_next  = __builtin_shufflevector (vx_next , vx_next , 1,2,3,-1) ;
      vy_next  = __builtin_shufflevector (vy_next , vy_next , 1,2,3,-1) ;
      vz_next  = __builtin_shufflevector (vz_next , vz_next , 1,2,3,-1) ;
      r_next  = __builtin_shufflevector (r_next , r_next , 1,2,3,-1) ;
      nnt_next = __builtin_shufflevector (nnt_next, nnt_next, 1,2,3,-1) ;

      vx_next  [3] = svx ;
      vy_next  [3] = svy ;
      vz_next  [3] = svz ;
      r_next   [3] = sr  ;
      nnt_next [3] = snnt ;
    }
  }
  {
    constexpr int k = L::Q - 1 ;
    constexpr int ik = L::inv (k) ;

    VT nnt = nnt_next ;
    V  vx  =  vx_next ;
    V  vy  =  vy_next ;
    V  vz  =  vz_next ;
    V  r_  =   r_next ;

    vx += fx ;
    vy += fy ;
    vz += fz ;

    V vel = (T)(L::ex(ik)) * vx + (T)(L::ey(ik)) * vy + (T)(L::ez(ik)) * vz ;
    V f_nghbr = (T)(L::w(ik))*r_*(one- (vx*vx+vy*vy+vz*vz)*oph+vel*three+vel*vel*fpf);
    V f_current = (T)(L::w(ik)) * r_current ;

    nnt &= 1 ;
    VI fluidMask =  __builtin_convertvector(nnt, VI) ;
    fluidMask -= 1 ;
    V f = (V)(fluidMask & (VI)(f_nghbr)) + (V)(~ fluidMask & (VI)(f_current)) ;

    vr  += f ;
    vvx += (T)(L::ex(ik))*f;
    vvy += (T)(L::ey(ik))*f;
    vvz += (T)(L::ez(ik))*f;
  }

  V vinvR = (T)(1) / vr ;
  vvx *= vinvR ;
  vvy *= vinvR ;
  vvz *= vinvR ;

#ifdef FILL_SOLID_NODES_WITH_NAN
  {
    nnt_current &= 1 ;
    VI fluidMask =  __builtin_convertvector(nnt_current, VI) ;
    fluidMask -= 1 ;
    const V vnan = {NAN,NAN,NAN,NAN} ;
    
    vvx = (V)(fluidMask & (VI)(vvx)) + (V)(~ fluidMask & (VI)(vnan)) ;
    vvy = (V)(fluidMask & (VI)(vvy)) + (V)(~ fluidMask & (VI)(vnan)) ;
    vvz = (V)(fluidMask & (VI)(vvz)) + (V)(~ fluidMask & (VI)(vnan)) ;
    vr  = (V)(fluidMask & (VI)(vr )) + (V)(~ fluidMask & (VI)(vnan)) ;
  }
#endif

  __builtin_nontemporal_store (vvx, ptrDst + 0) ;
  __builtin_nontemporal_store (vvy, ptrDst + 1) ;
  __builtin_nontemporal_store (vvz, ptrDst + 2) ;
  __builtin_nontemporal_store (vr , ptrDst + 3) ;
} ;

template <class L, class T>
inline
void
LBMTau1Vector3D <L,T>::
computeFluidBatch 
(
  IndexType xp, IndexType y, IndexType z,
  IndexType W, IndexType H, IndexType D,
  size_t idx_x00_y00_z00,
  V * ptrSrc, V * ptrDst,
  int direction,
  T fx, T fy, T fz  
)
{
  constexpr unsigned NN = N_NODES_PER_BATCH ;
  constexpr T one = 1 ;
  constexpr T oph = 1.5 ;
  constexpr T three = 3. ;
  constexpr T fpf = 4.5 ;

  V vvx = {(T)(0)} ;
  V vvy = {(T)(0)} ;
  V vvz = {(T)(0)} ;
  V vr  = {(T)(0)} ;

  // BEWARE: Assumes that direction 0 is (0,0,0).
  V vx_next   = ptrSrc [idx_x00_y00_z00 + xp + 0] ;
  V vy_next   = ptrSrc [idx_x00_y00_z00 + xp + 1] ;
  V vz_next   = ptrSrc [idx_x00_y00_z00 + xp + 2] ;
  V  r_next   = ptrSrc [idx_x00_y00_z00 + xp + 3] ;

  #pragma unroll
  for(int k=0; k < (L::Q - 1) ; k++)
  {
    const int ik = L::inv (k) ;

    const int k_next = k+1 ;  
    const int jp = ( y  + L::ey(k_next) + H ) % H;
    const int mp = ( z  + L::ez(k_next) + D ) % D;

    IndexType idx_0 = jp * W + mp * W * H + xp ;

    V  vx  =  vx_next ;
    V  vy  =  vy_next ;
    V  vz  =  vz_next ;
    V  r_  =   r_next ;

     vx_next = ptrSrc [idx_0 + 0] ;
     vy_next = ptrSrc [idx_0 + 1] ;
     vz_next = ptrSrc [idx_0 + 2] ;
      r_next = ptrSrc [idx_0 + 3] ;

    T svx, svy, svz, sr ;

    if (-1 == L::ex(k_next))
    {
      int ip = (xp - 1 + W) % W ;

      svx = get_vx (ip,jp,mp, direction) ;
      svy = get_vy (ip,jp,mp, direction) ;
      svz = get_vz (ip,jp,mp, direction) ;
      sr  = get_r  (ip,jp,mp, direction) ;
    }
    if (1 == L::ex(k_next))
    {
      int ip = (xp + NN) % W ;

      svx = get_vx (ip,jp,mp, direction) ;
      svy = get_vy (ip,jp,mp, direction) ;
      svz = get_vz (ip,jp,mp, direction) ;
      sr  = get_r  (ip,jp,mp, direction) ;
    }

    vx += fx ;
    vy += fy ;
    vz += fz ;

    V vel = (T)(L::ex(ik)) * vx + (T)(L::ey(ik)) * vy + (T)(L::ez(ik)) * vz ;
    V f   = (T)(L::w(ik))*r_*(one- (vx*vx+vy*vy+vz*vz)*oph+vel*three+vel*vel*fpf);

    vr  += f ;
    vvx += (T)(L::ex(ik))*f;
    vvy += (T)(L::ey(ik))*f;
    vvz += (T)(L::ez(ik))*f;  
    
    if (-1 == L::ex(k_next))
    {
      vx_next  = __builtin_shufflevector ( vx_next,  vx_next, -1,0,1,2) ;
      vy_next  = __builtin_shufflevector ( vy_next,  vy_next, -1,0,1,2) ;
      vz_next  = __builtin_shufflevector ( vz_next,  vz_next, -1,0,1,2) ;
      r_next   = __builtin_shufflevector (  r_next,   r_next, -1,0,1,2) ;

      vx_next  [0] = svx ;
      vy_next  [0] = svy ;
      vz_next  [0] = svz ;
      r_next   [0] = sr  ;
    }
    if (1 == L::ex(k_next))
    {
      vx_next  = __builtin_shufflevector (vx_next , vx_next , 1,2,3,-1) ;
      vy_next  = __builtin_shufflevector (vy_next , vy_next , 1,2,3,-1) ;
      vz_next  = __builtin_shufflevector (vz_next , vz_next , 1,2,3,-1) ;
      r_next   = __builtin_shufflevector (r_next  , r_next  , 1,2,3,-1) ;

      vx_next  [3] = svx ;
      vy_next  [3] = svy ;
      vz_next  [3] = svz ;
      r_next   [3] = sr  ;
    }
  }
  {
    constexpr int k = L::Q - 1 ;
    constexpr int ik = L::inv (k) ;

    V  vx  =  vx_next ;
    V  vy  =  vy_next ;
    V  vz  =  vz_next ;
    V  r_  =   r_next ;

    vx += fx ;
    vy += fy ;
    vz += fz ;

    V vel = (T)(L::ex(ik)) * vx + (T)(L::ey(ik)) * vy + (T)(L::ez(ik)) * vz ;
    V f   = (T)(L::w(ik))*r_*(one- (vx*vx+vy*vy+vz*vz)*oph+vel*three+vel*vel*fpf);

    vr  += f ;
    vvx += (T)(L::ex(ik))*f;
    vvy += (T)(L::ey(ik))*f;
    vvz += (T)(L::ez(ik))*f;
  }

  V vinvR = (T)(1) / vr ;
  vvx *= vinvR ;
  vvy *= vinvR ;
  vvz *= vinvR ;

  __builtin_nontemporal_store (vvx, ptrDst + 0) ;
  __builtin_nontemporal_store (vvy, ptrDst + 1) ;
  __builtin_nontemporal_store (vvz, ptrDst + 2) ;
  __builtin_nontemporal_store (vr , ptrDst + 3) ;
} ;

template <class L, class T>
void
LBMTau1Vector3D <L,T>::
iterVectorizedLoop (int direction)
{
  const T fx = fx_, fy = fy_, fz = fz_ ;

  const IndexType W = w() ;
  const IndexType H = h() ;
  const IndexType D = d() ;

  const IndexType offset = secondCopyOffset() ;
  const IndexType offsetSrc =      direction  * offset ;
  const IndexType offsetDst = (1 - direction) * offset ;

  V * __restrict__ vPtr = dataPtr() ;

  constexpr unsigned TS = TILE_SIZE ;
  #pragma omp parallel for 
  for (IndexType zp=0 ; zp < D ; zp += TS) 
    for (IndexType yp=0 ; yp < H ; yp += TS)
      for (IndexType zo=0 ; zo < TS ; zo++) 
      for (IndexType yo=0 ; yo < TS ; yo++) 
      {
        IndexType z = zp + zo ;
        IndexType y = yp + yo ;

        constexpr unsigned NN = N_NODES_PER_BATCH ;

        // Node masks load for full row:
        NodeMask nodeMasks [nMasksPerRow_] ;
        IndexType idxRow = z * H + y ;
        IndexType idxNodeMask = nMasksPerRow_ * idxRow ;
        for (unsigned ii=0 ; ii < nMasksPerRow_ ; ii++)
        {
          nodeMasks [ii] = nodeMask_ [idxNodeMask + ii] ;
        }

        V * ptrSrc = vPtr + offsetSrc ;
        IndexType idx_x00_y00_z00 = y * W + z * W * H ;  

        unsigned nextBatchType = getBatchType (0,y,z) ;
        unsigned prevBatchType = getBatchType (W-NN,y,z) ;
        unsigned batchType = prevBatchType ; // For the first iteration of the loop.
        for (unsigned mIdx=0 ; mIdx < nMasksPerRow_ ; mIdx ++)
        {
          const IndexType xBegin = (IndexType)(mIdx) * N_NODES_PER_MASK ;
          // TODO: Replace std::min with separate code for the last iteration?
          const IndexType xEnd = std::min (xBegin + N_NODES_PER_MASK, W) ;

          for (IndexType xp = xBegin ; xp < xEnd ; xp += NN)    
          {
            prevBatchType = batchType ;
            batchType = nextBatchType ;
            IndexType xNext = (xp + NN) % W ; // First node after current batch.
            nextBatchType = getBatchType (xNext,y,z) ;

            if (isBatchFluid (batchType))
            {
              computeBoundaryBatch 
              (
                xp,y,z, W,H,D, idx_x00_y00_z00, 
                ptrSrc, vPtr + offsetDst + idx_x00_y00_z00 + xp,
                direction, fx,fy,fz
              ) ;
            }
            else if (isBatchMixed (batchType)) // Process nodes separately.
            {
              computeBoundaryBatch 
              (
                xp,y,z, W,H,D, idx_x00_y00_z00, 
                ptrSrc, vPtr + offsetDst + idx_x00_y00_z00 + xp,
                direction, fx,fy,fz
              ) ;
            }
          }
        }
      }
}

//
//  Manually optimized kernels.
//
template <class L, class T>
void
LBMTau1Vector3D <L,T>::
iterD3Q19 (int direction)
{
  if constexpr (3 == L::D  &&  19 == L::Q)
  {
  using std::cout ;

  constexpr T one = 1 ;
  constexpr T oph = 1.5 ;
  constexpr T three = 3. ;
  constexpr T fpf = 4.5 ;

  constexpr V vone   = {one,one,one,one} ;
  constexpr V voph   = {oph,oph,oph,oph} ;
  constexpr V vthree = {three,three,three,three} ;
  constexpr V vfpf   = {fpf,fpf,fpf,fpf} ;
  constexpr V va     = {L::a,L::a,L::a,L::a} ;
  constexpr V vb     = {L::b,L::b,L::b,L::b} ;
  constexpr V vc     = {L::c,L::c,L::c,L::c} ;

  const V vfx    = {fx_,fx_,fx_,fx_} ;
  const V vfy    = {fy_,fy_,fy_,fy_} ;
  const V vfz    = {fz_,fz_,fz_,fz_} ;

  const IndexType W = w() ;
  const IndexType H = h() ;
  const IndexType D = d() ;

  const IndexType offset = secondCopyOffset() ;
  const IndexType offsetSrc =      direction  * offset ;
  const IndexType offsetDst = (1 - direction) * offset ;

  V * __restrict__ vPtr = dataPtr() ;

  constexpr unsigned TS = TILE_SIZE ;
  #pragma omp parallel for
  for (IndexType zp=0 ; zp < D ; zp += TS) 
    for (IndexType yp=0 ; yp < H ; yp += TS)
      for (IndexType zo=0 ; zo < TS ; zo++) 
      for (IndexType yo=0 ; yo < TS ; yo++) 
      {
        IndexType z = zp + zo ;
        IndexType y = yp + yo ;

        constexpr unsigned NN = N_NODES_PER_BATCH ;

        // Node masks load for full row:
        NodeMask nodeMasks [nMasksPerRow_] ;
        IndexType idxRow = z * H + y ;
        IndexType idxNodeMask = nMasksPerRow_ * idxRow ;
        for (unsigned ii=0 ; ii < nMasksPerRow_ ; ii++)
        {
          nodeMasks [ii] = nodeMask_ [idxNodeMask + ii] ;
        }

        IndexType ym1 = (y - 1 + H) % H ;
        IndexType zm1 = (z - 1 + D) % D ;
        IndexType yp1 = (y + 1) % H ;
        IndexType zp1 = (z + 1) % D ;

        // Row-major indices of first elements in rows.
        IndexType idx_x00_ym1_zm1 = ym1 * W + zm1 * W * H ;
        IndexType idx_x00_y00_zm1 = y   * W + zm1 * W * H ;
        IndexType idx_x00_yp1_zm1 = yp1 * W + zm1 * W * H ;

        IndexType idx_x00_ym1_z00 = ym1 * W + z   * W * H ;
        IndexType idx_x00_y00_z00 = y   * W + z   * W * H ;  
        IndexType idx_x00_yp1_z00 = yp1 * W + z   * W * H ;  

        IndexType idx_x00_ym1_zp1 = ym1 * W + zp1 * W * H ;
        IndexType idx_x00_y00_zp1 = y   * W + zp1 * W * H ;
        IndexType idx_x00_yp1_zp1 = yp1 * W + zp1 * W * H ;

        V * ptrSrc = vPtr + offsetSrc ;

        V vf12_x00 ;
        {
          // Edge at ym1,zm1.

          V vx00 = ptrSrc [idx_x00_ym1_zm1 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_ym1_zm1 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_ym1_zm1 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_ym1_zm1 + 0 + 3] * 1   ;
        
          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;

          V v12_x00 = + vy00 + vz00 ;

          vf12_x00 = r00 * (vfCommon_x00 + v12_x00*vthree + v12_x00*v12_x00*vfpf) * vc ;
        }
        V vf6_x00, vf10_x00, vf15_x00 ;
        T f10_xm1 ;
        {
          // Edge at y00,zm1.

          V vx00 = ptrSrc [idx_x00_y00_zm1 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_y00_zm1 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_y00_zm1 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_y00_zm1 + 0 + 3] * 1   ;
        
          V vxm1 = ptrSrc [idx_x00_y00_zm1 + W - NN + 0] + vfx ;
          V vym1 = ptrSrc [idx_x00_y00_zm1 + W - NN + 1] + vfy ;
          V vzm1 = ptrSrc [idx_x00_y00_zm1 + W - NN + 2] + vfz ;
          V  rm1 = ptrSrc [idx_x00_y00_zm1 + W - NN + 3] * 1 ;
        
          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
          V vfCommon_xm1 = vone - (vxm1*vxm1 + vym1*vym1 + vzm1*vzm1) * voph ;

          V  v6_x00 = + vz00 ;
          V v10_x00 = + vx00 + vz00 ;
          V v10_xm1 = + vxm1 + vzm1 ;
          V v15_x00 = - vx00 + vz00 ;

           vf6_x00 = r00 * (vfCommon_x00 +  v6_x00*vthree +  v6_x00* v6_x00*vfpf) * vb ;
          vf10_x00 = r00 * (vfCommon_x00 + v10_x00*vthree + v10_x00*v10_x00*vfpf) * vc ;
          V
          vf10_xm1 = rm1 * (vfCommon_xm1 + v10_xm1*vthree + v10_xm1*v10_xm1*vfpf) * vc ;
          f10_xm1 = vf10_xm1 [VSize - 1] ;
          vf15_x00 = r00 * (vfCommon_x00 + v15_x00*vthree + v15_x00*v15_x00*vfpf) * vc ;
        }
        V vf17_x00 ;
        {
          // Edge at yp1,zm1.

          V vx00 = ptrSrc [idx_x00_yp1_zm1 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_yp1_zm1 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_yp1_zm1 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_yp1_zm1 + 0 + 3] * 1   ;
        
          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;

          V v17_x00 = - vy00 + vz00 ;

          vf17_x00 = r00 * (vfCommon_x00 + v17_x00*vthree + v17_x00*v17_x00*vfpf) * vc ;
        }
        V vf4_x00, vf8_x00, vf13_x00 ;
        T f8_xm1 ;
        {
          // Edge at ym1,z00.

          V vx00 = ptrSrc [idx_x00_ym1_z00 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_ym1_z00 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_ym1_z00 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_ym1_z00 + 0 + 3] * 1   ;

          V vxm1 = ptrSrc [idx_x00_ym1_z00 + W - NN + 0] + vfx ;
          V vym1 = ptrSrc [idx_x00_ym1_z00 + W - NN + 1] + vfy ;
          V vzm1 = ptrSrc [idx_x00_ym1_z00 + W - NN + 2] + vfz ;
          V  rm1 = ptrSrc [idx_x00_ym1_z00 + W - NN + 3] * 1 ;
        
          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
          V vfCommon_xm1 = vone - (vxm1*vxm1 + vym1*vym1 + vzm1*vzm1) * voph ;

          V  v4_x00 = + vy00 ;
          V  v8_x00 = + vx00 + vy00 ;
          V  v8_xm1 = + vxm1 + vym1 ;
          V v13_x00 = - vx00 + vy00 ;

          vf4_x00 = r00 * (vfCommon_x00 + v4_x00*vthree + v4_x00*v4_x00*vfpf) * vb ;
          vf8_x00 = r00 * (vfCommon_x00 + v8_x00*vthree + v8_x00*v8_x00*vfpf) * vc ;
          V
          vf8_xm1 = rm1 * (vfCommon_xm1 + v8_xm1*vthree + v8_xm1*v8_xm1*vfpf) * vc ;
          f8_xm1 = vf8_xm1 [VSize - 1] ;
          vf13_x00 = r00 * (vfCommon_x00 + v13_x00*vthree + v13_x00*v13_x00*vfpf) * vc ;
        }

        V vf0_x00, vf1_x00, vf2_x00 ;
        T f2_xm1 ;
        {
          // Edge at y0,z0.
        
          V vxm1 = ptrSrc [idx_x00_y00_z00 + W - NN + 0] + vfx ;
          V vym1 = ptrSrc [idx_x00_y00_z00 + W - NN + 1] + vfy ;
          V vzm1 = ptrSrc [idx_x00_y00_z00 + W - NN + 2] + vfz ;
          V  rm1 = ptrSrc [idx_x00_y00_z00 + W - NN + 3] * 1  ;
        
          V vx00 = ptrSrc [idx_x00_y00_z00 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_y00_z00 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_y00_z00 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_y00_z00 + 0 + 3] * 1   ;
        
          V vfCommon_xm1 = vone - (vxm1*vxm1 + vym1*vym1 + vzm1*vzm1) * voph ;
          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
        
          V v0_x00 = {0,0,0,0} ;
          V v1_x00 = - vx00 ;
          V v2_x00 = + vx00 ;
          V v2_xm1 = + vxm1 ;

          vf0_x00 = r00 * (vfCommon_x00 + v0_x00*vthree + v0_x00*v0_x00*vfpf) * va ;
          vf1_x00 = r00 * (vfCommon_x00 + v1_x00*vthree + v1_x00*v1_x00*vfpf) * vb ;
          vf2_x00 = r00 * (vfCommon_x00 + v2_x00*vthree + v2_x00*v2_x00*vfpf) * vb ;
          V 
          vf2_xm1 = rm1 * (vfCommon_xm1 + v2_xm1*vthree + v2_xm1*v2_xm1*vfpf) * vb ;
          f2_xm1 = vf2_xm1 [VSize -1] ;
        }

        V vf3_x00, vf7_x00, vf14_x00 ;
        T f14_xm1 ;
        {
          // Edge at yp1,z00.

          V vx00 = ptrSrc [idx_x00_yp1_z00 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_yp1_z00 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_yp1_z00 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_yp1_z00 + 0 + 3] * 1   ;

          V vxm1 = ptrSrc [idx_x00_yp1_z00 + W - NN + 0] + vfx ;
          V vym1 = ptrSrc [idx_x00_yp1_z00 + W - NN + 1] + vfy ;
          V vzm1 = ptrSrc [idx_x00_yp1_z00 + W - NN + 2] + vfz ;
          V  rm1 = ptrSrc [idx_x00_yp1_z00 + W - NN + 3] * 1  ;

          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
          V vfCommon_xm1 = vone - (vxm1*vxm1 + vym1*vym1 + vzm1*vzm1) * voph ;

          V  v3_x00 = - vy00 ;
          V  v7_x00 = - vx00 - vy00 ;
          V v14_x00 = + vx00 - vy00 ;
          V v14_xm1 = + vxm1 - vym1 ;

           vf3_x00 = r00 * (vfCommon_x00 + v3_x00*vthree + v3_x00*v3_x00*vfpf) * vb ;
           vf7_x00 = r00 * (vfCommon_x00 + v7_x00*vthree + v7_x00*v7_x00*vfpf) * vc ;
          vf14_x00 = r00 * (vfCommon_x00 + v14_x00*vthree + v14_x00*v14_x00*vfpf) * vc ;
          V
          vf14_xm1 = rm1 * (vfCommon_xm1 + v14_xm1*vthree + v14_xm1*v14_xm1*vfpf) * vc ;
          f14_xm1 = vf14_xm1 [VSize - 1] ;
        }
        V vf18_x00 ;
        {
          // Edge at ym1,zp1.

          V vx00 = ptrSrc [idx_x00_ym1_zp1 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_ym1_zp1 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_ym1_zp1 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_ym1_zp1 + 0 + 3] * 1   ;
        
          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;

          V v18_x00 = + vy00 - vz00 ;

          vf18_x00 = r00 * (vfCommon_x00 + v18_x00*vthree + v18_x00*v18_x00*vfpf) * vc ;
        }
        V vf5_x00, vf9_x00, vf16_x00 ;
        T f16_xm1 ;
        {
          // Edge at y00,zp1.

          V vx00 = ptrSrc [idx_x00_y00_zp1 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_y00_zp1 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_y00_zp1 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_y00_zp1 + 0 + 3] * 1   ;
        
          V vxm1 = ptrSrc [idx_x00_y00_zp1 + W - NN + 0] + vfx ;
          V vym1 = ptrSrc [idx_x00_y00_zp1 + W - NN + 1] + vfy ;
          V vzm1 = ptrSrc [idx_x00_y00_zp1 + W - NN + 2] + vfz ;
          V  rm1 = ptrSrc [idx_x00_y00_zp1 + W - NN + 3] * 1  ;

          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
          V vfCommon_xm1 = vone - (vxm1*vxm1 + vym1*vym1 + vzm1*vzm1) * voph ;

          V  v5_x00 = - vz00 ;
          V  v9_x00 = - vx00 - vz00 ;
          V v16_x00 = + vx00 - vz00 ;
          V v16_xm1 = + vxm1 - vzm1 ;

          vf5_x00 = r00 * (vfCommon_x00 + v5_x00*vthree + v5_x00*v5_x00*vfpf) * vb ;
          vf9_x00 = r00 * (vfCommon_x00 + v9_x00*vthree + v9_x00*v9_x00*vfpf) * vc ;
          vf16_x00 = r00 * (vfCommon_x00 + v16_x00*vthree + v16_x00*v16_x00*vfpf) * vc ;
          V
          vf16_xm1 = rm1 * (vfCommon_xm1 + v16_xm1*vthree + v16_xm1*v16_xm1*vfpf) * vc ;
          f16_xm1 = vf16_xm1 [VSize - 1] ;
        }

        V vf11_x00 ;
        {
          // Edge at yp1,zp1.

          V vx00 = ptrSrc [idx_x00_yp1_zp1 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_yp1_zp1 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_yp1_zp1 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_yp1_zp1 + 0 + 3] * 1   ;
        
          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;

          V v11_x00 = - vy00 - vz00 ;

          vf11_x00 = r00 * (vfCommon_x00 + v11_x00*vthree + v11_x00*v11_x00*vfpf) * vc ;
        }

        for (unsigned mIdx=0 ; mIdx < nMasksPerRow_ ; mIdx ++)
        {
          const IndexType xBegin = (IndexType)(mIdx) * N_NODES_PER_MASK ;
          // TODO: Replace std::min with separate code for the last iteration.
          const IndexType xEnd = std::min (xBegin + N_NODES_PER_MASK, W) ;

          for (IndexType xp = xBegin ; xp < xEnd ; xp += NN)    
          {
            if (isBatchFluid (xp,y,z)) // Vector processing of full batch
            {
              V vvx = {0,0,0,0} ;
              V vvy = {0,0,0,0} ;
              V vvz = {0,0,0,0} ;
              V vr  = {0,0,0,0} ;

              IndexType xNext = (xp + NN) % W ; // First node after current batch.

              {
                // Edge at ym1,zm1.
              
                V vxp1 = ptrSrc [idx_x00_ym1_zm1 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_ym1_zm1 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_ym1_zm1 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_ym1_zm1 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;
              
                V v12_xp1 =        + vyp1 + vzp1 ;
              
                V vf12_xp1 = rp1 * (vfCommon_xp1 + v12_xp1*vthree + v12_xp1*v12_xp1*vfpf) * vc ;
              
                V vf12 = vf12_x00 ;
                vf12_x00 = vf12_xp1 ;

                vr  +=   vf12 ;
                vvy += + vf12 ;
                vvz += + vf12 ;
              }
              {
                // Edge at y00,zm1.
              
                V vxp1 = ptrSrc [idx_x00_y00_zm1 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_y00_zm1 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_y00_zm1 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_y00_zm1 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;
              
                V  v6_xp1 = + vzp1 ;
                V v10_xp1 = + vxp1 + vzp1 ;
                V v15_xp1 = - vxp1 + vzp1 ;
              
                V  vf6_xp1 = rp1 * (vfCommon_xp1 +  v6_xp1*vthree +  v6_xp1* v6_xp1*vfpf) * vb ;
                V vf10_xp1 = rp1 * (vfCommon_xp1 + v10_xp1*vthree + v10_xp1*v10_xp1*vfpf) * vc ;
                V vf15_xp1 = rp1 * (vfCommon_xp1 + v15_xp1*vthree + v15_xp1*v15_xp1*vfpf) * vc ;
              
                V vf6 = vf6_x00 ;
                vf6_x00 = vf6_xp1 ;

                V vf10 =  __builtin_shufflevector (vf10_x00,  vf10_x00, -1,0,1,2) ;
                vf10 [0] = f10_xm1 ;
                f10_xm1 = vf10_x00 [VSize -1] ;
                vf10_x00 = vf10_xp1 ;

                V vf15 = __builtin_shufflevector (vf15_x00, vf15_xp1, 1,2,3,4) ;
                vf15_x00 = vf15_xp1 ;

                vr  +=   vf6 + vf10 + vf15 ;
                vvx += + vf10 - vf15 ;
                vvz += + vf6 + vf10 + vf15 ;
              }
              {
                // Edge at yp1,zm1.

                V vxp1 = ptrSrc [idx_x00_yp1_zm1 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_yp1_zm1 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_yp1_zm1 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_yp1_zm1 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                V v17_xp1 =        - vyp1 + vzp1 ;

                V vf17_xp1 = rp1 * (vfCommon_xp1 + v17_xp1*vthree + v17_xp1*v17_xp1*vfpf) * vc ;

                V vf17 = vf17_x00 ;
                vf17_x00 = vf17_xp1 ;

                vr  +=   vf17 ;
                vvy += - vf17 ;
                vvz += + vf17 ;
              }
              {
                // Edge at ym1,z00.
              
                V vxp1 = ptrSrc [idx_x00_ym1_z00 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_ym1_z00 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_ym1_z00 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_ym1_z00 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;
              
                V  v4_xp1 = + vyp1 ;
                V  v8_xp1 = + vxp1 + vyp1 ;
                V v13_xp1 = - vxp1 + vyp1 ;
              
                V  vf4_xp1 = rp1 * (vfCommon_xp1 + v4_xp1*vthree + v4_xp1*v4_xp1*vfpf) * vb ;
                V  vf8_xp1 = rp1 * (vfCommon_xp1 + v8_xp1*vthree + v8_xp1*v8_xp1*vfpf) * vc ;
                V vf13_xp1 = rp1 * (vfCommon_xp1 + v13_xp1*vthree + v13_xp1*v13_xp1*vfpf) * vc ;
              
                V vf4 = vf4_x00 ;
                vf4_x00 = vf4_xp1 ;

                V vf8 =  __builtin_shufflevector (vf8_x00,  vf8_x00, -1,0,1,2) ;
                vf8 [0] = f8_xm1 ;
                f8_xm1 = vf8_x00 [VSize -1] ;
                vf8_x00 = vf8_xp1 ;

                V vf13 = __builtin_shufflevector (vf13_x00, vf13_xp1, 1,2,3,4) ;
                vf13_x00 = vf13_xp1 ;
              
                vr  +=   vf4 + vf8 + vf13 ;
                vvx += + vf8 - vf13 ;
                vvy += + vf4 + vf8 + vf13 ;
              }
              {
                // Edge at y00,z00.

                V vxp1 = ptrSrc [idx_x00_y00_z00 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_y00_z00 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_y00_z00 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_y00_z00 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                V v0_xp1 = {0,0,0,0} ;
                V v1_xp1 = - vxp1 ;
                V v2_xp1 = + vxp1 ;

                V  vf0_xp1 = rp1 * (vfCommon_xp1 + v0_xp1*vthree + v0_xp1*v0_xp1*vfpf) * va ;
                V  vf1_xp1 = rp1 * (vfCommon_xp1 + v1_xp1*vthree + v1_xp1*v1_xp1*vfpf) * vb ;
                V  vf2_xp1 = rp1 * (vfCommon_xp1 + v2_xp1*vthree + v2_xp1*v2_xp1*vfpf) * vb ;

                V vf0 = vf0_x00 ;
                vf0_x00 = vf0_xp1 ;

                V vf1 = __builtin_shufflevector (vf1_x00, vf1_xp1, 1,2,3,4) ;
                vf1_x00 = vf1_xp1 ;
              
                V vf2 =  __builtin_shufflevector (vf2_x00,  vf2_x00, -1,0,1,2) ;
                vf2 [0] = f2_xm1 ;
                f2_xm1 = vf2_x00 [VSize -1] ;
                vf2_x00 = vf2_xp1 ;

                vr  +=   vf0 + vf1 + vf2 ;
                vvx += - vf1 + vf2;
              }
              {
                // Edge at yp1,z00.

                V vxp1 = ptrSrc [idx_x00_yp1_z00 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_yp1_z00 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_yp1_z00 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_yp1_z00 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                V  v3_xp1 = - vyp1 ;
                V  v7_xp1 = - vxp1 - vyp1 ;
                V v14_xp1 = + vxp1 - vyp1 ;

                V  vf3_xp1 = rp1 * (vfCommon_xp1 + v3_xp1*vthree + v3_xp1*v3_xp1*vfpf) * vb ;
                V  vf7_xp1 = rp1 * (vfCommon_xp1 + v7_xp1*vthree + v7_xp1*v7_xp1*vfpf) * vc ;
                V vf14_xp1 = rp1 * (vfCommon_xp1 + v14_xp1*vthree + v14_xp1*v14_xp1*vfpf) * vc ;

                V vf3 = vf3_x00 ;
                vf3_x00 = vf3_xp1 ;
                V vf7 = __builtin_shufflevector (vf7_x00, vf7_xp1, 1,2,3,4) ;
                vf7_x00 = vf7_xp1 ;

                V vf14 =  __builtin_shufflevector (vf14_x00,  vf14_x00, -1,0,1,2) ;
                vf14 [0] = f14_xm1 ;
                f14_xm1 = vf14_x00 [VSize -1] ;
                vf14_x00 = vf14_xp1 ;

                vr  +=   vf3 + vf7 + vf14 ;
                vvx += - vf7 + vf14 ;
                vvy += - vf3 - vf7 - vf14 ;
              }
              {
                // Edge at ym1,zp1.

                V vxp1 = ptrSrc [idx_x00_ym1_zp1 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_ym1_zp1 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_ym1_zp1 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_ym1_zp1 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                V v18_xp1 = + vyp1 - vzp1 ;

                V vf18_xp1 = rp1 * (vfCommon_xp1 + v18_xp1*vthree + v18_xp1*v18_xp1*vfpf) * vc ;

                V vf18 = vf18_x00 ;
                vf18_x00 = vf18_xp1 ;

                vr  +=   vf18 ;
                vvy += + vf18 ;
                vvz += - vf18 ;
              }
              {
                // Edge at y00,zp1.

                V vxp1 = ptrSrc [idx_x00_y00_zp1 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_y00_zp1 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_y00_zp1 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_y00_zp1 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                V  v5_xp1 = - vzp1 ;
                V  v9_xp1 = - vxp1 - vzp1 ;
                V v16_xp1 = + vxp1 - vzp1 ;

                V vf5_xp1 = rp1 * (vfCommon_xp1 + v5_xp1*vthree + v5_xp1*v5_xp1*vfpf) * vb ;
                V vf9_xp1 = rp1 * (vfCommon_xp1 + v9_xp1*vthree + v9_xp1*v9_xp1*vfpf) * vc ;
                V vf16_xp1 = rp1 * (vfCommon_xp1 + v16_xp1*vthree + v16_xp1*v16_xp1*vfpf) * vc ;

                V vf5 = vf5_x00 ;
                vf5_x00 = vf5_xp1 ;
                V vf9 = __builtin_shufflevector (vf9_x00, vf9_xp1, 1,2,3,4) ;
                vf9_x00 = vf9_xp1 ;
                V vf16 =  __builtin_shufflevector (vf16_x00,  vf16_x00, -1,0,1,2) ;
                vf16 [0] = f16_xm1 ;
                f16_xm1 = vf16_x00 [VSize -1] ;
                vf16_x00 = vf16_xp1 ;

                vr  +=   vf5 + vf9 + vf16 ;
                vvx +=       - vf9 + vf16 ;
                vvz += - vf5 - vf9 - vf16 ;
              }
              {
                // Edge at yp1,zp1.

                V vxp1 = ptrSrc [idx_x00_yp1_zp1 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_yp1_zp1 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_yp1_zp1 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_yp1_zp1 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                V v11_xp1 = - vyp1 - vzp1 ;

                V vf11_xp1 = rp1 * (vfCommon_xp1 + v11_xp1*vthree + v11_xp1*v11_xp1*vfpf) * vc ;

                V vf11 = vf11_x00 ;
                vf11_x00 = vf11_xp1 ;

                vr  +=   vf11 ;
                vvy += - vf11 ;
                vvz += - vf11 ;
              }

              V vinvR = vone / vr ;
              vvx *= vinvR ;
              vvy *= vinvR ;
              vvz *= vinvR ;

              {
                V * ptrDst = vPtr + offsetDst + idx_x00_y00_z00 + xp ;

                __builtin_nontemporal_store (vvx, ptrDst + 0) ;
                __builtin_nontemporal_store (vvy, ptrDst + 1) ;
                __builtin_nontemporal_store (vvz, ptrDst + 2) ;
                __builtin_nontemporal_store ( vr, ptrDst + 3) ;
              }
            }
            else if (isBatchMixed (xp,y,z)) // Process nodes separately.
            {
              IndexType xNext = (xp + NN) % W ; // First node after current batch.
              IndexType xPrev = (xp - NN + W) % W ;
              
              //WARNING: Must precompute values for the next fluid-only batch.
              //         We can be sure that fluid-only batch is surrounded by 
              //         either fluid-only or mixed batches.
              if (isBatchFluid (xNext, y,z))
              { 
                if (isBatchFluid (xPrev, y,z))
                {
                  // Can use values precomputed in the previous batch.
                }
                else
                {
                  // Must compute values for the previous batch, too.
                  // FIXME: use values computed in the for loop placed below.
                  {
                    // Edge at y0,zm1.

                    V vx00 = ptrSrc [idx_x00_y00_zm1 + xp + 0] + vfx ;
                    V vy00 = ptrSrc [idx_x00_y00_zm1 + xp + 1] + vfy ;
                    V vz00 = ptrSrc [idx_x00_y00_zm1 + xp + 2] + vfz ;
                    V  r00 = ptrSrc [idx_x00_y00_zm1 + xp + 3] * 1   ;

                    V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
                    V v10_x00 = + vx00 + vz00 ;
                    vf10_x00 = r00 * (vfCommon_x00 + v10_x00*vthree + v10_x00*v10_x00*vfpf) * vc ;
                  }
                  {
                    // Edge at ym1,z0.

                    V vx00 = ptrSrc [idx_x00_ym1_z00 + xp + 0] + vfx ;
                    V vy00 = ptrSrc [idx_x00_ym1_z00 + xp + 1] + vfy ;
                    V vz00 = ptrSrc [idx_x00_ym1_z00 + xp + 2] + vfz ;
                    V  r00 = ptrSrc [idx_x00_ym1_z00 + xp + 3] * 1   ;

                    V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
                    V v8_x00 = + vx00 + vy00 ;
                    vf8_x00 = r00 * (vfCommon_x00 + v8_x00*vthree + v8_x00*v8_x00*vfpf) * vc ;
                  }
                  {
                    // Edge at y0,z0.

                    V vx00 = ptrSrc [idx_x00_y00_z00 + xp + 0] + vfx ;
                    V vy00 = ptrSrc [idx_x00_y00_z00 + xp + 1] + vfy ;
                    V vz00 = ptrSrc [idx_x00_y00_z00 + xp + 2] + vfz ;
                    V  r00 = ptrSrc [idx_x00_y00_z00 + xp + 3] * 1   ;

                    V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
                    V v2_x00 = + vx00 ;
                    vf2_x00 = r00 * (vfCommon_x00 + v2_x00*vthree + v2_x00*v2_x00*vfpf) * vb ;
                  }
                  {
                    // Edge at yp1,z0.

                    V vx00 = ptrSrc [idx_x00_yp1_z00 + xp + 0] + vfx ;
                    V vy00 = ptrSrc [idx_x00_yp1_z00 + xp + 1] + vfy ;
                    V vz00 = ptrSrc [idx_x00_yp1_z00 + xp + 2] + vfz ;
                    V  r00 = ptrSrc [idx_x00_yp1_z00 + xp + 3] * 1   ;

                    V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
                    V v14_x00 = + vx00 - vy00 ;
                    vf14_x00 = r00 * (vfCommon_x00 + v14_x00*vthree + v14_x00*v14_x00*vfpf) * vc ;
                  }
                  {
                    // Edge at y0,zp1.

                    V vx00 = ptrSrc [idx_x00_y00_zp1 + xp + 0] + vfx ;
                    V vy00 = ptrSrc [idx_x00_y00_zp1 + xp + 1] + vfy ;
                    V vz00 = ptrSrc [idx_x00_y00_zp1 + xp + 2] + vfz ;
                    V  r00 = ptrSrc [idx_x00_y00_zp1 + xp + 3] * 1   ;

                    V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
                    V v16_x00 = + vx00 - vz00 ;
                    vf16_x00 = r00 * (vfCommon_x00 + v16_x00*vthree + v16_x00*v16_x00*vfpf) * vc ;
                  }
                }

                {
                  // Edge at ym1,zm1.

                  V vxp1 = ptrSrc [idx_x00_ym1_zm1 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_ym1_zm1 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_ym1_zm1 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_ym1_zm1 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V v12_xp1 =        + vyp1 + vzp1 ;

                  V vf12_xp1 = rp1 * (vfCommon_xp1 + v12_xp1*vthree + v12_xp1*v12_xp1*vfpf) * vc ;
                  vf12_x00 = vf12_xp1 ;
                }
                {
                  // Edge at y00,zm1.

                  V vxp1 = ptrSrc [idx_x00_y00_zm1 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_y00_zm1 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_y00_zm1 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_y00_zm1 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V  v6_xp1 = + vzp1 ;
                  V v10_xp1 = + vxp1 + vzp1 ;
                  V v15_xp1 = - vxp1 + vzp1 ;

                  V vf6_xp1 = rp1 * (vfCommon_xp1 + v6_xp1*vthree + v6_xp1*v6_xp1*vfpf) * vb ;
                  vf6_x00 = vf6_xp1 ;
                  V vf10_xp1 = rp1 * (vfCommon_xp1 + v10_xp1*vthree + v10_xp1*v10_xp1*vfpf) * vc ;
                  f10_xm1 = vf10_x00 [VSize - 1] ;
                  vf10_x00 = vf10_xp1 ;
                  V vf15_xp1 = rp1 * (vfCommon_xp1 + v15_xp1*vthree + v15_xp1*v15_xp1*vfpf) * vc ;
                  vf15_x00 = vf15_xp1 ;
                }
                {
                  // Edge at yp1,zm1.

                  V vxp1 = ptrSrc [idx_x00_yp1_zm1 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_yp1_zm1 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_yp1_zm1 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_yp1_zm1 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V v17_xp1 =        - vyp1 + vzp1 ;

                  V vf17_xp1 = rp1 * (vfCommon_xp1 + v17_xp1*vthree + v17_xp1*v17_xp1*vfpf) * vc ;

                  vf17_x00 = vf17_xp1 ;
                }
                {
                  // Edge at ym1,z00.

                  V vxp1 = ptrSrc [idx_x00_ym1_z00 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_ym1_z00 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_ym1_z00 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_ym1_z00 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V  v4_xp1 = + vyp1 ;
                  V  v8_xp1 = + vxp1 + vyp1 ;
                  V v13_xp1 = - vxp1 + vyp1 ;

                  V vf4_xp1 = rp1 * (vfCommon_xp1 + v4_xp1*vthree + v4_xp1*v4_xp1*vfpf) * vb ;
                  vf4_x00 = vf4_xp1 ;
                  V vf8_xp1 = rp1 * (vfCommon_xp1 + v8_xp1*vthree + v8_xp1*v8_xp1*vfpf) * vc ;
                  f8_xm1 = vf8_x00 [VSize -1] ;
                  vf8_x00 = vf8_xp1 ;
                  V vf13_xp1 = rp1 * (vfCommon_xp1 + v13_xp1*vthree + v13_xp1*v13_xp1*vfpf) * vc ;
                  vf13_x00 = vf13_xp1 ;
                }
                {
                  // Edge at y0,z0.

                  V vxp1 = ptrSrc [idx_x00_y00_z00 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_y00_z00 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_y00_z00 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_y00_z00 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V v0_xp1 = {0,0,0,0} ;
                  V v1_xp1 = - vxp1 ;
                  V v2_xp1 = + vxp1 ;

                  V  vf0_xp1 = rp1 * (vfCommon_xp1 + v0_xp1*vthree + v0_xp1*v0_xp1*vfpf) * va ;
                  vf0_x00 = vf0_xp1 ;

                  V  vf1_xp1 = rp1 * (vfCommon_xp1 + v1_xp1*vthree + v1_xp1*v1_xp1*vfpf) * vb ;
                  vf1_x00 = vf1_xp1 ;

                  V  vf2_xp1 = rp1 * (vfCommon_xp1 + v2_xp1*vthree + v2_xp1*v2_xp1*vfpf) * vb ;
                  f2_xm1 = vf2_x00 [VSize -1] ;
                  vf2_x00 = vf2_xp1 ;
                }
                {
                  // Edge at yp1,z00.

                  V vxp1 = ptrSrc [idx_x00_yp1_z00 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_yp1_z00 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_yp1_z00 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_yp1_z00 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V  v3_xp1 = - vyp1 ;
                  V  v7_xp1 = - vxp1 - vyp1 ;
                  V v14_xp1 = + vxp1 - vyp1 ;

                  V vf3_xp1 = rp1 * (vfCommon_xp1 + v3_xp1*vthree + v3_xp1*v3_xp1*vfpf) * vb ;
                  V vf7_xp1 = rp1 * (vfCommon_xp1 + v7_xp1*vthree + v7_xp1*v7_xp1*vfpf) * vc ;
                  V vf14_xp1 = rp1 * (vfCommon_xp1 + v14_xp1*vthree + v14_xp1*v14_xp1*vfpf) * vc ;
                  vf3_x00 = vf3_xp1 ;
                  vf7_x00 = vf7_xp1 ;
                  f14_xm1 = vf14_x00 [VSize - 1] ;
                  vf14_x00 = vf14_xp1 ;
                }
                {
                  // Edge at ym1,zp1.

                  V vxp1 = ptrSrc [idx_x00_ym1_zp1 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_ym1_zp1 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_ym1_zp1 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_ym1_zp1 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V v18_xp1 = + vyp1 - vzp1 ;

                  V vf18_xp1 = rp1 * (vfCommon_xp1 + v18_xp1*vthree + v18_xp1*v18_xp1*vfpf) * vc ;

                  vf18_x00 = vf18_xp1 ;
                }
                {
                  // Edge at y00,zp1.

                  V vxp1 = ptrSrc [idx_x00_y00_zp1 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_y00_zp1 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_y00_zp1 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_y00_zp1 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V  v5_xp1 = - vzp1 ;
                  V  v9_xp1 = - vxp1 - vzp1 ;
                  V v16_xp1 = + vxp1 - vzp1 ;

                  V vf5_xp1 = rp1 * (vfCommon_xp1 + v5_xp1*vthree + v5_xp1*v5_xp1*vfpf) * vb ;
                  V vf9_xp1 = rp1 * (vfCommon_xp1 + v9_xp1*vthree + v9_xp1*v9_xp1*vfpf) * vc ;
                  V vf16_xp1 = rp1 * (vfCommon_xp1 + v16_xp1*vthree + v16_xp1*v16_xp1*vfpf) * vc ;

                  vf5_x00 = vf5_xp1 ;
                  vf9_x00 = vf9_xp1 ;
                  f16_xm1 = vf16_x00 [VSize - 1] ;
                  vf16_x00 = vf16_xp1 ;
                }
                {
                  // Edge at yp1,zp1.

                  V vxp1 = ptrSrc [idx_x00_yp1_zp1 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_yp1_zp1 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_yp1_zp1 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_yp1_zp1 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V v11_xp1 = - vyp1 - vzp1 ;

                  V vf11_xp1 = rp1 * (vfCommon_xp1 + v11_xp1*vthree + v11_xp1*v11_xp1*vfpf) * vc ;

                  vf11_x00 = vf11_xp1 ;
                }
              }

              computeBoundaryBatch 
              (
                xp,y,z, W,H,D, idx_x00_y00_z00, 
                ptrSrc, vPtr + offsetDst + idx_x00_y00_z00 + xp,
                direction, vfx[0],vfy[0],vfz[0]
              ) ;
            }
          }
        }
      }
  }
}

template <class L, class T>
void
LBMTau1Vector3D <L,T>::
iterD3Q27 (int direction)
{
  if constexpr (3 == L::D  &&  27 == L::Q)
  {
  using std::cout ;

  constexpr T one = 1 ;
  constexpr T oph = 1.5 ;
  constexpr T three = 3. ;
  constexpr T fpf = 4.5 ;

  constexpr V vone   = {one,one,one,one} ;
  constexpr V voph   = {oph,oph,oph,oph} ;
  constexpr V vthree = {three,three,three,three} ;
  constexpr V vfpf   = {fpf,fpf,fpf,fpf} ;
  constexpr V va     = {L::a,L::a,L::a,L::a} ;
  constexpr V vb     = {L::b,L::b,L::b,L::b} ;
  constexpr V vc     = {L::c,L::c,L::c,L::c} ;
  constexpr V vd     = {L::d,L::d,L::d,L::d} ;

  const V vfx    = {fx_,fx_,fx_,fx_} ;
  const V vfy    = {fy_,fy_,fy_,fy_} ;
  const V vfz    = {fz_,fz_,fz_,fz_} ;

  const IndexType W = w() ;
  const IndexType H = h() ;
  const IndexType D = d() ;

  const IndexType offset = secondCopyOffset() ;
  const IndexType offsetSrc =      direction  * offset ;
  const IndexType offsetDst = (1 - direction) * offset ;

  V * __restrict__ vPtr = dataPtr() ;

constexpr unsigned TS = TILE_SIZE ;
#pragma omp parallel for
  for (IndexType zp=0 ; zp < D ; zp += TS) 
    for (IndexType yp=0 ; yp < H ; yp += TS)
      for (IndexType zo=0 ; zo < TS ; zo++) 
      for (IndexType yo=0 ; yo < TS ; yo++) 
      {
        IndexType z = zp + zo ;
        IndexType y = yp + yo ;

        constexpr unsigned NN = N_NODES_PER_BATCH ;

        // Node masks load for full row:
        NodeMask nodeMasks [nMasksPerRow_] ;
        IndexType idxRow = z * H + y ;
        IndexType idxNodeMask = nMasksPerRow_ * idxRow ;
        for (unsigned ii=0 ; ii < nMasksPerRow_ ; ii++)
        {
          nodeMasks [ii] = nodeMask_ [idxNodeMask + ii] ;
        }

        IndexType ym1 = (y - 1 + H) % H ;
        IndexType zm1 = (z - 1 + D) % D ;
        IndexType yp1 = (y + 1) % H ;
        IndexType zp1 = (z + 1) % D ;

        // Row-major indices of first elements in rows.
        IndexType idx_x00_ym1_zm1 = ym1 * W + zm1 * W * H ;
        IndexType idx_x00_y00_zm1 = y   * W + zm1 * W * H ;
        IndexType idx_x00_yp1_zm1 = yp1 * W + zm1 * W * H ;

        IndexType idx_x00_ym1_z00 = ym1 * W + z   * W * H ;
        IndexType idx_x00_y00_z00 = y   * W + z   * W * H ;  
        IndexType idx_x00_yp1_z00 = yp1 * W + z   * W * H ;  

        IndexType idx_x00_ym1_zp1 = ym1 * W + zp1 * W * H ;
        IndexType idx_x00_y00_zp1 = y   * W + zp1 * W * H ;
        IndexType idx_x00_yp1_zp1 = yp1 * W + zp1 * W * H ;

        V * ptrSrc = vPtr + offsetSrc ;

        //TODO: Decrease the number of local variables through storing velocity
        //      and density of nodes (12 values for 3 neighbouring batches)
        //      instead of all distribution functions. 
        //      In this case, distribution functions can be accumulated in-fly.
        V vf12_x00, vf20_x00, vf26_x00 ;
        T f20_xm1 ;
        {
          // Edge at ym1,zm1.

          V vx00 = ptrSrc [idx_x00_ym1_zm1 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_ym1_zm1 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_ym1_zm1 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_ym1_zm1 + 0 + 3] * 1   ;
        
          V vxm1 = ptrSrc [idx_x00_ym1_zm1 + W - NN + 0] + vfx ;
          V vym1 = ptrSrc [idx_x00_ym1_zm1 + W - NN + 1] + vfy ;
          V vzm1 = ptrSrc [idx_x00_ym1_zm1 + W - NN + 2] + vfz ;
          V  rm1 = ptrSrc [idx_x00_ym1_zm1 + W - NN + 3] * 1 ;
        
          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
          V vfCommon_xm1 = vone - (vxm1*vxm1 + vym1*vym1 + vzm1*vzm1) * voph ;

          V v12_x00 =        + vy00 + vz00 ;
          V v20_x00 = + vx00 + vy00 + vz00 ;
          V v20_xm1 = + vxm1 + vym1 + vzm1 ;
          V v26_x00 = - vx00 + vy00 + vz00 ;

          vf12_x00 = r00 * (vfCommon_x00 + v12_x00*vthree + v12_x00*v12_x00*vfpf) * vc ;
          vf20_x00 = r00 * (vfCommon_x00 + v20_x00*vthree + v20_x00*v20_x00*vfpf) * vd ;
          V
          vf20_xm1 = rm1 * (vfCommon_xm1 + v20_xm1*vthree + v20_xm1*v20_xm1*vfpf) * vd ;
          f20_xm1 = vf20_xm1 [VSize - 1] ;
          vf26_x00 = r00 * (vfCommon_x00 + v26_x00*vthree + v26_x00*v26_x00*vfpf) * vd ;
        }
        V vf6_x00, vf10_x00, vf15_x00 ;
        T f10_xm1 ;
        {
          // Edge at y00,zm1.

          V vx00 = ptrSrc [idx_x00_y00_zm1 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_y00_zm1 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_y00_zm1 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_y00_zm1 + 0 + 3] * 1   ;
        
          V vxm1 = ptrSrc [idx_x00_y00_zm1 + W - NN + 0] + vfx ;
          V vym1 = ptrSrc [idx_x00_y00_zm1 + W - NN + 1] + vfy ;
          V vzm1 = ptrSrc [idx_x00_y00_zm1 + W - NN + 2] + vfz ;
          V  rm1 = ptrSrc [idx_x00_y00_zm1 + W - NN + 3] * 1 ;
        
          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
          V vfCommon_xm1 = vone - (vxm1*vxm1 + vym1*vym1 + vzm1*vzm1) * voph ;

          V  v6_x00 = + vz00 ;
          V v10_x00 = + vx00 + vz00 ;
          V v10_xm1 = + vxm1 + vzm1 ;
          V v15_x00 = - vx00 + vz00 ;

           vf6_x00 = r00 * (vfCommon_x00 +  v6_x00*vthree +  v6_x00* v6_x00*vfpf) * vb ;
          vf10_x00 = r00 * (vfCommon_x00 + v10_x00*vthree + v10_x00*v10_x00*vfpf) * vc ;
          V
          vf10_xm1 = rm1 * (vfCommon_xm1 + v10_xm1*vthree + v10_xm1*v10_xm1*vfpf) * vc ;
          f10_xm1 = vf10_xm1 [VSize - 1] ;
          vf15_x00 = r00 * (vfCommon_x00 + v15_x00*vthree + v15_x00*v15_x00*vfpf) * vc ;
        }
        V vf17_x00, vf21_x00, vf24_x00 ;
        T f24_xm1 ;
        {
          // Edge at yp1,zm1.

          V vx00 = ptrSrc [idx_x00_yp1_zm1 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_yp1_zm1 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_yp1_zm1 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_yp1_zm1 + 0 + 3] * 1   ;
        
          V vxm1 = ptrSrc [idx_x00_yp1_zm1 + W - NN + 0] + vfx ;
          V vym1 = ptrSrc [idx_x00_yp1_zm1 + W - NN + 1] + vfy ;
          V vzm1 = ptrSrc [idx_x00_yp1_zm1 + W - NN + 2] + vfz ;
          V  rm1 = ptrSrc [idx_x00_yp1_zm1 + W - NN + 3] * 1 ;
        
          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
          V vfCommon_xm1 = vone - (vxm1*vxm1 + vym1*vym1 + vzm1*vzm1) * voph ;

          V v17_x00 =        - vy00 + vz00 ;
          V v21_x00 = - vx00 - vy00 + vz00 ;
          V v24_x00 = + vx00 - vy00 + vz00 ;
          V v24_xm1 = + vxm1 - vym1 + vzm1 ;

          vf17_x00 = r00 * (vfCommon_x00 + v17_x00*vthree + v17_x00*v17_x00*vfpf) * vc ;
          vf21_x00 = r00 * (vfCommon_x00 + v21_x00*vthree + v21_x00*v21_x00*vfpf) * vd ;
          vf24_x00 = r00 * (vfCommon_x00 + v24_x00*vthree + v24_x00*v24_x00*vfpf) * vd ;
          V
          vf24_xm1 = rm1 * (vfCommon_xm1 + v24_xm1*vthree + v24_xm1*v24_xm1*vfpf) * vd ;
          f24_xm1 = vf24_xm1 [VSize - 1] ;
        }
        V vf4_x00, vf8_x00, vf13_x00 ;
        T f8_xm1 ;
        {
          // Edge at ym1,z00.

          V vx00 = ptrSrc [idx_x00_ym1_z00 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_ym1_z00 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_ym1_z00 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_ym1_z00 + 0 + 3] * 1   ;

          V vxm1 = ptrSrc [idx_x00_ym1_z00 + W - NN + 0] + vfx ;
          V vym1 = ptrSrc [idx_x00_ym1_z00 + W - NN + 1] + vfy ;
          V vzm1 = ptrSrc [idx_x00_ym1_z00 + W - NN + 2] + vfz ;
          V  rm1 = ptrSrc [idx_x00_ym1_z00 + W - NN + 3] * 1 ;
        
          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
          V vfCommon_xm1 = vone - (vxm1*vxm1 + vym1*vym1 + vzm1*vzm1) * voph ;

          V  v4_x00 = + vy00 ;
          V  v8_x00 = + vx00 + vy00 ;
          V  v8_xm1 = + vxm1 + vym1 ;
          V v13_x00 = - vx00 + vy00 ;

          vf4_x00 = r00 * (vfCommon_x00 + v4_x00*vthree + v4_x00*v4_x00*vfpf) * vb ;
          vf8_x00 = r00 * (vfCommon_x00 + v8_x00*vthree + v8_x00*v8_x00*vfpf) * vc ;
          V
          vf8_xm1 = rm1 * (vfCommon_xm1 + v8_xm1*vthree + v8_xm1*v8_xm1*vfpf) * vc ;
          f8_xm1 = vf8_xm1 [VSize - 1] ;
          vf13_x00 = r00 * (vfCommon_x00 + v13_x00*vthree + v13_x00*v13_x00*vfpf) * vc ;
        }

        V vf0_x00, vf1_x00, vf2_x00 ;
        T f2_xm1 ;
        {
          // Edge at y0,z0.
        
          V vxm1 = ptrSrc [idx_x00_y00_z00 + W - NN + 0] + vfx ;
          V vym1 = ptrSrc [idx_x00_y00_z00 + W - NN + 1] + vfy ;
          V vzm1 = ptrSrc [idx_x00_y00_z00 + W - NN + 2] + vfz ;
          V  rm1 = ptrSrc [idx_x00_y00_z00 + W - NN + 3] * 1  ;
        
          V vx00 = ptrSrc [idx_x00_y00_z00 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_y00_z00 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_y00_z00 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_y00_z00 + 0 + 3] * 1   ;
        
          V vfCommon_xm1 = vone - (vxm1*vxm1 + vym1*vym1 + vzm1*vzm1) * voph ;
          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
        
          V v0_x00 = {0,0,0,0} ;
          V v1_x00 = - vx00 ;
          V v2_x00 = + vx00 ;
          V v2_xm1 = + vxm1 ;

          vf0_x00 = r00 * (vfCommon_x00 + v0_x00*vthree + v0_x00*v0_x00*vfpf) * va ;
          vf1_x00 = r00 * (vfCommon_x00 + v1_x00*vthree + v1_x00*v1_x00*vfpf) * vb ;
          vf2_x00 = r00 * (vfCommon_x00 + v2_x00*vthree + v2_x00*v2_x00*vfpf) * vb ;
          V 
          vf2_xm1 = rm1 * (vfCommon_xm1 + v2_xm1*vthree + v2_xm1*v2_xm1*vfpf) * vb ;
          f2_xm1 = vf2_xm1 [VSize -1] ;
        }

        V vf3_x00, vf7_x00, vf14_x00 ;
        T f14_xm1 ;
        {
          // Edge at yp1,z00.

          V vx00 = ptrSrc [idx_x00_yp1_z00 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_yp1_z00 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_yp1_z00 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_yp1_z00 + 0 + 3] * 1   ;

          V vxm1 = ptrSrc [idx_x00_yp1_z00 + W - NN + 0] + vfx ;
          V vym1 = ptrSrc [idx_x00_yp1_z00 + W - NN + 1] + vfy ;
          V vzm1 = ptrSrc [idx_x00_yp1_z00 + W - NN + 2] + vfz ;
          V  rm1 = ptrSrc [idx_x00_yp1_z00 + W - NN + 3] * 1  ;

          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
          V vfCommon_xm1 = vone - (vxm1*vxm1 + vym1*vym1 + vzm1*vzm1) * voph ;

          V  v3_x00 = - vy00 ;
          V  v7_x00 = - vx00 - vy00 ;
          V v14_x00 = + vx00 - vy00 ;
          V v14_xm1 = + vxm1 - vym1 ;

           vf3_x00 = r00 * (vfCommon_x00 + v3_x00*vthree + v3_x00*v3_x00*vfpf) * vb ;
           vf7_x00 = r00 * (vfCommon_x00 + v7_x00*vthree + v7_x00*v7_x00*vfpf) * vc ;
          vf14_x00 = r00 * (vfCommon_x00 + v14_x00*vthree + v14_x00*v14_x00*vfpf) * vc ;
          V
          vf14_xm1 = rm1 * (vfCommon_xm1 + v14_xm1*vthree + v14_xm1*v14_xm1*vfpf) * vc ;
          f14_xm1 = vf14_xm1 [VSize - 1] ;
        }
        V vf18_x00, vf22_x00, vf23_x00 ;
        T f22_xm1 ;
        {
          // Edge at ym1,zp1.

          V vx00 = ptrSrc [idx_x00_ym1_zp1 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_ym1_zp1 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_ym1_zp1 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_ym1_zp1 + 0 + 3] * 1   ;
        
          V vxm1 = ptrSrc [idx_x00_ym1_zp1 + W - NN + 0] + vfx ;
          V vym1 = ptrSrc [idx_x00_ym1_zp1 + W - NN + 1] + vfy ;
          V vzm1 = ptrSrc [idx_x00_ym1_zp1 + W - NN + 2] + vfz ;
          V  rm1 = ptrSrc [idx_x00_ym1_zp1 + W - NN + 3] * 1  ;
        
          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
          V vfCommon_xm1 = vone - (vxm1*vxm1 + vym1*vym1 + vzm1*vzm1) * voph ;

          V v18_x00 = + vy00 - vz00 ;
          V v22_x00 = + vx00 + vy00 - vz00 ;
          V v22_xm1 = + vxm1 + vym1 - vzm1 ;
          V v23_x00 = - vx00 + vy00 - vz00 ;

          vf18_x00 = r00 * (vfCommon_x00 + v18_x00*vthree + v18_x00*v18_x00*vfpf) * vc ;
          vf22_x00 = r00 * (vfCommon_x00 + v22_x00*vthree + v22_x00*v22_x00*vfpf) * vd ;
          V
          vf22_xm1 = rm1 * (vfCommon_xm1 + v22_xm1*vthree + v22_xm1*v22_xm1*vfpf) * vd ;
          f22_xm1 = vf22_xm1 [VSize - 1] ;
          vf23_x00 = r00 * (vfCommon_x00 + v23_x00*vthree + v23_x00*v23_x00*vfpf) * vd ;
        }
        V vf5_x00, vf9_x00, vf16_x00 ;
        T f16_xm1 ;
        {
          // Edge at y00,zp1.

          V vx00 = ptrSrc [idx_x00_y00_zp1 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_y00_zp1 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_y00_zp1 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_y00_zp1 + 0 + 3] * 1   ;
        
          V vxm1 = ptrSrc [idx_x00_y00_zp1 + W - NN + 0] + vfx ;
          V vym1 = ptrSrc [idx_x00_y00_zp1 + W - NN + 1] + vfy ;
          V vzm1 = ptrSrc [idx_x00_y00_zp1 + W - NN + 2] + vfz ;
          V  rm1 = ptrSrc [idx_x00_y00_zp1 + W - NN + 3] * 1  ;

          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
          V vfCommon_xm1 = vone - (vxm1*vxm1 + vym1*vym1 + vzm1*vzm1) * voph ;

          V  v5_x00 = - vz00 ;
          V  v9_x00 = - vx00 - vz00 ;
          V v16_x00 = + vx00 - vz00 ;
          V v16_xm1 = + vxm1 - vzm1 ;

          vf5_x00 = r00 * (vfCommon_x00 + v5_x00*vthree + v5_x00*v5_x00*vfpf) * vb ;
          vf9_x00 = r00 * (vfCommon_x00 + v9_x00*vthree + v9_x00*v9_x00*vfpf) * vc ;
          vf16_x00 = r00 * (vfCommon_x00 + v16_x00*vthree + v16_x00*v16_x00*vfpf) * vc ;
          V
          vf16_xm1 = rm1 * (vfCommon_xm1 + v16_xm1*vthree + v16_xm1*v16_xm1*vfpf) * vc ;
          f16_xm1 = vf16_xm1 [VSize - 1] ;
        }

        V vf11_x00, vf19_x00, vf25_x00 ;
        T f25_xm1 ;
        {
          // Edge at yp1,zp1.

          V vx00 = ptrSrc [idx_x00_yp1_zp1 + 0 + 0] + vfx ;
          V vy00 = ptrSrc [idx_x00_yp1_zp1 + 0 + 1] + vfy ;
          V vz00 = ptrSrc [idx_x00_yp1_zp1 + 0 + 2] + vfz ;
          V  r00 = ptrSrc [idx_x00_yp1_zp1 + 0 + 3] * 1   ;
        
          V vxm1 = ptrSrc [idx_x00_yp1_zp1 + W - NN + 0] + vfx ;
          V vym1 = ptrSrc [idx_x00_yp1_zp1 + W - NN + 1] + vfy ;
          V vzm1 = ptrSrc [idx_x00_yp1_zp1 + W - NN + 2] + vfz ;
          V  rm1 = ptrSrc [idx_x00_yp1_zp1 + W - NN + 3] * 1  ;
        
          V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
          V vfCommon_xm1 = vone - (vxm1*vxm1 + vym1*vym1 + vzm1*vzm1) * voph ;

          V v11_x00 = - vy00 - vz00 ;
          V v19_x00 = - vx00 - vy00 - vz00 ;
          V v25_x00 = + vx00 - vy00 - vz00 ;
          V v25_xm1 = + vxm1 - vym1 - vzm1 ;

          vf11_x00 = r00 * (vfCommon_x00 + v11_x00*vthree + v11_x00*v11_x00*vfpf) * vc ;
          vf19_x00 = r00 * (vfCommon_x00 + v19_x00*vthree + v19_x00*v19_x00*vfpf) * vd ;
          vf25_x00 = r00 * (vfCommon_x00 + v25_x00*vthree + v25_x00*v25_x00*vfpf) * vd ;
          V
          vf25_xm1 = rm1 * (vfCommon_xm1 + v25_xm1*vthree + v25_xm1*v25_xm1*vfpf) * vd ;
          f25_xm1 = vf25_xm1 [VSize - 1] ;
        }

        for (unsigned mIdx=0 ; mIdx < nMasksPerRow_ ; mIdx ++)
        {
          NodeMask mask = nodeMasks [mIdx] ;
          const IndexType xBegin = (IndexType)(mIdx) * N_NODES_PER_MASK ;
          // TODO: Replace std::min with separate code for the last iteration.
          const IndexType xEnd = std::min (xBegin + N_NODES_PER_MASK, W) ;

          for (IndexType xp = xBegin ; xp < xEnd ; xp += NN)    
          {
            //TODO: For all batches containing at least one fluid node, load
            //      FULL vectors to use them later in fluid-only batches.

            if (isBatchFluid (xp,y,z)) // Vector processing of full batch
            {
              V vvx = {0,0,0,0} ;
              V vvy = {0,0,0,0} ;
              V vvz = {0,0,0,0} ;
              V vr  = {0,0,0,0} ;

              IndexType xNext = (xp + NN) % W ; // First node after current batch.

              {
                // Edge at ym1,zm1.
              
                V vxp1 = ptrSrc [idx_x00_ym1_zm1 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_ym1_zm1 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_ym1_zm1 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_ym1_zm1 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;
              
                V v12_xp1 =        + vyp1 + vzp1 ;
                V v20_xp1 = + vxp1 + vyp1 + vzp1 ;
                V v26_xp1 = - vxp1 + vyp1 + vzp1 ;
              
                V vf12_xp1 = rp1 * (vfCommon_xp1 + v12_xp1*vthree + v12_xp1*v12_xp1*vfpf) * vc ;
                V vf20_xp1 = rp1 * (vfCommon_xp1 + v20_xp1*vthree + v20_xp1*v20_xp1*vfpf) * vd ;
                V vf26_xp1 = rp1 * (vfCommon_xp1 + v26_xp1*vthree + v26_xp1*v26_xp1*vfpf) * vd ;
              
                V vf12 = vf12_x00 ;
                vf12_x00 = vf12_xp1 ;

                V vf20 =  __builtin_shufflevector (vf20_x00,  vf20_x00, -1,0,1,2) ;
                vf20 [0] = f20_xm1 ;
                f20_xm1 = vf20_x00 [VSize -1] ;
                vf20_x00 = vf20_xp1 ;

                V vf26 = __builtin_shufflevector (vf26_x00, vf26_xp1, 1,2,3,4) ;
                vf26_x00 = vf26_xp1 ;

                vr  +=   vf12 + vf20 + vf26 ;
                vvx +=        + vf20 - vf26 ;
                vvy += + vf12 + vf20 + vf26 ;
                vvz += + vf12 + vf20 + vf26 ;
              }
              {
                // Edge at y00,zm1.
              
                V vxp1 = ptrSrc [idx_x00_y00_zm1 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_y00_zm1 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_y00_zm1 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_y00_zm1 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;
              
                V  v6_xp1 = + vzp1 ;
                V v10_xp1 = + vxp1 + vzp1 ;
                V v15_xp1 = - vxp1 + vzp1 ;
              
                V  vf6_xp1 = rp1 * (vfCommon_xp1 +  v6_xp1*vthree +  v6_xp1* v6_xp1*vfpf) * vb ;
                V vf10_xp1 = rp1 * (vfCommon_xp1 + v10_xp1*vthree + v10_xp1*v10_xp1*vfpf) * vc ;
                V vf15_xp1 = rp1 * (vfCommon_xp1 + v15_xp1*vthree + v15_xp1*v15_xp1*vfpf) * vc ;
              
                V vf6 = vf6_x00 ;
                vf6_x00 = vf6_xp1 ;

                V vf10 =  __builtin_shufflevector (vf10_x00,  vf10_x00, -1,0,1,2) ;
                vf10 [0] = f10_xm1 ;
                f10_xm1 = vf10_x00 [VSize -1] ;
                vf10_x00 = vf10_xp1 ;

                V vf15 = __builtin_shufflevector (vf15_x00, vf15_xp1, 1,2,3,4) ;
                vf15_x00 = vf15_xp1 ;

                vr  +=   vf6 + vf10 + vf15 ;
                vvx += + vf10 - vf15 ;
                vvz += + vf6 + vf10 + vf15 ;
              }
              {
                // Edge at yp1,zm1.

                V vxp1 = ptrSrc [idx_x00_yp1_zm1 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_yp1_zm1 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_yp1_zm1 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_yp1_zm1 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                V v17_xp1 =        - vyp1 + vzp1 ;
                V v21_xp1 = - vxp1 - vyp1 + vzp1 ;
                V v24_xp1 = + vxp1 - vyp1 + vzp1 ;

                V vf17_xp1 = rp1 * (vfCommon_xp1 + v17_xp1*vthree + v17_xp1*v17_xp1*vfpf) * vc ;
                V vf21_xp1 = rp1 * (vfCommon_xp1 + v21_xp1*vthree + v21_xp1*v21_xp1*vfpf) * vd ;
                V vf24_xp1 = rp1 * (vfCommon_xp1 + v24_xp1*vthree + v24_xp1*v24_xp1*vfpf) * vd ;

                V vf17 = vf17_x00 ;
                vf17_x00 = vf17_xp1 ;

                V vf21 = __builtin_shufflevector (vf21_x00, vf21_xp1, 1,2,3,4) ;
                vf21_x00 = vf21_xp1 ;

                V vf24 =  __builtin_shufflevector (vf24_x00,  vf24_x00, -1,0,1,2) ;
                vf24 [0] = f24_xm1 ;
                f24_xm1 = vf24_x00 [VSize -1] ;
                vf24_x00 = vf24_xp1 ;

                vr  +=   vf17 + vf21 + vf24 ;
                vvx +=        - vf21 + vf24 ;
                vvy += - vf17 - vf21 - vf24 ;
                vvz += + vf17 + vf21 + vf24 ;
              }
              {
                // Edge at ym1,z00.
              
                V vxp1 = ptrSrc [idx_x00_ym1_z00 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_ym1_z00 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_ym1_z00 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_ym1_z00 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;
              
                V  v4_xp1 = + vyp1 ;
                V  v8_xp1 = + vxp1 + vyp1 ;
                V v13_xp1 = - vxp1 + vyp1 ;
              
                V  vf4_xp1 = rp1 * (vfCommon_xp1 + v4_xp1*vthree + v4_xp1*v4_xp1*vfpf) * vb ;
                V  vf8_xp1 = rp1 * (vfCommon_xp1 + v8_xp1*vthree + v8_xp1*v8_xp1*vfpf) * vc ;
                V vf13_xp1 = rp1 * (vfCommon_xp1 + v13_xp1*vthree + v13_xp1*v13_xp1*vfpf) * vc ;
              
                V vf4 = vf4_x00 ;
                vf4_x00 = vf4_xp1 ;

                V vf8 =  __builtin_shufflevector (vf8_x00,  vf8_x00, -1,0,1,2) ;
                vf8 [0] = f8_xm1 ;
                f8_xm1 = vf8_x00 [VSize -1] ;
                vf8_x00 = vf8_xp1 ;

                V vf13 = __builtin_shufflevector (vf13_x00, vf13_xp1, 1,2,3,4) ;
                vf13_x00 = vf13_xp1 ;
              
                vr  +=   vf4 + vf8 + vf13 ;
                vvx += + vf8 - vf13 ;
                vvy += + vf4 + vf8 + vf13 ;
              }
              {
                // Edge at y00,z00.

                V vxp1 = ptrSrc [idx_x00_y00_z00 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_y00_z00 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_y00_z00 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_y00_z00 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                V v0_xp1 = {0,0,0,0} ;
                V v1_xp1 = - vxp1 ;
                V v2_xp1 = + vxp1 ;

                V  vf0_xp1 = rp1 * (vfCommon_xp1 + v0_xp1*vthree + v0_xp1*v0_xp1*vfpf) * va ;
                V  vf1_xp1 = rp1 * (vfCommon_xp1 + v1_xp1*vthree + v1_xp1*v1_xp1*vfpf) * vb ;
                V  vf2_xp1 = rp1 * (vfCommon_xp1 + v2_xp1*vthree + v2_xp1*v2_xp1*vfpf) * vb ;

                V vf0 = vf0_x00 ;
                vf0_x00 = vf0_xp1 ;

                V vf1 = __builtin_shufflevector (vf1_x00, vf1_xp1, 1,2,3,4) ;
                vf1_x00 = vf1_xp1 ;
              
                V vf2 =  __builtin_shufflevector (vf2_x00,  vf2_x00, -1,0,1,2) ;
                vf2 [0] = f2_xm1 ;
                f2_xm1 = vf2_x00 [VSize -1] ;
                vf2_x00 = vf2_xp1 ;

                vr  +=   vf0 + vf1 + vf2 ;
                vvx += - vf1 + vf2;
              }
              {
                // Edge at yp1,z00.

                V vxp1 = ptrSrc [idx_x00_yp1_z00 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_yp1_z00 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_yp1_z00 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_yp1_z00 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                V  v3_xp1 = - vyp1 ;
                V  v7_xp1 = - vxp1 - vyp1 ;
                V v14_xp1 = + vxp1 - vyp1 ;

                V  vf3_xp1 = rp1 * (vfCommon_xp1 + v3_xp1*vthree + v3_xp1*v3_xp1*vfpf) * vb ;
                V  vf7_xp1 = rp1 * (vfCommon_xp1 + v7_xp1*vthree + v7_xp1*v7_xp1*vfpf) * vc ;
                V vf14_xp1 = rp1 * (vfCommon_xp1 + v14_xp1*vthree + v14_xp1*v14_xp1*vfpf) * vc ;

                V vf3 = vf3_x00 ;
                vf3_x00 = vf3_xp1 ;
                V vf7 = __builtin_shufflevector (vf7_x00, vf7_xp1, 1,2,3,4) ;
                vf7_x00 = vf7_xp1 ;

                V vf14 =  __builtin_shufflevector (vf14_x00,  vf14_x00, -1,0,1,2) ;
                vf14 [0] = f14_xm1 ;
                f14_xm1 = vf14_x00 [VSize -1] ;
                vf14_x00 = vf14_xp1 ;

                vr  +=   vf3 + vf7 + vf14 ;
                vvx += - vf7 + vf14 ;
                vvy += - vf3 - vf7 - vf14 ;
              }
              {
                // Edge at ym1,zp1.

                V vxp1 = ptrSrc [idx_x00_ym1_zp1 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_ym1_zp1 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_ym1_zp1 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_ym1_zp1 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                V v18_xp1 = + vyp1 - vzp1 ;
                V v22_xp1 = + vxp1 + vyp1 - vzp1 ;
                V v23_xp1 = - vxp1 + vyp1 - vzp1 ;

                V vf18_xp1 = rp1 * (vfCommon_xp1 + v18_xp1*vthree + v18_xp1*v18_xp1*vfpf) * vc ;
                V vf22_xp1 = rp1 * (vfCommon_xp1 + v22_xp1*vthree + v22_xp1*v22_xp1*vfpf) * vd ;
                V vf23_xp1 = rp1 * (vfCommon_xp1 + v23_xp1*vthree + v23_xp1*v23_xp1*vfpf) * vd ;

                V vf18 = vf18_x00 ;
                vf18_x00 = vf18_xp1 ;

                V vf22 =  __builtin_shufflevector (vf22_x00,  vf22_x00, -1,0,1,2) ;
                vf22 [0] = f22_xm1 ;
                f22_xm1 = vf22_x00 [VSize -1] ;
                vf22_x00 = vf22_xp1 ;

                V vf23 = __builtin_shufflevector (vf23_x00, vf23_xp1, 1,2,3,4) ;
                vf23_x00 = vf23_xp1 ;

                vr  +=   vf18 + vf22 + vf23 ;
                vvx += + vf22 - vf23 ;
                vvy += + vf18 + vf22 + vf23 ;
                vvz += - vf18 - vf22 - vf23 ;
              }
              {
                // Edge at y00,zp1.

                V vxp1 = ptrSrc [idx_x00_y00_zp1 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_y00_zp1 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_y00_zp1 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_y00_zp1 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                V  v5_xp1 = - vzp1 ;
                V  v9_xp1 = - vxp1 - vzp1 ;
                V v16_xp1 = + vxp1 - vzp1 ;

                V vf5_xp1 = rp1 * (vfCommon_xp1 + v5_xp1*vthree + v5_xp1*v5_xp1*vfpf) * vb ;
                V vf9_xp1 = rp1 * (vfCommon_xp1 + v9_xp1*vthree + v9_xp1*v9_xp1*vfpf) * vc ;
                V vf16_xp1 = rp1 * (vfCommon_xp1 + v16_xp1*vthree + v16_xp1*v16_xp1*vfpf) * vc ;

                V vf5 = vf5_x00 ;
                vf5_x00 = vf5_xp1 ;
                V vf9 = __builtin_shufflevector (vf9_x00, vf9_xp1, 1,2,3,4) ;
                vf9_x00 = vf9_xp1 ;
                V vf16 =  __builtin_shufflevector (vf16_x00,  vf16_x00, -1,0,1,2) ;
                vf16 [0] = f16_xm1 ;
                f16_xm1 = vf16_x00 [VSize -1] ;
                vf16_x00 = vf16_xp1 ;

                vr  +=   vf5 + vf9 + vf16 ;
                vvx +=       - vf9 + vf16 ;
                vvz += - vf5 - vf9 - vf16 ;
              }
              {
                // Edge at yp1,zp1.

                V vxp1 = ptrSrc [idx_x00_yp1_zp1 + xNext + 0] + vfx ;
                V vyp1 = ptrSrc [idx_x00_yp1_zp1 + xNext + 1] + vfy ;
                V vzp1 = ptrSrc [idx_x00_yp1_zp1 + xNext + 2] + vfz ;
                V  rp1 = ptrSrc [idx_x00_yp1_zp1 + xNext + 3] * 1   ;
              
                V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                V v11_xp1 = - vyp1 - vzp1 ;
                V v19_xp1 = - vxp1 - vyp1 - vzp1 ;
                V v25_xp1 = + vxp1 - vyp1 - vzp1 ;

                V vf11_xp1 = rp1 * (vfCommon_xp1 + v11_xp1*vthree + v11_xp1*v11_xp1*vfpf) * vc ;
                V vf19_xp1 = rp1 * (vfCommon_xp1 + v19_xp1*vthree + v19_xp1*v19_xp1*vfpf) * vd ;
                V vf25_xp1 = rp1 * (vfCommon_xp1 + v25_xp1*vthree + v25_xp1*v25_xp1*vfpf) * vd ;

                V vf11 = vf11_x00 ;
                vf11_x00 = vf11_xp1 ;

                V vf19 = __builtin_shufflevector (vf19_x00, vf19_xp1, 1,2,3,4) ;
                vf19_x00 = vf19_xp1 ;

                V vf25 =  __builtin_shufflevector (vf25_x00,  vf25_x00, -1,0,1,2) ;
                vf25 [0] = f25_xm1 ;
                f25_xm1 = vf25_x00 [VSize -1] ;
                vf25_x00 = vf25_xp1 ;

                vr  +=   vf11 + vf19 + vf25 ;
                vvx +=        - vf19 + vf25 ;
                vvy += - vf11 - vf19 - vf25 ;
                vvz += - vf11 - vf19 - vf25 ;
              }

              V vinvR = vone / vr ;
              vvx *= vinvR ;
              vvy *= vinvR ;
              vvz *= vinvR ;

              {
                V * ptrDst = vPtr + offsetDst + idx_x00_y00_z00 + xp ;

                __builtin_nontemporal_store (vvx, ptrDst + 0) ;
                __builtin_nontemporal_store (vvy, ptrDst + 1) ;
                __builtin_nontemporal_store (vvz, ptrDst + 2) ;
                __builtin_nontemporal_store ( vr, ptrDst + 3) ;
              }
            }
            else if (isBatchMixed (xp,y,z)) // Process nodes separately.
            {
              IndexType xNext = (xp + NN) % W ; // First node after current batch.
              IndexType xPrev = (xp - NN + W) % W ;
              
              //WARNING: Must precompute values for the next fluid-only batch.
              //         We can be sure that fluid-only batch is surrounded by 
              //         either fluid-only or mixed batches.
              if (isBatchFluid (xNext, y,z))
              { 
                if (isBatchFluid (xPrev, y,z))
                {
                  // Can use values precomputed in the previous batch.
                }
                else
                {
                  // Must compute values for the previous batch, too.
                  // FIXME: use values computed in the for loop placed below.
                  {
                    // Edge at ym1,zm1.

                    V vx00 = ptrSrc [idx_x00_ym1_zm1 + xp + 0] + vfx ;
                    V vy00 = ptrSrc [idx_x00_ym1_zm1 + xp + 1] + vfy ;
                    V vz00 = ptrSrc [idx_x00_ym1_zm1 + xp + 2] + vfz ;
                    V  r00 = ptrSrc [idx_x00_ym1_zm1 + xp + 3] * 1   ;

                    V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
                    V v20_x00 = + vx00 + +vy00 + vz00 ;
                    vf20_x00 = r00 * (vfCommon_x00 + v20_x00*vthree + v20_x00*v20_x00*vfpf) * vd ;
                  }
                  {
                    // Edge at y0,zm1.

                    V vx00 = ptrSrc [idx_x00_y00_zm1 + xp + 0] + vfx ;
                    V vy00 = ptrSrc [idx_x00_y00_zm1 + xp + 1] + vfy ;
                    V vz00 = ptrSrc [idx_x00_y00_zm1 + xp + 2] + vfz ;
                    V  r00 = ptrSrc [idx_x00_y00_zm1 + xp + 3] * 1   ;

                    V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
                    V v10_x00 = + vx00 + vz00 ;
                    vf10_x00 = r00 * (vfCommon_x00 + v10_x00*vthree + v10_x00*v10_x00*vfpf) * vc ;
                  }
                  {
                    // Edge at yp1,zm1.

                    V vx00 = ptrSrc [idx_x00_yp1_zm1 + xp + 0] + vfx ;
                    V vy00 = ptrSrc [idx_x00_yp1_zm1 + xp + 1] + vfy ;
                    V vz00 = ptrSrc [idx_x00_yp1_zm1 + xp + 2] + vfz ;
                    V  r00 = ptrSrc [idx_x00_yp1_zm1 + xp + 3] * 1   ;

                    V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
                    V v24_x00 = + vx00 - vy00 + vz00 ;
                    vf24_x00 = r00 * (vfCommon_x00 + v24_x00*vthree + v24_x00*v24_x00*vfpf) * vd ;
                  }
                  {
                    // Edge at ym1,z0.

                    V vx00 = ptrSrc [idx_x00_ym1_z00 + xp + 0] + vfx ;
                    V vy00 = ptrSrc [idx_x00_ym1_z00 + xp + 1] + vfy ;
                    V vz00 = ptrSrc [idx_x00_ym1_z00 + xp + 2] + vfz ;
                    V  r00 = ptrSrc [idx_x00_ym1_z00 + xp + 3] * 1   ;

                    V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
                    V v8_x00 = + vx00 + vy00 ;
                    vf8_x00 = r00 * (vfCommon_x00 + v8_x00*vthree + v8_x00*v8_x00*vfpf) * vc ;
                  }
                  {
                    // Edge at y0,z0.

                    V vx00 = ptrSrc [idx_x00_y00_z00 + xp + 0] + vfx ;
                    V vy00 = ptrSrc [idx_x00_y00_z00 + xp + 1] + vfy ;
                    V vz00 = ptrSrc [idx_x00_y00_z00 + xp + 2] + vfz ;
                    V  r00 = ptrSrc [idx_x00_y00_z00 + xp + 3] * 1   ;

                    V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
                    V v2_x00 = + vx00 ;
                    vf2_x00 = r00 * (vfCommon_x00 + v2_x00*vthree + v2_x00*v2_x00*vfpf) * vb ;
                  }
                  {
                    // Edge at yp1,z0.

                    V vx00 = ptrSrc [idx_x00_yp1_z00 + xp + 0] + vfx ;
                    V vy00 = ptrSrc [idx_x00_yp1_z00 + xp + 1] + vfy ;
                    V vz00 = ptrSrc [idx_x00_yp1_z00 + xp + 2] + vfz ;
                    V  r00 = ptrSrc [idx_x00_yp1_z00 + xp + 3] * 1   ;

                    V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
                    V v14_x00 = + vx00 - vy00 ;
                    vf14_x00 = r00 * (vfCommon_x00 + v14_x00*vthree + v14_x00*v14_x00*vfpf) * vc ;
                  }
                  {
                    // Edge at ym1,zp1.

                    V vx00 = ptrSrc [idx_x00_ym1_zp1 + xp + 0] + vfx ;
                    V vy00 = ptrSrc [idx_x00_ym1_zp1 + xp + 1] + vfy ;
                    V vz00 = ptrSrc [idx_x00_ym1_zp1 + xp + 2] + vfz ;
                    V  r00 = ptrSrc [idx_x00_ym1_zp1 + xp + 3] * 1   ;

                    V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
                    V v22_x00 = + vx00 + vy00 - vz00 ;
                    vf22_x00 = r00 * (vfCommon_x00 + v22_x00*vthree + v22_x00*v22_x00*vfpf) * vd ;
                  }
                  {
                    // Edge at y0,zp1.

                    V vx00 = ptrSrc [idx_x00_y00_zp1 + xp + 0] + vfx ;
                    V vy00 = ptrSrc [idx_x00_y00_zp1 + xp + 1] + vfy ;
                    V vz00 = ptrSrc [idx_x00_y00_zp1 + xp + 2] + vfz ;
                    V  r00 = ptrSrc [idx_x00_y00_zp1 + xp + 3] * 1   ;

                    V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
                    V v16_x00 = + vx00 - vz00 ;
                    vf16_x00 = r00 * (vfCommon_x00 + v16_x00*vthree + v16_x00*v16_x00*vfpf) * vc ;
                  }
                  {
                    // Edge at yp1,zp1.

                    V vx00 = ptrSrc [idx_x00_yp1_zp1 + xp + 0] + vfx ;
                    V vy00 = ptrSrc [idx_x00_yp1_zp1 + xp + 1] + vfy ;
                    V vz00 = ptrSrc [idx_x00_yp1_zp1 + xp + 2] + vfz ;
                    V  r00 = ptrSrc [idx_x00_yp1_zp1 + xp + 3] * 1   ;

                    V vfCommon_x00 = vone - (vx00*vx00 + vy00*vy00 + vz00*vz00) * voph ;
                    V v25_x00 = + vx00 - vy00 - vz00 ;
                    vf25_x00 = r00 * (vfCommon_x00 + v25_x00*vthree + v25_x00*v25_x00*vfpf) * vd ;
                  }
                }

                {
                  // Edge at ym1,zm1.

                  V vxp1 = ptrSrc [idx_x00_ym1_zm1 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_ym1_zm1 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_ym1_zm1 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_ym1_zm1 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V v12_xp1 =        + vyp1 + vzp1 ;
                  V v20_xp1 = + vxp1 + vyp1 + vzp1 ;
                  V v26_xp1 = - vxp1 + vyp1 + vzp1 ;

                  V vf12_xp1 = rp1 * (vfCommon_xp1 + v12_xp1*vthree + v12_xp1*v12_xp1*vfpf) * vc ;
                  vf12_x00 = vf12_xp1 ;
                  V vf20_xp1 = rp1 * (vfCommon_xp1 + v20_xp1*vthree + v20_xp1*v20_xp1*vfpf) * vd ;
                  f20_xm1 = vf20_x00 [VSize - 1] ;
                  vf20_x00 = vf20_xp1 ;
                  V vf26_xp1 = rp1 * (vfCommon_xp1 + v26_xp1*vthree + v26_xp1*v26_xp1*vfpf) * vd ;
                  vf26_x00 = vf26_xp1 ;
                }
                {
                  // Edge at y00,zm1.

                  V vxp1 = ptrSrc [idx_x00_y00_zm1 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_y00_zm1 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_y00_zm1 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_y00_zm1 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V  v6_xp1 = + vzp1 ;
                  V v10_xp1 = + vxp1 + vzp1 ;
                  V v15_xp1 = - vxp1 + vzp1 ;

                  V vf6_xp1 = rp1 * (vfCommon_xp1 + v6_xp1*vthree + v6_xp1*v6_xp1*vfpf) * vb ;
                  vf6_x00 = vf6_xp1 ;
                  V vf10_xp1 = rp1 * (vfCommon_xp1 + v10_xp1*vthree + v10_xp1*v10_xp1*vfpf) * vc ;
                  f10_xm1 = vf10_x00 [VSize - 1] ;
                  vf10_x00 = vf10_xp1 ;
                  V vf15_xp1 = rp1 * (vfCommon_xp1 + v15_xp1*vthree + v15_xp1*v15_xp1*vfpf) * vc ;
                  vf15_x00 = vf15_xp1 ;
                }
                {
                  // Edge at yp1,zm1.

                  V vxp1 = ptrSrc [idx_x00_yp1_zm1 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_yp1_zm1 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_yp1_zm1 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_yp1_zm1 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V v17_xp1 =        - vyp1 + vzp1 ;
                  V v21_xp1 = - vxp1 - vyp1 + vzp1 ;
                  V v24_xp1 = + vxp1 - vyp1 + vzp1 ;

                  V vf17_xp1 = rp1 * (vfCommon_xp1 + v17_xp1*vthree + v17_xp1*v17_xp1*vfpf) * vc ;
                  V vf21_xp1 = rp1 * (vfCommon_xp1 + v21_xp1*vthree + v21_xp1*v21_xp1*vfpf) * vd ;
                  V vf24_xp1 = rp1 * (vfCommon_xp1 + v24_xp1*vthree + v24_xp1*v24_xp1*vfpf) * vd ;

                  vf17_x00 = vf17_xp1 ;
                  vf21_x00 = vf21_xp1 ;
                  f24_xm1 = vf24_x00 [VSize -1] ;
                  vf24_x00 = vf24_xp1 ;
                }
                {
                  // Edge at ym1,z00.

                  V vxp1 = ptrSrc [idx_x00_ym1_z00 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_ym1_z00 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_ym1_z00 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_ym1_z00 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V  v4_xp1 = + vyp1 ;
                  V  v8_xp1 = + vxp1 + vyp1 ;
                  V v13_xp1 = - vxp1 + vyp1 ;

                  V vf4_xp1 = rp1 * (vfCommon_xp1 + v4_xp1*vthree + v4_xp1*v4_xp1*vfpf) * vb ;
                  vf4_x00 = vf4_xp1 ;
                  V vf8_xp1 = rp1 * (vfCommon_xp1 + v8_xp1*vthree + v8_xp1*v8_xp1*vfpf) * vc ;
                  f8_xm1 = vf8_x00 [VSize -1] ;
                  vf8_x00 = vf8_xp1 ;
                  V vf13_xp1 = rp1 * (vfCommon_xp1 + v13_xp1*vthree + v13_xp1*v13_xp1*vfpf) * vc ;
                  vf13_x00 = vf13_xp1 ;
                }
                {
                  // Edge at y0,z0.

                  V vxp1 = ptrSrc [idx_x00_y00_z00 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_y00_z00 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_y00_z00 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_y00_z00 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V v0_xp1 = {0,0,0,0} ;
                  V v1_xp1 = - vxp1 ;
                  V v2_xp1 = + vxp1 ;

                  V  vf0_xp1 = rp1 * (vfCommon_xp1 + v0_xp1*vthree + v0_xp1*v0_xp1*vfpf) * va ;
                  vf0_x00 = vf0_xp1 ;

                  V  vf1_xp1 = rp1 * (vfCommon_xp1 + v1_xp1*vthree + v1_xp1*v1_xp1*vfpf) * vb ;
                  vf1_x00 = vf1_xp1 ;

                  V  vf2_xp1 = rp1 * (vfCommon_xp1 + v2_xp1*vthree + v2_xp1*v2_xp1*vfpf) * vb ;
                  f2_xm1 = vf2_x00 [VSize -1] ;
                  vf2_x00 = vf2_xp1 ;
                }
                {
                  // Edge at yp1,z00.

                  V vxp1 = ptrSrc [idx_x00_yp1_z00 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_yp1_z00 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_yp1_z00 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_yp1_z00 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V  v3_xp1 = - vyp1 ;
                  V  v7_xp1 = - vxp1 - vyp1 ;
                  V v14_xp1 = + vxp1 - vyp1 ;

                  V vf3_xp1 = rp1 * (vfCommon_xp1 + v3_xp1*vthree + v3_xp1*v3_xp1*vfpf) * vb ;
                  V vf7_xp1 = rp1 * (vfCommon_xp1 + v7_xp1*vthree + v7_xp1*v7_xp1*vfpf) * vc ;
                  V vf14_xp1 = rp1 * (vfCommon_xp1 + v14_xp1*vthree + v14_xp1*v14_xp1*vfpf) * vc ;
                  vf3_x00 = vf3_xp1 ;
                  vf7_x00 = vf7_xp1 ;
                  f14_xm1 = vf14_x00 [VSize - 1] ;
                  vf14_x00 = vf14_xp1 ;
                }
                {
                  // Edge at ym1,zp1.

                  V vxp1 = ptrSrc [idx_x00_ym1_zp1 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_ym1_zp1 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_ym1_zp1 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_ym1_zp1 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V v18_xp1 = + vyp1 - vzp1 ;
                  V v22_xp1 = + vxp1 + vyp1 - vzp1 ;
                  V v23_xp1 = - vxp1 + vyp1 - vzp1 ;

                  V vf18_xp1 = rp1 * (vfCommon_xp1 + v18_xp1*vthree + v18_xp1*v18_xp1*vfpf) * vc ;
                  V vf22_xp1 = rp1 * (vfCommon_xp1 + v22_xp1*vthree + v22_xp1*v22_xp1*vfpf) * vd ;
                  V vf23_xp1 = rp1 * (vfCommon_xp1 + v23_xp1*vthree + v23_xp1*v23_xp1*vfpf) * vd ;

                  vf18_x00 = vf18_xp1 ;
                  f22_xm1 = vf22_x00 [VSize - 1] ;
                  vf22_x00 = vf22_xp1 ;
                  vf23_x00 = vf23_xp1 ;
                }
                {
                  // Edge at y00,zp1.

                  V vxp1 = ptrSrc [idx_x00_y00_zp1 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_y00_zp1 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_y00_zp1 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_y00_zp1 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V  v5_xp1 = - vzp1 ;
                  V  v9_xp1 = - vxp1 - vzp1 ;
                  V v16_xp1 = + vxp1 - vzp1 ;

                  V vf5_xp1 = rp1 * (vfCommon_xp1 + v5_xp1*vthree + v5_xp1*v5_xp1*vfpf) * vb ;
                  V vf9_xp1 = rp1 * (vfCommon_xp1 + v9_xp1*vthree + v9_xp1*v9_xp1*vfpf) * vc ;
                  V vf16_xp1 = rp1 * (vfCommon_xp1 + v16_xp1*vthree + v16_xp1*v16_xp1*vfpf) * vc ;

                  vf5_x00 = vf5_xp1 ;
                  vf9_x00 = vf9_xp1 ;
                  f16_xm1 = vf16_x00 [VSize - 1] ;
                  vf16_x00 = vf16_xp1 ;
                }
                {
                  // Edge at yp1,zp1.

                  V vxp1 = ptrSrc [idx_x00_yp1_zp1 + xNext + 0] + vfx ;
                  V vyp1 = ptrSrc [idx_x00_yp1_zp1 + xNext + 1] + vfy ;
                  V vzp1 = ptrSrc [idx_x00_yp1_zp1 + xNext + 2] + vfz ;
                  V  rp1 = ptrSrc [idx_x00_yp1_zp1 + xNext + 3] * 1   ;

                  V vfCommon_xp1 = vone - (vxp1*vxp1 + vyp1*vyp1 + vzp1*vzp1) * voph ;

                  V v11_xp1 = - vyp1 - vzp1 ;
                  V v19_xp1 = - vxp1 - vyp1 - vzp1 ;
                  V v25_xp1 = + vxp1 - vyp1 - vzp1 ;

                  V vf11_xp1 = rp1 * (vfCommon_xp1 + v11_xp1*vthree + v11_xp1*v11_xp1*vfpf) * vc ;
                  V vf19_xp1 = rp1 * (vfCommon_xp1 + v19_xp1*vthree + v19_xp1*v19_xp1*vfpf) * vd ;
                  V vf25_xp1 = rp1 * (vfCommon_xp1 + v25_xp1*vthree + v25_xp1*v25_xp1*vfpf) * vd ;

                  vf11_x00 = vf11_xp1 ;
                  vf19_x00 = vf19_xp1 ;
                  f25_xm1 = vf25_x00 [VSize -1] ;
                  vf25_x00 = vf25_xp1 ;
                }
              }

              computeBoundaryBatch 
              (
                xp,y,z, W,H,D, idx_x00_y00_z00, 
                ptrSrc, vPtr + offsetDst + idx_x00_y00_z00 + xp,
                direction, vfx[0],vfy[0],vfz[0]
              ) ;
            }
          }
        }
      }
  }
}


