#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <cstring>
#include <fstream>
using namespace std;

#include "LBMTau1Vector3D.hpp"


template <class Domain>
void exportVtk (Domain & d, unsigned stepIdx)
{
/*
1. # vtk DataFile Version 2.0
2. Comments
3. ASCII
4. DATASET STRUCTURED_POINTS
5. DIMENSIONS 4 2 1
6. ORIGIN 0 0 0
7. SPACING 1 1 1
8. POINT_DATA 8
10. VECTORS PolePredkosci double
11. 1.0 0.0 0.0*/

	unsigned width  = d.w() ;
	unsigned height = d.h() ;
	unsigned depth  = d.d() ;

	stringstream fName ;
	fName << "velocity_" << setw(10) << setfill('0') << stepIdx << ".vtk" ;
	ofstream file(fName.str());
	cout << "Saving " << fName.str() << "..." << flush ;
	
	file << "# vtk DataFile Version 2.0\nLBM, Symulacje Komputerowe w Fizyce 2, Maciej Matyka 2019\n";
	file << "ASCII\nDATASET STRUCTURED_POINTS\n";
	file << "DIMENSIONS " << width << " " << height << " " << depth << " \n";
	file << "ORIGIN 0 0 0\nSPACING 1 1 1\n";
	file << "POINT_DATA " << width*height*depth << "\n";
	file << "VECTORS Velocity double\n";
  file <<  std::setprecision(10);

	for (unsigned z = 0 ; z < depth ; z++)
		for (unsigned y = 0 ; y < height ; y++)
			for (unsigned x = 0 ; x < width ; x++)
			{
				if (d.t (x,y,z).isSolid())
				{
        	file << "0 0 0" << endl;
				}
				else
				{
					file << d.get_vx (x,y,z) << " "
							 << d.get_vy (x,y,z) << " " 
							 << d.get_vz (x,y,z) << endl ;
				}
			}

	file.close();
	cout << "  OK\n" ;

	fName.str("") ;
	fName << "density_" << setw(10) << setfill('0') << stepIdx << ".vtk" ;
	ofstream file2(fName.str());
	cout << "Saving " << fName.str() << "..." << flush ;
	
	file2 << "# vtk DataFile Version 2.0\nLBM, Symulacje Komputerowe w Fizyce 2, Maciej Matyka 2019\n";
	file2 << "ASCII\nDATASET STRUCTURED_POINTS\n";
	file << "DIMENSIONS " << width << " " << height << " " << depth << " \n";
	file << "ORIGIN 0 0 0\nSPACING 1 1 1\n";
	file << "POINT_DATA " << width*height*depth << "\n";
	file2 << "SCALARS Density double\n";
	file2 << "LOOKUP_TABLE default\n";

	for (unsigned z = 0 ; z < depth ; z++)
		for (unsigned y = 0 ; y < height ; y++)
			for (unsigned x = 0 ; x < width ; x++)
			{
				if (d.t (x,y,z).isSolid())
				{
        	file2 << 0 << endl;
				}
				else
				{
					file2 << ((d.get_r (x,y,z) - 1.0) * 1e4) << endl ;
				}
			}

	file2.close();
	cout << "  OK\n" ;
}

template <class Domain>
void buildEmptyChannel (Domain & d)
{
	d.resetAllNodesToFluid() ;

	for (unsigned z = 0 ; z < d.d() ; z++)
		for (unsigned x = 0 ; x < d.w() ; x++)
		{
			d.setSolid (x, 0    , z) ;
			d.setSolid (x, d.h()-1, z) ;
		}
	for (unsigned y = 0 ; y < d.h() ; y++)
		for (unsigned x = 0 ; x < d.w() ; x++)
		{
			d.setSolid (x, y, 0) ;
			d.setSolid (x, y, d.d()-1) ;
		}
}

template <class Domain>
void placeSphericalObstacle 
(
	Domain & d, 
	unsigned x0, unsigned y0, unsigned z0, // sphere center
	double radius
)
{
	const auto W = d.w() ;
	const auto H = d.h() ;
	const auto D = d.d() ;

	const unsigned long long sqrR = round (radius * radius) ;

	unsigned minX, maxX, minY, maxY, minZ, maxZ ;
	const unsigned rp2 = ceil (radius) + 2 ;

	maxX = std::min (W, x0 + rp2) ;
	maxY = std::min (H, y0 + rp2) ;
	maxZ = std::min (D, z0 + rp2) ;

	if (x0 > rp2) minX = x0 - rp2 ; else minX = 0 ;
	if (y0 > rp2) minY = y0 - rp2 ; else minY = 0 ;
	if (z0 > rp2) minZ = z0 - rp2 ; else minZ = 0 ;
 
	for (unsigned x=minX ; x < maxX ; x++)   
		for (unsigned y=minY ; y < maxY ; y++)   	   	
			for (unsigned z=minZ ; z < maxZ ; z++)   	   	
			{
				unsigned long long r2 = 
				 ((x-x0))*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0) ;
				if (r2 < sqrR) 
				{
					d.setSolid (x,y,z) ;
				}
			}		  
}

template <class Domain>
void buildChannelWithObstacles (Domain & d, int nObstacles = 20)
{
	buildEmptyChannel (d) ;

	const auto   W = d.w() ;
	const auto   H = d.h() ;
	const auto   D = d.d() ;
	const double r = H / 5.0 ;
 
	srand(1000);
  for (int i=0; i < nObstacles ; i++)
	{
		unsigned x0 = round (((float)W)*(rand()/(float)RAND_MAX)) ;
		unsigned y0 = round (((float)H)*(rand()/(float)RAND_MAX)) ;
		unsigned z0 = round (((float)D)*(rand()/(float)RAND_MAX)) ;

		placeSphericalObstacle (d, x0,y0,z0, r) ;
	}
}

template <class Domain>
void buildChannelWithMixedChunks (Domain & d)
{
	buildEmptyChannel (d) ;

	const auto   W = d.w() ;
	const auto   H = d.h() ;
	const auto   D = d.d() ;

	const size_t NN = d.N_NODES_PER_BATCH ;

	for (IndexType z=0 ; z < D ; z++)
		for (IndexType y=0 ; y < H ; y++)
			for (IndexType xb=0 ; xb < W ; xb += NN)
			{
				IndexType xo = 0 ;
				if (0 != z%2) xo = NN-2 ;
				if (0 != y%2) xo = NN-2 ;

				d.setSolid (xb + xo, y,z) ;
			}
}



int main (int argc, char ** argv)
{
	unsigned N_ITER = 1000 ;
	if (argc > 1) 
	{
		N_ITER = atoi (argv[1]) ;
	}
	cout << "N_ITER = " << N_ITER << "\n" ;

	unsigned N_MEASURES = 10 ;
	if (argc > 2) 
	{
		N_MEASURES = atoi (argv[2]) ;
	}
	cout << "N_MEASURES = " << N_MEASURES << "\n" ;

	size_t N = 40 ;
	if (argc > 3)
	{
		N = atol (argv[3]) ;
	}
	size_t width  = N * 2 ;
	size_t height = N ;
	size_t depth  = N ;
	cout << "width x height x depth = " 
			 << width << " x " << height << " x " 
			 << depth << "\n" ;

	int nObstacles = 20 ;
	if (argc > 4)
	{
		nObstacles = atol (argv[4]) ;
	}
	cout << "nObstacles = " << nObstacles << "\n" ;
	
	using Clock = std::chrono::high_resolution_clock ;

	#ifdef DEF_DATA_TYPE
		using T = DEF_DATA_TYPE ;
	#else
		using T = double ;
	#endif
	#ifdef DEF_Q
		constexpr int Q = DEF_Q ;
	#else
		constexpr int Q = 19 ;
	#endif
	LBMTau1Vector3D <LatticeArrangement<3,Q,T>,T> d (width,height,depth) ;

	using L = decltype(d)::LatticeArrangement ;
	cout << "sizeof (DataType)   = " << sizeof (T) << "\n" ;
	cout << "latticeArrangement  = D" << L::D << "Q" << L::Q << "\n" ;

	cout << "Building geometry...\n" ;
	buildChannelWithObstacles (d, nObstacles) ;
	//buildChannelWithMixedChunks (d) ;

	d.setF (-0.000065, 0, 0) ;

	// Slow, but useful in case of problems with geometry.
	// Uncomment if needed.
	//d.checkNodeMasks() ;

	const auto s = d.computeStatistics() ;
	double porosity = (double)(s.nComputationalNodes)/s.nNodes ;
	double pFlCh = (double)(s.nFluidChunks) / s.nComputationalChunks ;
	double pMxCh = (double)(s.nMixedChunks) / s.nComputationalChunks ;
	double pCpCh = (double)(s.nComputationalChunks) / s.nChunks ;
	double pSlCh = (double)(s.nSolidChunks) / s.nChunks ;
	double chunkPorosity = (double)(s.nComputationalNodes) / (s.nComputationalChunks * d.N_NODES_PER_BATCH) ;
	cout << "Geometry statistics:\n"
			 << "nNodes               = " << s.nNodes               << "\n"
			 << "nFluidNodes          = " << s.nFluidNodes          << "\n"
			 << "nBoundaryNodes       = " << s.nBoundaryNodes       << "\n"
			 << "nSolidNodes          = " << s.nSolidNodes          << "\n"
			 << "nComputationalNodes  = " << s.nComputationalNodes  << "\n"
			 << "porosity             = " << porosity               << "\n"
			 << "nChunks              = " << s.nChunks              << "\n"
			 << "nFluidChunks         = " << s.nFluidChunks 
			 << "  (" << pFlCh << " of computational chunks)\n"
			 << "nMixedChunks         = " << s.nMixedChunks        
			 << "  (" << pMxCh << " of computational chunks)\n"
			 << "nSolidChunks         = " << s.nSolidChunks         
			 << "  (" << pSlCh << " of all chunks)\n"
			 << "nComputationalChunks = " << s.nComputationalChunks
			 << "  (" << pCpCh << " of all chunks)\n"
			 << "nFluidChunksAfterMixed = " << s.nFluidChunksAfterMixed << "\n"
			 << "nMixedChunksBetweenFluid = " << s.nMixedChunksBetweenFluid << "\n"
			 << "chunk porosity       = " << chunkPorosity          << "\n"
			 ;

	cout << "Running simulation...\n" ;

	for (unsigned im=0 ; im < N_MEASURES ; im++)
	{
		unsigned timeStepCounter = 0 ;
		size_t nNodes = s.nComputationalNodes ;

		auto timeBegin = Clock::now() ;
		for (unsigned sIdx = 0 ; sIdx < N_ITER ; sIdx++)
		{
			d.iter() ;
			timeStepCounter += 1 ;
		}
		auto timeEnd = Clock::now() ;

		auto totalDuration = timeEnd - timeBegin ;
		double totalDuration_ns = 
			chrono::duration_cast<chrono::nanoseconds>(totalDuration).count() ;

		const auto nTimeSteps = timeStepCounter ;
		uint64_t LU = nNodes * nTimeSteps ;
		double MLUPS = (LU / 1e6) / (totalDuration_ns / 1e9) ;
		double GBs   = MLUPS * 2 * sizeof (T) * 4 / 1e3 ;

		cout << "nTimeSteps = " << nTimeSteps << ", "
			<< "duration = " << totalDuration_ns / 1e6 << " ms, " 
			<< totalDuration_ns / 1e6 / (nTimeSteps) << " ms per kernel, "
			<< totalDuration_ns / LU << " ns per lattice node "
			<< MLUPS << " MLUPS, " << GBs << " GB/s\n" ;
	}

	exportVtk (d, N_MEASURES * N_ITER) ;

	return 0 ;
}
