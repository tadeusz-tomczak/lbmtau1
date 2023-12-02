#!/bin/bash

measure_performance()
{
  local nThreads=$1
  local geometrySize=$2
  local nObstacles=$3

  echo -e "\n\n################  Begin of measurements ################\n\n"

  date
  local startTime=`date +%s`
  echo ""

  cmd="OMP_DISPLAY_AFFINITY=TRUE OMP_PROC_BIND=true OMP_PLACES={0}:$nThreads:1 OMP_NUM_THREADS=$nThreads ./main 1000 3 $geometrySize $nObstacles"

  echo $cmd
  eval $cmd
  
  echo ""
  date

  echo -e "\n\n################  End of measurements ################\n\n"
}

measure_performance_and_log()
{
  local nThreads=$1
  local geometrySize=$2
  local nObstacles=$3

  local logfile=log_${nThreads}_threads_${geometrySize}_geometrySize_${nObstacles}_nObstacles.txt

  echo $logfile

  measure_performance $nThreads $geometrySize $nObstacles > $logfile 2>&1
}

compile()
{
  local Q=$1
  local dataType=$2

  date
  make clean
  EXTRA_OPTIONS="-DDEF_DATA_TYPE=$dataType -DDEF_Q=$Q" make -j -B
  date
}

compile_and_log()
{
  local Q=$1
  local dataType=$2

  echo $Q $dataType

  local compilation_logfile="log_compilation.txt"
  echo $compilation_logfile

  compile $Q $dataType > $compilation_logfile 2>&1
}

compile_and_measure_vs_numThreads()
{
  local Q=$1
  local dataType=$2
  local geometrySize=$3
  local nObstacles=$4

  compile_and_log $Q $dataType

  for nThreads in {1..32}
  do
    echo "nThreads = " $nThreads
    measure_performance_and_log $nThreads $geometrySize $nObstacles
  done

  local packedFile=q_${Q}_${dataType}_geometrySize_${geometrySize}_nObstacles_${nObstacles}.tar.bz2

  tar -cjf $packedFile log*.txt
  rm -rf log*.txt
}

compile_and_measure_vs_gridSize()
{
  local Q=$1
  local dataType=$2
  local nThreads=$3
  local nObstacles=$4

  compile_and_log $Q $dataType

  for geometrySize in {8..256..8}
  do
    echo "geometrySize = " $geometrySize
    measure_performance_and_log $nThreads $geometrySize $nObstacles
  done

  local packedFile=q_${Q}_${dataType}_nThreads_${nThreads}_nObstacles_${nObstacles}.tar.bz2

  tar -cjf $packedFile log*.txt
  rm -rf log*.txt
}

compile_and_measure_vs_porosity()
{
  local Q=$1
  local dataType=$2
  local nThreads=$3
  local geometrySize=$4

  compile_and_log $Q $dataType

  for nObstacles in 0 2 5 10 15 20 30 {40..300..20}
  do
    echo "nObstacles = " $nObstacles
    measure_performance_and_log $nThreads $geometrySize $nObstacles
  done

  local packedFile=q_${Q}_${dataType}_nThreads_${nThreads}_geometrySize_${geometrySize}_vs_nObstacles.tar.bz2

  tar -cjf $packedFile log*.txt
  rm -rf log*.txt
}


compile_and_measure_vs_numThreads 19 float  256 0
compile_and_measure_vs_numThreads 19 double 256 0
compile_and_measure_vs_numThreads 27 float  256 0
compile_and_measure_vs_numThreads 27 double 256 0

compile_and_measure_vs_gridSize 19 double  5 0
compile_and_measure_vs_gridSize 19 float   6 0
compile_and_measure_vs_gridSize 27 double  5 0
compile_and_measure_vs_gridSize 27 float  11 0

compile_and_measure_vs_porosity 19 float  16 256
compile_and_measure_vs_porosity 19 double 16 256
compile_and_measure_vs_porosity 27 float  16 256
compile_and_measure_vs_porosity 27 double 16 256
