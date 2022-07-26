#ifndef RAPIDNJ_WRAPPER_H
#define RAPIDNJ_WRAPPER_H

#include "stdinclude.h"
#include "dataLoaderPointer.hpp"
#include "rapidNJ.h"
#include "rapidNJ.h"
#include "rapidNJDisk.h"
#include "rapidNJMem.hpp"
#include "simpleNJ.h"
#include "JCdistance.hpp"
#include "KimuraDistance.hpp"
#include "ProgressBar.hpp"
#include <iomanip>
#include <sstream>

#define BUILDING_DLL 1
#define PTW32_STATIC_LIB 1

#if defined _WIN32 || defined __CYGWIN__ || defined __MINGW32__
#ifdef BUILDING_DLL
#ifdef __GNUC__
#define DLL_PUBLIC __attribute__ ((dllexport))
#else
#define DLL_PUBLIC __declspec(dllexport) // Note: actually gcc seems to also supports this syntax.
#endif
#else
#ifdef __GNUC__
#define DLL_PUBLIC __attribute__ ((dllimport))
#else
#define DLL_PUBLIC __declspec(dllimport) // Note: actually gcc seems to also supports this syntax.
#endif
#endif
#define DLL_LOCAL
#else
#if __GNUC__ >= 4
#define DLL_PUBLIC __attribute__ ((visibility ("default")))
#define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
#else
#define DLL_PUBLIC
#define DLL_LOCAL
#endif
#endif


#endif

typedef void (*return_callback)(size_t, const char*);

extern "C"
{
	DLL_PUBLIC void BuildTreeFromAlignment(int maxMemory, int distance, int numCores, int bootstrapReplicates, int inputType, bool allowNegativeBranches, int inputSequenceCount, int inputSequenceLength, int* inputSequenceNamesLengths, char** inputSequenceNames, char** inputSequenceData, progress_callback callback, return_callback returnCallback, bool verbose);
	DLL_PUBLIC void BuildDistanceMatrixFromAlignment(int maxMemory, int distance, int numCores, int inputType, int inputSequenceCount, int inputSequenceLength, int* inputSequenceNamesLengths, char** inputSequenceNames, char** inputSequenceData, distType** outputMatrix, bool verbose);
	DLL_PUBLIC void BuildTreeFromDistanceMatrix(int maxMemory, int numCores, bool allowNegativeBranches, int inputSequenceCount, int* inputSequenceNamesLengths, char** inputSequenceNames, bool halfMatrix, distType** distMatrix, progress_callback callback, return_callback returnCallback, bool verbose);
}