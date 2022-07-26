#include "rapidNJWrapper.h"

using namespace std;

bool distanceMatrixInput = true;
bool distanceMatrixFromPointer = false;
int matrixSize = -1;
int numCores = 1;

namespace options {
	bool verbose;
	string fileName;
	int memSize;
	int cores;
	string cacheDir;
	string percentageMemoryUsage;
	string distMethod;
	string inputFormat;
	string outputFormat;
	bool fastdist;
	int replicates;
	string inputtype;
	bool rapidNJ;
	bool simpleNJ;
	bool gpu;
	bool negative_branches;
	string outputFile;
};

class distMatrixData {
public:
	distType** matrix;
	vector<string>* sequenceNames;
	diskMatrix* dm;
};

// Helper functions

void configureNumberOfCores() {
	// Configure number of cores to use

	if (options::cores > 0) {
		numCores = options::cores;
	}
	if (numCores < 1) {
		numCores = 1;
	}
	if (options::verbose) {
		cerr << "Using " << numCores << " core(s) for distance estimation" << endl;
	}
}

double getMemSize() {
	return options::memSize * 1024.0 * 1024.0;
}

void printDistanceMatrix(ostream& out, distMatrixData* data) {
	out << "\t" << matrixSize << endl;
	for (int i = 0; i < matrixSize; i++) {
		out << (*data->sequenceNames)[i] << "\t";
		for (int j = 0; j < matrixSize; j++) {
			out << setprecision(6) << fixed << data->matrix[i][j] << " ";
		}
		out << endl;
	}
}

void printDistanceMatrixDisk(ostream& out, distMatrixData* data) {
	out << "\t" << matrixSize << endl;
	distType* row = new distType[matrixSize];
	for (int i = 0; i < matrixSize; i++) {
		out << (*data->sequenceNames)[i] << "\t";
		for (int j = 0; j < matrixSize; j++) {
			data->dm->readArray(row, i, matrixSize);
			out << setprecision(6) << fixed << row[j] << " ";
		}
		out << endl;
	}
	delete[] row;
}

distMatrixData* computeDistanceMatrix(bool useDiskMatrix, ostream& out, bool printMatrix, dataloader* dl) {
	/*if(options::gpu){
	computeDistanceMatrixGPU();
	return;
	}*/

	if (options::fastdist && options::verbose) {
		cerr << "Fastdist is enabled" << endl;
	}

	distMatrixData* retVal = new distMatrixData();
	if (useDiskMatrix) {
		retVal->dm = new diskMatrix(options::cacheDir, matrixSize);
	}

	retVal->sequenceNames = dl->getSequenceNames();
	matrixSize = dl->getSequenceCount();
	// process data
	if (options::distMethod == "jc") {
		if (options::verbose) {
			cerr << "Using JC algorithm to calculate distances" << endl;
		}
		JCdistance* alg = new JCdistance(options::verbose, options::fastdist, dl, retVal->dm);
		alg->computeDistanceMatrix(numCores);
		retVal->matrix = alg->getDistanceMatrix();
		delete alg;
	}
	else if (options::distMethod == "kim" || options::distMethod == "") {
		if (options::verbose) {
			cerr << "Using Kimura algorithm to calculate distances" << endl;
		}
		KimuraDistance* alg = new KimuraDistance(options::verbose, options::fastdist, dl, retVal->dm);
		alg->computeDistances(numCores);
		retVal->matrix = alg->getDistanceMatrix();
		delete alg;
	}
	else {
		cerr << "ERROR: Unknown sequence evolution model" << endl;
		exit(1);
	}

	if (printMatrix) {
		//Print the matrix
		if (retVal->dm == NULL) {
			printDistanceMatrix(out, retVal);
		}
		else {
			printDistanceMatrixDisk(out, retVal);
			delete retVal->dm;
		}
	}
	return retVal;
}

distMatrixData* computeDistanceMatrix(bool useDiskMatrix, ostream& out, bool printMatrix, dataloader* dl, distType** distMatrix) {
	/*if(options::gpu){
	computeDistanceMatrixGPU();
	return;
	}*/

	if (options::fastdist && options::verbose) {
		cerr << "Fastdist is enabled" << endl;
	}

	distMatrixData* retVal = new distMatrixData();
	if (useDiskMatrix) {
		retVal->dm = new diskMatrix(options::cacheDir, matrixSize);
	}

	retVal->sequenceNames = dl->getSequenceNames();
	matrixSize = dl->getSequenceCount();
	// process data
	if (options::distMethod == "jc") {
		if (options::verbose) {
			cerr << "Using JC algorithm to calculate distances" << endl;
		}
		JCdistance* alg = new JCdistance(options::verbose, options::fastdist, dl, retVal->dm, distMatrix);
		alg->computeDistanceMatrix(numCores);
		retVal->matrix = alg->getDistanceMatrix();
		delete alg;
	}
	else if (options::distMethod == "kim" || options::distMethod == "") {
		if (options::verbose) {
			cerr << "Using Kimura algorithm to calculate distances" << endl;
		}
		KimuraDistance* alg = new KimuraDistance(options::verbose, options::fastdist, dl, distMatrix, retVal->dm);
		alg->computeDistances(numCores);
		retVal->matrix = alg->getDistanceMatrix();
		delete alg;
	}
	else {
		cerr << "ERROR: Unknown sequence evolution model" << endl;
		exit(1);
	}

	if (printMatrix) {
		//Print the matrix
		if (retVal->dm == NULL) {
			printDistanceMatrix(out, retVal);
		}
		else {
			printDistanceMatrixDisk(out, retVal);
			delete retVal->dm;
		}
	}
	return retVal;
}


distMatrixReader* getDistanceMatrixData(ostream& out, bool halfMatrix, dataloader* dl, vector<string>* sequenceNames, distType** distanceMatrix) {
	distMatrixReader* reader;
	if (distanceMatrixInput) {
		reader = new distMatrixReader(options::verbose, options::fileName, matrixSize, halfMatrix);
		if (options::verbose) {
			cerr << "Reading distance matrix... \n";
		}
		reader->read_data(NULL);
	}
	else {
		if (!distanceMatrixFromPointer)
		{
			if (options::verbose) {
				cerr << "Computing distance matrix... \n";
			}
			distMatrixData* matrixData = computeDistanceMatrix(false, out, false, dl);
			reader = new distMatrixReader(options::verbose, matrixSize, halfMatrix, matrixData->sequenceNames, matrixData->matrix);
			reader->initializeData();
			delete matrixData;
		}
		else
		{
			reader = new distMatrixReader(options::verbose, matrixSize, halfMatrix, sequenceNames, distanceMatrix);
		}
	}
	return reader;
}

polytree* runRapidNJ(int matrixSize, distMatrixReader* reader, ProgressBar* pb, bool deleteAfterwards) {
	if (options::verbose) {
		cerr << "Computing phylogetic tree... \n";
	}
	rapidNJ* sorted = new rapidNJ(reader, matrixSize, options::negative_branches, pb);
	polytree* tree = sorted->run();

	//TODO delete rd

	if (deleteAfterwards)
	{
		delete sorted;
	}

	return tree;
}

polytree* runSimpleNJ(distMatrixReader* reader, ProgressBar* pb) {
	simpleNJ* njs = new simpleNJ(reader, matrixSize, options::negative_branches, pb);
	polytree* tree = njs->run();
	delete njs;
	return tree;
}

polytree* runRapidMNJ(int sortedMatrixSize, distMatrixReader* reader, ProgressBar* pb, bool deleteAfterwards) {
	if (options::verbose) {
		cerr << "Computing phylogetic tree... \n";
	}
	rapidNJMem* nj = new rapidNJMem(reader, matrixSize, sortedMatrixSize, options::verbose, options::negative_branches, pb);
	polytree* tree = nj->run();

	if (deleteAfterwards)
	{
		delete nj;
	}

	return tree;
}

polytree* runDiskNJ(ostream& out, int datastructureSize, dataloader* dl, ProgressBar* pb) {
	if (options::verbose) {
		cerr << "Reading data... \n";
	}
	rdDataInitialiser* reader;

	if (distanceMatrixInput) {
		reader = new rdDataInitialiser(options::verbose, datastructureSize, options::cacheDir, options::fileName);
		bool status = reader->read_data();
		if (!status) {
			cerr << "Could not read distance matrix in file " << options::fileName << endl;
			exit(1);
		}
	}
	else {
		distMatrixData* matrixData = computeDistanceMatrix(true, out, false, dl);
		reader = new rdDataInitialiser(options::verbose, datastructureSize, options::cacheDir, matrixSize);
		reader->initializeFromExistingMatrix(matrixData->sequenceNames, matrixData->dm);
		delete matrixData;
	}

	if (options::verbose) {
		cerr << "Computing phylogetic tree... \n";
	}
	rapidNJDisk* rd = new rapidNJDisk(reader, options::verbose, options::negative_branches, pb);
	polytree* tree = rd->run();
	delete rd;
	delete reader;
	return tree;
}

polytree* computeTree(ostream& out, dataloader* dl, ProgressBar* pb, vector<string>* sequenceNames, distType** distanceMatrix, bool halfMatrix) {
	// Compute the memory requirements for the three different algorithms.
	bool autoDecide = true;
	double matrixSized = (double)matrixSize;
	double systemMemory = getMemSize();
	double matrixMemUsage = ((double)sizeof(distType)) * matrixSized * matrixSized;
	double sortedMatrixMemUsage = matrixSized * matrixSized * ((double)sizeof(cluster_pair));
	int sortedMatrixSize = (int)((systemMemory - (matrixMemUsage / 2.0)) / (matrixSized * sizeof(cluster_pair)));
	sortedMatrixSize = min(sortedMatrixSize, matrixSize);
	int diskSortedMatrixSize = (int)(systemMemory / (distType)(matrixSize * (sizeof(cluster_pair) + sizeof(distType))));
	diskSortedMatrixSize = min(diskSortedMatrixSize, matrixSize);
	diskSortedMatrixSize = max(diskSortedMatrixSize, min(5, matrixSize));
	polytree* tree;

	// select which algorithm to use based on either parameters or memory requirements
	if (options::rapidNJ || options::cacheDir != "" || options::percentageMemoryUsage != "" || options::simpleNJ) {
		autoDecide = false;
	}
	if (options::verbose) {
		cerr << "Matrix size: " << matrixSize << endl;
		cerr << (systemMemory / 1024 / 1024 / 0.8) << " MB of memory is available" << endl;
	}

	if (options::rapidNJ || (autoDecide && sortedMatrixMemUsage + matrixMemUsage <= systemMemory)) {
		if (options::verbose) {
			cerr << "Using RapidNJ \n";
			cerr << "Using " << (matrixMemUsage / 1024 / 1024) << " MB for distance matrix" << endl;
			cerr << "Using " << (sortedMatrixMemUsage / 1024 / 1024) << " MB for sortedMatrix" << endl;
			cerr << "Total memory consumption is " << (matrixMemUsage + sortedMatrixMemUsage) / 1024 / 1024 << " MB" << endl;
		}
		if (sortedMatrixMemUsage > systemMemory - matrixMemUsage) {
			cerr << "WARNING: There's not enough memory to use RapidNJ. Consider using another algorithm." << endl;
		}
		distMatrixReader* reader = getDistanceMatrixData(out, (distanceMatrix != NULL && halfMatrix), dl, sequenceNames, distanceMatrix);

		if (distanceMatrix != NULL && halfMatrix)
		{
			tree = runRapidMNJ(matrixSize, reader, pb, !distanceMatrixFromPointer);
		}
		else
		{
			tree = runRapidNJ(matrixSize, reader, pb, !distanceMatrixFromPointer);
		}
	}
	else if (options::percentageMemoryUsage != "" || (autoDecide && sortedMatrixSize >= matrixSize * MIN_SORTED_MATRIX_SIZE)) {
		if (options::verbose) {
			cerr << "Using Memory efficient RapidNJ \n";
		}
		if (options::percentageMemoryUsage != "") {
			// try to use the user supplied argument
			int percentage;
			percentage = atoi(options::percentageMemoryUsage.data());
			if (percentage < 0 || percentage > 100) {
				cerr << "The memory use percentage must be >=0 and <=100 " << endl;
				exit(1);
			}
			int tempSize = (int)(matrixSize * (percentage / 100.0));
			if (tempSize > sortedMatrixSize) {
				cerr << "WARNING: Not enough memory for " << percentage << "% of the sorted matrix. Reduce the size of the sorted matrix or use RapidDiskNJ." << endl;
			}
			sortedMatrixSize = tempSize;
		}
		if (sortedMatrixSize < matrixSize * MIN_SORTED_MATRIX_SIZE) {
			cerr << "WARNING: the amount of available memory is too low for the memory efficient RapidNJ algorithm to run efficiently. Consider using RapidDiskNJ." << endl;
		}
		if (sortedMatrixSize < 1) {
			sortedMatrixSize = 1;
		}
		if (options::verbose) {
			cerr << "Sorted matrix has " << sortedMatrixSize << " columns" << endl;
		}
		distMatrixReader* reader = getDistanceMatrixData(out, true, dl, sequenceNames, distanceMatrix);
		tree = runRapidMNJ(sortedMatrixSize, reader, pb, !distanceMatrixFromPointer);
	}
	else if (options::simpleNJ) {
		if (options::verbose) {
			cerr << "Using naive NJ \n";
		}
		distMatrixReader* reader = getDistanceMatrixData(out, (distanceMatrix != NULL && halfMatrix), dl, sequenceNames, distanceMatrix);
		tree = runSimpleNJ(reader, pb);
	}
	else {
		if (options::verbose) {
			cerr << "Using RapidDiskNJ algorithm\n";
			cerr << "Sorted matrix has " << diskSortedMatrixSize << " columns" << endl;
		}
		tree = runDiskNJ(out, diskSortedMatrixSize, dl, pb);
	}

	return tree;
}

void bootstrapTree(ostream& out, polytree* tree, dataloader* dl, ProgressBar* pb) {
	for (int i = 0; i < options::replicates; i++) {
		dl->sample_sequences();
		pb->childProgress(1.0 / (options::replicates + 1.0));
		polytree* replicate = computeTree(out, dl, pb, NULL, NULL, false);
		if (options::verbose) {
			cerr << "Comparing trees..." << endl;
		}
		//cout << "---------------------" << i << "-------------------------" << endl;
		tree->compareTreeBootstrap(replicate);
		delete replicate;
	}
	//cout << endl;
}

// Actual methods

extern "C"
{

	DLL_PUBLIC void BuildTreeFromAlignment(int maxMemory, int distance, int numCores, int bootstrapReplicates, int inputType, bool allowNegativeBranches, int inputSequenceCount, int inputSequenceLength, int* inputSequenceNamesLengths, char** inputSequenceNames, char** inputSequenceData, progress_callback callback, return_callback returnCallback, bool verbose)
	{
		options::fastdist = true;
		distanceMatrixInput = false;
		distanceMatrixFromPointer = false;
		options::verbose = verbose;

		options::memSize = maxMemory;

		if (distance == 0)
		{
			options::distMethod = "jc";
		}
		else
		{
			options::distMethod = "kim";
		}

		
		options::cores = numCores;
		options::replicates = bootstrapReplicates;

		InputType type;

		switch (inputType)
		{
		case 0:
			type = DNA;
			break;
		case 1:
			type = PROTEIN;
			break;
		default:
			type = UNKNOWN;
			break;
		}

		options::negative_branches = !allowNegativeBranches;

		configureNumberOfCores();

		dataloaderPointer* pointerDL = new dataloaderPointer(type, inputSequenceCount, inputSequenceLength, inputSequenceNamesLengths, inputSequenceNames, inputSequenceData);

		matrixSize = pointerDL->getSequenceCount();

		ProgressBar* myPB = new ProgressBar(callback);

		if (options::replicates > -1)
		{
			myPB->childProgress(1.0 / (options::replicates + 1.0));
		}

		ostringstream myOut;

		polytree* myTree = computeTree(myOut, pointerDL, myPB, NULL, NULL, false);

		if (options::replicates > -1)
		{
			bootstrapTree(myOut, myTree, pointerDL, myPB);
			myTree->serialize_tree(myOut);
		}
		else
		{
			myTree->serialize_tree(myOut);
		}

		delete myTree;
		delete myPB;
		delete pointerDL;

		string treeString = myOut.str();

		returnCallback(treeString.length(), treeString.c_str());
	}


	DLL_PUBLIC void BuildDistanceMatrixFromAlignment(int maxMemory, int distance, int numCores, int inputType, int inputSequenceCount, int inputSequenceLength, int* inputSequenceNamesLengths, char** inputSequenceNames, char** inputSequenceData, distType** outputMatrix, bool verbose)
	{
		options::fastdist = true;
		distanceMatrixInput = false;
		distanceMatrixFromPointer = false;
		options::verbose = verbose;

		options::memSize = maxMemory;

		if (distance == 0)
		{
			options::distMethod = "jc";
		}
		else
		{
			options::distMethod = "kim";
		}


		options::cores = numCores;
		options::replicates = -1;

		InputType type;

		switch (inputType)
		{
		case 0:
			type = DNA;
			break;
		case 1:
			type = PROTEIN;
			break;
		default:
			type = UNKNOWN;
			break;
		}

		configureNumberOfCores();

		dataloaderPointer* pointerDL = new dataloaderPointer(type, inputSequenceCount, inputSequenceLength, inputSequenceNamesLengths, inputSequenceNames, inputSequenceData);

		matrixSize = pointerDL->getSequenceCount();

		ostream out(cout.rdbuf());

		distMatrixData* computedMatrix = computeDistanceMatrix(false, out, false, pointerDL, outputMatrix);

		delete computedMatrix;
	}

	DLL_PUBLIC void BuildTreeFromDistanceMatrix(int maxMemory, int numCores, bool allowNegativeBranches, int inputSequenceCount, int* inputSequenceNamesLengths, char** inputSequenceNames, bool halfMatrix, distType** distMatrix, progress_callback callback, return_callback returnCallback, bool verbose)
	{
		options::fastdist = true;
		distanceMatrixInput = false;
		distanceMatrixFromPointer = true;
		options::verbose = verbose;

		options::memSize = maxMemory;

		options::cores = numCores;
		options::replicates = -1;

		options::negative_branches = !allowNegativeBranches;

		configureNumberOfCores();

		matrixSize = inputSequenceCount;

		vector<string>* sequenceNames = new vector<string>();

		for (int i = 0; i < inputSequenceCount; i++)
		{
			sequenceNames->push_back(string(inputSequenceNames[i], inputSequenceNamesLengths[i]));
		}

		ProgressBar* myPB = new ProgressBar(callback);

		ostringstream myOut;

		polytree* myTree = computeTree(myOut, NULL, myPB, sequenceNames, distMatrix, halfMatrix);

		myTree->serialize_tree(myOut);

		delete myTree;
		delete myPB;
		delete sequenceNames;

		string treeString = myOut.str();

		returnCallback(treeString.length(), treeString.c_str());
	}
}



