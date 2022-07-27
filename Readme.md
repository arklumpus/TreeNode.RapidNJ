# TreeNode.RapidNJ: Multiplatform .NET library to build neighbour-joining trees using RapidNJ

<img src="Logo.svg" width="256" align="right">

[![License: AGPL v3](https://img.shields.io/badge/License-GPL_v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Version](https://img.shields.io/nuget/v/TreeNode.RapidNJ)](https://nuget.org/packages/TreeNode.RapidNJ)

__TreeNode.RapidNJ__ is a set of multiplatform .NET bindings for [RapidNJ](https://birc.au.dk/software/rapidnj), making it possible to create neighbour-joining trees in C#. The trees are returned as [`TreeNode`]() objects for easy manipulation.

The library is released under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html) licence.

## Getting started

The TreeNode.RapidNJ library targets .NET Standard 2.1, thus it can be used in projects that target .NET Standard 2.1+, .NET Core 3.0+, .NET 5.0+, and possibly others. MuPDFCore includes a pre-compiled native library, which currently supports the following platforms:

* Windows x64 (64 bit)
* Linux x64 (64 bit)
* Linux arm64/aarch64 (ARM 64 bit)
* macOS Intel x86_64 (64 bit)
* macOS Apple silicon (ARM 64 bit)

To use the library in your project, you should install the [TreeNode.RapidNJ NuGet package](https://www.nuget.org/packages/TreeNode.RapidNJ/). You will probably also want to install the [TreeNode NuGet package](https://www.nuget.org/packages/TreeNode/). When you publish a program that uses TreeNode.RapidNJ, the correct native library for the target architecture will automatically be copied to the build folder.

**Note**: you should make sure that end users on Windows install the [Microsoft Visual C++ Redistributable for Visual Studio 2015, 2017, 2019 and 2022](https://docs.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-160#visual-studio-2015-2017-2019-and-2022) for their platform, otherwise they will get an error message stating that `rapidNJWrapper.dll` could not be loaded because a module was not found.

## Usage

The library provides a static class `PhyloTree.RapidNJ` with some methods to build neighbour-joining trees using RapidNJ. These can be used to build a tree starting from a sequence alignment or from a distance matrix, or to build a distance matrix starting from a sequence alignment.

### Quick example

The following example shows the simplest way to create a tree from a sequence alignment:

```CSharp
// A sequence alignment obtained somehow - e.g., by reading a FASTA file.
Dictionary<string, string> sequenceAlignment = ...

// Create the neighbour-joining tree.
TreeNode tree = PhyloTree.RapidNJ.BuildTreeFromAlignment(sequenceAlignment);

// Print the tree to the console in Newick format.
Console.WriteLine(tree.ToString());
```

The following sections provide more details into the various methods contained in the library and their parameters.

### Building a tree from a sequence alignment

The `BuildTreeFromAlignment` method can be used to build a tree from a sequence alignment. There are two overloads of this method: one that requires the sequence alignment to be provided as a `Dictionary<string, string>` (where the keys are the sequence names, and the values are the aligned sequences), and another that requires the sequence names and sequences to be provided as two separate lists of strings (elements with the same index from the two lists are assumed to belong together). Aside from this, the two methods have the same optional parameters:

```CSharp
public static TreeNode BuildTreeFromAlignment(Dictionary<string, string> alignment, EvolutionModel evolutionModel = EvolutionModel.Kimura, int bootstrapReplicates = 0, AlignmentType alignmentType = AlignmentType.Autodetect, bool allowNegativeBranches = true, bool verbose = false, int numCores = 0, Action<double> progress = null);

public static TreeNode BuildTreeFromAlignment(IReadOnlyList<string> sequenceNames, IReadOnlyList<string> sequences, EvolutionModel evolutionModel = EvolutionModel.Kimura, int bootstrapReplicates = 0, AlignmentType alignmentType = AlignmentType.Autodetect, bool allowNegativeBranches = true, bool verbose = false, int numCores = 0, Action<double> progress = null);
```

The alignment can consist of either DNA or amino acid (protein) sequences, which must all have the same length. By default, the type of alignment is detected automatically, but it can also be specified using an optional parameter. Note that the library provides no facility to read a sequence alignment from a file - you should have your own way e.g. of parsing a FASTA file.

The two methods have the following optional parameters:

* `EvolutionModel evolutionModel`: this parameter determines the evolutionary model used to compute the distance matrix. This is necessary because just measuring the proportion $p$ of sites that differ between the two sequences underestimates the actual evolutionary distance (because multiple substitutions in the same site become more and more probable as time passes). The `EvolutionModel` enum has two values, whose meaning depends on whether the alignment contains DNA or protein sequences. In all cases the default is `EvolutionModel.Kimura`.

  On a DNA alignment, `EvolutionModel.JukesCantor` represents the 1969 Jukes and Cantor model, assuming equal base frequencies and equal rates for every substitution. With this model, the distance between two sequences is computed as $d = -\frac{3}{4} \ln \left(1-\frac{4}{3}p \right)$. If two sequences differ in more than 75% of the sites ($p \geq \frac{3}{4}$), the distance is returned as -1. The `EvolutionModel.Kimura` value, instead, represents the 1980 Kimura model, assuming equal base frequencies, but different rates for transitions and transversions. If $t$ is the proportion of observed transitions and $v$ is the proportion of observed transversions, with this model the distance between two sequences is computed as $d = -\frac{1}{2}\ln \left( \left( 1 - 2t - v \right)\sqrt{1 - 2v} \right) $. If $1 - 2t - v \leq 0$ or $v \geq 0.5$, the distance is returned as -1.

  On a protein alignment, `EvolutionModel.JukesCantor` represents again a model with equal amino acid frequences and equal rates for every substitution, with the distance being computed as $d = -\frac{19}{20} \ln \left(1-\frac{20}{19}p \right)$ (again, if two sequences differ in more than 95% of the sites ($p \geq \frac{19}{20}$), the distance is returned as -1). `EvolutionModel.Kimura` instead represents the 1983 Kimura protein distance formula, which should approximate the Dayhoff model; for sequences differing in less than 75% of the sites ($p \lt 0.75$), the distance is computed as $d = -\ln \left (1 - p - \frac{1}{5} p^2 \right)$. For sequences differing more than this but less than 93% ($0.75 \leq p \leq 0.93$), the distance is extracted from a [lookup table](https://github.com/arklumpus/rapidNJ-library/blob/master/src/distanceCalculation/KimuraDistance.cpp#L34). For sequences differing more than 93%, the distance is returned as 10.

* `int bootstrapReplicates`: if this is greater than 0 (which is the default), the requested number of bootstrap replicates are performed, and the resulting support values are annotated in the `Support` attribute of the nodes of the tree that is returned by the method.

* `AlignmentType alignmentType`: this parameter specifies the alignment type (`AlignmentType.DNA` or `AlignmentType.Protein`). If this `AlignmentType.Autodetect` (the default), the alignment type is inferred by scanning through all the sequences until a letter that is not `A`, `a`, `C`, `c`, `G`, `g`, `T`, `t`, `U`, `u`, `N`, or `n` is found. If such a character can be found, the alignment is assumed to be a protein alignment, otherwise it is assumed to be a DNA alignment.

* `bool allowNegativeBranches`: during the neighbour-joining algorithm computations, it is possible that some branch lengths may have negative values. If this parameter is `true` (the default), nothing is done about this; otherwise, the algorithm is altered to avoid these negative branch lengths.

* `bool verbose`: if this is `true`, the `rapidNJ` library will print some debug information to the console's standard output. Otherwise (the default), nothing is printed.

* `int numCores`: indicates the number of cores/threads to use for the distance matrix computation. If this is 0 (the default), a number equal to the number of processors in the computer (as returned by `Environment.ProcessorCount`) is used. Setting this to a value greater than the number of physical cores on the computer may slow down the computation, rather than speeding it up.

* `Action<double> progress`: a callback that will be used by the library to provide progress updates on the tree-building process. This will be called by the library with a parameter ranging fron 0 to 1. Note that no progress is reported while computing the distance matrix; therefore, some time may pass between when you invoke the `BuildTreeFromAlignmentMethod` and the first time that the callback is invoked.

### Building a distance matrix from a sequence alignment

The `BuildDistanceMatrixFromAlignment` method can be used to create a distance matrix from a sequence alignment. Similar to the `BuildTreeFromAlignment` method, the alignment can be supplied either as a `Dictionary<string, string>`, or as two separate lists of sequence names and sequence data. The two overloads of this method also have some optional parameters:

```CSharp
public static float[][] BuildDistanceMatrixFromAlignment(Dictionary<string, string> alignment, EvolutionModel evolutionModel = EvolutionModel.Kimura, AlignmentType alignmentType = AlignmentType.Autodetect, bool verbose = false, int numCores = 0);

public static float[][] BuildDistanceMatrixFromAlignment(IReadOnlyList<string> sequenceNames, IReadOnlyList<string> sequences, EvolutionModel evolutionModel = EvolutionModel.Kimura, AlignmentType alignmentType = AlignmentType.Autodetect, bool verbose = false, int numCores = 0);
```

The meaning of the optional parameters is exactly the same as for the parameters of the `BuildTreeFromAlignment` method. The distance matrix is returned as a square symmetrical jagged array.

### Building a tree from a distance matrix

Finally, the `BuildTreeFromDistanceMatrix` method can be used to build a tree from a distance matrix. The distance matrix can be supplied in multiple formats:

* As a `float[][]` jagged array (in this case, the matrix should either be square and symmetrical, or a lower triangular matrix);
* As a `float[,]` 2D array (in this case, the matrix should be sequare and symmetrical);
* As a `ReadOnlySpan<float>` or `Span<float>` containing all the elements in the matrix one after the other, row by row (in this case, the matrix should either be square and symmetrical, or a lower triangular matrix). Note that there is an implicit conversion between `T[]` and `Span<T>`, thus this overload can also be used with a matrix stored in a 1D `float[]` array. This overload is also useful if you have a matrix stored in unmanaged memory.

As a result, the `BuildTreeFromDistanceMatrix` method has four overloads:

```CSharp
BuildTreeFromDistanceMatrix(float[][] distanceMatrix, IReadOnlyList<string> sequenceNames, bool copyMatrix = true, bool allowNegativeBranches = true, bool verbose = false, int numCores = 0, Action<double> progress = null);

BuildTreeFromDistanceMatrix(float[,] distanceMatrix, IReadOnlyList<string> sequenceNames, bool copyMatrix = true, bool allowNegativeBranches = true, bool verbose = false, int numCores = 0, Action<double> progress = null);

BuildTreeFromDistanceMatrix(Span<float> distanceMatrix, IReadOnlyList<string> sequenceNames, bool copyMatrix = true, bool allowNegativeBranches = true, bool verbose = false, int numCores = 0, Action<double> progress = null);

BuildTreeFromDistanceMatrix(ReadOnlySpan<float> distanceMatrix, IReadOnlyList<string> sequenceNames, bool allowNegativeBranches = true, bool verbose = false, int numCores = 0, Action<double> progress = null);
```

In addition to the distance matrix, a list of strings that represent the sequence names is also required. This should have the same size as the distance matrix. The optional parameters are the same for all the overloads:

* `bool copyMatrix`: if this is `true` (the default), the distance matrix is copied to unmanaged memory and then passed to the `rapidNJ` library. If this is `false`, the matrix is instead pinned in place (if necessary) and its address is passed to the library directly. As a result, this second case is a bit faster and uses less memory; however, note that the distance matrix will be modified by the library. This means that you should only set `copyMatrix` to `false` if you do not need to do anything else with the distance matrix. The overload where the distance matrix is provided as a `ReadOnlySpan<float>` lacks this parameter, because a `ReadOnlySpan` should not be modified; in this overload, the matrix is always copied.

* `bool allowNegativeBranches`, `bool verbose`, `int numCores` and `Action<double> progress` have the same meaning and effect as the optional parameters with the same name in the `BuildTreeFromAlignment` method.

## Building from source

Building the TreeNode.RapidNJ library from source requires the following steps:

1. Building the `rapidnj` and `rapidNJWrapper` native libraries
2. Creating the `TreeNode.RapidNJ` library NuGet package

Step 1 needs to be performed on all of Windows, macOS and Linux, and on the various possible architectures (x64 for Windows, x64/Intel and arm64/Apple for macOS, and x64 and arm64 for Linux - no cross-compiling)! Otherwise, some native assets will be missing and it will not be possible to build the NuGet package.

### 1. Building `rapidnj` and `rapidNJWrapper`

First of all, you should clone this repository; to do this, you need to run:

```
git clone https://github.com/arklumpus/TreeNode.RapidNJ.git
cd TreeNode.RapidNJ
```

On Linux and macOS you can run this command from the terminal (after ensuring that you have installed git); on Windows you should be able to use the WSL bash terminal or the git command line.

A [fork of rapidNJ](https://github.com/arklumpus/rapidNJ-library) modified to produce a library instead of an executable is included as a submodule of this repository. This is not downloaded automatically when you clone the repository; to download it you need to run:

```
git submodule update --init
```

To compile `rapidNJWrapper`, you will need [CMake](https://cmake.org/) (version 3.8 or higher) and (on Windows) [Ninja](https://ninja-build.org/).

On Windows, the easiest way to get all the required tools is probably to install [Visual Studio](https://visualstudio.microsoft.com/it/). By selecting the "Desktop development with C++" workload you should get everything you need.

On macOS, you will need to install at least the Command-Line Tools for Xcode (if necessary, you should be prompted to do this while you perform the following steps) and CMake.

Once you have everything at the ready, you will have to build the library on the five platforms.

#### Windows

1. Assuming you have installed Visual Studio, you should open the `x64 Native Tools Command Prompt for VS` (you should be able to find this in the Start menu - a normal command prompt will not work).
2. `CD` to the directory where you have downloaded the TreeNode.RapidNJ source code.
3. `CD` into the `native` directory.
4. Type `build`. This will start the `build.cmd` batch script that will delete any previous build and compile `rapidnj` and `rapidNJWrapper`.

After this finishes, you should find a file named `rapidNJWrapper.dll` in the `native/out/build/win-x64/rapidNJWrapper/` directory. Leave it there.

#### macOS and Linux

1. Assuming you have everything ready, open a terminal in the folder where you have downloaded the TreeNode.RapidNJ source code.
2. `cd` into the `native` directory.
3. Type `chmod +x build.sh`.
4. Type `./build.sh`. This will delete any previous build and compile `rapidnj` and `rapidNJWrapper`.

After this finishes, you should find a file named `librapidNJWrapper.dylib` in the `native/out/build/mac-x64/rapidNJWrapper/` directory (on macOS running on an Intel x64 processor) or in the `native/out/build/mac-arm64/rapidNJWrapper/` directory (on macOS running on an Apple silicon arm64 processor), and a file named `librapidNJWrapper.so` in the `native/out/build/linux-x64/rapidNJWrapper/` directory (on Linux running on an x64 processor) or in the `native/out/build/linux-arm64/rapidNJWrapper/` directory (on Linux running on an ARM processor). Leave it there.

### 2. Creating the TreeNode.RapidNJ NuGet package

Once you have the `rapidNJWrapper.dll`, `librapidNJWrapper.dylib` and `librapidNJWrapper.so` files, make sure they are in the correct folders (`native/out/build/xxx-yyy/rapidNJWrapper/`), __all on the same machine__.

To create the  NuGet package, you will need the [.NET Core 3.0 SDK or higher](https://dotnet.microsoft.com/download/dotnet/current) for your platform. Once you have installed it and have everything ready, open a terminal in the folder where you have downloaded the TreeNode.RapidNJ source code and type:

```
cd TreeNode.RapidNJ
dotnet pack -c Release
```

This will create a NuGet package in `TreeNode.RapidNJ/bin/Release`. You can install this package on your projects by adding a local NuGet source.

### 3. Running tests

To verify that everything is working correctly, you should build the TreeNode.RapidNJ test suite and run it on all platforms. To build the test suite, you will need the [.NET 6 SDK or higher](https://dotnet.microsoft.com/download/dotnet/current). You will also need to have enabled the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install).

To build the test suite:

1. Make sure that you have changed the version of the TreeNode.RapidNJ NuGet package so that it is higher than the latest version of TreeNode.RapidNJ in the NuGet repository (you should use a pre-release suffix, e.g. `1.1.0-a1` to avoid future headaches with new versions of TreeNode.RapidNJ). This is set in line 9 of the `TreeNode.RapidNJ/TreeNode.RapidNJ.csproj` file.
2. Add the `TreeNode.RapidNJ/bin/Release` folder to your local NuGet repositories (you can do this e.g. in Visual Studio).
3. If you have not done so already, create the TreeNode.RapidNJ NuGet package following step 2 above.
4. Update line 28 of the `Tests/Tests.csproj` project file so that it refers to the version of the TreeNode.RapidNJ package you have just created.

These steps ensure that you are testing the right version of TreeNode.RapidNJ (i.e. your freshly built copy) and not something else that may have been cached.

Now, open a windows command line in the folder where you have downloaded the TreeNode.RapidNJ source code, type `BuildTests` and press `Enter`. This will create a number of files in the `Release\RapidNJTests` folder, where each file is an archive containing the tests for a certain platform and architecture:

* `RapidNJTests-linux-x64.tar.gz` contains the tests for Linux environments on x64 processors.
* `RapidNJTests-linux-arm64.tar.gz` contains the tests for Linux environments on arm64 processors.
* `RapidNJTests-mac-x64.tar.gz` contains the tests for macOS environments on Intel processors.
* `RapidNJTests-mac-arm64.tar.gz` contains the tests for macOS environments on Apple silicon processors.
* `RapidNJTests-win-x64.tar.gz` contains the tests for Windows environments on x64 processors.

To run the tests, copy each archive to a machine running the corresponding operating system, and extract it. Then:

#### Windows
* Open a command prompt and `CD` into the folder where you have extracted the contents of the test archive.
* Enter the command `RapidNJTestHost` (this will run the test program).

#### macOS and Linux
* Open a terminal and `cd` into the folder where you have extracted the contents of the test archive.
* Enter the command `chmod +x RapidNJTestHost` (this will add the executable flag to the test program).
* Enter the command `./RapidNJTestHost` (this will run the test program).
* On macOS, depending on your security settings, you may get a message saying `zsh: killed` when you try to run the program. To address this, you need to sign the executable, e.g. by running `codesign --timestamp --sign <certificate> RapidNJTestHost`, where `<certificate>` is the name of a code signing certificate in your keychain (e.g. `Developer ID Application: John Smith`). After this, you can try again to run the test program with `./RapidNJTestHost`. Additionally, you may need to open System Preferences, go to `Security and Privacy`, and click on `Allow anyways` _multiple times_. A lot of dialog windows may also show up, and you will have to click on `Open` every time. Essentially, you will have to run the test program multiple times until you have allowed everything and no scary dialog shows up.

The test suite will start; it will print the name of each test, followed by a green `  Succeeded  ` or a red `  Failed  ` depending on the test result. If everything went correctly, all tests should succeed.

When all the tests have been run, the program will print a summary showing how many tests have succeeded (if any) and how many have failed (if any). If any tests have failed, a list of these will be printed, and then they will be run again one at a time, waiting for a key press before running each test (this makes it easier to follow what is going on). If you wish to kill the test process early, you can do so with `CTRL+C`.
