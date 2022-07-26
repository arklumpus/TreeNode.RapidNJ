using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace PhyloTree
{
    /// <summary>
    /// Contains method to compute neighbour-joining trees using the rapidNJ algorithm (doi:10.1007/978-3-540-87361-7_10). Methods in this class are not thread-safe!
    /// </summary>
    public static class RapidNJ
    {
        /// <summary>
        /// The sequence evolution model used to compute the distance matrix from the alignment.
        /// </summary>
        public enum EvolutionModel
        {
            /// <summary>
            /// Jukes-Cantor model (assuming equal nucleotide/amino acid frequencies and substitution rates).
            /// </summary>
            JukesCantor = 0,

            /// <summary>
            /// For DNA alignments, this represents the Kimura 1980 model (assuming equal base frequencies, and unequal transition/transversion rates). For protein alignments, this represents the Kimura 1983 model that approximates PAM distances based only on the fraction of differing amino acids.
            /// </summary>
            Kimura = 1
        }

        /// <summary>
        /// The kind of sequences in the alignment.
        /// </summary>
        public enum AlignmentType
        {
            /// <summary>
            /// DNA sequences.
            /// </summary>
            DNA = 0,

            /// <summary>
            /// Protein sequences.
            /// </summary>
            Protein = 1,

            /// <summary>
            /// The kind of sequences will be determined based on the first sequence.
            /// </summary>
            Autodetect = 2
        }

        private delegate void ProgressCallback(double progress);
        private delegate void ReturnCallback(long length, IntPtr byteData);

        [DllImport("rapidNJWrapper", CallingConvention = CallingConvention.Cdecl)]
        private static extern void BuildTreeFromAlignment(int maxMemory, int distance, int numCores, int bootstrapReplicates, int inputType, bool allowNegativeBranches, int inputSequenceCount, int inputSequenceLength, IntPtr inputSequenceNamesLengths, IntPtr inputSequenceNames, IntPtr inputSequenceData, [MarshalAs(UnmanagedType.FunctionPtr)] ProgressCallback callback, [MarshalAs(UnmanagedType.FunctionPtr)] ReturnCallback returnCallback, bool verbose);

        [DllImport("rapidNJWrapper", CallingConvention = CallingConvention.Cdecl)]
        private static extern void BuildDistanceMatrixFromAlignment(int maxMemory, int distance, int numCores, int inputType, int inputSequenceCount, int inputSequenceLength, IntPtr inputSequenceNamesLengths, IntPtr inputSequenceNames, IntPtr inputSequenceData, IntPtr outputMatrix, bool verbose);

        [DllImport("rapidNJWrapper", CallingConvention = CallingConvention.Cdecl)]
        private static extern void BuildTreeFromDistanceMatrix(int maxMemory, int numCores, bool allowNegativeBranches, int inputSequenceCount, IntPtr inputSequenceNamesLengths, IntPtr inputSequenceNames, bool halfMatrix, IntPtr distMatrix, [MarshalAs(UnmanagedType.FunctionPtr)] ProgressCallback callback, [MarshalAs(UnmanagedType.FunctionPtr)] ReturnCallback returnCallback, bool verbose);

        /// <summary>
        /// Builds a neighbour-joining tree from a sequence alignment.
        /// </summary>
        /// <param name="alignment">The sequence alignment.</param>
        /// <param name="evolutionModel">The sequence evolution model to use to compute the distance matrix.</param>
        /// <param name="bootstrapReplicates">The number of bootstrap replicates to use for assessing branch support. If this is ≤ 0, no bootstrap is performed.</param>
        /// <param name="alignmentType">The type of sequences in the alignment.</param>
        /// <param name="allowNegativeBranches">If this is <see langword="true"/> (the default), the tree may contain branches with negative length. If this is <see langword="false"/>, negative branch lengths are prevented.</param>
        /// <param name="verbose">If this is <see langword="true"/>, diagnostic information is printed to the standard error.</param>
        /// <param name="numCores">The number of concurrent threads to use to compute the distance matrix.</param>
        /// <param name="progress">A callback used to report progress while building the tree.</param>
        /// <returns>A neighbour-joining tree created from the supplied <paramref name="alignment"/>.</returns>
        /// <exception cref="RapidNJException">Thrown if the sequences in the alignment do not all have the same length.</exception>
        public static TreeNode BuildTreeFromAlignment(Dictionary<string, string> alignment, EvolutionModel evolutionModel = EvolutionModel.Kimura, int bootstrapReplicates = 0, AlignmentType alignmentType = AlignmentType.Autodetect, bool allowNegativeBranches = true, bool verbose = false, int numCores = 0, Action<double> progress = null)
        {
            if (bootstrapReplicates <= 0)
            {
                bootstrapReplicates = -1;
            }

            if (numCores <= 0)
            {
                numCores = Environment.ProcessorCount;
            }

            if (alignmentType == AlignmentType.Autodetect)
            {
                alignmentType = AlignmentType.DNA;

                foreach (KeyValuePair<string, string> kvp in alignment)
                {
                    foreach (char c in kvp.Value)
                    {
                        if (!(c < 65 ||
                              c == 'A' ||
                              c == 'a' ||
                              c == 'C' ||
                              c == 'c' ||
                              c == 'G' ||
                              c == 'g' ||
                              c == 'T' ||
                              c == 't' ||
                              c == 'U' ||
                              c == 'u' ||
                              c == 'N' ||
                              c == 'n'))
                        {
                            alignmentType = AlignmentType.Protein;
                            break;
                        }
                    }

                    break;
                }
            }

            int inputSequenceLength = -1;
            int inputSequenceCount = alignment.Count;
            int[] inputSequenceNamesLengths = new int[inputSequenceCount];
            byte[][] inputSequenceNames = new byte[inputSequenceCount][];
            byte[][] inputSequenceData = new byte[inputSequenceCount][];


            int index = 0;
            foreach (KeyValuePair<string, string> kvp in alignment)
            {
                if (inputSequenceLength < 0)
                {
                    inputSequenceLength = kvp.Value.Length;
                }
                else
                {
                    if (kvp.Value.Length != inputSequenceLength)
                    {
                        throw new RapidNJException("The sequences in the alignment have different lenghts!");
                    }
                }

                inputSequenceNamesLengths[index] = kvp.Key.Length;
                inputSequenceNames[index] = Encoding.UTF8.GetBytes(kvp.Key);
                inputSequenceData[index] = Encoding.UTF8.GetBytes(kvp.Value);

                index++;
            }

            GCHandle[] inputSequenceNamesHandles = new GCHandle[inputSequenceCount];
            GCHandle[] inputSequenceDataHandles = new GCHandle[inputSequenceCount];

            IntPtr[] inputSequenceNamesAddresses = new IntPtr[inputSequenceCount];
            IntPtr[] inputSequenceDataAddresses = new IntPtr[inputSequenceCount];

            for (int i = 0; i < inputSequenceCount; i++)
            {
                inputSequenceNamesHandles[i] = GCHandle.Alloc(inputSequenceNames[i], GCHandleType.Pinned);
                inputSequenceDataHandles[i] = GCHandle.Alloc(inputSequenceData[i], GCHandleType.Pinned);

                inputSequenceNamesAddresses[i] = inputSequenceNamesHandles[i].AddrOfPinnedObject();
                inputSequenceDataAddresses[i] = inputSequenceDataHandles[i].AddrOfPinnedObject();
            }

            GCHandle inputSequenceNamesLengthsHandle = GCHandle.Alloc(inputSequenceNamesLengths, GCHandleType.Pinned);
            GCHandle inputSequenceNamesHandle = GCHandle.Alloc(inputSequenceNamesAddresses, GCHandleType.Pinned);
            GCHandle inputSequenceDataHandle = GCHandle.Alloc(inputSequenceDataAddresses, GCHandleType.Pinned);

            IntPtr inputSequenceNamesLengthsAddress = inputSequenceNamesLengthsHandle.AddrOfPinnedObject();
            IntPtr inputSequenceNamesAddress = inputSequenceNamesHandle.AddrOfPinnedObject();
            IntPtr inputSequenceDataAddress = inputSequenceDataHandle.AddrOfPinnedObject();

            try
            {
                string tree = null;

                unsafe
                {
                    BuildTreeFromAlignment(int.MaxValue, (int)evolutionModel, numCores, bootstrapReplicates, (int)alignmentType, allowNegativeBranches, inputSequenceCount, inputSequenceLength, inputSequenceNamesLengthsAddress, inputSequenceNamesAddress, inputSequenceDataAddress, prog =>
                    {
                        progress?.Invoke(prog / 100);
                    },

                    (length, address) =>
                    {
                        tree = Encoding.UTF8.GetString((byte*)address, (int)length);
                    }, verbose);
                }

                return Formats.NWKA.ParseTree(tree);
            }
            finally
            {
                inputSequenceDataHandle.Free();
                inputSequenceNamesHandle.Free();
                inputSequenceNamesLengthsHandle.Free();

                for (int i = 0; i < inputSequenceCount; i++)
                {
                    inputSequenceNamesHandles[i].Free();
                    inputSequenceDataHandles[i].Free();
                }

            }
        }

        /// <summary>
        /// Builds a neighbour-joining tree from a sequence alignment.
        /// </summary>
        /// <param name="sequenceNames">The names of the sequences/taxa in the alignment.</param>
        /// <param name="sequences">The sequences in the alignment.</param>
        /// <param name="evolutionModel">The sequence evolution model to use to compute the distance matrix.</param>
        /// <param name="bootstrapReplicates">The number of bootstrap replicates to use for assessing branch support. If this is ≤ 0, no bootstrap is performed.</param>
        /// <param name="alignmentType">The type of sequences in the alignment.</param>
        /// <param name="allowNegativeBranches">If this is <see langword="true"/> (the default), the tree may contain branches with negative length. If this is <see langword="false"/>, negative branch lengths are prevented.</param>
        /// <param name="verbose">If this is <see langword="true"/>, diagnostic information is printed to the standard error.</param>
        /// <param name="numCores">The number of concurrent threads to use to compute the distance matrix.</param>
        /// <param name="progress">A callback used to report progress while building the tree.</param>
        /// <returns>A neighbour-joining tree created from the supplied sequence alignment.</returns>
        /// <exception cref="IndexOutOfRangeException">Thrown if the number of sequences and the number of sequence names are different.</exception>
        /// <exception cref="RapidNJException">Thrown if the sequences in the alignment do not all have the same length.</exception>
        public static TreeNode BuildTreeFromAlignment(IReadOnlyList<string> sequenceNames, IReadOnlyList<string> sequences, EvolutionModel evolutionModel = EvolutionModel.Kimura, int bootstrapReplicates = 0, AlignmentType alignmentType = AlignmentType.Autodetect, bool allowNegativeBranches = true, bool verbose = false, int numCores = 0, Action<double> progress = null)
        {
            if (sequences.Count != sequenceNames.Count)
            {
                throw new IndexOutOfRangeException("The number of sequences is not equal to the number of sequence names!");
            }

            if (bootstrapReplicates <= 0)
            {
                bootstrapReplicates = -1;
            }

            if (numCores <= 0)
            {
                numCores = Environment.ProcessorCount;
            }

            if (alignmentType == AlignmentType.Autodetect)
            {
                alignmentType = AlignmentType.DNA;

                foreach (char c in sequences[0])
                {
                    if (!(c < 65 ||
                          c == 'A' ||
                          c == 'a' ||
                          c == 'C' ||
                          c == 'c' ||
                          c == 'G' ||
                          c == 'g' ||
                          c == 'T' ||
                          c == 't' ||
                          c == 'U' ||
                          c == 'u' ||
                          c == 'N' ||
                          c == 'n'))
                    {
                        alignmentType = AlignmentType.Protein;
                        break;
                    }
                }
            }

            int inputSequenceLength = -1;
            int inputSequenceCount = sequences.Count;
            int[] inputSequenceNamesLengths = new int[inputSequenceCount];
            byte[][] inputSequenceNames = new byte[inputSequenceCount][];
            byte[][] inputSequenceData = new byte[inputSequenceCount][];

            for (int i = 0; i < inputSequenceCount; i++)
            {
                if (inputSequenceLength < 0)
                {
                    inputSequenceLength = sequences[i].Length;
                }
                else
                {
                    if (sequences[i].Length != inputSequenceLength)
                    {
                        throw new RapidNJException("The sequences in the alignment have different lenghts!");
                    }
                }

                inputSequenceNamesLengths[i] = sequenceNames[i].Length;
                inputSequenceNames[i] = Encoding.UTF8.GetBytes(sequenceNames[i]);
                inputSequenceData[i] = Encoding.UTF8.GetBytes(sequences[i]);
            }

            GCHandle[] inputSequenceNamesHandles = new GCHandle[inputSequenceCount];
            GCHandle[] inputSequenceDataHandles = new GCHandle[inputSequenceCount];

            IntPtr[] inputSequenceNamesAddresses = new IntPtr[inputSequenceCount];
            IntPtr[] inputSequenceDataAddresses = new IntPtr[inputSequenceCount];

            for (int i = 0; i < inputSequenceCount; i++)
            {
                inputSequenceNamesHandles[i] = GCHandle.Alloc(inputSequenceNames[i], GCHandleType.Pinned);
                inputSequenceDataHandles[i] = GCHandle.Alloc(inputSequenceData[i], GCHandleType.Pinned);

                inputSequenceNamesAddresses[i] = inputSequenceNamesHandles[i].AddrOfPinnedObject();
                inputSequenceDataAddresses[i] = inputSequenceDataHandles[i].AddrOfPinnedObject();
            }

            GCHandle inputSequenceNamesLengthsHandle = GCHandle.Alloc(inputSequenceNamesLengths, GCHandleType.Pinned);
            GCHandle inputSequenceNamesHandle = GCHandle.Alloc(inputSequenceNamesAddresses, GCHandleType.Pinned);
            GCHandle inputSequenceDataHandle = GCHandle.Alloc(inputSequenceDataAddresses, GCHandleType.Pinned);

            IntPtr inputSequenceNamesLengthsAddress = inputSequenceNamesLengthsHandle.AddrOfPinnedObject();
            IntPtr inputSequenceNamesAddress = inputSequenceNamesHandle.AddrOfPinnedObject();
            IntPtr inputSequenceDataAddress = inputSequenceDataHandle.AddrOfPinnedObject();

            try
            {
                string tree = null;

                unsafe
                {
                    BuildTreeFromAlignment(int.MaxValue, (int)evolutionModel, numCores, bootstrapReplicates, (int)alignmentType, allowNegativeBranches, inputSequenceCount, inputSequenceLength, inputSequenceNamesLengthsAddress, inputSequenceNamesAddress, inputSequenceDataAddress, prog =>
                    {
                        progress?.Invoke(prog / 100);
                    },

                    (length, address) =>
                    {
                        tree = Encoding.UTF8.GetString((byte*)address, (int)length);
                    }, verbose);
                }

                return Formats.NWKA.ParseTree(tree);
            }
            finally
            {
                inputSequenceDataHandle.Free();
                inputSequenceNamesHandle.Free();
                inputSequenceNamesLengthsHandle.Free();

                for (int i = 0; i < inputSequenceCount; i++)
                {
                    inputSequenceNamesHandles[i].Free();
                    inputSequenceDataHandles[i].Free();
                }

            }
        }

        /// <summary>
        /// Computes a distance matrix from a sequence alignment.
        /// </summary>
        /// <param name="alignment">The sequence alignment.</param>
        /// <param name="evolutionModel">The sequence evolution model to use to compute the distance matrix.</param>
        /// <param name="alignmentType">The type of sequences in the alignment.</param>
        /// <param name="verbose">If this is <see langword="true"/>, diagnostic information is printed to the standard error.</param>
        /// <param name="numCores">The number of concurrent threads to use to compute the distance matrix.</param>
        /// <returns>A distance matrix built from the <paramref name="alignment"/> using the specified <paramref name="evolutionModel"/>. This will be a full distance matrix (i.e., <c>n</c> rows of <c>n</c> elements, where <c>n</c> is the number of sequences in the alignment).</returns>
        /// <exception cref="RapidNJException">Thrown if the sequences in the alignment do not all have the same length.</exception>
        public static float[][] BuildDistanceMatrixFromAlignment(Dictionary<string, string> alignment, EvolutionModel evolutionModel = EvolutionModel.Kimura, AlignmentType alignmentType = AlignmentType.Autodetect, bool verbose = false, int numCores = 0)
        {
            if (numCores <= 0)
            {
                numCores = Environment.ProcessorCount;
            }

            if (alignmentType == AlignmentType.Autodetect)
            {
                alignmentType = AlignmentType.DNA;

                foreach (KeyValuePair<string, string> kvp in alignment)
                {
                    foreach (char c in kvp.Value)
                    {
                        if (!(c < 65 ||
                              c == 'A' ||
                              c == 'a' ||
                              c == 'C' ||
                              c == 'c' ||
                              c == 'G' ||
                              c == 'g' ||
                              c == 'T' ||
                              c == 't' ||
                              c == 'U' ||
                              c == 'u' ||
                              c == 'N' ||
                              c == 'n'))
                        {
                            alignmentType = AlignmentType.Protein;
                            break;
                        }
                    }

                    break;
                }
            }

            int inputSequenceLength = -1;
            int inputSequenceCount = alignment.Count;
            int[] inputSequenceNamesLengths = new int[inputSequenceCount];
            byte[][] inputSequenceNames = new byte[inputSequenceCount][];
            byte[][] inputSequenceData = new byte[inputSequenceCount][];


            int index = 0;
            foreach (KeyValuePair<string, string> kvp in alignment)
            {
                if (inputSequenceLength < 0)
                {
                    inputSequenceLength = kvp.Value.Length;
                }
                else
                {
                    if (kvp.Value.Length != inputSequenceLength)
                    {
                        throw new RapidNJException("The sequences in the alignment have different lenghts!");
                    }
                }

                inputSequenceNamesLengths[index] = kvp.Key.Length;
                inputSequenceNames[index] = Encoding.UTF8.GetBytes(kvp.Key);
                inputSequenceData[index] = Encoding.UTF8.GetBytes(kvp.Value);

                index++;
            }

            GCHandle[] inputSequenceNamesHandles = new GCHandle[inputSequenceCount];
            GCHandle[] inputSequenceDataHandles = new GCHandle[inputSequenceCount];

            IntPtr[] inputSequenceNamesAddresses = new IntPtr[inputSequenceCount];
            IntPtr[] inputSequenceDataAddresses = new IntPtr[inputSequenceCount];

            for (int i = 0; i < inputSequenceCount; i++)
            {
                inputSequenceNamesHandles[i] = GCHandle.Alloc(inputSequenceNames[i], GCHandleType.Pinned);
                inputSequenceDataHandles[i] = GCHandle.Alloc(inputSequenceData[i], GCHandleType.Pinned);

                inputSequenceNamesAddresses[i] = inputSequenceNamesHandles[i].AddrOfPinnedObject();
                inputSequenceDataAddresses[i] = inputSequenceDataHandles[i].AddrOfPinnedObject();
            }

            GCHandle inputSequenceNamesLengthsHandle = GCHandle.Alloc(inputSequenceNamesLengths, GCHandleType.Pinned);
            GCHandle inputSequenceNamesHandle = GCHandle.Alloc(inputSequenceNamesAddresses, GCHandleType.Pinned);
            GCHandle inputSequenceDataHandle = GCHandle.Alloc(inputSequenceDataAddresses, GCHandleType.Pinned);

            IntPtr inputSequenceNamesLengthsAddress = inputSequenceNamesLengthsHandle.AddrOfPinnedObject();
            IntPtr inputSequenceNamesAddress = inputSequenceNamesHandle.AddrOfPinnedObject();
            IntPtr inputSequenceDataAddress = inputSequenceDataHandle.AddrOfPinnedObject();

            float[][] tbr = new float[inputSequenceCount][];

            GCHandle[] matrixRowHandles = new GCHandle[inputSequenceCount];
            IntPtr[] matrixRowAddresses = new IntPtr[inputSequenceCount];

            for (int i = 0; i < tbr.Length; i++)
            {
                tbr[i] = new float[inputSequenceCount];

                matrixRowHandles[i] = GCHandle.Alloc(tbr[i], GCHandleType.Pinned);
                matrixRowAddresses[i] = matrixRowHandles[i].AddrOfPinnedObject();
            }

            GCHandle matrixHandle = GCHandle.Alloc(matrixRowAddresses, GCHandleType.Pinned);
            IntPtr matrixAddress = matrixHandle.AddrOfPinnedObject();


            try
            {
                BuildDistanceMatrixFromAlignment(int.MaxValue, (int)evolutionModel, numCores, (int)alignmentType, inputSequenceCount, inputSequenceLength, inputSequenceNamesLengthsAddress, inputSequenceNamesAddress, inputSequenceDataAddress, matrixAddress, verbose);

                return tbr;
            }
            finally
            {
                matrixHandle.Free();

                inputSequenceDataHandle.Free();
                inputSequenceNamesHandle.Free();
                inputSequenceNamesLengthsHandle.Free();

                for (int i = 0; i < inputSequenceCount; i++)
                {
                    matrixRowHandles[i].Free();
                    inputSequenceNamesHandles[i].Free();
                    inputSequenceDataHandles[i].Free();
                }

            }
        }

        /// <summary>
        /// Computes a distance matrix from a sequence alignment.
        /// </summary>
        /// <param name="sequenceNames">The names of the sequences/taxa in the alignment.</param>
        /// <param name="sequences">The sequences in the alignment.</param>
        /// <param name="evolutionModel">The sequence evolution model to use to compute the distance matrix.</param>
        /// <param name="alignmentType">The type of sequences in the alignment.</param>
        /// <param name="verbose">If this is <see langword="true"/>, diagnostic information is printed to the standard error.</param>
        /// <param name="numCores">The number of concurrent threads to use to compute the distance matrix.</param>
        /// <returns>A distance matrix built from the alignment using the specified <paramref name="evolutionModel"/>. This will be a full distance matrix (i.e., <c>n</c> rows of <c>n</c> elements, where <c>n</c> is the number of sequences in the alignment).</returns>
        /// <exception cref="IndexOutOfRangeException">Thrown if the number of sequences and the number of sequence names are different.</exception>
        /// <exception cref="RapidNJException">Thrown if the sequences in the alignment do not all have the same length.</exception>
        public static float[][] BuildDistanceMatrixFromAlignment(IReadOnlyList<string> sequenceNames, IReadOnlyList<string> sequences, EvolutionModel evolutionModel = EvolutionModel.Kimura, AlignmentType alignmentType = AlignmentType.Autodetect, bool verbose = false, int numCores = 0)
        {
            if (sequences.Count != sequenceNames.Count)
            {
                throw new IndexOutOfRangeException("The number of sequences is not equal to the number of sequence names!");
            }

            if (numCores <= 0)
            {
                numCores = Environment.ProcessorCount;
            }

            if (alignmentType == AlignmentType.Autodetect)
            {
                alignmentType = AlignmentType.DNA;

                foreach (char c in sequences[0])
                {
                    if (!(c < 65 ||
                          c == 'A' ||
                          c == 'a' ||
                          c == 'C' ||
                          c == 'c' ||
                          c == 'G' ||
                          c == 'g' ||
                          c == 'T' ||
                          c == 't' ||
                          c == 'U' ||
                          c == 'u' ||
                          c == 'N' ||
                          c == 'n'))
                    {
                        alignmentType = AlignmentType.Protein;
                        break;
                    }
                }
            }

            int inputSequenceLength = -1;
            int inputSequenceCount = sequences.Count;
            int[] inputSequenceNamesLengths = new int[inputSequenceCount];
            byte[][] inputSequenceNames = new byte[inputSequenceCount][];
            byte[][] inputSequenceData = new byte[inputSequenceCount][];

            for (int i = 0; i < inputSequenceCount; i++)
            {
                if (inputSequenceLength < 0)
                {
                    inputSequenceLength = sequences[i].Length;
                }
                else
                {
                    if (sequences[i].Length != inputSequenceLength)
                    {
                        throw new RapidNJException("The sequences in the alignment have different lenghts!");
                    }
                }

                inputSequenceNamesLengths[i] = sequenceNames[i].Length;
                inputSequenceNames[i] = Encoding.UTF8.GetBytes(sequenceNames[i]);
                inputSequenceData[i] = Encoding.UTF8.GetBytes(sequences[i]);
            }

            GCHandle[] inputSequenceNamesHandles = new GCHandle[inputSequenceCount];
            GCHandle[] inputSequenceDataHandles = new GCHandle[inputSequenceCount];

            IntPtr[] inputSequenceNamesAddresses = new IntPtr[inputSequenceCount];
            IntPtr[] inputSequenceDataAddresses = new IntPtr[inputSequenceCount];

            for (int i = 0; i < inputSequenceCount; i++)
            {
                inputSequenceNamesHandles[i] = GCHandle.Alloc(inputSequenceNames[i], GCHandleType.Pinned);
                inputSequenceDataHandles[i] = GCHandle.Alloc(inputSequenceData[i], GCHandleType.Pinned);

                inputSequenceNamesAddresses[i] = inputSequenceNamesHandles[i].AddrOfPinnedObject();
                inputSequenceDataAddresses[i] = inputSequenceDataHandles[i].AddrOfPinnedObject();
            }

            GCHandle inputSequenceNamesLengthsHandle = GCHandle.Alloc(inputSequenceNamesLengths, GCHandleType.Pinned);
            GCHandle inputSequenceNamesHandle = GCHandle.Alloc(inputSequenceNamesAddresses, GCHandleType.Pinned);
            GCHandle inputSequenceDataHandle = GCHandle.Alloc(inputSequenceDataAddresses, GCHandleType.Pinned);

            IntPtr inputSequenceNamesLengthsAddress = inputSequenceNamesLengthsHandle.AddrOfPinnedObject();
            IntPtr inputSequenceNamesAddress = inputSequenceNamesHandle.AddrOfPinnedObject();
            IntPtr inputSequenceDataAddress = inputSequenceDataHandle.AddrOfPinnedObject();

            float[][] tbr = new float[inputSequenceCount][];

            GCHandle[] matrixRowHandles = new GCHandle[inputSequenceCount];
            IntPtr[] matrixRowAddresses = new IntPtr[inputSequenceCount];

            for (int i = 0; i < tbr.Length; i++)
            {
                tbr[i] = new float[inputSequenceCount];

                matrixRowHandles[i] = GCHandle.Alloc(tbr[i], GCHandleType.Pinned);
                matrixRowAddresses[i] = matrixRowHandles[i].AddrOfPinnedObject();
            }

            GCHandle matrixHandle = GCHandle.Alloc(matrixRowAddresses, GCHandleType.Pinned);
            IntPtr matrixAddress = matrixHandle.AddrOfPinnedObject();


            try
            {
                BuildDistanceMatrixFromAlignment(int.MaxValue, (int)evolutionModel, numCores, (int)alignmentType, inputSequenceCount, inputSequenceLength, inputSequenceNamesLengthsAddress, inputSequenceNamesAddress, inputSequenceDataAddress, matrixAddress, verbose);

                return tbr;
            }
            finally
            {
                matrixHandle.Free();

                inputSequenceDataHandle.Free();
                inputSequenceNamesHandle.Free();
                inputSequenceNamesLengthsHandle.Free();

                for (int i = 0; i < inputSequenceCount; i++)
                {
                    matrixRowHandles[i].Free();
                    inputSequenceNamesHandles[i].Free();
                    inputSequenceDataHandles[i].Free();
                }

            }
        }

        /// <summary>
        /// Builds a neighbour-joining tree from a distance matrix.
        /// </summary>
        /// <param name="distanceMatrix">The distance matrix used to build the tree. This can either be a full distance matrix (i.e., <c>n</c> rows of <c>n</c> elements, where <c>n</c> is the number of sequences in the alignment) or a lower triangular matrix (i.e., <c>n</c> rows, each of which contains <c>i + 1</c> elements, where <c>i</c> is the 0-based index of the row).</param>
        /// <param name="sequenceNames">The names of the taxa that will appear in the tree.</param>
        /// <param name="copyMatrix">If this is <see langword="true"/> (the default), the distance matrix will be copied before building the tree. If this is <see langword="false"/>, the matrix will be pinned and used in-place. This is faster and requires less memory, but note that that the matrix will be modified!</param>
        /// <param name="allowNegativeBranches">If this is <see langword="true"/> (the default), the tree may contain branches with negative length. If this is <see langword="false"/>, negative branch lengths are prevented.</param>
        /// <param name="verbose">If this is <see langword="true"/>, diagnostic information is printed to the standard error.</param>
        /// <param name="numCores">The number of concurrent threads to use to compute the distance matrix.</param>
        /// <param name="progress">A callback used to report progress while building the tree.</param>
        /// <returns>A neighbour-joining tree built from the supplied <paramref name="distanceMatrix"/>.</returns>
        /// <exception cref="RapidNJException">Thrown if the distance matrix is neither a full matrix nor a lower triangular matrix.</exception>
        /// <exception cref="IndexOutOfRangeException">Thrown if the number of sequence names does not correspond to the size of the distance matrix.</exception>
        public static TreeNode BuildTreeFromDistanceMatrix(float[][] distanceMatrix, IReadOnlyList<string> sequenceNames, bool copyMatrix = true, bool allowNegativeBranches = true, bool verbose = false, int numCores = 0, Action<double> progress = null)
        {
            if (numCores <= 0)
            {
                numCores = Environment.ProcessorCount;
            }

            int inputSequenceCount = sequenceNames.Count;
            int[] inputSequenceNamesLengths = new int[inputSequenceCount];
            byte[][] inputSequenceNames = new byte[inputSequenceCount][];


            bool isFullMatrix = true;
            bool isHalfMatrix = true;

            for (int i = 0; i < sequenceNames.Count; i++)
            {
                inputSequenceNamesLengths[i] = sequenceNames[i].Length;
                inputSequenceNames[i] = Encoding.UTF8.GetBytes(sequenceNames[i]);

                if (distanceMatrix[i].Length != inputSequenceCount)
                {
                    isFullMatrix = false;
                }

                if (distanceMatrix[i].Length != i + 1)
                {
                    isHalfMatrix = false;
                }
            }

            if (isFullMatrix == isHalfMatrix)
            {
                throw new RapidNJException("The distance matrix is neither a full matrix nor a lower triangular matrix!");
            }

            if (distanceMatrix.Length != sequenceNames.Count)
            {
                throw new RapidNJException("The number of sequence names does not correspond to the size of the distance matrix!");
            }

            GCHandle[] inputSequenceNamesHandles = new GCHandle[inputSequenceCount];
            IntPtr[] inputSequenceNamesAddresses = new IntPtr[inputSequenceCount];

            for (int i = 0; i < inputSequenceCount; i++)
            {
                inputSequenceNamesHandles[i] = GCHandle.Alloc(inputSequenceNames[i], GCHandleType.Pinned);
                inputSequenceNamesAddresses[i] = inputSequenceNamesHandles[i].AddrOfPinnedObject();
            }

            GCHandle inputSequenceNamesLengthsHandle = GCHandle.Alloc(inputSequenceNamesLengths, GCHandleType.Pinned);
            GCHandle inputSequenceNamesHandle = GCHandle.Alloc(inputSequenceNamesAddresses, GCHandleType.Pinned);

            IntPtr inputSequenceNamesLengthsAddress = inputSequenceNamesLengthsHandle.AddrOfPinnedObject();
            IntPtr inputSequenceNamesAddress = inputSequenceNamesHandle.AddrOfPinnedObject();

            GCHandle matrixHandle = default;
            GCHandle[] matrixRowHandles = null;
            IntPtr matrixAddress;
            IntPtr matrixStorage = IntPtr.Zero;

            if (!copyMatrix)
            {
                matrixRowHandles = new GCHandle[inputSequenceCount];
                IntPtr[] matrixRowAddresses = new IntPtr[inputSequenceCount];

                for (int i = 0; i < distanceMatrix.Length; i++)
                {
                    matrixRowHandles[i] = GCHandle.Alloc(distanceMatrix[i], GCHandleType.Pinned);
                    matrixRowAddresses[i] = matrixRowHandles[i].AddrOfPinnedObject();
                }

                matrixHandle = GCHandle.Alloc(matrixRowAddresses, GCHandleType.Pinned);
                matrixAddress = matrixHandle.AddrOfPinnedObject();
            }
            else
            {
                unsafe
                {
                    int size;

                    int matrixSize = isFullMatrix ? (inputSequenceCount * inputSequenceCount) : (inputSequenceCount * (inputSequenceCount + 1) / 2);

                    size = sizeof(float) * matrixSize + sizeof(IntPtr) * inputSequenceCount;

                    matrixStorage = Marshal.AllocHGlobal(size);
                    matrixAddress = IntPtr.Add(matrixStorage, sizeof(float) * matrixSize);

                    int index = 0;
                    float* matrixPointer = (float*)matrixStorage;
                    IntPtr* rowPointers = (IntPtr*)matrixAddress;

                    if (isFullMatrix)
                    {
                        for (int i = 0; i < inputSequenceCount; i++)
                        {
                            rowPointers[i] = IntPtr.Add(matrixStorage, sizeof(float) * index);
                            for (int j = 0; j < inputSequenceCount; j++)
                            {
                                matrixPointer[index] = distanceMatrix[i][j];
                                index++;
                            }
                        }
                    }
                    else
                    {
                        for (int i = 0; i < inputSequenceCount; i++)
                        {
                            rowPointers[i] = IntPtr.Add(matrixStorage, sizeof(float) * index);
                            for (int j = 0; j <= i; j++)
                            {
                                matrixPointer[index] = distanceMatrix[i][j];

                                index++;
                            }
                        }
                    }
                }
            }

            try
            {
                string tree = null;

                unsafe
                {
                    BuildTreeFromDistanceMatrix(int.MaxValue, numCores, allowNegativeBranches, inputSequenceCount, inputSequenceNamesLengthsAddress, inputSequenceNamesAddress, isHalfMatrix, matrixAddress, prog =>
                    {
                        progress?.Invoke(prog / 100);
                    },

                    (length, address) =>
                    {
                        tree = Encoding.UTF8.GetString((byte*)address, (int)length);
                    }, verbose);
                }

                return Formats.NWKA.ParseTree(tree);
            }
            finally
            {
                inputSequenceNamesHandle.Free();
                inputSequenceNamesLengthsHandle.Free();

                if (!copyMatrix)
                {
                    matrixHandle.Free();
                }
                else
                {
                    Marshal.FreeHGlobal(matrixStorage);
                }

                for (int i = 0; i < inputSequenceCount; i++)
                {
                    if (!copyMatrix)
                    {
                        matrixRowHandles[i].Free();
                    }
                    inputSequenceNamesHandles[i].Free();
                }

            }
        }

        /// <summary>
        /// Builds a neighbour-joining tree from a distance matrix.
        /// </summary>
        /// <param name="distanceMatrix">The distance matrix used to build the tree. This must be a full, symmetrical distance matrix (i.e., <c>n</c> rows of <c>n</c> elements, where <c>n</c> is the number of sequences in the alignment). If the matrix is not symmetrical, no exception will be thrown, but the reconstructed tree may not be correct.</param>
        /// <param name="sequenceNames">The names of the taxa that will appear in the tree.</param>
        /// <param name="copyMatrix">If this is <see langword="true"/> (the default), the distance matrix will be copied before building the tree. If this is <see langword="false"/>, the matrix will be pinned and used in-place. This is faster and requires less memory, but note that that the matrix will be modified!</param>
        /// <param name="allowNegativeBranches">If this is <see langword="true"/> (the default), the tree may contain branches with negative length. If this is <see langword="false"/>, negative branch lengths are prevented.</param>
        /// <param name="verbose">If this is <see langword="true"/>, diagnostic information is printed to the standard error.</param>
        /// <param name="numCores">The number of concurrent threads to use to compute the distance matrix.</param>
        /// <param name="progress">A callback used to report progress while building the tree.</param>
        /// <returns>A neighbour-joining tree built from the supplied <paramref name="distanceMatrix"/>.</returns>
        /// <exception cref="RapidNJException">Thrown if the distance matrix is neither a full matrix nor a lower triangular matrix.</exception>
        /// <exception cref="IndexOutOfRangeException">Thrown if the number of sequence names does not correspond to the size of the distance matrix.</exception>
        public static TreeNode BuildTreeFromDistanceMatrix(float[,] distanceMatrix, IReadOnlyList<string> sequenceNames, bool copyMatrix = true, bool allowNegativeBranches = true, bool verbose = false, int numCores = 0, Action<double> progress = null)
        {
            if (numCores <= 0)
            {
                numCores = Environment.ProcessorCount;
            }

            int inputSequenceCount = sequenceNames.Count;
            int[] inputSequenceNamesLengths = new int[inputSequenceCount];
            byte[][] inputSequenceNames = new byte[inputSequenceCount][];

            if (distanceMatrix.GetLength(0) != distanceMatrix.GetLength(1))
            {
                throw new RapidNJException("The distance matrix is not square!");
            }

            if (distanceMatrix.GetLength(0) != sequenceNames.Count)
            {
                throw new RapidNJException("The number of sequence names does not correspond to the size of the distance matrix!");
            }

            bool isHalfMatrix = false;

            for (int i = 0; i < sequenceNames.Count; i++)
            {
                inputSequenceNamesLengths[i] = sequenceNames[i].Length;
                inputSequenceNames[i] = Encoding.UTF8.GetBytes(sequenceNames[i]);
            }

            GCHandle[] inputSequenceNamesHandles = new GCHandle[inputSequenceCount];
            IntPtr[] inputSequenceNamesAddresses = new IntPtr[inputSequenceCount];

            for (int i = 0; i < inputSequenceCount; i++)
            {
                inputSequenceNamesHandles[i] = GCHandle.Alloc(inputSequenceNames[i], GCHandleType.Pinned);
                inputSequenceNamesAddresses[i] = inputSequenceNamesHandles[i].AddrOfPinnedObject();
            }

            GCHandle inputSequenceNamesLengthsHandle = GCHandle.Alloc(inputSequenceNamesLengths, GCHandleType.Pinned);
            GCHandle inputSequenceNamesHandle = GCHandle.Alloc(inputSequenceNamesAddresses, GCHandleType.Pinned);

            IntPtr inputSequenceNamesLengthsAddress = inputSequenceNamesLengthsHandle.AddrOfPinnedObject();
            IntPtr inputSequenceNamesAddress = inputSequenceNamesHandle.AddrOfPinnedObject();

            GCHandle matrixHandle = default;
            GCHandle actualMatrixHandle = default;
            IntPtr matrixAddress;
            IntPtr matrixStorage = IntPtr.Zero;

            if (!copyMatrix)
            {
                actualMatrixHandle = GCHandle.Alloc(distanceMatrix, GCHandleType.Pinned);
                IntPtr baseAddress = actualMatrixHandle.AddrOfPinnedObject();

                IntPtr[] matrixRowAddresses = new IntPtr[inputSequenceCount];

                for (int i = 0; i < inputSequenceCount; i++)
                {
                    matrixRowAddresses[i] = IntPtr.Add(baseAddress, sizeof(float) * inputSequenceCount * i);
                }

                matrixHandle = GCHandle.Alloc(matrixRowAddresses, GCHandleType.Pinned);
                matrixAddress = matrixHandle.AddrOfPinnedObject();
            }
            else
            {
                unsafe
                {
                    int size;

                    int matrixSize = inputSequenceCount * inputSequenceCount;

                    size = sizeof(float) * matrixSize + sizeof(IntPtr) * inputSequenceCount;

                    matrixStorage = Marshal.AllocHGlobal(size);
                    matrixAddress = IntPtr.Add(matrixStorage, sizeof(float) * matrixSize);

                    int index = 0;
                    float* matrixPointer = (float*)matrixStorage;
                    IntPtr* rowPointers = (IntPtr*)matrixAddress;

                    for (int i = 0; i < inputSequenceCount; i++)
                    {
                        rowPointers[i] = IntPtr.Add(matrixStorage, sizeof(float) * index);
                        for (int j = 0; j < inputSequenceCount; j++)
                        {
                            matrixPointer[index] = distanceMatrix[i, j];
                            index++;
                        }
                    }
                }
            }

            try
            {
                string tree = null;

                unsafe
                {
                    BuildTreeFromDistanceMatrix(int.MaxValue, numCores, allowNegativeBranches, inputSequenceCount, inputSequenceNamesLengthsAddress, inputSequenceNamesAddress, isHalfMatrix, matrixAddress, prog =>
                    {
                        progress?.Invoke(prog / 100);
                    },

                    (length, address) =>
                    {
                        tree = Encoding.UTF8.GetString((byte*)address, (int)length);
                    }, verbose);
                }

                return Formats.NWKA.ParseTree(tree);
            }
            finally
            {
                inputSequenceNamesHandle.Free();
                inputSequenceNamesLengthsHandle.Free();

                if (!copyMatrix)
                {
                    matrixHandle.Free();
                    actualMatrixHandle.Free();
                }
                else
                {
                    Marshal.FreeHGlobal(matrixStorage);
                }

                for (int i = 0; i < inputSequenceCount; i++)
                {
                    inputSequenceNamesHandles[i].Free();
                }

            }
        }

        /// <summary>
        /// Builds a neighbour-joining tree from a distance matrix.
        /// </summary>
        /// <param name="distanceMatrix">The distance matrix used to build the tree. This can either be a full distance matrix (i.e., <c>n</c> rows of <c>n</c> elements, where <c>n</c> is the number of sequences in the alignment) or a lower triangular matrix (i.e., <c>n</c> rows, each of which contains <c>i + 1</c> elements, where <c>i</c> is the 0-based index of the row).</param>
        /// <param name="sequenceNames">The names of the taxa that will appear in the tree.</param>
        /// <param name="copyMatrix">If this is <see langword="true"/> (the default), the distance matrix will be copied before building the tree. If this is <see langword="false"/>, the matrix will be pinned and used in-place. This is faster and requires less memory, but note that that the matrix will be modified!</param>
        /// <param name="allowNegativeBranches">If this is <see langword="true"/> (the default), the tree may contain branches with negative length. If this is <see langword="false"/>, negative branch lengths are prevented.</param>
        /// <param name="verbose">If this is <see langword="true"/>, diagnostic information is printed to the standard error.</param>
        /// <param name="numCores">The number of concurrent threads to use to compute the distance matrix.</param>
        /// <param name="progress">A callback used to report progress while building the tree.</param>
        /// <returns>A neighbour-joining tree built from the supplied <paramref name="distanceMatrix"/>.</returns>
        /// <exception cref="RapidNJException">Thrown if the distance matrix is neither a full matrix nor a lower triangular matrix, or if the number of entries in the matrix does not correspond to the sequence names.</exception>
        public static TreeNode BuildTreeFromDistanceMatrix(Span<float> distanceMatrix, IReadOnlyList<string> sequenceNames, bool copyMatrix = true, bool allowNegativeBranches = true, bool verbose = false, int numCores = 0, Action<double> progress = null)
        {
            if (numCores <= 0)
            {
                numCores = Environment.ProcessorCount;
            }

            int inputSequenceCount = sequenceNames.Count;
            int[] inputSequenceNamesLengths = new int[inputSequenceCount];
            byte[][] inputSequenceNames = new byte[inputSequenceCount][];


            bool isFullMatrix = sequenceNames.Count * sequenceNames.Count == distanceMatrix.Length;
            bool isHalfMatrix = sequenceNames.Count * (sequenceNames.Count + 1) / 2 == distanceMatrix.Length;

            if (isFullMatrix == isHalfMatrix)
            {
                throw new RapidNJException("The distance matrix is neither a full matrix nor a lower triangular matrix, or the number of entries in the matrix does not correspond to the number of sequence names!");
            }

            for (int i = 0; i < sequenceNames.Count; i++)
            {
                inputSequenceNamesLengths[i] = sequenceNames[i].Length;
                inputSequenceNames[i] = Encoding.UTF8.GetBytes(sequenceNames[i]);
            }

            GCHandle[] inputSequenceNamesHandles = new GCHandle[inputSequenceCount];
            IntPtr[] inputSequenceNamesAddresses = new IntPtr[inputSequenceCount];

            for (int i = 0; i < inputSequenceCount; i++)
            {
                inputSequenceNamesHandles[i] = GCHandle.Alloc(inputSequenceNames[i], GCHandleType.Pinned);
                inputSequenceNamesAddresses[i] = inputSequenceNamesHandles[i].AddrOfPinnedObject();
            }

            GCHandle inputSequenceNamesLengthsHandle = GCHandle.Alloc(inputSequenceNamesLengths, GCHandleType.Pinned);
            GCHandle inputSequenceNamesHandle = GCHandle.Alloc(inputSequenceNamesAddresses, GCHandleType.Pinned);

            IntPtr inputSequenceNamesLengthsAddress = inputSequenceNamesLengthsHandle.AddrOfPinnedObject();
            IntPtr inputSequenceNamesAddress = inputSequenceNamesHandle.AddrOfPinnedObject();

            GCHandle matrixHandle = default;
            IntPtr matrixAddress;
            IntPtr matrixStorage = IntPtr.Zero;

            if (!copyMatrix)
            {
                unsafe
                {
                    fixed (float* source = &MemoryMarshal.GetReference(distanceMatrix))
                    {
                        IntPtr baseAddress = (IntPtr)source;

                        IntPtr[] matrixRowAddresses = new IntPtr[inputSequenceCount];

                        if (isFullMatrix)
                        {
                            for (int i = 0; i < inputSequenceCount; i++)
                            {
                                matrixRowAddresses[i] = IntPtr.Add(baseAddress, sizeof(float) * inputSequenceCount * i);
                            }
                        }
                        else
                        {
                            for (int i = 0; i < inputSequenceCount; i++)
                            {
                                matrixRowAddresses[i] = IntPtr.Add(baseAddress, sizeof(float) * i * (i + 1) / 2);
                            }
                        }

                        matrixHandle = GCHandle.Alloc(matrixRowAddresses, GCHandleType.Pinned);
                        matrixAddress = matrixHandle.AddrOfPinnedObject();

                        try
                        {
                            string tree = null;

                            unsafe
                            {
                                BuildTreeFromDistanceMatrix(int.MaxValue, numCores, allowNegativeBranches, inputSequenceCount, inputSequenceNamesLengthsAddress, inputSequenceNamesAddress, isHalfMatrix, matrixAddress, prog =>
                                {
                                    progress?.Invoke(prog / 100);
                                },

                                (length, address) =>
                                {
                                    tree = Encoding.UTF8.GetString((byte*)address, (int)length);
                                }, verbose);
                            }

                            return Formats.NWKA.ParseTree(tree);
                        }
                        finally
                        {
                            inputSequenceNamesHandle.Free();
                            inputSequenceNamesLengthsHandle.Free();

                            matrixHandle.Free();

                            for (int i = 0; i < inputSequenceCount; i++)
                            {
                                inputSequenceNamesHandles[i].Free();
                            }

                        }
                    }
                }
            }
            else
            {
                unsafe
                {
                    int size;

                    int matrixSize = isFullMatrix ? (inputSequenceCount * inputSequenceCount) : (inputSequenceCount * (inputSequenceCount + 1) / 2);

                    size = sizeof(float) * matrixSize + sizeof(IntPtr) * inputSequenceCount;

                    matrixStorage = Marshal.AllocHGlobal(size);
                    matrixAddress = IntPtr.Add(matrixStorage, sizeof(float) * matrixSize);

                    fixed (float* source = &MemoryMarshal.GetReference(distanceMatrix))
                    {
                        Buffer.MemoryCopy(source, (void*)matrixStorage, sizeof(float) * matrixSize, sizeof(float) * matrixSize);
                    }

                    IntPtr* rowPointers = (IntPtr*)matrixAddress;

                    if (isFullMatrix)
                    {
                        for (int i = 0; i < inputSequenceCount; i++)
                        {
                            rowPointers[i] = IntPtr.Add(matrixStorage, sizeof(float) * i * inputSequenceCount);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < inputSequenceCount; i++)
                        {
                            rowPointers[i] = IntPtr.Add(matrixStorage, sizeof(float) * i * (i + 1) / 2);
                        }
                    }
                }

                try
                {
                    string tree = null;

                    unsafe
                    {
                        BuildTreeFromDistanceMatrix(int.MaxValue, numCores, allowNegativeBranches, inputSequenceCount, inputSequenceNamesLengthsAddress, inputSequenceNamesAddress, isHalfMatrix, matrixAddress, prog =>
                        {
                            progress?.Invoke(prog / 100);
                        },

                        (length, address) =>
                        {
                            tree = Encoding.UTF8.GetString((byte*)address, (int)length);
                        }, verbose);
                    }

                    return Formats.NWKA.ParseTree(tree);
                }
                finally
                {
                    inputSequenceNamesHandle.Free();
                    inputSequenceNamesLengthsHandle.Free();

                    Marshal.FreeHGlobal(matrixStorage);

                    for (int i = 0; i < inputSequenceCount; i++)
                    {
                        inputSequenceNamesHandles[i].Free();
                    }
                }
            }
        }

        /// <summary>
        /// Builds a neighbour-joining tree from a distance matrix.
        /// </summary>
        /// <param name="distanceMatrix">The distance matrix used to build the tree. This can either be a full distance matrix (i.e., <c>n</c> rows of <c>n</c> elements, where <c>n</c> is the number of sequences in the alignment) or a lower triangular matrix (i.e., <c>n</c> rows, each of which contains <c>i + 1</c> elements, where <c>i</c> is the 0-based index of the row). The matrix will be copied before building the tree.</param>
        /// <param name="sequenceNames">The names of the taxa that will appear in the tree.</param>
        /// <param name="allowNegativeBranches">If this is <see langword="true"/> (the default), the tree may contain branches with negative length. If this is <see langword="false"/>, negative branch lengths are prevented.</param>
        /// <param name="verbose">If this is <see langword="true"/>, diagnostic information is printed to the standard error.</param>
        /// <param name="numCores">The number of concurrent threads to use to compute the distance matrix.</param>
        /// <param name="progress">A callback used to report progress while building the tree.</param>
        /// <returns>A neighbour-joining tree built from the supplied <paramref name="distanceMatrix"/>.</returns>
        /// <exception cref="RapidNJException">Thrown if the distance matrix is neither a full matrix nor a lower triangular matrix.</exception>
        /// <exception cref="RapidNJException">Thrown if the distance matrix is neither a full matrix nor a lower triangular matrix, or if the number of entries in the matrix does not correspond to the sequence names.</exception>
        public static TreeNode BuildTreeFromDistanceMatrix(ReadOnlySpan<float> distanceMatrix, IReadOnlyList<string> sequenceNames, bool allowNegativeBranches = true, bool verbose = false, int numCores = 0, Action<double> progress = null)
        {
            if (numCores <= 0)
            {
                numCores = Environment.ProcessorCount;
            }

            int inputSequenceCount = sequenceNames.Count;
            int[] inputSequenceNamesLengths = new int[inputSequenceCount];
            byte[][] inputSequenceNames = new byte[inputSequenceCount][];


            bool isFullMatrix = sequenceNames.Count * sequenceNames.Count == distanceMatrix.Length;
            bool isHalfMatrix = sequenceNames.Count * (sequenceNames.Count + 1) / 2 == distanceMatrix.Length;

            if (isFullMatrix == isHalfMatrix)
            {
                throw new RapidNJException("The distance matrix is neither a full matrix nor a lower triangular matrix, or the number of entries in the matrix does not correspond to the number of sequence names!");
            }

            for (int i = 0; i < sequenceNames.Count; i++)
            {
                inputSequenceNamesLengths[i] = sequenceNames[i].Length;
                inputSequenceNames[i] = Encoding.UTF8.GetBytes(sequenceNames[i]);
            }

            GCHandle[] inputSequenceNamesHandles = new GCHandle[inputSequenceCount];
            IntPtr[] inputSequenceNamesAddresses = new IntPtr[inputSequenceCount];

            for (int i = 0; i < inputSequenceCount; i++)
            {
                inputSequenceNamesHandles[i] = GCHandle.Alloc(inputSequenceNames[i], GCHandleType.Pinned);
                inputSequenceNamesAddresses[i] = inputSequenceNamesHandles[i].AddrOfPinnedObject();
            }

            GCHandle inputSequenceNamesLengthsHandle = GCHandle.Alloc(inputSequenceNamesLengths, GCHandleType.Pinned);
            GCHandle inputSequenceNamesHandle = GCHandle.Alloc(inputSequenceNamesAddresses, GCHandleType.Pinned);

            IntPtr inputSequenceNamesLengthsAddress = inputSequenceNamesLengthsHandle.AddrOfPinnedObject();
            IntPtr inputSequenceNamesAddress = inputSequenceNamesHandle.AddrOfPinnedObject();

            IntPtr matrixAddress;
            IntPtr matrixStorage = IntPtr.Zero;

            unsafe
            {
                int size;

                int matrixSize = isFullMatrix ? (inputSequenceCount * inputSequenceCount) : (inputSequenceCount * (inputSequenceCount + 1) / 2);

                size = sizeof(float) * matrixSize + sizeof(IntPtr) * inputSequenceCount;

                matrixStorage = Marshal.AllocHGlobal(size);
                matrixAddress = IntPtr.Add(matrixStorage, sizeof(float) * matrixSize);

                fixed (float* source = &MemoryMarshal.GetReference(distanceMatrix))
                {
                    Buffer.MemoryCopy(source, (void*)matrixStorage, sizeof(float) * matrixSize, sizeof(float) * matrixSize);
                }

                IntPtr* rowPointers = (IntPtr*)matrixAddress;

                if (isFullMatrix)
                {
                    for (int i = 0; i < inputSequenceCount; i++)
                    {
                        rowPointers[i] = IntPtr.Add(matrixStorage, sizeof(float) * i * inputSequenceCount);
                    }
                }
                else
                {
                    for (int i = 0; i < inputSequenceCount; i++)
                    {
                        rowPointers[i] = IntPtr.Add(matrixStorage, sizeof(float) * i * (i + 1) / 2);
                    }
                }
            }

            try
            {
                string tree = null;

                unsafe
                {
                    BuildTreeFromDistanceMatrix(int.MaxValue, numCores, allowNegativeBranches, inputSequenceCount, inputSequenceNamesLengthsAddress, inputSequenceNamesAddress, isHalfMatrix, matrixAddress, prog =>
                    {
                        progress?.Invoke(prog / 100);
                    },

                    (length, address) =>
                    {
                        tree = Encoding.UTF8.GetString((byte*)address, (int)length);
                    }, verbose);
                }

                return Formats.NWKA.ParseTree(tree);
            }
            finally
            {
                inputSequenceNamesHandle.Free();
                inputSequenceNamesLengthsHandle.Free();

                Marshal.FreeHGlobal(matrixStorage);

                for (int i = 0; i < inputSequenceCount; i++)
                {
                    inputSequenceNamesHandles[i].Free();
                }
            }
        }

        /// <summary>
        /// Represents errors that occur while preparing the data for computing neighbour-joining trees.
        /// </summary>
        public class RapidNJException : Exception
        {
            /// <summary>
            /// Creates a new <see cref="RapidNJException"/>.
            /// </summary>
            /// <param name="message">The message that describes the error.</param>
            public RapidNJException(string message) : base(message) { }
        }

    }
}
