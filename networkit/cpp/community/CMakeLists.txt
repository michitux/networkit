networkit_add_module(community
    AdjustedRandMeasure.cpp
    ClusteringGenerator.cpp
    CommunityDetectionAlgorithm.cpp
    Conductance.cpp
    ConductanceCommunityCleanup.cpp
    CoverHubDominance.cpp
    CoverF1Similarity.cpp
    OverlappingNMIDistance.cpp
    Coverage.cpp
    CutClustering.cpp
    DissimilarityMeasure.cpp
    DynamicNMIDistance.cpp
    EdgeCut.cpp
    GraphClusteringTools.cpp
    GraphStructuralRandMeasure.cpp
    HubDominance.cpp
    IntrapartitionDensity.cpp
    IsolatedInterpartitionConductance.cpp
    IsolatedInterpartitionExpansion.cpp
    JaccardMeasure.cpp
    LPDegreeOrdered.cpp
    LPPotts.cpp
    LocalCoverEvaluation.cpp
    LocalPartitionEvaluation.cpp
    Modularity.cpp
    NMIDistance.cpp
    NodeStructuralRandMeasure.cpp
    PLM.cpp
    PLP.cpp
    ParallelAgglomerativeClusterer.cpp
    PartitionFragmentation.cpp
    PartitionHubDominance.cpp
    PartitionIntersection.cpp
    SampledGraphStructuralRandMeasure.cpp
    SampledNodeStructuralRandMeasure.cpp
    StablePartitionNodes.cpp
	LouvainMapEquation.cpp
    cleanup/MergeCommunities.cpp
    cleanup/SignificanceCommunityCleanUp.cpp
    cleanup/SingleCommunityCleanUp.cpp
    cleanup/StochasticDistributionCalculator.cpp
    cleanup/SignificanceCalculator.cpp
    EgoSplitting.cpp
    )

networkit_module_link_modules(community
    auxiliary base coarsening components flow graph matching structures)

add_subdirectory(test)
add_subdirectory(cleanup/test)

# Enable math-optimizations and thus vectorized log for LouvainMapEquation
if (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
    ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR
    ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang"))
  set_source_files_properties(LouvainMapEquation.cpp DIRECTORY ../../.. PROPERTIES COMPILE_OPTIONS "-ffast-math")
endif()
