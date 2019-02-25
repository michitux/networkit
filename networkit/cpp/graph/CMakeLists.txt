networkit_add_module(graph
    DynamicForest.cpp
    Graph.cpp
    GraphBuilder.cpp
    GraphTools.cpp
    KruskalMSF.cpp
    RandomMaximumSpanningForest.cpp
    SpanningForest.cpp
    UnionMaximumSpanningForest.cpp
    )

networkit_module_link_modules(graph
    auxiliary base structures)

add_subdirectory(test)

