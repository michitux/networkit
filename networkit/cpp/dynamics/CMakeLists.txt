networkit_add_module(dynamics
    DGSStreamParser.cpp
    DGSWriter.cpp
    GraphDifference.cpp
    GraphEvent.cpp
    CommunityEvent.cpp
    GraphEventHandler.cpp
    GraphEventProxy.cpp
    GraphUpdater.cpp
    CoverUpdater.cpp
    )

networkit_module_link_modules(dynamics
        auxiliary base graph structures)

add_subdirectory(test)

