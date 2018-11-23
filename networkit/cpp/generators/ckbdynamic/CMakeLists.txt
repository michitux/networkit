networkit_add_module(generators_ckbdynamic
  CKBDynamicImpl.cpp
  Community.cpp
  CommunityChangeEvent.cpp
  CommunityBirthEvent.cpp
  CommunityMergeEvent.cpp
  CommunityDeathEvent.cpp
  CommunitySplitEvent.cpp
  EventStreamGenerator.cpp
  PowerlawCommunitySizeDistribution.cpp
)

networkit_module_link_modules(generators_ckbdynamic
  auxiliary base dynamics)
