networkit_add_module(generators_ckbdynamic
  CKBDynamicImpl.cpp
  Community.cpp
  CommunityEventListener.cpp
  CommunityChangeEvent.cpp
  CommunityBirthEvent.cpp
  CommunityMergeEvent.cpp
  CommunityDeathEvent.cpp
  CommunitySplitEvent.cpp
  PowerlawCommunitySizeDistribution.cpp
  BucketSampling.cpp
)

networkit_module_link_modules(generators_ckbdynamic
  auxiliary base dynamics)
