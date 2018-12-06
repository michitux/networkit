networkit_add_module(generators_ckbdynamic
  CKBDynamicImpl.cpp
  Community.cpp
  CommunityChangeEvent.cpp
  CommunityBirthEvent.cpp
  CommunityMergeEvent.cpp
  CommunityDeathEvent.cpp
  CommunitySplitEvent.cpp
  CustomCommunityMembershipDistribution.cpp
  CustomCommunitySizeDistribution.cpp
  EventStreamGenerator.cpp
  PowerlawCommunityMembershipDistribution.cpp
  PowerlawCommunitySizeDistribution.cpp
)

networkit_module_link_modules(generators_ckbdynamic
  auxiliary base dynamics)
