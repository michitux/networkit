#ifndef COMMUNITY_EVENT_H_
#define COMMUNITY_EVENT_H_

#include <string>
#include "../Globals.h"

/**
 * @ingroup dynamics
 */
namespace NetworKit {
	class CommunityEvent {
	public:
		enum Type {
			NODE_JOINS_COMMUNITY,
			NODE_LEAVES_COMMUNITY,
			TIME_STEP
		};

		Type type; //!< type of community event
		node u; //!< node that joins or leaves a community
		index community; //!< the community that is joined or left

		CommunityEvent(Type type = TIME_STEP, node u = none, index community = none);

		std::string toString() const;

		friend bool operator==(const CommunityEvent &a, const CommunityEvent &b);
	};
}

#endif
