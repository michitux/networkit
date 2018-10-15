#include "CommunityEvent.h"

namespace NetworKit {
	CommunityEvent::CommunityEvent(Type type, node u, index community) : type(type), u(u), community(community) {}

	std::string CommunityEvent::toString() const {
		std::stringstream ss;

		if (type == NODE_JOINS_COMMUNITY) {
			ss << "njc(" << u << "," << community << ")";
		} else if (type == NODE_LEAVES_COMMUNITY) {
			ss << "nlc(" << u << "," << community << ")";
		} else {
			if (type != TIME_STEP) {
				throw std::runtime_error("Invalid type of community event.");
			}

			ss << "st";
		}

		return ss.str();
	}

	bool operator==(const CommunityEvent &a, const CommunityEvent &b) {
		if (a.type != b.type) return false;
		if (a.type == CommunityEvent::TIME_STEP) return true;
		return a.u == b.u && a.community == b.community;
	}
}
