#include "CoverUpdater.h"

namespace NetworKit {


CoverUpdater::CoverUpdater(Cover& C) : C(C) {
}

void CoverUpdater::update(const std::vector<CommunityEvent>& stream) {
	for (CommunityEvent ev : stream) {
		switch (ev.type) {
			case CommunityEvent::NODE_JOINS_COMMUNITY : {
				if (ev.community >= C.upperBound()) {
					C.setUpperBound(ev.community + 1);
				}

				while (ev.u >= C.numberOfElements()) {
					C.extend();
				}

				C.addToSubset(ev.community, ev.u);
				break;
			}
			case CommunityEvent::NODE_LEAVES_COMMUNITY : {
				C.removeFromSubset(ev.community, ev.u);
				break;
			}
			case CommunityEvent::TIME_STEP : {
				break;
			}
			default: {
				throw std::runtime_error("unknown event type");
			}
		}
	}
}

} /* namespace NetworKit */
