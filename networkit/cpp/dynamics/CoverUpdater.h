#ifndef COVERUPDATER_H_
#define COVERUPDATER_H_

#include "../structures/Cover.h"
#include "../dynamics/CommunityEvent.h"

namespace NetworKit {

/**
 * @ingroup dynamics
 */
class CoverUpdater {

public:

	CoverUpdater(Cover& C);

	void update(const std::vector<CommunityEvent>& stream);
private:
	Cover& C;
};

} /* namespace NetworKit */

#endif /* GRAPHUPDATER_H_ */
