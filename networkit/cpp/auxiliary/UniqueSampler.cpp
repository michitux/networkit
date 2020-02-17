#include "UniqueSampler.h"

namespace Aux {
	UniqueSampler::UniqueSampler(count max) : dist(0, max-1), max(max), i(0) {}
	count UniqueSampler::draw() {
		dist.param(std::uniform_int_distribution<count>::param_type(i, max-1));
		const count r = dist(Aux::Random::getURNG());
		count result = r;

		auto rit = replace.find(r);
		if (rit != replace.end()) {
			result = rit->second;
		}

		auto iit = replace.find(i);
		if (iit == replace.end()) {
			replace[r] = i;
		} else {
			// Inserts invalidate all iterators - cache value to avoid invalid read
			count v = iit->second;
			replace[r] = v;
		}

		++i;

		return result;
	}
}
