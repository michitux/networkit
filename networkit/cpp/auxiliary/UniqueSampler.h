#ifndef UNIQUE_SAMPLER_H_
#define UNIQUE_SAMPLER_H_

#include "../Globals.h"
#include <random>
#include "Random.h"
#include <robin_hood.h>

namespace Aux {
	using NetworKit::count;

		class UniqueSampler {
		public:
			UniqueSampler(count max);
			count draw();

		private:
			std::uniform_int_distribution<count> dist;
			count max, i;
			robin_hood::unordered_flat_map<count, count> replace;
		};

}

#endif
