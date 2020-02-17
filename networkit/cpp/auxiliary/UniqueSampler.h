#ifndef UNIQUE_SAMPLER_H_
#define UNIQUE_SAMPLER_H_

#include "../Globals.h"
#include <random>
#include "Random.h"
#include <tsl/robin_map.h>

namespace Aux {
	using NetworKit::count;

		class UniqueSampler {
		public:
			UniqueSampler(count max);
			count draw();

		private:
			std::uniform_int_distribution<count> dist;
			count max, i;
			tsl::robin_map<count, count> replace;
		};

}

#endif
