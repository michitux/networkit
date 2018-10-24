#ifndef UNIQUE_SAMPLER_H_
#define UNIQUE_SAMPLER_H_

#include "../Globals.h"
#include <random>
#include "Random.h"
#include <unordered_map>

namespace Aux {
	using NetworKit::count;

		class UniqueSampler {
		public:
			UniqueSampler(count max);
			count draw();

		private:
			std::uniform_int_distribution<count> dist;
			count max, i;
			std::unordered_map<count, count> replace;
		};

}

#endif
