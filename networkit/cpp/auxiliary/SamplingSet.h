#ifndef SAMPLING_SET_H_
#define SAMPLING_SET_H_

#include <vector>
#include <tsl/robin_map.h>
#include "Random.h"

namespace Aux {
	template <class Key, class Hash = std::hash<Key>>
	class SamplingSet {
	private:
		tsl::robin_map<Key, size_t, Hash> positions;
		std::vector<Key> elements;

	public:
		using const_iterator = typename std::vector<Key>::const_iterator;

		const_iterator begin() const { return elements.cbegin(); }
		const_iterator end() const { return elements.cend(); }

		size_t insert(const Key& e) {
			bool inserted = positions.try_emplace(e, elements.size()).second;

			if (inserted) {
				elements.push_back(e);
			}

			return inserted;
		}

		size_t erase(const Key& e) {
			auto it = positions.find(e);
			if (it == positions.end()) return 0;

			size_t pos = it->second;

			if (pos != elements.size() - 1) {
				elements[pos] = elements.back();
				positions[elements[pos]] = pos;
			}

			elements.pop_back();
			positions.erase(it);

			return 1;
		}

		bool contains(const Key& e) const {
			return (positions.find(e) != positions.end());
		}

		size_t size() const {
			return elements.size();
		}

		bool empty() const {
			return elements.empty();
		}

		Key at(size_t pos) const {
			return elements[pos];
		}

		void clear() {
			positions.clear();
			elements.clear();
		}

		void reserve(size_t n) {
			positions.reserve(n);
			elements.reserve(n);
		}

		/**
		 * Selects a random sample of @a k elements that are then stored at positions 0..k-1.
		 *
		 * @param k The number of samples
		 */
		void random_sample(size_t k) {
			auto &urng = Aux::Random::getURNG();
			std::uniform_int_distribution<size_t> dist;
			using p_t = std::uniform_int_distribution<size_t>::param_type;
			const size_t n = elements.size();

			assert(k <= n);

			auto swap_items = [&](size_t i, size_t j) {
				if (i != j) {
					std::swap(elements[i], elements[j]);
					positions[elements[i]] = i;
					positions[elements[j]] = j;
				}
			};

			if (k <= n/2 || k < 2) {
				for (size_t i = 0; i < k; ++i) {
					p_t p(i, n - 1);
					size_t j = dist(urng, p);
					swap_items(i, j);
				}
			} else {
				for (size_t i = n - 1; i >= k; --i) {
					p_t p(0, i);
					size_t j = dist(urng, p);
					swap_items(i, j);
				}
			}
		}
	};
}

#endif
