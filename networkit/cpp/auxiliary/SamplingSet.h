#ifndef SAMPLING_SET_H_
#define SAMPLING_SET_H_

#include <vector>
#include <robin_hood.h>
#include "Random.h"

namespace Aux {
	template <class Key, class Hash = robin_hood::hash<Key>>
	class SamplingSet {
	private:
		robin_hood::unordered_flat_map<Key, size_t, Hash> positions;
		std::vector<Key> elements;

	public:
		//SamplingSet() { positions.min_load_factor(0.05f); }

		using const_iterator = typename std::vector<Key>::const_iterator;

		const_iterator begin() const { return elements.cbegin(); }
		const_iterator end() const { return elements.cend(); }

		size_t insert(const Key& e) {
			bool inserted = positions.insert({e, elements.size()}).second;

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

		const Key& at(size_t pos) const {
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

		/**
		 * Sample a random item out of the items i..size()-1 that is then stored at position i.
		 *
		 * @param i The position at which the item shall be stored.
		 * @return The sampled item at position i.
		 */
		const Key& sample_item(size_t i) {
			size_t j = Aux::Random::integer(i, elements.size() - 1);
			swap_items(i, j);
			return elements[i];
		}
	private:
		void swap_items(size_t i, size_t j) {
			if (i != j) {
				using std::swap;
				swap(elements[i], elements[j]);
				positions[elements[i]] = i;
				positions[elements[j]] = j;
			}
		}
	};
}

#endif
