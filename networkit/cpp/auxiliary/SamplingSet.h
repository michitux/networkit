#ifndef SAMPLING_SET_H_
#define SAMPLING_SET_H_

#include <vector>
#include <tsl/robin_map.h>

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
			if (contains(e)) return 0;

			elements.push_back(e);
			try {
				positions[e] = elements.size() - 1;
			} catch (...) {
				elements.pop_back();
				throw;
			}

			return 1;
		}

		size_t erase(const Key& e) {
			auto it = positions.find(e);
			if (it == positions.end()) return 0;

			size_t pos = it->second;
			elements[pos] = elements.back();
			elements.pop_back();
			positions[elements[pos]] = pos;

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
	};
}

#endif
