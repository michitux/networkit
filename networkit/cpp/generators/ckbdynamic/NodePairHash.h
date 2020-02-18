#ifndef CKBDYNAMIC_NODE_PAIR_HASH_H_
#define CKBDYNAMIC_NODE_PAIR_HASH_H_

#include <utility>
#include <functional>
#include <robin_hood.h>
#include "../../Globals.h"

namespace NetworKit {
  namespace CKBDynamicImpl {
    struct NodePairHash {
      std::size_t operator()(const std::pair<node, node>& k) const {
	std::size_t lhs = robin_hood::hash<node>()(k.first);
	std::size_t rhs = robin_hood::hash<node>()(k.second);
	lhs^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
	return lhs;
      }
    };
  }
}

#endif
