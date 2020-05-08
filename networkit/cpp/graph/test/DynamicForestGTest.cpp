/*
 *
 */

#include "DynamicForestGTest.h"
#include <networkit/graph/DynamicForest.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/components/RandomSpanningForest.hpp>
#include <networkit/community/QuasiThresholdEditingLinear.hpp>
#include <networkit/graph/GraphTools.hpp>


namespace NetworKit {

TEST_F(DynamicForestGTest, testConstruction) {
	Graph G(3, false, true);
	G.addEdge(1, 0);
	G.addEdge(2, 0);

	DynamicForest dynForest(G);

	EXPECT_EQ(0, dynForest.parent(1));
	EXPECT_EQ(0, dynForest.parent(2));
	EXPECT_EQ(none, dynForest.parent(0));
	auto c = dynForest.children(0);
	EXPECT_TRUE((c == std::vector<node>{1, 2} || c == std::vector<node>{2, 1}));
	EXPECT_EQ(0, dynForest.children(1).size());
	EXPECT_EQ(0, dynForest.children(2).size());
	ASSERT_EQ(std::vector<node>{0}, dynForest.children(none));
}

TEST_F(DynamicForestGTest, testNextChild) {
	Graph G(4, false, true);
	for (node c = 1; c < 4; ++c) {
		G.addEdge(c, 0);
	}

	DynamicForest dynForest(G);

	EXPECT_EQ(0, dynForest.nextDFSNodeOnEnter(none, none));
	EXPECT_EQ(1, dynForest.nextDFSNodeOnEnter(0, none));
	EXPECT_EQ(1, dynForest.nextDFSNodeOnEnter(0, 0));

	for (node u = 1; u < 3; ++u) {
		EXPECT_EQ(u + 1, dynForest.nextChild(u, 0));
		EXPECT_EQ(u + 1, dynForest.nextDFSNodeOnEnter(u, 0));
		EXPECT_EQ(u + 1, dynForest.nextDFSNodeOnEnter(u, none));
	}

	EXPECT_EQ(0, dynForest.nextChild(3, 0));

	EXPECT_EQ(1, dynForest.nextDFSNodeOnEnter(1, 1));
}

TEST_F(DynamicForestGTest, testIsolate) {
	Graph G(5, false, true);
	G.addEdge(1, 0);
	G.addEdge(2, 0);
	G.addEdge(3, 1);
	G.addEdge(4, 1);

	DynamicForest dynForest(G);

	dynForest.isolate(1);

	EXPECT_EQ(0, dynForest.parent(3));
	EXPECT_EQ(0, dynForest.parent(4));

	EXPECT_EQ(none, dynForest.parent(1));
	EXPECT_EQ(std::vector<node>(), dynForest.children(1));

	EXPECT_EQ(0, dynForest.parent(2));

	EXPECT_EQ(3, dynForest.children(0).size());
}

TEST_F(DynamicForestGTest, testPathStructures) {
	Graph karate = METISGraphReader().read("input/karate.graph");
	
	karate.indexEdges();
	QuasiThresholdEditingLinear editing(karate);
	editing.run();
	std::vector<node> parent = editing.getParents();
	
	Graph forest = Graph(GraphTools::copyNodes(karate), false, true);
	karate.forNodes([&](node u) {
		if (parent[u] != none) {
			forest.addEdge(u, parent[u]);
		}
	});
	DynamicForest dynForest(forest);
	forest.forNodesInRandomOrder ([&](node nodeToMove){
		TRACE("==================nodeToMove: ", nodeToMove,"=============================");
		forest.forNodesInRandomOrder([&](node parent){
			if(nodeToMove != parent){
				TRACE("==================parent: ", parent,"=============================");
				std::vector<node> children = dynForest.children(parent);
				for(int i = 0; i < children.size(); i++){
					if(children[i] == nodeToMove){
						children.erase(children.begin()+i);
						break;
					}
				}
				//TRACE("==================children: ", children,"=============================");
				count pow_set_size = pow(2, children.size());
				count counter, j;
				std::vector<node> adopted;
				for(counter = 0; counter < pow_set_size; counter++) {
					adopted.clear(); 
					for(j = 0; j < children.size(); j++) {
						if(counter & (1 << j)){
							adopted.push_back(children[j]);
						}
					}
					//TRACE("==================adopt: ", adopted,"=============================");
					dynForest.isolate(nodeToMove);
					dynForest.moveToPosition(nodeToMove, parent, adopted);
			}
			
			}
		});
	});
}
}