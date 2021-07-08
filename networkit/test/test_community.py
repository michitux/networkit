#!/usr/bin/env python3
import unittest
import os

from networkit import *


class Test_CommunityDetection(unittest.TestCase):

	def setUp(self):
		self.G = readGraph("input/PGPgiantcompo.graph",Format.METIS)

	def test_PLM(self):
		comms = community.detectCommunities(self.G, community.PLM(self.G))

		# each node must be assigned
		for v in self.G.iterNodes():
			self.assertTrue(comms.contains(v))

	def test_CutClustering(self):
		jazz = readGraph("input/jazz.graph",Format.METIS)

		cc = community.CutClustering(jazz,0.5076142131979697)

		comms = community.detectCommunities(jazz, cc)
		self.assertEqual(193, comms.numberOfSubsets())

		hierarchy = community.CutClustering.getClusterHierarchy(jazz)
		self.assertEqual(3, len(hierarchy))

	def test_EgoSplitting(self):
		lfrGraph = readGraph("input/lfr_small.graph", Format.METIS)

		egoSplitting = community.EgoSplitting(lfrGraph)
		egoSplitting.run()
		cover = egoSplitting.getCover()
		self.assertGreater(cover.numberOfSubsets(), 5)
		for size in cover.subsetSizes():
			self.assertGreater(size, 4)
			self.assertLess(size, 40)

		plm = community.PLM(lfrGraph, True, 1.0, "none randomized")
		mapeq = community.LouvainMapEquation(lfrGraph, True, 16, "relaxmap")
		egoSplitting = community.EgoSplitting(lfrGraph, plm, mapeq)
		egoSplitting.run()
		cover = egoSplitting.getCover()
		self.assertGreater(cover.numberOfSubsets(), 5)
		for size in cover.subsetSizes():
			self.assertGreater(size, 4)
			self.assertLess(size, 40)

if __name__ == "__main__":
	unittest.main()
