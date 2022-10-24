import unittest
from src import node, graph, components


# we can make a testing osw file, and then have function read it
# class TestInitialization(unittest.TestCase):
#     def test_(self):
#
#         example_graph = ppi.Graph()
#         accession_List = [
#             "1", "2", "3", "4"
#         ]
#
#         self.assertEqual(True, False)
#
#



class TestCollapse(unittest.TestCase):
    """this class contains test for method used in collapse
    """

    def setUp(self) -> None:
        # make graph
        self.graph_1 = graph.Graph()


        # make proteins objects
        self.pro1 = node.Protein(['1'], '1', 0)
        self.pro2 = node.Protein(['2'], '2', 0)
        self.pro3 = node.Protein(['3'], '3', 0)
        self.pro4 = node.Protein(['4'], '4', 0)
        self.pro5 = node.Protein(['5'], '5', 0)
        self.pro6 = node.Protein(['6'], '6', 0)
        self.pro7 = node.Protein(['7'], '7', 0)
        self.pro8 = node.Protein(['8'], '8', 0)
        self.pro9 = node.Protein(['9'], '9', 0)

        # make peptide objects
        self.pep1 = node.Peptide(['1'], 0.1, 0)
        self.pep2 = node.Peptide(['2'], 0.2, 0)
        self.pep3 = node.Peptide(['3'], 0.3, 0)
        self.pep4 = node.Peptide(['4'], 0.4, 0)
        self.pep5 = node.Peptide(['5'], 0.5, 0)
        self.pep6 = node.Peptide(['6'], 0.6, 0)
        self.pep7 = node.Peptide(['7'], 0.7, 0)
        self.pep8 = node.Peptide(['8'], 0.8, 0)
        self.pep9 = node.Peptide(['9'], 0.9, 0)
        self.pep10 = node.Peptide(['10'], 1.0, 0)

        # fill graph with proteins and their possible constituent peptides
        self.graph_1.node_dict[self.pro1] = [self.pep4, self.pep3, self.pep7,
                                             self.pep9, self.pep8]
        self.graph_1.node_dict[self.pro2] = [self.pep8]
        self.graph_1.node_dict[self.pro3] = [self.pep6]
        self.graph_1.node_dict[self.pro4] = [self.pep2, self.pep1]
        self.graph_1.node_dict[self.pro5] = [self.pep4, self.pep8]
        self.graph_1.node_dict[self.pro6] = [self.pep2, self.pep6]
        self.graph_1.node_dict[self.pro7] = [self.pep1, self.pep5]
        self.graph_1.node_dict[self.pro8] = [self.pep8]
        self.graph_1.node_dict[self.pro9] = [self.pep2, self.pep1]

        # fill graph with peptide and their possible protein
        self.graph_1.node_dict[self.pep1] = [self.pro7]
        self.graph_1.node_dict[self.pep2] = [self.pro4, self.pro6, self.pro9]
        self.graph_1.node_dict[self.pep3] = [self.pro1]
        self.graph_1.node_dict[self.pep4] = [self.pro1, self.pro5]
        self.graph_1.node_dict[self.pep5] = [self.pro7]
        self.graph_1.node_dict[self.pep6] = [self.pro3, self.pro6]
        self.graph_1.node_dict[self.pep7] = [self.pro1]
        self.graph_1.node_dict[self.pep8] = [self.pro1, self.pro5, self.pro2,
                                             self.pro8]
        self.graph_1.node_dict[self.pep9] = [self.pro1]
        self.graph_1.node_dict[self.pep10] = [self.pro4, self.pro9]

        # neighbours of pep8 (all node with degree
        self.protein_list = [self.pro1, self.pro5, self.pro2, self.pro8]
        # neighbours of pro1
        self.peptide_list = [self.pep4, self.pep3, self.pep7, self.pep9,
                             self.pep8]

    def test_categorize_node_degree(self):
        """
        reorder neighbours take a list of peptide and proteins
        and then
        :return:
        """
        """testing categorize_node_degree and peptides
        important thing to check about reorder node: that each node is in the 
        right category"""

        # run the function
        categorized_nodes = self.graph_1.categorize_node_degree(
            list(self.graph_1.node_dict))

        # then test
        # skip 0, since there are never nodes that have 0 degree (those were
        # filtered out when graph was initialized)
        for i in range(1, len(categorized_nodes)):
            node_degree_i_list = categorized_nodes[i]
            for a_node in node_degree_i_list:
                # for every node, check that their degree is 'i'
                num_neighbours = len(self.graph_1.node_dict[a_node])
                self.assertEqual(num_neighbours, i)

    def test_grouping_recursion(self):
        # run the previous function

        # this part is identical to the actual code
        categorized_nodes = self.graph_1.categorize_node_degree(
            list(self.graph_1.node_dict))

        max_degree = len(categorized_nodes)

        for i in range(1, len(categorized_nodes)):
            node_degree_i_list = categorized_nodes[i]

            if len(node_degree_i_list) == 0:
                print(i, max_degree, "no node")
                continue

            specific_dict = {}
            self.graph_1.build_specific_node_dict(node_degree_i_list, specific_dict)

            node_degree_i_list.sort()
            self.graph_1.grouping_recursion(node_degree_i_list, specific_dict)
        # identical part ends

        # check that the correct nodes were grouped together

        # TODO: i need to check this
        # I think making change to protein object in the node dict
        # also make change to the protein object
        # because I fill the node dict with that object

        # either protein 2 has id '2' and '8' or
        # protein 8 has id '2' and '8'
        # also I want to use XOR

        # checking all multi ones merged correctly,
        # thought this only check one of them is deleted, test should be rewritten so that
        # it test one of them is merged, and all other ones are deleted
        self.assertTrue(self.pro2.get_id().sort() == ['2', '8'] or
                        self.pro8.get_id().sort() == ['2', '8'])

        # and that the other one is deleted
        self.assertTrue(
            (self.pro2.get_first_id() + self.pro2.get_target_decoy()) in self.graph_1.node_to_delete or
            (self.pro8.get_first_id() + self.pro8.get_target_decoy()) in self.graph_1.node_to_delete
        )

        self.assertTrue(self.pro4.get_id().sort() == ['4', '9'] or
                        self.pro9.get_id().sort() == ['4', '9'])

        self.assertTrue(
            (self.pro4.get_first_id() + self.pro4.get_target_decoy()) in self.graph_1.node_to_delete or
            (self.pro9.get_first_id() + self.pro9.get_target_decoy()) in self.graph_1.node_to_delete
        )

        self.assertTrue(self.pep1.get_id().sort() == ['1', '5'] or
                        self.pep5.get_id().sort() == ['1', '5'])

        self.assertTrue(
            (self.pro1.get_first_id() + self.pro1.get_target_decoy()) in self.graph_1.node_to_delete or
            (self.pro5.get_first_id() + self.pro5.get_target_decoy()) in self.graph_1.node_to_delete
        )

        self.assertTrue(self.pep3.get_id().sort() == ['3', '7', '9'] or
                        self.pep7.get_id().sort() == ['3', '7', '9'] or
                        self.pep9.get_id().sort() == ['3', '7', '9'])

        self.assertTrue(
            (self.pro3.get_first_id() + self.pro3.get_target_decoy()) in self.graph_1.node_to_delete or
            (self.pro7.get_first_id() + self.pro7.get_target_decoy()) in self.graph_1.node_to_delete or
            (self.pro9.get_first_id() + self.pro9.get_target_decoy()) in self.graph_1.node_to_delete
        )

        # checking all singles stayed single
        self.assertEqual(self.pro1.get_id().sort(), ['1'])
        self.assertEqual(self.pro3.get_id().sort(), ['3'])
        self.assertEqual(self.pro5.get_id().sort(), ['5'])
        self.assertEqual(self.pro6.get_id().sort(), ['7'])
        self.assertEqual(self.pro7.get_id().sort(), ['7'])

        self.assertEqual(self.pep2.get_id().sort(), ['2'])
        self.assertEqual(self.pep4.get_id().sort(), ['4'])
        self.assertEqual(self.pep6.get_id().sort(), ['6'])
        self.assertEqual(self.pep8.get_id().sort(), ['8'])
        self.assertEqual(self.pep10.get_id().sort(), ['10s'])


        # TODO: I need to write a test so that whenever any method go uses
        #  the delete node, it does not work


class TestSeparate(unittest.TestCase):
    """this class contains test for method used in Separate
    """

    def setUp(self) -> None:
        self.graph_1 = graph.Graph()
        self.pro_1 = graph.Protein(['1'], '1', 0)
        self.pro_28 = graph.Protein(['2', '8'], '28', 0)
        self.pro_3 = graph.Protein(['3'], '3', 0)
        self.pro_49 = graph.Protein(['4', '9'], '49', 0)
        self.pro_5 = graph.Protein(['5'], '5', 0)
        self.pro_6 = graph.Protein(['6'], '6', 0)
        self.pro_7 = graph.Protein(['7'], '7', 0)

        self.pep_3_7_9 = graph.Peptide(['3', '7', '9'], 0.1, 0)
        self.pep_4 = graph.Peptide(['4'], 0.2, 0)
        self.pep_8 = graph.Peptide(['8'], 0.3, 0)
        self.pep_2 = graph.Peptide(['2'], 0.4, 0)
        self.pep_6 = graph.Peptide(['6'], 0.5, 0)
        self.pep_10 = graph.Peptide(['10'], 0.6, 0)
        self.pep_1_5 = graph.Peptide(['1', '5'], 0.7, 0)

        self.graph_1.node_dict[self.pro_1] = [self.pep_3_7_9, self.pep_4,
                                              self.pep_8]
        self.graph_1.node_dict[self.pro_28] = [self.pep_8]
        self.graph_1.node_dict[self.pro_3] = [self.pep_6]
        self.graph_1.node_dict[self.pro_49] = [self.pep_2, self.pep_10]
        self.graph_1.node_dict[self.pro_5] = [self.pep_4, self.pep_8]
        self.graph_1.node_dict[self.pro_6] = [self.pep_2, self.pep_6]
        self.graph_1.node_dict[self.pro_7] = [self.pep_1_5]

        self.graph_1.node_dict[self.pep_3_7_9] = [self.pro_1]
        self.graph_1.node_dict[self.pep_4] = [self.pro_1, self.pro_5]
        self.graph_1.node_dict[self.pep_8] = [self.pro_1, self.pro_5,
                                              self.pro_28]
        self.graph_1.node_dict[self.pep_2] = [self.pro_6, self.pro_49]
        self.graph_1.node_dict[self.pep_6] = [self.pro_3, self.pro_6]
        self.graph_1.node_dict[self.pep_10] = [self.pro_49]
        self.graph_1.node_dict[self.pep_1_5] = [self.pro_7]

        self.com_1 = graph.Component()
        self.com_2 = graph.Component()
        self.com_3 = graph.Component()



    def test_DFS(self):
        """test whether dfs put the correct protein and peptide into the
        correct connected components
        """
        # run dfs
        self.graph_1.dfs(self.pro_1, self.com_1)
        self.graph_1.dfs(self.pep_2, self.com_2)
        self.graph_1.dfs(self.pro_7, self.com_3)

        # rename for ease of calling
        self.com_1_pro = self.com_1._protein_dict.keys()
        self.com_1_pep = self.com_1._peptide_dict.keys()
        self.com_2_pro = self.com_2._protein_dict.keys()
        self.com_2_pep = self.com_2._peptide_dict.keys()
        self.com_3_pro = self.com_3._protein_dict.keys()
        self.com_3_pep = self.com_3._peptide_dict.keys()

        self.assertEqual(len(self.com_1_pro), 3)
        self.assertIn(self.pro_1, self.com_1_pro)
        self.assertIn(self.pro_28, self.com_1_pro)  # failed
        self.assertIn(self.pro_5, self.com_1_pro)

        self.assertEqual(len(self.com_1_pep), 3)
        self.assertIn(self.pep_3_7_9, self.com_1_pep)
        self.assertIn(self.pep_4, self.com_1_pep)
        self.assertIn(self.pep_8, self.com_1_pep)

        self.assertEqual(len(self.com_2_pro), 3)
        self.assertIn(self.pro_3, self.com_2_pro)
        self.assertIn(self.pro_49, self.com_2_pro)
        # failed, because forgot between pro49 and pep2
        self.assertIn(self.pro_6, self.com_2_pro)

        self.assertEqual(len(self.com_2_pep), 3)
        self.assertIn(self.pep_2, self.com_2_pep)
        self.assertIn(self.pep_6, self.com_2_pep)
        self.assertIn(self.pep_10, self.com_2_pep)
        # failed, because forgot between pro49 and pep2

        self.assertEqual(len(self.com_3_pro), 1)
        self.assertIn(self.pro_7, self.com_3_pro)

        self.assertEqual(len(self.com_3_pep), 1)
        self.assertIn(self.pep_1_5, self.com_3_pep)


class TestReduce(unittest.TestCase):

    def test_reduce(self):
        component_1 = components.Component()
        component_2 = components.Component()
        component_3 = components.Component()

        pro_1 = graph.Protein(['1'], '1', 0)
        pro_28 = graph.Protein(['2', '8'], '28', 0)
        pro_3 = graph.Protein(['3'], '3', 0)
        pro_49 = graph.Protein(['4', '9'], '49', 0)
        pro_5 = graph.Protein(['5'], '5', 0)
        pro_6 = graph.Protein(['6'], '6', 0)
        pro_7 = graph.Protein(['7'], '7', 0)

        pep_3_7_9 = graph.Peptide(['3', '7', '9'], 0.1, 0)
        pep_4 = graph.Peptide(['4'], 0.2, 0)
        pep_8 = graph.Peptide(['8'], 0.3, 0)
        pep_2 = graph.Peptide(['2'], 0.4, 0)
        pep_6 = graph.Peptide(['6'], 0.5, 0)
        pep_10 = graph.Peptide(['10'], 0.6, 0)
        pep_1_5 = graph.Peptide(['1', '5'], 0.7, 0)

        # not done through dfs, since I want to only test reduce
        component_1.add_protein(pro_1, [pep_3_7_9, pep_4, pep_8])
        component_1.add_protein(pro_28, [pep_8])
        component_1.add_protein(pro_5, [pep_8, pep_4])
        component_1.add_peptide(pep_3_7_9, [pro_1])
        component_1.add_peptide(pep_4, [pro_1, pro_5])
        component_1.add_peptide(pep_8, [pro_1, pro_5, pro_28])

        component_2.add_protein(pro_3, [pep_6])
        component_2.add_protein(pro_49, [pep_2, pep_10])
        component_2.add_protein(pro_6, [pep_2, pep_6])
        component_2.add_peptide(pep_2, [pro_49, pro_6])
        component_2.add_peptide(pep_6, [pro_3, pro_6])
        component_2.add_peptide(pep_10, [pro_49])

        component_3.add_protein(pro_7, [pep_1_5])
        component_3.add_peptide(pep_1_5, [pro_7])

        """a_component 1"""
        self.assertEqual(component_1.find_num_uncovered_peptides(pro_1), 3)
        self.assertEqual(component_1.find_most_uncovered_protein(), pro_1)
        self.assertEqual(component_1.make_protein_list(),
                         [pro_1.get_sqlite_id()])

        # now the peptide neighbours of pro1 should be set as covered
        self.assertTrue(pep_3_7_9.is_covered())
        self.assertTrue(pep_4.is_covered())
        self.assertTrue(pep_8.is_covered())

        """a_component 2"""
        self.assertEqual(component_2.find_num_uncovered_peptides(pro_49), 2)

        self.assertEqual(component_2.find_num_uncovered_peptides(pro_6), 2)

        self.assertTrue(
            component_2.find_most_uncovered_protein() == pro_49 or
            component_2.find_most_uncovered_protein() == pro_6
        )
        # TODO protein 6 returns

        # since pro3 only has 1 peptide neighbor, it should have the same
        # number of uncovered peptide before and after making protein list
        self.assertEqual(component_2.find_num_uncovered_peptides(pro_3), 1)

        self.assertEqual(
            component_2.make_protein_list().sort(),
            [pro_49.get_sqlite_id(), pro_6.get_sqlite_id()].sort()
        )

        self.assertTrue(pep_2.is_covered())
        self.assertTrue(pep_6.is_covered())
        self.assertTrue(pep_10.is_covered())

        """a_component 3"""
        self.assertEqual(component_3.find_num_uncovered_peptides(pro_7), 1)
        self.assertEqual(component_3.find_most_uncovered_protein(), pro_7)

        self.assertEqual(component_3.make_protein_list(),
                         [pro_7.get_sqlite_id()])

        self.assertTrue(pep_1_5.is_covered())


if __name__ == '__main__':
    unittest.main()
