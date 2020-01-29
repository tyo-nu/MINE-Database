__author__ = 'JGJeffryes'

from pymongo import MongoClient
import sys
import platform
import time
from optparse import OptionParser
from minedatabase.queries import quick_search
from minedatabase.databases import establish_db_client

"""
def write_pathway_html(db, path, outfile):
    html_printer = Printer(db)
    for _id in path:
        if _id[0] == 'C':
            html_printer.print_compound_html(_id, outfile, labels=['Names', 'KEGG_code'])

        if _id[0] == 'R':
            html_printer.print_rxn_html(_id, outfile)"""


class PathwaySearch():
    def __init__(self, options):
        client = establish_db_client()
        self.paths = []
        self.db = client[options.db]
        self.options = options
        self.native_comps = set()
        self.discovered = set()

    def __compound_tests(self, compound):
        if not compound:
                if self.options.verbose:
                    print "Could not find " + compound['_id']
                return False
        try:
            if float(compound['NP_likeness']) <= self.options.np_min:
                if self.options.verbose:
                    print "%s excluded for insufficient natural product likeness" % compound['_id']
                return False
        except KeyError:
            pass
        return True

    def __reaction_checks(self, reaction):
        if not reaction:
            if options.verbose:
                print "Could not find " + reaction['_id']
            return False
        try:
            if float(reaction['Energy']) >= self.options.gibbs_cap:
                if options.verbose:
                    print "%s excluded for thermodynamics" % reaction['_id']
                return False
        except (KeyError, ValueError):
            pass
        return True

    def dfs(self):
        self._DFS_ittr(self.options.start_comp, [])
        return sorted(self.paths, key=lambda x: len(x))

    def _DFS_ittr(self, comp_id, path):
        if comp_id == self.options.end_comp:
            self.paths.append(path+[comp_id])
            if not self.options.all_paths:
                return self.paths
        elif len(path)/2 < self.options.len_limit:
            compound = self.db.compounds.find_one({"_id": comp_id}, {'NP_likeness': 1, 'Reactant_in': 1})
            if self.__compound_tests(compound):
                try:
                    for rxn_id in compound["Reactant_in"]:
                        reaction = self.db.reactions.find_one({"_id": rxn_id}, {'Products': 1, 'Energy': 1})
                        if self.__reaction_checks(reaction):
                            for prod_tuple in reaction['Products']:
                                if prod_tuple[1][0] == "C":
                                    self._DFS_ittr(prod_tuple[1], path + [comp_id, rxn_id])
                except KeyError:
                    if self.options.verbose:
                        print "Dead End"
        else:
            if self.options.verbose:
                print "Max Depth"

    def bfs(self):
        self.__load_queue(self.options.start_comp)
        while len(self.queue) > 0:
            #remove a compound from the queue and add add to discovered. if the path length is at the limit we break
            compound = self.queue.pop(0)
            #the depth is half the path length because we only count reactions
            depth = len(compound['Path']) / 2
            if self.options.verbose:
                print "Depth: %s\tQueue:%s" % (depth, len(self.queue))
            if compound['_id'] == self.options.end_comp:
                #save the path
                self.paths.append(compound['Path'] + [compound['_id']])
                #continue or exit depending on the multi flag
                if not self.options.all_paths:
                    return self.paths
            if depth < self.options.len_limit:
                #look at the products for every reaction in which the start compound is a reactant
                try:
                    for rxn_id in compound["Reactant_in"]:
                        for product_id in self.__get_products(rxn_id):
                            product = self.db.compounds.find_one({"_id": product_id}, {'NP_likeness': 1, 'Reactant_in': 1})
                            #if the product passes, record the path taken to reach the new node and enqueue
                            if self.__compound_tests(product):
                                product['Path'] = compound['Path'] + [compound['_id'], rxn_id]
                                self.queue.append(product)
                except KeyError:
                    if self.options.verbose:
                        print "Dead End"

        return self.paths

    def __load_queue(self, start_id):
        self.queue = []
        #look at the products for every reaction in which the start compound is a reactant and add all the non-cofactors
        # to the queue.
        self.discovered.add(start_id)
        start_dict = self.db.compounds.find_one({"_id": start_id})
        for rxn_id in start_dict["Reactant_in"]:
            reaction = self.db.reactions.find_one({"_id": rxn_id})
            for comp_tuple in reaction["Products"]:
                #check if compound is a cofactor by changing the first part of its _id
                if comp_tuple[1][0] == "C":
                    #the second part of the tuple has the id, thus the 1 subscript
                    compound = self.db.compounds.find_one({"_id": str(comp_tuple[1])})
                    #Start a path field to keep track of all the steps that brought us to the compound
                    compound['Path'] = [start_id, rxn_id]
                    self.discovered.add(compound['_id'])
                    self.queue.append(compound)


    def __get_products(self, rxn_id):
        products = []
        reaction = self.db.reactions.find_one({"_id": rxn_id}, {'Products': 1, 'Energy': 1})
        #reaction passes tests, return all products that haven't been discovered and are not cofactors
        if self.__reaction_checks(reaction):
            for comp_tuple in reaction["Products"]:
                if (comp_tuple[1] not in self.discovered) and (comp_tuple[1][0] == "C"):
                    self.discovered.add(comp_tuple[1])
                    products.append(comp_tuple[1])
        return products


if __name__ == '__main__':
    start_time = time.time()
    #This block handles user flags and arguments. For more information see the optparse API documentation

    usage = "usage: %prog [options] run_name"
    parser = OptionParser(usage)
    parser.add_option("-d", "--database", dest="db", default="1GenEcoCyc",
                      help="The name of the database to search")
    parser.add_option("-s", "--start", dest="start_comp", default="Cb5b3273ab083d77ed29fbef8f7e464929af29c13",
                      help="Specify the starting compound for a search, defaults to glucose")
    parser.add_option("-e", "--end", dest="end_comp", default="C89b394fd02e5e5e60ae1e167780ea7ab3276288e",
                      help="Specify the ending compound for a search, defaults to ethanol")
    parser.add_option("-l", "--length", dest="len_limit", type="int", default=3,
                      help="Maximum length of pathways to return, defaults to 3 reactions")
    parser.add_option("-n", "--NP", dest="np_min", type="float", default=-3,
                      help="Set a floor on the minimum natural product likeness of compounds in a pathway")
    parser.add_option("-g", "--GibbsCap", dest="gibbs_cap", type="float", default=1000,
                      help="Set a cap on the gibbs free energy of a reaction in a pathway")
    parser.add_option("-a", "--all_paths", dest="all_paths", action="store_true", default=False,
                      help="Find all valid paths shorter than length cap")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False,
                      help="Print verbose output")
    (options, args) = parser.parse_args()

    client = establish_db_client()
    db = client[options.db]

    options.start_comp = quick_search(db, options.start_comp)[0]['_id']
    options.end_comp = quick_search(db, options.end_comp)[0]['_id']

    search = PathwaySearch(options)
    if options.all_paths:
        paths = search.dfs()
    else:
        paths = search.bfs()

    for i, path in enumerate(paths):
        """if args:
            with open(args[0]+'%s.html' % i, 'w') as outfile:
                write_pathway_html(db, path, outfile)
        else:"""
        print(', '.join(str(x) for x in path))

    print 'Found %s paths shorter than %s reactions' % (len(paths), options.len_limit)
    print 'Script complete, execution took %s sec' %(time.time() - start_time)