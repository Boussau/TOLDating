# Copyright or Copr. Centre National de la Recherche Scientifique (CNRS) (2018)
# Contributors:
# - Vincent Lanore <vincent.lanore@gmail.com>
# - Carine Rey <carine.rey@ens-lyon.org>
# - Bastien Boussau <bastien.boussau@univ-lyon1.fr

# This software is a computer program whose purpose is to provide a set of scripts for pre and post processing of data for
# convergence detection programs.

# This software is governed by the CeCILL-C license under French law and abiding by the rules of distribution of free software.
# You can use, modify and/ or redistribute the software under the terms of the CeCILL-C license as circulated by CEA, CNRS and
# INRIA at the following URL "http://www.cecill.info".

# As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users
# are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive
# licensors have only limited liability.

# In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or
# reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated
# to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth
# computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements
# in conditions enabling the security of their systems and/or data to be ensured and, more generally, to use and operate it in
# the same conditions as regards security.

# The fact that you are presently reading this means that you have had knowledge of the CeCILL-C license and that you accept
# its terms.


import importlib.util
spec = importlib.util.spec_from_file_location("diffsel_script_utils", "/home/boussau/Programming/PhyloDraw/convergence-scripts/script/diffsel_script_utils.py")
utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(utils)
# foo.MyClass()

import pandas as pd
import random

from ete3 import Tree, NodeStyle, TreeStyle, TextFace, faces
# from convergence-scripts/script/diffsel_script_utils import *

#===================================================================================================
utils.STEP("Parsing command line arguments")

from argparse import ArgumentParser, FileType
parser = ArgumentParser(description='Selects tips from a phylogenetic tree.')
parser.add_argument('inputFile', metavar="input", type=FileType('r'), nargs=1, help='the tree file (newick format)')
parser.add_argument('-d', '--dataFile', metavar="dataFile", type=FileType('r'), nargs=1, help="Information about tips")
parser.add_argument('-l', '--level', metavar="level", type=int, nargs=1, help="Level of selection")

#parser.add_argument('-t', '--transition', dest="add_transition", action='store_true', help="add the tag Transition where you put a convergent event.")

args = parser.parse_args()

tree_file = args.inputFile[0]
utils.MESSAGE("Tree file is "+utils.param(tree_file.name))


# sister = args.sister
# utils.MESSAGE("Sister branch condition: "+param(sister))
# add_transition = args.add_transition
# utils.MESSAGE("Add transtion condition: "+param(add_transition))
# out_file = tree_file.name+".annotated"
# utils.MESSAGE("Output file is "+utils.data(out_file))


data_file = args.dataFile[0]
utils.MESSAGE("Data file is "+utils.param(data_file.name)) # bac_metadata_r86.tsv

level = args.level[0]
utils.MESSAGE("Level is " + utils.param(level)) # bac_metadata_r86.tsv


#===================================================================================================

utils.MESSAGE("Reading tree from file")
t = Tree(tree_file.name)

utils.MESSAGE("tree read")


def get_genome_accession(s):
    #print (s.split("%")[0])
    return (s.split("%")[0])

def get_level (i, s):
#    print(s)
#    print(s.split("%")[i])
    return (s.split("%")[i])

def get_level_3 (s):
    return (s.split("%")[3])

group_to_tips = dict()
tip_to_dist = dict()
# Annotate inner nodes
for node in t.traverse("postorder"):
    if node.is_leaf():
#        group = get_level_3 (node.name)
        group = get_level (level, node.name)
        #print(group)
        node.add_feature("group", group)
        node.add_feature("dist_to_root", -1.0)
        if not group_to_tips.__contains__(group):
            group_to_tips[group] = list()
        group_to_tips[group].append(node.name)
    else:
        if (node.children[0].group == node.children[1].group):
            node.add_feature("group", node.children[0].group)
            node.add_feature("dist_to_root", -1.0)
        else:
            node.children[0].add_feature("dist_to_root", 0.0)
            node.children[1].add_feature("dist_to_root", 0.0)
            node.add_feature("dist_to_root",0.0)
            node.add_feature("group","NA")

print("Number of groups: " + str(len(group_to_tips)))

#print(group_to_tips.keys())

for node in t.traverse("preorder"):
    if (not node.is_root()) and node.dist_to_root==-1:
        node.dist_to_root = node.dist + node.up.dist_to_root
    if node.is_leaf():
        tip_to_dist[get_genome_accession(node.name)] = node.dist_to_root

#print (len(tip_to_dist.keys()))

utils.STEP("Getting table")
#df = pd.read_csv('bac_metadata_r86.tsv', sep='\t')
df = open(data_file.name,"r").readlines()

utils.MESSAGE("Table contains "+str(len(df)) + " lines.")

groups={}
accession_to_data = dict()
for line in df[1:]:
    words = line.split("\t")
    group=words[-5].split(";")[level-1]
    if not group in groups.keys():
        groups[group] = {"GC":[],"Size":[],"dist_to_root":[],"completeness":[],"contamination":[],"accession":[]}
    accession = words[0]
    if accession in tip_to_dist.keys():
        gc = float(words[4])
        siz = float(words[6])
        dist = tip_to_dist[accession]
        completeness = float(words[19])
        contamination=float(words[20])
        groups[group]["GC"].append([accession,gc])
        groups[group]["Size"].append([accession,siz])
        groups[group]["dist_to_root"].append([accession,dist])
        groups[group]["completeness"].append([accession,completeness])
        groups[group]["contamination"].append([accession,contamination])
        groups[group]["accession"].append(accession)
        accession_to_data[accession] = {"GC":gc,"Size":siz,"dist_to_root":dist,"completeness":completeness,"contamination":contamination}
# print(len(groups.keys()))
#
# print(list(groups.keys())[0])
#
# test = list(groups.keys())[0]

#print(groups[test]["GC"])

all_GC = []
group_mean_GCs = dict()
for k in groups.keys():
    #print(str(k) + " : "+str(len(groups[k]["GC"])))
    group_GCs = []
    if len(groups[k]["GC"]) > 0:
        for gc in groups[k]["GC"]:
            all_GC.append(gc[1])
            group_GCs.append(gc[1])
        group_mean_GCs[k] = float(sum(group_GCs))/len(group_GCs)
mean_GC = float(sum(all_GC))/len(all_GC)

print ("Mean GC: " + str(mean_GC))
representant = dict()
representantRandom = dict()
representantBestGC = dict()
representantBestCompleteness = dict()
representantBestContamination = dict()
representantBestDistance = dict()
representantBestSize = dict()
representant_noGC = dict()

nsample = 10

for k in groups.keys():
    #print (groups[k])
    score = dict()
    score_noGC = dict()
    tabgc=groups[k]["GC"]
    if (len(tabgc) > 0):
    # If we want to compare to the overall average    tabgc.sort(key=lambda x: abs(mean_GC-x[1]))
        print("Group "+str(k) + " ; Number of Genomes: " + str(len(tabgc)))
        # for i in range(min(len(tabgc), 10)):
        #     print("BEFORE GC : Group "+str(k) + " : " + str(abs(group_mean_GCs[k]-tabgc[i][1])))
        # print("MEAN : " + str(group_mean_GCs[k]))

        tabgc.sort(key=lambda x: abs(group_mean_GCs[k]-x[1]))
        for t in range(len(tabgc)):
            score[tabgc[t][0]] = t
        tabsize=groups[k]["Size"]
        tabsize.sort(key=lambda x: -x[1])
        for t in range(len(tabsize)):
            score[tabsize[t][0]] += t
            score_noGC[tabsize[t][0]] = t
        tabcomplet=groups[k]["completeness"]
        # for i in range(min(len(tabgc), 10)):
        #     print("BEFORE INCOMPLETENESS : Group "+str(k) + " : " + str(tabcomplet[i][1]))
        tabcomplet.sort(key=lambda x: -x[1])
        # for i in range(min(len(tabgc), 10)):
        #     print("AFTER INCOMPLETENESS : Group "+str(k) + " : " + str(tabcomplet[i][1]))
        for t in range(len(tabcomplet)):
            score[tabcomplet[t][0]] += t
            score_noGC[tabsize[t][0]] += t
            # score[tabcomplet[t][0]] += t
        tabconta=groups[k]["contamination"]
        tabconta.sort(key=lambda x: x[1])
        for t in range(len(tabconta)):
            score[tabconta[t][0]] += t
            score_noGC[tabsize[t][0]] += t
        tabdist=groups[k]["dist_to_root"]
        tabdist.sort(key=lambda x: x[1])
        for t in range(len(tabdist)):
            score[tabdist[t][0]] += t
            score_noGC[tabsize[t][0]] += t
        print ("\tBest GC: " + str(tabgc[0][1]) + "; Worst GC: " +  str(tabgc[-1][1]))
        print ("\tBest size: " + str(tabsize[0][1]) + "; Worst size: " +  str(tabsize[-1][1]))
        print ("\tBest completeness: " + str(tabcomplet[0][1]) + "; Worst completeness: " +  str(tabcomplet[-1][1]))
        print ("\tBest contamination: " + str(tabconta[0][1]) + "; Worst contamination: " +  str(tabconta[-1][1]))
        print ("\tBest distance: " + str(tabdist[0][1]) + "; Worst distance: " +  str(tabdist[-1][1]))
        keys = list(score.keys())
        keys.sort(key=lambda x: score[x])
        keys_noGC = list(score_noGC.keys())
        keys_noGC.sort(key=lambda x: score_noGC[x])
        if len(keys) > 0:
            representant[k] = keys[0]
            representant_noGC[k] = keys_noGC[0]
            representantRandom[k] = [random.choice(keys) for _ in range(nsample)]
            representantBestGC[k] = tabgc[0][0]
            representantBestCompleteness[k] = tabcomplet[0][0]
            representantBestContamination[k] = tabconta[0][0]
            representantBestDistance[k] = tabdist[0][0]
            representantBestSize[k] = tabsize[0][0]
            print("\tSelected representant genome: GC: " + str(accession_to_data[keys[0]]["GC"]) + " ; size: " + str(accession_to_data[keys[0]]["Size"]) + " ; completeness: " + str(accession_to_data[keys[0]]["completeness"]) + " ; contamination: " + str(accession_to_data[keys[0]]["contamination"]) +" ; distance:" + str(accession_to_data[keys[0]]["dist_to_root"]))
output = open("representant_"+str(level),"w")
for k in representant.keys():
    output.write(representant[k]+"\n")
output.close()

output = open("representantNoGC_"+str(level),"w")
for k in representant_noGC.keys():
    output.write(representant_noGC[k]+"\n")
output.close()

output = open("representantRandom_"+str(level),"w")
for k in representant.keys():
    for i in range(nsample):
        output.write(representantRandom[k][i]+"\n")
output.close()

output = open("representantBestSize_"+str(level),"w")
for k in representant.keys():
    output.write(representantBestSize[k]+"\n")
output.close()

output = open("representantBestGC_"+str(level),"w")
for k in representant.keys():
    output.write(representantBestGC[k]+"\n")
output.close()

output = open("representantBestCompleteness_"+str(level),"w")
for k in representant.keys():
    output.write(representantBestCompleteness[k]+"\n")
output.close()

output = open("representantBestContamination_"+str(level),"w")
for k in representant.keys():
    output.write(representantBestContamination[k]+"\n")
output.close()

output = open("representantBestDistance_"+str(level),"w")
for k in representant.keys():
    output.write(representantBestDistance[k]+"\n")
output.close()

#print (len(representant.keys()))



#
# condi_color_dic = {"0":"#E6E6FA", "1":"#ff0000", "2":"#90EE90"}  ##ADD8E6
#
# utils.MESSAGE("Setting node style")
# nstyle = NodeStyle()
# nstyle["fgcolor"] = "black"
# nstyle["size"] = 1
#
# utils.MESSAGE("Setting tree style")
# tree_style = TreeStyle()
# tree_style.show_leaf_name = False
# tree_style.show_branch_length = False
# tree_style.draw_guiding_lines = True
# tree_style.complete_branch_lines_when_necessary = True
# tree_style.legend_position = 1
# tree_style.mode = "c"
# tree_style.arc_start = -180 # 0 degrees = 3 o'clock
# tree_style.arc_span = 180
#
#
# # utils.MESSAGE("Setting legend with condition numbers and colors")
# # for condi_i in sorted(condi_color_dic.keys()):
# #     tf = TextFace("Condition      " + condi_i)
# #     tf.background.color = condi_color_dic[condi_i]
# #     tf.margin_right = 2
# #     tf.margin_top = 1
# #     tf.margin_left = 2
# #     tf.margin_bottom = 1
# #     tf.border.width = 1
# #     tree_style.legend.add_face(tf, column=1)
#
# # if add_transition:
# #     utils.MESSAGE("Setting transition style")
# #     tf = TextFace("Transition -> x")
# #     tf.background.color = "white"
# #     tf.margin_right = 2
# #     tf.margin_top = 1
# #     tf.margin_left = 2
# #     tf.margin_bottom = 1
# #     tf.border.color = "red"
# #     tf.border.width = 2
# #     tree_style.legend.add_face(tf, column=1)
#
# #===================================================================================================
# utils.STEP("Tree retrieval and preparation")
#
# def set_tag(node, tag, value):
#     if hasattr(node, tag):
#         setattr(node, tag, value)
#     else:
#         node.add_feature(tag, value)
#
# def set_if_no_tag(node, tag, value):
#     if not hasattr(node, tag):
#         node.add_feature(tag, value)
#
# utils.MESSAGE("Detect existing tags")
# features = []
# for n in t.traverse("postorder"): # get all features:
#     features.extend(list(set(dir(n)) - set(dir(Tree()))))
# features = list(set(features)) # list(set(*)) = remove duplicates
#
# utils.SUBMESSAGE("No detected tag" if not features else "Detected tags: "+", ".join([utils.data(f) for f in features]))
#
# if "i" in features:
#     utils.WARNING("\"i\" is in the detected tags but it will be removed by the programm")
#     features.remove("i")
#
# if not "Condition" in features:
#     utils.MESSAGE("Setting all nodes to Condition = "+utils.data(0))
#     features.append("Condition")
# else:
#     utils.MESSAGE("Setting all nodes without tag Condition to "+utils.data(0))
# for n in t.traverse("postorder"):
#     set_if_no_tag(n, "Condition", 0)
#
# #if add_transition and not "Transition" in features:
# #    features.append("Transition")
#
# utils.MESSAGE("Numberings nodes")
# i = 0
# for n in t.traverse("postorder"):
#     set_tag(n, "i", str(i))
#     i+=1
#
# #===================================================================================================
# utils.STEP("Convergent branch selection")
#
#
# def remove_redundancy(duplicate):
#     final_list = []
#     for num in duplicate:
#         if num not in final_list:
#             final_list.append(num)
#     return final_list
#
# def draw_tree(tree, sample_names):
#     tree_copy = tree.copy("newick-extended")
#     tree_copy.add_feature("i", tree.i)
#     tree_copy.add_feature("Condition", tree.Condition)
#
#     ancestor = tree_copy.get_common_ancestor(sample_names)
#     paths = []
#     for sa in sample_names:
#         node = tree_copy.get_leaves_by_name(name=sa)[0]
#         while node.up:
#             if node == ancestor:
#                 break
#             paths.append(node)
#             node = node.up
#     pruned_path = remove_redundancy(paths)
#
#     for n in pruned_path:
#         mark_node(n, 1)
#
# #    nameFace = faces.AttrFace("name", fgcolor=condi_color_dic[str(n.Condition)])  #fsize=20,
#
#     for n in tree_copy.traverse():
#         if n.is_leaf():
#                 n.set_style(nstyle)
#                 n.add_face(TextFace(str(n.name), fgcolor=condi_color_dic[str(n.Condition)]), column=0, position="aligned")
#         else:
#             n.set_style(nstyle)
#         nd = TextFace(str(n.i))
#         nd.background.color = condi_color_dic[str(n.Condition)]
#         nd.margin_right = 2
#         nd.margin_top = 1
#         nd.margin_left = 2
#         nd.margin_bottom = 1
#         nd.border.width = 1
#
#         style = NodeStyle()
#         style["fgcolor"] = "#0f0f0f"
#         style["size"] = 0
#         #style["vt_line_color"] = condi_color_dic[str(n.Condition)]  #"#ff0000"
#         style["hz_line_color"] = condi_color_dic[str(n.Condition)]  #"#ff0000"
#         style["vt_line_width"] = 8
#         style["hz_line_width"] = 8
#         style["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
#         style["hz_line_type"] = 0
#         n.img_style = style
#         #n.add_face(nameFace, n, column=0)
#         # if add_transition:
#         #     if hasattr(n, "Transition"):
#         #         nd.border.color = "red"
#         #         nd.border.width = 2
#         n.add_face(nd, column=0, position="float")
#         n.add_face(TextFace("       "), column=0, position="branch-bottom")
#
#
#     return tree_copy
#
# def mark_node(node, condition):
#     utils.SUBMESSAGE("adding tag Condition = " + utils.data(condition) + " to node " + utils.data(node.i))
#     node.Condition = condition
#
#
# def mark_subtree(node, condition):
#     utils.SUBMESSAGE("adding tag Condition = " + utils.data(condition) + " to the subtree rooted at node " + utils.data(node.i))
#     node.Condition = condition
#     for child in node.get_descendants():
#         child.Condition = condition
#
#     # if add_transition:
#     #     utils.SUBMESSAGE("adding tag Transition = " + utils.data(condition) + " at node " + utils.data(node.i))
#     #     set_tag(node, "Transition", condition)
#
# draw_tree(t, sample_names)
#
# utils.MESSAGE("Starting subtree selection")
# pdf_file = tree_file.name+".pdf"
# while True:
#     # asking user input
#     #utils.MESSAGE("Please look at "+utils.data(pdf_file)+" to see node numbers")
#     #draw_tree(t).render(pdf_file, tree_style=tree_style)
#     utils.MESSAGE("Please look at tree to see node numbers")
#     draw_tree(t, sample_names).show(tree_style=tree_style)
#
#     user_input = input(utils.ask_input("Please enter start of convergent subtree (type "+utils.green("s")+" to save and quit):"))
#
#     #processing input
#     if user_input.isdigit(): #if input is an int
#         n_i = t.search_nodes(i=user_input) # locating subtree whose root is at nb
#         if n_i: # if found
#             n_i = n_i[0] # we expect only one result as indices are supposed to be unique
#
#             mark_subtree(n_i, 1)
#
#             # if sister: # if sister is specified, handle sister trees
#             #     n_s = n_i.get_sisters()
#             #     for n_s_i in n_s:
#             #         mark_subtree(n_s_i, 2)
#
#         draw_tree(t, sample_names)
#
#     elif user_input == "s":
#         break
#     elif user_input == "q":
#         print(utils.MESSAGE("Exit"))
#         exit(1)
#     else:
#         print(utils.failure("Input was not an integer; try again"))
#
# #===================================================================================================
# utils.STEP("Writing result to file: ")
# utils.MESSAGE("-- Output file is " + utils.data(out_file))
# t.write(format=1, features=features, outfile = out_file)
