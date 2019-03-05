
python ../../Scripts/RenameLeavesAccordingToRefFile.py  arc122_r86.1.tree arc122_r86.1.treeRN arc_metadata_r86.tsv

python ../../Scripts/removeQuotesFromTree.py arc122_r86.1.tree arc122_r86.1.treeNoQuotes


# Selection at level 1:
python ../../Scripts/select_leaves.py arc122_r86.1.treeRN -d arc_metadata_r86.tsv -l 1 > outputSelectionLevel1
# Selection at level 2:
python ../../Scripts/select_leaves.py arc122_r86.1.treeRN -d arc_metadata_r86.tsv -l 2 > outputSelectionLevel2
# Selection at level 3:
python ../../Scripts/select_leaves.py arc122_r86.1.treeRN -d arc_metadata_r86.tsv -l 3 > outputSelectionLevel3


# extract a subtree with sampled genomes
python ~/Documents/Utils/scripts/pruneTreeWithListOfLeaves.py arc122_r86.1.treeRN representant_1 sampled_level1.dnd
python ~/Documents/Utils/scripts/pruneTreeWithListOfLeaves.py arc122_r86.1.treeRN representant_2 sampled_level2.dnd
python ~/Documents/Utils/scripts/pruneTreeWithListOfLeaves.py arc122_r86.1.treeRN representant_3 sampled_level3.dnd

