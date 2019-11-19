# Protocol

Preparing marker genes for analysis (Dated TOL project)

Use markers from Williams et al. (2019) [51]. 
Collect from all of the chosen genomes in a per-domain fashion.
Run phylter 0.9, dump the incongruent gene/species pairs (per-family).

Align:

ls *fa | parallel 'mafft --localpair --maxiterate 1000 {} > {.}_linsi.aln'

Mask:

ls *aln | parallel 'java -jar ~/bin/BMGE-1.12/BMGE.jar -t AA -i {} -m BLOSUM30 -of {.}_bmge.fas'

Concat:

~/bin/catfasta2phyml/catfasta2phyml.pl -f -c *fas > concat/euk_fams_linsi_bmge_concat.fas

Initial tree:

iqtree -nt 10 -s euk_fams_linsi_bmge_concat.fas -m C60+F -bb 1000 &


