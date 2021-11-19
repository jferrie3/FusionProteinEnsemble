#! /bin/bash
source activate lion
python ../FastFloppyLinker8.py -inpdb Halo-7xIg-SNAPf_WorkingFusion.pdb -linkers Halo-7xIg-SNAPf_WorkingFusion.linkers -loops Halo-7xIg-SNAPf_WorkingFusion.loops -t_frag Halo-7xIg-SNAPf_WorkingFusion.200.3mers -clusterid $1