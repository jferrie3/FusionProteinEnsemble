#! /bin/bash
source activate lion
python ../FastFloppyLinker8.py -inpdb Halo-4xIg-SNAPf_WorkingFusion.pdb -linkers Halo-4xIg-SNAPf_WorkingFusion.linkers -loops Halo-4xIg-SNAPf_WorkingFusion.loops -t_frag Halo-4xIg-SNAPf_WorkingFusion.200.3mers -clusterid $1