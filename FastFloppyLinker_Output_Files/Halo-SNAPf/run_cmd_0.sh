#! /bin/bash
source activate lion
python ../FastFloppyLinker8.py -inpdb Halo-SNAPf_WorkingFusion.pdb -linkers Halo-SNAPf_WorkingFusion.linkers -loops Halo-SNAPf_WorkingFusion.loops -t_frag Halo-SNAPf_WorkingFusion.200.3mers -clusterid $1