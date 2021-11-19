import sys
from optparse import OptionParser

def initialize_options() :
    parser = OptionParser();
    parser.add_option( "-p", "--pdb", dest="pdbname", type="string", help="The input pdb that should be renumbered" );
    parser.add_option( "-o", "--output", dest="output", type="string", help="The output file to write the renumbered pdb to; stdout if not given", default="" );
    parser.add_option( "-a", "--chainA_begin", dest="chAbegin", type="int", help="The first residue in chain A", default=1 );
    parser.add_option( "-b", "--chainB_begin", dest="chBbegin", type="int", help="The first residue in chain B", default=1 );
    parser.add_option( "-c", "--chainC_begin", dest="chCbegin", type="int", help="The first residue in chain C", default=1 );
    parser.add_option( "-d", "--chainD_begin", dest="chDbegin", type="int", help="The first residue in chain D", default=1 );
    parser.add_option( "-e", "--chainE_begin", dest="chEbegin", type="int", help="The first residue in chain E", default=1 );
    return parser

def remap_resnum_for_line( line, mapping ):
    if line[0:4] != "ATOM": return
    chain = line[21];
    if chain not in mapping: return line
    resstring = line[22:27]
    resnum, lastresstring = mapping[ chain ]
    if lastresstring == "":
    	stopidx = list(mapping.keys()).index(chain)
    	for chainkey in list(mapping.keys())[:stopidx]:
    		chainresnum, chainlastresstring = mapping[ chainkey ]
    		if chainresnum > 1:
    			mapping[ chain ] = (chainresnum + 1, resstring)
    			resnum, lastresstring = mapping[ chain ]
    if lastresstring == "" or resstring != lastresstring :
        if lastresstring != "" : resnum += 1
        mapping[ chain ] = (resnum, resstring )
    newresstring = str(resnum) + " "
    if len(newresstring) == 2: newresstring = "   " + newresstring
    elif len(newresstring) == 3: newresstring = "  " + newresstring
    elif len(newresstring) == 4: newresstring = " " + newresstring
    return line[0:22] + newresstring + line[27:]

if __name__ == "__main__":

    parser = initialize_options()
    (options, args) = parser.parse_args()

    pdblines = open( options.pdbname ).readlines()
    mapping = {}
    mapping["A"] = (options.chAbegin, "");
    mapping["B"] = (options.chBbegin, "");
    mapping["C"] = (options.chCbegin, "");
    mapping["D"] = (options.chDbegin, "");
    mapping["E"] = (options.chEbegin, "");
    newlines = []
    for line in pdblines:
    	lineout = remap_resnum_for_line( line, mapping )
    	if lineout != None:
    		newlines.append( lineout )
    if options.output == "" :
        for line in newlines :
            print(line,)
    else :
   		open( options.output, "w").writelines(newlines)
