#!/bin/sh
#
#


parmed protein.parm7 <<_EOF
loadRestrt protein.rst7
setOverwrite True
tiMerge :1-394 :395-788 :91 :485
outparm merged_protein.parm7 merged_protein.rst7
quit
_EOF

parmed complex.parm7 <<_EOF
loadRestrt complex.rst7
setOverwrite True
tiMerge :1-394 :395-788 :91 :485
outparm merged_complex.parm7 merged_complex.rst7
quit
_EOF
