Meld2Nexus
====================
Lightweight program to merge alignment files, preserving partition information in Nexus CHARSET format. Assumes one gene per file (although this seems to be a reasonable extension to pursue). Does not assume the same taxon sampling across files; however, for taxa present across files, make sure names do not vary (or they will be interpreted as distinct). Currently takes in (vanilla) Nexus-formatted alignments as input, although fasta support is on the way. As always, line returns are assumed to be of unix flavour.

Compile
---------------

In a  terminal prompt in the src directory, type:

	make

Usage
---------------

Call as either:

	./Meld2Nexus -c file1.NEX file2.NEX file3.NEX ... [-o outFileName]

or

	./Meld2Nexus -f fileContainingNexusFileNames [-o outFileName]

where 'outFileName' is optional (default = "Merged.Nex").

UPDATE
--------------
This is now implemented in the [phyx](https://github.com/FePhyFoFum/phyx) program `pxcat`.

