# MC-Annotate

**Description**
-----------
MC-Annotate is an annotation tool for PDB files. 
It generates information about nucleotides conformation and their spatial interactions.


**Dependency**
--------------
- MC-Core libraries
- CMake
- C++ compiler
- courage

**Documentation**
-----------------
To get used to the concepts, please consult the references.
The actual code documentation is on its way.

**Installation**
----------------
1. mkdir build && cd build
2. ccmake .. (and the rest of the configuration stuff)
3. make (inside the same folder)
- will probably go fubar, so iterate on step 2 and 3 until it works
- check out your environment variables with "env"

**Authors**
-----------
The general concepts were developped by Sebastien Lemieux and Patrick Gendron (see References).
The software was mainly developed by Martin Larose and Philippe Thibault.

**License**
-----------
Read LICENSE.txt

**References**
--------------
Gendron P, Lemieux S, Major F (2001) Quantitative analysis of nucleic acid three-dimensional structures. J Mol Biol 308:919–936.

Lemieux S, Major F (2002) RNA canonical and noncanonical base pairing types: A recognition method and complete repertoire. Nucleic Acids Res 30:4250–4263.



