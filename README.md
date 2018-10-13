# SMSD
Small Molecule Subgraph Detector (SMSD) is a Java based software library for finding Maximum Common Subgraph (MCS)/ Substructure between small molecules.
This enables help us to find similarity/distance between two molecules. MCS is also used for screening drug like compounds by hitting molecules, which share common subgraph (substructure).


The SMSD code:

The present code is part of the Small Molecule Subgraph Detector 
(SMSD http://www.ebi.ac.uk/thornton-srv/software/SMSD ) library.

Encourage this project by <b>citing the paper and the (SMSD) URL</b> 

S. A. Rahman, M. Bashton, G. L. Holliday, R. Schrader and J. M. Thornton, Small Molecule Subgraph Detector (SMSD) toolkit, Journal of Cheminformatics 2009, 1:12. DOI:10.1186/1758-2946-1-12

Wish you a happy coding!

-----------------------
Maven POM configuration
-----------------------

```xml

<dependencies>
    <dependency>
        <groupId>uk.ac.ebi.smsd</groupId>
        <artifactId>smsd-core</artifactId>
        <version>LATEST</version>
    </dependency>
    <dependency>
        <groupId>uk.ac.ebi.smsd</groupId>
        <artifactId>smsd-exec</artifactId>
        <version>LATEST</version>
    </dependency>
</dependencies>
```
--------------------
Compile (compiled jar might be available under exec/target/ folder)
--------------------
```
-compile core modules
mvn clean install

-compile with dependencies
mvn install
```

--------------------
Command Line Options
--------------------
```
java -Xms500M -Xmx512M -cp smsd-2.2.0.jar: uk.ac.ebi.smsd.cmd.SMSDcmd
```
===================

Windows Platform (DOS): SMSD.bat file
```
java -Xms500M -Xmx512M -cp smsd-2.2.0.jar: uk.ac.ebi.smsd.cmd.SMSDcmd %*
```

Unix/Mac: SMSD.sh file
```
java -Xms500M -Xmx512M -cp smsd-2.2.0.jar: uk.ac.ebi.smsd.cmd.SMSDcmd $@
```
===================

--------------------
Options
--------------------

```
Single SMILES query vs single SMILES target:

java -Xms500M -Xmx512M -cp smsd.jar: uk.ac.ebi.smsd.cmd.SMSDcmd -Q SMI -q "CCN" -T SMI -t "CCCNC"

Single SMILES query vs single MOL target:

java -Xms500M -Xmx512M -cp smsd.jar: uk.ac.ebi.smsd.cmd.SMSDcmd -Q SMI -q "CCN" -T MOL -t Data/ATP.mol

Signature query vs smiles target, outputting the subgraph to stdout:

java -Xms500M -Xmx512M -cp smsd.jar: uk.ac.ebi.smsd.cmd.SMSDcmd -Q SIG -q "[C]([C][C])" -T SMI -t "C(C)CC" -O SMI -o -

Multiway N-MCS, outputing the subgraph as a smiles to stdout:

java -Xms500M -Xmx512M -cp smsd.jar: uk.ac.ebi.smsd.cmd.SMSDcmd -T SDF -t Data/arom.sdf -N -O SMI -o -

Just use ./SMSD -I to list all the image options.

Few more options:

++++++++++++++++++++++++++++++++++++++++++++++
NOTE: The graph matching is performed by removing
the Hydrogens
++++++++++++++++++++++++++++++++++++++++++++++

usage:
       
 -A                  Appends output to existing files, else creates new
                     files
 -a                  Add Hydrogen
 -b                  Match Bond types (Single, Double etc)
 -d <WIDTHxHEIGHT>   Dimension of the image in pixels
 -f <number>         Default: 0, Stereo: 1, Stereo+Fragment: 2,
                     Stereo+Fragment+Energy: 3
 -g                  create png of the mapping
 -h,--help           Help page for command usage
 -I <option=value>   Image options
 -m                  Report all Mappings
 -N                  Do N-way MCS on the target SD file
 -o <filename>       Output the substructure to a file
 -O <type>           Output type
 -Q <type>           Query type (MOL, SMI, etc)
 -q <filepath>       Query filename
 -r                  Remove Hydrogen
 -s                  SubStructure detection
 -S <suffix>         Add suffix to the files
 -T <type>           Target type (MOL, SMI, etc)
 -t <filepath>       Target filename
 -z                  Ring Match

++++++++++++++++++++++++++++++++++++++++++++++

Allowed types for single-molecules (query or target):
MOL	MDL V2000 format
ML2	MOL2 Tripos format
PDB	Protein Databank Format
CML	Chemical Markup Language
SMI	SMILES string format
SIG	Signature string format

Allowed types for multiple-molecules (targets only):
SDF	SD file format
```

THIRD PARTY TOOL
---------------
You need the CDK (https://github.com/cdk/cdk) as SMSD depends on the CDK for processing the Chemical information.
