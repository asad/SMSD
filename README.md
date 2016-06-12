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
Command Line Options
--------------------
```
java -Xms500M -Xmx512M -cp bin/smsd.jar: uk.ac.ebi.smsd.cmd.SMSDcmd
```
===================
Windows Platform (DOS)
```
java -Xms500M -Xmx512M -cp bin/smsd.jar: uk.ac.ebi.smsd.cmd.SMSDcmd %*
```

Unix/Mac
```
java -Xms500M -Xmx512M -cp bin/smsd.jar: uk.ac.ebi.smsd.cmd.SMSDcmd $@
```
===================

THIRD PARTY TOOL
---------------
You need the CDK (https://github.com/cdk/cdk) as SMSD depends on the CDK for processing the Chemical information.
