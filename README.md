
## totallySAM.py 1.0.0 
  ### SAM file analysis tool               

<img align=right src="https://img.shields.io/pypi/pyversions/3?color=caa6f7&logo=Python&logoColor=white&style=flat-square">

#### Author :
Anaïs Prud'homme (prudhomme.anais.12@outlook.fr) <br>
M1 SNS - BCD (Bioinformatics) <br>
HMIN113M - Systèmes <br>
WOOHP coorporation <br>
  
#### Date :
 - creation date : 12-nov-2019
 - first version : 09-jan-2020
  
#### Requierements :
To accomplish its mission, totallySAM requires Python3. <br>
You can install Python3 by running : `sudo apt-get install python3`
    
#### Description :
totallySAM.py is a bioinformatic tool for SAM file analysis. It runs only for paired-end sequences SAM file. <br>
The standard programm will return :
 - Short file description
 - number of perfectly mapped reads, number of partially mapped reads and number of unmapped reads
 - number of mapped/unmapped pairs, number of mapped/mapped reads, number of partially_mapped/mapped reads and the number of partially_mapped/unmapped reads
    
You can use it with an option to get, for partially mapped reads, the number of substitutions. <br>
The script has been coded under Python3.
 
#### If you have a mission for totallySAM.py : How to run totallySAM.py ?
 - Download totallySAM.py program on the GitHub
 -  Open your terminal
 -  Go to the folder where you saved it
 -  Use the command line : <br>
    Standard version
    ` Python3 ./totallySAM.py samfile.sam ` <br>
    With variant calling option
    ` Python3 ./totallySAM.py samfile.sam -v `
    
 The program need SAM format file to run.
	
#### Results :
The program print SAM file information :
 - File Name
 - File Description (Format Version, Reference sequence name, Reference sequence length, Program Identifier, Read group identifier)
 - Sequences Description (Number of header lines, Number of reads, Number of corrupted reads (ignored because it doesn't have the right column number, some information are missing))
			
The program print analysis results : 
 - number of reads which are : perfectly mapped (without substitution), partially mapped (one or more nucléotides aren't mapped), unmapped (none nucleotide is mapped)
 - pairs :
	 - one read is partially mapped, the mate is unmapped
	 - one read is perfectly mapped, the mate is unmapped
	 -  one read is partially mapped, the mate is perfectly mapped
	 - both reads are unmapped
			
If the user run optionnal version, the program wil also print, for each read, the read name and the number of substitution (the nucleotide differs from reference).
	
    
#### Issues : 
 - The program doesn't compute flag score so it doesn't take into account uncommon flag (it only use common flag that you can find on https://www.samformat.info/sam-format-flag)
 - The program doesn't export results in a file 
  
#### Commits :
 - 23-nov-2019 : the program only return pairs and reads information
 - 17-dec-2019 : the program check some file information to specify if the file is corrupted
 -  8-jan-2019 : the program can now (optionnaly) count and return variant number per reads
  
#### Licence :

[<img align="right" src="https://img.shields.io/badge/License-GNU GPLv3-caa6f7?style=flat-square">](https://www.gnu.org/licenses/gpl-3.0)

The program is available on : https://github.com/NaisPrudhomme/totallySAM
