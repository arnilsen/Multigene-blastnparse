# Multigene-blastnparse

multigene_blastnparse.py is a simple script for retrieving associated sequences from blastn matches to a query sequence. This script is intended to be useful for multigene fungal phylogenies, but has also been successfully tested on small plant and animal datasets. The script is written in python3.x and makes extensive use of biopython. It takes a fasta file of user query sequences and returns a csv file of blast matches and all available markers for specimens. Additionally, a fasta file containing all sequences for the different makers is created.

**System requirements**

multigene_blastnparse.py requires that python 3.x and biopython is installed. Installation instructions for python are available [here](https://www.python.org/downloads/), and biopython [here](https://biopython.org/wiki/Download).

Use git clone to obtain the script, or copy and paste the script into a text editor.

**Running the script**
```python
python multigene_blastnparse.py -i <user supplied fasta file> -e <email required by NCBI>
```
Compulsory arguments are -i, and -e. The email address is a requirement from NCBI so you can be contacted in case of excessive usage. The user can also specify the optional argument -v, --verbose. This will retain all temporary files. The script can take sometime to run, depending on the load at NCBI.

**Output**

multigene_blastnparse.py outputs two files, a csv and a fasta file. The csv file lists species-vouchers with additional available markers. If no identifying voucher can be ascertained, then the species is listed with unknown appended to it. It also includes blast match information for each specimen: the pairwise identity (pairwise), length of the alignment (align_length), and which query sequence the match is to, numbered in the order that the sequences are in the user provided fasta file (query). Several conveniences have been integrated into the csv file: identification of type specimens, identification of joint ITS and LSU sequences, and concatenation of gene names. The script can identify if a matching sequence is from a type specimen, appending ‘-TYPE’ to the marker. Additionally, if the sequence is from the ITS region and is longer than 700 bp long, the sequence is labelled as ‘ITS-LSU?’. Sequences that cannot be identified are binned as unknown. Sequences that contain multiple genes have their unique marker names concatenated e.g. trnL-trnF.
The fasta file contains all the sequences in the csv file with the same species-voucher identifier in the first position of the header. The second position contains the accession and third the marker name. This allows the user to use tools like grep to retrieve all makers from the generated fasta file i.e.

```
grep -A 1 ‘>Cortinarius_acidophilus__O:T.-E.-Brandrud-61-79’ <fasta file> | sed '/^--$/d' > <file name.fasta>
```

**Known limitations**

There are two main limitations of this script: inconstancies in species names/specimen vouchers in GenBank entries and whole genomes. Occasionally GenBank entries contain discrepancies in their species names and specimen vouchers in different markers of the same individual. This causes the script to create separate rows for non-identical species and vouchers. These differences can mostly be detected by the user as the csv is ordered alphabetically and the different entries are normally grouped together.
The other limitation stems from blast matches to genomes, in particular plastid genomes. Blast matches to genomes results in the genomes being retrieved from NCBI. If the genome is annotated with multiple genes then all the unique gene names are concatenated in the csv and fasta file. The separate regions are not extracted. I intend to add functionality that will extract gene regions from genomes in the future.
