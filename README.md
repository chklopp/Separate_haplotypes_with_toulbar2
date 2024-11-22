# Separate haplotypes with toulbar2

Given a genome assembly in fasta format and a set of proteins from the same or a closely related species also in fasta format and the awaited ploitdy, our python script will align the proteins onto the genome, produce a protein to contig link file, transform the link file in a contraint file, resolve the contraints with toolbar2 and extract the haplotyped contigs lists in separate haplotype fasta files.

Dependencies 
- miniprot (https://github.com/lh3/miniprot)
- toulbar2 (https://github.com/toulbar2/toulbar2)
- xz

Running the script 

<pre>
python shtb2.py --proteins proteins.fasta --assembly assembly.fasta --ploidy 4
    
</pre>

Script help 

<pre>
python shtb2.py --help
usage: shtb2.py [-h] --assembly ASSEMBLY --proteins PROTEINS --ploidy PLOIDY [--output OUTPUT] [--mpthreads MPTHREADS] [--optime OPTIME]

split haplotypes wiht toulbar2.

options:
  -h, --help            show this help message and exit
  --assembly ASSEMBLY   input assembly fasta file
  --proteins PROTEINS   input proteins fasta file
  --ploidy PLOIDY       input genome ploïdy
  --output OUTPUT       haplotype prefix
  --mpthreads MPTHREADS threads to run miniprot
  --optime OPTIME       optimisation time in seconds
</pre>
