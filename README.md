# Separate haplotypes with toulbar2

Given a genome assembly in fasta format and a set of proteins from the same or a closely related species also in fasta format and the awaited ploitdy, our python script will align the proteins onto the genome, produce a protein to contig link file, transform the link file in a contraint file, resolve the contraints with toolbar2 and extract the haplotyped contigs lists in separate haplotype fasta files. One can also add a second set of constraints which indicates for each contig in which haplotype it should be placed. The protein based constraint weight is square the number of links. The group constraint weight is the number of proteins in the contig.

The preprint is available on biorxiv [Improving hifiasm haplotypes for autopolyploid genome assemblies using constraint programming](https://www.biorxiv.org/content/10.1101/2025.04.01.646355v1.article-metrics)
[doi : https://doi.org/10.1101/2025.04.01.646355](https://doi.org/10.1101/2025.04.01.646355)

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
usage: shtb2.py [-h] --assembly ASSEMBLY --proteins PROTEINS --ploidy PLOIDY [--output OUTPUT] [--mpthreads MPTHREADS] [--skip_align SKIP_ALIGN] [--optime OPTIME] [--groups GROUPS]
                    [--groups_weight GROUPS_WEIGHT]

split haplotypes wiht toulbar2.

options:
  -h, --help            show this help message and exit
  --assembly ASSEMBLY   input assembly fasta file
  --proteins PROTEINS   input proteins fasta file
  --ploidy PLOIDY       input genome ploïdy
  --output OUTPUT       haplotype prefix
  --mpthreads MPTHREADS
                        threads to run miniprot
  --skip_align SKIP_ALIGN
                        if the alignment is already performed
  --optime OPTIME       optimisation time in seconds
  --groups GROUPS       groups produced in a previous step (two columns, one for contigs and one for group, start at 0

</pre>
