#Copyright (c) 2022-2024, INRAE - MIAT
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#
#The python library pytoulbar2 is provided by the toulbar2 team at
#https://github.com/toulbar2/toulbar2
#


import argparse
import itertools as it
import pytoulbar2
import numpy as np
import sys
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import copy



pytoulbar2.tb2.option.verbose = 0

def write_containts(ploidy) :
    print("Writing contraint file")
    ctm = {} # dict mapping contigs to sets of transcripts
    cnm = {} # dict mapping contig names to indices
    ind = 0
    top = 100000000
    mycfn = pytoulbar2.CFN(top)
    
    with open("prot.links") as f:
        for l in f:
            l = l.split()
            #print(l)
            if (len(l) != 2):
                print("Malformed line")
                exit(1)
            ts = l[1].split(',')
            ctm[l[0]] = set()
            cnm[l[0]] = ind
            for t in ts:
                ctm[l[0]].add(t)
            ind += 1
            mycfn.AddVariable(l[0],[i for i in range(ploidy)])

    shift = 0
    for c1,c2 in it.combinations(ctm.keys(),2):
        inter = len(ctm[c1].intersection(ctm[c2]))
        if (inter > 0):
            #print(cnm[c1],cnm[c2],inter)
            shift += inter
            table = np.identity(ploidy)*inter # maxcut: pay inter if different. Minimization, pay -inter if different or inter if equal  
            mycfn.AddFunction([c1,c2],list(table.flatten()))

    # this block has been transformed in double loop 
    #mycfn.AddFunction([0],[0,top,top,top])
    #mycfn.AddFunction([1],[0,0,top,top])
    #mycfn.AddFunction([2],[0,0,0,top])

    l1 = [0 for i in range(ploidy)]
    for i in range(ploidy - 1) :
        l2 = copy.deepcopy(l1)
        #print("l2=",l2)
        for j in range(i+1,ploidy) :
            l2[j]=top
        mycfn.AddFunction([i],l2)
    #print([i],l2)
    
    #print("Cost shift =",shift)
    mycfn.Dump("constraints.cfn")

def align_proteins(proteins, genome, ploidy, threads):
    try:
        print("protein alignments started...this can take some time")
        command = "miniprot -N"+str(ploidy -1)+" -t"+str(threads)+" "+genome+" "+proteins+" > prot_align.paf"
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if result.returncode == 0 :
            print("proteins aligned ")
        else:
            print("proteins not aligned")
    except FileNotFoundError:
        print("protein alignment file could not be created")
        return False


    
def run_contraints(time):
    print("Solving constraints, this will run "+str(time)+" seconds")
    try:
        command = "xz -fz constraints.cfn"
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        command = "toulbar2 -vns constraints.cfn.xz -s=3 -w=solution -timer="+str(time)
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if result.returncode == 0 :
            print("constraint solved")
        else:
            print("constraint not solved")
    except Exception as e:
        print("contraint solving did not succeed")
        return False
      
def generate_link_file():
    haplo = []
    try:
        print("generate_link_file")
        d = {}
        f = open("prot_align.paf", "r")
        for l in f :
            ll = l.split("\t")
            if ll[5] not in d :
                d[ll[5]]=ll[0]
            else :
                d[ll[5]] = d[ll[5]] + "," + ll[0]
        f.close()
        fo = open("prot.links", "w")
        for k,v in d.items() :
            haplo.append(k)
            fo.write(k+"\t"+v+"\n")
        fo.close()
    except FileNotFoundError:
        print("protein alignment file could not be created")
        return False 
    return haplo
    
def get_haplotypes(names) :
    try :
        with open("solution") as f:
            hap = {}
            #print(names)
            for i,l in enumerate(f):
                if i == 0 :
                    ll = l.split(" ")
                    print("number of contigs to class = "+str(len(names)))
                    print("solution size = "+str(len(ll)))
                    for j,lll in enumerate(ll) :
            	        hap[names[j]] = int(lll)+1 
    except FileNotFoundError:
        print("contraint solution file not found")
        return False
    return(hap)             
                      
# Function to check if all specified software can be run on the system
def check_softwares(softwares):
    for software in softwares:
        #print(software)
        try:
            # Run the command with no arguments
            result = subprocess.run([software], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            # print(result.returncode)       
            # Check the return code
            if result.returncode == 0 or result.returncode == 1 :
                print("The "+software+" command is accessible.")
            else:
                print("The "+software+" command is not accessible.")
                sys.exit(1) 
                return False
        except FileNotFoundError:
            print("The "+software+" command is not found in the environment.")
            return False
        except Exception as e:
            print(f"An error occurred: {e}")
            return False

class SequenceManipulator:
    def __init__(self):
        self.sequences = {}  # Sequence dict 
        
    def load(self, file_path):
        """load sequences from FASTA file"""
        print("loading : "+file_path)
        try:
            records = list(SeqIO.parse(file_path, "fasta"))
            for record in records:
                self.sequences[record.id] = record
                #print(record.id)
            print(f"{len(records)} sequences loaded from {file_path}")
        except FileNotFoundError:
            print("File not found.")

    def print_haplotypes(self, ploidy, prefix, haplotype_dict):
        """create and fill haplotype files"""
        #print(ploidy, prefix, haplotype_dict)
        print("Writing haplotype files")
        f = [open(prefix+str(i)+".fasta",'w') for i in range(1, ploidy+1)]
 
        for k,v in haplotype_dict.items() :
            #print(k, v)
            f[v-1].write(">"+k+"\n")
            f[v-1].write(str(self.sequences[k].seq)+"\n")
            
        for i in range(1, ploidy + 1) :
            f[i - 1].close()
            
def main():
    parser = argparse.ArgumentParser(description="split haplotypes wiht toulbar2.")
    parser.add_argument("--assembly", type=str, required=True, help="input assembly fasta file")
    parser.add_argument("--proteins", type=str, required=True, help="input proteins fasta file")
    parser.add_argument("--ploidy", type=int, required=True, help="input genome ploïdy")
    parser.add_argument("--output", type=str, default="hap_", required=False, help="haplotype prefix")
    parser.add_argument("--mpthreads", type=int, default=4, required=False, help="threads to run miniprot")
    parser.add_argument("--optime", type=int, default=900, required=False, help="optimisation time in seconds")

    args = parser.parse_args()

    # check software 
    softwares = ("xz -h","miniprot -h","toulbar2 -h")
    check_softwares(softwares)

    # create sequence manipulator
    manipulator = SequenceManipulator()

    # load initial fasta file 
    manipulator.load(args.assembly)

    # compare protein file to genome file 
    align_proteins(args.proteins, args.assembly, args.ploidy, args.mpthreads)
    
    # generate link file 
    linked_contig_names = generate_link_file()
    print("Number of contigs having protein links :"+str(len(linked_contig_names)))
    
    # generate constraint file 
    write_containts(args.ploidy)
    
    # run constraint file with toulbar2
    run_contraints(args.optime)

    # run constraint file with toulbar2
    haplo_dict = get_haplotypes(linked_contig_names)
    
    # Sauvegarder les résultats dans le fichier de sortie
    manipulator.print_haplotypes(args.ploidy, args.output, haplo_dict)

if __name__ == "__main__":
    main()
