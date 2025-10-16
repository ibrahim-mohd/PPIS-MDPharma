from modeller import *
from modeller.automodel import *
import os
import argparse
# Define help message
help_message = """\
First, write an alignment file, i.e., find out which residues are missing by reading the PDB and FASTA sequences.

Modify the FASTA sequence from the PDB site as follows: use "/" to separate different chains and "*" to mark the end.

Example:  
>P1;4oc7
sequence:hmx:::::::0.00:0.00
GSHMTSSANEGSHMTSSANEDMPVERILEAELAVEPKTETYVEANMGLNPSSPNDPVTNICQAADKQLFTLVEWAKRIPHFSELPLDDQVILLRAGWNELLIASFSHRSIAVKDGILLATGLHVHRNSAHSAGVGAIFDRVLTELVSKMRDMQMDKTELGCLRAIVLFNPDSKGLSNPAEVEALREKVYASLEAYCKHKYPEQPGRFAKLLLRLPALRSIGLKCLEHLFFFKLIGDTPIDTFLMEMLEAPHQMT/
KHKILHRLLQDSS*
"""

# Argument parser
parser = argparse.ArgumentParser(
    description=help_message, formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument("-p", "--path", type=str, default=os.getcwd(), help="Path to find all files")
parser.add_argument("-f", "--pdb_file", type=str, default="test.pdb", help="PDB file")
parser.add_argument("-s", "--fasta_file", type=str, default="modeller.fasta", help="Modeller FASTA file")

# Parse arguments
args = parser.parse_args()


# input
path    = args.path
pdb_file  = args.pdb_file
fasta_file = args.fasta_file
## protein code
f = open  (path + f"/{fasta_file}", "r+").readlines()
code_ = f[0].split (";")[1][:-1]
#print ("\n\n",code_, "Here0", "\n\n")
###################################################


#
os.chdir(path)  

env = Environ()
aln = Alignment(env)
env.rand_seed = -1
#mdl = Model(env, file='./6i2i/6i2i_pro.pdb')
mdl = Model(env, file=path +f'/{pdb_file}')
aln.append_model(mdl, align_codes='initial_structure', atom_files=path +f'/{pdb_file}')
#aln.append(file=path+f'/{fasta_file}', align_codes='xxxx')
#aln.salign()
#aln.write(file=path+'/alignment.ali', alignment_format='PIR')
 #
env = Environ()
aln = Alignment(env)
 
mdl = Model(env, file=path +f'/{pdb_file}')
aln.append_model(mdl, align_codes='initial_structure', atom_files=path +f'/{pdb_file}')
aln.append(file=path+f'/{fasta_file}', align_codes=code_)
aln.salign()
aln.write(file=path+'/alignment.ali', alignment_format='PIR')



a = DOPELoopModel(env, alnfile = path+'/alignment.ali',knowns = 'initial_structure', sequence = code_)

a.starting_model= 1 # index of the first model
a.ending_model = 1 # index of the last model
# (determines how many models to calculate)
a.md_level = None # No refinement of model
a.loop.starting_model = 1 # First loop model
a.loop.ending_model = 4 # Last loop model
a.loop.md_level =  refine.very_slow #refine.fast # Loop model refinement level

a.make()

