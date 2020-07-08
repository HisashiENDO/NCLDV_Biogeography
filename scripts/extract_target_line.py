import multiprocessing
import itertools
import sys
import os

# sepicify query ID path (Specific family or total NCLDVs generated with parse_jplace.py)
target_path = sys.argv[1] # Output of the script "parse_jplace"

# input_path and output_path
input_path="file_of_abundance_profile" # Original file in https://www.ocean-microbiome.org
output_path="output.table" #Specify output file

with open(target_path,"r") as target_file:
    target_set=set(line.strip() for line in target_file)

def Is_target(line):
    return line.split("\t")[0] in target_set

def Split_input_file(input_path, size=10000000):
    input_file_b = open(input_path, "rb")
    file_size = os.path.getsize(input_path)
    pos_p = 0
    while pos_p + size < file_size:
        input_file_b.seek(pos_p + size)
        line = input_file_b.readline()
        pos = input_file_b.tell()
        yield (pos_p, pos)
        pos_p = pos
    yield (pos_p, file_size)

def Gen_res(pos_tup):
    pos_p, pos = pos_tup
    input_file = open(input_path, "r")
    input_file.seek(pos_p)
    result_line_arr=[line for line in input_file.readlines(pos - pos_p) if Is_target(line)]
    return result_line_arr

file_segment_gen=Split_input_file(input_path)
pool = multiprocessing.Pool(8)
res_gen = pool.imap_unordered(Gen_res, file_segment_gen)
with open(output_path,"w") as output_file:
    output_file.write("".join(itertools.chain.from_iterable(res_gen)))
