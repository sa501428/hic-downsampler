import sys

import numpy as np

from just_extract import write_short_format_line

## python3 split_up.py mnd.txt new_mnd.txt threshold

cutoff = float(sys.argv[3])


def get_new_num_counts(init_num):
    return int(np.sum(np.random.uniform(0, 1, size=(init_num, 1)) < cutoff))


total = 0
outfile = open(sys.argv[2], 'w')

if sys.argv[1].endswith(".bn"):
    with open(sys.argv[1]) as f:
        for line in f:
            line_vals = line.split()
            new_num_counts = get_new_num_counts(int(line_vals[4]))
            if new_num_counts > 0:
                total = total + new_num_counts
                write_short_format_line(outfile, line_vals[0], line_vals[1],
                                        line_vals[2], line_vals[3], new_num_counts)
else:
    with open(sys.argv[1]) as f:
        for line in f:
            line_vals = line.split()
            new_num_counts = get_new_num_counts(int(line_vals[4]))
            if new_num_counts > 0:
                total = total + new_num_counts
                write_short_format_line(outfile, line_vals[0], line_vals[1],
                                        line_vals[2], line_vals[3], new_num_counts)
outfile.close()

print(sys.argv[2], 'will have', total, 'contacts')
