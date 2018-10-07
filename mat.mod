param num_matrices >= 1, integer; # Number of matrices in the data file to be read
param num_rows >= 1, integer;     # Number of rows
param num_cols >= 1, integer;     # Number of columns 

set MATS    := 1 .. num_matrices; # set of matrices
set ROWS    := 1 .. num_rows;	  # set of rows
set COLUMNS := 1 .. num_cols;	  # set of columns

param req_dose >= 0; # values for entries of each matrix
param all_dose >= 0; 
param beam_dose {MATS, ROWS, COLUMNS} >= 0;
param tumor {ROWS, COLUMNS} >= 0;
param crit {ROWS, COLUMNS} >= 0;
param a >= 0; # weight factor

var I {MATS} >= 0; # variable capturing the maximum value at each 
var good >=0;
var bad >=0;				   # index across all matrices given

# Pushing all variables to the maximum value of their corresponding indices
minimize Dosage: sum{k in MATS} I[k];

# Each variable at an index is >= to the maximum value at 
# the index across all matrices given.
subject to Good: good = sum {i in ROWS, j in COLUMNS, k in MATS} I[k]*tumor[i,j]*beam_dose[k,i,j];
subject to Bad: bad = sum {i in ROWS, j in COLUMNS, k in MATS} I[k]*crit[i,j]*beam_dose[k,i,j];
subject to Satisfied {i in ROWS, j in COLUMNS}: req_dose * tumor[i,j] <= sum {k in MATS} I[k]*beam_dose[k,i,j];
subject to Allowed {i in ROWS, j in COLUMNS}: all_dose >= crit[i,j] * sum {k in MATS} I[k]*beam_dose[k,i,j];
