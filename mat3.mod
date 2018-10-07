param num_matrices >= 1, integer; # Number of matrices in the data file to be read
param num_rows >= 1, integer;     # Number of rows
param num_cols >= 1, integer;     # Number of columns 

set MATS    := 1 .. num_matrices; # set of matrices
set R := 0 .. (num_rows+1);	  # set of rows
set C := 0 .. (num_cols+1);	  # set of columns
set ROWS    := 1 .. num_rows;
set COLUMNS := 1 .. num_cols;

param req_dose >= 0; # values for entries of each matrix
param all_dose >= 0; 
param beam_dose {MATS, R, C} >= 0, default 0;
param tumor {ROWS, COLUMNS} >= 0;
param crit {ROWS, COLUMNS} >= 0;
param max_buf {ROWS, COLUMNS} >= 0;
param min_buf {ROWS, COLUMNS} >= 0;
param nabla;

var I {MATS} >= 0; # variable capturing the maximum value at each 
				   # index across all matrices given

# Pushing all variables to the maximum value of their corresponding indices
minimize Dosage: (sum {i in ROWS, j in COLUMNS, k in MATS, i2 in (i-1)..(i+1), j2 in (j-1)..(j+1)} crit[i,j]*I[k]*beam_dose[k,i2,j2]);

# Each variable at an index is >= to the maximum value at 
# the index across all matrices given.

subject to Satisfied {i in ROWS, j in COLUMNS}: (req_dose - max_buf[i,j]) * tumor[i,j] <= sum {k in MATS} I[k]*beam_dose[k,i,j];
subject to Allowed {i in ROWS, j in COLUMNS}: (all_dose + min_buf[i,j]) >= crit[i,j] * (sum {k in MATS} I[k]*beam_dose[k,i,j]);
