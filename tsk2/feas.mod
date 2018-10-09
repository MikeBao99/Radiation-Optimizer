param num_matrices >= 1, integer; # Number of matrices in the data file to be read
param num_rows >= 1, integer;     # Number of rows
param num_cols >= 1, integer;     # Number of columns

set MATS    := 1 .. num_matrices; # set of matrices
set ROWS    := 1 .. num_rows; # set of rows
set COLUMNS := 1 .. num_cols; # set of columns

param req_dose >= 0; # required dosage for tumor cells
param all_dose >= 0; # allowed dosage for critical cells
param beam_dose {MATS, R, C} >= 0, default 0; # beam base dosage
param tumor {ROWS, COLUMNS} >= 0; # tumor matrix
param crit {ROWS, COLUMNS} >= 0; # critical cell matrix

var I {MATS} >= 0; # intensity of beams
var max_buf {ROWS, COLUMNS} >=0; # feasibility buffer
var min_buf {ROWS, COLUMNS} >= 0;	# feasibility buffer

minimize Dosage: sum {i in ROWS, j in COLUMNS} (max_buf[i,j] + min_buf[i,j]);

subject to Satisfied {i in ROWS, j in COLUMNS}: (req_dose - max_buf[i,j]) * tumor[i,j] <= sum {k in MATS} I[k]*beam_dose[k,i,j];
subject to Allowed {i in ROWS, j in COLUMNS}: (all_dose + min_buf[i,j]) >= crit[i,j] * sum {k in MATS} I[k]*beam_dose[k,i,j];
