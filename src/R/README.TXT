On the command prompt, run
Rscript LLS_runv6.R L dt N ns T nosim test_name dir_name mu
where 
nosim=100 (100 runs)
test_name = simple_complex (for example)
dir_name = directory name. 
For example, the following code runs the medium simple/complex cases:
* Rscript LLS_runv6.R 50 10 1000 100 310 100 "medium_complex" "../wfsim/data/" 0.0001
* Rscript LLS_runv6.R 50 10 1000 100 1000 100 "medium_simple" "../wfsim/data/" 0.0001

The output is the CSV file containing the runtimes and selection coefficients.


