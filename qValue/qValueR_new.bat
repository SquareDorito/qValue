@echo off
title Testing R XML Parsing

set R_path="C:\Users\knoh1\Documents\qValue_lf_test\qValue_slave3.r"
set pXML_path="C:\Users\knoh1\Documents\qValue_lf_test\phospho.xml"
set npXML_path="C:\Users\knoh1\Documents\qValue_lf_test\nonphospho.xml"
set counter_path="C:\Users\knoh1\Documents\qValue_lf_test\counters.txt"
set wait_path="C:\Users\knoh1\Documents\qValue_lf_test\wait.txt"

Rscript %R_path% %pXML_path% %npXML_path% %counter_path% %wait_path% timepoints=1 replicates=5 label_free=3 silac_numerator=2 silac_denominator=3 min_replicates=3 unpaired=1 min min_peak_area=0

pause