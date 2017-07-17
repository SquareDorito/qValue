@echo off
title Testing R XML Parsing

set R_path="C:\Users\knoh1\Documents\qValue_product\qValue.r"
set pXML_path="C:\Users\knoh1\Documents\qValue_product\phospho_data.xml"
set npXML_path="C:\Users\knoh1\Documents\qValue_product\unphospho_data.xml"
set wait_path="C:\Users\knoh1\Documents\qValue_product\wait.txt"

Rscript %R_path% %pXML_path% %npXML_path% %wait_path% timepoints=12 replicates=5 label_free=2 silac_numerator=3 silac_denominator=2 min_replicates=3 unpaired=1 min min_peak_area=1000

pause