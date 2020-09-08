call activate timstof3
python timstof/src/extract_msn_nopd.py --skip-ms1 %1 --output-path %2
#python calculate_xics  -i %1 
call conda deactivate
