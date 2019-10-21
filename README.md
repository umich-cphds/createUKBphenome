# createUKBphenome

## Basic concepts
1. ICD code / PheWAS code mapping from phewascatalog (https://phewascatalog.org/phecodes and https://phewascatalog.org/phecodes_icd10)  
2. Collection of information about PheWAS codes and their inclusion / exclusion filters  
3. Collection and harmonization of ICD codes from UKB
4. Extraction of all ICD codes from the available fields in your UKB baskets  
5. Generatation of a phenome: case control study for each phecode

## Required R libraries
- data.table
- tidyr
- parallel
- intervals
- XML
- RCurl
- rlist
- bitops

## Step 1: Describe your data
Add all your TAB-delimited baskets in a text file here (one basket per line) --> `./data/baskets.txt`, e.g. ukb####.tab

## Step 2: Create Phenome
`cd createUKBphenome`  
`Rscript ./scripts/function.createUKBphenome.r`

## Output
1. Full ICD / PheWAS code tables with descriptions (what's the underlying ICD code for each phecode)
2. UKB phenome with exclusion criteria applied to controls 
3. UKB phenome without applying exclusion criteria to controls
4. Overview of all phecodes, their categories and general descriptions
5. Output of all ICD codes that were NOT mapped to phecodes (incl. sample sizes)
6. Output of all individuals that had sex-specific diagnose codes that did not match their sex

## Notes:
- This script requires a ton of memory (~20-30 GB), because it reads and collects a lot of data into memory.
- This script requires ICD data of the UK Biobank (ideally the most comprehensive list), `Genetic Sex` and `Sex`
- Only samples with `Genetic Sex` equals `Sex` are kept, because it's unclear why it should be different (potential sources for mismatch: gender identity, bone marrow transplant, sample swap)
