# createUKBphenome
Rscript that will create a PheWAS code based phenome using ICD9 and ICD10 data from baskets of the UK biobank (http://biobank.ndph.ox.ac.uk/showcase/).  
It uses the PheWAS R package (https://github.com/PheWAS/PheWAS) and its underlying mapping tables in a hopefully transparent way.

## Basic concepts

1. Collect all ICD code / PheWAS code mapping tables (https://phewascatalog.org/phecodes and https://phewascatalog.org/phecodes_icd10)  
2. Collect information about PheWAS codes and their inclusion / exclusion filters  
3. Collect and harmonize all official ICD codes from UKB and map them to phecodes  
4. Extract all ICD codes from the available fields in your UKB baskets  
5. Generate a phenome: case control study for each phecode  


## Step 1: Describe your data
Add all your TAB-delimited baskets in a text file here (one basket per line) --> `./data/baskets.txt`, e.g. ukb####.tab

## Step 2: Create Phenome
`Rscript /net/junglebook/home/larsf/Projects/createUKBphenome/scripts/function.createUKBphenome.r`

## Output
1. Full ICD / PheWAS code tables with descriptions (what's the underlying ICD code for each phecode)
2. UKB phenome with exclusion criteria applied to controls 
3. UKB phenome without applying exclusion criteria to controls
4. Overview of all phecodes, their categories and general descriptions

## Notes:
- This script might require a ton of memory (~20 GB), because it reads and collects a lot of data into memory.
- This script requires ICD data (ideally the most comprehensive list), `genetic sex` and sex data
- Only samples where sex == gender are kept, because it's unclear why sex differs from gender (potential sources for mismatch: gender identity, bone marrow transplant, sample swap)
