## GENERAL INFORMATION  
This README file was updated on 2026-11-10.  

A. Manuscript associated with this archive:  
Mutz, J. & S. Kalisz. 2026. Multiple herbivore-induced defenses vary in timing but confer sustained resistance in experimental plant populations. Ecology.

B. Author details:  
Jessie Mutz, Dartmouth College (jessie.mutz@dartmouth.edu; corresponding author)  
Susan Kalisz, University of Tennessee, Knoxville

C. Funding sources:  
This work was funded by USDA-NIFA Award 2021-67034-35232.  


## ACCESS INFORMATION 
Licenses/restrictions placed on the code: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication


## OVERVIEW OF REPOSITORY CONTENTS  
This repository consists of raw and cleaned data files, R scripts, output files from statistical analyses, and figures.

### Raw data files ("raw-data/."; all .CSV)

Guide to columns that are shared across many data files (these are also column names for **plant-guide.csv**).
| Column name | Data type | Description |    
| ----------- | ----------- | ----------- |  
| Cage_pos | alphanumeric | position of cage in the field (A-G as "columns", 1-10 as "rows") |
| Treatment | alphanumeric | treatment combination of timepoint (T1-T5), herbivore density treatment (L = low, M = mid, H = high, C = no-herbivore control), and replicate (1-5) |
| Timepoint | numeric | timepoint of leaf sampling (1-5) |
| Herb_dens | character | herbivore density treatment (low, mid, high, control) |
| Rep | numeric | replicate within treatment |
| Plant_ID | numeric | unique plant ID from tag (6001-6700) |
| Geno | numeric | plant genotype |


* **bioassay-data.csv**

Data associated with no-choice bioassays using 2 naive *Spodoptera* caterpillars per plant.
| Column name | Data type | Description |    
| ----------- | ----------- | ----------- |
| Group | alphanumeric | sampling group within timepoint (1A, 1B, 2A, 2B, etc.) |
| Start_date_X | date | start date of no-choice bioassay X (of 2), MM/DD/YYYY |
| Start_time_X | time | start time of no-choice bioassay X (of 2), HH:MM |
| Initial_mass_X | numeric | initial mass of bioassay caterpillar X (of 2) (g) |
| End_date_X | date | end date of no-choice bioassay X (of 2), MM/DD/YYYY |
| End_time_X | time | end time of no-choice bioassay X (of 2), HH:MM |
| Final_mass_X | numeric | final mass of bioassay caterpillar X (of 2) (g) |
| Dead_X | Boolean | 0 = caterpillar alive, 1 = caterpillar dead at end of bioassay X (of 2) |
| Damage_X | categorical | damage (> 0.5%) observed at end of bioassay X (of 2) (yes, some, unsure) |
| Days_elapsed_X | numeric | length of time between start and end of bioassay X (of 2), in days |
| RGR_X | numeric | relative growth rate of bioassay caterpillar X (of 2) |
| Initial_area | numeric | initial leaf area (both leaf halves) in sq. cm |
| Final_damage_X | numeric | leaf area removed by bioassay caterpillar X (of 2) in sq. cm |
| Area_consumed | numeric | total leaf area removed by both bioassay caterpillars  in sq. cm |

* **bioassay-RGR-data.csv**

bioassay-data.csv in long format with separate bioassay caterpillars as separate rows.

* **bioassay-leaves-initial.csv** & **bioassay-leaves-final.csv**

Estimated initial and final leaf areas and leaf area consumed from LeafByte for the leaves used in bioassays. The files have the same structure:
| Column name | Data type | Description |    
| ----------- | ----------- | ----------- |  
| Date_(year.month.day) | date | date of image analysis as YYYY.MM.DD |
| Time | time | time of image analysis as HH:MM:SS |
| Latitude_(degrees) | BLANK | not used |
| Longitude_(degrees) | BLANK | not used |
| Barcode | BLANK | not used |
| Sample_Number | numeric | unique image ID |
| Total_Leaf_Area_(cm2) | numeric | estimated leaf area in sq. cm |
| Consumed_Leaf_Area_(cm2) | numeric | estimated leaf area missing in sq. cm |
| Percent_Consumed| numeric | percentage of leaf area consumed |
| Notes | numeric | unique leaf ID (plant ID only for initial, plant ID-leaf half 1 or 2 for final) | 
| Scale_Length_(cm) | numeric | length of scale used in leaf images |

* **chem-data-formatted.csv**

Hemigossypolone sample masses and relative abundances. Structure is wide-format, where each column designates a sample ([Cage_pos][Plant_ID]). Second row (sample_mass) provides mass of leaf sample analyzed, third row (Hemigossypolone_3.6) provides counts of hemigossypolone with 3.6 minute retention time, fourth row (Hemigossypolone_10.5) provides counts of hemigossypolone with 10.5 minute retention time.

* **field-obs-data.csv**

Data from field observations of plants before *Spodoptera* introduction, at the time of leaf sampling, and at harvest.
| Column name | Data type | Description |    
| ----------- | ----------- | ----------- |  
| Date | date | date of observation as MM/DD/YYYY |
| Census_ID | categorical | timing of data collection (initial = before *Spodoptera* introduction, sampling = at the time of leaf sampling, final = at harvest) |
| Plot_ID | alphanumeric | position of cage within field |
| Obs | character | initials of observer |
| Rec | character | initials of recorder |
| Stem_lengths_raw | alphanumeric | length of main stem and all secondary stems > 4cm (in cm), provided separately |
| Stem_lengths | numeric | length of main stem and all secondary stems > 4cm (in cm), summed |
| Percent_damage_raw | alphanumeric | percent damage estimates for 10 leaves, provided separately |
| Percent_damage | numeric | average percent damage across 10 leaves |
| Bud count | numeric | number of flower buds present |
| Flw_count | numeric | number of open flowers present |
| Fru_count | numeric | number of fruits present |
| Aphids | categorical | number of aphids present (none = none observed, low = 1-50, mid = 50-100, high = 100+) |

* **leaf-trait-midvein-lengths.csv**

Midvein lengths of leaves used to quantify trichome and terpenoid gland densities.
| Column name | Data type | Description |    
| ----------- | ----------- | ----------- |  
| Rep | numeric | leaf half (1 or 2) |
| Plant_ID | numeric | unique plant ID from tag (6001-6700) |
| Midvein_length | numeric | measured length of midvein (cm) |

* **phys-trait-data.csv**

Counts of trichomes and terpenoid glands from leaf images.
| Column name | Data type | Description |    
| ----------- | ----------- | ----------- | 
| Plant_ID | numeric | unique plant ID from tag (6001-6700) |
| Photo_ID | numeric | photo used for image analysis |
| Obs | character | initials of person doing image analysis |
| Folder_name | alphanumeric | name of folder containing leaf image used |
| Leaf_surface | categorical | whether leaf image is of the upper or lower leaf surface |
| Glands | numeric | number of glands counted within specified 2.5mm-radius disc |
| Trichomes | numeric | number of trichomes counted within specified 2.5mm-radius disc |

* **phys-traits-by-observer.csv**

Counts of trichomes and terpenoid glands for a subset of leaves to standardize across observers. Structure same as above, with additional columns of JM_glands and JM_trichomes for JM counts of terpenoid glands and trichomes, respectively.

* **plants-attacked.csv**

Number of plants per cage with ant tending at 24 h after *Spodoptera* introduction.

| Column name | Data type | Description |    
| ----------- | ----------- | ----------- | 
| Cage_pos | alphanumeric | position of cage in the field |
| No_plants_attacked | numeric | number of plants with ant-tending |
| No_plants_observed | numeric | number of plants observed in the cage |

* **terpenoid-gland-areas.csv**
Output of ImageJ macro to estimate area of leaf containing terpenoid glands.

| Column name | Data type | Description |    
| ----------- | ----------- | ----------- | 
| Folder | alphanumeric | name of folder containing leaf image |
| Image_file | alphanumeric | file name of image analyzed |
| Total_area | numeric | estimated area of leaf surface containing terpenoid glands in sq. mm |

* **final-biomass.csv** & **stem-lengths-to-biomass.csv**

Not analyzed or reported in this manuscript.

### Cleaned data files ("cleaned-data/."; all .CSV)
These files contain data formatted and cleaned according to the code in scripts/bioassay-data-setup, scripts/chem-data-setup, scripts/field-obs-data-setup, scripts/full-data-setup, and scripts/phys-trait-data-setup.

### R scripts ("scripts/.", .R and .txt)

Code for cleaning and formatting data: scripts/bioassay-data-setup, scripts/chem-data-setup, scripts/field-obs-data-setup, scripts/full-data-setup, and scripts/phys-trait-data-setup.

Code for statistical analyses:
* Ant tending: **ant-attack.R** and **ant-attack-bin.txt**
* Population mean defense expression and resistance as a function of *Spodoptera* per plant across sampling timepoints: **induced-response-Bayesian-models.R** specifies .txt models (TimexDens-LinThr-Geno.txt, TimexDens-LinThr-Geno-bin.txt), **helper-funcs.R** specifies helper functions for organizing output, **full-analyses-vary-pop-threshold.R** implements analyses.
* KCV analyses: **kcv-models.R** specifies .txt model (NoTimexDens-KCV-beta.txt), **kcv-analyses.R** implements analyses.
* Relationship between defense expression and resistance (plant level): **plant-resistance-defenses.R** specifies .txt model (plant-resistance-defense.txt) and implements analysis.
* Relationship between population mean defense expression and population mean resistance: **pop-resistance-defenses.R** specifies .txt model (pop-resis-NoTimexDens-1.txt) and implements analysis.

Code for figures: figures.R

### Output files ("output/.", .csv, .pdf, and .txt files)
Output files as described in R code for statistical analyses. The files summarize output of statistical analyses and provide plots for determining model convergence.

### Figures ("figures/.", .pdf and .png)
Figures produced for the manuscript, as specified in "scripts/figures".
