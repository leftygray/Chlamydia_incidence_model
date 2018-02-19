## Australian Chlamydia Incidence Model

This project contains the code for a model to estimate trends in chlamydia incidence in Australia. The code implements a Bayesian statistical method based on a decision-pathway model producing estimates for chlamydia incidence in Australia since 2001 for male and female 16–29-year-olds. The original code was developed for the paper: 

* Hammad Ali, Ewan Cameron, Christopher C Drovandi, James M McCaw, Rebecca J Guy, Melanie Middleton, Carol El-Hayek, et al. “A New Approach to Estimating Trends in Chlamydia Incidence.” Sexually Transmitted Infections, 2015;91:513–519. doi:10.1136/sextrans-2014-051631

For the paper, estimates were generated for the 2001-2013 period. The original code has been modified slightly to produce annual estimates since 2013 for the HIV, viral hepatitis and sexually transmissible infections in [Australia Annual Surveillance Report](https://kirby.unsw.edu.au/report-type/annual-surveillance-reports) (ASR) produced annually by the [_The Kirby Institute_](https://kirby.unsw.edu.au/). The model is written in R (currently version 3.3.2) with inputs to produce estimates for the period 2001-2016.

**Original Model Developer:** Ewan Cameron     
**Email:** dr.ewan.cameron@gmail.com  
**Affiliations:** School of Mathematical Sciences, Queensland University of Technology, Brisbane, Australia; Spatial Ecology & Epidemiology Group, University of Oxford, Oxford, United Kingdom

**Repository owner:** Richard T. Gray  
**ORCID ID:** orcid.org/0000-0002-2885-0483    
**Affiliation:** [_The Kirby Institute_](https://kirby.unsw.edu.au/), UNSW Sydney, Sydney NSW 2052, Australia

### Aims ###

Directly measuring disease incidence in a population is difficult and not feasible to do  routinely. The aim of this project was to develop and apply a new method for estimating at a population level the number of incident genital chlamydia infections, and the corresponding incidence rates, by age and sex using routine surveillance data. 

### Code organization ###

All the model scripts, input data and outputs are stored in the main directory and 3 sub-directories. The main directory also contains an optional Rstudio project file `chlamydia_model.Rproj` (this is not required for running the model). The user is required to generate the input files in the appropriate format and specify the specific outputs to process. 

_Main directory files_

The following three R scripts run the model using the the `.csv` input files from the `data/` directory and generate the chlamydia incidence results and figures which are stored in the `output\` directory. Instructions for running the code are provided in the next section. 

- **abc.run.abc.R:** Contains the main algorithm calling other routines in `code\` to read in the data files in `data\` and runs the ABC fitting routine. Results are stored in the `output\` directory in a folder with a name corresponding to the date and time the script was sourced started. Output files are generated for each simulation with name `theta.test.N.dat` (where N is the simulation number). SMC-ABC control parameters are specified on lines 29-33. 
  
- **abc.posterior.processing.R:** Generates the posterior estimates of chlamydia notifications and test counts, incidence counts and proportions, positivity, and prevalence. The user must specify the date and time directory produced by `abc.run.abc.R` and the simulation number for the final output file. The output is saved in the same directory as `posterior.dat`. 

- **abc.plot.results.R:** Produces the final results and PDF plots corresponding to the results in `posterior.dat` in a user specified directory in `output\`. All results are stored in the `output/figures/` directory. The incidence estimates for each age and sex are stored in a separate `Chlamydia_Incidence.csv` file. WARNING: sourcing this file will overwrite previously generated output files. 

#### code ####

Contains all the specific R functions used in the analysis and by the main directory scripts.

#### data ####

Contains all the input data files used by the `abc.run.abc.R` script. These files are stored as `.csv` files containing data since 2001 to the year of analysis (currently 2016). Four input files are required for the model to run:

- `priors.csv`:  Specifications for prior distributions describing constant parameters for all sub-populations.
- `population.csv`: Estimated population size for males, females and overall for each age from 0 to 100 years for each year since 2001.  
- `notifications.csv`: Annual notifications for males and females aged < 15, 15-19, 20-24, 25-29, 30-33, and > 34 for each year since 2001. 
- `tests.csv`:  Annual number of chlamydia tests conducted for males and females aged < 15, 15-19, 20-24, 25-29, 30-33, and > 34 for each year since 2001. NOTE: data for 2006-2007 is missing for Australia (this is specifically accounted for in the scripts). 

#### output ####

Contains output directories storing the results produced by the main directory scripts and associated figures. Each run/source of `abc.run.abc.R` produces a folder with the date and time of running. These output directories need to be managed by the user. 

### Instructions for running the model ###

Running the model and its input/output requires some management by the user. Requiring some familiarity with the R language and the use of CSV files. The code is not so sophisticated that it can take care of everything automatically and the user should perform sanity checks along the way to test whether the contents of any given file have actually been read in or not.

To run the model perform the following steps:

1. Update the input files in the `data/` directory to contain required data. These files must have the same cloumn format but the number of rows can vary (this number must be consistent across the files though). 
2. Source `abc.run.abc.R`. This runs the  SMC-ABC algorithm. The fits will eventually stop by themselves once the code judges any further improvements to be outweighed by the time required in computing, but the job may also be killed manually from outside if required. Nothing is to be gained by running the fitting code for longer than 12 hours. Simulation output will be stored in a created directory in `output/` named with the date and time. The output files have the fomart `theta.test.N.dat` (where N is the simulation number).    
3. Edit lines 6-7 `abc.run.abc.R` to have the correct output date and time and simulation number (final run of `abc.run.abc.R` saved) and source. This produces the posterior estimates for chlamydia notifications and test counts, incidence counts and proportions, positivity, and prevalence and produces a `posteriors.dat` output file in associated output directory. This script should run until it stops (probably taking an hour). 
4. Source `abc.plot.results.R` to generate the plots and incidence results. 

Please contact Richard Gray (Rgray@kirby.unsw.edu.au) or Ewan Cameron (dr.ewan.cameron@gmail.com) if you have trouble running or deciphering the code. 
  
