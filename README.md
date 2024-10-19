# Harper_et_al_2024_sampling_approaches_great_crested_newt

Data processing workflow and supplementary data for:

Harper *et al.* (2024) Extending sampling approaches for great crested newt (*Triturus cristatus*) eDNA monitoring. *bioRxiv*


## Background

Paired water samples for ethanol precipitation and filtration were collected from ponds once a month from April to October to examine whether the great crested newt (*Triturus cristatus*) eDNA survey window could be extended beyond 30 June, and whether filtration could be used for great crested newt eDNA capture. 

The water sampling protocol set forth by Biggs *et al.* (2014) was used for both eDNA capture methods with minor modifications for filtration, where 20 x 125 mL subsamples were collected at equidistant intervals around the pond perimeter and pooled into a single sampling bag for homogenisation, following which as much water as possible was filtered. 

We aimed to include 17 positive and 3 negative ponds for great crested newt each month, and field blanks (mineral water) were included from May onwards. Conventional population estimates and Habitat Suitability Index assessments for each pond were also performed for comparison to eDNA results. 

The resulting data were analysed to examine:

1. The number of positive ponds produced each month by ethanol precipitation and filtration
2. The eDNA score (number of positive qPCR replicates) for each pond in each month produced by ethanol precipitation and filtration.
3. The number of positive ponds in-season (April to June) vs. out-of-season (July to October) produced by ethanol precipitation and filtration.
4. The eDNA score for each pond in-season (April to June) vs. out-of-season (July to October) produced by ethanol precipitation and filtration
5. Whether volume of water filtered influenced the eDNA scores for ponds (and thus detection).
6. Whether there was a relationship between population size estimates and eDNA scores.


## Contents

R scripts used to analyse metabarcoding data and produce figures [(here)](https://github.com/lrharper1/Harper_et_al_2024_sampling_approaches_great_crested_newt/tree/master/R%20scripts)

Data needed to run analyses in R [(here)](https://github.com/lrharper1/Harper_et_al_2024_sampling_approaches_great_crested_newt/tree/master/Data/)


## Setting up the environment

In order to retrieve scripts and associated data, start by cloning this repository to your current directory:

```
git clone --recursive https://github.com/lrharper1/Harper_et_al_2024_sampling_approaches_great_crested_newt.git
```

In order to make use of our scripts and data, you will have to install R on your computer. R is compatible with all major operating systems, but see CRAN (https://cran.r-project.org/) for details.

Once R is installed, open our R script in the R console and simply run the code!
