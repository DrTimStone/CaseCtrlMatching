# CaseCtrlMatching

This script takes all the available cancers in our most recent study (either aedinocarcinoma or intramucosal dysplasia samples) and matches them with each of our different classes of control (people with no positive diagnosis, healthy volunteers, high-grade barrett's or high-grade dysplasia).

The source of the information is a database download, which has not only the clinical status of the individuals, but also various lifestyle factors. These include age, sex and bmi, but also data about alcohol consumption, smoking, PPI consumption (PPIs are a class of heartburn drug) and heartburn information.

The lifestyle and clinical information is processed. Estimates are made of cigarette pack-years and analagous measures for alcohol, PPIs and heartburn. 

After this has happened, ALL these main covariates are normalised to have the same upper bound of 100. This is necessary because otherwise, the cigarette vector would have different length (and dominiate greatly) if, for example, we changed from packets of cigarettes to individual cigarette consumption.

These normalised vectors and then multiplied by a weighting factor, it is more important that sexes and ages are correctly matched than lifestyle ones, the greater weighting of age and sex gives them a matching priority.

The matching itself is done by calculating the Euclidean distance (via the base diff() function) between the normalised and weighted samples.

