Step1- Download Data in this folder from repository
https://data.mendeley.com/datasets/pxpd8vc3fh/draft?a=bcbf2750-774f-4a56-8afd-50f852bed4a8

Description of Data

This data contains the result of experiments on ball bearings of the series 6204_10.
For details on the procedure of how to prepare the bearings please refer to the section "Experimental Setup" of the article.

The conditions of the experiments are:

Experiment 1-6
	- Dry
	- Dry + interference
	- Lubrication 5%
	- Lubrication 100%

Experiment 7-12
	- Dry
	- Lubrication 100%


The data is tagged using the bearing ID. Each bearing ID was associated with an experiment number:

Experiment # = Bearing_ID
E1 = B71
E2 = B7
E3 = B8
E4 = B9
E5 = B10
E6 = B11
E7 = B26
E8 = B27
E9 = B28
E10 = B31
E11 = B32
E12 = B34

-For the fault indicators comparison data from experiments 1-6 was used.
-For the fast-kurtogram and envelope analysis the data from experiments 1-12 was used.
-For the run-to-failure test the Bearing_ID=B4 was used.

The filename contains tags to indicate modifications to the raw data. The meaning is as follows:
HF - The data was pre-filtered using a Highpass 1Khz, Bandstop 15Khz and Lowpass 20Khz.
P1 - Using shim of thickness of 5mils to increase interference between outer ring and bearing housing.
min - Lubrication 5%
full -  Lubrication 100%
Dry - Lubrication 0%

For details on the selection of filters please refer to the section "Experimental Setup" of the article.

The filename for run-to-failure (B4) contains tags to identify the conditions for interference. 
The shim thicknes used to increase the interference is tagged as follows:
I5 - 5 mils 
I10 - 10 mils
I15 - 15 mils







