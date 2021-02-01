# Covid-19 SEIRD implementation

In this repository, an R-implementation of an SEIRD compartment model for modeling the Covid-19 epidemic is stored. The main file is 'SEIRD_Analysis.R' which uses the script 'Model_core.R', which contains the core of the iterative implementation, the auxiliary scripts 'general_functions.R' and 'general_settings.R' and contact rate data 'ContactRates.txt'.

To run the main script, copy all files into a single folder, and set the working directory to the corresponding folder path.




## Scope

The implementation is intended to enable forecasting of hypothetical scenarios under different parameter specifications. Since an extended SEIRD model is used, the same limitations, which generally hold for compartment models, apply here. Although the considered modeling can be effectively used in investigating possible outcomes of the epidemic, we emphasize that it is not our primary goal of establishing a forecasting model with factual validity.



## Scenarios

In the main script, one can choose from 4 different hypothetical scenarios, which are studied in a related research article by Achim Doerre and Gabriele Doblhammer. These relate to different variations of partial lockdown measures. Since they are hypothetical in nature, the scenarios do not necessarily represent the actual course of the epidemic.




## Changing the parameters

All code is edited in a way that enables modification of specific parameters. In particular, if other specifications of the epidemiological parameters are desired, they can be changed in the main file. Similarly, if for the uncertain parameters, such as the manifestation index, other ranges seem plausible, they can be modified in the script. Finally, if the implementation is intended for a different setup of groups (e.g., different partition of age groups and gender), the matrix of contact rates has to be adjusted as well as related steps in the main code.
