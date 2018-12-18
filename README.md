# Welcome to the clear-sky irradiance model library
This Github repository contains many clear-sky irradiance models as coded in R, though occasionally Matlab.
This repository was created to conicide with our research publication titled "Rigorous worldwide performance assessment of 75 global clear-sky irradiance models using Principal Component Analysis" published in the Journal of Renewable and Sustainable Energy Reviews and authored by Xixi Sun, Jamie M. Bright, Christian A. Gueymard, Brendan Acord, Peng Wang and Nicholas A. Engerer.

The R code available in this repository was written by Xixi Sun.

## About the models
The reader is referred to our publication and its accompanying supplementary material in order to fully understand the workings of these models. There you will find our interpretation and links to the original work.
PLACEHOLDER URL TO JOURNAL PAPER

The models all require different input data


|No. |Clear-sky Model|Eo|zen|h|alb|p|T|TL|aod|alp|beta|O3|NO2 |H2O|tau|Tot 
|---|------|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:
|1 | TJ || •| | | | | | | | | | | | | 1 
|2 | Schulze || •| | | | | | | | | | | | | 1 
|3 | DPP || •| | | | | | | | | | | | | 1 
|4 | Adnot || •| | | | | | | | | | | | | 1 
|5 | Biga ||•| | | | | | | | | | | | | 1 
|6 | ASHRAE ||•| | | | | | | | | | | | | 1 
|7 | Sharma | •|•| | | | | | | | | | | | | 2 
|8 | El Mghouchi | •|•| | | | | | | | | | | | | 2 
|9 | Yang & Walsh | •|•| | | | | | | | | | | | | 2 
|10 | HLJ |•|•|• ||||||| | | ||| 3 
|11 | Kumar |•|• |||• ||||||| | || 3 
|12 | Campbell |•|• |||• ||||||| | || 3 
|13 | Fu \| Rich | • |• |•|| | |||| | | ||| 3 
|14 | Atwater & Ball-1 |•|•|| |•| | | | | | | |• | |4 
|15 | KASM |•|•| | |•| | | | | | | |• | |4 
|16 | Capderou |•|•| •| |•| | | | | | | | | |4 
|17-22 | Kasten | •|• |• ||||• || | | | | | | 4 
|23-28 | Heliosat-1 | •|•|||• ||• ||| || | | | 4 
|29-34 | ESRA | •|• |•| | | |• | | | | | | | | 4 
|35-40 | Heliosat-2 |•|• |• ||||• ||| | | || | 4 
|41-46 | Ineichen & Perez | •|•|• ||||• || | | ||| | 4 
|47 | CLS |•|•||•|• ||| |||||• | |5 
|48 | King |•|•||•|• ||| ||•||| | |5 
|49 | Josefsson | •|•||•|• ||| |||||• | |5 
|50 | Badescu | • |• |||• ||| |||•||• | |5 
|51 | Simplified Solis | •|•|||• |||• |||||• | |5 
|52 | Advanced Solis | •|•|||• |||• |||||• | |5 
|53 | Perrin | •|•|||• ||| ||•|•||• | |6 
|54 | CEM | •|•||•|• ||| • | ||||• | |6 
|55 | Atwater & Ball-2 | •|• ||• |•|| |•|||||•| |6 
|56 | RSC | •|•||•|•|| | ||•|||•| |6 
|57 | PSIM | •|•||•|•|| | ||•|||•| |6 
|58 | Bashahu |•|•| | |• | | | |•|•| ||•| |6 
|59 | MMAC | •|•||•|•|||•|||||•| |6 
|60 | Yang | •|•|||•|||||•|•||•| |6 
|61 | Calinoiu |•|•||||||||•|•|•|• | |6 
|62 | Hoyt | •|•||•|•|||||•|•||•| | 7 
|63 | MAC | •|•||•|•||||•|•|||•||7 
|64 | METSTAT | •|•| |•|•|||•| | |•| |•| | 7 
|65 | PSI-REST | •|•| |•|•| | | | |•|•| |•| | 7 
|66 | Paulescu & Schlett | •|•|•| |•| | | | |•|•| |•| | 7 
|67 | MRM v5 | •|•||•|•|||||•|•||• | |7 
|68 | MRM v7 | •|•||•|•|||•| | |•||•|| 7 
|69 | Janjai |•|•|•| |||||•|• |• ||•| | 7 
|70 | Bird | •|•||• |•||| |•|• |• ||• | |8 
|71 | Iqbal-c | •|•||•|•||| |•|• |• ||• | |8 
|72 | Modified Iqbal-c |•|•|•|•| ||||•|•|•||•| |8 
|73 | REST2v5 | •|•||•|•||||• |• |• |• |• || 9 
|74 | REST2v9 |•|•||•|•||||•|•|•| |•|•| 9 
|75 | McClear | •|•|•|•|•|•|| •|•| | •| |•| | 10 


We obtained these input data from the MERRA2 reanalysis database. 
The reader is again referred to the data section of our paper in order to learn how to access these inputs.

