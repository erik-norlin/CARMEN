# CARMEN: A Cellular Automaton Rain Flow Model with Evaporation and Infiltration

### Aim of this project ### 
This has been a smaller scientific group project where we were given the task to implement a mathematical model to simulate a complex system. The aim of this project was to implement a cellular automatan to simulate water flow dynamics over a real world region from gathered satellite data (hydrological modelling), to ultimately allow investigation of potential flooding for hypotheical scenarios given a data set of a region. This was done by implementing the mathematical model from the paper [*A two scale cellular automaton for flow dynamics modeling (2CAFDYM)*](https://www.sciencedirect.com/science/article/pii/S0307904X16305492) (H. Kassogué et al., 2017) in python using object oriented programming.

### What the project does ### 
The program allows one to observe how the water flow dynamics over a particular real world region in Morocco is being directly affected by the landuse and soil type. With this program, one can artificially switch the landuse and soil type of the region to all infrastructure, all forest and real world data to see how the water dynamics change with respect to the landuse and soil type.
  
### How to run the project ### 
There's no installation required to run this project. Just connect to the repository and run the main program. To choose between different soil types over the studied region, simply uncomment the lines of code in both the main program and the CellObject class where it specifies landuse and soil type. Then, run the program and observe the magic (plots are found in the Plots folder).

![Water running down landscape](https://github.com/erik-norlin/CARMEN/blob/master/Plots/Qps/Qps_forest/forest_t%3D360000.jpeg?raw=true)

### Creators of this project ### 
* Thibault Desjonquères
* Paulina Essunger
* Hannes Nilsson
* Erik Norlin

For more details, please read the project report.
