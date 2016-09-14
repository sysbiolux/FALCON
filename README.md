# FALCON
FAst Logical COntextualization  of Networks


Contact info: thomas.sauter@uni.lu  
Please cite this software as follows: De Landtsheer et al., XXX  
Copyright: This software is freely available for non-commercial users under the GPLv3 license
  

FALCON is a toolbox for the efficient contextualization of logical network models. 
Contextualized network models provide important qualitative and quantitative information about the system being modelled.
Specifically, FALCON is well suited to assess the relative contributions of different signal transduction mechanisms to the behavior of the system at steady-state. 
The logical formulation of the interaction between different molecules is intuitive, and often is the only available information to systems biologists. 
A Bayesian interpretation of the logical ‘gates’ allows for an algebraic formulation of the system and an efficient calculation of the long-term steady-state of the system given the specified inputs. 
A gradient-descent optimizer (fmincon) is used to minimize the error function calculated as the sum of squared residuals between the simulated and experimentally measured nodes of interest. 
FALCON offers multiple types of analysis including identifiability analysis, systematic knock-out and differential regulation which can be applied to various biological applications.
