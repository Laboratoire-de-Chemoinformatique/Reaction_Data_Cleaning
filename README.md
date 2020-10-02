# Reaction_Data_Cleaning
This repository contains chemical reactions data curation best practices.
#### In 'scripts' package you will find:
* standardization protocol;
* a script to run the standardization protocol in parallel mode using chunks.
#### In 'data' directory you will find:
* Our golden dataset zip archive curated and mapped manually;
* USPTO dataset curated by the standardization protocol and mapped by RXNMapper.

### Corresponding Authors:
    Alexandre Varnek (varnek@unistra.fr)
    Timur Madzhidov (tmadzhidov@gmail.com)
    
### Contributors:
    Arkadii Lin (arkadiyl18@gmail.com)
    Ramil Nugmanov (nougmanoff@hotmail.com)
    Natalia Duybankova (NDyubank@its.jnj.com)
    Jonas Verhoeven (jverhoe9@its.jnj.com)
    Timur Madzhidov (tmadzhidov@gmail.com)
    Alexandre Varnek (varnek@unistra.fr)
    Joerg Wegner (jwegner@its.jnj.com)
   
### Copyright:
    Copyright 2020, MaDeSmart, Machine Design of Small Molecules by AI VLAIO project HBC.2018.2287
    
### Credits:
    Kazan Federal University, Russia
    University of Strasbourg, France
    University of Linz, Austria
    University of Leuven, Belgium
    Janssen Pharmaceutica N.V., Beerse, Belgium
    Rail Suleymanov, Arcadia, St. Petersburg, Russia

# Reference
Please, cite the paper when you use the data or the scripts:
  Lin, Arkadii; Dyubankova, Natalia; Madzhidov, Timur; Nugmanov, Ramil; Rakhimbekova, Assima; Ibragimova, Zarina; Akhmetshin, Tagir; Gimadiev, Timur; Suleymanov, Rail; Verhoeven, Jonas; Wegner, Jörg Kurt; Ceulemans, Hugo; Varnek, Alexandre. (2020): Atom-to-Atom Mapping: A Benchmarking Study of Popular Mapping Algorithms and Consensus Strategies. ChemRxiv. Preprint. https://doi.org/10.26434/chemrxiv.13012679.v1 

### Dependencies
  - python: 3.7
  - CGRtools: 4.0.36
  - ordered-set: 4.0.2
  - pyjnius: 1.3.0
  - JChemSuite package from ChemAxon: 19.9.0
  
