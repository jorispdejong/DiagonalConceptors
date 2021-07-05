# Controlling Recurrent Neural Networks by Diagonal Conceptors
This repository is complementary to a master research project conducted at the Mathematics (statistics and big data) department of the RijksUniversiteit Groningen. 


## Introduction
The subject of this research was diagonal conceptors, a practical alternative to conceptors. The conceptors architecture is a neuro-computational mechanism, which was first introduced in the technical report "Controlling Recurrent Neural Networks by Conceptors", written by H. Jaeger. Among other things, conceptors provide a method for learning, self-generating, and combining a collection of temporal patterns in a Recurrent Neural Network (RNN). Intuitively, the idea of conceptors is that driving an RNN with patterns 1,2,...,p excites different areas A1,A2,...,Ap in the neural state space. After pattern j is learned by the RNN, it can be retrieved by projecting the neural state onto Aj via a projection matrix Cj, which is the conceptor matrix associated with pattern j. The conceptor matrices are square matrices, the size of the number of neurons in the RNN. Therefore, conceptors quickly become impractical. 

In this research, a alternative approach is proposed, called _diagonal conceptors_, which offers a practical alternative for conceptors. Diagonal conceptors are diagonal matrices, hence can be written as vectors. They are shown to yield equally good results as conceptors in most cases.

For more details, see the report.

## Repository Structure
The _scripts_ folder contains the main scripts for diagonal conceptors and conceptors for three examples (four periodic patterns, chaotic attractors, human motions). The scripts in this folder usually source a number of functions, which are defined in the scripts in the _source_ folder, and they import data from the _data_ folder. The outputs are saved in the _output_ folder. The output folder is divided into plots and the reservoirs, where a reservoir is a list containing everything that is required to set up the RNN.

The scripts are commented, so they should provide you with everything you need when running a simulation. The required packages should be in the _source/Libraries.R_ file, so before running any script, it may be required to install those packages. Besides that, the scripts can be run without any setup. After reading the theory of the report, the scripts should be intuitive and self-explanatory. The scripts do not take long to execute (mostly under 1 min) and the progress can be shown by setting the parameter _verbose=T_. Note that these scripts are not finished and can surely by optimized. However, they are a great introduction to diagonal conceptors and conceptors.


