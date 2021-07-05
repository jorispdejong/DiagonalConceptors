# Controlling Recurrent Neural Networks by Diagonal Conceptors
This repository is complementary to a master research project conducted at the Mathematics (statistics and big data) department of the RijksUniversiteit Groningen. 


## Introduction
The subject of this research was diagonal conceptors, a practical alternative to conceptors. The conceptors architecture is a neuro-computational mechanism, which was first introduced in the technical report "Controlling Recurrent Neural Networks by Conceptors", written by H. Jaeger. Among other things, conceptors provide a method for learning, self-generating, and combining a collection of temporal patterns in a Recurrent Neural Network (RNN). Intuitively, the idea of conceptors is that driving an RNN with patterns 1,2,...,p excites different areas A1,A2,...,Ap in the neural state space. After pattern j is learned by the RNN, it can be retrieved by projecting the neural state onto Aj via a projection matrix Cj, which is the conceptor matrix associated with pattern j. The conceptor matrices are square matrices, the size of the number of neurons in the RNN. Therefore, conceptors quickly become impractical. 

In this research, a alternative approach is proposed, called diagonal conceptors, which offers a practical alternative for conceptors. Diagonal conceptors are diagonal matrices, hence can be written as vectors. They are shown to yield equally good results as conceptors in most cases.

For more details, see the report.

## Repository Structure
