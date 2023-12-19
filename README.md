This is an example code of a flutter analysis for a 2D Section, with 3 degrees of freedom (Heave, Pitch, Control Surface). 

The goal of ths project is to calculate the open loop flutter speed for a 2D section from the PhD thesis of Prof. Moti Karpel and verify the speeds and flutter modes match. 

In order to incorporate the unsteady aerodynamics into the model, we utilize Theodorsen's Reduced frequency in order to model the continuous changing motion, phase shift and the additional lag states of the degrees of freedom.  

Once these equations were built, they are setup into a State Space Model, which reduces the order of the system to solve. 

Further details and equations are present in the theory file.
