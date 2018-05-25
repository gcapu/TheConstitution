# The Constitution
Basic library describing some constitutive relationships of materials in solid mechanics.
This is a work in progress and all function interfaces could--and very likely will--change. 

## Undecided aspects
Since each material has its own requirements, making a uniform interface for all of them is challenging. The easiest option would be to use a single function signature that possess all possible inputs, but the user would need to calculate unnecessary things to feed to that function. 

An alternative is to give each material its own signature (maybe chosen from a predefined set of variants) and have enable_if-like sections in the user code calculate the required inputs and feed them to the material functions. This approach makes the library difficult to use and the user code quite ugly.

We could also have both options: A set of material functions with different function signatures, and, for each one of them, a wrapper with the common signature. Users can then use either option. 

I haven't decided yet how to approach the problem and I'm just going to generate some particular materials with the variables they need. For the moment I will give them the deformation tensor F and a few other parameters. Once I reach a reasonable number of materials in the library I will decide how to proceed.