# Boolean.jl

This code is used to simulate networks using an Ising Boolean formalism, previously used by Font-Clos et al. Here is how to use it:

**1.** Clone the GitHub Repository
```
git clone https://github.com/askhari139/Boolean.jl
```
**2.** Install Julia (tested on 1.4 and 1.6 versions of Julia). 

**3.** Run the script ``dependencyInstaller.jl``. This will install all required packages.
```
julia dependencyInstaller.jl
```
**4.** Copy the ``script.jl`` to the folder that has the topo files. The topo files are space-delimited files that list out all the edges in a network. Each topo file has 3 columns: Source, Target and Type (1 for activation, 2 for inhibition).

**5.** Edit the first line in ``script.jl`` to include the path to the cloned repository.

**6.** Run ``script.jl``
```
julia script.jl
```

**7.** Voila! All topo files in your folder are set to be simulated!
