[![CI](https://github.com/kalmarek/PermutationGroups.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/kalmarek/PermutationGroups.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/kalmarek/PermutationGroups.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/kalmarek/PermutationGroups.jl)

# PermutationGroups.jl

This project is mostly for self-educational purposes, as my preparation for teaching a course on computational group theory.
Some basic performance considerations were taken into the account while designing the package, but no extensive performance tuning was performed.

The package is mostly aimed at working with permutation groups of small to medium degree (up to degree of a million lets say). The basic functionalities include computing the stabilizer chains via Schreier-Sims algorithm which gives access to compuitng the order of a group and membership tests.

The code should be well documented. See also `test` folder for example usage.
