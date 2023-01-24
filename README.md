# Compression using Binary Contours

The code in this repository is a prototype for the BinCon compression method
described in the article. It is not intended for operational or production
use, but rather as an example of how the compression method can be implemented.

The code is written in the [Julia](https://julialang.org/) programming language.
It shouldn't be too difficult to re-implement the BinCon method in other
languages. Note: When I wrote the code, I was learning Julia. So bear this in
mind if you are using it.

The most important file is [BinaryContour.jl](./BinaryContour.jl).

Two example scripts are provided to demonstrate how to use BinCon.
[GetGFS.jl](./GetGFS.jl) downloads recent NCEP GFS data and compresses it using
BinCon. [plot.jl](./plot.jl) plots compressed data files.

## Documentation

- Documentation on downloading and compressing GFS data: [GetGFS.md](./GetGFS.md)
- Documentation on plotting compressed data: [plot.md](./plot.md)
- Documentation on the predefined map regions: [regions.md](./regions.md)

## Installation

The dependencies of the software include Julia and [The Generic Mapping
Tools](https://www.generic-mapping-tools.org/), including the
[Julia wrappers](https://github.com/GenericMappingTools/GMT.jl/releases).

More installation instructions to be written.

## Contact

I'm happy to hear how people are using this code, or if people have questions.
My Email is tim@nomuka.com


Tim Hume
January 2023
