# Installation instructions for Windows

The following instructions explain how to install the BinCon software on a Windows
computer.

## Install the Generic Mapping Tools

Install the latest version of the Generic Mapping Tools. Currently this is:
[Latest GMT](https://github.com/GenericMappingTools/gmt/releases/tag/6.4.0)

Choose the appropriate installer, and run it. The default installation settings
are fine.

## Install Julia

Julia can be downloaded and installed from [here](https://julialang.org/downloads).
When installing Julia, make sure to click the "Add Julia to Path" box.

Now open a Windows Powershell terminal (you may need to install this if you don't
already have it).

Run Julia from the command line:

![Julia REPL](images/julia.png "Julia REPL")

Now install the Julia GMT wrappers (from within the Julia REPL):

```
]
add GMT
```

When this has finished you can exit from the REPL.

## Download the files you need from GitHub

Copy the files from [this GitHub repository](https://github.com/timo007/bincon)
into a folder on your computer. You can either download the files using a
browser, or if you have git installed on your computer, you could clone the
repository to somewhere local. You don't need the files in the test_data,
article or images subdirectories on GitHub.

Navigate to the folder you've just created in the Windows Powershell terminal.
Now you can run commands like this (you'll need to change the date and time
of the data as appropriate):

```
julia.exe .\GetGFS.jl -t 2023012412 -v prmslmsl --reg NZ --tol 0.25 --cnt 200
julia.exe .\plot.jl -i GFS_NZ_000-003-001_MSL_2023012412_048.bin 
```

To test downloading and plotting a compressed data file from a remote server, try this:

```
julia.exe .\plot.jl -i https://github.com/timo007/bincon/blob/main/test_data/GFS_NZ_000-003-001_MSL_2023010412_024.bin?raw=true -o example.pdf 
```

The resultant plot should look like this
![Example plot](images/example.png "Example plot")
