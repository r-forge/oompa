Notes for building the newest packages. Perform steps 1-3 for each
package separately. After you are finished, perform steps 4 and 5 once
to take care of all packages.

-------------------------------------------------------------------
[1] Update the 'DESCRIPTION' file to have the correct version number,
date, and 'depends' line.

[2] From an operating system command line, run:
    Rcmd check [PACKAGE]
on Windows (or "R CMD check [PACKAGE]"on UNIX). After everything
passes without warnings, move the resulting [PACKAGE].Rcheck directory
into the 'RCheck" folder. 

[3] From an operating system command line, run both 
    Rcmd build [PACKAGE]
and
    Rcmd build --binary [PACKAGE]
which will produce a .tar.gz file containing source and a .zip file
containing a Windows binary. Move these files into the correct
subdirectory of the 'Releases" directory.

NOTE: Some packages include binary R data objects that may need to be
rebuilt for newer versions of R. It is usually easiest to build and
install the package, run the code from the current version of R that
regenerates the data objects, put them in the right place, and then go
back and repeat steps 1-3 for the data-bearing package.

NOTE: Some packages include sample code that has dependencies other
than those included in the 'DESCRIPTION' file. These are primarily
related to examples showing that the methoids work with ExpressionSet
objects, and they use the data fromthe TailRank package for this
purpose. These should also be retested after installation. If changes
need to be made, then steps 1-3 must be repeated.

-------------------------------------------------------------------
[4] Copy all of the .tar.gz files into the repository at (*)
    v:\bioinformatics\OOMPA\src\contrib
Inside R, run the function
    library(tools)
    write_PACKAGES("V:/bioinformatics/OOMPA/src/contrib")

[5] Copy all of the .zip files into the repository at (*)
    V:\bioinformatics\OOMPA\bin\windows\contrib\[x].[y]
where [x].[y] represents the major version of R for which the binaries
were compiled. Then, inside R, run the function
    library(tools)
    write_PACKAGES("V:/bioinformatics/OOMPA/bin/windows/contrib/[x].[y]")
substituting appropriately for [x] and [y].

TODO: Find someone who knows how to produce mac binaries, and get them
to help make them and put them into the repository....

-------------------------------------------------------------------
(*) Note that drive V: on my Windows machine expands to the network
path: 
    \\mdadqsfs01\bio_web

-------------------------------------------------------------------
People can now download the latest versions from the repository using
the R command

install.packages(repos="http://bioinformatics.mdanderson.org/OOMPA")
