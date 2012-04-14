# Introduction

Original source code retrieved from: http://creckmodeling.chem.polimi.it/openfoam.html

A bit more information about these modifications: http://www.cfd-online.com/Forums/blogs/wyldckat/1280-libopensmoke-making-work-once-again-openfoam.html

# Installation
## Prerequisites
Firstly, you need the GNU Scientific Library:

*   On Ubuntu, install GSL by simply running:

        sudo apt-get install gsl-bin libgsl0-dev

*   To build your own GSL, read the section "Compilation and Installation" on the [doc/UserGuide.pdf](libOpenSMOKE/blob/master/doc/UserGuide.pdf?raw=true "User Guide") document.

    Keep in mind that after you build your own GSL, doing `git clone` and before running `Allwmake`, edit that file `Allwmake` and change the path `$HOME/gsl` to your own installation path of GSL:

        export GSL_INST_DIR="$HOME/gsl"

## Installing `libOpenSMOKE`
To use the build with OpenFOAM 1.7.x (maybe works with 1.7.0 and 1.7.1), simply run:

    git clone https://github.com/wyldckat/libOpenSMOKE.git
    cd libOpenSMOKE
    ./Allwmake

To build with OpenFOAM 2.0.x or 2.1.x or above, before running "Allwmake", check which versions exist by running:

    git branch -a

For example, if it shows `remotes/origin/v2.0.x`, then run these steps:

    git checkout v2.0.x
    ./Allwmake

The binaries will be installed in your own user folders, so please DO NOT INSTALL AS ROOT! To find out which are your own user folders for binaries, run the following commands:

    echo $FOAM_USER_LIBBIN
    echo $FOAM_USER_APPBIN

# Notes
* The files on the folder "utilities/exe" are static binaries. I don't know where the original source code is, so you better ask the authors at http://creckmodeling.chem.polimi.it/openfoam.html
