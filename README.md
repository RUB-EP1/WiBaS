# WiBaS

WiBaS ("Williams' Background Suppression") is a lightweight and 
easy-to-use C++ implementation of the Q-value method published by Mike Williams:

http://arxiv.org/abs/0809.2548  
http://iopscience.iop.org/article/10.1088/1748-0221/4/10/P10003

WiBaS is written by Julian Pychy as a member of the PANDA collaboration (https://panda.gsi.de) 
and is distributed under the GPLv3.
The package allows a simple definition of the phasespace metric and configuration of the fitting processes using C++ objects. 
Several analyses making use of it are in progress and first results have been published:

http://epjc.epj.org/articles/epjc/abs/2015/03/10052_2015_Article_3341/10052_2015_Article_3341.html

## Usage
The files  in the *src* folder can easily be included in any C++ project. 
The second folder contains two commented example analyses
in which the background below an omega meson signal is removed.

### ACLiC

In this case the package is used within an ACLiC compiled ROOT script (https://root.cern.ch/).
In an Unix environment, start by executing ``root -l root_startscript.cc`` in the corresponding folder.

### Standalone

This folder contains an example how to use the package in a gcc compiled application, which requires linking against 
installed ROOT and RooFit libraries. Just type ``make`` to build. 
