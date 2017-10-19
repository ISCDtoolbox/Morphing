# Shape Morphing [![Build Status](https://travis-ci.org/ISCDtoolbox/Morphing.svg?branch=master)](https://travis-ci.org/ISCDtoolbox/Morphing)

This repository contains the code morphing for deforming a template mesh onto a target shape associated to the journal article
http://www.sciencedirect.com/science/article/pii/S1631073X16300802.

### Installation

You can grab the sources by cloning this repository or downloading a .zip archive of the sources. In order to build the project, navigate to the created directory and in a command prompt type:

```
mkdir build
cd build
cmake ..
make
(sudo) make install
```


### Usage

In a terminal, run:

```
morphing [-h] [-nit n] [-dref nref [refs]] [-bref nref [refs]] [-elref nref [refs]] target_file[.mesh] template_file[.mesh]  
```

The different parameters correspond to:
* **dref**  : fixed surfaces , *inside* the mesh to be deformed.
* **bref**  : "follower" elements on the exterior surface, which must not be morphed (ears for instance)
* **elref** : tetrahedral elements inside the fixed surfaces

Together with the target_file[.mesh] a scalar field target_file.sol should be provided, corresponding to the signed distance function.

Please run:
```
morphing -h
```
to see the default parameters.

For the examples provided in the demo/ folder, the correct command to use would be:
```
morphing -dref 1 2 -elref 1 1 -bref 1 1 target.mesh template.mesh
```


### Authors and contributors

morphing has been initiated by Maya de Buhan (Université Paris Descartes) and Chiara Nardoni (Université Pierre et Marie Curie).
Contributors to this project are warmly welcomed.

### License

morphing is given under the [terms of the GNU Lesser General Public License] (LICENSE.md).

If you use morphing in your work, please refer to the journal article:

An optimization method for elastic shape matching, M. De Buhan, C. Dapogny, P. Frey, C. Nardoni, C.R. Acd. Sci., Paris, Sèrie I, 2016.
