# Source code for EPITA Ani3D labs

## Download

_Beware: Do not download the code via the github "Download ZIP" GUI (cgp submodule are not included)_

### Via Clone

```
git clone --recurse-submodules https://github.com/drohmer/epita_ani3d_lab_code.git
```

or 

```
git clone https://github.com/drohmer/epita_ani3d_lab_code.git
cd epita_ani3d_lab_code/
git submodule init
git submodule update
```

### Archive download

* [epita_ani3d_lab_code.zip](https://imagecomputing.net/course/2023_2024/epita_ani3d/lab_code/epita_ani3d_lab_code.zip)
* [epita_ani3d_lab_code.tar.gz](https://imagecomputing.net/course/2023_2024/epita_ani3d/lab_code/epita_ani3d_lab_code.tar.gz)



## Compilation

Scenes are in a directory _scenes_epita_ani3d/xx_Name/_


Each scene is provided with its own CMakeLists.txt and Makefile

Compilation via CMake
```
cd scenes_epita_ani3d/01_introduction/
cmake -B build/
cd build/
make -j$(nproc)
./01_introduction
```

Compilation via Makefile
```
cd scenes_epita_ani3d/01_introduction/
make -j$(nproc)
./01_introduction
```

More details are provided in the dedicated [Installation setup and Compilation page](https://imagecomputing.net/cgp/compilation) if needed (Windows, Linux, Mac).

## Requirement

The code requires the GLFW library in Linux.

ex. on Ubuntu

```
sudo apt-get install libglfw3-dev
```



## Links

* [Lab class exercises](https://imagecomputing.net/course/2023_2024/epita_ani3d/lab/index.html)
* [Documentation on CGP library](https://imagecomputing.net/cgp/index.html)

