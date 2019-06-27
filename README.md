# Nebula
*If you don't know what this is, just ignore it®™*

## Getting Nebula
Via the command-line with git installed:
```sh
git clone --recurse-submodules https://github.com/lvkessel/nebula.git
```
If you forgot the `--recurse-submodules` parameter, you can fix it by doing
```sh
git submodule init
git submodule update
```

Via a web browser, without git installed:
* Download the [zip file](https://github.com/lvkessel/nebula/archive/master.zip) from github.
* Download nvidia's [cub](https://nvlabs.github.io/cub/download_cub.html) library, and place it in the `3rdparty/cub/` folder. (So there is a `3rdparty/cub/cub/cub.cuh` file).

## Compiling
Nebula has been tested with the following compilers and versions:
* CUDA 8.0 / gcc 4.8
* CUDA 9.1 / Visual Studio 2017 15.4

Later compilers are pretty much guaranteed to work. Earlier versions of these compilers might also work, but they haven't been tested. In addition to these compilers, we require the following:
* cmake 3.8 or greater
* hdf5 1.8.13 or greater, with the high-level library.

There are three main programs at the moment:
* `nebula_gpu`, which uses CUDA for fast simulations
* `nebula_cpu_mt`, a multithreaded CPU version.
* `nebula_cpu_edep`, a CPU version that outputs information on energy deposits.

They are all compiled with cmake:
```sh
mkdir build
cd build
cmake ..
make
```
The executables can then be found in the `build/bin` folder.

## Running
Similar to `e-scatter`, though with fewer options:
```sh
./nebula tri.tri pri.pri mat1.mat [..] matN.mat
```
The output, which is fully compatible with e-scatter, is written to the standard output for the GPU and the multi-threaded CPU version. The "inspect" version writes the same output to a file `tmp.bin` (sorry), while writing interesting cascade data to stdout. New-style `.hdf5` material files are recommended. Old-style `.mat` files are also supported for now, but beware that this will be removed soon, and it does not work with the Penn inelastic model.

## Changing physical models
All programs include the `physics_config.h` file, and expect it to define two data types.
* `scatter_physics`, which is a `std::tuple`-like type enumerating the scattering processes present in the simulation.
* `intersect_t`, which defines the material-material boundary intersection physics.

For most purposes, the file should speak for itself. If not, contact me. You can flip some switches in the template parameter lists for the inelastic, elastic and intersection physics. You can use the Kieft inelastic model instead of the Penn inelastic model by commenting the relevant lines.

## Inspecting the cascade
There is one default way of inspecting the cascade. The `nebula_cpu_edep` program runs on the CPU, and outputs each event where energy was lost. By default, it outputs the usual detected electron data to a file `detected.bin`, while energy-deposition data is sent to stdout.

The file format is 5 floats and 2 ints. The floats are scattering position (x, y, z), the electron's energy before the event, and energy lost. The two ints are the pixel coordinates belonging to the primary electron. This data is only output for scattering events where energy was lost: pure elastic scattering events are not output. Be aware that HUGE amounts of data are generated.
