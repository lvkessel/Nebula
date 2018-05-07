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
* hdf5 1.8.12 or greater, with the C++ and high-level interfaces.

There are three main programs at the moment:
* `gpu_main.cu`, which, incidentally, also supports a simple CPU version with a preprocessor flag.
* `cpu_mt_main.cpp`, a multithreaded CPU version.
* `cpu_inspect_main.cpp`, a single-threaded CPU version that uses the `tracking_cpu_particle_manager` to store and print some cascade information.

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
If you want to extract information about the cascade, you need to tell the simulator what data to extract. Because of the large number of potentially interesting things, the design is deliberately "hacky": you have to dive into the simulator code, where you have access to all variables. Knowledge about the internals of the simulator is also required. If you get stuck, ask me. Documentation is on its way.

Use the `cpu_inspect_main.cpp` for compilation. It uses the "tracking particle manager" by default, in which you can do your hacking. It can be found in `drivers/cpu/tracking_cpu_particle_manager.h`. It overrides the `create_secondary` and `set_scatter_event` functions from the normal particle manager. These overrided functions store the positions and kinetic energies for all scattering events.

When the `flush_detected` function is run, all parent particles of each electron are guaranteed to still be in memory. The `data` member (defined in the base class, `drivers/cpu/cpu_particle_manager.h`) is an array of all particles in memory. Let's say you have an element from this `data` array, you can use:
* `data[i].primary_tag`: for identifying if two particles belong to the same cascade. This is an integer and it is not unique.
* `data[i].unique_tag`: a globally unique tag, can be used to identify a primary or secondary electron. Be warned that a primary electron's `primary_tag` is not the same as its `unique_tag`!
* `data[i].parent_unique_tag`: this particle's parent's unique tag. Can be used to find the parent. Equal to zero if this is a primary particle.
* `data[i].parent_create_event`: the index to the parent's `events` array corresponding to the event in which this particle was created.
* `data[i].events`: an array of `event_info`, containing the following:
  * `data[i].events[j].type`: by default, one for inelastic and two for elastic scattering
  * `data[i].events[j].position.x/y/z`: position at which the scatter event took place
  * `data[i].events[j].energy`: kinetic energy BEFORE the event

By default, the `tracking_cpu_particle_manager` finds the deepest position (i.e. the lowest z coordinate) at which the electron, including all its parents, has been. It then prints this depth, along with the kinetic energy at that depth and the kinetic energy at detection, to stdout. This should serve as an example you can use to understand the system.

Because the data is printed to stdout, it can be redirected to a file by running the simulation as `./nebula_inspect tri.tri pri.pri mat.mat > text_file.txt`.
