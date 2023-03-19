# A quick installation

Requirements:
- GNU make.
- Git version control system.
- A C programming compiler, for instance, `gcc`.
- A compatible `openmp` library with the compiler.
- The `UMFPACK` library (`suitesparse`).

`$ git clone --depth=1 git@github.com:rashti-alireza/Elliptica.git`

`$ make MyConfig`

Activate the projects of interest and link UMFPACK library in the `MyConfig` file.
You should also activate the project's repositories if you need to clone them.

Clone the added projects:

`$ make git_clone`

Install:

`$ make -j4`

Run:

`$ ./Exe/elliptica /path/to/the/parameter/file`




