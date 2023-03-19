# A quick installation

`$ git clone --depth=1 git@github.com:rashti-alireza/Elliptica.git`

`$ make MyConfig`

Activate the projects of interest and link UMFPACK library in the `MyConfig` file.
You should also activate the project's repositories if you need to clone them.

Clone the added projects:

`$ make git_clone`

Install:

`$ make -j4`

