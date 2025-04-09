Welcome to

```ascii
EEEEEEEEEEEE  ll  ll
EE            ll  ll
EE            ll  ll  O                       O
EE            ll  ll                   tt
EEEEEEEEEEEE  ll  ll  ii  pppppppp. tttttttt  ii    .ccccc   .aaaaa.
EE            ll  ll  ii  pp      p.   tt     ii   .c       .a     a.
EE            ll  ll  ii  pp       p.  tt     ii  .c       .a       a.
EE            ll  ll  ii  pp      p.   tt     ii   .c       .a     .a
EEEEEEEEEEEE  ll  ll  ii  pppppppp.    tt     ii    .ccccc   .aaaaa.aa
                          pp                                          aa
                          pp
                          pp
                          pp

```

# Installation

Requirements:
- GNU make.
- Git version control system.
- A C programming compiler, for instance, `gcc`.
- A compatible `openmp` library with the compiler.
- The `UMFPACK` library (`suitesparse`).

`$ git clone --depth=1 git@github.com:rashti-alireza/Elliptica.git`

`$ make MyConfig`

`$ make git_clone`

`$ make -j4`

You can see the list of available projects after issuing `make git_clone` in `Src/Projects`.

# Run

In each project there is a directory named `XXXX_ParFiles`, where `XXXX` denotes the project name, 
e.g., `NSNS_ParFiles`.
One can find various parameter files for different gravitational systems. 
Having set the parameters of interest, one invoke the code:

`$ ./Exe/elliptica /path/to/the/parameter/file`


# Citation

inspirehep.net citation key = `{Rashti:2021ihv,Rashti:2024drr}`


