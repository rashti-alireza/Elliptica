### interpreting analytic derivative functions 
### produced by Sympy (Derivative) or Mathematica (D):

### Derivative:

### d^3/d?^3

# (x,y,z)
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?x ?\, ?y ?\, ?z ?)\)/\(bbn_ks_ddd\1_D0D1D2 KS_func_pass_args_macro \)/g;

# (x,(y,2)) | ((y,2),x)
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?x ?\, ?\(y ?\, ?2 ?\) ?)\)/\(bbn_ks_ddd\1_D0D1D1 KS_func_pass_args_macro \)/g;
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\(y ?\,? ?2 ?\) ?\, ?x ?)\)/\(bbn_ks_ddd\1_D0D1D1 KS_func_pass_args_macro \)/g;

# (x,(z,2)) | ((z,2),x)
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?x ?\, ?\(z ?\, ?2 ?\) ?)\)/\(bbn_ks_ddd\1_D0D2D2 KS_func_pass_args_macro \)/g;
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\(z ?\,? ?2 ?\) ?\, ?x ?)\)/\(bbn_ks_ddd\1_D0D2D2 KS_func_pass_args_macro \)/g;

# (y,(x,2)) | ((x,2),y)
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?y ?\, ?\(x ?\, ?2 ?\) ?)\)/\(bbn_ks_ddd\1_D0D0D1 KS_func_pass_args_macro \)/g;
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\(x ?\,? ?2 ?\) ?\, ?y ?)\)/\(bbn_ks_ddd\1_D0D0D1 KS_func_pass_args_macro \)/g;

# (y,(z,2)) | ((z,2),y)
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?y ?\, ?\(z ?\, ?2 ?\) ?)\)/\(bbn_ks_ddd\1_D1D2D2 KS_func_pass_args_macro \)/g;
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\(z ?\,? ?2 ?\) ?\, ?y ?)\)/\(bbn_ks_ddd\1_D1D2D2 KS_func_pass_args_macro \)/g;


# (z,(x,2)) | ((x,2),z)
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?z ?\, ?\(x ?\, ?2 ?\) ?)\)/\(bbn_ks_ddd\1_D0D0D2 KS_func_pass_args_macro \)/g;
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\(x ?\,? ?2 ?\) ?\, ?z ?)\)/\(bbn_ks_ddd\1_D0D0D2 KS_func_pass_args_macro \)/g;


# (z,(y,2)) | ((y,2),z)
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?z ?\, ?\(y ?\, ?2 ?\) ?)\)/\(bbn_ks_ddd\1_D1D1D2 KS_func_pass_args_macro \)/g;
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\(y ?\,? ?2 ?\) ?\, ?z ?)\)/\(bbn_ks_ddd\1_D1D1D2 KS_func_pass_args_macro \)/g;

# (x,3) | (y,3) | (z,3)
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\( ?x ?\, ?3 ?\) ?)\)/\(bbn_ks_ddd\1_D0D0D0 KS_func_pass_args_macro \)/g;
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\( ?y ?\, ?3 ?\) ?)\)/\(bbn_ks_ddd\1_D1D1D1 KS_func_pass_args_macro \)/g;
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\( ?z ?\, ?3 ?\) ?)\)/\(bbn_ks_ddd\1_D2D2D2 KS_func_pass_args_macro \)/g;

### d^2/d?^2
# (x,2) | (y,2) | (z,2)
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\( ?x ?\, ?2 ?\) ?)\)/\(bbn_ks_dd\1_D0D0 KS_func_pass_args_macro \)/g;
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\( ?y ?\, ?2 ?\) ?)\)/\(bbn_ks_dd\1_D1D1 KS_func_pass_args_macro \)/g;
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\( ?z ?\, ?2 ?\) ?)\)/\(bbn_ks_dd\1_D2D2 KS_func_pass_args_macro \)/g;

# (x,y) | (x,z) | (y,z)
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?x ?\, y ?) ?\)/\(bbn_ks_dd\1_D0D1 KS_func_pass_args_macro \)/g;
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?y ?\, x ?) ?\)/\(bbn_ks_dd\1_D0D1 KS_func_pass_args_macro \)/g;

s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?x ?\, z ?) ?\)/\(bbn_ks_dd\1_D0D2 KS_func_pass_args_macro \)/g;
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?z ?\, x ?) ?\)/\(bbn_ks_dd\1_D0D2 KS_func_pass_args_macro \)/g;

s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?y ?\, z ?) ?\)/\(bbn_ks_dd\1_D1D2 KS_func_pass_args_macro \)/g;
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?z ?\, y ?) ?\)/\(bbn_ks_dd\1_D1D2 KS_func_pass_args_macro \)/g;

# x | y | z
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\, ?(x) ?\)/\(bbn_ks_d\1_D0 KS_func_pass_args_macro \)/g;
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\, ?(y) ?\)/\(bbn_ks_d\1_D1 KS_func_pass_args_macro \)/g;
s/\bDerivative\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\, ?(z) ?\)/\(bbn_ks_d\1_D2 KS_func_pass_args_macro \)/g;


### D:

# (x,y,z)
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?x ?\, ?y ?\, ?z ?)\)/\(bbn_ks_ddd\1_D0D1D2 KS_func_pass_args_macro \)/g;

# (x,(y,2)) | ((y,2),x)
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?x ?\, ?\(y ?\, ?2 ?\) ?)\)/\(bbn_ks_ddd\1_D0D1D1 KS_func_pass_args_macro \)/g;
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\(y ?\,? ?2 ?\) ?\, ?x ?)\)/\(bbn_ks_ddd\1_D0D1D1 KS_func_pass_args_macro \)/g;

# (x,(z,2)) | ((z,2),x)
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?x ?\, ?\(z ?\, ?2 ?\) ?)\)/\(bbn_ks_ddd\1_D0D2D2 KS_func_pass_args_macro \)/g;
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\(z ?\,? ?2 ?\) ?\, ?x ?)\)/\(bbn_ks_ddd\1_D0D2D2 KS_func_pass_args_macro \)/g;

# (y,(x,2)) | ((x,2),y)
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?y ?\, ?\(x ?\, ?2 ?\) ?)\)/\(bbn_ks_ddd\1_D0D0D1 KS_func_pass_args_macro \)/g;
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\(x ?\,? ?2 ?\) ?\, ?y ?)\)/\(bbn_ks_ddd\1_D0D0D1 KS_func_pass_args_macro \)/g;

# (y,(z,2)) | ((z,2),y)
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?y ?\, ?\(z ?\, ?2 ?\) ?)\)/\(bbn_ks_ddd\1_D1D2D2 KS_func_pass_args_macro \)/g;
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\(z ?\,? ?2 ?\) ?\, ?y ?)\)/\(bbn_ks_ddd\1_D1D2D2 KS_func_pass_args_macro \)/g;


# (z,(x,2)) | ((x,2),z)
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?z ?\, ?\(x ?\, ?2 ?\) ?)\)/\(bbn_ks_ddd\1_D0D0D2 KS_func_pass_args_macro \)/g;
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\(x ?\,? ?2 ?\) ?\, ?z ?)\)/\(bbn_ks_ddd\1_D0D0D2 KS_func_pass_args_macro \)/g;


# (z,(y,2)) | ((y,2),z)
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?z ?\, ?\(y ?\, ?2 ?\) ?)\)/\(bbn_ks_ddd\1_D1D1D2 KS_func_pass_args_macro \)/g;
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\(y ?\,? ?2 ?\) ?\, ?z ?)\)/\(bbn_ks_ddd\1_D1D1D2 KS_func_pass_args_macro \)/g;

# (x,3) | (y,3) | (z,3)
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\( ?x ?\, ?3 ?\) ?)\)/\(bbn_ks_ddd\1_D0D0D0 KS_func_pass_args_macro \)/g;
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\( ?y ?\, ?3 ?\) ?)\)/\(bbn_ks_ddd\1_D1D1D1 KS_func_pass_args_macro \)/g;
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\( ?z ?\, ?3 ?\) ?)\)/\(bbn_ks_ddd\1_D2D2D2 KS_func_pass_args_macro \)/g;

### d^2/d?^2
# (x,2) | (y,2) | (z,2)
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\( ?x ?\, ?2 ?\) ?)\)/\(bbn_ks_dd\1_D0D0 KS_func_pass_args_macro \)/g;
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\( ?y ?\, ?2 ?\) ?)\)/\(bbn_ks_dd\1_D1D1 KS_func_pass_args_macro \)/g;
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?\( ?z ?\, ?2 ?\) ?)\)/\(bbn_ks_dd\1_D2D2 KS_func_pass_args_macro \)/g;

# (x,y) | (x,z) | (y,z)
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?x ?\, ?y ?) ?\)/\(bbn_ks_dd\1_D0D1 KS_func_pass_args_macro \)/g;
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?y ?\, ?x ?) ?\)/\(bbn_ks_dd\1_D0D1 KS_func_pass_args_macro \)/g;

s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?x ?\, ?z ?) ?\)/\(bbn_ks_dd\1_D0D2 KS_func_pass_args_macro \)/g;
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?z ?\, ?x ?) ?\)/\(bbn_ks_dd\1_D0D2 KS_func_pass_args_macro \)/g;

s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?y ?\, ?z ?) ?\)/\(bbn_ks_dd\1_D1D2 KS_func_pass_args_macro \)/g;
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\,( ?z ?\, ?y ?) ?\)/\(bbn_ks_dd\1_D1D2 KS_func_pass_args_macro \)/g;

# x | y | z
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\, ?(x) ?\)/\(bbn_ks_d\1_D0 KS_func_pass_args_macro \)/g;
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\, ?(y) ?\)/\(bbn_ks_d\1_D1 KS_func_pass_args_macro \)/g;
s/\bD\b\(bbn_ks_([[:alnum:]_]+)\(([xyz\, ]+)\)\, ?(z) ?\)/\(bbn_ks_d\1_D2 KS_func_pass_args_macro \)/g;


## change to variable for optimization
s/(bbn_ks_d[XYZ]_D[012]) +KS_func_pass_args_macro/\(\1_\)/g;


















