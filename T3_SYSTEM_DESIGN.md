# Implementation details

Since we plan on testing this feature using an isometric contraction of parallel fibred muscle, we can modify the parameter file `IC_parameters.prm` in the `examples` folder.

## Parameter files

In `IC_parameters.prm`, we add the following parameters to the `Materials` subsection:
```
  # --------------- Fat properties ------------------

  # Bulk modulus fat [Pa]
  set Bulk modulus fat = 1.0e+07

  # Fat "fudge" factor [nondimensional].
  # Ideally this should be 1.
  set Fat factor = 1.0

  # Constants in Yeoh strain-energy function
  set Fat constant 1 = 323.91
  set Fat constant 2 = 5163.1
  set Fat constant 3 = -3872.9

  # Fat fraction
  set Fat fraction = 0.02
```

After testing, all other `.prm` files in the Flexodeal Lite repository should be updated to reflect the changes. Otherwise, the rest of the examples may not run since the code now expect new parameters that other files did not define.

## Main C++ files

In `flexodeal.cc`:

1. Modify the parser, i.e. the `MuscleProperties` struct. Add six new variables, including functions to _declare_ and _get_ the parameters.

2. Modify the class that defines the material properties of muscle according to the provided equations, i.e. `Muscle_Tissues_Three_Field`: 
    - add 6 new inputs to the constructor, 5 of them will become new class members. 
    - The input `kappa_fat` is used solely to homogenize the value of `kappa_muscle` and is not needed as a new class member.
    - The fictitious Kirchhoff stress, defined in `get_tau_bar()`, needs to be homogenized as explained in the accompanying document. This requires to define the new fictitious Kirchhoff stress for fat, call it `get_tau_fat_bar()`.
    - The fictitious elasticity tensor, as defined in `get_c_bar()` needs to be homogenized in the same way. Similarly, this also requires to define the new fictitious elasticity tensor for fat, call it `get_c_fat_bar()`.

3. Since the constructor for `Muscle_Tissues_Three_Field` was modified, we need to update the definition of any object of this class that is created in the code. Luckily, this is only needed once in `PointHistory::setup_lqp()` when the `material` pointer is created.

