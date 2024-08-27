# Periodic Boundary Condition (PBC)

## Setting up PBC

You can enable the periodic boundary condition by adding `[PBC_box]` block in the input file, as shown below.

```
[PBC_box]
size = [500.0, 500.0, 500.0]
```

In this example, a periodic bounding box with one side of 500 Å is set up. The rectangular box is placed so that its centre is at the coordinates (0, 0, 0). This means that in this example, the range of each dimension is between -250 Å and 250 Å. At present, only rectangular boxes are supported.

Currently, there is no option to set up a box with physical boundaries (i.e. without periodic boundaries). This may be added in the future if required.

### PBC on restart

Box size information is written to the `rst` file (from September 2024 version onwards).

If `[PBC_box]` is present in the input file and the size information is present in the restart file (`.rst`), the size loaded from the restart file will be used. To prevent this (e.g. if you want to use a different dimension for the restart simulation), you can set `ignore_rst = true`, which will cause the program to ignore the size information in the restart file.

```
[PBC_box]
size = [500.0, 500.0, 500.0]
ignore_rst = true
           ## false by default
```

If there is no `[PBC_box]` in the input file, there will be no PBC set even if the restart file contains the size information, so just remove the `[PBC_box]` from the input file if you need to switch simulations from PBC to no PBC for some reason.

### Deprecated format

The format below for specifying the box size is deprecated and will be removed in a future release.

```
[PBC_box]
x = 500.0
y = 500.0
z = 500.0
```

## Resize option

"Resize" option can be used to gradually compress or expand the periodic boundary box during the simulation.

```
[PBC_box]
size = [500.0, 500.0, 500.0]

resize_step = 1000
resize_change = [ -0.1,  -0.1,  -0.1]
resize_target = [400.0, 400.0, 400.0]
```

In the example above, the box size starts from [500, 500, 500] (unless a restart file is used and has PBC), shrinking by [-0.1, -0.1, -0.1] for every 1000 steps, until it reaches the target size [400, 400, 400].

### Deprecated format for resize

The following is an old format using `[Vaiable_box]`. This option is deprecated and will be removed in a future version. `target` cannot be set in this way. Use `resize` in `[PBC_box]` as above.

```
[Variable_box]
step = 1000
change_x = -0.1
change_y = -0.1
change_z = -0.1
```