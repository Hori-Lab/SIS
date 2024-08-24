# Restraint option

Restraint can be added by specifying a `restraint` filename in `[Files.In]`. 

```
[Files]
    [Files.In]
    restraint = "restraint.txt"
```

There is no other control parameters required in the input file. 

## Restraint file format

The restraint file is a sepeparete text file, containing lines, each of which specifies individual restraint. Blank lines are allowed. Any text preceded by # is treated as a comment. See below for examples.

## Sigmoid potential

You can apply a restraint to a pair of particles using a sigmoid-type function, 

$$\begin{split}
U(r) = -\frac{1}{2}\varepsilon\left[1 - \tanh{\left(\frac{r - d}{s}\right)}\right]
\end{split}
$$

where $r$ is the distance between the particles, $\varepsilon$ is the energy scale (i.e. strength of the restraint), $d$ is the parameter that detemines the range, $s$ is the parameter that controls the smoothness of the potential.

In the restraint file, this should be specified with a keyword `Sigmoid` with the parameters.


```
#                particle 1 ID  particle 2 ID    epsilon     d         s     r_cut
Sigmoid                1            47            5.0       5.0       1.8     12.0
Sigmoid               48             2            5.0       5.0       1.8     12.0
```

`r_cut` is the cut-off distance where the potential decays to zero. 
 

The next example, `Sigmoid-to-bead` is a variation of the Sigmoid restraint. The only difference is that the target particle (subject ID) is restrained around the reference particle, while the reference particle does not feel the force. In other words, the position of the reference particle is used to constrain the subject particle, but the reference particle itself is not affected at all. The potential function and parameters are the same as for `Sigmoid`.

```
#                 subject ID     reference ID    epsilon     d         s     r_cut
Sigmoid-to-bead        1            47            5.0       5.0       1.8     12.0
Sigmoid-to-bead       48             2            5.0       5.0       1.8     12.0
```
