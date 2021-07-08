```@meta
DocTestSetup = :(using PYMethod)
```

# Finite Element Analysis

Simulate lateral behavior of pile based on FEM using beam elements.

## STEP 1: Create pile object

If the bottom and top coordinates are ``0\,\mathrm{m}`` and ``1\,\mathrm{m}``, respectively,
then the object can be constructed as

```jldoctest pile
julia> pile = FEPileModel(0, 1, 20);
```

The last argument specifies the number of beam elements of the pile.
Thus, every element has the length of ``0.05\,\mathrm{m}``.
Check `pile.coordinates` for coordinates of all nodes.

## STEP 2: Initialize parameters

* `pile.E`: Young's modulus of pile
* `pile.I`: Second moment of area of pile
* `pile.D`: Diameter of pile

Note that each parameter is actually a vector with `length` of number of nodes.
For example, `pile.E[1]` and `pile.E[end]` are the Young's moduli at top and bottom of pile, respectively.
If the values change from one node to the other, they will be linearly interpolated for the inner elements.

In the following example, Young's modulus is ``2.0\times10^{11}\,\mathrm{N/m^2}``,
second moment of area is ``3.07\times10^{-7}\,\mathrm{m^4}``,
and diameter is ``0.05\,\mathrm{m}`` for all nodes.

```jldoctest pile
julia> pile.E .= 2.0e11;

julia> pile.I .= 3.07e-7;

julia> pile.D .= 0.05;
```

## STEP 3: Setup boundary conditions

* `pile.u`: Lateral displacement
* `pile.θ`: Angle of deflection (where `θ` can be typed by `\theta<tab>`)
* `pile.Fext`: Lateral force
* `pile.Mext`: Moment

Similar to parameters, above values are also vectors with `length` of number of nodes.
To apply the lateral force with ``100\,\mathrm{N}`` at the top of the pile, simply set the value as

```jldoctest pile
julia> pile.Fext[1] = 100;
```

## STEP 4: Setup p-y curves

* `pile.pycurves`: `pycurve(y, z) -> p`

P-y curves also needs to be setup on each nodes.
The object must be the function which has the lateral displacement `y` and
vertical coordinate `z` as arguments and returns the earth pressure `p`.
`pycurve(y, z) = 0` is used by defalut.

## STEP 5: Solve the problem

```jldoctest pile
julia> solve!(pile);

julia> pile.u
21-element view(::PYMethod.FEMVector{Float64}, 1:2:41) with eltype Float64:
 0.0005428881650347354
 0.0005022054831674393
 0.00046172638436203214
 0.00042165445168040216
 0.00038219326818443735
 0.00034354641693602555
 0.00030591748099705464
 0.00026951004342941267
 0.00023452768729498742
 0.00020117399565566698
 ⋮
 0.00011292073832721299
 8.811753528719052e-5
 6.596091205171245e-5
 4.665445168266666e-5
 3.0401737241941133e-5
 1.7406351791423797e-5
 7.871878393002527e-6
 2.0019001085652914e-6
 0.0
```
