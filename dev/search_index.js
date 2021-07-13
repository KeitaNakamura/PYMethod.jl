var documenterSearchIndex = {"docs":
[{"location":"api/ChangEquation/#Chang's-equation","page":"Chang's equation","title":"Chang's equation","text":"","category":"section"},{"location":"api/ChangEquation/","page":"Chang's equation","title":"Chang's equation","text":"Modules = [PYMethod]\nOrder = [:type, :function, :macro]\nPages = [\"ChangEquation.md\"]","category":"page"},{"location":"api/ChangEquation/","page":"Chang's equation","title":"Chang's equation","text":"Modules = [PYMethod]\nOrder = [:type, :function, :macro]\nPages = [\"ChangEquation.jl\"]","category":"page"},{"location":"api/ChangEquation/#PYMethod.ChangEquation","page":"Chang's equation","title":"PYMethod.ChangEquation","text":"ChangEquation(bottom, top; parameters...)\n\nParameters\n\nz_0: height of ground surface (z_0 = 0 by default)\nF_t: Lateral load at pile head\nM_t: Moment at pile head (M_t = 0 by default)\nD: Diameter of pile\nE: Young's modulus of pile\nI: Second moment of area of pile\nk: Modulus of subgrade reaction\n\nExamples\n\njulia> eq = ChangEquation(-19, 1; F_t = 10, D = 0.6, E = 2e8, I = 0.0002507, k = 3750);\n\n\n\n\n\n","category":"type"},{"location":"api/ChangEquation/#PYMethod.calculate_deflection-Tuple{ChangEquation, Real}","page":"Chang's equation","title":"PYMethod.calculate_deflection","text":"calculate_deflection(::ChangEquation, z)\n\nCalculate deflection at height z.\n\n\n\n\n\n","category":"method"},{"location":"api/ChangEquation/#PYMethod.calculate_moment-Tuple{ChangEquation, Real}","page":"Chang's equation","title":"PYMethod.calculate_moment","text":"calculate_moment(::ChangEquation, z)\n\nCalculate moment at height z.\n\n\n\n\n\n","category":"method"},{"location":"api/ChangEquation/#PYMethod.calculate_shearforce-Tuple{ChangEquation, Real}","page":"Chang's equation","title":"PYMethod.calculate_shearforce","text":"calculate_shearforce(::ChangEquation, z)\n\nCalculate shear force at height z.\n\n\n\n\n\n","category":"method"},{"location":"api/FEM/#FEM","page":"FEM","title":"FEM","text":"","category":"section"},{"location":"api/FEM/","page":"FEM","title":"FEM","text":"Modules = [PYMethod]\nOrder = [:type, :function, :macro]\nPages = [\"FEM.md\"]","category":"page"},{"location":"api/FEM/","page":"FEM","title":"FEM","text":"Modules = [PYMethod]\nOrder = [:type, :function, :macro]\nPages = [\"fem.jl\"]","category":"page"},{"location":"api/FEM/#PYMethod.FEPileModel-Tuple{Real, Real, Int64}","page":"FEM","title":"PYMethod.FEPileModel","text":"FEPileModel(bottom::Real, top::Real, nelements::Int) -> pile\n\nConstruct an object of the finite element model to simulate lateral behavior of pile. The ith Beam element can be accessed by pile[i]. The following values are vectors storing each nodal value.\n\nParameters\n\npile.E: Young's modulus of pile\npile.I: Second moment of area of pile\npile.D: Diameter of pile\n\nBoundary conditions\n\npile.u: Lateral displacement\npile.θ: Angle of deflection (where θ can be typed by \\theta<tab>)\npile.Fext: External lateral force\npile.Mext: External moment\n\nThe above vectors will be updated after solve! the problem.\n\nNodal variables\n\npile.coordinates: Coordinates of nodes\npile.F: Internal lateral force\npile.M: Internal moment\n\nThe internal lateral force and moment vectors will be updated after solve! the problem.\n\n\n\n\n\n","category":"method"},{"location":"api/FEM/#PYMethod.clear_boundary_conditions!-Tuple{FEPileModel}","page":"FEM","title":"PYMethod.clear_boundary_conditions!","text":"clear_boundary_conditions!(::FEPileModel)\n\nClear boundary conditions. The parameters and p-y curves are remained.\n\n\n\n\n\n","category":"method"},{"location":"api/FEM/#PYMethod.mass_matrix-Tuple{PYMethod.Beam}","page":"FEM","title":"PYMethod.mass_matrix","text":"mass_matrix(::Beam)\n\nConstruct mass matrix from Beam.\n\n\n\n\n\n","category":"method"},{"location":"api/FEM/#PYMethod.solve!-Union{Tuple{FEPileModel{T}}, Tuple{T}} where T","page":"FEM","title":"PYMethod.solve!","text":"solve!(::FEPileModel)\n\nSolve the finite element problem.\n\n\n\n\n\n","category":"method"},{"location":"api/FEM/#PYMethod.stiffness_matrix-Tuple{PYMethod.Beam}","page":"FEM","title":"PYMethod.stiffness_matrix","text":"stiffness_matrix(::Beam)\n\nConstruct element stiffness matrix from Beam.\n\n\n\n\n\n","category":"method"},{"location":"manual/FEM/#Finite-Element-Analysis","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"","category":"section"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"Simulate lateral behavior of pile based on FEM using beam elements.","category":"page"},{"location":"manual/FEM/#Overview-of-finite-element-formulation","page":"Finite Element Analysis","title":"Overview of finite element formulation","text":"","category":"section"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"The differential equation for the deflection curve of a beam on an elastic foundation is","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"EIv(x) + p(xv) = 0","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"where E is the Young's modulus, I is the second moment of area, v is the deflection, p is the pressure from the foundation. Multiplying the test function w and integral the equation over the domain reduces the form","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"int w left EIv(x) + p(xv) right mathrmdx = 0","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"Applying the divergence theorem two times leads following weak form","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"int w EIv(x) mathrmdx + int w p(xv) mathrmdx =\n    int_Gamma_S w barS mathrmdGamma + int_Gamma_M w barM mathrmdGamma","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"where barS and barM are the shear force and the bending moment on the boundary. Using the shape function N, we interpolate v, w and p from nodal values v_i, w_i and p_i, respectively, as","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"v = sum_i N_i v_iquad\nw = sum_i N_i w_iquad\np = sum_i N_i p_i","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"Thus the weak form can be discretized as","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"K_ij v_j(x) + M_ij p_j(xv) = barF_i","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"where","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"K_ij = int EI N_i N_j mathrmdxquad M_ij = int N_i N_j mathrmdxquad\n    barF_i = int_Gamma_S N_i barS mathrmdGamma + int_Gamma_M N_i barM mathrmdGamma","category":"page"},{"location":"manual/FEM/#How-to-simulate","page":"Finite Element Analysis","title":"How to simulate","text":"","category":"section"},{"location":"manual/FEM/#STEP-1:-Create-pile-object","page":"Finite Element Analysis","title":"STEP 1: Create pile object","text":"","category":"section"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"If the bottom and top coordinates are 0mathrmm and 20mathrmm, respectively, then the object can be constructed as","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"using PYMethod # hide\npile = FEPileModel(0, 20, 200)\nnothing # hide","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"The last argument specifies the number of beam elements of the pile. Thus, every element has the length of 01mathrmm. Check pile.coordinates for coordinates of all nodes.","category":"page"},{"location":"manual/FEM/#STEP-2:-Initialize-parameters","page":"Finite Element Analysis","title":"STEP 2: Initialize parameters","text":"","category":"section"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"pile.E: Young's modulus of pile\npile.I: Second moment of area of pile\npile.D: Diameter of pile","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"Note that each parameter is actually a vector with length of number of nodes. For example, pile.E[1] and pile.E[end] are the Young's moduli at top and bottom of pile, respectively. If the values change from one node to the other, they will be linearly interpolated for the inner elements.","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"In the following example, Young's modulus is 20times10^11mathrmNm^2, second moment of area is 307times10^-7mathrmm^4, and diameter is 005mathrmm for all nodes.","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"pile.E .= 2.0e11\npile.I .= 3.07e-7\npile.D .= 0.05\nnothing # hide","category":"page"},{"location":"manual/FEM/#STEP-3:-Setup-boundary-conditions","page":"Finite Element Analysis","title":"STEP 3: Setup boundary conditions","text":"","category":"section"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"pile.u: Lateral displacement\npile.θ: Angle of deflection (where θ can be typed by \\theta<tab>)\npile.Fext: Lateral force\npile.Mext: Moment","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"Similar to parameters, above values are also vectors with length of number of nodes. To apply the lateral force with 100mathrmN at the top of the pile, simply set the value as","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"pile.Fext[1] = 100\nnothing # hide","category":"page"},{"location":"manual/FEM/#STEP-4:-Setup-p-y-curves","page":"Finite Element Analysis","title":"STEP 4: Setup p-y curves","text":"","category":"section"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"pile.pycurves: pycurve(y, z) -> p","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"P-y curves also needs to be setup on each nodes. The object must be the function which has the lateral displacement y and vertical coordinate z as arguments and returns the earth pressure p. pycurve(y, z) = 0 is used by default.","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"k = 50e3\npile.pycurves .= pycurve(y, z) = z > 19 ? 0 : k*y\nnothing # hide","category":"page"},{"location":"manual/FEM/#STEP-5:-Solve-the-problem","page":"Finite Element Analysis","title":"STEP 5: Solve the problem","text":"","category":"section"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"solve!(pile)\npile.u","category":"page"},{"location":"manual/FEM/#STEP-6:-Visualize-the-results","page":"Finite Element Analysis","title":"STEP 6: Visualize the results","text":"","category":"section"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"using Plots\neq = ChangEquation(0, 20; z_0 = 19, E = 2e11, I = 3.07e-7, D = 0.05, F_t = 100, k = 50e3)\nplot(pile; label = \"FEM\")\nplot!(eq; label = \"Chang's equation\")\nplot!(legend = :bottomleft); nothing # hide\nsavefig(\"fem_chang.svg\"); nothing # hide","category":"page"},{"location":"manual/FEM/","page":"Finite Element Analysis","title":"Finite Element Analysis","text":"(Image: )","category":"page"},{"location":"#PYMethod","page":"Home","title":"PYMethod","text":"","category":"section"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"pkg> add https://github.com/KeitaNakamura/PYMethod.jl.git","category":"page"}]
}
