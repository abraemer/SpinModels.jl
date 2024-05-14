var documenterSearchIndex = {"docs":
[{"location":"quick-ref/#Quick-reference","page":"Quick reference","title":"Quick reference","text":"","category":"section"},{"location":"quick-ref/#Terms","page":"Quick reference","title":"Terms","text":"","category":"section"},{"location":"quick-ref/","page":"Quick reference","title":"Quick reference","text":"Add terms together to define the structure of your Hamiltonian and multiply with constants and coupling matrices (only once!) to set up a concrete realization.","category":"page"},{"location":"quick-ref/","page":"Quick reference","title":"Quick reference","text":"X(h), Y(h) and Z(h) stand for sum_i h_i sigma_alpha^(i) where alpha is x, y or z.\nXX(J), YY(J) and ZZ(J) stand for sum_ij J_ij sigma_alpha^(i)sigma_alpha^(j) where alpha is x, y or z.\nFlipFlop(J) (or equivalently Hopp(J)) is the same as 0.5*J*(XX()+YY())\nXXZ(Δ, J) is the same as 2*Hopp(J)+Δ*ZZ(J) which is the same as XX(J)+YY(J)+Δ*ZZ(J)","category":"page"},{"location":"quick-ref/#Geometries","page":"Quick reference","title":"Geometries","text":"","category":"section"},{"location":"quick-ref/","page":"Quick reference","title":"Quick reference","text":"Current implemented are (N always denotes the total number of spins):","category":"page"},{"location":"quick-ref/","page":"Quick reference","title":"Quick reference","text":"Chain(N)\nBox(N, dims)\nNoisyChain(N, σ; spacing=1)\nPartiallyFilledChain(numspins, numsites; spacing=1)\nBlockaded(geometry; retries=1000, blockade=1.0)\nPBC(geometry)\nNN(geometry, k=1)","category":"page"},{"location":"quick-ref/","page":"Quick reference","title":"Quick reference","text":"note: Note\nThe order of Blockaded, PBC and NN is irrelevant. They sort them themselves.","category":"page"},{"location":"quick-ref/#Interactions","page":"Quick reference","title":"Interactions","text":"","category":"section"},{"location":"quick-ref/","page":"Quick reference","title":"Quick reference","text":"Currently implemented are","category":"page"},{"location":"quick-ref/","page":"Quick reference","title":"Quick reference","text":"ConstantInteraction(value)\nPowerLaw(α)\nNN(interaction, k=1)","category":"page"},{"location":"quick-ref/","page":"Quick reference","title":"Quick reference","text":"note: Note\nThere are QoL overloads on Geometry that apply ConstantInteraction(1) automatically. Thus you should probably not need to use it directly. Please give feedback if that does unexpected things.","category":"page"},{"location":"quick-ref/#Note-about-nearest-neighbor-(NN)","page":"Quick reference","title":"Note about nearest neighbor (NN)","text":"","category":"section"},{"location":"quick-ref/","page":"Quick reference","title":"Quick reference","text":"For isotropic interactions (i.e. all currently implemented ones) it does not matter whether you apply NN to the interaction or the geometry. In principle for anisotropic interaction there will be a subtle difference:","category":"page"},{"location":"quick-ref/","page":"Quick reference","title":"Quick reference","text":"NN on geometry is based on distances and will remove all but the ksmallest distances\nNN on interactions is based on coupling strength and will remove all but the kstrongest couplings","category":"page"},{"location":"full-docs/geometry/#Full-docs-for-Geometry","page":"Full docs for Geometry","title":"Full docs for Geometry","text":"","category":"section"},{"location":"full-docs/geometry/#Types","page":"Full docs for Geometry","title":"Types","text":"","category":"section"},{"location":"full-docs/geometry/","page":"Full docs for Geometry","title":"Full docs for Geometry","text":"Chain\nBox\nNoisyChain\nPartiallyFilledChain\nBlockaded\nPBC\nNN(::SpinModels.Geometry)","category":"page"},{"location":"full-docs/geometry/#SpinModels.Chain","page":"Full docs for Geometry","title":"SpinModels.Chain","text":"Chain(L; spacing=1)\n\nRepresents a 1D chain of L spins with regular spacing. Spins are at x_i = spacing*i\n\n\n\n\n\n","category":"type"},{"location":"full-docs/geometry/#SpinModels.Box","page":"Full docs for Geometry","title":"SpinModels.Box","text":"Box(N, dims)\n\nRepresents a box with side-lengths given by dims. N spins are randomly positioned within that box.\n\n\n\n\n\n","category":"type"},{"location":"full-docs/geometry/#SpinModels.NoisyChain","page":"Full docs for Geometry","title":"SpinModels.NoisyChain","text":"NoisyChain(L, σ; spacing=1)\n\nRepresents a 1D chain of L spins with random noise on top of regular spaced positions. Spins are at x_i = spacing*i + delta where delta sim mathcalU(-fracsigma2fracsigma2)\n\n\n\n\n\n","category":"type"},{"location":"full-docs/geometry/#SpinModels.PartiallyFilledChain","page":"Full docs for Geometry","title":"SpinModels.PartiallyFilledChain","text":"PartiallyFilledChain(numspins, numsites; spacing=1)\n\nRepresents a 1D chain of L sites where N spins are placed randomly. Sites are spacing units apart.\n\n\n\n\n\n","category":"type"},{"location":"full-docs/geometry/#SpinModels.Blockaded","page":"Full docs for Geometry","title":"SpinModels.Blockaded","text":"Blockaded(geometry; retries=1000, blockade=1.0)\n\nModifies the underlying geometry to respect the blockade condition namely all spins need to be at least blockade units apart. Positions are resampled at most retries times.\n\nnote: Note\nIf you set retries=0, resampling will not stop until a valid configuration is found.\n\n\n\n\n\n","category":"type"},{"location":"full-docs/geometry/#SpinModels.PBC","page":"Full docs for Geometry","title":"SpinModels.PBC","text":"PBC(geometry)\n\nModify the underlying geometry to respect periodic boundary conditions.\n\n\n\n\n\n","category":"type"},{"location":"full-docs/geometry/#SpinModels.NN-Tuple{SpinModels.Geometry}","page":"Full docs for Geometry","title":"SpinModels.NN","text":"NN(geometry, k=1)\n\nModify the geometry to only keep the k closest connections from each spin. Please see nearest_neighbor_from_distances! for further information.\n\nnote: Note\nThis is not identical to NN(::SpinModels.Interaction) if the interaction is anisotropic.\n\n\n\n\n\n","category":"method"},{"location":"full-docs/geometry/#Functions","page":"Full docs for Geometry","title":"Functions","text":"","category":"section"},{"location":"full-docs/geometry/","page":"Full docs for Geometry","title":"Full docs for Geometry","text":"nspins\npositions\ndistance_matrix\nnearest_neighbor_from_distances\nnearest_neighbor_from_distances!","category":"page"},{"location":"full-docs/geometry/#SpinModels.nspins","page":"Full docs for Geometry","title":"SpinModels.nspins","text":"nspins(term)\nnspins(sum_of_terms)\n\nReturn the number of spins in the underlying geometry if that information was already provided.\n\n\n\n\n\nnspins(geometry)\n\nReturn number of spins the geometry.\n\n\n\n\n\n","category":"function"},{"location":"full-docs/geometry/#SpinModels.positions","page":"Full docs for Geometry","title":"SpinModels.positions","text":"positions(geometry)\npositions(geometry; rng)\n\nGenerate the positions of the spin in the geometry. Disordered geometries allow for setting an RNG.\n\n\n\n\n\n","category":"function"},{"location":"full-docs/geometry/#SpinModels.distance_matrix","page":"Full docs for Geometry","title":"SpinModels.distance_matrix","text":"distance_matrix(geometry[, positions]; rng)\n\nCompute the spin-spin distances. Positions are generated by the geometry via positions or provided directly. Disordered geometries allow for setting an RNG.\n\n\n\n\n\n","category":"function"},{"location":"full-docs/geometry/#SpinModels.nearest_neighbor_from_distances","page":"Full docs for Geometry","title":"SpinModels.nearest_neighbor_from_distances","text":"nearest_neighbor_from_distances(distance_matrix, k=1)\n\nNon-destructive form of nearest_neighbor_from_distances!. See also: nearest_neighbor_from_interactions!\n\n\n\n\n\n","category":"function"},{"location":"full-docs/geometry/#SpinModels.nearest_neighbor_from_distances!","page":"Full docs for Geometry","title":"SpinModels.nearest_neighbor_from_distances!","text":"nearest_neighbor_from_distances!(distance_matrix, k=1)\n\nSet all but the k smallest distances from each column to Inf. k only counts different values! Thus k=1 keeps only the closest distance per spin but allows to have the same distance multiple times.\n\nnote: Note\nThe distance to itself is always ignored.\n\nnote: Note\nThe resulting matrix won't necessarily be symmetric if nearest neighbors are ill-defined. See also: nearest_neighbor_from_interactions!\n\n\n\n\n\n","category":"function"},{"location":"full-docs/interaction/#Full-docs-for-Interaction","page":"Full docs for Interaction","title":"Full docs for Interaction","text":"","category":"section"},{"location":"full-docs/interaction/#Types","page":"Full docs for Interaction","title":"Types","text":"","category":"section"},{"location":"full-docs/interaction/","page":"Full docs for Interaction","title":"Full docs for Interaction","text":"ConstantInteraction\nPowerLaw\nNN(::SpinModels.Interaction)\nNN","category":"page"},{"location":"full-docs/interaction/#SpinModels.ConstantInteraction","page":"Full docs for Interaction","title":"SpinModels.ConstantInteraction","text":"ConstantInteraction(value)\n\nRepresents a simple constant valued interaction.\n\nThere are quality of life overloads, that removes the need to use ConstantInteraction(1)(geometry) such that you can simply use the geometry like an interaction matrix in most scenarios. Please complain if this produces unintuitive behavior.\n\n\n\n\n\n","category":"type"},{"location":"full-docs/interaction/#SpinModels.PowerLaw","page":"Full docs for Interaction","title":"SpinModels.PowerLaw","text":"PowerLaw(α)\n\nRepresents couplings decaying as x_i-x_j^-alpha\n\n\n\n\n\n","category":"type"},{"location":"full-docs/interaction/#SpinModels.NN-Tuple{SpinModels.Interaction}","page":"Full docs for Interaction","title":"SpinModels.NN","text":"NN(interaction, k=1)\n\nModify the interaction to only keep the k strongest couplings from each spin. Please see nearest_neighbor_from_interactions! for further information.\n\nnote: Note\nThis is not identical to NN(::SpinModels.Geometry) if the interaction is anisotropic.\n\n\n\n\n\n","category":"method"},{"location":"full-docs/interaction/#SpinModels.NN","page":"Full docs for Interaction","title":"SpinModels.NN","text":"NN(geometry, k=1)\n\nModify the geometry to only keep the k closest connections from each spin. Please see nearest_neighbor_from_distances! for further information.\n\nnote: Note\nThis is not identical to NN(::SpinModels.Interaction) if the interaction is anisotropic.\n\n\n\n\n\nNN(interaction, k=1)\n\nModify the interaction to only keep the k strongest couplings from each spin. Please see nearest_neighbor_from_interactions! for further information.\n\nnote: Note\nThis is not identical to NN(::SpinModels.Geometry) if the interaction is anisotropic.\n\n\n\n\n\n","category":"function"},{"location":"full-docs/interaction/#Functions","page":"Full docs for Interaction","title":"Functions","text":"","category":"section"},{"location":"full-docs/interaction/","page":"Full docs for Interaction","title":"Full docs for Interaction","text":"interaction_matrix\nnearest_neighbor_from_interactions\nnearest_neighbor_from_interactions!","category":"page"},{"location":"full-docs/interaction/#SpinModels.interaction_matrix","page":"Full docs for Interaction","title":"SpinModels.interaction_matrix","text":"interaction_matrix(interaction, geometry[, positions]; rng)\n(::Interaction)(geometry[, positions]; rng)\n\nCompute the corresponding interaction matrix from interaction and geometry. Positions are generated by the geometry via positions or provided directly. Disordered geometries allow for setting an RNG.\n\n\n\n\n\n","category":"function"},{"location":"full-docs/interaction/#SpinModels.nearest_neighbor_from_interactions","page":"Full docs for Interaction","title":"SpinModels.nearest_neighbor_from_interactions","text":"nearest_neighbor_from_interactions(interaction_matrix, k=1)\n\nNon-destructive form of nearest_neighbor_from_interactions!. See also: nearest_neighbor_from_distances!\n\n\n\n\n\n","category":"function"},{"location":"full-docs/interaction/#SpinModels.nearest_neighbor_from_interactions!","page":"Full docs for Interaction","title":"SpinModels.nearest_neighbor_from_interactions!","text":"nearest_neighbor_from_interactions!(interaction_matrix, k=1)\n\nSet all but the k strongest couplings from each column to 0. k only counts different values! Thus k=1 keeps only the strongest couplings per spin but allows to have the same value multiple times.\n\nnote: Note\nThe resulting matrix won't necessarily be symmetric if nearest neighbors are ill-defined. See also: nearest_neighbor_from_distances!\n\n\n\n\n\n","category":"function"},{"location":"full-docs/hamiltonian/#Full-docs-for-Hamiltonian","page":"Full docs for Hamiltonian","title":"Full docs for Hamiltonian","text":"","category":"section"},{"location":"full-docs/hamiltonian/#Types","page":"Full docs for Hamiltonian","title":"Types","text":"","category":"section"},{"location":"full-docs/hamiltonian/","page":"Full docs for Hamiltonian","title":"Full docs for Hamiltonian","text":"These Term types may be added together freely. Multiplication with scalars works as expected. You can only impose spatial information (via multiplication with a vector or matrix respective of the type of term) once! You can only add together compatible terms in terms of spatial structure.","category":"page"},{"location":"full-docs/hamiltonian/","page":"Full docs for Hamiltonian","title":"Full docs for Hamiltonian","text":"X\nY\nZ\nXX\nYY\nZZ\nHopp\nFlipFlop\nXXZ","category":"page"},{"location":"full-docs/hamiltonian/#SpinModels.X","page":"Full docs for Hamiltonian","title":"SpinModels.X","text":"X(J=1)\n\nRepresents a sigma_x^(i) term. Multiply with a vector h to get the term sum_i h_i sigma_x^(i) \n\n\n\n\n\n","category":"function"},{"location":"full-docs/hamiltonian/#SpinModels.Y","page":"Full docs for Hamiltonian","title":"SpinModels.Y","text":"Y(J=1)\n\nRepresents a sigma_y^(i) term. Multiply with a vector h to get the term sum_i h_i sigma_y^(i) \n\n\n\n\n\n","category":"function"},{"location":"full-docs/hamiltonian/#SpinModels.Z","page":"Full docs for Hamiltonian","title":"SpinModels.Z","text":"Z(J=1)\n\nRepresents a sigma_z^(i) term. Multiply with a vector h to get the term sum_i h_i sigma_z^(i) \n\n\n\n\n\n","category":"function"},{"location":"full-docs/hamiltonian/#SpinModels.XX","page":"Full docs for Hamiltonian","title":"SpinModels.XX","text":"XX(J=1)\n\nRepresents a sigma_x^(i) sigma_x^(j) term. Multiply with a matrix J to get the term sum_ineq j J_ij sigma_x^(i) sigma_x^(j) \n\nnote: Note\nTwo-body terms do 'double-counting' of the couplings meaning the sum runs over i,j independently and not i < j.\n\n\n\n\n\n","category":"function"},{"location":"full-docs/hamiltonian/#SpinModels.YY","page":"Full docs for Hamiltonian","title":"SpinModels.YY","text":"YY(J=1)\n\nRepresents a sigma_y^(i) sigma_y^(j) term. Multiply with a matrix J to get the term sum_ineq j J_ij sigma_y^(i) sigma_y^(j) \n\nnote: Note\nTwo-body terms do 'double-counting' of the couplings meaning the sum runs over i,j independently and not i < j.\n\n\n\n\n\n","category":"function"},{"location":"full-docs/hamiltonian/#SpinModels.ZZ","page":"Full docs for Hamiltonian","title":"SpinModels.ZZ","text":"ZZ(J=1)\n\nRepresents a sigma_z^(i) sigma_z^(j) term. Multiply with a matrix J to get the term sum_ineq j J_ij sigma_z^(i) sigma_z^(j) \n\nnote: Note\nTwo-body terms do 'double-counting' of the couplings meaning the sum runs over i,j independently and not i < j.\n\n\n\n\n\n","category":"function"},{"location":"full-docs/hamiltonian/#SpinModels.Hopp","page":"Full docs for Hamiltonian","title":"SpinModels.Hopp","text":"Hopp(J=1)\n\nRepresents a sigma_+^(i) sigma_-^(j) + mathrmhc term. Multiply with a matrix J to get the term sum_ineq j J_ij sigma_+^(i) sigma_-^(j) + mathrmhc \n\nnote: Note\nTwo-body terms do 'double-counting' of the couplings meaning the sum runs over i,j independently and not i < j.\n\n\n\n\n\n","category":"function"},{"location":"full-docs/hamiltonian/#SpinModels.FlipFlop","page":"Full docs for Hamiltonian","title":"SpinModels.FlipFlop","text":"FlipFlop(J=1)\n\nSame as Hopp.\n\n\n\n\n\n","category":"function"},{"location":"full-docs/hamiltonian/#SpinModels.XXZ","page":"Full docs for Hamiltonian","title":"SpinModels.XXZ","text":"XXZ(Δ, J=1)\n\nRepresents the terms: Hopp(2*J)+Δ*ZZ(J)\n\n\n\n\n\n","category":"function"},{"location":"full-docs/hamiltonian/#Functions","page":"Full docs for Hamiltonian","title":"Functions","text":"","category":"section"},{"location":"full-docs/hamiltonian/","page":"Full docs for Hamiltonian","title":"Full docs for Hamiltonian","text":"nspins(::SpinModels.Hamiltonian)","category":"page"},{"location":"full-docs/hamiltonian/#SpinModels.nspins-Tuple{SpinModels.Hamiltonian}","page":"Full docs for Hamiltonian","title":"SpinModels.nspins","text":"nspins(term)\nnspins(sum_of_terms)\n\nReturn the number of spins in the underlying geometry if that information was already provided.\n\n\n\n\n\n","category":"method"},{"location":"full-docs/internal/#Non-exported-symbols","page":"Non exported symbols","title":"Non exported symbols","text":"","category":"section"},{"location":"full-docs/internal/","page":"Non exported symbols","title":"Non exported symbols","text":"SpinModels.NNGeometry\nSpinModels.NNInteraction\nSpinModels._nbodyness\nSpinModels._hopping_coordinates!\nSpinModels._z_field!\nSpinModels._zcorrelator_values!\nSpinModels._check_add_compatible\nSpinModels._check_array","category":"page"},{"location":"full-docs/internal/#SpinModels.NNGeometry","page":"Non exported symbols","title":"SpinModels.NNGeometry","text":"NNGeometry(geometry, k)\n\nSee NN(::SpinModels.Geometry) and nearest_neighbor_from_distances!.\n\n\n\n\n\n","category":"type"},{"location":"full-docs/internal/#SpinModels.NNInteraction","page":"Non exported symbols","title":"SpinModels.NNInteraction","text":"NNInteraction(interaction, k)\n\nSee NN(::SpinModels.Interaction) and nearest_neighbor_from_interactions!.\n\n\n\n\n\n","category":"type"},{"location":"full-docs/internal/#SpinModels._nbodyness","page":"Non exported symbols","title":"SpinModels._nbodyness","text":"_nbodyness(term)\n\nNumber of spins the term acts on.\n\n\n\n\n\n","category":"function"},{"location":"full-docs/internal/#SpinModels._hopping_coordinates!","page":"Non exported symbols","title":"SpinModels._hopping_coordinates!","text":"_hopping_coordinates!(A, B, N, i, j)\n\nCompute the coordinates of the matrix entries of the hopping operator between site i and J in a spin chain. The values are 1 everywhere, so there is no input/output for these in this function.\n\nHopping operator: σ₊ⁱσ₋ʲ + σ₋ⁱσ₊ʲ = 𝟙⊗…⊗𝟙⊗σ₊⊗𝟙⊗…⊗𝟙⊗σ₋⊗𝟙⊗…⊗𝟙 + transpose\n\nParameters\n\nA,B coordinate vectors. These should have length 2^(N-2) will be overwritten!\nN number of sites in the system\ni,j sites the excitation can hop between. 1 ≤ i,j ≤ N\n\n\n\n\n\n","category":"function"},{"location":"full-docs/internal/#SpinModels._z_field!","page":"Non exported symbols","title":"SpinModels._z_field!","text":"_z_field!(V, N, hz, i)\n\nADD to V the diagonal of hz σzⁱ.\n\nParameters\n\nV array of length 2^N. Will be mutated!\nN number of sites in the system\nhz field strength\ni site to act on\n\n\n\n\n\n","category":"function"},{"location":"full-docs/internal/#SpinModels._zcorrelator_values!","page":"Non exported symbols","title":"SpinModels._zcorrelator_values!","text":"    _zcorrelator_values!(V, N, Jij, i, j)\n\nADD to V the diagonal of Jᵢⱼ σzⁱσzʲ.\n\nParameters\n\nV array of length 2^N. Will be mutated!\nN number of sites in the system\nJij interaction strength\ni,j sites to couple\n\n\n\n\n\n","category":"function"},{"location":"full-docs/internal/#SpinModels._check_add_compatible","page":"Non exported symbols","title":"SpinModels._check_add_compatible","text":"checkadd_compatible(term1, term2)\n\nTwo terms can be added, if either:\n\nboth have no spatial information\nboth carry spatial information for the same number of spins\n\nOtherwise error.\n\n\n\n\n\n","category":"function"},{"location":"full-docs/internal/#SpinModels._check_array","page":"Non exported symbols","title":"SpinModels._check_array","text":"_check_array(dims, nbody)\n\nBad name. Raise an error if an array with size dims is not suitable as coefficient for an nbody term. Generally a n-body term needs a n-dimensional array which axes are all equal.\n\n\n\n\n\n","category":"function"},{"location":"couplings/#Coupling-matrices","page":"Coupling matrices","title":"Coupling matrices","text":"","category":"section"},{"location":"couplings/","page":"Coupling matrices","title":"Coupling matrices","text":"The task of creating a coupling matrix is divided into two parts:","category":"page"},{"location":"couplings/","page":"Coupling matrices","title":"Coupling matrices","text":"The Geometry describes how many spins in what spatial configuration there are\nThe Interaction knows how to translate the positions into coupling strengths","category":"page"},{"location":"couplings/","page":"Coupling matrices","title":"Coupling matrices","text":"So generally you want to first build up the desired geometry and then apply the correct interaction to that. Of course you can always directly provide a coupling matrix if these standards do not cover your needs. Just keep in mind the index order.","category":"page"},{"location":"couplings/#Geometry","page":"Coupling matrices","title":"Geometry","text":"","category":"section"},{"location":"couplings/","page":"Coupling matrices","title":"Coupling matrices","text":"There are several simple base geometries: Chain, NoisyChain and Box. The latter two don't admit to fixed positions and instead draw them randomly each time positions is called on them.","category":"page"},{"location":"couplings/","page":"Coupling matrices","title":"Coupling matrices","text":"Then there are several modifiers, you can apply:","category":"page"},{"location":"couplings/","page":"Coupling matrices","title":"Coupling matrices","text":"PBC to enforce periodic boundary conditions. This basically changes how distances are computed within the geometry.\nBlockaded to enforce a minimal distance between spins. This only makes sense to apply to disordered geometries.\nNN to only keep nearest neighbour distances.","category":"page"},{"location":"couplings/","page":"Coupling matrices","title":"Coupling matrices","text":"You may use function like positions, distance_matrix and nearest_neighbor_from_distances on geometry objects.","category":"page"},{"location":"couplings/#Interaction","page":"Coupling matrices","title":"Interaction","text":"","category":"section"},{"location":"couplings/","page":"Coupling matrices","title":"Coupling matrices","text":"The next step is to define the desired interaction type. Currently implemented are ConstantInteraction and PowerLaw. You can also modify these to be nearest-neighbor only with NN. There are quality of life overloads, that removes the need to use ConstantInteraction(1)(geometry) such that you can simply use the geometry like an interaction matrix in most scenarios. Please complain if this produces unintuitive behavior.","category":"page"},{"location":"couplings/","page":"Coupling matrices","title":"Coupling matrices","text":"Once you have defined your desired geometry and interaction you can get the final interaction matrix via interaction_matrix (or equivalently by calling the interaction object with the geometry).","category":"page"},{"location":"#SpinModels.jl","page":"Home","title":"SpinModels.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Construct all of your favorite spin-frac12 models with ease.","category":"page"},{"location":"#Recipe","page":"Home","title":"Recipe","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Define a coupling matrix like J = PowerLaw(6)(Chain(10))\nAdd the relevant terms of your Hamiltonian and apply couplings like H = J*(XX() + YY() + Δ*ZZ())\nConvert to sparse/dense with SparseArrays.sparse(H)/LinearAlgebra.Matrix(H)","category":"page"},{"location":"#Usage-examples","page":"Home","title":"Usage examples","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The (standard) MBL Hamiltonian hatH = fracJ4sum_i left(sigma_+^(i)sigma_-^(j) + sigma_-^(i)sigma_+^(j) + 2Delta sigma_z^(i)sigma_z^(j)right) + sum_i h_i sigma_z^(i) could be constructed like so:","category":"page"},{"location":"","page":"Home","title":"Home","text":"N = 10 # number of spins\nh = 5.0 # strength of disorder\nΔ = 1.0 # strength of Ising coupling, here XXX model\ngeometry = NN(Chain(N)) # nearest neighbor chain...\nJ_matrix = interaction_matrix(ConstantInteraction(J), geometry) # ... with constant interactions\n# could also use ConstantInteraction(J)(geometry)\nh_vec = h*(rand(N) .- 1) # random fields\nH = J_matrix/4 * (Hopping() + 2*Δ*ZZ()) + h_vec*Z()","category":"page"},{"location":"","page":"Home","title":"Home","text":"This object H can now be converted to either sparse or dense matrices using SparseArrays.sparse(H) or LinearAlgebra.Matrix(H) respectively.","category":"page"},{"location":"#Structure-of-this-package","page":"Home","title":"Structure of this package","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A spin model's Hamiltonian consists of the sum of different terms, such as XX or Z, that denote different Pauli operator combinations (here: sigma^(i)_x sigma^(j)_x and sigma_z^(i)). To be able to realize the Hamiltonian's matrix, each term needs to have information about the coupling strengths. For one-body terms this is a vector and for two-body terms a matrix. The sizes of these coupling arrays correspond to the amount of spins in the model (and as such need to be equal across all terms).","category":"page"},{"location":"","page":"Home","title":"Home","text":"For ease of use, there are several helpful structs and functions implemented to make constructing commonly used coupling types easy. They are organized into Geometrys and Interactions. See Coupling matrices.","category":"page"},{"location":"#index_order","page":"Home","title":"A word on ordering the basis states","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The basis of a spin-frac12 model can be easily represented by numbers in binary, where each 0 denotes a uparrowrangle and each 1 denotes a downarrowrangle. The 3-spin state uparrowdownarrowuparrowrangle for example gets the number 010_2 = 2_10 and as such is the third basis state (1-based counting).","category":"page"},{"location":"","page":"Home","title":"Home","text":"There are now two ways of organizing the basis with respect to the coupling matrices: The coupling J_ij couples the ith spin with the jth spin, where the ith spin is","category":"page"},{"location":"","page":"Home","title":"Home","text":"the spin at position i in the state vector (order like 1rangleotimes2rangleotimes3rangleldots)\nthe spin corresponding to the ith digit in the binary expansion with value 2^i-1 (order like ldots3rangleotimes2rangleotimes1rangle)","category":"page"},{"location":"","page":"Home","title":"Home","text":"I'd argue that for the sake of doing bit manipulation on the state numbers the latter is more convenient as one does not need to know the total amount of spins. Currently this library counts its spins from the LEFT corresponding to the FORMER (1.) option (due to legacy code reasons).","category":"page"},{"location":"","page":"Home","title":"Home","text":"Ultimately this issue boils down to the inconsistency in writing generally from left to right except when writing numbers, where the lowest valued place is on the right and increasing to the left.","category":"page"},{"location":"#Full-documentation","page":"Home","title":"Full documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"full-docs/hamiltonian.md\",\"full-docs/geometry.md\",\"full-docs/interaction.md\"]\nDepth = 1","category":"page"}]
}
