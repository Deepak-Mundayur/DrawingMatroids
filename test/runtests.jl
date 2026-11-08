# Maximal matroid degenerations (src/maximat, a Julia port of Maximat, github.com/rprebet/maximal_matroids)
# checked against the results of Liwski, Mohammadi, Prébet, "Efficient algorithms for maximal
# matroid degenerations and irreducible decompositions of circuit varieties", J. Algebra 691
# (2026). The cases mirror upstream's tests/test_maximat.py.
using Test
using DrawingMatroids, Oscar

const MX = DrawingMatroids.Maximat

S(v...) = Set{Int}(v)
L(sets...) = Set{Int}[sets...]

const FANO = Any[L(), L(), L(S(1,4,5), S(1,2,3), S(1,6,7), S(3,4,7), S(2,5,7), S(2,4,6), S(3,5,6)), L(S(1:7...))]
const FANO_DUAL = Any[L(), L(), L(), L(S(4,5,6,7), S(2,3,5,6), S(2,3,4,7), S(1,3,5,7), S(1,3,4,6), S(1,2,4,5), S(1,2,6,7))]
const VAMOS = Any[L(), L(), L(), L(S(1,2,3,4), S(3,4,5,6), S(5,6,7,8), S(1,2,7,8), S(3,4,7,8))]
const STEINER = Any[L(), L(), L(), L(S(1,2,4,8), S(2,3,5,8), S(3,4,6,8), S(4,5,7,8), S(1,5,6,8), S(2,6,7,8), S(1,3,7,8),
                                    S(3,5,6,7), S(1,4,6,7), S(1,2,5,7), S(1,2,3,6), S(2,3,4,7), S(1,3,4,5), S(2,4,5,6))]
const GRID = L(S(1,2,3), S(4,5,6), S(7,8,9), S(1,4,7), S(2,5,8), S(3,6,9))
const M33 = Any[L(), L(), GRID, L(S(1,2,3,4,7), S(1,2,3,5,8), S(1,2,3,6,9), S(1,4,5,6,7), S(2,4,5,6,8),
                                  S(3,4,5,6,9), S(1,4,7,8,9), S(2,5,7,8,9), S(3,6,7,8,9))]
const EX56 = Any[L(), L(), L(), L(S(1,2,3,4), S(3,4,5,6), S(1,2,5,6))]
const AG23 = Any[L(), L(), L(S(1,2,3), S(4,5,6), S(7,8,9), S(1,4,7), S(2,5,8), S(3,6,9),
                             S(1,5,9), S(2,6,7), S(3,4,8), S(1,6,8), S(2,4,9), S(3,5,7)), L(S(1:9...))]

# (name, cyclic flats by rank 0..3, ground set size, number of maximal degenerations in the paper)
const EXAMPLES = [
    ("Fano (Example 4.7)", FANO, 7, 22),
    ("Fano dual (Section 6.3)", FANO_DUAL, 7, 22),
    ("Vámos (Section 6.1)", VAMOS, 8, 21),
    ("Vámos + free point", VAMOS, 9, 62),
    ("S(3,4,8) (Section 6.2)", STEINER, 8, 31),
    ("M(K33)* (Example 1.2)", M33, 9, 34),
    ("Example 5.6", EX56, 6, 10),
    ("AG(2,3) (Theorem 7.4)", AG23, 9, 31),
    ("U44 + 5||4", Any[L(), L(S(4,5)), L(), L()], 5, 3),
    ("Vámos + loop (Lemma 3.22)", Any[L(S(9)), L(), L(), VAMOS[4]], 9, 21),
    ("Vámos + 8||9 (Lemma 3.25)", Any[L(), L(S(8,9)), L(), VAMOS[4]], 9, 23),
]

@testset "maximal degenerations: paper examples" begin
    for (name, CF, d, expected) in EXAMPLES
        @testset "$name" begin
            result = MX.maximal_degenerations(CF, d)
            @test length(result) == expected
            # Every output N is a strict degeneration of M, and no two outputs are comparable.
            @test all(N -> MX.sup_cyclic(CF, N, d) && !MX.sup_cyclic(N, CF, d), result)
            @test all(MX.sup_cyclic(result[i], result[j], d) == false
                      for i in eachindex(result), j in eachindex(result) if i != j)
        end
    end
end

@testset "Fano degeneration types (Example 3.7)" begin
    result = MX.maximal_degenerations(FANO, 7)
    @test sort([minimum(only(N[1])) for N in result if !isempty(N[1])]) == 1:7   # M(i)
    @test count(N -> length(N[2]) == 3, result) == 7                             # B'_k
    @test count(N -> length(N[2]) == 1, result) == 7                             # A'_j
    @test any(N -> isempty(N[1]) && isempty(N[2]) && N[3] == [S(1:7...)] && isempty(N[4]), result)  # U(2,7)
end

canon(result) = sort([string([sort(sort.(collect.(level))) for level in N]) for N in result])

@testset "input formats, strata, parallelism (as upstream)" begin
    # Vámos with 8||9 and loop 10, given as true cyclic flats or in reduced form.
    full = Any[L(S(10)), L(S(8,9,10)), L(), L(S(1,2,3,4,10), S(3,4,5,6,10), S(5,6,7,8,9,10), S(1,2,7,8,9,10), S(3,4,7,8,9,10))]
    reduced = Any[L(S(10)), L(S(8,9)), L(), VAMOS[4]]
    @test canon(MX.maximal_degenerations(full, 10)) == canon(MX.maximal_degenerations(reduced, 10))
    # Computing all strata at once gives a subset of the union of the strata computed separately.
    parts = reduce(vcat, [MX.maximal_degenerations(VAMOS, 8; S=[i]) for i in 1:4])
    @test issubset(Set(canon(MX.maximal_degenerations(VAMOS, 8))), Set(canon(parts)))
    # Parallel and serial runs agree (start Julia with -t 2 or more to exercise threads).
    for (CF, d) in ((FANO, 7), (VAMOS, 9), (M33, 9))
        @test canon(MX.maximal_degenerations(CF, d; n_procs=2)) == canon(MX.maximal_degenerations(CF, d))
    end
end

@testset "former bugs fixed upstream" begin
    @test MX.rem_sub(Any[L(S(1,2,3), S(1,2,3)), L()]) == [L(S(1,2,3)), L()]
    @test MX.Hypergraph.diff_size(S(1,2,3,4), L(S(1,2,3,4,5)), 1)
end

@testset "input checking (check_cyclic_flats)" begin
    bad(CF, d) = @test_throws ArgumentError MX.maximal_degenerations(CF, d)
    bad(FANO, 6)                                                        # element outside [1..d]
    bad(FANO[1:3], 7)                                                   # not 4 levels
    bad(Any[L(S(1), S(2)), L(), L(), L()], 4)                           # two loop sets
    bad(Any[L(S(1)), L(S(1,2)), L(), L()], 4)                           # only part of the loops
    bad(Any[L(), L(S(1,2), S(2,3)), L(), L()], 4)                       # intersecting parallel classes
    bad(Any[L(), L(), L(S(1,2,3), S(1,2,4)), L()], 5)                   # lines sharing two points
    bad(Any[L(), L(), L(S(1,2,3), S(3,4,5)), L()], 5)                   # meeting lines not in a plane
    bad(Any[L(), L(), L(S(1,2,3)), L(S(1,2,3,4))], 5)                   # plane not cyclic
    @test_throws ArgumentError MX.maximal_degenerations(FANO, 7; S=[5])
    @test !hasmethod(Base.replace, Tuple{Vector{Any}, Vector{Any}})      # no methods added to Base
end

@testset "drawing adapter" begin
    # Oscar's cyclic flats give the same input as the paper's Fano data.
    @test length(MX.maximal_degenerations(DrawingMatroids.cyclic_flats_by_rank(Oscar.fano_matroid()), 7; S=1:3)) == 22
    # Every separating basis is a basis of M that becomes dependent in some maximal degeneration.
    M = Oscar.non_fano_matroid()
    separating_bases, degenerations = DrawingMatroids.maximal_degeneration_basis_transversal(M)
    bases = Set(sort(Int.(collect(B))) for B in Oscar.bases(M))
    @test all(B -> B in bases, separating_bases)
    @test all(N -> any(B -> DrawingMatroids.degeneration_rank(B, N) < 3, separating_bases), degenerations)
end

import RealizationSpaces as RSp
quiet(f) = redirect_stdout(f, devnull)
same_matroid(M1, M2) = Oscar.rank(M1) == Oscar.rank(M2) && Set(Oscar.nonbases(M1)) == Set(Oscar.nonbases(M2))

@testset "draw_component" begin
    # A rank-2 component (all points on a line) of a rank-3 space is drawn as a 2×n matrix.
    Z = quiet(() -> RSp.nonbasis_variety(RSp.non_fano_matroid()))
    C2 = only(filter(C -> Oscar.rank(RSp.matroid(C)) == 2, RSp.components(Z)))
    r = quiet(() -> draw_component(C2; filename=nothing, seed=1))
    @test size(r.matrix) == (2, 7) && RSp.satisfies_matroid(r.matrix, RSp.matroid(C2); tol=1e-7)
    # A wrong stored matroid is replaced by the one read off the witness points.
    M5 = matroid_from_nonbases([[1, 2, 3]], 5)
    C = only(RSp.components(quiet(() -> RSp.realization_space(M5))))
    RSp._set_component_matroid!(C, matroid_from_nonbases([[1, 2, 3], [2, 3, 4], [2, 3, 5]], 5))
    r = quiet(() -> draw_component(C; filename=nothing, seed=1))
    @test same_matroid(RSp.matroid(C), M5) && RSp.satisfies_matroid(r.matrix, M5; tol=1e-7)
    # Without a rerun that agrees with the witness points, the component is not drawn.
    RSp._set_component_matroid!(C, matroid_from_nonbases([[1, 2, 3], [2, 3, 4], [2, 3, 5]], 5))
    @test_throws ErrorException DrawingMatroids._component_matroid!(C; tol=1e-7, attempts=0)
end
