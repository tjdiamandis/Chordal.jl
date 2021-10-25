using Random

@testset "Max det completion" begin
    # -- Test 1: matrix with known completion --
    # From VA Figure 4.2
    n = 17
    ijs = [(3,1), (3,2), (4,1), (4,2), (4,3), (5,1), (5,3), (5,4), (8,7),
           (9,5), (9,6), (9,7), (9,8), (11,10), (13,10), (13,11), (13,12),
           (14,10), (14,11), (14,12), (14,13),
           (15, 1), (15,3), (15,4), (15,5), (15,7), (15,8), (15,9),
           (16,5), (16,6), (16,9), (16,12), (16,13), (16,14), (16,15),
           (17,10), (17,11), (17,12), (17,13), (17,14), (17,15), (17,16)]
    append!(ijs, [(i,i) for i in 1:n])
    II, JJ = Chordal.unzip(ijs)
    A = sparse(II, JJ, ones(length(II)))
    A = (A + A') / 2
    sp = sparsity_pattern(A)

    # True solution
    #=
    X = Semidefinite(n)
    p = maximize(logdet(X), [X[II[x],JJ[x]] == A[II[x],JJ[x]] for x in 1:length(II)])
    Convex.solve!(p, () -> Mosek.Optimizer())
    Xstar = X.value
    =#
    Xstar = [1.0 0.3333370598035662 0.5 0.5 0.5 0.22222287873041388 0.27778132246744036 0.2777811910239065 0.33333530246961796 0.1712975423472041 0.17129709225708564 0.20370126720855525 0.2037019856666862 0.20370171663098346 0.5 0.33333088023046814 0.27777638852035724; 0.3333370598035662 1.0 0.5 0.5 0.33333305125898244 0.14814707689070564 0.1851867428728771 0.1851869926336214 0.22222250153752335 0.11419886412281514 0.11419870887817024 0.135800056328276 0.13580093153877904 0.13580062045920008 0.33333343444912716 0.22221930265929504 0.18518415224568624; 0.5 0.5 1.0 0.5 0.5 0.22222196163744942 0.27778076043784233 0.2777809794496988 0.33333475836820464 0.1712991354296111 0.1712988996359925 0.2037013235821753 0.2037028777637407 0.2037025691566877 0.5 0.33333017624316663 0.2777768088531805; 0.5 0.5 0.5 1.0 0.5 0.2222219730839787 0.27778081718375297 0.27778099169886683 0.3333347721005158 0.171299355067919 0.1712991188490877 0.20370132886164274 0.20370290806545718 0.20370259032235122 0.5 0.3333303500560772 0.27777693617839855; 0.5 0.33333305125898244 0.5 0.5 1.0 0.33333571187398225 0.33333659387305287 0.33333587588587565 0.5 0.2222299016007586 0.222228389459611 0.27777554232420304 0.2777775941596431 0.27777738847461536 0.5 0.5 0.3333356606726846; 0.22222287873041388 0.14814707689070564 0.22222196163744942 0.2222219730839787 0.33333571187398225 1.0 0.2777807210172411 0.27777952955889934 0.5 0.19907461348157693 0.1990723384116564 0.25925704288495866 0.25925593341309217 0.2592559179558187 0.3333352130700529 0.5 0.2777786944527559; 0.27778132246744036 0.1851867428728771 0.27778076043784233 0.27778081718375297 0.33333659387305287 0.2777807210172411 1.0 0.5 0.5 0.1712997337016747 0.17129908722574652 0.2037044400165765 0.2037058562144702 0.20370577454105926 0.5 0.3333317798520487 0.2777766860348415; 0.2777811910239065 0.1851869926336214 0.2777809794496988 0.27778099169886683 0.33333587588587565 0.27777952955889934 0.5 1.0 0.5 0.17129896845287737 0.1712983380154316 0.20370409144796173 0.2037052452025856 0.2037051577978674 0.5 0.3333308462064361 0.27777608804056464; 0.33333530246961796 0.22222250153752335 0.33333475836820464 0.3333347721005158 0.5 0.5 0.5 0.5 1.0 0.22222720796785625 0.22222501887881926 0.2777758212711955 0.2777770394826086 0.27777696102252825 0.5 0.5 0.3333344233472463; 0.1712975423472041 0.11419886412281514 0.1712991354296111 0.171299355067919 0.2222299016007586 0.19907461348157693 0.1712997337016747 0.17129896845287737 0.22222720796785625 1.0 0.5 0.3750025210974436 0.5 0.5 0.2916717850653076 0.37500249564668897 0.5; 0.17129709225708564 0.11419870887817024 0.1712988996359925 0.1712991188490877 0.222228389459611 0.1990723384116564 0.17129908722574652 0.1712983380154316 0.22222501887881926 0.5 1.0 0.3750032717610717 0.5 0.5 0.29167218644030213 0.3750017383599886 0.5; 0.20370126720855525 0.135800056328276 0.2037013235821753 0.20370132886164274 0.27777554232420304 0.25925704288495866 0.2037044400165765 0.20370409144796173 0.2777758212711955 0.3750025210974436 0.3750032717610717 1.0 0.5 0.5 0.33333681129144077 0.5 0.5; 0.2037019856666862 0.13580093153877904 0.2037028777637407 0.20370290806545718 0.2777775941596431 0.25925593341309217 0.2037058562144702 0.2037052452025856 0.2777770394826086 0.5 0.5 0.5 1.0 0.5 0.333338503723327 0.5 0.5; 0.20370171663098346 0.13580062045920008 0.2037025691566877 0.20370259032235122 0.27777738847461536 0.2592559179558187 0.20370577454105926 0.2037051577978674 0.27777696102252825 0.5 0.5 0.5 0.5 1.0 0.3333385764072932 0.5 0.5; 0.5 0.33333343444912716 0.5 0.5 0.5 0.3333352130700529 0.5 0.5 0.5 0.2916717850653076 0.29167218644030213 0.33333681129144077 0.333338503723327 0.3333385764072932 1.0 0.5 0.5; 0.33333088023046814 0.22221930265929504 0.33333017624316663 0.3333303500560772 0.5 0.5 0.3333317798520487 0.3333308462064361 0.5 0.37500249564668897 0.3750017383599886 0.5 0.5 0.5 0.5 1.0 0.5; 0.27777638852035724 0.18518415224568624 0.2777768088531805 0.27777693617839855 0.3333356606726846 0.2777786944527559 0.2777766860348415 0.27777608804056464 0.3333344233472463 0.5 0.5 0.5 0.5 0.5 0.5 0.5 1.0]
    pstar = det(Xstar)

    W = CD.maxdet_completion(A)
    @test ≈(A[CartesianIndex.(ijs)], W[CartesianIndex.(ijs)])
    @test det(W) ≈ pstar

    L, D = Chordal.maxdet_completion_etree(A)
    W = inv(D) * inv(LowerTriangular(Matrix(L)))
    W = L' \ W
    @test ≈(A[CartesianIndex.(ijs)], W[CartesianIndex.(ijs)])
    @test det(W) ≈ pstar

    L, D = Chordal.maxdet_completion_factors(A)
    W = inv(Matrix(D)) * inv(LowerTriangular(Matrix(L)))
    W = L' \ W
    @test ≈(A[CartesianIndex.(ijs)], W[CartesianIndex.(ijs)])
    @test det(W) ≈ pstar


    # -- Test 2: ensure pattern matches on random matrices --
    n = 50
    r = 5
    Random.seed!(0)
    M = rand(n,r)
    M = M*M'

    sp = spzeros(n, n)
    for i in 1:n
        block_size = randn() < 1.2 ? 2 : 10
        sp[i:min(i+block_size,n), i:min(i+block_size,n)] .= 1
    end

    M_ = sparse(sp .* M)
    M_comp = CD.maxdet_completion(M_)
    @test ≈(M_, sp.*M_comp, atol=1e-6)

    L, D = CD.maxdet_completion_etree(M_)
    W = inv(D) * inv(LowerTriangular(Matrix(L)))
    W = L' \ W
    @test ≈(M_, sp.*W, atol=1e-6)

    # L, D = Chordal.maxdet_completion_factors(M_)
    # W = inv(Matrix(D)) * inv(LowerTriangular(Matrix(L)))
    # W = L' \ W
    # @test ≈(M_, sp.*W, atol=1e-6)
end


@testset "Min rank completion" begin
   ## Build a test matrix
    # From Xin Jiang's masters thesis
    n = 9
    nonzero_ind_val = [
        (1,1,923.0),
        (1,4,464.0),
        (1,8,1215.0),
        (2,2,1721.0),
        (2,3,1494.0),
        (2,5,1182.0),
        (2,8,1638.0),
        (3,3,1300.0),
        (3,4,726.0),
        (3,5,1016.0),
        (3,8,1432.0),
        (4,4,418.0),
        (4,5,580.0),
        (4,8,770.0),
        (5,5,921.0),
        (5,8,1209.0),
        (5,9,959.0),
        (6,6,677.0),
        (6,7,421.0),
        (6,9,746.0),
        (7,7,270.0),
        (7,8,661.0),
        (7,9,499.0),
        (8,8,1773.0),
        (8,9,1195.0),
        (9,9,1021.0)
    ]
    sp = sparse(CD.unzip(nonzero_ind_val)...)
    sp = sp + triu(sp, 1)'
    s_pat = sparsity_pattern([sp])


    Y = CD.minrank_completion(sp)
    A_comp = Y*Y'
    matched = s_pat .* A_comp
    @test matched ≈ sp

    function bigger_mrc_test(n=50, r=4)
        function sparsity_structure(n)
            sp = spzeros(n, n)
            for i in 1:n
                block_size = randn() < 1 ? 2 : 10
                sp[i:min(i+block_size,n), i:min(i+block_size,n)] .= 1
            end
            return Int.(sp)
        end

        M = rand(n,r)
        M = M*M'
    
        sp = sparsity_structure(n)
    
        M_ = sparse(sp .* M)
        U = CD.minrank_completion(M_)
        M_comp = U*U'
        err = abs.(sp.*M_comp - M_)
        rmse = norm(err) / sqrt(length(M_))
        # @show sp
        # println(sp[end-4:end, end-4:end] == ones(Int, 5,5), rmse > 1e-10)
        return rmse
    end
    
    for N in [10, 20, 50, 200], r in [2, 5, 10, 20]
        r >= N && continue
        @test maximum([bigger_mrc_test(N, r) for i in 1:10]) < sqrt(eps())
    end
    
end
