@testset "linalg" begin
    @testset "depth_vector" begin
        @test depth_vector([0, 0, 1, 0, 2]) == (3, 1)
        @test depth_vector([1.0, 0.0, 3.0, 4.1]) == (1, 1.0) 
        @test depth_vector(zeros(Int, 3)) == (0, 0)

        F = GF(5)
        @test depth_vector(zeros(F, 3)) == (0, zero(F))
        @test depth_vector(F.([5, 10, 7, 3])) == (3, F(2))
    end
    @testset "vector_in_subspace"begin
        W = SubSpace{Int}([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1,0]], [1, 2, 3])
        @test !(vector_in_subspace(W, [0, 0, 0, 1]))
        @test vector_in_subspace(W, [3, 1, 2, 0]) == (true, [3, 1, 2])
        @test vector_in_subspace(W, [1, 2, 3, 4]; add = true) == (false, [1, 2, 3, 4])

        @test_throws InexactError SubSpace([[1, 2, 3, 0], [1, 3, 2, 0], [4, 4, 1, 1], [2, 5, 5, 0]]) 
        F = GF(7)
        V = SubSpace([F.([1, 2, 3, 0]), F.([1, 3, 2, 0]), F.([4, 4, 1, 1]), F.([2, 5, 5, 0])])
        @test V.positions == [1, 2, 3]
        @test V.basis[1] == F.([1, 2, 3, 0])
        @test V.basis[2] == F.([0, 1, 6, 0])
        @test V.basis[3] == F.([0, 0, 1, 6])

    end

    function is_identity(I, n)
        (d, l) = size(I)
        @test d == l == n
        @test all([isone(I[i, i]) for i = 1:n])
        @test all([iszero(I[i,j]) for i = 1:n for j = 1:n if i!=j])
    end        
    @testset "invert_matrix" begin
        F = GF(7)
        A = F.([1 2 3; 4 5 6; 7 8 1])
        Ainv = invert_matrix(A)
        is_identity(A*Ainv, 3)
        is_identity(Ainv*A, 3)

        F = GF(31)
        A = F.([1 2 4 3; 5 6 7 8; 11 12 9 10; 15 16 13 14])
        is_identity(A*Ainv, 4)
        is_identity(Anv*A, 4)

    end 
    
    function is_nullspace(A, W, dim)
        @test length(W.positions) == dim
        for i = 1:dim
            @test all(iszero.(transpose(A)*W.basis[i]))
        end
    end
    @testset "nullspace" begin
        F = GF(31)
        A = F.([1 2 3 4; 5 6 7 8; 11 12 9 10; 15 16 13 14])
        W = nullspace(A)
        is_nullspace(A, W, 1)

        A = F.([1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16])
        W = nullspace(A)
        is_nullspace(A, W, 2)

        A = F.([1 2 3 4; 13 14 15 16; 5 6 7 8; 9 10 11 12])
        W = nullspace(A)
        is_nullspace(A, W, 2)


        A = F.([1 0 0; 2 0 0; 0 1 0; 1 2 0; 0 0 1; 1 2 2])
        W = nullspace(A)
        is_nullspace(A, W, 3)

    end


end

