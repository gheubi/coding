using Test, LinearAlgebra

include("../src/main.jl")

@testset "Polar Transform Test" begin 

    println("Bit Reveral Operator Test")
    b1 = [1, 2]
    for b in b1
        @test b == bit_rev_op(b, 1)
    end
    b2 = [1, 2, 3, 4]
    b2ex = [1, 3, 2, 4]
    for i in eachindex(b2)
        @test b2ex[i] == bit_rev_op(b2[i], 2)
    end
    b3ex = [1,5,3,7,2,6,4,8]
    for i in 1:8
        @test b3ex[i] == bit_rev_op(i, 3)
    end

    println("Polar Transform Test")
    println("   N=2 Polar Transform Test")
    F = [1 0; 1 1]
    F_n = [1 0; 1 1]
    for i in [0, 1], j in [0, 1]
        @test vec([i j] * F) .% 2 == polar_transform!(UInt8.([i, j]))
    end
    for n in 2:13
        N = 1 << n
        println("   Test N=$N polar transform")
        v = rand([0,1], N)
        F_n = kron(F_n, F)
        
        vB_n = map(i -> v[bit_rev_op(i, n)], 1:N)
        @test vec(transpose(vB_n)*F_n) .% 2 == polar_transform!(UInt8.(v))
    end

    println("Polar Transform benchmark")
    F = [1 0; 1 1]
    F_n = [1 0; 1 1]
    for n in 2:10
        N = 1 << n
        println("   Test N=$N polar transform")
        v = rand([0,1], N)
        F_n = kron(F_n, F)
        print("       Time to compute matrix product")
        @time vec(transpose(map(i -> v[bit_rev_op(i, n)], 1:N))*F_n) .% 2
        print("       Time to compute polar_transform!")
        @time polar_transform!(UInt8.(v))
    end
end