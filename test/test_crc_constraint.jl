include("../src/main.jl")


using Test

# Convert a byte to MSB-first bits
function byte_to_bits(x::UInt8)
    return BIT[(x >> shift) & 1 for shift in 7:-1:0]
end

@testset "CRC computation tests" begin

    @testset "CRC-8 / ATM known test vector" begin
        # CRC-8 polynomial: x^8 + x^2 + x + 1
        width = NAT(8)
        poly = BIT[1, 0, 0, 0, 0, 0, 1, 1, 1]

        message = reduce(
            vcat,
            [byte_to_bits(b) for b in codeunits("123456789")]
        )

        work = vcat(message, zeros(BIT, width))
        dst = zeros(BIT, width)

        result = crc_compute!(width, poly, work, dst)

        # Expected CRC-8/ATM = 0xF4
        expected = byte_to_bits(0xF4)

        @test result == expected
        @test dst == expected
    end

    @testset "Zero message" begin
        width = NAT(8)
        poly = BIT[1, 0, 0, 0, 0, 0, 1, 1, 1]

        message = zeros(BIT, 16)
        work = vcat(message, zeros(BIT, width))
        dst = zeros(BIT, width)

        result = crc_compute!(width, poly, work, dst)

        @test result == zeros(BIT, width)
    end

    @testset "Single bit message" begin
        width = NAT(8)
        poly = BIT[1, 0, 0, 0, 0, 0, 1, 1, 1]

        message = BIT[1]
        work = vcat(message, zeros(BIT, width))
        dst = zeros(BIT, width)

        result = crc_compute!(width, poly, work, dst)

        @test length(result) == width
        @test result == dst
    end

    @testset "Message beginning with zero" begin
        width = NAT(8)
        poly = BIT[1, 0, 0, 0, 0, 0, 1, 1, 1]

        message = BIT[0, 0, 1, 1, 0, 1, 0, 1]
        work = vcat(message, zeros(BIT, width))
        dst = zeros(BIT, width)

        result = crc_compute!(width, poly, work, dst)

        @test length(result) == width
        @test result == dst
    end

    @testset "Different message lengths" begin
        width = NAT(8)
        poly = BIT[1, 0, 0, 0, 0, 0, 1, 1, 1]

        for len in (8, 16, 24, 32, 64, 128)
            message = BIT[isodd(i) for i in 1:len]
            work = vcat(message, zeros(BIT, width))
            dst = zeros(BIT, width)

            result = crc_compute!(width, poly, work, dst)

            @test length(result) == width
            @test result == dst
        end
    end
end
