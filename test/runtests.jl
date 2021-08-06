using DigitalNets, SpecialFunctions, Statistics, Test

@testset "Digital Nets" begin

    @testset "Constructor Digital Nets" begin



    @testset "DigitalNet32(s)" begin
        digital_net=DigitalNet32(10)
        @test ndims(digital_net) == 10
    end

    @testset "DigitalNet64(s)" begin
        digital_net=DigitalNet64(10)
        @test ndims(digital_net) == 10
    end

    @testset "DigitalNet64InterlacedTwo(s)" begin
        digital_net=DigitalNet64InterlacedTwo(10)
        @test ndims(digital_net) == 10
    end

    @testset "DigitalNet64InterlacedThree(s)" begin
        digital_net=DigitalNet64InterlacedThree(10)
        @test ndims(digital_net) == 10
    end

    @testset "DigitalNet64InterlacedThree64(s)" begin
        digital_net=DigitalNet64InterlacedThree64(10)
        @test ndims(digital_net) == 10
    end


    @testset "DigitalNet64InterlacedFour(s)" begin
        digital_net=DigitalNet64InterlacedFour(10)
        @test ndims(digital_net) == 10
    end

end
end
