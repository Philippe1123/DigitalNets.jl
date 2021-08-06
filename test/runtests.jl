using DigitalNets, SpecialFunctions, Statistics, Test, PyPlot

@testset "Digital Nets" begin
    QMCPointsNoInterlacing=[0
    0.5
    0.75
    0.25
    0.375
    0.875
    0.625
    0.125
    0.1875
    0.6875
    0.9375
    0.4375
    0.3125
    0.8125
    0.5625
    0.0625
    0.09375
    0.59375
    0.84375
    0.34375
    0.46875
    0.96875
    0.71875
    0.21875
    0.15625
    0.65625
    0.90625
    0.40625
    0.28125
    0.78125
    0.53125
    0.03125
    0.046875
    0.546875
    0.796875
    0.296875
    0.421875
    0.921875
    0.671875
    0.171875
    0.234375
    0.734375
    0.984375
    0.484375
    0.359375
    0.859375
    0.609375
    0.109375
    0.078125
    0.578125
    0.828125
    0.328125
    0.453125
    0.953125
    0.703125
    0.203125
    0.140625
    0.640625
    0.890625
    0.390625
    0.265625
    0.765625
    0.515625
    0.015625
    0.0234375
    0.5234375
    0.7734375
    0.2734375
    0.3984375
    0.8984375
    0.6484375
    0.1484375
    0.2109375
    0.7109375
    0.9609375
    0.4609375
    0.3359375
    0.8359375
    0.5859375
    0.0859375
    0.1171875
    0.6171875
    0.8671875
    0.3671875
    0.4921875
    0.9921875
    0.7421875
    0.2421875
    0.1796875
    0.6796875
    0.9296875
    0.4296875
    0.3046875
    0.8046875
    0.5546875
    0.0546875
    0.0390625
    0.5390625
    0.7890625
    0.2890625
    0.4140625
    0.9140625
    0.6640625
    0.1640625
    0.2265625
    0.7265625
    0.9765625
    0.4765625
    0.3515625
    0.8515625
    0.6015625
    0.1015625
    0.0703125
    0.5703125
    0.8203125
    0.3203125
    0.4453125
    0.9453125
    0.6953125
    0.1953125
    0.1328125
    0.6328125
    0.8828125
    0.3828125
    0.2578125
    0.7578125
    0.5078125
    0.0078125
    0.01171875
    0.51171875
    0.76171875
    0.26171875
    0.38671875
    0.88671875
    0.63671875
    0.13671875
    0.19921875
    0.69921875
    0.94921875
    0.44921875
    0.32421875
    0.82421875
    0.57421875
    0.07421875
    0.10546875
    0.60546875
    0.85546875
    0.35546875
    0.48046875
    0.98046875
    0.73046875
    0.23046875
    0.16796875
    0.66796875
    0.91796875
    0.41796875
    0.29296875
    0.79296875
    0.54296875
    0.04296875
    0.05859375
    0.55859375
    0.80859375
    0.30859375
    0.43359375
    0.93359375
    0.68359375
    0.18359375
    0.24609375
    0.74609375
    0.99609375
    0.49609375
    0.37109375
    0.87109375
    0.62109375
    0.12109375
    0.08984375
    0.58984375
    0.83984375
    0.33984375
    0.46484375
    0.96484375
    0.71484375
    0.21484375
    0.15234375
    0.65234375
    0.90234375
    0.40234375
    0.27734375
    0.77734375
    0.52734375
    0.02734375
    0.01953125
    0.51953125
    0.76953125
    0.26953125
    0.39453125
    0.89453125
    0.64453125
    0.14453125
    0.20703125
    0.70703125
    0.95703125
    0.45703125
    0.33203125
    0.83203125
    0.58203125
    0.08203125
    0.11328125
    0.61328125
    0.86328125
    0.36328125
    0.48828125
    0.98828125
    0.73828125
    0.23828125
    0.17578125
    0.67578125
    0.92578125
    0.42578125
    0.30078125
    0.80078125
    0.55078125
    0.05078125
    0.03515625
    0.53515625
    0.78515625
    0.28515625
    0.41015625
    0.91015625
    0.66015625
    0.16015625
    0.22265625
    0.72265625
    0.97265625
    0.47265625
    0.34765625
    0.84765625
    0.59765625
    0.09765625
    0.06640625
    0.56640625
    0.81640625
    0.31640625
    0.44140625
    0.94140625
    0.69140625
    0.19140625
    0.12890625
    0.62890625
    0.87890625
    0.37890625
    0.25390625
    0.75390625
    0.50390625
    0.00390625 ]


    QMCPointsInterlacingFactorTwo=[0
    0.75
    0.6875
    0.4375
    0.234375
    0.984375
    0.546875
    0.296875
    0.10546875
    0.85546875
    0.66796875
    0.41796875
    0.15234375
    0.90234375
    0.58984375
    0.33984375
    0.0927734375
    0.8427734375
    0.6552734375
    0.4052734375
    0.1708984375
    0.9208984375
    0.6083984375
    0.3583984375
    0.0498046875
    0.7998046875
    0.7373046875
    0.4873046875
    0.1904296875
    0.9404296875
    0.5029296875
    0.2529296875
    0.065185546875
    0.815185546875
    0.627685546875
    0.377685546875
    0.174560546875
    0.924560546875
    0.612060546875
    0.362060546875
    0.045654296875
    0.795654296875
    0.733154296875
    0.483154296875
    0.217529296875
    0.967529296875
    0.530029296875
    0.280029296875
    0.029052734375
    0.779052734375
    0.716552734375
    0.466552734375
    0.232177734375
    0.982177734375
    0.544677734375
    0.294677734375
    0.111083984375
    0.861083984375
    0.673583984375
    0.423583984375
    0.126708984375
    0.876708984375
    0.564208984375
    0.314208984375 ]


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


    @testset "ComputePoints" begin
        digital_net=DigitalNet64(10)
        out=digital_net[0:100]
        point1=out[55]
        point2=digital_net[55-1]
        @test point1==point2

    end

    #Test if QMC points are correctly generated QMCPoints obtained from https://people.cs.kuleuven.be/~dirk.nuyens/qmc-generators/
    @testset "CheckCorrectQMCPointsNoInterlacing" begin


        len=length(QMCPointsNoInterlacing)
        digital_net64=DigitalNet64(1)
        digital_net32=DigitalNet32(1)
        QMCPointsFromdigital_net64=digital_net64[0:len-1]
        QMCPointsFromdigital_net32=digital_net64[0:len-1]

        for id=1:len
            @test  QMCPointsNoInterlacing[id] == QMCPointsFromdigital_net64[id][1]
            @test  QMCPointsNoInterlacing[id] == QMCPointsFromdigital_net32[id][1]
        end
    end



    @testset "CheckCorrectQMCPointsInterlacingFactorTwo" begin



        len=length(QMCPointsInterlacingFactorTwo)
        digital_net64_interlaced_two=DigitalNet64InterlacedTwo(1)
        QMCPointsFromdigital_net64=digital_net64_interlaced_two[0:len-1]

        for id=1:len
            @test  QMCPointsInterlacingFactorTwo[id] == QMCPointsFromdigital_net64[id][1]
        end

    end


    @testset "CheckRestartFromBegining" begin

        digital_net64=DigitalNet64(1)
        QMCPointsFromdigital_net64_set1=digital_net64[0:10]
        QMCPointsFromdigital_net64_set2=digital_net64[0:10]
        @test QMCPointsFromdigital_net64_set1==QMCPointsFromdigital_net64_set2

    end


    @testset "CheckPointGenerationFromOtherPointThanZero" begin

        digital_net64=DigitalNet64(1)
        StartPoint = 10
        EndPoint = 20
        QMCPointsFromdigital_net64_set1=digital_net64[StartPoint:EndPoint]
        for id=StartPoint+1:EndPoint+1
            @test digital_net64[id-1][1] == QMCPointsNoInterlacing[id]
        end
    end


end
