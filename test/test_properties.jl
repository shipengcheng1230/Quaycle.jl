using Test

@testset "properties" begin

    @testset "frictional properties" begin
        tmpfile = tempname()
        fp = FrictionalProperties([rand(5, 3) for _ in 1: 4]..., RForm(), DieterichStateLaw())
        save_properties(tmpfile, fp; option="w")
        fp_recover = read_properties(tmpfile)[description(fp)]
        for field in fieldnames(fp)
            @test getfield(fp, Symbol(field)) == getfield(fp_recover, Symbol(field))
        end
    end

    @testset "fault properties" begin
        tmpfile = tempname()
        sp = HomoFaultProperties(rand(6)...)
        save_properties(tmpfile, sp; option="w")
        sp_recover = read_properties(tmpfile)[description(sp)]
        for field in fieldnames(sp)
            @test getfield(sp, Symbol(field)) == getfield(sp_recover, Symbol(field))
        end
    end

    @testset "multi-properites" begin
        tmpfile = tempname()
        sp = HomoFaultProperties(rand(6)...)
        fp = FrictionalProperties([rand(5, 3) for _ in 1: 4]..., RForm(), DieterichStateLaw())
        save_properties(tmpfile, [sp, fp])
        p_recover = read_properties(tmpfile)
        sp_recover = p_recover[description(sp)]
        fp_recover = p_recover[description(fp)]
        for field in fieldnames(sp)
            @test getfield(sp, Symbol(field)) == getfield(sp_recover, Symbol(field))
        end
        for field in fieldnames(fp)
            @test getfield(fp, Symbol(field)) == getfield(fp_recover, Symbol(field))
        end
    end

end
