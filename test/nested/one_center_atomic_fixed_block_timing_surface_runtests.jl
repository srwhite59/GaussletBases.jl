# Integration/slow test. Do not include in default nested runner.

@testset "One-center atomic fixed-block timing surface" begin
    function _timing_labels(report::GaussletBases.TimeG.TimingReport)
        labels = String[]
        function _visit(node)
            push!(labels, node.label)
            foreach(_visit, node.children)
        end
        foreach(_visit, report.roots)
        return labels
    end

    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)

    timed_full = one_center_atomic_full_parent_fixed_block(
        basis;
        exponents = expansion.exponents,
        nside = 5,
        timing = true,
    )
    timed_legacy = one_center_atomic_legacy_profile_fixed_block(
        basis;
        exponents = expansion.exponents,
        working_box = 2:12,
        nside = 5,
        timing = true,
    )

    @test timed_full isa GaussletBases.TimedNestedFixedBlockBuild
    @test timed_legacy isa GaussletBases.TimedNestedFixedBlockBuild
    @test timed_full.timings isa GaussletBases.TimeG.TimingReport
    @test timed_legacy.timings isa GaussletBases.TimeG.TimingReport
    @test timed_full.fixed_block.shell.working_box == (1:13, 1:13, 1:13)
    @test timed_legacy.fixed_block.shell.working_box == (2:12, 2:12, 2:12)
    @test !hasproperty(timed_full.fixed_block, :gaussian_terms)
    @test !hasproperty(timed_full.fixed_block, :pair_terms)
    @test !hasproperty(timed_full.fixed_block, :term_storage)
    @test !isnothing(timed_full.fixed_block.gaussian_sum)
    @test !isnothing(timed_full.fixed_block.pair_sum)
    @test !hasproperty(timed_legacy.fixed_block, :gaussian_terms)
    @test !hasproperty(timed_legacy.fixed_block, :pair_terms)
    @test !hasproperty(timed_legacy.fixed_block, :term_storage)
    @test !isnothing(timed_legacy.fixed_block.gaussian_sum)
    @test !isnothing(timed_legacy.fixed_block.pair_sum)
    @test norm(timed_full.fixed_block.overlap - I, Inf) < 1.0e-10
    @test norm(timed_legacy.fixed_block.overlap - I, Inf) < 1.0e-10
    full_labels = _timing_labels(timed_full.timings)
    legacy_labels = _timing_labels(timed_legacy.timings)
    @test "fixed_block.total" in full_labels
    @test "fixed_block.parent_bundle" in full_labels
    @test "fixed_block.sequence_build" in full_labels
    @test "fixed_block.adapter" in full_labels
    @test "diatomic.packet.total" in full_labels
    @test "diatomic.packet.gaussian_terms" in full_labels
    @test "diatomic.packet.pair_terms" in full_labels
    @test "diatomic.packet.total" in legacy_labels

    full_report = nested_fixed_block_timing_report(timed_full)
    legacy_report = nested_fixed_block_timing_report(timed_legacy.timings)
    @test occursin("fixed_block.total", full_report)
    @test occursin("shell_layer.nonpacket", full_report)
    @test occursin("sequence_merge.nonpacket", full_report)
    @test occursin("diatomic.packet.gaussian_terms", full_report)
    @test occursin("diatomic.packet.pair_terms", full_report)
    @test occursin("diatomic.packet.total", legacy_report)
end
