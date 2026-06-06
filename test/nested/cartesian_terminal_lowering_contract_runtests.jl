using Test
using GaussletBases

const CSHForLowering = GaussletBases.CartesianShellification
const CTL = GaussletBases.CartesianTerminalLowering
const CPBForLowering = GaussletBases.CartesianCPB

function _terminal_lowering_module_source_files()
    module_dir =
        normpath(joinpath(@__DIR__, "..", "..", "src", "cartesian_terminal_lowering"))
    return sort(collect(
        joinpath(module_dir, file)
        for file in readdir(module_dir)
        if endswith(file, ".jl")
    ))
end

function _terminal_lowering_shell_plan()
    parent_axes = (collect(1:9), collect(1:9), collect(1:13))
    return CSHForLowering.shellify(
        parent_axes,
        ((5, 5, 4), (5, 5, 10)),
        CSHForLowering.AtomOutwardShellification(;
            core_side = 3,
            q = 3,
            bond_axis = :z,
        ),
    )
end

function _terminal_lowering_distorted_shell_plan()
    parent_axes = (collect(1:21), collect(1:7), collect(1:7))
    return CSHForLowering.shellify(
        parent_axes,
        ((3, 4, 4), (19, 4, 4)),
        CSHForLowering.AtomOutwardShellification(;
            core_side = 3,
            q = 3,
            bond_axis = :x,
        ),
    )
end

function _lowering_contract_count(plan, kind::Symbol)
    return count(contract -> CTL.lowering_kind(contract) == kind, CTL.contracts(plan))
end

@testset "CartesianTerminalLowering dependency direction" begin
    source_files = _terminal_lowering_module_source_files()
    @test !isempty(source_files)
    for path in source_files
        text = read(path, String)
        @test !occursin("parentmodule(@__MODULE__)", text)
        @test !occursin("GB._", text)
        @test !occursin("_cartesian_", text)
        @test !occursin("_nested_", text)
    end
end

@testset "CartesianTerminalLowering contracts" begin
    shell_plan = _terminal_lowering_shell_plan()
    terminal_region_count = length(CSHForLowering.terminal_regions(shell_plan))

    wl_plan = CTL.lower_terminal_regions(shell_plan, CTL.WhiteLindseyLowering())
    @test wl_plan.policy isa CTL.WhiteLindseyLowering
    @test length(CTL.contracts(wl_plan)) == terminal_region_count
    @test CTL.summary(wl_plan).all_terminal_regions_have_selected_contract
    @test !CTL.summary(wl_plan).materialized
    @test !CTL.summary(wl_plan).retained_spaces_materialized
    @test !CTL.summary(wl_plan).operator_blocks_materialized

    wl_shell_contract = first(
        contract for contract in CTL.contracts(wl_plan)
        if CTL.lowering_kind(contract) == :white_lindsey_boundary_strata
    )
    @test length(CTL.source_cpbs(wl_shell_contract)) == 26
    @test wl_shell_contract.metadata.facet_count == 6
    @test wl_shell_contract.metadata.edge_count == 12
    @test wl_shell_contract.metadata.corner_count == 8
    @test sum(CPBForLowering.support_count, CTL.source_cpbs(wl_shell_contract); init = 0) ==
          wl_shell_contract.metadata.shell_support_count
    wl_available_pqs = first(
        contract for contract in CTL.available_contracts(wl_plan)
        if CTL.lowering_kind(contract) == :pqs_filled_source_cpb
    )
    @test isnothing(wl_available_pqs.metadata.q)
    @test wl_available_pqs.metadata.parameter_status ==
          :available_but_unparameterized

    pqs_plan = CTL.lower_terminal_regions(shell_plan, CTL.PQSLowering(q = 3))
    @test pqs_plan.policy isa CTL.PQSLowering
    @test length(CTL.contracts(pqs_plan)) == terminal_region_count
    @test _lowering_contract_count(pqs_plan, :white_lindsey_boundary_strata) == 0
    @test _lowering_contract_count(pqs_plan, :pqs_filled_source_cpb) >= 1

    pqs_shell_contract = first(
        contract for contract in CTL.contracts(pqs_plan)
        if CTL.lowering_kind(contract) == :pqs_filled_source_cpb
    )
    @test length(CTL.source_cpbs(pqs_shell_contract)) == 1
    @test CPBForLowering.codimension(only(CTL.source_cpbs(pqs_shell_contract))) == 0
    @test pqs_shell_contract.retained_rule == :pqs_boundary_comx_product_modes
    @test pqs_shell_contract.realization_rule == :shell_projection_lowdin
    @test pqs_shell_contract.metadata.q == 3

    distorted_plan = CTL.lower_terminal_regions(
        _terminal_lowering_distorted_shell_plan(),
        CTL.WhiteLindseyLowering(),
    )
    @test _lowering_contract_count(distorted_plan, :distorted_product_box_comx) == 1
    distorted_contract = first(
        contract for contract in CTL.contracts(distorted_plan)
        if CTL.lowering_kind(contract) == :distorted_product_box_comx
    )
    @test length(CTL.source_cpbs(distorted_contract)) == 1
    @test distorted_contract.retained_rule == :distorted_product_comx_all_axes
    @test distorted_contract.metadata.source_mode_shape == (7, 3, 3)
end
