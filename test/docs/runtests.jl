@testset "Documentation consistency" begin
    read_doc(parts...) = read(joinpath(_PROJECT_ROOT, parts...), String)
    contains_all(text, phrases...) = all(phrase -> occursin(phrase, text), phrases)
    contains_all_lower(text, phrases...) = begin
        lowered = lowercase(text)
        all(phrase -> occursin(lowercase(phrase), lowered), phrases)
    end

    design = read_doc("DESIGN.md")
    readme = read_doc("README.md")
    status = read_doc("STATUS.md")
    roadmap = read_doc("ROADMAP.md")

    docs_architecture = read_doc("docs", "architecture.md")
    current_atomic_branch = read_doc("docs", "current_atomic_branch.md")
    current_ordinary_branch = read_doc("docs", "current_ordinary_branch.md")
    atomic_mean_field_supporting = read_doc("docs", "atomic_mean_field_supporting_notes.md")
    ordinary_pgdg_supporting = read_doc("docs", "ordinary_pgdg_supporting_notes.md")
    atomic_direct_note = read_doc("docs", "atomic_ida_direct.md")
    ordinary_pgdg_note = read_doc("docs", "ordinary_pgdg_decision.md")
    global_map_note = read_doc("docs", "global_map_local_contraction.md")
    leaf_pgdg_note = read_doc("docs", "leaf_pgdg_1d.md")

    docs_project = read_doc("docs", "Project.toml")
    docs_make = read_doc("docs", "make.jl")
    docs_workflow = read_doc(".github", "workflows", "docs.yml")

    docs_site_index = read_doc("docs", "src", "index.md")
    docs_site_manual = read_doc("docs", "src", "manual", "index.md")
    docs_site_examples = read_doc("docs", "src", "howto", "example_guide.md")
    docs_site_reference_index = read_doc("docs", "src", "reference", "index.md")
    docs_site_reference_bases = read_doc("docs", "src", "reference", "bases_and_mappings.md")
    docs_site_reference_ops = read_doc("docs", "src", "reference", "operators_and_diagnostics.md")
    docs_site_reference_atomic = read_doc("docs", "src", "reference", "atomic_and_ordinary.md")
    docs_site_reference_export = read_doc("docs", "src", "reference", "export.md")
    docs_site_developer = read_doc("docs", "src", "developer", "index.md")
    docs_site_developer_notes = read_doc("docs", "src", "developer", "supporting_notes.md")
    docs_site_architecture = read_doc("docs", "src", "developer", "architecture.md")
    docs_site_atomic = read_doc("docs", "src", "explanations", "current_atomic_branch.md")
    docs_site_ordinary = read_doc("docs", "src", "explanations", "current_ordinary_branch.md")
    docs_site_angular_track = read_doc("docs", "src", "explanations", "angular_research_track.md")

    @testset "Root Docs Authority And Story" begin
        @test contains_all(
            readme,
            "# GaussletBases.jl",
            "mature **radial / atomic workflow**",
            "experimental **angular and advanced research track**",
            "integral-diagonal approximation (IDA)",
            "## Best first path through the repository",
            "## Documentation map",
            "https://srwhite59.github.io/GaussletBases.jl/dev/manual/",
            "https://srwhite59.github.io/GaussletBases.jl/dev/reference/",
            "https://srwhite59.github.io/GaussletBases.jl/dev/developer/",
        )
        @test contains_all_lower(
            readme,
            "ordinary cartesian workflow",
            "bond-aligned diatomic nested fixed-source / fixed-block",
            "geometry payload support",
        )
        @test !occursin("Gausslets.jl", readme)

        @test contains_all(
            status,
            "# Current Status",
            "## Mature",
            "## Real but experimental",
            "## Legacy / quarantined",
            "### Radial / atomic line",
            "### Bond-aligned diatomic workflow",
            "### Angular research track",
            "### Flat supporting-note history in `docs/`",
        )
        @test contains_all_lower(
            status,
            "exact cartesian overlap / projector / transfer primitives",
            "supported one-build source reuse",
            "experimental chain / square-lattice nested producers",
            "old 1d comx-cleaned hybrid ordinary route",
        )

        @test contains_all(
            roadmap,
            "# Roadmap",
            "## Current center of gravity",
            "## Highest-value next work",
            "### 2. Deepen the one-center Cartesian line",
            "### 3. Mature the bond-aligned diatomic workflow",
        )
        @test contains_all_lower(
            roadmap,
            "one-center nested cartesian and bond-aligned diatomic workflow support",
            "supported source reuse",
            "compact-only cleanup",
            "experimental chain / square-lattice promotion",
        )

        @test contains_all(
            docs_architecture,
            "# Architecture",
            "This flat `docs/` file is no longer the current authority.",
            "`docs/src/developer/architecture.md`",
            "`docs/src/explanations/current_atomic_branch.md`",
            "`docs/src/explanations/current_ordinary_branch.md`",
        )

        @test contains_all(
            current_atomic_branch,
            "# Current Atomic Branch",
            "This flat `docs/` file is no longer the current authority.",
            "`docs/src/explanations/current_atomic_branch.md`",
            "`docs/src/developer/architecture.md`",
            "`docs/src/reference/index.md`",
        )

        @test contains_all(
            current_ordinary_branch,
            "# Current Ordinary Branch",
            "This flat `docs/` file is no longer the current authority.",
            "`docs/src/explanations/current_ordinary_branch.md`",
            "`docs/src/developer/architecture.md`",
            "`docs/src/reference/index.md`",
        )

        @test !occursin("primitive_kinetic_matrix", design)
        @test !occursin("CombinedMapping", design)
        @test !occursin("ScaledMapping", design)
        @test !occursin("NoDiagonalApproximation", design)
    end

    @testset "Rendered Docs Navigation And Authority" begin
        @test contains_all(
            docs_site_index,
            "# GaussletBases.jl",
            "## Start here",
            "## Primary documents",
            "## Manual first, Reference second",
            "## Current scope",
            "[Manual](manual/index.md)",
            "[Algorithms](algorithms/index.md)",
            "[Examples](examples/index.md)",
            "[Reference](reference/index.md)",
            "[Developer Notes](developer/index.md)",
        )
        @test contains_all_lower(
            docs_site_index,
            "mature radial / atomic workflow",
            "ordinary cartesian mapped/hybrid workflow",
            "advanced/research line",
        )

        @test contains_all(
            docs_site_manual,
            "# Manual",
            "## Recommended reading order",
            "## Branch-specific paths",
            "## If you want more depth later",
            "[Current atomic branch](../explanations/current_atomic_branch.md)",
            "[Current ordinary branch](../explanations/current_ordinary_branch.md)",
            "[Angular research track](../explanations/angular_research_track.md)",
            "[Developer Notes](../developer/index.md)",
        )
        @test contains_all_lower(
            docs_site_manual,
            "main user-facing manual",
            "recommended starting point",
            "ordinary cartesian branch",
        )

        @test contains_all(
            docs_site_examples,
            "# Example guide",
            "## Core starting sequence",
            "## Radial and atomic sequence",
            "## Ordinary Cartesian sequence",
            "## Primitive and hierarchy sequence",
            "[Current atomic branch](../explanations/current_atomic_branch.md)",
            "[Current ordinary branch](../explanations/current_ordinary_branch.md)",
            "[Developer Notes](../developer/index.md)",
            "38_qiu_white_reference_vee.jl",
        )
        @test contains_all_lower(
            docs_site_examples,
            "legacy/internal experimental regressions",
            "producer-side only",
            "explicitly experimental",
        )

        @test contains_all(
            docs_site_reference_index,
            "# Reference overview",
            "[Bases and mappings](bases_and_mappings.md)",
            "[Operators and diagnostics](operators_and_diagnostics.md)",
            "[Atomic and ordinary workflows](atomic_and_ordinary.md)",
            "[Export layer](export.md)",
            "[Manual](../manual/index.md)",
        )
        @test contains_all_lower(
            docs_site_reference_index,
            "first curated api-reference slice",
            "workflow and interpretation",
            "api questions",
        )

        @test contains_all(
            docs_site_developer,
            "# Developer Notes",
            "## Main developer-facing pages",
            "[Architecture and current direction](architecture.md)",
            "[PGDG Cartesian efficiency contract](pgdg_cartesian_efficiency_contract.md)",
            "[Numerical contracts](numerical_contracts.md)",
            "[Supporting note map](supporting_notes.md)",
        )
        @test contains_all_lower(
            docs_site_developer,
            "lower-priority development",
            "package-shape and architecture view",
            "supporting note chains behind the current branch interpretations",
        )

        @test contains_all(
            docs_site_developer_notes,
            "# Supporting note map",
            "## Current grouped supporting-note entry points",
            "docs/atomic_mean_field_supporting_notes.md",
            "docs/ordinary_pgdg_supporting_notes.md",
            "docs/documenter_transition_plan.md",
            "[Manual](../manual/index.md)",
            "[Reference](../reference/index.md)",
        )
        @test contains_all_lower(
            docs_site_developer_notes,
            "supporting material rather than as the main user manual",
            "recommended reading order",
            "scientific development record",
        )

        @test contains_all(
            docs_site_architecture,
            "# Architecture and current direction",
            "## Broad foundation",
            "## Mature public-facing path",
            "## Atomic line",
            "## Ordinary line",
            "## Current bottom line",
        )
        @test contains_all_lower(
            docs_site_architecture,
            "ordinary gausslets are the broad foundation",
            "radial gausslets are the mature current workflow",
            "primitive layers and contraction are the structural bridge to later work",
        )
    end

    @testset "Current Branch Contracts" begin
        @test contains_all(
            docs_site_atomic,
            "# Current atomic branch",
            "## What the atomic branch is today",
            "## Angular research track",
            "## Documentation authority for this branch",
            "[Developer Notes](../developer/index.md)",
            "[Angular research track](angular_research_track.md)",
        )
        @test contains_all_lower(
            docs_site_atomic,
            "shortest current user-facing status read",
            "density-density / ida",
            "solver-facing export is already supported",
            "supporting history",
        )
        @test occursin("HFDMRG", docs_site_atomic)

        @test contains_all(
            docs_site_ordinary,
            "# Current ordinary branch",
            "## What the ordinary branch is today",
            "## Documentation authority for this branch",
            "[Reference](../reference/index.md)",
            "[Algorithms](../algorithms/index.md)",
            "[Developer Notes](../developer/index.md)",
        )
        @test contains_all(
            docs_site_ordinary,
            ":numerical_reference",
            ":pgdg_localized_experimental",
            "AsinhMapping",
        )
        @test contains_all_lower(
            docs_site_ordinary,
            "pgdg-style analytic route is good enough on the mapped ordinary backbone",
            "bond-aligned diatomic molecular supplement direct-product and prebuilt",
            "nested fixed-block qiu-white routes are now also pgdg-capable",
            "nested source-building front doors",
            "experimental chain/square nested qw source wrappers",
        )

        @test contains_all(
            docs_site_angular_track,
            "# Angular Research Track",
            "HFDMRG",
            "sphere_point_set_orders",
            "build_atomic_injected_angular_hfdmrg_payload",
            "build_atomic_injected_angular_small_ed_benchmark",
            "write_angular_benchmark_exact_hamv6_jld2",
        )
        @test contains_all_lower(
            docs_site_angular_track,
            "manuscript-facing and experimental",
            "shell-local experimental construction",
            "exact common low",
            "not a frozen public api",
        )
    end

    @testset "Docs Build And Reference Surfaces" begin
        @test occursin("Documenter", docs_project)

        @test contains_all(
            docs_make,
            "makedocs",
            "doctest = true",
            "checkdocs = :none",
            "deploydocs(",
            "prettyurls = DOCS_CI",
            "\"Manual\"",
            "\"Examples\"",
            "\"Reference\"",
            "\"Developer Notes\"",
        )

        @test contains_all(
            docs_workflow,
            "name: Docs",
            "julia --project=docs docs/make.jl",
            "Pkg.develop(PackageSpec(path=pwd()))",
            "GITHUB_TOKEN",
        )

        @test contains_all(
            docs_site_reference_bases,
            "jldoctest",
            "UniformBasisSpec",
            "MappedUniformBasisSpec",
            "white_lindsey_atomic_mapping",
            "CombinedInvsqrtMapping",
        )

        @test contains_all(
            docs_site_reference_ops,
            "jldoctest",
            "basis_diagnostics",
            "radial_quadrature",
            "kinetic_matrix",
            "atomic_operators",
        )

        @test contains_all(
            docs_site_reference_atomic,
            "jldoctest",
            "AtomicIDAOperators",
            "mapped_ordinary_one_body_operators",
            "ordinary_cartesian_qiu_white_operators",
            "ordinary_cartesian_vee_expectation",
        )

        @test contains_all(
            docs_site_reference_export,
            "atomic_ida_density_interaction_matrix",
            "fullida_dense_payload",
            "sliced_ham_payload",
            "atomic_hamv6_payload",
            "angular_benchmark_exact_hamv6_payload",
            "experimental_homonuclear_chain_nested_dense_payload",
        )
        @test contains_all_lower(
            docs_site_reference_export,
            "producer-side only",
            "hfdmrg.solve_hfdmrg(...)",
            "experimental rather than a settled public standard",
        )
    end

    @testset "Representative Supporting Note Markers" begin
        @test contains_all_lower(
            atomic_mean_field_supporting,
            "recommended supporting-note order",
            "atomic_ida_direct.md",
            "atomic_ida_spin_fock.md",
        )
        @test contains_all_lower(
            ordinary_pgdg_supporting,
            "recommended supporting-note order",
            "ordinary_pgdg_decision.md",
            "ordinary_pgdg_backend_pivot.md",
            "asinhmapping",
        )

        @test startswith(atomic_direct_note, "> **Status:** supporting note.")
        @test contains_all_lower(
            atomic_direct_note,
            "direct/hartree",
            "radial-diagonal",
            "dense `angular_kernel`",
        )

        @test startswith(ordinary_pgdg_note, "> **Status:** supporting development note.")
        @test contains_all_lower(
            ordinary_pgdg_note,
            "primitive/contraction-level analytic prototype",
            "comx/localization",
            "135",
            "115",
        )

        @test startswith(global_map_note, "> **Note for new users:**")
        @test contains_all(
            global_map_note,
            "GlobalMappedPrimitiveLayer1D",
            "global map + local shell or box contraction",
        )

        @test startswith(leaf_pgdg_note, "> **Note for new users:**")
        @test contains_all(
            leaf_pgdg_note,
            "LeafLocalPGDG1D",
            "hierarchy-driven local basis generation",
            "no historical nested-driver port",
            "augment_leaf_pgdg",
        )
    end
end
