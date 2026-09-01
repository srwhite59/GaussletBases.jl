@testset "Documentation consistency" begin
    include(joinpath(_PROJECT_ROOT, "docs", "check_manager_log.jl"))
    include(joinpath(_PROJECT_ROOT, "docs", "make.jl"))

    read_doc(parts...) = read(joinpath(_PROJECT_ROOT, parts...), String)
    contains_all(text, phrases...) = all(phrase -> occursin(phrase, text), phrases)
    contains_all_lower(text, phrases...) = begin
        lowered = lowercase(text)
        all(phrase -> occursin(lowercase(phrase), lowered), phrases)
    end

    design = read_doc("DESIGN.md")
    readme = read_doc("README.md")
    root_project = read_doc("Project.toml")
    changelog = read_doc("CHANGELOG.md")
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
    ci_workflow = read_doc(".github", "workflows", "ci.yml")
    authority_check = read_doc("docs", "check_cartesian_authority.jl")

    docs_site_index = read_doc("docs", "src", "index.md")
    docs_site_manual = read_doc("docs", "src", "manual", "index.md")
    docs_site_examples = read_doc("docs", "src", "howto", "example_guide.md")
    docs_flat_examples = read_doc("docs", "example_guide.md")
    docs_site_examples_landing = read_doc("docs", "src", "examples", "index.md")
    docs_site_pqs = read_doc("docs", "src", "manual", "projected_q_shells.md")
    docs_site_external_gto = read_doc("docs", "src", "manual", "external_cartesian_gto_transfer.md")
    docs_site_diatomic_box = read_doc(
        "docs", "src", "algorithms", "cartesian_nested_diatomic_box_policy.md")
    docs_site_diatomic_distortion = read_doc(
        "docs", "src", "algorithms", "cartesian_nested_diatomic_coordinate_distortion.md")
    pqs_example = read_doc("examples", "39_pqs_h2plus.jl")
    docs_site_reference_index = read_doc("docs", "src", "reference", "index.md")
    docs_site_reference_bases = read_doc("docs", "src", "reference", "bases_and_mappings.md")
    docs_site_reference_ops = read_doc("docs", "src", "reference", "operators_and_diagnostics.md")
    docs_site_reference_atomic = read_doc("docs", "src", "reference", "atomic_and_ordinary.md")
    docs_site_reference_export = read_doc("docs", "src", "reference", "export.md")
    v0_2_exports = (
        :PQSH2PlusRow, :PQSH2PlusComparison, :pqs_h2plus_comparison,
        :ExactRepresentedHartreeField, :FittedReferenceHartreeField,
        :ScreenedHartreeCorrection, :screened_hartree_correction,
        :screened_hartree_delta_one_body, :screened_hartree_energy_constant,
        :screened_hartree_consistency_error, :screened_hartree_field_kind,
        :cartesian_residual_gto_mwg_system,
    )
    external_gto_exports = (
        :ExternalGTOOrbitalSpinBlock, :ExternalGTOOrbitalPacket,
        :ExternalGTOOrbitalSpinImport, :ExternalGTOOrbitalImportResult,
        :import_external_gto_orbitals, :read_external_cartesian_gto_packet,
        :closest_external_gto_determinant,
    )
    foundation_exports = (
        :uofx, :xofu, :dudx, :du2dx2, :basis_spec,
        :family, :mapping, :centers, :reference_centers, :integral_weights,
    )
    function_stencil_exports = (
        :value, :direct_value, :derivative, :center, :reference_center,
        :integral_weight, :stencil, :coefficients, :terms,
    )
    partition_leaf_exports = (
        :boxes, :leaf_boxes, :box_indices, :box_level, :box_parent, :box_children,
        :box_block, :box_coupling, :leaf_primitive_indices, :primitive_origins,
        :primitive_leaf_boxes, :leaf_contractions)
    partition_leaf_context = (
        :BasisBox1D, :BasisPartition1D, :basis_partition, :HierarchicalBasisBox1D,
        :HierarchicalBasisPartition1D, :hierarchical_partition, :refine_partition,
        :LeafGaussianSpec1D, :LeafLocalPGDG1D, :build_leaf_pgdg, :augment_leaf_pgdg,
        :GlobalMappedPrimitiveLayer1D, :build_global_mapped_primitive_layer,
        :LeafBoxContraction1D, :LeafBoxContractionLayer1D, :contract_leaf_boxes)
    atomic_ida_reference_exports = (
        :orbitals, :two_electron_states, :radial_multipole, :gaunt_coefficient,
        :gaunt_tensor, :angular_kernel, :apply_overlap, :apply_hamiltonian,
        :ground_state_energy, :lanczos_ground_state)
    atomic_ida_reference_context = (
        :AtomicOrbital, :AtomicIDAOperators, :atomic_ida_operators,
        :AtomicIDATwoElectronState, :AtomicIDATwoElectronProblem,
        :atomic_ida_two_electron_problem)
    supported_surface_exports = (
        :BondAlignedDiatomicQWBasis3D, :CoulombGaussianExpansion, :basis_metadata,
        :cartesian_base_hamiltonian, :external_gto_ordering_fingerprint,
        :external_gto_overlap_fingerprint)
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
            "https://srwhite59.github.io/GaussletBases.jl/stable/manual/",
            "https://srwhite59.github.io/GaussletBases.jl/stable/reference/",
            "https://srwhite59.github.io/GaussletBases.jl/stable/developer/",
            "[Projected q-shells](https://srwhite59.github.io/GaussletBases.jl/stable/manual/projected_q_shells/)",
            "[Reference-density Hartree screening](https://srwhite59.github.io/GaussletBases.jl/stable/manual/reference_density_hartree_screening/)",
            "[External Cartesian GTO transfer](https://srwhite59.github.io/GaussletBases.jl/stable/manual/external_cartesian_gto_transfer/)",
        )
        @test contains_all_lower(
            readme,
            "ordinary cartesian workflow",
            "bond-aligned diatomic nested fixed-source / fixed-block",
            "geometry payload support",
            "checkpoint-only exporter",
            "pyscf and numpy are optional dependencies",
        )
        @test !occursin("Gausslets.jl", readme)
        @test occursin("rev = \"v0.2.0\"", readme)
        @test !occursin("https://srwhite59.github.io/GaussletBases.jl/dev/", readme)
        @test occursin("version = \"0.2.0\"", root_project)
        @test startswith(changelog, "# Changelog\n\n## v0.2.0")
        @test contains_all(changelog, "cartesian_residual_gto_mwg_system",
            "strict reader", "checkpoint-only PySCF exporter",
            "caller-thresholded closest-determinant", "Supported-floor, PQS, and Screening",
            "invalid package exports", "support-local shell-seed construction",
            "exact-order four-element", "path-aware CI fail closed",
            "all unchanged release assertions", "preserve accepted public numerics",
            "## v0.2.0-rc2", "## v0.2.0-rc1")

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
        @test contains_all(docs_site_pqs, "# Projected q-shells (PQS)",
            "PQS shell construction", "examples/39_pqs_h2plus.jl")
        @test occursin("\"External Cartesian GTO transfer\" => \"manual/external_cartesian_gto_transfer.md\"", docs_make)
        @test occursin("[External Cartesian GTO transfer](external_cartesian_gto_transfer.md)", docs_site_manual)
        @test all(name -> name in names(GaussletBases), external_gto_exports)
        @test all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_export)\n"),
            external_gto_exports)
        @test !occursin("\ngto_overlap_matrix\n", "\n$(docs_site_reference_export)\n")
        @test contains_all(docs_site_external_gto, "checkpoint-only exporter",
            "unchanged raw projection", "minimum_gram_eigenvalue",
            "producer attestation", "cannot recover", "missing direction")
        @test all(term -> !occursin(term, docs_site_external_gto),
            ("/Users/", "Dropbox", "CloudStorage", "REQ-", "C2", "GaussletBases._"))
        @test contains_all(docs_site_examples_landing, "39_pqs_h2plus.jl",
            "Projected q-shells (PQS)")
        @test contains_all(pqs_example, "cartesian_base_hamiltonian",
            "CartesianIDAHamiltonian", "one_body_hamiltonian", "nuclear_repulsion")
        @test !occursin("GaussletBases._", pqs_example)
        @test all(name -> name in names(GaussletBases), v0_2_exports)
        @test all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_export)\n"),
            v0_2_exports)
        for page in (docs_site_diatomic_box, docs_site_diatomic_distortion)
            @test all(path -> !occursin(path, page),
                ("/Users/srw", "Dropbox", "CloudStorage"))
        end
        for guide in (docs_site_examples, docs_flat_examples)
            @test contains_all_lower(guide, "28_ordinary_one_body_fidelity.jl",
                "internal", "non-contractual")
        end
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
            "check_manager_log.jl",
            "ManagerLogPolicy.check_live_log()",
            "\"Manual\"",
            "\"Examples\"",
            "\"Reference\"",
            "\"Developer Notes\"",
        )

        @test contains_all(
            docs_workflow,
            "name: Docs",
            "tags:",
            "- \"v*\"",
            "cartesian-authority:",
            "needs: cartesian-authority",
            "startsWith(github.ref, 'refs/tags/v')",
            "julia --project=docs docs/make.jl",
            "Pkg.develop(PackageSpec(path=pwd()))",
            "docs/check_cartesian_authority.jl --check",
            "docs/check_cartesian_authority.jl --self-test",
        )
        @test _documentation_target("pull_request", "refs/pull/1/merge") ==
              (:pull_request, "dev")
        @test _documentation_target("push", "refs/heads/main") == (:dev, "dev")
        @test _documentation_target("push", "refs/tags/v0.2.0-rc1") ==
              (:tag, "v0.2.0-rc1")
        @test _documentation_target("push", "refs/tags/v0.2.0-rc2") ==
              (:tag, "v0.2.0-rc2")
        @test _documentation_target("push", "refs/tags/v0.2.0") ==
              (:tag, "v0.2.0")
        @test "$(_DOCS_SITE)/$(_documentation_target("push", "refs/heads/main")[2])/" ==
              "https://srwhite59.github.io/GaussletBases.jl/dev/"
        @test "$(_DOCS_SITE)/$(_documentation_target("push", "refs/tags/v0.2.0-rc1")[2])/" ==
              "https://srwhite59.github.io/GaussletBases.jl/v0.2.0-rc1/"
        @test "$(_DOCS_SITE)/$(_documentation_target("push", "refs/tags/v0.2.0-rc2")[2])/" ==
              "https://srwhite59.github.io/GaussletBases.jl/v0.2.0-rc2/"
        @test "$(_DOCS_SITE)/$(_documentation_target("push", "refs/tags/v0.2.0")[2])/" ==
              "https://srwhite59.github.io/GaussletBases.jl/v0.2.0/"
        for ref in ("refs/tags/v0.2", "refs/tags/release-v0.2.0",
                    "refs/tags/v0.2.0+build", "refs/tags/v01.2.3")
            @test_throws ErrorException _documentation_target("push", ref)
        end
        @test_throws ErrorException _documentation_target("schedule", "refs/heads/main")
        @test !isempty(VersionNumber("0.2.0-rc1").prerelease)
        @test isempty(VersionNumber("0.2.0").prerelease)
        @test contains_all(docs_make,
            "DOCS_DEPLOY && DOCS_TARGET[1] ∉ (:dev, :tag)",
            "canonical = DOCS_CANONICAL")
        @test length(findall("GITHUB_TOKEN:", docs_workflow)) == 1
        @test length(findall("GAUSSLETBASES_DOCS_DEPLOY", docs_workflow)) == 1
        @test length(findall("contents: write", docs_workflow)) == 1
        @test length(findall("contents: read", docs_workflow)) == 2
        docs_only_path = path -> path == "AGENTS.md" || startswith(path, "docs/") ||
                                   startswith(path, "test/docs/") ||
                                   path == ".github/workflows/docs.yml"
        @test all(docs_only_path, ("AGENTS.md", "docs/page.md", "test/docs/runtests.jl",
                                  ".github/workflows/docs.yml"))
        @test !any(docs_only_path, ("Project.toml", "README.md", "CHANGELOG.md",
                                   ".github/workflows/ci.yml", "src/GaussletBases.jl"))
        @test contains_all(ci_workflow, "pull_request:", "fetch-depth: 0",
            "github.event.before", "git cat-file -e \"\${before}^{commit}\"",
            "git diff --name-only --no-renames", "scope=full", "scope=docs_only",
            "AGENTS\\.md|docs/.+|test/docs/.+|\\.github/workflows/docs\\.yml")
        @test contains_all(ci_workflow, "always()", "needs.classify.result != 'success'",
            "needs.classify.outputs.scope != 'docs_only'",
            "needs.classify.result == 'success'",
            "needs.classify.outputs.scope == 'docs_only'")
        @test occursin("if: \${{ always() && !startsWith(github.ref, 'refs/tags/') }}",
            ci_workflow)
        full_step_condition =
            "if: \${{ needs.classify.result != 'success' || needs.classify.outputs.scope != 'docs_only' }}"
        @test length(findall(full_step_condition, ci_workflow)) == 6
        @test contains_all(ci_workflow, "Documentation-only marker",
            "Documentation-only change; numerical gate skipped.")
        @test occursin("julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'",
            ci_workflow)
        @test all(gate -> length(findall("gate: $gate", ci_workflow)) == 1,
            ("Supported floor", "PQS paper", "Screening paper"))
        @test occursin("groups: core,ida,cartesian,examples,radial,misc,angular_public", ci_workflow)
        @test !occursin(r"groups:[^\n]*[, ]angular(?:,|\n)", ci_workflow)
        @test contains_all(ci_workflow, "tag-identity-install:",
            "refs/tags/\${tag}:refs/tags/\${tag}", "refs/tags/\${tag}^{tag}",
            "refs/tags/\${tag}^{commit}", "Pkg.add(url=", "GITHUB_REF_NAME")
        @test !occursin("workflow_run", ci_workflow)
        @test _DOCS_VERSIONS == Any[
            "stable" => "v^",
            "v#.#",
            "v0.2.0-rc2" => "v0.2.0-rc2",
            "v0.2.0-rc1" => "v0.2.0-rc1",
            "dev" => "dev",
        ]
        @test occursin("versions = _DOCS_VERSIONS", docs_make)
        mktempdir() do directory
            mkpath(joinpath(directory, "dev"))
            mkpath(joinpath(directory, "v0.2.0-rc2"))
            mkpath(joinpath(directory, "v0.2.0-rc1"))
            script = """
                using Documenter
                entries, symlinks = Documenter.Writers.HTMLWriter.expand_versions(
                    $(repr(directory)), $(repr(_DOCS_VERSIONS)))
                @assert entries == ["v0.2.0-rc2", "v0.2.0-rc1", "dev"]
                @assert !any(pair -> pair.first == pair.second, symlinks)
                @assert !any(pair -> pair.first == "stable", symlinks)
                """
            docs_environment = joinpath(_PROJECT_ROOT, "docs")
            command = `$(Base.julia_cmd()) --project=$docs_environment --startup-file=no -e $script`
            @test success(command)
        end
        mktempdir() do directory
            foreach(entry -> mkpath(joinpath(directory, entry)),
                ("dev", "v0.2.0", "v0.2.0-rc2", "v0.2.0-rc1"))
            script = """
                using Documenter
                entries, symlinks = Documenter.Writers.HTMLWriter.expand_versions(
                    $(repr(directory)), $(repr(_DOCS_VERSIONS)))
                @assert entries == ["stable", "v0.2", "v0.2.0-rc2", "v0.2.0-rc1", "dev"]
                @assert symlinks == ["stable" => "v0.2.0", "v0.2" => "v0.2.0"]
                """
            command = `$(Base.julia_cmd()) --project=$(joinpath(_PROJECT_ROOT, "docs")) --startup-file=no -e $script`
            @test success(command)
        end
        @test !occursin("check_cartesian_authority", docs_make)

        @test ManagerLogPolicy.check_live_log() <= ManagerLogPolicy.LIVE_LOG_MAX_LINES
        mktempdir() do directory
            fixture = joinpath(directory, "manager_log.md")
            open(fixture, "w") do io
                for line in 1:ManagerLogPolicy.LIVE_LOG_MAX_LINES
                    println(io, "line $(line)")
                end
            end
            @test ManagerLogPolicy.check_live_log(fixture) ==
                  ManagerLogPolicy.LIVE_LOG_MAX_LINES
            open(fixture, "a") do io
                println(io, "over limit")
            end
            @test_throws ErrorException ManagerLogPolicy.check_live_log(fixture)
        end

        @test contains_all(
            authority_check,
            "module CartesianAuthority",
            "authority must set authoritative=true",
            "generated registry is stale or partially edited",
            "legacy authority artifact remains live",
            "--render",
        )

        @test contains_all(
            docs_site_reference_bases,
            "jldoctest",
            "UniformBasisSpec",
            "MappedUniformBasisSpec",
            "white_lindsey_atomic_mapping",
            "CombinedInvsqrtMapping",
        )
        @test all(name -> name in names(GaussletBases), foundation_exports)
        @test all(name -> Base.Docs.hasdoc(GaussletBases, name), foundation_exports)
        @test all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_bases)\n"),
            foundation_exports)
        @test all(name -> name in names(GaussletBases), function_stencil_exports)
        @test all(name -> Base.Docs.hasdoc(GaussletBases, name), function_stencil_exports)
        @test all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_bases)\n"),
            function_stencil_exports)
        @test all(name -> name in names(GaussletBases) && Base.Docs.hasdoc(GaussletBases, name),
            partition_leaf_exports)
        @test occursin("## Partitions and leaf-local layers", docs_site_reference_bases) &&
              all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_bases)\n"),
            (partition_leaf_exports..., partition_leaf_context...))

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
        @test all(name -> name in names(GaussletBases) && Base.Docs.hasdoc(GaussletBases, name),
            atomic_ida_reference_exports)
        @test all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_atomic)\n"),
            (atomic_ida_reference_exports..., atomic_ida_reference_context...))
        @test contains_all_lower(docs_site_reference_atomic, "read-only", "construction order",
            "rebuilt on demand", "not cached", "stored dense matrices", "fully reorthogonalized",
            "krylovdim", "value", "one-up/one-down", "not a production solver")
        @test all(name -> name in names(GaussletBases) && Base.Docs.hasdoc(GaussletBases, name),
            supported_surface_exports)
        @test !Base.Docs.hasdoc(GaussletBases, :AbstractBondAlignedOrdinaryQWBasis3D)
        @test occursin("\nbasis_metadata\n", "\n$(docs_site_reference_bases)\n")
        @test all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_ops)\n"),
            (:BondAlignedDiatomicQWBasis3D, :CoulombGaussianExpansion))
        @test all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_export)\n"),
            supported_surface_exports[4:6])
        @test contains_all_lower(docs_site_reference_bases, "owner-specific",
            "no universal field schema", "copy behavior")
        @test contains_all_lower(docs_site_reference_ops, "bond-aligned mapped-product", "homonuclear",
            "heteronuclear", "finite gaussian expansion", "not an exact operator", "universal-interval")
        @test contains_all_lower(docs_site_reference_export, "origin-centered hydrogen", "z axis",
            "strict packet-integrity", "not numerical", "permutation", "tolerance", "convention")

        geometry_inspection_exports = (
            :bond_aligned_diatomic_geometry_payload, :BondAlignedDiatomicGeometryPoint3D,
            :BondAlignedDiatomicGeometryNucleus3D, :BondAlignedDiatomicGeometryBox3D,
            :BondAlignedDiatomicGeometryPayload3D, :BondAlignedDiatomicGeometryPlaneSlice3D,
            :bond_aligned_diatomic_source_geometry_payload, :bond_aligned_diatomic_plane_slice)
        @test all(name -> name in names(GaussletBases) && Base.Docs.hasdoc(GaussletBases, name), geometry_inspection_exports)
        @test occursin("## Expert/experimental bond-aligned geometry inspection", docs_site_reference_bases) &&
              all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_bases)\n"), geometry_inspection_exports)
        @test contains_all_lower(docs_site_reference_bases, "inspection and visualization", "does not construct or mutate",
            "compressed fixed centers", "raw source-region points", "same basis geometry", "read-only inspection data",
            "not as a general molecular-geometry api", "plotting framework", "stable serialization format")
        @test length(filter(!=(:GaussletBases), Docs.undocumented_names(GaussletBases))) == 23

        @test contains_all(
            docs_site_reference_export,
            "atomic_ida_density_interaction_matrix",
            "fullida_dense_payload",
            "sliced_ham_payload",
            "atomic_hamv6_payload",
            "angular_benchmark_exact_hamv6_payload",
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

include(joinpath(@__DIR__, "cartesian_ham_builder_policy_runtests.jl"))
