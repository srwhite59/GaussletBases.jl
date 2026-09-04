@testset "Documentation consistency" begin
    include(joinpath(_PROJECT_ROOT, "docs", "check_manager_log.jl"))
    include(joinpath(_PROJECT_ROOT, "docs", "make.jl"))

    read_doc(parts...) = read(joinpath(_PROJECT_ROOT, parts...), String)
    contains_all(text, phrases...) = all(phrase -> occursin(phrase, text), phrases)

    design = read_doc("DESIGN.md")
    readme = read_doc("README.md")
    root_project = read_doc("Project.toml")
    changelog = read_doc("CHANGELOG.md")
    agents = read_doc("AGENTS.md")

    docs_architecture = read_doc("docs", "architecture.md")
    current_atomic_branch = read_doc("docs", "current_atomic_branch.md")
    current_ordinary_branch = read_doc("docs", "current_ordinary_branch.md")
    atomic_direct_note = read_doc("docs", "atomic_ida_direct.md")
    ordinary_pgdg_note = read_doc("docs", "ordinary_pgdg_decision.md")
    global_map_note = read_doc("docs", "global_map_local_contraction.md")
    leaf_pgdg_note = read_doc("docs", "leaf_pgdg_1d.md")

    docs_project = read_doc("docs", "Project.toml")
    docs_make = read_doc("docs", "make.jl")
    docs_workflow = read_doc(".github", "workflows", "docs.yml")
    ci_workflow = read_doc(".github", "workflows", "ci.yml")
    maintenance_workflow = read_doc(".github", "workflows", "cartesian-internal-maintenance.yml")
    test_runner = read_doc("test", "runtests.jl")
    authority_check = read_doc("docs", "check_cartesian_authority.jl")
    execution_whitelist = read_doc("docs", "src", "developer", "designs",
        "cartesian_hamiltonian_producer", "execution_whitelist.md")

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
    partition_leaf_context = (
        :BasisBox1D, :BasisPartition1D, :basis_partition, :HierarchicalBasisBox1D,
        :HierarchicalBasisPartition1D, :hierarchical_partition, :refine_partition,
        :LeafGaussianSpec1D, :LeafLocalPGDG1D, :build_leaf_pgdg, :augment_leaf_pgdg,
        :GlobalMappedPrimitiveLayer1D, :build_global_mapped_primitive_layer,
        :LeafBoxContraction1D, :LeafBoxContractionLayer1D, :contract_leaf_boxes)
    atomic_ida_reference_context = (
        :AtomicOrbital, :AtomicIDAOperators, :atomic_ida_operators,
        :AtomicIDATwoElectronState, :AtomicIDATwoElectronProblem,
        :atomic_ida_two_electron_problem)
    docs_site_developer = read_doc("docs", "src", "developer", "index.md")
    docs_site_developer_notes = read_doc("docs", "src", "developer", "supporting_notes.md")
    docs_site_atomic = read_doc("docs", "src", "explanations", "current_atomic_branch.md")
    docs_site_ordinary = read_doc("docs", "src", "explanations", "current_ordinary_branch.md")
    docs_site_angular_track = read_doc("docs", "src", "explanations", "angular_research_track.md")
    recommended_atomic_setup = read_doc("docs", "src", "howto", "recommended_atomic_setup.md")

    @testset "Root Docs Authority And Story" begin
        @test contains_all(
            readme,
            "# GaussletBases.jl",
            "https://srwhite59.github.io/GaussletBases.jl/stable/manual/",
            "https://srwhite59.github.io/GaussletBases.jl/stable/reference/",
            "https://srwhite59.github.io/GaussletBases.jl/stable/developer/",
            "[Projected q-shells](https://srwhite59.github.io/GaussletBases.jl/stable/manual/projected_q_shells/)",
            "[Reference-density Hartree screening](https://srwhite59.github.io/GaussletBases.jl/stable/manual/reference_density_hartree_screening/)",
            "[External Cartesian GTO transfer](https://srwhite59.github.io/GaussletBases.jl/stable/manual/external_cartesian_gto_transfer/)",
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
            docs_architecture,
            "`docs/src/developer/architecture.md`",
            "`docs/src/explanations/current_atomic_branch.md`",
            "`docs/src/explanations/current_ordinary_branch.md`",
        )

        @test contains_all(
            current_atomic_branch,
            "# Current Atomic Branch",
            "`docs/src/explanations/current_atomic_branch.md`",
            "`docs/src/developer/architecture.md`",
            "`docs/src/reference/index.md`",
        )

        @test contains_all(
            current_ordinary_branch,
            "# Current Ordinary Branch",
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
            "[Manual](manual/index.md)",
            "[Algorithms](algorithms/index.md)",
            "[Examples](examples/index.md)",
            "[Reference](reference/index.md)",
            "[Developer Notes](developer/index.md)",
        )
        @test contains_all(
            docs_site_manual,
            "# Manual",
            "[Current atomic branch](../explanations/current_atomic_branch.md)",
            "[Current ordinary branch](../explanations/current_ordinary_branch.md)",
            "[Angular research track](../explanations/angular_research_track.md)",
            "[Developer Notes](../developer/index.md)",
        )
        @test contains_all(
            docs_site_examples,
            "# Example guide",
            "[Current atomic branch](../explanations/current_atomic_branch.md)",
            "[Current ordinary branch](../explanations/current_ordinary_branch.md)",
            "[Developer Notes](../developer/index.md)",
            "38_qiu_white_reference_vee.jl",
        )
        @test contains_all(docs_site_pqs, "# Projected q-shells (PQS)",
            "PQS shell construction", "examples/39_pqs_h2plus.jl")
        @test occursin("\"External Cartesian GTO transfer\" => \"manual/external_cartesian_gto_transfer.md\"", docs_make)
        @test occursin("[External Cartesian GTO transfer](external_cartesian_gto_transfer.md)", docs_site_manual)
        @test all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_export)\n"),
            _DOCUMENTED_PUBLIC_SURFACE.external_gto)
        @test !occursin("\ngto_overlap_matrix\n", "\n$(docs_site_reference_export)\n")
        @test occursin(
            "# External Cartesian GTO Transfer",
            docs_site_external_gto,
        )
        @test all(term -> !occursin(term, docs_site_external_gto),
            ("/Users/", "Dropbox", "CloudStorage", "REQ-", "C2", "GaussletBases._"))
        @test contains_all(docs_site_examples_landing, "39_pqs_h2plus.jl",
            "Projected q-shells (PQS)")
        @test contains_all(pqs_example, "cartesian_base_hamiltonian",
            "CartesianIDAHamiltonian", "one_body_hamiltonian", "nuclear_repulsion")
        @test !occursin("GaussletBases._", pqs_example)
        @test all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_export)\n"),
            _DOCUMENTED_PUBLIC_SURFACE.v0_2)
        for page in (docs_site_diatomic_box, docs_site_diatomic_distortion)
            @test all(path -> !occursin(path, page),
                ("/Users/srw", "Dropbox", "CloudStorage"))
        end
        for guide in (docs_site_examples, docs_flat_examples)
            @test contains_all(lowercase(guide), "28_ordinary_one_body_fidelity.jl",
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
        @test contains_all(
            docs_site_developer,
            "[Architecture and current direction](architecture.md)",
            "[PGDG Cartesian efficiency contract](pgdg_cartesian_efficiency_contract.md)",
            "[Numerical contracts](numerical_contracts.md)",
            "[Supporting note map](supporting_notes.md)",
        )
        @test contains_all(
            docs_site_developer_notes,
            "docs/atomic_mean_field_supporting_notes.md",
            "docs/ordinary_pgdg_supporting_notes.md",
            "docs/documenter_transition_plan.md",
            "[Manual](../manual/index.md)",
            "[Reference](../reference/index.md)",
        )
    end

    @testset "Current Branch Contracts" begin
        @test contains_all(
            docs_site_atomic,
            "# Current atomic branch",
            "[Developer Notes](../developer/index.md)",
            "[Angular research track](angular_research_track.md)",
        )
        @test contains_all(
            docs_site_ordinary,
            "# Current ordinary branch",
            "[Reference](../reference/index.md)",
            "[Algorithms](../algorithms/index.md)",
            "[Developer Notes](../developer/index.md)",
        )
        @test contains_all(
            docs_site_angular_track,
            "# Angular Research Track",
            "sphere_point_set_orders",
            "build_atomic_injected_angular_hfdmrg_payload",
            "build_atomic_injected_angular_small_ed_benchmark",
            "write_angular_benchmark_exact_hamv6_jld2",
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
        @test occursin("groups: core,ida,cartesian,examples,radial,misc,angular_public,docs_fast", ci_workflow)
        @test !occursin(r"groups:[^\n]*[, ]angular(?:,|\n)", ci_workflow)
        @test contains_all(test_runner, ":docs_fast,",
            "_test_group_enabled(:docs_fast) || _test_group_enabled(:docs)",
            "docs\", \"public_surface_runtests.jl")
        @test length(findall("public_surface_runtests.jl", test_runner)) == 1
        @test occursin("const _FAST_TEST_GROUPS = Set((\n    :radial,\n    :core,\n    :ida,\n    :docs,\n))",
            test_runner)
        @test contains_all(ci_workflow, "tag-identity-install:",
            "refs/tags/\${tag}:refs/tags/\${tag}", "refs/tags/\${tag}^{tag}",
            "refs/tags/\${tag}^{commit}", "Pkg.add(url=", "GITHUB_REF_NAME")
        @test !occursin("workflow_run", ci_workflow)
        maintenance_paths = (
            "test/nested/cartesian_atomic_hf_reference_packet_runtests.jl", "test/nested/cartesian_occupied_first_injection_runtests.jl",
            "test/nested/cartesian_external_gto_import_runtests.jl", "test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl",
            "test/nested/cartesian_screened_hartree_correction_runtests.jl")
        @test issorted(first.(findfirst.(maintenance_paths, Ref(maintenance_workflow))))
        @test contains_all(maintenance_workflow, "workflow_dispatch:", "cron: \"17 10 * * 3\"",
            "timeout-minutes: 15", "contents: read")
        @test !occursin(r"(?m)^  (push|pull_request):", maintenance_workflow)
        @test length(findall("julia --project=. test/nested/", maintenance_workflow)) == 5
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
            "generated execution whitelist is stale or partially edited",
            "legacy authority artifact remains live",
            "--render",
        )
        @test startswith(execution_whitelist,
            "# Cartesian Hamiltonian Producer Execution Whitelist\n\n")
        @test contains_all(execution_whitelist, "Generated authority view. Do not edit.",
            "SHA-256", "HP-AUTHORITY-EXECUTION-WHITELIST-TEST-01")
        @test contains_all(agents, "authority.toml", "execution_whitelist.md",
            "grants nothing independently")
        @test all(needle -> !occursin(needle, agents),
            ("BEGIN CARTESIAN HAMILTONIAN PRODUCER EXECUTION WHITELIST",
             "Authority SHA-256:", "- `HP-"))
        @test all(needle -> !occursin(needle, authority_check),
            ("marked_whitelist_block", "_standalone_marker", "agents_whitelist_block.md"))

        @test contains_all(
            docs_site_reference_bases,
            "jldoctest",
            "UniformBasisSpec",
            "MappedUniformBasisSpec",
            "white_lindsey_atomic_mapping",
            "CombinedInvsqrtMapping",
        )
        @test all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_bases)\n"),
            _DOCUMENTED_PUBLIC_SURFACE.foundation)
        @test all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_bases)\n"),
            _DOCUMENTED_PUBLIC_SURFACE.function_stencil)
        @test occursin("## Partitions and leaf-local layers", docs_site_reference_bases) &&
              all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_bases)\n"),
            (_DOCUMENTED_PUBLIC_SURFACE.partition_leaf..., partition_leaf_context...))

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
        @test all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_atomic)\n"),
            (_DOCUMENTED_PUBLIC_SURFACE.atomic_ida..., atomic_ida_reference_context...))
        @test occursin(
            "## Atomic IDA inspection and tiny two-electron reference",
            docs_site_reference_atomic,
        )
        @test !Base.Docs.hasdoc(GaussletBases, :AbstractBondAlignedOrdinaryQWBasis3D)
        @test occursin("\nbasis_metadata\n", "\n$(docs_site_reference_bases)\n")
        @test all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_ops)\n"),
            (:BondAlignedDiatomicQWBasis3D, :CoulombGaussianExpansion))
        @test all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_export)\n"),
            _DOCUMENTED_PUBLIC_SURFACE.supported[4:6])

        @test occursin("## Expert/experimental bond-aligned geometry inspection", docs_site_reference_bases) &&
              all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_bases)\n"),
            _DOCUMENTED_PUBLIC_SURFACE.geometry)
        @test occursin("## Experimental Angular Profile And Sequence Producers",
            docs_site_reference_export,
        ) && all(
            name -> occursin("\n$(name)\n", "\n$(docs_site_reference_export)\n"),
            _DOCUMENTED_PUBLIC_SURFACE.angular,
        )
        @test !occursin("\nShellLocalAngularProfileKey\n", "\n$(docs_site_reference_export)\n")
        @test occursin("## Expert radial paper-parity prototype", docs_site_reference_bases) &&
              all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_bases)\n"),
            _DOCUMENTED_PUBLIC_SURFACE.radial_parity)
        @test contains_all(recommended_atomic_setup, "RadialBasisSpec", "0.09358986806", "0.02357750369", "0.0936", "0.0236")
        @test occursin("## Experimental QW geometry diagnostics", docs_site_reference_ops) &&
              all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_ops)\n"),
            _DOCUMENTED_PUBLIC_SURFACE.qw_geometry)
        @test occursin("## Experimental sliced hydrogen-chain operators", docs_site_reference_atomic) &&
              all(name -> occursin("\n$(name)\n", "\n$(docs_site_reference_atomic)\n"),
            _DOCUMENTED_PUBLIC_SURFACE.sliced_chain)
        @test occursin("\ncartesian_base_working_basis\n", "\n$(docs_site_reference_export)\n")
        @test occursin(
            "### Expert staged Cartesian construction",
            docs_site_reference_export,
        )
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
        @test startswith(atomic_direct_note, "> **Status:** supporting note.")
        @test startswith(ordinary_pgdg_note, "> **Status:** supporting development note.")
        @test startswith(global_map_note, "> **Note for new users:**")
        @test startswith(leaf_pgdg_note, "> **Note for new users:**")
    end
end

include(joinpath(@__DIR__, "cartesian_ham_builder_policy_runtests.jl"))
