const _PROTECTED_LOCALIZED_LADDER_BUNDLE_KIND = :protected_localized_ladder_bundle
const _PROTECTED_LOCALIZED_LADDER_BUNDLE_CONVENTION =
    :protected_localized_inherited_site_ladder_v1

_plb_clean(value) = replace(string(value), '\t' => ' ', '\n' => ' ')
_plb_sym(matrix) = Matrix{Float64}(0.5 .* (matrix .+ transpose(matrix)))

function _plb_write_tsv(path::AbstractString, rows, columns)
    open(path, "w") do io
        println(io, join(string.(columns), '\t'))
        for row in rows
            println(io, join((_plb_clean(getproperty(row, col)) for col in columns), '\t'))
        end
    end
    return path
end

function _plb_current_commit()
    try
        return readchomp(`git rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function _plb_source_commit(path::AbstractString)
    try
        return jldopen(path, "r") do file
            for key in (
                "producer_provenance/current_commit",
                "producer_provenance/git_commit",
                "producer_provenance/source_commit",
                "provenance/current_commit",
            )
                haskey(file, key) && return String(file[key])
            end
            return "unknown"
        end
    catch
        return "unknown"
    end
end

function _plb_stage!(rows, label::AbstractString, f)
    timed = @timed f()
    push!(rows, (; label, seconds = timed.time,
        alloc_mib = timed.bytes / 2.0^20, gc_s = timed.gctime))
    return timed.value
end

function _plb_recipe_from_artifact(path::AbstractString)
    jldopen(path, "r") do file
        expansion = _cartesian_read_coulomb_expansion_summary(file)
        locations = file["recipe_provenance/atom_locations"]
        recipe = (;
            source_artifact = String(path),
            source_commit = _plb_source_commit(path),
            atom_symbols = String.(file["recipe_provenance/atom_symbols"]),
            nuclear_charges = Float64.(file["recipe_provenance/nuclear_charges"]),
            atom_locations = [Tuple(Float64.(locations[i, :])) for i in axes(locations, 1)],
            nup = Int(file["recipe_provenance/nup"]),
            ndn = Int(file["recipe_provenance/ndn"]),
            ns = Int(file["recipe_provenance/ns"]),
            nesting = Symbol(file["recipe_provenance/nesting"]),
            core_spacing = Float64(file["recipe_provenance/core_spacing"]),
            s_factor = haskey(file, "recipe_provenance/s_factor") ?
                Float64(file["recipe_provenance/s_factor"]) : 1.0,
            xmax_parallel = Float64(file["recipe_provenance/xmax_parallel"]),
            xmax_transverse = Float64(file["recipe_provenance/xmax_transverse"]),
            basisname = String(file["recipe_provenance/basisname"]),
            basisfile = String(file["recipe_provenance/basisfile"]),
            lmax = Int(file["recipe_provenance/lmax"]),
            uncontracted = Bool(file["recipe_provenance/uncontracted"]),
            width_filtering = file["recipe_provenance/width_filtering"],
            residual_occupation_cutoff =
                Float64(file["supplement_provenance/occupation_cutoff"]),
            tau_neg_abs = Float64(file["supplement_provenance/tau_neg_abs"]),
            tau_neg_rel = Float64(file["supplement_provenance/tau_neg_rel"]),
            tau_merge_abs = Float64(file["supplement_provenance/tau_merge_abs"]),
            tau_merge_rel = Float64(file["supplement_provenance/tau_merge_rel"]))
        return isnothing(expansion) ? recipe :
            merge(recipe, (; coulomb_accuracy = expansion.policy))
    end
end

_plb_with_ns(recipe, ns::Integer) = merge(recipe, (; ns = Int(ns)))
_plb_get(recipe, key::Symbol, default) = hasproperty(recipe, key) ? getproperty(recipe, key) : default

function _plb_system_spec(recipe)
    return (; atom_symbols = recipe.atom_symbols,
        nuclear_charges = recipe.nuclear_charges,
        atom_locations = recipe.atom_locations,
        nup = recipe.nup,
        ndn = recipe.ndn)
end

function _plb_basis_spec(recipe)
    if hasproperty(recipe, :radius)
        return (; q = recipe.ns, core_spacing = recipe.core_spacing,
            s_factor = _plb_get(recipe, :s_factor, 1.0),
            coulomb_accuracy = _plb_get(recipe, :coulomb_accuracy, :compact),
            radius = recipe.radius,
            reference_spacing = _plb_get(recipe, :reference_spacing, 1.0),
            d = recipe.core_spacing)
    end
    return (; ns = recipe.ns, q = recipe.ns,
        nesting = _plb_get(recipe, :nesting, :pqs),
        core_spacing = recipe.core_spacing,
        s_factor = _plb_get(recipe, :s_factor, 1.0),
        coulomb_accuracy = _plb_get(recipe, :coulomb_accuracy, :compact),
        xmax_parallel = recipe.xmax_parallel,
        xmax_transverse = recipe.xmax_transverse,
        parent_axis_family = _plb_get(recipe, :parent_axis_family, :G10))
end

function _plb_supplement_spec(recipe)
    basis_by_center = hasproperty(recipe, :basis_by_center) ?
        recipe.basis_by_center :
        fill(recipe.basisname, length(recipe.atom_symbols))
    return (; basis_by_center,
        lmax = recipe.lmax,
        uncontracted = recipe.uncontracted,
        width_filtering = _plb_get(recipe, :width_filtering, nothing),
        basisfile = _plb_get(recipe, :basisfile, nothing))
end

function _plb_build_inputs(recipe, stages)
    base = _plb_stage!(stages, "ns$(recipe.ns) base working basis", () ->
        cartesian_base_working_basis(_plb_system_spec(recipe);
            basis = _plb_basis_spec(recipe), supplemented = true))
    supplement_basis = _plb_stage!(stages, "ns$(recipe.ns) supplement basis", () ->
        cartesian_residual_gto_supplement_basis(base, _plb_supplement_spec(recipe)))
    supplement = supplement_basis.basis
    X = _plb_stage!(stages, "ns$(recipe.ns) mixed terminal-supplement overlap", () ->
        CartesianFinalBasisRealization._terminal_residual_mixed_overlap(
            base.terminal_basis, base.parent.parent_axis_bundle_object, supplement))
    S_AA = _plb_stage!(stages, "ns$(recipe.ns) supplement overlap", () ->
        Matrix{Float64}(_cartesian_supplement_cross_overlap(supplement, supplement)))
    products = _plb_stage!(stages, "ns$(recipe.ns) base products", () ->
        cartesian_base_products(base))
    unit_nuclear = _plb_stage!(stages, "ns$(recipe.ns) base unit nuclear", () ->
        cartesian_base_unit_nuclear(base))
    vee = _plb_stage!(stages, "ns$(recipe.ns) base Vee", () ->
        cartesian_base_vee(base))
    ham = _plb_stage!(stages, "ns$(recipe.ns) base Hamiltonian", () ->
        cartesian_base_hamiltonian_assembly(base, products, unit_nuclear, vee))
    labels = CartesianResidualGaussians.residual_gaussian_candidate_labels(supplement)
    centers = CartesianResidualGaussians.residual_gaussian_candidate_centers(supplement)
    nuclei = CartesianResidualGaussians.residual_gaussian_float_centers(base.input.locations)
    owners = Int[CartesianResidualGaussians.residual_candidate_owner(center, nuclei)
        for center in centers]
    return (; recipe, base, supplement_basis, supplement, X, S_AA, ham,
        labels, centers, owners)
end

function _plb_compactness(inputs)
    return (; metric = :midpoint_weighted_tail, cutoff = 0.2,
        supplement = inputs.supplement, selector = :ordered_compact_first_mgs)
end

function _plb_compact_residual(inputs, stages)
    recipe = inputs.recipe
    return _plb_stage!(stages, "ns$(recipe.ns) compact-first residual basis", () ->
        CartesianResidualGaussians.build_residual_gaussian_basis(
            inputs.base.terminal_basis.final_dimension,
            inputs.X, inputs.S_AA, inputs.labels, inputs.centers, inputs.owners;
            residual_occupation_cutoff =
                _plb_get(recipe, :residual_occupation_cutoff, 1.0e-6),
            tau_neg_abs = _plb_get(recipe, :tau_neg_abs, 1.0e-12),
            tau_neg_rel = _plb_get(recipe, :tau_neg_rel, 1.0e-12),
            tau_merge_abs = _plb_get(recipe, :tau_merge_abs, 1.0e-12),
            tau_merge_rel = _plb_get(recipe, :tau_merge_rel, 1.0e-12),
            residual_compactness = _plb_compactness(inputs)))
end

function _plb_protected_geometry(inputs, residual, stages; occupied_blocks = Matrix{Float64}[])
    recipe = inputs.recipe
    return _plb_stage!(stages, "ns$(recipe.ns) staged protected-original geometry", () ->
        CartesianResidualGaussians.staged_protected_original_injection_geometry(
            inputs.base.terminal_basis.final_dimension,
            inputs.X, inputs.S_AA, inputs.labels, inputs.centers, inputs.owners,
            residual; occupied_blocks,
            s_cut = _plb_get(recipe, :protected_s_cut, 0.95),
            occ_cut = _plb_get(recipe, :protected_occ_cut, 0.003),
            candidate_overlap_atol = 1.0e-12,
            candidate_overlap_rtol = 1.0e-8))
end

function _plb_augmented(inputs, residual, stages)
    recipe = inputs.recipe
    products = _plb_stage!(stages, "ns$(recipe.ns) augmented products", () ->
        cartesian_residual_gto_augmented_products(
            inputs.base, inputs.supplement_basis, residual;
            base_kinetic = inputs.ham.kinetic))
    unit_nuclear = _plb_stage!(stages, "ns$(recipe.ns) augmented unit nuclear", () ->
        cartesian_residual_gto_augmented_unit_nuclear(
            inputs.base, residual, products;
            base_unit_nuclear = inputs.ham.nuclear_attraction_unit_by_center))
    vee = _plb_stage!(stages, "ns$(recipe.ns) augmented Vee", () ->
        cartesian_residual_gto_augmented_vee(
            inputs.base, inputs.ham, residual, products, unit_nuclear))
    return (; products, unit_nuclear, vee)
end

function _plb_protected_one_body(inputs, geometry, augmented)
    blocks = augmented.products.supplement_blocks
    raw_K = (; GG = inputs.ham.kinetic, GA = blocks.mixed.kinetic,
        AA = blocks.self.kinetic)
    raw_U = [(; GG, GA, AA) for (GG, GA, AA) in zip(
        inputs.ham.nuclear_attraction_unit_by_center,
        blocks.mixed.nuclear, blocks.self.nuclear)]
    return CartesianResidualGaussians.transform_protected_original_fixed_sector_one_body(
        raw_K, raw_U, inputs.base.input.charges, geometry).one_body_hamiltonian
end

function _plb_localized_raw_coefficients(geometry, loc)
    components = CartesianResidualGaussians.protected_original_fixed_sector_components(geometry)
    nG = size(geometry.T_G, 1)
    nZ = size(components.Z, 2)
    G_F = hcat(zeros(Float64, nG, nZ), components.G_perp)
    A_F = hcat(components.Z, components.A_perp)
    return (; G_L = G_F * loc.W, A_L = A_F * loc.W)
end

function _plb_member_diagnostics(case)
    geometry = case.geometry
    localized = case.localized.diagnostics
    return (;
        B_min = geometry.b_min,
        B_median = geometry.b_median,
        B_max = geometry.b_max,
        B_lt_0p999 = geometry.b_lt_0p999,
        B_lt_0p99 = geometry.b_lt_0p99,
        B_lt_0p98 = geometry.b_lt_0p98,
        B_lt_0p95 = geometry.b_lt_0p95,
        B_lt_0p9 = geometry.b_lt_0p9,
        F_S_F_identity_error = geometry.f_s_f_identity_block_max,
        Z_S_M_Qperp_error = geometry.z_m_qperp_max,
        Qperp_identity_error = geometry.qperp_identity_sample_max,
        protected_span_min_sv = geometry.protected_span_min_sv,
        L_identity_error = localized.L_identity_error,
        M_L_diag_delta_max = localized.M_L_diag_delta_max,
        M_L_offdiag_max = localized.M_L_offdiag_max,
        M_L_fro_delta = localized.M_L_fro_delta,
        H1_L_symmetry_error = localized.H1_L_symmetry_error,
        Vee_L_symmetry_error = localized.Vee_L_symmetry_error)
end

function _plb_finish_member(recipe, inputs, residual, geometry, stages)
    augmented = _plb_augmented(inputs, residual, stages)
    H_F = _plb_stage!(stages, "ns$(recipe.ns) protected fixed one-body", () ->
        _plb_protected_one_body(inputs, geometry, augmented))
    localized = _plb_stage!(stages, "ns$(recipe.ns) protected-localized H1/Vee", () ->
        CartesianResidualGaussians.protected_localized_inherited_site_hamiltonian(
            H_F, augmented.vee, geometry))
    row_locality = _plb_stage!(stages, "ns$(recipe.ns) row locality", () ->
        CartesianResidualGaussians.protected_localized_row_locality(
            localized.transform, augmented.products.position;
            sector_counts = (; base = inputs.base.terminal_basis.final_dimension,
                compact_R = residual.residual_dimension),
            x2 = augmented.products.x2))
    raw = _plb_stage!(stages, "ns$(recipe.ns) localized raw coefficients", () ->
        _plb_localized_raw_coefficients(geometry, localized.transform))
    return (; recipe, inputs, residual, geometry, augmented, localized,
        row_locality, raw, H = localized.H1_L, V = localized.Vee_L)
end

function _plb_build_member(recipe, stages)
    inputs = _plb_build_inputs(recipe, stages)
    residual = _plb_compact_residual(inputs, stages)
    geometry = _plb_protected_geometry(inputs, residual, stages)
    return _plb_finish_member(recipe, inputs, residual, geometry, stages)
end
_plb_reference_embeddings(inputs, placements, stages) = [_plb_stage!(stages,
        "ns$(inputs.recipe.ns) embed reference owner $(placement.owner_index)", () ->
        CartesianReferenceDensity.atomic_reference_packet_occupied_embedding(
            placement.packet, inputs.supplement, inputs.S_AA, inputs.owners;
            owner_index = placement.owner_index, center = placement.center,
            supplement_indices = placement.supplement_indices)
    ) for placement in placements]
function _plb_additive_reference_raw_blocks(inputs, embeddings, stages)
    raw_blocks = [_plb_stage!(stages,
            "ns$(inputs.recipe.ns) fitted reference field owner $(embedding.owner_index)", () ->
                CartesianReferenceDensity.atomic_reference_packet_fitted_potential_raw_blocks(
                    inputs.base, inputs.supplement, embedding.packet;
                    center = embedding.center)) for embedding in embeddings]
    GG, GA, AA = copy(first(raw_blocks).GG), copy(first(raw_blocks).GA),
        copy(first(raw_blocks).AA)
    for raw in raw_blocks[2:end]
        GG .+= raw.GG
        GA .+= raw.GA
        AA .+= raw.AA
    end
    return (; GG, GA, AA)
end
function _plb_build_additive_reference_member(
    recipe, stages, placements; correction_options = (;))
    isempty(placements) && return (; member = _plb_build_member(recipe, stages),
        correction = nothing, reference = nothing)
    inputs = _plb_build_inputs(recipe, stages)
    residual = _plb_compact_residual(inputs, stages)
    embeddings = _plb_reference_embeddings(inputs, placements, stages)
    occupied_blocks = [entry.Y for entry in embeddings]
    geometry = _plb_protected_geometry(inputs, residual, stages;
        occupied_blocks)
    member = _plb_finish_member(recipe, inputs, residual, geometry, stages)
    represented = _plb_stage!(stages, "ns$(recipe.ns) represented reference in L", () ->
        CartesianResidualGaussians.protected_original_reference_blocks_in_local_basis(
            member.raw.G_L, member.raw.A_L, inputs.X, inputs.S_AA,
            member.localized.transform,
            occupied_blocks))
    raw = _plb_additive_reference_raw_blocks(inputs, embeddings, stages)
    hartree = _plb_stage!(stages, "ns$(recipe.ns) additive reference J0 in F/L", () ->
        CartesianResidualGaussians.transform_protected_original_localized_exact_hartree(
            raw, geometry, member.localized.transform))
    energy = _plb_stage!(stages, "ns$(recipe.ns) additive reference density energy", () ->
        CartesianReferenceDensity.atomic_reference_packet_additive_density_energy(embeddings))
    correction = _plb_stage!(stages, "ns$(recipe.ns) additive screened-Hartree correction", () ->
        CartesianReferenceDensity.build_additive_screened_hartree_correction(
            member.V, hartree.J0_L, energy.total, represented.coefficient_blocks,
            [entry.occupations for entry in embeddings];
            packets = [entry.packet for entry in embeddings], correction_options...))
    reference = (; represented = represented.diagnostics,
        hartree = hartree.diagnostics, energy,
        geometry = geometry.additive_reference)
    return (; member, correction, reference)
end

function _plb_write_member_artifact(case, path::AbstractString, source_artifact, source_commit, current_commit)
    mkpath(dirname(path))
    recipe = case.recipe
    write_protected_localized_ida_hamiltonian(path;
        H1_L = case.H,
        Vee_L = case.V,
        nup = recipe.nup,
        ndn = recipe.ndn,
        nuclear_charges = recipe.nuclear_charges,
        nuclear_positions = recipe.atom_locations,
        sector_counts = (; base = case.inputs.base.terminal_basis.final_dimension,
            compact_R = case.residual.residual_dimension,
            protected_Z = case.geometry.protected_original_count,
            broad_Z = case.geometry.z_dimension - case.geometry.protected_original_count),
        diagnostics = _plb_member_diagnostics(case),
        provenance = (; source_artifact = String(source_artifact),
            source_commit = String(source_commit), current_commit = String(current_commit)),
        basis_controls = (; nesting = _plb_get(recipe, :nesting, :pqs),
            ns = recipe.ns, core_spacing = recipe.core_spacing,
            s_factor = _plb_get(recipe, :s_factor, 1.0),
            basisname = recipe.basisname, lmax = recipe.lmax),
        geometry_inputs = (; atom_symbols = recipe.atom_symbols,
            nuclear_charges = recipe.nuclear_charges,
            atom_locations = recipe.atom_locations),
        coulomb_expansion = _cartesian_coulomb_expansion_summary(
            case.inputs.base.input.coulomb_accuracy,
            case.inputs.base.coulomb_expansion),
        localized_ordering = (; matrix_order = :native,
            convention = _PROTECTED_LOCALIZED_LADDER_BUNDLE_CONVENTION),
        row_locality = case.row_locality)
    return path
end

function _plb_axis_centers(case, axis)
    return Float64.(_nested_axis_pgdg(
        case.inputs.base.parent.parent_axis_bundle_object, axis).centers)
end

function _plb_parent_lattice_rows(source, target; atol::Real = 1.0e-12)
    rows = NamedTuple[]
    ok = true
    for axis in (:x, :y, :z)
        source_centers = _plb_axis_centers(source, axis)
        target_centers = _plb_axis_centers(target, axis)
        same_count = length(source_centers) == length(target_centers)
        delta = same_count ? norm(source_centers - target_centers, Inf) : Inf
        axis_ok = same_count && delta <= atol
        ok &= axis_ok
        push!(rows, (; axis, source_count = length(source_centers),
            target_count = length(target_centers), same_count,
            center_max_abs_delta = delta, ok = axis_ok))
    end
    return ok, rows
end

function _plb_same_supplement(source, target; atol::Real = 0.0)
    source.inputs.labels == target.inputs.labels || return false
    length(source.inputs.centers) == length(target.inputs.centers) || return false
    for (a, b) in zip(source.inputs.centers, target.inputs.centers)
        maximum(abs(Float64(a[i]) - Float64(b[i])) for i in 1:3) <= atol || return false
    end
    return true
end

function _plb_support_overlap_matrix(left_block, right_block, overlaps)
    out = Matrix{Float64}(undef, length(left_block.support_states),
        length(right_block.support_states))
    CartesianFinalBasisRealization._support_cross!(
        out, left_block.support_states, right_block.support_states, overlaps)
    return out
end

function _plb_terminal_cross_action(left_basis, right_basis, right_coeffs, bundles)
    overlaps = (_nested_axis_pgdg(bundles, :x).overlap,
        _nested_axis_pgdg(bundles, :y).overlap,
        _nested_axis_pgdg(bundles, :z).overlap)
    result = zeros(Float64, left_basis.final_dimension, size(right_coeffs, 2))
    for left_block in left_basis.blocks
        acc = zeros(Float64, length(left_block.support_states), size(right_coeffs, 2))
        for right_block in right_basis.blocks
            cross = _plb_support_overlap_matrix(left_block, right_block, overlaps)
            right_slice = @view right_coeffs[right_block.column_range, :]
            if isnothing(right_block.coefficients)
                mul!(acc, cross, right_slice, 1.0, 1.0)
            else
                mul!(acc, cross, right_block.coefficients * right_slice, 1.0, 1.0)
            end
        end
        if isnothing(left_block.coefficients)
            result[left_block.column_range, :] .+= acc
        else
            result[left_block.column_range, :] .+= transpose(left_block.coefficients) * acc
        end
    end
    return result
end

function _plb_cross_overlap(target, source, stages)
    ok, _ = _plb_parent_lattice_rows(source, target)
    ok || throw(ArgumentError("protected ladder transfer requires shared parent lattice"))
    _plb_same_supplement(source, target) ||
        throw(ArgumentError("protected ladder transfer requires identical supplement candidates"))
    Sgg_GL = _plb_stage!(stages,
        "terminal cross action ns$(target.recipe.ns)<-ns$(source.recipe.ns)", () ->
            _plb_terminal_cross_action(
                target.inputs.base.terminal_basis,
                source.inputs.base.terminal_basis,
                source.raw.G_L,
                target.inputs.base.parent.parent_axis_bundle_object))
    X_target_source = _plb_stage!(stages, "mixed overlap G_target A_source", () ->
        CartesianFinalBasisRealization._terminal_residual_mixed_overlap(
            target.inputs.base.terminal_basis,
            target.inputs.base.parent.parent_axis_bundle_object,
            source.inputs.supplement))
    X_source_target = _plb_stage!(stages, "mixed overlap G_source A_target", () ->
        CartesianFinalBasisRealization._terminal_residual_mixed_overlap(
            source.inputs.base.terminal_basis,
            source.inputs.base.parent.parent_axis_bundle_object,
            target.inputs.supplement))
    S_AA = _plb_stage!(stages, "supplement cross overlap A_target A_source", () ->
        Matrix{Float64}(_cartesian_supplement_cross_overlap(
            target.inputs.supplement, source.inputs.supplement)))
    return _plb_stage!(stages, "assemble protected-localized cross overlap", () ->
        transpose(target.raw.G_L) * Sgg_GL +
        transpose(target.raw.G_L) * X_target_source * source.raw.A_L +
        transpose(target.raw.A_L) * transpose(X_source_target) * source.raw.G_L +
        transpose(target.raw.A_L) * S_AA * source.raw.A_L)
end

function _plb_singular_summary(S)
    values = sqrt.(max.(0.0, eigvals(Symmetric(_plb_sym(transpose(S) * S)))))
    sort!(values)
    return (; singular_min = first(values),
        singular_median = values[cld(length(values), 2)],
        singular_max = last(values),
        sv_lt_0p999 = count(<(0.999), values),
        sv_lt_0p99 = count(<(0.99), values),
        sv_lt_0p95 = count(<(0.95), values),
        sv_lt_0p9 = count(<(0.9), values))
end

function _plb_read_orbitals(path::AbstractString, row_locality, dim::Int)
    jldopen(path, "r") do file
        up = Matrix{Float64}(file["psiup"])
        dn = Matrix{Float64}(file["psidn"])
        size(up, 1) == dim && size(dn, 1) == dim ||
            throw(DimensionMismatch("saved orbital dimension mismatch"))
        order = haskey(file, "site_order_kind") ? Symbol(String(file["site_order_kind"])) :
            :native
        if order in (:z_order, :z, :bond_axis_z)
            perm = Int.(row_locality.native_to_z_order)
            up = up[perm, :]
            dn = dn[perm, :]
        elseif !(order in (:native, :inherited_M_site_order))
            throw(ArgumentError("unknown saved orbital site_order_kind $(order)"))
        end
        return (; up, dn, order, nup = size(up, 2), ndn = size(dn, 2))
    end
end

function _plb_orbital_path(saved_orbitals, ns)
    isnothing(saved_orbitals) && return nothing
    saved_orbitals isa AbstractString && return String(saved_orbitals)
    if saved_orbitals isa AbstractDict
        haskey(saved_orbitals, ns) && return String(saved_orbitals[ns])
        haskey(saved_orbitals, string(ns)) && return String(saved_orbitals[string(ns)])
    end
    hasproperty(saved_orbitals, Symbol("ns$(ns)")) &&
        return String(getproperty(saved_orbitals, Symbol("ns$(ns)")))
    return nothing
end

function _plb_transfer_diagnostics(source_orb, up_target, dn_target)
    up_overlap = transpose(up_target) * up_target
    dn_overlap = transpose(dn_target) * dn_target
    return (; source_alpha_trace = tr(transpose(source_orb.up) * source_orb.up),
        source_beta_trace = tr(transpose(source_orb.dn) * source_orb.dn),
        target_alpha_trace = tr(up_overlap),
        target_beta_trace = tr(dn_overlap),
        alpha_trace_loss = source_orb.nup - tr(up_overlap),
        beta_trace_loss = source_orb.ndn - tr(dn_overlap),
        alpha_orth_inf = norm(up_overlap - I, Inf),
        beta_orth_inf = norm(dn_overlap - I, Inf))
end

function _plb_density(C)
    return Matrix{Float64}(C * transpose(C))
end

function _plb_trace_product(A, B)
    value = 0.0
    @inbounds for j in axes(A, 2), i in axes(A, 1)
        value += A[i, j] * B[i, j]
    end
    return value
end

function _plb_exchange_energy(V, rho)
    value = 0.0
    @inbounds for j in axes(V, 2), i in axes(V, 1)
        value += V[i, j] * rho[i, j] * rho[i, j]
    end
    return -0.5 * value
end

function _plb_density_energy(case, up, dn)
    rho_up = _plb_density(up)
    rho_dn = _plb_density(dn)
    rho = rho_up + rho_dn
    q = diag(rho)
    one = _plb_trace_product(case.H, rho)
    hartree = 0.5 * dot(q, case.V * q)
    ex_up = _plb_exchange_energy(case.V, rho_up)
    ex_dn = _plb_exchange_energy(case.V, rho_dn)
    electronic = one + hartree + ex_up + ex_dn
    nuclear = _cartesian_ida_nuclear_repulsion(
        Float64.(case.recipe.nuclear_charges),
        _cartesian_nuclear_position_matrix(case.recipe.atom_locations))
    return (; one_body = one, hartree, exchange_alpha = ex_up,
        exchange_beta = ex_dn, interaction = hartree + ex_up + ex_dn,
        electronic, nuclear, total = electronic + nuclear)
end

function _plb_member_row(case, artifact_path::AbstractString)
    geometry = case.geometry
    return (; ns = case.recipe.ns,
        artifact_path,
        base_dimension = case.inputs.base.terminal_basis.final_dimension,
        compact_R = case.residual.residual_dimension,
        final_dimension = size(case.H, 1),
        protected_originals = geometry.protected_original_count,
        broad_Z = geometry.z_dimension - geometry.protected_original_count,
        Z_dimension = geometry.z_dimension,
        B_min = geometry.b_min,
        B_median = geometry.b_median,
        B_max = geometry.b_max,
        B_lt_0p999 = geometry.b_lt_0p999,
        B_lt_0p99 = geometry.b_lt_0p99,
        H1_symmetry = norm(case.H - transpose(case.H), Inf),
        Vee_symmetry = norm(case.V - transpose(case.V), Inf),
        H1_finite = all(isfinite, case.H),
        Vee_finite = all(isfinite, case.V))
end

function _plb_write_manifest(path, rows; coulomb_expansion)
    jldopen(path, "w") do file
        file["artifact_kind"] = _PROTECTED_LOCALIZED_LADDER_BUNDLE_KIND
        file["format_version"] = 1
        file["convention_id"] = _PROTECTED_LOCALIZED_LADDER_BUNDLE_CONVENTION
        file["convention_version"] = 1
        for (name, value) in pairs(rows)
            file[String(name)] = value
        end
        _cartesian_write_coulomb_expansion_summary!(file, coulomb_expansion)
    end
    return path
end

function _plb_write_restart(path, source_ns, target_ns, source_path, up, dn, diagnostics, energy)
    mkpath(dirname(path))
    jldopen(path, "w") do file
        file["artifact_kind"] = :protected_localized_ladder_transferred_orbitals
        file["format_version"] = 1
        file["convention_id"] = _PROTECTED_LOCALIZED_LADDER_BUNDLE_CONVENTION
        file["source_ns"] = Int(source_ns)
        file["target_ns"] = Int(target_ns)
        file["source_orbitals"] = String(source_path)
        file["site_order_kind"] = "native"
        file["psiup"] = Matrix{Float64}(up)
        file["psidn"] = Matrix{Float64}(dn)
        for (key, value) in pairs(diagnostics)
            file["diagnostics/$(key)"] = value
        end
        for (key, value) in pairs(energy)
            file["fixed_density_energy/$(key)"] = value
        end
    end
    return path
end

function build_protected_localized_ladder_bundle(
    recipe_or_artifact;
    ns_values,
    output_dir,
    saved_orbitals = nothing,
)
    recipe0 = recipe_or_artifact isa AbstractString ?
        _plb_recipe_from_artifact(String(recipe_or_artifact)) :
        recipe_or_artifact
    source_artifact = _plb_get(recipe0, :source_artifact, "")
    source_commit = _plb_get(recipe0, :source_commit, "unknown")
    current_commit = _plb_current_commit()
    ns_list = Int.(collect(ns_values))
    length(ns_list) >= 1 || throw(ArgumentError("ns_values must not be empty"))
    mkpath(output_dir)
    mkpath(joinpath(output_dir, "members"))
    mkpath(joinpath(output_dir, "transfers"))
    mkpath(joinpath(output_dir, "restarts"))
    mkpath(joinpath(output_dir, "summaries"))

    stages = NamedTuple[]
    cases = Dict{Int,Any}()
    member_rows = NamedTuple[]
    member_paths = Dict{Int,String}()
    coulomb_expansion = nothing
    for ns in ns_list
        recipe = _plb_with_ns(recipe0, ns)
        case = _plb_stage!(stages, "build ns$(ns) protected-localized member", () ->
            _plb_build_member(recipe, stages))
        cases[ns] = case
        member_expansion = _cartesian_coulomb_expansion_summary(
            case.inputs.base.input.coulomb_accuracy,
            case.inputs.base.coulomb_expansion)
        if isnothing(coulomb_expansion)
            coulomb_expansion = member_expansion
        else
            member_expansion == coulomb_expansion || throw(ArgumentError(
                "protected ladder members use different Coulomb expansions"))
        end
        member_path = joinpath(output_dir, "members", "ns$(ns)",
            "protected_localized_hamiltonian.jld2")
        _plb_write_member_artifact(case, member_path, source_artifact,
            source_commit, current_commit)
        member_paths[ns] = member_path
        push!(member_rows, _plb_member_row(case, member_path))
    end

    parent_rows = NamedTuple[]
    transfer_rows = NamedTuple[]
    restart_rows = NamedTuple[]
    transfer_paths = String[]
    restart_paths = String[]
    for pair in 2:length(ns_list)
        source_ns = ns_list[pair - 1]
        target_ns = ns_list[pair]
        source = cases[source_ns]
        target = cases[target_ns]
        ok, rows = _plb_parent_lattice_rows(source, target)
        append!(parent_rows, [merge(row, (; source_ns, target_ns)) for row in rows])
        ok || throw(ArgumentError("parent lattice mismatch for ns$(source_ns)->ns$(target_ns)"))
        same_supplement = _plb_same_supplement(source, target)
        same_supplement ||
            throw(ArgumentError("supplement mismatch for ns$(source_ns)->ns$(target_ns)"))
        S = _plb_cross_overlap(target, source, stages)
        transfer_path = joinpath(output_dir, "transfers",
            "S_ns$(target_ns)_ns$(source_ns).jld2")
        jldsave(transfer_path;
            artifact_kind = :protected_localized_ladder_cross_overlap,
            format_version = 1,
            convention_id = _PROTECTED_LOCALIZED_LADDER_BUNDLE_CONVENTION,
            source_ns,
            target_ns,
            S_BA = S)
        push!(transfer_paths, transfer_path)
        singular = _plb_singular_summary(S)
        row = merge((; source_ns, target_ns, transfer_path,
            cross_rows = size(S, 1), cross_cols = size(S, 2), same_supplement),
            singular)

        orbital_path = _plb_orbital_path(saved_orbitals, source_ns)
        if !isnothing(orbital_path)
            source_orb = _plb_read_orbitals(orbital_path, source.row_locality, size(S, 2))
            up_target = S * source_orb.up
            dn_target = S * source_orb.dn
            transfer_diag = _plb_transfer_diagnostics(source_orb, up_target, dn_target)
            energy = _plb_density_energy(target, up_target, dn_target)
            restart_path = joinpath(output_dir, "restarts",
                "ns$(target_ns)_from_ns$(source_ns)_occupied_orbitals.jld2")
            _plb_write_restart(restart_path, source_ns, target_ns, orbital_path,
                up_target, dn_target, transfer_diag, energy)
            push!(restart_paths, restart_path)
            push!(restart_rows, merge((; source_ns, target_ns, restart_path,
                source_orbitals = orbital_path), transfer_diag, energy))
            row = merge(row, transfer_diag, (; fixed_density_total = energy.total,
                restart_path))
        end
        push!(transfer_rows, row)
    end

    _plb_write_tsv(joinpath(output_dir, "summaries", "ladder_members.tsv"),
        member_rows, [:ns, :artifact_path, :base_dimension, :compact_R,
            :final_dimension, :protected_originals, :broad_Z, :Z_dimension,
            :B_min, :B_median, :B_max, :B_lt_0p999, :B_lt_0p99,
            :H1_symmetry, :Vee_symmetry, :H1_finite, :Vee_finite])
    _plb_write_tsv(joinpath(output_dir, "summaries", "parent_lattice.tsv"),
        parent_rows, [:source_ns, :target_ns, :axis, :source_count,
            :target_count, :same_count, :center_max_abs_delta, :ok])
    _plb_write_tsv(joinpath(output_dir, "summaries", "transfers.tsv"),
        transfer_rows, [:source_ns, :target_ns, :transfer_path, :cross_rows,
            :cross_cols, :same_supplement, :singular_min, :singular_median,
            :singular_max, :sv_lt_0p999, :sv_lt_0p99, :sv_lt_0p95, :sv_lt_0p9])
    if !isempty(restart_rows)
        _plb_write_tsv(joinpath(output_dir, "summaries", "restarts.tsv"),
            restart_rows, [:source_ns, :target_ns, :restart_path, :source_orbitals,
                :source_alpha_trace, :source_beta_trace, :target_alpha_trace,
                :target_beta_trace, :alpha_trace_loss, :beta_trace_loss,
                :alpha_orth_inf, :beta_orth_inf, :total])
    end
    _plb_write_tsv(joinpath(output_dir, "summaries", "stages.tsv"),
        stages, [:label, :seconds, :alloc_mib, :gc_s])

    manifest_path = joinpath(output_dir, "manifest.jld2")
    _plb_write_manifest(manifest_path, (;
        git_head = current_commit,
        source_artifact = String(source_artifact),
        source_commit = String(source_commit),
        ns_values = ns_list,
        member_ns = Int.(keys(member_paths) |> collect |> sort),
        member_paths = [member_paths[ns] for ns in sort(collect(keys(member_paths)))],
        transfer_paths,
        restart_paths,
        atom_symbols = String.(recipe0.atom_symbols),
        nuclear_charges = Float64.(recipe0.nuclear_charges),
        atom_locations = reduce(vcat, [reshape(collect(Float64.(loc)), 1, :)
            for loc in recipe0.atom_locations]),
        nup = Int(recipe0.nup),
        ndn = Int(recipe0.ndn),
        core_spacing = Float64(recipe0.core_spacing),
        s_factor = Float64(_plb_get(recipe0, :s_factor, 1.0)),
        basisname = String(recipe0.basisname),
        lmax = Int(recipe0.lmax),
        parent_lattice_ok = all(row -> row.ok, parent_rows),
        summary_dir = joinpath(output_dir, "summaries"));
        coulomb_expansion)
    return read_protected_localized_ladder_bundle(output_dir)
end

function read_protected_localized_ladder_bundle(path; load_members::Bool = false,
    load_transfers::Bool = false, load_restarts::Bool = false)
    root = isdir(path) ? String(path) : dirname(String(path))
    manifest_path = isdir(path) ? joinpath(root, "manifest.jld2") : String(path)
    manifest = jldopen(manifest_path, "r") do file
        Symbol(file["artifact_kind"]) == _PROTECTED_LOCALIZED_LADDER_BUNDLE_KIND ||
            throw(ArgumentError("not a protected-localized ladder bundle manifest"))
        Int(file["format_version"]) == 1 ||
            throw(ArgumentError("unsupported protected ladder bundle format"))
        (;
            artifact_kind = Symbol(file["artifact_kind"]),
            format_version = Int(file["format_version"]),
            convention_id = Symbol(file["convention_id"]),
            convention_version = Int(file["convention_version"]),
            git_head = String(file["git_head"]),
            source_artifact = String(file["source_artifact"]),
            source_commit = String(file["source_commit"]),
            ns_values = Int.(file["ns_values"]),
            member_ns = Int.(file["member_ns"]),
            member_paths = String.(file["member_paths"]),
            transfer_paths = String.(file["transfer_paths"]),
            restart_paths = String.(file["restart_paths"]),
            coulomb_expansion = _cartesian_read_coulomb_expansion_summary(file),
            parent_lattice_ok = Bool(file["parent_lattice_ok"]),
            summary_dir = String(file["summary_dir"]))
    end
    members = load_members ?
        [read_protected_localized_ida_hamiltonian(path) for path in manifest.member_paths] :
        nothing
    transfers = load_transfers ?
        [jldopen(path, "r") do file
            (; source_ns = Int(file["source_ns"]), target_ns = Int(file["target_ns"]),
                S_BA = Matrix{Float64}(file["S_BA"]))
        end for path in manifest.transfer_paths] :
        nothing
    restarts = load_restarts ?
        [jldopen(path, "r") do file
            (; source_ns = Int(file["source_ns"]), target_ns = Int(file["target_ns"]),
                psiup = Matrix{Float64}(file["psiup"]),
                psidn = Matrix{Float64}(file["psidn"]))
        end for path in manifest.restart_paths] :
        nothing
    return (; root, manifest_path, manifest, members, transfers, restarts)
end
