# Private RHF input contract helpers for complete core/shell PQS diagnostics.

const _PQS_COMPLETE_CORE_SHELL_RHF_FIXTURE_ROLES =
    (:route_smoke, :physics_endpoint)

function _pqs_multilayer_rhf_contract_property(object, key::Symbol, default = nothing)
    isnothing(object) && return default
    hasproperty(object, key) && return getproperty(object, key)
    return default
end

function _pqs_multilayer_rhf_contract_summary_property(
    object,
    key::Symbol,
    default = nothing,
)
    summary = _pqs_multilayer_rhf_contract_property(object, :summary)
    isnothing(summary) && return default
    return _pqs_multilayer_rhf_contract_property(summary, key, default)
end

function _pqs_multilayer_complete_core_shell_rhf_occupation(electron_count)
    if isnothing(electron_count)
        return nothing, :missing_electron_count
    end
    if !(electron_count isa Integer) || electron_count isa Bool || electron_count <= 0
        return nothing, :invalid_electron_count
    end
    iseven(electron_count) || return nothing, :unsupported_open_shell_rhf_input
    nocc = electron_count ÷ 2
    return (;
        electron_count,
        electron_count_source = :explicit,
        spin_policy = :closed_shell_rhf,
        nocc,
        occupancy = 2,
        occupation_vector = ntuple(_ -> 2, nocc),
        fractional_occupation_supported = false,
        open_shell_supported = false,
    ), nothing
end

function _pqs_multilayer_complete_core_shell_rhf_missing_inputs(;
    source_plan,
    final_basis,
    h1_payload,
    density_inputs,
    coulomb_expansion,
    electron_count,
    fixture_role,
)
    missing = Symbol[]
    source_plan_status =
        _pqs_multilayer_rhf_contract_property(source_plan, :status)
    source_plan_kind =
        _pqs_multilayer_rhf_contract_property(source_plan, :object_kind)
    if source_plan_kind !== :pqs_multilayer_shell_source_plan ||
       source_plan_status !== :available_pqs_multilayer_shell_source_plan
        push!(missing, :source_plan)
    end

    final_basis_status =
        _pqs_multilayer_rhf_contract_property(final_basis, :status)
    final_basis_kind =
        _pqs_multilayer_rhf_contract_property(final_basis, :object_kind)
    if final_basis_kind !== :pqs_complete_core_shell_final_basis ||
       final_basis_status !== :available_pqs_complete_core_shell_final_basis
        push!(missing, :final_basis)
    end

    h1_payload_status =
        _pqs_multilayer_rhf_contract_property(h1_payload, :status)
    h1_payload_kind =
        _pqs_multilayer_rhf_contract_property(h1_payload, :object_kind)
    final_hamiltonian =
        _pqs_multilayer_rhf_contract_property(h1_payload, :final_hamiltonian)
    final_h1_matrix =
        _pqs_multilayer_rhf_contract_property(
            final_hamiltonian,
            :hamiltonian_matrix,
        )
    if h1_payload_kind !== :pqs_multilayer_complete_core_shell_h1_payload ||
       h1_payload_status !==
       :materialized_pqs_multilayer_complete_core_shell_h1_payload ||
       isnothing(final_h1_matrix)
        push!(missing, :h1_payload)
    end

    density_inputs_status =
        _pqs_multilayer_rhf_contract_property(density_inputs, :status)
    axis_weights =
        _pqs_multilayer_rhf_contract_property(density_inputs, :axis_weights)
    raw_pair_factor_terms =
        _pqs_multilayer_rhf_contract_property(
            density_inputs,
            :raw_pair_factor_terms,
        )
    if density_inputs_status !== :available_complete_core_shell_density_inputs ||
       isnothing(axis_weights) ||
       isnothing(raw_pair_factor_terms)
        push!(missing, :density_inputs)
    end

    isnothing(coulomb_expansion) && push!(missing, :coulomb_expansion)
    isnothing(electron_count) && push!(missing, :electron_count)
    isnothing(fixture_role) && push!(missing, :fixture_role)
    return Tuple(missing)
end

function _pqs_multilayer_complete_core_shell_rhf_missing_blocker(missing_inputs)
    isempty(missing_inputs) && return nothing
    length(missing_inputs) > 1 && return :missing_rhf_input_contract_inputs
    only_missing = only(missing_inputs)
    only_missing === :source_plan && return :missing_source_plan
    only_missing === :final_basis && return :missing_final_basis
    only_missing === :h1_payload && return :missing_h1_payload
    only_missing === :density_inputs && return :missing_density_inputs
    only_missing === :coulomb_expansion && return :missing_coulomb_expansion
    only_missing === :electron_count && return :missing_electron_count
    only_missing === :fixture_role && return :missing_fixture_role
    return :missing_rhf_input_contract_inputs
end

function _pqs_multilayer_complete_core_shell_rhf_input_contract(;
    source_plan = nothing,
    final_basis = nothing,
    h1_payload = nothing,
    density_inputs = nothing,
    coulomb_expansion = nothing,
    electron_count = nothing,
    fixture_role = nothing,
    metadata = (;),
)
    object_kind = :pqs_multilayer_complete_core_shell_rhf_input_contract
    source_plan_status =
        _pqs_multilayer_rhf_contract_property(source_plan, :status)
    final_basis_status =
        _pqs_multilayer_rhf_contract_property(final_basis, :status)
    h1_payload_status =
        _pqs_multilayer_rhf_contract_property(h1_payload, :status)
    density_inputs_status =
        _pqs_multilayer_rhf_contract_property(density_inputs, :status)
    coulomb_coefficients =
        _pqs_multilayer_rhf_contract_property(coulomb_expansion, :coefficients)
    coulomb_term_count = isnothing(coulomb_coefficients) ?
                         nothing : length(coulomb_coefficients)
    density_term_count =
        _pqs_multilayer_rhf_contract_summary_property(
            density_inputs,
            :term_count,
        )
    final_dimension =
        _pqs_multilayer_rhf_contract_property(final_basis, :final_retained_count)

    missing_inputs =
        _pqs_multilayer_complete_core_shell_rhf_missing_inputs(
            ;
            source_plan,
            final_basis,
            h1_payload,
            density_inputs,
            coulomb_expansion,
            electron_count,
            fixture_role,
        )
    occupation, occupation_blocker =
        _pqs_multilayer_complete_core_shell_rhf_occupation(electron_count)
    fixture_role_blocker =
        isnothing(fixture_role) ? :missing_fixture_role :
        fixture_role in _PQS_COMPLETE_CORE_SHELL_RHF_FIXTURE_ROLES ? nothing :
        :unsupported_fixture_role
    term_count_mismatch =
        !isnothing(coulomb_term_count) &&
        !isnothing(density_term_count) &&
        coulomb_term_count != density_term_count

    blocker =
        !isempty(missing_inputs) ?
        _pqs_multilayer_complete_core_shell_rhf_missing_blocker(missing_inputs) :
        !isnothing(occupation_blocker) ? occupation_blocker :
        !isnothing(fixture_role_blocker) ? fixture_role_blocker :
        term_count_mismatch ? :coulomb_expansion_term_count_mismatch :
        nothing
    status =
        isnothing(blocker) ?
        :available_pqs_multilayer_complete_core_shell_rhf_input_contract :
        :blocked_pqs_multilayer_complete_core_shell_rhf_input_contract
    summary = (;
        status,
        blocker,
        source_plan_status,
        final_basis_status,
        h1_payload_status,
        density_inputs_status,
        final_dimension,
        coulomb_term_count,
        density_term_count,
        coulomb_term_count_compatible =
            !isnothing(coulomb_term_count) &&
            !isnothing(density_term_count) &&
            coulomb_term_count == density_term_count,
        electron_count,
        occupation_status =
            isnothing(occupation_blocker) ? :available_closed_shell_rhf :
            occupation_blocker,
        fixture_role,
        fixture_role_status =
            isnothing(fixture_role_blocker) ? :available_fixture_role :
            fixture_role_blocker,
        fock_materialized = false,
        scf_materialized = false,
        rhf_energy_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        private_diagnostic_only = true,
    )
    return (;
        object_kind,
        status,
        blocker,
        missing_inputs,
        electron_count,
        occupation,
        fixture_role,
        summary,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source = :pqs_multilayer_complete_core_shell_rhf_input_contract,
                closed_shell_rhf_only = true,
                electron_count_source = :explicit,
                fock_materialized = false,
                scf_materialized = false,
                rhf_energy_materialized = false,
                driver_route_materialized = false,
                exports_materialized = false,
                artifacts_materialized = false,
            ),
        ),
    )
end
