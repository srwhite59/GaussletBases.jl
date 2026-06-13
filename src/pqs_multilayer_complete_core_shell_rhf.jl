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

function _pqs_multilayer_complete_core_shell_rhf_one_step_missing_inputs(;
    input_contract,
    h1_payload,
    density_interaction,
    final_density,
)
    missing = Symbol[]
    contract_status =
        _pqs_multilayer_rhf_contract_property(input_contract, :status)
    contract_kind =
        _pqs_multilayer_rhf_contract_property(input_contract, :object_kind)
    if contract_kind !==
       :pqs_multilayer_complete_core_shell_rhf_input_contract ||
       contract_status !==
       :available_pqs_multilayer_complete_core_shell_rhf_input_contract
        push!(missing, :rhf_input_contract)
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

    density_status =
        _pqs_multilayer_rhf_contract_property(density_interaction, :status)
    density_kind =
        _pqs_multilayer_rhf_contract_property(density_interaction, :object_kind)
    pair_matrix =
        _pqs_multilayer_rhf_contract_property(
            density_interaction,
            :pre_final_pair_matrix,
        )
    final_to_pre_final =
        _pqs_multilayer_rhf_contract_property(
            density_interaction,
            :final_to_pre_final_coefficients,
        )
    if density_kind !== :pqs_complete_core_shell_pre_final_density_interaction ||
       density_status !==
       :materialized_pqs_complete_core_shell_pre_final_density_interaction ||
       isnothing(pair_matrix) ||
       isnothing(final_to_pre_final)
        push!(missing, :density_interaction)
    end

    isnothing(final_density) && push!(missing, :final_density)
    return Tuple(missing)
end

function _pqs_multilayer_complete_core_shell_rhf_one_step_missing_blocker(
    missing_inputs,
)
    isempty(missing_inputs) && return nothing
    length(missing_inputs) > 1 && return :missing_rhf_one_step_inputs
    only_missing = only(missing_inputs)
    only_missing === :rhf_input_contract && return :missing_rhf_input_contract
    only_missing === :h1_payload && return :missing_h1_payload
    only_missing === :density_interaction && return :missing_density_interaction
    only_missing === :final_density && return :missing_final_density
    return :missing_rhf_one_step_inputs
end

function _pqs_multilayer_complete_core_shell_rhf_one_step_blocked_payload(;
    blocker,
    missing_inputs = (),
    input_contract = nothing,
    h1_payload = nothing,
    density_interaction = nothing,
    final_density = nothing,
    final_density_trace = nothing,
    final_density_symmetry_error = nothing,
    metadata = (;),
)
    contract_summary =
        _pqs_multilayer_rhf_contract_property(input_contract, :summary, (;))
    final_dimension =
        _pqs_multilayer_rhf_contract_property(contract_summary, :final_dimension)
    electron_count =
        _pqs_multilayer_rhf_contract_property(input_contract, :electron_count)
    fixture_role =
        _pqs_multilayer_rhf_contract_property(input_contract, :fixture_role)
    summary = (;
        status = :blocked_pqs_multilayer_complete_core_shell_rhf_one_step_payload,
        blocker,
        missing_inputs,
        input_contract_status =
            _pqs_multilayer_rhf_contract_property(input_contract, :status),
        h1_payload_status =
            _pqs_multilayer_rhf_contract_property(h1_payload, :status),
        density_interaction_status =
            _pqs_multilayer_rhf_contract_property(density_interaction, :status),
        final_dimension,
        electron_count,
        fixture_role,
        final_density_trace,
        final_density_symmetry_error,
        fock_materialized = false,
        one_step_energy_materialized = false,
        scf_materialized = false,
        rhf_converged = false,
        rhf_energy_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        private_diagnostic_only = true,
    )
    return (;
        object_kind =
            :pqs_multilayer_complete_core_shell_rhf_one_step_payload,
        status = summary.status,
        blocker,
        missing_inputs,
        final_density,
        effective_fock_matrix = nothing,
        fock_matrix = nothing,
        one_body_energy = nothing,
        two_body_energy = nothing,
        total_energy = nothing,
        summary,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source =
                    :pqs_multilayer_complete_core_shell_rhf_one_step_payload,
                blocker,
                fock_materialized = false,
                one_step_energy_materialized = false,
                scf_materialized = false,
                rhf_converged = false,
                rhf_energy_materialized = false,
                driver_route_materialized = false,
                exports_materialized = false,
                artifacts_materialized = false,
            ),
        ),
    )
end

function _pqs_multilayer_complete_core_shell_rhf_density_interaction_payload(
    density_interaction,
    h1_j_payload,
)
    !isnothing(density_interaction) && return density_interaction
    direct =
        _pqs_multilayer_rhf_contract_property(h1_j_payload, :density_interaction)
    !isnothing(direct) && return direct
    nested =
        _pqs_multilayer_rhf_contract_property(h1_j_payload, :h1_j_payload)
    return _pqs_multilayer_rhf_contract_property(
        nested,
        :density_interaction,
    )
end

function _pqs_multilayer_complete_core_shell_rhf_h1_payload(
    h1_payload,
    h1_j_payload,
)
    !isnothing(h1_payload) && return h1_payload
    direct =
        _pqs_multilayer_rhf_contract_property(h1_j_payload, :h1_payload)
    !isnothing(direct) && return direct
    nested =
        _pqs_multilayer_rhf_contract_property(
            h1_j_payload,
            :complete_core_shell_diagnostic_route_payload,
        )
    return _pqs_multilayer_rhf_contract_property(nested, :h1_payload)
end

function _pqs_multilayer_complete_core_shell_rhf_initial_density_blocked_payload(;
    blocker,
    missing_inputs = (),
    input_contract = nothing,
    h1_payload = nothing,
    metadata = (;),
)
    contract_summary =
        _pqs_multilayer_rhf_contract_property(input_contract, :summary, (;))
    final_dimension =
        _pqs_multilayer_rhf_contract_property(contract_summary, :final_dimension)
    electron_count =
        _pqs_multilayer_rhf_contract_property(input_contract, :electron_count)
    fixture_role =
        _pqs_multilayer_rhf_contract_property(input_contract, :fixture_role)
    summary = (;
        status =
            :blocked_pqs_multilayer_complete_core_shell_rhf_initial_density_payload,
        blocker,
        missing_inputs,
        input_contract_status =
            _pqs_multilayer_rhf_contract_property(input_contract, :status),
        h1_payload_status =
            _pqs_multilayer_rhf_contract_property(h1_payload, :status),
        final_dimension,
        electron_count,
        fixture_role,
        initial_density_source = :h1_aufbau,
        initial_density_materialized = false,
        scf_materialized = false,
        rhf_converged = false,
        rhf_energy_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        private_diagnostic_only = true,
    )
    return (;
        object_kind =
            :pqs_multilayer_complete_core_shell_rhf_initial_density_payload,
        status = summary.status,
        blocker,
        missing_inputs,
        occupied_orbital_coefficients = nothing,
        final_density = nothing,
        eigenvalues = nothing,
        occupied_eigenvalues = nothing,
        electron_trace = nothing,
        summary,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source =
                    :pqs_multilayer_complete_core_shell_rhf_initial_density_payload,
                blocker,
                initial_density_source = :h1_aufbau,
                initial_density_materialized = false,
                scf_materialized = false,
                rhf_converged = false,
                rhf_energy_materialized = false,
                driver_route_materialized = false,
                exports_materialized = false,
                artifacts_materialized = false,
            ),
        ),
    )
end

function _pqs_multilayer_complete_core_shell_rhf_initial_density_payload(;
    input_contract = nothing,
    h1_payload = nothing,
    h1_j_payload = nothing,
    metadata = (;),
)
    h1_payload =
        _pqs_multilayer_complete_core_shell_rhf_h1_payload(
            h1_payload,
            h1_j_payload,
        )
    missing_inputs = Symbol[]
    contract_kind =
        _pqs_multilayer_rhf_contract_property(input_contract, :object_kind)
    contract_status =
        _pqs_multilayer_rhf_contract_property(input_contract, :status)
    if contract_kind !==
       :pqs_multilayer_complete_core_shell_rhf_input_contract ||
       contract_status !==
       :available_pqs_multilayer_complete_core_shell_rhf_input_contract
        push!(missing_inputs, :rhf_input_contract)
    end

    h1_payload_kind =
        _pqs_multilayer_rhf_contract_property(h1_payload, :object_kind)
    h1_payload_status =
        _pqs_multilayer_rhf_contract_property(h1_payload, :status)
    final_hamiltonian =
        _pqs_multilayer_rhf_contract_property(h1_payload, :final_hamiltonian)
    h1_matrix =
        _pqs_multilayer_rhf_contract_property(
            final_hamiltonian,
            :hamiltonian_matrix,
        )
    if h1_payload_kind !== :pqs_multilayer_complete_core_shell_h1_payload ||
       h1_payload_status !==
       :materialized_pqs_multilayer_complete_core_shell_h1_payload ||
       isnothing(h1_matrix)
        push!(missing_inputs, :h1_payload)
    end
    if !isempty(missing_inputs)
        blocker =
            length(missing_inputs) == 1 && only(missing_inputs) ===
            :rhf_input_contract ? :missing_rhf_input_contract :
            length(missing_inputs) == 1 && only(missing_inputs) ===
            :h1_payload ? :missing_h1_payload :
            :missing_rhf_initial_density_inputs
        return _pqs_multilayer_complete_core_shell_rhf_initial_density_blocked_payload(
            ;
            blocker,
            missing_inputs = Tuple(missing_inputs),
            input_contract,
            h1_payload,
            metadata,
        )
    end

    final_dimension = input_contract.summary.final_dimension
    matrix = Matrix{Float64}(h1_matrix)
    if size(matrix, 1) != size(matrix, 2) ||
       size(matrix, 1) != final_dimension
        return _pqs_multilayer_complete_core_shell_rhf_initial_density_blocked_payload(
            ;
            blocker = :h1_dimension_mismatch,
            input_contract,
            h1_payload,
            metadata,
        )
    end
    if !all(isfinite, matrix)
        return _pqs_multilayer_complete_core_shell_rhf_initial_density_blocked_payload(
            ;
            blocker = :nonfinite_h1_matrix,
            input_contract,
            h1_payload,
            metadata,
        )
    end

    occupation = input_contract.occupation
    nocc = occupation.nocc
    occupancy = occupation.occupancy
    if nocc > final_dimension
        return _pqs_multilayer_complete_core_shell_rhf_initial_density_blocked_payload(
            ;
            blocker = :insufficient_final_dimension_for_occupation,
            input_contract,
            h1_payload,
            metadata,
        )
    end

    symmetric_h1 = 0.5 .* (matrix .+ transpose(matrix))
    eig = eigen(Symmetric(symmetric_h1))
    occupied_orbital_coefficients = Matrix(eig.vectors[:, 1:nocc])
    final_density =
        Float64(occupancy) .*
        (occupied_orbital_coefficients * transpose(occupied_orbital_coefficients))
    final_density = 0.5 .* (final_density .+ transpose(final_density))
    eigenvalues = Vector{Float64}(eig.values)
    occupied_eigenvalues = eigenvalues[1:nocc]
    electron_trace = tr(final_density)

    summary = (;
        status =
            :materialized_pqs_multilayer_complete_core_shell_rhf_initial_density_payload,
        blocker = nothing,
        input_contract_status = input_contract.status,
        h1_payload_status = h1_payload.status,
        final_dimension,
        electron_count = input_contract.electron_count,
        nocc,
        occupancy,
        electron_trace,
        fixture_role = input_contract.fixture_role,
        initial_density_source = :h1_aufbau,
        lowest_eigenvalue = first(eigenvalues),
        highest_eigenvalue = last(eigenvalues),
        lowest_occupied_eigenvalue = first(occupied_eigenvalues),
        highest_occupied_eigenvalue = last(occupied_eigenvalues),
        initial_density_materialized = true,
        scf_materialized = false,
        rhf_converged = false,
        rhf_energy_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        private_diagnostic_only = true,
    )
    return (;
        object_kind =
            :pqs_multilayer_complete_core_shell_rhf_initial_density_payload,
        status = summary.status,
        blocker = nothing,
        missing_inputs = (),
        occupied_orbital_coefficients,
        final_density,
        eigenvalues,
        occupied_eigenvalues,
        electron_trace,
        summary,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source =
                    :pqs_multilayer_complete_core_shell_rhf_initial_density_payload,
                initial_density_source = :h1_aufbau,
                initial_density_materialized = true,
                scf_materialized = false,
                rhf_converged = false,
                rhf_energy_materialized = false,
                driver_route_materialized = false,
                exports_materialized = false,
                artifacts_materialized = false,
            ),
        ),
    )
end

function _pqs_multilayer_complete_core_shell_rhf_scf_controls_blocker(
    max_iterations,
    density_atol,
    energy_atol,
)
    if !(max_iterations isa Integer) ||
       max_iterations isa Bool ||
       max_iterations <= 0
        return :invalid_scf_controls
    end
    if !(density_atol isa Real) ||
       !isfinite(density_atol) ||
       density_atol <= 0
        return :invalid_scf_controls
    end
    if !(energy_atol isa Real) ||
       !isfinite(energy_atol) ||
       energy_atol <= 0
        return :invalid_scf_controls
    end
    return nothing
end

function _pqs_multilayer_complete_core_shell_rhf_scf_blocked_payload(;
    blocker,
    missing_inputs = (),
    input_contract = nothing,
    h1_payload = nothing,
    density_interaction = nothing,
    initial_density_payload = nothing,
    final_density = nothing,
    occupied_orbital_coefficients = nothing,
    final_one_step_payload = nothing,
    iteration_records = (),
    max_iterations = nothing,
    density_atol = nothing,
    energy_atol = nothing,
    metadata = (;),
)
    contract_summary =
        _pqs_multilayer_rhf_contract_property(input_contract, :summary, (;))
    final_dimension =
        _pqs_multilayer_rhf_contract_property(contract_summary, :final_dimension)
    electron_count =
        _pqs_multilayer_rhf_contract_property(input_contract, :electron_count)
    fixture_role =
        _pqs_multilayer_rhf_contract_property(input_contract, :fixture_role)
    summary = (;
        status = :blocked_pqs_multilayer_complete_core_shell_rhf_scf_payload,
        blocker,
        missing_inputs,
        input_contract_status =
            _pqs_multilayer_rhf_contract_property(input_contract, :status),
        h1_payload_status =
            _pqs_multilayer_rhf_contract_property(h1_payload, :status),
        density_interaction_status =
            _pqs_multilayer_rhf_contract_property(density_interaction, :status),
        initial_density_status =
            _pqs_multilayer_rhf_contract_property(
                initial_density_payload,
                :status,
            ),
        final_dimension,
        electron_count,
        fixture_role,
        iteration_count = length(iteration_records),
        max_iterations,
        density_atol,
        energy_atol,
        rhf_materialized = false,
        rhf_converged = false,
        driver_route_materialized = false,
        route_report_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        public_api = false,
        private_diagnostic_only = true,
    )
    return (;
        object_kind = :pqs_multilayer_complete_core_shell_rhf_scf_payload,
        status = summary.status,
        blocker,
        missing_inputs,
        final_density,
        occupied_orbital_coefficients,
        final_one_step_payload,
        iteration_records,
        summary,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source = :pqs_multilayer_complete_core_shell_rhf_scf_payload,
                blocker,
                rhf_materialized = false,
                rhf_converged = false,
                driver_route_materialized = false,
                route_report_materialized = false,
                exports_materialized = false,
                artifacts_materialized = false,
                public_api = false,
            ),
        ),
    )
end

function _pqs_multilayer_complete_core_shell_rhf_scf_payload(;
    input_contract = nothing,
    h1_payload = nothing,
    h1_j_payload = nothing,
    density_interaction = nothing,
    initial_density_payload = nothing,
    max_iterations::Integer = 25,
    density_atol::Real = 1.0e-8,
    energy_atol::Real = 1.0e-10,
    metadata = (;),
)
    controls_blocker =
        _pqs_multilayer_complete_core_shell_rhf_scf_controls_blocker(
            max_iterations,
            density_atol,
            energy_atol,
        )
    if !isnothing(controls_blocker)
        return _pqs_multilayer_complete_core_shell_rhf_scf_blocked_payload(
            ;
            blocker = controls_blocker,
            max_iterations,
            density_atol,
            energy_atol,
            metadata,
        )
    end

    h1_payload =
        _pqs_multilayer_complete_core_shell_rhf_h1_payload(
            h1_payload,
            h1_j_payload,
        )
    density_interaction =
        _pqs_multilayer_complete_core_shell_rhf_density_interaction_payload(
            density_interaction,
            h1_j_payload,
        )
    missing_inputs = Symbol[]
    input_contract_status =
        _pqs_multilayer_rhf_contract_property(input_contract, :status)
    input_contract_kind =
        _pqs_multilayer_rhf_contract_property(input_contract, :object_kind)
    if input_contract_kind !==
       :pqs_multilayer_complete_core_shell_rhf_input_contract ||
       input_contract_status !==
       :available_pqs_multilayer_complete_core_shell_rhf_input_contract
        push!(missing_inputs, :rhf_input_contract)
    end
    h1_payload_status =
        _pqs_multilayer_rhf_contract_property(h1_payload, :status)
    h1_payload_kind =
        _pqs_multilayer_rhf_contract_property(h1_payload, :object_kind)
    if h1_payload_kind !== :pqs_multilayer_complete_core_shell_h1_payload ||
       h1_payload_status !==
       :materialized_pqs_multilayer_complete_core_shell_h1_payload
        push!(missing_inputs, :h1_payload)
    end
    density_status =
        _pqs_multilayer_rhf_contract_property(density_interaction, :status)
    density_kind =
        _pqs_multilayer_rhf_contract_property(density_interaction, :object_kind)
    if density_kind !== :pqs_complete_core_shell_pre_final_density_interaction ||
       density_status !==
       :materialized_pqs_complete_core_shell_pre_final_density_interaction
        push!(missing_inputs, :density_interaction)
    end
    if !isempty(missing_inputs)
        blocker =
            length(missing_inputs) == 1 && only(missing_inputs) ===
            :rhf_input_contract ? :missing_rhf_input_contract :
            length(missing_inputs) == 1 && only(missing_inputs) ===
            :h1_payload ? :missing_h1_payload :
            length(missing_inputs) == 1 && only(missing_inputs) ===
            :density_interaction ? :missing_density_interaction :
            :missing_rhf_scf_inputs
        return _pqs_multilayer_complete_core_shell_rhf_scf_blocked_payload(
            ;
            blocker,
            missing_inputs = Tuple(missing_inputs),
            input_contract,
            h1_payload,
            density_interaction,
            max_iterations,
            density_atol,
            energy_atol,
            metadata,
        )
    end

    if isnothing(initial_density_payload)
        initial_density_payload =
            _pqs_multilayer_complete_core_shell_rhf_initial_density_payload(
                ;
                input_contract,
                h1_payload,
                metadata,
            )
    end
    if _pqs_multilayer_rhf_contract_property(
        initial_density_payload,
        :status,
    ) !==
       :materialized_pqs_multilayer_complete_core_shell_rhf_initial_density_payload
        blocker =
            _pqs_multilayer_rhf_contract_property(
                initial_density_payload,
                :blocker,
                :missing_initial_density,
            )
        return _pqs_multilayer_complete_core_shell_rhf_scf_blocked_payload(
            ;
            blocker,
            missing_inputs = (:initial_density,),
            input_contract,
            h1_payload,
            density_interaction,
            initial_density_payload,
            max_iterations,
            density_atol,
            energy_atol,
            metadata,
        )
    end

    density = Matrix{Float64}(initial_density_payload.final_density)
    occupation = input_contract.occupation
    nocc = occupation.nocc
    occupancy = occupation.occupancy
    final_dimension = input_contract.summary.final_dimension
    previous_energy = nothing
    iteration_records = NamedTuple[]
    final_one_step_payload = nothing
    occupied_orbital_coefficients =
        initial_density_payload.occupied_orbital_coefficients

    for iteration in 1:max_iterations
        one_step =
            _pqs_multilayer_complete_core_shell_rhf_one_step_payload(
                ;
                input_contract,
                h1_payload,
                density_interaction,
                final_density = density,
                metadata,
            )
        if one_step.status !==
           :materialized_pqs_multilayer_complete_core_shell_rhf_one_step_payload
            return _pqs_multilayer_complete_core_shell_rhf_scf_blocked_payload(
                ;
                blocker = isnothing(one_step.blocker) ?
                    :missing_density_interaction : one_step.blocker,
                input_contract,
                h1_payload,
                density_interaction,
                initial_density_payload,
                final_density = density,
                occupied_orbital_coefficients,
                final_one_step_payload = one_step,
                iteration_records = Tuple(iteration_records),
                max_iterations,
                density_atol,
                energy_atol,
                metadata,
            )
        end

        fock_matrix =
            0.5 .* (
                one_step.effective_fock_matrix .+
                transpose(one_step.effective_fock_matrix)
            )
        eig = eigen(Symmetric(fock_matrix))
        occupied_orbital_coefficients = Matrix(eig.vectors[:, 1:nocc])
        next_density =
            Float64(occupancy) .*
            (
                occupied_orbital_coefficients *
                transpose(occupied_orbital_coefficients)
            )
        next_density = 0.5 .* (next_density .+ transpose(next_density))
        density_change = norm(next_density - density, Inf)
        energy_change =
            isnothing(previous_energy) ? nothing :
            abs(one_step.total_energy - previous_energy)
        density_converged = density_change <= Float64(density_atol)
        energy_converged =
            isnothing(energy_change) || energy_change <= Float64(energy_atol)
        converged = density_converged && energy_converged
        push!(
            iteration_records,
            (;
                iteration,
                total_energy = one_step.total_energy,
                density_change,
                energy_change,
                density_converged,
                energy_converged,
                converged,
            ),
        )

        final_one_step_payload = one_step
        density = next_density
        previous_energy = one_step.total_energy
        if converged
            summary = (;
                status =
                    :materialized_pqs_multilayer_complete_core_shell_rhf_scf_payload,
                blocker = nothing,
                input_contract_status = input_contract.status,
                h1_payload_status = h1_payload.status,
                density_interaction_status = density_interaction.status,
                initial_density_status = initial_density_payload.status,
                final_dimension,
                electron_count = input_contract.electron_count,
                electron_trace = tr(density),
                fixture_role = input_contract.fixture_role,
                iteration_count = length(iteration_records),
                converged_iteration = iteration,
                max_iterations,
                density_atol,
                energy_atol,
                first_iteration_energy_change_rule =
                    :energy_converged_when_previous_energy_missing,
                final_total_energy = one_step.total_energy,
                final_density_change = density_change,
                final_energy_change = energy_change,
                rhf_materialized = true,
                rhf_converged = true,
                driver_route_materialized = false,
                route_report_materialized = false,
                exports_materialized = false,
                artifacts_materialized = false,
                public_api = false,
                private_diagnostic_only = true,
            )
            return (;
                object_kind =
                    :pqs_multilayer_complete_core_shell_rhf_scf_payload,
                status = summary.status,
                blocker = nothing,
                missing_inputs = (),
                final_density = density,
                occupied_orbital_coefficients,
                final_one_step_payload,
                iteration_records = Tuple(iteration_records),
                summary,
                metadata = merge(
                    NamedTuple(metadata),
                    (;
                        source =
                            :pqs_multilayer_complete_core_shell_rhf_scf_payload,
                        rhf_materialized = true,
                        rhf_converged = true,
                        driver_route_materialized = false,
                        route_report_materialized = false,
                        exports_materialized = false,
                        artifacts_materialized = false,
                        public_api = false,
                    ),
                ),
            )
        end
    end

    return _pqs_multilayer_complete_core_shell_rhf_scf_blocked_payload(
        ;
        blocker = :scf_not_converged,
        input_contract,
        h1_payload,
        density_interaction,
        initial_density_payload,
        final_density = density,
        occupied_orbital_coefficients,
        final_one_step_payload,
        iteration_records = Tuple(iteration_records),
        max_iterations,
        density_atol,
        energy_atol,
        metadata,
    )
end

function _pqs_multilayer_complete_core_shell_rhf_one_step_payload(;
    input_contract = nothing,
    h1_payload = nothing,
    density_interaction = nothing,
    h1_j_payload = nothing,
    final_density = nothing,
    density_trace_atol::Real = 1.0e-8,
    symmetry_atol::Real = 1.0e-8,
    metadata = (;),
)
    h1_payload =
        _pqs_multilayer_complete_core_shell_rhf_h1_payload(
            h1_payload,
            h1_j_payload,
        )
    density_interaction =
        _pqs_multilayer_complete_core_shell_rhf_density_interaction_payload(
            density_interaction,
            h1_j_payload,
        )
    missing_inputs =
        _pqs_multilayer_complete_core_shell_rhf_one_step_missing_inputs(
            ;
            input_contract,
            h1_payload,
            density_interaction,
            final_density,
        )
    if !isempty(missing_inputs)
        return _pqs_multilayer_complete_core_shell_rhf_one_step_blocked_payload(
            ;
            blocker =
                _pqs_multilayer_complete_core_shell_rhf_one_step_missing_blocker(
                    missing_inputs,
                ),
            missing_inputs,
            input_contract,
            h1_payload,
            density_interaction,
            final_density,
            metadata,
        )
    end

    contract_summary = input_contract.summary
    final_dimension = contract_summary.final_dimension
    density_matrix = Matrix{Float64}(final_density)
    if size(density_matrix, 1) != size(density_matrix, 2)
        return _pqs_multilayer_complete_core_shell_rhf_one_step_blocked_payload(
            ;
            blocker = :nonsquare_final_density,
            input_contract,
            h1_payload,
            density_interaction,
            final_density = density_matrix,
            metadata,
        )
    end
    if size(density_matrix, 1) != final_dimension
        return _pqs_multilayer_complete_core_shell_rhf_one_step_blocked_payload(
            ;
            blocker = :final_density_dimension_mismatch,
            input_contract,
            h1_payload,
            density_interaction,
            final_density = density_matrix,
            metadata,
        )
    end
    if !all(isfinite, density_matrix)
        return _pqs_multilayer_complete_core_shell_rhf_one_step_blocked_payload(
            ;
            blocker = :nonfinite_final_density,
            input_contract,
            h1_payload,
            density_interaction,
            final_density = density_matrix,
            metadata,
        )
    end

    density_symmetry_error =
        norm(density_matrix - transpose(density_matrix), Inf)
    if density_symmetry_error > Float64(symmetry_atol)
        return _pqs_multilayer_complete_core_shell_rhf_one_step_blocked_payload(
            ;
            blocker = :nonsymmetric_final_density,
            input_contract,
            h1_payload,
            density_interaction,
            final_density = density_matrix,
            final_density_symmetry_error = density_symmetry_error,
            metadata,
        )
    end

    density_trace = tr(density_matrix)
    electron_count = input_contract.electron_count
    if abs(density_trace - electron_count) > Float64(density_trace_atol)
        return _pqs_multilayer_complete_core_shell_rhf_one_step_blocked_payload(
            ;
            blocker = :electron_trace_mismatch,
            input_contract,
            h1_payload,
            density_interaction,
            final_density = density_matrix,
            final_density_trace = density_trace,
            final_density_symmetry_error = density_symmetry_error,
            metadata,
        )
    end

    h1_matrix = Matrix{Float64}(h1_payload.final_hamiltonian.hamiltonian_matrix)
    size(h1_matrix) == (final_dimension, final_dimension) ||
        return _pqs_multilayer_complete_core_shell_rhf_one_step_blocked_payload(
            ;
            blocker = :h1_dimension_mismatch,
            input_contract,
            h1_payload,
            density_interaction,
            final_density = density_matrix,
            final_density_trace = density_trace,
            final_density_symmetry_error = density_symmetry_error,
            metadata,
        )
    all(isfinite, h1_matrix) ||
        return _pqs_multilayer_complete_core_shell_rhf_one_step_blocked_payload(
            ;
            blocker = :nonfinite_h1_matrix,
            input_contract,
            h1_payload,
            density_interaction,
            final_density = density_matrix,
            final_density_trace = density_trace,
            final_density_symmetry_error = density_symmetry_error,
            metadata,
        )

    final_to_pre_final =
        Matrix{Float64}(density_interaction.final_to_pre_final_coefficients)
    pair_matrix = Matrix{Float64}(density_interaction.pre_final_pair_matrix)
    size(final_to_pre_final, 2) == final_dimension ||
        return _pqs_multilayer_complete_core_shell_rhf_one_step_blocked_payload(
            ;
            blocker = :final_to_pre_final_dimension_mismatch,
            input_contract,
            h1_payload,
            density_interaction,
            final_density = density_matrix,
            final_density_trace = density_trace,
            final_density_symmetry_error = density_symmetry_error,
            metadata,
        )
    pre_final_dimension = size(final_to_pre_final, 1)
    size(pair_matrix) == (pre_final_dimension, pre_final_dimension) ||
        return _pqs_multilayer_complete_core_shell_rhf_one_step_blocked_payload(
            ;
            blocker = :density_interaction_dimension_mismatch,
            input_contract,
            h1_payload,
            density_interaction,
            final_density = density_matrix,
            final_density_trace = density_trace,
            final_density_symmetry_error = density_symmetry_error,
            metadata,
        )
    all(isfinite, final_to_pre_final) && all(isfinite, pair_matrix) ||
        return _pqs_multilayer_complete_core_shell_rhf_one_step_blocked_payload(
            ;
            blocker = :nonfinite_density_interaction,
            input_contract,
            h1_payload,
            density_interaction,
            final_density = density_matrix,
            final_density_trace = density_trace,
            final_density_symmetry_error = density_symmetry_error,
            metadata,
        )

    occupancy = input_contract.occupation.occupancy
    orbital_density_final = density_matrix ./ Float64(occupancy)
    orbital_density_pre =
        final_to_pre_final * orbital_density_final * transpose(final_to_pre_final)
    orbital_density_pre =
        0.5 .* (orbital_density_pre .+ transpose(orbital_density_pre))
    interaction_matrix = 0.5 .* (pair_matrix .+ transpose(pair_matrix))
    occupations = vec(diag(orbital_density_pre))
    coulomb_pre =
        2.0 .* Diagonal(interaction_matrix * occupations) .-
        orbital_density_pre .* interaction_matrix
    coulomb_final =
        transpose(final_to_pre_final) * Matrix(coulomb_pre) * final_to_pre_final
    effective_fock_matrix = h1_matrix + coulomb_final
    one_body_energy = tr(density_matrix * h1_matrix)
    two_body_energy = 0.5 * tr(density_matrix * coulomb_final)
    total_energy = one_body_energy + two_body_energy

    summary = (;
        status =
            :materialized_pqs_multilayer_complete_core_shell_rhf_one_step_payload,
        blocker = nothing,
        input_contract_status = input_contract.status,
        h1_payload_status = h1_payload.status,
        density_interaction_status = density_interaction.status,
        final_dimension,
        pre_final_dimension,
        electron_count,
        density_trace,
        final_density_symmetry_error = density_symmetry_error,
        fixture_role = input_contract.fixture_role,
        density_convention = :spin_summed_closed_shell_final_density,
        contraction_rule =
            :pre_final_restricted_direct_minus_exchange_from_orbital_density,
        h1_matrix_materialized = true,
        fock_materialized = true,
        one_step_energy_materialized = true,
        scf_materialized = false,
        rhf_converged = false,
        rhf_energy_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        private_diagnostic_only = true,
    )
    return (;
        object_kind =
            :pqs_multilayer_complete_core_shell_rhf_one_step_payload,
        status = summary.status,
        blocker = nothing,
        missing_inputs = (),
        final_density = density_matrix,
        effective_fock_matrix,
        fock_matrix = effective_fock_matrix,
        one_body_energy,
        two_body_energy,
        total_energy,
        summary,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source =
                    :pqs_multilayer_complete_core_shell_rhf_one_step_payload,
                density_convention = summary.density_convention,
                contraction_rule = summary.contraction_rule,
                final_to_pre_final_coefficients_source =
                    :combined_lowdin_cleanup_times_final_coefficients,
                h1_matrix_source =
                    :pqs_multilayer_complete_core_shell_h1_payload,
                fock_materialized = true,
                one_step_energy_materialized = true,
                scf_materialized = false,
                rhf_converged = false,
                rhf_energy_materialized = false,
                driver_route_materialized = false,
                exports_materialized = false,
                artifacts_materialized = false,
            ),
        ),
    )
end
