if !isdefined(@__MODULE__, :_CARTESIAN_PAIR_STAGE_FINGERPRINT_HELPERS_INCLUDED)
const _CARTESIAN_PAIR_STAGE_FINGERPRINT_HELPERS_INCLUDED = true

function _cartesian_pair_stage_property(record, field, default = nothing)
    return hasproperty(record, field) ? getproperty(record, field) : default
end

function _cartesian_pair_stage_fingerprint(record, spec)
    names = Tuple(item[1] for item in spec)
    values = Tuple(
        _cartesian_pair_stage_property(record, item[2], item[3])
        for item in spec
    )
    return NamedTuple{names}(values)
end

function _cartesian_pair_stage_contract_counts_fingerprint(entries)
    return Tuple(
        (;
            unit_key = _cartesian_pair_stage_property(entry, :unit_key),
            lowering_contract_count = _cartesian_pair_stage_property(
                entry,
                :lowering_contract_count,
            ),
            selected_contract_count = _cartesian_pair_stage_property(
                entry,
                :selected_contract_count,
            ),
        )
        for entry in entries
    )
end

const _CARTESIAN_PAIR_STAGE_TERMINAL_UNIT_FIELDS = (
    (:available, :terminal_shellification_unit_inventory_available, false),
    (:region_count, :terminal_shellification_region_count, 0),
    (:unit_count, :terminal_shellification_unit_count, 0),
    (:unit_keys, :terminal_shellification_unit_keys, ()),
    (:unit_roles, :terminal_shellification_unit_roles, ()),
    (:unit_kinds, :terminal_shellification_unit_kinds, ()),
    (:support_counts, :terminal_shellification_unit_support_counts, ()),
    (
        :central_gap_region_count,
        :terminal_shellification_central_gap_region_count,
        0,
    ),
    (
        :central_midpoint_slab_count,
        :terminal_shellification_central_midpoint_slab_count,
        0,
    ),
    (
        :central_distorted_product_box_count,
        :terminal_shellification_central_distorted_product_box_count,
        0,
    ),
)

function _cartesian_pair_stage_terminal_unit_inventory_fingerprint(stage)
    return _cartesian_pair_stage_fingerprint(
        stage,
        _CARTESIAN_PAIR_STAGE_TERMINAL_UNIT_FIELDS,
    )
end

const _CARTESIAN_PAIR_STAGE_LOWERING_CONTRACT_FIELDS = (
    (
        :available,
        :terminal_shellification_lowering_contract_inventory_available,
        false,
    ),
    (:status, :terminal_shellification_lowering_contract_inventory_status, nothing),
    (:count, :terminal_shellification_lowering_contract_count, 0),
    (:kinds, :terminal_shellification_lowering_contract_kinds, ()),
    (:kind_counts, :terminal_shellification_lowering_contract_kind_counts, ()),
    (
        :lw_complete_shell_cpb_count,
        :terminal_shellification_lw_complete_shell_cpb_count,
        0,
    ),
    (
        :lw_complete_shell_cpb_family_counts,
        :terminal_shellification_lw_complete_shell_cpb_family_counts,
        (),
    ),
)

function _cartesian_pair_stage_lowering_contract_fingerprint(stage)
    base = _cartesian_pair_stage_fingerprint(
        stage,
        _CARTESIAN_PAIR_STAGE_LOWERING_CONTRACT_FIELDS,
    )
    return merge(
        base,
        (;
            counts_by_unit = _cartesian_pair_stage_contract_counts_fingerprint(
                _cartesian_pair_stage_property(
                    stage,
                    :terminal_shellification_contract_counts_by_unit,
                    (),
                ),
            ),
        ),
    )
end

const _CARTESIAN_PAIR_STAGE_SELECTED_CONTRACT_FIELDS = (
    (
        :available,
        :terminal_shellification_selected_lowering_contract_inventory_available,
        false,
    ),
    (
        :status,
        :terminal_shellification_selected_lowering_contract_inventory_status,
        nothing,
    ),
    (:family, :terminal_shellification_selected_lowering_family, nothing),
    (:count, :terminal_shellification_selected_contract_count, 0),
    (:kinds, :terminal_shellification_selected_contract_kinds, ()),
    (:kind_counts, :terminal_shellification_selected_contract_kind_counts, ()),
    (
        :all_units_have_exactly_one_selected_contract,
        :terminal_shellification_all_units_have_exactly_one_selected_contract,
        false,
    ),
    (:unselected_count, :terminal_shellification_unselected_contract_count, 0),
    (:unselected_kinds, :terminal_shellification_unselected_contract_kinds, ()),
)

function _cartesian_pair_stage_selected_contract_fingerprint(stage)
    base = _cartesian_pair_stage_fingerprint(
        stage,
        _CARTESIAN_PAIR_STAGE_SELECTED_CONTRACT_FIELDS,
    )
    return merge(
        base,
        (;
            counts_by_unit = _cartesian_pair_stage_contract_counts_fingerprint(
                _cartesian_pair_stage_property(
                    stage,
                    :terminal_shellification_selected_contract_counts_by_unit,
                    (),
                ),
            ),
        ),
    )
end

const _CARTESIAN_PAIR_STAGE_CRC_SIDECAR_FIELDS = (
    (:object_kind, :object_kind, nothing),
    (:status, :status, nothing),
    (:selected_contract_count, :selected_contract_count, 0),
    (:sidecar_available_count, :sidecar_available_count, 0),
    (:sidecar_missing_count, :sidecar_missing_count, 0),
    (:sidecar_inventory_complete, :sidecar_inventory_complete, false),
    (:missing_sidecar_reasons, :missing_sidecar_reasons, ()),
    (:missing_sidecar_kinds, :missing_sidecar_kinds, ()),
    (
        :final_retained_unit_inventory_available,
        :final_retained_unit_inventory_available,
        false,
    ),
    (:pair_inventory_available, :pair_inventory_available, false),
    (:pair_inventory_status, :pair_inventory_status, nothing),
    (:operator_blocks_materialized, :operator_blocks_materialized, false),
    (:pair_operator_blocks_materialized, :pair_operator_blocks_materialized, false),
    (:hamiltonian_data_materialized, :hamiltonian_data_materialized, false),
    (:artifacts_materialized, :artifacts_materialized, false),
)

function _cartesian_pair_stage_crc_sidecar_summary_fingerprint(summary)
    return _cartesian_pair_stage_fingerprint(
        summary,
        _CARTESIAN_PAIR_STAGE_CRC_SIDECAR_FIELDS,
    )
end

const _CARTESIAN_PAIR_STAGE_PAIR_METADATA_FIELDS = (
    (:pair_inventory_known, :pair_inventory_known, false),
    (:pair_inventory_source, :pair_inventory_source, nothing),
    (:pair_family_counts, :pair_family_counts, ()),
    (:helper_by_pair_family, :helper_by_pair_family, ()),
    (:pair_operator_helper_by_family, :pair_operator_helper_by_family, ()),
    (:pair_helper_status_by_family, :pair_helper_status_by_family, ()),
    (:pair_operator_blocks_materialized, :pair_operator_blocks_materialized, false),
    (:operator_pairs_materialized, :operator_pairs_materialized, false),
    (
        :route_core_pair_inventory_available,
        :route_core_pair_inventory_available,
        false,
    ),
    (:route_core_pair_inventory_status, :route_core_pair_inventory_status, nothing),
    (:route_core_pair_count, :route_core_pair_count, 0),
    (:route_core_pair_operator_ready, :route_core_pair_operator_ready, false),
    (
        :route_core_pair_operator_readiness_status,
        :route_core_pair_operator_readiness_status,
        nothing,
    ),
    (:route_core_pair_operator_blocker, :route_core_pair_operator_blocker, nothing),
)

function _cartesian_pair_stage_pair_metadata_fingerprint(stage)
    pair_entries = _cartesian_pair_stage_property(stage, :pair_entries, ())
    route_core_pair_keys =
        _cartesian_pair_stage_property(stage, :route_core_pair_keys, ())
    base = _cartesian_pair_stage_fingerprint(
        stage,
        _CARTESIAN_PAIR_STAGE_PAIR_METADATA_FIELDS,
    )
    return merge(
        base,
        (;
            pair_entry_count = length(pair_entries),
            pair_count = _cartesian_pair_stage_property(
                stage,
                :pair_count,
                length(pair_entries),
            ),
            route_core_pair_key_count = length(route_core_pair_keys),
            route_core_pair_keys,
        ),
    )
end

end
