# Text report and optional artifact helpers for `bin/cartesian_ham_builder.jl`.
#
# Keep these after report construction so the core route helper file can focus
# on setup, skeletons, diagnostics, and dry-run report assembly.


# Text report helpers. These are intentionally simple print utilities, not a
# general logging framework.

function _pqs_route_driver_print_section(title)
    println()
    println("[", title, "]")
    return nothing
end

function _pqs_route_driver_print_kv(key, value)
    println(rpad(String(key), 42), "  ", value)
    return nothing
end

function _pqs_route_driver_print_named_tuple(title, values)
    _pqs_route_driver_print_section(title)
    for field in keys(values)
        _pqs_route_driver_print_kv(field, getproperty(values, field))
    end
    return nothing
end

function _pqs_source_box_route_driver_materialization_status_fields()
    return (
        :status,
        :materialized_report_kind,
        :retained_dimension,
        :pqs_materialization_status,
    )
end

function _pqs_source_box_route_driver_print_materialization(materialization)
    _pqs_route_driver_print_section("route_materialization")
    for field in _pqs_source_box_route_driver_materialization_status_fields()
        _pqs_route_driver_print_kv(field, getproperty(materialization, field))
    end
    return nothing
end


# Optional artifact helpers for the private driver.

function _pqs_route_driver_write_tsv_row(io, section, key, value)
    println(io, section, '\t', key, '\t', repr(value))
    return nothing
end

function _pqs_route_driver_write_named_tuple_tsv(io, section, values)
    for field in keys(values)
        _pqs_route_driver_write_tsv_row(io, section, field, getproperty(values, field))
    end
    return nothing
end


function _pqs_source_box_route_driver_durable_materializer_payload(payload)
    isnothing(payload) && return nothing
    return (;
        object_kind = payload.object_kind,
        private_development_only = payload.private_development_only,
        transient_only = payload.transient_only,
        durable_report_serialization = :sanitized_before_save,
        source = payload.source,
        parent_basis_object_available = payload.parent_basis_object_available,
        parent_qw_basis_object_available = payload.parent_qw_basis_object_available,
        parent_axis_bundle_object_available =
            payload.parent_axis_bundle_object_available,
        parent_basis_object_type_label = payload.parent_basis_object_type_label,
        parent_qw_basis_object_type_label = payload.parent_qw_basis_object_type_label,
        parent_axis_bundle_object_type_label =
            payload.parent_axis_bundle_object_type_label,
        axis_bundle_backend = payload.axis_bundle_backend,
        axis_bundle_backend_available = payload.axis_bundle_backend_available,
        heavy_objects_elided = true,
    )
end

function _pqs_source_box_route_driver_durable_report(report)
    hasproperty(report, :route_materializer_payload) || return report
    return (;
        (
            field => field == :route_materializer_payload ?
                     _pqs_source_box_route_driver_durable_materializer_payload(
                         report.route_materializer_payload,
                     ) :
                     getproperty(report, field) for field in keys(report)
        )...,
    )
end

function _pqs_source_box_route_driver_durable_diatomic_materialization(materialization)
    isnothing(materialization) && return nothing
    hasproperty(materialization, :source) || return materialization
    sanitized = (;
        (
            field => field == :source ?
                     nothing :
                     getproperty(materialization, field) for
            field in keys(materialization)
        )...,
    )
    return merge(sanitized, (; source_elided_for_durable_report = true))
end

function _pqs_source_box_route_driver_durable_diatomic_probe(probe)
    isnothing(probe) && return nothing
    hasproperty(probe, :materialization) || return probe
    return (;
        (
            field => field == :materialization ?
                     _pqs_source_box_route_driver_durable_diatomic_materialization(
                         probe.materialization,
                     ) :
                     getproperty(probe, field) for field in keys(probe)
        )...,
    )
end

function _pqs_source_box_route_driver_durable_materialization(materialization)
    isnothing(materialization) && return nothing
    hasproperty(materialization, :route_configured_diatomic_materializer_probe) ||
        return materialization
    return (;
        (
            field => field == :route_configured_diatomic_materializer_probe ?
                     _pqs_source_box_route_driver_durable_diatomic_probe(
                         materialization.route_configured_diatomic_materializer_probe,
                     ) :
                     getproperty(materialization, field) for field in keys(materialization)
        )...,
    )
end

function _pqs_source_box_route_driver_save(
    report;
    save_artifact, save_tsv, outfile, tsvfile, materialization = nothing,
    input_path = nothing,
)
    durable_report = _pqs_source_box_route_driver_durable_report(report)
    durable_materialization =
        _pqs_source_box_route_driver_durable_materialization(materialization)

    if save_artifact
        println("saving JLD2 report ", outfile)
        jldopen(outfile, "w") do file
            file["report"] = durable_report
            isnothing(durable_materialization) ||
                (file["materialization"] = durable_materialization)
        end
    end

    if save_tsv
        println("saving TSV report ", tsvfile)
        open(tsvfile, "w") do io
            println(io, "section\tkey\tvalue")
            _pqs_route_driver_write_named_tuple_tsv(io, "system_metadata", report.system_metadata)
            _pqs_route_driver_write_named_tuple_tsv(io, "recipe_metadata", report.recipe_metadata)
            _pqs_route_driver_write_named_tuple_tsv(
                io, "parent_contract", report.parent_contract)
            _pqs_route_driver_write_named_tuple_tsv(
                io, "standard_setup_diagnostics", report.standard_setup.diagnostics,
            )
            _pqs_route_driver_write_named_tuple_tsv(
                io, "parent_axis_readiness_diagnostics",
                report.parent_axis_readiness.diagnostics,
            )
            if !isnothing(report.parent_axis_probe)
                _pqs_route_driver_write_named_tuple_tsv(
                    io, "parent_axis_probe_diagnostics",
                    report.parent_axis_probe.diagnostics,
                )
            end
            _pqs_route_driver_write_named_tuple_tsv(
                io, "route_axis_counts_diagnostics", report.route_axis_counts.diagnostics,
            )
            if !isnothing(report.raw_product_box_probe)
                _pqs_route_driver_write_named_tuple_tsv(
                    io, "raw_product_box_probe_diagnostics",
                    report.raw_product_box_probe.diagnostics,
                )
            end
            _pqs_route_driver_write_named_tuple_tsv(
                io, "standard_unit_inventory", report.standard_unit_inventory)
            for unit in report.retained_units
                _pqs_route_driver_write_tsv_row(io, "retained_unit", unit.unit_key, unit)
            end
            for entry in report.pair_entries
                _pqs_route_driver_write_tsv_row(io, "pair_entry", entry.pair_key, entry)
            end
            if !isnothing(materialization)
                for field in _pqs_source_box_route_driver_materialization_status_fields()
                    _pqs_route_driver_write_tsv_row(
                        io,
                        "route_materialization",
                        field,
                        getproperty(materialization, field),
                    )
                end
            end
            _pqs_route_driver_write_named_tuple_tsv(io, "diagnostics", report.diagnostics)
        end
    end
    return nothing
end
