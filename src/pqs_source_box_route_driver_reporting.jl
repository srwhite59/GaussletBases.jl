# Text report and optional artifact helpers for `bin/cartesian_ham_builder.jl`.
#
# Keep these after report construction so the core route helper file can focus
# on setup and construction.


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

function _pqs_source_box_route_driver_print_materialization(materialization)
    _pqs_route_driver_print_section("route_materialization")
    for field in (
        :route_family,
        :result_kind,
        :requested,
        :materialized,
        :retained_dimension,
        :final_dimension,
        :h1_lowest,
        :h1_symmetry_error,
        :overlap_identity_error,
        :provider_blocks_included,
        :residual_rank,
        :augmented_dimension,
        :augmented_h1_lowest,
        :augmented_h1_symmetry_error,
        :augmented_density_dimension,
        :augmented_pair_matrix_symmetry_error,
        :augmented_h1_j_self_coulomb,
        :private_augmented_rhf_converged,
        :private_augmented_rhf_iterations,
        :private_augmented_rhf_electronic_energy,
        :private_augmented_rhf_total_with_nuclear_repulsion,
        :private_augmented_rhf_commutator_residual,
        :basisfile,
        :hamfile,
    )
        hasproperty(materialization, field) ||
            continue
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


function _pqs_source_box_route_driver_durable_report(report)
    return (;
        (
            field =>
                field in (
                    :parent_basis_object,
                    :parent_axis_bundle_object,
                    :pqs_gto_sidecar_inputs,
                ) ?
                     nothing :
                     getproperty(report, field) for field in keys(report)
        )...,
        heavy_objects_elided_for_durable_report = true,
    )
end

function _pqs_source_box_route_driver_save(
    report;
    save_artifact, save_tsv, outfile, tsvfile, materialization = nothing,
    input_path = nothing,
)
    durable_report = _pqs_source_box_route_driver_durable_report(report)

    if save_artifact
        println("saving JLD2 report ", outfile)
        jldopen(outfile, "w") do file
            file["report"] = durable_report
            isnothing(materialization) ||
                (file["materialization"] = materialization)
        end
    end

    if save_tsv
        println("saving TSV report ", tsvfile)
        open(tsvfile, "w") do io
            println(io, "section\tkey\tvalue")
            _pqs_route_driver_write_named_tuple_tsv(io, "system_metadata", report.system_metadata)
            _pqs_route_driver_write_named_tuple_tsv(io, "recipe_metadata", report.recipe_metadata)
            _pqs_route_driver_write_named_tuple_tsv(io, "parent_summary", report.parent_summary)
            _pqs_route_driver_write_named_tuple_tsv(io, "route_summary", report.route_summary)
            _pqs_route_driver_write_named_tuple_tsv(io, "pair_summary", report.pair_summary)
            for unit in report.retained_units
                _pqs_route_driver_write_tsv_row(io, "retained_unit", unit.unit_key, unit)
            end
            for entry in report.pair_entries
                _pqs_route_driver_write_tsv_row(io, "pair_entry", entry.pair_key, entry)
            end
            if !isnothing(materialization)
                for field in keys(materialization)
                    _pqs_route_driver_write_tsv_row(
                        io,
                        "route_materialization",
                        field,
                        getproperty(materialization, field),
                    )
                end
            end
        end
    end
    return nothing
end
