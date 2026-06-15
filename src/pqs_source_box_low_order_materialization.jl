function _pqs_source_box_route_driver_materialization(
    report;
    materialize_route::Bool = false,
    save_basis_artifact::Bool = false,
    save_ham_artifact::Bool = false,
    basisfile = nothing,
    hamfile = nothing,
    kwargs...,
)
    requested = materialize_route || save_basis_artifact || save_ham_artifact
    retained_dimension =
        hasproperty(report, :retained_dimension) ? report.retained_dimension : nothing
    return (;
        object_kind = :cartesian_nesting_route_driver_materialization,
        route_family = hasproperty(report, :route_family) ? report.route_family : nothing,
        status =
            requested ?
            :blocked_materialization_after_route_scaffold_demolition :
            :not_requested,
        blocker =
            requested ?
            :route_configured_low_order_materializer_removed :
            nothing,
        materialized_report = nothing,
        materialized_report_kind = nothing,
        retained_dimension,
        pqs_materialization_status =
            requested ?
            :blocked_materialization_after_route_scaffold_demolition :
            :not_requested,
        save_basis_artifact_requested = save_basis_artifact,
        save_ham_artifact_requested = save_ham_artifact,
        basisfile,
        hamfile,
        ignored_legacy_keyword_count = length(kwargs),
    )
end
