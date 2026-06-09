# Private/nonproduction override file for bin/cartesian_ham_builder.jl.
#
# Usage:
#
#     julia --project=. bin/cartesian_ham_builder.jl examples/private_global_overlap_option.jl
#
# This demonstrates the option shape only. The tiny numerical factors below are
# test-style inputs for a 2 x 2 final-local overlap fixture; they are not a
# production Be2 route configuration and should not be assumed to match the
# default driver route. A real run must provide overlap factors and a global
# dimension consistent with its actual pair-block materialization plan.

private_global_overlap_requested = true
private_global_overlap_global_dimension = 2
private_global_overlap_inputs = (;
    parent_axis_counts = (2, 2, 2),
    overlap_1d = (;
        x = [1.0 0.2; 0.2 1.1],
        y = [1.2 0.3; 0.3 1.3],
        z = [1.4 0.4; 0.4 1.5],
    ),
)

# The driver returns structured private_global_overlap_result and
# private_global_overlap_summary fields in the materialization result. This is
# not production driver wiring, not Hamiltonian assembly, not term summing, and
# not Coulomb, IDA/MWG, PQS Lowdin/projection, export, or artifact work.
