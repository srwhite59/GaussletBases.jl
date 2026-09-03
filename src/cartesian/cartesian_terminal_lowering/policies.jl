# Terminal lowering policies.

abstract type TerminalLoweringPolicy end

"""
    WhiteLindseyLowering()

Select White--Lindsey-style lowering for complete terminal shells. Direct
terminal regions and central distorted product boxes keep their direct/product
contracts.
"""
struct WhiteLindseyLowering <: TerminalLoweringPolicy end

"""
    PQSLowering(; q)

Select PQS lowering for complete terminal shells. `q` records the product-mode
order used by the source-box retained-rule contract.
"""
struct PQSLowering <: TerminalLoweringPolicy
    q::Int

    function PQSLowering(q::Integer)
        q > 0 || throw(ArgumentError("PQS lowering q must be positive"))
        return new(Int(q))
    end
end

PQSLowering(; q::Integer) = PQSLowering(q)

policy_kind(::WhiteLindseyLowering) = :white_lindsey_terminal_lowering
policy_kind(::PQSLowering) = :pqs_terminal_lowering
