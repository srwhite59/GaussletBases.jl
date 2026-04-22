function _quick_radial_operator_fixture()
    return _cached_fixture(:quick_radial_operator_fixture, () -> begin
        _radial_operator_fixture()
    end)
end

function _quick_radial_atomic_fixture()
    return _cached_fixture(:quick_radial_atomic_fixture, () -> begin
        rb, grid = _quick_radial_operator_fixture()
        radial_ops = atomic_operators(rb, grid; Z = 2.0, lmax = 2)
        channels = ylm_channels(2)
        atom = atomic_one_body_operators(radial_ops, channels)
        ida = atomic_ida_operators(radial_ops, channels)
        (rb, grid, radial_ops, channels, atom, ida)
    end)
end
