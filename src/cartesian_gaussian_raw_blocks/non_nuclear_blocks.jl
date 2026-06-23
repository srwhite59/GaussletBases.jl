function gaussian_non_nuclear_raw_blocks(proxy, supplement, expansion)
    cross = getfield(_GB_PARENT, :_qwrg_cartesian_shell_cross_moment_blocks_3d)(
        (x = proxy.x, y = proxy.y, z = proxy.z),
        supplement,
        expansion,
        proxy.ncart;
        include_factor_terms = false,
    )
    self = getfield(_GB_PARENT, :_qwrg_cartesian_shell_self_moment_blocks_3d)(
        supplement,
        expansion;
        include_factor_terms = false,
    )
    return (;
        ga = (;
            overlap = cross.overlap_ga,
            kinetic = cross.kinetic_ga,
            position = (x = cross.position_x_ga, y = cross.position_y_ga,
                z = cross.position_z_ga),
            x2 = (x = cross.x2_x_ga, y = cross.x2_y_ga, z = cross.x2_z_ga),
        ),
        aa = (;
            overlap = _symmetrize_raw_block(self.overlap_aa),
            kinetic = _symmetrize_raw_block(self.kinetic_aa),
            position = (x = _symmetrize_raw_block(self.position_x_aa),
                y = _symmetrize_raw_block(self.position_y_aa),
                z = _symmetrize_raw_block(self.position_z_aa)),
            x2 = (x = _symmetrize_raw_block(self.x2_x_aa),
                y = _symmetrize_raw_block(self.x2_y_aa),
                z = _symmetrize_raw_block(self.x2_z_aa)),
        ),
    )
end
