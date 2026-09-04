# Cartesian/PQS Source Relocation

`HP-SOURCE-LAYOUT-CARTESIAN-FN-01` and
`HP-SOURCE-LAYOUT-CARTESIAN-TEST-01` authorize only the Cartesian/PQS half of
source-layout Step 4 from baseline
`96802b00cc380707306106b529ff22a09ccb415e`. The external follow-up review is
audit evidence, not authority. Independent inspection fixes the closure at
`76` tracked files, `40,736` lines, and `1,671,542` bytes: `28` flat owners and
`48` files contained in `12` tracked submodule directories.

One current docstring prevents a truthful baseline-to-final claim that every
moved blob is unchanged. Before moving files, repo-manager must make one local
prerequisite commit that changes only this sentence in
`src/cartesian_terminal_shellification_geometry.jl`:

```text
The implementation lives in `src/cartesian_shellification/`.
```

to the path-neutral sentence:

```text
The implementation is provided by the internal `CartesianShellification` module.
```

That file's baseline blob/SHA-256 are
`ff1ba22edc86b1b594c90a75b5157ba86bcbf073` and
`a3c8ba89840acf58df4cda47885f5796dba1c5bbcf90cc4933d5d78866ed73a5`;
the required post-edit blob/SHA-256 are
`3c83d034a01ebd46f5c9fe49fdc09b9369a8092f` and
`cf40eb13810a7b3da9e5bd36cd907adc1c10ea23657d8ba4bc4cbefd331d52be`.
No executable line or other source body may change. The following relocation
commit must then use `git mv` and report all `76` files as `100%` renames
relative to that prerequisite commit. Both local commits are pushed together
so the source-bearing transaction incurs one remote full-matrix cycle.

The frozen baseline map is:

| Pre-move path | Authorized target | Git blob | SHA-256 |
| --- | --- | --- | --- |
| `src/CartesianParentAxisFactors.jl` | `src/cartesian/CartesianParentAxisFactors.jl` | `d797e5c1facd400c64f5291e4d6198dc865fa55a` | `5c70d7923668ff9d89e7ad38002ba617f31c39fdff71fc66c45912f2ac788acc` |
| `src/CartesianParentGaussletBases.jl` | `src/cartesian/CartesianParentGaussletBases.jl` | `f0854f87180e81f5e3806a25b203012d359802ce` | `0cc2839cead709c0789d74b9643a04307f90ac33000e748e163d68e06ba63d12` |
| `src/bond_aligned_diatomic_geometry.jl` | `src/cartesian/bond_aligned_diatomic_geometry.jl` | `ad0609f1ea70949f608f5659e3a1a461ffe36a28` | `1a53d0ac3a94ec2384a540c9089bfeb522fc02026fdb819c27d6e9532a332512` |
| `src/cartesian_base_hamiltonian.jl` | `src/cartesian/cartesian_base_hamiltonian.jl` | `97e9494a6d6f815a97dcf4f6bcbd161558b6ea1c` | `cfd23c492f2422080f6c4a681a6f3595a4280a61485bab3b29291b33f4101b9a` |
| `src/cartesian_basis_representation.jl` | `src/cartesian/cartesian_basis_representation.jl` | `a95c6be6236a57e4fe50173564e1df43752699dd` | `8aed35aa4ac55238ebfc2934c0e5fe048bdd68257f7c81c8e308ffc1a4ada24b` |
| `src/cartesian_cpb/CartesianCPB.jl` | `src/cartesian/cartesian_cpb/CartesianCPB.jl` | `564d026e3b93461d588e29ec6029f9404c37f575` | `2495e9806a13d3432056f55e57a129bc00e08018820c00494a61ba80feb45634` |
| `src/cartesian_cpb/coordinate_product_boxes.jl` | `src/cartesian/cartesian_cpb/coordinate_product_boxes.jl` | `4f44c85a5bc994390f720404a85d73b509187194` | `12071157ed9844b2c0f3f4c5d3c1a3c3569f3a926550bd5df1003ae7b673046e` |
| `src/cartesian_cross_overlap.jl` | `src/cartesian/cartesian_cross_overlap.jl` | `2ddb2b0608e5ce168e057c08b91e3c96fbc71841` | `ba6b897886007016068cfd8108ca1b838fc3463be293343a01d0a69e61372255` |
| `src/cartesian_external_gto_import.jl` | `src/cartesian/cartesian_external_gto_import.jl` | `bbf61e324a638b620ccaeb0d3a8661e59c308250` | `062120391dc3c14a95b48978b6b311582aac6f75f36cee4212e89c1af5ae3704` |
| `src/cartesian_external_gto_interchange.jl` | `src/cartesian/cartesian_external_gto_interchange.jl` | `68d7d94664d9f1997ed0180ac7aa06f81c63af58` | `fda3cf43c40bc227f2a3777200ed7157e9725cb4a250f435469cd55cb3eddb7e` |
| `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl` | `src/cartesian/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl` | `92b0b722f874e94b73ec4bf429a2d20f7d5828e7` | `e95fa5d86c6965f3a3d455a877897a869bc93e8a2ae151c79f18bbc1b6ef4ea2` |
| `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl` | `src/cartesian/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl` | `7eb5b6e02f8eab33e7f58d957bd44611ce29597e` | `8223c9105881032204a8aa7e40099c5f8c614a68e6857d80f2659729bb6a8a34` |
| `src/cartesian_final_basis_realization/pqs_terminal_ida.jl` | `src/cartesian/cartesian_final_basis_realization/pqs_terminal_ida.jl` | `cbb2d44f9d784968b6a9ff6baec2b71cc8645f66` | `1b11dc51d8052430b20c20e31af73513a58770a52b194fe54d209a29a42811da` |
| `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl` | `src/cartesian/cartesian_final_basis_realization/pqs_terminal_one_body.jl` | `d7884a2b56bc8e43444fdf5073e2bd0e1eb62e4a` | `3ab3ebf48053191b29829e25fa378397ce0d41e8bdf3aacb2efbe6f22554b314` |
| `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` | `src/cartesian/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` | `3902ff1462b6f066f2c3ac01a37e08a2fcac5a70` | `8c48eb3086b47ac1a5116fe914fed9cefbb65a081e8da8fffe81fd83baed5c44` |
| `src/cartesian_final_basis_realization/terminal_face_product_blocks.jl` | `src/cartesian/cartesian_final_basis_realization/terminal_face_product_blocks.jl` | `304d70e31c1721e20c808c4dc32472557ab535fc` | `0c5a4c72ce17cae738d20e36296174b06f397c2786cf3e1c02097a33e48826c8` |
| `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl` | `src/cartesian/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl` | `c86d395db497c83a37cecc5afe04d3078b7dfca2` | `ef450ce51c3bb631ba87d6f3fbd3441c3fd3040ce46a8be26ec4adefce4f3ade` |
| `src/cartesian_gaussian_axis_integrals.jl` | `src/cartesian/cartesian_gaussian_axis_integrals.jl` | `63778d8a2be3aee973625b15a43f8507f0b20bb7` | `b8b0c1cc46324f7e09601c7cacff4cbcf1b96a790dc3fc3b15c450a3dd1d3fed` |
| `src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl` | `src/cartesian/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl` | `7644f29eab7788993ab5b05f49ce3f7086b0f16f` | `15eb073d38106e934131c1234508973ad99feaf62d730a9cd1e970fc926d701a` |
| `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl` | `src/cartesian/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl` | `3f9e4ba9e337c04b1d8b832b5494f54c3dd5d03f` | `5909ce6f0bfaca8b9ee4a75597ae17cf93eeec5da16417392b79c18b47eba5d6` |
| `src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl` | `src/cartesian/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl` | `f3474f0a1748337d14d938b03498df8314e3e199` | `cfb1610e2e8677a5d754156acdd5cd598a16ca39135929888940094a2783c0c2` |
| `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl` | `src/cartesian/cartesian_gaussian_raw_blocks/nuclear_blocks.jl` | `b2414ae994143b7b963935f43c5f2986504ac888` | `2b1cac9f5af79c0d2237eb82c97d850eeeeffc8f6821d1800c9c9548394090ce` |
| `src/cartesian_gto_probes.jl` | `src/cartesian/cartesian_gto_probes.jl` | `9de967865c1a1ef806bfdc8b5e2ffb82a0a0fc5d` | `5410aef19c4d134f3319aabd967be19aaaacd0dd476eb149726d17186d43e94c` |
| `src/cartesian_ida_hamiltonian.jl` | `src/cartesian/cartesian_ida_hamiltonian.jl` | `932d800912502604a941c6baa21ea6e46944bb41` | `8fb4046f2e82522ea644418d44fc23342f1d41727c9eda9c55094f85ed9b63c4` |
| `src/cartesian_nested_atomic.jl` | `src/cartesian/cartesian_nested_atomic.jl` | `e64e4fc51c928c826b92fd7f053633172fccd01f` | `b4db7df0ce26f94a6afa2dda0e839034ce1f6d84fad66997f6a36087d9c80f8a` |
| `src/cartesian_nested_diatomic.jl` | `src/cartesian/cartesian_nested_diatomic.jl` | `6b3432949b421a022efd4984c9beb456f5388993` | `fdac52901f6ace1de99520e870cda3f1eec6f125d26503af792f3b85b36fd140` |
| `src/cartesian_nested_experimental_geometries.jl` | `src/cartesian/cartesian_nested_experimental_geometries.jl` | `1cfe9023c71eccb615f6ec529b29e543e2f17de5` | `26aba50390ce4ecb66516b700cc1c1097270d1fc20f95e5537d5744daf452c68` |
| `src/cartesian_nested_faces.jl` | `src/cartesian/cartesian_nested_faces.jl` | `e2dd1827d2f71be602eb73bcd9f4db7214abe33e` | `91d6f81f6c2530fb86bacc94b7f0d21974b3f37f53e51319e7d34430eba76d92` |
| `src/cartesian_nested_owned_units.jl` | `src/cartesian/cartesian_nested_owned_units.jl` | `6ed0f927f593be33730d3e0282c1466b97ade80c` | `3f2b86bc359eb284351a54cd8efe7b77a187dc8924109ac10f8539e62a49eb80` |
| `src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl` | `src/cartesian/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl` | `e48a5f320314ac225a7d1c94dbf51d58cc94b14f` | `0255490a871433c99c2e65f245299763caa3097602acde18f7cdd596b8c92835` |
| `src/cartesian_pair_block_materialization/pqs_source_axis_transforms.jl` | `src/cartesian/cartesian_pair_block_materialization/pqs_source_axis_transforms.jl` | `00e4e8cbe2cf360859c18faf9a4738c9b580215f` | `23f5f622251a33285826c8df32098bfac910997a742da3376832fc80fdc6be66` |
| `src/cartesian_protected_ladder_bundle.jl` | `src/cartesian/cartesian_protected_ladder_bundle.jl` | `2200dee2d358bbd296e345f02c9c5a612535d2a7` | `521741d7fe3de3762c4b5cf74b61c099606349a73aedef014d9bb35868248c75` |
| `src/cartesian_qw_hybrid_representation.jl` | `src/cartesian/cartesian_qw_hybrid_representation.jl` | `db118cfafbdcdb3009093aa93a0a091cc33b337d` | `8663b9d3b0592b4346c5ffc0fccdf34a2db7e3cc656e46a3a6f96c2077e6f65b` |
| `src/cartesian_raw_product_sources/CartesianRawProductSources.jl` | `src/cartesian/cartesian_raw_product_sources/CartesianRawProductSources.jl` | `2afe574816d91319045408dece53820dd5d88f8a` | `9de2a226586a6235a7bbd16607ba4ed338077ccb69448ab578e69a44192bdc25` |
| `src/cartesian_raw_product_sources/axis_transform_facts.jl` | `src/cartesian/cartesian_raw_product_sources/axis_transform_facts.jl` | `66e35a50cc868612369217844bec666673d331d4` | `350185507dde8c8ffea41a26203fabbf6f3aa63bceb1c287339b9225b95a062a` |
| `src/cartesian_raw_product_sources/records.jl` | `src/cartesian/cartesian_raw_product_sources/records.jl` | `e8db8bb183750e9fdb15c8e9b653cb3a8f03448e` | `b26c1cb4234e1f5866b43158955ad0f0eddb7b0412853bdcd07b1aa6ca481ca2` |
| `src/cartesian_raw_product_sources/source_mode_indices.jl` | `src/cartesian/cartesian_raw_product_sources/source_mode_indices.jl` | `3ff9af60c6499155e1f07ce7135386e29dbf1b4e` | `795994dbe360e4c298f0afc86efd403e309c1dab0bd2035235d8ac65e1a3e4f9` |
| `src/cartesian_raw_product_sources/summaries.jl` | `src/cartesian/cartesian_raw_product_sources/summaries.jl` | `959a8b87c424ded8a3e9a3926d55a8f623cebb7b` | `f83a30ff5755e7393da71daf26977229f5f60a1a9e73f5ea262974d62ac38e63` |
| `src/cartesian_reference_density/CartesianReferenceDensity.jl` | `src/cartesian/cartesian_reference_density/CartesianReferenceDensity.jl` | `e87422f8b7aefe43447121f567083955b15f4f91` | `de67c6dc92d1c4399bb893fc46b24bdb0c1be41db4adc81f446e25df785f2a5f` |
| `src/cartesian_reference_density/atomic_hf_reference_packets.jl` | `src/cartesian/cartesian_reference_density/atomic_hf_reference_packets.jl` | `9e49f686ce81808a482c259c8acde51f13ab82c7` | `a3d679058fa53cf5907d1ff57610f9a96cb2415fa5e506b55947e8d30f63b903` |
| `src/cartesian_reference_density/represented_molecular_hartree.jl` | `src/cartesian/cartesian_reference_density/represented_molecular_hartree.jl` | `8bab7bf8353f4c1c5906a3a366857f694e22be5a` | `e4f8f0b54b73f3cc3696ce600a267012cd76b17f7e6f19aee11571fe48680ac1` |
| `src/cartesian_reference_density/screened_hartree_correction.jl` | `src/cartesian/cartesian_reference_density/screened_hartree_correction.jl` | `009ae11378514ace6495080941b1b75a0508c135` | `517e5bf800a638f2c8e63e66f2080a7c28128c8e36eaf40010424aab8d4ae306` |
| `src/cartesian_representation_constructors.jl` | `src/cartesian/cartesian_representation_constructors.jl` | `29a98fb5d632f56a658257d423fadc686aceed9c` | `3d5fe0634243b2ff8f731290625d5c118e72bf5b2416646f5480cbf305c53a2f` |
| `src/cartesian_representation_transfer.jl` | `src/cartesian/cartesian_representation_transfer.jl` | `cda545dc4ccd0cd4b70d89521a8a1a121a1a089e` | `ba1733d9e756ae4887b50127811ea793c9589268ab40f9d6bdc36dc109f485a6` |
| `src/cartesian_residual_gaussians/CartesianResidualGaussians.jl` | `src/cartesian/cartesian_residual_gaussians/CartesianResidualGaussians.jl` | `eb02ef25fd5c394e83abdfe150e0f1b4e2cd81b3` | `535df854b26ccd61ec4752b1a021f2a1a52ee8d71b5c3f71f404e8b1c7c78011` |
| `src/cartesian_residual_gaussians/augmented_operators.jl` | `src/cartesian/cartesian_residual_gaussians/augmented_operators.jl` | `253b48f0a2ed6a604e688a830071957673ec95fd` | `b0cfa6b120ebb46eeb0659591a917ba4cfd60df22559b916b7cc6cf0fd06d30b` |
| `src/cartesian_residual_gaussians/mwg_interaction.jl` | `src/cartesian/cartesian_residual_gaussians/mwg_interaction.jl` | `35397f52094657a1862cdf82abd4c2d0e80d91c2` | `dfe94391052da0ce3ae2f2e266304d80d8c7aa5941b5b6805c2d5145d2028b08` |
| `src/cartesian_residual_gaussians/residual_basis.jl` | `src/cartesian/cartesian_residual_gaussians/residual_basis.jl` | `7c81c2370e5228f0750ce0af6f0fac3f7b8f5f35` | `390d8a10831f65f8f61eeb67784ce6b61a5ec314e713f7de471a779b975716d3` |
| `src/cartesian_retained_unit_transform_contracts/CartesianRetainedUnitTransformContracts.jl` | `src/cartesian/cartesian_retained_unit_transform_contracts/CartesianRetainedUnitTransformContracts.jl` | `aa90211cfc104b345e1df41a29aab80b7b6c05a4` | `b1bb9e046dc442396960256d911c97dd73bde51f0ee3fde1a68ac1b7584d659c` |
| `src/cartesian_retained_unit_transform_contracts/records.jl` | `src/cartesian/cartesian_retained_unit_transform_contracts/records.jl` | `2f312ea9876b62697eab2911fee4b266454dc7c4` | `141d6c467a7873a40490ea1224335f4d87d6f28573ba5d632a889448377c3dcb` |
| `src/cartesian_retained_unit_transform_contracts/summaries.jl` | `src/cartesian/cartesian_retained_unit_transform_contracts/summaries.jl` | `54f8b721c692c27a7a1c30ce7570532b4d18c701` | `251a25a815b13bb20c3f3bf59cb7622d360e78758326b33270355e9e07fb6a68` |
| `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl` | `src/cartesian/cartesian_retained_unit_transform_contracts/unit_contracts.jl` | `ee1fd22d4bde6d08d290a891bc65fa61bcb56bec` | `ec95ee327931a01b2df18f27e1849a7cdf5951c6eb1c6c6f9d6871d263971f67` |
| `src/cartesian_retained_units/CartesianRetainedUnits.jl` | `src/cartesian/cartesian_retained_units/CartesianRetainedUnits.jl` | `889c10a96fe8b7a1dec125b2b79bf94c08817088` | `a7930e34cfaf16c5340075f9ef2fa91808ca6051d0f7ad2f6a8f34da0af3a383` |
| `src/cartesian_retained_units/lower_contract_units.jl` | `src/cartesian/cartesian_retained_units/lower_contract_units.jl` | `cecf4e5e56c98a915427088fc3237bb728b83615` | `08355a1b80db3ed0a69173e3053fdfeb71c5cf4f6f82f41e81d1a57f2016058d` |
| `src/cartesian_retained_units/records.jl` | `src/cartesian/cartesian_retained_units/records.jl` | `494348f3f4d112ae3b6fec97923cd62f1fa925fe` | `f071d6663879389976949d9f565b6b1215a5b9370cbb470ed6f972cf1af58777` |
| `src/cartesian_retained_units/summaries.jl` | `src/cartesian/cartesian_retained_units/summaries.jl` | `4ce1c3f8bdbc201809f12b44962f7aa7b397e2e2` | `1a2264cd0e4e244d500a647ab0ea10d5f3446615a0085fd35379ae1f8ae931ed` |
| `src/cartesian_route_core/CartesianRouteCore.jl` | `src/cartesian/cartesian_route_core/CartesianRouteCore.jl` | `67c97fb1377b2eb009a96d9c5fa6b244f48f3253` | `006435b7b40682e882377da00578cc88d59b9d248894c6bdc94c780d4b3b7f09` |
| `src/cartesian_route_core/lowering_sources.jl` | `src/cartesian/cartesian_route_core/lowering_sources.jl` | `127987b2f1cced43112fa55f0b584adbade9aff3` | `ebf41b63db95f7d2e856cbb96d8985e0022600f20fc589b558b011fb80dea73b` |
| `src/cartesian_route_core/retained_spaces.jl` | `src/cartesian/cartesian_route_core/retained_spaces.jl` | `73024f30b13a07ce73406b3523ea05f959270499` | `3ca0764870601708bf284903d786464bbaf58cdcdba7d84233d4e78cb5d01bf2` |
| `src/cartesian_route_core/shellification_regions.jl` | `src/cartesian/cartesian_route_core/shellification_regions.jl` | `c5311e9c9e46ebe90fa064b07a4143f5f6d62fcb` | `0b48addf184604ac8ec6fefc2e3fdd604db8963273a022e6d3efda212a6b8473` |
| `src/cartesian_shellification/CartesianShellification.jl` | `src/cartesian/cartesian_shellification/CartesianShellification.jl` | `87c4b3a78a2952fe46f688d22153203958a11b6d` | `406e39767f22b9d70bb19472bfe32f2bf9a933cc75e763bf3e3a18965c472b24` |
| `src/cartesian_shellification/terminal_geometry.jl` | `src/cartesian/cartesian_shellification/terminal_geometry.jl` | `7fd076d64bbbecd34ec007797bb62a71345c2739` | `516cb56b4f4842dec52154a68d4dc4520ef637207113016e684acefea5f5c4f7` |
| `src/cartesian_terminal_lowering/CartesianTerminalLowering.jl` | `src/cartesian/cartesian_terminal_lowering/CartesianTerminalLowering.jl` | `9c8a4a3d9f734ada30cae09d98030a334a1a8f32` | `64d04272d9649bdb5a014b0aba7865a37ffc0a027e853122990b1fa0fb549863` |
| `src/cartesian_terminal_lowering/contracts.jl` | `src/cartesian/cartesian_terminal_lowering/contracts.jl` | `95b5f607558c87154f8f47cfa06fab2e1c811624` | `e36267648d77eca10bbf4cfc60e937ebea7439ca3abac62601a752ea9a858ae1` |
| `src/cartesian_terminal_lowering/policies.jl` | `src/cartesian/cartesian_terminal_lowering/policies.jl` | `6a53736bdefd2a72d45b22215dc5d7f3afea0914` | `561c5bd833a54f80030ad8399f5e0cccf6a0f09ebc9044ec0890aef1ff145dae` |
| `src/cartesian_terminal_lowering/region_contracts.jl` | `src/cartesian/cartesian_terminal_lowering/region_contracts.jl` | `bdc44de72792c6853fbd0bbeb0fc566ce53e5a95` | `ec235409e237e08f46e04e1bb2ef1036f3060781efa8d637bd46b2eb81b487c1` |
| `src/cartesian_terminal_lowering/selection.jl` | `src/cartesian/cartesian_terminal_lowering/selection.jl` | `bdaa6ea8dc23b36d0f33579daa62e48f0b3cc28c` | `881736b188e1182c103ed73d17b8c6bcb3d8da5188e3b7b43139d99d0426e0f5` |
| `src/cartesian_terminal_lowering/summaries.jl` | `src/cartesian/cartesian_terminal_lowering/summaries.jl` | `bff840819a2a4eb166a4122ac6b4c30d8bb7db03` | `e69250d55a3b332ed4e9aa188811d1d0235e07061cd3e38539aa240be2529454` |
| `src/cartesian_terminal_shellification_geometry.jl` | `src/cartesian/cartesian_terminal_shellification_geometry.jl` | `ff1ba22edc86b1b594c90a75b5157ba86bcbf073` | `a3c8ba89840acf58df4cda47885f5796dba1c5bbcf90cc4933d5d78866ed73a5` |
| `src/gaussian_coulomb_reference.jl` | `src/cartesian/gaussian_coulomb_reference.jl` | `a315467fd6de5458ffade5b0c49f17712e46b7d4` | `ed600bfda69527c925710a82abec82e6ff4cc47d94f7f235bfcc2b42cdac94fc` |
| `src/hamiltonian_corrections.jl` | `src/cartesian/hamiltonian_corrections.jl` | `c58854b7156060e776a1fa4b30bfe7118cc36ef5` | `3cdde68ebb068c30381e1a3ab83ea2a3b1af7bec8c13692c38da24f54ed1ffaa` |
| `src/pqs_matched_h2plus.jl` | `src/cartesian/pqs_matched_h2plus.jl` | `f7fc4e7663993e23dd6e4bb842219ab83e004370` | `322d7fbe3eb7a0050cb045210a8bdb786f18f0b7d511ed26cdda71b17f572c8d` |
| `src/pqs_source_box_diatomic_complete_core_shell.jl` | `src/cartesian/pqs_source_box_diatomic_complete_core_shell.jl` | `b5838e49703ae81702e4a8a12135d17c6d51291f` | `5b118a3f9ca218df7ed7bf3f17fc283c2a1b2ce448bf8b3d4f964d870f518bf6` |
| `src/pqs_source_box_low_order_materialization.jl` | `src/cartesian/pqs_source_box_low_order_materialization.jl` | `dd153e71745469f2f26ffbcaaf37cd8be509ed0f` | `bd7b60c06d681a870ec646dbb06706e86fe6f70da96eb8389e1b8207dc781379` |
| `src/pqs_source_box_route_driver_helpers.jl` | `src/cartesian/pqs_source_box_route_driver_helpers.jl` | `69f83fb552aacd03c7488efc3b5e6d73414239c2` | `b2f68d09fff5e4622173597f915d3278b7eafc07d090ea2fdd5cd8096be973f8` |
| `src/pqs_source_box_route_driver_skeletons.jl` | `src/cartesian/pqs_source_box_route_driver_skeletons.jl` | `aa1644ebc46574508b8c40f583e9b56fd148856e` | `6dc13c1790c4e33ba04767ee6095e5aafb17743cd232a100da20eccd3a764a3b` |

In `src/GaussletBases.jl`, only the corresponding `39` quoted include paths
may change, in place. `CartesianParentAxisFactors.jl` is the sole mapped flat
owner not included there; `CartesianParentGaussletBases.jl` includes it as one
of the `37` byte-identical nested sibling includes. The complete `88`-entry
include order is fixed; after removing the new `cartesian/` prefix its SHA-256
must remain
`8ab555c688804b36f384f0516d75407e618b9fe78b875060a685e0aa454f5f7e`.
The `37` nested sibling includes and `17` module-relative imports in the moved
owners remain byte-identical. None of these owners uses `@__DIR__`, `@__FILE__`,
the package data helper, or a source-relative resource path, and no source,
test, example, script, tool, or workflow outside the root loader directly
includes one of them.

That 88-entry digest remains the exact Step 4B relocation identity. The
accepted `HP-FOUNDATION-LANCZOS-FN-01` transaction inserted only
`foundation/lanczos.jl` immediately after `foundation/primitive_sets.jl`,
without changing the relative order of any existing entry. Its normalized
89-entry digest is
`46cc0f341af7c41663359b317e780b4b9fa7bc2e2a61292b0d460762af588016`;
this does not reopen any Cartesian/PQS move or body.

Current-pointer reconciliation is selective and mechanical. Update `340`
affected `records.paths` entries, including the one still-planned represented-
Hartree path, the `7` existing authority scope path literals, and the `3`
current `repo_path` evidence literals owned by `HP-RHO0-FAPP-FN-01` and
`HP-RHO0-JANCHOR-FN-01`. Update the `264` current code-pointer occurrences
across the `53` hand-edited current-pointer documents identified by the
preflight. The additional listed canonical map remains frozen. Regenerate
`registry.md` and `execution_whitelist.md` as whole files from `authority.toml`.
The pre-move column above, `324` exact-path occurrences in `23` archived or
historical documents, and broader retired-path evidence remain frozen. Do not
perform a repository-wide replacement. The anti-bloat policy in `AGENTS.md`
was made path-neutral in the authorizing commit and requires no implementation
edit.

Implementation target card:

```text
Target: relocate the exact 76-file Cartesian/PQS inventory under src/cartesian/
Physics endpoint: none changed; existing Cartesian, PQS, and Screening owners remain the oracles
Allowed files: one path-neutral docstring commit; 28 file and 12 directory git mv operations; 39 root include substitutions; classified current authority/documentation paths
Forbidden files/surfaces: every executable body, test, workflow, API, numerical contract, release, and follow-on refactor
Must delete or simplify: remove the 28 flat files and 12 flat submodule directories; delete no implementation
Forbidden additions: modules, namespaces, renamed bindings/files, wrappers, shims, helpers, caches, fixtures, dependencies, tests, or files beyond the destination directory
Success condition: exact prerequisite hash, then 76 R100 files, exact path substitutions, unchanged include order and behavior
Validation: hashes/renames, include/import order, path closure, package load, bounded public and maintenance owners, docs/authority/Documenter/diff, normal full CI/Docs
Line budget: prerequisite docstring +1/-1; moved bodies +0/-0; src/GaussletBases.jl +39/-39; current authority/documentation substitutions only; tests/workflows +0/-0
Failure rule: make no relocation commit if any hash, rename score, include/import position, executable body, test, path count, package load, or unchanged owner differs
```

Local acceptance runs package load; the existing
`core,ida,cartesian,examples,docs_fast` selection; `pqs_release`;
`screening_release`; and the four scheduled Cartesian maintenance owners once,
each in its own Julia process. It also runs full docs tests, authority check and
self-test, two independent generated-view renders, Documenter, manager-log
bound, and diff checks. The source-bearing push must classify as full and pass
the unchanged Supported floor, PQS paper, and Screening paper jobs plus the
separate Docs workflow. Do not run the full angular owner, occupied-first,
represented-Hartree, or unrelated nested suites.

This authority adds no umbrella module, namespace, binding, export, dispatch,
dependency, test policy, workflow, numerical policy, cache, or release work.
`cartesian_nested_faces.jl` decomposition, further Lanczos work, Foundation-to-
Ordinary ownership changes, orphan-test repair, conditioning work, empty-
directory cleanup, and every later layout or architecture change remain
separate.
