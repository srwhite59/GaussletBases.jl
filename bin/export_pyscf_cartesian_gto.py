#!/usr/bin/env python3
"""Export one attested standard-Hamiltonian PySCF HF checkpoint as XGTO v1."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pyscf
from pyscf import gto, scf

KIND = "gaussletbases.external_cartesian_gto"
VERSION = 1
EXPORTER_VERSION = "1.0.0"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def quote(value: str) -> str:
    return json.dumps(str(value), ensure_ascii=True)


def atom_array(values) -> str:
    def one(value):
        if isinstance(value, str):
            return quote(value)
        if isinstance(value, (bool, np.bool_)):
            return "true" if value else "false"
        if isinstance(value, (int, np.integer)):
            return str(int(value))
        value = float(value)
        if not np.isfinite(value):
            raise ValueError("manifest values must be finite")
        return repr(value)
    return "[" + ", ".join(one(x) for x in values) + "]"


def cartesian_powers(l: int):
    return [(lx, ly, l - lx - ly)
            for lx in range(l, -1, -1)
            for ly in range(l - lx, -1, -1)]


def cartesian_molecule(mol: gto.Mole) -> gto.Mole:
    cart = gto.loads(mol.dumps())
    cart.cart = True
    cart.symmetry = False
    cart.output = None
    cart.verbose = 0
    cart.build()
    return cart


def checkpoint_state(path: Path):
    mol, state = scf.chkfile.load_scf(str(path))
    if not isinstance(mol, gto.Mole) or hasattr(mol, "a"):
        raise ValueError("only finite molecular PySCF checkpoints are supported")
    if any(mol.atom_nelec_core(i) for i in range(mol.natm)) or getattr(mol, "pseudo", None):
        raise ValueError("ECP and pseudopotential checkpoints are unsupported")
    coeff = np.asarray(state["mo_coeff"])
    occ = np.asarray(state["mo_occ"])
    if np.iscomplexobj(coeff) or not np.isrealobj(coeff):
        raise ValueError("complex, spinor, and GHF coefficients are unsupported")
    if not np.isfinite(coeff).all() or not np.isfinite(occ).all():
        raise ValueError("checkpoint orbitals and occupations must be finite")
    if np.any(occ < 0):
        raise ValueError("checkpoint occupations must be nonnegative")
    if coeff.ndim == 2 and occ.ndim == 1:
        if mol.spin != 0:
            raise ValueError("single-block checkpoints must be spin-zero RHF; ROHF is rejected")
        blocks = [("restricted", coeff, occ)]
        state_kind = "rhf"
    elif coeff.ndim == 3 and coeff.shape[0] == 2 and occ.shape[0] == 2:
        blocks = [("alpha", coeff[0], occ[0])]
        if np.any(occ[1] > 0):
            blocks.append(("beta", coeff[1], occ[1]))
            state_kind = "uhf"
        else:
            state_kind = "alpha_only"
    else:
        raise ValueError("checkpoint is not real RHF or collinear UHF")
    expected = ((mol.nelectron,) if state_kind == "rhf" else
                (mol.nelec[0],) if state_kind == "alpha_only" else tuple(mol.nelec))
    actual = tuple(float(np.sum(block[2])) for block in blocks)
    if len(actual) != len(expected) or not np.allclose(actual, expected, atol=1e-10, rtol=0):
        raise ValueError(f"checkpoint occupations {actual} disagree with electron sector {expected}")
    return mol, state_kind, blocks


def explicit_aos(mol: gto.Mole):
    labels = gto.cart_labels(mol, fmt=True)
    records = []
    row = 0
    for ib in range(mol.nbas):
        atom = int(mol.bas_atom(ib))
        l = int(mol.bas_angular(ib))
        exponents = np.asarray(mol.bas_exp(ib), dtype=float)
        coefficients = np.asarray(mol.bas_ctr_coeff(ib), dtype=float)
        for ictr in range(int(mol.bas_nctr(ib))):
            for powers in cartesian_powers(l):
                records.append({
                    "ao_index_1based": row + 1,
                    "label": labels[row],
                    "atom_index_1based": atom + 1,
                    "shell_index_1based": ib + 1,
                    "contraction_index_1based": ictr + 1,
                    "angular_powers": powers,
                    "center_bohr": mol.atom_coord(atom, unit="Bohr"),
                    "exponents": exponents,
                    "coefficients": coefficients[:, ictr],
                })
                row += 1
    if row != mol.nao_nr():
        raise RuntimeError(f"Cartesian AO inventory mismatch: {row} != {mol.nao_nr()}")
    return records


def occupied_blocks(mol, cart, source_kind, source_blocks, overlap):
    X = None if source_kind == "cartesian" else mol.cart2sph_coeff(normalized="sp")
    result = []
    for spin, coeff, occupations in source_blocks:
        keep = np.flatnonzero(occupations > 0)
        if not len(keep):
            raise ValueError(f"{spin} block has no occupied orbitals")
        C = np.asarray(coeff[:, keep], dtype=float)
        if X is not None:
            C = X @ C
        if not np.isfinite(C).all():
            raise ValueError(f"{spin} Cartesian occupied coefficients are not finite")
        error = np.linalg.norm(C.T @ overlap @ C - np.eye(len(keep)), ord=np.inf)
        if not np.isfinite(error) or error > 1e-10:
            raise ValueError(f"{spin} occupied metric error {error:.3e} exceeds 1e-10")
        result.append({"spin": spin, "coefficients": C,
                       "occupations": np.asarray(occupations[keep], dtype=float),
                       "mo_indices": keep + 1, "metric_error": float(error)})
    return result


def payload_arrays(overlap, blocks):
    arrays = [("source_overlap", overlap)]
    for block in blocks:
        prefix = "alpha" if block["spin"] in ("restricted", "alpha") else "beta"
        arrays.extend([(prefix + "_coefficients", block["coefficients"]),
                       (prefix + "_occupations", block["occupations"])])
    payload = bytearray()
    records = []
    for name, value in arrays:
        array = np.asarray(value, dtype="<f8", order="F")
        raw = array.tobytes(order="F")
        records.append({"name": name, "shape": array.shape,
                        "byte_offset": len(payload), "element_count": array.size,
                        "sha256": sha256(raw)})
        payload.extend(raw)
    return bytes(payload), records


def write_manifest(path, payload_name, mol, cart, state_kind, source_kind,
                   aos, blocks, arrays, payload_hash, checkpoint_hash):
    conversion = ("mol.cart2sph_coeff(normalized=sp)"
                  if source_kind == "spherical" else "not_applied")
    lines = [
        "[format]", f"kind = {quote(KIND)}", f"version = {VERSION}", "",
        "[encoding]", 'scalar_type = "ieee754_float64"', 'byte_order = "little_endian"',
        'array_order = "column_major"', "",
        "[units]", 'length = "bohr"', 'exponent = "bohr_inverse_square"', "",
        "[producer]", 'name = "PySCF"', f"version = {quote(pyscf.__version__)}",
        f"python_version = {quote(platform.python_version())}",
        f"numpy_version = {quote(np.__version__)}",
        'exporter = "export_pyscf_cartesian_gto.py"',
        f"exporter_version = {quote(EXPORTER_VERSION)}", "",
        "[state]", f"kind = {quote(state_kind)}",
        f"source_ao_representation = {quote(source_kind)}",
        'exported_ao_representation = "cartesian"',
        'cartesian_normalization = "pyscf_libcint_sp"',
        f"spherical_conversion = {quote(conversion)}",
        'hamiltonian_kind = "all_electron_nuclear_coulomb"',
        f"molecular_charge = {int(mol.charge)}", f"spin_2s = {int(mol.spin)}",
        f"electron_count = {int(mol.nelectron)}", f"nalpha = {int(mol.nelec[0])}",
        f"nbeta = {int(mol.nelec[1])}", f"ao_count = {cart.nao_nr()}", "",
        "[source]", f"checkpoint_sha256 = {quote(checkpoint_hash)}", "",
        "[payload]", f"basename = {quote(payload_name)}",
        f"byte_count = {sum(8 * x['element_count'] for x in arrays)}",
        f"sha256 = {quote(payload_hash)}", "",
    ]
    for ia in range(mol.natm):
        lines += ["[[nuclei]]", f"atom_index_1based = {ia+1}",
                  f"symbol = {quote(mol.atom_symbol(ia))}",
                  f"nuclear_charge = {repr(float(mol.atom_charge(ia)))}",
                  f"center_bohr = {atom_array(mol.atom_coord(ia, unit='Bohr'))}", ""]
    for ao in aos:
        lines += ["[[aos]]", f"ao_index_1based = {ao['ao_index_1based']}",
                  f"label = {quote(ao['label'])}",
                  f"atom_index_1based = {ao['atom_index_1based']}",
                  f"shell_index_1based = {ao['shell_index_1based']}",
                  f"contraction_index_1based = {ao['contraction_index_1based']}",
                  f"angular_powers = {atom_array(ao['angular_powers'])}",
                  f"center_bohr = {atom_array(ao['center_bohr'])}",
                  f"exponents_bohr_inverse_square = {atom_array(ao['exponents'])}",
                  f"contraction_coefficients = {atom_array(ao['coefficients'])}", ""]
    for block in blocks:
        lines += ["[[orbital_blocks]]", f"spin = {quote(block['spin'])}",
                  f"mo_indices_1based = {atom_array(block['mo_indices'])}",
                  f"metric_error_inf = {repr(block['metric_error'])}", ""]
    for array in arrays:
        lines += ["[[arrays]]", f"name = {quote(array['name'])}",
                  f"shape = {atom_array(array['shape'])}",
                  f"byte_offset = {array['byte_offset']}",
                  f"element_count = {array['element_count']}",
                  f"sha256 = {quote(array['sha256'])}", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=(
        "Export occupied orbitals from an existing real molecular PySCF HF checkpoint. "
        "Invoking this command attests that the checkpoint used the ordinary all-electron "
        "nuclear-Coulomb Hamiltonian; checkpoint fields cannot prove that condition."))
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.suffix != ".toml":
        raise ValueError("output must have .toml extension")
    payload_path = args.output.with_suffix(".f64")
    if args.output.exists() or payload_path.exists():
        raise FileExistsError("refusing to overwrite output bundle")
    mol, state_kind, source_blocks = checkpoint_state(args.checkpoint)
    source_kind = "cartesian" if mol.cart else "spherical"
    cart = cartesian_molecule(mol)
    overlap = cart.intor_symmetric("int1e_ovlp")
    aos = explicit_aos(cart)
    blocks = occupied_blocks(mol, cart, source_kind, source_blocks, overlap)
    payload, arrays = payload_arrays(overlap, blocks)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload_path.write_bytes(payload)
    write_manifest(args.output, payload_path.name, mol, cart, state_kind, source_kind,
                   aos, blocks, arrays, sha256(payload), sha256(args.checkpoint.read_bytes()))


if __name__ == "__main__":
    main()
