"""Tests for the workflow layer.

``test_plan_reproduces_the_example`` checks that, from the cif2cell output and the TOML of every reference example, the new
implementation alone produces all nine input files, byte for byte.
"""

import sys
import types
import warnings

import numpy as np
import pytest

from cif2qewan.config import Config
from cif2qewan.exceptions import (
    ConfigError,
    InputModelError,
    PseudopotentialError,
    StructureError,
)
from cif2qewan.qe import kpoints
from cif2qewan.qe.model import KPathPoint, KPointsAutomatic, KPointsPath, RawValue
from cif2qewan.qe.pseudopotential import PseudopotentialEntry, PseudopotentialTable
from cif2qewan.qe.writer import render_namelist_input
from cif2qewan.structure.compare import compare_structures
from cif2qewan.structure.model import AtomicSite, NormalizedStructure
from cif2qewan.structure.readers import Cif2cellReader, PymatgenReader
from cif2qewan.workflow.builder import (
    STARTING_MAGNETIZATION,
    StructureHints,
    WannierCounts,
    WorkflowBuilder,
    build_plan,
    plan_from_cif2cell_output,
)
from cif2qewan.workflow.render import render_plan
from conftest import EXAMPLE_CASES, PACKAGE, example_config

FE = EXAMPLE_CASES[
    0
].reference  # PSLibrary/Fe: cif_scf.in and CIF shared by all Fe examples
PSL = PACKAGE / "pp_psl_rrkj.csv"


# --------------------------------------------------------------------------
# Config
# --------------------------------------------------------------------------


def test_config_from_the_example_toml():
    config = Config.from_toml(FE / "cif2qewan.toml", so=True)
    assert config.pseudo_dir == "/path/to/pslibrary"
    assert config.scf_k_resolution == pytest.approx(0.15)
    assert config.degauss == pytest.approx(0.01)
    assert config.pw2wan == {"write_unk": ".true."}
    assert config.so and not config.mag and config.spinor


def test_config_errors(tmp_path):
    with pytest.raises(ConfigError, match="not found"):
        Config.from_toml(tmp_path / "missing.toml")
    with pytest.raises(ConfigError, match="missing configuration keys"):
        Config.from_dict({"pseudo_dir": "x"})
    base = {
        "cif2cell_path": "c",
        "pseudo_dir": "p",
        "pp_list_path": "l",
        "scf_k_resolution": 0.15,
        "degauss": 0.01,
    }
    with pytest.raises(ConfigError, match=r"\[pw2wan\]"):
        Config.from_dict(base)
    with pytest.raises(ConfigError, match="write_unk"):
        Config.from_dict({**base, "pw2wan": {}})
    with pytest.raises(ConfigError, match="degauss"):
        Config.from_dict({**base, "degauss": 0.0, "pw2wan": {"write_unk": ".true."}})
    with pytest.raises(ConfigError, match="invalid configuration value"):
        Config.from_dict({**base, "degauss": "abc", "pw2wan": {"write_unk": ".true."}})
    with pytest.raises(ConfigError, match="pw2wan.write_unk"):
        Config.from_dict({**base, "pw2wan": {"write_unk": "sometimes"}})
    bad = tmp_path / "bad.toml"
    bad.write_text("this = [is not toml\n")
    with pytest.raises(ConfigError, match="cannot read"):
        Config.from_toml(bad)


# --------------------------------------------------------------------------
# Pseudopotential table
# --------------------------------------------------------------------------


def test_table_lookup_and_relativistic_names():
    table = PseudopotentialTable.from_csv(PSL)
    fe = table.lookup("Fe")
    assert fe == PseudopotentialEntry(
        "Fe", "Fe.pbe-spn-rrkjus_psl.0.2.1.UPF", 4, "spd", 64.0, 782.0
    )
    assert fe.num_wann == 9
    assert fe.projection_orbitals == ("s", "p", "d")
    assert fe.pseudo_file(relativistic=True) == "Fe.rel-pbe-spn-rrkjus_psl.0.2.1.UPF"
    assert fe.pseudo_file(relativistic=False) == fe.file_name
    assert "Fe" in table and "Fr" not in table

    dojo = PseudopotentialTable.from_csv(PACKAGE / "nc-sr-05_pbe_standard_upf.csv")
    assert dojo.lookup("Fe").pseudo_file(relativistic=True) == "Fe_fr.UPF"
    assert dojo.lookup("He").num_wann == 0  # empty orbitals: no projection


def test_table_errors(tmp_path):
    table = PseudopotentialTable.from_csv(PSL)
    with pytest.raises(PseudopotentialError, match="no pseudopotential"):
        table.lookup("Fr")
    with pytest.raises(PseudopotentialError, match="not listed"):
        table.lookup("Xx")
    with pytest.raises(PseudopotentialError, match="not found"):
        PseudopotentialTable.from_csv(tmp_path / "missing.csv")
    wrong = tmp_path / "wrong.csv"
    wrong.write_text("atom,file\nFe,x\n")
    with pytest.raises(PseudopotentialError, match="expected columns"):
        PseudopotentialTable.from_csv(wrong)
    with pytest.raises(PseudopotentialError, match="unknown orbital"):
        PseudopotentialEntry("Fe", "Fe.UPF", 0, "spdg", 1.0, 4.0)


# --------------------------------------------------------------------------
# k-point policies
# --------------------------------------------------------------------------


def fe_structure():
    return Cif2cellReader.read_output(FE / "cif_scf.in").structure


def test_meshes():
    assert kpoints.scf_mesh(fe_structure(), 0.15) == (21, 21, 21)
    assert kpoints.scf_mesh(fe_structure(), 0.30) == (10, 10, 10)
    assert kpoints.nscf_mesh((21, 21, 21)) == (8, 8, 8)
    assert kpoints.nscf_mesh((2, 6, 40)) == (4, 6, 8)
    points = kpoints.mesh_points((2, 1, 3))
    assert len(points) == 6
    assert points[:3] == ((0.0, 0.0, 0.0), (0.0, 0.0, 1 / 3), (0.0, 0.0, 2 / 3))
    listed = kpoints.mesh_kpoints_list((2, 2, 2))
    assert listed.option == "crystal" and listed.num_points == 8
    assert sum(p[3] for p in listed.points) == pytest.approx(1.0)


def test_band_path_of_fe_matches_the_reference():
    path, used_seekpath = kpoints.band_path(fe_structure())
    assert used_seekpath
    assert [(p.label, p.ndiv) for p in path.points] == [
        ("G", 20),
        ("H", 14),
        ("N", 14),
        ("G", 17),
        ("P", 17),
        ("H", 0),
        ("P", 10),
        ("N", 10),
    ]
    segments = kpoints.win_kpoint_path(path)
    assert [(s.start_label, s.end_label) for s in segments] == [
        ("G", "H"),
        ("H", "N"),
        ("N", "G"),
        ("G", "P"),
        ("P", "H"),
        ("P", "N"),
    ]


def test_band_path_falls_back_without_seekpath(monkeypatch):
    monkeypatch.setattr(kpoints, "_seekpath", lambda: None)
    path, used_seekpath = kpoints.band_path(fe_structure())
    assert not used_seekpath
    assert [p.label for p in path.points] == ["R", "G", "X", "M", "G"]
    assert {p.ndiv for p in path.points} == {20}


def test_old_seekpath_is_rejected(monkeypatch):
    fake = types.ModuleType("seekpath")
    fake.__version__ = "2.0.5"
    monkeypatch.setitem(sys.modules, "seekpath", fake)
    assert kpoints._seekpath() is None
    fake.__version__ = "2.2.1"
    assert kpoints._seekpath() is fake


def test_win_kpoint_path_skips_jumps():
    path = KPointsPath(
        "crystal_b",
        (
            KPathPoint((0, 0, 0), 5, "G"),
            KPathPoint((0.5, 0, 0), 0, "X"),
            KPathPoint((0.5, 0.5, 0), 5, "M"),
        ),
    )
    assert [(s.start_label, s.end_label) for s in kpoints.win_kpoint_path(path)] == [
        ("G", "X")
    ]


# --------------------------------------------------------------------------
# WorkflowBuilder: the reference examples
# --------------------------------------------------------------------------


@pytest.fixture(params=EXAMPLE_CASES, ids=lambda c: c.example)
def case(request):
    return request.param


def test_plan_reproduces_the_example(case):
    output = Cif2cellReader.read_output(case.reference / "cif_scf.in")
    plan = plan_from_cif2cell_output(output, example_config(case))
    files = render_plan(plan)

    assert set(files) == {
        "scf.in",
        "nscf.in",
        "pw2wan.in",
        "pwscf.win",
        "check_wannier/nscf.in",
        "band/nscf.in",
        "band/band.in",
        "band/proj.in",
        "band/pp.in",
    }
    for path, text in files.items():
        assert text == (case.reference / path).read_text(), f"{case.example}/{path}"


def test_pymatgen_path_gives_an_equivalent_plan(case):
    """Without cif2cell hints the plan differs only in the cell representation."""
    pytest.importorskip("pymatgen")
    config = example_config(case)
    from_cif2cell = plan_from_cif2cell_output(
        Cif2cellReader.read_output(case.reference / "cif_scf.in"), config
    )
    structure = PymatgenReader().read(case.reference / "mp-13_Fe.cif")
    from_pymatgen = build_plan(structure, config)

    assert compare_structures(from_cif2cell.structure, from_pymatgen.structure) == []
    assert (
        from_pymatgen.scf.kpoints == from_cif2cell.scf.kpoints
    )  # 21 21 21 from the resolution
    for name in ("scf", "nscf", "check_wannier", "bands_nscf"):
        a, b = getattr(from_cif2cell, name), getattr(from_pymatgen, name)
        assert a.control == b.control and a.electrons == b.electrons
        assert a.species == b.species  # same masses and pseudopotential files
        assert {k: v for k, v in a.system.entries.items() if k != "A"} == {
            k: v for k, v in b.system.entries.items() if k != "A"
        }
        volume_a = (
            abs(np.linalg.det(a.cell.matrix()))
            * float(a.system.entries["A"].literal) ** 3
        )
        volume_b = (
            abs(np.linalg.det(b.cell.matrix()))
            * float(b.system.entries["A"].literal) ** 3
        )
        assert volume_a == pytest.approx(volume_b, rel=1e-4)
    w_a, w_b = from_cif2cell.wannier90, from_pymatgen.wannier90
    assert (
        w_a.num_wann,
        w_a.num_bands,
        w_a.exclude_bands,
        w_a.spinors,
        w_a.mp_grid,
    ) == (w_b.num_wann, w_b.num_bands, w_b.exclude_bands, w_b.spinors, w_b.mp_grid)
    assert w_a.projections == w_b.projections
    assert (
        render_plan(from_pymatgen)["pw2wan.in"]
        == render_plan(from_cif2cell)["pw2wan.in"]
    )


# --------------------------------------------------------------------------
# WorkflowBuilder: policies
# --------------------------------------------------------------------------


def test_wannier_counts():
    counts = WannierCounts(num_wann=9, nexclude=4, spinor=False)
    assert (counts.num_wann_total, counts.nexclude_total, counts.num_bands) == (
        9,
        4,
        27,
    )
    assert (counts.nbnd_nscf, counts.nbnd_check) == (31, 17)
    spinor = WannierCounts(num_wann=9, nexclude=4, spinor=True)
    assert (spinor.num_wann_total, spinor.nexclude_total, spinor.num_bands) == (
        18,
        8,
        54,
    )
    assert (spinor.nbnd_nscf, spinor.nbnd_check) == (62, 34)


def test_mag_policy_is_collinear_scf_then_noncollinear_nscf():
    case = EXAMPLE_CASES[0]
    plan = build_plan(fe_structure(), example_config(case, so=False, mag=True))
    scf, nscf = plan.scf.system.entries, plan.nscf.system.entries
    assert (
        scf["nspin"] == 2 and scf["starting_magnetization(1)"] == STARTING_MAGNETIZATION
    )
    assert "noncolin" not in scf
    assert nscf["noncolin"] is True and nscf["lforcet"] is True
    assert nscf["lspinorb"] is False and (nscf["angle1"], nscf["angle2"]) == (0, 0)
    assert nscf["starting_magnetization(1)"] == STARTING_MAGNETIZATION
    assert (
        plan.scf.species[0].pseudo_file == plan.nscf.species[0].pseudo_file
    )  # no --so: scalar rel.
    assert plan.wannier90.spinors


def test_cutoff_warning_and_missing_projection(tmp_path):
    table = tmp_path / "table.csv"
    table.write_text(
        "atom,pp_file_name,nexclude,orbitals,ecutwfc,ecutrho\n"
        "Fe,Fe.pbe,4,,60.0,200.0\n"
    )
    config = Config.from_dict(
        {
            "cif2cell_path": "c",
            "pseudo_dir": "p",
            "pp_list_path": str(table),
            "scf_k_resolution": 0.15,
            "degauss": 0.01,
            "pw2wan": {"write_unk": ".false.", "wannier_plot_supercell": 3},
        }
    )
    builder = WorkflowBuilder(config)
    with pytest.raises(InputModelError, match="no Wannier projections"):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            builder.build(fe_structure())
    assert any("ecut_rho should be bigger" in str(w.message) for w in caught)

    text = render_namelist_input(builder._pw2wan_input())
    assert " wannier_plot_supercell = 3\n" in text
    assert " write_unk = .false.\n" in text


@pytest.mark.parametrize(
    ("write_unk", "enabled"),
    [(".true.", True), (".false.", False), (True, True), (False, False)],
)
def test_wannier_plot_requires_unk_files(write_unk, enabled):
    config = Config.from_dict(
        {
            "cif2cell_path": "c",
            "pseudo_dir": "/pp",
            "pp_list_path": str(PSL),
            "scf_k_resolution": 0.15,
            "degauss": 0.01,
            "pw2wan": {"write_unk": write_unk},
        }
    )
    plan = build_plan(fe_structure(), config)

    assert config.write_unk is enabled
    assert plan.wannier90.parameters["wannier_plot"] is enabled
    assert (
        f" write_unk = {'.true.' if enabled else '.false.'}\n"
        in render_namelist_input(plan.pw2wan)
    )


def test_hints_must_match_the_structure():
    builder = WorkflowBuilder(example_config(EXAMPLE_CASES[1]))
    with pytest.raises(InputModelError, match="do not match"):
        builder.build(fe_structure(), StructureHints(masses={"Co": 58.9}))
    plan = builder.build(fe_structure(), StructureHints(alat=1.0, scf_kmesh=(3, 3, 3)))
    assert plan.scf.system.entries["A"] == RawValue(f"{1.0:10.5f}")
    assert plan.scf.kpoints == KPointsAutomatic((3, 3, 3))
    assert plan.wannier90.mp_grid == (4, 4, 4)  # clamped nscf mesh


def test_builder_rejects_partial_occupancy():
    partial = NormalizedStructure(
        np.eye(3) * 4.0,
        (AtomicSite("Na", (0, 0, 0), occupancy=0.5),),
    )
    with pytest.raises(StructureError, match="partial occupancy"):
        build_plan(partial, example_config(EXAMPLE_CASES[1]))


def test_two_species_use_the_table_order_and_sum_the_counts():
    structure = NormalizedStructure(
        np.eye(3) * 4.0, (AtomicSite("O", (0.5, 0.5, 0.5)), AtomicSite("Fe", (0, 0, 0)))
    )
    plan = build_plan(structure, example_config(EXAMPLE_CASES[1]))
    assert [s.label for s in plan.scf.species] == ["O", "Fe"]
    fe, o = (PseudopotentialTable.from_csv(PSL).lookup(el) for el in ("Fe", "O"))
    assert plan.wannier90.num_wann == fe.num_wann + o.num_wann
    assert (
        plan.wannier90.exclude_bands == (1, fe.nexclude + o.nexclude)
        if fe.nexclude + o.nexclude
        else None
    )
    assert plan.scf.system.entries["ecutwfc"] == max(fe.ecutwfc, o.ecutwfc)
