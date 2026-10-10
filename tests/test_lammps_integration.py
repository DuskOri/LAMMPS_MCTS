"""Optional execution checks; set LAMMPS_TEST_EXE and LAMMPS_TEST_DATA."""

from pathlib import Path
from tempfile import TemporaryDirectory
import os
import re
import unittest

from generator.polymer_builder import build_polymer_from_sequence
from generator.packmol_runner import prepare_initial_system
from md_engine.lammps_runner import run_lammps_input
from md_engine.process_control import reset_cancellation
from md_engine.thermal_conductivity_writer import write_direct_nemd_input, write_adaptive_nemd_input, _render_lammps_header, _normalise_params
from mcts import configure_fragment_registry


@unittest.skipUnless(os.environ.get("LAMMPS_TEST_EXE") and os.environ.get("LAMMPS_TEST_DATA"), "LAMMPS integration paths not configured")
class LammpsExecutionTests(unittest.TestCase):
    def setUp(self):
        reset_cancellation()
        self.temp = TemporaryDirectory()
        self.root = Path(os.environ.get("LAMMPS_TEST_OUTPUT", self.temp.name))
        self.root.mkdir(parents=True, exist_ok=True)
        self.executable = os.environ["LAMMPS_TEST_EXE"]
        self.data = Path(os.environ["LAMMPS_TEST_DATA"])
        self.params = dict(density_equilibration=False, nvt_steps=1, npt_steps=1, timestep=0.1,
                           nemd_steady_steps=2, production_steps=30, nemd_production_block_steps=10,
                           nemd_min_production_blocks=3, nemd_profile_nevery=1, nemd_profile_repeat=1,
                           dump_every=10, thermo_every=10, nemd_heat_flux_tolerance=1e-6)

    def tearDown(self):
        self.temp.cleanup()
        reset_cancellation()

    def run_script(self, name):
        return run_lammps_input(self.root / (name + ".in"), self.executable,
                                log_path=self.root / (name + ".log"), timeout=60,
                                work_dir=self.root, stream_output=False)

    def check_measurement(self, writer, name):
        result = writer(self.data, self.root / (name + ".in"), self.root / name, self.params)
        self.assertTrue(result.success, result.message)
        run = self.run_script(name)
        log = Path(run.log_path).read_text(encoding="utf-8")
        self.assertNotIn("ERROR", log)
        self.assertIn("Heat flux convergence check", log)
        self.assertIn(run.returncode, (0, 2))
        if "WARNING: Heat flux did not converge" in log:
            self.assertFalse(run.success)
            self.assertEqual(run.returncode, 2)
        for suffix in ("_heat_flux.dat", "_temperature_profile.dat", "_bath_energy.dat", "_final.data"):
            self.assertTrue((self.root / (name + suffix)).exists(), suffix)

    def test_adaptive_template_executes_and_rejects_nonconvergence(self):
        self.check_measurement(write_adaptive_nemd_input, "adaptive")

    def test_combined_template_has_all_required_fixes(self):
        self.check_measurement(write_direct_nemd_input, "combined")

    def test_flux_convergence_success_path_executes_cleanup(self):
        self.params["nemd_heat_flux_tolerance"] = 1.0
        write_adaptive_nemd_input(self.data, self.root / "converged.in", self.root / "converged", self.params)
        result = self.run_script("converged")
        self.assertTrue(result.success, result.message)
        self.assertIn("Heat flux converged;", Path(result.log_path).read_text(encoding="utf-8"))

    def test_per_atom_flux_matches_lammps_and_mask_updates(self):
        header = _render_lammps_header(self.data.as_posix(), _normalise_params(self.params))
        script = header + '''
velocity all create 300.0 87287 mom yes rot yes dist gaussian
region flux block EDGE EDGE EDGE EDGE EDGE EDGE units box
compute myKE all ke/atom
compute myPE all pe/atom
compute stressPair all stress/atom NULL pair kspace
compute stressBond all centroid/stress/atom NULL bond angle dihedral improper
compute fluxPair all heat/flux myKE myPE stressPair
compute fluxBond all heat/flux myKE myPE stressBond
variable local_Jx atom "rmask(flux)*((c_myKE+c_myPE)*vx-((c_stressPair[1]+c_stressBond[1])*vx+(c_stressPair[4]+c_stressBond[4])*vy+(c_stressPair[5]+c_stressBond[5])*vz)/68568.415)"
compute regional_flux all reduce sum v_local_Jx
variable reference equal c_fluxPair[1]+c_fluxBond[1]-c_fluxPair[4]
variable difference equal c_regional_flux-v_reference
variable occupants atom rmask(flux)
compute selected_count all reduce sum v_occupants
thermo_style custom step v_difference c_selected_count
run 0
print "FLUX_DIFFERENCE=$(v_difference:%.12g)"
region flux delete
region flux block $(xlo) $(0.5*(xlo+xhi)) EDGE EDGE EDGE EDGE units box
group probe id 1
set group probe x $(xlo+0.1*lx)
run 0
print "COUNT_BEFORE=$(c_selected_count:%.0f)"
set group probe x $(xlo+0.9*lx)
run 0
print "COUNT_AFTER=$(c_selected_count:%.0f)"
'''
        (self.root / "flux_identity.in").write_text(script, encoding="utf-8")
        result = self.run_script("flux_identity")
        self.assertTrue(result.success, result.message)
        log = Path(result.log_path).read_text(encoding="utf-8")
        difference = float(re.search(r"^FLUX_DIFFERENCE=([^\s]+)", log, re.M)[1])
        self.assertLess(abs(difference), 1e-10)
        before = int(re.search(r"^COUNT_BEFORE=(\d+)", log, re.M)[1])
        after = int(re.search(r"^COUNT_AFTER=(\d+)", log, re.M)[1])
        self.assertEqual(before - after, 1)


@unittest.skipUnless(os.environ.get("PACKMOL_TEST_EXE"), "Packmol integration path not configured")
class PackmolExecutionTests(unittest.TestCase):
    def test_real_packmol_and_uff_paths_with_spaces(self):
        reset_cancellation()
        configure_fragment_registry({"active_families": "polyimide"})
        with TemporaryDirectory(prefix="packmol space ") as tmp:
            root = Path(tmp)
            build = build_polymer_from_sequence(["A", "B", "End"], dp=1, output_dir=root)
            self.assertTrue(build.success, build.message)
            packed = prepare_initial_system(build.pdb_path, root / "packed", molecule_count=2,
                box_size=30, packmol_executable=os.environ["PACKMOL_TEST_EXE"], template_mol=build.mol_path, timeout=30)
            self.assertTrue(packed.success, packed.message)
            self.assertTrue(packed.force_field_ready)
            self.assertTrue(Path(packed.system_data_path).exists())


if __name__ == "__main__":
    unittest.main()
