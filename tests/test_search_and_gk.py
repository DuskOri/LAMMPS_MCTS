"""重复单元状态机和快速 Green-Kubo 输入的基础检查。"""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from generator import (
    RepeatUnitLengthEstimate,
    estimate_repeat_unit_contour_length,
    resolve_degree_of_polymerization,
)
from generator.packmol_runner import _resolve_packmol
from md_engine import (
    write_direct_nemd_input,
    write_rapid_gk_input,
    write_uff_lammps_data_from_template,
)
from md_engine.lammps_runner import _should_stream_line
from mcts import configure_fragment_registry
from mcts.Poly_Build import build_poly_chain
from mcts.mcts_engine import MCTSEngine
from mcts.state import PolymerState
from post_process.thermal_analyzer import (
    analyze_direct_nemd_outputs,
    read_plateau_conductivity,
)
from rdkit import Chem
from rdkit.Chem import AllChem


class RepeatUnitStateTests(unittest.TestCase):
    def setUp(self):
        configure_fragment_registry({"active_families": "polyimide"})

    def test_end_does_not_count_as_repeat_unit_fragment(self):
        state = PolymerState(max_steps=5)
        for action in ("A", "B", "A", "B", "A"):
            state = state.take_action(action)

        self.assertFalse(state.is_terminal())
        self.assertEqual(len(state.repeat_unit_fragments), 5)
        self.assertEqual(state.legal_actions(), ["End"])

        terminal = state.take_action("End")
        self.assertTrue(terminal.is_terminal())
        self.assertEqual(len(terminal.repeat_unit_fragments), 5)
        self.assertEqual(terminal.to_mcts_output()[-1], "End")

    def test_repeat_unit_may_end_before_limit(self):
        state = PolymerState(max_steps=5).take_action("A").take_action("B")
        self.assertIn("End", state.legal_actions())
        terminal = state.take_action("End")
        self.assertEqual(terminal.to_mcts_output(), ["A", "B", "End"])

    def test_custom_fragments_only_load_in_custom_space(self):
        custom = {
            "key": "CUSTOM_TEST",
            "smiles": "[1*]CC[2*]",
            "family": "custom",
        }
        registry = configure_fragment_registry(
            {"active_families": "polyimide", "custom_fragments": [custom]}
        )
        self.assertNotIn("CUSTOM_TEST", registry.keys())
        self.assertNotIn("CUSTOM_TEST", registry.allowed_fragments("Start"))

        registry = configure_fragment_registry(
            {"active_families": "custom", "custom_fragments": [custom]}
        )
        self.assertIn("CUSTOM_TEST", registry.keys())
        self.assertIn("CUSTOM_TEST", registry.allowed_fragments("Start"))
        configure_fragment_registry({"active_families": "polyimide"})


class AutomaticChainLengthTests(unittest.TestCase):
    def setUp(self):
        configure_fragment_registry({"active_families": "polyimide"})

    def test_dp_is_selected_nearest_to_target_length(self):
        estimate = RepeatUnitLengthEstimate(
            internal_length_angstrom=11.2,
            connection_length_angstrom=1.52,
        )
        dp = resolve_degree_of_polymerization(estimate, target_length_angstrom=120.0)

        self.assertEqual(dp, 10)
        current_error = abs(estimate.chain_length(dp) - 120.0)
        self.assertLessEqual(current_error, abs(estimate.chain_length(dp - 1) - 120.0))
        self.assertLessEqual(current_error, abs(estimate.chain_length(dp + 1) - 120.0))

    def test_real_repeat_unit_has_positive_contour_length(self):
        repeat_unit = build_poly_chain(["A", "B", "End"], dp=1)
        estimate = estimate_repeat_unit_contour_length(repeat_unit)
        dp = resolve_degree_of_polymerization(estimate, target_length_angstrom=120.0)

        self.assertGreater(estimate.repeat_increment_angstrom, 0.0)
        self.assertGreaterEqual(dp, 1)
        self.assertLess(abs(estimate.chain_length(dp) - 120.0), estimate.repeat_increment_angstrom)


class SearchProgressTests(unittest.TestCase):
    def setUp(self):
        configure_fragment_registry({"active_families": "polyimide"})

    def test_ranked_search_reports_iterations_and_live_candidates(self):
        events = []

        def collect(event, message, details):
            events.append((event, message, details))

        engine = MCTSEngine(
            max_steps=2,
            reward_fn=lambda sequence: 1.0 / max(len(sequence), 1),
            random_seed=7,
            progress_callback=collect,
        )
        engine.ranked_candidates(iterations=4, top_k=2)

        iteration_events = [item for item in events if item[0] == "search"]
        candidate_events = [item for item in events if item[0] == "candidate"]
        self.assertEqual([item[2]["iteration"] for item in iteration_events], [1, 2, 3, 4])
        self.assertTrue(candidate_events)
        self.assertIn("sequence", candidate_events[-1][2]["candidate"])
        self.assertIn("score", candidate_events[-1][2]["candidate"])


class RapidGreenKuboTests(unittest.TestCase):
    def test_lammps_console_output_keeps_progress_lines(self):
        self.assertTrue(_should_stream_line("LAMMPS (22 Jul 2025)"))
        self.assertTrue(_should_stream_line("Step Temp Press v_kappa"))
        self.assertTrue(_should_stream_line("3000 301.2 -20.0 0.18"))
        self.assertTrue(_should_stream_line("WARNING: test warning"))
        self.assertTrue(_should_stream_line("ERROR: test error"))
        self.assertTrue(_should_stream_line("Density plateau check: rho=1.02"))
        self.assertFalse(_should_stream_line("Pair | 18.5 | 60.8"))

    def test_packmol_resolution_prefers_current_python_environment(self):
        with TemporaryDirectory() as directory:
            environment_dir = Path(directory)
            environment_packmol = environment_dir / "Scripts" / "packmol.exe"
            environment_packmol.parent.mkdir(parents=True)
            environment_packmol.write_text("test", encoding="utf-8")

            with patch(
                "generator.packmol_runner.sys.executable",
                str(environment_dir / "python.exe"),
            ), patch(
                "generator.packmol_runner.shutil.which",
                return_value="fallback-packmol",
            ):
                resolved = _resolve_packmol("packmol")

            self.assertEqual(resolved, str(environment_packmol))

    def test_uff_writer_generates_force_field_sections(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            molecule = Chem.AddHs(Chem.MolFromSmiles("CC"))
            AllChem.EmbedMolecule(molecule, randomSeed=7)
            AllChem.UFFOptimizeMolecule(molecule)
            mol_file = root / "ethane.mol"
            pdb_file = root / "ethane.pdb"
            data_file = root / "ethane.data"
            Chem.MolToMolFile(molecule, str(mol_file))
            Chem.MolToPDBFile(molecule, str(pdb_file))

            write_uff_lammps_data_from_template(
                mol_file,
                pdb_file,
                data_file,
                molecule_count="1",
                box_size=30.0,
            )

            data_text = data_file.read_text(encoding="utf-8")
            self.assertIn("Pair Coeffs", data_text)
            self.assertIn("Bond Coeffs", data_text)
            self.assertIn("Angle Coeffs", data_text)
            result = write_rapid_gk_input(data_file, root / "gk.in", root / "gk")
            self.assertTrue(result.success, result.message)

    def test_topology_only_data_is_rejected(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            data_file = root / "topology.data"
            data_file.write_text(
                "1 atoms\n1 atom types\n1 bonds\n1 bond types\n\nMasses\n\n1 12.011\n",
                encoding="utf-8",
            )
            result = write_rapid_gk_input(
                data_file,
                root / "gk.in",
                root / "gk",
            )

            self.assertFalse(result.success)
            self.assertIn("Pair Coeffs", result.message)
            self.assertIn("Bond Coeffs", result.message)
            self.assertFalse((root / "gk.in").exists())

    def test_parameterized_data_writes_intramolecular_hfacf(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            data_file = root / "parameterized.data"
            data_file.write_text(
                "\n".join(
                    [
                        "2 atoms",
                        "1 bonds",
                        "1 atom types",
                        "1 bond types",
                        "",
                        "Pair Coeffs",
                        "",
                        "1 0.1 3.5",
                        "",
                        "Bond Coeffs",
                        "",
                        "1 300.0 1.4",
                    ]
                ),
                encoding="utf-8",
            )
            result = write_rapid_gk_input(
                data_file,
                root / "gk.in",
                root / "gk",
                params={"molecule_count": 2, "force_field_include": None},
            )

            self.assertTrue(result.success)
            script = (root / "gk.in").read_text(encoding="utf-8")
            self.assertIn("group           mol_1 molecule 1", script)
            self.assertIn("group           mol_2 molecule 2", script)
            self.assertIn("type auto", script)
            self.assertNotIn("pair_style      zero", script)
            self.assertNotIn("exclude molecule/inter", script)
            self.assertNotIn("include         None", script)

    def test_production_thermo_matches_correlation_window(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            data_file = root / "parameterized.data"
            data_file.write_text(
                "\n".join(
                    [
                        "1 atoms",
                        "1 atom types",
                        "",
                        "Pair Coeffs",
                        "",
                        "1 0.1 3.5",
                    ]
                ),
                encoding="utf-8",
            )
            result = write_rapid_gk_input(
                data_file,
                root / "gk.in",
                root / "gk",
                params={
                    "sample_nevery": 10,
                    "correlation_samples": 300,
                    "production_steps": 3500,
                    "thermo_every": 500,
                },
            )

            self.assertTrue(result.success)
            script = (root / "gk.in").read_text(encoding="utf-8")
            self.assertIn("thermo          ${d}", script)
            self.assertIn("run             6000", script)

    def test_density_equilibration_uses_plateau_instead_of_final_density(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            data_file = root / "parameterized.data"
            data_file.write_text(
                "\n".join(
                    [
                        "1 atoms",
                        "1 atom types",
                        "",
                        "Pair Coeffs",
                        "",
                        "1 0.1 3.5",
                    ]
                ),
                encoding="utf-8",
            )
            result = write_rapid_gk_input(
                data_file,
                root / "gk.in",
                root / "gk",
                params={
                    "density_equilibration": True,
                    "precompression_density": 0.7,
                    "density_block_steps": 1000,
                    "density_sample_every": 100,
                    "density_max_blocks": 5,
                    "density_plateau_tolerance": 0.02,
                },
            )

            self.assertTrue(result.success)
            script = (root / "gk.in").read_text(encoding="utf-8")
            self.assertIn("variable        rho_pre equal 0.7", script)
            self.assertIn("variable        density_loop loop 3", script)
            self.assertIn("Density plateau check", script)
            self.assertIn("jump SELF density_plateau_done", script)
            self.assertIn("Density plateau reached within tolerance", script)
            self.assertIn("Density plateau was not reached", script)
            self.assertNotIn("target density = 1.1", script)

    def test_conductivity_uses_tail_average(self):
        with TemporaryDirectory() as directory:
            output = Path(directory) / "kappa.dat"
            output.write_text("# step kappa\n10 0.10\n20 0.20\n30 0.30\n40 0.40\n", encoding="utf-8")
            value = read_plateau_conductivity(output, tail_fraction=0.50)
            self.assertAlmostEqual(value, 0.35)


class DirectNemdTests(unittest.TestCase):
    def test_writer_keeps_reference_hot_baths_and_central_heat_flux(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            data_file = root / "parameterized.data"
            data_file.write_text(
                "1 atoms\n1 atom types\n\nPair Coeffs\n\n1 0.1 3.5\n",
                encoding="utf-8",
            )
            result = write_direct_nemd_input(
                data_file,
                root / "nemd.in",
                root / "nemd",
                params={"nemd_steady_steps": 1000, "production_steps": 2000},
            )

            self.assertTrue(result.success, result.message)
            script = (root / "nemd.in").read_text(encoding="utf-8")
            self.assertIn("fix             fhot all langevin", script)
            self.assertIn("fix_modify      fhot temp Thot", script)
            self.assertIn("compute         heat_flux flux heat/flux", script)
            self.assertIn("compute         layers all chunk/atom bin/1d x", script)
            self.assertNotIn("thermal/conductivity", script)

    def test_analyzer_uses_heat_flux_over_linear_temperature_gradient(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            profile = root / "temp.profile"
            flux = root / "flux.profile"
            rows = ["# profile", "1000 20 1000"]
            for index in range(20):
                x = 0.25 + 0.5 * index
                temperature = 450.0 - 20.0 * x
                rows.append(f"{index + 1} {x} 50 {temperature}")
            profile.write_text("\n".join(rows), encoding="utf-8")
            flux.write_text("# step Jx\n1000 2.0e-6\n2000 2.0e-6\n", encoding="utf-8")

            result = analyze_direct_nemd_outputs(
                root / "nemd.in",
                profile,
                flux,
                min_gradient_r2=0.99,
                min_temperature_span=20.0,
            )

            expected_flux = 2.0e-6 * (4184.0 / 6.02214076e23) / 1.0e-35
            expected_kappa = expected_flux / (20.0 * 1.0e10)
            self.assertTrue(result.success, result.message)
            self.assertAlmostEqual(result.gradient_k_per_a, 20.0)
            self.assertAlmostEqual(result.gradient_r2, 1.0)
            self.assertAlmostEqual(result.conductivity_w_mk, expected_kappa)

    def test_analyzer_rejects_profile_without_linear_gradient(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            profile = root / "temp.profile"
            flux = root / "flux.profile"
            rows = ["# profile", "1000 20 1000"]
            temperatures = [300, 450, 280, 430, 260] * 4
            for index, temperature in enumerate(temperatures):
                rows.append(f"{index + 1} {0.25 + 0.5 * index} 50 {temperature}")
            profile.write_text("\n".join(rows), encoding="utf-8")
            flux.write_text("# step Jx\n1000 2.0e-6\n", encoding="utf-8")

            result = analyze_direct_nemd_outputs(
                root / "nemd.in",
                profile,
                flux,
                min_gradient_r2=0.70,
            )

            self.assertFalse(result.success)
            self.assertIn("R2", result.message)


if __name__ == "__main__":
    unittest.main()
