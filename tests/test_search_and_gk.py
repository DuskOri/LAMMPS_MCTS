"""重复单元状态机和快速 Green-Kubo 输入的基础检查。"""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from generator import (
    RepeatUnitLengthEstimate,
    estimate_repeat_unit_contour_length,
    resolve_degree_of_polymerization,
)
from md_engine import write_rapid_gk_input, write_uff_lammps_data_from_template
from mcts import configure_fragment_registry
from mcts.Poly_Build import build_poly_chain
from mcts.mcts_engine import MCTSEngine
from mcts.state import PolymerState
from post_process.thermal_analyzer import read_plateau_conductivity
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

    def test_conductivity_uses_tail_average(self):
        with TemporaryDirectory() as directory:
            output = Path(directory) / "kappa.dat"
            output.write_text("# step kappa\n10 0.10\n20 0.20\n30 0.30\n40 0.40\n", encoding="utf-8")
            value = read_plateau_conductivity(output, tail_fraction=0.50)
            self.assertAlmostEqual(value, 0.35)


if __name__ == "__main__":
    unittest.main()
