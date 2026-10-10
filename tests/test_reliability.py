"""Regression coverage for scientific validation and run lifecycle."""

from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch
import http.client
import json
import math
import subprocess
import sys
import threading
import unittest

from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers, AllChem, rdMolTransforms
import main
import web_server
from generator.polymer_builder import build_polymer_from_sequence
from generator.packmol_runner import _run_packmol
from md_engine.process_control import PipelineCancelled, reset_cancellation, cancel_processes
from md_engine.lammps_runner import run_lammps_input
from md_engine.thermal_conductivity_writer import _normalise_params
from md_engine.uff_data_writer import _dihedral_parameters
from mcts import configure_fragment_registry
from mcts.mcts_engine import MCTSEngine
from mcts.state import PolymerState
from mcts.thermal_feedback import ThermalFeedbackEvaluator
from post_process.thermal_analyzer import analyze_direct_nemd_outputs


class ReliabilityTests(unittest.TestCase):
    def setUp(self):
        reset_cancellation()
        configure_fragment_registry({"active_families": "polyimide"})
        self.old_state = deepcopy(web_server.RUN_STATE)
        web_server.RUN_STATE["running"] = False
        web_server.RUN_CANCEL_EVENT.clear()

    def tearDown(self):
        reset_cancellation()
        web_server.RUN_STATE.clear()
        web_server.RUN_STATE.update(self.old_state)
        web_server.RUN_CANCEL_EVENT.clear()
        configure_fragment_registry({"active_families": "polyimide"})

    def config(self, root):
        config = deepcopy(main.DEFAULT_CONFIG)
        config["output"]["directory"] = str(root / "results")
        config["mcts"]["feedback_output_dir"] = str(root / "data")
        config["chemistry"]["custom_fragment_file"] = ""
        return config

    def test_nested_config_roundtrip_without_yaml(self):
        config = deepcopy(main.DEFAULT_CONFIG)
        config["chemistry"].update(custom_fragments=[{"key": "X", "smiles": "[1*]CC#CC[2*]"}], custom_transitions={"Start": ["X"], "X": []})
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.yaml"
            web_server.write_config(config, path)
            with patch.object(main, "yaml", None):
                loaded = main.load_config(path)
        self.assertEqual(loaded["chemistry"], config["chemistry"])

    def test_unchanged_profile_preserves_hidden_parameters(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = self.config(root)
            config["thermal_conductivity"].update(profile="standard", nemd_thermostat_damp=42, production_steps=1000000)
            web_server.write_config(config, root / "config.yaml")
            with patch.object(web_server, "PROJECT_ROOT", root):
                saved = web_server.update_config({"thermal_profile": "standard", "iterations": 2})
            self.assertEqual(saved["thermal_conductivity"]["nemd_thermostat_damp"], 42)
            self.assertEqual(saved["thermal_conductivity"]["production_steps"], 1000000)

    def test_insufficient_production_limit_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "at least"):
            _normalise_params({"production_steps": 1000})
        params = _normalise_params({"production_steps": 101001})
        self.assertLessEqual(params["production_steps"], 101001)

    def test_custom_start_empty_edges_and_step_limit_are_respected(self):
        registry = configure_fragment_registry({"active_families": "custom", "custom_fragments": [{"key": "X", "smiles": "[1*]CC[2*]"}, {"key": "Y", "smiles": "[1*]COC[2*]"}], "custom_transitions": {"Start": ["X"], "X": ["Y"], "Y": []}})
        self.assertEqual(registry.allowed_fragments("Start"), ["X"])
        self.assertEqual(registry.allowed_fragments("Y"), [])
        dead = PolymerState(max_steps=1).take_action("X")
        self.assertEqual(dead.legal_actions(), [])
        with self.assertRaises(ValueError):
            dead.to_mcts_output()
        engine = MCTSEngine(max_steps=1, reward_fn=Mock(side_effect=AssertionError("must not evaluate dead ends")))
        self.assertEqual(engine.ranked_candidates(4), [])

    def test_invalid_custom_fragments_rejected_before_saving(self):
        with TemporaryDirectory() as tmp, patch.object(web_server, "PROJECT_ROOT", Path(tmp)):
            for key, smiles in [("Start", "[1*]CC[2*]"), ("X", "C"), ("X", "not smiles")]:
                with self.assertRaises(ValueError):
                    web_server.save_custom_graph({"fragments": [{"key": key, "smiles": smiles}], "transitions": {}})
            self.assertFalse((Path(tmp) / "custom_fragments.json").exists())

    def test_bad_build_returns_candidate_failure(self):
        with patch("generator.polymer_builder._load_poly_build_tools", return_value=(Mock(side_effect=ValueError("bad valence")), Mock(), Mock())):
            result = build_polymer_from_sequence(["A", "B", "End"])
        self.assertFalse(result.success)
        self.assertIn("bad valence", result.message)

    def test_evaluations_are_isolated_and_persisted_before_callback(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            original = root / "data" / "eval_0001" / "nemd_heat_flux.dat"
            original.parent.mkdir(parents=True)
            original.write_text("1000 999\n", encoding="utf-8")
            first = main.prepare_run_config(self.config(root))
            second = main.prepare_run_config(self.config(root))
            self.assertNotEqual(first["runtime"]["run_id"], second["runtime"]["run_id"])
            evaluator = ThermalFeedbackEvaluator(first, first["mcts"]["feedback_output_dir"], run_lammps=False)
            with patch.object(evaluator, "_evaluate", side_effect=ValueError("candidate invalid")):
                evaluator(["A", "B", "End"])
            record = evaluator.records[0]
            self.assertFalse(record.success)
            self.assertTrue(Path(record.metadata_path).exists())
            rows = web_server.read_csv(Path(first["output"]["directory"]) / "feedback_records.csv")
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["run_id"], first["runtime"]["run_id"])
            self.assertEqual(original.read_text(), "1000 999\n")

    def test_search_checkpoint_survives_callback_cancel(self):
        with TemporaryDirectory() as tmp:
            config = main.prepare_run_config(self.config(Path(tmp)))
            config["mcts"].update(feedback_mode="heuristic", max_steps=2, iterations=3)
            def stop(event, message, details):
                if event == "candidate":
                    raise PipelineCancelled("stop after first candidate")
            with self.assertRaises(PipelineCancelled):
                main.run_search(config, progress_callback=stop)
            directory = Path(config["output"]["directory"])
            self.assertEqual(len(web_server.read_csv(directory / "candidates.csv")), 1)
            self.assertEqual(json.loads((directory / "search_checkpoint.json").read_text())["iteration"], 1)

    def test_custom_output_names_and_thermal_api_records(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = self.config(root)
            config["output"]["feedback_records_csv"] = "measurements.csv"
            web_server.write_config(config, root / "config.yaml")
            prepared = main.prepare_run_config(config)
            directory = Path(prepared["output"]["directory"])
            web_server.write_csv_rows(directory / "measurements.csv", [{"success": "True", "conductivity_w_mk": "0.2"}], ["success", "conductivity_w_mk"])
            with patch.object(web_server, "PROJECT_ROOT", root):
                result = web_server.load_results()
            self.assertEqual(len(result["thermal_results"]), 1)
            self.assertEqual(result["thermal_results"][0]["measurement_status"], "未通过新校验")

    def test_exports_failing_do_not_erase_search_records(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = self.config(root)
            config["mcts"].update(feedback_mode="heuristic", iterations=1, max_steps=2)
            web_server.write_config(config, root / "config.yaml")
            with patch.object(web_server, "PROJECT_ROOT", root), patch("main.build_candidates", side_effect=ValueError("export failure")):
                web_server.run_pipeline()
            pointer = json.loads((root / "results" / "latest_run.json").read_text())
            self.assertEqual(pointer["status"], "failed")
            self.assertTrue(web_server.read_csv(Path(pointer["directory"]) / "candidates.csv"))

    def test_packmol_timeout_closes_real_child(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / "input.py"
            source.write_text("import time\nprint('started', flush=True)\ntime.sleep(30)\n", encoding="utf-8")
            result = _run_packmol(sys.executable, source, timeout=0.3)
            self.assertEqual(result.returncode, 124)
            self.assertIn("timeout", result.stdout)

    def test_callback_exception_closes_lammps_child(self):
        process = Mock()
        process.poll.return_value = None
        process.stdout.readline.side_effect = ["LAMMPS (test)\n", ""]
        with TemporaryDirectory() as tmp, patch("md_engine.lammps_runner.subprocess.Popen", return_value=process):
            with self.assertRaisesRegex(ValueError, "callback"):
                run_lammps_input(Path(tmp) / "job.in", line_callback=Mock(side_effect=ValueError("callback failed")), stream_output=False)
        process.terminate.assert_called_once()
        process.wait.assert_called()
        process.stdout.close.assert_called()

    def test_ethane_torsion_barrier_is_distributed(self):
        molecule = Chem.AddHs(Chem.MolFromSmiles("CC"))
        rows, types = _dihedral_parameters(molecule)
        self.assertEqual(len(rows), 9)
        expected = rdForceFieldHelpers.GetUFFTorsionParams(molecule, 2, 0, 1, 5)
        self.assertAlmostEqual(sum(types[row[0]-1][0] * 2 for row in rows), expected)

    def test_ethane_torsional_energy_scan_matches_rdkit(self):
        molecule = Chem.AddHs(Chem.MolFromSmiles("CC"))
        self.assertEqual(AllChem.EmbedMolecule(molecule, randomSeed=7), 0)
        rows, types = _dihedral_parameters(molecule)
        differences = []
        for angle in (0, 20, 40, 60, 80, 100, 120):
            rdMolTransforms.SetDihedralDeg(molecule.GetConformer(), 2, 0, 1, 5, angle)
            # Exclude nonbonded terms to compare the torsional energy change alone.
            force_field = AllChem.UFFGetMoleculeForceField(molecule, vdwThresh=0)
            reference = force_field.CalcEnergy()
            exported = 0.0
            for type_id, *atoms in rows:
                coefficient, sign, period = types[type_id - 1]
                phi = rdMolTransforms.GetDihedralDeg(molecule.GetConformer(), *(atom - 1 for atom in atoms))
                exported += coefficient * (1 + sign * math.cos(period * math.radians(phi)))
            differences.append(reference - exported)
        self.assertLess(max(differences) - min(differences), 1e-6)

    def test_offline_input_preparation_never_analyzes_old_outputs(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = main.prepare_run_config(self.config(root))
            data = root / "prepared.data"
            data.write_text("1 atoms\n1 atom types\n\nPair Coeffs\n\n1 0.1 3.5\n", encoding="utf-8")
            evaluator = ThermalFeedbackEvaluator(config, config["mcts"]["feedback_output_dir"], run_lammps=False)
            build = Mock(success=True, pdb_path="candidate.pdb", mol_path="candidate.mol")
            system = Mock(success=True, system_data_path=str(data))
            with patch("mcts.thermal_feedback.build_polymer_from_sequence", return_value=build), patch("mcts.thermal_feedback.prepare_initial_system", return_value=system), patch("mcts.thermal_feedback.analyze_direct_nemd_outputs") as analyze:
                evaluator(["A", "B", "End"])
            analyze.assert_not_called()
            self.assertFalse(evaluator.records[0].success)
            self.assertIn("not executed", evaluator.records[0].message)

    def test_torsion_periods_include_mixed_and_group6_exceptions(self):
        propene = Chem.AddHs(Chem.MolFromSmiles("CC=C"))
        _, types = _dihedral_parameters(propene)
        self.assertIn((1, 3), [row[1:] for row in types])
        self.assertIn((-1, 6), [row[1:] for row in types])
        peroxide = Chem.AddHs(Chem.MolFromSmiles("COOC"))
        rows, types = _dihedral_parameters(peroxide)
        oxygens = [row for row in rows if peroxide.GetAtomWithIdx(row[2]-1).GetAtomicNum() == 8 and peroxide.GetAtomWithIdx(row[3]-1).GetAtomicNum() == 8]
        self.assertTrue(oxygens)
        self.assertTrue(all(types[row[0]-1][1:] == (1, 2) for row in oxygens))

    def test_cancel_stops_packmol_child(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / "input.py"
            source.write_text("import time\nprint('started', flush=True)\ntime.sleep(30)\n", encoding="utf-8")
            launched = threading.Event()
            real_popen = subprocess.Popen
            def launch(*args, **kwargs):
                process = real_popen(*args, **kwargs)
                launched.set()
                return process
            timer = threading.Timer(0.4, cancel_processes)
            with patch("generator.packmol_runner.subprocess.Popen", side_effect=launch):
                timer.start()
                try:
                    with self.assertRaises(PipelineCancelled):
                        _run_packmol(sys.executable, source, timeout=10)
                finally:
                    timer.cancel()
            self.assertTrue(launched.is_set())


class ThermalQualityTests(unittest.TestCase):
    def analyze(self, root, values=None, balance=1.0, agreement=1.0, reverse=False):
        values = values or [2e-6] * 6
        profile = []
        for step in range(100, 601, 100):
            profile.append(f"{step} 20 1000")
            for index in range(20):
                x = index * 0.5 + 0.25
                profile.append(f"{index+1} {x} 50 {450 - 20*x}")
        (root / "temperature.dat").write_text("\n".join(profile), encoding="utf-8")
        (root / "flux.dat").write_text("\n".join(f"{(i+1)*100} {-v if reverse else v}" for i, v in enumerate(values)), encoding="utf-8")
        (root / "baths.dat").write_text("\n".join(f"{step} {-step*0.5*2e-4*agreement} {step*0.5*2e-4*agreement*balance} 100" for step in range(100, 601, 100)), encoding="utf-8")
        return analyze_direct_nemd_outputs("job.in", root / "temperature.dat", root / "flux.dat", bath_energy_path=root / "baths.dat", timestep_fs=0.5, block_steps=100)

    def test_balanced_stable_measurement_passes(self):
        with TemporaryDirectory() as tmp:
            result = self.analyze(Path(tmp))
        self.assertTrue(result.success, result.message)
        self.assertAlmostEqual(result.conductivity_w_mk, 0.006947695457055375)
        self.assertLess(result.energy_balance_error, 1e-12)

    def test_energy_imbalance_rejected(self):
        with TemporaryDirectory() as tmp:
            result = self.analyze(Path(tmp), balance=0.2)
        self.assertFalse(result.success)
        self.assertIn("balance", result.message)

    def test_microscopic_bath_disagreement_rejected(self):
        with TemporaryDirectory() as tmp:
            result = self.analyze(Path(tmp), agreement=3)
        self.assertFalse(result.success)
        self.assertIn("disagree", result.message)

    def test_noisy_flux_and_reversed_flux_rejected(self):
        with TemporaryDirectory() as tmp:
            noisy = self.analyze(Path(tmp), values=[2e-6, 2e-6, 2e-6, 1e-7, 8e-6, 1e-7])
            reversed_result = self.analyze(Path(tmp), reverse=True)
        self.assertFalse(noisy.success)
        self.assertIn("error", noisy.message)
        self.assertFalse(reversed_result.success)


class LocalHttpTests(unittest.TestCase):
    def test_same_origin_token_and_static_file_boundary(self):
        server = web_server.ThreadingHTTPServer(("127.0.0.1", 0), web_server.LammpsMctsHandler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        connection = http.client.HTTPConnection("127.0.0.1", server.server_port)
        try:
            connection.request("GET", "/api/session")
            response = connection.getresponse()
            token = json.loads(response.read())["token"]
            connection.request("POST", "/api/config", body="{}", headers={"Content-Type": "application/json"})
            response = connection.getresponse()
            self.assertEqual(response.status, 403)
            response.read()
            with patch("web_server.update_config", return_value=deepcopy(main.DEFAULT_CONFIG)) as update:
                connection.request("POST", "/api/config", body="{}", headers={"Content-Type": "application/json", "X-CSRF-Token": token})
                response = connection.getresponse()
                self.assertEqual(response.status, 200)
                self.assertTrue(json.loads(response.read())["ok"])
                update.assert_called_once()
            for path in ("/config.yaml", "/main.py", "/frontend/%2e%2e/config.yaml"):
                connection.request("GET", path)
                response = connection.getresponse()
                self.assertEqual(response.status, 404, path)
                response.read()
            connection.request("GET", "/api/config", headers={"Origin": "https://evil.example"})
            response = connection.getresponse()
            self.assertEqual(response.status, 403)
            self.assertIsNone(response.getheader("Access-Control-Allow-Origin"))
            response.read()
            connection.request("GET", "/frontend/")
            response = connection.getresponse()
            self.assertEqual(response.status, 200)
            response.read()
        finally:
            connection.close()
            server.shutdown()
            server.server_close()
            thread.join(timeout=3)


if __name__ == "__main__":
    unittest.main()
