import os
import subprocess
import sys
import unittest
from pathlib import Path

from rdkit import Chem

from dimorphite_dl import dimorphite_dl


class MolConversionTests(unittest.TestCase):
    @unittest.skipIf(os.name == "nt", "Requires POSIX file descriptors")
    def test_conversion_keeps_stderr_usable_with_concurrent_file_open(self):
        script = """
import os
import sys
import threading
from unittest.mock import patch

from dimorphite_dl import dimorphite_dl

original_close = os.close
saved_stderr = os.dup(2)
unrelated_files = []

def close_and_open_in_another_thread(fd):
    original_close(fd)
    if fd == 2:
        # Force another thread to take the descriptor before stderr is restored.
        thread = threading.Thread(
            target=lambda: unrelated_files.append(open(os.devnull, "w"))
        )
        thread.start()
        thread.join(timeout=5)
        assert not thread.is_alive()
        assert unrelated_files[-1].fileno() == 2

try:
    with patch.object(os, "close", close_and_open_in_another_thread):
        assert dimorphite_dl.UtilFuncs.convert_smiles_str_to_mol("CCN") is not None
    for stream in unrelated_files:
        stream.close()
    sys.stderr.write("worker startup flush probe")
    sys.stderr.flush()
finally:
    os.dup2(saved_stderr, 2)
    original_close(saved_stderr)
"""
        result = subprocess.run(
            [sys.executable, "-c", script],
            cwd=Path(__file__).resolve().parents[1],
            capture_output=True,
            text=True,
            check=False,
            timeout=15,
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_physiological_charges(self):
        for smiles, expected in (
            ("CCN", 1),
            ("CC(=O)O", -1),
            ("CCO", 0),
            ("C[N+](C)(C)C", 1),
        ):
            with self.subTest(smiles=smiles):
                mol = dimorphite_dl.run_single_mol(smiles, ph=7.35)
                self.assertEqual(Chem.GetFormalCharge(mol), expected)
