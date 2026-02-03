import json
import pathlib
import shutil
import subprocess
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]


class TestCpSmoothing2412(unittest.TestCase):
    def test_cp_smoothing_le_spike(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS smoothing checks")

        script = pathlib.Path(__file__).with_name("compare_cp_smoothing_2412.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            text=True,
            capture_output=True,
            check=True,
        )
        data = json.loads(proc.stdout)["results"]

        self.assertTrue(data.get("converged", False), msg="Viscous solve did not converge")

        for side in ("upper", "lower"):
            before = data["before"][side]["spike"]
            after = data["after"][side]["spike"]
            if before is None or after is None:
                self.fail(f"Not enough points to compute LE spike for {side} side")
            # Smoothing should not amplify the LE spike.
            self.assertLessEqual(
                after,
                before * 1.05 + 1.0e-6,
                msg=f"LE spike grew on {side} side: before={before}, after={after}",
            )


if __name__ == "__main__":
    unittest.main()
