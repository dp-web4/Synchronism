"""Exercise the actual CI PDF step with disposable files and a fake builder."""

import re
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


class PdfCiContract(unittest.TestCase):
    @staticmethod
    def step_script():
        root = Path(__file__).resolve().parents[2]
        workflow = (root / ".github/workflows/build_whitepaper.yml").read_text()
        match = re.search(
            r"(?m)^      - name: Build PDF Version\n"
            r"        working-directory: \./whitepaper\n"
            r"        run: \|\n((?: {10}[^\n]*\n)+)",
            workflow,
        )
        if not match:
            raise AssertionError("PDF build step is missing or has an unexpected shape")
        return textwrap.dedent(match.group(1))

    def run_step(self, builder):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "build").mkdir()
            artifact = root / "build/Synchronism_Whitepaper.pdf"
            artifact.write_bytes(b"stale artifact")
            (root / "make-pdf.sh").write_text("#!/bin/bash\n" + builder + "\n")
            result = subprocess.run(
                ["bash", "-e", "-o", "pipefail", "-c", self.step_script()],
                cwd=root, capture_output=True, text=True, check=False,
            )
            return result.returncode, artifact.read_bytes() if artifact.exists() else None

    def test_failed_builder_is_not_reported_as_success(self):
        code, artifact = self.run_step("exit 7")
        self.assertNotEqual(code, 0)
        self.assertIsNone(artifact)

    def test_success_without_new_pdf_does_not_reuse_stale_output(self):
        code, artifact = self.run_step("exit 0")
        self.assertNotEqual(code, 0)
        self.assertIsNone(artifact)

    def test_empty_new_pdf_is_not_success(self):
        code, _ = self.run_step(": > build/Synchronism_Whitepaper.pdf")
        self.assertNotEqual(code, 0)

    def test_nonempty_new_pdf_succeeds(self):
        code, artifact = self.run_step("printf fresh > build/Synchronism_Whitepaper.pdf")
        self.assertEqual(code, 0)
        self.assertEqual(artifact, b"fresh")


if __name__ == "__main__":
    unittest.main()
