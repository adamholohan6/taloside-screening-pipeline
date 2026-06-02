import subprocess
import sys
from pathlib import Path


def run_phase(module_name: str) -> None:
    print("--- Starting {} ---".format(module_name))

    # Ensure the `src/` layout is discoverable for subprocesses.
    repo_root = Path(__file__).resolve().parent
    src_path = str(repo_root / "src")

    env = {}
    env["PYTHONPATH"] = src_path

    subprocess.run(
        [sys.executable, "-m", module_name],
        check=True,
        env=env,
    )


if __name__ == "__main__":
    try:
        run_phase("taloside_pipeline.glycolibrary_generator")
        run_phase("taloside_pipeline.phase2_integration")
        run_phase("taloside_pipeline.phase3_docking")
        print("PIPELINE SUCCESSFUL")
    except Exception as e:
        print("PIPELINE FAILED: {}".format(e))
