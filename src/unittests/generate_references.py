#!/usr/bin/env python3
"""Freeze minimal full-system trajectories; never run automatically by tests."""

import argparse
import copy
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time

from reference_cases import reference_cases, stage_options


ROOT = Path(__file__).resolve().parent
SECTIONS = {"ir": "INFRARED", "raman": "RAMAN", "sfg": "SFG", "power": "POWER"}


def first_frames(path, count=2):
    """Preserve atom records and metadata; select frames only, never molecules."""
    lines = path.read_text().splitlines(keepends=True)
    position = 0
    xyz = path.suffix == ".xyz"
    for _ in range(count):
        header_size = 2 if xyz else 9
        atoms = int(lines[position if xyz else position + 3])
        end = position + header_size + atoms
        if end > len(lines):
            raise ValueError(f"Incomplete input: {path}")
        position = end
    return "".join(lines[:position])


def add_unused_dipole_columns(text):
    """Satisfy the SFG reader when PS replaces all trajectory dipoles."""
    lines = text.splitlines()
    names = ("v_dipx_perm", "v_dipy_perm", "v_dipz_perm",
             "v_dipx_ind", "v_dipy_ind", "v_dipz_ind")
    position = 0
    while position < len(lines):
        count = int(lines[position + 3])
        columns = lines[position + 8].split()[2:]
        missing = [name for name in names if name not in columns]
        if missing:
            lines[position + 8] += " " + " ".join(missing)
            for atom in range(position + 9, position + 9 + count):
                lines[atom] += " 0" * len(missing)
        position += 9 + count
    return "\n".join(lines) + "\n"


def numeric_table(text):
    lines = text.splitlines()
    columns = lines[0].lstrip("# ").split()
    rows = [[float(value) for value in line.split()] for line in lines[1:] if line.strip()]
    return {"columns": columns, "rows": rows}


def xyz_properties(binary, trajectory, config, directory, environment, section, columns):
    """Evaluate variable-cell XYZ frames using each frame's JSON cell."""
    lines = trajectory.read_text().splitlines(keepends=True)
    position = 0
    rows = []
    while position < len(lines):
        count = int(lines[position])
        comment = lines[position + 1]
        lattice = [float(value) for value in comment.split('Lattice="', 1)[1].split('"', 1)[0].split()]
        frame_config = copy.deepcopy(config)
        frame_config[section]["box"] = lattice
        frame_input = directory / "frame.xyz"
        frame_input.write_text("".join(lines[position:position + 2 + count]))
        config_path = directory / "frame.json"
        config_path.write_text(json.dumps(frame_config, indent=2) + "\n")
        with (directory / (binary.name + "-frame-" + str(len(rows)) + ".log")).open("w") as stderr:
            result = subprocess.run([str(binary), str(frame_input), str(config_path)],
                                    cwd=directory, env=environment, stdout=subprocess.PIPE,
                                    stderr=stderr, check=True, timeout=180, text=True)
        values = [[float(value) for value in line.split()]
                  for line in result.stdout.splitlines() if line.strip()]
        if len(values) != 1 or len(values[0]) != columns:
            raise ValueError("Expected one XYZ property sample per frame")
        values[0][0] = len(rows) * config[section]["dump_frequency"]
        rows.append(values[0])
        position += 2 + count
    return "".join(" ".join(format(value, ".17g") for value in row) + "\n" for row in rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replace", action="store_true", help="intentionally replace an existing baseline")
    parser.add_argument("--suite", choices=SECTIONS, help="generate references for one suite only")
    parser.add_argument("--bindir", type=Path, default=ROOT.parent,
                        help="compiled src directory; defaults to the source tree")
    args = parser.parse_args()
    binary_root = args.bindir.resolve()
    selected_suites = [args.suite] if args.suite else list(SECTIONS)
    target = ROOT / "references" / "ir.json"
    if any((target.parent / (suite + ".json")).exists() for suite in selected_suites) and not args.replace:
        parser.error("Reference exists; review changes before using --replace")

    specifications = [case for case in reference_cases(ROOT) if case["suite"] in selected_suites]
    environment = dict(os.environ, OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE", MKL_NUM_THREADS="1")
    run_root = ROOT / ".build" / "reference-generation"
    run_root.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="baseline-", dir=run_root))
    cases = []
    started = time.monotonic()
    stage_cache = {}
    cache_root = ROOT / ".build" / "reference-cache"
    cache_root.mkdir(exist_ok=True)
    for index, spec in enumerate(specifications):
        suite, label, filename = spec["suite"], spec["name"], spec["input"]
        case_start = time.monotonic()
        directory = work / str(index)
        directory.mkdir()
        source = ROOT / "inputs" / filename
        frame_count = spec["frames"]
        config = spec["config"]
        section = SECTIONS[suite]
        options = config[section]
        trajectory = directory / filename
        selected = first_frames(source, frame_count)
        if suite == "sfg" and options["use_ps_dipole"] and options["dipole_mb"] == "1b":
            selected = add_unused_dipole_columns(selected)
        if spec["cell"] == "matrix":
            box = options["box"]
            options["box"] = [box[1]-box[0],0,0,0,box[3]-box[2],0,0,0,box[5]-box[4]]
            if source.suffix == ".lammpstrj":
                lines = selected.splitlines(keepends=True)
                position = 0
                while position < len(lines):
                    count = int(lines[position+3])
                    lines[position+4] = "ITEM: BOX BOUNDS xy xz yz pp pp pp\n"
                    for axis in range(3):
                        lines[position+5+axis] = lines[position+5+axis].rstrip() + " 0\n"
                    position += 9+count
                selected = "".join(lines)
        trajectory.write_text(selected)
        config_path = directory / "config.json"
        config_path.write_text(json.dumps(config, indent=2) + "\n")
        programs = {
            "ir": [("infrared/dipole",4), ("infrared/response",2), ("infrared/spectrum",2)],
            "raman": [("raman/polarizability",7), ("raman/response",2), ("raman/spectrum",2)],
            "sfg": [("sfg/response",4), ("sfg/spectrum",4)],
            "power": [("power/response",2), ("power/spectrum",2)],
        }[suite]
        stages, reuse = [], {}
        current = trajectory
        for stage_index, (executable, columns) in enumerate(programs):
            stage_name = executable.split("/")[-1]
            output_name = stage_name + ".dat"
            output = directory / output_name
            binary = binary_root / executable
            signature = dict(program=executable, binary=hashlib.sha256(binary.read_bytes()).hexdigest(),
                             input=hashlib.sha256(current.read_bytes()).hexdigest(),
                             options=stage_options(suite,stage_name,options))
            per_frame_xyz = stage_name in ("dipole", "polarizability") and spec["cell"] == "per-frame-xyz"
            if per_frame_xyz:
                signature["per_frame_xyz"] = True
            key = hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()
            # Separate tiers keep ordinary cases independent of skipped extended cases.
            tier_key = (key, spec["extended"])
            cached = cache_root / key
            if tier_key in stage_cache:
                previous, data = stage_cache[tier_key]
                reuse[output_name] = previous
                output.write_text(data)
            elif (cached / "output.dat").exists():
                output.write_bytes((cached / "output.dat").read_bytes())
                sidecar = "molecular_dipoles.dat" if suite == "ir" else "molecular_polarizabilities.dat"
                if (cached / sidecar).exists():
                    (directory / sidecar).write_bytes((cached / sidecar).read_bytes())
            else:
                if per_frame_xyz:
                    output.write_text(xyz_properties(binary, current, config, directory, environment,
                                                     section, columns))
                else:
                    with output.open("w") as stdout, (directory / (stage_name+".log")).open("w") as stderr:
                        subprocess.run([str(binary),str(current),str(config_path)],cwd=directory,
                                       env=environment,stdout=stdout,stderr=stderr,check=True,timeout=180)
            rows = [[float(value) for value in line.split()] for line in output.read_text().splitlines()]
            if not rows or any(len(row) != columns or not all(map(math.isfinite,row)) for row in rows):
                raise ValueError(f"Invalid reference output: {output}")
            cached.mkdir(exist_ok=True)
            (cached / "output.dat").write_bytes(output.read_bytes())
            (cached / "settings.json").write_text(json.dumps(signature,indent=2)+"\n")
            stage_cache.setdefault(tier_key,(label+"/"+output_name,output.read_text()))
            stages.append(dict(executable=executable,expected=rows))
            current = output
        molecular = None
        if options.get("write_molecular_properties"):
            sidecar = "molecular_dipoles.dat" if suite == "ir" else "molecular_polarizabilities.dat"
            path = directory / sidecar
            # A same-run reused property stage keeps its sidecar in the stage cache.
            if not path.exists():
                for stage in stages:
                    if stage["executable"].endswith(("dipole","polarizability")):
                        break
                raise RuntimeError(f"Missing molecular output: {path}")
            molecular = (sidecar,path.read_text())
            # The property stage's persistent entry is distinct when output is enabled.
            property_binary = binary_root / programs[0][0]
            signature = dict(program=programs[0][0],binary=hashlib.sha256(property_binary.read_bytes()).hexdigest(),
                             input=hashlib.sha256(trajectory.read_bytes()).hexdigest(),
                             options=stage_options(suite,programs[0][0].split("/")[-1],options))
            property_key = hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()
            (cache_root / property_key / sidecar).write_text(molecular[1])
        cases.append(dict(suite=suite,label=label,input=filename,extended=spec["extended"],
                          frame_count=frame_count,config=config,
                          stages=stages,reuse=reuse,cell=spec["cell"],covers=spec["covers"],molecular=molecular))
        print(f"{label}: {time.monotonic()-case_start:.2f}s",flush=True)

    headers = {
        "infrared/dipole": "time_ps mu_x_D mu_y_D mu_z_D",
        "infrared/response": "time_ps C_mu_mu_C2_m2",
        "infrared/spectrum": "wavenumber_cm-1 n_alpha_cm-1",
        "raman/polarizability": "time_ps alpha_xx_A3 alpha_xy_A3 alpha_xz_A3 alpha_yy_A3 alpha_yz_A3 alpha_zz_A3",
        "raman/response": "time_ps C_alpha_alpha_A6",
        "raman/spectrum": "wavenumber_cm-1 intensity",
        "sfg/response": "time_ps response raw_sum time_origins",
        "sfg/spectrum": "wavenumber_cm-1 imaginary real magnitude_squared",
        "power/response": "time_ps velocity_correlation",
        "power/spectrum": "wavenumber_cm-1 normalized_intensity",
    }
    catalog = {suite: [] for suite in SECTIONS}
    generated = datetime.datetime.now(datetime.timezone.utc).isoformat()
    bundles = {suite: {"generated": generated, "cases": []} for suite in selected_suites}
    for suite in SECTIONS:
        path = target.parent / (suite + ".json")
        if suite not in selected_suites and path.exists():
            catalog[suite] = json.loads(path.read_text())["cases"]
    target.parent.mkdir(exist_ok=True)
    for case in cases:
        identifier = case["label"]
        outputs = []
        for stage in case["stages"]:
            filename = stage["executable"].split("/")[-1] + ".dat"
            outputs.append({"file": filename,
                            "columns": headers[stage["executable"]].split(),
                            "rows": stage["expected"]})
        entry = dict(name=identifier, input=case["input"], frames=case["frame_count"],
                     extended=case["extended"], outputs=outputs,
                     config=case["config"], reuse=case["reuse"], cell=case["cell"],
                     covers=case["covers"])
        if case["molecular"]:
            filename, data = case["molecular"]
            entry["molecular"] = {"filename": filename, "table": numeric_table(data)}
        catalog[case["suite"]].append(entry)
        bundles[case["suite"]]["cases"].append(entry)
    if "sfg" in selected_suites:
        directory = work / "fresnel"
        directory.mkdir()
        config = json.loads((ROOT / "inputs" / "fresnel.json").read_text())
        components = {"xxz": "sfg-air-1b/spectrum.dat",
                      "zzz": "sfg-zzz-both/spectrum.dat",
                      "xzx": "sfg-xzx-both/spectrum.dat"}
        settings = copy.deepcopy(config)
        for component, source in components.items():
            label = source.split("/")[0]
            case = next(entry for entry in catalog["sfg"] if entry["name"] == label)
            rows = next(output["rows"] for output in case["outputs"]
                        if output["file"] == "spectrum.dat")
            # Fresnel rejects zero infrared wavenumber; retain every positive sample.
            path = directory / (component + ".dat")
            path.write_text("".join(" ".join(format(value, ".17g") for value in row) + "\n"
                                    for row in rows if row[0] > 0))
            settings[component + "_file"] = str(path)
        for key in ("water_refractive_index_n_file", "water_refractive_index_k_file"):
            path = ROOT.parent / "sfg" / "fresnel_factors" / config[key]
            settings[key] = str(path)
        config_path = directory / "fresnel.json"
        config_path.write_text(json.dumps(settings, indent=2) + "\n")
        with (directory / "fresnel.log").open("w") as stderr:
            subprocess.run([str(binary_root / "sfg" / "fresnel"), str(config_path)],
                           cwd=directory, env=environment, stderr=stderr, check=True, timeout=180)
        outputs = []
        for polarization in ("ssp", "pss", "sps", "ppp"):
            filename = settings[polarization + "_output_file"]
            rows = [[float(value) for value in line.split()]
                    for line in (directory / filename).read_text().splitlines()]
            if not rows or any(len(row) != 2 or not all(map(math.isfinite, row)) or row[1] < 0
                               for row in rows) or not any(row[1] > 0 for row in rows):
                raise ValueError("Invalid or entirely zero Fresnel reference: " + filename)
            outputs.append(dict(file=filename, columns=["wavenumber_cm-1", "intensity"], rows=rows))
        bundles["sfg"]["fresnel"] = dict(name="sfg-fresnel-air-1b", components=components,
                                          positive_wavenumbers_only=True, config=config,
                                          outputs=outputs)
        print("sfg-fresnel-air-1b: saved SSP, PSS, SPS, and PPP intensities", flush=True)
    for suite, bundle in bundles.items():
        path = target.parent / (suite + ".json")
        with path.open("w" if args.replace else "x") as stream:
            json.dump(bundle, stream, indent=2)
            stream.write("\n")
    print(f"Saved {len(cases)} cases to {target.parent}; total {time.monotonic()-started:.2f}s")
    print(f"Generation logs: {work}")


if __name__ == "__main__":
    main()
