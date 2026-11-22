"""The settings used for saved trajectory results."""

import copy
import json

SECTIONS = {"ir": "INFRARED", "raman": "RAMAN", "sfg": "SFG", "power": "POWER"}


def reference_cases(root):
    cases = []

    def add(suite, name, filename, changes=None, *, frames=4, extended=False,
            cell=None, covers=()):
        config = json.loads((root / "inputs" / (suite + ".json")).read_text())
        options = config[SECTIONS[suite]]
        options.update(dump_frequency=.004, int_max=.016, spectrum_transform_time_ps=.016,
                       trajectory_duration_ps=0., skip_initial_atoms=0, tau=0.,
                       smoothing_fwhm_cm1=0.)
        if suite == "ir":
            options.update(DMS=True, dipole_mb="1b")
        elif suite == "raman":
            options.update(polarizability_mb="1b", raman_response="depolarized",
                           subtract_mean_polarizability=True)
        elif suite == "sfg":
            options.update(dipole_mb="1b", polarizability_mb="1b", ijk="xxz",
                           include_autocorrelation=True, include_cross_correlation=False,
                           spectrum_num_wavenumbers=33, cross_correlation_sampling_stride=1)
        else:
            options["skip_initial_atoms"] = 1
        options.update(changes or {})
        cases.append(dict(suite=suite, name=name, input=filename, config=config, frames=frames,
                          extended=extended, cell=cell, covers=list(covers)))

    # Raman NB modes are intentionally outside the test scope.
    add("ir", "ir-bulk-1b", "bulk-water.xyz", covers=["DMS=true", "dipole_mb=1b"])
    add("ir", "ir-internal-2b", "bulk-water.xyz", {"dipole_mb": "1b+2b"},
        extended=True, covers=["DMS=true", "dipole_mb=1b+2b"])
    modes = ["1b", "nb", "1b+nb", "1b+2b", "1b+2b+nb"]
    for mode, tag in zip(modes, ["1b", "nb", "1b-nb", "1b-2b", "mb"]):
        add("ir", "ir-bulk-mb" if tag == "mb" else "ir-external-" + tag, "bulk-water.lammpstrj",
            {"DMS": False, "dipole_mb": mode}, extended="2b" in mode,
            covers=["DMS=false", "dipole_mb=" + mode])
    add("raman", "raman-bulk-1b", "bulk-water.xyz", covers=["polarizability_mb=1b"])
    for mode, tag in [("1b","1b"),("1b+2b","1b-2b")]:
        add("raman", "raman-lammps-" + tag, "bulk-water.lammpstrj",
            {"polarizability_mb": mode}, extended=mode != "1b",
            covers=["polarizability_mb=" + mode])
    for dipole, tag in zip(modes, ["1b", "nb", "1b-nb", "1b-2b", "mb"]):
        for polar, suffix in [("1b", "1b"), ("1b+2b", "2b")]:
            name = "sfg-air-1b" if (tag,suffix) == ("1b","1b") else "sfg-air-mb" if (tag,suffix) == ("mb","2b") else "sfg-model-" + tag + "-" + suffix
            add("sfg", name, "air-water.lammpstrj",
                {"dipole_mb": dipole, "polarizability_mb": polar},
                extended="2b" in dipole or "2b" in polar,
                covers=["dipole_mb=" + dipole, "polarizability_mb=" + polar])
            if name == "sfg-air-mb":
                cases[-1]["config"]["SFG"].update(include_cross_correlation=True,
                                                cross_correlation_sampling_stride=1)
    add("sfg", "sfg-graphene", "graphene-water.lammpstrj", {"skip_initial_atoms":180},
        covers=["skip_initial_atoms=180", "carbon exclusion"])
    add("power", "power-xyz", "f-water.xyz", covers=["XYZ velocities", "skip_initial_atoms=1"])
    add("power", "power-lammps", "f-water.lammpstrj", covers=["LAMMPS velocities", "skip_initial_atoms=1"])

    sources = {"ir":"bulk-water.xyz", "raman":"bulk-water.xyz",
               "sfg":"air-water.lammpstrj", "power":"f-water.xyz"}
    # Numeric choices are varied one at a time, plus a damping/smoothing pair.
    for suite, filename in sources.items():
        damping_tau = .008
        smoothing_width = 15.0
        sampling_options = {"dump_frequency": .008}
        choices = [
            ("sampling", sampling_options,
             ["dump_frequency=0.008"]),
            ("transform", {"spectrum_transform_time_ps":.012}, ["spectrum_transform_time_ps=0.012"]),
            ("damping", {"tau":damping_tau}, ["tau=" + str(damping_tau)]),
            ("smoothing", {"smoothing_fwhm_cm1":smoothing_width},
             ["smoothing_fwhm_cm1=" + str(smoothing_width)]),
            ("window", {"tau":damping_tau, "smoothing_fwhm_cm1":smoothing_width},
             ["damping and smoothing"]),
        ]
        for tag, changes, coverage in choices:
            if suite == "sfg" and tag == "sampling":
                continue
            if suite in ("ir", "power") and tag in {"damping", "window"}:
                continue
            if suite == "raman" and tag in {"sampling", "damping", "window"}:
                continue
            add(suite, suite + "-" + tag, filename, changes, covers=coverage)
        add(suite, suite + "-duration", filename,
            {"trajectory_duration_ps":.016}, frames=6, covers=["trajectory_duration_ps=0.016"])
        add(suite, suite + "-lags", filename, {"int_max":.020}, frames=6,
            covers=["int_max=0.020"])
        for imcon in ([1,2,3] if suite in ("ir","raman","sfg") else []):
            if suite == "sfg" and imcon == 3:
                comment = (root / "inputs" / "alumina.xyz").read_text().splitlines()[1]
                lattice = [float(value) for value in comment.split('Lattice="', 1)[1].split('"', 1)[0].split()]
                add(suite, "sfg-cell-3", "alumina.lammpstrj",
                    {"imcon":3, "box":lattice, "skip_initial_atoms":800,
                     "dipole_mb":"1b", "polarizability_mb":"1b", "use_ps_dipole":True},
                    covers=["imcon=3", "nonzero triclinic tilt", "skip_initial_atoms=800",
                            "dipole_mb=1b", "polarizability_mb=1b", "use_ps_dipole=true"])
                continue
            if suite in ("ir", "raman") and imcon == 3:
                comment = (root / "inputs" / "alumina.xyz").read_text().splitlines()[1]
                lattice = [float(value) for value in comment.split('Lattice="', 1)[1].split('"', 1)[0].split()]
                settings = {"imcon": 3, "box": lattice, "skip_initial_atoms": 800}
                model_key = "dipole_mb" if suite == "ir" else "polarizability_mb"
                for mode, suffix in [("1b", ""), ("1b+2b", "-2b")]:
                    for extension, format_suffix in [("lammpstrj", ""), ("xyz", "-xyz")]:
                        add(suite, suite + "-cell-3" + format_suffix + suffix, "alumina." + extension,
                            dict(settings, **{model_key: mode}), extended=mode == "1b+2b",
                            cell="per-frame-xyz" if extension == "xyz" else None,
                            covers=["imcon=3", "nonzero triclinic tilt", "skip_initial_atoms=800",
                                    model_key + "=" + mode])
                continue
            source = "bulk-water.lammpstrj" if suite in ("ir","raman") else filename
            options = {"imcon":imcon}
            # A zero-tilt triclinic representation preserves the original cell.
            add(suite, suite + "-cell-" + str(imcon), source, options,
                cell="matrix" if imcon == 3 else None,
                covers=["imcon=" + str(imcon), "box matrix" if imcon == 3 else "box bounds"])
        if suite in ("ir", "sfg"):
            add(suite, suite + "-temperature", filename, {"temperature":600.}, covers=["temperature=600"])
        if suite in ("ir","raman"):
            add(suite, suite + "-skip-water", filename, {"skip_initial_atoms":3}, covers=["skip_initial_atoms=3"])
            add(suite, suite + "-molecules", filename, {"write_molecular_properties":True},
                covers=["write_molecular_properties=true"])
        if suite == "ir":
            add(suite, suite + "-volume", "bulk-water.lammpstrj", {"box":[0,30,0,20,0,20]}, covers=["box changed"])
            add(suite, suite + "-cutoffs", filename,
                {"dipole_mb":"1b+2b"} if suite == "ir" else {"polarizability_mb":"1b+2b"},
                extended=True, covers=["cutoff_1=2.5", "cutoff_2=4.5"])
            cases[-1]["config"][SECTIONS[suite]].update(cutoff_1=2.5, cutoff_2=4.5)

    for value in [False, True]:
        add("ir", "ir-mean-" + str(value).lower(), sources["ir"],
            {"subtract_mean_dipole":value}, covers=["subtract_mean_dipole="+str(value).lower()])
        add("power", "power-mean-" + str(value).lower(), sources["power"],
            {"subtract_mean_velocity":value}, covers=["subtract_mean_velocity="+str(value).lower()])
        for channel in ["isotropic","depolarized"]:
            add("raman", "raman-" + channel + "-mean-" + str(value).lower(), sources["raman"],
                {"raman_response":channel,"subtract_mean_polarizability":value},
                covers=["raman_response="+channel,"subtract_mean_polarizability="+str(value).lower()])
        add("raman", "raman-bose-" + str(value).lower(), sources["raman"],
            {"bose_einstein_correction":value}, covers=["bose_einstein_correction="+str(value).lower()])
        add("power", "power-fluoride-" + ("excluded" if value else "included"), "f-water.lammpstrj",
            {"skip_initial_atoms":int(value)}, covers=["skip_initial_atoms="+str(int(value))])

    add("power", "power-fluoride-included-mean-true", "f-water.lammpstrj",
        {"skip_initial_atoms":0, "subtract_mean_velocity":True},
        covers=["skip_initial_atoms=0", "subtract_mean_velocity=true"])

    for ijk in ["xxz","zzz","xzx","zxx"]:
        for surface in ["positive","negative","both"]:
            add("sfg", "sfg-" + ijk + "-" + surface, sources["sfg"],
                {"ijk":ijk,"surface":surface}, covers=["ijk="+ijk,"surface="+surface])
    for key in ["center_slab","shift_slab_by_half_cell","time_dependent_switch","use_ps_dipole",
                "subtract_initial_response"]:
        for value in [False,True]:
            add("sfg", "sfg-"+key.replace('_','-')+"-"+str(value).lower(), sources["sfg"],
                {key:value}, covers=[key+"="+str(value).lower()])
    for ijk in ["zzz", "xzx", "zxx"]:
        add("sfg", "sfg-ps-" + ijk, sources["sfg"],
            {"ijk":ijk, "use_ps_dipole":True},
            covers=["ijk=" + ijk, "use_ps_dipole=true"])
    for self_terms, neighbors in [(True,False),(False,True),(True,True)]:
        add("sfg", "sfg-correlation-" + ("self" if self_terms else "") + ("neighbors" if neighbors else ""),
            sources["sfg"], {"include_autocorrelation":self_terms,"include_cross_correlation":neighbors},
            covers=["include_autocorrelation="+str(self_terms).lower(),"include_cross_correlation="+str(neighbors).lower()])
    for tag, changes in [
        ("interface", {"zc1":4.,"zc2":9.}),
        ("neighbor-cutoff", {"include_cross_correlation":True,"cross_correlation_cutoff":3.}),
        ("stride-1", {"include_cross_correlation":True,"cross_correlation_sampling_stride":1}),
        ("stride-2", {"include_cross_correlation":True,"cross_correlation_sampling_stride":2}),
        ("range", {"spectrum_max_wavenumber_cm1":3000.}),
        ("grid", {"spectrum_num_wavenumbers":17}),
        ("cutoffs", {"dipole_mb":"1b+2b","cutoff_1":2.5,"cutoff_2":4.5}),
    ]:
        add("sfg", "sfg-"+tag, sources["sfg"], changes, extended=tag=="cutoffs",
            covers=[key+"="+str(value).lower() for key,value in changes.items()])
    unique_cases = []
    seen = {}
    for case in cases:
        key = json.dumps({key: value for key, value in case.items()
                          if key not in ("name", "covers")}, sort_keys=True)
        if key in seen:
            for coverage in case["covers"]:
                if coverage not in seen[key]["covers"]:
                    seen[key]["covers"].append(coverage)
            continue
        seen[key] = case
        unique_cases.append(case)
    cases = unique_cases
    assert len({case["name"] for case in cases}) == len(cases)
    return cases


def stage_options(suite, stage, options):
    """Only settings actually read by this stage belong in its reuse key."""
    property_keys = ["box","imcon","dump_frequency","skip_initial_atoms","trajectory_duration_ps",
                     "write_molecular_properties","cutoff_1","cutoff_2"]
    if stage == "dipole":
        keys = property_keys + ["DMS","dipole_mb"]
    elif stage == "polarizability":
        keys = property_keys + ["polarizability_mb"]
    elif stage == "response":
        if suite in ("ir","raman"):
            keys = ["int_max","trajectory_duration_ps","subtract_mean_dipole"] if suite == "ir" else \
                   ["int_max","trajectory_duration_ps","subtract_mean_polarizability","raman_response"]
        elif suite == "power":
            keys = ["box","imcon","dump_frequency","skip_initial_atoms","trajectory_duration_ps",
                    "int_max","subtract_mean_velocity"]
        else:
            keys = [key for key in options if key not in ["temperature","tau","smoothing_fwhm_cm1",
                    "spectrum_transform_time_ps","subtract_initial_response",
                    "spectrum_max_wavenumber_cm1","spectrum_num_wavenumbers"]]
    else:
        keys = ["int_max","spectrum_transform_time_ps","tau","smoothing_fwhm_cm1"]
        if suite in ("ir","raman"):
            keys += ["box","imcon","temperature"]
        if suite == "raman": keys += ["bose_einstein_correction"]
        if suite == "sfg": keys += ["temperature","subtract_initial_response",
                                    "spectrum_max_wavenumber_cm1","spectrum_num_wavenumbers"]
    return {key:copy.deepcopy(options[key]) for key in keys if key in options}

