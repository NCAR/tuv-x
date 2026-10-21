#!/usr/bin/env python3
"""Build and run the sibling TUV-x reference with either Homebrew toolchain.

Paths can be overridden for other installations; no reference dependency is
required to build the C library. Logs and copied runtime data stay in build/.
"""
import argparse
import os
from pathlib import Path
import subprocess
import shutil

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("toolchain", choices=["gcc", "llvm"])
    p.add_argument("--source", type=Path, default=ROOT.parent / "TUV-x")
    p.add_argument("--prefix", type=Path)
    p.add_argument("--netcdf", type=Path)
    p.add_argument("--configure-only", action="store_true")
    p.add_argument("--jobs", type=int, default=8)
    a = p.parse_args()
    home = Path.home()
    brew = Path(subprocess.check_output(["brew", "--prefix"], text=True).strip())
    llvm = a.toolchain == "llvm"
    version = "llvm-23.1.3" if llvm else "gcc-16.2.0"
    prefix = a.prefix or home / "software" / ("chempas-" + version) / "install"
    netcdf = a.netcdf or home / "miniforge3/envs/ceresmip-laptop"
    build = ROOT / "build" / ("reference-" + a.toolchain)
    build.mkdir(parents=True, exist_ok=True)
    # Keep portability corrections isolated from the upstream checkout and
    # its numerical library. LLVM needs the same legacy real-kind promotion
    # as GNU; BSD sed does not implement GNU's bare -i option.
    source = ROOT / "build/reference-source"
    if not source.exists():
        shutil.copytree(a.source, source, ignore=shutil.ignore_patterns(".git", "build*"))
        old = source / "test/oldtuv/CMakeLists.txt"
        old.write_text(old.read_text().replace('MATCHES "GNU"', 'MATCHES "GNU|LLVMFlang"'))
        for script in (source / "test/regression").rglob("*.sh"):
            lines = []
            for line in script.read_text().splitlines():
                if line.startswith("sed ") and " -i " in line:
                    line = line.replace(" -i ", " ")
                    target = line.split()[-1]
                    line += f' > {target}.tmp && mv {target}.tmp {target}'
                lines.append(line)
            script.write_text("\n".join(lines) + "\n")
    # Reopening an already connected input file has processor-dependent
    # positioning. Make the legacy reader's intended restart explicit.
    legacy = source / "test/oldtuv/rxn.f"
    old = "      OPEN(UNIT=kin,FILE='odat/DATAJ1/CH2O/CH2O_jpl11.abs'\n     $     ,STATUS='old')\n      do i = 1, 4"
    new = old.replace("\n      do i", "\n      REWIND(kin)\n      do i")
    legacy.write_text(legacy.read_text().replace(old, new))
    # Darwin/arm64 Flang does not use the C ABI for these large structs passed
    # and returned by value. A test-only pointer shim preserves their contents.
    bridge = source / "test/regression/solvers/delta_eddington.F90"
    text = (a.source / "test/regression/solvers/delta_eddington.F90").read_text()
    text = text.replace("function run_delta_eddington_solver_c( input ) bind(c, name='RunDeltaEddingtonSolver')",
                        "subroutine run_delta_eddington_solver_c( input, output ) bind(c, name='RunDeltaEddingtonSolverPointer')")
    text = text.replace("type(solver_input_t_c), value :: input", "type(solver_input_t_c), intent(in) :: input")
    text = text.replace("type(solver_output_t_c) run_delta_eddington_solver_c", "type(solver_output_t_c), intent(out) :: output")
    text = text.replace("end function run_delta_eddington_solver_c", "end subroutine run_delta_eddington_solver_c")
    text = text.replace("bind(c, name='FreeOutput')", "bind(c, name='FreeOutputPointer')")
    text = text.replace("type(solver_output_t_c), value :: output", "type(solver_output_t_c), intent(in) :: output")
    text = text.replace("output = run_delta_eddington_solver_c( input )", "call run_delta_eddington_solver_c( input, output )")
    bridge.write_text(text)
    cpp = source / "test/regression/solvers/delta_eddington.cpp"
    text = (a.source / "test/regression/solvers/delta_eddington.cpp").read_text()
    text += '\nextern "C" void RunDeltaEddingtonSolverPointer(const SolverInput* input, SolverOutput* output) {\n'
    text += '  *output = RunDeltaEddingtonSolver(*input);\n}\n'
    text += 'extern "C" void FreeOutputPointer(const SolverOutput* output) { FreeOutput(*output); }\n'
    cpp.write_text(text)
    env = dict(os.environ)
    env["PKG_CONFIG_PATH"] = os.pathsep.join(map(str, [prefix / "lib/pkgconfig", netcdf / "lib/pkgconfig"]))
    env["PATH"] = str(netcdf / "bin") + os.pathsep + env["PATH"]
    env["PYTHONPATH"] = str(ROOT / "build/python") + os.pathsep + env.get("PYTHONPATH", "")
    cc = brew / ("opt/llvm/bin/clang" if llvm else "bin/gcc-16")
    cxx = brew / ("opt/llvm/bin/clang++" if llvm else "bin/g++-16")
    fc = brew / ("bin/flang" if llvm else "bin/gfortran-16")
    yaml = prefix.parent / "build/musica/_deps/yaml-cpp-src"
    command = ["cmake", "-S", str(source), "-B", str(build), "-G", "Ninja",
               f"-DCMAKE_C_COMPILER={cc}", f"-DCMAKE_CXX_COMPILER={cxx}",
               f"-DCMAKE_Fortran_COMPILER={fc}", "-DCMAKE_BUILD_TYPE=Release",
               "-DCMAKE_C_FLAGS=-ffp-contract=off", "-DCMAKE_CXX_FLAGS=-ffp-contract=off",
               "-DCMAKE_Fortran_FLAGS=-ffp-contract=off", "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"]
    if yaml.exists():
        command.append(f"-DFETCHCONTENT_SOURCE_DIR_YAML-CPP={yaml}")
    for stage, cmd in [("configure", command),
                       ("build", ["cmake", "--build", str(build), "-j", str(a.jobs)]),
                       ("test", ["ctest", "--test-dir", str(build), "--output-on-failure", "-j", "1"])]:
        if a.configure_only and stage != "configure":
            break
        log = build / (stage + ".log")
        print(f"{a.toolchain}: {stage} ({log})", flush=True)
        with log.open("w") as out:
            result = subprocess.run(cmd, env=env, stdout=out, stderr=subprocess.STDOUT)
        if result.returncode:
            print(log.read_text()[-16000:], flush=True)
            raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
