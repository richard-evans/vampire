#!/usr/bin/env python3
#------------------------------------------------------------------------------
#
#   This file is part of the VAMPIRE open source package under the
#   Free BSD licence (see licence file for details).
#
#------------------------------------------------------------------------------
"""
VAMPIRE validation suite: runs reference cases with any vampire binary and
compares the results against analytic solutions or reference outputs.

   vtest.py run --exe ./vampire-cuda --label gpu-dp [--precision sp] [--cases a,b]
                [--extra "gpu:acceleration=off"] [--seed-offset N] [--seeds N] [--jobs N] [--reference]
   vtest.py compare --test gpu-dp --ref cpu [--cases a,b] [--verbose] [--summary] [--paired RTOL]
   vtest.py identical --test A --ref B [--cases a,b]
   vtest.py configs --test A --ref B [--cases a,b] [--magnetic-only]

See README.md for case.json format and tolerance policy.
"""
import argparse, concurrent.futures, datetime, json, os, re, shutil, socket, subprocess, sys
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
CASES = os.path.join(ROOT, "cases")
RESULTS = os.path.join(ROOT, "results")      # scratch results (git-ignored)
REFERENCE = os.path.join(ROOT, "reference")  # committed reference outputs

KB = 1.380649e-23        # Boltzmann constant (J/K)
MU_B = 9.2740100783e-24  # Bohr magneton (J/T)
GAMMA = 1.76085963023e11 # gyromagnetic ratio used by vampire (1/(T s))

#------------------------------------------------------------------------------
# helpers
#------------------------------------------------------------------------------
def load_case(name):
   with open(os.path.join(CASES, name, "case.json")) as f:
      return json.load(f)

def case_names(selection):
   names = sorted(d for d in os.listdir(CASES) if os.path.isfile(os.path.join(CASES, d, "case.json")))
   if not selection: return [n for n in names if load_case(n)["kind"] != "benchmark"]   # benchmarks only on request
   else:
      wanted = selection.split(",")
      missing = [w for w in wanted if w not in names]
      if missing: sys.exit("Unknown case(s): " + ", ".join(missing))
      names = wanted
   return names

def run_names(case, seeds=None):
   # seeded cases get one run per seed (default: the case's own seed count), deterministic cases a single run
   if case["kind"] in ("deterministic", "benchmark") or case.get("model") == "precession": return ["run"]
   return ["seed%d" % s for s in range(1, (seeds or case["seeds"]) + 1)]

def seed_runs(case, d):
   # all seeds present in a result directory (the case's own seed count, or more when run with --seeds)
   runs = [r for r in os.listdir(d) if r.startswith("seed") and r[4:].isdigit()]
   return sorted(runs, key=lambda r: int(r[4:])) if len(runs) > case["seeds"] else run_names(case)

def read_output(path):
   """Returns dict column name -> numpy array (requires output:column-headers)."""
   with open(path) as f:
      lines = [l for l in f if l.strip() and not l.startswith("#")]
   names = lines[0].split()
   data = np.array([[float(v) for v in l.split()] for l in lines[1:]], ndmin=2)
   if data.shape[1] != len(names): sys.exit("Column header mismatch in " + path)
   return {n: data[:, i] for i, n in enumerate(names)}

def quantities(case, path):
   cols = read_output(path)
   q = {q: np.asarray(eval(expr, {"__builtins__": {}}, cols), dtype=float) for q, expr in case["quantities"].items()}
   if case.get("grain"):
      # grain.txt (no column headers): time, temperature, then the grain columns, one row per output
      q["grain"] = np.loadtxt(os.path.join(os.path.dirname(path), "grain.txt"), ndmin=2)[:, 2:]
   return q

def result_dir(label):
   # look in scratch results first, then committed references
   for base in (RESULTS, REFERENCE):
      if os.path.isdir(os.path.join(base, label)): return os.path.join(base, label)
   sys.exit("No results found for label '%s'" % label)

#------------------------------------------------------------------------------
# run
#------------------------------------------------------------------------------
def run_one(cmd, case_name, run, out_dir, extra, seed_offset):
   src = os.path.join(CASES, case_name)
   if os.path.isdir(out_dir): shutil.rmtree(out_dir)
   shutil.copytree(src, out_dir, ignore=shutil.ignore_patterns("case.json"))
   with open(os.path.join(out_dir, "input"), "a") as f:
      f.write("\n")   # input files may lack a trailing newline
      for line in extra: f.write(line + "\n")
      if run.startswith("seed"): f.write("sim:integrator-random-seed = %d\n" % (int(run[4:]) + seed_offset))
   with open(os.path.join(out_dir, "stdout"), "w") as so:
      rc = subprocess.call(cmd, cwd=out_dir, stdout=so, stderr=subprocess.STDOUT)
   ok = rc == 0 and os.path.isfile(os.path.join(out_dir, "output"))
   return case_name, run, ok

def cmd_run(args):
   exe = os.path.abspath(args.exe)
   base = os.path.join(REFERENCE if args.reference else RESULTS, args.label)
   names = case_names(args.cases)
   version = subprocess.run([exe, "--version"], capture_output=True, text=True).stdout.strip()
   # MPI versions are started with mpirun (command from the MPIRUN environment variable)
   cmd = os.environ.get("MPIRUN", "mpirun").split() + ["-np", str(args.ranks), exe] if args.ranks else [exe]
   gpu = ""
   try: gpu = subprocess.run(["nvidia-smi", "--query-gpu=name,driver_version", "--format=csv,noheader"], capture_output=True, text=True).stdout.strip()
   except OSError: pass

   jobs = [(c, r) for c in names for r in run_names(load_case(c), args.seeds)]
   failed = []
   # runs of an earlier set with more seeds would otherwise remain and be compared
   for c in names:
      shutil.rmtree(os.path.join(base, c), ignore_errors=True)
      os.makedirs(os.path.join(base, c))   # before threads, avoids makedirs races
   with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
      futures = [pool.submit(run_one, cmd, c, r, os.path.join(base, c, r), args.extra, args.seed_offset) for c, r in jobs]
      for fut in concurrent.futures.as_completed(futures):
         c, r, ok = fut.result()
         print("%-30s %-8s %s" % (c, r, "done" if ok else "FAILED (see stdout in run directory)"), flush=True)
         if not ok: failed.append((c, r))

   if args.reference:
      # keep only what comparisons need so committed references stay small
      for c, r in jobs:
         d = os.path.join(base, c, r)
         for f in os.listdir(d):
            p = os.path.join(d, f)
            if f not in ("output", "grain.txt"):
               shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
            else:
               # the host name and directory in the file header describe the machine, not the result
               with open(p, errors="replace") as fh: lines = fh.readlines()
               with open(p, "w") as fh: fh.writelines(l for l in lines if not re.match(r"#\s+(hostname|path)\s+:", l))

   meta = {"label": args.label, "exe": os.path.basename(exe) if args.reference else exe, "ranks": args.ranks, "version": version,
           "precision": args.precision, "extra": args.extra, "seed_offset": args.seed_offset,
           "date": datetime.datetime.now().isoformat(timespec="seconds"), "gpu": gpu}
   if not args.reference: meta["host"] = socket.gethostname()
   with open(os.path.join(base, "meta.json"), "w") as f: json.dump(meta, f, indent=3)
   return 1 if failed else 0

#------------------------------------------------------------------------------
# analytic models
#------------------------------------------------------------------------------
def precession(t, p):
   # single spin initially along x, field H along z (Landau-Lifshitz-Gilbert, damping alpha)
   arg = t * GAMMA * p["H"] / (1.0 + p["alpha"]**2)
   larg = p["alpha"] * arg
   return {"mx": np.cos(arg) / np.cosh(larg), "my": np.sin(arg) / np.cosh(larg), "mz": np.tanh(larg)}

def thermal_mz(p):
   # <cos theta> for non-interacting spins, field and uniaxial anisotropy along z (Boltzmann average);
   # fourth and sixth order terms as in vampire (anisotropy/uniaxial_order_4.cpp, _6.cpp)
   c = np.linspace(-1.0, 1.0, 400001)
   e = (p["mu_s"] * MU_B * p["H"] * c + p["ku"] * c * c + p.get("ku4", 0.0) * (c**4 - 6.0 / 7.0 * c**2)
        + p.get("ku6", 0.0) * (c**6 - 15.0 / 11.0 * c**4 + 5.0 / 11.0 * c**2)) / (KB * p["T"])
   w = np.exp(e - e.max())
   return float((c * w).sum() / w.sum())

#------------------------------------------------------------------------------
# compare
#------------------------------------------------------------------------------
def seed_means(case, d):
   # time average (or per-row value) of each quantity per seed, then mean and standard error over seeds
   per = {}
   for r in seed_runs(case, d):
      for q, v in quantities(case, os.path.join(d, r, "output")).items():
         per.setdefault(q, []).append(v if case.get("per_row") else v[case["skip_rows"]:].mean())
   return {q: (np.mean(v, axis=0), np.std(v, ddof=1, axis=0) / np.sqrt(len(v))) for q, v in per.items()}

def noise_ok(se_test, se_ref):
   # a statistical match is meaningless if the test results are non-finite or far noisier than the reference
   se_test, se_ref = np.atleast_1d(se_test), np.atleast_1d(se_ref)
   return bool(np.all(np.isfinite(se_test)) and np.all(se_test <= 10.0 * se_ref + 1e-300))

Z_SCORES = []   # (case, quantity, signed z, noise ratio) of every statistical comparison, for --summary

def compare_paired(case, test_dir, ref_dir, rtol):
   """Every run of the test set against the same run of the reference set, row by row:
   abs(test-ref) <= rtol*max(abs(ref)) + rtol*abs(ref). For runs with the same dynamics
   whose output differs only by rounding (e.g. the order of a summation)."""
   runs = sorted(r for r in os.listdir(test_dir) if os.path.isfile(os.path.join(test_dir, r, "output")))
   if not runs: return [("-", False, "no output: run failed (see stdout in the run directory)")]
   worst = {}
   for run in runs:
      if not os.path.isfile(os.path.join(ref_dir, run, "output")): return [("-", False, "no reference output for " + run)]
      t, r = quantities(case, os.path.join(test_dir, run, "output")), quantities(case, os.path.join(ref_dir, run, "output"))
      for k in r:
         if t[k].shape != r[k].shape: return [(k, False, "%s: row count %d vs %d" % (run, len(t[k]), len(r[k])))]
         rel = np.abs(t[k] - r[k]) / np.maximum(np.abs(r[k]).max() + np.abs(r[k]), 1e-300)
         worst[k] = max(worst.get(k, 0.0), float(rel.max()))
   return [(k, w <= rtol, "identical" if w == 0 else "max scaled diff %.3e over %d run(s) (tol %.0e)" % (w, len(runs), rtol)) for k, w in worst.items()]

def compare_case(name, case, test_dir, ref_dir, precision, paired=None):
   """Returns list of (quantity, passed, description)."""
   if paired is not None and case["kind"] != "benchmark": return compare_paired(case, test_dir, ref_dir, paired)
   rows = []
   tol = case.get("tolerance", {}).get(precision, {})
   missing = [r for r in run_names(case) if not os.path.isfile(os.path.join(test_dir, r, "output"))]
   if missing:
      return [("-", False, "no output for %s: run failed (see stdout in the run directory)" % ", ".join(missing))]
   if case["kind"] == "analytic" and case["model"] == "precession":
      cols = read_output(os.path.join(test_dir, "run", "output"))
      q = quantities(case, os.path.join(test_dir, "run", "output"))
      exact = precession(cols["Real_time"], case["params"])
      for k in ("mx", "my", "mz"):
         err = np.abs(q[k] - exact[k]).max()
         rows.append((k, err <= tol["atol"], "max |err| = %.3e (tol %.1e)" % (err, tol["atol"])))
   elif case["kind"] == "analytic" and case["model"] == "thermal":
      exact = thermal_mz(case["params"])
      m, se = seed_means(case, test_dir)["mz"]
      mr, sr = seed_means(case, ref_dir)["mz"]
      Z_SCORES.append((name, "mz", np.atleast_1d((m - mr) / max(np.hypot(se, sr), 1e-300)), np.atleast_1d(se / max(sr, 1e-300))))
      lim = case["k"] * se + tol["atol"]
      ok = abs(m - exact) <= lim and noise_ok(se, sr)
      rows.append(("mz", ok, "%.5f +/- %.5f vs exact %.5f (limit %.1e, reference +/- %.5f)" % (m, se, exact, lim, sr)))
   elif case["kind"] == "deterministic":
      t = quantities(case, os.path.join(test_dir, "run", "output"))
      r = quantities(case, os.path.join(ref_dir, "run", "output"))
      n = tol.get("rows")   # optionally compare only the first n rows (before round-off amplification dominates)
      if n: t, r = ({k: v[:n] for k, v in x.items()} for x in (t, r))
      for k in r:
         if t[k].shape != r[k].shape:
            rows.append((k, False, "row count %d vs %d" % (len(t[k]), len(r[k])))); continue
         # atol is relative to the largest reference magnitude, so it applies equally to unit vectors and energies (J);
         # the optional absolute floor "abs" covers quantities that are close to zero in all compared rows
         diff = np.abs(t[k] - r[k])
         ok = bool(np.all(diff <= tol.get("abs", 0.0) + tol["atol"] * np.abs(r[k]).max() + tol["rtol"] * np.abs(r[k])))
         desc = "identical" if not diff.any() else "max |diff| = %.3e, max rel = %.3e" % (diff.max(), (diff / np.maximum(np.abs(r[k]), 1e-300)).max())
         if n: desc += " (first %d rows)" % n
         rows.append((k, ok, desc))
   elif case["kind"] == "benchmark":
      with open(os.path.join(test_dir, "run", "stdout")) as f:
         times = [l.split()[-1] for l in f if l.startswith("Simulation run time")]
      rows.append(("time", bool(times), "simulation run time %s s" % (times[0] if times else "not found")))
   elif case["kind"] == "statistical":
      t, r = seed_means(case, test_dir), seed_means(case, ref_dir)
      for k in r:
         (mt, st), (mr, sr) = (np.atleast_1d(x).ravel() for x in t[k]), (np.atleast_1d(x).ravel() for x in r[k])
         if mt.shape != mr.shape:
            rows.append((k, False, "row count %d vs %d" % (len(mt), len(mr)))); continue
         Z_SCORES.append((name, k, (mt - mr) / np.maximum(np.hypot(st, sr), 1e-300), st / np.maximum(sr, 1e-300)))
         z = np.abs(mt - mr) / np.maximum(np.hypot(st, sr), 1e-300)
         i = int(np.argmax(z))   # report the worst row
         desc = "%.6g +/- %.2g vs ref %.6g +/- %.2g" % (mt[i], st[i], mr[i], sr[i])
         if len(mt) > 1: desc = "worst of %d rows (row %d): " % (len(mt), i) + desc
         rows.append((k, bool(z[i] <= case["k"]) and noise_ok(st, sr), desc))
   return rows

def cmd_compare(args):
   test_base, ref_base = result_dir(args.test), result_dir(args.ref)
   with open(os.path.join(test_base, "meta.json")) as f: precision = json.load(f)["precision"]
   failures = 0
   print("Comparing '%s' (%s) against '%s'" % (args.test, precision, args.ref))
   for name in case_names(args.cases):
      case = load_case(name)
      td, rd = os.path.join(test_base, name), os.path.join(ref_base, name)
      if not os.path.isdir(td):
         print("%-30s MISSING results" % name); failures += 1; continue
      try:
         rows = compare_case(name, case, td, rd, precision, args.paired)
      except (OSError, ValueError, KeyError) as e:
         rows = [("-", False, "error: %s" % e)]
      if not args.verbose and all(ok for _, ok, _ in rows):
         same = all(d == "identical" for _, _, d in rows)
         print("%-30s %-12s OK   %d quantities%s" % (name, case["kind"], len(rows), ", all identical" if same else ""))
         continue
      for q, ok, desc in rows:
         print("%-30s %-12s %-4s %s" % (name, q, "OK" if ok else "FAIL", desc))
         failures += 0 if ok else 1
   if args.summary and Z_SCORES: print_summary()
   print("%s: %d failure(s)" % ("OVERALL PASS" if failures == 0 else "OVERALL FAIL", failures))
   return 1 if failures else 0

def print_summary():
   # If test and reference sample the same distribution, the signed z-scores of all statistical
   # comparisons are approximately standard normal: mean ~0, rms ~1, ~4.6 % beyond 2 and ~0.3 % beyond 3,
   # and the ratio of test to reference standard errors is ~1 (for equal numbers of seeds).
   print("Summary of statistical comparisons (signed z = (test - ref) / combined standard error):")
   print("   %-26s %-14s %5s %7s %6s %7s %11s" % ("case", "quantity", "n", "mean z", "rms z", "max |z|", "noise ratio"))
   for name, q, z, ratio in Z_SCORES:
      print("   %-26s %-14s %5d %7.2f %6.2f %7.2f %11.2f" % (name, q, len(z), z.mean(), np.sqrt((z**2).mean()), np.abs(z).max(), np.median(ratio)))
   z = np.concatenate([x[2] for x in Z_SCORES]); ratio = np.concatenate([x[3] for x in Z_SCORES])
   print("   %-41s %5d %7.2f %6.2f %7.2f %11.2f" % ("all", len(z), z.mean(), np.sqrt((z**2).mean()), np.abs(z).max(), np.median(ratio)))
   print("   fraction |z| > 2: %.1f %% (expected 4.6 %%), |z| > 3: %.1f %% (expected 0.3 %%)" % (100 * np.mean(np.abs(z) > 2), 100 * np.mean(np.abs(z) > 3)))

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# configs: final spin configurations of two result sets, atom by atom
#------------------------------------------------------------------------------
def read_config(d):
   # binary config files (config:output-format=binary): atoms-coords.data holds
   # n, n types, n categories, 3n coordinates; spins-<counter>.data holds n, 3n spins;
   # with config:output-mode=file-per-node each output group writes its own files,
   # atoms-coords-<group>.data and spins-<counter>-<group>.data, read here in group order
   def read(path, offset):
      with open(path, "rb") as f:
         n = int(np.frombuffer(f.read(8), dtype=np.uint64)[0])
         f.seek(8 + offset * n)
         return n, np.frombuffer(f.read(24 * n), dtype=np.float64).reshape(n, 3)
   files = sorted(os.listdir(d))
   coord_files = [f for f in files if f.startswith("atoms-coords") and f.endswith(".data")]
   counter = sorted({f[6:14] for f in files if f.startswith("spins-") and f.endswith(".data")})[-1]
   spin_files = [f for f in files if f.startswith("spins-" + counter) and f.endswith(".data")]
   if len(coord_files) != len(spin_files): sys.exit("Coordinate and spin files of %s have different groups" % d)
   coords = [read(os.path.join(d, f), 8) for f in coord_files]
   spins = [read(os.path.join(d, f), 0) for f in spin_files]
   if [n for n, _ in coords] != [n for n, _ in spins]: sys.exit("Coordinate and spin files of %s have different atom counts" % d)
   coords, spins = np.concatenate([c for _, c in coords]), np.concatenate([s for _, s in spins])
   order = np.lexsort((coords[:, 2], coords[:, 1], coords[:, 0]))   # files list atoms in rank order
   return coords[order], spins[order], "spins-" + counter

def keep_materials(case):
   # 0-based indices of the case's non-magnetic=keep materials (as in the configuration files)
   keep = set()
   for f in os.listdir(os.path.join(CASES, case)):
      if not f.endswith(".mat"): continue
      for l in open(os.path.join(CASES, case, f)):
         m = re.match(r"\s*material\[(\d+)\]:non-magnetic\s*=\s*keep\b", l)
         if m: keep.add(int(m.group(1)) - 1)
   return keep

def read_types(d):
   # material of every atom, in the order of read_config
   files = sorted(f for f in os.listdir(d) if f.startswith("atoms-coords") and f.endswith(".data"))
   types, coords = [], []
   for f in files:
      with open(os.path.join(d, f), "rb") as fh:
         n = int(np.frombuffer(fh.read(8), dtype=np.uint64)[0])
         types.append(np.frombuffer(fh.read(4 * n), dtype=np.int32))
         fh.seek(8 + 8 * n)
         coords.append(np.frombuffer(fh.read(24 * n), dtype=np.float64).reshape(n, 3))
   coords = np.concatenate(coords)
   return np.concatenate(types)[np.lexsort((coords[:, 2], coords[:, 1], coords[:, 0]))]

def cmd_configs(args):
   test_dir, ref_dir = result_dir(args.test), result_dir(args.ref)
   failures = 0
   for c in case_names(args.cases):
      keep = keep_materials(c) if args.magnetic_only else set()
      for r in sorted(os.listdir(os.path.join(test_dir, c))) if os.path.isdir(os.path.join(test_dir, c)) else []:
         t, f = os.path.join(test_dir, c, r), os.path.join(ref_dir, c, r)
         if not os.path.isdir(t): continue
         missing = [d for d in (t, f) if not any(x.startswith("atoms-coords") and x.endswith(".data") for x in os.listdir(d))] if os.path.isdir(f) else [f]
         if missing:
            print("%-30s %-8s FAIL no binary configuration in %s" % (c, r, ", ".join(missing))); failures += 1; continue
         tc, ts, tname = read_config(t)
         rc, rs, rname = read_config(f)
         if tc.shape != rc.shape or not np.array_equal(tc, rc):
            print("%-30s %-8s FAIL atoms differ (%d vs %d atoms or coordinates)" % (c, r, len(tc), len(rc))); failures += 1; continue
         same = np.all(ts.view(np.uint64) == rs.view(np.uint64), axis=1)
         if keep:
            # only the magnetic atoms: atoms of non-magnetic=keep materials may differ
            magnetic = ~np.isin(read_types(t), sorted(keep))
            ts, rs, same = ts[magnetic], rs[magnetic], same[magnetic]
         status = "IDENTICAL" if same.all() else "DIFFERENT"
         if not same.all(): failures += 1
         print("%-30s %-8s %-9s %s: %d %satoms, %d bit-identical spins, max |diff| = %.3e" % (c, r, status, tname, len(ts), "magnetic " if keep else "", same.sum(), np.abs(ts - rs).max()))
   print("OVERALL", "PASS" if failures == 0 else "FAIL: %d run(s) differ" % failures)
   return 1 if failures else 0

#------------------------------------------------------------------------------
# identical: output and grain.txt of two result sets, byte for byte (comment lines ignored)
#------------------------------------------------------------------------------
def data_lines(path):
   with open(path) as f: return [l for l in f if not l.startswith("#")]

def cmd_identical(args):
   test_dir, ref_dir = result_dir(args.test), result_dir(args.ref)
   failures = 0
   for c in case_names(args.cases):
      if not os.path.isdir(os.path.join(test_dir, c)): continue
      for r in sorted(os.listdir(os.path.join(test_dir, c))):
         t, f = os.path.join(test_dir, c, r), os.path.join(ref_dir, c, r)
         if not os.path.isdir(t): continue
         files = [x for x in ("output", "grain.txt") if os.path.isfile(os.path.join(f, x))]
         missing = [x for x in files if not os.path.isfile(os.path.join(t, x))]
         differ = [x for x in files if x not in missing and data_lines(os.path.join(t, x)) != data_lines(os.path.join(f, x))]
         ok = files and not missing and not differ
         if not ok: failures += 1
         detail = ", ".join(["%s missing" % x for x in missing] + ["%s differs" % x for x in differ]) if files else "no output in reference"
         print("%-30s %-8s %-9s %s" % (c, r, "IDENTICAL" if ok else "DIFFERENT", ", ".join(files) if ok else detail))
   print("OVERALL", "PASS" if failures == 0 else "FAIL: %d run(s) differ" % failures)
   return 1 if failures else 0

def main():
   p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
   sub = p.add_subparsers(dest="cmd", required=True)
   r = sub.add_parser("run")
   r.add_argument("--exe", required=True)
   r.add_argument("--label", required=True)
   r.add_argument("--precision", choices=["dp", "sp"], default="dp")
   r.add_argument("--cases")
   r.add_argument("--extra", action="append", default=[], help="extra input line (repeatable)")
   r.add_argument("--seed-offset", type=int, default=0, help="shift seeds, e.g. to re-run a statistical failure with fresh seeds")
   r.add_argument("--seeds", type=int, help="number of seeds for statistical cases (default: each case's own)")
   r.add_argument("--jobs", type=int, default=1, help="concurrent runs (keep 1 for GPU)")
   r.add_argument("--ranks", type=int, help="run with mpirun -np RANKS (MPI versions)")
   r.add_argument("--reference", action="store_true", help="store as committed reference")
   c = sub.add_parser("compare")
   c.add_argument("--test", required=True)
   c.add_argument("--ref", required=True)
   c.add_argument("--cases")
   c.add_argument("--verbose", action="store_true", help="show every quantity, not just failures")
   c.add_argument("--summary", action="store_true", help="pooled z-score statistics of all statistical comparisons")
   c.add_argument("--paired", type=float, metavar="RTOL", help="compare every run with the same run of the reference row by row (same dynamics, rounding differences only)")
   k = sub.add_parser("configs", help="compare the final binary spin configurations of two result sets atom by atom (bit for bit)")
   k.add_argument("--test", required=True)
   k.add_argument("--ref", required=True)
   k.add_argument("--cases")
   k.add_argument("--magnetic-only", action="store_true", help="ignore atoms of non-magnetic=keep materials (read from the case's .mat files)")
   i = sub.add_parser("identical", help="compare output and grain.txt of two result sets byte for byte, ignoring comment lines")
   i.add_argument("--test", required=True)
   i.add_argument("--ref", required=True)
   i.add_argument("--cases")
   args = p.parse_args()
   sys.exit({"run": cmd_run, "compare": cmd_compare, "configs": cmd_configs, "identical": cmd_identical}[args.cmd](args))

if __name__ == "__main__":
   main()
