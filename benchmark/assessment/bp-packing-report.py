# The packed layout against the gather. Every figure is a median over the
# repeats of one case, which is the mistake the previous analysis made: a dict
# keyed by case kept whichever repeat came last.
import csv, glob, os, sys, statistics, collections

rows = collections.defaultdict(list)   # (scene, bp, ver, case, type) -> [prep+broad]
for f in sorted(glob.glob(os.path.join(sys.argv[1], "*.csv"))):
    scene, bp, ver, _ = os.path.basename(f).split(".")
    with open(f) as fh:
        for r in csv.DictReader(fh):
            rows[(scene, bp, ver, r["case"], r["type"])].append(
                float(r["prep_ms"]) + float(r["broad_ms"]))

med = {k: statistics.median(v) for k, v in rows.items()}
scenes = sorted({k[0] for k in med})
print("Whole broad phase (prep + queries), median over 3 repeats per case, ms")
print("%-20s %-10s %-4s %5s %10s %10s %8s" % (
    "scene", "strategy", "type", "cases", "gather", "packed", "ratio"))
tot = collections.defaultdict(lambda: [0.0, 0.0])
for scene in scenes:
    for bp in ("cell2dmin", "cell2d"):
        for ty in ("ee", "vf"):
            a = {k[3]: v for k, v in med.items() if k[:3] == (scene, bp, "seg") and k[4] == ty}
            b = {k[3]: v for k, v in med.items() if k[:3] == (scene, bp, "pack") and k[4] == ty}
            common = sorted(set(a) & set(b))
            if not common:
                continue
            A, B = sum(a[c] for c in common), sum(b[c] for c in common)
            tot[(bp, ty)][0] += A
            tot[(bp, ty)][1] += B
            print("%-20s %-10s %-4s %5d %10.1f %10.1f %7.3fx" % (
                scene, bp, ty, len(common), A, B, A / B if B else 0))
print()
for (bp, ty), (A, B) in sorted(tot.items()):
    print("%-20s %-10s %-4s %5s %10.1f %10.1f %7.3fx" % ("ALL SCENES", bp, ty, "", A, B, A / B if B else 0))
