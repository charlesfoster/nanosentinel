"""Build a fake MinKNOW run (fastq_pass/barcode01..02) from any ONT .fastq.gz, for trying nanosentinel.

Usage: python dev/make_fake_run.py reads.fastq.gz OUT_DIR   (then use --minknow_data OUT_DIR --run_name exp1)
"""
import gzip, os, sys, datetime
src, root = sys.argv[1], sys.argv[2]
run = os.path.join(root, "exp1", "sample1", "20261109_1000_P2S-00001-A_FLOabc_deadbeef")
t0 = datetime.datetime(2026, 11, 9, 10, 0, 0)
recs = []
with gzip.open(src, "rt") as fh:
    while True:
        h = fh.readline()
        if not h: break
        recs.append((h.rstrip().split(" ")[0], fh.readline(), fh.readline(), fh.readline()))
n = min(len(recs), 3000)
for i, (name, seq, plus, qual) in enumerate(recs[:n]):
    bc = "barcode01" if i % 2 == 0 else "barcode02"
    chunk = i // 500
    d = os.path.join(run, "fastq_pass", bc); os.makedirs(d, exist_ok=True)
    st = (t0 + datetime.timedelta(seconds=i * 2)).strftime("%Y-%m-%dT%H:%M:%SZ")
    f = os.path.join(d, f"FLOabc_pass_{bc}_deadbeef_{chunk}.fastq.gz")
    with gzip.open(f, "at") as out:
        out.write(f"{name} runid=deadbeef read={i} ch={i%512} start_time={st} barcode={bc}\n{seq}+\n{qual}")
print(run)
