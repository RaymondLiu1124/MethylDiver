"""Durable checkpoints for chunked CLI runs. Only the parent writes the database."""
import hashlib
import json
from pathlib import Path
import sqlite3
import sys


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


class Checkpoint:
    def __init__(self, output, resume_from, targets):
        self.resuming = resume_from is not None
        self.path = Path(resume_from or (str(output) + ".checkpoint")).absolute()
        self.targets = targets
        self.lock = self.db = None

    def __enter__(self):
        return self

    def _open(self):
        if self.db is not None:
            return
        if self.resuming:
            if not all((self.path / name).is_file() for name in ("owner.sqlite", "chunks.sqlite")):
                raise ValueError(f"No MethylDiver checkpoint at {self.path}.")
        else:
            try:
                self.path.mkdir(mode=0o700)
            except FileExistsError:
                raise ValueError(f"Checkpoint already exists: {self.path}. Use --resume-from or choose a new output.") from None
        # An OS-released SQLite lock covers the full run, while the separate
        # chunk database can commit each finished chunk independently.
        self.lock = sqlite3.connect((self.path / "owner.sqlite").as_uri() + "?mode=" +
                                    ("rw" if self.resuming else "rwc"), uri=True, timeout=0)
        try:
            self.lock.execute("BEGIN EXCLUSIVE")
        except sqlite3.OperationalError as exc:
            raise ValueError(f"Cannot lock checkpoint {self.path}: {exc}") from exc
        self.db = sqlite3.connect((self.path / "chunks.sqlite").as_uri() + "?mode=" +
                                  ("rw" if self.resuming else "rwc"), uri=True)
        if not self.resuming:
            with self.db:
                self.db.execute("CREATE TABLE signatures (target TEXT PRIMARY KEY, config TEXT NOT NULL)")
                self.db.execute("""CREATE TABLE chunks (target TEXT, id TEXT, summary TEXT,
                    checksum TEXT, complete_data BLOB, incomplete_data BLOB, manifest_data BLOB,
                    PRIMARY KEY (target, id))""")
        print(f"methyldiver: checkpoint {self.path}", file=sys.stderr, flush=True)

    def bind(self, target, signature):
        self._open()
        config = _json({"format": 1, "targets": self.targets, **signature})
        saved = self.db.execute("SELECT config FROM signatures WHERE target=?", (target,)).fetchone()
        if saved:
            if saved[0] != config:
                raise ValueError("Incompatible checkpoint: inputs, analysis options, software or chunk plan changed.")
        else:
            # A second target may be new after a failure in the first target.
            other = self.db.execute("SELECT config FROM signatures LIMIT 1").fetchone()
            if other and json.loads(other[0])["targets"] != self.targets:
                raise ValueError("Incompatible checkpoint: modification selection changed.")
            with self.db:
                self.db.execute("INSERT INTO signatures VALUES (?,?)", (target, config))

    def save(self, target, stage, summary, limits):
        name = summary["id"]
        complete = stage / "complete_patterns" / (name + ".parquet")
        incomplete = stage / "incomplete_patterns" / (name + ".parquet")
        meta = _json({"summary": summary, "limits": limits})
        with self.db:
            self.db.execute("INSERT OR REPLACE INTO chunks VALUES (?,?,?,?,?,?,?)",
                (target, name, meta, hashlib.sha256(meta.encode()).hexdigest(), complete.read_bytes(),
                 incomplete.read_bytes() if summary["outputs"][f"incomplete_patterns/{name}.parquet"]["retained"] else None,
                 (stage / summary["manifest"]).read_bytes() if summary["manifest"] else None))

    def restore(self, target, job, stage, keep_patterns, chunk_manifests, limits):
        name = job["id"]
        row = self.db.execute("SELECT summary,checksum,complete_data,incomplete_data,manifest_data "
                              "FROM chunks WHERE target=? AND id=?", (target, name)).fetchone()
        if row is None:
            return None
        meta, checksum, complete, incomplete, manifest = row
        if hashlib.sha256(meta.encode()).hexdigest() != checksum:
            return None
        record = json.loads(meta)
        summary = record["summary"]
        # A stricter guard must be applied by the estimator for this chunk.
        if any(limits[k] < v for k, v in record["limits"].items()):
            return None
        if summary["id"] != name or summary["n_units"] != job["n_units"]:
            return None
        relative = f"complete_patterns/{name}.parquet"
        files = [(relative, complete, summary["outputs"][relative]["sha256"])]
        if keep_patterns:
            relative = f"incomplete_patterns/{name}.parquet"
            files.append((relative, incomplete, summary["outputs"][relative]["sha256"]))
        if chunk_manifests:
            files.append((f"chunk_manifests/{name}.json.gz", manifest, summary["manifest_sha256"]))
        if any(data is None or hashlib.sha256(data).hexdigest() != digest for _, data, digest in files):
            return None
        for relative, data, _ in files:
            (stage / relative).write_bytes(data)
        summary["reused_with_original_limits"] = record["limits"]
        return summary

    def __exit__(self, exc_type, exc, traceback):
        # Only delete our own two files after final output publication. Leave
        # user files and incompatible checkpoints intact.
        opened = self.db is not None
        try:
            if self.db is not None:
                self.db.close()
            if opened and exc_type is None:
                (self.path / "chunks.sqlite").unlink()
        finally:
            if self.lock is not None:
                self.lock.close()
        if opened and exc_type is None:
            (self.path / "owner.sqlite").unlink(missing_ok=True)
            if not any(self.path.iterdir()):
                self.path.rmdir()
        elif opened:
            print(f"methyldiver: resume with --resume-from {self.path}", file=sys.stderr)
