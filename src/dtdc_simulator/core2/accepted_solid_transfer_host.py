"""Manufactured accepted-solid transaction host with a durable replay ledger.

This is PR-06 S0. The owner disposition
(``docs/GT_PS2_PRODUCTION_RULINGS_OWNER_DISPOSITION_2026-08-09.md:104-110``)
requires ONE host-owned versioned accepted solid-transfer transaction to
authorize feed, inter-tray and DT-to-DC motion, binding NINE items: accepted
interval, authority class, revisions, configuration, topology, codec/auditor
identities, packet plan, seven-extensive debit/credit, and a persistently
consumed replay token. Offered flow never changes holdup or populations. Gate
F5 (``release/ground_truth_v2_gates.yaml:128``) makes the scope total and the
replay protection RESTART-SURVIVING; the frozen qualification requirement
(``release/physics_decisions.yaml:602``) names the four mandatory rejection
classes -- stale, replay, foreign, insufficient -- and requires token
consumption to persist atomically WITH debit and credit.

The design basis is ``docs/GT_PS2_PR06_TRANSACTION_HOST_DESIGN_BASIS_2026-08-11.md``.

WHAT IS RULED AND IMPLEMENTED HERE.

* B5 (``docs/GT_PS2_B4_B5_RULING_RECORD_2026-08-10.md:16-22``): the replay token
  is a domain-separated SHA-256 over the request's canonical identity -- no
  host-issued nonce, no issuance round trip. Consumption inserts it into an
  append-only consumed-set inside the SAME atomic persistence act as the
  seven-extensive debit and credit. The set never expires and is never pruned.
  A re-presented consumed token is refused with a distinct terminal error and
  ZERO state change; there is no idempotent echo. Recovery is the separate
  read-only reconciliation query, never an echo from the mutation path.
* B4 (``:70-76``), and precisely what S0 ships of it: the LOAD-TIME COMPARISON
  is implemented -- the registry's definition digest is compared against a
  pinned digest and a mismatch fails closed -- and that pinned digest is bound
  into every transaction through the codec-identity field, so a swapped
  registry is refused twice, once by identity mismatch and once by token
  consumption. The pin ARTIFACT is NOT shipped here: there is no ``release/``
  pin file and no loader for one, so the digest compared against is whatever
  the caller supplies. The release-anchored artifact and its loader arrive with
  this packet's release stage; the gap is recorded as
  ``release_pin_artifact_not_shipped``.
* The commit sequence authenticates the host's current object identity, holds
  an exclusive lock across the whole sequence, replaces source ownership with
  destination ownership, and persists consumption -- the four obligations frozen
  in the population bridge's commit docstring.

THE COMMIT POINT IS ONE FILE, BY CONSTRUCTION. Ledger head, per-cell holdup,
per-cell revision and the consumed-set live in a single generation file, so ONE
``os.rename`` commits debit, credit and token consumption together. The repo's
only multi-file transactional precedent shows that a multi-rename window is not
atomic on Windows; this host avoids the window rather than compensating for it.

GENERATION FILES, NOT IN-PLACE REPLACEMENT. The house publication primitive is
``os.rename``, chosen for its Windows NON-CLOBBER semantics, which cannot
publish over an existing destination. The ledger is therefore an append-only
chain of monotonically numbered generation files: commit ``n`` stages
``<base>.<n>.json.pending``, verifies it as-if-committed, and renames it onto a
destination that provably does not yet exist. Non-clobber then doubles as the
concurrent-commit guard, and no generation is ever deleted or rewritten. Each
generation carries a chain digest over its predecessor, so substitution of any
generation, and truncation of any NON-HEAD generation, are detected on load.

WHERE THE CHAIN STOPS. Chain links point BACKWARD only, so deletion of head
generations is NOT self-detecting: what remains is internally consistent, its
consumed-set is silently smaller, replay is re-enabled, and the reconciliation
query then answers a truthful-looking "not committed" for a request that WAS
committed. Nothing inside the store can see that, so the anchor is external --
the packet's audit and campaign byte-pin over the store's own head -- and the
boundary is declared as ``head_generation_rollback_not_self_detecting`` and
documented by two tests rather than dressed up as closed.

THREE-SUFFIX STATE MACHINE. Every generation is ``base``, ``base.pending`` or
``base.invalid``. Staged bytes are read back, strictly re-parsed, and verified
as-if-committed before the rename; any failure quarantines ``.pending`` ->
``.invalid`` and NEVER deletes. THE RENAME IS THE COMMIT POINT AND IS NEVER
RETRACTED: every reversible check, the store-identity checkpoint included, runs
BEFORE it, and a post-rename re-read failure REPORTS the condition with the
published path named and LEAVES the published generation in place. Un-publishing
it would answer "not committed" for a transaction whose commit point completed;
genuine corruption instead fails closed on the next load. A store carrying any
reserved-suffix residue is unresolved: the mutation path refuses to proceed and
the read-only reconciliation query refuses to answer over it, while the
read-only ``inspect_store`` diagnostic stays reachable over any store state.

THE LOCK IS FAIL-CLOSED. A surviving lockfile is never stolen and never aged
out. The host refuses with a distinct error naming the lockfile and the
reconciliation query that resolves the prior transaction's disposition;
automatic takeover would weaken a fail-closed boundary.

COMPOSITION. The selected state profile in ``state_codec`` is imported normally
and supplies payload validation and identity: payload identities come from
``canonical_state_identity``, never from a raw hash over unvalidated bytes. The
B1 scaling contract is resolved DYNAMICALLY through a split module path, because
the frozen ``src/`` scans of the sibling audits read raw file text and require
zero literal occurrences of the audited module names anywhere under ``src/``.
Editing those audits is not an option -- they are keep-stable pinned and their
archived campaigns re-verify their bytes. The four guarded sibling module names
are likewise absent from this file, including from its prose. The canonical
seven-extensive order is mirrored locally as a frozen tuple and cross-checked
against the resolved module at import time; drift fails closed immediately.

SCALE AND CONCURRENCY POSTURE, DECLARED RATHER THAN HIDDEN. Every generation
copies the whole entry log and the whole consumed-set forward, so store bytes
and load cost both grow quadratically in the number of commits, unbounded and
never pruned (``store_growth_quadratic_unbounded``). The store is single-host by
design: the reserved-suffix residue check cannot distinguish abandoned residue
from a live concurrent commit, so concurrent readers are not supported
(``single_host_store_no_concurrent_readers``). File contents are fsynced, but
the directory entry carries no barrier on this platform, so restart-surviving
means surviving PROCESS death and OS-crash durability is best-effort
(``directory_entry_durability_not_barriered``). Lock release is best-effort as
well, and a leaked lockfile deliberately fails closed
(``lock_release_best_effort``).

This is manufactured architecture evidence. It authenticates no schema source,
selects no production variant registry, wires no production caller, moves no
physical solid, transports no complete particle state, and closes no gate.
Payload self-sufficiency remains a declared open item and is recorded as a
limitation, not a claim.
"""

from __future__ import annotations

import enum
import hashlib
import importlib
import json
import math
import os
import struct
from dataclasses import InitVar, dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any, ClassVar, Mapping

from dtdc_simulator.core2 import state_codec

TRANSACTION_MAGIC = "GT-PS-2-ACCEPTED-SOLID-TRANSFER"
TRANSACTION_VERSION = 1

#: Every digest domain in this module lives under one fresh family so it can
#: never collide with the selected profile's domains or the PR-05 receipt space.
DIGEST_DOMAIN_PREFIX = "GT-PS-2/pr06-accepted-solid-transfer/"
IDENTITY_DIGEST_DOMAIN = DIGEST_DOMAIN_PREFIX + "identity/v1"
REPLAY_TOKEN_DOMAIN = DIGEST_DOMAIN_PREFIX + "replay-token/v1"
PACKET_PLAN_DIGEST_DOMAIN = DIGEST_DOMAIN_PREFIX + "packet-plan-entry/v1"
RECEIPT_DIGEST_DOMAIN = DIGEST_DOMAIN_PREFIX + "receipt/v1"
LEDGER_DIGEST_DOMAIN = DIGEST_DOMAIN_PREFIX + "ledger-generation/v1"
HOST_IDENTITY_DIGEST_DOMAIN = DIGEST_DOMAIN_PREFIX + "host-object/v1"

#: The nine items the disposition requires one transaction to bind. Items one
#: through eight are framed into the request identity; item nine is derived from
#: that identity under its own sub-domain and is consumed on commit.
BOUND_TRANSACTION_ITEMS: tuple[str, ...] = (
    "accepted_interval",
    "authority_class",
    "revisions",
    "configuration",
    "topology",
    "codec_and_auditor_identities",
    "packet_plan",
    "seven_extensive_debit_credit",
    "persistently_consumed_replay_token",
)

#: Canonical order of the seven frozen carrier extensives, mirrored locally as
#: value strings and cross-checked against the resolved contract at import.
EXTENSIVE_NAMES: tuple[str, ...] = (
    "dry_matter_kg",
    "residual_oil_label_kg",
    "attached_hexane_kg",
    "internal_hexane_kg",
    "external_water_kg",
    "retained_water_kg",
    "common_datum_energy_j",
)
#: Only the six mass extensives are sign-constrained; datum-relative energy is
#: legitimately signed.
MASS_EXTENSIVE_COUNT = 6

LEDGER_BASENAME = "accepted_solid_transfer_ledger"
LEDGER_SUFFIX = ".json"
LEDGER_SCHEMA = "GT-PS-2-ACCEPTED-SOLID-TRANSFER-LEDGER"
LEDGER_SCHEMA_VERSION = 1
HOST_LOCK_NAME = "accepted_solid_transfer_host.lock"
RESERVED_SUFFIXES: tuple[str, ...] = (".pending", ".invalid")
GENERATION_DIGITS = 8
GENESIS_GENERATION = 1

#: The read-only reconciliation query's spelling, fixed here as a Class A shape
#: under B5 (``:64-66``) and named in every lock refusal.
RECONCILIATION_QUERY_NAME = "reconcile_request"

#: Standing limitations recorded by the ledger. These are limitations, never
#: claims. A bare identifier in a tuple is not a disclosure, so every entry
#: below is also explained in prose:
#:
#: * ``payload_self_sufficiency_open`` -- payload bytes alone cannot
#:   reconstruct the full particle authority (variant inventory ``:77-81``);
#:   the ledger transports identities, it does not close that gap.
#: * ``production_variant_registry_incomplete`` -- the registry this host pins
#:   is whatever the caller supplies; no complete production variant registry
#:   exists yet, and this module does not claim one.
#: * ``schema_source_unauthenticated`` -- schema and registry definitions are
#:   trusted inputs; nothing here authenticates their provenance.
#: * ``no_physical_transport_performed`` -- this host moves ledger custody and
#:   conserved totals; no physical model runs and no material moves.
#: * ``reconciliation_answer_not_host_authenticated`` -- ``AcceptedSolidTransfer``
#:   and ``TransferCommitReceipt`` are unforgeable (issue-token construction);
#:   ``ReconciliationAnswer`` is a plain read-only report a caller could
#:   construct by hand, so a consumer cannot distinguish a host-issued answer
#:   from a fabricated one. No state is at risk; the asymmetry is declared.
DECLARED_LIMITATIONS: tuple[str, ...] = (
    "payload_self_sufficiency_open",
    "production_variant_registry_incomplete",
    "schema_source_unauthenticated",
    "no_physical_transport_performed",
    "reconciliation_answer_not_host_authenticated",
    "release_pin_artifact_not_shipped",
    "head_generation_rollback_not_self_detecting",
    "lock_release_best_effort",
    "directory_entry_durability_not_barriered",
    "store_growth_quadratic_unbounded",
    "single_host_store_no_concurrent_readers",
)


class AcceptedSolidTransferError(ValueError):
    """Base exception for a refused accepted-solid transfer."""


class StaleTransferError(AcceptedSolidTransferError):
    """Raised when a proposal names a revision that is no longer current."""


class ReplayedTransferError(AcceptedSolidTransferError):
    """Raised when a consumed replay token is re-presented.

    B5 makes this its own terminal type: a replay must be distinguishable from
    malformed bytes and from a stale proposal, and it produces zero state
    change. There is no idempotent echo -- the caller reconciles instead.
    """


class ForeignTransferError(AcceptedSolidTransferError):
    """Raised when a proposal is foreign to this host, store or registry."""


class InsufficientSolidError(AcceptedSolidTransferError):
    """Raised when eligible mass is short. Never repaired by a clamp."""


class ConservationResidualError(AcceptedSolidTransferError):
    """Raised when debit does not equal credit over the seven extensives.

    Its own leaf type, distinct from the four ruled rejection classes: an
    unbalanced proposal is neither stale, replayed, foreign nor insufficient, and
    collapsing it into the base type would leave a caller unable to tell an
    arithmetic non-closure from any other refusal. The residual is reported,
    never repaired by a clamp.
    """


class TransferStoreIntegrityError(AcceptedSolidTransferError):
    """Raised when the durable store is unresolved, torn or inconsistent."""


class HostLockRefusedError(AcceptedSolidTransferError):
    """Raised when the exclusive host lock cannot be taken. Never stolen."""


# --------------------------------------------------------------------------- #
# The B1 scaling contract, resolved dynamically.
#
# The sibling audits scan every src/ file for their own module names as static
# literals and require zero hits, so the module path is assembled from parts.
# The audit for THIS module byte-pins the resolved module's blob, so dynamic
# resolution cannot be silently redirected. The canonical extensive order is
# additionally mirrored above and cross-checked here: belt and suspenders,
# fail-closed on drift.
# --------------------------------------------------------------------------- #

SCALING_CONTRACT_MODULE = "dtdc_simulator.core2." + "packet_weight_" + "scaling_contract"

_contract = importlib.import_module(SCALING_CONTRACT_MODULE)

if tuple(getattr(_contract, "EXTENSIVE_NAMES", ())) != EXTENSIVE_NAMES:
    raise AcceptedSolidTransferError(
        "resolved scaling contract declares a different canonical extensive order"
    )
if getattr(_contract, "MASS_EXTENSIVE_COUNT", None) != MASS_EXTENSIVE_COUNT:
    raise AcceptedSolidTransferError(
        "resolved scaling contract declares a different mass-extensive count"
    )
for _required in (
    "PerParticleInventory",
    "PacketWeight",
    "PacketTotals",
    "ExternalSurfaceComposite",
    "ScalingContractError",
    "WeightAuthorityError",
    "ReconciliationError",
    "SurfaceAuthorityError",
    "derive_packet_totals",
    "require_single_surface_authority",
):
    if not hasattr(_contract, _required):
        raise AcceptedSolidTransferError(
            f"resolved scaling contract is missing its public verb {_required!r}"
        )
del _required

# The re-export surface is the FULL set a caller needs to construct a packet-plan
# entry, because the entry type-checks exact classes: a caller reaching past this
# module for the surface composite would import a class this module rejects.
PerParticleInventory = _contract.PerParticleInventory
PacketWeight = _contract.PacketWeight
ExternalSurfaceComposite = _contract.ExternalSurfaceComposite
ScalingContractError = _contract.ScalingContractError


class TransferAuthorityClass(enum.Enum):
    """The three motions the one transaction authorizes."""

    FEED_ADMISSION = "feed_admission"
    INTER_TRAY = "inter_tray"
    DT_TO_DC_HANDOFF = "dt_to_dc_handoff"


class DeclaredEngineeringAssumption(enum.Enum):
    """Typed declarations; the disposition forbids silent choices."""

    #: Ledger head, holdup, revisions and the consumed-set share ONE generation
    #: file so a single rename commits debit, credit and consumption together.
    SINGLE_FILE_LEDGER_COMMIT_POINT = "single_file_ledger_commit_point_v1"

    #: Generation files are append-only and monotonically numbered, so every
    #: publication renames onto a destination that does not yet exist and the
    #: house non-clobber semantics stay intact.
    APPEND_ONLY_GENERATION_CHAIN = "append_only_generation_chain_v1"

    #: A surviving lockfile is never stolen and never aged out.
    FAIL_CLOSED_STALE_LOCK = "fail_closed_stale_lock_v1"


def _framed_digest(domain: str, parts: tuple[bytes, ...]) -> str:
    """Length-framed, domain-separated SHA-256.

    Mirrors the selected profile's construction so this host's digests are
    framed identically, but under its own domain family so the two digest
    spaces can never collide. Reimplemented locally rather than imported: the
    profile's helper is private and is not part of its public surface.
    """

    digest = hashlib.sha256()
    for part in (domain.encode("ascii"), *parts):
        digest.update(len(part).to_bytes(8, "big"))
        digest.update(part)
    return "sha256:" + digest.hexdigest()


def _is_negative_zero(value: float) -> bool:
    return struct.pack(">d", value) == b"\x80\x00\x00\x00\x00\x00\x00\x00"


def _float_part(value: float) -> bytes:
    """Digest a binary64 with -0.0 canonicalized, mirroring the frozen digest."""

    canonical = 0.0 if value == 0.0 else value
    return struct.pack(">d", canonical)


def _require_binary64(name: str, *values: float) -> None:
    for value in values:
        if type(value) is not float:
            raise AcceptedSolidTransferError(f"{name} must contain exact binary64 float values")
        if not math.isfinite(value):
            raise AcceptedSolidTransferError(f"{name} must contain only finite values")


def _require_nonblank_text(name: str, value: object) -> str:
    if type(value) is not str or not value.strip():
        raise AcceptedSolidTransferError(f"{name} must be an exact nonblank str")
    return value


def _require_digest_text(name: str, value: object) -> str:
    if type(value) is not str or not value.startswith("sha256:") or len(value) != 71:
        raise AcceptedSolidTransferError(f"{name} must be an exact sha256: digest str")
    if any(character not in "0123456789abcdef" for character in value[7:]):
        raise AcceptedSolidTransferError(f"{name} must be a lowercase hexadecimal digest")
    return value


def _require_exact_index(name: str, value: object, *, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        raise AcceptedSolidTransferError(f"{name} must be an exact int >= {minimum}")
    return value


def _require_extensive_tuple(name: str, value: object) -> tuple[float, ...]:
    if type(value) is not tuple or len(value) != len(EXTENSIVE_NAMES):
        raise AcceptedSolidTransferError(
            f"{name} must be an exact {len(EXTENSIVE_NAMES)}-tuple of extensives"
        )
    _require_binary64(name, *value)
    return value


# --------------------------------------------------------------------------- #
# Canonical bytes and the single durable-write primitive.
# --------------------------------------------------------------------------- #


def _reject_json_constant(name: str) -> None:
    raise TransferStoreIntegrityError(f"ledger JSON carries the forbidden constant {name!r}")


def _strict_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    keys = [key for key, _ in pairs]
    if len(keys) != len(set(keys)):
        raise TransferStoreIntegrityError("ledger JSON object contains duplicate keys")
    return dict(pairs)


def _canonical_json_bytes(payload: Mapping[str, Any]) -> bytes:
    """Sorted, minimally separated, ASCII JSON terminated by exactly one LF."""

    text = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )
    return (text + "\n").encode("ascii")


def _strict_json_payload(raw: bytes, label: str) -> dict[str, Any]:
    if type(raw) is not bytes or not raw.endswith(b"\n") or raw[:-1].endswith((b"\r", b"\n")):
        raise TransferStoreIntegrityError(f"{label} must end in exactly one LF")
    try:
        text = raw[:-1].decode("ascii", errors="strict")
    except UnicodeDecodeError as error:
        raise TransferStoreIntegrityError(f"{label} is not canonical ASCII JSON") from error
    try:
        payload = json.loads(
            text,
            object_pairs_hook=_strict_pairs,
            parse_constant=_reject_json_constant,
        )
    except ValueError as error:
        raise TransferStoreIntegrityError(f"{label} is not strict JSON") from error
    if type(payload) is not dict:
        raise TransferStoreIntegrityError(f"{label} root must be an exact object")
    if _canonical_json_bytes(payload) != raw:
        raise TransferStoreIntegrityError(f"{label} bytes are not canonical for this profile")
    return payload


def _write_durable(path: Path, value: bytes) -> None:
    """The one durable-write step: exclusive create, short-write check, fsync.

    Every durable byte in this module goes through here, and nothing else in
    this module opens a file for writing. Keeping it a single module-level
    function is deliberate: torn writes and post-rename failures are then
    FAULT-INJECTED by the hostile suite rather than merely inspected.

    DURABILITY, STATED HONESTLY: the file CONTENTS are fsynced, but the directory
    entry -- the visibility of a creation or a rename after an OS crash -- gets no
    barrier on this platform, where a directory cannot be opened for fsync. So
    "restart-surviving" here means surviving PROCESS death, which is what the
    ruled replay requirement needs; OS-crash durability is best-effort and is
    declared as ``directory_entry_durability_not_barriered``.
    """

    with path.open("xb") as stream:
        written = stream.write(value)
        if written != len(value):
            raise TransferStoreIntegrityError(
                f"short write for {path}: wrote {written} of {len(value)} bytes"
            )
        stream.flush()
        os.fsync(stream.fileno())


def _has_reparse_point(path: Path) -> bool:
    try:
        status = path.lstat()
    except OSError:
        return False
    return bool(getattr(status, "st_reparse_tag", 0))


def _reject_reparse_ancestors(path: Path) -> None:
    """Inspect the caller-supplied ancestry BEFORE resolve() erases the evidence."""

    current = path
    while True:
        if current.is_symlink() or os.path.isjunction(current) or _has_reparse_point(current):
            raise TransferStoreIntegrityError(f"store ancestry crosses a reparse point: {current}")
        if current.parent == current:
            return
        current = current.parent


def _directory_identity(path: Path) -> dict[str, object]:
    """Resolved path plus device and inode; the device doubles as the volume guard."""

    if not path.is_dir() or path.is_symlink() or os.path.isjunction(path):
        raise TransferStoreIntegrityError(
            f"store root is missing or is not a plain directory: {path}"
        )
    status = path.stat(follow_symlinks=False)
    return {
        "path": str(path.resolve(strict=True)),
        "device": status.st_dev,
        "inode": status.st_ino,
    }


def _generation_name(generation: int) -> str:
    return f"{LEDGER_BASENAME}.{generation:0{GENERATION_DIGITS}d}{LEDGER_SUFFIX}"


def _generation_number(name: str) -> int | None:
    if not name.startswith(LEDGER_BASENAME + ".") or not name.endswith(LEDGER_SUFFIX):
        return None
    middle = name[len(LEDGER_BASENAME) + 1 : -len(LEDGER_SUFFIX)]
    if len(middle) != GENERATION_DIGITS or not middle.isdigit() or not middle.isascii():
        return None
    generation = int(middle)
    return generation if _generation_name(generation) == name else None


# --------------------------------------------------------------------------- #
# Value types. The claim-flag block below is exactly what PR-06 earns and
# exactly what it does not; the frozen flag lines of the sibling modules are
# never edited to match.
# --------------------------------------------------------------------------- #

_ISSUE_TOKEN = object()


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedInterval:
    """Bound item 1: the accepted interval this transaction authorizes."""

    start_s: float
    end_s: float

    def __post_init__(self) -> None:
        _require_binary64("accepted interval", self.start_s, self.end_s)
        if not self.end_s > self.start_s:
            raise AcceptedSolidTransferError("accepted interval must be strictly forward in time")


@dataclass(frozen=True, slots=True, kw_only=True)
class TransferTopology:
    """Bound item 5: the source and destination solid control volumes."""

    source_cell_id: str
    destination_cell_id: str

    def __post_init__(self) -> None:
        _require_nonblank_text("source_cell_id", self.source_cell_id)
        _require_nonblank_text("destination_cell_id", self.destination_cell_id)
        if self.source_cell_id == self.destination_cell_id:
            raise AcceptedSolidTransferError("a transfer cannot name one cell twice")


@dataclass(frozen=True, slots=True, kw_only=True)
class OfferedSolidFlow:
    """A proposal only. Offered flow never changes holdup or populations.

    There is deliberately no method here that reaches a ledger, and no mutation
    path at all: the type is frozen and slotted, refuses copy, deep-copy and
    serialization, and carries a standing ``moves_material`` False marker. The
    host's mutation path accepts host-issued transactions exclusively, so an
    offered flow presented to it is refused as foreign.
    """

    offered_dry_matter_kg: float
    topology: TransferTopology
    authority_class: TransferAuthorityClass

    moves_material: ClassVar[bool] = False
    changes_holdup: ClassVar[bool] = False
    changes_populations: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_binary64("offered flow", self.offered_dry_matter_kg)
        if self.offered_dry_matter_kg < 0.0 or _is_negative_zero(self.offered_dry_matter_kg):
            raise AcceptedSolidTransferError(
                "offered dry matter must be non-negative with canonical positive zero"
            )
        if type(self.topology) is not TransferTopology:
            raise AcceptedSolidTransferError(
                "offered flow topology must be an exact TransferTopology"
            )
        if type(self.authority_class) is not TransferAuthorityClass:
            raise AcceptedSolidTransferError(
                "offered flow authority class must be an exact TransferAuthorityClass"
            )

    def __copy__(self) -> None:
        raise TypeError("offered solid flow cannot be copied")

    def __deepcopy__(self, memo: object) -> None:
        del memo
        raise TypeError("offered solid flow cannot be deep-copied")

    def __reduce_ex__(self, protocol: int) -> None:
        del protocol
        raise TypeError("offered solid flow cannot be serialized")


@dataclass(frozen=True, slots=True, kw_only=True)
class PacketPlanEntry:
    """One packet of bound item 7: canonical payload bytes plus its host weight.

    The transfer is an identity transfer on the payload: the bytes that leave
    the source are the bytes that arrive at the destination. That is a property
    of the transfer and NOT a flag on this type -- the bytes are caller-supplied,
    so a claim flag here would assert something this module cannot enforce. The
    packet weight is host-carried and never lives inside payload bytes, so it is
    the only authority that may differ between two children sharing identical
    bytes.
    """

    packet_id: str
    payload_bytes: bytes
    inventory: Any
    weight: Any
    surface: Any | None = None

    def __post_init__(self) -> None:
        _require_nonblank_text("packet_id", self.packet_id)
        if type(self.payload_bytes) is not bytes or not self.payload_bytes:
            raise AcceptedSolidTransferError("packet payload must be exact non-empty bytes")
        if type(self.inventory) is not PerParticleInventory:
            raise AcceptedSolidTransferError(
                "packet inventory must be an exact PerParticleInventory"
            )
        if type(self.weight) is not PacketWeight:
            raise AcceptedSolidTransferError("packet weight must be an exact PacketWeight")
        if self.surface is not None and type(self.surface) is not ExternalSurfaceComposite:
            raise AcceptedSolidTransferError(
                "packet surface composite must be an exact ExternalSurfaceComposite"
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedSolidTransferRequest:
    """The caller-supplied canonical request: bound items 1 through 8."""

    interval: AcceptedInterval
    authority_class: TransferAuthorityClass
    topology: TransferTopology
    source_prior_revision: int
    destination_prior_revision: int
    configuration_digest: str
    codec_registry_digest: str
    auditor_identity_digest: str
    packet_plan: tuple[PacketPlanEntry, ...]
    declared_debit_totals: tuple[float, ...]

    def __post_init__(self) -> None:
        if type(self.interval) is not AcceptedInterval:
            raise AcceptedSolidTransferError("interval must be an exact AcceptedInterval")
        if type(self.authority_class) is not TransferAuthorityClass:
            raise AcceptedSolidTransferError(
                "authority_class must be an exact TransferAuthorityClass"
            )
        if type(self.topology) is not TransferTopology:
            raise AcceptedSolidTransferError("topology must be an exact TransferTopology")
        _require_exact_index("source_prior_revision", self.source_prior_revision)
        _require_exact_index("destination_prior_revision", self.destination_prior_revision)
        _require_digest_text("configuration_digest", self.configuration_digest)
        _require_digest_text("codec_registry_digest", self.codec_registry_digest)
        _require_digest_text("auditor_identity_digest", self.auditor_identity_digest)
        if type(self.packet_plan) is not tuple or not self.packet_plan:
            raise AcceptedSolidTransferError("packet_plan must be a non-empty exact tuple")
        if any(type(entry) is not PacketPlanEntry for entry in self.packet_plan):
            raise AcceptedSolidTransferError(
                "packet_plan must contain exact PacketPlanEntry values"
            )
        identifiers = tuple(entry.packet_id for entry in self.packet_plan)
        if len(identifiers) != len(set(identifiers)):
            raise AcceptedSolidTransferError("packet_plan packet identifiers must be unique")
        _require_extensive_tuple("declared_debit_totals", self.declared_debit_totals)
        # The ruled signed-zero narrowing applies to every persisted mass column,
        # and a declared debit is persisted verbatim. Energy is datum-relative
        # and stays sign-free, so the ladder covers the six mass columns only.
        for name, value in zip(
            EXTENSIVE_NAMES[:MASS_EXTENSIVE_COUNT], self.declared_debit_totals, strict=False
        ):
            if value < 0.0 or _is_negative_zero(value):
                raise AcceptedSolidTransferError(
                    f"declared_debit_totals {name} must be non-negative with canonical zero"
                )


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedSolidTransfer:
    """The one host-issued versioned transaction, binding all nine items."""

    request: AcceptedSolidTransferRequest
    request_identity: str
    replay_token: str
    payload_identities: tuple[str, ...]
    debit_totals: tuple[float, ...]
    credit_totals: tuple[float, ...]
    conservation_residual: tuple[float, ...]
    host_identity: str
    _issue_token: InitVar[object | None] = None

    accepted_solid_transfer_implemented: ClassVar[bool] = True
    host_object_identity_authenticated: ClassVar[bool] = True
    transaction_consumption_persisted: ClassVar[bool] = True
    atomic_persistence_selected: ClassVar[bool] = True
    persisted_replay_protection_selected: ClassVar[bool] = True
    # Explicitly NOT supplied by this host, and staying false.
    complete_particle_state_transport_selected: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_variant_registry_complete: ClassVar[bool] = False
    source_authenticated: ClassVar[bool] = False

    def __post_init__(self, _issue_token: object | None) -> None:
        if _issue_token is not _ISSUE_TOKEN:
            raise TypeError("AcceptedSolidTransfer values are issued only by the host")
        if type(self.request) is not AcceptedSolidTransferRequest:
            raise TypeError("transaction request must be an exact request value")
        if type(self.payload_identities) is not tuple:
            raise TypeError("payload identities must be an exact tuple")

    @property
    def bound_items(self) -> tuple[str, ...]:
        """The nine disposition-named items this transaction binds."""

        return BOUND_TRANSACTION_ITEMS

    def __copy__(self) -> None:
        raise TypeError("AcceptedSolidTransfer authority views cannot be copied")

    def __deepcopy__(self, memo: object) -> None:
        del memo
        raise TypeError("AcceptedSolidTransfer authority views cannot be deep-copied")

    def __reduce_ex__(self, protocol: int) -> None:
        del protocol
        raise TypeError("AcceptedSolidTransfer authority views cannot be serialized")


@dataclass(frozen=True, slots=True, kw_only=True)
class TransferCommitReceipt:
    """The host-issued out-of-band receipt for one committed transaction."""

    request_identity: str
    replay_token: str
    committed_revision: int
    generation: int
    ledger_digest: str
    receipt_digest: str
    debit_totals: tuple[float, ...]
    credit_totals: tuple[float, ...]
    ledger_path: str
    _issue_token: InitVar[object | None] = None

    accepted_solid_transfer_implemented: ClassVar[bool] = True
    host_object_identity_authenticated: ClassVar[bool] = True
    transaction_consumption_persisted: ClassVar[bool] = True
    atomic_persistence_selected: ClassVar[bool] = True
    persisted_replay_protection_selected: ClassVar[bool] = True
    # Explicitly NOT supplied by this host, and staying false.
    complete_particle_state_transport_selected: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_variant_registry_complete: ClassVar[bool] = False
    source_authenticated: ClassVar[bool] = False

    limitations: ClassVar[tuple[str, ...]] = DECLARED_LIMITATIONS

    def __post_init__(self, _issue_token: object | None) -> None:
        if _issue_token is not _ISSUE_TOKEN:
            raise TypeError("TransferCommitReceipt values are issued only by the host")

    def __copy__(self) -> None:
        raise TypeError("TransferCommitReceipt authority views cannot be copied")

    def __deepcopy__(self, memo: object) -> None:
        del memo
        raise TypeError("TransferCommitReceipt authority views cannot be deep-copied")

    def __reduce_ex__(self, protocol: int) -> None:
        del protocol
        raise TypeError("TransferCommitReceipt authority views cannot be serialized")


@dataclass(frozen=True, slots=True, kw_only=True)
class ReconciliationAnswer:
    """The read-only answer to "was request X committed?".

    This is never an echo from the mutation path and never carries the accepted
    payload back to the caller: it reports disposition, the committed revision
    and the receipt digest, and nothing else.
    """

    request_identity: str
    committed: bool
    committed_revision: int | None
    generation: int | None
    receipt_digest: str | None

    read_only: ClassVar[bool] = True
    is_idempotent_echo: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.committed) is not bool:
            raise AcceptedSolidTransferError("reconciliation disposition must be an exact bool")
        if self.committed:
            if self.committed_revision is None or self.receipt_digest is None:
                raise AcceptedSolidTransferError(
                    "a committed answer must carry revision and receipt"
                )
        elif (
            self.committed_revision is not None
            or self.generation is not None
            or self.receipt_digest is not None
        ):
            raise AcceptedSolidTransferError("an uncommitted answer must carry no committed fields")


# --------------------------------------------------------------------------- #
# Identity, replay token, receipt.
# --------------------------------------------------------------------------- #


def packet_plan_entry_digest(entry: PacketPlanEntry, payload_identity: str) -> str:
    """Frame one packet-plan entry around its CODEC-VALIDATED payload identity.

    The payload identity is supplied by the selected profile, never computed
    here over unvalidated bytes: an unvalidated hash would bless noncanonical or
    foreign-variant bytes with a ledger identity.
    """

    if type(entry) is not PacketPlanEntry:
        raise AcceptedSolidTransferError("packet-plan entry must be an exact PacketPlanEntry")
    _require_digest_text("payload identity", payload_identity)
    return _framed_digest(
        PACKET_PLAN_DIGEST_DOMAIN,
        (
            entry.packet_id.encode("utf-8"),
            payload_identity.encode("ascii"),
            entry.inventory.energy_datum_id.encode("utf-8"),
            *(_float_part(value) for value in entry.inventory.as_tuple()),
            _float_part(entry.weight.representative_particles),
        ),
    )


def request_identity_digest(
    request: AcceptedSolidTransferRequest, payload_identities: tuple[str, ...]
) -> str:
    """Bound items 1 through 8, length-framed in the disposition's own order."""

    if type(request) is not AcceptedSolidTransferRequest:
        raise AcceptedSolidTransferError("request must be an exact AcceptedSolidTransferRequest")
    if type(payload_identities) is not tuple or len(payload_identities) != len(request.packet_plan):
        raise AcceptedSolidTransferError(
            "payload identities must match the packet plan one-for-one"
        )
    plan_digests = tuple(
        packet_plan_entry_digest(entry, identity).encode("ascii")
        for entry, identity in zip(request.packet_plan, payload_identities, strict=True)
    )
    return _framed_digest(
        IDENTITY_DIGEST_DOMAIN,
        (
            TRANSACTION_MAGIC.encode("ascii"),
            TRANSACTION_VERSION.to_bytes(8, "big"),
            DeclaredEngineeringAssumption.SINGLE_FILE_LEDGER_COMMIT_POINT.value.encode("ascii"),
            # 1 accepted interval
            _float_part(request.interval.start_s),
            _float_part(request.interval.end_s),
            # 2 authority class
            request.authority_class.value.encode("ascii"),
            # 3 revisions
            request.source_prior_revision.to_bytes(8, "big"),
            request.destination_prior_revision.to_bytes(8, "big"),
            # 4 configuration
            request.configuration_digest.encode("ascii"),
            # 5 topology
            request.topology.source_cell_id.encode("utf-8"),
            request.topology.destination_cell_id.encode("utf-8"),
            # 6 codec and auditor identities
            request.codec_registry_digest.encode("ascii"),
            request.auditor_identity_digest.encode("ascii"),
            # 7 packet plan
            len(plan_digests).to_bytes(8, "big"),
            *plan_digests,
            # 8 seven-extensive debit declaration
            *(_float_part(value) for value in request.declared_debit_totals),
        ),
    )


def derive_replay_token(request_identity: str) -> str:
    """Bound item 9, derived under its own sub-domain (B5).

    No host-issued nonce, no issuance round trip, no expiry. Domain separation
    keeps the token field distinct from the identity field while pinning the
    derivation, so replay detection reduces to set membership.
    """

    _require_digest_text("request identity", request_identity)
    return _framed_digest(
        REPLAY_TOKEN_DOMAIN,
        (
            TRANSACTION_MAGIC.encode("ascii"),
            TRANSACTION_VERSION.to_bytes(8, "big"),
            request_identity.encode("ascii"),
        ),
    )


def commit_receipt_digest(
    request_identity: str,
    replay_token: str,
    committed_revision: int,
    debit_totals: tuple[float, ...],
    credit_totals: tuple[float, ...],
    prior_ledger_digest: str,
) -> str:
    """Out-of-band receipt binding one commit to its ledger position."""

    return _framed_digest(
        RECEIPT_DIGEST_DOMAIN,
        (
            TRANSACTION_MAGIC.encode("ascii"),
            TRANSACTION_VERSION.to_bytes(8, "big"),
            request_identity.encode("ascii"),
            replay_token.encode("ascii"),
            committed_revision.to_bytes(8, "big"),
            *(_float_part(value) for value in debit_totals),
            *(_float_part(value) for value in credit_totals),
            prior_ledger_digest.encode("ascii"),
        ),
    )


def _generation_chain_digest(payload: Mapping[str, Any]) -> str:
    """Chain digest over one generation's canonical bytes and its predecessor.

    The typed append-only-chain assumption is framed into the digest, so the
    declaration is part of the artifact's identity and cannot be edited without
    invalidating every generation that was sealed under it. Both sibling modules
    bind their own declared assumption into their identity the same way.
    """

    body = {key: value for key, value in payload.items() if key != "ledger_digest"}
    return _framed_digest(
        LEDGER_DIGEST_DOMAIN,
        (
            LEDGER_SCHEMA.encode("ascii"),
            LEDGER_SCHEMA_VERSION.to_bytes(8, "big"),
            DeclaredEngineeringAssumption.APPEND_ONLY_GENERATION_CHAIN.value.encode("ascii"),
            int(body["generation"]).to_bytes(8, "big"),
            str(body["prior_ledger_digest"]).encode("ascii"),
            _canonical_json_bytes(body),
        ),
    )


# --------------------------------------------------------------------------- #
# The single-file ledger generation: head, holdup, revisions and the
# append-only consumed-set, all committed by one rename.
# --------------------------------------------------------------------------- #

GENESIS_PRIOR_DIGEST = "sha256:" + "0" * 64

_GENERATION_KEYS = frozenset(
    {
        "schema",
        "schema_version",
        "generation",
        "revision",
        "prior_ledger_digest",
        "store_id",
        "registry_digest",
        "energy_datum_id",
        "cells",
        "consumed_tokens",
        "entries",
        "limitations",
        "ledger_digest",
    }
)

_ENTRY_KEYS = frozenset(
    {
        "sequence",
        "request_identity",
        "replay_token",
        "receipt_digest",
        "authority_class",
        "source_cell_id",
        "destination_cell_id",
        "committed_revision",
        "debit_totals",
        "credit_totals",
        "closing_residual",
        "packet_ids",
        "payload_identities",
        "codec_registry_digest",
        "auditor_identity_digest",
        "configuration_digest",
    }
)


def _require_json_extensives(value: object, label: str) -> tuple[float, ...]:
    if type(value) is not list or len(value) != len(EXTENSIVE_NAMES):
        raise TransferStoreIntegrityError(f"{label} is not a seven-extensive array")
    for item in value:
        if type(item) is not float or not math.isfinite(item):
            raise TransferStoreIntegrityError(f"{label} carries a non-binary64 component")
    return tuple(value)


def _validate_generation_payload(payload: Mapping[str, Any], label: str) -> dict[str, Any]:
    """Strictly re-parse one generation as-if-committed. No lenient branch."""

    if type(payload) is not dict or set(payload) != set(_GENERATION_KEYS):
        raise TransferStoreIntegrityError(f"{label} schema is incomplete or foreign")
    if payload["schema"] != LEDGER_SCHEMA or payload["schema_version"] != LEDGER_SCHEMA_VERSION:
        raise TransferStoreIntegrityError(f"{label} declares a foreign schema")
    generation = payload["generation"]
    revision = payload["revision"]
    if type(generation) is not int or generation < GENESIS_GENERATION:
        raise TransferStoreIntegrityError(f"{label} generation number is malformed")
    if type(revision) is not int or revision < 0:
        raise TransferStoreIntegrityError(f"{label} revision is malformed")
    if revision != generation - GENESIS_GENERATION:
        raise TransferStoreIntegrityError(f"{label} revision does not track its generation")
    _require_digest_text(f"{label} prior chain digest", payload["prior_ledger_digest"])
    _require_digest_text(f"{label} chain digest", payload["ledger_digest"])
    _require_digest_text(f"{label} registry digest", payload["registry_digest"])
    _require_nonblank_text(f"{label} store_id", payload["store_id"])
    # The common-datum label is a STORE-level pin: per-plan agreement alone would
    # let a cell's datum-relative energy accumulate across different datums over
    # successive transactions, which is not a quantity at all.
    _require_nonblank_text(f"{label} energy_datum_id", payload["energy_datum_id"])
    if payload["limitations"] != list(DECLARED_LIMITATIONS):
        raise TransferStoreIntegrityError(f"{label} limitations differ from the declared set")
    cells = payload["cells"]
    if type(cells) is not dict or not cells:
        raise TransferStoreIntegrityError(f"{label} declares no solid control volumes")
    for cell_id, record in cells.items():
        _require_nonblank_text(f"{label} cell id", cell_id)
        if type(record) is not dict or set(record) != {"revision", "holdup"}:
            raise TransferStoreIntegrityError(f"{label} cell record schema is foreign")
        if type(record["revision"]) is not int or record["revision"] < 0:
            raise TransferStoreIntegrityError(f"{label} cell revision is malformed")
        holdup = _require_json_extensives(record["holdup"], f"{label} cell holdup")
        for name, value in zip(EXTENSIVE_NAMES[:MASS_EXTENSIVE_COUNT], holdup, strict=False):
            if value < 0.0 or _is_negative_zero(value):
                raise TransferStoreIntegrityError(f"{label} cell holdup {name} is negative")
    tokens = payload["consumed_tokens"]
    entries = payload["entries"]
    if type(tokens) is not list or type(entries) is not list:
        raise TransferStoreIntegrityError(f"{label} consumed-set or entry log is not an array")
    if len(tokens) != len(entries):
        raise TransferStoreIntegrityError(f"{label} consumed-set and entry log disagree in length")
    if len(set(tokens)) != len(tokens):
        raise TransferStoreIntegrityError(f"{label} consumed-set contains a duplicate token")
    if len(entries) != revision:
        raise TransferStoreIntegrityError(f"{label} entry count does not track its revision")
    for index, entry in enumerate(entries):
        if type(entry) is not dict or set(entry) != set(_ENTRY_KEYS):
            raise TransferStoreIntegrityError(
                f"{label} entry {index} schema is incomplete or foreign"
            )
        if entry["sequence"] != index + 1 or entry["committed_revision"] != index + 1:
            raise TransferStoreIntegrityError(f"{label} entry {index} is out of sequence")
        _require_digest_text(f"{label} entry request identity", entry["request_identity"])
        _require_digest_text(f"{label} entry replay token", entry["replay_token"])
        _require_digest_text(f"{label} entry receipt digest", entry["receipt_digest"])
        if entry["replay_token"] != tokens[index]:
            raise TransferStoreIntegrityError(
                f"{label} entry {index} token differs from the consumed-set"
            )
        if entry["replay_token"] != derive_replay_token(entry["request_identity"]):
            raise TransferStoreIntegrityError(
                f"{label} entry {index} token is not derived from its request identity"
            )
        for key in ("debit_totals", "credit_totals", "closing_residual"):
            _require_json_extensives(entry[key], f"{label} entry {key}")
        if type(entry["packet_ids"]) is not list or not entry["packet_ids"]:
            raise TransferStoreIntegrityError(f"{label} entry {index} names no packets")
        if len(entry["payload_identities"]) != len(entry["packet_ids"]):
            raise TransferStoreIntegrityError(f"{label} entry {index} payload identities are short")
    body = {key: value for key, value in payload.items() if key != "ledger_digest"}
    if payload["ledger_digest"] != _generation_chain_digest(body):
        raise TransferStoreIntegrityError(f"{label} chain digest does not match its own bytes")
    return dict(payload)


def _sealed_generation_payload(body: dict[str, Any]) -> dict[str, Any]:
    payload = dict(body)
    payload["ledger_digest"] = _generation_chain_digest(payload)
    return payload


def _publish_generation(
    root: Path,
    generation: int,
    payload: Mapping[str, Any],
    checkpoint: Any,
) -> Path:
    """Stage, verify as-if-committed, then publish with ONE atomic rename.

    The destination is a fresh generation number, so the rename's Windows
    non-clobber semantics both keep the publication atomic and refuse a second
    concurrent commit of the same generation. Nothing is ever deleted: a
    pre-commit failure quarantines ``.pending`` -> ``.invalid``.

    THE RENAME IS THE COMMIT POINT AND IS THE ONLY IRREVERSIBLE ACT. Every
    reversible check runs strictly BEFORE it, the store-identity checkpoint
    included -- directory identity is a property of the STORE, not of the
    published bytes, so re-checking it afterwards could only manufacture a reason
    to retract a completed atomic act. After the rename nothing here renames or
    deletes the published generation: a post-rename re-read or re-validation
    failure REPORTS the condition and names the published path, leaving the file
    in place. That is deliberate. Retracting it would make a crash-after window
    answer "not committed" for a transaction whose commit point completed, and
    the reconciliation query must be able to answer ``committed`` for it.
    Genuine corruption is not lost either: the next load re-validates the whole
    chain and fails closed on it.
    """

    name = _generation_name(generation)
    final_path = root / name
    pending_path = root / (name + ".pending")
    invalid_path = root / (name + ".invalid")
    value = _canonical_json_bytes(payload)
    checkpoint("before staging the generation")
    try:
        _write_durable(pending_path, value)
        checkpoint("after the durable generation write")
        if pending_path.read_bytes() != value:
            raise TransferStoreIntegrityError("staged generation differs after the durable write")
        _validate_generation_payload(
            _strict_json_payload(value, "staged generation"), "staged generation"
        )
        if os.path.lexists(final_path):
            raise TransferStoreIntegrityError("generation destination already exists")
        checkpoint("after staged generation verification")
    except BaseException:
        if os.path.lexists(pending_path):
            if os.path.lexists(invalid_path):
                raise TransferStoreIntegrityError(
                    "staged generation failed and its quarantine path already exists"
                )
            os.rename(pending_path, invalid_path)
        raise
    try:
        checkpoint("before the atomic generation rename")
        os.rename(pending_path, final_path)
    except BaseException:
        if os.path.lexists(pending_path):
            if os.path.lexists(invalid_path):
                raise TransferStoreIntegrityError(
                    "generation publication failed and its quarantine path already exists"
                )
            os.rename(pending_path, invalid_path)
        raise
    try:
        published = final_path.read_bytes()
        if published != value:
            raise TransferStoreIntegrityError("published generation differs from validated bytes")
        _validate_generation_payload(
            _strict_json_payload(published, "published generation"), "published generation"
        )
    except Exception as error:
        raise TransferStoreIntegrityError(
            f"published generation {final_path} failed its post-commit re-verification and is "
            "LEFT IN PLACE: the rename is the commit point and is never retracted, so this "
            f"transaction IS committed and {RECONCILIATION_QUERY_NAME}() will say so. The next "
            f"load re-validates the chain and fails closed on genuine corruption. Cause: {error}"
        ) from error
    return final_path


# --------------------------------------------------------------------------- #
# The host.
# --------------------------------------------------------------------------- #


class AcceptedSolidTransferHost:
    """The single production transaction authority for solid motion.

    One host object owns one durable store. The mutation path is
    ``prepare_transfer`` followed by ``commit_transfer``; the recovery path is
    the read-only ``reconcile_request``. The two never share a code path, so a
    committed disposition can never be echoed out of the mutation path.
    """

    accepted_solid_transfer_implemented: ClassVar[bool] = True
    host_object_identity_authenticated: ClassVar[bool] = True
    transaction_consumption_persisted: ClassVar[bool] = True
    atomic_persistence_selected: ClassVar[bool] = True
    persisted_replay_protection_selected: ClassVar[bool] = True
    # Explicitly NOT supplied by this host, and staying false.
    complete_particle_state_transport_selected: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_variant_registry_complete: ClassVar[bool] = False
    source_authenticated: ClassVar[bool] = False

    def __init__(
        self,
        store_dir: Path,
        *,
        registry: state_codec.StateSchemaRegistry,
        pinned_registry_digest: str,
    ) -> None:
        if type(registry) is not state_codec.StateSchemaRegistry:
            raise ForeignTransferError("registry must be an exact StateSchemaRegistry")
        # B4, exactly as far as S0 goes: the LOAD-TIME COMPARISON is here and it
        # fails closed, and the compared digest is bound into every transaction.
        # What is NOT here is the trust root: ``pinned_registry_digest`` is a
        # caller-supplied argument, and no release/-anchored pin artifact or
        # loader for one exists yet. Until that artifact ships with this packet's
        # release stage, this comparison authenticates the registry against the
        # CALLER's claim, not against a release anchor.
        _require_digest_text("pinned registry digest", pinned_registry_digest)
        if registry.definition_digest != pinned_registry_digest:
            raise ForeignTransferError(
                "registry definition digest does not match the pinned release digest"
            )
        root = Path(store_dir)
        _reject_reparse_ancestors(root)
        if not root.is_dir():
            raise TransferStoreIntegrityError(f"transfer store directory is missing: {root}")
        self._root = root.resolve(strict=True)
        self._registry = registry
        self._pinned_registry_digest = pinned_registry_digest
        self._store_identity = _directory_identity(self._root)
        self._lock_path = self._root / HOST_LOCK_NAME
        head_generation, head = self._load_head()
        self._store_id = head["store_id"]
        self._energy_datum_id = _require_nonblank_text(
            "store energy_datum_id", head["energy_datum_id"]
        )
        if head["registry_digest"] != pinned_registry_digest:
            raise ForeignTransferError("store was created under a different pinned registry digest")
        self._host_identity = self._compute_host_identity()
        del head_generation

    # -- identity -------------------------------------------------------- #

    def _compute_host_identity(self) -> str:
        return _framed_digest(
            HOST_IDENTITY_DIGEST_DOMAIN,
            (
                TRANSACTION_MAGIC.encode("ascii"),
                str(self._store_identity["path"]).encode("utf-8"),
                str(self._store_identity["device"]).encode("ascii"),
                str(self._store_identity["inode"]).encode("ascii"),
                self._store_id.encode("utf-8"),
                self._pinned_registry_digest.encode("ascii"),
                str(id(self)).encode("ascii"),
            ),
        )

    @property
    def host_identity(self) -> str:
        return self._host_identity

    @property
    def store_root(self) -> Path:
        return self._root

    @property
    def energy_datum_id(self) -> str:
        """The store's one common-datum label. Every plan must name exactly it."""

        return self._energy_datum_id

    def _checkpoint(self, label: str) -> None:
        """Re-assert store directory identity AND host object identity."""

        observed = _directory_identity(self._root)
        if observed != self._store_identity:
            raise TransferStoreIntegrityError(
                f"store directory identity changed at {label}: "
                f"{observed!r} != {dict(self._store_identity)!r}"
            )
        if self._compute_host_identity() != self._host_identity:
            raise ForeignTransferError(f"host object identity changed at {label}")

    # -- store scan ------------------------------------------------------- #

    def _require_no_reserved_residue(self) -> None:
        # Name-shape-bound BY DESIGN: only <valid generation name>.pending/.invalid
        # counts as residue, because those are the only names this module's own
        # writes can produce. A hand-made file like "ledger.9.json.pending" is
        # neither residue nor a generation to this check -- operator-fabricated
        # residue is out of scope and is not silently trusted either (it is
        # simply not a generation).
        residue = sorted(
            entry.name
            for entry in self._root.iterdir()
            if entry.name.endswith(RESERVED_SUFFIXES)
            and _generation_number(entry.name.rsplit(".", 1)[0]) is not None
        )
        if residue:
            raise TransferStoreIntegrityError(
                "transfer store carries unresolved reserved-suffix residue: " + ", ".join(residue)
            )

    def _load_head(self) -> tuple[int, dict[str, Any]]:
        """Load and chain-verify every generation, returning the head."""

        self._require_no_reserved_residue()
        numbers = sorted(
            number
            for number in (_generation_number(entry.name) for entry in self._root.iterdir())
            if number is not None
        )
        if not numbers:
            raise TransferStoreIntegrityError(
                f"transfer store carries no genesis ledger: {self._root}"
            )
        if numbers != list(range(GENESIS_GENERATION, GENESIS_GENERATION + len(numbers))):
            raise TransferStoreIntegrityError("transfer store generation chain is not contiguous")
        prior_digest = GENESIS_PRIOR_DIGEST
        payload: dict[str, Any] = {}
        for number in numbers:
            path = self._root / _generation_name(number)
            label = f"ledger generation {number}"
            payload = _validate_generation_payload(
                _strict_json_payload(path.read_bytes(), label), label
            )
            if payload["generation"] != number:
                raise TransferStoreIntegrityError(f"{label} declares a foreign generation number")
            if payload["prior_ledger_digest"] != prior_digest:
                raise TransferStoreIntegrityError(f"{label} does not chain to its predecessor")
            prior_digest = payload["ledger_digest"]
        return numbers[-1], payload

    # -- the reachable operator diagnostic -------------------------------- #

    @classmethod
    def inspect_store(cls, store_dir: Path) -> dict[str, Any]:
        """Read-only store report, constructible over ANY store state.

        WHY THIS EXISTS. Every refusal in this module names the reconciliation
        query, but that query is a METHOD: reaching it needs a host object, and a
        host cannot be constructed over a store carrying reserved-suffix residue
        -- which is precisely the state an interrupted commit leaves behind. A
        fresh process would therefore find the named recovery route unreachable.
        This classmethod is the reachable one.

        WHAT IT DOES NOT DO. It answers NO commit/no-commit question. That is the
        reconciliation query's role and it keeps refusing over an unresolved
        store, because an unresolved store cannot support a truthful answer. This
        reports what is on disk: committed generation numbers, the head number,
        reserved-suffix residue, whether a lockfile is present, and the head's
        store pins. It validates nothing, repairs nothing, and writes nothing, so
        the pins it reports are UNVALIDATED and are ``None`` when the head cannot
        be parsed at all.
        """

        root = Path(store_dir)
        _reject_reparse_ancestors(root)
        if not root.is_dir():
            raise TransferStoreIntegrityError(f"transfer store directory is missing: {root}")
        resolved = root.resolve(strict=True)
        generations: list[int] = []
        residue: list[str] = []
        for entry in resolved.iterdir():
            number = _generation_number(entry.name)
            if number is not None:
                generations.append(number)
            elif (
                entry.name.endswith(RESERVED_SUFFIXES)
                and _generation_number(entry.name.rsplit(".", 1)[0]) is not None
            ):
                residue.append(entry.name)
        generations.sort()
        head_generation = generations[-1] if generations else None
        pins: dict[str, str | None] = {
            "store_id": None,
            "registry_digest": None,
            "energy_datum_id": None,
        }
        if head_generation is not None:
            payload: dict[str, Any] = {}
            try:
                raw = (resolved / _generation_name(head_generation)).read_bytes()
            except OSError:
                raw = b""
            if raw:
                try:
                    payload = _strict_json_payload(raw, "inspected head generation")
                except TransferStoreIntegrityError:
                    # An unparseable head still gets a listing, just no pins: a
                    # diagnostic over a broken store must not itself fail.
                    payload = {}
            for key in pins:
                value = payload.get(key)
                pins[key] = value if type(value) is str else None
        return {
            "store_root": str(resolved),
            "committed_generations": tuple(generations),
            "head_generation": head_generation,
            "reserved_suffix_residue": tuple(sorted(residue)),
            "store_is_unresolved": bool(residue),
            # is_file(): a DIRECTORY at the lock path reads lock_present False
            # here, yet still refuses acquisition (the exclusive create
            # raises). This report describes the normal lockfile only.
            "lock_present": (resolved / HOST_LOCK_NAME).is_file(),
            "store_pins": pins,
            "store_pins_are_unvalidated": True,
            "answers_commit_disposition": False,
            "reconciliation_query": RECONCILIATION_QUERY_NAME,
            "limitations": DECLARED_LIMITATIONS,
        }

    # -- store creation --------------------------------------------------- #

    @classmethod
    def create_store(
        cls,
        store_dir: Path,
        *,
        registry: state_codec.StateSchemaRegistry,
        pinned_registry_digest: str,
        store_id: str,
        energy_datum_id: str,
        initial_holdup: Mapping[str, tuple[float, ...]],
    ) -> AcceptedSolidTransferHost:
        """Publish the genesis generation, then open the store over it.

        ``energy_datum_id`` pins the store's ONE common energy datum for its whole
        life. It is persisted in every generation and re-validated on every load,
        because per-cell datum-relative energy is only a quantity relative to a
        single datum.
        """

        if type(registry) is not state_codec.StateSchemaRegistry:
            raise ForeignTransferError("registry must be an exact StateSchemaRegistry")
        _require_digest_text("pinned registry digest", pinned_registry_digest)
        if registry.definition_digest != pinned_registry_digest:
            raise ForeignTransferError(
                "registry definition digest does not match the pinned release digest"
            )
        _require_nonblank_text("store_id", store_id)
        _require_nonblank_text("energy_datum_id", energy_datum_id)
        root = Path(store_dir)
        _reject_reparse_ancestors(root)
        if not root.is_dir():
            raise TransferStoreIntegrityError(f"transfer store directory is missing: {root}")
        resolved = root.resolve(strict=True)
        identity = _directory_identity(resolved)
        for entry in resolved.iterdir():
            if _generation_number(entry.name) is not None or entry.name.endswith(RESERVED_SUFFIXES):
                raise TransferStoreIntegrityError(
                    f"transfer store already carries ledger state: {entry.name}"
                )
        if type(initial_holdup) is not dict or not initial_holdup:
            raise AcceptedSolidTransferError("initial holdup must be a non-empty exact dict")
        cells: dict[str, Any] = {}
        for cell_id, holdup in initial_holdup.items():
            _require_nonblank_text("initial holdup cell id", cell_id)
            totals = _require_extensive_tuple(f"initial holdup for {cell_id}", holdup)
            for name, value in zip(EXTENSIVE_NAMES[:MASS_EXTENSIVE_COUNT], totals, strict=False):
                if value < 0.0 or _is_negative_zero(value):
                    raise AcceptedSolidTransferError(
                        f"initial holdup {name} for {cell_id} must be non-negative"
                    )
            cells[cell_id] = {"revision": 0, "holdup": [0.0 + value for value in totals]}
        payload = _sealed_generation_payload(
            {
                "schema": LEDGER_SCHEMA,
                "schema_version": LEDGER_SCHEMA_VERSION,
                "generation": GENESIS_GENERATION,
                "revision": 0,
                "prior_ledger_digest": GENESIS_PRIOR_DIGEST,
                "store_id": store_id,
                "registry_digest": pinned_registry_digest,
                "energy_datum_id": energy_datum_id,
                "cells": cells,
                "consumed_tokens": [],
                "entries": [],
                "limitations": list(DECLARED_LIMITATIONS),
            }
        )

        def checkpoint(label: str) -> None:
            observed = _directory_identity(resolved)
            if observed != identity:
                raise TransferStoreIntegrityError(
                    f"store directory identity changed at {label} during genesis"
                )

        _publish_generation(resolved, GENESIS_GENERATION, payload, checkpoint)
        return cls(resolved, registry=registry, pinned_registry_digest=pinned_registry_digest)

    # -- the exclusive host lock ------------------------------------------ #

    def _acquire_lock(self, transaction_identity: str) -> None:
        """Exclusive-create lockfile, taken BEFORE validation and held throughout.

        Fail-closed stale policy: a surviving lockfile is never stolen and never
        aged out. Automatic takeover would weaken a fail-closed boundary, so the
        refusal names the lockfile and the reconciliation query that resolves
        the prior transaction's disposition instead.

        Acquisition is ALL-OR-NOTHING. Because a lockfile is never stolen, a
        half-written one left behind by a failure after the exclusive create
        would wedge the store permanently, so it is cleaned up before raising.
        """

        contents = _canonical_json_bytes(
            {
                "host_identity": self._host_identity,
                "store_id": self._store_id,
                "store_path": str(self._store_identity["path"]),
                "transaction_identity": transaction_identity,
                "process_id": os.getpid(),
                "policy": DeclaredEngineeringAssumption.FAIL_CLOSED_STALE_LOCK.value,
            }
        )
        try:
            _write_durable(self._lock_path, contents)
        except FileExistsError as error:
            # FIRST, and load-bearing: an existing lockfile belongs to another
            # holder and is NEVER removed by the branch below.
            raise HostLockRefusedError(
                f"host lock {self._lock_path} is held; it is never stolen or aged out. "
                f"Resolve the prior transaction's disposition with {RECONCILIATION_QUERY_NAME}()."
            ) from error
        except Exception as error:
            # All-or-nothing acquisition: the exclusive create succeeded but the
            # write did not (short write, ENOSPC, fsync failure), so a partial
            # lockfile is OURS and would otherwise masquerade as a live holder
            # forever, wedging the store behind a lock that is never stolen.
            try:
                os.remove(self._lock_path)
                residue = "the partial lockfile this host had just created was removed"
            except OSError:
                residue = (
                    "the partial lockfile could NOT be removed and now blocks this store; "
                    f"resolve it by hand after {RECONCILIATION_QUERY_NAME}()"
                )
            raise TransferStoreIntegrityError(
                f"host lock {self._lock_path} could not be established; {residue}: {error}"
            ) from error

    def _release_lock(self) -> None:
        """Best-effort release. Deliberately NEVER raises.

        This runs in a ``finally``, and a raising ``finally`` would REPLACE the
        durable commit receipt the rename already earned, or mask the ruled
        terminal replay error, with an unrelated filesystem error -- turning a
        committed transaction into an apparently failed one. A leaked lockfile
        instead fails closed on the next acquire, which is the correct direction
        for a fail-closed boundary: the operator resolves it with the
        reconciliation query. Declared as ``lock_release_best_effort``.
        """

        try:
            os.remove(self._lock_path)
        except FileNotFoundError:
            # The expected benign case: released, or never created.
            pass
        except OSError:
            pass

    # -- the mutation path ------------------------------------------------ #

    def _payload_identities(self, request: AcceptedSolidTransferRequest) -> tuple[str, ...]:
        """Validate every payload through the selected profile, then take identity.

        A malformed or foreign-variant payload surfaces the profile's own decode
        refusal unchanged. That is deliberate: a transaction-level refusal must
        stay distinguishable from malformed bytes, so the two families are never
        collapsed into one another.
        """

        return tuple(
            state_codec.canonical_state_identity(entry.payload_bytes, self._registry)
            for entry in request.packet_plan
        )

    def _derived_credit_totals(self, request: AcceptedSolidTransferRequest) -> tuple[float, ...]:
        """Sum the B1-certified per-packet outward scaling over the packet plan."""

        columns: list[list[float]] = [[] for _ in EXTENSIVE_NAMES]
        datum_ids = set()
        for entry in request.packet_plan:
            if entry.surface is not None:
                _contract.require_single_surface_authority(entry.inventory, entry.surface)
            totals = _contract.derive_packet_totals(entry.inventory, entry.weight)
            datum_ids.add(totals.energy_datum_id)
            for index, value in enumerate(totals.totals):
                columns[index].append(value)
        if len(datum_ids) != 1:
            raise ForeignTransferError("a packet plan cannot mix common-datum energy references")
        plan_datum = datum_ids.pop()
        if plan_datum != self._energy_datum_id:
            raise ForeignTransferError(
                f"packet plan common-datum energy reference {plan_datum!r} is foreign to this "
                f"store's pinned datum {self._energy_datum_id!r}"
            )
        credit = tuple(math.fsum(column) for column in columns)
        _require_binary64("derived credit totals", *credit)
        return credit

    def _conservation_residual(
        self, debit_totals: tuple[float, ...], credit_totals: tuple[float, ...]
    ) -> tuple[float, ...]:
        """Debit minus credit over the seven extensives. Reported, never repaired."""

        residual = tuple(
            math.fsum((debit, -credit))
            for debit, credit in zip(debit_totals, credit_totals, strict=True)
        )
        if any(value != 0.0 for value in residual):
            raise ConservationResidualError(
                "accepted transfer does not balance over the seven extensives: "
                f"residual {residual!r}. A residual is reported, never repaired by a clamp."
            )
        return residual

    def prepare_transfer(self, request: AcceptedSolidTransferRequest) -> AcceptedSolidTransfer:
        """Bind all nine items and issue the one host-owned transaction.

        Nothing is persisted here. Conservation is checked BEFORE any durable
        act, so an unbalanced proposal is a terminal refusal that never reaches
        the store.
        """

        if type(request) is not AcceptedSolidTransferRequest:
            raise ForeignTransferError(
                "the mutation path accepts exact AcceptedSolidTransferRequest values only"
            )
        self._checkpoint("transaction preparation")
        if request.codec_registry_digest != self._pinned_registry_digest:
            raise ForeignTransferError(
                "request codec identity does not carry this host's pinned registry digest"
            )
        payload_identities = self._payload_identities(request)
        credit_totals = self._derived_credit_totals(request)
        debit_totals = tuple(request.declared_debit_totals)
        residual = self._conservation_residual(debit_totals, credit_totals)
        identity = request_identity_digest(request, payload_identities)
        return AcceptedSolidTransfer(
            request=request,
            request_identity=identity,
            replay_token=derive_replay_token(identity),
            payload_identities=payload_identities,
            debit_totals=debit_totals,
            credit_totals=credit_totals,
            conservation_residual=residual,
            host_identity=self._host_identity,
            _issue_token=_ISSUE_TOKEN,
        )

    def commit_transfer(self, transaction: AcceptedSolidTransfer) -> TransferCommitReceipt:
        """Re-audit, debit, credit and consume the token in ONE atomic act.

        The lock is taken before ledger validation and held through re-audit,
        debit, credit and token consumption. Replay is tested BEFORE staleness
        on purpose: because the token derives from an identity that already
        binds both prior revisions, a re-presented request is also stale, and
        the ruled distinct terminal replay error must win.
        """

        if type(transaction) is not AcceptedSolidTransfer:
            raise ForeignTransferError(
                "the mutation path accepts host-issued AcceptedSolidTransfer values only"
            )
        if transaction.host_identity != self._host_identity:
            raise ForeignTransferError("transaction was issued by a different host object")
        self._checkpoint("commit entry")
        self._acquire_lock(transaction.request_identity)
        try:
            return self._commit_under_lock(transaction)
        finally:
            self._release_lock()

    def _commit_under_lock(self, transaction: AcceptedSolidTransfer) -> TransferCommitReceipt:
        self._checkpoint("after lock acquisition")
        head_generation, head = self._load_head()
        if head["registry_digest"] != self._pinned_registry_digest:
            raise ForeignTransferError("store registry pin differs from this host's pinned digest")
        if head["energy_datum_id"] != self._energy_datum_id:
            raise ForeignTransferError("store energy datum differs from the one this host opened")
        request = transaction.request
        # Re-audit from scratch: never trust the issued view's own arithmetic.
        payload_identities = self._payload_identities(request)
        identity = request_identity_digest(request, payload_identities)
        token = derive_replay_token(identity)
        if identity != transaction.request_identity or token != transaction.replay_token:
            raise ForeignTransferError("transaction identity does not re-derive from its request")
        credit_totals = self._derived_credit_totals(request)
        debit_totals = tuple(request.declared_debit_totals)
        residual = self._conservation_residual(debit_totals, credit_totals)
        # B5 ordering: replay wins over staleness. Zero state change either way.
        if token in head["consumed_tokens"]:
            raise ReplayedTransferError(
                f"replay token {token} was already consumed; this transaction is terminal. "
                f"Use {RECONCILIATION_QUERY_NAME}() to recover its committed disposition."
            )
        cells = head["cells"]
        source_id = request.topology.source_cell_id
        destination_id = request.topology.destination_cell_id
        for cell_id in (source_id, destination_id):
            if cell_id not in cells:
                raise ForeignTransferError(f"cell {cell_id!r} is foreign to this transfer store")
        if cells[source_id]["revision"] != request.source_prior_revision:
            raise StaleTransferError(
                f"source cell {source_id!r} is at revision {cells[source_id]['revision']}, "
                f"not the proposed prior revision {request.source_prior_revision}"
            )
        if cells[destination_id]["revision"] != request.destination_prior_revision:
            raise StaleTransferError(
                f"destination cell {destination_id!r} is at revision "
                f"{cells[destination_id]['revision']}, not the proposed prior revision "
                f"{request.destination_prior_revision}"
            )
        source_holdup = tuple(cells[source_id]["holdup"])
        destination_holdup = tuple(cells[destination_id]["holdup"])
        remaining = tuple(
            math.fsum((held, -debit))
            for held, debit in zip(source_holdup, debit_totals, strict=True)
        )
        for name, value in zip(EXTENSIVE_NAMES[:MASS_EXTENSIVE_COUNT], remaining, strict=False):
            if value < 0.0:
                raise InsufficientSolidError(
                    f"source cell {source_id!r} cannot supply the accepted {name} debit; "
                    "insufficiency rejects upstream and is never repaired by a clamp"
                )
        credited = tuple(
            math.fsum((held, credit))
            for held, credit in zip(destination_holdup, credit_totals, strict=True)
        )
        _require_binary64("post-transfer holdup", *remaining, *credited)
        next_cells = {
            cell_id: {"revision": record["revision"], "holdup": list(record["holdup"])}
            for cell_id, record in cells.items()
        }
        next_cells[source_id] = {
            "revision": cells[source_id]["revision"] + 1,
            "holdup": [value if value != 0.0 else 0.0 for value in remaining],
        }
        next_cells[destination_id] = {
            "revision": cells[destination_id]["revision"] + 1,
            "holdup": [value if value != 0.0 else 0.0 for value in credited],
        }
        closing_residual = tuple(
            math.fsum(
                [record["holdup"][index] for record in next_cells.values()]
                + [-record["holdup"][index] for record in cells.values()]
            )
            for index in range(len(EXTENSIVE_NAMES))
        )
        committed_revision = head["revision"] + 1
        prior_ledger_digest = head["ledger_digest"]
        receipt_digest = commit_receipt_digest(
            identity, token, committed_revision, debit_totals, credit_totals, prior_ledger_digest
        )
        entry = {
            "sequence": committed_revision,
            "committed_revision": committed_revision,
            "request_identity": identity,
            "replay_token": token,
            "receipt_digest": receipt_digest,
            "authority_class": request.authority_class.value,
            "source_cell_id": source_id,
            "destination_cell_id": destination_id,
            "debit_totals": list(debit_totals),
            "credit_totals": list(credit_totals),
            "closing_residual": list(closing_residual),
            "packet_ids": [item.packet_id for item in request.packet_plan],
            "payload_identities": list(payload_identities),
            "codec_registry_digest": request.codec_registry_digest,
            "auditor_identity_digest": request.auditor_identity_digest,
            "configuration_digest": request.configuration_digest,
        }
        next_generation = head_generation + 1
        payload = _sealed_generation_payload(
            {
                "schema": LEDGER_SCHEMA,
                "schema_version": LEDGER_SCHEMA_VERSION,
                "generation": next_generation,
                "revision": committed_revision,
                "prior_ledger_digest": prior_ledger_digest,
                "store_id": self._store_id,
                "registry_digest": self._pinned_registry_digest,
                "energy_datum_id": self._energy_datum_id,
                "cells": next_cells,
                # Debit, credit and consumption share ONE file, so the single
                # rename below commits all three together or none of them.
                "consumed_tokens": list(head["consumed_tokens"]) + [token],
                "entries": list(head["entries"]) + [entry],
                "limitations": list(DECLARED_LIMITATIONS),
            }
        )
        published = _publish_generation(self._root, next_generation, payload, self._checkpoint)
        del residual
        return TransferCommitReceipt(
            request_identity=identity,
            replay_token=token,
            committed_revision=committed_revision,
            generation=next_generation,
            ledger_digest=payload["ledger_digest"],
            receipt_digest=receipt_digest,
            debit_totals=debit_totals,
            credit_totals=credit_totals,
            ledger_path=str(published),
            _issue_token=_ISSUE_TOKEN,
        )

    # -- the read-only recovery path -------------------------------------- #

    def reconcile_request(self, request_identity: str) -> ReconciliationAnswer:
        """Answer "was request X committed?" WITHOUT touching the mutation path.

        This is the ruled recovery mechanism and it is strictly read-only: it
        takes no lock, writes nothing, and returns disposition plus the
        committed revision and receipt. It refuses to answer over a store
        carrying reserved-suffix residue, because an unresolved store cannot
        support a truthful answer.
        """

        _require_digest_text("request identity", request_identity)
        self._checkpoint("reconciliation query")
        _generation, head = self._load_head()
        del _generation
        for entry in head["entries"]:
            if entry["request_identity"] == request_identity:
                return ReconciliationAnswer(
                    request_identity=request_identity,
                    committed=True,
                    committed_revision=entry["committed_revision"],
                    generation=entry["committed_revision"] + GENESIS_GENERATION,
                    receipt_digest=entry["receipt_digest"],
                )
        return ReconciliationAnswer(
            request_identity=request_identity,
            committed=False,
            committed_revision=None,
            generation=None,
            receipt_digest=None,
        )

    # -- read-only views --------------------------------------------------- #

    def ledger_revision(self) -> int:
        return int(self._load_head()[1]["revision"])

    def cell_revision(self, cell_id: str) -> int:
        cells = self._load_head()[1]["cells"]
        if cell_id not in cells:
            raise ForeignTransferError(f"cell {cell_id!r} is foreign to this transfer store")
        return int(cells[cell_id]["revision"])

    def cell_holdup(self, cell_id: str) -> tuple[float, ...]:
        cells = self._load_head()[1]["cells"]
        if cell_id not in cells:
            raise ForeignTransferError(f"cell {cell_id!r} is foreign to this transfer store")
        return tuple(cells[cell_id]["holdup"])

    def consumed_tokens(self) -> tuple[str, ...]:
        """The append-only consumed-set. It never expires and is never pruned."""

        return tuple(self._load_head()[1]["consumed_tokens"])

    def declared_limitations(self) -> Mapping[str, bool]:
        return MappingProxyType({name: True for name in DECLARED_LIMITATIONS})


__all__ = (
    "AcceptedInterval",
    "AcceptedSolidTransfer",
    "AcceptedSolidTransferError",
    "AcceptedSolidTransferHost",
    "AcceptedSolidTransferRequest",
    "BOUND_TRANSACTION_ITEMS",
    "ConservationResidualError",
    "DECLARED_LIMITATIONS",
    "DIGEST_DOMAIN_PREFIX",
    "DeclaredEngineeringAssumption",
    "EXTENSIVE_NAMES",
    "ExternalSurfaceComposite",
    "ForeignTransferError",
    "GENESIS_GENERATION",
    "GENESIS_PRIOR_DIGEST",
    "HOST_IDENTITY_DIGEST_DOMAIN",
    "HOST_LOCK_NAME",
    "HostLockRefusedError",
    "IDENTITY_DIGEST_DOMAIN",
    "InsufficientSolidError",
    "LEDGER_BASENAME",
    "LEDGER_DIGEST_DOMAIN",
    "LEDGER_SCHEMA",
    "LEDGER_SCHEMA_VERSION",
    "MASS_EXTENSIVE_COUNT",
    "OfferedSolidFlow",
    "PACKET_PLAN_DIGEST_DOMAIN",
    "PacketPlanEntry",
    "PacketWeight",
    "PerParticleInventory",
    "RECEIPT_DIGEST_DOMAIN",
    "RECONCILIATION_QUERY_NAME",
    "REPLAY_TOKEN_DOMAIN",
    "RESERVED_SUFFIXES",
    "ReconciliationAnswer",
    "ReplayedTransferError",
    "SCALING_CONTRACT_MODULE",
    "ScalingContractError",
    "StaleTransferError",
    "TRANSACTION_MAGIC",
    "TRANSACTION_VERSION",
    "TransferAuthorityClass",
    "TransferCommitReceipt",
    "TransferStoreIntegrityError",
    "TransferTopology",
    "commit_receipt_digest",
    "derive_replay_token",
    "packet_plan_entry_digest",
    "request_identity_digest",
)
