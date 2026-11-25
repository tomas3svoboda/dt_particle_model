"""A compiled IAPWS-95 water residual, bit-identical to the Python one by construction (2026-11-03).

Record: docs/GT_PS2_DTDC_KERNEL_IAPWS_COMPILED_BIT_IDENTITY_2026-11-03.md (hint C of the external
parallelization study, briefed by the lead; a measurement and a prototype behind the kernel-envelope
knob ``DTDC_KERNEL_ENVELOPE_IAPWS_COMPILED=1``; not adopted).

What it is.  ``dtdc_simulator.core2.props.water._residual_uncached`` (the 56-term reduced residual
Helmholtz energy with its eight derivatives, and the ideal part ``_phi0``) compiled to native code
through llvmlite, which the project venv already carries (this host has no C compiler).  The machine
code is generated from the SAME Python source text: the two functions' source is parsed, two declared
rewrites are applied (below), and the result is executed once with traced operands, so every IEEE
operation is emitted in the order and association Python evaluates it.  Python-level constant
sub-expressions (``2.0 * C``, ``1.0 / (2.0 * be)``, ``_G0[i] ** 2``, ``t * (t - 1.0)``) are evaluated by
this same interpreter at trace time, i.e. with the very operations the interpreted residual performs
on every call.

The operations reaching the C runtime are exactly the interpreter's: ``math.exp`` and ``math.log``
call the C runtime's ``exp`` and ``log``; ``float ** float`` calls its ``pow`` after CPython's own
special-case handling (``float_pow``: exponent zero, a zero, negative or unit base), which is mirrored
here operation for operation.  The emitted calls go to neutrally named symbols bound BY ADDRESS to
``ucrtbase.dll``'s ``exp``, ``log`` and ``pow``, the functions the interpreter's DLL imports through
the ``api-ms-win-crt-math-l1-1-0`` forwarder; the binding is verified at build time by reading back
the addresses the generated code itself holds.

Code generation: no IR optimisation pass is run; no fast-math flag is set on any instruction; the
target machine is the generic ``x86-64`` CPU with an empty feature string (no FMA instructions exist
for it) and LLVM's default ``AllowFPOpFusion=Standard``, which fuses only the explicit ``fmuladd``
intrinsic (never emitted here); ``fmul`` then ``fadd`` stay two rounded operations.  A known
FMA-sensitive expression is compiled with the same machine and checked at build time.

Every situation in which the interpreted residual would raise, or would leave the plain-double
domain, sets a flag instead and the call falls back to the Python residual itself: a division by
zero; ``math.exp`` or ``math.log`` outside their finite non-raising domain; a ``pow`` whose base is
zero, non-finite, or negative with a non-integer exponent, or whose result is not finite; any
non-finite output.  A fallback returns the interpreter's own result (or its own exception).

The two declared rewrites of the traced source: (1) the ridge guard ``if dm1 == 0.0: dm1 = 1.0e-12``
becomes the branch-free select of the same values (it is a data-dependent branch); (2) ``math`` and
``_phi0`` resolve to the traced versions.  Anything else in the source that cannot be traced
(a data-dependent branch, an unsupported operation) refuses the build.

physically_qualifying: false.
"""

from __future__ import annotations

import ast
import ctypes
import inspect
import math
import os
import struct
import textwrap
import time

RECORD = "docs/GT_PS2_DTDC_KERNEL_IAPWS_COMPILED_BIT_IDENTITY_2026-11-03.md"
SYMBOL_PREFIX = "dtdc_ucrt_"
RUNTIME_FUNCTIONS = ("exp", "log", "pow")
FUNCTION_NAME = "dtdc_iapws95_residual"
TARGET_CPU = "x86-64"
TARGET_FEATURES = ""
CODEGEN_OPT_LEVEL = 2
OUTPUT_FIELDS = (
    "phir",
    "phir_d",
    "phir_dd",
    "phir_t",
    "phir_tt",
    "phir_dt",
    "phir_dtt",
    "phir_dttt",
    "phi0",
    "phi0_t",
    "phi0_tt",
)
MUTATIONS = (None, "reversed_term_order")
#: The water module's coefficient tuples the trace bakes in; checked by identity at every call.
TRACED_CONSTANTS = (
    "_N0",
    "_G0",
    "_C",
    "_D",
    "_T",
    "_NR",
    "_G_D",
    "_G_T",
    "_G_N",
    "_G_ALPHA",
    "_G_BETA",
    "_G_GAMMA",
    "_G_EPS",
    "_NA_A",
    "_NA_B",
    "_NA_BB",
    "_NA_N",
    "_NA_C",
    "_NA_D",
    "_NA_AA",
    "_NA_BETA",
)


def _constants_guard(constants):
    """A function of the water module's ``__dict__``: True while every traced tuple is still the
    module's own object (a rebinding, such as a datum shift of ``_N0``, makes it False)."""

    names = tuple(constants)
    objects = tuple(constants[n] for n in names)
    code = "lambda d: " + " and ".join(f"d[{n!r}] is c[{i}]" for i, n in enumerate(names))
    return eval(code, {"c": objects})  # noqa: S307 - generated from the fixed name tuple


class TraceRefused(RuntimeError):
    """The residual's source carries something this tracer does not reproduce exactly."""


# ---------------------------------------------------------------------------------------------
# The tracer: one straight-line basic block, every operation an unflagged IEEE instruction.
# ---------------------------------------------------------------------------------------------
class _Tracer:
    def __init__(self, ir, module, reuse=True):
        self.ir = ir
        # Exact reuse: an operation already emitted on the same operands (the same SSA values or
        # constants of the same bits) is not emitted twice; its first result is used again.  Every
        # operation here is a deterministic function of its operands' bits (an IEEE instruction in
        # the default rounding mode, or ucrtbase's exp, log, pow in this process), so the reused
        # value is the value the repeated operation would return.  The interpreter repeats them
        # (``delta**2`` once per term); the generated code computes each once.
        self.reuse = reuse
        self.memo = {}
        self.reused = 0
        self.double = ir.DoubleType()
        self.module = module
        fnty = ir.FunctionType(self.double, [self.double])
        powty = ir.FunctionType(self.double, [self.double, self.double])
        self.exp = ir.Function(module, fnty, SYMBOL_PREFIX + "exp")
        self.log = ir.Function(module, fnty, SYMBOL_PREFIX + "log")
        self.pow = ir.Function(module, powty, SYMBOL_PREFIX + "pow")
        out_ty = ir.PointerType(self.double)
        sig = ir.FunctionType(ir.IntType(32), [self.double, self.double, out_ty])
        self.fn = ir.Function(module, sig, FUNCTION_NAME)
        self.builder = ir.IRBuilder(self.fn.append_basic_block("entry"))
        self.bad = ir.Constant(ir.IntType(1), 0)
        self.counts = {"fadd": 0, "fsub": 0, "fmul": 0, "fdiv": 0, "fneg": 0, "exp": 0, "log": 0}
        self.counts.update({"pow_calls": 0, "pow_constant_one": 0, "select": 0})

    def const(self, value):
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise TraceRefused(f"unsupported constant {value!r}")
        return self.ir.Constant(self.double, float(value))

    def flag(self, condition):
        self.bad = self.builder.or_(self.bad, condition)

    def finite(self, value):
        b = self.builder
        lo = b.fcmp_ordered(">", value, self.const(-math.inf))
        hi = b.fcmp_ordered("<", value, self.const(math.inf))
        return b.and_(lo, hi)

    def not_(self, condition):
        return self.builder.xor(condition, self.ir.Constant(self.ir.IntType(1), 1))

    def _key(self, *operands):
        out = []
        for v in operands:
            if isinstance(v, self.ir.Constant):
                out.append(("c", struct.pack("<d", float(v.constant))))
            else:
                out.append(("v", id(v)))
        return tuple(out)

    def _memo(self, name, operands, emit):
        if not self.reuse:
            return emit()
        key = (name, self._key(*operands))
        hit = self.memo.get(key)
        if hit is not None:
            self.reused += 1
            return hit[0]
        value = emit()
        self.memo[key] = (value, operands)  # the operands kept alive: their ids stay unique
        return value

    def binary(self, op, left, right):
        return self._memo(op, (left, right), lambda: self._binary(op, left, right))

    def _binary(self, op, left, right):
        self.counts[op] += 1
        b = self.builder
        if op == "fdiv":
            # float_div raises ZeroDivisionError on a zero divisor (either sign)
            self.flag(b.fcmp_ordered("==", right, self.const(0.0)))
        return getattr(b, op)(left, right)

    def neg(self, value):
        return self._memo("fneg", (value,), lambda: self._neg(value))

    def _neg(self, value):
        self.counts["fneg"] += 1
        return self.builder.fneg(value)

    def exp_(self, value):
        return self._memo("exp", (value,), lambda: self._exp(value))

    def _exp(self, value):
        self.counts["exp"] += 1
        r = self.builder.call(self.exp, [value])
        self.flag(self.not_(self.finite(value)))
        self.flag(self.not_(self.finite(r)))
        return r

    def log_(self, value):
        return self._memo("log", (value,), lambda: self._log(value))

    def _log(self, value):
        self.counts["log"] += 1
        b = self.builder
        self.flag(
            self.not_(b.and_(b.fcmp_ordered(">", value, self.const(0.0)), self.finite(value)))
        )
        return b.call(self.log, [value])

    def power(self, base, exponent):
        """CPython float_pow with a trace-time constant exponent and a traced base."""
        if isinstance(exponent, bool) or not isinstance(exponent, (int, float)):
            raise TraceRefused("a traced exponent")
        iw = float(exponent)
        if iw == 0.0:  # v**0 is 1.0, even 0**0 and nan**0
            self.counts["pow_constant_one"] += 1
            return self.const(1.0)
        if not math.isfinite(iw):
            raise TraceRefused("a non-finite exponent")
        return self._memo("pow", (base, self.const(iw)), lambda: self._power(base, iw))

    def _power(self, base, iw):
        b = self.builder
        integer = iw == math.floor(iw)
        odd = integer and math.fmod(abs(iw), 2.0) == 1.0  # DOUBLE_IS_ODD_INTEGER
        # non-finite or zero base: the special-case branches of float_pow; fall back
        self.flag(self.not_(self.finite(base)))
        self.flag(b.fcmp_ordered("==", base, self.const(0.0)))
        negative = b.fcmp_ordered("<", base, self.const(0.0))
        if integer:
            absolute = b.select(negative, b.fneg(base), base)
            self.counts["select"] += 1
        else:
            self.flag(negative)  # a complex result in Python
            absolute = base
        self.counts["pow_calls"] += 1
        raw = b.call(self.pow, [absolute, self.const(iw)])
        unit = b.fcmp_ordered("==", absolute, self.const(1.0))
        result = b.select(unit, self.const(1.0), raw)  # 1**w is 1 without calling pow
        self.counts["select"] += 1
        self.flag(self.not_(self.finite(result)))  # _Py_ADJUST_ERANGE1: an infinite result raises
        if odd:
            result = b.select(negative, b.fneg(result), result)
            self.counts["select"] += 1
        return result

    def select_if_equal(self, value, probe, replacement):
        self.counts["select"] += 1
        b = self.builder
        hit = b.fcmp_ordered("==", value.v, self.const(probe))
        return _V(self, b.select(hit, self.const(replacement), value.v))


class _V:
    """A traced double.  Python operator semantics map one to one onto IEEE instructions."""

    __slots__ = ("t", "v")

    def __init__(self, tracer, value):
        self.t, self.v = tracer, value

    def _other(self, other):
        if isinstance(other, _V):
            return other.v
        return self.t.const(other)

    def __add__(self, o):
        return _V(self.t, self.t.binary("fadd", self.v, self._other(o)))

    def __radd__(self, o):
        return _V(self.t, self.t.binary("fadd", self._other(o), self.v))

    def __sub__(self, o):
        return _V(self.t, self.t.binary("fsub", self.v, self._other(o)))

    def __rsub__(self, o):
        return _V(self.t, self.t.binary("fsub", self._other(o), self.v))

    def __mul__(self, o):
        return _V(self.t, self.t.binary("fmul", self.v, self._other(o)))

    def __rmul__(self, o):
        return _V(self.t, self.t.binary("fmul", self._other(o), self.v))

    def __truediv__(self, o):
        return _V(self.t, self.t.binary("fdiv", self.v, self._other(o)))

    def __rtruediv__(self, o):
        return _V(self.t, self.t.binary("fdiv", self._other(o), self.v))

    def __neg__(self):
        return _V(self.t, self.t.neg(self.v))

    def __pow__(self, o):
        if isinstance(o, _V):
            raise TraceRefused("a traced exponent")
        return _V(self.t, self.t.power(self.v, o))

    def __rpow__(self, o):
        raise TraceRefused("a traced exponent")

    def __bool__(self):
        raise TraceRefused("a data-dependent branch")

    def __eq__(self, o):
        raise TraceRefused("a data-dependent comparison")

    __ne__ = __lt__ = __le__ = __gt__ = __ge__ = __eq__
    __hash__ = None

    def __abs__(self):
        raise TraceRefused("abs")

    def __float__(self):
        raise TraceRefused("float() of a traced value")

    __pos__ = __mod__ = __floordiv__ = __abs__


class _MathShim:
    def __init__(self, tracer):
        self._t = tracer

    def exp(self, x):
        return _V(self._t, self._t.exp_(x.v)) if isinstance(x, _V) else math.exp(x)

    def log(self, x):
        return _V(self._t, self._t.log_(x.v)) if isinstance(x, _V) else math.log(x)

    def __getattr__(self, name):
        raise TraceRefused(f"math.{name}")


# ---------------------------------------------------------------------------------------------
# The source rewrite: the water module's own two functions, two declared changes.
# ---------------------------------------------------------------------------------------------
def _function_source(water, name):
    """The function's text from the water module's own source file (never through the module
    attribute, which may be the installed wrapper or a probe's pass-through)."""
    text = inspect.getsource(water)
    found = [
        node
        for node in ast.parse(text).body
        if isinstance(node, ast.FunctionDef) and node.name == name
    ]
    if len(found) != 1:
        raise TraceRefused(f"expected one module-level {name} in the water module")
    return textwrap.dedent(ast.get_source_segment(text, found[0]))


def _rewrite(tree):
    """Replace the single ridge guard ``if dm1 == 0.0: dm1 = <c>`` by a select of the same values."""

    rewritten = 0

    class _Guard(ast.NodeTransformer):
        def visit_If(self, node):
            nonlocal rewritten
            self.generic_visit(node)
            test = node.test
            if (
                isinstance(test, ast.Compare)
                and isinstance(test.left, ast.Name)
                and test.left.id == "dm1"
                and len(test.ops) == 1
                and isinstance(test.ops[0], ast.Eq)
                and isinstance(test.comparators[0], ast.Constant)
                and not node.orelse
                and len(node.body) == 1
                and isinstance(node.body[0], ast.Assign)
                and len(node.body[0].targets) == 1
                and isinstance(node.body[0].targets[0], ast.Name)
                and node.body[0].targets[0].id == "dm1"
                and isinstance(node.body[0].value, ast.Constant)
            ):
                rewritten += 1
                call = ast.Call(
                    func=ast.Name(id="__select_if_equal__", ctx=ast.Load()),
                    args=[
                        ast.Name(id="dm1", ctx=ast.Load()),
                        test.comparators[0],
                        node.body[0].value,
                    ],
                    keywords=[],
                )
                return ast.copy_location(
                    ast.Assign(targets=[ast.Name(id="dm1", ctx=ast.Store())], value=call), node
                )
            return node

    tree = ast.fix_missing_locations(_Guard().visit(tree))
    if rewritten != 1:
        raise TraceRefused(f"expected exactly one ridge guard, found {rewritten}")
    return tree


def _collector(**fields):
    return fields


def build_ir(mutation=None, reuse=True):
    """Trace the water module's residual into an llvmlite IR module.  Returns (ir_text, meta)."""

    if mutation not in MUTATIONS:
        raise ValueError(f"unknown mutation {mutation!r}")
    from llvmlite import ir

    from dtdc_simulator.core2.props import water

    module = ir.Module(name="dtdc_iapws95_compiled")
    tracer = _Tracer(ir, module, reuse=reuse)
    sources = {name: _function_source(water, name) for name in ("_phi0", "_residual_uncached")}
    tree = ast.parse(sources["_phi0"] + "\n\n" + sources["_residual_uncached"])
    tree = _rewrite(tree)
    namespace = {name: getattr(water, name) for name in TRACED_CONSTANTS}
    traced_constants = dict(namespace)  # the module's own objects, before any mutation
    if mutation == "reversed_term_order":
        # the mutation of the test: the same 51 terms of groups 1 and 2 summed in reverse order
        for name in ("_NR", "_C", "_D", "_T"):
            namespace[name] = tuple(reversed(namespace[name]))
    namespace.update(
        math=_MathShim(tracer),
        HelmholtzResult=_collector,
        HelmholtzResult_fields=OUTPUT_FIELDS,
        __select_if_equal__=lambda v, probe, rep: tracer.select_if_equal(v, probe, rep),
        __builtins__={"zip": zip, "range": range, "float": float, "tuple": tuple},
    )
    namespace["__annotations__"] = {}
    exec(compile(tree, "<traced water residual>", "exec"), namespace)  # noqa: S102
    delta = _V(tracer, tracer.fn.args[0])
    tau = _V(tracer, tracer.fn.args[1])
    fields = namespace["_residual_uncached"](delta, tau)
    if fields.get("delta") is not delta or fields.get("tau") is not tau:
        raise TraceRefused("the result does not carry the input delta and tau")
    if set(fields) != {"delta", "tau", *OUTPUT_FIELDS}:
        raise TraceRefused(f"unexpected result fields {sorted(fields)}")
    b = tracer.builder
    out = tracer.fn.args[2]
    for index, name in enumerate(OUTPUT_FIELDS):
        value = fields[name]
        value = value.v if isinstance(value, _V) else tracer.const(value)
        tracer.flag(tracer.not_(tracer.finite(value)))
        b.store(value, b.gep(out, [ir.Constant(ir.IntType(32), index)]))
    b.ret(b.zext(tracer.bad, ir.IntType(32)))
    # the function-pointer readback of the runtime binding
    for name in RUNTIME_FUNCTIONS:
        fn = getattr(tracer, name)
        g = ir.GlobalVariable(module, fn.type, "dtdc_binding_" + name)
        g.initializer = fn
        g.global_constant = True
    # the FMA probe: a * b + c as two separate unflagged instructions
    sig = ir.FunctionType(tracer.double, [tracer.double] * 3)
    probe = ir.Function(module, sig, "dtdc_fma_probe")
    pb = ir.IRBuilder(probe.append_basic_block("entry"))
    pb.ret(pb.fadd(pb.fmul(probe.args[0], probe.args[1]), probe.args[2]))
    text = str(module)
    for forbidden in (
        "fast",
        "contract",
        "reassoc",
        "nnan",
        "ninf",
        "nsz",
        "arcp",
        "afn",
        "fmuladd",
    ):
        if f" {forbidden} " in text or f"llvm.{forbidden}" in text:
            raise TraceRefused(f"the emitted IR carries {forbidden!r}")
    meta = {
        "operations": dict(tracer.counts),
        "exact_reuse": reuse,
        "operations_reused": tracer.reused,
        "source_sha256": {k: _sha(v) for k, v in sources.items()},
        "mutation": mutation,
    }
    meta["_traced_constants"] = traced_constants
    return text, meta


def _sha(text):
    import hashlib

    return hashlib.sha256(text.encode("utf-8")).hexdigest()


# ---------------------------------------------------------------------------------------------
# The runtime binding and the JIT.
# ---------------------------------------------------------------------------------------------
def runtime_addresses():
    """The C runtime's exp, log, pow: ucrtbase.dll's exports and the api-set forwarder the
    interpreter's DLL imports them through must resolve to the same address."""

    ucrt = ctypes.CDLL("ucrtbase")
    forwarder = ctypes.CDLL("api-ms-win-crt-math-l1-1-0")
    out = {}
    for name in RUNTIME_FUNCTIONS:
        a = ctypes.cast(getattr(ucrt, name), ctypes.c_void_p).value
        f = ctypes.cast(getattr(forwarder, name), ctypes.c_void_p).value
        if a != f:
            raise RuntimeError(f"{name}: ucrtbase {a:#x} != api-set forwarder {f:#x}")
        out[name] = a
    return out


class Compiled:
    """One JIT-compiled residual.  Keeps the engine alive for the function's lifetime."""

    def __init__(self, mutation=None, opt=CODEGEN_OPT_LEVEL, reuse=True):
        import llvmlite.binding as llvm

        started = time.perf_counter()
        if os.name != "nt":
            raise RuntimeError(
                "the UCRT binding is Windows-only; see the record's portability section"
            )
        llvm.initialize_native_target()
        llvm.initialize_native_asmprinter()
        addresses = runtime_addresses()
        for name, address in addresses.items():
            llvm.add_symbol(SYMBOL_PREFIX + name, address)
        text, meta = build_ir(mutation, reuse=reuse)
        self.constants_current = _constants_guard(meta.pop("_traced_constants"))
        traced = time.perf_counter()
        target = llvm.Target.from_triple(llvm.get_process_triple())
        machine = target.create_target_machine(
            cpu=TARGET_CPU, features=TARGET_FEATURES, opt=opt, jit=True
        )
        module = llvm.parse_assembly(text)
        module.verify()
        self.assembly = machine.emit_assembly(llvm.parse_assembly(text))
        self.engine = llvm.create_mcjit_compiler(module, machine)
        self.engine.finalize_object()
        self.llvm_module = module
        # the binding read back from the code's own function-pointer globals
        bound = {}
        for name in RUNTIME_FUNCTIONS:
            cell = self.engine.get_global_value_address("dtdc_binding_" + name)
            bound[name] = ctypes.c_void_p.from_address(cell).value
            if bound[name] != addresses[name]:
                raise RuntimeError(
                    f"the generated code binds {name} to {bound[name]:#x}, not ucrtbase {addresses[name]:#x}"
                )
        self._fn = ctypes.CFUNCTYPE(
            ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double)
        )(self.engine.get_function_address(FUNCTION_NAME))
        self._probe = ctypes.CFUNCTYPE(
            ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double
        )(self.engine.get_function_address("dtdc_fma_probe"))
        self._buffer = (ctypes.c_double * len(OUTPUT_FIELDS))()
        self.fma = fma_check(self._probe, self.assembly)
        if not self.fma["unfused"]:
            raise RuntimeError(f"the target machine contracts a*b+c: {self.fma}")
        self.meta = {
            **meta,
            "target_triple": machine.triple,
            "target_cpu": TARGET_CPU,
            "target_features": TARGET_FEATURES,
            "codegen_opt_level": opt,
            "ir_optimisation_passes": "none",
            "runtime_addresses": {k: f"{v:#x}" for k, v in addresses.items()},
            "bound_addresses": {k: f"{v:#x}" for k, v in bound.items()},
            "fma_check": self.fma,
            "ir_sha256": _sha(text),
            "assembly_sha256": _sha(self.assembly),
            "trace_s": traced - started,
            "build_s": time.perf_counter() - started,
            "llvmlite": __import__("llvmlite").__version__,
            "llvm_version": ".".join(map(str, llvm.llvm_version_info)),
        }

    def values(self, delta, tau):
        """The eleven doubles, or None where the Python residual must answer (fallback)."""

        buffer = self._buffer
        if self._fn(delta, tau, buffer):
            return None
        return buffer[:]


def fma_check(probe, assembly):
    """a*b+c with a = 1+2^-30, b = 1-2^-30, c = -1: rounded twice it is 0.0; fused it is -2^-60."""

    a, b, c = 1.0 + 2.0**-30, 1.0 - 2.0**-30, -1.0
    got = probe(a, b, c)
    lowered = assembly.lower()
    fma_mnemonics = [m for m in ("vfmadd", "vfmsub", "vfnmadd", "vfnmsub") if m in lowered]
    return {
        "probe_result_hex": struct.pack(">d", got).hex(),
        "python_result_hex": struct.pack(">d", a * b + c).hex(),
        "fused_result_hex": struct.pack(">d", -(2.0**-60)).hex(),
        "fma_mnemonics_in_assembly": fma_mnemonics,
        "unfused": got == 0.0
        and struct.pack(">d", got) == struct.pack(">d", a * b + c)
        and not fma_mnemonics,
    }


# ---------------------------------------------------------------------------------------------
# The switch: water._residual_uncached replaced in-process, the cache and its keys unchanged.
# ---------------------------------------------------------------------------------------------
_STATE = {"compiled": None, "python": None, "calls": 0, "fallbacks": 0}


def compiled_residual_uncached(delta, tau):
    """Drop-in for water._residual_uncached; falls back to the Python residual on any flag."""

    state = _STATE
    compiled = state["compiled"]
    if (
        type(delta) is float  # noqa: E721 - exact float only
        and type(tau) is float  # noqa: E721
        and compiled.constants_current(state["water_dict"])
    ):
        values = compiled.values(delta, tau)
        state["calls"] += 1
        if values is not None:
            (
                phir,
                phir_d,
                phir_dd,
                phir_t,
                phir_tt,
                phir_dt,
                phir_dtt,
                phir_dttt,
                p0,
                p0_t,
                p0_tt,
            ) = values
            return state["result_type"](
                delta=delta,
                tau=tau,
                phir=phir,
                phir_d=phir_d,
                phir_dd=phir_dd,
                phir_t=phir_t,
                phir_tt=phir_tt,
                phir_dt=phir_dt,
                phir_dtt=phir_dtt,
                phir_dttt=phir_dttt,
                phi0=p0,
                phi0_t=p0_t,
                phi0_tt=p0_tt,
            )
    state["fallbacks"] += 1
    return state["python"](delta, tau)


def _wrapper_chain(fn):
    """The attribute and every function it wraps (``__wrapped__``), outermost first."""
    seen = []
    while fn is not None and fn not in seen and len(seen) < 16:
        seen.append(fn)
        fn = getattr(fn, "__wrapped__", None)
    return seen


def install(compiled=None):
    """Arm the compiled residual in this process (idempotent).  Returns the build metadata."""

    from dtdc_simulator.core2.props import water

    if _STATE["compiled"] is None:
        _STATE["compiled"] = compiled or Compiled()
        _STATE["python"] = water._residual_uncached
        _STATE["result_type"] = water.HelmholtzResult
        _STATE["water_dict"] = water.__dict__
    if compiled_residual_uncached not in _wrapper_chain(water._residual_uncached):
        if water._residual_uncached is not _STATE["python"]:
            raise RuntimeError("water._residual_uncached was replaced by someone else")
        # inspect.getsource(water._residual_uncached) must keep returning the module's own source:
        # native_water_path_bounds.interval_namespace re-executes it on Interval operands
        compiled_residual_uncached.__wrapped__ = _STATE["python"]
        water._residual_uncached = compiled_residual_uncached
    return _STATE["compiled"].meta


def install_if_armed(environ=None):
    """Kernel-envelope IAPWS_COMPILED (armed by exactly "1"): install for the rest of the process
    and return its stamp; None when the knob is unset (nothing installed, nothing imported)."""

    import particle_kernel_envelope as envelope

    if not envelope.armed(envelope.IAPWS_COMPILED, environ):
        return None
    install()
    return {**stamp(), "build": _STATE["compiled"].meta}


def uninstall():
    from dtdc_simulator.core2.props import water

    if _STATE["python"] is not None and water._residual_uncached is compiled_residual_uncached:
        water._residual_uncached = _STATE["python"]


def installed():
    from dtdc_simulator.core2.props import water

    return compiled_residual_uncached in _wrapper_chain(water._residual_uncached)


def stamp():
    meta = _STATE["compiled"].meta if _STATE["compiled"] is not None else None
    return {
        "knob": "DTDC_KERNEL_ENVELOPE_IAPWS_COMPILED",
        "installed": installed(),
        "calls": _STATE["calls"],
        "fallbacks": _STATE["fallbacks"],
        "ir_sha256": None if meta is None else meta["ir_sha256"],
        "record": RECORD,
    }


__all__ = [
    "Compiled",
    "TraceRefused",
    "build_ir",
    "compiled_residual_uncached",
    "fma_check",
    "install",
    "installed",
    "runtime_addresses",
    "stamp",
    "uninstall",
]
