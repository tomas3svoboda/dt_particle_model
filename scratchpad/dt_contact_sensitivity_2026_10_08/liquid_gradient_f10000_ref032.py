"""Algebraically identical Luikov log-activity difference without exp/ratio loss."""
import ast
import hashlib
import inspect
import math


def stable_log_activity_gradient(w_left, p_left, w_right, p_right, distance):
    if not all(math.isfinite(w) and w > 0 for w in (w_left, w_right)):
        raise ValueError("liquid log-activity gradient requires positive finite water")
    if distance <= 0 or not math.isfinite(distance):
        raise ValueError("liquid activity distance must be positive and finite")
    if p_left.A1 != p_right.A1 or p_left.A2 != p_right.A2:
        raise ValueError("stable difference requires the same Luikov law on both sides")
    # ln(a_right)-ln(a_left) = A1/A2 * (W_right-W_left)/(W_left*W_right).
    # Native activity evaluation in the caller still enforces both domains.
    return (p_left.A1/p_left.A2) * ((w_right-w_left)/w_left/w_right) / distance


def installed_function(module):
    """Replace exactly the liquid-gradient expression, preserving all other AST."""
    original = module.liquid_connected_fast_path
    source_function = inspect.unwrap(original)
    source = inspect.getsource(source_function)
    old = "gradient = math.log(aw_right / aw_left) / args[4]"
    new = ("gradient = _stable_log_activity_gradient(left.retained_water_loading, args[6].luikov, "
           "right.retained_water_loading, args[7].luikov, args[4])")
    if source.count(old) != 1:
        raise ValueError("source drift: expected exactly one liquid gradient expression")
    rewritten = source.replace(old, new)
    old_ast, new_ast = ast.parse(source), ast.parse(rewritten)
    class Restore(ast.NodeTransformer):
        changes = 0
        def visit_Assign(self, node):
            if isinstance(node.targets[0], ast.Name) and node.targets[0].id == "gradient":
                self.changes += 1
                return ast.parse(old).body[0]
            return self.generic_visit(node)
    restore = Restore()
    restored = restore.visit(new_ast)
    assert restore.changes == 1
    assert ast.dump(restored, include_attributes=False) == ast.dump(old_ast, include_attributes=False)
    namespace = dict(source_function.__globals__)
    namespace["_stable_log_activity_gradient"] = stable_log_activity_gradient
    exec(compile(rewritten, inspect.getsourcefile(source_function), "exec"), namespace)
    replacement = namespace[original.__name__]
    return replacement, {"original_function_sha256": hashlib.sha256(source.encode()).hexdigest(),
                         "rewritten_function_sha256": hashlib.sha256(rewritten.encode()).hexdigest(),
                         "only_liquid_gradient_AST_expression_changed": True,
                         "force": "(A1/A2)*(W_right-W_left)/(W_left*W_right)/distance",
                         "physical_law_changed": False}
