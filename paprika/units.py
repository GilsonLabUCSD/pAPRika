from simtk import unit
import logging

logger = logging.getLogger(__name__)


def _ast_eval(node):
    """
    Performs an algebraic syntax tree evaluation of a unit.
    Parameters
    ----------
    node : An ast parsing tree node
    """
    import ast
    import operator as op

    operators = {ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul,
        ast.Div: op.truediv, ast.Pow: op.pow, ast.BitXor: op.xor,
        ast.USub: op.neg}

    if isinstance(node, ast.Num): # <number>
        return node.n
    elif isinstance(node, ast.BinOp): # <left> <operator> <right>
        return operators[type(node.op)](_ast_eval(node.left), _ast_eval(node.right))
    elif isinstance(node, ast.UnaryOp): # <operator> <operand> e.g., -1
        return operators[type(node.op)](_ast_eval(node.operand))
    elif isinstance(node, ast.Name):
        # see if this is a simtk unit
        b = getattr(unit, node.id)
        return b
    # TODO: This was a quick hack that surprisingly worked. We should validate this further.
    elif isinstance(node, ast.List):
        return ast.literal_eval(node)
    else:
        raise TypeError(node)

def string_to_quantity(quantity_string):
    """
    Takes a string representation of a quantity and returns a simtk.unit.Quantity
    Parameters
    ----------
    quantity_string : str
        The quantity to deserialize
    Returns
    -------
    output_quantity : simtk.unit.Quantity
        The deserialized quantity
    """
    if quantity_string is None:
        return None
    # This can be the exact same as string_to_unit
    import ast
    output_quantity = _ast_eval(ast.parse(quantity_string, mode='eval').body)
    return output_quantity
