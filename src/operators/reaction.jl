# Intentionally minimal for now.
# Concrete reaction implementations should subtype AbstractReactionOperator
# and implement one or more of:
#   rhs!
#   implicit_rhs!
#   step!
#   residual!
#   jacobian!
