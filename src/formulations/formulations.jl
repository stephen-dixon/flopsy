struct UnsplitFormulation <: AbstractFormulation end
struct IMEXFormulation <: AbstractFormulation end

struct LieSplit end
struct StrangSplit end

struct SplitFormulation{S} <: AbstractFormulation
    scheme::S
end

struct ResidualFormulation <: AbstractFormulation end
