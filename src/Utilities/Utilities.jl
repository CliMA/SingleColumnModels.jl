module Utilities

export @unpack_fields

"""
    unpack_fields(_struct, syms...)

Unpack struct properties `syms`
from struct `_struct`

# Example
```julia
julia> struct Foo;a;b;c;end

julia> f = Foo(1,2,3)
Foo(1, 2, 3)

julia> @unpack_fields f a c; @show a c
a = 1
c = 3
```
"""
macro unpack_fields(_struct, syms...)
    thunk = Expr(:block)
    for sym in syms
        push!(
            thunk.args,
            :($(esc(sym)) = getproperty($(esc(_struct)), $(QuoteNode(sym)))),
        )
    end
    push!(thunk.args, nothing)
    return thunk
end

end