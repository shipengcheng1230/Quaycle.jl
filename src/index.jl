## some conventions

_ϵσindex(::Val{:xx}) = 1
_ϵσindex(::Val{:xy}) = 2
_ϵσindex(::Val{:xz}) = 3
_ϵσindex(::Val{:yy}) = 4
_ϵσindex(::Val{:yz}) = 5
_ϵσindex(::Val{:zz}) = 6

const _diagcomponent = (:xx, :yy, :zz)
