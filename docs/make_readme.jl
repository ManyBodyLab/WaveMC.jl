using Literate: Literate
using WaveMC

Literate.markdown(
    joinpath(pkgdir(WaveMC), "docs", "files", "README.jl"),
    joinpath(pkgdir(WaveMC));
    flavor = Literate.CommonMarkFlavor(),
    name = "README",
)
