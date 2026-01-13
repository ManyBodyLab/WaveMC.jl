using Literate: Literate
using WaveMC

Literate.markdown(
    joinpath(pkgdir(WaveMC), "docs", "files", "README.jl"),
    joinpath(pkgdir(WaveMC), "docs", "src");
    flavor = Literate.DocumenterFlavor(),
    name = "index",
)
