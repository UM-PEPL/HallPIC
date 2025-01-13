using Documenter
using HallPIC

makedocs(
    sitename = "HallPIC",
    format = Documenter.HTML(),
    modules = [HallPIC]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/UM-PEPL/HallPIC.git"
)
