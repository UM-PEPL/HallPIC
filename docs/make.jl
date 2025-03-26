using Documenter
using HallPIC

makedocs(
    sitename = "HallPIC",
    pages = [
        "Design Overview" => "index.md",
        "Reactions Overview" => "reactions.md"
    ],
    format = Documenter.HTML(),
    modules = [HallPIC]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/UM-PEPL/HallPIC.git"
)
