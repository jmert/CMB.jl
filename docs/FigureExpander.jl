module FigureExpander

import Markdown
import Documenter: Expanders, Documents
import Documenter.Utilities: Utilities, Selectors

struct FigureNode
    img::Markdown.Image
    caption::Vector{Any}
end

abstract type FigureBlock <: Expanders.ExpanderPipeline end

Selectors.order(::Type{FigureBlock}) = 50.0
Selectors.matcher(::Type{FigureBlock}, node, page, doc) = Expanders.iscode(node, r"^@figure")

function Selectors.runner(::Type{FigureBlock}, x, page, doc)
    matched = match(r"^@figure\s+\"(?<filename>.+?)(?<!\\)\"\s*(?<params>;.*)?$", x.language)
    matched === nothing && error("invalid '@figure' syntax: $(x.language))")
    filename, params = matched.captures
    alt = ""
    if params !== nothing
        if occursin(r"\balt\b", params)
            alt = match(r"\balt\s*=\s*\"(.+?)(?<!\\)\"", params).captures[1]
        end
    end
    page.mapping[x] = FigureNode(Markdown.Image(filename, alt),
                                 Utilities.mdparse(x.code, mode=:blocks))
end

Documents.walk(f, meta, block::FigureNode) = f(block) ? begin
        Documents.walk(f, meta, block.img)
        Documents.walk(f, meta, block.caption)
    end : nothing

#########################################################################################

import Documenter.Writers: HTMLWriter
import Documenter.Utilities.DOM: @tags

function HTMLWriter.domify(ctx, navnode, node::FigureNode)
    @tags figure figcaption
    return figure(
            HTMLWriter.domify(ctx, navnode, node.img),
            figcaption(HTMLWriter.domify(ctx, navnode, node.caption))
    )
end

#########################################################################################

import Documenter: Utilities
import Documenter.Writers: LaTeXWriter
import .LaTeXWriter: _print, _println, wrapinline, wrapblock

# near-copy of LaTeXWriter.latexinline(io::IO, md::Markdown.Image)
function LaTeXWriter.latex(io::IO, node::FigureNode)
    img = node.img
    caption = node.caption
    wrapblock(io, "figure") do
        _println(io, "\\centering")
        url = if Utilities.isabsurl(img.url)
            @warn "images with absolute URLs not supported in LaTeX output in $(Utilities.locrepr(io.filename))" url = img.url
            # We nevertheless output an \includegraphics with the URL. The LaTeX build will
            # then give an error, indicating to the user that something wrong. Only the
            # warning would be drowned by all the output from LaTeX.
            img.url
        elseif startswith(img.url, '/')
            # URLs starting with a / are assumed to be relative to the document's root
            normpath(lstrip(img.url, '/'))
        else
            normpath(joinpath(dirname(io.filename), img.url))
        end
        url = replace(url, "\\" => "/") # use / on Windows too.
        url, _ = splitext(url) # let LaTeX figure out the file extension
        wrapinline(io, "includegraphics[max width=\\linewidth]") do
            _print(io, url)
        end
        _println(io)
        wrapinline(io, "caption[$(img.alt)]") do
            LaTeXWriter.latex(io, caption)
        end
        _println(io)
    end
end

end # FigureExpander
