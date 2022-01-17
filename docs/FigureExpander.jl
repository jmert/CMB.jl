module FigureExpander

import Markdown
import Documenter: Expanders, Documents
import Documenter.Utilities: Utilities, Selectors

struct FigureNode
    img::Markdown.Image
    params::Dict{Symbol,Any}
    caption::Vector{Any}
end

"""
Expands to a figure block, appropriate for both HTML and LaTeX output.

## Example

````md
```@figure; htmllabel = "Figure 1", latexlabel = "fig:key_name"
![alt text for accesibility](./source_file.ext)

Any _Markdown_ content can then be written as the rest of the "code" block,
which will be passed through to the `<figcaption>` / `\\caption` blocks in
HTML / ``\\LaTeX``, respectively.  The preceding Markdown image is parsed and
used to define the `<img>` / `\\includegraphics` elements, so that the image and
caption are included in the document as they should be when written natively.

Following the opening ` ```@figure`, optional parameters may be included
as a key-value list comma-separated values. The only two recognized parameters
at this time are:

- htmllabel

  A name which is prepended to the caption included in HTML output. Intended
  to give explicit text for use in referring to the figure.

- latexlabel

  The label name to emit in a `\\label` command after the caption, for use
  in references elsewhere in the LaTeX document.
```
````
"""
abstract type FigureBlock <: Expanders.ExpanderPipeline end

Selectors.order(::Type{FigureBlock}) = 50.0
Selectors.matcher(::Type{FigureBlock}, node, page, doc) = Expanders.iscode(node, r"^@figure")

let
    param_regex = r"""
        (?<name>\w+) # word-like name
        \s*=\s*
        (?:
            # a quoted string expression (with escapes) --- https://stackoverflow.com/a/10786066
            "(?<strval>[^"\\]*(?:\\.[^"\\]*)*)"
        |
            # a simple literal expression -- no whitespace, not starting with "
            (?<litval>[^"][^\s]*)
        )
        \s*,?\s*  # eat comma separator and spaces, as necessary
        """x

    function Selectors.runner(::Type{FigureBlock}, x, page, doc)
        matched = match(r"^@figure\s*(?<params>;.*)?$", x.language)
        matched === nothing && error("invalid '@figure' syntax: $(x.language))")
        content = Utilities.mdparse(x.code, mode = :blocks)::Vector{Any}
        isempty(content) && error("invalid `@figure` block: block is empty")

        # expect to find image as the first element in the block
        para1 = content[1]::Markdown.Paragraph
        para1content = para1.content::Vector{Any}
        if !(para1content[1] isa Markdown.Image)
            error("invalid `@figure` block: expected first contents to be image, none found")
        end
        # remove the image node from the rest of the caption
        img = popfirst!(para1content)::Markdown.Image
        # avoid an empty paragraph if the image node was the only thing in the paragraph
        if isempty(para1content)
            popfirst!(content)
        end

        # parse parameters and store
        if matched[:params] !== nothing
            params = Dict{Symbol,String}()
            paramstr = strip(lstrip(isequal(';'), matched[:params]))
            idx = firstindex(paramstr)
            while true
                m = match(param_regex, paramstr, idx)
                m === nothing && break
                params[Symbol(m[:name])] = something(m[:strval], m[:litval])
                idx = nextind(paramstr, idx, length(m.match))
            end
            if idx <= lastindex(paramstr)
                error("invalid `@figure` syntax: trailing string in `$paramstr`")
            end
        end
        page.mapping[x] = FigureNode(img, params, content)
    end
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
    caption = node.caption
    if haskey(node.params, :htmllabel)
        pushfirst!(caption[1].content, node.params[:htmllabel] * " ")
    end
    return figure(
            HTMLWriter.domify(ctx, navnode, node.img),
            figcaption(HTMLWriter.domify(ctx, navnode, caption))
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
        if haskey(node.params, :latexlabel)
            wrapinline(io, "label") do
                _print(io, node.params[:latexlabel])
            end
        end
        _println(io)
    end
end

end # FigureExpander
