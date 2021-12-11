using Documenter, Logging

# Outputs are Julia version specific
if VERSION >= v"1.6"
    # Disable Documeter's Info logging
    oldlvl = Logging.min_enabled_level(current_logger())
    disable_logging(Logging.Info)
    try
        DocMeta.setdocmeta!(CMB, :DocTestSetup, :(using CMB); recursive=true)
        doctest(CMB, testset="Doc Tests")
    finally
        disable_logging(oldlvl - 1)
    end
end
