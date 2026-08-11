# Deletes the freeze so every document re-executes under the `full` profile.
# Refuses when the freeze has uncommitted changes that git could not restore.
# Set FORCE_RESET_FREEZE=1 to override.

const ROOT = normpath(joinpath(@__DIR__, ".."))
const TARGETS = ["_freeze", "slides/_freeze"]

cd(ROOT)

function tracked_changes()
    out = try
        read(`git status --porcelain -- $TARGETS`, String)
    catch
        error("cannot run git to verify the freeze is safe to delete; " *
              "set FORCE_RESET_FREEZE=1 to delete anyway")
    end
    # Untracked entries ("?? ") are regenerated output, not work to protect.
    filter(line -> !isempty(line) && !startswith(line, "??"), split(out, '\n'))
end

if get(ENV, "FORCE_RESET_FREEZE", "") != "1"
    changes = tracked_changes()
    if !isempty(changes)
        error("""
              refusing to delete the freeze: $(length(changes)) uncommitted change(s).
              $(join(changes, '\n'))
              Commit them, or set FORCE_RESET_FREEZE=1 to discard them.
              """)
    end
end

for target in TARGETS
    rm(joinpath(ROOT, target); recursive = true, force = true)
end
println("freeze reset: ", join(TARGETS, ", "))
