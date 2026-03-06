#!/bin/bash

# Default paths (allow override via environment)
DEFAULT_CONFIG="${SPECIATOR_CONFIG:-/speciator/config.toml}"
DEFAULT_LIBRARY="${SPECIATOR_LIBRARY:-/speciator/library}"

# Check if config-file option is present (handles both forms)
if [[ ! " $* " =~ " --config-file " ]] && [[ ! " $* " =~ " -c " ]]; then
    set -- --config-file "$DEFAULT_CONFIG" "$@"
fi

# Check if library-location option is present (handles both forms)
if [[ ! " $* " =~ " --library-location " ]] && [[ ! " $* " =~ " -l " ]]; then
    set -- --library-location "$DEFAULT_LIBRARY" "$@"
fi

# Execute speciator with the processed arguments
exec speciator "$@"
